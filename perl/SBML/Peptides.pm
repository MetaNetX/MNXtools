package Peptides;
#Filename Peptides.pm

use strict;
use warnings;
use diagnostics;

use List::Util qw( uniq );
use FindBin qw( $RealBin );
use lib "$RealBin/../perl/SBML";
use Constants;
use Formaters;
use Annotation;


our $peptides;

sub parse_GeneProducts_Info {
    my ($SBML_model, $SBML_object) = @_;

    #Gene products info available in  SBML 3+ / fbc  only!
    if ( $SBML_object->getLevel() >= 3 ){
        my $fbcModel = $SBML_model->getPlugin('fbc');
        if ( $fbcModel ){
            for ( my $i=0; $i < $fbcModel->getListOfGeneProducts->getNumGeneProducts(); $i++ ){
                my $geneProduct = $fbcModel->getListOfGeneProducts->get($i);
                my $Sid       = $geneProduct->getId();
                my $clean_Sid = Formaters::format_GeneProduct($Sid);

                $peptides->{ $clean_Sid }->{'Sid'}               = $clean_Sid;
                $peptides->{ $clean_Sid }->{'source'}            = $Sid;
                $peptides->{ $clean_Sid }->{'name'}              = Formaters::trim( $geneProduct->getName()              // '' );
                $peptides->{ $clean_Sid }->{'label'}             = Formaters::trim( $geneProduct->getLabel()             // '' );
                $peptides->{ $clean_Sid }->{'associatedSpecies'} = Formaters::trim( $geneProduct->getAssociatedSpecies() // '' );
                $peptides->{ $clean_Sid }->{'xrefs'}             = [];

                # Parse annotations
                if ( $geneProduct->isSetAnnotation() ){
                    my $xrefs;
                    $xrefs = Annotation::parse_annotation($geneProduct, $xrefs, 'pep');
                    push @{ $peptides->{ $clean_Sid }->{'xrefs'} }, map { Annotation::clean_pep_annotation($_) }
                                                                    map { @{ $xrefs->{$_} } }
                                                                    grep { /_pep$/ }
                                                                    keys %$xrefs;
                }
            }
        }
    }

    #NOTE In SBML2 no info about gene, so no need to fill peptides.tsv

    return;
}

sub update_GeneProducts {
    my ($pept_id, $mnet) = @_;

    my $gpr_id = $pept_id->getId();
    my ($gpr_desc, $gpr_xrefs, $gpr_gene) = $mnet->get_pept_info($gpr_id);
    if ( $gpr_desc ){
        $pept_id->setName($gpr_desc);
    }
    if ( $gpr_xrefs || $gpr_gene ){
        my @labels;
        push @labels, $gpr_xrefs  if ( $gpr_xrefs );
        push @labels, $gpr_gene   if ( $gpr_gene  );
        my $gpr_label = join(' -- ', @labels);
        $pept_id->setLabel($gpr_label);
    }

    #NOTE skip notes, this is SBML level 3

    # Add annotations
#TODO add missing pept prefixes
    #sort xrefs, to have uniprot first to set it as BQB_IS
    my @gpr_xrefs = sort {
                            ($b =~ /^uniprotkb:/) <=> ($a =~ /^uniprotkb:/)
                         || ($b =~ /^uniprot:/)   <=> ($a =~ /^uniprot:/)
                         || $a cmp $b
                         }
                    split(';', $gpr_xrefs);
    if ( $gpr_xrefs ){
        if ( !$pept_id->isSetMetaId() ){
            $pept_id->setMetaId( $gpr_id );
        }
        my $CV = new LibSBML::CVTerm();
        $CV->setQualifierType($LibSBML::BIOLOGICAL_QUALIFIER);
        #is(_a) -> reference
        $CV->setBiologicalQualifierType($LibSBML::BQB_IS);
        my $annotation_ref = Formaters::guess_annotation_link('pept', $gpr_xrefs[0]);
        if ( $annotation_ref ){
            $CV->addResource($annotation_ref);
            $pept_id->addCVTerm($CV);
        }
        #isRelatedTo -> other xrefs
        if ( exists $gpr_xrefs[1] ){
            my $CV2 = new LibSBML::CVTerm();
            $CV2->setQualifierType($LibSBML::BIOLOGICAL_QUALIFIER);
            $CV2->setBiologicalQualifierType($LibSBML::BQB_IS_HOMOLOG_TO); #BQB_IS_RELATED_TO looks better but deprecated now!
            #NOTE SBML looks to keep only uniq annotations, so remove identical URIs
            for my $xref ( grep { !/:G_/ } splice(@gpr_xrefs, 1) ){
                my $annotation = Formaters::guess_annotation_link('pept', $xref);
                if ( $annotation ){
                    $CV2->addResource($annotation);
                }
            }
            $pept_id->addCVTerm($CV2);
        }
    }

    return;
}

1;

