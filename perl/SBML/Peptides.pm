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

1;

