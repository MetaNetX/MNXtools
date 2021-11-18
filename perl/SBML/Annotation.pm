package Annotation;
#Filename Annotation.pm

use strict;
use warnings;
use diagnostics;

use FindBin qw( $RealBin );
use lib "$RealBin/../perl/SBML";
use Formaters;
use Constants;

my $id_res = qr/<rdf:li rdf:resource="https?:\/\/identifiers\.org/;

sub parse_annotation {
    my ($element, $element_info, $category) = @_;
    die "No category defined\n"  if ( !$category );

    #NOTE Assume one rdf resource per line!
    #TODO Try to get the qualifier attribute here
    for my $annotation ( grep { /rdf:resource=/ } split(/\n/, $element->getAnnotationString()) ){
        if ( $annotation =~ m|$id_res/([^\/]+?):(\S+?)\s*/?"| ){
            push @{ $element_info->{"$1_$category"} }, "$1:$2";
        }
        elsif ( $annotation =~ m|$id_res/(.+?)/(\S+?)\s*/?"| ){
            push @{ $element_info->{"$1_$category"} }, "$1:$2";
        }
        else {
            $Constants::logger->debug("annotMISS: $annotation");
        }
    }

    return $element_info;
}


sub clean_comp_annotation {
    my ($annotation) = @_;
    my ($prefix, $id) = split(/:/, $annotation, 2);

    if ( $prefix eq 'go' ){
        $annotation = 'GO:'.$id;
    }
    elsif ( $prefix eq 'bigg.compartment' ){
        $annotation = 'biggC:'.$id;
    }
    elsif ( $prefix eq 'metanetx.compartment' ) {
        $annotation = 'mnx:'.$id;
    }
    elsif ( $prefix eq 'cl' ){
        $annotation = 'CL:'.$id;
    }

    return $annotation;
}

sub clean_chem_annotation {
    my ($annotation) = @_;
    my ($prefix, $id) = split(/:/, $annotation, 2);

    if ( $prefix eq 'bigg.metabolite' ){
        $annotation = 'biggM:'.$id;
    }
    elsif ( $prefix eq 'biocyc' ){
        $id =~ s{^META:}{};
        $annotation = 'metacycM:'.$id;
    }
    elsif ( $prefix eq 'metacyc.compound' ){
        $annotation = 'metacycM:'.$id;
    }
    elsif ( $prefix eq 'chebi' ){
        $id =~ s{^CHEBI:}{};
        $annotation = 'CHEBI:'.$id;
    }
    elsif ( $prefix eq 'kegg.compound' ){
        $annotation = 'keggC:'.$id;
    }
    elsif ( $prefix eq 'kegg.drug' ){
        $annotation = 'keggD:'.$id;
    }
    elsif ( $prefix eq 'kegg.glycan' ){
        $annotation = 'keggG:'.$id;
    }
    elsif ( $prefix eq 'kegg.environ' ){
        $annotation = 'keggE:'.$id;
    }
    elsif ( $prefix eq 'metanetx.chemical' ){
        $annotation = 'mnx:'.$id;
    }
    elsif ( $prefix eq 'reactome.compound' || $prefix eq 'reactome' ){
        $annotation = 'reactomeM:'.$id;
    }
    elsif ( $prefix eq 'sabiork.compound' ){
        $annotation = 'sabiorkM:'.$id;
    }
    elsif ( $prefix eq 'seed.compound' ){
        $annotation = 'seedM:'.$id;
    }
    elsif ( $prefix eq 'swisslipid' ){
        $annotation = lc $id;
    }
    elsif ( $prefix eq 'slm' ){
        $annotation = 'SLM:'.$id;
    }
    elsif ( $prefix eq 'lipidmaps' ){
        $annotation = 'lipidmapsM:'.$id;
    }
    elsif ( $prefix eq 'envipath' ){
        $annotation = 'envipathM:'.$id;
    }
    elsif ( $prefix eq 'umbbd.compound' ){
        $annotation = 'umbbd:'.$id;
    }
    elsif ( $prefix eq 'unipathway.compound' ){
        $annotation = 'upa:'.$id;
    }
    # Chem structure
    elsif ( $prefix eq 'inchi' ){
        $annotation = $id;
    }
    elsif ( $prefix eq 'inchikey' || $prefix eq 'inchi_key' ){
        $annotation = $id;
    }
    # Extra
    elsif ( $prefix =~ /^pubchem\./ ){
        $annotation = 'pubchem:'.$id;
    }

    #Nothing to change for hmdb
    #Nothing to change for some extra also
    return $annotation;
}

sub clean_reac_annotation {
    my ($annotation) = @_;
    my ($prefix, $id) = split(/:/, $annotation, 2);

    if ( $prefix eq 'bigg.reaction' ){
        $annotation = 'biggR:'.$id;
    }
    elsif ( $prefix eq 'rhea' ){
        $annotation = 'rheaR:'.$id;
    }
    elsif ( $prefix eq 'kegg.reaction' ){
        $annotation = 'keggR:'.$id;
    }
    elsif ( $prefix eq 'biocyc' ){
        $id =~ s{^META:}{};
        $annotation = 'metacycR:'.$id;
    }
    elsif ( $prefix eq 'metacyc.reaction' ){
        $annotation = 'metacycR:'.$id;
    }
    elsif ( $prefix eq 'metanetx.reaction' ){
        $annotation = 'mnx:'.$id;
    }
    elsif ( $prefix eq 'reactome.reaction' || $prefix eq 'reactome' ){
        $annotation = 'reactomeR:'.$id;
    }
    elsif ( $prefix eq 'sabiork.reaction' ){
        $annotation = 'sabiorkR:'.$id;
    }
    elsif ( $prefix eq 'seed.reaction' ){
        $annotation = 'seedR:'.$id;
    }
    elsif ( $prefix eq 'unipathway.reaction' ){
        $annotation = 'upa:'.$id;
    }
    #EC
    elsif ( $prefix eq 'ec-code' ){
        $annotation = 'EC:'.$id;
    }

    #Nothing to change for some extra also
    return $annotation;
}

sub clean_pep_annotation {
    my ($annotation) = @_;
    my ($prefix, $id) = split(/:/, $annotation, 2);

    if ( $prefix eq 'ncbigi' ){
        $id =~ s{^GI:}{};
        $annotation = 'ncbigi:'.$id;
    }
    elsif ( $prefix =~ /^Ensembl/i ){
        $annotation = 'ensembl:'.$id;
    }
    elsif ( $prefix eq 'biocyc' && $id =~ /^ECOLI:(.+)$/ ){
        $annotation = 'ecocyc_gene:'.$1;
    }

    #Nothing to change for uniprot
    #Nothing to change for some extra also
    return $annotation;
}

1;

