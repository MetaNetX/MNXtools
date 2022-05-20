package Compartments;
#Filename Compartments.pm

use strict;
use warnings;
use diagnostics;

use List::Util qw( uniq );
use FindBin qw( $RealBin );
use lib "$RealBin/../perl/SBML";
use Constants;
use Formaters;
use Notes;
use Annotation;

our $compartments;


sub parse_Compartments_Info {
    my ($SBML_model) = @_;

    for my $comp ( $SBML_model->getListOfCompartments() ){
        my $Sid = $comp->getId();
        $compartments->{ $Sid }->{'source'}  = $Sid;
        $compartments->{ $Sid }->{'name'}    = Formaters::trim( $comp->getName() // '' );
        $compartments->{ $Sid }->{'Sid'}     = Formaters::format_Compartment($Sid);
        $compartments->{ $Sid }->{'xrefs'}   = [];
        $compartments->{ $Sid }->{'sboTerm'} = $comp->getSBOTermID();

        # Parse notes
        if ( $comp->isSetNotes() ){
            my $notes;
            $notes = Notes::parse_notes($comp, $notes);
            for my $note (keys %$notes){
                $compartments->{ $Sid }->{lc $note} = join(';', @{ $notes->{$note} });
            }
        }

        # Parse annotations
        #NOTE annotations overload notes!
        if ( $comp->isSetAnnotation() ){
            my $xrefs;
            $xrefs = Annotation::parse_annotation($comp, $xrefs, 'comp');
            push @{ $compartments->{ $Sid }->{'xrefs'} }, map { Annotation::clean_comp_annotation($_) }
                                                          map { @{ $xrefs->{$_} } }
                                                          grep { /_comp$/ }
                                                          keys %$xrefs;
        }
    }

    return;
}

sub create_compartment {
    my ($SBML_model, $compartment_Sid, $compartment_name) = @_;

    # Add boundary compartment if not already there
    if ( ! $SBML_model->getCompartment($compartment_Sid) ){
        my $compartment_obj = $SBML_model->createCompartment();
        $compartment_obj->setId($compartment_Sid);
        $compartment_obj->setName($compartment_name);
        $SBML_model->addCompartment($compartment_obj);
    }

    return;
}

# For example in case of a compartment all as boundaryCondition=true
# So it is nerver used, replaced by BOUNDARY
sub remove_compartment {
    my ($SBML_model, $compartment_Sid) = @_;

    if ( $SBML_model->getCompartment($compartment_Sid) ){
        $SBML_model->removeCompartment($compartment_Sid);
    }

    return;
}


sub create_SBML_compartment {
    my ($MetNet, $mnet_id, $SBML_model, $use_notes, $comp_id) = @_;

    my $comp_id_fixed = $comp_id;
    $comp_id_fixed =~ s{^UNK:}{};
    my $comp = $SBML_model->createCompartment();
    $comp->setId($comp_id_fixed);
    $comp->setConstant(1); # True
    my ($comp_name, $comp_xrefs) = $MetNet->get_comp_info($comp_id);
    $comp->setName($comp_name);
    $comp->setSBOTerm( $comp_id eq $Constants::boundary_comp_id ? $Constants::boundary_comp_sbo : $Constants::default_comp_sbo );
    #notes
    if ( $use_notes ){
        my $notes = '';
        if ( $MetNet->get_comp_source($mnet_id, $comp_id) ){
            $notes .= '<html:p>SOURCE: '.$MetNet->get_comp_source($mnet_id, $comp_id).'</html:p>';
        }
        if ( $comp_xrefs ){
            $notes .= '<html:p>REFERENCE: '.$comp_xrefs.'</html:p>';
        }
        $comp->setNotes($notes)  if ( $notes );
    }
    #TODO use the prefix from the prefix file
    #TODO annotations: fix is/isRelatedTo and different identifiers.org... when we will have all mapped xrefs
    #TODO add other possible xref sources (e.g., CCO, ...)
    if ( $comp_xrefs ){
        if ( !$comp->isSetMetaId() ){
            $comp->setMetaId( $comp->getId() );
        }
        my $CV = new LibSBML::CVTerm();
        $CV->setQualifierType($LibSBML::BIOLOGICAL_QUALIFIER);
        $CV->setBiologicalQualifierType($LibSBML::BQB_IS);
        if ( $comp_xrefs =~ /^go:\d+$/i ){
            $CV->addResource($Constants::identifiers_go.uc($comp_xrefs));
            $comp->addCVTerm($CV);
        }
        elsif ( $comp_xrefs =~ /^bigg:(..?)$/i ){
            $CV->addResource($Constants::identifiers_biggc.$1);
            $comp->addCVTerm($CV);
        }
        elsif ( $comp_xrefs =~ /^mnx:(.+)$/i ){
            $CV->addResource($Constants::identifiers_mnxc.$1);
            $comp->addCVTerm($CV);
        }
    }

    return;
}

1;

