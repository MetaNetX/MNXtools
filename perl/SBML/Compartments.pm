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

1;

