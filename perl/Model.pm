package Model;
#Filename Model.pm

use strict;
use warnings;
use diagnostics;

use Data::Dumper;

use File::Basename;
use List::Util qw( uniq );
use FindBin qw( $RealBin );
use lib "$RealBin/../perl";
use Constants;
use Formaters;
use Notes;
use Annotation;
use Reactions;


# Initialize fields
my $ID               = ''; # the identifier for the model (mandatory)
my $Description      = ''; # a short description for the model (model name usually)
#my $MNXref_Version   = ''; # MNXRef version used at mapping time #NOTE Not defined for unmapped
my $Source_ID        = ''; # the original name of the model, before its mapping to the MNXRef namespace
my $Processed_Date   = ''; # the date at which the model has been converted
my $Taxid            = ''; # the NCBI taxonomic ID of the organism related to this model
my $LB               = -1000; # Global (min) Lower Bound
my $UB               =  1000; # Global (max) Upper Bound
# Extra
my $extra; # e.g. $Source_Reference = ''; # model reference. E.g. PMID  |  my $Source_Link      = ''; # model external source


sub get_Processed_Date {
    return  sprintf('%d/%02d/%02d', ((localtime())[5] + 1900), ((localtime())[4] + 1), ((localtime())[3]));
}


sub set_ID {
    my ($new_ID) = @_;
    $ID = Formaters::convert_to_Sid($new_ID);
    return;
}
sub get_ID {
    return $ID;
}

sub set_Description {
    my ($new_Description) = @_;
    $Description = $new_Description;
    return;
}
sub get_Description {
    return $Description;
}

sub set_Source_ID {
    my ($new_Source_ID) = @_;
    $Source_ID = Formaters::convert_to_Sid($new_Source_ID);
    return;
}
sub get_Source_ID {
    return $Source_ID;
}

sub set_Taxid {
    my ($new_Taxid) = @_;
    $Taxid = $new_Taxid;
    return;
}
sub get_Taxid {
    return $Taxid;
}

sub set_LB {
    my ($new_LB) = @_;
    $LB = $new_LB;
    return;
}
sub get_LB {
    return $LB;
}

sub set_UB {
    my ($new_UB) = @_;
    $UB = $new_UB;
    return;
}
sub get_UB {
    return $UB;
}

sub get_extra {
    return $extra;
}
sub set_extra {
    my ($new_extra) = @_;
    $extra = $new_extra;
    return;
}


sub parse_Model_Info {
    my ($SBML_model, $SBML_input, $SBML_object) = @_;

    my $SBML_filename = basename($SBML_input);
    $SBML_filename =~ s{\.(sb|x)ml\d?.*?$}{};
    Model::set_ID( $SBML_filename );
    Model::set_Description( $SBML_model->getName() || '' );
    Model::set_Source_ID( $SBML_model->getId()     || '' );

    #Parse notes
    if ( $SBML_model->isSetNotes() ){
        Model::set_extra( Notes::parse_notes($SBML_model, Model::get_extra()) );
    }

    #Parse annotations
    #NOTE annotations overload notes!
    if ( $SBML_model->isSetAnnotation() ){
        Model::set_extra( Annotation::parse_annotation($SBML_model, Model::get_extra(), 'model') );
    }

    #Get LB/UB from reaction parameters or model parameters
    Reactions::get_min_max_bounds($SBML_model, $SBML_object);


    my $extra = Model::get_extra();
    $Taxid = exists $extra->{'taxonomy_model'}->[0] ? join(';', uniq sort {$a <=> $b} map { s/^.+?://; $_ } @{ $extra->{'taxonomy_model'} }) : $Taxid;
    $UB    = $UB    // $Constants::MAX_BOUND_VAL;
    $LB    = $LB    // $Constants::MIN_BOUND_VAL;
    Model::set_Taxid($Taxid);
    Model::set_UB($UB);
    Model::set_LB($LB);

    return;
}

1;

