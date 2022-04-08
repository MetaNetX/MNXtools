#!/usr/bin/env perl

# Convert MetaNetX TSV files format to SBML

use strict;
use warnings;
use diagnostics;

use Data::Dumper;

use File::Basename;
use Getopt::Long;
use List::Util qw( uniq );
use LibSBML;
use Sort::Naturally;
use FindBin qw( $RealBin );
use lib "$RealBin/../perl";
use lib "$RealBin/../perl/SBML";
use Constants;
use SBML;
use Model;
use Compartments;
use Chemicals;
use Reactions;
use Peptides;
use MetNet;
#use Toolbox;
#my $tb = Toolbox->new();

my $VERSION = '0.0.1';
my $soft_name = basename $0;


# Initialize and parse command line arguments
my ($debug, $help, $version)     = (0, 0, 0);
my ($TSV_directory)              = ('');
my ($check_SBML, $validate_SBML) = (1, 0);
my ($use_notes)                  = (1);
GetOptions ('debug'          => \$debug,
            'help|H|?'       => sub { help(0); },
            'version|V'      => sub { print "$soft_name\n$Constants::copyright\nVersion: $VERSION\n"; },
            'tsv|metanetx=s' => \$TSV_directory,
            'notes|n!'       => \$use_notes,
            'check!'         => \$check_SBML,
            'validate'       => sub {$validate_SBML = 1; $check_SBML = 1; },
           )
           or help(1, 'Error in command line arguments');

if ( $TSV_directory && (!-e "$TSV_directory" || !-d "$TSV_directory" || !-r "$TSV_directory") ){
    help(1, "Cannot read files in [$TSV_directory] directory: directory does not exist or is not readable");
}

Constants::set_debug($debug);


##############################################################################
# MetaNetX TSV to SBML
if ( $TSV_directory ){
    my $MetNet = MetNet->new();
    my $mnet_id = $MetNet->load( "$TSV_directory" );

    # Create SBML document
#FIXME SBML 3.1 FBC 2 only now!
    my $SBML_ns = new LibSBML::SBMLNamespaces($Constants::SBML_level, $Constants::SBML_version);
    $SBML_ns->addPackageNamespace('fbc', $Constants::FBC_version)   if ( $SBML_ns->getLevel() >= 3 );
    #Create an empty SBML document with additional XML namespace for XML notes
    $SBML_ns->addNamespace('http://www.w3.org/1999/xhtml', 'html')  if ( $use_notes );
    my $SBML_doc = new LibSBML::SBMLDocument($SBML_ns);
    $SBML_doc->setPackageRequired('fbc', 0)  if ( $SBML_ns->getLevel() >= 3 );


    # Set model attributes
    #sid
    my ($model_id) = map { $_->[1] } grep { $_->[0] =~ /Source ID/ } $MetNet->get_mnet_info($mnet_id);
    my $SBML_model = $SBML_doc->createModel($model_id // $mnet_id);
    #fbc in model
    if ( $SBML_ns->getLevel() >= 3 ){
        my $SBML_model_fbc = $SBML_model->getPlugin('fbc');
        $SBML_model_fbc->setStrict(1);
    }
    #name
    my ($model_name) = map { $_->[1] } grep { $_->[0] =~ /Description/ } $MetNet->get_mnet_info($mnet_id);
    if ( $model_name ){
        $SBML_model->setName($model_name);
    }
    #fbc model SBO
    $SBML_model->setSBOTerm($Constants::SBO_FLUXBALANCEFRMW);
    #TODO metaid???
    #notes
    if ( $use_notes ){
        my $notes = '';
        for my $info ( $MetNet->get_mnet_info($mnet_id) ){
            $notes .= '<html:p>'.$info->[0].': '.$info->[1].'</html:p>'  if ( $info->[1] );
        }
        if ( $MetNet->get_LU_bounds($mnet_id) ){
            my ($LB, $UB) = $MetNet->get_LU_bounds($mnet_id);
            $notes .= '<html:p>LB: '.$LB.'</html:p>';
            $notes .= '<html:p>UB: '.$UB.'</html:p>';
        }
        $SBML_model->setNotes($notes);
    }
    #annotations
    my ($taxid) = map { $_->[1] } grep { $_->[0] =~ /Taxid/ } $MetNet->get_mnet_info($mnet_id);
    if ( $taxid && $taxid =~ /^\d+$/ ){
        my $CV = new LibSBML::CVTerm();
        $CV->setQualifierType($LibSBML::BIOLOGICAL_QUALIFIER);
        $CV->setBiologicalQualifierType($LibSBML::BQB_HAS_TAXON);
        $CV->addResource($Constants::identifiers_taxid.$taxid);
        $SBML_model->addCVTerm($CV);
        $SBML_model->setAnnotation($CV);
        $SBML_model->appendAnnotation();
#FIXME swig issue that does not make setAnnotation available???
    }


    # UnitDefinitions / Parameters
    my $udef = $SBML_model->createUnitDefinition();
    $udef->setId($Constants::unit);
    # kind, exponent, multiplier, scale
    for my $u ( [sprintf("%d", $LibSBML::UNIT_KIND_MOLE), 1, 1, -3], [sprintf("%d", $LibSBML::UNIT_KIND_GRAM), -1, 1, 0], [sprintf("%d", $LibSBML::UNIT_KIND_SECOND), -1, 3600, 0] ){
        my $unit = $udef->createUnit();
        $unit->setKind($u->[0]);
        $unit->setExponent($u->[1]);
        $unit->setScale($u->[3]);
        $unit->setMultiplier($u->[2]);
    }

    my ($LB, $UB) = ($Constants::MIN_BOUND_VAL, $Constants::MAX_BOUND_VAL);
    my $other_bounds;
    if ( $MetNet->get_reac_bounds($mnet_id) ){
        my ($LB_all, $UB_all);
        ($LB_all, $UB_all, $LB, $UB) = $MetNet->get_reac_bounds($mnet_id);
        push @{ $other_bounds }, uniq sort grep { $_ != $LB && $_ != $UB } (values %$LB_all, values %$UB_all);
    }
    my $all_bounds;
    push @{ $all_bounds }, [$LB, $Constants::MIN_BOUND_ID, 'default'];
    push @{ $all_bounds }, [$UB, $Constants::MAX_BOUND_ID, 'default'];
    my $extra_bound_count = 0;
    for my $other_bound ( @{ $other_bounds } ){
        my $id  = $other_bound == $Constants::ZERO_BOUND_VAL ? $Constants::ZERO_BOUND_ID : 'extrabound'.++$extra_bound_count;
        my $sbo = $other_bound == $Constants::ZERO_BOUND_VAL ? 'default' : '';
        push @{ $all_bounds }, [$other_bound, $id, $sbo];
    }
    for my $B ( @$all_bounds ){
        my $param = $SBML_model->createParameter();
        $param->setId($B->[1]);
        $param->setValue($B->[0]);
        $param->setConstant(1); # True
        if ( $B->[2] eq 'default' ){
            $param->setSBOTerm($Constants::SBO_FLUXDEFAULT);
        }
        else {
            $param->setSBOTerm($Constants::SBO_FLUX);
        }
        $param->setUnits($Constants::unit);
    }


    # <listOfCompartments>
    for my $comp_id ( uniq nsort $MetNet->select_comp_ids() ){
        my $comp = $SBML_model->createCompartment();
        $comp->setId($comp_id);
        $comp->setConstant(1); # True
        my ($comp_name, $comp_xref) = $MetNet->get_comp_info($comp_id);
        $comp->setName($comp_name);
        $comp->setSBOTerm( $comp_id eq $Constants::boundary_comp_id ? $Constants::boundary_comp_sbo : $Constants::default_comp_sbo );
        #notes
        if ( $use_notes ){
            my $notes = '';
            if ( $MetNet->get_comp_source($mnet_id, $comp_id) ){
                $notes .= '<html:p>SOURCE: '.$MetNet->get_comp_source($mnet_id, $comp_id).'</html:p>';
            }
            if ( $comp_xref ){
                $notes .= '<html:p>REFERENCE: '.$comp_xref.'</html:p>';
            }
            $comp->setNotes($notes)  if ( $notes );
        }
        #TODO annotations: fix is/isRelatedTo and different identifiers.org...
        if ( $comp_xref ){
            my $CV = new LibSBML::CVTerm();
            $CV->setQualifierType($LibSBML::BIOLOGICAL_QUALIFIER);
            $CV->setBiologicalQualifierType($LibSBML::BQB_IS);
            $CV->addResource($Constants::identifiers_go.$comp_xref);
            $comp->addCVTerm($CV);
            #$comp->setAnnotation($CV);
            #$comp->appendAnnotation();
        }
    }


#TODO <listOfSpecies>
#TODO <listOfReactions>
#TODO <fbc:listOfObjectives>  -> <fbc:listOfFluxObjectives>
#TODO <fbc:listOfGeneProducts>


    # Write SBML
    my $SBML_writer = new LibSBML::SBMLWriter();
    $SBML_writer->setProgramName($soft_name);
    $SBML_writer->setProgramVersion($VERSION);
    print $SBML_writer->writeSBMLToString($SBML_doc);
}
else {
    help(0);
}

exit 0;



sub help {
    my ($exit_code, $msg) = @_;

    if ( $exit_code != 0 ){
        *STDOUT = *STDERR;
    }

    if ( $msg ){
        print "\n\t$msg\n\n";
    }
    print <<"EOF";
$Constants::copyright
Version: $VERSION
$soft_name converts MetaNetX TSV files to SBML (by default SBML $Constants::SBML_level.$Constants::SBML_version - FBCv$Constants::FBC_version)

$soft_name [options] -tsv <TSV repository with MetaNetX TSV files to convert to SBML>

         -tsv|t    <dir>   TSV repository with MetaNetX TSV files to convert to SBML

MetaNetX TSV to SBML options:
         -notes|n          Use SBML notes (enabled by default)

General options:
         -check|c          Check SBML (enabled by default)
         -validate         Validate SBML (check consistency)
         -debug|d          Show debug information
         -h                print this help
         -V                print version information
EOF
    exit $exit_code;
#TODO Add other TSV2SBML options
}

