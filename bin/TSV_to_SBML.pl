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
#FIXME id must be SBML compliant: do not contain ':', ...
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
        if ( !$SBML_model->isSetMetaId() ){
            #NOTE its vital to set a metaId first, otherwise the cvterms can't be added
            #NOTE Annotations can only be added to elements, that have a metaId set
            $SBML_model->setMetaId( $SBML_model->getId() );
        }
        my $CV = new LibSBML::CVTerm();
        $CV->setQualifierType($LibSBML::BIOLOGICAL_QUALIFIER);
        $CV->setBiologicalQualifierType($LibSBML::BQB_HAS_TAXON);
        my $annotation = Formaters::guess_annotation_link('other', "taxon:$taxid");
        if ( $annotation ){
            $CV->addResource($annotation);
            $SBML_model->addCVTerm($CV);
        }
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


#TODO allow a model to be kegg-oriented for example:
#     without conflicts the id will be kegg ids

#FIXME from reac list, list chem and comp (and geneproducts) AND then warn if defined in reac but not in chem/comp TSVs!
    # <listOfReactions>
    my (%reac_chem, %reac_comp, %reac_gpr);#TODO
    for my $reac_id ( uniq nsort $MetNet->select_reac_ids() ){
        Reactions::create_SBML_reaction($MetNet, $mnet_id, $SBML_model, $use_notes, $reac_id, $all_bounds);
    }


    # <listOfCompartments>
    for my $comp_id ( uniq nsort $MetNet->select_comp_ids() ){
#        Compartments::create_SBML_compartment($MetNet, $mnet_id, $SBML_model, $use_notes, $comp_id);
    }


    # <listOfSpecies>
    for my $chem_id ( uniq nsort $MetNet->select_chem_ids() ){
#        Chemicals::create_SBML_chemical($MetNet, $mnet_id, $SBML_model, $use_notes, $chem_id);
    }


    # <fbc:listOfGeneProducts> setLabel/setName
    if ( $SBML_ns->getLevel() >= 3 ){
        my $gpr_fbc = $SBML_model->getPlugin('fbc');
        if ( $gpr_fbc ){
            my $gpr = $gpr_fbc->getListOfGeneProducts();
            for (my $i = 0; $i < $gpr->getNumGeneProducts(); $i++){
                my $gpr_id = $gpr->get($i)->getId();
                my ($gpr_desc, $gpr_xrefs, $gpr_gene) = $MetNet->get_pept_info($gpr_id);
                if ( $gpr_desc ){
                    $gpr->get($i)->setName($gpr_desc);
                }
                if ( $gpr_xrefs ){
                    $gpr->get($i)->setLabel($gpr_xrefs);
                }
            }
        }
    }


#TODO <fbc:listOfObjectives>  -> <fbc:listOfFluxObjectives>


    # Check model
    SBML::check_SBML($SBML_doc, $check_SBML, $validate_SBML);
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

