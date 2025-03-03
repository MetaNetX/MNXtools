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


#TODO SBML validator
#     SBO term 'SBO:0000289' on the <compartment> is not in the appropriate branch
#     SBO term 'SBO:0000649' on the <species>     is not in the appropriate branch
#
#TODO -validate         Validate SBML (check consistency)
#     option returns *Segmentation fault*
#
#TODO MEMOTE to do
#
#TODO allow a model to be kegg-oriented for example:
#     without conflicts the id will be kegg ids

    # <listOfReactions>
    # Add reactions in the model, and subsequently the linked chemical species and compartments
    my @reac_obj = $MetNet->get_growth_reac_ids($mnet_id);
    for my $reac_id ( uniq nsort $MetNet->select_reac_ids() ){
        Reactions::create_SBML_reaction($MetNet, $mnet_id, $SBML_model, $use_notes, $reac_id, $all_bounds, \@reac_obj);
    }
    #check reaction lists between SBML and MetNet
    my %diff_reactions = map { $_->getId() => -1 } $SBML_model->getListOfReactions();
    map { $diff_reactions{ $_ }++ } $MetNet->select_reac_ids();
    for my $missed_reac ( grep { $diff_reactions{ $_ } != 0 } keys %diff_reactions ){
        if ( $diff_reactions{ $missed_reac } < 0 ){
            warn "[$missed_reac] reaction in SBML but not in TSV\n";
        }
        elsif ( $diff_reactions{ $missed_reac } > 0 ){
            warn "[$missed_reac] reaction in TSV but not in SBML\n";
        }
    }

    # <listOfCompartments>
    #check compartment lists between SBML and MetNet
    my %diff_compartments = map { $_->getId() => -1 } $SBML_model->getListOfCompartments();
    map { $diff_compartments{ $_ }++ } $MetNet->select_comp_ids();
    for my $missed_comp ( grep { $diff_compartments{ $_ } != 0 } keys %diff_compartments ){
        if ( $diff_compartments{ $missed_comp } < 0 ){
            warn "[$missed_comp] compartment in SBML but not in TSV\n";
        }
        elsif ( $diff_compartments{ $missed_comp } > 0 ){
            warn "[$missed_comp] compartment in TSV but not in SBML\n";
        }
    }


    # <listOfSpecies>
    #check chemical species lists between SBML and MetNet
    my %diff_chemicals = map { my $cid = $_->getId(); $cid =~ s{__64__.+$}{}; $cid => -1 } $SBML_model->getListOfSpecies();
    map { $diff_chemicals{ $_ }++ } $MetNet->select_chem_ids();
    for my $missed_chem ( grep { $diff_chemicals{ $_ } != 0 } keys %diff_chemicals ){
        if ( $diff_chemicals{ $missed_chem } < 0 ){
            warn "[$missed_chem] chemical in SBML but not in TSV\n";
        }
        elsif ( $diff_chemicals{ $missed_chem } > 0 ){
            warn "[$missed_chem] chemical in TSV but not in SBML\n";
        }
    }


    # <fbc:listOfGeneProducts> setLabel/setName
    if ( $SBML_ns->getLevel() >= 3 ){
        my $gpr_fbc = $SBML_model->getPlugin('fbc');
        if ( $gpr_fbc ){
            my $gpr = $gpr_fbc->getListOfGeneProducts();
            for (my $i = 0; $i < $gpr->getNumGeneProducts(); $i++){
                Peptides::update_GeneProducts($gpr->get($i), $MetNet);
            }

            #check gpr lists between SBML and MetNet
            my $model_gpr_nbr = $gpr->getNumGeneProducts() - 1;
            my %diff_gprs = map { $gpr->get($_)->getId() => -1 } 0..$model_gpr_nbr;
            map { $diff_gprs{ $_ }++ } $MetNet->select_pept_ids();
            for my $missed_gpr ( keys %diff_gprs ){
                if ( $diff_gprs{ $missed_gpr } < 0 ){
                    warn "[$missed_gpr] gpr in SBML but not in TSV\n";
                }
                elsif ( $diff_gprs{ $missed_gpr } > 0 ){
                    warn "[$missed_gpr] gpr in TSV but not in SBML\n";
                }
            }
        }
    }
    #TODO add peptide notes only for SBML level 2 ???


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

