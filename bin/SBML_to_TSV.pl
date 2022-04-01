#!/usr/bin/env perl

# Convert SBML format to MetaNetX TSV files format

use strict;
use warnings;
use diagnostics;

use Data::Dumper;

use File::Basename;
use Getopt::Long;
use List::Util qw( uniq );
use LibSBML;
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

my $VERSION = '0.1.4';
my $soft_name = basename $0;


# Initialize and parse command line arguments
my ($debug, $help, $version)     = (0, 0, 0);
my ($SBML_input)                 = ('');
my ($check_SBML, $validate_SBML) = (1, 0);
my ($outdir)                     = ('.');
GetOptions ('debug'          => \$debug,
            'help|H|?'       => sub { help(0); },
            'version|V'      => sub { print "$soft_name\n$Constants::copyright\nVersion: $VERSION\n"; },
            'sbml=s'         => \$SBML_input,
            'outdir=s'       => \$outdir,
            'check!'         => \$check_SBML,
            'validate'       => sub {$validate_SBML = 1; $check_SBML = 1; },
           )
           or help(1, 'Error in command line arguments');

if ( $SBML_input && (!-e "$SBML_input" || !-r "$SBML_input" || !-s "$SBML_input") ){
    help(1, "Cannot read [$SBML_input] file: it does not exist, is not readable or empty");
}

Constants::set_debug($debug);


##############################################################################
# SBML to MetaNetX
if ( $SBML_input ){
    my $SBML_reader = new LibSBML::SBMLReader();
    my $SBML_object = $SBML_reader->readSBML($SBML_input);

    if ( !-e "$outdir" ){
        $Constants::logger->fatal("Directory does not exist: [$outdir]");
        exit 1;
    }
#FIXME not compatible the way import works on the web site
#    system( "find $outdir/ -mindepth 1 -maxdepth 1 | xargs rm -rf" ); # Remove everything in previously existing $outdir folder
    # Check model
    SBML::check_SBML($SBML_object, $check_SBML, $validate_SBML);


    # Get general model information
    #NOTE Order of parsing is important because chemicals with oundaryCondition=true will create the boundary compartment
    Model::parse_Model_Info($SBML_object->getModel(), $SBML_input, $SBML_object);
    Chemicals::parse_Chemicals_Info($SBML_object->getModel(), $SBML_object);
    Reactions::parse_Reactions_Info($SBML_object->getModel(), $SBML_object);
    Compartments::parse_Compartments_Info($SBML_object->getModel());
    Peptides::parse_GeneProducts_Info($SBML_object->getModel(), $SBML_object);

    # Find biomass reaction(s) if any
    Reactions::find_growth_reaction($Reactions::reactions);

    # Write TSV files
    my $MetNet = MetNet->new();
    my $model_id = Model::get_ID();
    $MetNet->add_mnet( $model_id, Model::get_LB(), Model::get_UB() );
    $MetNet->set_mnet_desc( $model_id, Model::get_Description() || '' );
    $MetNet->push_mnet_info( $model_id,
        [ 'Source ID',      Model::get_Source_ID()       || '' ],
        [ 'Processed Date', Model::get_Processed_Date()        ],
        [ 'Taxid',          Model::get_Taxid()           || '' ],
    );
    my %is_objective = (); # to keep track of objective under the new ID
    my $counter = 0;
    for my $source ( sort keys %$Reactions::reactions ){
        my $reaction = $Reactions::reactions->{$source};
        my $reac_id = $source;
        unless( $reac_id=~ /^[a-zA-Z]\w{1,31}$/ ){ # because some "old" models are still using invalid SMBL identifers
            $reac_id = MetNet::compute_reac_id( $reaction->{equation}) . '_' . ++$counter;
        }
        $is_objective{$source} = 1  if $reaction->{'is_objective'};
        my $dir = $reaction->{'reversible'} ? 'B' : 'LR' ; # SBML does not know 'RL'
        my $LB = $reaction->{'LB'} // Model::get_LB();
        $LB = 0  if ( $LB < 0 && $dir eq 'LR' );
        $MetNet->add_reac_add_enzy($model_id,
                                   $reac_id,
                                   $reaction->{'equation'},
                                   '',
                                   ($reaction->{'peptides'} || ''),
                                   $LB,
                                   ($reaction->{'UB'} // Model::get_UB()),
                                   $dir);
        $MetNet->set_reac_source( $model_id, $reac_id, $source);
        $MetNet->set_reac_info( $reac_id,
                ( $reaction->{'EC'}       ? join( ';', uniq sort @{$reaction->{'EC'}} )       : '' ),
                ( $reaction->{'pathways'} ? join( ';', uniq sort @{$reaction->{'pathways'}} ) : '' ),
                ( $reaction->{'xrefs'}    ? join( ';', uniq sort @{$reaction->{'xrefs'}} )    : '' ),
            );
    }
#    print Dumper $Peptides::peptides;
    for my $source ( keys %$Peptides::peptides ){
        my $peptide = $Peptides::peptides->{$source};
        $MetNet->set_pept_info($source, $peptide->{'name'}, join(';', uniq sort @{ $peptide->{'xrefs'} }), $peptide->{'label'});
    }
#    print Dumper $Chemicals::chemicals;
    for my $source ( keys %$Chemicals::chemicals ){
        unless ( exists $MetNet->{'look'}{'chem'}{$source} ){
            $Constants::logger->warn("Declared chemical [$source] not used");
            next;
        }
        my $chemical = $Chemicals::chemicals->{$source};
        $MetNet->set_chem_info($source, $chemical->{'name'}, $chemical->{'formula'}, $chemical->{'mass'}, $chemical->{'charge'}, join(';', uniq sort @{ $chemical->{'xrefs'} }));
        $MetNet->set_chem_source($model_id, $source, join(';', uniq sort @{ $chemical->{'source'} }));
    }
#    print Dumper $Compartments::compartments;
    for my $source ( keys %$Compartments::compartments ){
        unless ( exists $MetNet->{'look'}{'comp'}{$source} ){
            $Constants::logger->warn("Declared compartment [$source] not used");
            next;
        }
        my $compartment = $Compartments::compartments->{$source};
        eval { $MetNet->_get_comp_num( $source ) };
        if( $@ ){
            $Constants::logger->warn('Ignoring compartment not used in reaction: ' . $source);
            #$tb->warn( 'Ignoring compartment not used in reaction: ' . $source );
        }
        else{
            $MetNet->set_comp_info(
                $source,
                $compartment->{'name'},
                ( $compartment->{'xrefs'}    ? join( ';', uniq sort @{$compartment->{'xrefs'}} ) : '' ),
            );
            $MetNet->set_comp_source( $model_id, $source, $compartment->{'source'});
        }
    }


    # try hard to fix ambiguities with model objectives:
    my @growth = $MetNet->get_growth_reac_ids( $model_id ); # already smartly sorted: the first one is the best candidate
    if( @growth > 1 ){
        my %growth_info = ();
        my $obj_count   = 0;
        foreach( @growth ){
            $growth_info{$_} = [ $MetNet->get_enzy_info( $model_id, $_ ) ];
            $obj_count++  if exists $is_objective{$_};
        }
        if( $obj_count == 0 ){ # no objective found
            $is_objective{$growth[0]} = 1;
        }
        elsif( $obj_count == 1 ){
            # already ok
        }
        else{
            # $tb->warn( "Multiple objective... chosing arbitrarily: $growth[0]!!!\n";
            $Constants::logger->warn( "Multiple objective... chosing arbitrarily: $growth[0]!!!" );
            %is_objective = (); # reset
            $is_objective{$growth[0]} = 1;
        }
        foreach( @growth ){
            $growth_info{$_}[1] = 0;
            $growth_info{$_}[2] = $is_objective{$_} ? 'NA' : 0;
        }
        $MetNet->rename( $model_id, '#tmp' );
        $MetNet->add_mnet( $model_id, $MetNet->get_LU_bounds( '#tmp' ));
        $MetNet->set_mnet_desc( $model_id, $MetNet->get_mnet_desc( '#tmp' ));
        $MetNet->push_mnet_info( $model_id, $MetNet->get_mnet_info( '#tmp' ));
        foreach my $reac_id ( $MetNet->select_reac_ids( mnet => '#tmp' )){
            if( exists $growth_info{$reac_id} ){
                $MetNet->copy_reac_add_enzy( $model_id, $reac_id, @{$growth_info{$reac_id}} );
            }
            else{
                $MetNet->copy_reac_copy_enzy( '#tmp', $model_id, $reac_id );
            }
        }
        foreach my $chem_id ( $MetNet->select_chem_ids( mnet => '#tmp' )){
            $MetNet->set_chem_source( $model_id, $chem_id, $MetNet->get_chem_source( '#tmp', $chem_id ));
        }
        foreach my $comp_id ( $MetNet->select_comp_ids( mnet => '#tmp' )){
            $MetNet->set_comp_source( $model_id, $comp_id, $MetNet->get_comp_source( '#tmp', $comp_id ));
        }
        foreach my $reac_id ( $MetNet->select_reac_ids( mnet => '#tmp' )){
            $MetNet->set_reac_source( $model_id, $reac_id, $MetNet->get_reac_source( '#tmp', $reac_id ));
        }
    }
    $MetNet->write( $outdir, $model_id );
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
$soft_name converts an SBML file to the MetaNetX TSV format

$soft_name [options] -sbml <SBML file to convert to MetaNetX TSV format> -outdir <dir>

         -sbml|s   <file>  SBML file to convert to MetaNetX TSV format

SBML to MetaNetX TSV options:
         -outdir   <dir>   Where to store produced TSV files

General options:
         -check|c          Check SBML (enabled by default)
         -validate         Validate SBML (check consistency)
         -debug|d          Show debug information
         -h                print this help
         -V                print version information
EOF
    exit $exit_code;
}

