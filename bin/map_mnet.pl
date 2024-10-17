#!/usr/bin/env perl

# ------------------------------------------------ #
#
# ------------------------------------------------ #

use strict;

use POSIX qw( _exit );
use Data::Dumper;
use Getopt::Std;
use FindBin qw( $RealBin );
use Cwd qw( abs_path );

use lib "$RealBin/../perl";
# use MNXref qw( read_taxo_space );
use ChemSpace;
use GeneSpace;
use Convert;

use lib "$RealBin/../util";
use Toolbox;
my $tb = Toolbox->new();

# ------------------------------------------------------
# Read command line and check its arguments
# ------------------------------------------------------

my $usage = "$0 [options] <indir> <outdir>

<indir>   directory with a mnet in MetaNetX TSV format
<outdir>  directory to write the mapped mnet and the yaml mapping
          report.
          Nota Bene: the content of <outdir> is completely erased!

options:

    -h            this help

    -V <dumpfile> read namespace cache file
                  (/usr/local/MNXtools/cache/ChemSpace.bindump by default)
    -p <dir>      path to dir with id_map and peptide.tsv files

    # -L            use late merge data

    -c <list>     a comma/slash separated list of prefixes that specify which
                  metabolite cross-references should be considered.
                  Comma separated list are evaluated alltogether thus possibly
                  generating conflicts.
                  Slash separated list are evaluated successively, thus
                  possibly generating discrepency

                  In addition the following keywords are defined:

                  ID    Utilize the chemichal identifer (e.g. the SBML identifer)
                  XREFS Utilize all cross references (not recommened)

                  For example:

                  -c ID                (this is the default when -c is not set)
                  -c ID,XREF           (use all information available, possibly generating many conflicts)
                  -c chebi/kegg/ID     (it should speak for itself)
                  -c chebi/slm/ID      (ditto)

    -y <regexp>   use pept xref for mapping
                  ('.' means alli, which is the default)
    -G            produce reactions with MetaNetX generic compartments,
                  destroying the connectivity of compartimentalized models.
";

my %opt;  # GLOBAL: to store the options
Getopt::Std::getopts( 'hV:p:Gc:', \%opt );
if( $opt{h} or @ARGV != 2 ){
    print "$usage\n";
    exit 1;
}

$opt{V} = '/usr/local/MNXtools/cache/ChemSpace.bindump' unless $opt{V};
$tb->die( 'Option -V <*.bindump> is mandatory!' ) unless $opt{V};
$opt{c} = 'ID' unless $opt{c};

my( $unmapped_dir, $mapped_dir ) = @ARGV;
foreach( $unmapped_dir, $mapped_dir ){
    $tb->die( "Dir does not exists: $_" )  unless -d $_;
}
$tb->system( "rm -rf $mapped_dir/*" );

# ------------------------------------------------ #
# Prepare and/or retrieve namespace object
# ------------------------------------------------ #

$tb->report( 'retrieve', $opt{V} );
my $ns = ChemSpace->new( $tb )->retrieve_from_file( $opt{V} );

$ns->enable_late_merge()    if $opt{L};
# $ns->enable_un_merge()       if $opt{U};
# $ns->disable_normalize_coef  if $opt{C};

my $option = {};
$tb->die( 'Options -x <regexp> and -X <regexp> are mutually exclusive!' ) if $opt{x} and $opt{X};
if( $opt{x} ){
    $option->{use_chem_xref} = $opt{x};
}
else{
    $option->{use_chem_ID} = 1;
    $option->{use_chem_xref} = $opt{X} if $opt{X};
}

# $option->{max_range} = $opt{R}  if $opt{R};
# $option->{prefix}    = $opt{P}  if $opt{P};

# ------------------------------------------------ #
# Do the job
# ------------------------------------------------ #

my $gene_space = GeneSpace->new( $opt{p} );
my $convertor  = Convert->new( $ns, $gene_space );
my $metnet     = MetNet->new();
$tb->report( 'load GSMN', $unmapped_dir );
my $model_name = $metnet->load( $unmapped_dir );
my $metnet2 = MetNet->new();
$convertor->convert( $metnet, $model_name, $metnet2, $model_name, $opt{c}, $option );
$tb->report( 'write GSMN', $mapped_dir );
$metnet2->write( $mapped_dir, $model_name );
my $yaml = $tb->open( "> $mapped_dir/convert_log.yaml" );
print $yaml join( "\n", @{$convertor->{log4yaml}} ) . "\n";
close $yaml;

_exit 0;

