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

my $usage = "$0 [options] <in-dir> <out-dir>

options: -h          this help

Nota bene: the content of <out-dir> is erased!

         -V <file-bindump> read namespace cache file ( ../cache/ChemSpace.bindump by default )

         -p <dir>          path to dir with id_map and peptide.tsv files

         -L          use late merge data
         -x <regexp> use chem xref for mapping ('.' means all, all by default)
         -y <regexp> use pept xref for mapping ('.' means all, all by default)
         -G          use generic compartment
";

my %opt;  # GLOBAL: to store the options
Getopt::Std::getopts( 'hV:p:Lx:y:G', \%opt ); # FIXME: implement -D options
if( $opt{h} or @ARGV != 2 ){
    print "$usage\n";
    exit 1;
}
$opt{V} = '../cache/ChemSpace.bindump' unless $opt{V};

$tb->die( 'Option -V <*.bindump> is mandatory!' ) unless $opt{V};
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
$convertor->convert( $metnet, $model_name, $metnet2, $model_name, $option );
$tb->report( 'write GSMN', $mapped_dir );
$metnet2->write( $mapped_dir, $model_name );
my $fh = $tb->open( "> $mapped_dir/convert.log" );
print $fh join( "\n", sort $convertor->get_logs() ) . "\n";
close $fh;
my $yaml = $tb->open( "> $mapped_dir/convert_log.yaml" );
print $yaml join( "\n", @{$convertor->{log4yaml}} ) . "\n";
close $yaml;


exit 0;

