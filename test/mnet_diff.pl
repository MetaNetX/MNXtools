#!/bin/env perl

use strict;

use Data::Dumper;
use Getopt::Std;
use FindBin qw( $RealBin );

use lib "$RealBin/../../perl/util";
use Toolbox;
my $tb = Toolbox->new();

my $usage = "$0 [option] col_format <dir_1> <dir_2>

options: -h          this help
";

my %opt;  # GLOBAL: to store the options
Getopt::Std::getopts( 'h', \%opt );
if( $opt{'h'} or @ARGV != 2 ){
    print "$usage\n";
    exit 0;
}
my( $dir_1, $dir_2 ) = @ARGV;
$tb->system( "$RealBin/tab_diff.pl ss     $dir_1/compartments.tsv  $dir_2/compartments.tsv" );
$tb->system( "$RealBin/tab_diff.pl ss.sss $dir_1/chemicals.tsv     $dir_2/chemicals.tsv"    );
$tb->system( "$RealBin/tab_diff.pl ss.s   $dir_1/reactions.tsv     $dir_2/reactions.tsv"    );
$tb->system( "$RealBin/tab_diff.pl sssss  $dir_1/enzymes.tsv       $dir_2/enzymes.tsv"      );
$tb->system( "$RealBin/tab_diff.pl ss     $dir_1/peptides.tsv      $dir_2/peptides.tsv"      );

