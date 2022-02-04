#!/usr/bin/env perl

use strict;

use Data::Dumper;
use Getopt::Std;
use FindBin qw( $RealBin );

use lib "$RealBin/../perl";
use Toolbox;
my $tb = Toolbox->new();

my $usage = "$0 [option] col_format <file_1> <file_2>

options: -h          this help

         -t <real>   tolerance (absolute value) for real number
                     comparison (1e-6 by default)
         -z          Zombie mode: do not die on error!

In addition the environment variable

    COPY_RESULT_AS_NEW_REF

can be set to 1 to create new reference dataset
";

my %opt;  # GLOBAL: to store the options
Getopt::Std::getopts( 'ht:z', \%opt );
if( $opt{'h'} or @ARGV != 3 ){
    print "$usage\n";
    exit 0;
}
$opt{t} = 1e-6 unless $opt{t};

my( $format, $file_1, $file_2 ) = @ARGV;

if( $ENV{COPY_RESULT_AS_NEW_REF} == 1 ){
    $tb->warn( "COPYING RESULT AS NEW REF!" );
    $tb->system( "cp $file_1 $file_2" );
}

$tb->report( 'format', $format );
my @format = split //, $format;
my $tsv_1 = $tb->scan_tsv( $file_1 );
my $tsv_2 = $tb->scan_tsv( $file_2 );
my $L = @$tsv_1;
$tb->die( "The two input files have a different number of lines!" ) unless $L == @$tsv_2;
foreach my $i ( 0 .. @format -1 ){
    if( $format[$i] eq 's' ){ # string comparison
        foreach( 0 .. $L - 1 ){
            $tsv_1->[$_][$i] =~ s/_(B|LR|RL)$//;
            $tsv_2->[$_][$i] =~ s/_(B|LR|RL)$//;
            unless( $tsv_1->[$_][$i] eq $tsv_2->[$_][$i] ){
                msg(
                    "Strings differ: $tsv_1->[$_][$i] vs $tsv_2->[$_][$i] at line "
                    . ( $_ + 1 )
                    . ', col '
                    . ($i + 1 )
                );
            }
        }
    }
    elsif( $format[$i] eq 'i' ){ # integer comparison
        foreach( 0 .. $L - 1 ){
            unless( $tsv_1->[$_][$i] == $tsv_2->[$_][$i] ){
                msg(
                    "Numbers not equal: $tsv_1->[$_][$i] vs $tsv_2->[$_][$i] at line "
                    . ( $_ + 1 )
                    . ', col '
                    . ( $i + 1 )
                );
            }
        }
    }
    elsif( $format[$i] eq 'r' ){ # real comparison (with some tolerance)
        foreach my $j ( 0 .. $L - 1 ){
            if( abs( $tsv_1->[$j][$i] - $tsv_2->[$j][$i] ) > $opt{t} ) {
                msg(
                    "Numbers differ too much: $tsv_1->[$j][$i] vs $tsv_2->[$j][$i] at line "
                    . ( $j + 1 )
                    . ', col '
                    . ( $i + 1 )
                );
            }
        }
    }
    elsif( $format[$i] eq '.' ){ # skip that column

    }
    else{
        $tb->die( "Unsupported column format: $format[$i]" );
    }
}

sub msg{
    my $msg = shift;
    if( $opt{z} ){
        $tb->warn( $msg );
    }
    else{
        $tb->die( $msg );
    }
}
