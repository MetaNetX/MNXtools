#!/usr/bin/env perl

use strict;

use Data::Dumper;
use Getopt::Std;
use FindBin qw( $RealBin );

use lib "$RealBin/../perl";
use Toolbox;
my $tb = Toolbox->new();

my $usage = "$0 [option] test_file...

options: -h           this help
         -t <regexp>  execute the test with matching name(s) (all by default)
";

my %opt;  # GLOBAL: to store the options
Getopt::Std::getopts( 'ht:', \%opt );
if( $opt{'h'} or @ARGV == 0 ){
    print "$usage\n";
    exit 0;
}
$opt{t} = '.'  unless $opt{t};

my %test_info = ();
my %elapsed   = ();
my $rank = 0;
my $test_id = '';
while( <> ){
    next if /^\#/;
    next if /^\s*$/;
    if( /^TEST=(.+)/ ){
        $rank++;
		$test_id = $1;
		$tb->die( "Duplicated test: $test_id" ) if exists $test_info{$test_id};
		$test_info{$test_id} = {
			rank   => $rank,
			cmd    => '',
		}
	}
    else{
        $tb->die( "Not test defined yet: $_" ) unless $test_id;
        while( /\$(\w+)/g ){
            $tb->die( "Environment variable '$1' not set in line: $_" ) unless exists $ENV{$1};
        }
        $test_info{$test_id}{cmd} .= $_;
	}
}

foreach my $test_id ( sort { $test_info{$a}{rank} <=> $test_info{$b}{rank} } keys %test_info ){
    unless( $test_id =~ /$opt{t}/o ){
        $tb->report( 'skip', $test_id );
        next;
    }
    my @cmd = (
        'TEST_NAME=' . $test_id,
        'SOFT='  . '../..',
        'TMP='   . $ENV{TMP},
        'MNET='   . $ENV{MNET},
        'NCORE=' . $ENV{NCORE},
    );
    push @cmd, $ENV{INIT} if $ENV{INIT};
    push @cmd, $test_info{$test_id}{cmd};
    my $t0 = time();
    warn "\n";
    $tb->system( join "\n", @cmd );
    $test_info{$test_id}{elapsed} = time() - $t0;
    $tb->report( $test_id, $test_info{$test_id}{elapsed} . ' s');
}

warn "\n";
$tb->report( '', '### Test Summary ###' );
foreach my $test_id ( sort { $test_info{$a}{rank} <=> $test_info{$b}{rank} } keys %test_info ){
    $tb->report( $test_id, $test_info{$test_id}{elapsed} . ' s');
}

