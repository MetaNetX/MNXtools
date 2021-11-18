package SBML;
#Filename SBML.pm

use strict;
use warnings;
use diagnostics;

use FindBin qw( $RealBin );
use lib "$RealBin/../perl";
use Constants;


sub check_SBML {
    my ($SBML_object, $check_SBML, $validate_SBML) = @_;

    if ( $check_SBML ){
        if ( $validate_SBML ){
            $SBML_object->checkConsistency();
        }

        my @SBML_errors   = grep { / \[Error\]\)/ }   split(/\n\n/, $SBML_object->getErrorLog()->toString() );
        my @SBML_warnings = grep { / \[Warning\]\)/ } split(/\n\n/, $SBML_object->getErrorLog()->toString() );
        push @SBML_errors,  grep { / \[Fatal\]\)/ }   split(/\n\n/, $SBML_object->getErrorLog()->toString() ); #Fatal are super errors!

        my $SBML_errors   = scalar @SBML_errors   // 0;
        my $SBML_warnings = scalar @SBML_warnings // 0;

        if ( $SBML_warnings > 0 && Constants::get_debug() ){
            warn join("\n\n", @SBML_warnings), "\n\n";
        }
        if ( $SBML_errors > 0 ){
            warn join("\n\n", @SBML_errors), "\n\n";
            $Constants::logger->fatal("\t$SBML_errors error(s), $SBML_warnings warning(s)");
            exit 1;
        }
    }

    return;
}

1;

