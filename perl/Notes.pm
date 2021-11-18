package Notes;
#Filename Notes.pm

use strict;
use warnings;
use diagnostics;

use FindBin qw( $RealBin );
use lib "$RealBin/../perl";
use Formaters;

sub parse_notes {
    my ($element, $element_info) = @_;

    for my $note ( grep { !/<[A-Za-z]+/ } map { s/<\/(html:)?p>.*$//s; $_ } grep { /<\/(html:)?p>/ } split(/<(html:p|p)>/, $element->getNotesString(), -1) ){
        if ( $note =~ /^((.+?):.+)$/ ){
            push @{ $element_info->{$2} }, Formaters::clean_Note($1);
        }
    }

    return $element_info;
}

1;

