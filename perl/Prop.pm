package Prop;

use strict;

sub new{
    my( $package, $prop_id, $desc, $type, $data_ref ) = @_;
    my $self = bless {
        id       => $prop_id,   # /^[\w\-\.]+$/
        desc     => $desc,      # any string
        type     => $type,      # one of MetNet type: 'react', 'spec', 'pept' ...
        data     => $data_ref,  # a hash table where
                                #    keys are the internal accession number of react, spec, pept ... in MetNet
                                #    values can be anything
        is_num   => 0,          # not numeric by default
    }, $package;
    return $self;
}
sub get_id{
    shift->{id};
}
sub get_desc{
    shift->{desc};
}
sub get_type{
    shift->{type};
}
sub get_data_ref{
    shift->{data};
}
sub is_num{
    shift->{is_num};
}

sub get_counts{
    my( $self ) = @_;
    die "Abstract method called: Prop::get_counts\n";
}
sub select_items{
    my( $self, $condition) = @_;
    die "Abstract method called: Prop::get_select_items\n";
}
sub get_html_table{ # This is a placeholder sub, only useful for the metaNetX Web site
    my( $self, $mnx, $mnet_id, $analid, $checkbox ) = @_;
    die "Abstract method called: Prop::get_html_table()\n";
}

1;

