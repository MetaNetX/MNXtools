package LinkGroup;

use strict;

use Data::Dumper;

sub cluster{
    my( $double_hash_ref ) = @_; # any data structure that "starts" with a double HASH
    # create a bidirectional data structure:
    my %ref = ();
    foreach my $key ( keys %$double_hash_ref ){
        foreach( keys %{$double_hash_ref->{$key}} ){
            $ref{$key}{$_} = 1;
            $ref{$_}{$key} = 1;
        }
    }
    # collect clusters
    my %done    = ();
    my @cluster = ();
    foreach my $key ( sort keys %ref ){
        next if $done{$key};
        my $grp = { $key => 1 };
        _collect( \%ref, $grp );
        my %buf = ();
        foreach( keys %$grp ){
            $done{$_} = 1;
            $buf{$_} = $double_hash_ref->{$_} if exists $double_hash_ref->{$_};
        }
        push @cluster, \%buf;
    }
    return @cluster;
}

sub _collect{
    my( $ref, $grp ) = @_;
    my $done = 1;
    foreach my $key ( keys %{$grp} ){
        foreach( keys %{$ref->{$key}} ){
            if( exists $grp->{$_} ){
                next;
            }
            else{
                $grp->{$_} = 1;
                $done      = 0;
            }
        }
    }
    _collect( $ref, $grp ) unless $done; # recursion
}

# https://en.wikipedia.org/wiki/Transitive_reduction
# https://stackoverflow.com/questions/1690953/transitive-reduction-algorithm-pseudocode
# $list can be any structure that "starts" with a double HASH
# $list is modified by this routine
# Should once implement the oposite algo: transitive closure

# The basic gist of the transitive reduction algorithm I used is

# foreach x in graph.vertices
#    foreach y in graph.vertices
#      foreach z in graph.vertices
#         delete edge xz if edges xy and yz exist

# The transitive closure algorithm I used in the same script is very similar but the last line is
#         add edge xz if edges xy and yz OR edge xz exist

# The implementation below should be much faster on large graph
# For example: head => y must be removed
#        {
#            head => {
#                x => 1,
#                y => 1,
#            },
#            x => {
#                y => 1,
#            },
#        };

sub transitive_reduction{
    my( $list ) = @_;
    my $done = 0;
    while( ! $done ){ # is this loop necessary?
        $done = 1;
        foreach my $head ( keys %$list ){
            my @child = sort keys %{$list->{$head}};
            next if @child == 1; # skip
            foreach my $x ( 0 .. @child - 2 ){ # @child > 1
                foreach my $y ( $x + 1 .. @child - 1 ){
                    if( exists $list->{$child[$x]} and exists $list->{$child[$x]}{$child[$y]} ){
                        delete $list->{$head}{$child[$y]};
                        $done = 0;
                    }
                    elsif( exists $list->{$child[$y]} and exists $list->{$child[$y]}{$child[$x]} ){
                        delete $list->{$head}{$child[$x]};
                        $done = 0;
                    }
                }
            }
        }
    }
}
# foreach x in graph.vertices
# #    foreach y in graph.vertices
# #      foreach z in graph.vertices
# #         delete edge xz if edges xy and yz exist
#
# # The transitive closure algorithm I used in the same script is very similar but the last line is
# #         add edge xz if edges xy and yz OR edge xz exist
sub transitive_closure{
    my( $list ) = @_;
    my @node = all_nodes( $list );
    foreach my $x ( @node ){
        foreach my $y ( @node ){
            foreach my $z ( @node ){
                if( exists $list->{$x} and exists  $list->{$x}{$y} and
                    exists $list->{$y} and exists  $list->{$y}{$z} ){
                    $list->{$x}{$z} = 1;
                }
            }
        }
    }
}

sub has_linear_topology{ # distinct "non-connexe" graphs can pass the this test
    my $g = shift;
    my %rev = ();
    foreach my $head ( keys %$g ){
        return 0 if keys %{$g->{$head}} > 1;
        $rev{$_}{$head} = 1 foreach keys %{$g->{$head}};
    }
    foreach my $tail ( keys %rev ){
        return 0 if keys %{$rev{$tail}} > 1;
    }
    return 1;
}

sub all_nodes{
    my $g = shift;
    my %node = ();
    foreach my $head ( keys %$g ){
        $node{$head} = 1;
        $node{$_}    = 1 foreach keys %{$g->{$head}};
    }
    return keys %node;
}

1;

