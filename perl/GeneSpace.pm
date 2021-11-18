package GeneSpace;

use strict;

use Exporter;
our @ISA = qw( Exporter );
use List::Util qw( uniq );
use Data::Dumper;

use Toolbox;
my $tb = Toolbox->new();

sub new{
    my( $package, $pept_dir, $xref_filter ) = @_;
    my $self = bless {
        id2pept       => {},
        pept_info     => {},
        filter_regexp => $xref_filter ? qr/$xref_filter/ : qr/./,
    }, $package;
    $self->add_pept_dir( $pept_dir ) if $pept_dir;
    return $self;
}
sub add_pept_dir{
    my( $self, $pept_dir ) = @_;
    my $tsv = $tb->scan_tsv( $pept_dir . '/id_map.tsv' );
    $_->[1] =~ s/^(arch|bact|euk):// foreach @$tsv;  # FIXME: this is a temporary patch
    $self->{id2pept}{$_->[0]} = $_->[1]  foreach @$tsv;
    $tsv = $tb->scan_tsv(  $pept_dir . '/peptides.tsv' );
    foreach( @$tsv ){ # FIXME: these are temporary patches
        $_->[0] =~ s/^(arch|bact|euk)://;
        $_->[2] =~ s/uniprot:/uniprotkb:/;
        $_->[3] =~ s/ /;/g;
        $self->{pept_info}{$_->[0]} = [ $_->[1], $_->[2], $_->[3] ];
    }
}
sub map_pept{
    my( $self, $metnet, $pept_id ) = @_;
    return $self->{id2pept}{$pept_id} if exists $self->{id2pept}{$pept_id};
    return $pept_id if $pept_id eq 'SPONTANEOUS';
    my %buf = ();
    foreach my $id ( sort split /;/, $metnet->get_pept_xrefs( $pept_id )){
        next unless $id =~ $self->{filter_regexp};
        if( exists $self->{id2pept}{$id} ){
            push @{$buf{$self->{id2pept}{$id}}}, $id;
        }
    }
    unless( keys %buf ){
        foreach my $id ( sort split /;/, $metnet->get_pept_xrefs( $pept_id )){
            next unless $id =~ $self->{filter_regexp};
            if( $id =~ /(.+):(.+)/ ){
                next if $1 eq 'refseq_synonym'; # FIXME: don't hardcode stuff like this
                next if $2 =~ /^\d+$/;          #        really
                push @{$buf{$self->{id2pept}{$2}}}, $2 if exists $self->{id2pept}{$2};
            }
        }
    }
    unless( keys %buf ){
        foreach my $id ( sort split /;/, $metnet->get_pept_gene( $pept_id )){
            next unless $id =~ $self->{filter_regexp};
            push @{$buf{$self->{id2pept}{$id}}}, $id if exists $self->{id2pept}{$id};
        }
    }
    my @id = keys %buf;
    if( @id > 1){
        $tb->warn( "Multiple mapping for $pept_id: " . Dumper \%buf );
    }
    else{
        return $id[0] if @id;
    }
    return $pept_id;
#
    return 'UNK:' . $pept_id;

#    if( $id =~ /^gene:(.+)/ ){ # FIXME: this is a temporary patch
#        return $self->{id2pept}{$1}  if exists $self->{id2pept}{$1};
#    }
#    if( $id =~ /^gene:G_(.+)/ ){ # mmmh stuff? # FIXME: this is a temporary patch
#        return $self->{id2pept}{$1}  if exists $self->{id2pept}{$1};
#    }
#    return $i;

}
# uniq was added to fix some problem in source SBML
sub map_complex{
    my( $self, $metnet, $complex ) = @_;
    my @prot = ();
    foreach my $prot ( split /;/, $complex ){
        my @pept = ();
        foreach( split /\+/, $prot ){
            if( $_ ){
                push @pept, $self->map_pept( $metnet, $_ );
            }
            else{
                $tb->warn( "Empty pept found in prot: $prot" ) unless $_;
            }
        }
        if( @pept != uniq @pept ){
            $tb->warn( "Pept not unique ($complex): " . join ' ', @pept );
            @pept = uniq @pept;
        }
        push @prot, join '+', sort @pept;
    }
    if( @prot != uniq @prot ){
        $tb->warn( "Prot not unique ($complex): " . join ' ', @prot );
        @prot = uniq @prot;
    }
    return join ';', sort @prot;
}
sub exists_pept{
    my( $self, $pept_id ) = @_;
    return exists $self->{pept_info}{$pept_id};
}
sub get_pept_info{
    my( $self, $pept_id ) = @_;
    return @{$self->{pept_info}{$pept_id}};
}

1;

