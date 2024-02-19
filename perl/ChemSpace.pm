package ChemSpace;

# ------------------------------------------------ #
# A helper class to help compute new version of
# MNXref and to convert models. It is optimized for
# readability and not really for speed!
# ------------------------------------------------ #

use strict;

use Sort::Naturally;
use Data::Dumper;
use Storable;
use Digest::CRC;
use FindBin qw( $RealBin );

use lib "$RealBin/../util";
use MetNet;

# ------------------------------------------------ #
# Constructors
# ------------------------------------------------ #

sub new{
    my( $package, $tb ) = @_;
    my $self = {
        tb                => $tb,
        option            => { PMF => 1 },
        chem_prop         => {},
        chem_xref         => {},
        chem_depr         => {},
        chem_late         => {},
        chem_isom         => {},
        id_to_chem        => {},
        reac_prop         => {},
        reac_xref         => {},
        reac_depr         => {},
        comp_prop         => {},
        comp_xref         => {},
        comp_depr         => {},
        id_to_comp        => {},
        eq_to_mnxr        => {},
    };
    return bless $self, $package;
}
sub load{
    my( $self, $chemistry_dir ) = @_;
    $self->_load_chem_prop( $chemistry_dir . '/chem_prop.tsv' );
    $self->_load_chem_xref( $chemistry_dir . '/chem_xref.tsv' );
    $self->_load_chem_depr( $chemistry_dir . '/chem_depr.tsv' );
    $self->_load_chem_isom( $chemistry_dir . '/chem_isom.tsv' );
    $self->_load_comp_prop( $chemistry_dir . '/comp_prop.tsv' );
    $self->_load_comp_xref( $chemistry_dir . '/comp_xref.tsv' );
    $self->_load_comp_depr( $chemistry_dir . '/comp_depr.tsv' );
    $self->_load_reac_prop( $chemistry_dir . '/reac_prop.tsv' );
    $self->_load_reac_xref( $chemistry_dir . '/reac_xref.tsv' );
    $self->_load_reac_depr( $chemistry_dir . '/reac_depr.tsv' );
    $self->_load_reac_isom( $chemistry_dir . '/reac_isom.tsv' );
    return $self;
}

# ------------------------------------------------ #
# Save itself to disk
# ------------------------------------------------ #

sub write{
    my( $self, $chemistry_dir, $header ) = @_;
    $self->_write_chem_prop( $chemistry_dir . '/chem_prop.tsv', $header );
    $self->_write_chem_xref( $chemistry_dir . '/chem_xref.tsv', $header );
    $self->_write_chem_depr( $chemistry_dir . '/chem_depr.tsv', $header );
    # $self->_write_chem_late( $chemistry_dir . '/Compounds_LateMerge.tsv',     $header );
    $self->_write_chem_isom( $chemistry_dir . '/chem_isom.tsv', $header );
    $self->_write_comp_prop( $chemistry_dir . '/comp_prop.tsv', $header );
    $self->_write_comp_xref( $chemistry_dir . '/comp_xref.tsv', $header );
    $self->_write_comp_depr( $chemistry_dir . '/comp_depr.tsv', $header );
    $self->_write_reac_prop( $chemistry_dir . '/reac_prop.tsv', $header );
    $self->_write_reac_xref( $chemistry_dir . '/reac_xref.tsv', $header );
    $self->_write_reac_depr( $chemistry_dir . '/reac_depr.tsv', $header );
    $self->_write_reac_isom( $chemistry_dir . '/reac_isom.tsv', $header );
}

# ------------------------------------------------ #
# The two subs below to be removed ?
# ------------------------------------------------ #

sub enable_late_merge{ # Cannot be undone, currently
    my $self = shift;
    my %replace = ();
    foreach my $id_1 ( keys %{$self->{chem_late}} ){
        my @buf = keys %{ $self->{id_to_chem}{$id_1}};
        $self->{tb}->die( "Missing or multiple mapping for chem ID: $id_1 (" .join( ',', @buf ) . ')' ) unless @buf == 1;
        my $mnxm_1 = $buf[0];
        my $id_2 = $self->{chem_late}{$id_1};
        @buf = keys %{ $self->{id_to_chem}{$id_2}};
        $self->{tb}->die( "Missing or multiple mapping for chem ID: $id_2 (" .join( ',', @buf ) . ')' ) unless @buf == 1;
        my $mnxm_2 = $buf[0];
        warn "$mnxm_1 <- $mnxm_2 ( $id_1 <- $id_2)\n";
        $replace{$mnxm_1} = $mnxm_2;
    }
    $self->{id_to_chem} = {};
    foreach my $mnxm_old ( keys %{$self->{chem_xref}} ){
        my $mnxm_new = exists $replace{$mnxm_old} ? $replace{$mnxm_old} : $mnxm_old ;
        $self->{id_to_chem}{$mnxm_old}{$mnxm_new} = 1;
        foreach my $source ( keys %{$self->{chem_xref}{$mnxm_old}} ){
            $self->{id_to_chem}{$source}{$mnxm_new} = 1;
            $self->{id_to_chem}{$1}{$mnxm_new}      = 1 if $source =~ /:(\S+)/;
        }
    }
}
sub enable_un_merge{
    my $self = shift;
}
sub disable_PMF{
    my $self = shift;
    $self->{option}{PMF} = 0;
}
sub disable_normalize_coef{
    my $self = shift;
    no strict;
    *{_normalize_coef} = sub{ return shift };
    use strict;
}
sub add_chem_alt_id{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_txt( $file, 2 );
    foreach( @$tsv ){
        my( $mnxm_id, $msg ) = $self->map_chem( $_->[1] );
        $self->{id_to_chem}{$_->[0]}{$mnxm_id} = 1;
        warn "*** $_->[0] ===> $mnxm_id *** \n";
    }
}

# ------------------------------------------------------------
# store/retrieve itself to/from cache file
# ------------------------------------------------------------

sub retrieve_from_file{ # a constructor
    my( $self, $cache_file ) = @_;
    my $tb = $self->{tb};
    $self = Storable::retrieve( $cache_file ) or $self->{tb}->confess( "Cannot retrieve from cache file: $cache_file" );
    $self->{tb} = $tb;
    return $self;
}
sub store_to_file{ # Nota bene: nstore explodes memory consumption (FIXME: verify)!
    my( $self, $cache_file ) = @_;
    my $tb = $self->{tb};
    $self->{tb} = undef;
    Storable::store( $self, $cache_file )  or $self->{tb}->confess( "Cannot store to cache file: $cache_file" );
    $self->{tb} = $tb;
}

# ------------------------------------------------ #
# Chem prop related routines
# ------------------------------------------------ #

my @chem_prop_key = (
    'name',
    'reference',
    'formula',
    'charge',
    'mass',
    'InChI',
    'InChIKey',
    'SMILES',
);
my %chem_prop_default  = (); # set default value
$chem_prop_default{$_} = '' foreach @chem_prop_key;

sub add_chem_prop{
    my( $self, $ID, %arg ) = @_;
    $self->{chem_prop}{$ID} = {};
    $self->set_chem_prop( $ID, %arg ); # update default
    $self->{id_to_chem}{$ID}{$ID} = 1;
}
sub set_chem_prop{
    my( $self, $ID, %arg ) = @_;
    $self->{tb}->confess( "Unknown chem ID: $ID" ) unless exists $self->{chem_prop}{$ID};
    foreach( keys %arg ){
        $self->{tb}->confess( "Invalid prop key passed: ID=$ID key=$_ value=$arg{$_}" ) unless exists $chem_prop_default{$_};
        $self->{chem_prop}{$ID}{$_} = $arg{$_};
    }
}
sub get_chem_prop{
    my( $self, $ID ) = @_;
    $self->{tb}->confess( "Unknown chem ID: $ID" ) unless exists $self->{chem_prop}{$ID};
    return %{$self->{chem_prop}{$ID}};
}
sub _load_chem_prop{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->add_chem_prop(
            $_->[0],
            name       => $_->[1],
            reference  => $_->[2],
            formula    => $_->[3],
            charge     => $_->[4],
            mass       => $_->[5],
            InChI      => $_->[6],
            InChIKey   => $_->[7],
            SMILES     => $_->[8],
        );
    }
}
sub _write_chem_prop{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#ID', @chem_prop_key ), "\n";
    foreach my $id ( sort keys %{$self->{chem_prop}} ){
        my @token = ( $id );
        foreach( @chem_prop_key ){
            push @token, $self->{chem_prop}{$id}{$_};
        }
        print $tsv join( "\t", @token ) . "\n";
    }
    close $tsv;
}

sub get_chem_ids{
    my $self = shift;
    return sort keys %{$self->{chem_prop}};
}

# ------------------------------------------------ #
# Chem xref related routines
# ------------------------------------------------ #

my @chem_xref_key = (
    'description',
);
my %chem_xref_default  = (); # set default value
$chem_xref_default{$_} = '' foreach @chem_xref_key;

sub add_chem_xref{
    my( $self, $ID, $source, %arg ) = @_;
    $self->{chem_xref}{$ID}{$source} = {};
    $self->set_chem_xref( $ID, $source, %arg ); # update default
    $self->{id_to_chem}{$source}{$ID} = 1;
    $self->{id_to_chem}{$1}{$ID}      = 1 if $source =~ /:(\S+)/;
}
sub set_chem_xref{
    my( $self, $ID, $source, %arg ) = @_;
    $self->{tb}->confess( "Unknown chem ID: $ID" ) unless exists $self->{chem_xref}{$ID};
    foreach( keys %arg ){
        $self->{tb}->confess( "Invalid key passed: ID=$ID source=$source key=$_ value=$arg{$_}" ) unless exists $chem_xref_default{$_};
        $self->{chem_xref}{$ID}{$source}{$_} = $arg{$_};
    }
}
sub get_chem_xref{
    my( $self, $ID ) = @_;
    $self->{tb}->confess( "Unknown chem ID: $ID" ) unless exists $self->{chem_xref}{$ID};
    return $self->{chem_xref}{$ID};
}
sub _load_chem_xref{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->add_chem_xref(
            $_->[1],
            $_->[0],
            description => $_->[2],
        );
    }
}
sub _write_chem_xref{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#source', 'ID', @chem_xref_key ), "\n";
    foreach my $id ( sort keys %{$self->{chem_xref}} ){
        foreach my $source ( sort keys %{$self->{chem_xref}{$id}} ){
            my @token = ( $source, $id );
            foreach( @chem_xref_key ){
                push @token, $self->{chem_xref}{$id}{$source}{$_};
            }
            print $tsv join( "\t", @token), "\n";
        }
    }
    close $tsv;
}

# ------------------------------------------------ #
# dealing with deprecated identifiers
# ------------------------------------------------ #

sub add_chem_depr{
    my( $self, $old_mnx_id, $new_mnx_id, $label ) = @_;
    $self->{chem_depr}{$old_mnx_id}{$new_mnx_id} = $label;
}
sub add_reac_depr{
    my( $self, $old_mnx_id, $new_mnx_id, $label ) = @_;
    $self->{reac_depr}{$old_mnx_id}{$new_mnx_id} = $label;
}
sub add_comp_depr{
    my( $self, $old_mnx_id, $new_mnx_id, $label ) = @_;
    $self->{comp_depr}{$old_mnx_id}{$new_mnx_id} = $label;
}
sub _write_chem_depr{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#deprecated_ID', 'ID', 'version' ), "\n";
    foreach my $old_id ( sort keys %{$self->{chem_depr}} ){
        foreach ( sort keys %{$self->{chem_depr}{$old_id}} ){
            print $tsv join "\t", $old_id, $_, $self->{chem_depr}{$old_id}{$_} . "\n";
        }
    }
    close $tsv;
}
sub _write_reac_depr{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#deprecated_ID', 'ID', 'version' ), "\n";
    foreach my $old_id ( sort keys %{$self->{reac_depr}} ){
        foreach ( sort keys %{$self->{reac_depr}{$old_id}} ){
            print $tsv join "\t", $old_id, $_, $self->{reac_depr}{$old_id}{$_} . "\n";
        }
    }
    close $tsv;
}
sub _write_comp_depr{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#deprecated_ID', 'ID', 'version' ), "\n";
    foreach my $old_id ( sort keys %{$self->{comp_depr}} ){
        foreach ( sort keys %{$self->{comp_depr}{$old_id}} ){
            print $tsv join "\t", $old_id, $_, $self->{comp_depr}{$old_id}{$_} . "\n";
        }
    }
    close $tsv;
}
sub _load_chem_depr{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->{chem_depr}{$_->[0]}{$_->[1]} = $_->[2];
    }
}
sub _load_comp_depr{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->{comp_depr}{$_->[0]}{$_->[1]} = $_->[2];
    }
}
sub _load_reac_depr{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->{reac_depr}{$_->[0]}{$_->[1]} = $_->[2];
    }
}

# ------------------------------------------------ #
# chem_late is still experimental
# ------------------------------------------------ #

sub _load_chem_late{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->{chem_late}{$_->[1]} = $_->[0];
    }
}
sub set_chem_late{
    my( $self, $from, $to ) = @_;
    $self->{chem_late}{$from} = $to;
}
sub _write_chem_late{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#ID1', 'ID2s' ), "\n";
    foreach( sort keys %{$self->{chem_late}} ){
        print $tsv join "\t", $_, $self->{chem_late}{$_} . "\n";
    }
    close $tsv;
}
# ------------------------------------------------ #
#
# ------------------------------------------------ #

sub _load_chem_isom{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        push @{$self->{chem_isom}{$_->[1]}}, $_->[0];
    }
}
sub set_chem_isom{
    my( $self, $from, $to ) = @_;
    push @{$self->{chem_isom}{$from}}, $to;
}
sub get_chem_isom{
    my( $self, $ID ) = @_;
    $self->{tb}->confess( "Unknown chem ID: $ID" ) unless exists $self->{chem_prop}{$ID};
    if( exists $self->{chem_isom}{$ID} ){
        return @{$self->{chem_isom}{$ID}};
    }
    else{
        return ();
    }
}
sub _write_chem_isom{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#parent', 'child', 'description'), "\n";
    foreach my $id ( sort keys %{$self->{chem_isom}} ){
        foreach( @{$self->{chem_isom}{$id}} ){
            my $comment = # $self->{chem_prop{$id}{name}
            print $tsv join "\t",
            $id,
            $_,
            $self->{chem_prop}{$id}{name} . ' -> ' . $self->{chem_prop}{$_}{name} . "\n";
        }
    }
    close $tsv;
}

# ------------------------------------------------ #
# Comp prop related routines
# ------------------------------------------------ #

my @comp_prop_key = (
        'name',
        'reference',
#         'comment',
);
my %comp_prop_default  = (); # set default value
$comp_prop_default{$_} = '' foreach @comp_prop_key;

sub add_comp_prop{
    my( $self, $ID, %arg ) = @_;
    $self->{comp_prop}{$ID} = {};
    $self->set_comp_prop( $ID, %arg ); # update default
    $self->{id_to_comp}{$ID}{$ID}     = 1;
}
sub set_comp_prop{
    my( $self, $ID, %arg ) = @_;
    $self->{tb}->confess( "Unknown comp ID: $ID" ) unless exists $self->{comp_prop}{$ID};
    foreach( keys %arg ){
        $self->{tb}->confess( "Invalid prop key passed: ID=$ID key=$_ value=$arg{$_}" ) unless exists $comp_prop_default{$_};
        $self->{comp_prop}{$ID}{$_} = $arg{$_};
    }
}
sub get_comp_prop{
    my( $self, $ID ) = @_;
    $self->{tb}->confess( "Unknown comp ID: $ID" ) unless exists $self->{comp_prop}{$ID};
    return %{$self->{comp_prop}{$ID}};
}
sub _load_comp_prop{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->add_comp_prop(
            $_->[0],
            name       => $_->[1],
            reference  => $_->[2],
        );
    }
}
sub _write_comp_prop{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#ID', @comp_prop_key ), "\n";
    foreach my $id ( sort keys %{$self->{comp_prop}} ){
        print $tsv $id;
        foreach( @comp_prop_key ){
            print $tsv "\t" . $self->{comp_prop}{$id}{$_};
        }
        print $tsv "\n";
    }
    close $tsv;
}
sub get_comp_ids{
    my $self = shift;
    return sort keys %{$self->{comp_prop}};
}

# ------------------------------------------------ #
# Comp xref related routines
# ------------------------------------------------ #

my @comp_xref_key = (
#    'name',
#    'desc',
    'description',
);
my %comp_xref_default  = (); # set default value
$comp_xref_default{$_} = '' foreach @comp_xref_key;

sub add_comp_xref{
    my( $self, $ID, $source, %arg ) = @_;
    $self->{comp_xref}{$ID}{$source} = {};
    $self->set_comp_xref( $ID, $source, %arg ); # update default
    $self->{id_to_comp}{$source}{$ID} = 1;
    $self->{id_to_comp}{$1}{$ID}      = 1 if $source =~ /:(\S+)/;
}
sub set_comp_xref{
    my( $self, $ID, $source, %arg ) = @_;
    $self->{tb}->confess( "Unknown comp ID: $ID" ) unless exists $self->{comp_xref}{$ID};
    foreach( keys %arg ){
        $self->{tb}->confess( "Invalid key passed: ID=$ID source=$source key=$_ value=$arg{$_}" ) unless exists $comp_xref_default{$_};
        $self->{comp_xref}{$ID}{$source}{$_} = $arg{$_};
    }
}
sub get_comp_xref{
    my( $self, $ID ) = @_;
    $self->{tb}->confess( "Unknown comp ID: $ID" ) unless exists $self->{comp_xref}{$ID};
    return $self->{comp_xref}{$ID};
}
sub _load_comp_xref{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->add_comp_xref(
            $_->[1],
            $_->[0],
            description   => $_->[2],
        );
    }
}
sub _write_comp_xref{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#source', 'ID', @comp_xref_key ), "\n";
    foreach my $id ( sort keys %{$self->{comp_xref}} ){
        foreach my $source ( sort keys %{$self->{comp_xref}{$id}} ){
            my @token = ( $source, $id );
            foreach( @comp_xref_key ){
                push @token, $self->{comp_xref}{$id}{$source}{$_};
            }
            print $tsv join( "\t", @token), "\n";
        }
    }
    close $tsv;
}

# ------------------------------------------------ #
# Reac prop related routines
# ------------------------------------------------ #

my @reac_prop_key = (
    'mnx_equation',
    'reference',
    'classifs',
    'is_balanced',
    'is_transport',
);
my %reac_prop_default  = (); # set default value
$reac_prop_default{$_} = '' foreach @reac_prop_key;

sub add_reac_prop{
    my( $self, $ID, %arg ) = @_;
    $self->{reac_prop}{$ID} = {};
    $self->set_reac_prop( $ID, %arg ); # update default
}
sub set_reac_prop{
    my( $self, $ID, %arg ) = @_;
    $self->{tb}->confess( "Unknown reac ID: $ID" ) unless exists $self->{reac_prop}{$ID};
    foreach( keys %arg ){
        $self->{tb}->confess( "Invalid prop key passed: ID=$ID key=$_ value=$arg{$_}" ) unless exists $reac_prop_default{$_};
        $self->{reac_prop}{$ID}{$_} = $arg{$_};
    }
    if( exists $arg{mnx_equation} ){
        $self->{eq_to_mnxr}{$arg{mnx_equation}} = $ID;
    }
}
sub get_reac_prop{
    my( $self, $ID ) = @_;
    $self->{tb}->confess( "Unknown reac ID: $ID" ) unless exists $self->{reac_prop}{$ID};
    return %{$self->{reac_prop}{$ID}};
}
sub get_all_reac_IDs{
    my( $self ) =  @_;
    return keys %{$self->{reac_prop}};
}
sub _load_reac_prop{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->add_reac_prop(
            $_->[0],
            mnx_equation            => $_->[1],
            reference               => $_->[2],
            classifs                => $_->[3],
            is_balanced             => $_->[4],
            is_transport            => $_->[5],
        );
    }
}
sub _write_reac_prop{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#ID', @reac_prop_key ), "\n";
    foreach my $id ( sort keys %{$self->{reac_prop}} ){
        print $tsv $id;
        foreach( @reac_prop_key ){
            print $tsv "\t" . $self->{reac_prop}{$id}{$_};
        }
        print $tsv "\n";
    }
    close $tsv;
}
sub get_reac_ids{
    my $self = shift;
    return sort keys %{$self->{reac_prop}};
}

# ------------------------------------------------ #
# Reac xref related routines
# ------------------------------------------------ #

my @reac_xref_key = (
#    'Source_equation',
#    'Description',
    'description',
);
my %reac_xref_default  = (); # set default value
$reac_xref_default{$_} = '' foreach @reac_xref_key;

sub add_reac_xref{
    my( $self, $ID, $source, %arg ) = @_;
    $self->{reac_xref}{$ID}{$source} = {};
    $self->set_reac_xref( $ID, $source, %arg ); # update default
#    $self->{prefix_id_to_reac}{$source} = $ID;
#    $source =~ /:(\S+)/;
#    push @{$self->{id_to_reac}{$1}}, $ID;
    $self->{id_to_reac}{$source}{$ID} = 1;
    $self->{id_to_reac}{$1}{$ID}      = 1 if $source =~ /:(\S+)/;
}
sub set_reac_xref{
    my( $self, $ID, $source, %arg ) = @_;
    $self->{tb}->confess( "Unknown reac ID: $ID" ) unless exists $self->{reac_xref}{$ID};
    foreach( keys %arg ){
        $self->{tb}->confess( "Invalid key passed: ID=$ID source=$source key=$_ value=$arg{$_}" ) unless exists $reac_xref_default{$_};
        $self->{reac_xref}{$ID}{$source}{$_} = $arg{$_};
    }
}
sub get_reac_xref{
    my( $self, $ID ) = @_;
    $self->{tb}->confess( "Unknown reac ID: $ID" ) unless exists $self->{reac_xref}{$ID};
    return $self->{reac_xref}{$ID};
}
sub _load_reac_xref{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        $self->add_reac_xref(
            $_->[1],
            $_->[0],
            description => $_->[2],
        );
    }
}
sub _write_reac_xref{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#source', 'ID', @reac_xref_key ), "\n";
    foreach my $id ( sort keys %{$self->{reac_xref}} ){
        foreach my $source ( sort keys %{$self->{reac_xref}{$id}} ){
            print $tsv "$source\t$id";
            foreach( @reac_xref_key ){
                print $tsv "\t" . $self->{reac_xref}{$id}{$source}{$_};
            }
            print $tsv "\n";
        }
    }
    close $tsv;
}
# ------------------------------------------------ #
# Reac isom related routines
# ------------------------------------------------ #

sub _load_reac_isom{
    my( $self, $file ) = @_;
    my $tsv = $self->{tb}->scan_tsv( $file );
    foreach( @$tsv ){
        push @{$self->{reac_isom}{$_->[1]}}, $_->[0];
    }
}
sub set_reac_isom{
    my( $self, $from, $to ) = @_;
    push @{$self->{reac_isom}{$from}}, $to;
}
sub get_reac_isom{
    my( $self, $ID ) = @_;
    $self->{tb}->confess( "Unknown reac ID: $ID" ) unless exists $self->{reac_prop}{$ID};
    if( exists $self->{reac_isom}{$ID} ){
        return @{$self->{reac_isom}{$ID}};
    }
    else{
        return ();
    }
}
sub _write_reac_isom{
    my( $self, $file, $header ) = @_;
    my $tsv = $self->{tb}->open( '> ' . $file );
    print $tsv $header;
    print $tsv join( "\t", '#parent', 'child' ), "\n";
    foreach my $id ( sort keys %{$self->{reac_isom}} ){
        foreach( @{$self->{reac_isom}{$id}} ){
            print $tsv join "\t",
            $id,
            $_ . "\n";
        }
    }
    close $tsv;
}

# ------------------------------------------------ #
# Equation-related routines
# ------------------------------------------------ #

sub map_equation{
    my( $self, $eq_str, $prefix ) = @_;
    my( $eq, $msg )                           = $self->_map_equation( $eq_str, $prefix );
    my( $is_balanced, $str_spec, $sign_spec ) = $self->_balance_equation( $eq );
    if( $str_spec eq ' = ' ){
        push @$msg, '- code: REAC_MAP_EMPTY_WARN'; # might be changed latter to REAC_MAP_EMPTY_MNXREF in Convert.pm
        return ( '', ' = ', 1, 1, $msg, 'EMPTY', ' = ', 1 );
    }
    # determine reac_id and mnxr_id below
    my @comp_old = keys %$eq;
    my %comp_new = ();
    $comp_new{"MNXD$_"} = 1 foreach 1 .. @comp_old;
    my $comp_map = _build_comp_map( \@comp_old, \%comp_new );
    my @eq_map = ();
    foreach my $label ( sort keys %$comp_map ){
        my %eq_new = ();
        foreach my $comp_old ( @comp_old ){
            my $comp_new = $comp_map->{$label}{$comp_old};
            $eq_new{$comp_new} = $eq->{$comp_old};
        }
        push @eq_map, [ _canonicalize( \%eq_new ), $comp_map->{$label} ];
    }
    # this 'sort' is critical to chose the canonical representation of the equation:
    @eq_map = sort { $a->[0] cmp $b->[0] } @eq_map;
    my $str_best = $eq_map[0][0];
    my %reac_id = ();
    foreach my $num ( 0 .. @eq_map - 1 ){
        next unless $eq_map[$num][0] eq $str_best;
        my( $str_gen, $sign_gen, $comp_map ) = @{$eq_map[$num]};
        next unless $self->{eq_to_mnxr}{$str_gen};
        my $reac_id = $self->{eq_to_mnxr}{$str_gen};
        my $ok = 1;
        foreach( sort { $comp_map->{$a} cmp $comp_map->{$b} } keys %$comp_map ){
            if( $comp_map->{$_} =~ /^MNXD[1-9]$/){
                if( /^BOUNDARY$/ ){
                    $reac_id .= 'B';
                }
                elsif( /^MNXC(\d+)$/ ){
                    $reac_id .= 'C' . $1;
                }
                elsif( /^MNXD(\d+)$/ ){
                    $reac_id .= 'D' . $1;
                }
                else{
                    $ok = 0;
                }
            }
            else{
                $ok = 0;
            }
            next unless $ok;
        }
        $reac_id =~ s/B(C\d+)$/$1B/; # move BOUNDARY to the last position
        $reac_id{ lc $reac_id } = $num;
    }
    my $reac_id = '';
    my $num     = 0;
    if( %reac_id ){
         # this 'sort' is critical to chose the canonical representation of the reaction identifier:
        my @reac_id = nsort keys %reac_id;
        $reac_id = $reac_id[0];
        $num = $reac_id{$reac_id};
    }
    if( ! $reac_id or 16 < length $reac_id or $reac_id eq 'mnxr01' ){
        $reac_id = MetNet::compute_reac_id( $str_spec );
    }
    my( $str_gen, $sign_gen, $comp_map ) = @{$eq_map[$num]};
    my $mnxr_id = $self->{eq_to_mnxr}{$str_gen} || '';
    return ( $reac_id, $str_spec, $sign_spec, $is_balanced, $msg, $mnxr_id, $str_gen, $sign_gen );
}

sub _map_equation{
    my( $self, $str, $prefix ) = @_;
    return( {}, [] ) if $str eq ' = '; # Empty equation
    my %msg = ();
    my @warn = ();
    my $eq_old = $self->_normalize_coef( _parse_equation( $str ));
#    if( exists $eq_old->{BOUNDARY} ){ # check (restricted) syntax of boundary reaction # FIXME: is this the right place?
#        my @comp_id = keys %$eq_old;
#        if( @comp_id != 2 ){
#            $msg{ERROR}{'Compartment salad at BOUNDARY'} = 1;
#        }
#        my @chem_id_1 = keys %{$eq_old->{$comp_id[0]}};
#        my @chem_id_2 = keys %{$eq_old->{$comp_id[1]}};
#        if( @chem_id_1 != 1 or @chem_id_2 != 1 or $chem_id_1[0] ne $chem_id_2[0] ){
#            $msg{ERROR}{'Chemical salad at BOUNDARY'} = 1;
#        }
#    }
    my %eq_new    = ();
    my %proton    = ();
    my %hydroxyde = ();
    my %new2old   = ();
    my $has_PMF   = 0;
    foreach my $comp_old ( keys %$eq_old ){
        my( $comp_new, $msg ) = $self->map_comp( $comp_old, $prefix );
        foreach my $chem_old ( keys %{$eq_old->{$comp_old}} ){
            my( $chem_new, $msg ) = $self->map_chem( $chem_old, $prefix );
            my $txt = $chem_old . ' (' . ( exists $self->{chem_xref}{$chem_old} ? $self->{chem_xref}{$chem_old}{name} : '' ) . ')';
            if( $chem_new eq 'MNXM1' ){
                $proton{$comp_new} += $eq_old->{$comp_old}{$chem_old};
            }
            elsif( $chem_new eq 'MNXM02' ){
                $hydroxyde{$comp_new} += $eq_old->{$comp_old}{$chem_old};
            }
            elsif( $chem_new eq 'MNXM03' ){
                $proton{$comp_new}        += $eq_old->{$comp_old}{$chem_old};
                $eq_new{$comp_new}{WATER} += $eq_old->{$comp_old}{$chem_old}; # add water
            }
            else{
                $has_PMF = 1 if $chem_new eq 'MNXM01';
                $eq_new{$comp_new}{$chem_new} += $eq_old->{$comp_old}{$chem_old};
                $new2old{$chem_new}{$chem_old} = 1;
            }
        }
    }
    foreach my $comp_new ( keys %eq_new ){
        foreach my $chem_new ( keys %{$eq_new{$comp_new}} ){
            if( $eq_new{$comp_new}{$chem_new} == 0 ){            # This is an VERY important test to find problem in the namepace
                delete $eq_new{$comp_new}{$chem_new};            # to permit debugging in caller context a WARNING is issued that
                foreach( sort keys %{$new2old{$chem_new}} ){
                    push @warn,
                         '- code: REAC_MAP_LOSS',
                         "  chem: $chem_new # " . $self->{chem_prop}{$chem_new}{name},
                         "  comp: $comp_new # " . $self->{comp_prop}{$comp_new}{name};
                }
            }
        }
    }

    if( exists $eq_old->{BOUNDARY} or $has_PMF ){
        foreach( keys %proton ){              # assume PMF already properly defined
            $eq_new{$_}{MNXM1} = $proton{$_}; # hence, restore protons
        }
        return ( \%eq_new, \@warn );           # ... and return
    }
    if( %hydroxyde ){ # replace HYDROXYDE with WATER minus PROTON in the right compartment(s)
        foreach my $comp ( keys %hydroxyde ){
            $proton{$comp}        -= $hydroxyde{$comp}; # remove a proton
            $eq_new{$comp}{WATER} += $hydroxyde{$comp}; # and add water
        }
    }
    my @comp = keys %proton; # the compartment(s) where proton are found
    if( @comp == 1 ){
        if( exists $eq_new{$comp[0]} ){ # likely to be a balance proton
            $eq_new{$comp[0]}{MNXM1} += $proton{$comp[0]};
        }
        else{ # this proton is in its own compartment => this is an implicit transport !!!
            my @buf = keys %eq_new;
            if( @buf == 1 ){
                $proton{$buf[0]} = - $proton{$comp[0]}; # balance equation w.r.t. PMF
                @comp = keys %proton; # update it
            }
            else{ # => there is three compartments!
                push @warn, '- code: REAC_MAP_PROTON_SALAD';
            }
        }
    }
    if( @comp == 2 ){ # introduce translocated protons
        my $count = 0; # how many protons are translocated
        if( keys %eq_new == 1 ){ # the whole reaction exepted proton, occurs in a single compartment
            $count = exists $eq_new{$comp[0]}
                   ? abs( $proton{$comp[1]} )
                   : abs( $proton{$comp[0]} );
        }
        elsif( abs( $proton{$comp[0]} ) == abs( $proton{$comp[1]} ) ){
            $count = abs( $proton{$comp[0]} );
        }
        else{ # this start to be more complicated: resolve ambiguity by arbitrily taking the smallest possible number
            $count = abs( $proton{$comp[0]} ) < abs( $proton{$comp[1]} )
                   ? abs( $proton{$comp[0]} )
                   : abs( $proton{$comp[1]} );
            #$msg{WARNING} = "Ambiguous identification of translocated protons";
            # warn "Ambiguous identification of translocated protons: $str\n";
        }
        my $sgn = $proton{$comp[0]} < 0 ? -1 : 1;
        $eq_new{$comp[0]}{MNXM01} =   $sgn * $count;
        $eq_new{$comp[1]}{MNXM01} = - $sgn * $count;
        if( ( abs( $proton{$comp[0]} ) - $count ) > 0 ){
            $eq_new{$comp[0]}{MNXM1} = $sgn * ( abs( $proton{$comp[0]} ) - $count );
        }
        elsif( ( abs( $proton{$comp[1]} ) - $count ) > 0 ){
            $eq_new{$comp[1]}{MNXM1} = - $sgn * ( abs( $proton{$comp[1]} ) - $count ),
        }
    }
    if( @comp > 2 ){
        push @warn, '- code: REAC_MAP_PROTON_SALAD';
    }
    return( \%eq_new, \@warn );
}
sub _chose_chem{
    my( $self, @chem ) = @_;
    my $warn = [];
    if( @chem > 1 ){
        @chem = nsort @chem;
        $warn = [
            '- code: CHEM_MAP_MULTIPLE',
            '  IDs_dst:',
        ];
        push @$warn, "      - $_ # " . $self->{chem_prop}{$_}{name} foreach @chem;
    }
    return( $chem[0], $warn );
}
sub _search_chem_depr{
    my( $self, $id ) = @_;
    my %buf = ();
    if( exists $self->{chem_depr}{$id} ){
        foreach my $id2 ( keys %{$self->{chem_depr}{$id}} ){
            if( exists $self->{chem_prop}{$id2} ){
                $buf{$id2} = 1;
            }
            else{
                $buf{$_} = 1 foreach $self->_search_chem_depr( $id2 ); # reccursion
            }
        }
    }
    return keys %buf;
}
sub map_chem{
    my( $self, $id_old ) = @_;
    return ( $id_old, [ '- code: CHEM_MAP_UNKNOWN' ] ) if $id_old =~ /^UNK:/;
    my $id_new = '';
    my $msg = [];
    if( exists $self->{id_to_chem}{$id_old} ){
        ( $id_new, $msg ) = $self->_chose_chem( keys %{$self->{id_to_chem}{$id_old}} );
        $msg = [ '- code: CHEM_MAP_OK' ] unless @$msg;
    }
    # Numeric identifier as in chebi or pubchem are possibly misleading
    if( ! $id_new and $id_old =~ /:(.+)/ and exists $self->{id_to_chem}{$1} ){
        ( $id_new, $msg ) = $self->_chose_chem( keys %{$self->{id_to_chem}{$1}} );
        $msg = [ '- code: CHEM_MAP_OK' ] unless @$msg;
    }
    if( ! $id_new and my @id  = $self->_search_chem_depr( $id_old ) ){
        ( $id_new, $msg ) = $self->_chose_chem( @id );
        unshift @$msg, '- code: CHEM_MAP_DEPRECATED';
    }
    if( $id_new ){
        return( $id_new, $msg );
    }
    else{
        return ( _get_unk_id( 'chem', $id_old ), [ '- code: CHEM_MAP_UNKNOWN' ] );
    }
}
sub _get_unk_id{
    my( $prefix, $id ) = @_;
    $id =~ s/[^\w\-\#\.]/_/g; # security
    if( length $id > 24 ){
        my $ctx = Digest::CRC->new( type => 'crc32', poly => 0x2D );
        $ctx->add( $id );
        $id = $prefix . uc $ctx->hexdigest();
    }
    return 'UNK:' . $id
}
sub search_chem_isom{ # return the list of all isom children + itself
    my( $self, $mnxm_id ) = @_;
    my %id = ( $mnxm_id => 1 );
    foreach my $id ( $self->get_chem_isom( $mnxm_id )){
        foreach( $self->search_chem_isom( $id )){
            $id{$_} = 1;
        }
    }
    return sort keys %id;
}
sub _chose_comp{
    my( $self, @comp ) = @_;
    my $warn = [];
    if( @comp > 1 ){
        @comp = nsort @comp;
        $warn = [
            '- code: COMP_MAP_MULTIPLE',
            '  IDs_dst:',
        ];
        push @$warn, "      - $_ # " . $self->{comp_prop}{$_}{name}  foreach @comp;
    }
    return ( $comp[0], $warn );
}
sub _search_comp_depr{
    my( $self, $id ) = @_;
    my %buf = ();
    if( exists $self->{comp_depr}{$id} ){
        foreach my $id2 ( keys %{$self->{comp_depr}{$id}} ){
            if( exists $self->{comp_prop}{$id2} ){
                $buf{$id2} = 1;
            }
            else{
                $buf{$_} = 1  foreach $self->_search_comp_depr( $id2 ); # reccursion
            }
        }
    }
    return keys %buf;
}
sub map_comp{ # same procedure as in map_chem
    my( $self, $id_old ) = @_; # , $prefix ) = @_;
    return ( $id_old, [ '- code: COMP_MAP_UNKNOWN' ] )  if $id_old =~ /^UNK:/;
    my $id_new = '';
    my $msg = [];
    if( exists $self->{id_to_comp}{$id_old} ){
        ( $id_new, $msg ) = $self->_chose_comp( keys %{$self->{id_to_comp}{$id_old}} );
        $msg = [ '- code: COMP_MAP_OK' ]  unless @$msg;
    }
    if( ! $id_new and $id_old =~ /:(.+)/ and exists $self->{id_to_comp}{$1} ){
        ( $id_new, $msg ) = $self->_chose_comp( keys %{$self->{id_to_comp}{$1}} );
        $msg = [ '- code: COMP_MAP_OK' ] unless @$msg;
    }
    if( ! $id_new and my @id  = $self->_search_comp_depr( $id_old ) ){
        ( $id_new, $msg ) = $self->_chose_comp( @id );
        unshift @$msg, '- code: COMP_MAP_DEPRECATED';
    }
    if( $id_new ){
        return( $id_new, $msg );
    }
    else{
        return ( _get_unk_id( 'comp', $id_old ), [ '- code: COMP_MAP_UNKNOWN' ]);
    }
}
sub _chose_reac{
    my( $self, @reac ) = @_;
    my $warn = [];
    if( @reac > 1 ){
        @reac = nsort @reac;
        $warn = [
            '- code: REAC_MAP_MULTIPLE',
            '  IDs_dst:',
        ];
        push @$warn, "      - $_ # " . $self->{reac_prop}{$_}{name} foreach @reac;
    }
    return ( $reac[0], $warn );
}
sub _search_reac_depr{
    my( $self, $id ) = @_;
    my %buf = ();
    if( exists $self->{reac_depr}{$id} ){
        foreach my $id2 ( keys %{$self->{reac_depr}{$id}} ){
            if( exists $self->{reac_prop}{$id2} ){
                $buf{$id2} = 1;
            }
            else{
                $buf{$_} = 1 foreach $self->_search_reac_depr( $id2 ); # reccursion
            }
        }
    }
    return keys %buf;
}
sub map_reac{ # this subs is rarely used, prefer to use map_equation() ALWAYS
    my( $self, $id_old ) = @_;
    return ( $id_old, [ '- code: REAC_MAP_UNKNOWN' ] ) if $id_old =~ /^UNK:/;
    my $id_new = '';
    my $msg = [];
    if( exists $self->{id_to_reac}{$id_old} ){
        ( $id_new, $msg ) = $self->_chose_reac( keys %{$self->{id_to_reac}{$id_old}} );
        $msg = [ '- code: REAC_MAP_OK' ] unless @$msg;
    }
    # Numeric identifier as in chebi or pubchem are possibly misleading
    if( ! $id_new and $id_old =~ /:(.+)/ and exists $self->{id_to_reac}{$1} ){
        ( $id_new, $msg ) = $self->_chose_reac( keys %{$self->{id_to_reac}{$1}} );
        $msg = [ '- code: REAC_MAP_OK' ] unless @$msg;
    }
    if( ! $id_new and my @id  = $self->_search_reac_depr( $id_old ) ){
        ( $id_new, $msg ) = $self->_chose_reac( @id );
        unshift @$msg, '- code: REAC_MAP_DEPRECATED';
    }
    if( $id_new ){
        return( $id_new, $msg );
    }
    else{
        return ( _get_unk_id( 'reac', $id_old ), [ '- code: REAC_MAP_UNKNOWN' ] );
    }
}

sub _balance_equation{ # Nota bene: MNXM1 and MNXM01 are treated differently
    my( $self, $ref ) = @_;
    my %equation      = %{ $ref };
    return( 1, _canonicalize( \%equation )) if exists $equation{BOUNDARY}; # save time and permit for MNXM1@BOUNDARY
    my %molec_count   = ();
    my $delta_charge  = 0;
    my %delta_formula = ();
    my %delta_molec   = (); # to count compounds without charge and formula
    foreach my $comp ( keys %equation ){
        foreach my $chem ( keys %{$equation{$comp}} ){
            my $coef =$equation{$comp}{$chem};
            if( $chem eq 'MNXM1' ){
                delete $equation{$comp}{$chem}; # as a consequence, unbalanced equation will have no proton
            }
            elsif( exists $self->{chem_prop}{$chem} and
                   $self->{chem_prop}{$chem}{charge} ne '' and
                   $self->{chem_prop}{$chem}{formula} ne '' and
                   $self->{chem_prop}{$chem}{formula} !~ /\*$/ ){
                $delta_charge += $coef * $self->{chem_prop}{$chem}{charge};
                while( $self->{chem_prop}{$chem}{formula} =~ /([A-Z][a-z]*)(\d*)/g ){
                    $delta_formula{$1} += $coef * ( $2 || 1 );
                }
                $molec_count{$comp} += abs( $coef );
            }
            else{
                $delta_molec{$chem} += $coef;
                $molec_count{$comp} += abs( $coef );
            }
        }
    }
    my $delta_H = $delta_formula{H} || 0;
    delete $delta_formula{H} if exists $delta_formula{H};
    my $is_balanced = $delta_H == $delta_charge; # first test if balance can be restored with protons
    if( $is_balanced ){
        foreach( keys %delta_formula ){ # are the the formula really balanced?
            $is_balanced = 0 if $delta_formula{$_} != 0;
        }
    }
    if( $is_balanced ){
        foreach( keys %delta_molec ){ # are the symbol really balanced?
            $is_balanced = 0 if  $delta_molec{$_} != 0;
        }
    }
    if( $is_balanced and $delta_H != 0 ){ # now fix the equation, the challenge is to chose the right compartment
        my @buf = sort keys %molec_count; # first sort alphabetically
        my @comp = sort { $molec_count{$b} <=> $molec_count{$a} } @buf; # then attempt to select the compartment with most molecules
        $equation{$comp[0]}{MNXM1} -= $delta_H; # chose $comp[0] to add/remove proton
    }
    return( $is_balanced, _canonicalize( \%equation ));
}

sub _canonicalize{ # FIXME: do not place proton at first place
    my( $eq ) = @_;
    my %L = ();
    my %R = ();
    foreach my $comp ( keys %$eq ){
        foreach( keys %{$eq->{$comp}} ){
            my $spec = $_ . '@' . $comp;
            if( $eq->{$comp}{$_} < 0 ){
                $L{$spec} -= $eq->{$comp}{$_};
            }
            elsif( $eq->{$comp}{$_} > 0 ){
                $R{$spec} += $eq->{$comp}{$_};
            }
            # else{} # ignore coef with a value of zero
        }
    }
    my @term_left  = ();
    my @term_right = ();
    foreach( sort keys %L ){ # this is a critical sort
        push @term_left, $L{$_} . ' ' . $_
    }
    foreach( sort keys %R ){  # this is a critical sort
        push @term_right, $R{$_} . ' ' . $_
    }
    my $left_str  = join( ' + ', @term_left );
    my $right_str = join( ' + ', @term_right );
    my $L = ' ' . $left_str  . ' ';
    my $R = ' ' . $right_str . ' ';
    if( $R =~ /\@BOUNDARY / or $L !~ / BIOMASS\@/ && $R =~ / BIOMASS\@/ ){
        return( "$left_str = $right_str", 1 );
    }
    elsif( $L =~ /\@BOUNDARY / or $L =~ / BIOMASS\@/ && $R !~ / BIOMASS\@/ ){
        return( "$right_str = $left_str", -1 );
    }
    elsif( join( ' ', sort keys %L ) lt join( ' ', sort keys %R )){
        return ( "$left_str = $right_str", 1 );
    }
    else{ # $left_str gt $right_str
        return ( "$right_str = $left_str", -1 );
    }
}

# ------------------------------------------------ #
# A helper sub
# ------------------------------------------------ #

sub _parse_equation{
    my( $left_str, $arrow_str, $right_str ) = $_[0] =~ /(.*) (<\-+|<=+>|\-+>|<\?>|=+) (.*)/;
    my %equation = ();
    foreach my $token ( split / \+ /, $left_str ){
        my( $coef, $chem, $comp ) = $token =~ /(\S+) (\S+)\@(\S+)/;
        $equation{$comp}{$chem} -= $coef;
    }
    foreach my $token ( split / \+ /, $right_str ){
        my( $coef, $chem, $comp ) = $token =~ /(\S+) (\S+)\@(\S+)/;
        $equation{$comp}{$chem} += $coef;
    }
    return \%equation;
}

# ------------------------------------------------ #
# Exhaustively enumerate all possible compartment
# mappings (grows as factorial!)
# ------------------------------------------------ #

sub _build_comp_map{
    my( $comp_mnx, $comp_new ) = @_;
    my @map = __build_comp_map( $comp_mnx, $comp_new );
    my %map = ();
    foreach my $i ( 1 .. @map ){
        my %tmp = @{$map[$i-1]}; # funky typecasting no?
        $map{"Map_$i"} = \%tmp;
    }
    return \%map;
}
sub __build_comp_map{
    my( $comp_mnx, $comp_new ) = @_;
    my @map = ();
    my @comp_mnx = @{$comp_mnx};
    my $mnx = shift @comp_mnx;
    foreach my $new ( sort keys %$comp_new ){
        my %new = %$comp_new;
        delete $new{$new};
        if( @comp_mnx ){
            foreach( __build_comp_map( \@comp_mnx, \%new )){ # recursion
                push @map, [ @$_, $mnx, $new ];
            }
        }
        else{
            push @map, [ $mnx, $new ];
        }
    }
    return @map;
}

# ------------------------------------------------ #
# This routine attempt to produce integer and
# canonical stoichiometric coefficients. However it
# may also gives up to avoid modifying the reaction
# too much.
# ------------------------------------------------ #

sub _normalize_coef{
    my $self = shift;
    my $eq   = shift;
    my %coef = ();
    my $non_integer = 0;
    my $min_float   = 1;
    my $is_empty    = 1;
    foreach my $comp ( keys %$eq ){
        foreach my $chem ( keys %{$eq->{$comp}} ){
            return $eq if $chem eq 'BIOMASS';
            if( $eq->{$comp}{$chem} == 0 ){
                next;
            }
            else{
                $is_empty = 0;
            }
            my $coef = abs $eq->{$comp}{$chem};
            $coef{$coef} = 1;
            if( $coef !~ /^\-?\d+$/ ){
                $non_integer = 1;
                $min_float = $coef if $min_float > $coef;
            }
        }
    }
    return $eq if $is_empty;
    if( $non_integer ){
        if( $min_float >= 1e-3 ){ # attempt to catch 0.5, 1/8 and % otherwise give up
            my @coef = ();
            push @coef, 1000 * $_ foreach keys %coef;
            my $gcd = _find_gcd( @coef );
            if( $gcd == 1000 ){
                # $self->{tb}->warn( "Small float (gcd failed): " . Dumper $eq );
                return $eq;
            }
            else{ # fix
                my $eq_new = $eq;
                my $ok = 1;
                foreach my $comp ( keys %$eq ){
                    foreach my $chem ( keys %{$eq->{$comp}} ){
                        $eq_new->{$comp}{$chem} = 1000 * $eq->{$comp}{$chem} / $gcd;
                        $ok = 0 unless $eq_new->{$comp}{$chem} =~ /^\-?\d+$/;
                    }
                }
                if( $ok ) {
                    # $self->{tb}->warn( "Small float (gcd=$gcd): " . Dumper $eq_new );
                    return $eq_new;
                }
                else{
                    # $self->{tb}->warn( "Small float (rounding failed): " . Dumper $eq );
                    return $eq;
                }
            }
        }
        else{ # $min_float < 1e-3
            # $self->{tb}->warn( "Very small float: " . Dumper $eq );
            return $eq;
        }
    }
    else{ # only integer coefficient
        my $gcd = _find_gcd( keys %coef );
        return $eq if $gcd == 1;
        my $eq_new = $eq;
        foreach my $comp ( keys %$eq ){
            foreach my $chem ( keys %{$eq->{$comp}} ){
                $eq_new->{$comp}{$chem} = $eq->{$comp}{$chem} / $gcd;
            }
        }
        # $self->{tb}->warn( "Fix integer (gcd=$gcd): " . Dumper $eq_new );
        return $eq_new;
    }
    $self->{tb}->unreachable();
}
sub _gcd{
    my( $a, $b ) = @_;
    return $b if $a == 0;
    return _gcd( $b % $a, $a );
}
sub _find_gcd{ # Inspired by https://www.geeksforgeeks.org/gcd-two-array-numbers/
    my @int = @_;
    my $res = $_[0];
    for( 1 .. @_ - 1 ){
        $res = _gcd( $_[$_], $res );
    }
    return $res;
}

sub chem_exists{
    my( $self, $chem_id ) = @_;
    return exists $self->{chem_prop}{$chem_id};
}
sub comp_exists{
    my( $self, $comp_id ) = @_;
    return exists $self->{comp_prop}{$comp_id};
}
sub reac_exists{
    my( $self, $reac_id ) = @_;
    return exists $self->{reac_prop}{$reac_id};
}
sub get_reac_info{  # return the same list of arguments as MetNet->get_reac_info()
    my( $self, $reac_id ) = @_;
    # xref will be sorted such as the first item is the reference one
    my $reference = $self->{reac_prop}{$reac_id}{reference};
    my @xref = ();
    foreach( keys %{$self->get_reac_xref( $reac_id )}){
        next if $_ eq $reference;
        push @xref, $_;
    }
    return(
        $self->{reac_prop}{$reac_id}{classifs},
        '', # $self->{reac_prop}{$reac_id}{pathway} does not exist actually,
        join( ';', $reference, sort @xref),
    );
}
sub get_chem_info{ # return the same list of five arguments as MetNet->get_chem_info()
    my( $self, $chem_id ) = @_;
    # xref will be sorted such as the first item is the reference one
    my $reference = $self->{chem_prop}{$chem_id}{reference};
    my @xref = ();
    foreach( keys %{$self->get_chem_xref( $chem_id )}){
        next if $_ eq $reference;
        push @xref, $_;
    }
    if( exists $self->{chem_prop}{$chem_id}{InChIKey} ){
        my $inchikey = $self->{chem_prop}{$chem_id}{InChIKey};
        $inchikey =~ s/^InChIKey=//;
        push @xref, 'inchikey:' . $inchikey if $inchikey;
    }
    return (
        $self->{chem_prop}{$chem_id}{name},
        $self->{chem_prop}{$chem_id}{formula},
        $self->{chem_prop}{$chem_id}{mass},
        $self->{chem_prop}{$chem_id}{charge},
        join( ';', $reference, sort @xref),
    );
}
sub get_comp_info{ # same method interface as MetNet->get_comp_info()
    my( $self, $comp_id ) = @_;
    # FIXME: implement xref as above
    return (
        $self->{comp_prop}{$comp_id}{name},
        $self->{comp_prop}{$comp_id}{reference},
    );
}
sub _search_reac_isom{
    my( $self, $equation ) = @_;
    if( my( $child_list ) = $equation =~ / (\w+\|[\w\|]+)\@/ ){
        my %mnxr = ();
        foreach my $item ( split /\|/, $child_list ){
            my $str = $equation;
            $str =~ s/ \Q$child_list\E\@/ $item\@/;
            my @mnxr = $self->_search_reac_isom( $str );
            foreach( @mnxr ){
                $mnxr{$_} = 1 if $_;
            }
        }
        return keys %mnxr;
    }
    else{
        my( $new_id, $eq_new, $sign, $is_balanced, $msg, $mnxr_id, $str_gen, $sign_gen ) = $self->map_equation( $equation );
        return '' if ! $mnxr_id;
        return $mnxr_id;
    }
}
sub search_reac_isom{
    my( $self, $eq_str ) = @_;
    my( $new_id, $eq_new, $sign, $is_balanced, $msg, $mnxr_id, $str_gen, $sign_gen ) = $self->map_equation( $eq_str );
    my $num = $eq_str =~ tr /@/@/;
    return $mnxr_id if $num > 4; # FIXME: this is a trick to avoid combinatorial explosion
    my $ok = 1;
    while( $eq_new =~ / (\w+)\@/g ){ # skip UNK:
        my $chem_id = $1;
        my @child = $self->search_chem_isom( $chem_id );
        if( @child > 1 ){
            my $child_list = join( '|', @child );
            $eq_new =~ s/ $chem_id\@/ $child_list\@/;
            $ok = 0;
        }
    }
    if( $ok ){ # nothing to search for
        return $mnxr_id;
    }
    else{
        return sort grep /\S/, $self->_search_reac_isom( $eq_new );
    }
}

1;


