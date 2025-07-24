package MetNet;

use strict;
no warnings;

use POSIX qw( strftime );
use Digest::CRC;
use Carp;
#use Data::Dumper;

sub compute_reac_id{ # A class method
    my( $equation, $prefix) = @_;
    $prefix = 'reac' unless $prefix;
    my $str = $equation;
    $str =~ s/(<\-+|<=+>|\-+>|\=+)/<?>/;
    my $ctx = Digest::CRC->new( type => 'crc32', poly => 0x2D );
    $ctx->add( $str );
    my $reac_id = $prefix . uc $ctx->hexdigest();
    # $reac_id =~ s/^R0+/R/; # FIXME this is to emulate the behaviour of Java (should be removed at some point)
    return $reac_id;
}

sub new{
    my( $package, %arg ) = @_;
    my $self = {
        comp      => {},   # Compartments  Cytoplasm
        spec      => {},   # Species       Glucos@Cytoplasm.
        chem      => {},   # Chemical compounds     Glucose
        reac      => {},   # Biochemical reactions (with species and direction)
        equa      => {},   # Equation
        enzy      => {},
        pept      => {},
        mnet      => {},   # Models        ..
        look      => {},   # For quick lookup
        counter   => 0,    # Internal counter to give a unique accession number to every object
        number    => 0,    # A mnet counter used by the GUI
    };
    bless $self, $package;
}

sub load{ # load a single mnet from disk
    my( $self, $path, $new_mnet_id ) = @_;

    my $filename = $path.'/model.tsv';
    my $entry = '';
    open(FILE, $filename)  or confess "Cannot read from file: $filename\n";
    while( <FILE> ){
        next  if /^\#/;
        $entry .= $_;
    }
    close FILE;
    my( $old_mnet_id, $desc, $LB, $UB, @info ) = $self->_parse_model_entry( $entry );
    my $mnet_id = $new_mnet_id ? $new_mnet_id : $old_mnet_id;

    $filename = $path.'/enzymes.tsv';
    open(FILE, $filename)  or confess "Cannot read from file: $filename\n";
    my %enzy = ();
    while( <FILE> ){
        next  if /^\#/;
        chop;
        my @token =  split /\t/;
        my $reac_id = shift @token;
        my $dir = undef;
        if( @token == 4 ){
            $dir = $token[3];
            $reac_id =~ s/_(LR|RL|B)$//;
        }
        elsif( $reac_id =~ s/_(LR|RL|B)$// ){
            $dir = $1;
        }
        else{
            confess "Cannot obtain reaction direction: $reac_id\t" . join( "\t", @token ) . "\n";
        }
        push @{$enzy{$reac_id}}, [ splice( @token, 0, 3 ), $dir ];
    }
    close FILE;
    unless( $LB and $UB ){
        foreach my $reac_id ( keys %enzy ){
            foreach( @{$enzy{$reac_id}} ){
                $LB = $_->[1]  if $LB > $_->[1];
                $UB = $_->[2]  if $UB < $_->[2];
            }
        }
        confess 'Cannot guess LB and UB from enzy!'  unless $LB and $UB; # i.e. non-zero
    }
    unless( -$LB == $UB ){
        warn "Changing bounds: $LB and $UB\n";
        if( -$LB < $UB ){
            $LB = -$UB;
        }
        else{
            $UB = -$LB;
        }
    }
    $self->add_mnet( $mnet_id, $LB, $UB );
    $self->set_mnet_desc( $mnet_id, $desc );
    $self->push_mnet_info( $mnet_id, @info );
    my $mnet_num = $self->_get_mnet_num( $mnet_id );

    $filename = $path . '/reactions.tsv';
    open(FILE, $filename)  or confess "Cannot read from file: $filename\n";
    while(<FILE>){
        next  if /^\#/;
        chop;
        my @token =  split /\t/;
        my $reac_id = shift @token;
        $reac_id =~ s/_(LR|RL|B)$//; # FIXME: remove this fossil statement
        my $source = splice @token, 1, 1;
        $self->_store_reac( $mnet_num, $reac_id, splice @token, 0, 2 );
        $self->set_reac_source( $mnet_id, $reac_id, $source );
        $self->set_reac_info( $reac_id, @token );
        foreach( @{$enzy{$reac_id}} ){
            if( $_->[1] == $LB ){
                $_->[1] = 'NA';
            }
            else{
                $_->[1] += 0; # enforce numeric context => reformat
            }
            if( $_->[2] == $UB ){
                $_->[2] = 'NA';
            }
            else{
                $_->[2] += 0; # enforce numeric context => reformat
            }
            $self->_store_enzy($mnet_num, $reac_id, @$_ );
        }
    }
    close FILE;

    $filename = $path.'/chemicals.tsv';
    open(FILE, $filename)  or confess "Cannot read from file: $filename\n";
    while( <FILE> ){
        next  if /^\#/;
        chop;
        my @token  = split /\t/;
        my $source = splice @token, 2, 1;
        $self->set_chem_source( $mnet_id, $token[0], $source )  if $source; # FIXME: this was a quick patch!!!
        $self->set_chem_info( @token );
    }
    close FILE;

    $filename = $path . '/compartments.tsv';
    open(FILE, $filename )  or confess "Cannot read from file: $filename\n";
    while( <FILE> ){
        next  if /^\#/;
        chop;
        my @token = split /\t/;
        my $source = splice @token, 2, 1;
        $self->set_comp_source( $mnet_id, $token[0], $source );
        $self->set_comp_info( @token );
    }
    close FILE;

    $filename = $path . '/peptides.tsv';
    if( -s $filename ){
        open(FILE, $filename )  or confess "Cannot read from file: $filename\n";
        while( <FILE> ){
            next  if /^\#/;
            chop;
            $self->set_pept_info( split /\t/ );
        }
        close FILE;
    }

    return $old_mnet_id; # even if the name has been changed => this permit to save the original name
}
sub mnet_store{ # Direct call is DEPRECATED and this method must become private
    my( $self, $mnet_id, $desc ) = @_;
    my $mnet_num;
    if(exists $self->{look}{mnet}{$mnet_id}){
        $mnet_num = $self->{look}{mnet}{$mnet_id};
    }
    else{
        $mnet_num = ++$self->{counter};
        $self->{look}{mnet}{$mnet_id} = $mnet_num;
    }
    $self->{mnet}{$mnet_num} = {
        id          => $mnet_id,
        desc        => $desc   || '',
        info        => [],
        number      => ++$self->{number},
        reac        => {},
        equa        => {},
        spec        => {},
        chem        => {},
        comp        => {},
        enzy        => {},
        pept        => {},
        chem_source => {},
        comp_source => {},
        reac_source => {},
        LB          => 0,
        UB          => 0,
    };
    $self->{mnet}{$mnet_num}{stamp} = time();
    return $mnet_num;
}

sub _check_LU_bounds{
    my( $self, $LB, $UB ) = @_;
    confess "LB is not negative: $LB"  unless $LB < 0;
    confess "UB is not positive: $UB"  unless $UB > 0;
    confess "LB is not equal to -UB: $LB vs $UB"  unless $LB == -$UB;
}
sub add_mnet{
    my( $self, $mnet_id, $LB, $UB ) = @_;
    confess "Mnet already exists: $mnet_id\n"  if exists $self->{look}{mnet}{$mnet_id};
    $self->_check_LU_bounds( $LB, $UB );
    my $mnet_num = $self->mnet_store( $mnet_id, 'No description' );
    $self->{mnet}{$mnet_num}{LB} = $LB;
    $self->{mnet}{$mnet_num}{UB} = $UB;
}
sub set_mnet_desc{
    my( $self, $mnet_id, $desc ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    $self->{mnet}{$mnet_num}{desc} = $desc;
}
sub get_mnet_desc{
    my( $self, $mnet_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    return $self->{mnet}{$mnet_num}{desc};
}
sub push_mnet_info{
    my $self = shift;
    my $mnet_num = $self->_get_mnet_num( shift );
    push @{$self->{mnet}{$mnet_num}{info}}, @_;
}
sub set_mnet_info{
    my( $self, $mnet_id, @info ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    $self->{mnet}{$mnet_num}{info} = \@info;
}
sub get_mnet_info{
    my $self = shift;
    my $mnet_num = $self->_get_mnet_num( shift );
    return @{$self->{mnet}{$mnet_num}{info}};
}
sub reac_store{
    my $self = shift;
    warn "Method reac_store() is deprecated, use _store_reac() instead!";
    $self->_store_reac( @_ );
}
sub _store_reac{ # if $reac_id is already stored then the following fields are ignored !!!
    my( $self, $mnet_num, $reac_id, $equation, $equa_id ) = @_;
    my $reac_num = 0;
    if( exists $self->{look}{reac}{$reac_id} ){
        $reac_num = $self->{look}{reac}{$reac_id};
        foreach( keys( %{$self->{reac}{$reac_num}{left}} ), keys( %{$self->{reac}{$reac_num}{right}} )){
            my $spec_id = $self->{spec}{$_}{id};
            $self->spec_store( $mnet_num, $reac_num, $spec_id ); # keep the internal data structure up to date
        }
        my $equa_num = $self->{reac}{$reac_num}{equa};
        $equa_id = $self->{equa}{$equa_num}{id};
        $self->equa_store( $mnet_num, $reac_num, $equa_id ); # keep the internal data structure up to date
    }
    else{
        #FIXME check here the values of $reac_id, $equation, $equa_id
        $reac_num = ++$self->{counter};
        $self->{look}{reac}{$reac_id} = $reac_num;
        my( $left, $arrow, $right ) = $equation =~ /(.*) (<\-+|<=+>|\-+>|<\?>|=+) (.*)/;
        confess "Cannot parse equation $reac_id: $equation\n"  unless $left and $arrow and $right;
        $equation =~ s/ (<\-+|<=+>|\-+>|=+) / <?> /;
        my( %left, %right, %comp );
        foreach(split / \+ /, $left){
            my ( $coef, $spec_id ) = /(\S+) (\S+)/;
            my $spec_num = $self->spec_store( $mnet_num, $reac_num, $spec_id );
            $left{$spec_num} = $coef;
            $comp{$self->{spec}{$spec_num}{comp}} = 1;
        }
        foreach(split / \+ /, $right){
            my ( $coef, $spec_id ) = /(\S+) (\S+)/;
            my $spec_num = $self->spec_store( $mnet_num, $reac_num, $spec_id );
            $right{$spec_num} = $coef;
            $comp{$self->{spec}{$spec_num}{comp}} = 1;
        }
        my $equa_num = $self->equa_store( $mnet_num, $reac_num, $equa_id );
        $self->{reac}{$reac_num} = {
            id       => $reac_id,
            equation => $equation,
            equa     => $equa_num,
            left     => \%left,
            right    => \%right,
            EC       => '',
            pathway  => '',
            xref     => '',
            comp     => \%comp,
            mnet     => {}
        };
    }
    $self->{reac}{$reac_num}{mnet}{$mnet_num} = 1;
    $self->{mnet}{$mnet_num}{reac}{$reac_num} = 1;
    return $reac_num;
}
sub set_reac_info{
    my( $self, $reac_id, $EC, $pathway, $xref ) = @_;
    my $reac_num = $self->_get_reac_num( $reac_id );
    $self->{reac}{$reac_num}{EC}      = $EC      || '',
    $self->{reac}{$reac_num}{pathway} = $pathway || '',
    $self->{reac}{$reac_num}{xref}    = $xref    || '',
}
sub get_reac_info{
    my( $self, $reac_id ) = @_;
    my $reac_num = $self->_get_reac_num( $reac_id );
    return (
        $self->{reac}{$reac_num}{EC},
        $self->{reac}{$reac_num}{pathway},
        $self->{reac}{$reac_num}{xref}
    );
}
sub equa_store{
    my($self, $mnet_num, $reac_num, $equa_id) = @_;
    $equa_id = 'NA' unless $equa_id;
    my $equa_num = 0;
    if(exists $self->{look}{equa}{$equa_id}){
        $equa_num = $self->{look}{equa}{$equa_id};
    }
    else{
        $equa_num = ++$self->{counter};
        $self->{look}{equa}{$equa_id} = $equa_num;
        $self->{equa}{$equa_num}{id}  = $equa_id;
    }
    $self->{equa}{$equa_num}{reac}{$reac_num} = 1;
    $self->{equa}{$equa_num}{mnet}{$mnet_num} = 1;
    $self->{mnet}{$mnet_num}{equa}{$equa_num} = 1;
    return $equa_num;
}
sub spec_store{
    my($self,$mnet_num,$reac_num,$spec_id) = @_;
    my $spec_num = 0;
    if( exists $self->{look}{spec}{$spec_id} ){
        $spec_num = $self->{look}{spec}{$spec_id};
        my( $chem_id, $comp_id ) = split '@',$spec_id;
        my $chem_num = $self->{look}{chem}{$chem_id};
        $self->{mnet}{$mnet_num}{chem}{$chem_num} = 1;
        my $comp_num = $self->{look}{comp}{$comp_id};
        $self->{mnet}{$mnet_num}{comp}{$comp_num} = 1;
    }
    else{
        $spec_num = ++$self->{counter};
        $self->{look}{spec}{$spec_id} = $spec_num;
        my( $chem_id, $comp_id ) = split '@', $spec_id;
        $self->{spec}{$spec_num} = {
            id   => $spec_id,
            chem => $self->_store_chem( $mnet_num, $chem_id ),
            comp => $self->_store_comp( $mnet_num, $comp_id )
        };
    }
    $self->{mnet}{$mnet_num}{spec}{$spec_num} = 1;
    $self->{spec}{$spec_num}{mnet}{$mnet_num} = 1;
    $self->{spec}{$spec_num}{reac}{$reac_num} = 1;
    return $spec_num;
}
sub _store_chem{
    my( $self, $mnet_num, $chem_id ) = @_;
    my $chem_num = 0;
    if( exists $self->{look}{chem}{$chem_id} ){
        $chem_num = $self->{look}{chem}{$chem_id};
    }
    else{
        $chem_num = ++$self->{counter};
        $self->{look}{chem}{$chem_id} = $chem_num;
        $self->{chem}{$chem_num} = {
            id      => $chem_id,
            desc    => '',
            source  => '',
            formula => '',
            mass    => '',
            charge  => '',
            xref    => '',
        };
    }
    $self->{mnet}{$mnet_num}{chem}{$chem_num} = 1;
    return $chem_num;
}
sub set_chem_info{
    my( $self, $chem_id, $desc, $formula, $mass, $charge, $xref  ) = @_;
    if( exists $self->{look}{chem}{$chem_id} ){
        my $chem_num = $self->{look}{chem}{$chem_id};
        $self->{chem}{$chem_num}{desc}    = $desc    || '';
        $self->{chem}{$chem_num}{formula} = $formula || '';
        $self->{chem}{$chem_num}{mass}    = defined $mass   ? $mass   : '';  # support for a charge of 0
        $self->{chem}{$chem_num}{charge}  = defined $charge ? $charge : '';  # support for a charge of 0
        $self->{chem}{$chem_num}{xref}    = $xref    || '';
    }
}
sub get_chem_info{
    my( $self, $chem_id ) = @_;
    my $chem_num = $self->_get_chem_num( $chem_id );
    return (
        $self->{chem}{$chem_num}{desc},
        $self->{chem}{$chem_num}{formula},
        $self->{chem}{$chem_num}{mass},
        $self->{chem}{$chem_num}{charge},
        $self->{chem}{$chem_num}{xref}
    );
}
sub get_chem_source{
    my( $self, $mnet_id, $chem_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my $chem_num = $self->_get_chem_num( $chem_id );
    return $self->{mnet}{$mnet_num}{chem_source}{$chem_num};
}
sub set_chem_source{
    my( $self, $mnet_id, $chem_id, $source) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my $chem_num = $self->_get_chem_num( $chem_id );
    $self->{mnet}{$mnet_num}{chem_source}{$chem_num} = $source;
}
sub _store_comp{
    my( $self, $mnet_num, $comp_id ) = @_;
    my $comp_num = 0;
    if( exists $self->{look}{comp}{$comp_id} ){
        $comp_num = $self->{look}{comp}{$comp_id};
    }
    else{
        $comp_num = ++$self->{counter};
        $self->{look}{comp}{$comp_id} = $comp_num;
        $self->{comp}{$comp_num} = {
            id     => $comp_id,
            desc   => '',
            source => ''
        }
    }
    $self->{mnet}{$mnet_num}{comp}{$comp_num} = 1;
    return $comp_num;
}
sub set_comp_info{
    my( $self, $comp_id, $desc, $xref ) = @_;
    if( exists $self->{look}{comp}{$comp_id} ){
        my $comp_num = $self->{look}{comp}{$comp_id};
        $self->{comp}{$comp_num}{desc} = $desc || '';
        $self->{comp}{$comp_num}{xref} = $xref || '';
    }
}
sub get_comp_info{
    my( $self, $comp_id ) = @_;
    my $comp_num = $self->_get_comp_num( $comp_id );
    return (
        $self->{comp}{$comp_num}{desc},
        $self->{comp}{$comp_num}{xref},
    );
}
sub get_comp_source{
    my( $self, $mnet_id, $comp_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my $comp_num = $self->_get_comp_num( $comp_id );
    return $self->{mnet}{$mnet_num}{comp_source}{$comp_num};
}
sub set_comp_source{
    my( $self, $mnet_id, $comp_id, $source) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my $comp_num = $self->_get_comp_num( $comp_id );
    $self->{mnet}{$mnet_num}{comp_source}{$comp_num} = $source;
}
sub _store_enzy{
    my( $self, $mnet_num, $reac_id, $complex, $lower_bound, $upper_bound, $dir) = @_;
    my $enzy_id  = join '|', $reac_id, $complex, $lower_bound, $upper_bound, $dir;
    my $reac_num = $self->{look}{reac}{$reac_id};
    my $enzy_num;
    if(exists $self->{look}{enzy}{$enzy_id}){
        $enzy_num = $self->{look}{enzy}{$enzy_id};
    }
    else{
        $enzy_num = ++$self->{counter};
        $self->{look}{enzy}{$enzy_id}             = $enzy_num;
        $self->{enzy}{$enzy_num}{id}              = $enzy_id; # required for garbage collection ?
        $self->{reac}{$reac_num}{enzy}{$enzy_num} = 1;
        $self->{enzy}{$enzy_num}{reac}            = $reac_num;
        $self->{enzy}{$enzy_num}{complex}         = $complex || '';
        $self->{enzy}{$enzy_num}{LB}              = $lower_bound;
        $self->{enzy}{$enzy_num}{UB}              = $upper_bound;
        $self->{enzy}{$enzy_num}{dir}             = $dir;
    }
    $self->{enzy}{$enzy_num}{mnet}{$mnet_num} = 1;
    $self->{mnet}{$mnet_num}{enzy}{$enzy_num} = 1;
    $complex =~s /\d+\*/ /g; # remove stoichiometric ceofficients FIXME: is this deprecated?
    my( @pept_id ) = $complex =~ /([^;,\+]+)/g; # hence valid separators are ';', '+' and ','
    $self->pept_store( $mnet_num, $reac_num, $_ ) foreach @pept_id;
    return $enzy_num;
}
sub pept_store{
    my( $self, $mnet_num, $reac_num, $pept_id ) = @_;
    my $pept_num;
    if( exists $self->{look}{pept}{$pept_id} ){
        $pept_num = $self->{look}{pept}{$pept_id};
    }
    else{
        $pept_num = ++$self->{counter};
        $self->{look}{pept}{$pept_id} = $pept_num;
        $self->{pept}{$pept_num} = {
            id    => $pept_id,
            desc  => '',
            xrefs => '',
            gene  => '',
        };
    }
    $self->{mnet}{$mnet_num}{pept}{$pept_num} = 1;
    $self->{pept}{$pept_num}{reac}{$reac_num} = 1;
    $self->{reac}{$reac_num}{pept}{$pept_num} = 1;
    $self->{pept}{$pept_num}{mnet}{$mnet_num} = 1;

    return $pept_num;
}
sub set_pept_info{
    my( $self, $pept_id, $desc, $xrefs, $gene ) = @_;
    if( exists $self->{look}{pept}{$pept_id} # store only if gene exists
        and $desc || $xrefs || $gene ){      # and one of the field is non-empty
        my $pept_num = $self->_get_pept_num( $pept_id );
        $self->{pept}{$pept_num}{desc}  = $desc  || '';
        $self->{pept}{$pept_num}{xrefs} = $xrefs || '';
        $self->{pept}{$pept_num}{gene}  = $gene  || '';
    }
}
sub get_pept_info{
    my( $self, $pept_id ) = @_;
    my $pept_num = $self->_get_pept_num( $pept_id );
    return (
        $self->{pept}{$pept_num}{desc},
        $self->{pept}{$pept_num}{xrefs},
        $self->{pept}{$pept_num}{gene},
    );
}

# -------------------------------------------------------- #
# A few methods to propagate info from one model to another
# one.
# -------------------------------------------------------- #

sub get_mnet_ids{
    my $self = shift;
    my @mnet_id = ();
    foreach( keys %{$self->{mnet}} ){
        push @mnet_id, $self->{mnet}{$_}{id};
    }
    return @mnet_id;
}
sub get_reac_id{ # DEPRECATED call select_reac_ids( mnet => $mnet_id ) instead
    my( $self, $mnet_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my @reac_id = ();
    foreach my $reac_num ( keys %{$self->{mnet}{$mnet_num}{reac}} ){
        push @reac_id, $self->{reac}{$reac_num}{id};
    }
    return @reac_id;
}
sub get_growth_reac_ids{ # return a list of reac ID(s) for the growth equation
                         # The list is temptatively sorted by reac source
                         # and possibly not blocked
    my( $self, $mnet_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my %buf = (); # to sort multiple growth equation
    foreach my $reac_num ( $self->_get_growth_reac_nums( $mnet_id )){
        my $reac_id = $self->{reac}{$reac_num}{id};
        my @enzy_info = $self->get_enzy_info( $mnet_id, $reac_id );
        my $source = exists $self->{mnet}{$mnet_num}{reac_source}{$reac_num}
                   ? $self->{mnet}{$mnet_num}{reac_source}{$reac_num}
                   :  '~' . $reac_id; # use '~' because it has ASCII rank 127
        $buf{$reac_id} = ( $enzy_info[2] eq 'NA' ? 'A' : 'B' ) . ' ' . $source; # Alive / Blocked then source
    }
    my @reac_id = ();
    if( keys %buf ){
        @reac_id = sort { $buf{$a} cmp $buf{$b} } keys %buf;
    }
    return @reac_id;
}
sub _get_mnet_num{ # written to be fast
    return $_[0]->{look}{mnet}{$_[1]} // confess "Cannot retrieve mnet: $_[1]\n";
}
sub _get_reac_num{ # written to be fast
    return $_[0]->{look}{reac}{$_[1]} // confess "Cannot retrieve reac: $_[1]\n";
}
sub _get_spec_num{ # written to be fast
    return $_[0]->{look}{spec}{$_[1]} // confess "Cannot retrieve spec: $_[1]\n";
}
sub _get_chem_num{ # written to be fast
    return $_[0]->{look}{chem}{$_[1]} // confess "Cannot retrieve chem: $_[1]\n";
}
sub _get_comp_num{ # written to be fast
    return $_[0]->{look}{comp}{$_[1]} // confess "Cannot retrieve comp: $_[1]\n";
}
sub _get_equa_num{ # written to be fast
    return $_[0]->{look}{equa}{$_[1]} // confess "Cannot retrieve equa: $_[1]\n";
}
sub _get_pept_num{ # written to be fast
    return $_[0]->{look}{pept}{$_[1]} // confess "Cannot retrieve pept: $_[1]\n";
}
sub get_reac_source{
    my( $self, $mnet_id, $reac_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my $reac_num = $self->_get_reac_num( $reac_id );
    return $self->{mnet}{$mnet_num}{reac_source}{$reac_num};
}
sub set_reac_source{
    my( $self, $mnet_id, $reac_id, $source ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my $reac_num = $self->_get_reac_num( $reac_id );
    $self->{mnet}{$mnet_num}{reac_source}{$reac_num} = $source;
}
sub get_reac_dir{                         # return one of 'LR', 'RL', 'B' or 'C' (i.e. composite)
    my( $self, $reac_id, $mnet_id ) = @_; # @mnet_id is optional, but none means all
    my $reac_num = $self->_get_reac_num( $reac_id );
    my @mnet_id = $mnet_id || $self->get_mnet_ids();
    my %dir = ();
    foreach( @mnet_id ){
        my $mnet_num = $self->_get_mnet_num( $_ );
        foreach my $enzy_num ( keys %{$self->{reac}{$reac_num}{enzy}} ){
            next unless exists $self->{enzy}{$enzy_num}{mnet}{$mnet_num};
            $dir{$self->{enzy}{$enzy_num}{dir}} = 1;
        }
    }
    my @dir = keys %dir;
    return $dir[0] if @dir == 1;
    return 'C';
}
sub get_reac_equation{
    my( $self, $reac_id ) = @_;
    my $reac_num = $self->_get_reac_num( $reac_id );
    return $self->{reac}{$reac_num}{equation};
}
sub get_reac_human{ # a helper sub
    my( $self, $reac_id ) = @_;
    my $reac_num = $self->_get_reac_num( $reac_id );
    my %left  = ();
    my %right = ();
    foreach( keys %{$self->{reac}{$reac_num}{left}} ){
        my $coef = $self->{reac}{$reac_num}{left}{$_} == 1 ? '' : "$self->{reac}{$reac_num}{left}{$_} ";
        push @{$left{$self->{comp}{$self->{spec}{$_}{comp}}{desc}}}, "$coef$self->{chem}{$self->{spec}{$_}{chem}}{desc}";
    }
    foreach( keys %{$self->{reac}{$reac_num}{right}} ){
        my $coef = $self->{reac}{$reac_num}{right}{$_} == 1 ? '' : "$self->{reac}{$reac_num}{right}{$_} ";
        push @{$right{$self->{comp}{$self->{spec}{$_}{comp}}{desc}}}, "$coef$self->{chem}{$self->{spec}{$_}{chem}}{desc}";
    }
    my @left_comp  = sort keys %left;
    my @right_comp = sort keys %right;
    if( @left_comp == 1 and @right_comp == 1 and $left_comp[0] eq $right_comp[0] ){
        return join ' ',
            '[',
            join( ' + ', sort @{$left{$left_comp[0]}} ),
            '=',
            join( ' + ', sort @{$right{$right_comp[0]}} ),
            ']@' . $left_comp[0];
    }
    my @left_token  = ();
    foreach( sort keys %left ){
        push @left_token, join ' ',
             '[',
             join( ' + ', sort @{$left{$_}} ),
             ']@' . $_;
    }
    my @right_token = ();
    foreach( sort keys %right ){
        push @right_token, join ' ',
             '[',
             join( ' + ', sort @{$right{$_}} ),
             ']@' . $_;
    }
    return join ' ',
        join( ' + ', @left_token ),
        '=',
        join( ' + ', @right_token );
}
sub get_reac_mnxr{
    my( $self, $reac_id ) = @_;
    my $reac_num = $self->_get_reac_num( $reac_id );
    my $equa_num = $self->{reac}{$reac_num}{equa};
    return $self->{equa}{$equa_num}{id} || '';
}
sub get_enzy_info{ # return a flat list which length is a multiple of 4
    my( $self, $mnet_id, $reac_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my $reac_num = $self->_get_reac_num( $reac_id );
    my @info = ();
    foreach my $enzy_num (keys %{$self->{reac}{$reac_num}{enzy}}){
        next unless exists $self->{enzy}{$enzy_num}{mnet}{$mnet_num};
        push @info,
            $self->{enzy}{$enzy_num}{complex},
            $self->{enzy}{$enzy_num}{LB},
            $self->{enzy}{$enzy_num}{UB},
            $self->{enzy}{$enzy_num}{dir};
    }
    return @info if @info;            # FIXME: remove this test
    return ( '', 'NA', 'NA', 'B');    #        and die here
}
sub get_pept_desc{
    my( $self, $pept_id ) = @_;
    my $pept_num = $self->_get_pept_num( $pept_id );
    return $self->{pept}{$pept_num}{desc};
}
sub get_pept_gene{
    my( $self, $pept_id ) = @_;
    my $pept_num = $self->_get_pept_num( $pept_id );
    return $self->{pept}{$pept_num}{gene};
}
sub get_pept_xrefs{
    my( $self, $pept_id ) = @_;
    my $pept_num = $self->_get_pept_num( $pept_id );
    return $self->{pept}{$pept_num}{xrefs};
}

sub merge_enzy_info{
    my( $self, @enzy_info ) = @_;
    return @enzy_info if @enzy_info == 4;
    confess "Incomplete enzy info: " . join( "\t", @enzy_info ) . "\n" if @enzy_info < 4;
    my %complex = ();
    %complex    = ( $enzy_info[0] => 1 ) if $enzy_info[0]; # case not empty (was a nasty bug to find)
    my $LB      = $enzy_info[1];
    my $UB      = $enzy_info[2];
    my %dir     = ( $enzy_info[3] => 1 );
    splice @enzy_info, 0 , 4;
    while( @enzy_info ){
        my( $complex, $lb, $ub, $dir ) = splice @enzy_info, 0 , 4;
        $complex{$complex} = 1 if $complex;
        $LB = $lb eq 'NA' ? 'NA'
            : $LB eq 'NA' ? 'NA'
            : $lb + $LB; # should work in most cases
        $UB = $ub eq 'NA' ? 'NA'
            : $UB eq 'NA' ? 'NA'
            : $ub + $UB; # should work in most cases
        $dir{$dir} = 1;
    }
    my @dir = keys %dir;
    return ( join( ';', sort keys %complex ), $LB, $UB, ( @dir == 1 ? $dir[0] : 'B' ));
}
sub get_chem_desc{
    my( $self, $chem_id ) = @_;
    my $chem_num = $self->_get_chem_num( $chem_id );
    return $self->{chem}{$chem_num}{desc};
}
sub _copy_reac_enzy{
    my( $self, $old_mnet_num, $new_mnet_num, $reac_num ) = @_;
    $self->_store_reac( $new_mnet_num, $self->{reac}{$reac_num}{id} );
    foreach my $enzy_num (keys %{$self->{reac}{$reac_num}{enzy}}){
        next unless exists $self->{enzy}{$enzy_num}{mnet}{$old_mnet_num};
        $self->_store_enzy( $new_mnet_num,
                            $self->{reac}{$reac_num}{id},
                            $self->{enzy}{$enzy_num}{complex},
                            $self->{enzy}{$enzy_num}{LB},
                            $self->{enzy}{$enzy_num}{UB},
                            $self->{enzy}{$enzy_num}{dir}
        );
    }
}

sub copy_reac_enzy{ # DEPRECATED
    shift->copy_reac_copy_enzy( @_ );
}
sub copy_reac_copy_enzy{ # ( $self, $old_mnet_id, $new_mnet_id, @reac_id )
    my $self         = shift;
    my $old_mnet_num = $self->_get_mnet_num( shift );
    my $new_mnet_num = $self->_get_mnet_num( shift );
    foreach( @_ ){
        $self->_copy_reac_enzy( $old_mnet_num, $new_mnet_num, $self->_get_reac_num( $_ ));
    }
}
sub copy_reac_add_enzy{
    my( $self, $new_mnet_id, $reac_id, @enzy_info ) = @_;
    my $new_mnet_num = $self->_get_mnet_num( $new_mnet_id );
    $self->_store_reac( $new_mnet_num, $reac_id );
    while(  @enzy_info ){
        my( $complex, $LB, $UB, $dir ) = splice( @enzy_info, 0, 4 );
        $self->_store_enzy( $new_mnet_num, $reac_id, $complex, $LB, $UB, $dir );
    }
}
sub add_reac_add_enzy{
    my( $self, $new_mnet_id, $reac_id, $equation, $mnxr_id, @enzy_info ) = @_;
    my $new_mnet_num = $self->_get_mnet_num( $new_mnet_id );
    $self->_store_reac( $new_mnet_num, $reac_id, $equation, $mnxr_id );
    while(  @enzy_info ){
        my( $complex, $LB, $UB, $dir ) = splice( @enzy_info, 0, 4 );
        $self->_store_enzy( $new_mnet_num, $reac_id, $complex, $LB, $UB, $dir );
    }
}

# -------------------------------------------------------- #
# How to delete a single model
# -------------------------------------------------------- #

sub remove{
    my( $self, $mnet_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    foreach my $enzy_num ( keys %{$self->{mnet}{$mnet_num}{enzy}} ){
    #     warn "ENZY_NUM: $enzy_num\n";
    #     foreach my $reac_num ( keys %{$self->{enzy}{$enzy_num}{reac}} ){
    #             # next unless exists $self->{reac}{$reac_num}{enzy}{$enzy_num};
    #         delete $self->{reac}{$reac_num}{enzy}{$enzy_num};
    #     }
        delete $self->{enzy}{$enzy_num}{mnet}{$mnet_num};
    }
    foreach my $reac_num (keys %{$self->{mnet}{$mnet_num}{reac}}){
         delete $self->{reac}{$reac_num}{mnet}{$mnet_num};
    }
    foreach my $spec_num (keys %{$self->{mnet}{$mnet_num}{spec}}){
         delete $self->{spec}{$spec_num}{mnet}{$mnet_num};
    }
    foreach my $equa_num (keys %{$self->{mnet}{$mnet_num}{equa}}){
         delete $self->{equa}{$equa_num}{mnet}{$mnet_num};
    }
    foreach my $pept_num (keys %{$self->{mnet}{$mnet_num}{pept}}){
         delete $self->{pept}{$pept_num}{mnet}{$mnet_num};
    }
    delete $self->{look}{mnet}{$mnet_id};
    delete $self->{mnet}{$mnet_num};

    # and now act as a reference-based garbage collector
    foreach my $reac_num ( keys %{$self->{reac}} ){
        if( 0 == keys %{$self->{reac}{$reac_num}{mnet}} ){
            delete $self->{look}{reac}{$self->{reac}{$reac_num}{id}};
            delete $self->{reac}{$reac_num};
        }
    }
    foreach my $enzy_num ( keys %{$self->{enzy}} ){
        if( 0 == keys %{$self->{enzy}{$enzy_num}{mnet} }){
            delete $self->{look}{enzy}{$self->{enzy}{$enzy_num}{id}};
            delete $self->{enzy}{$enzy_num};
        }
    }
    foreach my $equa_num ( keys %{$self->{equa}} ){
        if( 0 == keys %{$self->{equa}{$equa_num}{mnet}} ){
            delete $self->{look}{equa}{$self->{equa}{$equa_num}{id}};
            delete $self->{equa}{$equa_num};
        }
    }
    foreach my $pept_num ( keys %{$self->{pept}} ){
        if( 0 == keys %{$self->{pept}{$pept_num}{mnet}} ){
            delete $self->{look}{pept}{$self->{pept}{$pept_num}{id}};
            delete $self->{pept}{$pept_num};
        }
    }
    my(%chem_ok,%comp_ok);
    foreach my $spec_num (keys %{$self->{spec}}){
        if( 0 == keys %{$self->{spec}{$spec_num}{mnet}} ){
            delete $self->{look}{spec}{$self->{spec}{$spec_num}{id}};
            delete $self->{spec}{$spec_num};
        }
        else{
            $chem_ok{$self->{spec}{$spec_num}{chem}} = 1;
            $comp_ok{$self->{spec}{$spec_num}{comp}} = 1;
        }
    }
    foreach my $chem_num ( keys %{$self->{chem}} ){
        unless($chem_ok{$chem_num}){
            delete $self->{look}{chem}{$self->{chem}{$chem_num}{id}};
            delete $self->{chem}{$chem_num};
        }
    }
    foreach my $comp_num ( keys %{$self->{comp}} ){
        unless($comp_ok{$comp_num}){
            delete $self->{look}{comp}{$self->{comp}{$comp_num}{id}};
            delete $self->{comp}{$comp_num};
        }
    }
}

sub rename{
    my( $self, $mnet_id, $new_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    confess "mnet already exists: $new_id\n" if exists $self->{look}{mnet}{$new_id};
    delete $self->{look}{mnet}{$mnet_id};
    $self->{look}{mnet}{$new_id} = $mnet_num;
    $self->{mnet}{$mnet_num}{id} = $new_id;
}
sub describe{
    my($self,$mnet_id,$text) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    $self->{mnet}{$mnet_num}{desc} = $text;
}

#MetaNetX/TSV/model.tsv 1.1
#property   value
#MetaNetX/TSV/reactions.tsv 1.1
#ID    equation    source  mnxr_ID  classifs   pathways    xrefs
#MetaNetX/TSV/enzymes.tsv 1.1
#reac_ID    complex     LB   UB     direction
#MetaNetX/TSV/chemicals.tsv 1.1
#ID    name    source    formula    mass    charge    xrefs
#MetaNetX/TSV/compartments.tsv 1.1
#ID    name    source  xrefs
#MetaNetX/TSV/peptides.tsv 1.1
#ID    description    xrefs    gene_names

#MetaNetX/MNXref/chem_prop.tsv 1.1
#ID    name    reference  formula  charge  mass InChI  InChIKey SMILES
#MetaNetX/MNXref/chem_xref.tsv 1.1
#xref    chem_ID   description
#MetaNetX/MNXref/comp_prop.tsv 1.1
#ID    name    reference
#MetaNetX/MNXref/comp_xref.tsv 1.1
#xref   comp_ID     description
#MetaNetX/MNXref/reac_prop.tsv 1.1
#ID    equation     reference   classifs   is_balanced  is_transport
#MetaNetX/MNXref/reac_xref.tsv 1.1
#xref    mnxr_ID    description

# UNK:

sub write{
    my( $self, $path, $mnet_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my $filename = "$path/model.tsv";
    open( METNET, "> $filename" )  or confess "Cannot write to file: $filename\n";
    my $info = $self->{mnet}{$mnet_num};
    print METNET join "\n",
        "#MetaNetX/TSV/model.tsv 1.1",
        "#property\tvalue",
        "ID\t$info->{id}",
        "Description\t$info->{desc}",
        "LB\t$info->{LB}",
        "UB\t$info->{UB}" . "\n";
    foreach( @{$info->{info}} ){
        print METNET "$_->[0]\t$_->[1]\n";
    }
    close METNET;
    $filename = "$path/reactions.tsv";
    open( REAC, "> $filename")  or confess "Cannot write to file: $filename\n";
    print REAC "#MetaNetX/TSV/reactions.tsv 1.1\n";
    print REAC "#ID\tequation\tsource\tmnxr_ID\tclassifs\tpathways\txrefs\n";
    close REAC;
    open( REAC, "| sort - >> $filename") or confess "Cannot append to file: $filename\n";
    $filename = "$path/chemicals.tsv";
    open( CHEM, "> $filename")  or confess "Cannot write to file: $filename\n";
    print CHEM "#MetaNetX/TSV/chemicals.tsv 1.1\n";
    print CHEM "#ID    name    source    formula    mass    charge    xrefs\n";
    close CHEM;
    open( CHEM, "| sort - >> $filename") or confess "Cannot append to file: $filename\n";
    $filename = "$path/compartments.tsv";
    open( COMP, "> $filename") or confess "Cannot write to file: $filename\n";
    print COMP "#MetaNetX/TSV/compartments.tsv 1.1\n";
    print COMP "#ID\tname\tsource\txrefs\n";
    close COMP;
    open( COMP, "| sort - >> $filename") or confess "Cannot append to file: $filename\n";
    $filename = "$path/enzymes.tsv";
    open( ENZY, "> $filename")  or confess "Cannot write to file: $filename\n";
    print ENZY "#MetaNetX/TSV/enzymes.tsv 1.1\n";
    print ENZY "#reac_ID\tcomplex\tLB\tUB\tdirection\n";
    close ENZY;
    open( ENZY, "| sort - >> $filename")  or confess "Cannot appendto file: $filename\n";
    $filename = "$path/peptides.tsv";
    open( PEPT, "> $filename")  or confess "Cannot write to file: $filename\n";
    print PEPT "#MetaNetX/TSV/peptides.tsv 1.1\n";
    print PEPT "#ID\tdescription\txrefs\tgene_names\n";
    close PEPT;
    open( PEPT, "| sort - >> $filename")  or confess "Cannot append to file: $filename\n";
    $filename = "$path/stoichiometries.tsv";
    open( STOI, "> $filename")  or confess "Cannot write to file: $filename\n";
    print STOI "#MetaNetX/TSV/stoichiometries.tsv 1.1\n";
    print STOI "#reac_ID\tspecies_ID\tstoichiometry\tchemical_ID\tcompartment_ID\n";
    close STOI;
    open( STOI, "| sort - >> $filename")  or confess "Cannot append to file: $filename\n";
    my %seen = ();
    foreach my $reac_num (keys %{$self->{mnet}{$mnet_num}{reac}}){
        my $reac_info = $self->{reac}{$reac_num};
        my $reac_dir  = $self->get_reac_dir( $reac_info->{id}, $mnet_id );
        my $equation  = $reac_info->{equation};
        if( $reac_dir eq 'LR' ){
            $equation=~ s/ <\?> / --> /;
        }
        elsif( $reac_dir eq 'RL' ){
            $equation=~ s/ <\?> / <-- /;
        }
        else{
            $equation=~ s/ <\?> / <==> /;
        }
        print REAC join "\t",
            $reac_info->{id},
            $equation,
            $self->get_reac_source( $mnet_id, $reac_info->{id} ),
            $self->{equa}{$reac_info->{equa}}{id} ne 'NA' ? $self->{equa}{$reac_info->{equa}}{id} : '',
            $reac_info->{EC},
            $reac_info->{pathway},
            $reac_info->{xref} . "\n";
        foreach my $spec_num ( keys %{$reac_info->{left}},
                               keys %{$reac_info->{right}} ){
            unless( $seen{$spec_num} ){
                unless( $seen{$self->{spec}{$spec_num}{chem}} ){
                    my $chem_num  = $self->{spec}{$spec_num}{chem};
                    my $chem_info = $self->{chem}{$chem_num};
                    next  if ( !$chem_info->{id} );
                    print CHEM join( "\t",
                                     $chem_info->{id},
                                     $chem_info->{desc},
                                     $self->get_chem_source( $mnet_id, $chem_info->{id} ),
                                     $chem_info->{formula},
                                     $chem_info->{mass},
                                     $chem_info->{charge},
                                     $chem_info->{xref} ) ."\n";
                    $seen{$chem_num} = 1;
                }
                unless($seen{$self->{spec}{$spec_num}{comp}}){
                    my $comp_num  = $self->{spec}{$spec_num}{comp};
                    my $comp_info = $self->{comp}{$comp_num};
                    next  if ( !$comp_info->{id} );
                    print COMP join "\t",
                                    $comp_info->{id},
                                    $comp_info->{desc},
                                    $self->get_comp_source( $mnet_id, $comp_info->{id} ),
                                    $comp_info->{xref} . "\n";
                    $seen{$comp_num} = 1;
                }
                $seen{$spec_num} = 1;
            }
        }
    }
    foreach my $enzy_num ( keys %{$self->{mnet}{$mnet_num}{enzy}} ){
        my $reac_num = $self->{enzy}{$enzy_num}{reac};
        print ENZY join( "\t",
                         $self->{reac}{$reac_num}{id},
                         $self->{enzy}{$enzy_num}{complex},
                         $self->{enzy}{$enzy_num}{LB} eq 'NA' ? $self->{mnet}{$mnet_num}{LB} : $self->{enzy}{$enzy_num}{LB},
                         $self->{enzy}{$enzy_num}{UB} eq 'NA' ? $self->{mnet}{$mnet_num}{UB} : $self->{enzy}{$enzy_num}{UB},
                         $self->{enzy}{$enzy_num}{dir} ). "\n";
    }
    foreach my $pept_num ( keys %{$self->{mnet}{$mnet_num}{pept}} ){
        next unless $self->{pept}{$pept_num}{desc} || $self->{pept}{$pept_num}{xrefs} || $self->{pept}{$pept_num}{gene};
        print PEPT join "\t",
            $self->{pept}{$pept_num}{id},
            $self->{pept}{$pept_num}{desc},
            $self->{pept}{$pept_num}{xrefs},
            $self->{pept}{$pept_num}{gene} . "\n";
    }
    foreach my $reac_num (keys %{$self->{mnet}{$mnet_num}{reac}}){
        my $reac_info = $self->{reac}{$reac_num};
        foreach my $spec_num ( sort keys %{$reac_info->{left}} ){
            print STOI join "\t",
                $reac_info->{id},
                $self->{spec}{$spec_num}{id},
                - $reac_info->{left}{$spec_num},
                $self->{chem}{$self->{spec}{$spec_num}{chem}}{id},
                $self->{comp}{$self->{spec}{$spec_num}{comp}}{id} . "\n";
        }
        foreach my $spec_num ( sort keys %{$reac_info->{right}} ){
            print STOI join "\t",
                $reac_info->{id},
                $self->{spec}{$spec_num}{id},
                $reac_info->{right}{$spec_num},
                $self->{chem}{$self->{spec}{$spec_num}{chem}}{id},
                $self->{comp}{$self->{spec}{$spec_num}{comp}}{id} . "\n";
        }
    }
    close REAC;
    close CHEM;
    close COMP;
    close ENZY;
    close PEPT;
    close STOI;
}
sub _parse_model_entry{ # Internal cooking !
    my( $self, $entry ) = @_;
    my( $mnet_id, $desc ) = ( '', '' );
    my( $LB, $UB ) = ( 0, 0 );
    my @info = ();
    foreach( split /\n/, $entry ){
        next  if /^\#/;
        my($key, $value) = split /\t/;
        if( $key eq 'ID' ){
            $mnet_id = $value;
        }
        elsif( $key eq 'Description' ){
            $desc = $value;
        }
        elsif( $key eq 'LB' ){
            $LB = $value;
        }
        elsif( $key eq 'UB' ){
            $UB = $value;
        }
        # elsif( $key =~ /^(MNXref Version|Source (ID|Reference|Link|Date)|Processed Date|Taxid)$/ ){
        else{
            push @info, [ $key, $value ];
        }
    }
    return ( $mnet_id, $desc, $LB, $UB, @info );
}
sub logic{
    my( $self, $mnet_id_1, $mnet_id_2, $reac_op, $prot_op, $new_id ) = @_;
    my $mnet_1 = $self->_get_mnet_num( $mnet_id_1);
    my $mnet_2 = $self->_get_mnet_num( $mnet_id_2 );
    confess "mnet already exists: $new_id\n" if exists $self->{look}{mnet}{$new_id};
    my $mnet_num = $self->mnet_store( $new_id,
                                      "$mnet_id_1 $reac_op $mnet_id_2" );
                                        # strftime('%d %b %Y',localtime()) );
    my @reac_num;
    if($reac_op eq 'OR'){
        my %buf;
        $buf{$_} = 1 foreach keys %{$self->{mnet}{$mnet_1}{reac}};
        $buf{$_} = 1 foreach keys %{$self->{mnet}{$mnet_2}{reac}};
        @reac_num = keys %buf;
    }
    elsif($reac_op eq 'AND'){
        foreach ( keys %{$self->{mnet}{$mnet_1}{reac}} ){
            push @reac_num,$_ if exists $self->{mnet}{$mnet_2}{reac}{$_};
        }
    }
    elsif($reac_op eq 'M1'){
        @reac_num = keys %{$self->{mnet}{$mnet_1}{reac}};
    }
    elsif($reac_op eq 'M2'){
        @reac_num = keys %{$self->{mnet}{$mnet_2}{reac}};
    }
    elsif($reac_op eq 'NOT'){
        foreach ( keys %{$self->{mnet}{$mnet_1}{reac}} ){
            push @reac_num,$_ unless exists $self->{mnet}{$mnet_2}{reac}{$_};
        }
    }
    foreach(@reac_num){
        $self->_store_reac( $mnet_num, $self->{reac}{$_}{id} );
    }
    my( $LB1, $UB1 ) = $self->get_LU_bounds( $mnet_id_1 );
    my( $LB2, $UB2 ) = $self->get_LU_bounds( $mnet_id_2 );
    my $LB = $LB1 < $LB2 ? $LB1 : $LB2;
    my $UB = $UB1 > $UB2 ? $UB1 : $UB2;
    $self->set_LU_bounds( $new_id, $LB , $UB);
    if( $prot_op eq 'M1' ){
        foreach my $reac_num (@reac_num){
            foreach my $enzy_num (keys %{$self->{reac}{$reac_num}{enzy}}){
                next unless exists $self->{enzy}{$enzy_num}{mnet}{$mnet_1};
                $self->_store_enzy( $mnet_num,
                                    $self->{reac}{$reac_num}{id},
                                    $self->{enzy}{$enzy_num}{complex},
                                    $self->{enzy}{$enzy_num}{LB},
                                    $self->{enzy}{$enzy_num}{UB},
                                    $self->{enzy}{$enzy_num}{dir}
                );
            }
        }
    }
    elsif($prot_op eq 'M2' ){
        foreach my $reac_num (@reac_num){
            foreach my $enzy_num (keys %{$self->{reac}{$reac_num}{enzy}}){
                 next unless exists $self->{enzy}{$enzy_num}{mnet}{$mnet_2};
                 $self->_store_enzy( $mnet_num,
                                     $self->{reac}{$reac_num}{id},
                                     $self->{enzy}{$enzy_num}{complex},
                                     $self->{enzy}{$enzy_num}{LB},
                                     $self->{enzy}{$enzy_num}{UB},
                                     $self->{enzy}{$enzy_num}{dir}
                );
            }
        }
    }
    elsif( $prot_op eq 'OR'){
        foreach my $reac_num (@reac_num){
            foreach my $enzy_num (keys %{$self->{reac}{$reac_num}{enzy}}){
                if( exists $self->{mnet}{$mnet_1}{enzy}{$enzy_num}
                    or exists $self->{mnet}{$mnet_2}{enzy}{$enzy_num} ){
                    $self->_store_enzy( $mnet_num,
                                       $self->{reac}{$reac_num}{id},
                                       $self->{enzy}{$enzy_num}{complex},
                                       $self->{enzy}{$enzy_num}{LB},
                                       $self->{enzy}{$enzy_num}{UB},
                                       $self->{enzy}{$enzy_num}{dir}
                    );
                }
            }
        }
    }
    else{# $prot_op eq 'AND'
        foreach my $reac_num (@reac_num){
            foreach my $enzy_num ( keys %{$self->{reac}{$reac_num}{enzy}} ){
                if( exists $self->{mnet}{$mnet_1}{enzy}{$enzy_num}
                    and exists $self->{mnet}{$mnet_2}{enzy}{$enzy_num} ){
                    $self->_store_enzy( $mnet_num,
                                       $self->{reac}{$reac_num}{id},
                                       $self->{enzy}{$enzy_num}{complex},
                                       $self->{enzy}{$enzy_num}{LB},
                                       $self->{enzy}{$enzy_num}{UB},
                                       $self->{enzy}{$enzy_num}{dir}
                    );
                }
            }
        }
    }
}
sub get_mnets{
    my($self,%arg) = @_; # key = chem/comp/spec/reac
    return sort {$a<=>$b} keys %{$self->{mnet}} unless %arg;
    my %num;
    $num{$_} = 1 foreach keys %{$self->{mnet}};
    if(exists $arg{chem}){
        my $chem_num = $self->{look}{chem}{$arg{chem}};
        foreach(keys %num){
            delete $num{$_} unless exists $self->{mnet}{$_}{chem}{$chem_num};
        }
    }
    if(exists $arg{comp}){
        my $comp_num = $self->{look}{chem}{$arg{comp}};
        foreach(keys %num){
            delete $num{$_} unless exists $self->{mnet}{$_}{comp}{$comp_num};
        }
    }
    if(exists $arg{spec}){
        my $spec_num = $self->{look}{spec}{$arg{spec}};
        foreach(keys %num){
            delete $num{$_} unless exists $self->{mnet}{$_}{spec}{$spec_num};
        }
    }
    if(exists $arg{reac}){
        my $reac_num = $self->{look}{reac}{$arg{reac}};
        foreach(keys %num){
            delete $num{$_} unless exists $self->{mnet}{$_}{reac}{$reac_num};
        }
    }
    return sort {$a<=>$b} keys %num;
}
sub get_ready_mnet_nums{
    my $self = shift;
    my @mnet_num = ();
    foreach my $mnet_num ( keys %{$self->{mnet}} ){
        push @mnet_num, $mnet_num if keys %{$self->{mnet}{$mnet_num}{reac}};
    }
    return sort { $a <=> $b } @mnet_num;
}
sub get_mnet_id{
    my($self,$mnet_num) = @_;
    return $self->{mnet}{$mnet_num}{id};
}
sub get_mnet_code{
    my($self, $mnet_num) = @_;
    return '#'.$self->{mnet}{$mnet_num}{number};
}

# -------------------------------------------------------- #
# Internal query methods. All these are returning HASH ref
# -------------------------------------------------------- #

sub _intersect_num{ # a class method
    my( $num_1, $num_2 ) = @_;
    my %num = ();
    foreach( keys %$num_1 ){
        $num{$_} = 1 if exists $num_2->{$_};
    }
    return \%num;
}
sub _select_reac_num{
    my( $self, %arg ) = @_;
    return $self->{reac} unless %arg;
    foreach( keys %arg ){
        confess "Invalid filter key for reac: $_\n" unless /^(mnet|spec|chem|comp|pept|equa)$/; # FIXME: add chem and comp
    }
    my $num = undef; # evaluates to FALSE, while {} evaluates to TRUE
    if( exists $arg{pept} ){
         my $pept_num = $self->_get_pept_num( $arg{pept} );
         $num = $self->{pept}{$pept_num}{reac};
    }
    if( exists $arg{equa} ){
        my $equa_num = $self->_get_equa_num( $arg{equa} );
        my $num2     = $self->{equa}{$equa_num}{reac};
        $num         = $num ? _intersect_num( $num, $num2 ) : $num2;
    }
    if( exists $arg{spec} ){
        my $spec_num = $self->_get_spec_num( $arg{spec} );
        my $num2     = $self->{spec}{$spec_num}{reac};
        $num         = $num ? _intersect_num( $num, $num2 ) : $num2;
    }
    if( exists $arg{mnet} ){
        my $mnet_num = $self->_get_mnet_num( $arg{mnet} );
        my $num2     = $self->{mnet}{$mnet_num}{reac};
        $num         = $num ? _intersect_num( $num, $num2 ) : $num2;
    }
    if( exists $arg{chem} ){
        my %num2 = ();
        foreach my $spec_num ( keys %{$self->_select_spec_num( chem => $arg{chem} )}){
            foreach( keys %{$self->{spec}{$spec_num}{reac}} ){
                $num2{$_} = 1;
            }
        }
        $num = $num ? _intersect_num( $num, \%num2 ) : \%num2;
    }
    if( exists $arg{comp} ){
        my %num2 = ();
        foreach my $spec_num ( keys %{$self->_select_spec_num( comp => $arg{comp} )}){
            foreach( keys %{$self->{spec}{$spec_num}{reac}} ){
                $num2{$_} = 1;
            }
        }
        $num = $num ? _intersect_num( $num, \%num2 ) : \%num2;
    }
    return $num // {};
}
sub _select_spec_num{
    my( $self, %arg ) = @_;
    return $self->{spec} unless %arg;
    foreach( keys %arg ){
        confess "Invalid filter key for spec: $_\n" unless /^(mnet|comp|reac|chem)$/; # FIXME: add chem and comp
    }
    my $num = undef; # evaluates to FALSE, while {} evaluates to TRUE
    if( exists $arg{chem} ){
        my $chem_num = $self->_get_chem_num( $arg{chem} );
        my %num2 = ();
        foreach( keys %{$self->{spec}} ){
            $num2{$_} = 1 if $self->{spec}{$_}{chem} eq $chem_num;
        }
        $num = $num ? _intersect_num( $num, \%num2 ) : \%num2;
    }
    if( exists $arg{comp} ){
        my $comp_num = $self->_get_comp_num( $arg{comp} );
        my %num2 = ();
        foreach( keys %{$self->{spec}} ){
            $num2{$_} = 1 if $self->{spec}{$_}{comp} eq $comp_num;
        }
        $num = $num ? _intersect_num( $num, \%num2 ) : \%num2;
    }
    if( exists $arg{reac} ){
        my $reac_num = $self->_get_reac_num( $arg{reac} );
        my %num2 = ();
        foreach( keys %{$self->{reac}{$reac_num}{left}}, keys %{$self->{reac}{$reac_num}{right}} ){
            $num2{$_} = 1;
        }
        $num = $num ? _intersect_num( $num, \%num2 ) : \%num2;
    }
    if( exists $arg{mnet} ){
        my $mnet_num = $self->_get_mnet_num( $arg{mnet} );
        my $num2     = $self->{mnet}{$mnet_num}{spec};
        $num         = $num ? _intersect_num( $num, $num2 ) : $num2;
    }
    return $num // {};
}
sub _select_chem_num{
    my( $self, %arg ) = @_;
    return $self->{chem} unless %arg;
    foreach( keys %arg ){
        confess "Invalid filter key for _select_chem_num: $_\n" unless /^(mnet|reac|comp)$/; # FIXME: add more
    }
    my $num = undef; # evaluates to FALSE, while {} evaluates to TRUE
    if( exists $arg{reac} ){
        my %num2 = ();
        foreach( keys %{$self->_select_spec_num( reac => $arg{reac} )}){
            $num2{$self->{spec}{$_}{chem}} = 1;
        }
        $num = $num ? _intersect_num( $num, \%num2 ) : \%num2;
    }
    if( exists $arg{comp} ){
        my %num2 = ();
        foreach( keys %{$self->_select_spec_num( comp => $arg{comp} )}){
            $num2{$self->{spec}{$_}{chem}} = 1;
        }
        $num = $num ? _intersect_num( $num, \%num2 ) : \%num2;
    }
    if( exists $arg{mnet} ){
        my $mnet_num = $self->_get_mnet_num( $arg{mnet} );
        my $num2 = $self->{mnet}{$mnet_num}{chem};
        $num = $num ? _intersect_num( $num, $num2 ) : $num2;
    }
    return $num // {};
}
sub _select_comp_num{
    my( $self, %arg ) = @_;
    return $self->{comp} unless %arg;
    foreach( keys %arg ){
        confess "Invalid filter key for _select_comp_num: $_\n" unless /(mnet|reac)$/; # FIXME: add more
    }
    my $num = undef; # evaluates to FALSE, while {} evaluates to TRUE
    if( exists $arg{reac} ){
        my %num2 = ();
        foreach( keys %{$self->_select_spec_num( reac => $arg{reac} )}){
            $num2{$self->{spec}{$_}{comp}} = 1;
        }
        $num = $num ? _intersect_num( $num, \%num2 ) : \%num2;
    }
    if( exists $arg{mnet} ){
        my $mnet_num = $self->_get_mnet_num( $arg{mnet} );
        $num = $num ? _intersect_num( $num, $self->{mnet}{$mnet_num}{comp} ) : $self->{mnet}{$mnet_num}{comp};
    }
    return $num // {};
}
sub _select_pept_num{
    my( $self, %arg ) = @_;
    return $self->{pept} unless %arg;
    foreach( keys %arg ){
        confess "Invalid filter key for _select_pept_num: $_\n" unless /^(mnet|reac)$/; # FIXME: add chem and comp
    }
    my $num = undef; # evaluates to FALSE, while {} evaluates to TRUE
    if( exists $arg{mnet} ){
        my $mnet_num = $self->_get_mnet_num( $arg{mnet} );
        $num = $num ? _intersect_num( $num, $self->{mnet}{$mnet_num}{pept} ) : $self->{mnet}{$mnet_num}{pept};
    }
    if( exists $arg{reac} ){
        my $reac_num = $self->_get_reac_num( $arg{reac} );
        $num = $num ? _intersect_num( $num, $self->{reac}{$reac_num}{pept} ) : $self->{reac}{$reac_num}{pept};
    }
    return $num // {};
}

# -------------------------------------------------------- #
# Public query methods
# -------------------------------------------------------- #

sub _num_to_ids{
    my( $self, $type, $num ) = @_;
    my $buf = $self->{$type};
    my @id = ();
    foreach( keys %$num ){
        push @id, $buf->{$_}{id};
    }
    return @id;
}

sub select_reac_count{
    return scalar keys %{ shift->_select_reac_num( @_ ) };
}
sub select_reac_ids{
    my $self = shift;
    return( $self->_num_to_ids( 'reac', $self->_select_reac_num( @_ )));
}
sub select_spec_count{
    return scalar keys %{ shift->_select_spec_num( @_ ) };
}
sub select_spec_ids{
    my $self = shift;
    return( $self->_num_to_ids( 'spec', $self->_select_spec_num( @_ )));
}
sub select_chem_count{
    return scalar keys %{ shift->_select_chem_num( @_ ) };
}
sub select_chem_ids{
    my $self = shift;
    return( $self->_num_to_ids( 'chem', $self->_select_chem_num( @_ )));
}
sub select_comp_count{
    return scalar keys %{ shift->_select_comp_num( @_ ) };
}
sub select_comp_ids{
    my $self = shift;
    return( $self->_num_to_ids( 'comp', $self->_select_comp_num( @_ )));
}
sub select_pept_count{
    return scalar keys %{ shift->_select_pept_num( @_ ) };
}
sub select_pept_ids{
    my $self = shift;
    return( $self->_num_to_ids( 'pept', $self->_select_pept_num( @_ )));
}

# -------------------------------------------------------- #
# More methods
# -------------------------------------------------------- #

sub get_equas{ # DEPRECATED
    my($self,%arg) = @_; # Supported filter keys: mnet, chem, comp
    return keys %{$self->{equa}} unless %arg;
    my @equa_num;
    if( exists $arg{mnet} ){
        my $mnet_num = $self->_get_mnet_num( $arg{mnet} );
        return keys %{$self->{mnet}{$mnet_num}{equa}};
    }
}
sub get_specs{ # DEPRECATED
    my($self,%arg) = @_; # key = mnet OR reac|spec|chem|comp
    return keys %{ $self->_select_spec_num( %arg )};
}
sub get_chems{ # DEPRECATED
    my($self,%arg) = @_;
    return keys %{ $self->_select_chem_num( %arg )};
}
sub get_comps{
    my($self,%arg) = @_;
    return keys %{ $self->_select_comp_num( %arg )};
}
sub get_enzys{ # DEPRECATED
    my( $self, %arg ) = @_;
    if( exists $arg{mnet} ){
        my $mnet_num = $self->_get_mnet_num( $arg{mnet} );
        return keys %{$self->{mnet}{$mnet_num}{enzy}};
    }
    elsif( exists $arg{equa} ){ # really ???
        my $equa_num = $self->{look}{equa}{$arg{equa}};
        my %enzy_ok;
        foreach my $reac_num (keys %{$self->{equa}{$equa_num}{reac}}){
            foreach (keys %{$self->{reac}{$reac_num}{enzy}}){
                $enzy_ok{$_} = 1;
            }
        }
        return keys %enzy_ok;
    }
    else{
        return keys %{$self->{enzy}};
    }
}
sub get_pepts{ # DEPRECATED
    my($self,%arg) = @_;
    return keys %{ $self->_select_pept_num( %arg )};
}
sub _get_growth_reac_nums{ # return the unsorted list of all growth reactions
    my( $self , $mnet_id ) = @_;
    eval{ $self->_get_chem_num( 'BIOMASS' ) };
    return () if $@; # BIOMASS not found!!!
    my @reac_num = ();
    foreach my $reac_num ( keys %{$self->_select_reac_num( mnet => $mnet_id, chem => 'BIOMASS' )}){
        my $ok = 0;
        foreach my $spec_num ( keys %{$self->{reac}{$reac_num}{left}},
                               keys %{$self->{reac}{$reac_num}{right}} ){
            my $chem_num = $self->{spec}{$spec_num}{chem};
            $ok = 1 if $self->{chem}{$chem_num}{id} ne 'BIOMASS'; # at least another chemical compounds => GROWTH by definition!!!
        }
        push @reac_num, $reac_num if $ok;
    }
    return @reac_num;
}
# sub get_boundary_reacs{ # DEPRECATED
#    my( $self , $mnet_id ) = @_;
#    eval{ $self->_get_comp_num( 'BOUNDARY' ) };
#    return () if $@; # BOUNDARY not found!!!
#    return $self->get_reacs( mnet => $mnet_id, comp => 'BOUNDARY' );
# }
sub get_boundary_reac_ids{ # may return ()
    my( $self , $mnet_id ) = @_;
    eval{ $self->_get_comp_num( 'BOUNDARY' ) };
    return () if $@; # BOUNDARY not found!!!
    return $self->select_reac_ids( mnet => $mnet_id, comp => 'BOUNDARY' );
}
sub get_LU_bounds{
    my( $self, $mnet_id ) = @_;
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    return( $self->{mnet}{$mnet_num}{LB}, $self->{mnet}{$mnet_num}{UB} );
}
sub set_LU_bounds{
    my( $self, $mnet_id, $LB, $UB ) = @_;
    $self->_check_LU_bounds( $LB, $UB );
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    $self->{mnet}{$mnet_num}{LB} = $LB;
    $self->{mnet}{$mnet_num}{UB} = $UB;
}
sub get_reac_bounds{ # ensure that every reaction in the model has bounds
    my( $self, $mnet_id ) = @_;
    my( %LB, %UB );
    my $mnet_num = $self->_get_mnet_num( $mnet_id );
    my( $LB, $UB ) = $self->get_LU_bounds( $mnet_id );
    foreach my $enzy_num ( keys %{$self->{mnet}{$mnet_num}{enzy}} ){
        my $reac_num = $self->{enzy}{$enzy_num}{reac};
        $LB{$reac_num} += $self->{enzy}{$enzy_num}{LB} eq 'NA'
                        ? $self->{mnet}{$mnet_num}{LB}
                        : $self->{enzy}{$enzy_num}{LB};
        $UB{$reac_num} += $self->{enzy}{$enzy_num}{UB} eq 'NA'
                        ? $self->{mnet}{$mnet_num}{UB}
                        : $self->{enzy}{$enzy_num}{UB};
    }
    foreach my $reac_num ( keys %{$self->{mnet}{$mnet_num}{reac}} ){
        $LB{$reac_num} = $LB if $LB{$reac_num} < $LB;
        $UB{$reac_num} = $UB if $UB{$reac_num} > $UB;
    }
    return ( \%LB, \%UB, $LB, $UB );
}

# Return true if the reaction is balanced with respect to mass and charge.
# Return false if the reaction is not balanced or if some information is missing
# The is_balanced column of the Coupound_Props.tsv file rely on a better algorithm
# able to deal with some generic compounds
sub is_balanced{
    my( $self, $reac_id ) = @_;
    my $reac_num = $self->_get_reac_num( $reac_id );
    my $delta_mass   = 0;
    my $delta_charge = 0;
    foreach my $spec_num ( keys %{$self->{reac}{$reac_num}{left}} ){
        my $chem_num = $self->{spec}{$spec_num}{chem};
        return 0 if $self->{chem}{$chem_num}{mass} eq '';
        $delta_mass   += $self->{reac}{$reac_num}{left}{$spec_num} * $self->{chem}{$chem_num}{mass};
        return 0 if $self->{chem}{$chem_num}{charge} eq '';
        $delta_charge += $self->{reac}{$reac_num}{left}{$spec_num} * $self->{chem}{$chem_num}{charge};
    }
    foreach my $spec_num ( keys %{$self->{reac}{$reac_num}{right}} ){
        my $chem_num = $self->{spec}{$spec_num}{chem};
        return 0 if $self->{chem}{$chem_num}{mass} eq '';
        $delta_mass   -= $self->{reac}{$reac_num}{right}{$spec_num} * $self->{chem}{$chem_num}{mass};
        return 0 if $self->{chem}{$chem_num}{charge} eq '';
        $delta_charge -= $self->{reac}{$reac_num}{right}{$spec_num} * $self->{chem}{$chem_num}{charge};
    }
    return 1 if abs( $delta_mass ) < 0.1 and $delta_charge == 0;
    return 0;
}
sub unify_reac_enzy{
    my( $self, $old_mnet_id, $new_mnet_id ) = @_;
    $self->add_mnet( $new_mnet_id, $self->get_LU_bounds( $old_mnet_id ));
    foreach my $reac_id ( $self->select_reac_ids( mnet => $old_mnet_id )){
        $self->copy_reac_add_enzy(
            $new_mnet_id,
            $reac_id,
            $self->merge_enzy_info( $self->get_enzy_info( $old_mnet_id, $reac_id )),
        );
        $self->set_reac_source( $new_mnet_id, $reac_id, $self->get_reac_source( $old_mnet_id, $reac_id ));
    }
}

sub one_dir_reac{ # split every metabolite in two (substrate and product, except BIOMASS)
                  # and make every reaction unidirectional
    my( $self, $old_mnet_id, $new_mnet_id ) = @_;
    my $tmp_mnet_id = '#' . $new_mnet_id;
    $self->unify_reac_enzy( $old_mnet_id, $tmp_mnet_id );
    my $tmp_mnet_num = $self->_get_mnet_num( $tmp_mnet_id );
    my( $LB, $UB ) = $self->get_LU_bounds( $tmp_mnet_id );
    $self->add_mnet( $new_mnet_id, $LB, $UB );
    my $new_mnet_num = $self->_get_mnet_num( $new_mnet_id );
    foreach my $reac_id ( $self->select_reac_ids( mnet => $tmp_mnet_id )){
        my( $complex, $lb, $ub, $dir ) = $self->get_enzy_info( $tmp_mnet_id, $reac_id );
        my $equation = $self->get_reac_equation( $reac_id );
        my $equa_id  = $self->get_reac_mnxr( $reac_id );
        my $source = '';
        if( $dir eq 'LR' or $dir eq 'B' ){
            my( $left, $right ) = split / <\?> /, $equation;
            $left  =~ s/\@/_S\@/g; # substrate
            $right =~ s/\@/_P\@/g; # product
            $left  =~ s/BIOMASS_S/BIOMASS/;
            $right =~ s/BIOMASS_P/BIOMASS/;
            $self->_store_reac(  $new_mnet_num, $reac_id . '_SP', $left . ' <?> ' . $right, $equa_id );
            $self->_store_enzy( $new_mnet_num, $reac_id . '_SP', $complex, ( $lb > 0 ? $lb : 0 ), $ub, 'LR' );
        }
        if( $dir eq 'RL' or $dir eq 'B' ){
             my( $left, $right ) = split / <\?> /, $equation;
            $left  =~ s/\@/_P\@/g; # substrate
            $right =~ s/\@/_S\@/g; # product
            $left  =~ s/BIOMASS_P/BIOMASS/;
            $right =~ s/BIOMASS_S/BIOMASS/;
            $self->_store_reac( $new_mnet_num, $reac_id . '_PS', $left . ' <?> ' . $right, $equa_id );
            $self->_store_enzy( $new_mnet_num, $reac_id . '_PS', $complex, $lb, ( $ub < 0 ? $ub : 0 ), 'RL' );
        }
    }
    foreach my $spec_id ( $self->select_spec_ids( mnet => $tmp_mnet_id )){
        my( $chem_id, $comp_id ) = split /\@/, $spec_id;
        next if $chem_id eq 'BIOMASS';
        next if $comp_id eq 'BOUNDARY';
        my $left  = "$chem_id\_P\@$comp_id";
        my $right = "$chem_id\_S\@$comp_id";
        my $reac_id = 'PS_' . $chem_id . '__64__' . $comp_id;
        $self->_store_reac(  $new_mnet_num, $reac_id, "1 $left <?> 1 $right", '' );
        $self->_store_enzy( $new_mnet_num, $reac_id, '', 0, $UB, 'LR' );
    }
    foreach my $chem_id ( $self->select_chem_ids( mnet => $tmp_mnet_id )){
        next if $chem_id eq 'BIOMASS';
        $self->set_chem_info( $chem_id .'_S' , $self->get_chem_info( $chem_id ));
        $self->set_chem_info( $chem_id .'_P' , $self->get_chem_info( $chem_id ));
    }
    $self->remove( $tmp_mnet_id );
}

sub shave_boundary{ # remove boudary reactions that are disconnected from the rest of the network
    my( $self, $mnet_id, $new_id ) = @_;
    $self->add_mnet( $new_id, $self->get_LU_bounds( $mnet_id ));
    my %skip = ();
    foreach my $reac_id ( $self->get_boundary_reac_ids( $mnet_id )){
        my @spec_id = $self->select_spec_ids( reac => $reac_id );
        foreach( @spec_id ){
            next if /\@BOUNDARY/;
            $skip{$reac_id} = 1 if $self->select_reac_ids( mnet => $mnet_id, spec => $_ ) == 1;
        }
    }
    my @reac_id = ();
    foreach( $self->select_reac_ids( mnet => $mnet_id )){
        push @reac_id, $_ unless $skip{$_};
    }
    $self->copy_reac_copy_enzy( $mnet_id, $new_id, @reac_id );
}
sub _get_mnet_dict{
    my( $self, @mnet_id ) = @_;
    my %mnet_ok = ();
    foreach my $mnet_id ( @mnet_id ){
        my $mnet_num = $self->_get_mnet_num( $mnet_id );
        $mnet_ok{$mnet_num} = $mnet_id;
    }
    return %mnet_ok;
}

sub get_protnet_view{ # FIXME: handle SPONTANEOUS
    my( $self, $max_freq, @mnet_id ) = @_;
    # First establish the list of mnet, reac and enzy that fit with @mnet_id
    my %mnet_num_ok = $self->_get_mnet_dict( @mnet_id );
    my %reac_num_ok = ();
    my %enzy_num_ok = ();
    my %complex     = ();

    foreach my $mnet_num ( keys %mnet_num_ok ){
        foreach my $reac_num ( keys %{$self->{mnet}{$mnet_num}{reac}} ){
            foreach my $enzy_num ( keys %{$self->{reac}{$reac_num}{enzy}} ){
                next unless exists $self->{enzy}{$enzy_num}{mnet}{$mnet_num};
                if( $self->{enzy}{$enzy_num}{complex} ){
                    foreach( split /;/, $self->{enzy}{$enzy_num}{complex} ){
                        $complex{$reac_num}{$_} = 1;
                    }
                }
            }
            $reac_num_ok{$reac_num} = 1 if exists $complex{$reac_num};
        }
    }

    # Secondly, compute the frequency of occurence of every chemical species
    my %spec_freq = ();
    foreach my $reac_num ( keys %reac_num_ok ){
        foreach my $spec_num ( keys %{$self->{reac}{$reac_num}{left}},
                               keys %{$self->{reac}{$reac_num}{right}} ){
            $spec_freq{$spec_num}++;
        }
    }

   # Now create the network

    my %netw = ();
    foreach my $reac_num ( keys %reac_num_ok ){
        my $reac_id = $self->{reac}{$reac_num}{id};
        foreach my $spec_num ( keys %{$self->{reac}{$reac_num}{left}},
                               keys %{$self->{reac}{$reac_num}{right}} ){
            next if $spec_freq{$spec_num} == 1;        # skip dead end
            next if $spec_freq{$spec_num} > $max_freq; # skip frequent species
            foreach my $reac_num_2 ( keys %reac_num_ok ){
                next if $reac_num == $reac_num_2;
                next unless exists $self->{reac}{$reac_num_2}{left}{$spec_num} or exists $self->{reac}{$reac_num_2}{right}{$spec_num};
                my $reac_id_2 = $self->{reac}{$reac_num_2}{id};
                if( exists $netw{$reac_id}{reac}{$reac_id_2} ){
                    if( $netw{$reac_id}{reac}{$reac_id_2} > $spec_freq{$spec_num} ){
                        $netw{$reac_id}{reac}{$reac_id_2} = $spec_freq{$spec_num};
                    }
                }
                else{
                    $netw{$reac_id}{reac}{$reac_id_2} = $spec_freq{$spec_num};
                }
            }
        }
    }

    # Now document the network

    foreach my $reac_id ( keys %netw ){
        my $reac_num = $self->{look}{reac}{$reac_id};
        $netw{$reac_id}{equa} = $self->{reac}{$reac_num}{equation};
        $netw{$reac_id}{enzy} = $complex{$reac_num};
        # my $dir = $self->get_reac_dir( $reac_id, @mnet_id );
        # $dir = 'B' unless $dir =~ /^(B|LR|RL)$/;
        # $netw{$reac_id}{dir}  = $dir;
    }
    return \%netw;
}

#sub DESTROY{
#    my $self = shift;
#    warn "DESTROY: $self\n";
#    $self->SUPER::DESTROY;
#}

sub filter_remove_reac{
    my( $self, $mnet_id_src, $mnet_id_dst, %arg ) = @_; # %arg is passed to select_reac_ids()
    $self->add_mnet( $mnet_id_dst, $self->get_LU_bounds( $mnet_id_src ));
    my %skip = ();
    $skip{$_} = 1 foreach $self->select_reac_ids( %arg );
    foreach( $self->select_reac_ids( mnet => $mnet_id_src )){
        next if $skip{$_};
        $self->copy_reac_copy_enzy( $mnet_id_src, $mnet_id_dst, $_ );
    }
}

sub get_pept2reac{
    my( $self, $mnet_id ) = @_;
    my %essential = ();
    foreach my $reac_id ( $self->select_reac_ids( mnet => $mnet_id )){
        my @enzy_info = $self->get_enzy_info( $mnet_id, $reac_id );
        if( @enzy_info == 4 ){
            my %count = ();
            my @prot = split ';', $enzy_info[0];
            foreach my $prot ( @prot ){
                foreach my $pept_id ( split /\+/, $prot ){
                    $count{$pept_id}++;
                }
            }
            foreach my $pept_id ( keys %count ){
                if( $count{$pept_id} == @prot ){
                    push @{$essential{$pept_id}}, $reac_id;
                }
            }
        }
        else{
            die "Reac/enzy not unified!";
        }
    }
    my %pept2reac = ();
    foreach my $pept_id ( $self->select_pept_ids( mnet => $mnet_id )){
        $pept2reac{$pept_id} = exists $essential{$pept_id} ? $essential{$pept_id} : [];
    }
    return \%pept2reac;
}

1;

