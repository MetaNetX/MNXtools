package Convert;

use strict;

use Exporter;
our @ISA = qw(Exporter);

use Data::Dumper;

use FindBin qw( $RealBin );
use lib "$RealBin/../util";
use Toolbox;
my $tb = Toolbox->new();

use Prefix;
use MetNet;
use LinkGroup;

sub new{
    my( $package, $chem_space, $gene_space ) = @_;
    my $self = {
        ns           => $chem_space,
        gene_space   => $gene_space,
        prefix_space => Prefix->new(),
        log4yaml   => [],
        chem_log   => {}, # helper data structure to build the log
        comp_log   => {}, # ditto
        reac_log   => {}, # ditto
        count      => {},
        chem_dict  => {},
        option => {
            true_PMF     => 0, # the default is to blurr the distinction between H(+) and PMF
            comp_premap  => {},
        },
    };
    return bless $self, $package;
}
sub set_comp_premap{
    my( $self, $comp_premap ) = @_;
    $self->{option}{comp_premap} = $comp_premap;
}
sub _use_xref{
    my( $self, $metnet, $source_name, $xref_regexp, $use_chem_id ) = @_;
    my %chem_map = ();
    my %xref_unk = ();
    my %chem_seen = ();
    $chem_seen{$self->{chem_dict}{$_}} = 1 foreach keys %{$self->{chem_dict}}; # already attributed
    foreach my $chem_old ( sort $metnet->select_chem_ids( mnet => $source_name )){
        next if exists $self->{chem_dict}{$chem_old}; # already assigned
        my $source_old = $metnet->get_chem_source( $source_name, $chem_old ) || $chem_old;
        my( $desc, $formula, $mass, $charge, $xref ) =  $metnet->get_chem_info( $chem_old );
        my @chem_unk = ();
        my @chem_id = $use_chem_id ? ( $chem_old ) : ();
        if( defined $xref_regexp ){
            foreach( split /;/, $xref ){
                push @chem_id, $_ if $_ =~ $xref_regexp;
            }
        }
        foreach my $chem_id ( @chem_id ){
            my ( $chem_new, $msg ) = $self->{ns}->map_chem( $chem_id );
            if( $chem_new =~ /^UNK:/ and $chem_id =~ /^(\w+):(\S+)$/){
                my( $dbkey, $id) = ( $1, $2 );
                my %buf = ();
                foreach my $prefix ( $self->{prefix_space}->get_prefix_from_depr( 'chem', $dbkey )){
                    ( $chem_new, $msg ) = $self->{ns}->map_chem( $prefix . ':' . $id );
                    $buf{$chem_new} = 1 if $chem_new !~ /^UNK:/;
                }
                my @buf = keys %buf;
                if( @buf == 1 ){
                    $chem_new = $buf[0];
                }
            }
            if( $chem_new =~ /^UNK:/ ){
                push @{$xref_unk{$chem_old}}, $chem_id;
            }
            else{
                push @{$chem_map{$chem_old}{$chem_new}}, $chem_id;
            }
        }
    }
    my @cluster = LinkGroup::cluster( \%chem_map );
    foreach my $cluster ( @cluster ){
        my @chem_old = keys %$cluster;
        my %val = ();
        foreach my $chem_old ( @chem_old ){
            $val{$_} = 1 foreach keys %{$cluster->{$chem_old}};
        }
        my @chem_new = keys %val;
        my $ok = 1;
        foreach( @chem_new ){
            $ok = 0 if exists $chem_seen{$_}
        }
        if( ! $ok ){
            foreach my $chem_old ( @chem_old ){
                $self->{chem_dict}{$chem_old} = 'UNK:' . $chem_old;
                push @{$self->{chem_log}{$chem_old}{status}},'- code: CHEM_XREF_DISCREPANCY', "  mappings:";
                foreach my $chem_new ( sort keys %{$chem_map{$chem_old}} ){
                    my %prop = $self->{ns}->get_chem_prop( $chem_new );
                    foreach my $xref ( sort @{$chem_map{$chem_old}{$chem_new}} ){
                        push @{$self->{chem_log}{$chem_old}{status}}, "    $xref: $chem_new # " . $prop{name};
                    }
                }
            }
        }
        elsif( @chem_new == 1 ){ # unambiguous one-to-one or many-to-one mapping
             foreach my $chem_old ( @chem_old ){
                $self->{chem_dict}{$chem_old} = $chem_new[0];
                my %prop = $self->{ns}->get_chem_prop( $chem_new[0] );
                if( defined $xref_regexp ){ # else do not report
                    push @{$self->{chem_log}{$chem_old}{status}},'- code: CHEM_XREF_OK', "  mappings:";
                    foreach my $xref ( sort @{$chem_map{$chem_old}{$chem_new[0]}} ){
                        push @{$self->{chem_log}{$chem_old}{status}}, "    $xref: $chem_new[0] # " . $prop{name};
                    }
                }
            }
        }
        else{ # @chem_new > 1 ambiguous or conflicting mappings
            my $parent = ''; # looks for incests, i.e. parent/child relations
            foreach my $chem_new ( @chem_new ){
                next if $parent;
                my %ok = ();
                $ok{$_} = 1 foreach $self->{ns}->search_chem_isom( $chem_new );
                $parent = $chem_new;
                foreach( @chem_new ){
                    $parent = '' unless exists $ok{$_};
                }
            }
            if( $parent ){ # one of chem is a parent of all the others
                foreach my $chem_old ( @chem_old ){
                    $self->{chem_dict}{$chem_old} = $parent;
                    push @{$self->{chem_log}{$chem_old}{status}}, '- code: CHEM_XREF_AMBIGUOUS', "  mappings:";
                    foreach my $chem_new ( sort keys %{$chem_map{$chem_old}} ){
                        my %prop = $self->{ns}->get_chem_prop( $chem_new );
                        foreach my $xref ( sort @{$chem_map{$chem_old}{$chem_new}} ){
                            push @{$self->{chem_log}{$chem_old}{status}}, "    $xref: $chem_new # " . $prop{name};
                        }
                    }
                }
            }
            else{ # create UNK: and report conflict
                foreach my $chem_old ( @chem_old ){
                    $self->{chem_dict}{$chem_old} = 'UNK:' . $chem_old;
                    push @{$self->{chem_log}{$chem_old}{status}}, '- code: CHEM_XREF_CONFLICT', "  mappings:";
                    foreach my $chem_new ( sort keys %{$chem_map{$chem_old}} ){
                        my %prop = $self->{ns}->get_chem_prop( $chem_new );
                        foreach my $xref ( sort @{$chem_map{$chem_old}{$chem_new}} ){
                            push @{$self->{chem_log}{$chem_old}{status}}, "    $xref: $chem_new # " . $prop{name};
                        }
                    }
                }
            }
        }
    }
}
sub _premap_chem{
    my( $self, $chem_id ) = @_;
    return exists $self->{chem_dict}{$chem_id} ? $self->{chem_dict}{$chem_id} : $chem_id;
}
sub _premap_map_equation{
    my( $self, $old_eq_str ) = @_;
    $old_eq_str .= ' '; # add trailing ' '
    my $new_eq_str = $old_eq_str;
    while( $old_eq_str =~ / (\S+)\@(\S+) /g ){
        my $spec_old = $1;
        my( $chem_id, $comp_id ) = ($1, $2 );
        $spec_old = $1 . '@' . $2;
        my $spec_new = $self->_premap_chem( $chem_id ) . '@' . $comp_id;
        $new_eq_str =~ s/ $spec_old / $spec_new /;
    }
    chop $new_eq_str; # remove trailing ' '
    return $self->{ns}->map_equation( $new_eq_str );
}
sub _interpret_chem_map_rule{
    my $chem_map_rule = shift;
    $chem_map_rule =~ s/[^\w,\.\-\/]//; # security
    return ( '', 1 ) unless $chem_map_rule; # ( no xref, use ID ) by default
    my @rule = ();
    foreach my $rule ( split /\//, $chem_map_rule ){
        my @xref     = ();
        my $use_xref = 0;
        my $use_id   = 0; # false
        foreach( split ',', $rule ){
            if( /^ID$/ ){
                 $use_id = 1;
            }
            elsif( /^XREF$/ ){
                $use_xref = 1;
            }
            else{
                push @xref, $_;
            }
        }
        if( $use_xref ){
            push @rule, [ qr/./, $use_id ];
        }
        elsif( @xref ){
            my $str = '^(' . join ( '|', @xref ) . ')';
            push @rule, [ qr/$str/, $use_id ];
        }
        else{
            push @rule, [ undef, $use_id ];
        }
    }
    return @rule;
}
sub convert{
    my( $self, $metnet, $source_name, $metnet2, $dest_name, $chem_map_rule, $option ) = @_;
    foreach( sort keys %$option ){
        $tb->die( "Invalid option key: $_" ) unless /^prefix|use_chem_ID|use_chem_xref|generic_comp$/;
    }
    $self->{log4yaml}   = []; # reset
    $self->{chem_dict}  = {};
    $self->{chem_log}   = {};
    $self->{comp_log}   = {};
    $self->{reac_log}   = {};
    my @chem_map_rule = _interpret_chem_map_rule( $chem_map_rule );
    $self->_use_xref( $metnet, $source_name, @{$_} ) foreach @chem_map_rule;
    my %reac_info = (); # Reac-centric temporary data structure to prepare the new model
    foreach my $reac_id ( sort $metnet->select_reac_ids( mnet => $source_name ) ){
        my $eq_orig = $metnet->get_reac_equation( $reac_id );
        $eq_orig    = $self->_premap_comp( $eq_orig ) if $self->{option}{comp_premap};
        my $source  = $metnet->get_reac_source( $source_name, $reac_id );
        my( $new_id, # might be '' if reaction is empty
            $eq_new,
            $sign,
            $is_balanced,
            $msg,
            $mnxr_id,
            $str_gen,
            $sign_gen ) = $self->_premap_map_equation( $eq_orig );
        if( $eq_new eq ' = ' or $mnxr_id eq 'EMPTY' ){
            my( $id_new, $msg ) = $self->{ns}->map_reac( $reac_id );
            if( $id_new eq 'EMPTY' ){
                $self->{reac_log}{$reac_id}{status} = [ '- code: REAC_EMPTY_MNXREF' ];
            }
            else{
                $self->{reac_log}{$reac_id}{status} = [ '- code: REAC_EMPTY_UNKNOWN' ];
            }
            next;
        }
        push @{$self->{reac_log}{$reac_id}{status}}, @$msg;
        $self->{reac_log}{$reac_id}{ID_dst} = $new_id;
        my $spec_loss = 0;
        foreach( @$msg ){
            $spec_loss = 1 if /REAC_MAP_LOSS/;
        }
        if( ! $spec_loss ){
            if( $mnxr_id ){
                push @{$self->{reac_log}{$reac_id}{status}}, '- code: REAC_MAP_MNXREF', "  mnxr: " . $mnxr_id;
            }
            elsif( $str_gen =~ /UNK:/){
                push @{$self->{reac_log}{$reac_id}{status}}, '- code: REAC_MAP_UNKNOWN';
            }
            else{
                push @{$self->{reac_log}{$reac_id}{status}}, '- code: REAC_MAP_OK';
            }
        }
        if( ! exists $reac_info{$new_id} ){
            $reac_info{$new_id} = {
                eq_orig   => $eq_orig, # to ease debugging
                equation  => $option->{generic_comp} ? $str_gen : $eq_new,
                mnxr      => $mnxr_id,
                source    => {},
                enzy      => [],
                from      => [ $reac_id ],
            };
            $reac_info{$new_id}{source}{$_} = 1 foreach split /;/, $source;
        }
        else{ # Equations merge: just need to update source and evid
            $reac_info{$new_id}{source}{$_} = 1 foreach split /;/, $source;
            push @{$reac_info{$new_id}{from}}, $reac_id;
        }
        $sign = $option->{generic_comp} ? $sign_gen : $sign;
        my @enzy_info = $metnet->get_enzy_info( $source_name, $reac_id );
        while( my( $complex_old, $lb, $ub, $dir ) = splice @enzy_info, 0, 4 ){
            my $complex_new = $self->{gene_space}->map_complex( $metnet, $complex_old );
            if( $sign == 1 ){
                push @{$reac_info{$new_id}{enzy}}, $complex_new, $lb, $ub, $dir;
            }
            else{
                my $lb_new  = $ub eq 'NA' ? 'NA'
                            : $ub == 0    ? 0
                            : - $ub;
                my $ub_new  = $lb eq 'NA' ? 'NA'
                            : $lb == 0    ? 0
                            : - $lb;
                my $dir_new = $dir eq 'LR' ? 'RL'
                            : $dir eq 'RL' ? 'LR'
                            : $dir; # case for 'B' or 'C';
                push @{$reac_info{$new_id}{enzy}}, $complex_new, $lb_new, $ub_new, $dir_new;
            }
        }
    }

    # ------------------------------------------------
    # Assembling the final model
    # ------------------------------------------------

    $metnet2->add_mnet( $dest_name, $metnet->get_LU_bounds( $source_name ));
    $metnet2->set_mnet_desc( $dest_name, $metnet->get_mnet_desc( $source_name ));
    $metnet2->push_mnet_info( $dest_name, $metnet->get_mnet_info( $source_name ));
    my %chem_ok = ();
    my %comp_ok = ();
    foreach my $new_id ( sort keys %reac_info ){
        eval{
            $metnet2->add_reac_add_enzy(
                $dest_name,
                $new_id,
                $reac_info{$new_id}{equation},
                $reac_info{$new_id}{mnxr},
                $metnet2->merge_enzy_info( @{$reac_info{$new_id}{enzy}} ), # will unify the resulting models
            );
        };
        if( $@ ){
            warn Dumper $reac_info{$new_id}; # prints out orig equation, which is essential to debug
            die $@;
        }

        my @info = $metnet->get_reac_info( $reac_info{$new_id}{from}[0] ); # $EC, $pathway, $xref
        @info    = $self->{ns}->get_reac_info( $reac_info{$new_id}{mnxr} ) if $reac_info{$new_id}{mnxr};
        $info[1] = ''; # pathway became a placeholder from now on (Oct 2020)
        $metnet2->set_reac_info( $new_id, @info );
        $metnet2->set_reac_source( $dest_name, $new_id, join ';', sort keys %{$reac_info{$new_id}{source}} );
        if( @{$reac_info{$new_id}{from}} > 1 ){
            foreach my $reac_old ( @{$reac_info{$new_id}{from}} ){
                push @{$self->{reac_log}{$reac_old}{status}},
                    '- code: REAC_MNET_MERGE',
                    '  IDs_src:',
                    map { '    - ' . $_ } sort @{$reac_info{$new_id}{from}};
            }
        }
    }

    my %chem_source  = ();
    my %chem_new2old = ();
    foreach my $chem_old ( $metnet->select_chem_ids( mnet => $source_name ) ){
        my $source_old = $metnet->get_chem_source( $source_name, $chem_old ) || $chem_old;
        my @info = $metnet->get_chem_info( $chem_old );
        my( $chem_new, $msg ) = $self->{ns}->map_chem( $self->_premap_chem( $chem_old )); # not the first call!
        $self->{chem_log}{$chem_old}{ID_dst}   = $chem_new;
        $self->{chem_log}{$chem_old}{name_src} = $info[0];
        $self->{chem_log}{$chem_old}{sources}  = $source_old;
        push @{$chem_new2old{$chem_new}}, $chem_old;
        $chem_source{$chem_new}{$_} = 1 foreach split /;/, $source_old;
#        if( $chem_new =~ /^UNK:/ ){
#            push @{$self->{chem_log}{$chem_old}{status}}, @$msg;
#        }
#        els
        if( @$msg > 0 ){
            push @{$self->{chem_log}{$chem_old}{status}}, @$msg;
        }
        else{ # this should NOT happen
            push @{$self->{chem_log}{$chem_old}{status}}, [ '- code: CHEM_MAP_FIXME' ];
        }
    }
    foreach my $chem_id ( $metnet2->select_chem_ids( mnet => $dest_name ) ){
        my @info = ();
        if( $chem_id =~ /^UNK:(.+)/ ){
            @info = eval{ $metnet->get_chem_info( $1 )};
            @info = eval{ $metnet->get_chem_info( $chem_id )} if $@;
        }
        else{
            @info = $self->{ns}->get_chem_info( $chem_id );
        }
        $info[0] = $chem_id unless $info[0];
        $metnet2->set_chem_info( $chem_id, @info );
        $metnet2->set_chem_source( $dest_name, $chem_id, join ';', sort keys %{$chem_source{$chem_id}} );
    }
    my %comp_source  = ();
    my %comp_new2old = ();
    if( $option->{generic_comp} ){
        foreach my $comp_old ( $metnet->select_comp_ids( mnet => $source_name ) ){
            $self->{comp_log}{$comp_old}{ID_dst} = '';
            push @{$self->{comp_log}{$comp_old}{status}}, '- code: COMP_GENERIC';
        }
    }
    else{
        foreach my $comp_old ( $metnet->select_comp_ids( mnet => $source_name ) ){
            my $source_old = $metnet->get_comp_source( $source_name, $comp_old ) || $comp_old;
            my @info = $metnet->get_comp_info( $comp_old );
            my( $comp_new, $msg ) = $self->{ns}->map_comp( $comp_old );
            $self->{comp_log}{$comp_old}{ID_dst} = $comp_new;
            $comp_source{$comp_new}{$_} = 1 foreach split /;/, $source_old;
            push @{$comp_new2old{$comp_new}}, $comp_old;
            push @{$self->{comp_log}{$comp_old}{status}}, @$msg;
        }
    }
    foreach my $comp_id ( $metnet2->select_comp_ids( mnet => $dest_name ) ){
        my @info = ();
        if( $self->{ns}->comp_exists( $comp_id )){
            @info = $self->{ns}->get_comp_info( $comp_id );
        }
        else{
            if( $comp_id =~ /^(UNK|NOMAP):(.+)/ ){
                @info = eval{ $metnet->get_comp_info( $2 )};
            }
            unless( $info[0] ){
                $info[0] = 'No description'; # don't leave the field empty for the sake of display
            }
        }
        $metnet2->set_comp_info( $comp_id, @info );
        $metnet2->set_comp_source( $dest_name, $comp_id, join ';', sort keys %{$comp_source{$comp_id}} );
    }

    foreach my $pept_id ( $metnet2->select_pept_ids( mnet => $dest_name ) ){
        my @info = ();
        if( $self->{gene_space}->exists_pept( $pept_id )){
            @info = $self->{gene_space}->get_pept_info( $pept_id );
        }
        else{
            eval{ @info = $metnet->get_pept_info( $pept_id ) };
            if( $@ and $pept_id =~ /^(UNK|NOMAP):(\S+)/ ){
                eval{ @info = $metnet->get_pept_info( $2) };
            }
        }
        $info[0] = $pept_id unless $info[0]; # don't leave this field empty for the sake of display
        $metnet2->set_pept_info( $pept_id, @info );
    }

    # ------------------------------------------------ #
    # Add blue and red protons exchange reactions, if
    # needed. BOUNDARY deserve a special treatment
    # ------------------------------------------------ #

    my @red_proton  = eval{ $metnet2->select_spec_ids( mnet => $source_name, chem => 'MNXM01') };
    my @blue_proton = eval{ $metnet2->select_spec_ids( mnet => $source_name, chem => 'MNXM1') };
    if( @red_proton > 0 and @blue_proton > 0 ){ # most likely case is when MNXM01 does not exist
        foreach my $red_proton ( @red_proton ){
            my $blue_proton = $red_proton;
            $blue_proton =~ s/^MNXM01\@/MNXM1\@/;
            my( $comp ) = $red_proton =~ /\@(\S+)/;
            next if $comp eq 'BOUNDARY';
            my( $new_id, $eq_new, $sign, $is_balanced, $msg, $mnxr_id, $str_gen, $sign_gen ) = $self->{ns}->map_equation( "1 $red_proton = 1 $blue_proton" );
            next if exists $reac_info{$new_id};
            next if $new_id eq 'mnxr01d2'; # mnxr01d1 is enough
            $metnet2->add_reac_add_enzy( $dest_name, $new_id, $eq_new, $mnxr_id, 'SPONTANEOUS', 'NA', 'NA', 'B' );
        }
    }

    # ------------------------------------------------ #
    # Test for isomeric incest
    # ------------------------------------------------ #

    my %parent  = ();
    %chem_ok = ();
    foreach( $metnet2->select_chem_ids( mnet => $dest_name )){
        $chem_ok{$_} = 1;
        $parent{$_}  = 1 if exists $self->{ns}{chem_isom}{$_}; # else ignore to speed up
    }
    foreach my $parent ( sort keys %parent ){
        my @isom  = $self->{ns}->search_chem_isom( $parent );
        my @child = ();
        foreach my $child ( @isom ){
            next unless $chem_ok{$child};
            next if $child eq $parent;
            push @child, $child;
        }
        if( @child ){
            foreach my $parent_old ( @{$chem_new2old{$parent}} ){
                my @member = ();
                foreach my $chem_new ( @child ){
                    foreach my $chem_old ( @{$chem_new2old{$chem_new}} ){
                        push @member, "    $chem_old: $chem_new # " . $metnet2->get_chem_desc( $chem_new );
                    }
                }
                push @{$self->{chem_log}{$parent_old}{status}},
                    '- code: CHEM_MNET_ISOMERIC',
                    '  children:',
                    sort @member;
            }
        }
    }

    # ------------------------------------------------ #
    # Detect merged chems
    # ------------------------------------------------ #

    my %buf = ();
    foreach my $chem_old ( keys %{$self->{chem_log}} ){
        push @{$buf{$self->{chem_log}{$chem_old}{ID_dst}}}, $chem_old;
    }
    foreach my $chem_new ( keys %buf ){
        my @chem_old = @{$buf{$chem_new}};
        if( @chem_old > 1 ){
            my @member = ();
            foreach( @chem_old ){
                push @member, "    - $_ # $self->{chem_log}{$_}{name_src}";
            }
            foreach( @chem_old ){
                push @{$self->{chem_log}{$_}{status}},
                    '- code: CHEM_MNET_MERGE',
                    '  IDs_src:',
                    sort @member;
            }
        }
    }

    # ------------------------------------------------ #
    # Assemble YAML conver log
    # ------------------------------------------------ #

    push @{$self->{log4yaml}}, '---', '', "\nchem:\n";
    foreach my $chem_old ( sort $metnet->select_chem_ids( mnet => $source_name )){
        my @info_old = $metnet->get_chem_info( $chem_old );
        my $chem_new = $self->{chem_log}{$chem_old}{ID_dst};
        my @info_new = eval{ $metnet2->get_chem_info( $chem_new )};
        @info_new = ( 'NOT FOUND FIXME', '', '', '' )  unless @info_new;
        push @{$self->{log4yaml}},
             '',  # spacer
             "  - ID_src: $chem_old # $self->{chem_log}{$chem_old}{name_src}",
             "    ID_dst: $chem_new # $info_new[0]",
             '    status:',
             '      ' . join( "\n      ", @{$self->{chem_log}{$chem_old}{status}} );
    }
    push @{$self->{log4yaml}}, "\ncomp:";
    foreach my $comp_old ( sort $metnet->select_comp_ids( mnet => $source_name )){
        my @info_old = $metnet->get_comp_info( $comp_old );
        my $comp_new = $self->{comp_log}{$comp_old}{ID_dst} || '~';
        my $desc = $option->{generic_comp}
                 ? 'MNXD1 or MNXD2'
                 : ( $metnet2->get_comp_info( $comp_new ))[0];
        push @{$self->{log4yaml}},
             '',  # spacer
             "  - ID_src: $comp_old # " . ( $metnet->get_comp_info( $comp_old ))[0],
             "    ID_dst: $comp_new # $desc",
             '    status:',
             '      ' . join( "\n      ", @{$self->{comp_log}{$comp_old}{status}} );
    }
    push @{$self->{log4yaml}}, "\nreac:";
    foreach my $reac_old ( sort $metnet->select_reac_ids( mnet => $source_name )){
        my $reac_new = $self->{reac_log}{$reac_old}{ID_dst} || '~';
        my $desc     = $self->{reac_log}{$reac_old}{ID_dst}
                     ? $metnet2->get_reac_equation( $reac_new )
                     : 'EMPTY reaction';
        push @{$self->{log4yaml}},
            '',  # spacer
            "  - ID_src: $reac_old # " . $metnet->get_reac_equation( $reac_old ),
            "    ID_dst: $reac_new # $desc",
            '    status:',
            '      ' . join( "\n      ", @{$self->{reac_log}{$reac_old}{status}} );
    }
}
sub _premap_comp{
    my( $self, $str_old ) = @_;
    $str_old .= ' ';
    my $str_new = $str_old;
    my %seen = ();
    while( $str_old =~ /\@(\S+) /g ){
        next if $seen{$1};
        if( exists $self->{option}{comp_premap}{$1} ){
            my $old_comp = $1;
            $seen{$old_comp} = 1;
            my $new_comp = $self->{option}{comp_premap}{$old_comp};
            $str_new =~ s/\@$old_comp /\@$new_comp\@ /g;
        }
    }
    $str_new =~ s/\@ / /g;
    if( $str_old ne $str_new ){
        warn "$str_old\n$str_new\n\n";
    }
    chop $str_new;
    return $str_new;
}
sub yaml_escape{
    my $str = shift;
    $str =~ s/'/\\'/g;
    return $str;
}

1;

