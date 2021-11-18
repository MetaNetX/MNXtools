package Toolbox;

##  https://tldp.org/LDP/abs/html/exitcodes.html

use strict;

use Carp qw( longmess );
use Data::Dumper;
use IO::File;
use IO::Dir;
use Term::ANSIColor;

sub new{
    my( $package, $is_quiet ) = @_;
    my $self = bless {
        verbose => ! $is_quiet,
        log     => undef,
    }, $package;
    # $self->{verbose} = ! $is_quiet; # verbose by default: report actions on STDERR
    # $self->{indent}  = $ENV{MNX_INDENT_COUNTER} || 0;
    return $self;
}
sub verbose{ # turn verbosity on/off and query it
    my( $self, $verbose ) = @_;
    $self->{verbose} = $verbose  if defined $verbose;
    return $self->{verbose};
}
sub set_log{
    my( $self, $log_fh ) = @_;
    $self->{log} = $log_fh;
}
sub _format_msg{
    my( $self, $title, $arg, $color ) = @_;
    $arg = Dumper $arg if ref $arg;
    $arg = '' unless defined $arg;
    $color = 'black' unless $color;
    chomp $arg;
    $title = ' ' . $title  while 12 > length $title;
    my $txt = '#'; #  . '#' x ($ENV{'MNX_INDENT_COUNTER'} || 2);
    print {$self->{log}} $title . ' : ' . $arg . "\n" if  $self->{log};
    return join '',
        colored( "$txt $title : ", 'blue' ),
        colored( $arg, $color ),
        "\n";
}
sub report{ # the "two-arguments version of warn"
    my( $self, $title, $arg, $color) = @_; # $title is expected to be a single short word
    warn $self->_format_msg( $title, $arg, $color || 'black' );
}
sub warn{
    my( $self, $arg ) = @_;
    $self->report( 'WARN', $arg, 'red');
}
sub die{
    my( $self, $arg ) = @_;
    die $self->_format_msg( 'ERROR', $arg, 'red');
    exit 1;
}
sub confess{
    my( $self, $arg ) = @_;
    CORE::die $self->_format_msg( 'ERROR', $arg . longmess(), 'red');
}
sub try_system{ # returns ( <exit code>, <message> )
    my( $self, $cmd ) = @_;
    $self->report( 'system', $cmd, 'magenta' ) if $self->{verbose};
    if( system( $cmd ) == 0 ){
        return( 0, '' );
    }
    else{
        my $exit = $? >> 8;
        my $msg  = $!;
        $self->report( 'FAILED', "EXIT: $exit ($msg)", 'red' );
        return( $exit, $msg );
    }
}
sub system{
    my( $self, $cmd ) = @_;

    $self->report( 'system', $cmd, 'magenta' ) if $self->{verbose};
    $ENV{MNX_INDENT_COUNTER} += 1;
    unless( system( $cmd ) == 0 ){
        $self->warn( $cmd ) unless $self->{verbose};
        # the following code was adapted from `perldoc -f system`
        if( $? == -1 ){
            $self->die( "Failed to execute: $!" );
        }
        elsif( $? & 127 ) {
            $self->die( sprintf
                'Child died with signal %d, %s coredump',
                ($? & 127),  ($? & 128) ? 'with' : 'without' );
        }
        else{
           $self->die( sprintf 'Child exited with value %d', $? >> 8 );
           # $self->warn( sprintf 'Child exited with value %d', $? );
        }
    }
    $ENV{MNX_INDENT_COUNTER} -= 1;
}
sub open{
    my( $self, $file_arg ) = @_; # something like 'foo.txt'; '> foot.txt'; '| sort >> foo.txt'
    my $fh = new IO::File;
    $fh->open( $file_arg ) or $self->die( "Cannot open: $file_arg" );
    $self->report( 'open', $file_arg ) if $self->{'verbose'};
    return $fh;
}
sub slurp{
    my( $self, $filename ) = @_;
    my $fh = new IO::File;
    $fh->open( $filename )  or $self->die( "Cannot slurp file: $filename" );
    $self->report( 'slurp', $filename ) if $self->{'verbose'};
    my $buf = do{ local $/; <$fh>; }; # Faster file slurping method I am aware of
    close $fh;
    return $buf;
}
sub scan_tsv{ # somehow a memory greedy method: use only on "small" file
    my( $self, $filename, $delim ) = @_;
    $delim = "\t"  unless $delim;
    my $fh = new IO::File; # !!! Some operating systems may perform IO::File::new() or IO::File::open() on a directory without errors. (from perldoc IO::File )
    $fh->open( $filename ) or $self->die( "Cannot scan_tsv file: $filename" );
    $self->report( 'scan_tsv', $filename ) if $self->{'verbose'};
    my $buf = do{ local $/; <$fh>; }; # Maybe sysread is even faster!
    close $fh;
    chomp $buf;
    my @tab = ();
    foreach( split /\n/, $buf ){
        next if /^\#/;
        push @tab, [ split /$delim/o, $_, -1 ]; # arg -1 permits to catch the rightmost empty fields ( perldoc -f split )
    }
    return \@tab;
}
sub scan_ssv{ # same as above, much more tolerant to manual edition
    my( $self, $filename, $count ) = @_;
    my $fh = $self->open( $filename );
    my @tab = ();
    while( my $line = <$fh> ){
        chomp $line;
        $line =~ s/\#.*//;
        next if $line =~ /^\s*$/;
        $line =~ s/^\s+//;
        my @token = split /\s+/, $line;
        if( @token < $count ){
            $self->warn( 'Ignoring token: ' . join ' ', @token );
        }
        else{
            push @tab, [ splice @token, 0, $count ];
        }
    }
    return \@tab;
}

sub ls{ # return the file names found at $path
    my( $self, $path, $regexp ) = @_;
    my $dh = IO::Dir->new( $path ) or $self->die( "Invalid dir: $path" );
    my @filename = ();
    $regexp = '^[^\.]' unless $regexp; # to ignore hidden files by default
    while ( defined ( $_ = $dh->read )){
        push @filename, $_ if /$regexp/;
    }
    return @filename;
}

# -----------------------------------------------------------------
#
# These methods return new object
#
# -----------------------------------------------------------------

sub new_R_bridge{ # provides the exec and stop methods
    my $self = shift;
    return Toolbox::RBrigde->new( $self, @_ );
}
sub new_ssh_bridge{ # provides the exec and stop methods
    my( $self, $hostname ) = @_;
    return Toolbox::SshBridge->new( $self, $hostname);
}

1;

# ---------------------------------------------------- #
# R-interactive session
# ---------------------------------------------------- #

package Toolbox::RBrigde;

use strict;

use IPC::Run qw( start pump finish );

sub new{
    my( $package, $toolbox ) = @_;
    my $self = bless { toolbox => $toolbox }, $package;
    my $R_launch_cmd = 'R --slave --interactive'; # FIXME: don't hardcode this !!!
    $toolbox->report('start R', $R_launch_cmd ) if $self->{toolbox}->verbose();
    $self->{R_in}       = '';
    $self->{R_out}      = '';
    $self->{R_err}      = '';
    $self->{R_harness} = start [ split ' ', $R_launch_cmd ],
                                 '<' ,    \$self->{R_in},
                                 '>pty>', \$self->{R_out},
                                 '2>',    \$self->{R_err} # perldoc IPC::Run for details
        or die "Cannot launch R";
    return $self;
}
sub exec{
    my( $self, $cmd, $no_try_catch ) = @_;
    $self->{toolbox}->report( 'R', $cmd, 'green') if $self->{toolbox}->verbose();
    my $DEBUG = 0;
    unless( $no_try_catch ){
#        $cmd = "tryCatch({ " . $cmd . " },
#    error=function(e){ simpleError(e) },
#    warning=function(w){ message(w) }
#)";
   $cmd = "tryCatch({ " . $cmd . " },
    error=function(e){ simpleError(e) }
)";
    }
    $cmd .= "\n"  unless $cmd =~ /\n$/;
    $self->{toolbox}->report( 'R CMD', $cmd )  if $DEBUG;
    my $line_count = ( $cmd =~ tr/\n/\n/ );
    $self->{R_log_fh}->print($cmd)  if exists $self->{R_log_fh};
    $self->{R_in}  = $cmd . "print('__DONE__')\n";
    $self->{R_out} = '';
    $self->{R_err} = '';
    pump $self->{R_harness}  until $self->{R_out} =~ /"__DONE__"\r\n/;
    $self->{toolbox}->report( 'R OUT', $self->{R_out} ) if $DEBUG;
    $self->{toolbox}->report( 'R msg', $self->{R_err} ) if $self->{R_err};
    if( my($msg) = $self->{R_out} =~ /\<simpleError: (.+)/){
        $self->{R_log_fh}->close() if exists $self->{R_log_fh};
        $self->{toolbox}->report( 'Faulty R', $cmd );
        $self->{toolbox}->die( $msg );
    }
    my @out = split /\r\n/, $self->{R_out};
    splice @out, 0, $line_count;       # remove the echoed command
    pop @out; pop @out;   # remove the 2 lines with '__DONE__' at the end
    return join( "\n", @out ). "\n";
}
sub stop{
    my $self = shift;
    $self->{toolbox}->report('stop R', 'q()')  if $self->{toolbox}->verbose();
    $self->{R_in}  ="q()\n";
    pump $self->{R_harness};
    finish $self->{R_harness}  or die "Cannot finish: $?\n";
    if ( exists $self->{R_log_fh} ){
        $self->{R_log_fh}->close();
        delete $self->{R_log_fh};
    }
    delete $self->{$_}  foreach ('R_in', 'R_out', 'R_harness');
}

1;

package Toolbox::SshBridge;

use strict;

use IPC::Run qw( start pump finish );

sub new{
    my($package, $toolbox, $hostname) = @_;
    my $self = bless { toolbox    => $toolbox,
                       hostname => $hostname }, $package;
    $self->{toolbox}->report('start',"ssh -T $hostname") if $self->{toolbox}{verbose};
    $self->{SSH_in}       = '';
    $self->{SSH_out}      = '';
    $self->{SSH_err}      = '';
    $self->{SSH_harness} = start [ 'ssh','-T',$hostname ],
                                   '<' , \$self->{SSH_in},
                                   '>',  \$self->{SSH_out},
                                   '2>', \$self->{SSH_err}, or die "Cannot launch SSH";
    return $self,
}
sub exec{
    my( $self, $cmd, $no_echo ) = @_;
    $self->{toolbox}->report( 'ssh exec' , $cmd, 'magenta' ) if $self->{toolbox}{verbose} and ! $no_echo;
    $cmd .= "\n" unless $cmd =~ /\n$/;
    $self->{SSH_out} = '';
    $self->{SSH_err} = '';
    $self->{SSH_in}  = $cmd . "echo __DONE__\n";
    eval {
        pump $self->{SSH_harness} until $self->{SSH_out} =~ /__DONE__\n/;
    };
    $self->{SSH_out} =~ s/(.+)\n$//;     # remove the last line;
    return ($self->{SSH_out},$self->{SSH_err});
}
sub stop{
    my $self = shift;
    $self->{toolbox}->report('ssh stop');
    $self->{SSH_in}  ="exit\n";
    pump $self->{SSH_harness};
    finish $self->{SSH_harness} or die "Cannot finish: $?\n";
    delete $self->{$_} foreach ('SSH_in','SSH_out','SSH_err','SSH_harness');
}

1;
