package Analysis;

use strict;

# The constructor is supposed to make all checks to 'accept' a job
# and must call $mnx->user_error() if something goes wrong.

# Later when launch( $mnx ) is called no user_error() should
# be produced

# The update_analysis_status() sub is essential for piloting
# the analysis

sub html_form{ # static, inherited by "simple" method
    my( $mnx, $metnet, $anal_id ) = @_;
    return join "\n",
        '<h1>' . $mnx->get_title( $anal_id ) . '</h1>',
        '<div class="center">',
        "<table class='center'><tr><td class='form_box'>",
        "<form action='/cgi-bin/mnxweb/launch' id='form_id' method='post' accept-charset='utf-8'>",
        "<input type='hidden' name='analysis' value='$anal_id'>",
        "<div class='tcenter'>Select at least one model below &nbsp;",
        "<button type='button' onclick='fade2grey(\"form_id\");'>Analyze</button>",
        '</div>',
        $mnx->get_mnet_summary( $metnet, 'checkbox', anal => $anal_id ),
        "</form>",
        "</td></tr></table></div>",
        '<br>',
        "<table class='mnx_table' style='width:80%'>",
        '<caption>Instructions</caption>',
        "<tr><th style='width:20%'>What</th><th>Comment</th></tr>",
        "<tr><td>Principle</td><td class='justif'>" . $mnx->get_desc( $anal_id ) . "</td></tr>",
        '</table>';
}

sub html_form_2{
    my( $mnx, $metnet, $anal_id, $title, $desc ) = @_;
    return join "\n",
        '<h1>' . $title . '</h1>',
        '<div class="center">',
        "<table class='center'><tr><td class='form_box'>",
        "<form action='/cgi-bin/mnxweb/launch' id='form_id' method='post' accept-charset='utf-8'>",
        "<input type='hidden' name='analysis' value='$anal_id'>",
        "<div class='tcenter'>Select at least one model below &nbsp;",
        "<button type='button' onclick='fade2grey(\"form_id\");'>Analyze</button>",
        '</div>',
        $mnx->get_mnet_summary( $metnet, 'checkbox', anal => $anal_id ),
        "</form>",
        "</td></tr></table></div>",
        '<br>',
        "<table class='mnx_table' style='width:80%'>",
        '<caption>Instructions</caption>',
        "<tr><th style='width:20%'>What</th><th>Comment</th></tr>",
        "<tr><td>Principle</td><td class='justif'>" . $desc . "</td></tr>",
        '</table>';
}

sub new{
    my( $package, $anal_id ) = @_;
    my $self = bless { # not a good idea to store $mnx here!
        id       => $anal_id,
        task_id  => '',     # used to monitor job execution
        status   => 'INIT', # one of INIT, RUN, DONE, FAILED
        props    => [],     # to store the results of the analysis
    }, $package;
    return $self;
}
sub get_id{
    return shift->{id};
}
sub get_status{
    return shift->{status};
}
sub set_status{
    my( $self, $status ) = @_;
    $self->{status} = $status;
}
sub get_task_id{
    shift->{task_id};
}
sub get_chunk_count{
    my $self = shift;
    return 1; # default
}
sub lauch{ # must return the updated status
    my( $self, $mnx ) = @_;
    die "Abstract method Analysis::launch() called\n";
}
sub collect{
    my( $self, $mnx ) = @_;
    die "Abstract method Analysis::collect() called\n";
}
sub get_error_msg{
    return $_[0]->{error_msg} || 'unknown error';
}
sub get_html_summary{
    my( $self, $mnx, $metnet, $mnet_id ) = @_;
    my $anal_id = $self->{id};
    my $title   = $mnx->get_title( $anal_id );
    my $desc    = $mnx->get_desc( $anal_id );
    my @table = ();
    foreach my $prop ( @{$self->{props}} ){
        push @table, $prop->get_html_table( $mnx, $mnet_id );
    }
    my @html = (
        "<h1>$title of <em>$mnet_id</em></h1>",
        '<div>',
        "<table class='mnx_table' style='width:80%'>",
        "<tr><th style='width:50%'>Result</th><th>Comment</th></tr>",
        "<tr><td class='center'>",
        join( '<br>', @table ),
        "</td>\n<td class='justif'>$desc</td></tr>\n",
        "</table></div>\n",
        "<p class='notification'>",
        $mnx->get_title( 'results' ),
        "</p>\n",
    );
    return join "\n", @html;
}
sub get_html_error_report{
    my( $self, $mnx, $metnet, $mnet_id ) = @_;
    my $anal_id = $self->{id};
    my @table = ();
    my @html = (
        "<h1>$anal_id of $mnet_id failed: <span style='color:red'>unknown error</h1>"
    );
    return join "\n", @html;
}
sub set_props{
    my( $self, $props ) = @_;
    $self->{props} = $props;
}
sub get_props{
    return @{shift->{props}};
}

1;

