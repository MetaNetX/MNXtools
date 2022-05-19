package Constants;
#Filename Constants.pm

use strict;
use warnings;
use diagnostics;

use Log::Log4perl qw(:easy);

#Log::Log4perl->easy_init($WARN);
Log::Log4perl->init(\ <<'EOT');
  log4perl.category = WARN, Screen
  log4perl.appender.Screen = \
      Log::Log4perl::Appender::ScreenColoredLevels
  log4perl.appender.Screen.layout = \
      Log::Log4perl::Layout::PatternLayout
  log4perl.appender.Screen.layout.ConversionPattern = \
      %d %F{1} %L> %m %n
EOT
our $logger = get_logger();


my $debug = 0;
my $current_year = (localtime())[5] + 1900;
our $copyright = "MetaNetX copyright 2011-$current_year SystemsX, SIB Swiss Institute of Bioinformatics";

sub set_debug {
    my ($main_debug) = @_;

    if ( $main_debug ){
        $debug = $main_debug;
        Log::Log4perl->init(\ <<'EOTC');
  log4perl.category = DEBUG, Screen
  log4perl.appender.Screen = \
      Log::Log4perl::Appender::ScreenColoredLevels
  log4perl.appender.Screen.layout = \
      Log::Log4perl::Layout::PatternLayout
  log4perl.appender.Screen.layout.ConversionPattern = \
      %d %F{1} %L> %m %n
EOTC
    }

    return;
}

sub get_debug {
    return $debug;
}

# MetaNetX TSV file names
our $MetaNetX_model_filename       = 'model.tsv';
our $MetaNetX_compartment_filename = 'compartments.tsv';
our $MetaNetX_chemical_filename    = 'chemicals.tsv';
our $MetaNetX_reaction_filename    = 'reactions.tsv';
our $MetaNetX_enzyme_filename      = 'enzymes.tsv';
our $MetaNetX_peptide_filename     = 'peptides.tsv';


# SBML
##Default bounds
our $MIN_BOUND_VAL  = -1000;
our $MAX_BOUND_VAL  =  1000;
our $ZERO_BOUND_VAL =  0;
our $MIN_BOUND_ID   = 'cobra_default_lb';
our $MAX_BOUND_ID   = 'cobra_default_ub';
our $ZERO_BOUND_ID  = 'cobra_0_bound';

##Default SBML version and level
our $SBML_level   = 3; # Possible values: 1.1 2.1 2.2 2.3 2.4 2.5 3.1 3.2
our $SBML_version = 1;
##Default FBC version
our $FBC_version  = 2; # Possible values: 1 2
##SBO
our $SBO_FLUXBALANCEFRMW = 624; # eq SBO:0000624
our $SBO_FLUXDEFAULT     = 626;
our $SBO_FLUX            = 625;
##Identifiers
my $identifier = 'https://identifiers.org/';
our $identifiers_taxid = $identifier.'taxonomy:';
our $identifiers_mnxc  = $identifier.'metanetx.compartment:';
our $identifiers_go    = $identifier.''; #expect e.g. GO:0006915
our $identifiers_biggc = $identifier.'bigg.compartment:';
##Unit
our $unit = 'mmol_per_gDW_per_hr';


# MetaNetX/MNXref
##Private variables/names
our $boundary_comp_id     = 'BOUNDARY';
our $boundary_comp_name   = 'model boundary';
our $boundary_comp_sbo    = 'SBO:0000289';
our $default_comp_sbo     = 'SBO:0000290';

our $default_chem_sbo     = 'SBO:0000247';

our $transport_reac_sbo   = 'SBO:0000185';

our $biomass_chem_id      = 'BIOMASS';
our $biomass_reac_id      = 'BIOMASS_EXT';
our @biomass_alias        = qw( biomass M_biomass cpd11416 cpd29751 S_Biomass_biomass Biomass_biomass s_0450 );
our @are_not_biomass_reac = qw( PRISM_growth_room rxn12877 rxn13038 );
our $biomass_chem_sbo     = 'SBO:0000649'; #NOTE is most of the time used for ALL growth equation components, not for BIOMASS chemical species itself
our $biomass_reac_sbo     = 'SBO:0000629';

our $spontaneous_enz      = 'SPONTANEOUS';

1;

