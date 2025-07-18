package Prefix;

use strict;

use Data::Dumper;

use Carp;

# scope: /^(chem|reac|comp|pept|other)$/
# value: not empty base IRI (mandatory)
# ident: prefix from identifiers.org (optional )
#        NB: identifiers.org may be used in the primary prefix value, e.g. hmdb
# depr:  deprecated prefixes (optional)

my %prefix_data =(
    ### RDF base ###
    xsd => {
        scope => 'other',
        value => 'http://www.w3.org/2001/XMLSchema#',
    },
    rdf => {
        scope => 'other',
        value => 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    },
    rdfs => {
        scope => 'other',
        value => 'http://www.w3.org/2000/01/rdf-schema#',
    },
    owl => {
        scope => 'other',
        value => 'http://www.w3.org/2002/07/owl#',
    },
    sh => {
        scope => 'other',
        value => 'http://www.w3.org/ns/shacl#',
    },
    ### MetaNetX ###
    comp => {
        scope => 'comp',
        value => 'https://rdf.metanetx.org/comp/',
        ident => 'metanetx.compartment',
    },
    chem => {
        scope => 'chem',
        value => 'https://rdf.metanetx.org/chem/',
        ident => 'metanetx.chemical',
    },
    reac => {
        scope => 'reac',
        value => 'https://rdf.metanetx.org/reac/',
        ident => 'metanetx.reaction',
    },
    pept => {
        scope => 'pept',
        value => 'https://rdf.metanetx.org/pept/',
    },
    mnet => {
        scope => 'other',
        value => 'https://rdf.metanetx.org/mnet/',
    },
    mnx => {
        scope => 'other',
        value => 'https://rdf.metanetx.org/schema/',    # for predicates and special symbols like BIOMASS, PMF or SPONTANEOUS
    },
    ### Structural formats ###
#    SMILES => {
#        scope => 'chem',
#        value => ???,
#    },
    inchi => {
        scope => 'chem',
        value => 'https://identifiers.org/inchi:',
        ident => 'inchi',
    },
    inchikey => {
        scope => 'chem',
        value => 'https://identifiers.org/inchikey:',
        ident => 'inchikey',
    },
    ### BiGG ###
    biggM => {
        scope => 'chem',
        value => 'http://bigg.ucsd.edu/universal/metabolites/',
        ident => 'bigg.metabolite',
        depr  => [ 'bigg' ],
    },
    biggC => {
        scope => 'comp',
        value => 'http://bigg.ucsd.edu/compartments/',
        ident => 'bigg.compartment',
        depr  => [ 'bigg' ],
    },
    biggR => {
        scope => 'reac',
        value => 'http://bigg.ucsd.edu/universal/reactions/',
        ident => 'bigg.reaction',
        depr  => [ 'bigg' ],
    },
    ### VMH ###
    vmhM => {
        scope => 'chem',
        value => 'https://vmh.life/#metabolite/',
        ident => 'vmhmetabolite',
    },
    vmhR => {
        scope => 'reac',
        value => 'https://vmh.life/#reaction/',
        ident => 'vmhreaction'
    },
    vmhG => {
        scope => 'pept',
        value => 'https://vmh.life/#gene/',
        ident => 'vmhgene',
    },
#    vmhC => {
#        scope => 'comp',
#        value => 'http://bigg.ucsd.edu/compartments/',
#    },
    ### ChEBI (NB: for MetaNetX, CHEBI[:_] is not part of the identifier )
    chebi => {
        scope => 'chem',
        value => 'http://purl.obolibrary.org/obo/CHEBI_',
        ident => 'CHEBI',
    },
    ### HMDB ###
    hmdb => {
        scope => 'chem',
        value => 'https://identifiers.org/hmdb:',
        ident => 'hmdb',
    },
    # KEGG (the orignial IRI are not very consistent!)
    keggC => {
        scope => 'chem',
        value => 'https://www.genome.jp/kegg/',
        ident => 'kegg.compound',
        depr  => [ 'kegg' ],
    },
    keggD => {
        scope => 'chem',
        value => 'https://www.genome.jp/kegg/drug/',
        ident => 'kegg.drug',
        depr  => [ 'kegg' ],
    },
    keggG => {
        scope => 'chem',
        value => 'https://www.genome.jp/kegg/glycan/',
        ident => 'kegg.glycan',
        depr  => [ 'kegg' ],
    },
    keggE => {
        scope => 'chem',
        value => 'https://www.genome.jp/kegg/',
        ident => 'kegg.environ',
        depr  => [ 'kegg' ],
    },
    keggR => {
        scope => 'reac',
        value => 'https://www.genome.jp/kegg/reaction/',
        ident => 'kegg.reaction',
        depr  => [ 'kegg' ],
    },
    ### SEED ###
    seedM => {
        scope => 'chem',
        value => 'https://modelseed.org/biochem/compounds/',
        ident => 'seed.compound',
        depr  => [ 'seed' ],
    },
    seedR => {
        scope => 'reac',
        value => 'https://modelseed.org/biochem/reactions/',
        ident => 'seed.reaction',
        depr  => [ 'seed' ],
    },
    seedC => {
        scope => 'comp',
        value => 'https://a_placeholder_base_IRI_for_modelseed_compartments/', # invented IRI as far as I know !!!
    },
    ### EnviPath ###
    envipathM => {
        scope => 'chem',
        value => 'https://envipath.org/package/',
        ident  => 'envipath',
    },
    ### LIPIDMAPS ###
    lipidmapsM => {
        scope => 'chem',
        value => 'http://www.lipidmaps.org/data/LMSDRecord.php?LMID=',
        ident => 'lipidmaps',
    },
    ### MetaCyC ###
    metacycM => {
        scope => 'chem',
        value => 'https://biocyc.org/compound?id=',
        ident => 'metacyc.compound',
        depr  => [ 'metacyc' ],
    },
    metacycR => {
        scope => 'reac',
        value => 'https://metacyc.org/META/NEW-IMAGE?type=REACTION&object=',
        ident => 'metacyc.reaction',
        depr  => [ 'metacyc' ],
    },
    cco => {
        scope => 'comp',
        value => 'http://brg.ai.sri.com/CCO/downloads/cco.html#',
    },
    ### REACTOME ### ( no distinction between the chem/reac/comp )
    reactomeM => {
        scope => 'chem',
        value => 'https://reactome.org/content/detail/',
        ident => 'reactome',
        depr  => [ 'reactome' ],
    },
    reactomeR => {
        scope => 'reac',
        value => 'https://reactome.org/content/detail/', # same as above !
        ident => 'reactome', # reactome and identifiers.org utilise the same prefix for metabolites and reactions.
                             # Duplicated prefix declarations may cause problem in some implementation of SPARQL engine.
        depr  => [ 'reactome' ],
    },
    ### Sabio-RK ###
    sabiorkM => {
        scope => 'chem',
        value => 'https://sabiork.h-its.org/newSearch?q=',
        ident => 'sabiork.compound',
        depr  => [ 'sabiork' ],
    },
    sabiorkR => {
        scope => 'reac',
        value => 'https://sabiork.h-its.org/newSearch?q=sabioreactionid:',
        ident => 'sabiork.reaction',
        depr  => [ 'sabiork' ],
    },
    ### SwissLipids
    slm => {
        scope => 'chem',
        value => 'https://swisslipids.org/rdf/SLM_',
        ident => 'SLM',
    },
    glycosphingo => {
        scope => 'chem',
        value => 'http://slm.example.org/glycosphingo/', # invented IRI as far as I know !!!
    },
    rh => {
        scope => 'reac', # this is the dominant scope, but it is also used for compartment and vocabulary!
        value => 'http://rdf.rhea-db.org/',
        ident => 'rhea',
    },
    rheaP => {
        scope => 'chem',
        value => 'http://rhea.example.org/polymer/', # invented IRI as far as I know !!!
    },
    rheaG => {
        scope => 'chem',
        value => 'http://rhea.example.org/generic/', # invented IRI as far as I know !!!
    },
    rheaC => {
        scope => 'comp',
        value => 'http://rhea.example.org/compartment/', # invented IRI as far as I know !!!
    },
    elr => {
        scope => 'reac',
        value => 'http://rhea.example.org/elr/', # invented IRI as far as I know !!!
    },
    ### UniProt ###
    uniprotkb => {
        scope => 'pept',
        value => 'http://purl.uniprot.org/uniprot/',
        ident => 'uniprot',             # 'protein entries'
        depr  => [ 'up', 'euk', 'arch', 'bact' ],
    },
    ### Taxonomy ###
    taxon => {
        scope => 'other',
        value => 'http://purl.uniprot.org/taxonomy/',
        ident => 'taxonomy',
    },
    ### PubMed ###
    pmid => {
        scope => 'other',
        value => 'http://rdf.ncbi.nlm.nih.gov/pubmed/',
        ident => 'pubmed',
    },

    ### Miscelaneous vocabulary ( most of them are RDF predicates, hence they are not in identifiers.org )
    #rheaV => {
    #    scope => 'other',
    #    value => 'http://rdf.rhea-db.org/',
    #},
    up => {
        scope => 'other',
        value => 'http://purl.uniprot.org/core/',
    },
    sbo => {
        scope => 'other',
        value => 'https://identifiers.org/SBO:',
    },
    go => {
        scope => 'comp',
        value => 'http://purl.obolibrary.org/obo/GO_',
        ident => 'GO',
    },
    cl => {
        scope => 'comp',
        value => 'http://purl.obolibrary.org/obo/CL_',
        ident => 'CL',
    },
);

sub _validate_prefix_data{
    foreach my $prefix ( sort keys %prefix_data ){
        foreach( 'scope', 'value' ){
            die "Missing mandatory key $_ for prefix: $prefix\n" unless exists $prefix_data{$prefix}{$_};
        }
        foreach( sort keys %{$prefix_data{$prefix}} ){
            die "Invalid key for prefix $prefix: $_\n" unless /^(scope|value|ident|depr)$/;
        }
        unless( $prefix_data{$prefix}{scope} =~ /^(chem|reac|comp|pept|other)$/ ){
            die "Invalid scope for prefix $prefix: $prefix_data{$prefix}{scope}\n";
        }
        unless( $prefix_data{$prefix}{value} ){
            die "Empty value for prefix $prefix!\n";
        }
        if( exists $prefix_data{$prefix}{ident} and $prefix_data{$prefix}{ident} !~ /^[\w\.]+$/ ){
            die "Invalid ident syntax for prefix $prefix: $prefix_data{$prefix}{ident}\n";
        }
        if( exists $prefix_data{$prefix}{depr} ){
            foreach( @{$prefix_data{$prefix}{depr}} ){
                die "Invalid ident syntax for prefix $prefix: $_\n" unless /^\w+$/;
            }
        }
    }
}
sub new{
    my( $package ) = @_;
    _validate_prefix_data();
    my $self = bless {}, $package;
    foreach my $prefix ( sort keys %prefix_data ){
        $self->{prefix}{$prefix} = $prefix_data{$prefix}{value};
        if( exists $prefix_data{$prefix}{ident} and $prefix_data{$prefix}{ident} ne $prefix ){
            my $IRI = 'https://identifiers.org/' . $prefix_data{$prefix}{ident} . ':';
            $self->{prefix2}{$prefix_data{$prefix}{ident}} = $IRI;
            $self->{same_as}{$prefix} = $prefix_data{$prefix}{ident};
        }
        if( exists $prefix_data{$prefix}{ident} ){
            $self->{'fromSBML'}{$prefix_data{$prefix}{'scope'}}{$prefix_data{$prefix}{'ident'}} = $prefix;
            $self->{'toSBML'}{$prefix_data{$prefix}{'scope'}}{$prefix}                          = $prefix_data{$prefix}{'ident'};
        }
        if( exists $prefix_data{$prefix}{depr} ){
            foreach( @{$prefix_data{$prefix}{depr}} ){
                push @{$self->{depr}{$prefix_data{$prefix}{scope}}{$_}}, $prefix;
                $self->{'fromSBML'}{$prefix_data{$prefix}{'scope'}}{$_} = $prefix;
                $self->{'toSBML'}{ $prefix_data{$prefix}{'scope'} }{$_} = $prefix_data{$prefix}{'ident'};
            }
        }
    }
    return $self;
}
sub get_prefix_data{
    return \%prefix_data;
}

sub get_turtle_prefixes{
    my $self = shift;
    my @line = (
        '# This list of prefixes was automatically generated with a perl script.',
        '# !!! Manual modifications will be lost !!!',
        '# This is how to regenerate it:',
        '# clone https://github.com/MetaNetX/MNXtools.git MNXtools',
        '# cd MNXtools/perl',
        q|# perl -e 'use lib "."; use Prefix; print Prefix->new()->get_turtle_prefixes()' > prefixes.ttl|,
        '',
    );
    my %seen;
    foreach my $dbkey ( sort keys %{$self->{prefix}} ){
        $seen{$dbkey} = 1;
        push @line, '@prefix ' . $dbkey . ': <' . $self->{prefix}{$dbkey} . "> .";
        if( my $dbkey2 = $self->{same_as}{$dbkey} ){
            next  if ( $seen{$dbkey2} );
            $seen{$dbkey2} = 1;
            push @line, '@prefix ' . $dbkey2 . ': <' . $self->{prefix2}{$dbkey2} . "> . # proxy for $dbkey: via identifiers.org";
        }
    }
    return join( "\n", @line ) . "\n";
}

sub get_prefix_ontology{
    my $self = shift;
    my @line = (
        '# This output is automatically generated with a perl script.',
        '# !!! Manual modifications will be lost !!!',
        '# This is how to regenerate it:',
        '# clone https://github.com/MetaNetX/MNXtools.git MNXtools',
        '# cd MNXtools/perl',
        q|# perl -e 'use lib "."; use Prefix; print Prefix->new()->get_prefix_ontology()' > prefixes.ttl|,
        '',
        'prefix   sh: <http://www.w3.org/ns/shacl#>',
        'prefix  xsd: <http://www.w3.org/2001/XMLSchema#>',
        'prefix  owl: <http://www.w3.org/2002/07/owl#>',
        'prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>',
        '',
        '_:sparql_examples_prefixes',
        '   a owl:Ontology ;',
        '   rdfs:comment """
This is a collection of SPARQL prefixes used by MetaNetX.
Nota Bene: There are often two prefixes for the same entity. One corresponds to
the IRI used by the RDF community at SIB, and often starts with "http://purl.". It
is the recommended prefix. The other prefix corresponds to the "MIRIAM" prefixes
which were adopted by the Systems Biology community (https://sbml.org/documents/elaborations/miriam_annotation_syntax/),
and typically starts with "https://identifiers.org/". Very unfortunately, identifiers.org
has promoted the usage of the short form of IRIs in SBML annotations, and is maintaining
a list of "official" MIRIAM prefixes at https://registry.identifiers.org. To ensure
interoperability with the Systems Biology community and avoid namespace clashes, MetaNetX
has to respects the MIRIAM prefixe nomenclature, and has no other choice than to define
ad hoc prefixes for those not covered by MIRIAM.

""" ;',
        'owl:imports sh: .',
       ''
    );
    foreach my $dbkey ( sort keys %{$self->{prefix}} ){
        my $scope = $prefix_data{$dbkey}{scope} || 'other';
        push @line,
            "_:sparql_examples_prefixes sh:declare _:prefix_$dbkey .",
            "_:prefix_$dbkey ",
            "   rdfs:comment 'A primary prefix for $scope entries';",
            "   sh:prefix '$dbkey' ;",
            "   sh:namespace '$self->{prefix}{$dbkey}'^^xsd:anyURI .",
            '';
        if( my $dbkey2 = $self->{same_as}{$dbkey} ){
            push @line,
                "_:sparql_examples_prefixes sh:declare _:prefix_$dbkey2 .",
                "_:prefix_$dbkey2 ",
                "   rdfs:comment 'A secondary prefix for $scope entries, i.e. this is a proxy for $dbkey via identifiers.org';",
                "   rdfs:seeAlso _:prefix_$dbkey ;",
                "   sh:prefix '$dbkey2' ;",
                "   sh:namespace '$self->{prefix2}{$dbkey2}'^^xsd:anyURI .",
                '';
        }
    }
    return join( "\n", @line ) . "\n";
}

sub get_identifier{
    my( $self, $dbkey, $id ) =@_;
    if( exists $self->{prefix}{$dbkey} ){
        if( $id =~ /\W/ ){
            return '<' . $self->{prefix}{$dbkey} . $id .'>';
        }
        else{
            return $dbkey . ':' . $id;
        }
    }
    elsif( exists $self->{prefix2}{$dbkey} ){
        if( $id =~ /\W/ ){
            return '<' . $self->{prefix2}{$dbkey} . $id .'>';
        }
        else{
             return $dbkey . ':' . $id;
        }
    }
    else{
        confess "Unknown chem prefix: $dbkey:$id\n";
    }
}
sub get_chem_identifier{
    my( $self, $name ) = @_;
    if( $name =~ /:/ ){
        my( $dbkey, $id ) = $name =~ /^([\w\.]+):(.+)/;
#        if( $dbkey eq 'kegg' ){
#            if( $id =~ /^M_(\w)/ ){
#                $dbkey .= $1;
#            }
#            else{
#                $id =~ /^(\w)/;
#                $dbkey .= $1;
#            }
#        }
#        elsif( $dbkey =~ /^(bigg|seed|metacyc|sabiork)$/ ){
#            $dbkey .= 'M';
#        }
        return $self->get_identifier( $dbkey, $id );
    }
    elsif( $name =~ /^MNXM\d+$/ or $name =~ /^BIOMASS|PROTON|PMF|WATER|HYDROXYDE$/ ){
        return 'chem:' . $name;
    }
    else{
        return '_:' . $name;
    }
}
sub get_comp_identifier{
    my( $self, $name ) = @_;
    if( $name =~ /:/ ){
        my( $dbkey, $id ) = $name =~ /^([\w+\.]+):(.+)/;
#        if( $dbkey =~ /^(kegg|bigg|seed)$/ ){
#            $dbkey .= 'C';
#        }
#        elsif( $dbkey eq 'metacyc' ){
#            $dbkey = 'cco';
#        }
#        elsif( $dbkey eq 'rhea' ){
#            $dbkey = 'rheaC';
#        }
#        elsif( $dbkey eq 'sabiork' ){
#            $dbkey = 'sabiorkC';
#        }
        return $self->get_identifier( $dbkey, $id );
    }
    elsif( $name =~ /^MNX[CD]\d+$/ or $name =~ /^BOUNDARY$/ ){
        return 'comp:' . $name;
    }
    elsif( $name eq 'MNXDX' ){ # FIXME: this is forbidden !!!
        return 'comp:' . $name;
    }
    else{
        return '_:' . $name;
    }
}
sub get_reac_identifier{
    my( $self, $name ) = @_;
    if( $name =~ /:/ ){
        my( $dbkey, $id ) = $name =~ /^([\w\.]+):(.+)/;
#        if( $dbkey =~ /^(kegg|bigg|seed|metacyc|sabiork)$/ ){
#            $dbkey .= 'R';
#        }
        return $self->get_identifier( $dbkey, $id );
    }
    elsif( $name =~ /^MNXR\d+$/ or $name =~ /^EMPTY$/ ){
        return 'reac:' . $name;
    }
    elsif( $name eq 'BIOMASS_EXT' ){ # FIXME: this is forbidden !!!
        return 'reac:' . $name;
    }
    else{
        confess "Cannot guess prefix: $name";
    }
}
sub get_pept_identifier{
    my( $self, $name ) = @_;
    if( $name =~ /:/ ){
        my( $dbkey, $id ) = $name =~ /^([\w\.]+):(.+)/;
        if( $dbkey eq 'uniprot' ){
            $dbkey = 'uniprotkb';
        }
        if( exists $self->{prefix}{$dbkey} ){
            if( $id =~ /\W/ ){
                return '<' . $self->{prefix}{$dbkey} . $id .'>';
            }
            else{
                return $dbkey . ':' . $id;
            }
        }
        else{
             confess "Unknown pept prefix: $dbkey\n";
        }
    }
    elsif( $name =~ /^SPONTANEOUS$/ ){
        return 'mnx:' . $name;
    }
    else{
        confess "Cannot guess prefix: $name";
    }
}
sub prefix_exists{
    my( $self, $prefix ) = @_;
    return exists $self->{prefix}{$prefix} || exists $self->{prefix2}{$prefix};
}
sub get_prefix{
    return shift->{prefix};
}
sub get_prefix2{
    return shift->{prefix2};
}
sub get_same_as{
    return shift->{same_as};
}
sub has_secondary_prefix{
    my( $self, $id ) = @_;
    if( $id =~ /^([\w\.]+):/ ){
        return exists $self->{prefix2}{$1};
    }
    else{
        return 0;
    }
}
sub get_prefix_from_depr{
    my( $self, $scope, $depr ) = @_;
    if( exists $self->{depr}{$scope}{$depr} ){
        return @{$self->{depr}{$scope}{$depr}};
    }
    else{
        return ();
    }
}

1;

