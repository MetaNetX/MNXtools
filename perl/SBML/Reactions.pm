package Reactions;
#Filename Reactions.pm

use strict;
use warnings;
use diagnostics;

use XML::Parser;
use Data::Dumper;
use List::Util qw(min max uniq any);
use FindBin qw( $RealBin );
use lib "$RealBin/../perl/SBML";
use Constants;
use Formaters;
use Notes;
use Annotation;
use Compartments;
use Chemicals;


our $reactions;


sub get_min_max_bounds {
    my ($SBML_model, $SBML_object) = @_;

    #SBML level 2 (and 1?)
    if ( $SBML_object->getLevel() < 3 ){
        my %lfluxes;
        my %ufluxes;
        for my $reac ( $SBML_model->getListOfReactions() ){
            if ( eval {$reac->getKineticLaw()->getParameter('LOWER_BOUND')} ){
                $lfluxes{ $reac->getKineticLaw()->getParameter('LOWER_BOUND')->getValue() } = 1;
            }
            if ( eval {$reac->getKineticLaw()->getParameter('UPPER_BOUND')} ){
                $ufluxes{ $reac->getKineticLaw()->getParameter('UPPER_BOUND')->getValue() } = 1;
            }
        }
        Model::set_LB( (min keys %lfluxes) // Model::get_LB() );
        Model::set_UB( (max keys %ufluxes) // Model::get_UB() );
    }
    #SBML level 3
    else {
        my @param;
        # getNumParameters may contain values declared BUT not used in the model. So not really the min/max of the model!
        if ( $SBML_model->getNumParameters() ){
            @param = map { $SBML_model->getParameter($_-1)->getValue() } 1..$SBML_model->getNumParameters();
        }
        my %fluxes;
        # Get really used parameter/flux values over reactions
        for my $reac ( $SBML_model->getListOfReactions() ){
            my $fbcReaction = $reac->getPlugin('fbc');
            if ( $fbcReaction ){
                if ( $fbcReaction->isSetLowerFluxBound() ){
                    $fluxes{ $SBML_model->getParameter( $fbcReaction->getLowerFluxBound() )->getValue() } = 1;
                }
                if ( $fbcReaction->isSetUpperFluxBound() ){
                    $fluxes{ $SBML_model->getParameter( $fbcReaction->getUpperFluxBound() )->getValue() } = 1;
                }
            }
        }
        for my $param ( sort @param ){
            $Constants::logger->warn("Declared parameter/flux [$param] not used")  if ( !exists $fluxes{$param} );
        }
        Model::set_LB( (min keys %fluxes) // Model::get_LB() );
        Model::set_UB( (max keys %fluxes) // Model::get_UB() );
    }

    return;
}


sub parse_Reactions_Info {
    my ($SBML_model, $SBML_object) = @_;

   my $xml_parser = XML::Parser->new( Style => 'Tree' ); # to latter help parsing GPR

    #Get objective functions if any (SBML 3+)
    my $objectives;
    if ( $SBML_object->getLevel() >= 3 ){
        my $fbcModel = $SBML_model->getPlugin('fbc');
        if ( $fbcModel ){
            for ( my $i=0; $i < $fbcModel->getListOfObjectives()->getNumObjectives(); $i++ ){
                for ( my $j=0; $j < $fbcModel->getObjective($i)->getListOfFluxObjectives()->getNumFluxObjectives(); $j++ ){
                    push @{ $objectives }, $fbcModel->getObjective($i)->getListOfFluxObjectives()->get($j)->getReaction();
                }
            }
        }
    }

    for my $reaction ( $SBML_model->getListOfReactions() ){
        my $Sid = $reaction->getId();
        my $clean_Sid = Formaters::format_Reaction($Sid);

        if ( !exists $reactions->{ $clean_Sid } ){
            $reactions->{ $clean_Sid }->{'Sid'}        = $clean_Sid;
            $reactions->{ $clean_Sid }->{'name'}       = Formaters::trim( $reaction->getName() // '' );
            $reactions->{ $clean_Sid }->{'EC'}         = [];
            $reactions->{ $clean_Sid }->{'xrefs'}      = [];
            $reactions->{ $clean_Sid }->{'pathways'}   = [];
            $reactions->{ $clean_Sid }->{'source'}     = [$Sid];
            $reactions->{ $clean_Sid }->{'fast'}       = $reaction->getFast();
            $reactions->{ $clean_Sid }->{'reversible'} = $reaction->getReversible();
            # Get compartment if any
            if ( $SBML_object->getLevel() == 3 && $reaction->isSetCompartment() ){
                $reactions->{ $clean_Sid }->{'compartment'} = $reaction->getCompartment();
            }

            my $fbcReaction = $reaction->getPlugin('fbc');
            # Get bounds
            #SBML level 3
            if ( $SBML_object->getLevel() >= 3 && $fbcReaction ){
                if ( $fbcReaction->isSetLowerFluxBound() ){
                    $reactions->{ $clean_Sid }->{'LB'} = $SBML_model->getParameter( $fbcReaction->getLowerFluxBound() )->getValue();
                }
                if ( $fbcReaction->isSetUpperFluxBound() ){
                    $reactions->{ $clean_Sid }->{'UB'} = $SBML_model->getParameter( $fbcReaction->getUpperFluxBound() )->getValue();
                }
                if ( grep { $_ eq $reaction->getId() } @{$objectives} ){
                    $reactions->{ $clean_Sid }->{'is_objective'} = 1;
                }
            }
            #SBML level 2 (and 1?)
            elsif ( $SBML_object->getLevel() < 3 ){
                if ( eval {$reaction->getKineticLaw()->getParameter('LOWER_BOUND')} ){
                    $reactions->{ $clean_Sid }->{'LB'} = $reaction->getKineticLaw()->getParameter('LOWER_BOUND')->getValue();
                }
                if ( eval {$reaction->getKineticLaw()->getParameter('UPPER_BOUND')} ){
                    $reactions->{ $clean_Sid }->{'UB'} = $reaction->getKineticLaw()->getParameter('UPPER_BOUND')->getValue();
                }
                if ( eval {$reaction->getKineticLaw()->getParameter('OBJECTIVE_COEFFICIENT')} && $reaction->getKineticLaw()->getParameter('OBJECTIVE_COEFFICIENT')->getValue() > 0 ){
                    $reactions->{ $clean_Sid }->{'is_objective'} = 1;
                }
            }

            # Get sboTerm
            $reactions->{ $clean_Sid }->{'sboTerm'} = $reaction->getSBOTermID();

            # Get equation parts
            my $reactants = parse_Reactants_Info([ $reaction->getListOfReactants() ], $SBML_object);
            my $products  = parse_Products_Info( [ $reaction->getListOfProducts() ],  $SBML_object);
            my $modifiers = parse_Modifiers_Info([ $reaction->getListOfModifiers() ], $SBML_object);
            # If a SINGLE reactant but no product, assume this is an implicit boundary reaction
            ($reactants, $products) = find_implicit_boundary_reaction($SBML_model, $reactants, $products);
            #NOTE do not assume a SINGLE product but no reactant is an implicit boundary reaction
            #($products, $reactants) = find_implicit_boundary_reaction($SBML_model, $products, $reactants);
            $reactants = exists $reactants->[0] ? $reactants : [''];
            $products  = exists $products->[0]  ? $products  : [''];
            $modifiers = exists $modifiers->[0] ? $modifiers : [''];

            # Build equation
            #NOTE Keep the same order as in the original model
            $reactions->{ $clean_Sid }->{'equation'} = join(' + ', @$reactants) . ' = ' . join(' + ', @$products);


            # Parse notes
            if ( $reaction->isSetNotes() ){
                my $notes;
                $notes = Notes::parse_notes($reaction, $notes);
                for my $note (keys %$notes){
                    #NOTE skip (redundant and maybe not sync): REACTION, SOURCE
                    next  if ( $note eq 'REACTION' || $note eq 'SOURCE' );
                    next  if ( $note =~ /^PROP B/ );
                    next  if ( $note eq 'Recon2 formula' || $note eq 'Model Reaction Number' || $note eq 'Description' || $note eq 'MitoCarta' || $note eq 'MitoCarta score' ||
                               $note eq 'Directionality rationale' || $note eq 'HPA protein level in heart tissue' || $note eq 'HPA mRNA level in heart tissue' || $note eq 'Gene name' ||
                               $note eq 'HGNC' || $note eq '(HGNC' || $note eq 'GENE_LIST' || $note eq 'NOTES' || $note eq 'AUTHORS' ||
                               $note eq 'PROTEIN_ASSOCIATION' || $note eq 'TYPE' || $note eq 'STATUS' || $note eq 'HOLE' || $note eq 'SCORE' || $note eq 'GENERIC' );
                    next  if ( $notes->{$note}->[0] eq '' || $notes->{$note}->[0] eq 'S_' );
                    if ( $note eq 'SUBSYSTEM' || $note eq 'PATHWAYS' || $note eq 'KEGG_MAP' ){
                        push @{ $reactions->{ $clean_Sid }->{'pathways'} }, @{$notes->{$note}};
                    }
                    elsif ( $note eq 'PROTEIN_CLASS' || $note eq 'CLASSIFICATION' || $note eq 'EC_NUMBER' ||
                            $note eq 'EC NUMBER' || $note eq 'EC Number' || $note eq 'EC_Number' || $note eq 'ec-code' || $note eq 'sabiork.ec' || $note eq 'ec' ){
                        if ( $notes->{$note}->[0] ne 'No_Assignment' ){
                            push @{ $reactions->{ $clean_Sid }->{'EC'} }, @{$notes->{$note}};
                        }
                    }
                    elsif ( $note eq 'GENE_ASSOCIATION' || $note eq 'GENE ASSOCIATION' || $note eq 'GPR_ASSOCIATION' ){
                        if ( $notes->{$note}->[0] =~ /GAPFILLING/ || $notes->{$note}->[0] eq 'BOF' || $notes->{$note}->[0] eq 'UNIVERSAL' || $notes->{$note}->[0] eq 'AUTOCOMPLETION' ){
                            #Not real genes/peptides
                            next;
                        }
                        if ( $notes->{$note}->[0] ne '()' && $notes->{$note}->[0] ne 'Unknown' && $notes->{$note}->[0] ne 'N/A' && $notes->{$note}->[0] ne 'Non-enzymatic' && $notes->{$note}->[0] ne 'Non-Enzymatic' && $notes->{$note}->[0] ne 'NfC' && $notes->{$note}->[0] ne 'NoAssignment' ){
                            push @{ $reactions->{ $clean_Sid }->{'peptides_notes'} }, @{$notes->{$note}};
                        }
                    }
                    elsif ( $note eq 'REFERENCE' ){
                        push @{ $reactions->{ $clean_Sid }->{'xrefs'} }, map { Formaters::guess_prefix('reac', $_) } split(/;/, $notes->{$note}->[0]);
                    }
                    elsif ( $note eq 'XREFS' ){
                        push @{ $reactions->{ $clean_Sid }->{'xrefs'} }, map { Formaters::guess_prefix('reac', $_) } split(/;/, $notes->{$note}->[0]);
                    }
                    elsif ( $note eq 'KEGG_RID' || $note eq 'KEGG id' ){
                        push @{ $reactions->{ $clean_Sid }->{'xrefs'} }, map { 'keggR:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'BIOCYC' || $note eq 'METACYC' || $note eq 'ECOCYC' ){
                        push @{ $reactions->{ $clean_Sid }->{'xrefs'} }, map { 'metacycR:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'BIGG' || $note eq 'RECON2' || $note eq 'Recon2 id' ){
                        push @{ $reactions->{ $clean_Sid }->{'xrefs'} }, map { 'biggR:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note =~ /Confidence level/i ){
                        $reactions->{ $clean_Sid }->{'Confidence level'} = $notes->{$note}->[0];
                    }
                    else {
                        $Constants::logger->debug("reacMISS: $note\t\t$notes->{$note}->[0]");
                    }
                }
                #Cleaning
                @{ $reactions->{ $clean_Sid }->{'xrefs'} } = grep { !/:?N\/A$/ } @{ $reactions->{ $clean_Sid }->{'xrefs'} };
            }

            # Parse annotations
            #NOTE annotations overload notes!
            if ( $reaction->isSetAnnotation() ){
                my $xrefs;
                $xrefs = Annotation::parse_annotation($reaction, $xrefs, 'reac');
                push @{ $reactions->{ $clean_Sid }->{'xrefs'} }, map { Annotation::clean_reac_annotation($_) }
                                                                 grep { !/^ec-code:/ }
                                                                 map  { @{ $xrefs->{$_} } }
                                                                 grep { /_reac$/ }
                                                                 keys %$xrefs;
                push @{ $reactions->{ $clean_Sid }->{'EC'} },    map  { s{^ec-code:}{}; $_ }
                                                                 grep { /^ec-code:/ }
                                                                 map  { @{ $xrefs->{$_} } }
                                                                 grep { /_reac$/ }
                                                                 keys %$xrefs;
            }


            # Parse gene association
            ## SBML3 (fbc:geneProduct)
            if ( $SBML_object->getLevel() >= 3 && $fbcReaction ){
                if ( $fbcReaction->isSetGeneProductAssociation() ){
                    #NOTE does not look to have stoichiometric coefficients (yet) in geneProductAssociation
                    my $gpa = $fbcReaction->getGeneProductAssociation();
                    # $gpa->getAssociation() is currently broken in SBML/fbc Perl binding. Hence for the fix below
                    my $gpa_xml = $gpa->toSBML();
                    my $parse_tree = $xml_parser->parse( $gpa_xml );
                    my @cplx = parse_cplx_tree( $parse_tree, $gpa_xml );
                    $reactions->{ $clean_Sid }{'peptides'} = join ';', @cplx; # @cplx was already sorted in parse_cplx_tree
                }
            }
            ## then (notes)
            elsif ( exists $reactions->{ $clean_Sid }->{'peptides_notes'} && scalar @{ $reactions->{ $clean_Sid }->{'peptides_notes'} } >= 1 ){
                my $peptides_notes = join(' ', @{ $reactions->{ $clean_Sid }->{'peptides_notes'} });
                #NOTE Replace bad GENE_ASSOCIATION (e.g. SPLC1_S180040 and (cbiO (SPLC1_S040400)) and (cbiO (SPLC1_S170660)) ) that creates a huge memory leak
                # or (glpX-SEBP (SPLC1_S081020))
                $peptides_notes =~ s{\([\w\-\.]+ \((\w+)\)\)}{$1}g;
                $peptides_notes =~ s{\[}{\(}g; # In maranase iRS1597
                $peptides_notes =~ s{\]}{\)}g;
                $peptides_notes =~ s{\)\s+\(}{\)\(}g;
                $peptides_notes =~ s{\)\s+\)}{\)\)}g;
                $peptides_notes =~ s{\(\s+\(}{\(\(}g;
#                $Constants::logger->debug($peptides_notes);
                #Convert to our format
                $peptides_notes =~ s{spontaneous}{$Constants::spontaneous_enz}ig;
                $peptides_notes =~ s{\s+and\s+}{\+}g;
                $peptides_notes =~ s{\s+or\s+}{;}g;
                #NOTE assume they are gene/protein identifiers, without any spaces
                $peptides_notes =~ s{\(([^\(\)\s]+)\)}{$1}g;
                $peptides_notes =~ s{[\(\)]}{}g;
                $reactions->{ $clean_Sid }->{'peptides'} = join('', map { s{\s*(\S+)\s*}{$1}g; $_ }
                                                                    grep { /./ }
                                                                    split(/([\+;])/, $peptides_notes)
                                                               );
#                $Constants::logger->debug($reactions->{ $clean_Sid }->{'peptides'});
                #FIXME Issues with associations like this: (ENSG00000124406 or ENSG00000143515 or ENSG00000081923 or ENSG00000104043 or ENSG00000054793 or ENSG00000166377 or ENSG00000206190 or ENSG00000145246 or ENSG00000068650 or ENSG00000058063 or ENSG00000101974) and (ENSG00000112697 or ENSG00000182107)
                # or (SPLC1_S082520 and SPLC1_S082530 and SPLC1_S082540 and SPLC1_S230690 and SPLC1_S230690 and SPLC1_S170060 and SPLC1_S630130) and (SPLC1_S033260 and SPLC1_S033260 and SPLC1_S033270 and SPLC1_S033240 and SPLC1_S033230 and SPLC1_S082540 and SPLC1_S082540 and SPLC1_S033220 and SPLC1_S033210 and SPLC1_S411020 and SPLC1_S033250)
                #FIXME deal with stochiometric coefficient on enzymes!
            }
            ## then modifiers
            #NOTE discard because modifiers can also be co-factor compounds, not enzymes!
#            elsif ( scalar @$modifiers >= 1 ){
#                $reactions->{ $clean_Sid }->{'peptides'} = join(';', @$modifiers);
#            }

#TODO  get transport and boundary through their sbo terms
        }
        else {
            push @{ $reactions->{ $clean_Sid }->{'source'} }, $Sid;
        }
    }

    return;
}

sub check_elem{
    my( $hash_ref, $xml ) = @_;
    foreach( sort keys %$hash_ref ){
        die "Unexpected '$_': $xml\n"  unless /^(metaid)$/;
    }
}

sub parse_cplx_tree{ # return an ARRAY of complexes (OR) writen using + (AND)
    my( $tree, $xml ) = @_;
    my $action = $tree->[0];
    if( $action eq 'fbc:geneProductAssociation' ){ # this is a placeholder
        check_elem( $tree->[1][0], $xml );
        return parse_cplx_tree( [ $tree->[1][3], $tree->[1][4]], $xml ); # reccursion starts here
    }
    elsif( $action eq 'fbc:or' ){
        my @buf = @{$tree->[1]};
        check_elem( shift @buf, $xml );
        my @token = ();
        while( @buf ){
            my @tree = splice @buf, 0, 2;
            next  if $tree[0] eq '0'; # see XML::parse for details
            push @token, parse_cplx_tree( \@tree, $xml );
        }
        return sort @token;
    }
    elsif( $action eq 'fbc:and' ){
        my @buf = @{$tree->[1]};
        check_elem( shift @buf, $xml );
        my @token = ();
        while( @buf ){
            my @tree = splice @buf, 0, 2;
            next  if $tree[0] eq '0';
            my @res = parse_cplx_tree( \@tree, $xml ); # Can be longer than one, if there is an inner OR !
            if( @token == 0 ){
                @token = @res;
            }
            else{
                my @tmp = ();
                foreach my $outer ( @token ){
                    foreach( @res ){
                        push @tmp, join '+', sort( $_, split /\+/, $outer );
                    }
                }
                @token = @tmp;
            }
        }
        return sort @token;
    }
    elsif( $action eq 'fbc:geneProductRef' ){ # reccursion ends here
        if( exists $tree->[1][0]{'fbc:geneProduct'} ){
            die 'Invalid geneProduct syntax: ' . $tree->[1][0]{'fbc:geneProduct'}  if $tree->[1][0]{'fbc:geneProduct'} =~ /;\+/;
            $tree->[1][0]{'fbc:geneProduct'} =~ s/^G_//;
            return $tree->[1][0]{'fbc:geneProduct'};
        }
        else{
            die 'fbc:geneProduct not found!';
        }
    }
    else{
        die 'Cannot interpret "$action": ' . $xml;
    }
}

sub parse_Reactants_Info {
    my ($element_list, $SBML_object) = @_;

    return parse_Reaction_Elements_Info($element_list, $SBML_object);
}

sub parse_Products_Info {
    my ($element_list, $SBML_object) = @_;

    return parse_Reaction_Elements_Info($element_list,  $SBML_object);
}

sub parse_Modifiers_Info {
    my ($element_list, $SBML_object) = @_;

    return parse_Reaction_Elements_Info($element_list, $SBML_object);
}

sub parse_Reaction_Elements_Info {
    my ($element_list, $SBML_object) = @_;

    my $elements;
    for my $element ( @{ $element_list } ){
        my $species = Formaters::format_Chemical($element->getSpecies(), $SBML_object->getModel->getSpecies( $element->getSpecies() )->getCompartment(), $Chemicals::original_boundary_compartment);
        # Find BIOMASS compound if any
        if ( grep { uc $species eq uc $_ } @Constants::biomass_alias ){
            if ( $species ne $Constants::biomass_chem_id ){
                #Need to replace the chemicals
                my $biomass_chem_hash->{ $Constants::biomass_chem_id } = $Chemicals::chemicals->{ $species };
                $biomass_chem_hash->{ $Constants::biomass_chem_id }->{'Sid'} = $Constants::biomass_chem_id;
                Chemicals::replace_chemical($biomass_chem_hash, $species);
            }
            $species = $Constants::biomass_chem_id;
        }
        if ( ! $SBML_object->getModel->getSpecies( $element->getSpecies() ) ){
            $Constants::logger->fatal($element->getSpecies(). ' is not declared in the model');
            exit 1;
        }
        push @{ $elements }, ($element->getStoichiometry() || 1).' '.$species.'@'.$SBML_object->getModel->getSpecies( $element->getSpecies() )->getCompartment();
#        print $elements->[-1], "\n";
    }

    return $elements;
}


sub find_implicit_boundary_reaction {
    my ($SBML_model, $side1, $side2) = @_;

    #e.g. 1 reactant <==> no product
    if ( exists $side1->[0] && scalar @$side1 == 1 && ! exists $side2->[0] ){
        my $new_chem = $side1->[0];
        $new_chem =~ s{\@.+$}{\@$Constants::boundary_comp_id};
        $side2->[0] = $new_chem;
        Compartments::create_compartment($SBML_model, $Constants::boundary_comp_id, $Constants::boundary_comp_name);
    }
#TODO  use $Constants::boundary_comp_sbo

    return ($side1, $side2);
}


sub find_growth_reaction {
    my ($reactions) = @_;

    REAC:
    for my $reac ( keys %$reactions ){
        my $reaction = $reactions->{$reac};
        #Not real growth equation
        next REAC  if ( grep { $_ eq $reaction->{'Sid'} } @Constants::are_not_biomass_reac );
        for my $reac_src ( @{ $reaction->{'source'} } ){
            next REAC  if ( grep { $_ eq $reac_src } @Constants::are_not_biomass_reac );
        }

        #Exclude biomass export
        next REAC  if ( $reaction->{'equation'} =~ /\@$Constants::boundary_comp_id( |$)/ );


        my ($reactants, $products) = split(/ += +/, $reaction->{'equation'}, -1);
        my @reactants = map { m/ (.+)\@/; $1 } split(/\s+\+\s+/, $reactants);
        my @products  = map { m/ (.+)\@/; $1 } split(/\s+\+\s+/, $products);

        #Find by SBO term (and not biomass transport!)
        if ( $reaction->{'sboTerm'} eq $Constants::biomass_reac_sbo && scalar @reactants > 0 ){
            $Constants::logger->debug('Growth equation by SBO: ['. $reaction->{'Sid'}. '] ['. ( $reaction->{'name'} // '' ). ']');
        }
        #Find by biomass chemical
        elsif ( grep { $_ eq $Constants::biomass_chem_id } @products || grep { $_ eq $Constants::biomass_chem_id } @reactants ){
            $Constants::logger->debug('Growth equation by biomass chemical: ['. $reaction->{'Sid'}. '] ['. ( $reaction->{'name'} // '' ). ']');
        }
        #Find by name (growth or biomass in reac names/ids)
        elsif ( ($reaction->{'Sid'} =~ /.*(biomass|growth).*/i || $reaction->{'name'} =~ /.*(biomass|growth).*/i) &&
                (scalar @reactants > 5 || scalar @products > 5) ){
            $Constants::logger->debug('Growth equation by reaction name/id: ['. $reaction->{'Sid'}. '] ['. ( $reaction->{'name'} // '' ). ']');
        }
        #NOTE find through objective functions as lower priority because models could optimize something else than growth
        #Find with objective functions (SBML 3+ with fbc / SBML 2 (and 1?))
        elsif ( $reaction->{'is_objective'} ){
            $Constants::logger->debug('Growth equation by objective function: ['. $reaction->{'Sid'}. '] ['. ( $reaction->{'name'} // '' ). ']');
        }
#        #Find by reactant number
#        elsif ( scalar @reactants > 15 || scalar @products > 15 ){
#            #NOTE return too many false positives such as seed:rxn05296 "Protein synthesis" or seed:rxn11921 "Generic lipid content reaction"
#        }
        # Not a growth/biomass reaction
        else {
            next REAC;
        }


        # Is biomass compound already in the growth equation?
        my $has_biomass_compound = 0;
        my $local_biomass        = '';
        if ( grep { $_ eq $Constants::biomass_chem_id } @products ){
            ($local_biomass) = grep { / $Constants::biomass_chem_id\@/ } split(/\s+\+\s+/, $products);
            $has_biomass_compound = 1;
        }
        elsif ( grep { $_ eq $Constants::biomass_chem_id } @reactants ){
            ($local_biomass) = grep { / $Constants::biomass_chem_id\@/ } split(/\s+\+\s+/, $reactants);
            $has_biomass_compound = 1;
        }
        #If not, add BIOMASS compound (in the right compartment, by default in products)
        else {
            #Find the most frequent compartment in products
            my %compartments;
            map { m/\@(.+)$/; $compartments{$1}++ } split(/\s+\+\s+/, $products);
            my ($best_compartment) = sort { $compartments{$b} <=> $compartments{$a} || $a cmp $b } keys %compartments;
            $local_biomass = "1 $Constants::biomass_chem_id\@$best_compartment";
            #Add in chemicals
            my $biomass_chem;
            $biomass_chem->{ $Constants::biomass_chem_id }->{'Sid'}  = $Constants::biomass_chem_id;
            $biomass_chem->{ $Constants::biomass_chem_id }->{'name'} = 'Biomass';
            push @{ $biomass_chem->{ $Constants::biomass_chem_id }->{'source'} }, "$Constants::biomass_chem_id\@$best_compartment";
            $biomass_chem->{ $Constants::biomass_chem_id }->{'formula'} = '';
            $biomass_chem->{ $Constants::biomass_chem_id }->{'mass'}    = '';
            $biomass_chem->{ $Constants::biomass_chem_id }->{'charge'}  = '';
            push @{ $biomass_chem->{ $Constants::biomass_chem_id }->{'xrefs'} }, "mnx:$Constants::biomass_chem_id";
            Chemicals::add_chemical($biomass_chem);
            $Reactions::reactions->{ $reac }->{'equation'} .= " + $local_biomass";
            $has_biomass_compound = 1;
        }


        # Add biomass consumption (export to boundary) if not there
        if ( $has_biomass_compound && $local_biomass ){
            # A biomass consumption reaction already exists?
            for my $reac2 ( keys %$reactions ){
                my $export = $reactions->{$reac2};
                my @biomass_reactants = grep { /^$local_biomass$/ } split(/\s+\+\s+/, $reactants);
                my @biomass_products  = grep { /^$local_biomass$/ } split(/\s+\+\s+/, $products);
                #biomass consumption exists
                if (     exists $biomass_reactants[0] && !exists $biomass_products[0] && scalar @reactants == 1 && scalar @products <= 1 && $export->{'Sid'} ne $Constants::biomass_reac_id ){
                    delete $Reactions::reactions->{ Formaters::format_Reaction($export->{'Sid'}) };
                }
                elsif ( !exists $biomass_reactants[0] &&  exists $biomass_products[0] && scalar @reactants >= 1 && scalar @products == 1 && $export->{'Sid'} ne $Constants::biomass_reac_id ){
                    delete $Reactions::reactions->{ Formaters::format_Reaction($export->{'Sid'}) };
                }
            }
        }

        # (Re-)create biomass consumption
        $Reactions::reactions->{ $Constants::biomass_reac_id }->{'Sid'}      = $Constants::biomass_reac_id;
        $Reactions::reactions->{ $Constants::biomass_reac_id }->{'name'}     = $Constants::biomass_reac_id;
        $Reactions::reactions->{ $Constants::biomass_reac_id }->{'source'}   = $Constants::biomass_reac_id;
        my ($local_biomass_chem) = $local_biomass =~ /^(.+?)\@/;
        $Reactions::reactions->{ $Constants::biomass_reac_id }->{'equation'} = "$local_biomass = $local_biomass_chem\@$Constants::boundary_comp_id";
        $Reactions::reactions->{ $Constants::biomass_reac_id }->{'LB'}       = 0;
        $Reactions::reactions->{ $Constants::biomass_reac_id }->{'UB'}       = Model::get_UB();
    }

    return;
}

sub parse_reac_side {
    my ($side) = @_;

    my @side;
    for my $coef_species ( split(/ +\+ +/, $side) ){
        my ($coef, $species) = split(/ /,  $coef_species, 2);
        my ($chem, $comp)    = split(/\@/, $species,      2);
        push @side, [$coef, $chem, $comp];
    }

    return \@side;
}


sub create_SBML_reaction {
    my ($MetNet, $mnet_id, $SBML_model, $use_notes, $reac_id) = @_;

    my $reac_id_fixed = $reac_id;
    $reac_id_fixed =~ s{^UNK:}{};#FIXME should do it better, and more general
    my $reac = $SBML_model->createReaction();
    $reac->setId($reac_id_fixed);
    #ID           equation                               source    mnxr_ID    classifs    pathways    xrefs
    #mnxr02c2b    1 MNXM1@MNXC2 <==> 1 MNXM1@BOUNDARY    EX_h_e    MNXR02                             mnx:MNXR02;MNXR02;bigg.reaction:EX_h_e;...
    my $reac_equation = $MetNet->get_reac_equation($reac_id);
    my $reac_source   = eval {$MetNet->get_reac_source($reac_id)} || ''; #NOTE for our own reactions, no source in the model
    my $reac_mnxr     = $MetNet->get_reac_mnxr($reac_id);
    my ($reac_classifs, $reac_pathways, $reac_xrefs) = $MetNet->get_reac_info($reac_id);
    # Parse equation
    my ($reac_left, $reac_right) = split(/ +<\?> +/, $reac_equation);
    for my $reactant ( @{ Reactions::parse_reac_side($reac_left) } ){
        my $react = $reac->createReactant();
        #TODO need to be added through a species ref !
        #$react->addReactant("$reactant->[1]\@$reactant->[2]", $reactant->[0]);
    }
    my @lefts  = @{ Reactions::parse_reac_side($reac_left) };
    my @rights = @{ Reactions::parse_reac_side($reac_right) };
    #warn "[$reac_left] [$reac_right] [", join(':', map { "$_->[0] $_->[1] $_->[2]" } @lefts), "]\n";
#TODO get reaction direction: get_reac_dir($reac_id);
#TODO are growth reactions there?
#TODO metaid="MAR03905" sboTerm="SBO:0000176" reversible="false" fast="false" fbc:lowerFluxBound="FB2N0" fbc:upperFluxBound="FB3N1000"
#TODO list of reactants/products(/modifiers)

    return;
}

1;

