package Chemicals;
#Filename Chemicals.pm

use strict;
use warnings;
use diagnostics;

use List::Util qw( uniq );
use FindBin qw( $RealBin );
use lib "$RealBin/../perl/SBML";
use Constants;
use Formaters;
use Notes;
use Annotation;
use Compartments;


our $chemicals;
#my $chem_list;
our $original_boundary_compartment = $Constants::boundary_comp_id;


sub parse_Chemicals_Info{
    my ($SBML_model, $SBML_object) = @_;

    for my $species ( $SBML_model->getListOfSpecies() ){
        my $Sid = $species->getId();
        my $clean_Sid = Formaters::format_Chemical($Sid, $species->getCompartment(), '');
#        $chem_list->{$Sid} = $clean_Sid; #TODO Will be used to check reaction chemicals if exist!

        #Fix compartment in case of undeclared boundary compartment
        if ( $species->getBoundaryCondition() == 1 ){
            $original_boundary_compartment = $species->getCompartment();
            $species->setCompartment($Constants::boundary_comp_id);
            Compartments::create_compartment($SBML_model, $Constants::boundary_comp_id, $Constants::boundary_comp_name);
        }

        if ( !exists $chemicals->{ $clean_Sid } ){
            $chemicals->{ $clean_Sid }->{'Sid'}               = $clean_Sid;
            $chemicals->{ $clean_Sid }->{'name'}              = Formaters::trim( $species->getName() // '' );
            $chemicals->{ $clean_Sid }->{'xrefs'}             = [];
            $chemicals->{ $clean_Sid }->{'source'}            = [$Sid];
            $chemicals->{ $clean_Sid }->{'sboTerm'}           = $species->getSBOTermID();
            $chemicals->{ $clean_Sid }->{'boundaryCondition'} = $species->getBoundaryCondition();
            $chemicals->{ $clean_Sid }->{'compartment'}       = [$species->getCompartment()];


            my $fbcExtension = $species->getPlugin('fbc');
            # Get formula
            #NOTE For people storing formula with name description
            if ( $chemicals->{ $clean_Sid }->{'name'} =~ /^(.+?)_([A-Z][a-z]?\d*([A-Z][a-z]?\d*)*)$/ ){
                $chemicals->{ $clean_Sid }->{'name'}    = $1;
                $chemicals->{ $clean_Sid }->{'formula'} = $2;
            }
            if ( $SBML_object->getLevel() >= 3 && $fbcExtension ){
                if ( $fbcExtension->isSetChemicalFormula() ){
                    $chemicals->{ $clean_Sid }->{'formula'} = $fbcExtension->getChemicalFormula();
                }
            }

            # Get charge
            #SBML levels 2 and 1
            if ( $SBML_object->getLevel() <= 2 && $species->isSetCharge() ){
                $chemicals->{ $clean_Sid }->{'charge'} = $species->getCharge();
            }
            else {
                if ( $fbcExtension && $fbcExtension->isSetCharge() ){
                    $chemicals->{ $clean_Sid }->{'charge'} = $fbcExtension->getCharge();
                }
            }


            # Parse notes
            if ( $species->isSetNotes() ){
                my $notes;
                $notes = Notes::parse_notes($species, $notes);
                for my $note (keys %$notes){
                    #NOTE skip (redundant and maybe not sync): COMPOUND_ID, SOURCE
                    next  if ( $note eq 'COMPOUND_ID' || $note eq 'SOURCE' );
                    next  if ( $note eq 'NOTES' );
                    #NOTE skip (what is it?): STRINGCODE
                    if ( $note eq 'CHARGE' ){
                        $chemicals->{ $clean_Sid }->{'charge'} = $notes->{$note}->[0]  if ( !exists $chemicals->{ $clean_Sid }->{'charge'} );
                        $chemicals->{ $clean_Sid }->{'charge'} = ''                    if ( $chemicals->{ $clean_Sid }->{'charge'} eq 'NaN' );
                    }
                    elsif ( $note eq 'FORMULA' || $note eq 'formula' ){
                        $chemicals->{ $clean_Sid }->{'formula'} = $notes->{$note}->[0]  if ( !exists $chemicals->{ $clean_Sid }->{'formula'} );
                    }
                    elsif ( $note eq 'MASS' ){
                        $chemicals->{ $clean_Sid }->{'mass'} = $notes->{$note}->[0];
                    }
                    elsif ( $note eq 'REFERENCE' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { Formaters::guess_prefix('chem', $_) } split(/;/, $notes->{$note}->[0]);
                    }
                    elsif ( $note eq 'XREFS' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { Formaters::guess_prefix('chem', $_) } split(/;/, $notes->{$note}->[0]);
                    }
                    elsif ( $note eq 'CHEBI' || $note eq 'ChEBIID' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'CHEBI:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'HMDB' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'hmdb:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'LIPID-MAPS' || $note eq 'LIPIDMAPS' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'lipidmapsM:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'BIGG' || $note eq 'RECON2' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'biggM:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'BIOCYC' || $note eq 'METACYC' || $note eq 'ECOCYC' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'metacycM:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'BIOPATH' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'biopath:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'EAWAG-BBD-CPD' || $note eq 'UMBBDCPD' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'umbbd:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'KEGG' || $note eq 'KEGG ID' || $note eq 'KEGG-GLYCAN' || $note eq 'LIGANDCPD' || $note eq 'KEGGID' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { s/^(.)/kegg$1:$1/; $_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'MNXREF' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'mnx:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'REACTOME' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'reactomeM:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'SEED' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'seedM:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'UPA' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'upa:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'CAS' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'cas:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'CHEMSPIDER' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'chemspider:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'DRUGBANK' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'drugbank:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'KNAPSACK' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'knapsack:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'METABOLIGHTS' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'metabolights:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'NCI' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'nci:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'PUBCHEM' || $note eq 'PubChemID' ){
                        push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map { 'pubchem:'.$_ } @{$notes->{$note}};
                    }
                    elsif ( $note eq 'INCHI' || $note eq 'InChIString' ){
                        if ( $notes->{$note}->[0] !~ /^InChI=/ ){
                            $notes->{$note}->[0] = 'InChI='.$notes->{$note}->[0];
                        }
                        $chemicals->{ $clean_Sid }->{'inchi'}    = $notes->{$note}->[0];
                    }
                    elsif ( $note eq 'INCHIKEY' ){
                        $chemicals->{ $clean_Sid }->{'inchikey'} = $notes->{$note}->[0];
                    }
                    elsif ( $note eq 'SMILES' ){
                        $chemicals->{ $clean_Sid }->{'smiles'}   = $notes->{$note}->[0];
                    }
                    else {
                        $Constants::logger->debug("chemMISS: $note\t\t$notes->{$note}->[0]");
                    }
                }
                #Cleaning
                @{ $chemicals->{ $clean_Sid }->{'xrefs'} } = uniq sort grep { !/:N\/A$/ } @{ $chemicals->{ $clean_Sid }->{'xrefs'} };
                $chemicals->{ $clean_Sid }->{'formula'} = ''  if ( exists $chemicals->{ $clean_Sid }->{'formula'} && $chemicals->{ $clean_Sid }->{'formula'} eq 'N/A' );
            }

            # Parse annotations
            #NOTE annotations overload notes!
            if ( $species->isSetAnnotation() ){
                my $xrefs;
                $xrefs = Annotation::parse_annotation($species, $xrefs, 'chem');
                push @{ $chemicals->{ $clean_Sid }->{'xrefs'} }, map  { Annotation::clean_chem_annotation($_) }
                                                                 grep { !/(bio|meta|eco)cyc:.+?\-RXN/ &&
                                                                        !/kegg:R/ && !/kegg\.reaction/ &&
                                                                        !/metanetx\.reaction:/ &&
                                                                        !/mnx:MNXR/ &&
                                                                        !/ec-code/ &&
                                                                        !/bigg\.reaction/
                                                                      } # filter reactions in chemical annotations!
                                                                 map  { @{ $xrefs->{$_} } }
                                                                 grep { !/inchi/i }
                                                                 grep { /_chem$/ }
                                                                 keys %$xrefs;
                ($chemicals->{ $clean_Sid }->{'inchi'})        = map  { Annotation::clean_chem_annotation($_) }
                                                                 map  { $xrefs->{$_}->[0] }
                                                                 grep { $_ eq 'inchi_chem' }
                                                                 keys %$xrefs;
                ($chemicals->{ $clean_Sid }->{'inchikey'})     = map  { Annotation::clean_chem_annotation($_) }
                                                                 map  { $xrefs->{$_}->[0] }
                                                                 grep { $_ eq 'inchikey_chem' || $_ eq 'inchi_key_chem' }
                                                                 keys %$xrefs;
            }
        }
        else {
            push @{ $chemicals->{ $clean_Sid }->{'source'} },      $Sid;
            push @{ $chemicals->{ $clean_Sid }->{'compartment'} }, $species->getCompartment();
            @{ $chemicals->{ $clean_Sid }->{'xrefs'} } = uniq sort @{ $chemicals->{ $clean_Sid }->{'xrefs'} };
        }
    }

    return;
}


sub add_chemical {
    my ($hash) = @_;

    my ($Sid) = keys %$hash;
    if ( !exists $chemicals->{ $Sid } && 1 == scalar keys %$hash ){
        %$chemicals = (%$chemicals, %$hash);
        @{ $chemicals->{ $Sid }->{'xrefs'} } = uniq sort @{ $chemicals->{ $Sid }->{'xrefs'} };
    }

    return;
}

sub replace_chemical {
    my ($hash, $old_Sid) = @_;

    my ($new_Sid) = keys %$hash;
    if ( !exists $chemicals->{ $new_Sid } && 1 == scalar keys %$hash && exists $chemicals->{ $old_Sid } ){
        $chemicals->{ $new_Sid } = $hash->{ $new_Sid };
        delete $chemicals->{ $old_Sid };
        @{ $chemicals->{ $new_Sid }->{'xrefs'} } = uniq sort @{ $chemicals->{ $new_Sid }->{'xrefs'} };
    }

    return;
}

sub update_chemical {
    my ($hash) = @_;

    my ($Sid) = keys %$hash;
    if ( exists $chemicals->{ $Sid } && 1 == scalar keys %$hash ){
        $chemicals->{ $Sid } = $hash->{ $Sid };
        @{ $chemicals->{ $Sid }->{'xrefs'} } = uniq sort @{ $chemicals->{ $Sid }->{'xrefs'} };
    }

    return;
}


sub create_SBML_chemical {
    my ($MetNet, $mnet_id, $SBML_model, $use_notes, $chem_id) = @_;

    my $chem_id_fixed = $chem_id;
    my $comp_id = '';
    ($chem_id, $comp_id) = $chem_id =~ /^(.+?)__64__(.+?)$/;
    $chem_id_fixed =~ s{^UNK:}{};
    my $chem = $SBML_model->createSpecies();
    $chem->setId($chem_id_fixed);
    #beta-D-fructose C6H12O6 180.06339 0 chebi:28645;CHEBI:28645;bigg.metabolite:fru;biggM:M_fru;biggM:fru;chebi:10373;chebi:22766;chebi:42560;kegg.compound:C02336;...
    my ($chem_name, $chem_formula, $chem_mass, $chem_charge, $chem_xrefs) = $MetNet->get_chem_info($chem_id);
    $chem->setName($chem_name);
    $chem->setCompartment($comp_id);
    $chem->setHasOnlySubstanceUnits(0); #false
    $chem->setConstant(0);              #false
    $chem->setBoundaryCondition( $comp_id eq $Constants::boundary_comp_id ? 1 : 0 );
    $chem->setSBOTerm( $chem_id eq $Constants::biomass_chem_id ? $Constants::biomass_chem_sbo : $Constants::default_chem_sbo );
    my $chem_fbc = $chem->getPlugin('fbc');
    $chem_fbc->setCharge($chem_charge)            if ( $chem_charge  ne '' );
    $chem_fbc->setChemicalFormula($chem_formula)  if ( $chem_formula ne '' );
    my @chem_xrefs = split(';', $chem_xrefs);
    #notes
    if ( $use_notes ){
        my $notes = '';
        $notes .= '<html:p>COMPOUND_ID: '.$chem_id.                                    '</html:p>';
        $notes .= '<html:p>SOURCE: '.     $MetNet->get_chem_source($mnet_id, $chem_id).'</html:p>'  if ( $MetNet->get_chem_source($mnet_id, $chem_id) );
        $notes .= '<html:p>FORMULA: '.    $chem_formula.                               '</html:p>'  if ( $chem_formula ne '' );
        $notes .= '<html:p>MASS: '.       $chem_mass.                                  '</html:p>'  if ( $chem_mass    ne '' );
        $notes .= '<html:p>CHARGE: '.     $chem_charge.                                '</html:p>'  if ( $chem_charge  ne '' );
        $notes .= '<html:p>XREFS: '.      $chem_xrefs.                                 '</html:p>'  if ( $chem_xrefs );
        $notes .= '<html:p>REFERENCE: '.  $chem_xrefs[0].                              '</html:p>'  if ( exists $chem_xrefs[0] );
        $chem->setNotes($notes)  if ( $notes );
    }
    #annotations
    if ( $chem_xrefs ){
        if ( !$chem->isSetMetaId() ){
            $chem->setMetaId( $chem->getId() );
        }
        my $CV = new LibSBML::CVTerm();
        $CV->setQualifierType($LibSBML::BIOLOGICAL_QUALIFIER);
        #is(_a) -> reference
        $CV->setBiologicalQualifierType($LibSBML::BQB_IS);
        $CV->addResource($chem_xrefs[0]);#$Constants::identifiers_biggc.$1);
        $chem->addCVTerm($CV);
        #isRelatedTo -> other xrefs
        if ( exists $chem_xrefs[1] ){
            my $CV2 = new LibSBML::CVTerm();
            $CV2->setQualifierType($LibSBML::BIOLOGICAL_QUALIFIER);
            $CV2->setBiologicalQualifierType($LibSBML::BQB_IS_HOMOLOG_TO); #BQB_IS_RELATED_TO looks better but deprecated now!
#TODO skip duplicated ids: biggM:a, bigg.metabolite:a, biggM:M_a
            for my $xref ( splice(@chem_xrefs, 1) ){
                $CV2->addResource($xref);
            }
            $chem->addCVTerm($CV2);
        }
    }
#TODO use the prefix from the prefix file
#TODO inchikey are missing as xref

    return;
}

1;

