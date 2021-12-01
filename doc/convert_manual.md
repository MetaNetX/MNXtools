# Status codes (conversion to MNXref)

In the yaml output, the mapping reports are organized according to the `chem:`/`comp:`/`reac:` identifiers of the source model to permit the utilisation of the report without considering further the mapped model.

The mapping always provides the source `ID_src:` and the destination `ID_dst:` identifiers and a `status` list made of the codes explained below.
Some of these codes come with additional attributes.

The names of the different entities are supplied as comments in the yaml lof to facilitate its by human (as most of us are still part of them), some of these names are propagated from the source model, the others are taken from MNXref, dependig from the context.

_Suggested model improvements_ are given to help improve the formulation of genome-scale metabolic networks. However, it should be realized that these suggestions might be totally inadapted for other applications. The users are strongly advised not to follow them blindly, and report the problemeatic cases here or there. 

The following code are produced:

**CHEM_XREF_OK**

* This code is only produced when cross-refs are exploited
* There is a single cross-ref or there are multiple cross-refs that correspond to a same MNXref identifier

**CHEM_XREF_CONFLICT**

* This code is only produced when cross-refs are exploited
* There are multiple cross-refs that correspond to different unrelated MNXref identifiers
* _Suggested model improvement_: select which cross ref is the best one, keep it and move the others elsewhere (in comments, for example) to avoid recreating this conflict later.

**CHEM_XREF_AMBIGUOUS**

* This code is only produced when cross-refs are exploited
* There are multiple cross-refs that correspond to different MNXref identifiers, but these are related by stereochemical parent/child relationships.
* The parent MNXref identifier is selected automatically to preserve chemistry, at the expense of precision.
* _Suggested model improvement_: select which cross ref is the best one, possibly among those with the more detailed stereochemistry. It however cannot be excluded that the parent cross ref corresponds to a different metabolite, which precise stereochemistry could be specified from a different xref.

**CHEM_MAP_OK**

* A unique mapping returns a single MNXref identifier

**CHEM_MAP_MULTIPLE**

* The mapping returns more than one MNXref identifiers and one of them is **arbitrarily** selected in the mapped model. 
* This situation may arise because of not enough information is supplied in the source model or because of a deprecated MNXref identifiers which was remapped onto more than one identifiers in the latest MNXref releases.
* _Suggested model improvement_: update the identifier in the source model to remove the ambiguity.

**CHEM_MAP_DEPRECATED**

* The source identifier is a deprecated MNXref identifier.
* It might cause a CHEM_MAP_MULTIPLE code, or not.
* _Suggested model improvement_: update the identifier in the source model 

**CHEM_MAP_UNKNOWN**

* The metabolite cannot be mapped to an MNXref identifier.
* The prefix `UNK:` was added to the original metabolite identifier in the ad hoc MNX format. 
* _Suggested model improvement_: update the identifier in the source model if (and only if) the chemical exist somewhere in MNXref namespace 

**CHEM_MNET_MERGE**

* Two chemicals from the original model were merged into the same MNXref identifier.
* _Suggested model improvement_: merge the two original identifiers or correct one of them to realize the distinction or report a bug in MNXref

**CHEM_MNET_ISOMERIC**

* Two or more different MNXref identifiers were found related by isomeric parent/child relationships.
* Identification of parent/child is (unfortunately) not specified in the report.
* _Suggested model improvement_: precise the stereochemistry of the parent identifier, as it might be the same metabolite as the child or a different one.

**COMP_MAP_OK**

**COMP_MAP_DEPRECATED**

**COMP_MAP_UNKNOWN**

**COMP_MAP_MULTIPLE**

* The mapping returns more than one MNXref identifiers and one of them is **arbitrarily** selected in the mapped model.
* This situation may arise because of not enough information is supplied in the source model or because of a deprecated MNXref identifiers which was remapped onto more than one identifiers in the latest MNXref releases.
* _Suggested model improvement_: update the identifier in the source model to remove the ambiguity.

**COMP_GENERIC**

* The original compartments have been replaced by generic compartments MNXD1 and MNXD2 as in the MNXref distribution.
* The connectivity of the network is likely have been lost, if it contained more than one compartment.

**REAC_MAP_OK**

* The original equation was converted into an equation which reactants are all mapped to MNXref

**REAC_MAP_MNXREF**

* The original equation was converted into an equation which reactants are all mapped to MNXref
* In addition, The mapped equation correspond to an known equation int the MNXref repositry

**REAC_MAP_UNKNOWN**

* The original equation was converted into an equation which contain one or several reactant not found in MNXref.

**REAC_MAP_EMPTY**

* The original equation was converted into an EMPTY equation, _i.e._ ` = `
* It shlould be the case for any acid-base and/or tautomerization reaction. 
* If present, the accompanying code REAC_MAP_MNXREF indicates that the reaction belongs to the list of known empty reactions in the MNXref reconciliation. Otherwise the accompanying code is REAC_MAP_UNKNOWN 
* _Suggested model improvement_: empty reactions should get removed after merging the implied metabolite identifiers (see CHEM_MNET_MERGE code) 

**REAC_MAP_LOSS**

* A reactants has been lost from an equation, because it was present on both side of the equation with the same stoichiometric coeficient. As a consequence, the mapped reactions likely differ from the source reaction, which is possibly the most severe problem that can be encountered.
* This code is not reported for empty reactions.
* _Suggested model improvement_: Work on the metabolites to enforce the distinguish the one in the left and right term of the equatoin. It cannot be excluded that the problem arise from a mistake in the MNXref reconciliation, that shuld be reported here (thanks in advance)

**REAC_MAP_PROTON_SALAD**

* Failure to distinguish proton from PMF, for different reasons, for example if the equation implies three compartments.  

**REAC_MNET_MERGE**

* Two or more different original reaction were mapped onto a single one.
* _Suggested model improvement_: first, merge the implied metabolites into a single one; Secondly, merge the reactions into a single one.
* _Nota Bene_: On the contrary to SBML, MNXtools are agnostic here with respect to reaction directions. Directionality constraintsare placed on top of (undirected) equation, together with enzyme descriptions. This might explain part of the observed merges.

