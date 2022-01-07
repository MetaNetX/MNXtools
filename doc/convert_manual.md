# Status codes (conversion to MNXref)

In the yaml output, the mapping reports are organized according to the `chem:`/`comp:`/`reac:` identifiers of the source model to permit the utilisation of the report without considering further the mapped model.

The mapping always provides the source `ID_src:` and the destination `ID_dst:` identifiers and a `status` list made of the codes explained below.
Some of these codes come with additional attributes.

The names of the different entities are supplied as comments in the yaml log to facilitate it for humans (as most of us are still part of them), some of these names are propagated from the source model, the others are taken from MNXref, depending on the context.

_Suggested model improvements_ are given to help improve the formulation of genome-scale metabolic networks. However, it should be realized that these suggestions might be totally not adapted for other applications. The users are strongly advised not to follow them blindly, and to report the problematic cases here or there. You can report MetaNetX related issues via our dedicated github https://github.com/MetaNetX/MNXref/issues

The following code are produced:

**CHEM_XREF_OK**

* This code is only produced when cross-refs are exploited.
* There is a single cross-ref or there are multiple cross-refs that correspond to a single MNXref identifier.

**CHEM_XREF_CONFLICT**

* This code is only produced when cross-refs are exploited.
* There are multiple cross-refs that correspond to different unrelated MNXref identifiers.
* _Suggested model improvement_: select which cross-ref is the best one, keep it and move the others elsewhere (in comments, for example) to avoid recreating this conflict later.

**CHEM_XREF_AMBIGUOUS**

* This code is only produced when cross-refs are exploited.
* There are multiple cross-refs that correspond to different MNXref identifiers, but these are related by stereochemical parent/child relationships.
* The parent MNXref identifier is selected automatically to preserve chemistry, at the expense of precision.
* _Suggested model improvement_: select which cross-ref is the best one, possibly among those with the more detailed stereochemistry. It however cannot be excluded that the parent cross-ref corresponds to a different metabolite, which precise stereochemistry could be specified from a different xref.

**CHEM_MAP_OK**

* A unique mapping returns a single MNXref identifier.

**CHEM_MAP_MULTIPLE**

* The mapping returns more than one MNXref identifier and one of them is **arbitrarily** selected in the mapped model.
* This situation may arise because not enough information is supplied in the source model or because of a deprecated MNXref identifier which was remapped onto more than one identifier in the latest MNXref releases.
* _Suggested model improvement_: update the identifier in the source model to remove the ambiguity.

**CHEM_MAP_DEPRECATED**

* The source identifier is a deprecated MNXref identifier.
* It might cause a CHEM_MAP_MULTIPLE code, or not.
* _Suggested model improvement_: update the identifier in the source model.

**CHEM_MAP_UNKNOWN**

* The metabolite cannot be mapped to an MNXref identifier.
* The prefix `UNK:` was added to the original metabolite identifier in the ad hoc MNX format.
* _Suggested model improvement_: update the identifier in the source model if (and only if) the chemical exists somewhere in MNXref namespace.

**CHEM_MNET_MERGE**

* Two chemicals from the original model were merged into the same MNXref identifier.
* _Suggested model improvement_: merge the two original identifiers or correct one of them to realize the distinction or report a bug in MNXref.

**CHEM_MNET_ISOMERIC**

* Two or more different MNXref identifiers were found related by isomeric parent/child relationships.
* The code is used for the parent and list the children.
* _Suggested model improvement_: precise the stereochemistry of the parent identifier, as it might be the same metabolite as the child or a different one.

**COMP_MAP_OK**

**COMP_MAP_DEPRECATED**

**COMP_MAP_UNKNOWN**

**COMP_MAP_MULTIPLE**

* The mapping returns more than one MNXref identifier and one of them is **arbitrarily** selected in the mapped model.
* This situation may arise because not enough information is supplied in the source model or because of a deprecated MNXref identifier which was remapped onto more than one identifier in the latest MNXref releases.
* _Suggested model improvement_: update the identifier in the source model to remove the ambiguity.

**COMP_GENERIC**

* The original compartments have been replaced by generic compartments MNXD1 and MNXD2 as in the MNXref distribution.
* The connectivity of the network is likely to have been lost, if it contained more than one compartment.

**REAC_MAP_OK**

* The original equation was converted into an equation which reactants are all mapped to MNXref.

**REAC_MAP_MNXREF**

* The original equation was converted into an equation which reactants are all mapped to MNXref.
* In addition, the mapped equation corresponds to a known equation into the MNXref repository.

**REAC_MAP_UNKNOWN**

* The original equation was converted into an equation which contains one or several reactants not found in MNXref.

**REAC_EMPTY_MNXREF**

* The original equation was converted into an EMPTY equation, _i.e._ ` = `.
* It is an acid-base reaction and/or a tautomerization.
* The identifier of this reaction is registered in MNXref as an EMPTY reaction.
* _Suggested model improvement_: empty reactions should be removed, after merging the implied metabolites (see CHEM_MNET_MERGE code).

**REAC_EMPTY_UNKNOWN**

* The original equation was converted into an EMPTY equation, _i.e._ ` = `.
* The identifier of this reaction is NOT registered in MNXref as an EMPTY reaction.
* _Suggested model improvement_:
	* Best case scenario: it is an acid-base reaction and/or a tautomerization, the metabolites should be merged and the reaction removed.
	* Worst case scenario: all reactants are lost because of mapping errors, fix the chemistry of the original metabolites, to preserve the reaction in the model.

**REAC_MAP_LOSS**

* A reactant has been lost from an equation, because it was present on both sides of the equation with the same stoichiometric coefficient. As a consequence, the mapped reactions likely differ from the source reaction, which is possibly the most severe problem that can be encountered.
* This code is not reported for empty reactions.
* _Suggested model improvement_: Work on the metabolites to enforce the distinction between the one in the left and right terms of the equation. It cannot be excluded that the problem arises from a mistake in the MNXref reconciliation, that should be reported here (thanks in advance).

**REAC_MAP_PROTON_SALAD**

* Failure to distinguish proton from PMF, for different reasons, for example if the equation implies three compartments.

**REAC_MNET_MERGE**

* Two or more different original reactions were mapped onto a single one.
* _Suggested model improvement_: first, merge the implied metabolites into a single one; secondly, merge the reactions into a single one.
* _Nota Bene_: On the contrary to SBML, MNXtools are agnostic here with respect to reaction directions. Directionality constraints are placed on top of (undirected) equations, together with enzyme descriptions. This might explain part of the observed merges.

