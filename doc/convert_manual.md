# Status codes (conversion to MNXref)

In the yaml output, the mapping reports are organized according to the `chem:`/`comp:`/`reac:` identifiers of the source model to permit the utilisation of the report without considering further the mapped model.

The mapping always provides the source `ID_src:` and the destination `ID_dst:` identifiers and a `status` list made of the codes explained below.
Some of these codes come with additional attributes.
The names of the different entities are supplied as comments to facilitate the reading of the yaml log by a human, some of these names are propagated from the source model, the others are taken from MNXref.

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

**CHEM_MAP_WARN**

* The mapping returns a single MNXref identifier, with a warning message.

**CHEM_MAP_UNKNOWN**

* The metabolite cannot be mapped to an MNXref identifier.
* The prefix `UNK:` was added to the original metabolite identifier in the ad hoc MNX format. It should be converted into `UNK_` in SBML.

**CHEM_MERGE**

* Two chemicals from the original model were merged into the same MNXref identifier.
* _Suggested model improvement_: merge the two original identifiers or correct one of them to realize the distinction or report a bug in MNXref

**CHEM_ISOMERIC**

* Two or more different MNXref identifiers are related by isomeric parent/child relationships
* _Suggested model improvement_: precise the stereochemistry of the parent identifier, as it might be the same metabolite as the child or a different one.

**COMP_MAP_OK**

**COMP_MAP_WARN**

**COMP_UNKNOWN**

**COMP_MERGE**

**COMP_TO_GENERIC**

* The original compartments have been replaced by generic compartments MNXD1 and MNXD2 as in the MNXref distribution.
* The network connectivity is likely to be lost

**REAC_MAP_OK**

**REAC_MAP_NEW**

**REAC_MAP_EMPTY**

* The original equation was converted into an 'empty' equation.
* The MNXref reconciliation contains a list of known empty reactions, most of them are acid-base and/or tautomerization reactions.
* Valid empty reaction (e.g. acid-base) may be reported under code.

**REAC_MAP_WARN**

**REAC_MERGE**

* Two or more different original reaction identifiers were mapped onto a single MNXref identifier.
* _Suggested model improvement_: first, merge the implied metabolites into a single one; Secondly, merge the reactions into a single one.
* _Nota Bene_: MNXref is agnostic with respect to reaction directions and places directionality constraints on top of the (undirected) equation, together with enzyme descriptions.

