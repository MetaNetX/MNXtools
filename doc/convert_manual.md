Report of conversion to MNXref are produced in YAML format. 

The reported observations are centered on the identifiers given in the source model.



The following code are produced:

**CHEM_XREF_OK** 

* This code is only produced when cross-refs are exploited
* There is a single cross-ref or there are multiple cross-refs that correspond to a single MNXref identifier

**CHEM_XREF_CONFLICT**

* This code is only produced when cross-refs are exploited
* There are multiple cross-refs that correspont to different unrelated MNXref identifiers
* _Suggested model improvement_: select which cross ref is the best one, keep it and move the others elswhere (in comments, for example) to avoid recreating this conflict latter.

CHEM_XREF_AMBIGUOUS

* This code is only produced when cross-refs are exploitied
* There are multiple cross-refs that correspont to different MNXref identifiers, but these are related by stereochemical parent/child relationships.
* The parent MNXref identifier is selected automatically to preserve chemistry, at the expense of precision.
* _Suggested model improvement_: select which cross ref is the best one, possibly among those with the more detailed stereochemistry. It ihowever cannot be excluded that the parent cross ref correspond to a different metabolites, which precise stereochistry could be specified from a different xref.

COMP_MAP_OK

* The metabolite was succesfully mapped to a single identifier of MNXref.

CHEM_MAP_UNKNOWN

* The metabolite cannot be mapped to an MNXref identifier. 
* The prefix UNK: was added to the metabolite identifer.

CHEM_MAP_WARN

* 

CHEM_PROTECT


CHEM_MERGE

* Two or more different original identfiers were mapped onto a single MNXref identifier.
* _Suggested model improvement_: merge the original identifiers into a single metabolite. 

CHEM_ISOMERIC

* Two or more different MNXref identifers are related by isomeric parent/child relationships
* _Suggested model improvement_: precise the sterochistry of the parent identifier, as it might be the same metabolite as the child or a different one.

COMP_MAP_OK

COMP_MAP_WARN

COMP_UNKNOWN

COMP_TO_GENERIC

* The original compartments have been replaced by generic compartments MNXD1 and MNXD2 as in the MNXref distribution. 
* The network connectivity is likely to be lost

REAC_MAP_OK

REAC_MAP_NEW

REAC_MAP_EMPTY

* The original equation was converted into an 'empty' equation 
* The MNXref reconciliatoin contains a list of known empty reactions, most of them are acid-base and/or tauterization reactions. 
* Valid empty reaction (e.g. acid-base) may be reported under code 

REAC_MAP_WARN

REAC_MERGE

* Two or more different original reaction identifier were mapped onto a single MNXref identifier.
* _Suggested model improvement_: first, merge the original metabilte identifiers into a single metabolite (if any). Secondly, merge the reactions into a single one.
* (SBML) Nota Bene: MNXref is agnostic with respect to reaction directions and place direcionality constraints on top of (undirected) equation, together with enzyme descriptions.  

