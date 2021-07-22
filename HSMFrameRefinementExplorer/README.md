# HSMFrameRefinementExplorer
Construction and analysis of universal state space for HSM model using a model frame refinement approach.

The Frame Refinement Exploration tool produces a graph of symbolic states and their transitions, where a symbolic state comprises the mode (binding site association state) and the (partial) information on the binding site association/disassociation threshold ordering and possibly also substance affinity ordering with respect to the association/disassociation thresholds (where they are not recorded within the binding site association states); the transitions are marked by a binding site association or disassociation events.
Furthermore, the tool allows computing a factor graph of the initial state space graph by means of bisimulation: two symbolic states that can not be distinguished by their available transitions, are grouped together. 
The strongly connected components in the factor graph are computed, as well.
The prototype is implemented in Haskell. The main program Statespace.hs contains the algorithm implementation. The data model is also described as a collection of functions and is placed in auxiliary files to be imported during the main program execution.
For more details of the implementation and its running possibilities, consult the Frame Refinement Exploration Tool folder in this repository.
The results of running the algorithm with different model parameters are in the Results folder. 

## Running the program
Install Haskell from https://haskell.org/
The imports of Data.Map, Data.Set and Data.Graph are resolved within the standard Haskell installation.
The import of Data.Graph.SCC requires a Hackage module from https://hackage.haskell.org/package/GraphSCC-1.0.4/docs/Data-Graph-SCC.html

A copy of the Data.Graph.SCC module is provided within this repository for the convenience, it is a copy from https://hackage.haskell.org/package/GraphSCC-1.0.4/docs/Data-Graph-SCC.html
Open the main file Statespace.hs in the Haskell environment and type 'main' at the command prompt.

## Running options
The choices to be done for running the algorithm:
1.	change the import clause in the main program to import data from different models (e.g. Base_Red, Base_RedEssential, Base_BlueEssential).
2.	Within each model definition, a choice can be made between the zero-value and arbitrary initial state consideration. This is done by selecting one of the two definitions for the initlist functions in the model file.
3.	Choice of the function to be run: 
a.	Strongly connected components
b.	Bisimulation graph
c.	State list only

## The models
The “Red” model corresponds to functions P_ME and P_R, the proteins cI, cro and cII, and to the binding sites bCII-1, bOR1, bOR2 and bOR3. It can be demonstrated that in Lambda phage model the behavior of this fragment can not be influenced from outside.
The “Red Essential” model furthermore groups the binding sites bOR1 and bOR2 together into a single binding site, as these sites appear in functions determining protein growth always together.
The “Blue” model involves the “Red” model, plus the binding sites bOL3, bOL2 and bOL1, the function P_L and the protein N (although the growth/degradation of N is observed in the model, its binding to bN is not, as bN is not included into the fragment). As with the “Red” fragment, the “Blue” fragment can not be influenced from outside.
The “Blue Essential” model is obtained from the “Blue” model by:
1)	Grouping the binding sites bOR1 and bOR2 together (represented by bOR2)
2)	Grouping the binding sites bOL1 and bOL2 together (represented by bOL2)
3)	Eliminating the bOL3 binding site, as it does not influence any protein growth function.
The “Blue” model appears to be too large to analyze by a brute force algorithm considering all threshold orderings. The “Blue Essential” model, however, contains all relevant dependencies, as the “Blue” model, and it admits analysis possibilities already by the simple prototype, written in a functional programming language.

## The state graph encoding
A state is a vector of binding site “states”. The binding sites are enumerated by an index from 0 to 6 in the order:
bCII-1, bOR1, bOR2, bOR3, bOL1, bOL2, bOL3.
Some of the models use only a part of these binding sites.
A “state” of a binding site is a number:
-	For a binary binding site (e.g. bCII-1) 0 means a free site, 1 means an occupied site
-	For a binding site able to attract two proteins (e.g. cI and cro), values 0 and from 2 to 5 are used:
o	0 – a free binding site
o	2,3 – cI is attached (2 – the level of cro is below its attachment threshold, 3 – the level of cro is above its attachment threshold)
o	4,5 – cro is attached (4 – the level of cI is below its attachment threshold, 5 – the level of cI is above its attachment threshold)
An example initial state in the “Red” model is [(0,0),(1,0),(2,0),(3,0)].
Such a state is presented to the user as “0000+-+”, where “0000” is the sequence of the binding site states and “+-+” is a growth/degradation factor vector for the considered proteins (it is determined by the functions from the binding site state, however, is presented for the user perception convenience). The order of proteins in the growth/degradation sequence (as well as in the program code) comes from the pname function in the model code, and is cII, cI, cro, N.


