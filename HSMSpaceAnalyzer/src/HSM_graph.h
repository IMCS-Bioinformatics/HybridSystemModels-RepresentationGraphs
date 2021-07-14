//---------------------------------------------------------------------------

#ifndef HSM_graphH
#define HSM_graphH
//---------------------------------------------------------------------------

#include "my_bcb.h"


class Mode;
class HSMGraph;

// some simple useful functions

int exp2(int n); // computes 2^n, since c++ lacks integer power function
int exp3(int n); // computes 3^n, since c++ lacks integer power function
int factorial(int n); // computes n! (although used only in other modules)
void split_bits(int i,char *barray,int dim);  // split int i into bit array of size dim (min dim=1, max dim=31)
int join_bits(char *barray,int dim); // return integer corresponding to bit representation contained in bit array of size dim (min dim=1, max dim=31)
char comp_gf(int n,int argval,HSMGraph *G); // compute the argval value of n-th function from Defsgf
char comp_gf(int n,Mode *argM,HSMGraph *G); // compute the value of n-th function from Defsgf for mode argM

// these functions should be dynamic + also passed as parameters to constructors
/*
char f_pme(char x11,char x12,char x21,char x22,char x31,char x32,char y);
char f_pr(char x11,char x12,char x21,char x22);
*/



class Mode
{
public:
    Mode(int iId,HSMGraph *iG);
    ~Mode() {delete[] Genes; delete[] BS; delete[] BSwork; delete[] Tout_m; delete[] Tout_g; delete[] Tout_t; delete Tin_m;}
    void SetConsistency(); // set consistency status according to the current threshold oredering in G
    void ProcessTransitions(); // create transitions according to the current threshold ordewring
    int CreateTrans(int g,int s,int ss,int t); // a helper function - create new transition for gene g, next site s and sibling site ss (-1 for binary) and type t (0 - disactivate, 1 - activate)
    // variable declarations below
    int isConsistent; // mark whether the mode is consistent with the current ordering =1 if yes, =0 if no
    int Id; // id of this mode
    HSMGraph *G; // pointer to graph to which mode belongs
    int Ntrout; // number of outgoing transitions; to be computed on initialisation
    int Ntrin; // number of in-coming transitions; to be computed on initialisation
    char *Genes; // array showing growth or decreseasing of gene concentrations, size = Ngenes
    char *BS; // array with binding site states, size = Nbs+2*Nbs3
    char *BSwork ; /// work area for transition generation
    // max number of out transitions corresponds to the number of genes
    int *Tout_m; // transitions (outgoing) - destination mode Id
    int *Tout_g; // transitions (outgoing) - number of gene activated or disactivated (as in genes)
    int *Tout_t; // transitions (outgoing) - type (0 - activation, 1- disactivation)

    MyListS *Tin_m; // transitions (incoming) - source mode Id
};

class HSMGraph
{
public:
    HSMGraph(MyString iHSMName,int iNgenes,MyString *igNames,int iNtf,int iNb2,MyString *ib2Names,int iNb3,MyString *ib3Names,int *iBStfs,int *iNgf,int *iParamsgf,int *iDefsgf);
    ~HSMGraph();
    void ResetThresholds(char *iTHorder); // reset threshold values in THorder to the order given in iTHorder and recompute all Modes and transitions
    void ResetFuncDefs2(int f_def0,int f_def1/*,int f_def2*/); // reset 3 function values to that defined by inetgers f_def0,f_def1,f_def2  (this assumes that there are >= 3 genes)
    // service functions to skip explicit calculation of 2-dim array positions
    char getTHtype(int i,int j);
    char getTHorder(int i,int j);
    void setTHtype(int i,int j,char val);
    void setTHorder(int i,int j,char val);
    // variable declarations below
    MyString HSMName; // model name
    int Ngenes; // number of genes
    int Ntf; // number of TF-s (the Ntf of genes are also TFs)
    int Nb2; // number of normal binary binding sites
    int Nb3; // number of ternary bs represented by pairs - currently: bOR1[cI],bOR1[cro],bOR2[cI],bOR2[cro],bOR3[cI],bOR3[cro]
    MyString *gNames; // gene names
    MyString *b2Names; // binary bs names
    MyString *b3Names; // ternary bs names
    int *BStfs; // references to TFs linked to specific bs (the array of size Nb2+2*Nb3)
    int *Ngf; // number of regulators for each gene (size=Ngenes)
    int *Paramsgf; // indexes of bs in the right order that are parameters of the regulatory function (size=Ngenes, each subarray size 10)
    int *Defsgf; // definitions of regulatory functions (size=Ngenes, each subarray size 1024)
    char *THcount; // number of thresholds for each TF (number of times it is contained in BStfs *2)
    int THmax; // max # of thresholds per binding site
    char *THtype; // types of BS - 0 - dissociation, 1 - association (array of size Ntf*THmax)
    char *THorder; // threshold orders of BS
    int Nmodes; // number of modes; should be equal to 2^(Nbs+2*Nbs3)
    int NmodesA; // number of active modes, equal to 2^Nb2*3^Nb3  - these will be included in state space graph
    int NmodesC; // number of modes consistent with threshold ordering (subset of active ones)
    int NmodesZero; // consistent modes reachable from zero (set to zero here, redefined in SCC analysis)
    int NmodesIni; // consistent modes reachable from initial states with single "1' in bs (set to zero here, redefined in SCC analysis)
    MyListL *Modes; // should be TList or smth that allows both real and 0 pointers
    int Ntrans; // total number of transitions
    int Max_tout; // maximal number of outgoing transitions
    int Max_tin;  // maximal number of outgoing transitions
    // Reintitialisation number: =0 for default constructor, increase +1 by each re-in itialisation
    int IniCount;
};
#endif
