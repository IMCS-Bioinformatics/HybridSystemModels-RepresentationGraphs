//---------------------------------------------------------------------------


#pragma hdrstop
  
#include "HSM_graph.h"
#include "my_bcb.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

// some simple useful functions

// computes 2^n, since c++ lacks integer power function
int exp2(int n)
{
int exp2 = 1; for (int i = 0; i < n; i++) exp2 = exp2*2;
return exp2;
}

// computes 3^n, since c++ lacks integer power function
int exp3(int n)
{
int exp3 = 1; for (int i = 0; i < n; i++) exp3 = exp3*3;
return exp3;
}

// computes n! (although used only in other modules)
int factorial(int n)
{
if (n == 0) return 1;
return n*factorial(n-1);
}

// split int i into bit array of size dim (min dim=1, max dim=31)
void split_bits(int i,char *barray,int dim)
{
int res = i;
for (int i = dim-1; i>=0; i--) {barray[i] = res % 2; res = (res - (res % 2))/2;}
}

// return integer corresponding to bit representation contained in bit array of size dim (min dim=1, max dim=31)
int join_bits(char *barray,int dim)
{
int res = barray[0];
for (int i = 1; i < dim; i++) res = res*2+barray[i];
return res;
}

// compute the value of n-th function from Defsgf for mode argM
char comp_gf(int n,Mode *argM,HSMGraph *G)
{
/*
if (argM->Id == 1 && n == 0)
    {
    int x = argM->Id;
    }
*/
int dim = G->Ngf[n]; // number of arguments
char argbuf[100]; //character buffer for argument values (max 31)
for (int i = 0; i < dim; i++) argbuf[i] = argM->BS[G->Paramsgf[10*n+i]];
int argval = join_bits(argbuf,dim);
return (char) G->Defsgf[1024*n+argval];
}

HSMGraph::HSMGraph(MyString iHSMName,int iNgenes,MyString *igNames,int iNtf,int iNb2,MyString *ib2Names,int iNb3,MyString *ib3Names,int *iBStfs,int *iNgf,int *iParamsgf,int *iDefsgf)
{
// initialise model ************************************************************
IniCount = 0; // initialisation value (=0 for default constructure)
// general: genes, TFs and BS
HSMName = iHSMName; // model name
Ngenes = iNgenes; // genes
gNames = new MyString[Ngenes];
for (int i = 0; i < Ngenes; i++) gNames[i] = igNames[i];
Ntf = iNtf; // transcription factors
Nb2 = iNb2; // binary bs
b2Names = new MyString[Nb2];
for (int i = 0; i < Nb2; i++) b2Names[i] = ib2Names[i];
Nb3 = iNb3; // ternary bs
b3Names = new MyString[Nb3];
for (int i = 0; i < Nb3; i++) b3Names[i] = ib3Names[i];
// some test bits below
/*
Application->MessageBox(gNames[0].c_str(),"gene",0);
Application->MessageBox(gNames[1].c_str(),"gene",0);
Application->MessageBox(gNames[2].c_str(),"gene",0);
*/
// binding site <-> tf map
BStfs = new int[Nb2+2*Nb3];
for (int i = 0; i < Nb2+2*Nb3; i++) BStfs[i] = iBStfs[i];
// gene regulatory functions
Ngf = new int[Ngenes]; // number of regs for each gene (gene order as in gNames)
for (int i = 0; i < Ngenes; i++) Ngf[i] = iNgf[i];
Paramsgf = new int[Ngenes*10]; // lists of regulators for each of the genes
for (int i = 0; i < Ngenes*10; i++) Paramsgf[i] = iParamsgf[i];
Defsgf = new int[Ngenes*1024]; // definitions of reg functions for each of the genes
for (int i = 0; i < Ngenes*1024; i++) Defsgf[i] = iDefsgf[i];

// check function defs
/*
ofstream outf;
outf.open("out.txt");
char buf[5000];

for (int i = 0; i < 16; i++)
    {
    split_bits(i,buf,4);
    int val1 = f_pr(buf[0],buf[1],buf[2],buf[3]);
    int val2 = Defsgf[1024*1+i];
    int val3 = Defsgf[1024*2+i];
    outf << val1 << " " << val2 << " " << val3 << "\n";
    }
outf << "\n";
for (int i = 0; i < 128; i++)
    {
    split_bits(i,buf,7);
    int val1 = f_pme(buf[2],buf[3],buf[4],buf[5],buf[0],buf[1],buf[6]);
    int val2 = Defsgf[1024*0+i];
    outf << val1 << " " << val2 << " " << "\n";
    }
outf.close();
  */

// create arrays for state space ***********************************************
Nmodes = exp2(Nb2+2*Nb3); // number of modes in complete state space
NmodesA = exp2(Nb2)*exp3(Nb3); // number of active modes (just 3 modes per ternary site instead of 4) - these modes will be included in state space graph
NmodesZero = 0; // consistent modes reachable from zero (set to zero here, redefined in SCC analysis)
NmodesIni = 0; // consistent modes reachable from initial states with single "1' in bs (set to zero here, redefined in SCC analysis) 
Modes = new MyListL;

// set transcription factor ordering - create arrays and initialise to "default" order (as TFs appear in BStfs array)
THcount = new char[Ntf];
for (int i = 0; i < Ntf; i++) THcount[i] = 0; // number of thresholds for each TF
for (int i = 0; i < Nb2+2*Nb3; i++) {THcount[BStfs[i]]++; THcount[BStfs[i]]++;} // two thresholds for each TF binding site
// set maximal # of thersholds per BS
THmax = 0;
for (int i = 0; i < Ntf; i++) if (THcount[i] > THmax) THmax = THcount[i];
// create threshold ordering arrays
THtype = new char[Ntf*THmax]; // cannot directly create 2-dim array of non-constant size...
THorder = new char[Ntf*THmax]; // cannot directly create 2-dim array of non-constant size...
// set binding site types - 0 - dissociation, 1 - association (array of size Ntf*THmax)
// *** this part will remain fixed until we will allow overlapping as/dis intervals
for (int i = 0; i < Ntf*THmax; i++) THtype[i] = -1; // initialise all to -1 (could help for testing)
for (int i = 0; i < Ntf; i++) for (int j = 0; j < THcount[i]; j++)
    {
    // just set THtype=0 for odd and THtype=1 for even positions
    setTHtype(i,j,j % 2);
    }
// set binding site threshold order as they in BStfs array
// this might later be changed with HSMGraph::ResetThresholds(...) function
for (int i = 0; i < Ntf*THmax; i++) THorder[i] = -1; // initialise all to -1 (could help for testing)
// go through Bstfs and for TF i in pos j set THorder=j for positons 2*c_pos and 2*c_pos+1
// **** very careful with that **** assumes that bs are ordered the same as tfs + probably something else ** should be changed asap ***
for (int i = 0; i < Ntf; i++)
    {
    int c_pos = 0;
    for (int j = 0; j < Nb2+2*Nb3; j++) if (BStfs[j] == i)
        {
        setTHorder(i,2*c_pos,j);
        setTHorder(i,2*c_pos+1,j);
        c_pos++;
        }
    }
// **** very careful with that **** assumes that bs are ordered the same as tfs + probably something else ** should be changed asap ***


// go through mode array and set preliminary pointers to 1 to real modes and to 0 for placeholders for bogus ternary ones
// assumption that ternary occupies positions Nbs...Nbs+2*Nbs3-1
for (int i = 0; i < Nmodes; i++)
    {
    char barray[100]; // just make sufficiently large - probably 100 BS will be max...
    split_bits(i,barray,Nb2+2*Nb3);
    bool is_active = true;
    // set as active only these that does not contain both positions of ternary sites occupied (no two censecutive "1" for ternary ones)
    for (int i = 0; i < Nb3; i++) if (barray[Nb2+2*i] == 1 && barray[Nb2+2*i+1] == 1) is_active = false;
    if (is_active) Modes->Add((void*)1);
    else Modes->Add(0);
    }

// *** the part below is currently just 1-1 repeated in HSMGraph::ResetThresholds
// go through mode array once more and now create modes themselves, additionally count the consitent modes
NmodesC = 0;
for (int i = 0; i < Nmodes; i++) if (Modes->Items[i] != 0)
    {
    Modes->Items[i] = new Mode(i,this);
    Mode *M = (Mode*) Modes->Items[i];
    if (M->isConsistent) NmodesC++;
    }
Ntrans = 0; // total number of transitions
Max_tout = 0; // maximal number of outgoing transitions
Max_tin = 0;  // maximal number of outgoing transitions


// go through mode array once more and explicitly assign in-coming transitions
for (int i = 0; i < Nmodes; i++)  if (Modes->Items[i] != 0)
    {
    Mode *orig_Mode = (Mode*) Modes->Items[i];
    for (int j = 0; j < orig_Mode->Ntrout; j++)
        {
        int dest = (int) orig_Mode->Tout_m[j];
        Mode *dest_Mode = (Mode*) Modes->Items[dest];
        dest_Mode->Tin_m->Add((void*) orig_Mode->Id);
        dest_Mode->Ntrin++;
        }
    }
// go through once more and count transitions
for (int i = 0; i < Nmodes; i++)  if (Modes->Items[i] != 0)
    {
    Mode *M = (Mode*) Modes->Items[i];
    Ntrans = Ntrans+M->Ntrout;
    if (M->Ntrout > Max_tout) Max_tout = M->Ntrout;
    if (M->Ntrin > Max_tin) Max_tin = M->Ntrin;
    }
}

HSMGraph::~HSMGraph()
{
delete[] gNames; delete[] b2Names; delete[] b3Names;
delete[] BStfs;
delete[] Ngf; delete[] Paramsgf; delete[] Defsgf;
delete[] THcount; delete[] THtype; delete[] THorder;
for (int i = 0; i < Nmodes; i++) if (Modes->Items[i]) delete Modes->Items[i];
delete Modes;
}

// reset threshold values in THorder to the order given in iTHorder and recompute all Modes and transitions
void HSMGraph::ResetThresholds(char *iTHorder)
{
IniCount++; // initialisation value (+1 for each resetting)
NmodesZero = 0;
for (int i = 0; i < Ntf*THmax; i++)  {
THorder[i] = iTHorder[i];
char a = THorder[i];
char b = THorder[i];
}
// process only existing modes
//for (int i = 0; i < Nmodes; i++) if (Modes->Items[i] == 0) return;
// delete all existing modes and replace TList ppointers with 1
//for (int i = 0; i < Nmodes; i++) if (Modes->Items[i]) {delete Modes->Items[i]; Modes->Items[i] = (void*) 1;}
// the part below is currently 1-1 with initial default creation
NmodesC = 0;
for (int i = 0; i < Nmodes; i++) if (Modes->Items[i] != 0)
    {
    Mode *M = (Mode*) Modes->Items[i];
    // reset transitions
    M->Ntrout = 0; M->Ntrin = 0;
    //M->Tout_m->Clear(); M->Tout_g->Clear(); M->Tout_t->Clear();
    M->Tin_m->Clear();
    M->SetConsistency();
    M->ProcessTransitions();
    // if such, add to count of consistent modes
    if (M->isConsistent) NmodesC++;
    }
Ntrans = 0; // total number of transitions
Max_tout = 0; // maximal number of outgoing transitions
Max_tin = 0;  // maximal number of outgoing transitions


// go through mode array once more and explicitly assign in-coming transitions
for (int i = 0; i < Nmodes; i++)  if (Modes->Items[i] != 0)
    {
    Mode *orig_Mode = (Mode*) Modes->Items[i];
    for (int j = 0; j < orig_Mode->Ntrout; j++)
        {
        int dest = (int) orig_Mode->Tout_m[j];
        Mode *dest_Mode = (Mode*) Modes->Items[dest];
        dest_Mode->Tin_m->Add((void*) orig_Mode->Id);
        dest_Mode->Ntrin++;
        }
    }
// go through once more and count transitions
for (int i = 0; i < Nmodes; i++)  if (Modes->Items[i] != 0)
    {
    Mode *M = (Mode*) Modes->Items[i];
    Ntrans = Ntrans+M->Ntrout;
    if (M->Ntrout > Max_tout) Max_tout = M->Ntrout;
    if (M->Ntrin > Max_tin) Max_tin = M->Ntrin;
    }
}


void HSMGraph::ResetFuncDefs2(int f_def0,int f_def1/*,int f_def2*/) // reset 3 function values to that defined by inetgers f_def0,f_def1,f_def2
{
// currentlu mostly copy/paste from ResetThresholds
IniCount++; // initialisation value (+1 for each resetting)
NmodesZero = 0;
// just change the definition of the first 3 functions (this assumes Ngenes >= 3)
char f_def_buf[100]; // should be ok :)
int argN;
argN = exp2(Ngf[0]);
split_bits(f_def0,f_def_buf,argN);
for (int j = 0; j < argN; j++) Defsgf[1024*0+j] = (int) f_def_buf[j];
argN = exp2(Ngf[1]);
split_bits(f_def1,f_def_buf,argN);
for (int j = 0; j < argN; j++)
    {
    Defsgf[1024*1+j] = (int) f_def_buf[j];
    int x1 = Defsgf[1024*1+j];
      int x2 = Defsgf[1024*0+j]; 
      }
/*
argN = exp2(Ngf[2]);
split_bits(f_def2,f_def_buf,argN);
for (int j = 0; j < argN; j++)
    {
    Defsgf[1024*2+j] = (int) f_def_buf[j];
    }

    */
// process only existing modes
//for (int i = 0; i < Nmodes; i++) if (Modes->Items[i] == 0) return;
// delete all existing modes and replace TList ppointers with 1
//for (int i = 0; i < Nmodes; i++) if (Modes->Items[i]) {delete Modes->Items[i]; Modes->Items[i] = (void*) 1;}
// the part below is currently 1-1 with initial default creation
NmodesC = 0;
for (int i = 0; i < Nmodes; i++) if (Modes->Items[i] != 0)
    {
    Mode *M = (Mode*) Modes->Items[i];
    // reset transitions
    M->Ntrout = 0; M->Ntrin = 0;
    //M->Tout_m->Clear(); M->Tout_g->Clear(); M->Tout_t->Clear();
    M->Tin_m->Clear();
    M->SetConsistency();
    M->ProcessTransitions();
    // this is extra for function redefinition - change the activity states of all genes
    for (int ii = 0; ii < M->G->Ngenes; ii++) M->Genes[ii] = comp_gf(ii,M,M->G);
    // if such, add to count of consistent modes
    if (M->isConsistent) NmodesC++;
    }
Ntrans = 0; // total number of transitions
Max_tout = 0; // maximal number of outgoing transitions
Max_tin = 0;  // maximal number of outgoing transitions


// go through mode array once more and explicitly assign in-coming transitions
for (int i = 0; i < Nmodes; i++)  if (Modes->Items[i] != 0)
    {
    Mode *orig_Mode = (Mode*) Modes->Items[i];
    for (int j = 0; j < orig_Mode->Ntrout; j++)
        {
        int dest = (int) orig_Mode->Tout_m[j];
        Mode *dest_Mode = (Mode*) Modes->Items[dest];
        dest_Mode->Tin_m->Add((void*) orig_Mode->Id);
        dest_Mode->Ntrin++;
        }
    }
// go through once more and count transitions
for (int i = 0; i < Nmodes; i++)  if (Modes->Items[i] != 0)
    {
    Mode *M = (Mode*) Modes->Items[i];
    Ntrans = Ntrans+M->Ntrout;
    if (M->Ntrout > Max_tout) Max_tout = M->Ntrout;
    if (M->Ntrin > Max_tin) Max_tin = M->Ntrin;
    }
}

// service functions to skip explicit calculation of 2-dim array positions
char HSMGraph::getTHtype(int i,int j)
{
return (int) THtype[THmax*i+j];
}

char HSMGraph::getTHorder(int i,int j)
{
return (int) THorder[THmax*i+j];
}

void HSMGraph::setTHtype(int i,int j,char val)
{
THtype[THmax*i+j] = val;
}

void HSMGraph::setTHorder(int i,int j,char val)
{
THorder[THmax*i+j] = val;
}


Mode::Mode(int iId,HSMGraph *iG)
{
Id = iId;
G = iG;
Ntrout = 0; // start with 0, then see add the ones that are needed
Ntrin = 0;  // this will be optionally computed at processing stage
// set the values of binding sites: 0 or 1  - occuppied and non-occupied BS (this defined by bist of Id)
BS = new char[G->Nb2+(2*G->Nb3)]; // the values will be 0 or 1  - occuppied and non-occupied BS
BSwork = new char[G->Nb2+(2*G->Nb3)]; // additional work area for transition generation
split_bits(Id,BS,G->Nb2+2*(G->Nb3));
// set the values of genes (0 - inactive, 1 - active) depending from the current states of BS (defined by bits of Id)
Genes = new char[G->Ngenes]; // the values will be 0 or 1 indicated non-active and active genes
for (int i = 0; i < G->Ngenes; i++) Genes[i] = comp_gf(i,this,G);
// create TLists for transition information
Tout_m = new int[G->Ngenes]; // transitions (outgoing) - destination mode Id
Tout_g = new int[G->Ngenes]; // transitions (outgoing) - number of gene activated or disactivated (as in genes)
Tout_t = new int[G->Ngenes]; // transitions (outgoing) - type (0 - activation, 1- disactivation)
Tin_m = new MyListS; // transitions (incoming) - source mode Id
// call SetConsistency function
SetConsistency();
// call transition computing function
ProcessTransitions();
}

// set consistency status according to the current threshold ordering in G
void Mode::SetConsistency()
{
// set the state of mode consistency with threshold ordering
isConsistent = 1; // assume that yes and then check :)
// set isConsistent to no, if for some gene there is inactive bs such that:
// 1) bs is not a ternary site occupied by alternative TF
// 2) there is an active bs with higher threshold ordering
// *** breakpoint for testing specific mode
 /*
 if (Id == 84 && G->IniCount == 11)
    {
    int xx  = 0;
    }
 */
for (int i = 0; i < G->Ntf; i++)
    {
    for (int j = 1; j < G->THcount[i]; j++)
        {
        int check_pos = G->getTHorder(i,j-1);
        if (BS[check_pos] == 1) continue; // not inactive
        if (check_pos >= G->Nb2 && (((check_pos-G->Nb2) % 2 == 0 && BS[check_pos+1] == 1) || ((check_pos-G->Nb2) % 2 == 1 && BS[check_pos-1] == 1))) continue; // active as ternary
        for (int k = j; k < G->THcount[i]; k++)
            {
            int high_pos = G->getTHorder(i,k);
            if (BS[high_pos] == 1) isConsistent = 0; // higher and active found
            // *** the first inequelity seems extra bogus requirement for non-consistency
            //if (high_pos > check_pos && BS[high_pos] == 1) isConsistent = 0; // higher and active found
            }
        }
    }
}

// create transitions according to the current threshold ordering
void Mode::ProcessTransitions()
{
// compute all the transitions
for (int i = 0; i < G->Ntf; i++)
    {
    // leave this in place if additional behavior tests for specific genes needed
     /*
    if (Id == 38232 && G->IniCount == 1 && i == 2)
        {
        int aaa = 0;
        }
     */
    int site_pos = -1; // site position in BS which is the next in the ordering to be activated or disactivated (-1 - undefined)
    // gene is active and growing
    if (Genes[i] == 1)
        {
        for (int j = 0; j < G->THcount[i]; j++)
            if (G->getTHtype(i,j) == 1 && BS[G->getTHorder(i,j)] == 0)
                {
                site_pos = G->getTHorder(i,j);
                int x1 = G->THorder[24]; int x6 = G->THorder[25]; int x3 = G->THorder[26];
                // try to create transition in this position, if fails (for ternary sites) go to the next; if succeeds, the break
                int sibling_site = -1; // asumme binary and then check for ternary site, if this is the case
                if (site_pos >= G->Nb2) // ternary in this case
                    {
                    if ((site_pos-G->Nb2) % 2 == 0) sibling_site = site_pos+1; // sibling is the next
                    if ((site_pos-G->Nb2) % 2 == 1) sibling_site = site_pos-1; // sibling is the previous
                    }
                if (CreateTrans(i,site_pos,sibling_site,Genes[i]) == 1) break;
                }
        }

    // gene is inactive and decreasing
    if (Genes[i] == 0) // disactive and decreasing
        {
        for (int j = G->THcount[i]-1; j >= 0; j--)
            if (G->getTHtype(i,j) == 0 && BS[G->getTHorder(i,j)] == 1)
                {
                site_pos = G->getTHorder(i,j);
                // create transition in this position
                int sibling_site = -1; // asumme binary and then check for ternary site, if this is the case
                if (site_pos >= G->Nb2) // ternary in this case
                    {
                    if ((site_pos-G->Nb2) % 2 == 0) sibling_site = site_pos+1; // sibling is the next
                    if ((site_pos-G->Nb2) % 2 == 1) sibling_site = site_pos-1; // sibling is the previous
                    }
                if (CreateTrans(i,site_pos,sibling_site,Genes[i]) == 1) break; // this should always succeed for disactivation
                }
        }
    }
}

int Mode::CreateTrans(int g,int s,int ss,int t) // create new transition for gene g, next site s and sibling site ss (-1 for binary) and type t (0 - disactivate, 1 - activate)
{
// leave this in place if additional behavior tests for specific genes needed
/*
if (Id == 18 && g == 2)
    {
    int aaa = 0;
    }
*/
// copy activity array
for (int i = 0; i < G->Nb2+(2*G->Nb3); i++) BSwork[i] = BS[i];
// activate or diasactivate site s
BSwork[s] = t;
// binary case - just add transition
if (ss == -1)
    {
    int destId = join_bits(BSwork,G->Nb2+(2*G->Nb3));
    Tout_m[Ntrout] = destId;
    Tout_g[Ntrout] = g;
    Tout_t[Ntrout] = t;
    //Tout_m->Add((void*) destId);
    //Tout_g->Add((void*) g);
    //Tout_t->Add((void*) t);
    Ntrout++;
    }
// ternary case
else
    {
    if (BSwork[s] == 1 && BSwork[ss] == 1) return 0; // do nothing if both site and sibling will become activated (transition to impossible state)
    // for site disactivation check whether sibling needs to be activated
    if (t == 0 && BSwork[ss] == 0)
        {
        int sg = G->BStfs[ss]; // sibling gene
        // find sites order position
        int ss_pos = -1;
        for (int i = 0; i < G->THcount[sg]; i++)
        if (G->getTHorder(sg,i) == ss && G->getTHtype(sg,i) == 1) ss_pos = i;
        bool activate_sibling = false;
        for (int i = ss_pos+1; i < G->THcount[sg]; i++)
            {
            int bs_numb = G->getTHorder(sg,i);
            int xx = BS[bs_numb];
            if (BS[bs_numb] == 1) activate_sibling = true;
            }
        if (activate_sibling) BSwork[ss] = 1;
        }
    // add transition now
    int destId = join_bits(BSwork,G->Nb2+(2*G->Nb3));
    Tout_m[Ntrout] = destId;
    Tout_g[Ntrout] = g;
    Tout_t[Ntrout] = t;
    //Tout_m->Add((void*) destId);
    //Tout_g->Add((void*) g);
    //Tout_t->Add((void*) t);
    Ntrout++;
    }
return 1;
}





