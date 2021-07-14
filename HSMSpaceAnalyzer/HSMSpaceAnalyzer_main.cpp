//---------------------------------------------------------------------------

#include <iostream>
using namespace std;

#pragma hdrstop

//#include "HSMSpaceAnalyzer_main.h"
#include "HSM_graph.h"
#include "graph_SCC.h"
#include "HSM_io.h"
#include "my_bcb.h"

//---------------------------------------------------------------------------

#pragma argsused



// *** a range of global variables *********************************************

// specify model here by un-commenting the currently used value
MyString ModelFileName = MyString(""); // initialise from empty, use the first argument from commant line
// the visualisation level: if >= create GEW, if >= 2 save as picture
// does not have meaning without visual components
//int VisLevel = 2;


// *** the values below might be useful to edit by hand befor compilation ******
// "consistency" threshold for drawn part of state space: 1 - all consistent, 3 - consistent and reachable from 0
// 2 - consistent and reachable from ini states with single "1" at bs
int ConsThreshold = 1;
// use the latter for actual drawing (might be reset for non-zero Flag values)
int TempConsThreshold = 1;
int Flag = 0;  // use this for "forced" output on "interesting" state spaces
// input and output file locations
MyString ModelFileDirectory = MyString("HSM_Models/");
MyString ModelStateSpaceDirectory = MyString("HSM_StateSpaces/");

// *** if need this can be edited to generate state spaces with idxes from-to **
// *** could be needed or useful for large or numerous SSs *********************
// put Min and Max model numbers here for convenience
int MinMOD = 0;
int MaxMOD = 0; // 0 means no max constraints
int CheckORD = 0; // should be 1 for models with enforced constraints, 0 for others
int save_numb = -1; // order number of saved state spaces

// HSM graph

HSMGraph *G = 0;

// matrix graph representations corresponding to generated HSMgraph G (G_orig - the original graphs, G_trans - the transposed)
int *G_orig = 0;
int *G_trans = 0;
// for convenience - numbers of out edges from each vertex in G_orig and G_trans
int *GN_orig = 0;
int *GN_trans = 0;
// information about SCCs
int CNumb; // number of components
int *Cidx = 0; // component index (DFS tree root) in which concrete vertex is contained
int *Cidx2 = 0; // DFS tree root from which vertex has been reached first ( = 0 for values reached from 0 state)
int *Csize = 0; // number of vertices in component to which this concrete vertex belongs
MyListL *CompList = 0; // pointers to Cnumb int arrays containing lists of all components
MyListL *CLsizes; // sizes of component lists contained in CompList
MyListL *CLtypes; // types of components: 1 for attractors, 0 for temporary, 2 for attractors with outgoing edges (temporality can not be proved)

// *** declared functions ******************************************************

int findidx(int *pos, int size,int val);
MyString FormatSaveNumber(int n);

// generate all pn permutations outA of n integers from array inA, return number of them (should coinside with pn)
// n_max denotes the value of the first call and is needed to know the dimension of the outA array
// out_pos is the position in outA from which recursive calls start filling values (0 for initial call)
int AllPermutations(int n_max,int n,int pn,int out_pos,int *inA,int *outA);

// ranamed functions initially associated with button clicks and other events **

void GenerateModel();

// *** main function ***********************************************************

int main(int argc, char* argv[])
{
if (argc < 2) {cout << "Provide HSM file name as argument!!!"; return -1;}
ModelFileName = argv[1];
GenerateModel(); // all the work is being done here
return 0;
}
//---------------------------------------------------------------------------

// *** function implementations below ******************************************
void GenerateModel()
{
// *****************************************************************************
// initialise HSM model
if (G) delete G;

MyString HSMName;  // model name
int Ngenes; // number of genes
MyString *gNames; // list of gene names
int Ntf; // transcription factors - the first Ntf genes are also TFs (should have Ntf <= Ngenes)
int Nb2; // number of binary binding sites
MyString *b2Names; // list of binary binding sites
int Nb3; // number of ternary binding sites
MyString *b3Names; // list of ternary binding sites
// specify correspondence of binding sites to TFs
// the length of BStfs and array is Nb2+2*Nb3 and contains references of TFs that bind to each of the BS
int *BStfs;
// specify gene regulatory functions
int *Ngf; // number of regulators (# of bs) for each of genes, the size of Ngf should be Ngenes (gene order as in gNames)
// for each gene specify the indexes of bs (as in BStfs) that acts as its regulators
// for each gene g there should be Ngf[g] indexes
int *Paramsgf; // assume max 10 arguments - size Ngenes*10
// specify regulatory functions themselves (each described by array of 0 and 1 of size 2^Ngf[i])
int *Defsgf; // assume max 2^10 arguments
// read the model from file
int read_res = 0;


MyString ModelFileString =  ModelFileDirectory + ModelFileName;
// read the model form file
read_res = read_HSM(ModelFileString.c_str(),HSMName,Ngenes,gNames,Ntf,Nb2,b2Names,Nb3,b3Names,BStfs,Ngf,Paramsgf,Defsgf);
if (read_res) cout << "HSM model read succcessfully!";
else {cout << "Failed to read model!!!"; return;}


// create model
G = new HSMGraph(HSMName,Ngenes,gNames,Ntf,Nb2,b2Names,Nb3,b3Names,BStfs,Ngf,Paramsgf,Defsgf);

// *** non-GUI version specific only - enforce constrain checking for specific model names
if (G->HSMName == MyString("Lambda_Core_blue") || G->HSMName == MyString("Lambda_Complete") || G->HSMName == MyString("Lambda_Oppenheim") || G->HSMName == MyString("HK22_Complete"))
CheckORD = 1;


// clean up specification data (if needed)
delete[] Paramsgf; delete[] Defsgf;
// write out model initialisation data
cout << "HSM graph created!";
cout << "Model name:\t\t" << G->HSMName.c_str() << "\n";
cout << "Number of genes:\t\t" << G->Ngenes;
cout << "Gene list:\t\t\t";
for (int i = 0; i < G->Ngenes; i++) cout << G->gNames[i].c_str() << " ";
cout << "\n";
cout << "Number of TFs:\t\t" <<  G->Ntf << "\n";
cout << "TF list:\t\t\t";
for (int i = 0; i < G->Ntf; i++) cout << G->gNames[i].c_str() << " ";
cout << "\n";
cout << "Number of binary BS:\t" << G->Nb2 << "\n";
cout << "B2 list:\t\t\t";
for (int i = 0; i < G->Nb2; i++) cout << G->b2Names[i].c_str() << " ";
cout << "\n";
cout << "Number of ternary BS:\t" << G->Nb3 << "\n";
cout << "B3 list:\t\t\t";
for (int i = 0; i < G->Nb3; i++) cout << G->b3Names[i].c_str() << " ";
cout << "\n\n";
// *****************************************************************************

// create arrays and Tlists for components here sine these can be reused
// just re-initialise array values and empty TList for each new ordering
if (Cidx) delete[] Cidx; if (Cidx) delete[] Cidx2; if (Csize) delete[] Csize;
if (CompList) {for (int i = 0; i < CNumb; i++) delete[] (int*) CompList->Items[i]; delete CompList;}
CompList = new MyListL;
if (CLsizes) delete CLsizes; if (CLtypes) delete CLtypes;
CLsizes = new MyListL; CLtypes = new MyListL;
Cidx = new int[G->Nmodes]; // component index (DFS tree root) in which concrete vertex is contained
Cidx2 = new int[G->Nmodes]; // DFS tree root from which vertex has been reached first ( = 0 for values reached from 0 state)
Csize = new int[G->Nmodes]; // number of vertices in component to which this concrete vertex belongs


// now try to re-initialise according to all the possible BS orderings
// assume non-overlapping intervals as yet
// also permutations are hard coded for specific models
// assume max 10 BS for a gene (ini_cro,ini_cI,ini_cII sizes = 10) , and max 20 genes (THorder size 400 = Ngenes*2*#BS per gene)

int n_cI = -1; int n_cro = -1; int n_cII = -1;
int pos_cI = -1; int pos_cro = -1; int pos_cII = -1;
// allocate constraint arrays, assume 1024 elements will be enough
int *ini_cI = new int[1024];
int *ini_cro = new int[1024];
int *ini_cII = new int[1024];
int THsize = G->Ntf*G->THmax;
char *THorder = new char[THsize];
// lambda red
if (G->HSMName == MyString("Lambda_Core_red"))
    {
    n_cI = 3;
    n_cro = 3;
    n_cII = 1;
    int ini_cI_01[] = {1,3,5}; //cout << "\n aa1  " << ini_cI_01[0] << "\n\n";
    int ini_cro_01[] = {2,4,6};
    int ini_cII_01[] = {0};
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_01[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_01[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_01[i];
    // offset positions of cI and cro in THorder
    pos_cI = 0;
    pos_cro = 12;
    pos_cII = 6;
    char THorder_01[] = {1,1,3,3,5,5,0,0,-1,-1,-1,-1,6,6,4,4,2,2};
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_01[i];
    }
    //cout << "\n bb4 " << ini_cI[0] << "\n\n";
// lambda blue
if (G->HSMName == MyString("Lambda_Core_blue"))
    {
    n_cI = 6;
    n_cro = 6;
    n_cII = 1;
    int ini_cI_02[] = {1,3,5,7,9,11};
    int ini_cro_02[] = {2,4,6,8,10,12};
    int ini_cII_02[] = {0};
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_02[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_02[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_02[i];
    // offset positions of cI and cro in THorder
    pos_cI = 0;
    pos_cro = 24;
    pos_cII = 12;
    char THorder_02[] = {1,1,3,3,5,5,7,7,9,9,11,11,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,12,12,10,10,8,8,6,6,4,4,2,2};
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_02[i];
    }
    //cout << "\n bb5 " << ini_cI[0] << "\n\n";
if (G->HSMName == MyString("Lambda_Complete"))
    {
    n_cI = 6;
    n_cro = 6;
    n_cII = 2;
    int ini_cI_03[] = {4,6,8,10,12,14};
    int ini_cro_03[] = {5,7,9,11,13,15};
    int ini_cII_03[] = {0,1};
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_03[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_03[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_03[i];
    // offset positions of cI and cro in THorder
    pos_cI = 0;
    pos_cro = 24;
    pos_cII = 12;
    char THorder_03[] = {4,4,6,6,8,8,10,10,12,12,14,14,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,15,15,13,13,11,11,9,9,7,7,5,5,2,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_03[i];
    }
    //cout << "\n bb6 " << ini_cI[0] << "\n\n";
if (G->HSMName == MyString("Lambda_Oppenheim"))
    {
    n_cI = 6;
    n_cro = 6;
    n_cII = 3;
    int ini_cI_05[] = {5,7,9,11,13,15};
    int ini_cro_05[] = {6,8,10,12,14,16};
    int ini_cII_05[] = {0,1,2};
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_05[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_05[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_05[i];
    // offset positions of cI and cro in THorder
    pos_cI = 0;
    pos_cro = 24;
    pos_cII = 12;
    char THorder_05[] = {5,5,7,7,9,9,11,11,13,13,15,15,1,1,0,0,2,2,-1,-1,-1,-1,-1,-1,16,16,14,14,12,12,10,10,8,8,6,6,3,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_05[i];
    }
    //cout << "\n bb7 " << ini_cI[0] << "\n\n";
if (G->HSMName == MyString("HK22_Complete"))
    {
    n_cI = 6;
    n_cro = 6;
    n_cII = 2;
    int ini_cI_04[] = {3,5,7,9,11,13};
    int ini_cro_04[] = {4,6,8,10,12,14};
    int ini_cII_04[] = {0,1};
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_04[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_04[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_04[i];
    // offset positions of cI and cro in THorder
    pos_cI = 0;
    pos_cro = 24;
    pos_cII = 12;
    char THorder_04[] = {3,3,5,5,7,7,9,9,11,11,13,13,1,1,0,0,-1,-1,-1,-1,-1,-1,-1,-1,14,14,12,12,10,10,8,8,6,6,4,4,2,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_04[i];
    }
    //cout << "\n bb8 " << ini_cI[0] << "\n\n";
if (G->HSMName == MyString("Loops_2"))
    {
    // for loops test this is a dummy part....
    n_cI = 1;
    n_cro = 1;
    n_cII = 1;
    int ini_cI_15[] = {0};
    int ini_cro_15[] = {1};
    int ini_cII_15[] = {2};
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_15[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_15[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_15[i];
    // offset positions of cI and cro in THorder
    pos_cI = 0;
    pos_cro = 2;
    pos_cII = 4;
    char THorder_15[] = {0,0,1,1};
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_15[i];
    }
    //cout << "\n bb9 " << ini_cI[0] << "\n\n";
if (G->HSMName == MyString("Circadian_03") || G->HSMName == MyString("Circadian_03_2"))
    {
    // for circadian cycles this is a dummy part....
    n_cI = 1;
    n_cro = 1;
    n_cII = 1;
    int ini_cI_06[] = {0};
    int ini_cro_06[] = {1};
    int ini_cII_06[] = {2};
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_06[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_06[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_06[i];
    // offset positions of cI and cro in THorder
    pos_cI = 0;
    pos_cro = 2;
    pos_cII = 4;
    char THorder_06[] = {0,0,1,1,2,2}; // changed for 4 genes!
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_06[i];
    }
    //cout << "\n bb10 " << ini_cI[0] << "\n\n";
    if (G->HSMName == MyString("Circadian_05"))
    {
    // for circadian cycles this is a dummy part....
    n_cI = 1;
    n_cro = 1;
    n_cII = 1;
    int ini_cI_07[] = {0};
    int ini_cro_07[] = {1};
    int ini_cII_07[] = {2};
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_07[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_07[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_07[i];
    // offset positions of cI and cro in THorder
    pos_cI = 0;
    pos_cro = 2;
    pos_cII = 4;
    char THorder_07[] = {0,0,1,1,2,2,3,3,4,4};
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_07[i];
    }
    //cout << "\n bb11 " << ini_cI[0] << "\n\n";
    if (G->HSMName == MyString("Mu_3_01"))
    {
    // reuse cI, cII and cro for permutation count and constraints
    n_cI = 3; // refers to c
    n_cro = 2; // refers to ner
    n_cII = 1; // not used
    int ini_cI_08[] = {2,3,4};
    int ini_cro_08[] = {0,1};
    int ini_cII_08[] = {0}; // not used, but use dummy array for permutation generation
    for (int i = 0; i < n_cI; i++) ini_cI[i] = ini_cI_08[i];
    for (int i = 0; i < n_cro; i++) ini_cro[i] = ini_cro_08[i];
    for (int i = 0; i < n_cII; i++) ini_cII[i] = ini_cII_08[i];
    // offset positions of cI and cro in THorder
    pos_cI = 6; // c
    pos_cro = 0;  // ner
    pos_cII = 0; // not used
    pos_cII = 0; // not used
    char THorder_08[] = {0,0,1,1,-1,-1,2,2,3,3,4,4};
    for (int i = 0; i < THsize; i++) THorder[i] = THorder_08[i];
    }
    //cout << "\n bb12 " << ini_cI[0] << "\n\n";
int pn_cI = factorial(n_cI);
int pn_cro = factorial(n_cro);
int pn_cII = factorial(n_cII);
int *p_cI = new int[n_cI*pn_cI];
int *p_cro = new int[n_cro*pn_cro];
int *p_cII = new int[n_cII*pn_cII];
// cout << "\n cc3  " << ini_cI[0] << "\n\n";
AllPermutations(n_cI,n_cI,pn_cI,0,ini_cI,p_cI);

// *** test printout ***
/*
cout << "pn_cI *******************************************************************\n";
 cout << "\n dd4  " << ini_cI[0] << "\n\n";
cout << n_cI << " " << pn_cI << "\n\n";
for (int i = 0; i < n_cI; i++) cout << ini_cI[i] << " ";
cout << "\n\n";
for (int i = 0; i < pn_cI; i++)
{
    for (int j = 0; j < n_cI; j++) cout << p_cI[n_cI*i+j] << " ";
    cout << "\n";
}
cout << "pn_cI *******************************************************************\n";
*/

AllPermutations(n_cro,n_cro,pn_cro,0,ini_cro,p_cro);
AllPermutations(n_cII,n_cII,pn_cII,0,ini_cII,p_cII);
// for cycle tests assign dummy values of permutation count to match the number of functions
if (G->HSMName == MyString("Loops_2")) {pn_cI = 16; pn_cro = 16; pn_cII = 1;}
// for circadian cycles assign dummy values of permutation count to match the number of functions
if (G->HSMName == MyString("Circadian_03")) {pn_cI = 16; pn_cro = 16; pn_cII = 1;}
if (G->HSMName == MyString("Circadian_03_2")) {pn_cI = 1; pn_cro = 1; pn_cII = 1;}
cout << "Permutation count:" << pn_cI << " " << pn_cro << " " << pn_cII << "\n";   
save_numb = 0; // # of generated and saved model
for (int ii = 0; ii < pn_cI; ii++)
    {
    // check consistency for cI ordering
    bool cI_ok = false;
    //if (p_cI[n_cI*i+0] < p_cI[n_cI*i+1] &&  p_cI[n_cI*i+1] < p_cI[n_cI*i+2] && p_cI[n_cI*i+3] < p_cI[n_cI*i+4] && p_cI[n_cI*i+4] < p_cI[n_cI*i+5]) cI_ok = true;
    // currently quick and stupid workaround - find idx for each site position
    int cI_01,cI_02,cI_03,cI_04,cI_05,cI_06;
    int cro_01,cro_02,cro_03,cro_04,cro_05,cro_06;
    int cII_01,cII_02,cII_03;
    // red model
    if (G->HSMName == MyString("Lambda_Core_red"))
        {
        cI_01 = findidx(p_cI+n_cI*ii,n_cI,1);
        cI_02 = findidx(p_cI+n_cI*ii,n_cI,3);
        cI_03 = findidx(p_cI+n_cI*ii,n_cI,5);
        // just smth to pass the constraint check
        cI_04 = 3;
        cI_05 = 4;
        cI_06 = 5;
        }
    // blue model
    if (G->HSMName == MyString("Lambda_Core_blue"))
        {
        cI_01 = findidx(p_cI+n_cI*ii,n_cI,1);
        cI_02 = findidx(p_cI+n_cI*ii,n_cI,3);
        cI_03 = findidx(p_cI+n_cI*ii,n_cI,5);
        cI_04 = findidx(p_cI+n_cI*ii,n_cI,7);
        cI_05 = findidx(p_cI+n_cI*ii,n_cI,9);
        cI_06 = findidx(p_cI+n_cI*ii,n_cI,11);
        }
    // full model
    if (G->HSMName == MyString("Lambda_Complete"))
        {
        cI_01 = findidx(p_cI+n_cI*ii,n_cI,4);
        cI_02 = findidx(p_cI+n_cI*ii,n_cI,6);
        cI_03 = findidx(p_cI+n_cI*ii,n_cI,8);
        cI_04 = findidx(p_cI+n_cI*ii,n_cI,10);
        cI_05 = findidx(p_cI+n_cI*ii,n_cI,12);
        cI_06 = findidx(p_cI+n_cI*ii,n_cI,14);
        }
    // Oppenheim model
    if (G->HSMName == MyString("Lambda_Oppenheim"))
        {
        cI_01 = findidx(p_cI+n_cI*ii,n_cI,5);
        cI_02 = findidx(p_cI+n_cI*ii,n_cI,7);
        cI_03 = findidx(p_cI+n_cI*ii,n_cI,9);
        cI_04 = findidx(p_cI+n_cI*ii,n_cI,11);
        cI_05 = findidx(p_cI+n_cI*ii,n_cI,13);
        cI_06 = findidx(p_cI+n_cI*ii,n_cI,15);
        }
    // HK22 model
    if (G->HSMName == MyString("HK22_Complete"))
        {
        cI_01 = findidx(p_cI+n_cI*ii,n_cI,3);
        cI_02 = findidx(p_cI+n_cI*ii,n_cI,5);
        cI_03 = findidx(p_cI+n_cI*ii,n_cI,7);
        cI_04 = findidx(p_cI+n_cI*ii,n_cI,9);
        cI_05 = findidx(p_cI+n_cI*ii,n_cI,11);
        cI_06 = findidx(p_cI+n_cI*ii,n_cI,13);
        }
    // Mu_3 model
    if (G->HSMName == MyString("Mu_3_01"))
        {
        // cI refers to c, constraint order as for cI
        cI_01 = findidx(p_cI+n_cI*ii,n_cI,2);
        cI_02 = findidx(p_cI+n_cI*ii,n_cI,3);
        cI_03 = findidx(p_cI+n_cI*ii,n_cI,4);
        // just smth to pass the constraint check
        cI_04 = 5;
        cI_05 = 6;
        cI_06 = 7;
        }
    if (cI_01 < cI_02 && cI_02 < cI_03 && cI_04 < cI_05 && cI_05 < cI_06) cI_ok = true;
    // experimentally allow also 1<3<2 and 2<3<1 and 2<1<3 for cI (default 1<2<3)
    //if (cI_01 < cI_03 && cI_03 < cI_02 && cI_04 < cI_05 && cI_05 < cI_06) cI_ok = true;
    //if (cI_02 < cI_03 && cI_03 < cI_01 && cI_04 < cI_05 && cI_05 < cI_06) cI_ok = true;
    //if (cI_02 < cI_01 && cI_01 < cI_03 && cI_04 < cI_05 && cI_05 < cI_06) cI_ok = true;
    if ((CheckORD == 1) && !cI_ok)
    {
        //cout << "cI" << ii << " " << CheckORD << " " << cI_ok << " " << cI_01 << " " << cI_02 << " " << cI_03 << " " << cI_04 << " " << cI_05 << " " << cI_06 << "\n";
        continue;
    }
    for (int jj = 0; jj < pn_cro; jj++)
    {
    // check consistency for cro ordering
    bool cro_ok = false;
    //if (p_cro[n_cro*j+0] > p_cro[n_cro*j+1] && p_cro[n_cro*j+1] > p_cro[n_cro*j+2] && p_cro[n_cro*j+3] > p_cro[n_cro*j+4] && p_cro[n_cro*j+4] > p_cro[n_cro*j+5]) cro_ok = true;
        // currently quick and stupid workaround - find idx for each site position
        // red model
    if (G->HSMName == MyString("Lambda_Core_red"))
        {
        cro_01 = findidx(p_cro+n_cro*jj,n_cro,2);
        cro_02 = findidx(p_cro+n_cro*jj,n_cro,4);
        cro_03 = findidx(p_cro+n_cro*jj,n_cro,6);
        // just smth to pass the constraint check
        cro_04 = 12;
        cro_05 = 10;
        cro_06 = 8;
        }
    // blue model
    if (G->HSMName == MyString("Lambda_Core_blue"))
        {
        cro_01 = findidx(p_cro+n_cro*jj,n_cro,2);
        cro_02 = findidx(p_cro+n_cro*jj,n_cro,4);
        cro_03 = findidx(p_cro+n_cro*jj,n_cro,6);
        cro_04 = findidx(p_cro+n_cro*jj,n_cro,8);
        cro_05 = findidx(p_cro+n_cro*jj,n_cro,10);
        cro_06 = findidx(p_cro+n_cro*jj,n_cro,12);
        }
    // full model
    if (G->HSMName == MyString("Lambda_Complete"))
        {
        cro_01 = findidx(p_cro+n_cro*jj,n_cro,5);
        cro_02 = findidx(p_cro+n_cro*jj,n_cro,7);
        cro_03 = findidx(p_cro+n_cro*jj,n_cro,9);
        cro_04 = findidx(p_cro+n_cro*jj,n_cro,11);
        cro_05 = findidx(p_cro+n_cro*jj,n_cro,13);
        cro_06 = findidx(p_cro+n_cro*jj,n_cro,15);
        }
    // Oppenheim model
    if (G->HSMName == MyString("Lambda_Oppenheim"))
        {
        cro_01 = findidx(p_cro+n_cro*jj,n_cro,6);
        cro_02 = findidx(p_cro+n_cro*jj,n_cro,8);
        cro_03 = findidx(p_cro+n_cro*jj,n_cro,10);
        cro_04 = findidx(p_cro+n_cro*jj,n_cro,12);
        cro_05 = findidx(p_cro+n_cro*jj,n_cro,14);
        cro_06 = findidx(p_cro+n_cro*jj,n_cro,16);
        }
    // HK22 model
    if (G->HSMName == MyString("HK22_Complete"))
        {
        cro_01 = findidx(p_cro+n_cro*jj,n_cro,4);
        cro_02 = findidx(p_cro+n_cro*jj,n_cro,6);
        cro_03 = findidx(p_cro+n_cro*jj,n_cro,8);
        cro_04 = findidx(p_cro+n_cro*jj,n_cro,10);
        cro_05 = findidx(p_cro+n_cro*jj,n_cro,12);
        cro_06 = findidx(p_cro+n_cro*jj,n_cro,14);
        }
    // Mu_3 model
    if (G->HSMName == MyString("Mu_3_01"))
        {
        // cro refers to ner, no constraints for now
        // just smth to pass the constraint check
        cro_01 = 11;
        cro_02 = 10;
        cro_03 = 9;
        cro_04 = 8;
        cro_05 = 7;
        cro_06 = 6;
        }
    if (cro_01 > cro_02 && cro_02 > cro_03 && cro_04 > cro_05 && cro_05 > cro_06) cro_ok = true;
    // experimentally allow also 3<1<2 for cro (default 3<1<2)
    //if (cro_02 > cro_01 && cro_01 > cro_03 && cro_04 > cro_05 && cro_05 > cro_06) cro_ok = true;
    if ((CheckORD == 1) && !cro_ok) continue;
    for (int kk = 0; kk < pn_cII; kk++)
    {
    // check consistency for cro ordering
    bool cII_ok = false;
    // currently quick and stupid workaround - find idx for each site position
    // red model
    if (G->HSMName == MyString("Lambda_Core_red"))
        {
        // just smth to pass the constraint check
        cII_01 = 1;
        cII_02 = 0;
        }
    // blue model
    if (G->HSMName == MyString("Lambda_Core_blue"))
        {
        // just smth to pass the constraint check
        cII_01 = 1;
        cII_02 = 0;
        }
    // full model
    if (G->HSMName == MyString("Lambda_Complete"))
        {
        cII_01 = findidx(p_cII+n_cII*kk,n_cII,0);
        cII_02 = findidx(p_cII+n_cII*kk,n_cII,1);
        }
    // Oppenheim model
    if (G->HSMName == MyString("Lambda_Oppenheim"))
        {
        cII_01 = findidx(p_cII+n_cII*kk,n_cII,0);
        cII_02 = findidx(p_cII+n_cII*kk,n_cII,1);
        }
    // HK22 model
    if (G->HSMName == MyString("HK22_Complete"))
        {
        cII_01 = findidx(p_cII+n_cII*kk,n_cII,0);
        cII_02 = findidx(p_cII+n_cII*kk,n_cII,1);
        }
    // Mu_3 model
    if (G->HSMName == MyString("Mu_3_01"))
        {
        // just smth to pass the constraint check
        cII_01 = 1;
        cII_02 = 0;
        }
    if (cII_01 > cII_02) cII_ok = true;
    if ((CheckORD == 1) && !cII_ok) continue;
    {
    // *** any general processing starts below *********************************
    save_numb++;
    Flag  = 0; // assume "non-interesting" model at the start
    //Application->MessageBoxA("M",MyString(MyString(i)+ " " + MyString(j)).c_str(),0); continue;
    //    if (save_numb > 61) break;
    // fill in threshold values in THorder array - cI
    for (int k = 0; k < n_cI; k++) {
    int x = p_cI[n_cI*ii+k];
    THorder[pos_cI+2*k] = p_cI[n_cI*ii+k]; THorder[pos_cI+2*k+1] = p_cI[n_cI*ii+k];}
    // fill in threshold values in THorder array - cro
    for (int k = 0; k < n_cro; k++) {
     int y = p_cro[n_cro*jj+k];
    THorder[pos_cro+2*k] = p_cro[n_cro*jj+k]; THorder[pos_cro+2*k+1] = p_cro[n_cro*jj+k];}
    // fill in threshold values in THorder array - cII (full model only)
    for (int k = 0; k < n_cII; k++) {THorder[pos_cII+2*k] = p_cII[n_cII*kk+k]; THorder[pos_cII+2*k+1] = p_cII[n_cII*kk+k];}
    // if (save_numb != 418) continue; // use this to process a single specific ordering
    // if (save_numb < 970) continue; // use this to run in smaller batches
    // if (save_numb > 1) break;

    if (save_numb < MinMOD) continue; // only models with larger permutation numbers
    if (MaxMOD && save_numb > MaxMOD) break; // if max is defined, stop processing exceeding that trheshold
    // specifically for loop tests some artificial constraints on cycle variables to select only a subset of functions
    if (G->HSMName == "Loops_2" && ii > jj) continue; // for symmetric 2 loop models exclude models with simply swapped function pairs
// use ResetThresholds for real models and ResetFuncDefs for loop tests
// do not do anything for other models
    if (G->HSMName == MyString("Loops_2") || G->HSMName == MyString("Circadian_03")) G->ResetFuncDefs2(ii,jj);
    if (G->HSMName == MyString("Lambda_Core_red") || G->HSMName == MyString("Lambda_Core_blue") || G->HSMName == MyString("Lambda_Complete") ||
        G->HSMName == MyString("Lambda_Oppenheim") || G->HSMName == MyString("HK22_Complete") || G->HSMName == MyString("Mu_3_01"))
        G->ResetThresholds(THorder);


// write out model state space data
//Memo1->Lines->Add("INICOUNT:\t\t" + MyString(G->IniCount));
cout << "Number of modes:\t\t" << G->Nmodes << "\n";
cout << "Number of active modes:\t" << G->NmodesA << "\n";
cout << "Number of consistent modes:\t" << G->NmodesC << "\n";
cout << "Number of genes:\t\t" << G->Ngenes << "\n";
cout << "Total number of transitions:\t" << G->Ntrans << "\n";
cout << "Maximal # of out transitions:\t" << G->Max_tout << "\n";
cout << "Maximal # of in transitions:\t" << G->Max_tin << "\n";
// create 2 matrix representations - G_orig - the original digraph; G_trans - the transposed digraph
// delete arrays if already exist
if (G_orig) delete[] G_orig; if (G_trans) delete[] G_trans; if (GN_orig) delete[] GN_orig; if (GN_trans) delete[] GN_trans;
G_orig = new int[G->Nmodes*G->Max_tout];
for (int i = 0; i < G->Nmodes*G->Max_tout; i++) G_orig[i] = -1; // initialise to -1
G_trans = new int[G->Nmodes*G->Max_tin];
for (int i = 0; i < G->Nmodes*G->Max_tin; i++) G_trans[i] = -1; // initialise to -1
GN_orig = new int[G->Nmodes]; GN_trans = new int[G->Nmodes];
for (int i = 0; i < G->Nmodes; i++)
    {
    GN_orig[i] = -1; GN_trans[i] = -1; // set number of edges to -1 to identify non-existent and inconsistent modes
    Mode *M = (Mode*) G->Modes->Items[i];
    if (M == 0) continue;
    // here pass to SCC analysis all consistent ones (isCionsistent will be set > 1 only there)
    if (M->isConsistent  == 0) continue; // ignore the ones not consistent with threshold oredering
    GN_orig[i] = 0; GN_trans[i] = 0; // set number of edges to 0 and count afterwards
    for (int j = 0; j < M->Ntrout; j++)
        {
        Mode *M1 = (Mode*) G->Modes->Items[(int) M->Tout_m[j]];
        if (M1->isConsistent)
            {
            G_orig[i*G->Max_tout+GN_orig[i]] = (int)  M->Tout_m[j];
            GN_orig[i]++;
            }
        //else Application->MessageBox(MyString(j).c_str(),MyString(G->IniCount).c_str(),0);
        }
    for (int j = 0; j < M->Ntrin; j++)
        {
        Mode *M1 = (Mode*) G->Modes->Items[(int) M->Tin_m->Items[j]];
        if (M1->isConsistent)
            {
            G_trans[i*G->Max_tin+GN_trans[i]] = (int) M->Tin_m->Items[j];
            GN_trans[i]++;
            }
        }
    }
// create SCC from G_orig and G_trans pair
// reinitialise arrays ant TLists
if (CompList) for (int i = 0; i < CNumb; i++) delete[] (int*) CompList->Items[i];
CompList->Clear(); CLsizes->Clear(); CLtypes->Clear();
for (int i = 0; i < G->Nmodes; i++) {Cidx[i] = -1; Cidx2[i] = -1; Csize[i] = -1;}
// compute SCCs of state space
CNumb = graphSCC(G->Nmodes,G->Max_tout,G->Max_tin,G_orig,G_trans,GN_orig,GN_trans,Cidx,Cidx2,Csize,CompList,CLsizes,CLtypes);
// *** checking for SCC temoprality below **************************************
// strict checking -- see whether there is any +/- gene that eventually will reach the threshold
for (int i = 0; i < CNumb; i++)
    {
    int temporal_c = 0; // assume non-temporal for now
    for (int j = 0; j < G->Ngenes; j++)
        {
        int *v_list = (int*) CompList->Items[i];
        // check for gene j whether it is strictly +/- for all SCC elements
        int count_plus = 0; int count_minus = 0;
        int is_plus = 0; int is_minus = 0;
        for (int k = 0; k < (int) CLsizes->Items[i]; k++)
            {
            Mode *M3 = (Mode*) G->Modes->Items[v_list[k]];
            if (M3->Genes[j] == 0) count_minus++;
            if (M3->Genes[j] == 1) count_plus++;
            }
        if (count_minus == (int) CLsizes->Items[i]) is_minus = 1;
        if (count_plus == (int) CLsizes->Items[i]) is_plus = 1;
        // if all decreasing check for occupied site in an arbitrary mode (for SCC then should be occupied also in all others)
        // if all increasing check for non-occupied site in an arbitrary mode (for SCC then should be non-occupied also in all others)
        Mode *M4 = (Mode*) G->Modes->Items[Cidx[v_list[0]]]; // use the 1st component mode as a representative
        for (int k = 0; k < G->Nb2+2*G->Nb3; k++)
            {
            // for decreasing simple - check thus for occuppied site
            if (G->BStfs[k] == j && is_minus && M4->BS[k] == 1) temporal_c = 1;
            // for increasing additional check for ternary - sibling site must also be unoccupied
            if (G->BStfs[k] == j && is_plus && M4->BS[k] == 0)
                {
                if (k < G->Nb2) temporal_c = 1; // just binary site
                if (k >= G->Nb2 && ((k-G->Nb2) % 2 == 0) && M4->BS[k+1] == 0) temporal_c = 1; // ternary, sibling next
                if (k >= G->Nb2 && ((k-G->Nb2) % 2 == 1) && M4->BS[k-1] == 0) temporal_c = 1; // ternary, sibling previous
                }
            }
        }
    // SHOUT dialog if non-temporal component with outgoing edges is found
    // also set component type value to 2
    if (temporal_c == 0 && (int*) CLtypes->Items[i] == 0)
        {
        //Application->MessageBox(MyString(G->IniCount).c_str(),"Non-temporal SCC with outgoing edges!!!",0);
        // Application->MessageBox(MyString(i).c_str(),"Non-temporal SCC with outgoing edges!!!",0);
        Flag = 2; // set flag value to 1
        CLtypes->Items[i] = (void*) 2;
        }
    // set component ty
    }
// *** checking for SCC temoprality above **************************************
// find the modes reachable from 0, mark these with consistency 3
for (int i = 0; i < G->Nmodes; i++) if (G->Modes->Items[i]) if (Cidx2[i] == 0)
    {
    Mode *M2 = (Mode*) G->Modes->Items[i];
    M2->isConsistent = 3;
    G->NmodesZero++;
    G->NmodesIni++;
    }

// then check for attractors outside 0-reachibility area
int non_zero_attr = 0;
for (int i = 0; i < CNumb; i++)
    {
    Mode *M3 = (Mode*) G->Modes->Items[((int*)CompList->Items[i])[0]]; // chose the 1-st component vertex as representative
    if (M3->isConsistent == 1 && (int*) CLtypes->Items[i] == (int*) 1) non_zero_attr = 1; // attractor not reachable from 0 found
    }
Flag = Flag+non_zero_attr; // thus = 1 - temporal without out edges, = 2 - "proper" non-zero attractor, = 3 both of these separately
// write out component data
cout << "SCCs constructed!" << "\n";
cout << "Number of proper components:\t\t" << CNumb << "\n";
cout << "Number of modes reachable from 0:\t" << G->NmodesZero << "\n";
// write out model state space in file
MyString SaveNumber  = FormatSaveNumber(save_numb);
MyString SaveNumberFlag  = FormatSaveNumber(save_numb);
if (Flag) SaveNumberFlag = SaveNumberFlag +MyString( "_") + MyString(Flag); // append Flag to name, if non-zero
MyString ModelStateSpaceName = G->HSMName + MyString("/model_") +  SaveNumber + MyString(".txt");



// for function test use function indices instead of save number
if (G->HSMName == "Circadian_03" || G->HSMName == "Circadian_03_2")
    {
    ModelStateSpaceName = G->HSMName + MyString("/model_") +  MyString(ii) + MyString("_") + MyString(jj) + MyString(".txt");
    }
int write_res = write_HSM(G,ModelStateSpaceDirectory+ModelStateSpaceName,CompList,Cidx,Cidx2,CLsizes,CLtypes,ConsThreshold);
if (write_res)
    {
    cout << "State space written in file:\t\t" << ModelStateSpaceName.c_str() << "\n";
    cout << "Used consistency threshold:\t\t" << ConsThreshold << "\n";
    cout  << "\n";
    }
else {cout << "Failed to write model state space file!!!"; return;}

// *** drawing removed here from non-GUI version *******************************
// *** drawing removed here from non-GUI version *******************************

cout <<"\n";

// a batch of closing "}" to match # of order checking if-s above :)
}}}}

cout <<"\n\n\n";
cout << "All done!!!";
}

// *****************************************************************************
// additional helper functions (not auto-generated) are below ******************

// format file number to appropropriete length (1-6 chars)
MyString FormatSaveNumber(int n)
{
MyString ret_str("000000");
if (n < 10) sprintf(ret_str.sbuf+5, "%d", n);
if (n >= 10 && n < 100) sprintf(ret_str.sbuf+4, "%d", n);
if (n >= 100 && n < 1000) sprintf(ret_str.sbuf+3, "%d", n);
if (n >= 1000 && n < 10000) sprintf(ret_str.sbuf+2, "%d", n);
if (n >= 10000 && n < 100000) sprintf(ret_str.sbuf+1, "%d", n);
if (n >= 100000) sprintf(ret_str.sbuf+0, "%d", n);
return ret_str;
}

// generate all pn permutations outA of n integers from array inA, return number of them (should coinside with pn)
// inBuf is extra memory for copying inA values for subsequent calls
int AllPermutations(int n_max,int n,int pn,int out_pos,int *inA,int *outA)
{
if (n < 1) return 0; // smth wrong
if (n == 1) {outA[out_pos] = inA[0]; return 1;}
// part below assumes n > 1
int p_count = 0;
int pn2 = pn/n;
int p_count2 = 0;
// first make n copies in outA of the called permutations at positions i*pn2
for (int i = 0; i < n; i++)
    {
    // copy the whole array
    for (int j = 0 ; j < n; j++) outA[out_pos+n_max*(i*pn2)+j] = inA[j];
    // then swap (n-1-i)-th and (n-1)-th elements
    outA[out_pos+n_max*(i*pn2)+(n-1-i)] = inA[n-1]; outA[out_pos+n_max*(i*pn2)+(n-1)] = inA[n-1-i];
    // fill with outA[out_pos+n_max*(i*pn2)+(n-1)] elements n blocks of length pn2 in positions n_max*(i*pn2)
    for (int k = 0 ; k < pn2; k++)  outA[out_pos+n_max*(i*pn2+k)+(n-1)] = outA[out_pos+n_max*(i*pn2)+(n-1)];
    }        // return 0;
// now make n successive recursive calls
for (int i = 0; i < n; i++)
    {
    // create inA from the value saved in outA
    for (int j = 0 ; j < n; j++) inA[j] = outA[out_pos+n_max*(i*pn2)+j];
    // call recursively for n-1, put the result starting at n_max*(i*pn2)-th position
    p_count2 = AllPermutations(n_max,n-1,pn2,out_pos+n_max*(i*pn2),inA,outA);
    p_count = n*p_count2;
    }
return(p_count);
}

// a quick and dirty as for now :)
int findidx(int *pos, int size,int val)
{
for (int i = 0 ; i < size; i++) if (pos[i] == val) return(i);
return -1;
}




