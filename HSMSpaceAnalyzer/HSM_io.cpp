//---------------------------------------------------------------------------

#include <fstream>
#include <iostream>
using namespace std;

#pragma hdrstop

#include "HSM_io.h"
#include "HSM_graph.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

// write state space of HSM defined by G in file Fname
int write_HSM(HSMGraph *G,MyString Fname,MyListL *CompList,int *Cidx,int *Cidx2,MyListL *CLsizes,MyListL *CLtypes,int ConsThreshold)
{
if (G == 0) return 0; // nothing to write
ofstream outf;
outf.open(Fname.c_str(),ios::trunc);
if (!outf.is_open()) {cout << "Failed to open file for writing:\t\t" << Fname.c_str() << "\n"; return 0;}

// *** ConsThreshold not used just now, should be properly added later ***

// write out the threshold ordering array "as is" (should be later relocated in more appropriate place)
outf << G->Ntf*G->THmax << "\n";
for (int i = 0; i < G->Ntf*G->THmax; i++) outf << (int) G->THorder[i] << " ";
outf << "\n\n";

// general HSM model description
outf << G->Ngenes << "\n";
for (int i = 0; i < G->Ngenes; i++) outf << G->gNames[i].c_str() << " ";
outf << "\n";
outf << G->Nb2 << " " << G->Nb3 << "\n";
// BFs by names
for (int i = 0; i < G->Nb2; i++) outf << G->b2Names[i].c_str() << " ";
for (int i = 0; i < G->Nb3; i++) outf << G->b3Names[i].c_str() << " ";
outf << "\n";
// BFs by gene numbers
for (int i = 0; i < G->Nb2; i++) outf << G->BStfs[i] << " ";
for (int i = 0; i < 2*G->Nb3; i++) outf << G->BStfs[G->Nb2+i] << " ";
// NB! This is newly addded!
// Arguments and defintions of gene regulatory functions
outf << "\n";
outf << "\n";
for (int i = 0; i < G->Ngenes; i++) outf << G->Ngf[i] << " ";  // # of arguments
outf << "\n";
outf << "\n";
// argument defintions
for (int i = 0; i < G->Ngenes; i++)
    {
    for (int j = 0; j < G->Ngf[i]; j++) outf << G->Paramsgf[10*i+j] << " ";
    outf << "\n";
    }
outf << "\n";
// function defintions
for (int i = 0; i < G->Ngenes; i++)
    {
    int fargN = exp2(G->Ngf[i]);
    for (int j = 0; j < fargN; j++) outf << G->Defsgf[1024*i+j] << " ";
    outf << "\n";
    }
outf << "\n";
outf << "\n";
// description of HSM state space graph
outf << G->NmodesZero << " " << G->NmodesC << " " << G->Nmodes  << "\n";
outf << CompList->Count << "\n";
for (int i = 0; i < CompList->Count; i++)
    {
    int *v_list = (int*) CompList->Items[i];
    outf << Cidx[v_list[0]] << " " << (int) CLsizes->Items[i] << " " << (int) CLtypes->Items[i] << " ";
    }
outf << "\n";
for (int i = 0; i < CompList->Count; i++)
    {
    int *v_list = (int*) CompList->Items[i];
    for (int j = 0; j < (int) CLsizes->Items[i]; j++) outf << v_list[j] << " ";
    outf << "\n";
    }
outf << "\n";
// what follows are vN exactly blocks describing each of the vertices
for (int i = 0; i < G->Nmodes; i++) if (G->Modes->Items[i] != 0) // the mode is active
    {
    Mode *M = (Mode*) G->Modes->Items[i];
    if (M->isConsistent == 0) continue; // this mode is inconsistent with ordering
    outf << M->Id << "\n";
    outf << Cidx[M->Id] << " " << Cidx2[M->Id] << "\n";
    for (int j = 0; j < G->Ngenes; j++) outf << (int) M->Genes[j] << " ";
    outf << "\n";
    for (int j = 0; j < G->Nb2+2*G->Nb3; j++) outf << (int) M->BS[j] << " ";
    outf << "\n";
    // double work - check for number of edges to consistent vertices first, and print out that number of vertices afterwards
    int count_out = 0;
    for (int j = 0; j < M->Ntrout; j++)
        {
        Mode *M1 = (Mode*) G->Modes->Items[(int) M->Tout_m[j]];
        if (M1->isConsistent) count_out++;
        }
    outf << count_out << "\n";
    for (int j = 0; j < M->Ntrout; j++)
        {
        Mode *M1 = (Mode*) G->Modes->Items[(int) M->Tout_m[j]];
        if (M1->isConsistent) outf << (int) M->Tout_m[j] << " " << (int) M->Tout_g[j] << " ";
        }
    outf << "\n";
    outf << "\n";
    }
outf.close();
return 1; // writing completed
}

// read HSM model from file Fname
int read_HSM(MyString Fname,MyString &HSMName,int &Ngenes,MyString* &gNames,int &Ntf,int &Nb2,MyString* &b2Names,int &Nb3,MyString* &b3Names,int* &BStfs,int* &Ngf,int* &Paramsgf,int* &Defsgf)
{
ifstream inf;
inf.open(Fname.c_str());
if (!inf.is_open()) {cout << "Failed to open file for reading:\t\t" << Fname.c_str() << "\n"; return 0;}
char buf[1000]; // should be enough for everything :)
inf >> buf; HSMName = MyString(buf);
inf >> Ngenes;
gNames = new MyString[Ngenes];
for (int i = 0; i < Ngenes; i++) {inf >> buf; gNames[i] = MyString(buf);}
inf >> Ntf;
inf >> Nb2;
b2Names = new MyString[Nb2];
for (int i = 0; i < Nb2; i++) {inf >> buf; b2Names[i] = MyString(buf);}
inf >> Nb3;
b3Names = new MyString[Nb3];
for (int i = 0; i < Nb3; i++) {inf >> buf; b3Names[i] = MyString(buf);}
int bs_numb = Nb2+2*Nb3;
BStfs = new int[bs_numb];
for (int i = 0; i < bs_numb; i++) inf >> BStfs[i];
Ngf = new int[Ngenes];
for (int i = 0; i < Ngenes; i++) inf >> Ngf[i];
// create gene reg function parameter arrays
Paramsgf = new int[Ngenes*10]; // assume max 10 arguments
for (int i = 0; i < Ngenes*10; i++) Paramsgf[i] = -1; // initialise to -1 - just in case for clarity
// read parameter values themselves
for (int i = 0; i < Ngenes; i++) for (int j = 0; j < Ngf[i]; j++) inf >> Paramsgf[10*i+j];
// create reg function arrays
Defsgf = new int[Ngenes*1024]; // assume max 10 arguments
for (int i = 0; i < Ngenes; i++)
    {
    int argN = exp2(Ngf[i]);
    for (int j = 0; j < argN; j++) inf >> Defsgf[1024*i+j];
    }
inf.close();
return 1;
}
