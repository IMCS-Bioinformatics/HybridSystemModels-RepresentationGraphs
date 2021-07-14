//---------------------------------------------------------------------------

#ifndef HSM_ioH
#define HSM_ioH
//---------------------------------------------------------------------------
#endif

#include "HSM_graph.h"
#include "my_bcb.h"

// write state space of HSM defined by G in file Fname, SCC information is provided in TLists, Cidx contains identifying id (DFS root) for vertices component
// write only states at least of consistency threshold ConsThreshold
int write_HSM(HSMGraph *G,MyString Fname,MyListL *CompList,int *Cidx,int *Cidx2,MyListL *CLsizes,MyListL *CLtypes,int ConsThreshold);
// read HSM model from file Fname
int read_HSM(MyString Fname,MyString &HSMName,int &Ngenes,MyString* &gNames,int &Ntf,int &Nb2,MyString* &b2Names,int &Nb3,MyString* &b3Names,int* &BStfs,int* &Ngf,int* &Paramsgf,int* &Defsgf);
