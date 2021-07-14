//---------------------------------------------------------------------------

#ifndef graph_SCCH
#define graph_SCCH

#include "my_bcb.h"
//---------------------------------------------------------------------------
int graphSCC(int n,int m_out,int m_in,int *G_orig,int *G_trans,int *GN_orig,int *GN_trans,int *Cidx,int *Cidx2,int *Csize,MyListL *CompList,MyListL *CLsizes,MyListL *CLtypes);

// DFS functions below
void Visit(int u,int root,int *GN,int *G);
void VisitT(int u,int root,int *GN,int *G);

// check whether the vertex is non-zero intial (with sin gle "1" - i.e. a power of 2)
int CheckVertIni(int u);
#endif
