//---------------------------------------------------------------------------


#pragma hdrstop

#include "graph_SCC.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

//#include <stdio.h>
//using namespace std;

int max_out;
int max_in;
int *visited;
int *assignment;
int *f_list; // correspondence of vertices to finalizing times (from 0 to # of active vertices - 1)
int ftime; // "finalizing" time

int graphSCC(int n,int m_out,int m_in,int *G_orig,int *G_trans,int *GN_orig,int *GN_trans,int *Cidx,int *Cidx2,int *Csize,MyListL *CompList,MyListL *CLsizes,MyListL *CLtypes)
{
max_out = m_out; max_in = m_in;
int cn = 0; // number of componnets
int cn2 = 0; // number of components with at least 2 vertices
visited = new int[n]; assignment = new int[n]; f_list = new int[n]; // create temporary arrays
ftime = 0; // "finalizing" time
// initialise arrays
for (int i = 0; i < n; i++) {visited[i] = -1; assignment[i] = -1; f_list[i] = -1;}
// do 2xDFS
// direct DFS
for (int i = 0; i < n; i++) if (GN_orig[i] != -1) Visit(i,i,GN_orig,G_orig); // consider only active  ones
// transposed DFS
for (int i = ftime-1; i >= 0; i--)
    {
    if (assignment[f_list[i]] == -1) cn++;  // new component found
	VisitT(f_list[i],f_list[i],GN_trans,G_trans);
	}
// count the number of components
int *csizes = new int[n]; // size of component to which given vertex belongs
int *cidx2 = new int[n]; // the position of component in CompList
for (int i = 0; i < n; i++) {csizes[i] = 0; cidx2[i] = -1;}
for (int i = 0; i < n; i++) if (GN_orig[i] != -1) csizes[assignment[i]]++;
// count components of size 2 or more
for (int i = 0; i < n; i++) if (csizes[i] > 1) cn2++;
int *cpos = new int[cn2]; // next free index in component members array
for (int i = 0; i < cn2; i++) cpos[i] = 0;
// int a1 = cpos[0];
// fill component assignment and size arrays
for (int i = 0; i < n; i++) if (GN_orig[i] != -1)
    {
    Cidx[i] = assignment[i];
    Cidx2[i] = visited[i];
    Csize[i] = csizes[assignment[i]];
    }
// fill the CompList - lists of integer arrays containing members of each of the components (currently consider only components with sizes > 1)
// the first iteration - create arrays of component vertices and add them to CompList
// choose: 1) active, 2) single representative (with root id), with size > 1
int c2_count = 0;
for (int i = 0; i < n; i++) if (GN_orig[i] != -1) if (Cidx[i] == i) if (Csize[i] > 1)
    {
    int* comp_list = new int[csizes[i]]; // create new integer list of vertices belonging to the component
    CompList->Add((void*) comp_list); // add it to TList
    CLsizes->Add((void*) Csize[i]); // set component size
    CLtypes->Add((void*) 1); // set component type - currently set all types to 1 (attractors)
    cidx2[i] = c2_count;
    c2_count++;
    }
// the second iteration - add vertices to component lists
for (int i = 0; i < n; i++) if (GN_orig[i] != -1) if (Csize[i] > 1)
    {
    int* comp_list2 = (int*) CompList->Items[cidx2[Cidx[i]]];
    comp_list2[cpos[cidx2[Cidx[i]]]] = i;
    cpos[cidx2[Cidx[i]]]++;
    }
// *** checking for SCC temoprality below **************************************
// this is ersatz version for defining non-active components - include as nonactive all SCCs with outgoing edges
// strict checking in "main" file (for this access to G is needed)
for (int i = 0; i < cn2; i++) for (int j = 0; j < (int) CLsizes->Items[i]; j++)
    {
     /*
    if (i == 9)
        {
        int aaa = 12;
        }
      */
    int *v_list = (int*) CompList->Items[i];
    int v = v_list[j];
    for (int k = 0; k < GN_orig[v]; k++)  // check for outgoing destinations
            {
            // find transition destination
            int dest = G_orig[v*max_out+k];
            bool is_out = true; // destination leaves component
            for (int l = 0; l < (int) CLsizes->Items[i]; l++) if (dest == ((int*)CompList->Items[i])[l]) is_out = false;
            if (is_out) CLtypes->Items[i] = (int*) 0;
            }
    }
// *** checking for SCC temoprality above **************************************
// free memory and return
delete[] visited; delete[] assignment; delete[] f_list; delete[] csizes; delete[] cidx2; delete[] cpos;
return cn2;
}

void Visit(int u,int root,int *GN,int *G)
{
if (visited[u] != -1) return;
visited[u] = root;
for (int i = 0; i < GN[u]; i++) Visit(G[u*max_out+i],root,GN,G);
f_list[ftime] = u; ftime++;
}

void VisitT(int u,int root,int *GN,int *G)
{
if (assignment[u] != -1) return;
assignment[u] = root;
for (int i = 0; i < GN[u]; i++) VisitT(G[u*max_in+i],root,GN,G);
}


// check whether the vertex is non-zero intial (with sin gle "1" - i.e. a power of 2)
int CheckVertIni(int u)
{
if (u == 1) return 1; // just 1
if (u % 2 != 0) return 0; // not 1 and odd
return CheckVertIni(u/2); // divide by 2 and continue
}

                   