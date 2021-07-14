//---------------------------------------------------------------------------

#ifndef HSM_componentsH
#define HSM_componentsH
//---------------------------------------------------------------------------

#include "HSM_graph.h"

class HSMComponent;

int buildComponents(HSMGraph *G,int *Cidx,HSMComponent ***CList);

class HSMComponent
{
public:
    HSMComponent(int iCsize);
    ~HSMComponent() {delete[] Midx;}
    int addMidx(int idxM);
    int Cnumb;
    char isActive;
    int *Midx; // array of indexes to component modes
private:
    int ccount;
};

#endif
