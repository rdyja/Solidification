#include <Grid/Grid.h>
#include "solid_grid_field.hpp"
#include "solid_input_data.hpp"

using namespace TALYFEMLIB;

void SolidGridField::SetIC(int nsd) {
    for(int nodeID=1; nodeID<=pGrid->nodeno; nodeID++) {
        SolidNodeData* pData = &(Node(nodeID));
        if (nsd == 3) {
            pData->set_curr_temp(inputData_.initial_temperature());
        }
    }
    PrintStatusStream(std::cerr, "IC set ");
}

void SolidGridField::SetBndrIndicator(double width, double length, double
height)
{
    for (int elmID = 1; elmID < pGrid->elmno; elmID++) {
        ELEM* e = pGrid->GetElm(elmID);
        for (int i = 0; i < e->getSurfaceNo(); i++) {
            if (e->pElemIDArray[i] == 0) {
                for (int n = 0; n < e->nodeno; n++) {
                    if (i == n) continue;
                    int inds[] = { 1 };
                    int nInd = 1;
                    int nodeID = e->pNodeIDArray[n] + 1;
                    pGrid->GetNode(nodeID)->setIndicator(nInd, inds);
                }
            }
        }
    }

    pGrid->setCaredSurfaceIndicator();
    pGrid->genElmSurfaceIndicator();
    pGrid->CalcHalfBandWidth();
}

