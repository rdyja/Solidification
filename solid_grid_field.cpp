#include <Grid/node.h>
#include "solid_grid_field.hpp"
#include "solid_input_data.hpp"

using namespace TALYFEMLIB;

void SolidGridField::SetIC(int nsd) {
    for(int nodeID=0; nodeID <  p_grid_->n_nodes(); nodeID++) {
        SolidNodeData* pData = &(Node(nodeID));
        if (nsd == 3) {
            int elem_id = p_grid_->node_array_[nodeID]->elem_id_;
            ELEM *elem = p_grid_->elm_array_[elem_id];
            int mat_ind = elem->mat_ind();
            pData->set_curr_temp(inputData_.initial_temperature(mat_ind));
        }
    }
    PrintStatusStream(std::cerr, "IC set ");
}

void SolidGridField::SetBndrIndicator(double width, double length, double
height)
{
    /*for (int elmID = 0; elmID < pGrid->n_elements(); elmID++) {
        ELEM* e = pGrid->GetElm(elmID);
        for (int i = 0; i < e->getSurfaceNo(); i++) {
            if (e->pElemIDArray[i] == 0) {
                for (int n = 0; n < e->nodeno; n++) {
                    if (i == n) continue;
                    //int inds[] = { 1 };
                    //int nInd = 1;
                    int nodeID = e->pNodeIDArray[n];
                    ///???
                    pGrid->GetNode(nodeID)->addIndicatorNum(1);
                }
            }
        }
    }*/

    ///??? ???
    ///pGrid->setCaredSurfaceIndicator();
    ///pGrid->genElmSurfaceIndicator();    
    ///???
    //pGrid->CalcHalfBandWidth();
}

