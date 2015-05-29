#include <Grid/node.h>
#include "solid_grid_field.hpp"
#include "solid_input_data.hpp"

using namespace TALYFEMLIB;

void SolidGridField::SetIC(int nsd) {
    for(int nodeID=0; nodeID<pGrid->n_nodes(); nodeID++) {
        SolidNodeData* pData = &(Node(nodeID));
        if (nsd == 3) {
            pData->set_curr_temp(inputData_.initial_temperature());
        }
    }
    PrintStatusStream(std::cerr, "IC set ");
}


// temporary function to set proper boundary and surface indicators
// needed for contact boundary condition
void SolidGridField::SetBndrIndicator(const ContactBounds& pcb)
{
	PrintInfo("SolidGridField::SetBndrIndicator");
	for (int i = 0; i < pGrid->n_nodes() ; i++) {
		if (pcb.GetPeriodicPartner(i) != -1) {
			pGrid->GetNode(i)->addIndicatorNum(7); //we need to decide which indicators are for the 4th type BC
		}
	}

    pGrid->SetCaredSurfaceIndicator();
    pGrid->GenElmSurfaceIndicator();

    //pGrid->CalcHalfBandWidth();
}

