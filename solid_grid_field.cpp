#include <Grid/node.h>
#include "solid_grid_field.hpp"
#include "solid_input_data.hpp"

using namespace TALYFEMLIB;

void SolidGridField::SetIC(int nsd) {
	/*
    for(int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
        SolidNodeData* pData = &(Node(nodeID));
        if (nsd == 3) {
            int elem_id = p_grid_->node_array_[nodeID]->elem_id_;
            ELEM *elem = p_grid_->elm_array_[elem_id];
            int mat_ind = elem->mat_ind();
            pData->set_curr_temp(inputData_.initial_temperature(mat_ind));
            PrintInfo("SolidGridField::elem_id: ", elem_id);
        }
    }
    PrintStatusStream(std::cerr, "IC set ");
    */
    for(int elemID = 0; elemID < p_grid_->n_elements(); elemID++) {
//        SolidNodeData* pData = &(Node(nodeID));
//            int elem_id = p_grid_->node_array_[nodeID]->elem_id_;
		ELEM *elem = p_grid_->elm_array_[elemID];
		int mat_ind = elem->mat_ind();

		for(int i = 0; i < elem->n_nodes(); i++) {
			SolidNodeData* pData = &(GetNodeData(elem->node_id_array(i)));
			pData->set_curr_temp(inputData_.initial_temperature(mat_ind));
//			PrintInfo("SolidGridField::elem_id: ", elemID);
		}
	}


}


// temporary function to set proper boundary and surface indicators
// needed for contact boundary condition
void SolidGridField::SetBndrIndicator(const ContactBounds& pcb)
{
/*
	PrintInfo("SolidGridField::SetBndrIndicator");
	for (int i = 0; i < p_grid_->n_nodes() ; i++) {
		if (pcb.GetPeriodicPartner(i) != -1) {
			p_grid_->GetNode(i)->addIndicatorNum(7); //we need to decide which indicators are for the 4th type BC
		}
	}

    p_grid_->SetCaredSurfaceIndicator();
    p_grid_->GenElmSurfaceIndicator();
*/
    //pGrid->CalcHalfBandWidth();
}

