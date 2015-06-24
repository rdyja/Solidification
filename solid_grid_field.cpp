#include <Grid/node.h>
#include "solid_grid_field.hpp"
#include "solid_input_data.hpp"

using namespace TALYFEMLIB;

void SolidGridField::SetIC(int nsd) {

	for(int elemID = 0; elemID < p_grid_->n_elements(); elemID++) {
		ELEM *elem = p_grid_->elm_array_[elemID];
		int mat_ind = elem->mat_ind();

		for(int i = 0; i < elem->n_nodes(); i++) {
			SolidNodeData* pData = &(Node(elem->node_id_array(i)));
			pData->set_curr_temp(inputData_.initial_temperature(mat_ind));
		}
	}

}


