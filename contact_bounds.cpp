/*
 * contact_bounds.cpp
 *
 *  Created on: 08-05-2015
 *      Author: rdyja
 */

#include "contact_bounds.hpp"

using namespace TALYFEMLIB;

void ContactBounds::LoadContactBounds(GRID *pg,
		  const std::vector<PetscInt>& oldID, const std::vector<PetscInt>& newID) {

	  p_grid_ = pg;
	  is_periodic_ = true;

	  // set the partners for the boundaries
	  is_node_periodic_.redim(this->p_grid_->n_nodes());
	  is_node_periodic_.fill(false);
	  if (oldID.size() != newID.size()) {
		  throw TALYException() << "Unpaired number of periodic nodes";
	  }
	  for (size_t i = 0; i < oldID.size(); i++) {
//	    PetscInt newID = NewNodeID(oldID);
	    if (p_grid_->parallel_type_ == kWithDomainDecomp) {
	    	throw TALYException() << "Domain Decomposition unsupported in LoadPeriodicBounds";
	    } else {
			is_node_periodic_.set(oldID[i], true);
			pbc_partners_[oldID[i]] = newID[i];
			pbc_sol_partners_[oldID[i]] = newID[i];
	    }
	  }

}


