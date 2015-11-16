/*
 * contact_bounds.cpp
 *
 *  Created on: 08-05-2015
 *      Author: rdyja
 */

#include "contact_bounds.hpp"

using namespace TALYFEMLIB;
using namespace std;

void ContactBounds::LoadContactBounds(GRID *pg,
		  const std::vector<PetscInt>& oldID, const std::vector<PetscInt>& newID) {

	  p_grid_ = pg;
	  has_periodic_ = true;

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

void ContactBounds::ImportContactBounds(GRID *pg) {

	  if (pg->get_node_periodics().size() > 0) {
		  has_periodic_ = true;
//		  PrintInfo("TransferContactBounds");
	  } else {
		  has_periodic_ = false;
		  return;
	  }

	  p_grid_ = pg;

	  const map<int,int>& grid_periodics = pg->get_node_periodics();

	  // set the partners for the boundaries
	  is_node_periodic_.redim(p_grid_->n_nodes());
	  is_node_periodic_.fill(false);

//	  for (int i = 0; i < pg->n_nodes(); i++) {
//		  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
//		  				  "%d) physical_map=%d  solution_map=%d\n",
//						  i, pg->physical_map(i), pg->solution_map(i));
//	  }
//	  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

	  for (map<int,int>::const_iterator it = grid_periodics.begin();
			  it != grid_periodics.end(); ++it) {
		  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				  "grid_periodics: first=%d  second=%d\n", it->first, it->second);
	  }
	  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

	  for (map<int,int>::const_iterator it = grid_periodics.begin();
			  it != grid_periodics.end(); ++it) {

		if (p_grid_->parallel_type_ == kWithDomainDecomp) {
			pbc_partners_[it->first] = it->second;
			pbc_sol_partners_[it->first] = it->second;
			// is_node_periodic must be handled in special way
		} else {
			pbc_partners_[it->first] = it->second;
			pbc_sol_partners_[it->first] = it->second;
			is_node_periodic_.set(it->first, true);
		}

	  }

	  if (p_grid_->parallel_type_ == kWithDomainDecomp) {
		  for (int i = 0; i < pg->n_nodes(); i++) {
			  if (grid_periodics.find(p_grid_->solution_map(i)) != grid_periodics.end())
				  is_node_periodic_.set(i, true);
		  }
	  }

}

