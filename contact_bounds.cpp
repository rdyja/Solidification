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

	  p_grid_ = pg;

//	  PrintInfo("TransferContactBounds");


	  const map<int,int>& grid_periodics = pg->get_node_periodics();

	  // set the partners for the boundaries
	  is_node_periodic_.redim(this->p_grid_->n_nodes());
	  is_node_periodic_.fill(false);

	  for (map<int,int>::const_iterator it = grid_periodics.begin();
			  it != grid_periodics.end(); ++it) {
		  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				  "grid_periodics: first=%d  second=%d\n", it->first, it->second);
		  if (p_grid_->parallel_type_ == kWithDomainDecomp) {
			  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
					  "grid_periodics: solution_map(first)=%d solution_map(second)=%d\n",
					  p_grid_->solution_map(it->first), p_grid_->solution_map(it->second));
		  }
		  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	  }

	  for (map<int,int>::const_iterator it = grid_periodics.begin();
			  it != grid_periodics.end(); ++it) {

		if (p_grid_->parallel_type_ == kWithDomainDecomp) {
//			throw TALYException() << "Domain Decomposition unsupported in LoadPeriodicBounds";
			is_node_periodic_.set(it->first, true);
			pbc_partners_[it->first] = it->second;
			pbc_sol_partners_[it->first] = p_grid_->solution_map(it->second);
		} else {
			is_node_periodic_.set(it->first, true);
			pbc_partners_[it->first] = it->second;
			pbc_sol_partners_[it->first] = it->second;
		}

	  }

	  if (pg->get_node_periodics().size() > 0) {
		  has_periodic_ = true;
	  } else {
		  has_periodic_ = false;
	  }

}

