/*
 * contact_bounds.hpp
 *
 *  Created on: 07-05-2015
 *      Author: rdyja
 */

#ifndef CONTACT_BOUNDS_HPP_
#define CONTACT_BOUNDS_HPP_

#include <petsc.h>  // for PETSc objects and types

#include <map>  // for std::map
#include <vector>  // for std::vector

#include "Common/TALYException.h"  // for throwing exception
#include "DataStructure/zeroarray.h" //for ZEROARRAY
#include "Grid/grid_types/grid.h"  // for GRID class
#include "FEM/periodic_bounds.h" // for PeriodicMap

namespace TALYFEMLIB {

/*
 * ContactBounds is similar to PeriodicBounds class. The main difference
 * is not calculating periodic partners, but loading them from file. The other
 * difference is that it doesn't store data about periodic boundaries, only
 * about periodic nodes.
*/
class ContactBounds {
public:
	// Default constructor
	ContactBounds(): p_grid_(nullptr),  is_periodic_(false) { }
	virtual ~ContactBounds() { }

	inline bool is_periodic() { return is_periodic_; }
	inline bool IsNodePeriodic(int index) {
		return is_node_periodic_.get(index);
	}
	inline int GetPeriodicPartner(int index) const {
		PeriodicPhysicalMap::const_iterator item = pbc_partners_.find(index);
		return item != pbc_partners_.end() ? item->second : -1;
	}
	inline int GetPeriodicSolPartner(int index) const {
		PeriodicSolutionMap::const_iterator item = pbc_sol_partners_.find(index);
		return item != pbc_sol_partners_.end() ? item->second : -1;
	}
	inline const PeriodicPhysicalMap* pbc_partners() const {
		return &pbc_partners_;
	}
	inline const PeriodicSolutionMap* pbc_sol_partners() const {
		return &pbc_sol_partners_;
	}

	inline GRID* Grid() { return p_grid_; }


	void LoadContactBounds(GRID *pg,
			const std::vector<PetscInt>& oldID, const std::vector<PetscInt>& newID);
	void LoadContactBounds(GRID *pg);

private:
	GRID* p_grid_;
	bool is_periodic_;
	ZEROARRAY<bool> is_node_periodic_;
	PeriodicPhysicalMap pbc_partners_;
	PeriodicSolutionMap pbc_sol_partners_;
};

}

#endif /* CONTACT_BOUNDS_HPP_ */
