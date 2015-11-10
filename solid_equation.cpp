#include <Grid/femelm.h>
#include <DataStructure/zeroarray.h>
#include <DataStructure/zeromatrix.h>

#include "solid_equation.hpp"

using namespace TALYFEMLIB;

SolidEquation::SolidEquation(SolidInputData* id, TALYFEMLIB::ContactBounds* cb)
 :idata(id), contact_bounds_(NULL) {
    timers_[kTimerSolve].set_label("Solve");
    timers_[kTimerAssemble].set_label("Assemble");
    timers_[kTimerKSPSolve].set_label("KSPSolve");
    timers_[kTimerUpdate].set_label("Update");
}

SolidEquation::~SolidEquation() {
  timers_[kTimerSolve].PrintGlobalAverageSeconds();
  timers_[kTimerSolve].PrintGlobalTotalSeconds();
  timers_[kTimerAssemble].PrintGlobalAverageSeconds();
  timers_[kTimerAssemble].PrintGlobalTotalSeconds();
  timers_[kTimerKSPSolve].PrintGlobalAverageSeconds();
  timers_[kTimerKSPSolve].PrintGlobalTotalSeconds();
  timers_[kTimerUpdate].PrintGlobalAverageSeconds();
  timers_[kTimerUpdate].PrintGlobalTotalSeconds();
}


void SolidEquation::Solve(double t, double dt) {
	timers_[kTimerSolve].Start();  // we want to time the entire solve process
    this->t_ = t;
    this->dt_ = dt;

    fillEssBC();
    ApplyEssBCToSolution();
    timers_[kTimerAssemble].Start();  // we're timing just the assembly
    Assemble();
    timers_[kTimerAssemble].Stop();
    ApplyEssBC();

    timers_[kTimerKSPSolve].Start();  // time the KSPSolve step
    SolveKSP(solution_, 1, 0);
    timers_[kTimerKSPSolve].Stop();

    timers_[kTimerUpdate].Start();  // time the process of saving the solution
    for(int A=0; A < p_grid_->n_nodes(); A++){
        double newval=solution_(A);
        //pData->GetNodeData(A).UpdateDataStructures();
        p_data_->GetNodeData(A).set_curr_temp(newval);
    }
    timers_[kTimerUpdate].Stop();
    timers_[kTimerSolve].Stop();

}

double SolidEquation::compute_average_velocity(const TALYFEMLIB::FEMElm& fe, int nbf2) {
	const int nbf = fe.pElm->n_nodes();

    double Vsr = 0.0;
    for(int i = 0; i < nbf; ++i) {
        int J = fe.pElm->ElemToLocalNodeID(i);
        Vsr += p_data_->GetNodeData(J).get_velocity();
    }
    return Vsr /= nbf;
}

double SolidEquation::compute_average_temp(const TALYFEMLIB::FEMElm& fe, int nbf2) {
	const int nbf = fe.pElm->n_nodes();

    double Tsr = 0.0;
    for(int i = 0; i < nbf; ++i) {
        int J = fe.pElm->ElemToLocalNodeID(i);
        Tsr += p_data_->GetNodeData(J).get_prev_temp();
    }
    return Tsr /= nbf;
}

double SolidEquation::compute_average_temp_prev(const TALYFEMLIB::FEMElm& fe, int nbf2) {
	const int nbf = fe.pElm->n_nodes();

    double Tsr = 0.0;
    for(int i = 0; i < nbf; ++i) {
        int J = fe.pElm->ElemToLocalNodeID(i);
        Tsr += p_data_->GetNodeData(J).get_prev_minus_1_temp();
    }
    return Tsr /= nbf;
}

void SolidEquation::fillEssBC() {
    initEssBC();
}

void SolidEquation::Integrands4side(const TALYFEMLIB::FEMElm& fe, int sideInd,
		TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be) {

    NewtonBC* bc;
    if(idata->get_bc(sideInd, bc)) {
    	bc->calculate(fe, Ae, be, dt_);
    }

}

bool SolidEquation::Integrands4contact(TALYFEMLIB::FEMElm& fe, int sideInd,
		TALYFEMLIB::ZeroMatrix<double>& Ae1, TALYFEMLIB::ZeroMatrix<double>& Ae2,
		TALYFEMLIB::ZEROARRAY<double>& be) {

	ContactBC* bc;
	if (idata->get_bc(sideInd, bc)) {
		bc->calculate(fe, Ae1, Ae2, be, dt_);
		return true;
    } else {
    	return false;
    }

}

void SolidEquation::Integrands(const TALYFEMLIB::FEMElm& fe, TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be) {
    const int nsd = p_grid_->nsd();
    const int nbf = fe.pElm->n_nodes();
    const double detJxW = fe.detJxW();
    double Tsr = compute_average_temp(fe, nbf);
    double Tpsr = compute_average_temp_prev(fe, nbf);
    double Vsr = compute_average_velocity(fe, nbf);
    int mat_ind = fe.pElm->mat_ind();
    //TODO: Change to calculated values
    const SolidMaterial& material = idata->get_material(mat_ind);
    double lambda = material.conductivity(Tsr, Vsr);
    double capacity = material.heat_capacity(fe, p_data_, Tsr, Tpsr, Vsr);
//    if(capacity != capacity) {
//        material.heat_capacity(fe, p_data_, Tsr, Tpsr, Vsr);
//        PrintError("SolidEquation::Integrands -> capacity != capacity");
//    }

    for (int a = 0; a < nbf; a++) {
	    for (int b = 0; b < nbf; b++) {
			double M = capacity * fe.N(a)*fe.N(b)*detJxW;
            double N = 0;
            for (int k = 0; k < nsd; k++) {
				N +=  lambda * fe.dN(a,k)*fe.dN(b,k)*detJxW;
            }
            Ae(a,a) += M/dt_;
            Ae(a,b) += N;
            int J = fe.pElm->ElemToLocalNodeID(b);
            be(a) += M/dt_*p_data_->GetNodeData(J).get_prev_temp();
        }
	    p_data_->GetNodeData(fe.pElm->ElemToLocalNodeID(a)).set_capprox(capacity);
    }
}

void SolidEquation::InitializeContactBC(ContactBounds* cb) {
	contact_bounds_ = cb;

	has_contact_bc_.redim(n_total_dof_);
	has_contact_bc_.fill(false);

	for (int i = 0; i < p_grid_->n_nodes(); i++) {
		for (int k = 0; k < n_dof_; k++) { // TODO: not all DOFs have contact BC
			has_contact_bc_.set(i*n_dof_ + k, contact_bounds_->IsNodePeriodic(i));
		}
	}
}

void SolidEquation::Assemble(bool assemble_surface) {
//	PrintInfo("Assemble from SolidEquation");

	try {
		PetscErrorCode ierr;
		ierr = UpdateMatPreallocation(); CHKERRV(ierr);
		MatSetOption(Ag_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);  //Fixme: This has to be removed (when custom preallocator will be finished)
		// zero the stiffness matrix and load
		if (recalc_matrix_) { MatZeroEntries(Ag_); }
		// trying to do VecZeroEntries (bg_);
		PetscInt nlocal;
		double *array;
		VecGetLocalSize(bg_, &nlocal);
		VecGetArray(bg_, &array);
		memset(array, 0, sizeof(double)*nlocal);
		VecRestoreArray(bg_, &array);

		if (has_uniform_mesh_) {
		  AssembleVolumeUniformMesh(assemble_surface);
		} else {
		  AssembleVolume(assemble_surface);
		}
		if (assemble_surface) {
		  AssembleSurface();
		}
		if (recalc_matrix_) {
		  ierr = MatAssemblyBegin(*p_Ag_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
		  ierr = MatAssemblyEnd(*p_Ag_, MAT_FLUSH_ASSEMBLY); CHKERRV(ierr);
		}
		ierr = VecAssemblyBegin(*p_bg_); CHKERRV(ierr);
		ierr = VecAssemblyEnd(*p_bg_); CHKERRV(ierr);

    } catch (TALYException& e) {
      // if an exception was thrown (by the error handler or otherwise),
      // we need to pop the error handler on the way up the stack
      PetscPopErrorHandler();
      throw e;
    }

    PetscPopErrorHandler();
}

void SolidEquation::AssembleVolume(bool assemble_surface) {
//  PrintInfo("AssembleVolume from SolidEquation");
  FEMElm fe;
  for (int elmID = 0; elmID < p_grid_->n_elements(); elmID++) {
    if (!IsMyElement(elmID)) continue;
    fe.refill(p_grid_, elmID);
    // int order = 0;//p_grid_->basis_order();
    fe.setRelativeOrder(order_);
    int n = fe.pElm->n_nodes() * n_dof_;
    ZeroMatrix<double> Ae;
    if (recalc_matrix_) {
      Ae.redim(n, n);
      Ae.fill(0);
    }
    ZEROARRAY<double> be;
    be.redim(n);
    be.fill(0);
    AssembleElement(elmID, Ae, be, fe, assemble_surface);

    if (contact_bounds_->has_periodic()) {
        ZeroMatrix<double> Ae1, Ae2;
        if (recalc_matrix_) {
          Ae1.redim(n, n); Ae2.redim(n, n);
          Ae1.fill(0); Ae2.fill(0);
        }
        ZEROARRAY<double> be1;
        be1.redim(n);
        be1.fill(0);
    	AssembleElementContact(elmID, Ae1, Ae2, be1, fe, assemble_surface);
    }

  }
}

void SolidEquation::AssembleElementContact(int elmID,
		ZeroMatrix<double>& Ae1, ZeroMatrix<double>& Ae2,
		ZEROARRAY<double>& be1, FEMElm& fe, bool assemble_surface) {

//	PrintInfo("AssembleElementContact");

	bool flag = false; // flag to check if really there is a need to add contact BC to global matrix

//	PrintInfo("AssembleElementContact::fe.pElm->surface_indicator_.size() = ", fe.pElm->surface_indicator_.size());

	for (ELEM::SurfaceList_type::const_iterator it = fe.pElm->surface_indicator_.begin();
			it != fe.pElm->surface_indicator_.end(); it++) {

		fe.refill(p_grid_, elmID, &(*it));
		fe.setRelativeOrder(order_);
		fe.initNumItg();

		while (fe.moreItgPoints()) {
			fe.update4nextItgPt();
			for (unsigned int i = 0; i < SurfaceIndicator::MAX_SURFACE_INDICATORS; i++) {
//				PrintInfo("AssembleElementContact::it->has_indicator(", i, ")", it->has_indicator(i));
				if (it->has_indicator(i)) {
					flag = Integrands4contact(fe, i, Ae1, Ae2, be1) or flag;
				}
			}
		}
	}

	if (flag)
		AssembleAebeWithContact(fe, elmID, Ae1, Ae2, be1);

}

void SolidEquation::AssembleAebeWithContact(FEMElm& fe, int elmID,
		ZeroMatrix<double>& Ae1, ZeroMatrix<double>& Ae2, ZEROARRAY<double>& be) {

//	PrintInfo("AssembleAebeWithContact");

	ZEROARRAY<PetscInt> vertex_arr1, vertex_col_arr1; // global indices on first side
	ZEROARRAY<PetscInt> vertex_arr2, vertex_col_arr2; // global indices on second side

	CalcAebeIndicesWithContact(fe, vertex_arr1, vertex_col_arr1,
									vertex_arr2, vertex_col_arr2);

	AssembleAebeWithIndex(Ae1.data(), be.data(), vertex_arr1.size(), vertex_arr1.data(),
								vertex_col_arr1.size(), vertex_col_arr1.data());

	be.fill(0.0); //we should not assemble the same 'be' multiple times
	AssembleAebeWithIndex(Ae2.data(), be.data(), vertex_arr1.size(), vertex_arr1.data(),
								vertex_col_arr2.size(), vertex_col_arr2.data()); //this requires special handling in preallocate


	Ae1.ratio(-1.0); // values on the second side have opposite sign
	Ae2.ratio(-1.0); // values on the second side have opposite sign

	AssembleAebeWithIndex(Ae1.data(), be.data(), vertex_arr2.size(), vertex_arr2.data(),
								vertex_col_arr1.size(), vertex_col_arr1.data()); //this requires special handling in preallocate
	AssembleAebeWithIndex(Ae2.data(), be.data(), vertex_arr2.size(), vertex_arr2.data(),
								vertex_col_arr2.size(), vertex_col_arr2.data());

}


void SolidEquation::CalcAebeIndicesWithContact(FEMElm& fe,
		ZEROARRAY<PetscInt>& rows_out1, ZEROARRAY<PetscInt>& cols_out1,
		ZEROARRAY<PetscInt>& rows_out2, ZEROARRAY<PetscInt>& cols_out2) {

//	PrintInfo("CalcAebeIndicesWithContact");

	const int vertexn = n_dof_ * fe.pElm->n_nodes();
	rows_out1.redim(vertexn);
	cols_out1.redim(vertexn);
	PetscInt* rows_ptr1 = rows_out1.data();
	PetscInt* cols_ptr1 = cols_out1.data();
	rows_out2.redim(vertexn);
	cols_out2.redim(vertexn);
	PetscInt* rows_ptr2 = rows_out2.data();
	PetscInt* cols_ptr2 = cols_out2.data();

	if (p_grid_->parallel_type_ == kNoDomainDecomp) {
		for (int i = 0; i < fe.pElm->n_nodes(); i++) {
			for (int k = 0; k < n_dof_; k++) {
				LocalVarIdx lid = fe.pElm->node_id_array(i) * n_dof_ + k;
				int idx = i * n_dof_ + k;

				if (has_ess_bc_.get(lid)) {
					rows_ptr1[idx] = -1;
				} else {
					rows_ptr1[idx] = lid;
				}
				cols_ptr1[idx] = lid;

//				PrintInfo("gid = ", gid, "\thas_contact_bc_.get(gid) = ", has_contact_bc_.get(gid));
				if (has_contact_bc_.get(lid)) {
					int lclnodeID = fe.pElm->node_id_array(i);

					int newNode = contact_bounds_->GetPeriodicSolPartner(lclnodeID);
					// we don't want this to override an essential condition
					if (!has_ess_bc_.get(lid)) {
						rows_ptr2[idx] = newNode * n_dof_ + k;
					} else {
						rows_ptr2[idx] = -1;
					}
					cols_ptr2[idx] = newNode * n_dof_ + k;
//					PrintInfo("gid = ", gid, "\topposite gid = ", newNode);
				} else {
					rows_ptr2[idx] = -1;
					cols_ptr2[idx] = -1;
//					PrintInfo("gid = ", gid);
//					PrintError("CalcAebeIndicesWithContact called for surface without contact data");
				}
			}
		}
	} else {

//		PrintError("Contact boudary condition not implemented for ",
//														"Domain Decomposition");
//		exit(1);
		for (int i = 0; i < fe.pElm->n_nodes(); i++) {
			for (int k = 0; k < n_dof_; k++) {
				LocalVarIdx lid = fe.pElm->node_id_array(i) * n_dof_ + k;
		        GlobalVarIdx gid = (p_grid_->solution_map(fe.pElm->node_id_array(i)))
		                           * n_dof_ + k;
				int idx = i * n_dof_ + k;

				if (has_ess_bc_.get(lid)) {
					rows_ptr1[idx] = -1;
				} else {
					rows_ptr1[idx] = gid;
				}
				cols_ptr1[idx] = gid;

//				PrintInfo("gid = ", gid, "\thas_contact_bc_.get(gid) = ", has_contact_bc_.get(gid));
				if (has_contact_bc_.get(lid)) {
					LocalNodeID lclnodeID = fe.pElm->node_id_array(i);

					SolutionNodeID newNode = contact_bounds_->GetPeriodicSolPartner(lclnodeID);
					// we don't want this to override an essential condition
					if (!has_ess_bc_.get(lid)) {
						rows_ptr2[idx] = newNode * n_dof_ + k;
					} else {
						rows_ptr2[idx] = -1;
					}
					cols_ptr2[idx] = newNode * n_dof_ + k;
//					PrintInfo("gid = ", gid, "\topposite gid = ", newNode);
				} else {
					rows_ptr2[idx] = -1;
					cols_ptr2[idx] = -1;
//					PrintInfo("gid = ", gid);
//					PrintError("CalcAebeIndicesWithContact called for surface without contact data");
				}
			}
		}


	}
}


void SolidEquation::compute_additional_values() {
    compute_solid_fraction();
    compute_real_solidus_temperature();
    compute_grain_size();
    compute_heat_flux();
//    compute_capprox();
}

void SolidEquation::compute_solid_fraction() {
    for(int elemID = 0; elemID < p_grid_->n_elements(); elemID++) {
        ELEM *elem = p_grid_->elm_array_[elemID];
        int mat_ind = elem->mat_ind();

		for(int i = 0; i < elem->n_nodes(); i++) {
			SolidMaterial& solid_material = idata->get_material(mat_ind);
			SolidNodeData* pData = &p_data_->GetNodeData(elem->ElemToLocalNodeID(i));//&(Node(elem->node_id_array(i)));
			double v = pData->get_velocity();
			double Twe = pData->get_prev_temp();
			double fs = solid_material.get_solidification_model().solid_phase_fraction(Twe, v);
			pData->set_solid_fraction(fs);
		}
    }
}

void SolidEquation::compute_real_solidus_temperature() {
    for(int elemID = 0; elemID < p_grid_->n_elements(); elemID++) {
        ELEM *elem = p_grid_->elm_array_[elemID];
	int mat_ind = elem->mat_ind();

	for(int i = 0; i < elem->n_nodes(); i++) {
            SolidMaterial& solid_material = idata->get_material(mat_ind);
            SolidNodeData* pData = &p_data_->GetNodeData(elem->ElemToLocalNodeID(i));//&(Node(elem->node_id_array(i)));
            double v = pData->get_velocity();
            double ts = solid_material.get_solidification_model().real_solidus_temperature(v);
            pData->set_real_solidus_temperature(ts);
	}
    }
}

void SolidEquation::compute_grain_size() {
    for(int elemID = 0; elemID < p_grid_->n_elements(); elemID++) {
        ELEM *elem = p_grid_->elm_array_[elemID];
	int mat_ind = elem->mat_ind();

	for(int i = 0; i < elem->n_nodes(); i++) {
            SolidMaterial& solid_material = idata->get_material(mat_ind);
            SolidNodeData* pData = &p_data_->GetNodeData(elem->ElemToLocalNodeID(i));//&(Node(elem->node_id_array(i)));
            double v = pData->get_velocity();
            double rz = solid_material.grain_size(v);
            pData->set_grain_size(rz);
	}
    }
}

void SolidEquation::compute_capprox() {
    for(int elemID = 0; elemID < p_grid_->n_elements(); elemID++) {
        ELEM *elem = p_grid_->elm_array_[elemID];

        FEMElm fe;
        fe.refill(p_grid_, elemID);
        fe.setRelativeOrder(0);
        fe.initNumItg();

        const int nbf = fe.pElm->n_nodes();
        double Tsr = compute_average_temp(fe, nbf);
        double Tpsr = compute_average_temp_prev(fe, nbf);
        double Vsr = compute_average_velocity(fe, nbf);

        int mat_ind = fe.pElm->mat_ind();
//        int mat_ind = elem->mat_ind();
        const SolidMaterial& material = idata->get_material(mat_ind);

        double capacity = material.heat_capacity(fe, p_data_, Tsr, Tpsr, Vsr);

//        if (capacity != capacity) {
//            material.heat_capacity(fe, p_data_, Tsr, Tpsr, Vsr);
//            PrintError("SolidEquation::compute_capprox -> capacity != capacity");
//        }

		for(int i = 0; i < elem->n_nodes(); i++) {
			SolidNodeData* pData = &p_data_->GetNodeData(elem->ElemToLocalNodeID(i));
//			double v = pData->get_velocity();
//			double T = pData->get_prev_temp();
//			double Tprev = pData->get_prev_minus_1_temp();
//			double capacity = material.heat_capacity(fe, p_data_, T, Tprev, v);

			pData->set_capprox(capacity);
		}
    }
}

void SolidEquation::compute_heat_flux() {
    for(int elemID = 0; elemID < p_grid_->n_elements(); elemID++) {
    }
}


