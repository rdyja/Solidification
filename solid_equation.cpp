#include <Grid/femelm.h>
#include <DataStructure/zeroarray.h>
#include <DataStructure/zeromatrix.h>

#include "solid_equation.hpp"

using namespace TALYFEMLIB;

SolidEquation::SolidEquation(SolidInputData* id) :
	idata(id), contact_bounds_(NULL) {
}

void SolidEquation::Solve(double t, double dt) {
    this->t_ = t;
    this->dt_ = dt;

    fillEssBC();
    ApplyEssBCToSolution();
    Assemble();
    ApplyEssBC();

    SolveKSP(solution_, 1, 0);

    for(int A=0; A < p_grid_->n_nodes(); A++){
        double newval=solution_(A);
        //pData->Node(A).UpdateDataStructures();
        p_data_->Node(A).set_curr_temp(newval);
    }
}

double SolidEquation::compute_average_velocity(TALYFEMLIB::FEMElm& fe, int nbf) {
    double Vsr = 0.0;
    for(int i = 0; i < nbf; ++i) {
        int J = fe.pElm->glbNodeID(i);
        Vsr += p_data_->Node(J).get_velocity();
    }
    return Vsr /= nbf;
}

double SolidEquation::compute_average_temp(TALYFEMLIB::FEMElm& fe, int nbf) {
    double Tsr = 0.0;
    for(int i = 0; i < nbf; ++i) {
        int J = fe.pElm->glbNodeID(i);
        Tsr += p_data_->Node(J).get_prev_temp();
    }
    return Tsr /= nbf;
}

double SolidEquation::compute_average_temp_prev(TALYFEMLIB::FEMElm& fe, int nbf) {
    double Tsr = 0.0;
    for(int i = 0; i < nbf; ++i) {
        int J = fe.pElm->glbNodeID(i);
        Tsr += p_data_->Node(J).get_prev_minus_1_temp();
    }
    return Tsr /= nbf;
}

void SolidEquation::fillEssBC() {
    initEssBC();
}

void SolidEquation::Integrands4side(TALYFEMLIB::FEMElm& fe, int sideInd,
		TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be) {

    double alpha = idata->heat_exchange_coeff();
    double Tamb = idata->ambient_temperature();

    if (sideInd >= 1 && sideInd <= 6) { //we need to decide which indicators are for the 3rd type BC
        int nbf = fe.pElm->nodeno;
        double detSideJxW = fe.detJxW();
        for (int a = 0; a < nbf; ++a) {
            for(int b = 0; b < nbf; ++b) {
                double M = fe.N(a) * detSideJxW;
                Ae(a,b) += alpha * M/dt_;
                be(a) += alpha * Tamb * M/dt_;
            }
        }
    }
}

bool SolidEquation::Integrands4contact(TALYFEMLIB::FEMElm& fe, int sideInd,
		TALYFEMLIB::ZeroMatrix<double>& Ae1, TALYFEMLIB::ZeroMatrix<double>& Ae2,
		TALYFEMLIB::ZEROARRAY<double>& be) {

    if (sideInd >= 7 && sideInd <= 9) { //we need to decide which indicators are for the 4th type BC
        int nbf = fe.pElm->nodeno;
        double detSideJxW = fe.detJxW();
        double kappa = 0.0;
        for (int a = 0; a < nbf; ++a) {
            for (int b = 0; b < nbf; ++b) {
            	double M = fe.N(a) * detSideJxW;
            	Ae1(a,b) += -kappa * M/dt_;;
            	Ae2(a,b) +=  kappa * M/dt_;
            }
        }
        return true;
    } else {
    	return false;
    }

}

void SolidEquation::Integrands(TALYFEMLIB::FEMElm& fe, TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be) {
    const int nsd = p_grid_->nsd();
    const int nbf = fe.pElm->nodeno;
    const double detJxW = fe.detJxW();
    double Tsr = compute_average_temp(fe, nbf);
    double Tpsr = compute_average_temp_prev(fe, nbf);
    double Vsr = compute_average_velocity(fe, nbf);
    //TODO: Change to calculated values
    double lambda = idata->get_material().conductivity(Tsr, Vsr);
    double capacity = idata->get_material().heat_capacity(Tsr, Tpsr, Vsr);

    for(int a=0; a< nbf; a++){
	    for(int b=0; b < nbf; b++){
			double M = capacity * fe.N(a)*fe.N(b)*detJxW;
            double N = 0;
            for(int k=0; k < nsd; k++){
				N +=  lambda * fe.dN(a,k)*fe.dN(b,k)*detJxW;
            }
            Ae(a,b) += M/dt_ + N;
            int J = fe.pElm->glbNodeID (b);
            be(a) += M/dt_*p_data_->Node(J).get_prev_temp();
        }
    }
}

void SolidEquation::InitializeContactBC(ContactBounds* cb) {
	contact_bounds_ = cb;

	has_contact_bc_.redim(n_total_dof_);
	has_contact_bc_.fill(false);

	for (int i = 0; i < p_grid_->n_nodes(); i++) {
		for (int k = 0; k < n_dof_; k++) { // TODO: not all DOFs have to have contact BC
			has_contact_bc_.set(i*n_dof_ + k, contact_bounds_->IsNodePeriodic(i));
		}
	}
}

PetscErrorCode SolidEquation::Assemble(bool assemble_surface) {
	PrintInfo("Assemble from SolidEquation");

    PetscErrorCode ierr;
    ierr = UpdateMatPreallocation(); CHKERRQ(ierr);

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
      ierr = MatAssemblyBegin(*p_Ag_, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(*p_Ag_, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(*p_bg_); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*p_bg_); CHKERRQ(ierr);

	return 0;
}

void SolidEquation::AssembleVolume(bool assemble_surface) {
  PrintInfo("AssembleVolume from SolidEquation");
  FEMElm fe;
  for (int elmID = 0; elmID < p_grid_->n_elements(); elmID++) {
    if (!IsMyElement(elmID)) continue;
    fe.refill(p_grid_, elmID);
    // int order = 0;//p_grid_->basis_order();
    fe.setRelativeOrder(order_);
    int n = fe.pElm->nodeno * n_dof_;
    ZeroMatrix<double> Ae;
    if (recalc_matrix_) {
      Ae.redim(n, n);
      Ae.fill(0);
    }
    ZEROARRAY<double> be;
    be.redim(n);
    be.fill(0);
    AssembleElement(elmID, Ae, be, fe, assemble_surface);

    if (contact_bounds_ != NULL) {
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

	for (ELEM::SurfaceList_type::iterator it = fe.pElm->surfaceIndicator.begin();
			it != fe.pElm->surfaceIndicator.end(); it++) {

		fe.refill(p_grid_, elmID, it->surfaceID());
		fe.setRelativeOrder(order_);
		fe.initNumItg();

		while (fe.moreItgPoints()) {
			fe.update4nextItgPt();
			for (unsigned int i = 0; i < SurfaceIndicator::MAX_SURFACE_INDICATORS; i++) {
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

	AssembleAebeWithIndex(Ae1, be, vertex_arr1.size(), vertex_arr1.data(),
								vertex_col_arr1.size(), vertex_col_arr1.data());

	be.fill(0.0); //we should not assemble the same 'be' multiple times

	AssembleAebeWithIndex(Ae2, be, vertex_arr1.size(), vertex_arr1.data(),
								vertex_col_arr2.size(), vertex_col_arr2.data()); //this requires special handling in preallocate

	Ae1.ratio(-1.0); // values on the second side have opposite sign
	Ae2.ratio(-1.0); // values on the second side have opposite sign

	AssembleAebeWithIndex(Ae1, be, vertex_arr2.size(), vertex_arr2.data(),
								vertex_col_arr1.size(), vertex_col_arr1.data()); //this requires special handling in preallocate
	AssembleAebeWithIndex(Ae2, be, vertex_arr2.size(), vertex_arr2.data(),
								vertex_col_arr2.size(), vertex_col_arr2.data());

}


void SolidEquation::CalcAebeIndicesWithContact(FEMElm& fe,
		ZEROARRAY<PetscInt>& rows_out1, ZEROARRAY<PetscInt>& cols_out1,
		ZEROARRAY<PetscInt>& rows_out2, ZEROARRAY<PetscInt>& cols_out2) {

//	PrintInfo("CalcAebeIndicesWithContact");

	const int vertexn = n_dof_ * fe.pElm->nodeno;
	rows_out1.redim(vertexn);
	cols_out1.redim(vertexn);
	PetscInt* rows_ptr1 = rows_out1.data();
	PetscInt* cols_ptr1 = cols_out1.data();
	rows_out2.redim(vertexn);
	cols_out2.redim(vertexn);
	PetscInt* rows_ptr2 = rows_out2.data();
	PetscInt* cols_ptr2 = cols_out2.data();

	if (p_grid_->parallel_type_ == kNoDomainDecomp) {
		for (int i = 0; i < fe.pElm->nodeno; i++) {
			for (int k = 0; k < n_dof_; k++) {
				int gid = fe.pElm->pNodeIDArray[i] * n_dof_ + k;
				int idx = i * n_dof_ + k;

				if (has_ess_bc_.get(gid)) {
					rows_ptr1[idx] = -1;
				} else {
					rows_ptr1[idx] = gid;
				}
				cols_ptr1[idx] = gid;

				if (has_contact_bc_.get(gid)) {
					int lclnodeID = fe.pElm->pNodeIDArray[i];
					int newNode = contact_bounds_->GetPeriodicSolPartner(lclnodeID);
					// we don't want this to override an essential condition
					if (!has_ess_bc_.get(gid)) {
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

		PrintError("Contact boudary condition not implemented for ",
														"Domain Decomposition");
		exit(1);
	}
}


void SolidEquation::compute_additional_values() {
    compute_solid_fraction();
    compute_real_solidus_temperature();
    compute_grain_size();
    compute_heat_flux();
}

void SolidEquation::compute_solid_fraction() {

}

void SolidEquation::compute_real_solidus_temperature() {

}

void SolidEquation::compute_grain_size() {

}

void SolidEquation::compute_heat_flux() {

}


