#include "solid_equation.hpp"

SolidEquation::SolidEquation(SolidInputData* id) : idata(id) {
}

void SolidEquation::Solve(double t, double dt) {
    this->t_ = t;
    this->dt_ = dt;

    fillEssBC();
    ApplyEssBCToSolution();
    Assemble(1);
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

void SolidEquation::Integrands4side(TALYFEMLIB::FEMElm& fe,
                    int sideInd, TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be) {
    double alpha = idata->heat_exchange_coeff();
    double Tamb = idata->ambient_temperature();

    if (sideInd >= 1 && sideInd <= 6) {
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


