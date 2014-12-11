#include "solid_equation.hpp"

SolidEquation::SolidEquation(SolidInputData* id) : idata(id) {
}

void SolidEquation::Solve(double t, double dt, int bCalculateMatrix) {
    this->t = t;
    this->dt = dt;

    fillEssBC();
    ApplyEssBCToSolution();
    Assemble(1);
    ApplyEssBC(linear);

    SolveKSP(solution, 1, 0);

    for(int A=1; A<=pGrid->nodeno; A++){
        double newval=solution(A);
        //pData->Node(A).UpdateDataStructures();
        pData->Node(A).set_curr_temp(newval);
    }
}

double SolidEquation::compute_average_velocity(TALYFEMLIB::FEMElm& fe, int nbf) {
    double Vsr = 0.0;
    for(int i = 1; i <= nbf; ++i) {
        int J = fe.pElm->glbNodeID(i);
        Vsr += pData->Node(J).get_velocity();
    }
    return Vsr /= nbf;
}

double SolidEquation::compute_average_temp(TALYFEMLIB::FEMElm& fe, int nbf) {
    double Tsr = 0.0;
    for(int i = 1; i <= nbf; ++i) {
        int J = fe.pElm->glbNodeID(i);
        Tsr += pData->Node(J).get_prev_temp();
    }
    return Tsr /= nbf;
}

double SolidEquation::compute_average_temp_prev(TALYFEMLIB::FEMElm& fe, int nbf) {
    double Tsr = 0.0;
    for(int i = 1; i <= nbf; ++i) {
        int J = fe.pElm->glbNodeID(i);
        Tsr += pData->Node(J).get_prev_minus_1_temp();
    }
    return Tsr /= nbf;
}

void SolidEquation::fillEssBC() {
    initEssBC();
}

void SolidEquation::integrands4side(TALYFEMLIB::FEMElm& fe,
                    int sideInd, TALYFEMLIB::Matrix<double>& Ae, TALYFEMLIB::ARRAY<double>& be) {
    double alpha = idata->heat_exchange_coeff();
    double Tamb = idata->ambient_temperature();

    if (sideInd >= 1 && sideInd <= 6) {
        int nbf = fe.pElm->nodeno;
        double detSideJxW = fe.detJxW();
        for (int a = 1; a <= nbf; ++a) {
            for(int b = 1; b <= nbf; ++b) {
                double M = fe.N(a) * detSideJxW;
                Ae(a,b) += alpha * M/dt;
                be(a) += alpha * Tamb * M/dt;
            }
        }
    }
}

void SolidEquation::integrands(TALYFEMLIB::FEMElm& fe, TALYFEMLIB::Matrix<double>& Ae, TALYFEMLIB::ARRAY<double>& be) {
    const int nsd = pGrid->nsd;
    const int nbf = fe.pElm->nodeno;
    const double detJxW = fe.detJxW();
    double Tsr = compute_average_temp(fe, nbf);
    double Tpsr = compute_average_temp_prev(fe, nbf);
    double Vsr = compute_average_velocity(fe, nbf);
    //TODO: Change to calculated values
    double lambda = idata->get_material().conductivity(Tsr, Vsr);
    double capacity = idata->get_material().heat_capacity(Tsr, Tpsr, Vsr);

    for(int a=1; a<=nbf; a++){
	    for(int b=1; b<=nbf; b++){
			double M = capacity * fe.N(a)*fe.N(b)*detJxW;
            double N = 0;
            for(int k=1; k<=nsd; k++){
				N +=  lambda * fe.dN(a,k)*fe.dN(b,k)*detJxW;
            }
            Ae(a,b) += M/dt + N;
            int J = fe.pElm->glbNodeID (b);
            be(a) += M/dt*pData->Node(J).get_prev_temp();
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


