#include "newton_bc.hpp"

void NewtonBC::calculate(const TALYFEMLIB::FEMElm& fe, 
        TALYFEMLIB::ZeroMatrix<double>& Ae, 
        TALYFEMLIB::ZEROARRAY<double>& be, double dt) {
    
    double detSideJxW = fe.detJxW();
    int nbf = fe.pElm->n_nodes();
    
    for (int a = 0; a < nbf; ++a) {
        for(int b = 0; b < nbf; ++b) {
            double M = fe.N(a) * detSideJxW;
            Ae(a,b) += alpha_ * M/dt;
            be(a) += alpha_ * Tamb_ * M/dt;
        }
    }    
}
