#include "newton_bc_coeff.hpp"

void NewtonBC::calculate(TALYFEMLIB::FEMElm& fe, 
        TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be) {
    
    double detSideJxW = fe.detJxW();
    
    for (int a = 0; a < nbf; ++a) {
        for(int b = 0; b < nbf; ++b) {
            double M = fe.N(a) * detSideJxW;
            Ae(a,b) += alpha * M/dt_;
            be(a) += alpha * Tamb * M/dt_;
        }
    }    
}
