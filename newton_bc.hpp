#ifndef NEWTON_BC_COEFF_HPP
#define	NEWTON_BC_COEFF_HPP

#include <Grid/femelm.h>

class NewtonBC {
public:
    double heat_exchange_coeff() const {
        return alpha_;
    }
    
    double ambient_temperature() const {
        return Tamb_;    
    }
    void calculate(TALYFEMLIB::FEMElm& fe, TALYFEMLIB::ZeroMatrix<double>& Ae, 
        TALYFEMLIB::ZEROARRAY<double>& be); 
private:
    double alpha_, Tamb_;
};

#endif	/* NEWTON_BC_COEFF_HPP */

