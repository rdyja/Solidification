#ifndef NEWTON_BC_COEFF_HPP
#define	NEWTON_BC_COEFF_HPP

#include <Grid/femelm.h>

class NewtonBC {
public:
    NewtonBC(double alpha, double Tamb, std::string n)
    : alpha_(alpha),Tamb_(Tamb), name_(n) {   
    }
    double heat_exchange_coeff() const {
        return alpha_;
    }
    
    double ambient_temperature() const {
        return Tamb_;    
    }
    std::string name() const {
        return name_;
    }
    
    void calculate(const TALYFEMLIB::FEMElm& fe, TALYFEMLIB::ZeroMatrix<double>& Ae, 
        TALYFEMLIB::ZEROARRAY<double>& be, double dt); 
private:
    double alpha_, Tamb_;
    std::string name_;
};

#endif	/* NEWTON_BC_COEFF_HPP */

