/*
 * contact_bc.hpp
 *
 *  Created on: 18-06-2015
 *      Author: rdyja
 */

#ifndef CONTACT_BC_HPP_
#define CONTACT_BC_HPP_

#include <Grid/femelm.h>

class ContactBC {
public:
	ContactBC(double kappa, std::string n)
    : kappa_(kappa), name_(n) {
    }
    double heat_exchange_coeff() const {
        return kappa_;
    }

    std::string name() const {
        return name_;
    }

    void calculate(const TALYFEMLIB::FEMElm& fe,
    	TALYFEMLIB::ZeroMatrix<double>& Ae1, TALYFEMLIB::ZeroMatrix<double>& Ae2,
        TALYFEMLIB::ZEROARRAY<double>& be, double dt);
private:
    double kappa_;
    std::string name_;
};



#endif /* CONTACT_BC_HPP_ */
