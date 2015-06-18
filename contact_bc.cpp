/*
 * contact_bc.cpp
 *
 *  Created on: 18-06-2015
 *      Author: rdyja
 */

#include "contact_bc.hpp"

void ContactBC::calculate(const TALYFEMLIB::FEMElm& fe,
		TALYFEMLIB::ZeroMatrix<double>& Ae1, TALYFEMLIB::ZeroMatrix<double>& Ae2,
        TALYFEMLIB::ZEROARRAY<double>& be, double dt) {

	const int nbf = fe.pElm->n_nodes();
	const double detSideJxW = fe.detJxW();

	for (int a = 0; a < nbf; ++a) {
		for (int b = 0; b < nbf; ++b) {
			double M = fe.N(a) * fe.N(b) * detSideJxW;
			Ae1(a,b) += kappa_ * M/dt;
			Ae2(a,b) += -kappa_ * M/dt;
		}
	}

}

