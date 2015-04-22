#ifndef SOLID_EQUATION_H
#define SOLID_EQUATION_H

#include <FEM/cequation.h>
#include "solid_node_data.hpp"
#include "solid_input_data.hpp"

class SolidEquation
: public TALYFEMLIB::CEquation<SolidNodeData> {
    public:
        SolidEquation(SolidInputData* id);
        virtual void Solve(double t, double dt);
        virtual void Integrands(TALYFEMLIB::FEMElm& fe, 
                TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be);
        virtual void Integrands4side(TALYFEMLIB::FEMElm& fe, 
                int sideInd, TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be);            
        virtual void fillEssBC();
    private:                
        SolidInputData* idata;
        void compute_additional_values();
        void compute_solid_fraction();
        void compute_real_solidus_temperature();
        void compute_grain_size();
        void compute_heat_flux();
        double compute_average_temp(TALYFEMLIB::FEMElm& fe, int nbf);
        double compute_average_temp_prev(TALYFEMLIB::FEMElm& fe, int nbf);
        double compute_average_velocity(TALYFEMLIB::FEMElm& fe, int nbf);
};

#endif
