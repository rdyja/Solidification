#ifndef SOLID_EQUATION_H
#define SOLID_EQUATION_H

#include <FEM/FEMSolver.h>
#include "solid_node_data.hpp"
#include "solid_input_data.hpp"

class SolidEquation
: public TALYFEMLIB::CEquation<SolidNodeData, TALYFEMLIB::GPData, TALYFEMLIB::SurfaceGPData> {
    public:
        SolidEquation(SolidInputData* id);
        virtual void Solve(double t, double dt, int bCalculateMatrix = 1);
        virtual void integrands(TALYFEMLIB::FEMElm& fe, 
                TALYFEMLIB::Matrix<double>& Ae, TALYFEMLIB::ARRAY<double>& be);
        virtual void integrands4side(TALYFEMLIB::FEMElm& fe, 
                int sideInd, TALYFEMLIB::Matrix<double>& Ae, TALYFEMLIB::ARRAY<double>& be);            
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
