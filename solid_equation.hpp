#ifndef SOLID_EQUATION_H
#define SOLID_EQUATION_H

#include <FEM/cequation.h>
#include "solid_node_data.hpp"
#include "solid_input_data.hpp"
#include "contact_bounds.hpp"

class SolidEquation
: public TALYFEMLIB::CEquation<SolidNodeData> {
    public:
        SolidEquation(SolidInputData* id, TALYFEMLIB::ContactBounds* cb);
        virtual void Solve(double t, double dt);
        virtual void Integrands(TALYFEMLIB::FEMElm& fe,
                TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be);
        virtual void Integrands4side(TALYFEMLIB::FEMElm& fe, int sideInd,
        		TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be);
        virtual void fillEssBC();

        // methods adding contact BC to main equation
        virtual void Integrands4contact(TALYFEMLIB::FEMElm& fe, int sideInd,
        		TALYFEMLIB::ZeroMatrix<double>& Ae1, TALYFEMLIB::ZeroMatrix<double>& Ae2,
        		TALYFEMLIB::ZEROARRAY<double>& be);
        void AssembleElementContact(int elmID, TALYFEMLIB::ZeroMatrix<double>& Ae1,
        		TALYFEMLIB::ZeroMatrix<double>& Ae2, TALYFEMLIB::ZEROARRAY<double>& be,
				TALYFEMLIB::FEMElm& fe, bool assemble_surface = true);
        void AssembleAebeWithContact(TALYFEMLIB::FEMElm& fe, int elmID,
        		TALYFEMLIB::ZeroMatrix<double>& Ae1, TALYFEMLIB::ZeroMatrix<double>& Ae2,
				TALYFEMLIB::ZEROARRAY<double>& be);
        void CalcAebeIndicesWithContact(TALYFEMLIB::FEMElm& fe,
        		TALYFEMLIB::ZEROARRAY<PetscInt>& rows_out1, TALYFEMLIB::ZEROARRAY<PetscInt>& cols_out1,
				TALYFEMLIB::ZEROARRAY<PetscInt>& rows_out2, TALYFEMLIB::ZEROARRAY<PetscInt>& cols_out2);

        // redefined methods from base class in order to add contact BC
        virtual void Assemble(bool assemble_surface = true);
        virtual void AssembleVolume(bool assemble_surface = true);   
        
        
    private:
        SolidInputData* idata;        
        TALYFEMLIB::ContactBounds* contact_bounds_;           
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
