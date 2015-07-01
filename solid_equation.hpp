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
        void Solve(double t, double dt) override;
        void Integrands(const TALYFEMLIB::FEMElm& fe,
                TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be) override;
        void Integrands4side(const TALYFEMLIB::FEMElm& fe, int sideInd,
        		TALYFEMLIB::ZeroMatrix<double>& Ae, TALYFEMLIB::ZEROARRAY<double>& be) override;
        void fillEssBC() override;

        // additional methods to use contact BC
        void InitializeContactBC(TALYFEMLIB::ContactBounds* cb);

        // methods adding contact BC to main equation
        virtual bool Integrands4contact(TALYFEMLIB::FEMElm& fe, int sideInd,
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

        void compute_additional_values();

    private:
        SolidInputData* idata;
        TALYFEMLIB::ContactBounds* contact_bounds_;
        TALYFEMLIB::ZEROARRAY<bool> has_contact_bc_;

        void compute_solid_fraction();
        void compute_real_solidus_temperature();
        void compute_grain_size();
        void compute_heat_flux();
        void compute_capprox();
        double compute_average_temp(const TALYFEMLIB::FEMElm& fe, int nbf);
        double compute_average_temp_prev(const TALYFEMLIB::FEMElm& fe, int nbf);
        double compute_average_velocity(const TALYFEMLIB::FEMElm& fe, int nbf);
};

#endif
