// Copyright 2004- ZISIO ICIS, Czestochowa University of Technology
// http://icis.pcz.pl

#include <cassert>
#include "enthalpy_model.hpp"
#include "solid_model.hpp"
#include "solid_material.hpp"

EnthalpyModel::EnthalpyModel(SolidMaterial& heat_material,
                             SolidificationModel& solidification_model,
                             int model_id) :
    mat_(heat_material),
    enthalpy_model_id_(model_id),
    sol_(solidification_model)
{ }


double EnthalpyModel::enthalpy_in_mushy_zone(double T, double vT)
{

    double TL = mat_.liquidus_temperature();
    double TS = sol_.real_solidus_temperature(vT);
    double cs = mat_.specific_heat(TS-1, vT), rs = mat_.density(TS-1, vT);
    double cl = mat_.specific_heat(TL+1, vT), rl = mat_.density(TL+1, vT);
    double L = mat_.latent_heat();
    double HL = cs*rs*TS + rs*L + 0.5 * (cs*rs + cl*rl) * (TL - TS);
    double HS = cs*rs*TS;

    switch (enthalpy_model_id_) {

    case LINEAR_ENTHALPY_TYPE: {
        return HS + (T - TS)*(HL - HS)/(TL - TS);
    }

    case QUADRATIC_ENTHALPY_TYPE: {
        double a, b, c;
        abc(TL, TS, HL, HS, cs*rs,  a, b, c);
        return a*T*T + b*T + c;
    }

    case PIECEWISE_QUADRATIC_ENTHALPY_TYPE: {
        double a, b, c;
        abcabc(TL, TS, HL, HS, cl*rl, cs*rs, T, a, b, c);
        return a*T*T + b*T + c;
    }

    default: {
        assert(0); return 0.0;
    }

    }
}

double EnthalpyModel::heat_capacity_in_mushy_zone(double T, double vT)
{
    double TL = mat_.liquidus_temperature();
    double TS = sol_.real_solidus_temperature(vT);
    double cs = mat_.specific_heat(TS-1.0, vT), rs = mat_.density(TS-1.0, vT);
    double cl = mat_.specific_heat(TL+1.0, vT), rl = mat_.density(TL+1.0, vT);
    double L = mat_.latent_heat();
    double HL = cs*rs*TS + rs*L + 0.5 * (cs*rs + cl*rl) * (TL-TS);
    double HS = cs*rs*TS;

    switch (enthalpy_model_id_) {

    case LINEAR_ENTHALPY_TYPE: {
        return (HL - HS)/(TL - TS);
    }

    case QUADRATIC_ENTHALPY_TYPE: {
        double a, b, c;
        abc(TL, TS, HL, HS, cs*rs,  a, b, c);
        return 2.0*a*T + b;
    }

    case PIECEWISE_QUADRATIC_ENTHALPY_TYPE: {
        double a, b, c;
        abcabc(TL, TS, HL, HS, cl*rl, cs*rs, T, a, b, c);
        return 2.0*a*T + b;
    }

    default: {
        assert(0); return 0.0;
    }

    }
}

void EnthalpyModel::abc(double TL, double TS, double HL, double HS, double Ss,
                        double& a, double& b, double& c)
{
    double denom = (TL-TS)*(TL-TS);
    a = (HL - HS + Ss*(-TL + TS))/denom;
    b = (-2.0*HL*TS + 2.0*HS*TS + Ss*(TL*TL - TS*TS)) /denom;
    c = (-2.0*HS*TL*TS + HL*(TS*TS) + HS*(TL*TL)
        + Ss*(TL*(TS*TS) - (TL*TL)*TS)) /denom;
}

void EnthalpyModel::abcabc(double TL, double TS, double HL, double HS,
                           double  CL, double CS, double T,
                           double& a, double& b, double& c)
{
    const double alpha = 0.2; // how far is the interior point Tp from TL?

    assert(TL > TS);
    assert(HL > HS);

    double Tp = TL - alpha * (TL - TS);

    // solution of the 6x6 linear system directly copied from MuPAD
    if (T < Tp) {
        a = -(HL/(TL*TS - TL*Tp + TS*Tp - TS*TS)) + HS/(TL*TS - TL*Tp + TS*Tp - TS*TS) + (CL*(TL - Tp))/(2.0*TL*TS - 2.0*TL*Tp + 2.0*TS*Tp - (TS*TS)*2.0) + (CS*(-TL + TS*2.0 - Tp))/(-2.0*TL*TS + 2.0*TL*Tp - 2.0*TS*Tp + (TS*TS)*2.0);
        b = (2.0*HS*TS)/(-(TL*TS) + TL*Tp - TS*Tp + TS*TS) + (2.0*HL*TS)/(TL*TS - TL*Tp + TS*Tp - TS*TS) + (CL*(-(TL*TS) + TS*Tp))/(TL*TS - TL*Tp + TS*Tp - TS*TS) + (CS*(-(TL*Tp) + TS*TS))/(TL*TS - TL*Tp + TS*Tp - TS*TS);
        c =  -((HL*(TS*TS))/(TL*TS - TL*Tp + TS*Tp - TS*TS)) + (HS*(TL*TS - TL*Tp + TS*Tp))/(TL*TS - TL*Tp + TS*Tp - TS*TS) + (CL*(-(TL*(TS*TS)) + (TS*TS)*Tp))/(-2.0*TL*TS + 2.0*TL*Tp - 2.0*TS*Tp + (TS*TS)*2.0) + (CS*(2.0*TL*TS*Tp - TL*(TS*TS) - (TS*TS)*Tp))/(2.0*TL*TS - 2.0*TL*Tp + 2.0*TS*Tp - (TS*TS)*2.0);
    }
    else {
        a =  -(HL/(-(TL*TS) - TL*Tp + TS*Tp + TL*TL)) + HS/(-(TL*TS) - TL*Tp + TS*Tp + TL*TL) + (CS*(-TS + Tp))/(-2.0*TL*TS - 2.0*TL*Tp + 2.0*TS*Tp + (TL*TL)*2.0) + (CL*(TL*2.0 - TS - Tp))/(-2.0*TL*TS - 2.0*TL*Tp + 2.0*TS*Tp + (TL*TL)*2.0);
        b =  (2.0*HL*TL)/(-(TL*TS) - TL*Tp + TS*Tp + TL*TL) - (2.0*HS*TL)/(-(TL*TS) - TL*Tp + TS*Tp + TL*TL) + (CS*(TL*TS - TL*Tp))/(-(TL*TS) - TL*Tp + TS*Tp + TL*TL) + (CL*(TS*Tp - TL*TL))/(-(TL*TS) - TL*Tp + TS*Tp + TL*TL);
        c =   (HS*(TL*TL))/(-(TL*TS) - TL*Tp + TS*Tp + TL*TL) + (HL*(-(TL*TS) - TL*Tp + TS*Tp))/(-(TL*TS) - TL*Tp + TS*Tp + TL*TL) + (CS*(-((TL*TL)*TS) + (TL*TL)*Tp))/(-2.0*TL*TS - 2.0*TL*Tp + 2.0*TS*Tp + (TL*TL)*2.0) + (CL*(-2.0*TL*TS*Tp + (TL*TL)*TS + (TL*TL)*Tp))/(-2.0*TL*TS - 2.0*TL*Tp + 2.0*TS*Tp + (TL*TL)*2.0);
    }
}
