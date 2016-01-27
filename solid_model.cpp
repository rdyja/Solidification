// Copyright 2004- ZISIO ICIS, Czestochowa University of Technology
// http://icis.pcz.pl

#include "solid_model.hpp"
#include "solid_material.hpp"
#include <cmath>
#include <cassert>


SolidificationModel::SolidificationModel(const SolidMaterial& heat_material,
                                         int model_id) :
    epsT_(1.0),
    mat_(heat_material),
    solidification_model_id_(model_id)
{
}

double SolidificationModel::solid_phase_fraction(double T, double vT) const
{
    if (!mat_.is_casting()) return 1.0;

    const double TL = mat_.liquidus_temperature();
    const double TS = mat_.solidus_temperature();
    const double TE = mat_.eutectic_temperature();

    if (T > TL) return 0.0;

    switch(solidification_model_id_) {

    case EQUILIBRIUM_SOLIDIFICATION_TYPE:
        if (T < TS) return 1.0;
        else return equilibrium_solid_phase_fraction(T);

        break;

    case SCHEIL_SOLIDIFICATION_TYPE: {
        if (T < TE - epsT_) return 1.0;
        else {
            if (T > TE + epsT_) return scheil_solid_phase_fraction(T);
            else {  // eutectic range
                double fsE = scheil_solid_phase_fraction(TE);
                return 1.0 - (1.0 - fsE)*((mat_.enthalpy(T) - mat_.enthalpy(TE - epsT_))/
                                (mat_.enthalpy(TE + epsT_) - mat_.enthalpy(TE - epsT_)));
            }
        }
        break;
    }

    case INDIRECT_SOLIDIFICATION_TYPE: {
        double fsE = indirect_solid_phase_fraction(TE, vT);
        if (fsE > 1.0) {
            if (T < TS) return 1.0;
            else return indirect_solid_phase_fraction(T, vT);
        }
        else {
            double TSE = real_solidus_temperature(vT);
            if (T < TSE - epsT_) return 1.0;
            else {
                if (T > TSE + epsT_) return indirect_solid_phase_fraction(T, vT);
                else {  // eutectic range
                    return 1.0 - (1.0 - fsE)*((mat_.enthalpy(T, vT) - mat_.enthalpy(TE - epsT_, vT))/
                            (mat_.enthalpy(TE + epsT_, vT) - mat_.enthalpy(TE - epsT_, vT)));
                }
            }
        }
        break;
    }

    default:
        assert(0); return 0;
    }
}

double SolidificationModel::equilibrium_solid_phase_fraction(double T) const
{
    const double TL = mat_.liquidus_temperature();
    const double TM = mat_.pure_melting_temperature();
    const double k = mat_.get_k();

    return 1.0/(1.0 - k)*(TL - T)/(TM - T);
}

double SolidificationModel::scheil_solid_phase_fraction(double T) const
{
    const double TL = mat_.liquidus_temperature();
    const double TM = mat_.pure_melting_temperature();
    const double k = mat_.get_k();

    return 1.0 - pow((TM - T)/(TM - TL), 1.0/(k - 1.0));
}

double SolidificationModel::indirect_solid_phase_fraction(double T, double vT) const
{
    const double TL = mat_.liquidus_temperature();
    const double TM = mat_.pure_melting_temperature();
    const double k = mat_.get_k();
    const double alfa = mat_.coefficient_BF()/pow(mat_.grain_size(vT), 2.0);
    const double omega = alfa*(1.0 - exp(-1.0/alfa)) - 0.5*exp(-1.0/(2.0*alfa));
    const double n = mat_.coefficient_n();

    return 1.0/(1.0 - n*k*omega) * (1.0 - pow((TM - T)/(TM - TL), (1.0 - n*k*omega)/(k - 1.0)));
}

double SolidificationModel::real_solidus_temperature(double vT) const
{
    if (!mat_.is_casting()) return 0.0;

    switch (solidification_model_id_) {

    case EQUILIBRIUM_SOLIDIFICATION_TYPE:
        return mat_.solidus_temperature();

    case SCHEIL_SOLIDIFICATION_TYPE:
        return mat_.eutectic_temperature();

    case INDIRECT_SOLIDIFICATION_TYPE:
        assert(vT > -1.0);
        if (indirect_solid_phase_fraction(mat_.eutectic_temperature(), vT) > 1.0)
            return mat_.solidus_temperature();

    else {
		const double k = mat_.get_k();
		const double n = mat_.coefficient_n();
        double alfa = mat_.coefficient_BF()/pow(mat_.grain_size(vT), 2.0);
        double nkOmega = n*k*(alfa*(1.0 - exp(-1.0/alfa)) - 0.5*exp(-1.0/(2.0*alfa)));

        double TSE = mat_.pure_melting_temperature() - (mat_.pure_melting_temperature() - mat_.liquidus_temperature())*pow(nkOmega,(k-1.0)/(1.0-nkOmega));
        if (TSE < mat_.eutectic_temperature()) TSE = mat_.eutectic_temperature();
        if (TSE > mat_.solidus_temperature()) TSE = mat_.solidus_temperature();

        return TSE;
    }

    default:
        assert(0); return 0;
    }
}

bool SolidificationModel::is_in_eutectic_range(double T, double vT) const
{
    if (mat_.is_casting()) {
        double TS = real_solidus_temperature(vT);
        return T > (TS - epsT_) && T < (TS + epsT_) ? true : false;
    } else {
        return false;
    }
}
