// Copyright 2004- ZISIO ICIS, Czestochowa University of Technology
// http://icis.pcz.pl


#ifndef NSC_SOLIDIFICATION_MODEL_H
#define NSC_SOLIDIFICATION_MODEL_H

const int EQUILIBRIUM_SOLIDIFICATION_TYPE = 1;
const int SCHEIL_SOLIDIFICATION_TYPE = 2;
const int INDIRECT_SOLIDIFICATION_TYPE = 3;

class SolidMaterial;

class SolidificationModel
{
public:
    SolidificationModel(const SolidMaterial& heat_material, int model_id);
    double solid_phase_fraction(double T, double vT) const;
    double real_solidus_temperature(double vT) const;
	bool is_in_eutectic_range(double T, double vT) const;

private:
    double equilibrium_solid_phase_fraction(double T) const;
    double scheil_solid_phase_fraction(double T) const;
    double indirect_solid_phase_fraction(double T, double vT) const;

    double epsT_;
    const SolidMaterial& mat_;
    int solidification_model_id_;
	double k_;
};

#endif
