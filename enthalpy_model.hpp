// Copyright 2004- ZISIO ICIS, Czestochowa University of Technology
// http://icis.pcz.pl


#ifndef NSC_ENTHALPY_MODEL_H
#define NSC_ENTHALPY_MODEL_H

const int LINEAR_ENTHALPY_TYPE = 1;
const int QUADRATIC_ENTHALPY_TYPE = 2;
const int PIECEWISE_QUADRATIC_ENTHALPY_TYPE = 3;

class SolidMaterial;
class SolidificationModel;

class EnthalpyModel
{
public:
    EnthalpyModel(SolidMaterial& heat_material,
	              SolidificationModel& solidification_model,
	              int model_id);
	double enthalpy_in_mushy_zone(double T, double vT = -1.0);
	double heat_capacity_in_mushy_zone(double T, double vT = -1.0);

private:
	void abc(double TL, double TS, double HL, double HS, double Ss,
	         double& a, double& b, double& c);

	void abcabc(double TL, double TS, double HL, double HS,
                double  CL, double CS, double T,
                double& a, double& b, double& c);

	SolidMaterial& mat_;
	int enthalpy_model_id_;
	SolidificationModel& sol_;

};

#endif
