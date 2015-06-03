// Copyright 2004- ZISIO ICIS, Czestochowa University of Technology
// http://icis.pcz.pl


// author Arkadiusz Nagorka, nagorka@icis.pcz.pl

#include "solid_material.hpp"
#include "solid_model.hpp"
#include "enthalpy_model.hpp"

#include <string>
#include <stdexcept>
#include <cmath>

#include <cfloat>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>

using namespace std;

//static const char* MY_MATERIAL_TYPE = "NSC.Solidification.Material";

SolidMaterial::SolidMaterial(int index)
//: Material(MY_MATERIAL_TYPE, 16, m, 0)
{
    initialize_property_map();
}

void SolidMaterial::initialize_property_map() {
    properties["lambdaS"] = &lambdaS_;
    properties["cS"] = &cS_;
    properties["rhoS"] = &rhoS_;
    //properites["isCasting",
    properties["T0"] = &T0_;
    properties["sourceTerm"] = &sourceTerm_;
    properties["Ts"] = &Ts_;    
    properties["Tl"] = &Tl_;
    properties["Te"] = &Te_;
    properties["Tp"] = &Tp_;
    properties["lambdaL"] = &lambdaL_;
    properties["cL"] = &cL_;
    properties["rhoL"] = &rhoL_;
    properties["L"] = &L_;          
    properties["maxGrainSize"] = &maxGrainSize_;
    properties["coeffBF"] = &coeffBF_;  
}

double SolidMaterial::grain_size(double vT) const
{
  if (is_casting())
    return max_grain_size() * (1.0 - exp(-1.0/vT));
  else
    return 0.0;
}

void SolidMaterial::set_solidification_model(int model_id)
{
    solidification_model_ = unique_ptr<SolidificationModel>(new SolidificationModel(*this, model_id));
}

void SolidMaterial::set_enthalpy_model(int model_id)
{
    enthalpy_model_ = unique_ptr<EnthalpyModel>(new EnthalpyModel(*this, *solidification_model_, model_id));
}

double SolidMaterial::conductivity(double T, double vT) const
{
  if (!is_casting()) return get_property(0).evaluate(T);

  double TS = solidification_model_->real_solidus_temperature(vT), TL = liquidus_temperature();

  if (T < TS) {
    return get_property(0).evaluate(T);
  }
  else if (T > TL) {
    return get_property(10).evaluate(T);
  }
  else {
    // mushy zone
    double ks = get_property(0).evaluate(TS-1);
    double kl = get_property(10).evaluate(TL+1);
    double f = solidification_model_->solid_phase_fraction(T, vT);
    return f*ks + (1-f)*kl;
  }

  return 0;
}


double SolidMaterial::specific_heat(double T, double vT) const
{
  if (!is_casting()) return get_property(1).evaluate(T);

  double TS = solidification_model_->real_solidus_temperature(vT), TL = liquidus_temperature();

  if (T < TS) {
    return get_property(1).evaluate(T);
  }
  else if (T > TL) {
    return get_property(11).evaluate(T);
  }
  else {
    // mushy zone
    double cs =  get_property(1).evaluate(TS-1);
    double cl =  get_property(11).evaluate(TL+1);
    double f = solidification_model_->solid_phase_fraction(T, vT);
    return f*cs + (1-f)*cl;
  }

  return 0;
}

double SolidMaterial::density(double T, double vT) const
{
  if (!is_casting()) return get_property(2).evaluate(T);

  double TS = solidification_model_->real_solidus_temperature(vT), TL = liquidus_temperature();

  if (T < TS) {
    return get_property(2).evaluate(T);
  }
  else if (T > TL) {
    return get_property(12).evaluate(T);
  }
  else {
    // mushy zone
    double rs =  get_property(2).evaluate(TS-1);
    double rl =  get_property(12).evaluate(TL+1);
    double f = solidification_model_->solid_phase_fraction(T, vT);
    return f*rs + (1-f)*rl;
  }

  return 0;
}

double SolidMaterial::heat_capacity(double T, double T_p, double vT) const
{
  if (!is_casting()) return specific_heat(T)*density(T);

  double TL = liquidus_temperature(), TS = solidification_model_->real_solidus_temperature(vT);
  if (T < TS || T > TL) {
    return specific_heat(T, vT)*density(T, vT);
  }

  assert(enthalpy_model_.get() != NULL);
  return (enthalpy_model_->enthalpy_in_mushy_zone(T, vT) - enthalpy_model_->enthalpy_in_mushy_zone(T_p, vT))/(T - T_p);
}

double SolidMaterial::enthalpy(double T, double vT) const
{
  if (!is_casting()) return specific_heat(T) * density(T) * T;

  double TL = liquidus_temperature(), TS = solidification_model_->real_solidus_temperature(vT);
  double H = 0;

  if (T < TS - 1.0) {
    // solid
    H = specific_heat(T, vT) * density(T, vT) * T;
  }
  else {

    if (T < TL + 1.0) {
      // mushy zone

      H = enthalpy_model_->enthalpy_in_mushy_zone(T, vT);
    }
    else {
      // liquid
      double cs = specific_heat(TS-1, vT), rs = density(TS-1, vT);
      double cl = specific_heat(TL+1, vT), rl = density(TL+1, vT);
      double L = latent_heat();
      double HL = cs*rs*TS + rs*L + 0.5 * (cs*rs + cl*rl) * (TL-TS);

      H = HL + cl*rl*(T-TL);
    }
  }

  return H;
}

/*Material* SolidMaterialFactory::create(const char* type, const Mesh& mesh)
{
  if (string(type) == MY_MATERIAL_TYPE) {
    return new SolidMaterial(0, mesh);
  }
  else {
    // >>> here another material type can be handled
  }

  // unknown material type
  return 0;
}*/

MaterialProperty SolidMaterial::get_property(int num) const{
	switch(num) {
		case 0:
			return MaterialProperty(lambdaS_);
		case 1:
			return MaterialProperty(cS_);
		case 2:
			return MaterialProperty(rhoS_);
		case 3:
			return MaterialProperty(0.0);
		case 4:
			return MaterialProperty(T0_);
		case 5:
			return MaterialProperty(sourceTerm_);
		case 6:
			return MaterialProperty(Ts_);
		case 7:
			return MaterialProperty(Tl_);
		case 8:
			return MaterialProperty(Tp_);
		case 9:
			return MaterialProperty(Te_);
		case 10:
			return MaterialProperty(lambdaL_);
		case 11:
			return MaterialProperty(cL_);
		case 12:
			return MaterialProperty(rhoL_);
		case 13:
			return MaterialProperty(L_);
		case 14:
			return MaterialProperty(maxGrainSize_);
		case 15:
			return MaterialProperty(coeffBF_);
	}
	throw std::string("Reading wrong material property");	
}

void SolidMaterial::set_property(const string& propertyName, double value) {                    
                        
    if(properties.find(propertyName) == properties.end()) {
        throw std::string("Setting wrong material property " + propertyName);	
    }
    *properties[propertyName] = value;    
}

void SolidMaterial::update_k() {
    k_ = 1.0 - (liquidus_temperature() - solidus_temperature())/
               (pure_melting_temperature() - solidus_temperature());
}

double SolidMaterial::get_k() const {
	return k_;
}
