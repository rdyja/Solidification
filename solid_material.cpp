// Copyright 2004- ZISIO ICIS, Czestochowa University of Technology
// http://icis.pcz.pl


// author Arkadiusz Nagorka, nagorka@icis.pcz.pl

#include "solid_material.hpp"
#include "solid_model.hpp"
#include "enthalpy_model.hpp"
#include "solid_grid_field.hpp"

#include <string>
#include <stdexcept>
#include <cmath>

#include <cfloat>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <cmath>

using namespace std;

const double DOUBLE_VERY_SMALL = 1E-16;

SolidMaterial::SolidMaterial(int index)
//: Material(MY_MATERIAL_TYPE, 16, m, 0)
{
    initialize_property_map();
}

void SolidMaterial::initialize_property_map() {
    properties["lambdaS"] = 0.0;
    properties["cS"] = 0.0;
    properties["rhoS"] = 0.0;
    //properites["isCasting",
    properties["T0"] = 0.0;
    properties["sourceTerm"] = 0.0;
    properties["Ts"] = 0.0;    
    properties["Tl"] = 0.0;
    properties["Te"] = 0.0;
    properties["Tp"] = 0.0;
    properties["lambdaL"] = 0.0;
    properties["cL"] = 0.0;
    properties["rhoL"] = 0.0;
    properties["L"] = 0.0;          
    properties["maxGrainSize"] = 0.0;
    properties["coeffBF"] = 0.0;  
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

void SolidMaterial::set_heatcapacity_model(int model_id)
{
    heatcapacity_model_id_ = model_id;
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

double SolidMaterial::heat_capacity(const TALYFEMLIB::FEMElm& fe,
		TALYFEMLIB::GridField<SolidNodeData>* p_data, double T, double T_p, double vT) const
{
 if (!is_casting()) return specific_heat(T)*density(T);

  const double TL = liquidus_temperature();
  const double TS = solidification_model_->real_solidus_temperature(vT);

  if (T < TS || T > TL) {
    return specific_heat(T, vT)*density(T, vT);
  }

  assert(enthalpy_model_.get() != NULL);

  switch (heatcapacity_model_id_) {

  case MORGAN_HEAT_CAPACITY: {

	  const double H = enthalpy_model_->enthalpy_in_mushy_zone(T, vT);
	  const double H_p = enthalpy_model_->enthalpy_in_mushy_zone(T_p, vT);

  	  return (H - H_p)/(T - T_p);
  }

  case COMMINI_HEAT_CAPACITY: {
	  const int nsd = fe.pElm->nsd();
	  const int nnd = fe.pElm->n_nodes();

	  vector<double> Tprev(nnd), Vprev(nnd), Hprev(nnd);
	  for (ElemNodeID a = 0; a < nnd; a++) {
		  const int J = fe.pElm->ElemToLocalNodeID(a);
		  Tprev[a] = p_data->GetNodeData(J).get_prev_temp();
		  Vprev[a] = p_data->GetNodeData(J).get_velocity();
		  Hprev[a] = enthalpy_model_->enthalpy_in_mushy_zone(Tprev[a], Vprev[a]);
	  }

	  bool diff = true;
	  for (ElemNodeID a = 0; a < nnd; a++) {
		  if (fabs(Tprev[a] - Tprev[a == nnd - 1 ? 0 : a + 1]) < DOUBLE_VERY_SMALL) {
			  diff = false;
			  break;
		  }
	  }

	  if (diff) {
		double sum = 0.0;
		for (int dir = 0; dir < nsd; dir++) {
			double sum_n = 0.0, sum_d = 0.0;
			for (ElemNodeID a = 0; a < fe.pElm->n_nodes(); a++) {
				sum_n += (fe.dN(a, dir) * Hprev[a]);
				sum_d += (fe.dN(a, dir) * Tprev[a]);
			}
			sum += sum_n/sum_d;
		}
		return sum/nsd;

	  } else {
		if (fabs(Tprev[0]) < DOUBLE_VERY_SMALL) {
			return specific_heat(T, vT)*density(T, vT);
		} else {
			return Hprev[0]/Tprev[0] ;
		}
	  }



  }

  case LEMMON_HEAT_CAPACITY: {
	  const int nsd = fe.pElm->nsd();
	  const int nnd = fe.pElm->n_nodes();

	  vector<double> Tprev(nnd), Vprev(nnd), Hprev(nnd);
	  for (ElemNodeID a = 0; a < nnd; a++) {
		  const int J = fe.pElm->ElemToLocalNodeID(a);
		  Tprev[a] = p_data->GetNodeData(J).get_prev_temp();
		  Vprev[a] = p_data->GetNodeData(J).get_velocity();
		  Hprev[a] = enthalpy_model_->enthalpy_in_mushy_zone(Tprev[a], Vprev[a]);
	  }

	  bool diff = true;
	  for (ElemNodeID a = 0; a < nnd; a++) {
		  if (fabs(Tprev[a] - Tprev[a == nnd - 1 ? 0 : a + 1]) < DOUBLE_VERY_SMALL) {
			  diff = false;
			  break;
		  }
	  }
	  if (diff) {
		  double sum_n_2 = 0.0, sum_d_2 = 0.0;
		  for (int dir = 0; dir < nsd; dir++) {
			double sum_n = 0.0, sum_d = 0.0;
			for (ElemNodeID a = 0; a < fe.pElm->n_nodes(); a++) {
				sum_n += fe.dN(a, dir) * Hprev[a];
				sum_d += fe.dN(a, dir) * Tprev[a];
			}
			sum_n_2 += pow(sum_n, 2.0);
			sum_d_2 += pow(sum_d, 2.0);
		  }
		  return sqrt(sum_n_2/sum_d_2);

  	  } else {
  		if (fabs(Tprev[0]) < DOUBLE_VERY_SMALL) {
  			return specific_heat(T, vT)*density(T, vT);
  		} else {
  			return Hprev[0]/Tprev[0] ;
  		}
  	  }
  }

  case DELGUIDICE_HEAT_CAPACITY: {
	  const int nsd = fe.pElm->nsd();
	  const int nnd = fe.pElm->n_nodes();

	  vector<double> Tprev(nnd), Vprev(nnd), Hprev(nnd);
	  for (ElemNodeID a = 0; a < nnd; a++) {
		  const int J = fe.pElm->ElemToLocalNodeID(a);
		  Tprev[a] = p_data->GetNodeData(J).get_prev_temp();
		  Vprev[a] = p_data->GetNodeData(J).get_velocity();
		  Hprev[a] = enthalpy_model_->enthalpy_in_mushy_zone(Tprev[a], Vprev[a]);
	  }

	  bool diff = true;
	  for (ElemNodeID a = 0; a < nnd; a++) {
		  if (fabs(Tprev[a] - Tprev[a == nnd - 1 ? 0 : a + 1]) < DOUBLE_VERY_SMALL) {
			  diff = false;
			  break;
		  }
	  }

	  if (diff) {
		double sum_n_2 = 0.0, sum_d_2 = 0.0;
		for (int dir = 0; dir < nsd; dir++) {
			double sum_h = 0.0, sum_t = 0.0;
			for (ElemNodeID a = 0; a < nnd; a++) {
				sum_h += fe.dN(a, dir) * Hprev[a];// * fe.dN(a, dir) * Tprev[a];
				sum_t += fe.dN(a, dir) * Tprev[a];// * fe.dN(a, dir) * Tprev[a];
			}
			sum_n_2 += sum_h*sum_t;
			sum_d_2 += pow(sum_t, 2.0);
		}
		return sum_n_2/sum_d_2;

	  } else {
		if (fabs(Tprev[0]) < DOUBLE_VERY_SMALL) {
			return specific_heat(T, vT)*density(T, vT);
		} else {
			return Hprev[0]/Tprev[0] ;
		}
	  }
  }

  default:
	  throw(std::string("Unknown heat capacity approximation type"));
  }
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
                    return MaterialProperty(properties.at("lambdaS"));
		case 1:
			return MaterialProperty(properties.at("cS"));
		case 2:
			return MaterialProperty(properties.at("rhoS"));
		case 3:
			return MaterialProperty(0.0);
		case 4:
			return MaterialProperty(properties.at("T0"));
		case 5:
			return MaterialProperty(properties.at("sourceTerm"));
		case 6:
			return MaterialProperty(properties.at("Ts"));
		case 7:
			return MaterialProperty(properties.at("Tl"));
		case 8:
			return MaterialProperty(properties.at("Tp"));
		case 9:
			return MaterialProperty(properties.at("Te"));
		case 10:
			return MaterialProperty(properties.at("lambdaL"));
		case 11:
			return MaterialProperty(properties.at("cL"));
		case 12:
			return MaterialProperty(properties.at("rhoL"));
		case 13:
			return MaterialProperty(properties.at("L"));
		case 14:
			return MaterialProperty(properties.at("maxGrainSize"));
		case 15:
			return MaterialProperty(properties.at("coeffBF"));
	}
	throw std::string("Reading wrong material property");	
}

void SolidMaterial::set_property(const string& propertyName, double value) {                    
                        
    if(properties.find(propertyName) == properties.end()) {
        throw std::string("Setting wrong material property " + propertyName);	
    }
    properties[propertyName] = value;    
}

void SolidMaterial::update_k() {
    k_ = 1.0 - (liquidus_temperature() - solidus_temperature())/
               (pure_melting_temperature() - solidus_temperature());
}

double SolidMaterial::get_k() const {
	return k_;
}
