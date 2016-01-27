// Copyright 2004- ZISIO ICIS, Czestochowa University of Technology
// http://icis.pcz.pl


#ifndef NSC_HEAT_MATERIAL_H
#define NSC_HEAT_MATERIAL_H

// author Arkadiusz Nagorka, nagorka@icis.pcz.pl

//#include "solid_model.hpp"
//#include "enthalpy_model.hpp"
#include <memory>
#include <map>
#include <string>
#include <algorithm>

#include <Grid/femelm.h>
#include <Grid/gridfield.h>

#include "enthalpy_model.hpp"

const int MORGAN_HEAT_CAPACITY = 1;
const int DELGUIDICE_HEAT_CAPACITY = 2;
const int LEMMON_HEAT_CAPACITY = 3;
const int COMMINI_HEAT_CAPACITY = 4;

class MaterialProperty {
		double param;
	public:
	   explicit MaterialProperty(double p = 0.0) : param(p)
	   {}
	   double evaluate(double T = 0.0) const {
		   return param;
	   }
};

//~ class EnthalpyModel;
class SolidificationModel;
class SolidNodeData;

class SolidMaterial
{
public:
  explicit SolidMaterial(int index = 0);

  SolidMaterial(SolidMaterial&& sm):
		solidification_model_(std::move(sm.solidification_model_)),
		enthalpy_model_(std::move(sm.enthalpy_model_)),
		properties(std::move(sm.properties)),
		isCasting_(sm.isCasting_),
		name_(sm.name_),
		heatcapacity_model_id_(sm.heatcapacity_model_id_)
		{
		}
  SolidMaterial& operator=(SolidMaterial&& sm) {
	  solidification_model_=std::move(sm.solidification_model_);
	  enthalpy_model_ = std::move(sm.enthalpy_model_);
	  properties =std::move(sm.properties);
	  isCasting_ = sm.isCasting_;
	  name_ = sm.name_;
	  heatcapacity_model_id_ = sm.heatcapacity_model_id_;
	  return *this;
  }

  SolidMaterial(const SolidMaterial&) = delete;
  SolidMaterial& operator=(const SolidMaterial&) = delete;

  double conductivity(double T = 0.0, double vT = -1.0) const;
  double specific_heat(double T = 0.0, double vT = -1.0) const;
  double density(double T = 0.0, double vT = -1.0) const;
  double source_term(double T = 0) const;
  double initial_temperature() const;

  bool is_casting() const;
  double solidus_temperature() const;
  double liquidus_temperature() const;
  double pure_melting_temperature() const;
  double latent_heat() const;

  double heat_capacity(const TALYFEMLIB::FEMElm& fe, TALYFEMLIB::GridField<SolidNodeData>* p_data,
		  double T = 0.0, double T_p = 0.0, double vT = -1.0) const;
  double enthalpy(double T, double vT = -1.0) const;
  std::string name() const;

  const SolidificationModel& get_solidification_model() const;
  void set_solidification_model(int model_id);
  void set_enthalpy_model(int model_id);
  void set_heatcapacity_model(int id);
  void set_casting(bool casting);
  void set_name(std::string);

  void update_k();
  double get_k() const;

  double max_grain_size() const;
  double coefficient_BF() const;
  double coefficient_n() const;
  double grain_size(double vT) const;
  double eutectic_temperature() const;
  bool is_in_eutectic_range(double T, double vT) const;
  MaterialProperty get_property(int num) const;
  void set_property(const std::string&, double);
  void initialize_property_map();

private:
  std::unique_ptr<SolidificationModel> solidification_model_;
  std::unique_ptr<EnthalpyModel> enthalpy_model_;
  //std::map<std::string, MaterialProperty> properties;
  std::vector<MaterialProperty> properties;
  std::map<std::string, int> map_property_name_to_num;
  double k_;

  bool isCasting_;
  std::string name_;
  int heatcapacity_model_id_;
};

class SolidMaterialFactory
{
public:
  //virtual SolidMaterial* create(const char* material_type, const Mesh& mesh);
};

// ===== inline

inline void SolidMaterial::set_casting(bool casting) {
    isCasting_ = casting;
}

inline double SolidMaterial::source_term(double T) const
{
  return get_property(3).evaluate(T);
}

inline double SolidMaterial::initial_temperature() const
{
  return get_property(4).evaluate();
}

inline bool SolidMaterial::is_casting() const
{
  return isCasting_;
}

inline double SolidMaterial::solidus_temperature() const
{
  return get_property(6).evaluate();
}

inline double SolidMaterial::liquidus_temperature() const
{
  return get_property(7).evaluate();
}

inline double SolidMaterial::pure_melting_temperature() const
{
  return get_property(8).evaluate();
}

inline double SolidMaterial::eutectic_temperature() const
{
  return get_property(9).evaluate();
}

inline double SolidMaterial::latent_heat() const
{
  return get_property(13).evaluate();
}

inline double SolidMaterial::max_grain_size() const
{
  return get_property(14).evaluate();
}

inline double SolidMaterial::coefficient_BF() const
{
  return get_property(15).evaluate();
}

inline double SolidMaterial::coefficient_n() const
{
  return get_property(16).evaluate();
}

inline const SolidificationModel& SolidMaterial::get_solidification_model() const
{
    return *solidification_model_;
}

inline std::string SolidMaterial::name() const {
    return name_;
}

inline void SolidMaterial::set_name(std::string n) {
    name_ = n;
}

#endif
