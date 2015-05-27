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

class MaterialProperty {
		double param;
	public:
	   MaterialProperty(double p) : param(p)
	   {}
	   double evaluate(double T = 0.0) const {
		   return param;
	   }
};

class EnthalpyModel;
class SolidificationModel;

class SolidMaterial
{
public:
  explicit SolidMaterial(int index = 0);  
  
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

  double heat_capacity(double T = 0.0, double T_p = 0.0, double vT = -1.0) const;
  double enthalpy(double T, double vT = -1.0) const;
  std::string name() const;

  const SolidificationModel& get_solidification_model() const;
  void set_solidification_model(int model_id);
  void set_enthalpy_model(int model_id);
  void set_casting(bool casting);
  void set_name(std::string);

  void update_k();
  double get_k() const;

  double max_grain_size() const;
  double coefficient_BF() const;
  double grain_size(double vT) const;
  double eutectic_temperature() const;
  bool is_in_eutectic_range(double T, double vT) const;
  MaterialProperty get_property(int num) const;
  void set_property(const std::string&, double);
  
private:
  std::unique_ptr<SolidificationModel> solidification_model_;
  std::unique_ptr<EnthalpyModel> enthalpy_model_;
  std::map<std::string, double*> properties;
  double rhoL_, rhoS_, cL_, cS_, lambdaL_, lambdaS_, Tl_, Ts_, L_, T0_, k_;
  //Tp_ - pure melting temperature
  //Te_ - euthectic temperature
  double Tp_, Te_, maxGrainSize_, coeffBF_, sourceTerm_;
  bool isCasting_;
  std::string name_;
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
