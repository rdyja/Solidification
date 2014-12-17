#ifndef SOLID_INPUT_DATA_H
#define SOLID_INPUT_DATA_H

#include <string>
#include <memory>
#include <InputData/InputData.hpp>
#include "solid_model.hpp"
#include "enthalpy_model.hpp"
#include "solid_material.hpp"

class SolidInputData : public TALYFEMLIB::InputData {
	public:
	    bool ReadFromFile(const std::string& fileName = "config.txt");
	    int num_steps() const {
			return num_steps_;
		}
	    double time_step() const {
			return dt_;
		}
	    int save_each_step() const {
			return save_each_step_;
		}
	    int time_log_start() const {
			return time_log_start_;
		}
	    int time_log_stop() const {
			return time_log_stop_;
		}
	    double heat_exchange_coeff() const {
			return alpha_;
		}
	    double ambient_temperature() const {
			return Tamb_;
		}
            double initial_temperature() const {
                return solid_material_.initial_temperature();
            }
            SolidMaterial& get_material() {
                return solid_material_;
            }
	private:
            int recognize_solid_model(const string&);
            int recognize_enthalpy_model(const string& model);
	    SolidMaterial solid_material_;
	    double dt_;
	    int num_steps_;
        int save_each_step_;
        int time_log_start_, time_log_stop_;
        double alpha_;
        double Tamb_;
};

#endif

