#ifndef SOLID_INPUT_DATA_H
#define SOLID_INPUT_DATA_H

#include <string>
#include <memory>
#include <vector>
#include <InputData/InputData.h>
#include <Grid/grid_types/grid.h>
#include "solid_model.hpp"
#include "enthalpy_model.hpp"
#include "solid_material.hpp"
#include "newton_bc.hpp"
#include "extended_input.hpp"

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
            double initial_temperature(int ind = 0) const {
                int mat_ind = map_materials_.find(ind)->second;
                return solid_materials_[mat_ind].initial_temperature();
            }
            SolidMaterial& get_material(int ind = 0) {
                return solid_materials_[map_materials_[ind]];
            }            
            NewtonBC& get_bc(int ind = 0) {
                return boundary_conditions_coeff_[map_bc_[ind]];
            }  
            void find_maping_materials(const TALYFEMLIB::GRID& grid);
	private:
            int recognize_solid_model(const std::string&);
            int recognize_enthalpy_model(const std::string& model);
            bool recognize_newton_bc(const MapConf& conf);
            bool recognize_contact_bc(const MapConf& conf);
            bool recognize_material(const MapConf& grup);
            std::string recognize_name(const MapConf& map);
            bool set_casting_properties(const MapConf& conf, SolidMaterial& solid_material);
            bool set_conductivity_properties(const MapConf& conf, SolidMaterial& solid_material);            
	    std::vector<SolidMaterial> solid_materials_;
            std::vector<NewtonBC> boundary_conditions_coeff_;
            std::map<int, int> map_materials_;
            std::map<int, int> map_bc_;            
	    double dt_;
	    int num_steps_;
            int save_each_step_;
            int time_log_start_, time_log_stop_;
            double alpha_;
            double Tamb_;
};

#endif

