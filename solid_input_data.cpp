#include <cassert>
#include <map>
#include <Grid/material.h>
#include "solid_input_data.hpp"
#include "extended_input.hpp"

int SolidInputData::recognize_solid_model(const std::string& model) {
    std::string types[] = {"EQUILIBRIUM", "SCHEIL", "INDIRECT"};
    for(int i = 0; i < 3; ++i) {
        if(model == types[i])
            return i + 1;
    }
    throw std::string("Invalid solidification model");
}

int SolidInputData::recognize_enthalpy_model(const std::string& model) {
    std::string types[] = {"LINEAR_ENTHALPY_TYPE", "QUADRATIC_ENTHALPY_TYPE",
        "PIECEWISE_QUADRATIC_ENTHALPY_TYPE"};

    for(int i = 0; i < 3; ++i) {
        if(model == types[i])
            return i + 1;
    }
    throw std::string("Invalid enthalpy model");
}

bool SolidInputData::recognize_newton_bc(const MapConf& conf) {
    double alpha, Tamb;    
    std::string name = recognize_name(conf);
        
     if (!ReadValue(conf, "alpha", alpha)) {
        return false;
    }

    if (!ReadValue(conf, "Tamb", Tamb)) {
        return false;
    }
    
    NewtonBC bc(alpha, Tamb, name);
    boundary_conditions_coeff_.push_back(std::move(bc));
    
    return true;
}

std::string SolidInputData::recognize_name(const MapConf& conf) {
    std::string name;
    
    auto iter = conf.find("name");
    if(iter != conf.end()) {
        name = iter->second;
    } 
    
    if(name.empty())
        throw std::string("Invalid material name");
    
    return name;
}

bool SolidInputData::recognize_contact_bc(const MapConf& conf) {
    double kappa;
    
    std::string name = recognize_name(conf);
    
     if (!ReadValue(conf, "kappa", kappa)) {
        return false;
    }    
    
    return true;
}

bool SolidInputData::set_casting_properties(const MapConf& conf, SolidMaterial& solid_material) {
    double tmp;
    
    if (!ReadValue(conf, "lambdaS", tmp)) {
        return false;
    }

    solid_material.set_property("lambdaS", tmp);

    if (!ReadValue(conf, "lambdaL", tmp)) {
        return false;
    }

    solid_material.set_property("lambdaL", tmp);

    if (!ReadValue(conf, "rhoS", tmp)) {
        return false;
    }

    solid_material.set_property("rhoS", tmp);

    if (!ReadValue(conf, "rhoL", tmp)) {
        return false;
    }

    solid_material.set_property("rhoL", tmp);

    if (!ReadValue(conf, "cS", tmp)) {
        return false;
    }

    solid_material.set_property("cS", tmp);

    if (!ReadValue(conf, "cL", tmp)) {
        return false;
    }

    solid_material.set_property("cL", tmp);

    if (!ReadValue(conf, "Tp", tmp)) {
        return false;
    }

    solid_material.set_property("Tp", tmp);

    if (!ReadValue(conf, "Te", tmp)) {
        return false;
    }

    solid_material.set_property("Te", tmp);

    if (!ReadValue(conf, "Ts", tmp)) {
        return false;
    }

    solid_material.set_property("Ts", tmp);

    if (!ReadValue(conf, "Tl", tmp)) {
        return false;
    }

    solid_material.set_property("Tl", tmp);

    if (!ReadValue(conf, "L", tmp)) {
        return false;
    }

    solid_material.set_property("L", tmp); 

    if (!ReadValue(conf, "maxGrainSize", tmp)) {
        return false;
    }
    solid_material.set_property("maxGrainSize", tmp);

    if (!ReadValue(conf, "coeffBF", tmp)) {
        return false;
    }
    solid_material.set_property("coeffBF", tmp);
    
    return true;
}

bool SolidInputData::set_conductivity_properties(const MapConf& conf, SolidMaterial& solid_material) {
    double tmp;
    const double LARGE_NUMBER = 1e30;
    
    if (!ReadValue(conf, "lambdaS", tmp)) {
        return false;
    }

    solid_material.set_property("lambdaS", tmp);
    
    solid_material.set_property("lambdaL", 0.0);

    if (!ReadValue(conf, "rhoS", tmp)) {
        return false;
    }

    solid_material.set_property("rhoS", tmp);

    solid_material.set_property("rhoL", 0.0);

    if (!ReadValue(conf, "cS", tmp)) {
        return false;
    }

    solid_material.set_property("cS", tmp);

    solid_material.set_property("cL", 0.0); 
    solid_material.set_property("Tp", LARGE_NUMBER);
    solid_material.set_property("Te", LARGE_NUMBER);
    solid_material.set_property("Ts", LARGE_NUMBER);
    solid_material.set_property("Tl", LARGE_NUMBER);
    solid_material.set_property("L", 0.0); 
    solid_material.set_property("maxGrainSize", 0.0);
    solid_material.set_property("coeffBF", 0.0);
    
    return true;
}

bool SolidInputData::recognize_material(const MapConf& conf) {            
    solid_materials_.push_back(SolidMaterial());
    SolidMaterial& solid_material = solid_materials_[solid_materials_.size() - 1];
    
    solid_material.set_name(recognize_name(conf));
    solid_material.initialize_property_map();
    
    std::string model;
    if(!ReadValue(conf, "solid_model", model)) {
        return false;
    }
    
    solid_material.set_solidification_model(recognize_solid_model(model));

    if(!ReadValue(conf, "enthalpy_model", model)) {
        return false;
    }

    solid_material.set_enthalpy_model(recognize_enthalpy_model(model));

    int num;

    if(!ReadValue(conf, "isCasting", num)) {
        return false;
    }

    solid_material.set_casting(num);

    double tmp;

    if(!ReadValue(conf, "T0", tmp)) {
        return false;
    }

    solid_material.set_property("T0", tmp);

    if (solid_material.is_casting()) {
        set_casting_properties(conf, solid_material);
    } else {
        set_conductivity_properties(conf, solid_material);
    }

    solid_material.update_k();
    
    return true;
}

bool SolidInputData::ReadFromFile(const std::string& fileName) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::map<std::string, std::string> conf;
    std::vector<MapConf> groupConfs;
    jaz::ConfigFile cf;

    if (!read_config_file(fileName, conf, groupConfs)) {
        if(rank == 0)
            std::cerr << "Config file error ";// << cf.error() << std::endl;
        return false;
    }
    
    InputData::Initialize(conf, cf);

    if(!ReadValue(conf, "gridFile", inputFilenameGrid)) {
        return false;
    }

    if (!ReadValue(conf, "numberOfSteps", num_steps_)) {
        return false;
    }

    if (!ReadValue(conf, "dt", dt_)) {
        return false;
    }

    if (!ReadValue(conf, "saveEachStep", save_each_step_)) {
        return false;
    }
    if(save_each_step_ == 0) return false;

    if (!ReadValue(conf, "timeLogStart", time_log_start_)) {
        return false;
    }
    if (!ReadValue(conf, "timeLogStop", time_log_stop_)) {
        return false;
    }   
    
    for(auto group : groupConfs) {
        if(group.find("alpha") != group.end()) {
            if(!recognize_newton_bc(group)) {
                return false;
            }
        } else if(group.find("kappa") != group.end()) {
            if(!recognize_contact_bc(group))
                return false;
        } else {
            if(!recognize_material(group)) 
                return false;
        }
    }    

    return true;
}

void SolidInputData::find_maping_materials(const TALYFEMLIB::GRID& grid) {    
    std::map<std::string, int> reverse_map;
    
    for(unsigned i = 0; i < solid_materials_.size(); ++i) {
        reverse_map[solid_materials_[i].name()] = i;
    }
    
    //std::map<std::string, int> reverse_map_bc;
    
    for(unsigned i = 0; i < boundary_conditions_coeff_.size(); ++i) {
        reverse_map[boundary_conditions_coeff_[i].name()] = i;
    }
    
    for(auto iter = grid.material_desc_.begin(); iter != grid.material_desc_.end(); ++iter) {
        int real_ind = reverse_map[iter->second.tag];
        if(iter->second.type == TALYFEMLIB::TYPE_MATERIAL::VOLUME) {
            map_materials_[iter->first] = real_ind;
        } else {
            map_bc_[iter->first] = real_ind;
        }        
    }
}