#include <cassert>
#include <map>
#include "solid_input_data.hpp"

int SolidInputData::recognize_solid_model(const string& model) {
    string types[] = {"EQUILIBRIUM", "SCHEIL", "INDIRECT"};
    for(int i = 0; i < 3; ++i) {
        if(model == types[i])
            return i + 1;
    }
    throw std::string("Invalid solidification model");
}

int SolidInputData::recognize_enthalpy_model(const string& model) {
    string types[] = {"LINEAR_ENTHALPY_TYPE", "QUADRATIC_ENTHALPY_TYPE",
        "PIECEWISE_QUADRATIC_ENTHALPY_TYPE"};

    for(int i = 0; i < 3; ++i) {
        if(model == types[i])
            return i + 1;
    }
    throw std::string("Invalid enthalpy model");
}

bool SolidInputData::ReadFromFile(const std::string& fileName) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::map<std::string, std::string> conf;
    jaz::ConfigFile cf;

    if (!cf.read(fileName, conf)) {
        if(rank == 0)
            std::cerr << "Config file error " << cf.error() << std::endl;
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

    std::string model;
    if(!ReadValue(conf, "solid_model", model)) {
        return false;
    }

    solid_material_.set_solidification_model(recognize_solid_model(model));

    if(!ReadValue(conf, "enthalpy_model", model)) {
        return false;
    }

    solid_material_.set_enthalpy_model(recognize_enthalpy_model(model));

    int num;

    if(!ReadValue(conf, "isCasting", num)) {
        return false;
    }

    solid_material_.set_casting(num);

    double tmp;

    if(!ReadValue(conf, "T0", tmp)) {
        return false;
    }

    solid_material_.set_property("T0", tmp);

    if (!ReadValue(conf, "lambdaS", tmp)) {
        return false;
    }

    solid_material_.set_property("lambdaS", tmp);

    if (!ReadValue(conf, "lambdaL", tmp)) {
        return false;
    }

    solid_material_.set_property("lambdaL", tmp);

    if (!ReadValue(conf, "rhoS", tmp)) {
        return false;
    }

    solid_material_.set_property("rhoS", tmp);

    if (!ReadValue(conf, "rhoL", tmp)) {
        return false;
    }

    solid_material_.set_property("rhoL", tmp);

    if (!ReadValue(conf, "cS", tmp)) {
        return false;
    }

    solid_material_.set_property("cS", tmp);

    if (!ReadValue(conf, "cL", tmp)) {
        return false;
    }

    solid_material_.set_property("cL", tmp);

    if (!ReadValue(conf, "Tp", tmp)) {
        return false;
    }

    solid_material_.set_property("Tp", tmp);

    if (!ReadValue(conf, "Te", tmp)) {
        return false;
    }

    solid_material_.set_property("Te", tmp);

    if (!ReadValue(conf, "Ts", tmp)) {
        return false;
    }

    solid_material_.set_property("Ts", tmp);

    if (!ReadValue(conf, "Tl", tmp)) {
        return false;
    }

    solid_material_.set_property("Tl", tmp);

    if (!ReadValue(conf, "L", tmp)) {
        return false;
    }

    solid_material_.set_property("L", tmp);

    if (!ReadValue(conf, "alpha", alpha_)) {
        return false;
    }

    if (!ReadValue(conf, "Tamb", Tamb_)) {
        return false;
    }

    if (!ReadValue(conf, "maxGrainSize", tmp)) {
        return false;
    }
    solid_material_.set_property("maxGrainSize", tmp);

    if (!ReadValue(conf, "coeffBF", tmp)) {
        return false;
    }
    solid_material_.set_property("coeffBF", tmp);

    return true;
}
