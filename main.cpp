#include <TALYFEMLibInclude.h>
#include "solid_equation.hpp"
#include "solid_input_data.hpp"
#include "solid_grid_field.hpp"
#include <ExternalLibs/jaz/jaz/TimeLog.hpp>

namespace tf = TALYFEMLIB;

inline bool SetIC(SolidGridField& data, SolidInputData& idata) {
    data.SetIC(idata.nsd);
    return false;
}

void compute_cooling_velocity_before_liquidus(SolidInputData& inputData,
        SolidGridField& data, double dTime)
{
    int ne = data.pGrid->n_nodes();
    TALYFEMLIB::FEMElm fe;
    for (int i = 0; i < ne; ++i) {
        fe.refill(data.pGrid, i);
        //fe.setRelativeOrder(inputData.orderOfBF);

        const SolidMaterial& mat = inputData.get_material();

        if (mat.is_casting())
            for (int j = 0; j < fe.pElm->nodeno; ++j) {
                int globalNumNode = fe.pElm->glbNodeID(j);
                if (data.Node(globalNumNode).get_prev_temp() > mat.liquidus_temperature()) {
                    double vel = fabs(mat.initial_temperature() - data.Node(globalNumNode).get_prev_temp())/dTime;
                    data.Node(globalNumNode).set_velocity(vel);
                }
            }
        else
            for (int j = 0; j < fe.pElm->nodeno; ++j) {
                int globalNumNode = fe.pElm->glbNodeID(j);
                data.Node(globalNumNode).set_velocity(0.0);
            }

    }
}

void performCalculation(SolidInputData& inputData, SolidGridField& data,
    SolidEquation& solidEq, int rank) {
    double dt = inputData.time_step();
    double finalStep = inputData.num_steps();
    int logStart = inputData.time_log_start();
    int logStop = inputData.time_log_stop();
    int save = inputData.save_each_step() == 0 ? 1 : inputData.save_each_step();
    double	t = 0.0;
    string resultFileNamePrefix = "data_final", extension = ".plt";
    stringstream sfln, st;

    SetIC(data, inputData);
    jaz::TimeLog timeLogger("Total_solve", finalStep);

    for(int i = 0; i < finalStep; ++i) {
        if (i == logStart) timeLogger.start();
        t += dt;
        compute_cooling_velocity_before_liquidus(inputData, data, t);
        solidEq.Solve(t, dt);
        data.UpdateDataStructures();

        if (i >= logStart && i < logStop) {
            st << t;
            timeLogger.add(st.str());
            st.str(std::string());
        }

        PrintStatus("(", t, ") time stepping completed! ", rank);
        //cout << sfln.str().c_str() << endl;
        if (i%save == 0 && !((logStop - logStart) > 0)) {
            if (inputData.ifBoxGrid) {
                sfln << resultFileNamePrefix << i << extension;
                TALYFEMLIB::save_gf(&data, &inputData, sfln.str().c_str(), t);
            } else {
                sfln << resultFileNamePrefix << i << "_" << rank << extension;
                //data.printNodeData(sfln.str().c_str());                
                TALYFEMLIB::save_gf(&data, &inputData, sfln.str().c_str(), t);
            }
            sfln.str(std::string());
        }
        if (i == logStop - 1) timeLogger.stop();
    }
    
    if ((logStop - logStart) > 0 && (rank == 0)) {
        timeLogger.print(std::cout);
    }
}

void readGrid(SolidInputData& inputData, GRID *&pGrid) {
    /*if (inputData.ifBoxGrid) {
        CreateGrid(pGrid, &inputData);
    }    
    else {
        pGrid = new GRID;
        if (inputData.ifDD) {
            pGrid->loadFromALBERTAFileDD(inputData.inputFilenameGrid.c_str());
        } else {
            pGrid->loadFromALBERTAFile(inputData.inputFilenameGrid.c_str());
        }
    }*/
    CreateGrid(pGrid, &inputData);
}

void readConfigFile(SolidInputData& inputData, GRID *& pGrid) {
    if(!inputData.ReadFromFile()) {
        throw(std::string("Problem with reading input data!"));
    }

    if(!inputData.CheckInputData()) {
        throw(std::string("Problem with input data, check the config file"));
    }

    readGrid(inputData, pGrid);
}

int main(int argc, char** argv) {
    char help[] = "Solves a transient transfer problem!";
    int returnCode = 0;

    PetscInitialize(&argc,&argv,(char *)0, help);

    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        SolidInputData inputData;
        GRID *pGrid = nullptr;

        readConfigFile(inputData, pGrid);

	SolidEquation solidEq(&inputData);
	SolidGridField data(inputData);

        data.redimGrid(pGrid);
        data.redimNodeData(pGrid);
        if (!inputData.ifBoxGrid) data.SetBndrIndicator(1.0, 1.0, 1.0);

        int nOfDofPerNode = 1;	//< number of degree of freedom per node
	solidEq.redimSolver (pGrid, nOfDofPerNode);
	solidEq.setData( &data );

	performCalculation(inputData, data, solidEq, rank);
    }
    catch(const std::string& s) {
        std::cerr << s << std::endl;
	returnCode = -1;
    }
    catch(std::bad_alloc e) {
	std::cerr << "Problem with memory allocation " << e.what();
	returnCode = -1;
    }
    catch(...) {
        std::cerr << "Unknown error" << std::endl;
	returnCode = -1;
    }

    PetscFinalize();

    return returnCode;
}
