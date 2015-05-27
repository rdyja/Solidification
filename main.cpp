#include <TALYFEMLibInclude.h>
#include "solid_equation.hpp"
#include "solid_input_data.hpp"
#include "solid_grid_field.hpp"
#include "contact_bounds.hpp"
#include <ExternalLibs/jaz/jaz/TimeLog.hpp>

namespace tf = TALYFEMLIB;

inline bool SetIC(SolidGridField& data, SolidInputData& idata) {
    data.SetIC(idata.nsd);
    return false;
}

void compute_cooling_velocity_before_liquidus(SolidInputData& inputData,
        SolidGridField& data, double dTime)
{
    int ne = data.p_grid_->n_elements();
    TALYFEMLIB::FEMElm fe;
    for (int i = 0; i < ne; ++i) {
        fe.refill(data.p_grid_, i);
        //fe.setRelativeOrder(inputData.orderOfBF);

        const SolidMaterial& mat = inputData.get_material();

        if (mat.is_casting())
            for (int j = 0; j < fe.pElm->n_nodes(); ++j) {
                int globalNumNode = fe.pElm->ElemToLocalNodeID(j);
                if (data.Node(globalNumNode).get_prev_temp() > mat.liquidus_temperature()) {
                    double vel = fabs(mat.initial_temperature() - data.Node(globalNumNode).get_prev_temp())/dTime;
                    data.Node(globalNumNode).set_velocity(vel);
                }
            }
        else
            for (int j = 0; j < fe.pElm->n_nodes(); ++j) {
                int globalNumNode = fe.pElm->ElemToLocalNodeID(j);
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
    std::string resultFileNamePrefix = "data_final", extension = ".plt";
    std::stringstream sfln, st;

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

void readConfigFile(SolidInputData& inputData, GRID *& pGrid) {
    if(!inputData.ReadFromFile()) {
        throw(std::string("Problem with reading input data!"));
    }

    if(!inputData.CheckInputData()) {
        throw(std::string("Problem with input data, check the config file"));
    }

    CreateGrid(pGrid, &inputData);
    inputData.find_maping_materials(*pGrid);
}

ContactBounds* createContactBounds(SolidInputData& inputData, GRID *& pGrid) {

	// prepare data (normally should be read from gmsh file)
	std::vector<PetscInt> pbc_indices_master = {
			 3,  4,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 29, 30, 31, 32,
			45, 46, 50, 52, 51, 53, 56, 55, 54, 57, 60, 58, 59, 93, 94, 95
	};
	std::vector<PetscInt> pbc_indices_slave = {
			25, 26, 27, 28, 17, 18, 19, 20, 21, 22, 23, 24, 18, 20, 23, 24,
			91, 92, 65, 67, 66, 68, 71, 70, 69, 64, 74, 72, 73, 61, 62, 63
	};
//	std::vector<int> ind = { 1, 2, 3, 4, 8, 9 };
//	int indno = 9; // should be number of tags

    // set up periodic boundary object
    ContactBounds *pcb = new ContactBounds();
    pcb->LoadContactBounds(pGrid, pbc_indices_master, pbc_indices_slave);
    return pcb;
}

int main(int argc, char** argv) {
    char help[] = "Solves a solidification problem!";
    int returnCode = 0;

    PetscInitialize(&argc,&argv,(char *)0, help);

    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        SolidInputData inputData;
        GRID *pGrid = nullptr;
		{
			readConfigFile(inputData, pGrid);
			ContactBounds *pcb = NULL;//createContactBounds(inputData, pGrid);

			SolidEquation solidEq(&inputData, pcb);
			SolidGridField data(inputData);

			data.redimGrid(pGrid);
			data.redimNodeData();
			if (!inputData.ifBoxGrid)
				data.SetBndrIndicator(1.0, 1.0, 1.0);

			int nOfDofPerNode = 1;	//< number of degree of freedom per node
			solidEq.redimSolver (pGrid, nOfDofPerNode);
			solidEq.setData( &data );

			performCalculation(inputData, data, solidEq, rank);

			delete pcb;
		}
		DestroyGrid(pGrid); // to make sure that grid is destroyed after gridfield
    }
    catch(const std::string& s) {
        std::cerr << s << std::endl;
		returnCode = -1;
    }
    catch(std::bad_alloc& e) {
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
