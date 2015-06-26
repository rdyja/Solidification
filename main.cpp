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
    //int ne = data.pGrid->n_elements();

    int ne = data.p_grid_->n_elements();

    TALYFEMLIB::FEMElm fe;
    for (int i = 0; i < ne; ++i) {
        fe.refill(data.p_grid_, i);
        //fe.setRelativeOrder(inputData.orderOfBF);

        int mat_ind = fe.pElm->mat_ind();
        const SolidMaterial& mat = inputData.get_material(mat_ind);

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
//	std::vector<PetscInt> pbc_indices_master = {
//			1, 2, 5, 6, 20
//			5, 6, 20
//	};
//	std::vector<PetscInt> pbc_indices_slave = {
//			8, 11, 12, 15, 24
//			12, 15, 24
//	};

    // set up periodic boundary object
    ContactBounds *pcb = new ContactBounds();
//    pcb->LoadContactBounds(pGrid, pbc_indices_master, pbc_indices_slave);
    pcb->ImportContactBounds(pGrid);
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
			ContactBounds *pcb = createContactBounds(inputData, pGrid);

			SolidEquation solidEq(&inputData, pcb);
			SolidGridField data(inputData);

			data.redimGrid(pGrid);
			data.redimNodeData();

			int nOfDofPerNode = 1;	//< number of degree of freedom per node
			solidEq.redimSolver (pGrid, nOfDofPerNode);
			solidEq.setData( &data );
			if (pcb)
				solidEq.InitializeContactBC(pcb);

			if (!inputData.ifBoxGrid && pcb)
				data.SetBndrIndicator(*pcb);
//			pGrid->PrintElmSurfaceIndicator();

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
