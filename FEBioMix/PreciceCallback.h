#include <precice/SolverInterface.hpp>
#include <FECore/sdk.h>
#include <FECore/FECallBack.h>
#include <FECore/Callback.h>
#include <FECore/FENode.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEAnalysis.h>
#include <FECore/DumpMemStream.h>
#include <FEBioMix/FESolutesMaterialPoint.h>
#include <utility>

#define PARTICIPANT_NAME "FEBio"
#define ELEMENT_SET "CouplingDomain"
#define MESH_NAME "FEBioMesh"
#define READ_DATA "S_ext'"
#define READ_DATA2 "P_ext'"
#define READ_DATA3 "uec(PEX, P_ext)"
#define READ_DATA4 "uec(SIM, S_ext)"


#define WRITE_DATA "S_ext"
#define WRITE_DATAP "P_ext"
#define WRITE_DATA3 "volume"
//#define WRITE_DATA3 "vol_point"
//#define WRITE_DATA4 "f_fluid"
//#define WRITE_DATA5 "f_tissue"
//#define WRITE_DATA6 "norm_position"

//add other variables

class PreciceCallback : public FECallBack {
public:
   	PreciceCallback(FEModel *pfem) : FECallBack(pfem, CB_INIT | CB_UPDATE_TIME | CB_MAJOR_ITERS), dmp(*pfem) {}
    	void ReadData(FEModel *fem);
    	void WriteData(FEModel *fem);
    	void Init(FEModel *fem);
    	bool Execute(FEModel &fem, int nreason);
    	std::pair<int, vector<double>> getRelevantMaterialPoints(FEModel *fem, const std::string &elementName);
		template <typename T> void ReadScalarDataTemplate(FEModel *fem, T FESolutesMaterialPoint::*member, const std::string MACRO);
		template <typename T> void ReadVectorDataTemplate(FEModel *fem, std::vector<T> FESolutesMaterialPoint::*member, int index, const std::string MACRO);
		template <typename T> void WriteScalarDataTemplate(FEModel *fem, T FESolutesMaterialPoint::*member, const std::string MACRO);
		template <typename T> void WriteVectorDataTemplate(FEModel *fem, std::vector<T> FESolutesMaterialPoint::*member, int index, const std::string MACRO);

protected:
    	precice::SolverInterface *precice = NULL;
    	int dimensions; 		// precice dimensions
    	double dt; 			// solver timestep size
    	double precice_dt; 		// precice timestep size 
    	int numberOfVerticies; 		// number of vertices of muscle 
    	int meshID;			// precice ID of the muscle mesh
    	std::vector<int> vertexIDs;	// vertex IDs of the muscle mesh

	// Checkpoints used for implicit coupling
    	const std::string &coric = precice::constants::actionReadIterationCheckpoint();
    	const std::string &cowic = precice::constants::actionWriteIterationCheckpoint();
    	FEAnalysis *checkPointStep;
    	DumpMemStream dmp;
    	double checkpoint_time = 0;
    	FETimeStepController *checkpointTimeStepController = nullptr;
};


