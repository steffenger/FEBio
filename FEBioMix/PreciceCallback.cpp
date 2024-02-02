#include "PreciceCallback.h"
#include <FECore/log.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FETimeStepController.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include "FESolutesMaterialPoint.h"
#include <utility>

//try template function to simplify read and write function
template <typename T>
void PreciceCallback::ReadScalarDataTemplate(FEModel *fem, T FESolutesMaterialPoint::*member, const std::string dataName) {
// Read data from precice
    	    std::vector<double> data(this->numberOfVertices);
			double preciceDt = precice->getMaxTimeStepSize();
			double dt = min(preciceDt, fem->GetCurrentStep()->m_dt);
    	    precice->readData(MESH_NAME, dataName, this->numberOfVertices, this->vertexIDs.data(), dt, data.data());

    	    // Write data to febio
    	    int counter = 0;
    	    FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
    	    for (int i = 0; i < elementSet->Elements(); i++) {
    	    	    FEElement &element = elementSet->Element(i);
    	    	    for (int j = 0; j < element.GaussPoints(); j++) {
    	    	    	    FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
			auto solute = materialPoint->ExtractData<FESolutesMaterialPoint>();
			if (solute == nullptr) {
    	    			throw FEException("MaterialPoint is not an instance of FESolutesMaterialPoint");
				}
			solute->*member = data[counter];
    	    	    	    counter++;
    	    	    	}
    	    	}

}

template <typename T>
void PreciceCallback::ReadVectorDataTemplate(FEModel *fem, std::vector<T> FESolutesMaterialPoint::*member, int index, const std::string dataName) {
// Read data from precice
    	    std::vector<double> data(this->numberOfVertices);
			double preciceDt = precice->getMaxTimeStepSize();
			double dt = min(preciceDt, fem->GetCurrentStep()->m_dt);
    	    precice->readData(MESH_NAME, dataName, this->numberOfVertices, this->vertexIDs.data(), dt, data.data());

    	    // Write data to febio
    	    int counter = 0;
    	    FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
    	    for (int i = 0; i < elementSet->Elements(); i++) {
    	    	    FEElement &element = elementSet->Element(i);
    	    	    for (int j = 0; j < element.GaussPoints(); j++) {
    	    	    	    FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
			auto solute = materialPoint->ExtractData<FESolutesMaterialPoint>();
			if (solute == nullptr) {
    	    			throw FEException("MaterialPoint is not an instance of FESolutesMaterialPoint");
				}
			(solute->*member)[index] = data[counter];
    	    	    	    counter++;
    	    	    	}
    	    	}

}

template <typename T>
void PreciceCallback::WriteScalarDataTemplate(FEModel *fem, T FESolutesMaterialPoint::*member, const std::string dataName) {
        // Read data from febio
	    std::vector<double> data(this->numberOfVertices);
        int counter = 0;
        FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
        for (int i = 0; i < elementSet->Elements(); i++) {
            FEElement &element = elementSet->Element(i);
            for (int j = 0; j < element.GaussPoints(); j++) {
                FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
		auto solute = materialPoint->ExtractData<FESolutesMaterialPoint>();
			if (solute == nullptr) {
    	    			throw FEException("MaterialPoint is not an instance of FESolutesMaterialPoint");
				}
		data[counter] = solute->*member;
                counter++;
            }
        }

	// Write data to precice
    this->precice->writeData(MESH_NAME, dataName, this->numberOfVertices, this->vertexIDs.data(), data.data());

}	

template <typename T>
void PreciceCallback::WriteVectorDataTemplate(FEModel *fem, std::vector<T> FESolutesMaterialPoint::*member, int index,  const std::string dataName) {
        // Read data from febio
	    std::vector<double> data(this->numberOfVertices);
        int counter = 0;
        FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
        for (int i = 0; i < elementSet->Elements(); i++) {
            FEElement &element = elementSet->Element(i);
            for (int j = 0; j < element.GaussPoints(); j++) {
                FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
		auto solute = materialPoint->ExtractData<FESolutesMaterialPoint>();
			if (solute == nullptr) {
    	    			throw FEException("MaterialPoint is not an instance of FESolutesMaterialPoint");
				}
		data[counter] = (solute->*member)[index];
                counter++;
            }
        }

	// Write data to precice
    this->precice->writeData(MESH_NAME, dataName, this->numberOfVertices, this->vertexIDs.data(), data.data());

}	

// Get number of material points and their initial position
std::pair<int, vector<double>> PreciceCallback::getRelevantMaterialPoints(FEModel *fem, const std::string &elementName) {
    	vector<double> vertexPositions;
    	int counter = 0;
    	FEElementSet *elementSet = fem->GetMesh().FindElementSet(elementName);
    	if (!elementSet) {
    	    	feLogError((elementName + std::string("ElementSet not found")).c_str());
    	    	throw FEException("ElementSet not found");
    	}
    	for (int i = 0; i < elementSet->Elements(); i++) { FEElement &element = elementSet->Element(i); for (int j = 0; j < element.GaussPoints(); j++) {
    	        	// iterate over all materialpoints and add initial position to vector
    	        	FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
    	        	vec3d coord = materialPoint->m_r0; 
    	        	vertexPositions.push_back(coord.x);
    	        	vertexPositions.push_back(coord.y);
    	        	vertexPositions.push_back(coord.z);
    	        	counter++;
    	    	}
    	}
    	return std::pair<int, vector<double>> (counter, vertexPositions);
}

// Initialize the precice adapter
void PreciceCallback::Init(FEModel *fem) {
    	feLogInfo("PreciceCallback::Init");

	// Get config path from envrironment
	const char *config = getenv("BFP_CONFIG");
	if (!config) {
		config = "./precice-config.xml";
	}

    	// initialize precice
    	this->precice = new precice::Participant(PARTICIPANT_NAME, config, 0, 1);
    	this->dimensions = this->precice->getMeshDimensions(MESH_NAME);

    	// Get material point positions
    	FEMesh &femMesh = fem->GetMesh();
    	std::pair<int, vector<double>> vertexInfo = this->getRelevantMaterialPoints(fem, ELEMENT_SET);
    	this->numberOfVertices = vertexInfo.first;
    	vector<double> vertexPositions = vertexInfo.second;

    	// Initialize precice mesh
    	this->vertexIDs.resize(this->numberOfVertices);
    	this->precice->setMeshVertices(MESH_NAME, this->numberOfVertices, vertexPositions.data(), this->vertexIDs.data());

    	// Finish initializing precice
    	precice->initialize();     

		//WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::volume, WRITE_DATA3);
    	feLogInfo("Finished PreciceCallback::Init");
}

bool PreciceCallback::Execute(FEModel &fem, int nreason) {
    	feLogInfo("PreciceCallback::Execute");

    	if (nreason == CB_INIT) {
    	    	this->Init(&fem);
            //communicate the initialization variables here?
            //WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::volume, WRITE_DATA3);
			//WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::norm_position, WRITE_DATA);


    	} else if (nreason == CB_UPDATE_TIME) {
    	    	if (this->precice->requiresWritingCheckpoint()) {
    	    	    	feLogInfo("CB_UPDATE_TIME - Saving Checkpoint\n");
    	    	    	// Save
    	    	    	// this uses dmp.open(true,true) which leads to the time controller not beeing serialized
    	    	    	// Also setting dmp.open(true,false) leads to segfault dont know why yet
    	    	    	// Switch time controller
    	    	    	/*delete this->checkpointTimeStepController;
    	    	    	this->checkpointTimeStepController = fem.GetCurrentStep()->m_timeController;
    	    	    	FETimeStepController *newTimeController = new FETimeStepController(&fem);

    	    	    	newTimeController->SetAnalysis(fem.GetCurrentStep());

						//warning this seems not to work for complex timestepping
    	    	    	newTimeController->CopyFrom(this->checkpointTimeStepController);
    	    	    	//fem.GetCurrentStep()->m_timeController = newTimeController;
						
    	    	    	this->checkpoint_time = fem.GetTime().currentTime;
    	    	    	this->dmp.clear();
    	    	    	fem.Serialize(this->dmp);*/
    	    	}
    	    	// advance timestep
				double preciceDt = precice->getMaxTimeStepSize();
    	    	double dt = min(preciceDt, fem.GetCurrentStep()->m_dt);
    	    	feLogInfo("Current Simulation Time %f\n", fem.GetTime().currentTime);
    	    	feLogInfo("Timestep %f\n", dt);
    	    	fem.GetCurrentStep()->m_dt = dt;
    	} else if (nreason == CB_MAJOR_ITERS) {
    	    	if (this->precice->isCouplingOngoing()) {
    	    	    	// Read and write precice data
    	    	    	this->ReadData(&fem);
    	    	    	this->WriteData(&fem);
						double dt = this->precice->getMaxTimeStepSize();
    	    	    	this->precice->advance(dt);
    	    	    	if (this->precice->requiresReadingCheckpoint()) {
    	    	    	    	feLogInfo("CB_MAJOR_ITERS - Restoring Checkpoint\n");
    	    	    	    	// Restore
    	    	    	    	// taken from FEAnalysis.cpp Line 475 ff
    	    	    	    	// restore the previous state
    	    	    	    	/*this->dmp.Open(false, true); // This does not restore the time controller only if bshallow is false
    	    	    	    	fem.Serialize(this->dmp);
    	    	    	    	FETimeStepController *newTimeController = new FETimeStepController(&fem);
    	    	    	    	newTimeController->SetAnalysis(fem.GetCurrentStep());
    	    	    	    	newTimeController->CopyFrom(this->checkpointTimeStepController);
    	    	    	    	fem.GetCurrentStep()->m_timeController = newTimeController;
    	    	    	    	fem.GetTime().currentTime = this->checkpoint_time;
    	    	    	    	fem.GetCurrentStep()->m_ntimesteps--; // Decrease number of steps because it gets increased right after this*/
    	    	    	}
    	    	}
    	} else if (nreason == CB_SOLVED) {
    	    	this->precice->finalize();
    	    	delete precice;
    	}
    	feLogInfo("Finished PreciceCallback::Execute");
    	return true;
}

// Read data from precice to febio
void PreciceCallback::ReadData(FEModel *fem) {
    feLogInfo("PreciceCallback::ReadData");
                
	ReadScalarDataTemplate(fem, &FESolutesMaterialPoint::m_sourceterm, READ_DATA);
	ReadScalarDataTemplate(fem, &FESolutesMaterialPoint::m_sourceterm2, READ_DATA2);
	ReadScalarDataTemplate(fem, &FESolutesMaterialPoint::m_tangent1, READ_DATA3);
	ReadScalarDataTemplate(fem, &FESolutesMaterialPoint::m_tangent2, READ_DATA4);

    feLogInfo("Finished PreciceCallback::ReadData");
}

// Write data from precice to febio
void PreciceCallback::WriteData(FEModel *fem) {
    feLogInfo("PreciceCallback::WriteData");

	WriteVectorDataTemplate(fem, &FESolutesMaterialPoint::m_ca, 0, WRITE_DATA);
	WriteVectorDataTemplate(fem, &FESolutesMaterialPoint::m_ca, 1, WRITE_DATAP);

    feLogInfo("Finished PreciceCallback::WriteData");
}
