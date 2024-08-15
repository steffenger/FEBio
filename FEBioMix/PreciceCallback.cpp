#include "PreciceCallback.h"
#include <FECore/log.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FETimeStepController.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include "FESolutesMaterialPoint.h"
#include "FEMultiphasic.h"
#include <FECore/FELoadCurve.h>
#include "HelperProteinPosition.h"
#include <utility>


//change to receiving mesh
void PreciceCallback::ReadBoundaryConditionData (FEModel *fem, int index, const string otherMesh, const std::string dataName){
    /*
    const int otherMeshSize = precice.getMeshVertexSize(otherMesh);

    // Now finally get information about the mesh vertices.
    // First allocate memory for the IDs and coordinates
    std::vector<double> otherCoordinates(otherMeshSize * dim);
    std::vector<VertexID> otherVertexIDs(otherMeshSize);
    // ... and afterwards ask preCICE to fill the vectors
    precice.getMeshVertexIDsAndCoordinates(otherMesh,
                                           otherVertexIDs,
                                           otherCoordinates);

    // continue with time loop and write data directly to the mesh using
    // the received ids, which correspond to the vertices
    const int dataDim = precice.getDataDimensions(otherMesh, "OtherData");
    std::vector<double> data(dataDim * otherMeshSize);
    precice.writeData(otherMesh, "OtherData", otherVertexIDs, data);
    */
     const int dim = fem->participant->getMeshDimensions(otherMesh);
     const int otherMeshSize = fem->participant->getMeshVertexSize(otherMesh);
     
     std::vector<double> otherCoordinates(otherMeshSize * dim);
     std::vector<precice::VertexID> otherVertexIDs(otherMeshSize);
     // ... and afterwards ask preCICE to fill the vectors
     fem->participant->getMeshVertexIDsAndCoordinates(otherMesh,
                                           otherVertexIDs,
                                           otherCoordinates);

     const int dataDim = fem->participant->getDataDimensions(otherMesh, dataName);
     std::vector<double> data(dataDim * otherMeshSize);

     //std::vector<double> data(this->numberOfVertices);
                        double preciceDt = fem->participant->getMaxTimeStepSize();
                        double dt = min(preciceDt, fem->GetCurrentStep()->m_dt);
     fem->participant->readData("WholeBodyMesh", dataName, otherVertexIDs, dt, data);

     FELoadCurve* plc = dynamic_cast<FELoadCurve*>(fem->GetLoadController(index));
     double time =  fem->GetTime().currentTime;
     plc->Add(time, data[0]);
}

void PreciceCallback::WriteBoundaryConditionData (FEModel *fem, int index, const string dataSet, const std::string dataName) {

        std::vector<double> data(this->numberOfVertices);
        double f;
        double surface_jn;
        //probably need helper for outflux
        FEMesh &mesh = fem->GetMesh();
        FESurface *SurfaceOut = mesh.FindSurface(dataSet);
        FEElementSet* elementSetOutflow = mesh.FindElementSet(dataSet);
        for (int i = 0; i < SurfaceOut->Elements(); i++) {
            FEElement &element = SurfaceOut->Element(i);
            vec3d w(0,0,0);
            double c = 0;
            int nint = element.GaussPoints();
            for (int n=0; n<nint; ++n) {
            FESurfaceMaterialPoint *mp = dynamic_cast<FESurfaceMaterialPoint*>(element.GetMaterialPoint(n));
            //FESurfaceMaterialPoint *mp = element.GetSurfaceMaterialPoint(n);
            FEMaterialPoint *pt = element.GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(pt->ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(pt->ExtractData<FESolutesMaterialPoint>());
            w += pb.m_w;
            c += ps.m_ca[index];
            }
        w /= nint;
        c /= nint;
        
        // evaluate desired natural solute flux
        FESurfaceMaterialPoint *mp = dynamic_cast<FESurfaceMaterialPoint*>(element.GetMaterialPoint(1));
        vec3d dxt = mp->dxr ^ mp->dxs;
        double jn = c*(w*dxt);
        //if (flux->m_bshellb) jn = -jn;
        double dt = fem->GetCurrentStep()->m_dt;
        // molar flow rate
        f = jn* dt;
        surface_jn += f;
        }
        //double H_i = dof_a.shape;
        //fa[0] = H_i * f;
        std::cout << "surface_jn" << surface_jn << std::endl;
        data[0] = surface_jn;
        
        fem->participant->writeData(MESH_NAME, dataName, this->vertexIDs, data);
        
}

void PreciceCallback::UpdateCouplingData (FEModel *fem) {
       //create roadrunner template for every material point
       FEMesh &mesh = fem->GetMesh();
       FEElementSet* elementSetOutflow = mesh.FindElementSet("outflow");
       if (!elementSetOutflow) {
        feLogError((std::string("ElementSet not found")).c_str());
        throw FEException("ElementSet not found");
       }

       FEElementSet* elementSetInflow = mesh.FindElementSet("inflow");
       if (!elementSetInflow) {
        feLogError((std::string("ElementSet not found")).c_str());
        throw FEException("ElementSet not found");
       }



      FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
      for (int i = 0; i < elementSet->Elements(); i++) {
    	    	    FEElement &element = elementSet->Element(i);
    	    	    for (int j = 0; j < element.GaussPoints(); j++) {
    	    	    	FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);

						FESolutesMaterialPoint &ps = *(materialPoint->ExtractData<FESolutesMaterialPoint>());

		        //get Volume for microsimulation
                        FEMesh &m = fem->GetMesh();
                        double V = m.CurrentElementVolume(element);
                        //remove this for Center of element coupling
                        double V_gauss = V; //(element.GaussPoints());
						ps.volume = V_gauss;
                        
                        //double phiw = m_pMat->Porosity(mp);
                        

			//use helper function to get normalized position
                        double norm_position;
                        norm_position = get_normalized_position(*materialPoint, elementSetInflow, elementSetOutflow);
                        ps.norm_position = norm_position;


					}
	}
}



//try template function to simplify read and write function
template <typename T>
void PreciceCallback::ReadScalarDataTemplate(FEModel *fem, T FESolutesMaterialPoint::*member, const std::string dataName) {
// Read data from precice
    	    std::vector<double> data(this->numberOfVertices);
			double preciceDt = fem->participant->getMaxTimeStepSize();
			double dt = min(preciceDt, fem->GetCurrentStep()->m_dt);
    	    fem->participant->readData("FEBioMesh", dataName, this->vertexIDs, dt, data);

    	    // Write data to febio
    	    int counter = 0;
    	    FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
    	    for (int i = 0; i < elementSet->Elements(); i++) {
    	    	    FEElement &element = elementSet->Element(i);
                    double tmp_coe = data[counter];
    	    	    for (int j = 0; j < element.GaussPoints(); j++) {
    	    	    	    FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
			auto solute = materialPoint->ExtractData<FESolutesMaterialPoint>();
			if (solute == nullptr) {
    	    			throw FEException("MaterialPoint is not an instance of FESolutesMaterialPoint");
				}
			solute->*member = tmp_coe/element.GaussPoints();
    	    	    	}
                counter++;
    	    	}

}

template <typename T>
void PreciceCallback::ReadVectorDataTemplate(FEModel *fem, std::vector<T> FESolutesMaterialPoint::*member, int index, const std::string dataName) {
// Read data from precice
    	    std::vector<double> data(this->numberOfVertices);
			double preciceDt = fem->participant->getMaxTimeStepSize();
			double dt = min(preciceDt, fem->GetCurrentStep()->m_dt);
    	    fem->participant->readData(MESH_NAME, dataName, this->vertexIDs, dt, data);

    	    // Write data to febio
    	    int counter = 0;
    	    FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
    	    for (int i = 0; i < elementSet->Elements(); i++) {
    	    	    FEElement &element = elementSet->Element(i);
                    double tmp_coe = data[counter];
    	    	    for (int j = 0; j < element.GaussPoints(); j++) {
    	    	    	    FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
			auto solute = materialPoint->ExtractData<FESolutesMaterialPoint>();
			if (solute == nullptr) {
    	    			throw FEException("MaterialPoint is not an instance of FESolutesMaterialPoint");
				}
			(solute->*member)[index] = tmp_coe/element.GaussPoints();
    	    	    	}
            counter++;
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
            double tmp_coe = 0;
            for (int j = 0; j < element.GaussPoints(); j++) {
                FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
		auto solute = materialPoint->ExtractData<FESolutesMaterialPoint>();
			if (solute == nullptr) {
    	    			throw FEException("MaterialPoint is not an instance of FESolutesMaterialPoint");
				}
                tmp_coe += solute->*member;
            }
	data[counter] = tmp_coe/element.GaussPoints();
                counter++;
        }

	// Write data to precice
    //this->precice->writeData(MESH_NAME, dataName, this->vertexIDs, data);
    fem->participant->writeData(MESH_NAME, dataName, this->vertexIDs, data);

}

template <typename T>
void PreciceCallback::WriteVectorDataTemplate(FEModel *fem, std::vector<T> FESolutesMaterialPoint::*member, int index,  const std::string dataName) {
        // Read data from febio
	    std::vector<double> data(this->numberOfVertices);
        int counter = 0;
        FEElementSet *elementSet = fem->GetMesh().FindElementSet("CouplingDomain");//ELEMENT_SET);
        for (int i = 0; i < elementSet->Elements(); i++) {
            FEElement &element = elementSet->Element(i);
            double tmp_coe = 0;
            for (int j = 0; j < element.GaussPoints(); j++) {
                FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
		auto solute = materialPoint->ExtractData<FESolutesMaterialPoint>();
			if (solute == nullptr) {
    	    			throw FEException("MaterialPoint is not an instance of FESolutesMaterialPoint");
				}
                tmp_coe += (solute->*member)[index];
            }
		data[counter] = tmp_coe/element.GaussPoints();
                counter++;
        }

	// Write data to precice
    //this->precice->writeData(MESH_NAME, dataName, this->vertexIDs, data);
    fem->participant->writeData(MESH_NAME, dataName, this->vertexIDs, data);
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
    	for (int i = 0; i < elementSet->Elements(); i++) {
            FEElement &element = elementSet->Element(i);
            vec3d com = vec3d(0, 0, 0);
                for (int j = 0; j < element.GaussPoints(); j++) {
    	        	// iterate over all materialpoints and add initial position to vector
    	        	FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
                        com.x += materialPoint->m_r0.x/element.GaussPoints();
                        com.y += materialPoint->m_r0.y/element.GaussPoints();
                        com.z += materialPoint->m_r0.z/element.GaussPoints();
    	        	//vec3d coord = materialPoint->m_r0; 
    	        	/*vertexPositions.push_back(coord.x);
    	        	vertexPositions.push_back(coord.y);
    	        	vertexPositions.push_back(coord.z);
    	        	counter++;*/
    	    	}
            std::cout <<"x " << com.x << std::endl;
            std::cout <<"y " << com.y << std::endl;
            std::cout <<"z " << com.z << std::endl;
            vertexPositions.push_back(com.x);
            vertexPositions.push_back(com.y);
            vertexPositions.push_back(com.z);
            counter++;
    	}
    	return std::pair<int, vector<double>> (counter, vertexPositions);
}

// Initialize the precice adapter
void PreciceCallback::Init(FEModel *fem) {
    	feLogInfo("PreciceCallback::Init");
        
        //PARTICIPANT_NAME = "FEBio";
        //ELEMENT_SET = "CouplingDomain";
        /*
        //Define strings
        PARTICIPANT_NAME = "FEBio";
        const std::string MESH_NAME = "FEBioMesh";
        const std::string PARTICIPANT_NAME = "CouplingDomain";
        const char *config = "./precice-config.xml";
        */
	// Get config path from envrironment
	const char *config = getenv("BFP_CONFIG");
	if (!config) {
		config = "./precice-config.xml";
	}

    	// initialize precice
    	//this->precice = new precice::Participant(PARTICIPANT_NAME, config, 0, 1);
        fem->participant = new precice::Participant(PARTICIPANT_NAME, config, 0, 1);
    	this->dimensions = fem->participant->getMeshDimensions(MESH_NAME);

    	// Get material point positions
    	FEMesh &femMesh = fem->GetMesh();
    	std::pair<int, vector<double>> vertexInfo = this->getRelevantMaterialPoints(fem, ELEMENT_SET);
    	this->numberOfVertices = vertexInfo.first;
    	vector<double> vertexPositions = vertexInfo.second;

    	// Initialize precice mesh
    	this->vertexIDs.resize(this->numberOfVertices);
    	fem->participant->setMeshVertices(MESH_NAME, vertexPositions, this->vertexIDs);

    	// Finish initializing precice
    	//precice->initialize();
        fem->participant->initialize();
        
		//WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::volume, WRITE_DATA3);
    	feLogInfo("Finished PreciceCallback::Init");
}

bool PreciceCallback::Execute(FEModel &fem, int nreason) {
    	feLogInfo("PreciceCallback::Execute");

    	if (nreason == CB_INIT) {
    	    	this->Init(&fem);
            //communicate the initialization variables here?
			/*UpdateCouplingData(&fem);
			WriteVectorDataTemplate(&fem, &FESolutesMaterialPoint::m_ca, 0, WRITE_DATA);
	        WriteVectorDataTemplate(&fem, &FESolutesMaterialPoint::m_ca, 1, WRITE_DATAP);
            WriteScalarDataTemplate(&fem, &FESolutesMaterialPoint::volume, WRITE_DATAV);
			WriteScalarDataTemplate(&fem, &FESolutesMaterialPoint::norm_position, WRITE_DATANP);*/


    	} else if (nreason == CB_UPDATE_TIME) {
    	    	if (fem.participant->requiresWritingCheckpoint()) {
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
				double preciceDt = fem.participant->getMaxTimeStepSize();
    	    	double dt = min(preciceDt, fem.GetCurrentStep()->m_dt);
    	    	feLogInfo("Current Simulation Time %f\n", fem.GetTime().currentTime);
    	    	feLogInfo("Timestep %f\n", dt);
    	    	fem.GetCurrentStep()->m_dt = dt;
    	} else if (nreason == CB_MAJOR_ITERS) {
    	    	if (fem.participant->isCouplingOngoing()) {
    	    	    	// Read and write precice data
    	    	    	this->ReadData(&fem);
    	    	    	this->WriteData(&fem);
                        double preciceDt = fem.participant->getMaxTimeStepSize();
                        double dt = min(preciceDt, fem.GetCurrentStep()->m_dt);
		        //double dt = this->precice->getMaxTimeStepSize();
    	    	    	fem.participant->advance(dt);
    	    	    	if (fem.participant->requiresReadingCheckpoint()) {
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
    	    	fem.participant->finalize();
    	    	//delete precice;
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
    //void function like overwrite internal values
	UpdateCouplingData(fem);
	WriteVectorDataTemplate(fem, &FESolutesMaterialPoint::m_ca, 0, WRITE_DATA);
	WriteVectorDataTemplate(fem, &FESolutesMaterialPoint::m_ca, 1, WRITE_DATAP);
        WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::volume, WRITE_DATAV);
	//WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::f_fluid, WRITE_DATAP);
	//WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::f_tissue, WRITE_DATAP);
	WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::norm_position, WRITE_DATANP);



    feLogInfo("Finished PreciceCallback::WriteData");
}
