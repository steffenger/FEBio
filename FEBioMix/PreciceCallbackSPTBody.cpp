#include "PreciceCallbackSPTBody.h"
#include <FECore/log.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FETimeStepController.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include "FESolutesMaterialPoint.h"
#include <FEBioMech/FEElasticMaterialPoint.h>
#include "HelperProteinPosition.h"
#include <FECore/FELoadCurve.h>
#include "FEMultiphasic.h"
//#include <FECore/FEModel.h>
#include <FECore/FESurface.h>
#include <utility>

//function to update data which is used in coupling but not updated in the solver itself
void PreciceCallbackSPTBody::UpdateCouplingData (FEModel *fem) {
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


      double Full_Volume;
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
                        Full_Volume += V;
                        
			//use helper function to get normalized position
                        double norm_position;
                        norm_position = get_normalized_position(*materialPoint, elementSetInflow, elementSetOutflow);
                        ps.norm_position = norm_position;


					}
	}

std::cout << "Full_Volume " << Full_Volume << std::endl;
}


// Initialize the precice adapter
void PreciceCallbackSPTBody::Init(FEModel *fem) {
    	feLogInfo("PreciceCallback::Init");
        //PARTICIPANT_NAME = "FEBio";
        //ELEMENT_SET = "CouplingDomain";
        
        //Define strings
        PARTICIPANT_NAME = "FEBio";
        ELEMENT_SET = "CouplingDomain";
        MESH_NAME = "FEBioMesh";
        /*
        const std::string PARTICIPANT_NAME = "FEBio";
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
        
        otherMesh = "WholeBodyMesh";
        
        //test direct mesh access
        otherMesh = "WholeBodyMesh";
        std::cout << "test1" << std::endl;
        // Get the spacial dimensionality of the mesh
        const int dim = fem->participant->getMeshDimensions(otherMesh);

        // Allocate and fill the  'boundingBox' according to the interested region
        // with the desired bounds, in our example we use the unit cube.
        // Assuming dim == 3, means that the bounding box has dim * 2 == 6 elements.
        std::vector<double> boundingBox {
        -0.0005, 0.0005, -0.0005, // minimum corner
        0.0005, 0.000, 0.0002 // maximum corner
        };

       // Define region of interest, where we want to obtain the direct access.
       // See also the API documentation of this function for further notes.
        fem->participant->setMeshAccessRegion(otherMesh, boundingBox);
        


    	// Finish initializing precice
    	fem->participant->initialize();     

		//WriteScalarDataTemplate(fem, &FESolutesMaterialPoint::volume, WRITE_DATA3);
    	feLogInfo("Finished PreciceCallback::Init");
}

bool PreciceCallbackSPTBody::Execute(FEModel &fem, int nreason) {
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
                        //fem.PushState();
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

                                //fem.PopState();
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
void PreciceCallbackSPTBody::ReadData(FEModel *fem) {
    feLogInfo("PreciceCallback::ReadData");
                
	ReadScalarDataTemplate(fem, &FESolutesMaterialPoint::m_sourceterm, READ_DATA);
	ReadScalarDataTemplate(fem, &FESolutesMaterialPoint::m_sourceterm2, READ_DATA2);
	ReadScalarDataTemplate(fem, &FESolutesMaterialPoint::m_tangent1, READ_DATA3);
	ReadScalarDataTemplate(fem, &FESolutesMaterialPoint::m_tangent2, READ_DATA4);
        
        //fem->participant->getMeshDimensions(otherMesh);
        /*std::cout << "test1" << std::endl;
        const int dim = fem->participant->getMeshDimensions(otherMesh);
        const int otherMeshSize = fem->participant->getMeshVertexSize(otherMesh);
        std::cout << "test2" << std::endl;
        std::vector<double> otherCoordinates(otherMeshSize * dim);
        std::vector<precice::VertexID> otherVertexIDs(otherMeshSize);
        // ... and afterwards ask preCICE to fill the vectors
        fem->participant->getMeshVertexIDsAndCoordinates(otherMesh,
                                           otherVertexIDs,
                                           otherCoordinates);
        std::cout << "test3" << std::endl;
        const int dataDim = fem->participant->getDataDimensions(otherMesh, "S_ext_inflow");
        std::vector<double> data(dataDim * otherMeshSize);
        std::cout << "test4" << std::endl;
        //std::vector<double> data(this->numberOfVertices);
                        double preciceDt = fem->participant->getMaxTimeStepSize();
                        double dt = min(preciceDt, fem->GetCurrentStep()->m_dt);
        fem->participant->readData("WholeBodyMesh", "S_ext_inflow", otherVertexIDs, dt, data);
        std::cout << "test5" << std::endl;
        //FELoadCurve* plc = dynamic_cast<FELoadCurve*>(fem->GetLoadController(2));
        double time =  fem->GetTime().currentTime; //+ fem->GetCurrentStep()->m_dt;
        std::cout << time << std::endl;
        if (!data.empty()) {
           double test = data[0];
           //plc->Add(time, test);
           std::cout << "First element: " << test << std::endl;
        } else {
           std::cerr << "Vector is empty." << std::endl;
        }
        fem->EvaluateLoadControllers(time);*/
        //std::cout << "test6" << data.at(0) << std::endl;
        //double test = data[0];
        //plc->Add(time, test);
        //ReadBoundaryConditionData(fem, 0, otherMesh, "S_ext_inflow");
        //ReadBoundaryConditionData(fem, 1, otherMesh, "P_ext_inflow");
        
    feLogInfo("Finished PreciceCallback::ReadData");
}

// Write data from precice to febio
void PreciceCallbackSPTBody::WriteData(FEModel *fem) {
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

/*BEGIN_FECORE_CLASS(PreciceCallbackSPTBody, FECallback)
     ADD_PARAMETER(otherMesh, "otherMesh")
 END_FECORE_CLASS();*/

MyLoadController::MyLoadController(FEModel *fem) : FELoadController(fem)
{   index = 0.0;
    m_val0 = 0.0;
    m_val1 = 1.0;
    m_duration = 1.0;
}

bool MyLoadController::Init()
{   m_duration = 0.0;
    std::cout << "test1" << std::endl;
    // add initialization here
    /*std::string PARTICIPANT_NAME = "FEBioBC";
    this->precice = new precice::Participant(PARTICIPANT_NAME, "../precice-config.xml", 0, 1);
    std::string otherMesh = "WholeBodyMesh";

    std::vector<double> boundingBox {
         -0.0005, 0.0005, -0.0005, // minimum corner
         0.0005, 0.000, 0.0002 // maximum corner
         };
    std::cout << "test1" << std::endl;
    this->precice->setMeshAccessRegion(otherMesh, boundingBox);
    // Get the spacial dimensionality of the mesh
    const int dim = precice->getMeshDimensions(otherMesh);


    // Finish initializing precice
    this->precice->initialize();
    // call base class*/
    return FELoadController::Init();
}

double MyLoadController::GetValue(double currentTime)
{   
    if (currentTime > 0) {
    double test;
    //if (precice->isCouplingOngoing()) {
    FEModel *fem = GetFEModel();
    if (fem->participant->isCouplingOngoing()) {
    std::string otherMesh = "WholeBodyMesh";
    const int otherMeshSize = fem->participant->getMeshVertexSize(otherMesh);
    const int dim = fem->participant->getMeshDimensions(otherMesh);
    const int dataDim = fem->participant->getDataDimensions(otherMesh, r_data);
    std::vector<double> data(dataDim * otherMeshSize);
    std::vector<double> otherCoordinates(otherMeshSize * dim);
    std::vector<precice::VertexID> otherVertexIDs(otherMeshSize);
    fem->participant->getMeshVertexIDsAndCoordinates(otherMesh,
                                            otherVertexIDs,
                                            otherCoordinates);
    std::cout << "test4 "<< dataDim << std::endl;
     //std::vector<double> data(this->numberOfVertices);
    double preciceDt = fem->participant->getMaxTimeStepSize();
    double dt = min(preciceDt, fem->GetCurrentStep()->m_dt);
    std::cout << "dt " << dt << std::endl;
    fem->participant->readData("WholeBodyMesh", r_data, otherVertexIDs, dt, data);
    
    std::vector<double> write_data(dataDim * otherMeshSize);
    
    FEMesh &mesh = fem->GetMesh();
    FEElementSet* elementSetOutflow = mesh.FindElementSet("outflow");
    write_data[0] = getOutflow(elementSetOutflow, index); //data[0]*0.9;
    std::cout << "index " << index << std::endl; 
    fem->participant->writeData("WholeBodyMesh", w_data, otherVertexIDs, write_data);
    std::cout << w_data << r_data << std::endl;
    std::cout << "write_data "<< write_data[0] << std::endl;
    std::cout << "data factor " << write_data[0]/data[0] << std::endl;
    std::cout << "data size" << data.size() << std::endl;

    //std::cout << "test4.1 "<< dataDim << std::endl;
    //test
    FESurface *SurfaceIn = mesh.FindSurface("Soluteflux3");
    FESurface *SurfaceOut = mesh.FindSurface("SoluteNaturalFlux1");
    double A_in = 0.0;
    double A_out = 0.0;
    /*std::cout << "test5 "<< dataDim << std::endl;
        //! This function calculates the area of a surface element
        //double FESurface::CurrentFaceArea(FESurfaceElement& el)
    for (int i = 0; i < SurfaceIn->Elements(); i++) {

        A_in += SurfaceIn->CurrentFaceArea(SurfaceIn->Element(i));
        
        }

    for (int i = 0; i < SurfaceOut->Elements(); i++) {

        A_out += SurfaceOut->CurrentFaceArea(SurfaceOut->Element(i));

        }

       std::cout << "A_in " << A_in << std::endl;
       std::cout << "A_out " << A_out << std::endl;
    */
    //FELoadCurve* plc = dynamic_cast<FELoadCurve*>(fem->GetLoadController(2));
    double time =  fem->GetTime().currentTime; //+ fem->GetCurrentStep()->m_dt;

    if (!data.empty()) {
            test = data[0]; //0.01*time;
            //plc->Add(time, test);
            std::cout << "read_data " << test << std::endl;
         } else {
            test = 0.0;
            std::cerr << "Vector is empty." << std::endl;
         }

    std::cout << "update_load_with " << test << std::endl;
    return test; //(currentTime > m_duration ? m_val1 : m_val0);

    //const int dim = precice->getMeshDimensions(otherMesh);
    if(dt!=0.0){
    fem->participant->advance(dt);
    }
    }
    else {
    fem->participant->finalize();
    //delete precice;
    test = 0;
    }

    }

    else {
    std::cout << "return 0" << std::endl;
    return 0; //(currentTime > m_duration ? m_val1 : m_val0);

    }
    //return 0;
  }

BEGIN_FECORE_CLASS(MyLoadController, FELoadController)
    ADD_PARAMETER(index, "id");
    ADD_PARAMETER(r_data, "read_data");
    ADD_PARAMETER(w_data, "write_data");
    ADD_PARAMETER(m_duration, "duration");
END_FECORE_CLASS();

REGISTER_FECORE_CLASS(MyLoadController, "toggle");



MyReflowController::MyReflowController(FEModel *fem) : FELoadController(fem)
{   index = 0.0;
    m_duration = 0.0;
}

bool MyReflowController::Init()
{
     return FELoadController::Init();
}


double MyReflowController::GetValue(double currentTime)
{
  if (currentTime > 3) {
  double reflow = 0;
  double reflow_count = 0;
  std::cout << "femtest" << std::endl;
  FEModel *fem = GetFEModel();
  std::cout << "femtest" << std::endl;
  FEMesh &mesh = fem->GetMesh();
       FEElementSet* elementSetOutflow = mesh.FindElementSet("outflow");
       if (!elementSetOutflow) {
        feLogError((std::string("ElementSet not found")).c_str());
        throw FEException("ElementSet not found");
       }
  std::cout << "femtest 3" << std::endl;
  for (int i = 0; i < elementSetOutflow->Elements(); i++) {
       FEElement &element = elementSetOutflow->Element(i);
       for (int j = 0; j < element.GaussPoints(); j++) {
            std::cout << "femtest 4" << std::endl;
            FEMaterialPoint &materialPoint = *element.GetMaterialPoint(j);
            FESolutesMaterialPoint &ps = *materialPoint.ExtractData<FESolutesMaterialPoint>();
            std::cout << "femtest 5 " << currentTime << std::endl;
            reflow += ps.m_ca[0];
            //reflow += sqrt(ps.m_j[index].x*ps.m_j[index].x + ps.m_j[index].y*ps.m_j[index].y + ps.m_j[index].z*ps.m_j[index].z);
            
            std::cout << "reflow" << reflow << std::endl;
            std::cout << "femtest 6" << std::endl;
            reflow_count += 1;
            std::cout << "reflow_count" << reflow_count << std::endl;
        }
       } 

 
 }

  //return reflow/reflow_count;
  std::cout << "test" << index << std::endl;
  return 1*currentTime;
}

BEGIN_FECORE_CLASS(MyReflowController, FELoadController)
    ADD_PARAMETER(index, "id");
    ADD_PARAMETER(m_duration, "duration");
END_FECORE_CLASS();

REGISTER_FECORE_CLASS(MyReflowController, "reflow");


double getOutflow(FEElementSet* elementSetOutflow, int index){
  double outflow;
  double conc = 0;
  double f;
  double surface_jn;
  //FEModel& fem = *GetFEModel();
  //FEModel *fem = GetFEModel();
  /*FEMesh &mesh = fem->GetMesh();
  FEElementSet* elementSetOutflow = mesh.FindElementSet(outflow_set);
  if (!elementSetOutflow) {
    feLogError((std::string("ElementSet not found")).c_str());
    throw FEException("ElementSet not found");
  }*/
  for (int i = 0; i < elementSetOutflow->Elements(); i++) {
     FEElement &element = elementSetOutflow->Element(i);
     for (int j = 0; j < element.GaussPoints(); j++) {
     FEMaterialPoint &materialPoint = *element.GetMaterialPoint(j);
     FESolutesMaterialPoint &ps = *materialPoint.ExtractData<FESolutesMaterialPoint>();
     outflow += sqrt(ps.m_j[index].x*ps.m_j[index].x + ps.m_j[index].y*ps.m_j[index].y + ps.m_j[index].z*ps.m_j[index].z);
     conc += ps.m_ca[index];
  }
 }
 
 //
 //FESurface *SurfaceOut = mesh.FindSurface(dataSet);
        //FEElementSet* elementSetOutflow = mesh.FindElementSet(dataSet);
        /*for (int i = 0; i < elementSetOutflow->Elements(); i++) {
            FEElement &element = elementSetOutflow->Element(i);
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
        std::cout << "c und w " << w.x << w.y << w.z << c << nint << std::endl;
        // evaluate desired natural solute flux
        //FESurfaceMaterialPoint *mp = dynamic_cast<FESurfaceMaterialPoint*>(element.GetMaterialPoint(1));
        //vec3d dxt = mp->dxr ^ mp->dxs;
        vec3d dxt(1.0, 1.0, 1.0);
        double jn = c*(w*dxt);
        //if (flux->m_bshellb) jn = -jn;
        //FEModel *fem = GetFEModel();
        double dt = 0.1; //fem->GetCurrentStep()->m_dt;
        // molar flow rate
        f = jn* dt;
        surface_jn += f;
        }
        //double H_i = dof_a.shape;
        //fa[0] = H_i * f;
        std::cout << "surface_jn" << surface_jn << std::endl;
        outflow = surface_jn;*/


  //std::cout << "c und w " << w << c << std::endl;
  double outflow_return = outflow/(1);
  //std::cout << "outflow " << outflow_return << "conc " << conc << std::endl;
  
  return conc/8;
}

