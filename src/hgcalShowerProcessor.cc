#include "hgcalShowerProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
using namespace lcio ;
using namespace marlin ;

hgcalShowerProcessor ahgcalShowerProcessor ;

hgcalShowerProcessor::hgcalShowerProcessor() : Processor("hgcalShowerProcessor") {

  // modify processor description
  _description = "hgcalShowerProcessor displays events in HGCAL" ;
  
  
  std::vector<std::string> hgcalCollections;
  hgcalCollections.push_back(std::string("HGCALCalorimeterHit"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "CollectionName" , 
			    "HGCAL Collection Names"  ,
			    _hgcalCollections  ,
			    hgcalCollections);

  registerProcessorParameter( "OutputName" ,
			      "Name of the output root file " ,
			      _outName ,
			      std::string("toto.root") );

  /*------------caloobject::CaloGeom------------*/
  registerProcessorParameter( "Geometry::NLayers" ,
 			      "Number of layers",
 			      m_CaloGeomSetting.nLayers,
 			      (int) 28 ); 

  registerProcessorParameter( "Geometry::NPixelsPerLayer" ,
 			      "Number of pixels per layer (assume square geometry)",
 			      m_CaloGeomSetting.nPixelsPerLayer,
 			      (int) 64 ); 

  registerProcessorParameter( "Geometry::PixelSize" ,
 			      "Pixel size (assume square pixels)",
 			      m_CaloGeomSetting.pixelSize,
 			      (float) 10.0 ); 
  /*--------------------------------------------*/
  
  /*------------algorithm::InteractionFinder-----------*/
  registerProcessorParameter( "InteractionFinder::MinSize" ,
    			      "Minimum cluster size for to define an interaction point",
    			      m_InteractionFinderSetting.minSize,
    			      (int) 4 ); 

  registerProcessorParameter( "InteractionFinder::MaxRadius" ,
    			      "Maximum transversal distance to look for clusters",
    			      m_InteractionFinderSetting.maxRadius,
    			      (float) 50.0 ); 

  registerProcessorParameter( "InteractionFinder::MaxDepth" ,
    			      "Maximum depth (number of layers) to look for clusters",
    			      m_InteractionFinderSetting.maxDepth,
    			      (int) 4 ); 

  registerProcessorParameter( "InteractionFinder::MinNumberOfCluster" ,
    			      "Minimum number of found clusters (big enough) after the interaction point",
    			      m_InteractionFinderSetting.minNumberOfCluster,
    			      (int) 3 );
  
  registerProcessorParameter( "InteractionFinder::UseAnalogEnergy" ,
    			      "set true to use cluster energy rather than topology (siwcal)",
    			      m_InteractionFinderSetting.useAnalogEnergy,
    			      (bool) true );

  registerProcessorParameter( "InteractionFinder::MinEnergy" ,
    			      "minimum cluster energy to define interaction, default value corresponds to 1 MeV",
    			      m_InteractionFinderSetting.minEnergy,
    			      (float) 0.001 );  
  /*---------------------------------------------------*/

  /*------------algorithm::ShowerAnalyser------------*/
  registerProcessorParameter( "ShowerAnalyser::EnergyCalibrationOption" ,
 			      "Option to set energy calibration method",
 			      m_ShowerAnalyserSetting.energyCalibrationOption,
 			      std::string("SiWEcal") );
  
  registerProcessorParameter( "ShowerAnalyser::FirstSectionLastLayer" ,
 			      "Number of layer to define the end of first calorimeter section",
 			      m_ShowerAnalyserSetting.firstSectionLastLayer,
 			      (int) 10 );

  std::vector<float> vec;
  vec.push_back(1.10325471172506013e+02);
  registerProcessorParameter( "ShowerAnalyser::EnergyCalibrationFactors" ,
 			      "Calibration factors to reconstruct the shower energy",
 			      m_ShowerAnalyserSetting.energyCalibrationFactors,
 			      vec );
  /*-------------------------------------------------*/

  /*------------algorithm::Cluster------------*/
  registerProcessorParameter( "MaxTransversalCellID" ,
 			      "Maximum difference between two hits cellID (0 and 1) to build a cluster",
 			      m_ClusterParameterSetting.maxTransversal,
 			      (int) 1 ); 

  registerProcessorParameter( "MaxLongitudinalCellID" ,
 			      "Maximum difference between two hits cellID (2) to build a cluster",
 			      m_ClusterParameterSetting.maxLongitudinal,
 			      (int) 0 ); 

  registerProcessorParameter( "UseDistanceInsteadCellID" ,
 			      "Boolean to know if clustering algorithms uses distance instead of cellID to cluster hits together",
 			      m_ClusterParameterSetting.useDistanceInsteadCellID,
 			      (bool) false ); 

  registerProcessorParameter( "MaxTransversalDistance" ,
 			      "Maximum transversal distance (in mm) between two hits to gathered them in one cluster",
 			      m_ClusterParameterSetting.maxTransversalDistance,
 			      (float) 11.0 ); 

  registerProcessorParameter( "MaxLongitudinalDistance" ,
 			      "Maximum longitudinal distance (in mm) between two hits to gathered them in one cluster",
 			      m_ClusterParameterSetting.maxLongitudinalDistance,
 			      (float) 27.0 ); 
  /*------------------------------------------*/

  /*------------algorithm::ClusteringHelper------------*/
  registerProcessorParameter( "LongitudinalDistanceForIsolation" ,
    			      "Minimum longitudinal distance (in mm) between one hits and its neighbours to decide if it is isolated",
    			      m_ClusteringHelperParameterSetting.longitudinalDistance,
    			      (float) 100.0 ); 

  registerProcessorParameter( "TranversalDistanceForIsolation" ,
    			      "Minimum transversal distance (in mm) between one hits and its neighbours to decide if it is isolated",
    			      m_ClusteringHelperParameterSetting.transversalDistance,
    			      (float) 200.0 ); 
  /*---------------------------------------------------*/

}


void hgcalShowerProcessor::init()
{ 
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  /*--------------------Algorithms initialisation--------------------*/
  m_ShowerAnalyserSetting.interactionFinderParams = m_InteractionFinderSetting;
  algo_ShowerAnalyser=new algorithm::ShowerAnalyser();
  algo_ShowerAnalyser->SetShowerAnalyserParameterSetting(m_ShowerAnalyserSetting);
  
  algo_InteractionFinder=new algorithm::InteractionFinder();
  algo_InteractionFinder->SetInteractionFinderParameterSetting(m_InteractionFinderSetting);

  algo_Cluster=new algorithm::Cluster();
  algo_Cluster->SetClusterParameterSetting(m_ClusterParameterSetting);

  algo_ClusteringHelper=new algorithm::ClusteringHelper();
  algo_ClusteringHelper->SetClusteringHelperParameterSetting(m_ClusteringHelperParameterSetting);
  
  std::ostringstream os( std::ostringstream::ate );
  os.str(_outName);
  if( os.str().find(std::string(".root"))>os.str().size() )
    os << std::string(".root");
  outFile=new TFile(os.str().c_str(),"RECREATE");
  outTree = new TTree("tree","Hough efficiency tree");

  outTree->Branch("energy",&energy);
  outTree->Branch("edep",&edep);
  outTree->Branch("nlayer",&nlayer);
  outTree->Branch("reconstructedCosTheta",&reconstructedCosTheta);
  outTree->Branch("transverseRatio",&transverseRatio);
  outTree->Branch("eta",&eta);
  outTree->Branch("phi",&phi);
  outTree->Branch("f1",&f1);
  outTree->Branch("showerMax",&showerMax);
  outTree->Branch("edepAtMax",&edepAtMax);
  outTree->Branch("edepPerCell","std::vector<double>",&edepPerCell);
  outTree->Branch("beginX",&beginX);
  outTree->Branch("beginY",&beginY);
  outTree->Branch("beginZ",&beginZ);
  outTree->Branch("findInteraction",&findInteraction);
  outTree->Branch("longitudinal","std::vector<double>",&longitudinal);
  outTree->Branch("transverse","std::vector<double>",&transverse);
  outTree->Branch("distanceToAxis","std::vector<double>",&distanceToAxis);
  outTree->Branch("clustersEnergy","std::vector<double>",&clustersEnergy);

}

/*---------------------------------------------------------------------------*/

void hgcalShowerProcessor::DoShower()
{
  std::vector<caloobject::CaloCluster2D*> clusters;
  
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
    algo_Cluster->Run(it->second,clusters);
  }
  std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);

  caloobject::Shower* shower=new caloobject::Shower(clusters);
  algo_ShowerAnalyser->Run(shower);

  energy=shower->getEnergy();
  edep=shower->getEdep()*1000;
  nlayer=shower->getNlayer();
  reconstructedCosTheta=shower->getReconstructedCosTheta();
  transverseRatio=shower->getTransverseRatio();
  eta=shower->getEta();
  phi=shower->getPhi();
  f1=shower->getF1();
  showerMax=shower->getShowerMax();
  edepAtMax=shower->getEdepAtMax()*1000;
  edepPerCell=shower->getEdepPerCell();
  beginX=shower->getStartingPosition().x();
  beginY=shower->getStartingPosition().y();
  beginZ=shower->getStartingPosition().z();
  findInteraction=shower->findInteraction();
  longitudinal=shower->getLongitudinal();
  transverse=shower->getTransverse();
  distanceToAxis=shower->getDistancesToAxis();
  clustersEnergy=shower->getClustersEnergy();
  //*1000 -> MeV unit
  
  outTree->Fill();
  delete shower;
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete (*it);
  clusters.clear();

}

/*---------------------------------------------------------------------------*/

  void hgcalShowerProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void hgcalShowerProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HGCAL Collections of CalorimeterHits* 
  //
  CLHEP::Hep3Vector posShift(0.,0.,0.);
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  
  for (unsigned int i(0); i < _hgcalCollections.size(); ++i) {
    std::string colName =  _hgcalCollections[i] ;
    try{
      col = evt->getCollection( _hgcalCollections[i].c_str() ) ;
      numElements = col->getNumberOfElements();
      //if(numElements<10)continue;
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
	int cellID[]={IDdecoder(hit)["I"],IDdecoder(hit)["J"],IDdecoder(hit)["K-1"]};
	caloobject::CaloHit *aHit=new caloobject::CaloHit(cellID,vec,hit->getEnergy(),hit->getTime(),posShift);
	hitMap[cellID[2]].push_back(aHit);
      }
      DoShower();
      clearVec();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "DataNotAvailableException" << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

void hgcalShowerProcessor::clearVec()
{
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    for( std::vector<caloobject::CaloHit*>::iterator jt=(it->second).begin(); jt!=(it->second).end(); ++jt)
      delete *(jt);

  hitMap.clear();
}


void hgcalShowerProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void hgcalShowerProcessor::end(){ 
  outFile->Write();
  outFile->Close();
  delete algo_ShowerAnalyser;
  delete algo_InteractionFinder;
  delete algo_Cluster;
  delete algo_ClusteringHelper;
}
