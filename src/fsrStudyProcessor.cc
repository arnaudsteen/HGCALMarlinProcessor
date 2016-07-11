#include "fsrStudyProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
using namespace lcio ;
using namespace marlin ;

fsrStudyProcessor afsrStudyProcessor ;

fsrStudyProcessor::fsrStudyProcessor() : Processor("fsrStudyProcessor") {

  // modify processor description
  _description = "fsrStudyProcessor displays events in HGCAL" ;
  
  
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

  registerProcessorParameter( "MinDistanceToProjection" ,
			      "Minimum distance (in mm) to projection to consider reconstructed track as inpiut muon" ,
			      _minDistanceToProjection ,
			      (float) 15.0);

  registerProcessorParameter( "Do3DClustering" ,
			      "boolean to set at true to perform 3D clustering after the hough" ,
			      _do3DClustering ,
			      (bool) false);

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

  registerProcessorParameter( "Geometry::FirstLayerZ" ,
 			      "Z position of the first calorimeter layer",
 			      m_CaloGeomSetting.firstLayerZ,
 			      (float) -144.15); 

  registerProcessorParameter( "Geometry::FirstSectionLastLayer" ,
 			      "Number of layer to define the end of first calorimeter section",
 			      m_CaloGeomSetting.firstSectionLastLayer,
 			      (int) 10 );

  registerProcessorParameter( "Geometry::SecondSectionLastLayer" ,
 			      "Number of layer to define the end of second calorimeter section",
 			      m_CaloGeomSetting.firstSectionLastLayer,
 			      (int) 20 );
  
  registerProcessorParameter( "Geometry::ThirdSectionLastLayer" ,
			      "Number of layer to define the end of third calorimeter section",
 			      m_CaloGeomSetting.firstSectionLastLayer,
 			      (int) 30 );
  
  registerProcessorParameter( "Geometry::FourthSectionLastLayer" ,
 			      "Number of layer to define the end of fourth calorimeter section",
 			      m_CaloGeomSetting.firstSectionLastLayer,
 			      (int) 40 );
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

  /*------------algorithm::Tracking-----------*/
  registerProcessorParameter( "Tracking::ChiSquareLimit" ,
    			      "Maximum value of tracking fit chi2 to construct a track",
    			      m_TrackingParameterSetting.chiSquareLimit,
    			      (float) 100.0 ); 

  registerProcessorParameter( "Tracking::MaxTransverseRatio" ,
    			      "Maximum value of transverse ratio to construct a track",
    			      m_TrackingParameterSetting.maxTransverseRatio,
    			      (float) 0.05 ); 
  
  registerProcessorParameter( "Tracking::PrintDebug" ,
    			      "If true, Tracking algorithm will print some debug information",
    			      m_TrackingParameterSetting.printDebug,
    			      (bool) false ); 

  registerProcessorParameter( "Tracking::CosThetaLimit",
    			      "minimum value for cos theta to keep the track",
    			      m_TrackingParameterSetting.cosThetaLimit,
    			      (float) 0.0 );
  /*------------------------------------------*/

  /*------------algorithm::Hough-----------*/
  registerProcessorParameter( "Hough::NThetas" ,
    			      "Number of steps to loop over theta values",
    			      m_HoughParameterSetting.thetaSteps,
    			      (int) 100 );

  registerProcessorParameter( "Hough::MinimumNBins" ,
    			      "Minimum number of bins to consider HoughBin as a track",
    			      m_HoughParameterSetting.minimumNBins,
    			      (int) 6 );

  registerProcessorParameter( "Hough::MaximumClusterSizeForMip" ,
    			      "Maximum cluster size to be considered as MIP",
    			      m_HoughParameterSetting.maximumClusterSizeForMip,
    			      (int) 2 );

  registerProcessorParameter( "Hough::MaximumNumberOfNeighboursForMip" ,
    			      "Maximum number of neighbours to consider a cluster as a MIP",
    			      m_HoughParameterSetting.maximumNumberOfNeighboursForMip,
    			      (int) 2 );
  
  registerProcessorParameter( "Hough::MaximumNumberOfCoreNeighboursForMip" ,
    			      "Maximum number of core neighbours to consider a cluster as a MIP",
    			      m_HoughParameterSetting.maximumNumberOfCoreNeighboursForMip,
    			      (int) 0 );
  
  registerProcessorParameter( "Hough::TransversalDistance" ,
    			      "Maximum transversal distance between two clusters to consider them as neigbors",
    			      m_HoughParameterSetting.transversalDistance,
    			      (float) 50. );
  
  registerProcessorParameter( "Hough::IsolationDistance" ,
    			      "Maximum distance (in layer unit) for isolation criterion",
    			      m_HoughParameterSetting.isolationDistance,
    			      (int) 2 );

  registerProcessorParameter( "Hough::UseAnalogEnergy" ,
    			      "Set true to use 2D cluster deposited energy in active detectors (true for hgcal; false for sdhcal); default value=false",
    			      m_HoughParameterSetting.useAnalogEnergy,
    			      (bool) true );
  
  registerProcessorParameter( "Hough::MaxEnergy" ,
    			      "Maximum deposited energy in GeV in 2D cluster for mip candidate; default value = 1 MeV",
    			      m_HoughParameterSetting.maxEnergy,
    			      (float) 0.001 );

  registerProcessorParameter( "Hough::PrintDebug" ,
    			      "If true, Hough algorithm will print some debug information",
    			      m_HoughParameterSetting.printDebug,
    			      (bool) false );
  /*---------------------------------------*/
  
  /*------------algorithm::Cluster3D------------*/
  registerProcessorParameter( "Cluster3D::MaxLongitudinal" ,
 			      "Maximum difference between two hits cellID (2) to build a cluster",
 			      m_Cluster3DParameterSetting.maxLongitudinal,
 			      (int) 3 ); 

  registerProcessorParameter( "Cluster3D::MaxTransverseDistance" ,
 			      "Maximum transverse distance (in mm) between two hits to gathered them in one cluster",
 			      m_Cluster3DParameterSetting.maxTransverseDistance,
 			      (float) 30.0 ); 
  /*--------------------------------------------*/

}


void fsrStudyProcessor::init()
{ 
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  /*--------------------Algorithms initialisation--------------------*/
  m_ShowerAnalyserSetting.geometry=m_CaloGeomSetting;
  m_ShowerAnalyserSetting.interactionFinderParams = m_InteractionFinderSetting;
  algo_ShowerAnalyser=new algorithm::ShowerAnalyser();
  algo_ShowerAnalyser->SetShowerAnalyserParameterSetting(m_ShowerAnalyserSetting);
      
  algo_InteractionFinder=new algorithm::InteractionFinder();
  algo_InteractionFinder->SetInteractionFinderParameterSetting(m_InteractionFinderSetting);

  algo_Cluster=new algorithm::Cluster();
  algo_Cluster->SetClusterParameterSetting(m_ClusterParameterSetting);

  algo_ClusteringHelper=new algorithm::ClusteringHelper();
  algo_ClusteringHelper->SetClusteringHelperParameterSetting(m_ClusteringHelperParameterSetting);
  
  algo_Tracking=new algorithm::Tracking();
  algo_Tracking->SetTrackingParameterSetting(m_TrackingParameterSetting);

  m_HoughParameterSetting.geometry=m_CaloGeomSetting;
  algo_Hough=new algorithm::Hough();
  algo_Hough->SetHoughParameterSetting(m_HoughParameterSetting);

  m_Cluster3DParameterSetting.geometry=m_CaloGeomSetting;
  algo_Cluster3D=new algorithm::Cluster3D();
  algo_Cluster3D->SetCluster3DParameterSetting(m_Cluster3DParameterSetting);

  std::ostringstream os( std::ostringstream::ate );
  os.str(_outName);
  if( os.str().find(std::string(".root"))>os.str().size() )
    os << std::string(".root");
  outFile=new TFile(os.str().c_str(),"RECREATE");
  outTree = new TTree("tree","FSR study processor tree");

  outTree->Branch("energy",&energy);
  outTree->Branch("edep",&edep);
  outTree->Branch("meanEdep",&meanEdep);
  outTree->Branch("rmsEdep",&rmsEdep);
  outTree->Branch("nlayer",&nlayer);
  outTree->Branch("reconstructedCosTheta",&reconstructedCosTheta);
  outTree->Branch("transverseRatio",&transverseRatio);
  outTree->Branch("eta",&eta);
  outTree->Branch("phi",&phi);
  outTree->Branch("f1",&f1);
  outTree->Branch("f2",&f2);
  outTree->Branch("f3",&f3);
  outTree->Branch("f4",&f4);
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
  outTree->Branch("hitTimes","std::vector<double>",&hitTimes);

  outTree->Branch("muonGunPosition","std::vector<double>",&muonGunPosition);
  outTree->Branch("muonGunMomentum","std::vector<double>",&muonGunMomentum);
  outTree->Branch("muonClusterEnergy","std::vector<double>",&muonClusterEnergy);
  outTree->Branch("distanceToProjection",&distanceToProjection);
  outTree->Branch("ntrack",&ntrack);
  outTree->Branch("muonCosTheta",&muonCosTheta);

}

/*---------------------------------------------------------------------------*/

void fsrStudyProcessor::DoHough(std::vector<caloobject::CaloCluster2D*> &clusters, caloobject::CaloTrack* &track)
{
  std::vector< caloobject::CaloTrack* > tracks;
  algo_Hough->runHough(clusters,tracks,algo_Tracking);
  muonClusterEnergy.clear();
  CLHEP::Hep3Vector gunPos(muonGunPosition.at(0),
   			   muonGunPosition.at(1),
   			   muonGunPosition.at(2));
  float coeff= ( m_CaloGeomSetting.firstLayerZ - muonGunPosition.at(2) )/muonGunMomentum.at(2);
  
  CLHEP::Hep3Vector gunProj(muonGunPosition.at(0) + coeff*muonGunMomentum.at(0),
			    muonGunPosition.at(1) + coeff*muonGunMomentum.at(1),
			    muonGunPosition.at(2) + coeff*muonGunMomentum.at(2));
  
  algorithm::Distance<CLHEP::Hep3Vector,float> dist;
  distanceToProjection=std::numeric_limits<float>::max();
  muonCosTheta=-1;
  ntrack=tracks.size();

  if( ntrack>0 ){
    std::vector<caloobject::CaloTrack*>::iterator bestTrackIt;
    for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it){
      float dist=(float)(gunProj-(*it)->expectedTrackProjection(gunProj.z())).mag() ;
      if( dist<distanceToProjection ){
	distanceToProjection=dist;
	muonCosTheta=(*it)->getCosTheta();
	bestTrackIt=it;
      }
    }
    if( distanceToProjection < _minDistanceToProjection ){
      track=(*bestTrackIt);
      for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it)
	if( it != bestTrackIt )
	  delete (*it);
    }
    else 
      for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it)
	if( it != bestTrackIt )
	  delete (*it);
    tracks.clear();
  }
}

void fsrStudyProcessor::DoShower()
{
  std::vector<caloobject::CaloCluster2D*> clusters;
  
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
    algo_Cluster->Run(it->second,clusters);
  }
  std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);

  caloobject::CaloTrack* track = NULL;
  DoHough(clusters,track);
   
  if( NULL != track ){
    for(std::vector<caloobject::CaloCluster2D*>::iterator it=track->getClusters().begin(); it!=track->getClusters().end(); ++it){
      muonClusterEnergy.push_back( (*it)->getEnergy() );
      for(std::vector<caloobject::CaloHit*>::iterator jt=(*it)->getHits().begin(); jt!=(*it)->getHits().end(); ++jt)
	if( std::find( hitMap[ (*it)->getLayerID() ].begin(), hitMap[ (*it)->getLayerID() ].end(), (*jt) )!= hitMap[ (*it)->getLayerID() ].end() ){
	  hitMap[ (*it)->getLayerID() ].erase( std::find( hitMap[ (*it)->getLayerID() ].begin(), hitMap[ (*it)->getLayerID() ].end(), (*jt) ) );
	  delete (*jt);
	}
    }
  }
  std::vector<caloobject::CaloHit*> hitVec;
  for(std::map< int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    hitVec.insert(hitVec.begin(), it->second.begin(), it->second.end() );
  
  std::sort(hitVec.begin(), hitVec.end(), CalohitHelper::sortByLayer);

  caloobject::Shower* shower=NULL;
  std::vector<caloobject::CaloCluster3D*> clusters3D;
  if( _do3DClustering ){
    if( NULL != track ){
      algo_Cluster3D->Run( hitVec , clusters3D , track->expectedTrackProjection( m_CaloGeomSetting.firstLayerZ ) , track->orientationVector() );
      if(clusters3D.size()>0){
	caloobject::CaloCluster3D *cluster3D=(*clusters3D.begin());
	shower=new caloobject::Shower(cluster3D);
      }
    }
  }
  else
    shower=new caloobject::Shower(hitVec);

  if( NULL != shower){
    algo_ShowerAnalyser->Run(shower);
    energy=shower->getEnergy();
    edep=shower->getEdep()*1000;
    meanEdep=shower->getMeanEdep()*1000;
    rmsEdep=shower->getRMSEdep()*1000;
    nlayer=shower->getNlayer();
    reconstructedCosTheta=shower->getReconstructedCosTheta();
    transverseRatio=shower->getTransverseRatio();
    eta=shower->getEta();
    phi=shower->getPhi();
    f1=shower->getF1();
    f2=shower->getF2();
    f3=shower->getF3();
    f4=shower->getF4();
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
    hitTimes=shower->getHitTimes();
    //*1000 -> MeV unit
    
    outTree->Fill();
    delete shower;
  }
  for( std::vector<caloobject::CaloCluster3D*>::iterator it=clusters3D.begin(); it!=clusters3D.end(); ++it )
    delete (*it);
  clusters3D.clear();

  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete (*it);
  clusters.clear();
}

/*---------------------------------------------------------------------------*/

void fsrStudyProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void fsrStudyProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HGCAL Collections of CalorimeterHits* 
  //
  CLHEP::Hep3Vector posShift(0.,0.,0.);
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  
  std::vector<float> vecP;evt->parameters().getFloatVals(std::string("GunPosition"),vecP);
  std::vector<float> vecM;evt->parameters().getFloatVals(std::string("particleMomentum_0"),vecM);
  for(int i=0; i<3; i++){
    muonGunPosition.push_back(vecP.at(i));
    muonGunMomentum.push_back(vecM.at(i));
  }

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

void fsrStudyProcessor::clearVec()
{
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    for( std::vector<caloobject::CaloHit*>::iterator jt=(it->second).begin(); jt!=(it->second).end(); ++jt)
      delete *(jt);

  hitMap.clear();
}


void fsrStudyProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void fsrStudyProcessor::end(){ 
  outFile->Write();
  outFile->Close();
  delete algo_ShowerAnalyser;
  delete algo_InteractionFinder;
  delete algo_Cluster;
  delete algo_ClusteringHelper;
  delete algo_Hough;
  delete algo_Tracking;
}
