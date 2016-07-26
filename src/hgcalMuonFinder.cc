#include "hgcalMuonFinder.h"
#include <iostream>
#include <time.h>
#include <algorithm>
#include <math.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
using namespace lcio ;
using namespace marlin ;

hgcalMuonFinder ahgcalMuonFinder ;

hgcalMuonFinder::hgcalMuonFinder() : Processor("hgcalMuonFinder") {

  // modify processor description
  _description = "hgcalMuonFinder displays events in HGCAL" ;
  
  
  std::vector<std::string> hgcalCollections;
  hgcalCollections.push_back(std::string("HGCALCalorimeterHit"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "CollectionName" , 
			    "HGCAL Collection Names"  ,
			    _hgcalCollections  ,
			    hgcalCollections);

  registerProcessorParameter( "OutputName" ,
			      "Name of the output root file " ,
			      outputName ,
			      std::string("toto.root") );

  registerProcessorParameter( "EfficiencyDistance" ,
			      "Maximum distance between track expected projection and gun position" ,
			      efficiencyDistance,
			      (float)15.0 );

  std::vector<float> vec;
  registerProcessorParameter( "LayerZPosition" ,
			      "z position (in mm ) of each hgcal layers" ,
			      layerZPosition,
			      vec );

  registerProcessorParameter( "CylinderRadius" ,
			      "Cylinder radius (in mm ) arround the muon track to count hits and energy" ,
			      cylinderRadius,
			      float(20.0) );

  /*------------caloobject::CaloGeom------------*/
  registerProcessorParameter( "Geometry::NLayers" ,
 			      "Number of layers",
 			      m_CaloGeomSetting.nLayers,
 			      (int) 40 ); 
  registerProcessorParameter( "Geometry::NPixelsPerLayer" ,
 			      "Number of pixels per layer (assume square geometry)",
 			      m_CaloGeomSetting.nPixelsPerLayer,
 			      (int) 192 ); 
  registerProcessorParameter( "Geometry::PixelSize" ,
 			      "Pixel size (assume square pixels)",
 			      m_CaloGeomSetting.pixelSize,
 			      (float) 10.0 ); 
  registerProcessorParameter( "Geometry::FirstLayerZ" ,
 			      "Z position of the first calorimeter layer",
 			      m_CaloGeomSetting.firstLayerZ,
 			      (float) -209.60);  
  /*--------------------------------------------*/

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

  /*------------algorithm::ClusteringHelper------------*/
  registerProcessorParameter( "LongitudinalDistanceForIsolation" ,
    			      "Minimum longitudinal distance (in mm) between one hits and its neighbours to decide if it is isolated",
    			      m_ClusteringHelperParameterSetting.longitudinalDistance,
    			      (float) 100.0 ); 

  registerProcessorParameter( "TranversalDistanceForIsolation" ,
    			      "Minimum transversal distance (in mm) between one hits and its neighbours to decide if it is isolated",
    			      m_ClusteringHelperParameterSetting.transversalDistance,
    			      (float) 200.0 ); 

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

  /*------------algorithm::InteractionFinder-----------*/
  registerProcessorParameter( "InteractionFinder::MinSize" ,
    			      "Minimum cluster size for to define an interaction point",
    			      m_InteractionFinderParameterSetting.minSize,
    			      (int) 4 ); 

  registerProcessorParameter( "InteractionFinder::MaxRadius" ,
    			      "Maximum transversal distance to look for clusters",
    			      m_InteractionFinderParameterSetting.maxRadius,
    			      (float) 50.0 ); 

  registerProcessorParameter( "InteractionFinder::MaxDepth" ,
    			      "Maximum depth (number of layers) to look for clusters",
    			      m_InteractionFinderParameterSetting.maxDepth,
    			      (int) 4 ); 

  registerProcessorParameter( "InteractionFinder::MinNumberOfCluster" ,
    			      "Minimum number of found clusters (big enough) after the interaction point",
    			      m_InteractionFinderParameterSetting.minNumberOfCluster,
    			      (int) 3 ); 

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
}


void hgcalMuonFinder::init()
{ 
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  /*--------------------Geomttry initialisation--------------------*/
  m_HoughParameterSetting.geometry=m_CaloGeomSetting;
  /*---------------------------------------------------------------*/

  /*--------------------Algorithms initialisation--------------------*/
  algo_Cluster=new algorithm::Cluster();
  algo_Cluster->SetClusterParameterSetting(m_ClusterParameterSetting);

  algo_ClusteringHelper=new algorithm::ClusteringHelper();
  algo_ClusteringHelper->SetClusteringHelperParameterSetting(m_ClusteringHelperParameterSetting);
  
  algo_Tracking=new algorithm::Tracking();
  algo_Tracking->SetTrackingParameterSetting(m_TrackingParameterSetting);

  algo_InteractionFinder=new algorithm::InteractionFinder();
  algo_InteractionFinder->SetInteractionFinderParameterSetting(m_InteractionFinderParameterSetting);

  m_HoughParameterSetting.geometry=m_CaloGeomSetting;
  algo_Hough=new algorithm::Hough();
  algo_Hough->SetHoughParameterSetting(m_HoughParameterSetting);

  std::ostringstream os( std::ostringstream::ate );
  os.str(outputName);
  if( os.str().find(std::string(".root"))>os.str().size() )
    os << std::string(".root");
  outFile=new TFile(os.str().c_str(),"RECREATE");
  outTree = new TTree("tree","Hough efficiency tree");
  outTree->Branch("distanceToProjection",&distanceToProjection);
  outTree->Branch("ntrack",&ntrack);
  outTree->Branch("eta",&eta);
  outTree->Branch("phi",&phi);
  outTree->Branch("simEta",&simEta);
  outTree->Branch("simPhi",&simPhi);
  outTree->Branch("angleSimRec",&angleSimRec);
  outTree->Branch("muonClusterEnergy","std::vector<double>",&muonClusterEnergy);
  outTree->Branch("nhitInCylinder","std::vector<int>",&nhitInCylinder);
  outTree->Branch("energyInCylinder","std::vector<double>",&energyInCylinder);
  trackPosition=new TH2D("trackPosition","trackPosition",100,-500,500,100,-500,500);
  particlesPosition=new TH2D("particlesPosition","particlesPosition",100,-500,500,100,-500,500);
  particlesEta=new TH1D("particlesEta","particlesEta",1000,0,15);

  for(int i=0; i<m_CaloGeomSetting.nLayers; i++){
    nhitInCylinder.push_back(0.);
    energyInCylinder.push_back(0.);
  }
}

/*---------------------------------------------------------------------------*/

void hgcalMuonFinder::DoHough(std::vector<caloobject::CaloTrack*> &tracks)
{
  std::vector<caloobject::CaloCluster2D*> clusters;
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
    algo_Cluster->Run(it->second,clusters);
  }
  std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);
  muonClusterEnergy.clear();
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    muonClusterEnergy.push_back( (*it)->getEnergy() );
    if(algo_ClusteringHelper->IsIsolatedCluster(*it,clusters)){
      delete *it; 
      clusters.erase(it); 
      it--;
    }
  }
  algo_Hough->runHough(clusters,tracks,algo_Tracking);
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete (*it);
  clusters.clear();
}

/*---------------------------------------------------------------------------*/

void hgcalMuonFinder::tryToFindMuon()
{
  std::vector<caloobject::CaloTrack*> tracks;
  DoHough(tracks);
  ntrack = tracks.size();
  algorithm::Distance<CLHEP::Hep3Vector,float> dist;
  float minDist=std::numeric_limits<float>::max();
  std::vector<caloobject::CaloTrack*>::iterator bestIt;
  for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it){
    float dist = (float)(gunProjection-(*it)->expectedTrackProjection(gunProjection.z())).mag() ;
    if( dist<minDist ){
      minDist=dist;
      eta=(*it)->orientationVector().eta();
      phi=(*it)->orientationVector().phi();
      distanceToProjection=dist;
      bestIt=it;
    }
  }
  simEta=gunMomentum.eta();
  simPhi=gunMomentum.phi();

  if( ntrack>0 )
    angleSimRec =  (*bestIt)->orientationVector().angle( gunMomentum ) ;
  else{
    angleSimRec=-1.;
    std::cout << "angleSimRec = " << angleSimRec << std::endl;
  }
  
  
  if( ntrack>0 && distanceToProjection<efficiencyDistance )
    trackPosition->Fill( (*bestIt)->expectedTrackProjection( m_CaloGeomSetting.firstLayerZ ).x() ,
			 (*bestIt)->expectedTrackProjection( m_CaloGeomSetting.firstLayerZ ).y() );

  for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it)
    delete (*it);
  tracks.clear(); 
}

/*---------------------------------------------------------------------------*/

void hgcalMuonFinder::fillCylinder()
{
  for(int i=0; i<m_CaloGeomSetting.nLayers; i++){
    nhitInCylinder.at(i)=0.;
    energyInCylinder.at(i)=0.;
  }
  
  algorithm::Distance<caloobject::CaloHit,CLHEP::Hep3Vector> dist;
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
    float coeff=( layerZPosition.at(it->first) - gunPosition.z() )/gunMomentum.z();
    CLHEP::Hep3Vector trackPos( gunPosition.x() + coeff*gunMomentum.x() ,
				gunPosition.y() + coeff*gunMomentum.y() ,
				gunPosition.z() + coeff*gunMomentum.z() );
    for( std::vector<caloobject::CaloHit*>::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt){
      if( dist.getDistance( (*jt),trackPos )<cylinderRadius ){
	nhitInCylinder.at(it->first) += 1;
	energyInCylinder.at(it->first) += (*jt)->getEnergy();
      }
    }
  }
}

/*---------------------------------------------------------------------------*/

  void hgcalMuonFinder::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void hgcalMuonFinder::processEvent( LCEvent * evt )
{   
  //
  // * Reading HGCAL Collections of CalorimeterHits* 
  //
  CLHEP::Hep3Vector posShift(0.,0.,0.);
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  
  std::vector<float> vecP;evt->parameters().getFloatVals(std::string("GunPosition"),vecP);
  std::vector<float> vecM;evt->parameters().getFloatVals(std::string("particleMomentum_0"),vecM);

  float coeff= ( m_CaloGeomSetting.firstLayerZ - vecP.at(2) )/vecM.at(2);

  gunPosition = CLHEP::Hep3Vector( vecP.at(0), vecP.at(1), vecP.at(2) );
  gunMomentum = CLHEP::Hep3Vector( vecM.at(0), vecM.at(1), vecM.at(2) );
  gunProjection = CLHEP::Hep3Vector( vecP.at(0) + coeff*vecM.at(0) ,
				     vecP.at(1) + coeff*vecM.at(1) ,
				     vecP.at(2) + coeff*vecM.at(2) );
    
  for (unsigned int i(0); i < _hgcalCollections.size(); ++i) {
    std::string colName =  _hgcalCollections[i] ;
    try{
      col = evt->getCollection( _hgcalCollections[i].c_str() ) ;
      numElements = col->getNumberOfElements();
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
	int cellID[]={IDdecoder(hit)["I"],IDdecoder(hit)["J"],IDdecoder(hit)["K-1"]};
	caloobject::CaloHit *aHit=new caloobject::CaloHit(cellID,vec,hit->getEnergy(),hit->getTime(),posShift);
	hitMap[cellID[2]].push_back(aHit);
      }
      tryToFindMuon();
      fillCylinder();
      outTree->Fill();
      clearVec();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "DataNotAvailableException" << std::endl;
    }
  }
  _nEvt ++ ;
  eventProperties(evt);
  std::cout << "Event processed : " << _nEvt << std::endl;
}

void hgcalMuonFinder::eventProperties( LCEvent * evt )
{
  int nparticles=evt->parameters().getIntVal(std::string("NumberOfParticles"));
  std::vector<float> vecP;evt->parameters().getFloatVals(std::string("GunPosition"),vecP);
  
  std::ostringstream os( std::ostringstream::ate );
  os.str("");
  for(int i=0; i<nparticles; i++){
    os.str("particleMomentum_");
    os << i;
    std::vector<float> vecM;evt->parameters().getFloatVals(os.str(),vecM);
    float coeff= ( m_CaloGeomSetting.firstLayerZ - vecP.at(2) )/vecM.at(2);
    particlesPosition->Fill( vecP.at(0) + coeff*vecM.at(0), vecP.at(1) + coeff*vecM.at(1) );
    particlesEta->Fill( CLHEP::Hep3Vector(vecM.at(0),
					  vecM.at(1),
					  vecM.at(2)).eta() );
  }
}

void hgcalMuonFinder::clearVec()
{
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    for( std::vector<caloobject::CaloHit*>::iterator jt=(it->second).begin(); jt!=(it->second).end(); ++jt)
      delete *(jt);

  hitMap.clear();
}


void hgcalMuonFinder::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void hgcalMuonFinder::end(){
  outFile->Write();
  outFile->Close();
  delete algo_Cluster;
  delete algo_ClusteringHelper;
  delete algo_Tracking;
  delete algo_InteractionFinder;
}
