#include "hgcalEfficiencyProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
using namespace lcio ;
using namespace marlin ;

hgcalEfficiencyProcessor ahgcalEfficiencyProcessor ;

hgcalEfficiencyProcessor::hgcalEfficiencyProcessor() : Processor("hgcalEfficiencyProcessor") {

  // modify processor description
  _description = "hgcalEfficiencyProcessor calculates a HGCAL efficiency and pad multiplcity for each SHDCAL layer" ;
  
  
  std::vector<std::string> hgcalCollections;
  hgcalCollections.push_back(std::string("HGCALCalorimeterHit"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "CollectionName" , 
			    "HGCAL Collection Names"  ,
			    _hgcalCollections  ,
			    hgcalCollections);

  registerProcessorParameter( "RootFileName" ,
			      "File name for the root output",
			      outputRootName,
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

  std::vector<float> vec;
  vec.push_back(-320.0);
  vec.push_back(320.0);
  registerProcessorParameter( "Geometry::DetectorTransverseSize" ,
     			      "Define the detector transverse size used by efficiency algorithm (vector size must be 2 or 4; if 2 -> first value is min, second value is max; if 4 -> two first values define x edges , two last values define y edges) ",
     			      edges,
     			      vec ); 
  if( edges.size()==2 ){
     m_CaloGeomSetting.xmin=edges[0];
     m_CaloGeomSetting.ymin=edges[0];
     m_CaloGeomSetting.xmax=edges[1];
     m_CaloGeomSetting.ymax=edges[1];
   } 
   else if( edges.size()==4 ){
     m_CaloGeomSetting.xmin=edges[0];
     m_CaloGeomSetting.xmax=edges[1];
     m_CaloGeomSetting.ymin=edges[2];
     m_CaloGeomSetting.ymax=edges[3];
   }
  else{
    std::cout << "WARING : Wrong number of values in paramater Geometry::DetectorTransverseSize => will use default value -500.0, +500.0" << std::endl;
  }
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

  registerProcessorParameter( "LongitudinalDistanceForIsolation" ,
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
    			      (float) 0.9 ); 

  /*------------algorithm::Efficiency-----------*/
  registerProcessorParameter( "Efficiency::MaxRadius" ,
    			      "Maximum distance parameter to find a hit to consider the layer as efficient",
    			      m_EfficiencyParameterSetting.maxRadius,
    			      (float) 25.0 ); 

  registerProcessorParameter( "Efficiency::HGCALReadout" ,
    			      "Boolean to know if the detector used the semi digital readout",
    			      m_EfficiencyParameterSetting.semiDigitalReadout,
    			      (bool) false ); 

  registerProcessorParameter( "Efficiency::PrintDebug" ,
    			      "If true, Efficiency algorithm will print some debug information",
    			      m_EfficiencyParameterSetting.printDebug,
    			      (bool) false ); 

  m_EfficiencyParameterSetting.trackingParams=m_TrackingParameterSetting;

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

}


void hgcalEfficiencyProcessor::AlgorithmRegistrationParameters()
{
}

void hgcalEfficiencyProcessor::init()
{ 
  printParameters() ;

  file = new TFile(outputRootName.c_str(),"RECREATE");
  
  tree = (TTree*)file->Get("tree");
  if(!tree){
    std::cout << "tree creation" << std::endl; 
    tree = new TTree("tree","Shower variables");
  }
  tree->Branch("eventNum",&_nEvt);
  tree->Branch("eventChi2",&_eventChi2);
  tree->Branch("trackCosTheta",&_trackCosTheta);
  tree->Branch("trackX0",&_trackX0);
  tree->Branch("trackY0",&_trackY0);
  tree->Branch("trackZ0",&_trackZ0);
  tree->Branch("transverseRatio",&_transverseRatio);

  tree->Branch("Efficiency","std::vector<double>",&_efficiency);
  tree->Branch("EffEnergy","std::vector<double>",&_effEnergy);
  tree->Branch("Multiplicity","std::vector<double>",&_multiplicity);
  tree->Branch("Chi2","std::vector<double>",&_chi2);

  _nRun = 0 ;
  _nEvt = 0 ;

  /*--------------------Geometry initialisation--------------------*/
  m_EfficiencyParameterSetting.geometry = m_CaloGeomSetting;
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

  algo_Efficiency=new algorithm::Efficiency();
  algo_Efficiency->SetEfficiencyParameterSetting(m_EfficiencyParameterSetting);

  for(int i=0; i<m_CaloGeomSetting.nLayers; i++){
    caloobject::CaloLayer* aLayer=new caloobject::CaloLayer(i);
    layers.push_back(aLayer);
  }
}

void hgcalEfficiencyProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void hgcalEfficiencyProcessor::DoTracking()
{
  std::vector<caloobject::CaloCluster2D*> clusters;
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
    algo_Cluster->Run(it->second,clusters);
  }
  std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if(algo_ClusteringHelper->IsIsolatedCluster(*it,clusters)){
      delete *it; 
      clusters.erase(it); 
      it--;
    }
  }
  
  caloobject::CaloTrack* track=NULL;
  algo_Tracking->Run(clusters,track);
  if( track!=NULL ){
    _trackX0=track->getTrackStartingCluster()->getPosition().x();
    _trackY0=track->getTrackStartingCluster()->getPosition().y();
    _trackZ0=track->getTrackStartingCluster()->getPosition().z();
    algo_InteractionFinder->Run(clusters,track->getTrackParameters());
    if( algo_InteractionFinder->FindInteraction()==false ){
      _transverseRatio=algo_Tracking->getTransverseRatio();
      _eventChi2=track->getChi2();
      CLHEP::Hep3Vector nx(-1.,0.,track->getTrackParameters()[1]);
      CLHEP::Hep3Vector ny(0.,-1.,track->getTrackParameters()[3]);
      _trackCosTheta=(nx.cross(ny)).cosTheta();
      LayerProperties(clusters);
    }
  }
  file->cd();
  tree->Fill();

  delete track;
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete (*it);
}

void hgcalEfficiencyProcessor::LayerProperties(std::vector<caloobject::CaloCluster2D*> &clusters)
{
  int trackBegin= (*clusters.begin())->getLayerID();
  int trackEnd=(*(clusters.end()-1))->getLayerID();
  if(trackBegin==1) trackBegin=0;
  if(trackEnd==m_CaloGeomSetting.nLayers-2) trackEnd=m_CaloGeomSetting.nLayers-1;
  
  for(int i=0; i<m_CaloGeomSetting.nLayers; i++){
    _efficiency.push_back(-1); 
    _effEnergy.push_back(0); 
    _multiplicity.push_back(0); 
    _chi2.push_back(0); 
  }
  
  for(int K=trackBegin; K<=trackEnd; K++){
    //caloobject::CaloLayer* aLayer=new caloobject::CaloLayer(K);
    algo_Efficiency->Run(layers.at(K),clusters);

    if( layers.at(K)->getNTracks()!=0 ){
      _efficiency.at(K)=layers.at(K)->getEfficiency();
      _multiplicity.at(K)=layers.at(K)->getMultiplicity();
      _effEnergy.at(K)=layers.at(K)->getEfficiencyEnergy();
    }
    layers.at(K)->Reset();
  }
}

void hgcalEfficiencyProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HGCAL Collections of CalorimeterHits* 
  //
  //std::string initString;
  CLHEP::Hep3Vector posShift(0.,0.,-625.213);
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");

  for (unsigned int i(0); i < _hgcalCollections.size(); ++i) {
    std::string colName =  _hgcalCollections[i] ;
    try{
      col = evt->getCollection( _hgcalCollections[i].c_str() ) ;
      //initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      numElements = col->getNumberOfElements();
      //      UTIL::CellIDDecoder<CalorimeterHit*> idDecoder(col);
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
	int cellID[]={IDdecoder(hit)["I"],IDdecoder(hit)["J"],IDdecoder(hit)["K-1"]};
	caloobject::CaloHit *aHit=new caloobject::CaloHit(cellID,vec,hit->getEnergy(),hit->getTime(),posShift);
	hitMap[cellID[2]].push_back(aHit);
      }
      DoTracking();
      clearVec();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

void hgcalEfficiencyProcessor::clearVec()
{
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    for( std::vector<caloobject::CaloHit*>::iterator jt=(it->second).begin(); jt!=(it->second).end(); ++jt)
      delete *(jt);

  hitMap.clear();
  _efficiency.clear();
  _effEnergy.clear();
  _multiplicity.clear();
  _chi2.clear();
}


void hgcalEfficiencyProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void hgcalEfficiencyProcessor::end(){ 

  delete algo_Cluster;
  delete algo_ClusteringHelper;
  delete algo_Tracking;
  delete algo_InteractionFinder;
  delete algo_Efficiency;

  for(std::vector<caloobject::CaloLayer*>::iterator it=layers.begin(); it!=layers.end(); ++it)
    delete (*it);
  layers.clear();
  
  file->Write();
  file->Close();
}
