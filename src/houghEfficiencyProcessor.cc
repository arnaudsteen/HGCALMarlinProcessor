#include "houghEfficiencyProcessor.h"
#include <iostream>
#include <time.h>
#include <algorithm>

#include <EVENT/LCCollection.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
using namespace lcio ;
using namespace marlin ;

houghEfficiencyProcessor ahoughEfficiencyProcessor ;

houghEfficiencyProcessor::houghEfficiencyProcessor() : Processor("houghEfficiencyProcessor") {

  // modify processor description
  _description = "houghEfficiencyProcessor displays events in HGCAL" ;
  
  
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

  registerProcessorParameter( "EfficiencyDistance" ,
			      "Maximum distance between track expected projection and gun position" ,
			      _efficiencyDistance,
			      (float)15.0 );

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

  registerProcessorParameter( "Hough::PrintDebug" ,
    			      "If true, Hough algorithm will print some debug information",
    			      m_HoughParameterSetting.printDebug,
    			      (bool) false );
}


void houghEfficiencyProcessor::init()
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

  algo_Hough=new algorithm::Hough();
  algo_Hough->SetHoughParameterSetting(m_HoughParameterSetting);

  std::ostringstream os( std::ostringstream::ate );
  os.str(_outName);
  if( os.str().find(std::string(".root"))>os.str().size() )
    os << std::string(".root");
  outFile=new TFile(os.str().c_str(),"RECREATE");
  outTree = new TTree("tree","Hough efficiency tree");
  outTree->Branch("noiseRate",&noiseRate);
  outTree->Branch("gunPosition","std::vector<double>",&gunPosition);
  outTree->Branch("gunMomentum","std::vector<double>",&gunMomentum);
  outTree->Branch("distanceToVertex",&distanceToVertex);
  outTree->Branch("ntrack",&ntrack);
  outTree->Branch("cosTheta",&cosTheta);
  outTree->Branch("eta",&eta);
  outTree->Branch("theta",&theta);
  
  hDistance=new TH1D("hDistance","Distance",100,0,50);
  hEfficiency=new TH1D("hEfficiency","Efficiency",100,0,2);
  hCosThetaSim=new TH1D("hCosThetaSim","CosThetaSim",100,0,1.1);
  hCosThetaRec=new TH1D("hCosThetaRec","CosThetaRec",100,0,1.1);
  hNtrack=new TH1D("hNtrack","Ntrack",100,0,10);

}

/*---------------------------------------------------------------------------*/

void houghEfficiencyProcessor::DoHough()
{
  std::vector<caloobject::CaloCluster*> clusters;
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it){
    algo_Cluster->Run(it->second,clusters);
  }
  std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);
  for(std::vector<caloobject::CaloCluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if(algo_ClusteringHelper->IsIsolatedCluster(*it,clusters)){
      delete *it; 
      clusters.erase(it); 
      it--;
    }
  }
  
  std::vector< caloobject::CaloTrack* > tracks;
  algo_Hough->runHough(clusters,tracks,algo_Tracking);

  fillHistograms(tracks);
  
  for(std::vector<caloobject::CaloCluster*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete (*it);
  clusters.clear();
  for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it)
    delete (*it);
  tracks.clear();
  
}

/*---------------------------------------------------------------------------*/

void houghEfficiencyProcessor::fillHistograms(std::vector<caloobject::CaloTrack*> &tracks)
{
  CLHEP::Hep3Vector gunPos(gunPosition.at(0),
			   gunPosition.at(1),
			   gunPosition.at(2));
  algorithm::Distance<CLHEP::Hep3Vector,float> dist;
  float minDist=std::numeric_limits<float>::max();
  cosTheta=-1;
  for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it){
    float dist=(float)(gunPos-(*it)->expectedTrackProjection(gunPos.z())).mag() ;
    if( dist<minDist ){
      minDist=dist;
      cosTheta=(*it)->getCosTheta();
      eta=(*it)->orientationVector().eta();
      theta=(*it)->orientationVector().theta();
    }
  }
  
  hDistance->Fill(minDist);
  distanceToVertex=minDist;
  if(minDist<_efficiencyDistance)
    hEfficiency->Fill(1);
  else
    hEfficiency->Fill(0);
  
  hCosThetaRec->Fill(cosTheta);

  CLHEP::Hep3Vector gunMom(gunMomentum.at(0),
			   gunMomentum.at(1),
			   gunMomentum.at(2));
  hCosThetaSim->Fill(gunMom.cosTheta());
  hNtrack->Fill(tracks.size());
  ntrack=tracks.size();
}

/*---------------------------------------------------------------------------*/

  void houghEfficiencyProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void houghEfficiencyProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HGCAL Collections of CalorimeterHits* 
  //
  CLHEP::Hep3Vector posShift(0.,0.,0.);
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");
  
  std::vector<float> vecP;evt->parameters().getFloatVals(std::string("GunPosition"),vecP);
  std::vector<float> vecM;evt->parameters().getFloatVals(std::string("ParticleMomentum"),vecM);
  noiseRate=evt->getParameters().getFloatVal(std::string("noiseRate"));
  for(int i=0; i<3; i++){
    gunPosition.push_back(vecP.at(i));
    gunMomentum.push_back(vecM.at(i));
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
      DoHough();
      outTree->Fill();
      clearVec();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "DataNotAvailableException" << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

void houghEfficiencyProcessor::clearVec()
{
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    for( std::vector<caloobject::CaloHit*>::iterator jt=(it->second).begin(); jt!=(it->second).end(); ++jt)
      delete *(jt);

  hitMap.clear();
  gunPosition.clear();
  gunMomentum.clear();
}


void houghEfficiencyProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void houghEfficiencyProcessor::end(){ 
  outFile->Write();
  outFile->Close();
  delete algo_Cluster;
  delete algo_ClusteringHelper;
  delete algo_Tracking;
  delete algo_InteractionFinder;
}
