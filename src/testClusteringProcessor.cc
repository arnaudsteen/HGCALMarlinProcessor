#include "testClusteringProcessor.h"
#include <iostream>
#include <time.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>

#include <TStyle.h>
#include <TText.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <string.h>
using namespace lcio ;
using namespace marlin ;

testClusteringProcessor atestClusteringProcessor ;

testClusteringProcessor::testClusteringProcessor() : Processor("testClusteringProcessor") {

  // modify processor description
  _description = "testClusteringProcessor displays events in HGCAL" ;
  
  
  std::vector<std::string> hgcalCollections;
  hgcalCollections.push_back(std::string("HGCALCalorimeterHit"));
  registerInputCollections( LCIO::CALORIMETERHIT,
			    "CollectionName" , 
			    "HGCAL Collection Names"  ,
			    _hgcalCollections  ,
			    hgcalCollections);

  registerProcessorParameter( "PrefixPlotName" ,
			      "Prefix name for generated plots" ,
			      _prefixPlotName ,
			      std::string("muon") );

  registerProcessorParameter( "PauseAfterDraw" ,
			      "Boolean to set if a pause is wanted after diplaying event (using WaitPrimitive method of TCanvas)" ,
			      _pauseAfterDraw ,
			      (bool) false );
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
 			      "z position of the first hgcal layer",
 			      m_CaloGeomSetting.firstLayerZ,
 			      (float) -144.15 ); 
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

  registerProcessorParameter( "Hough::PadSize" ,
    			      "Transversal pad (pixel) size (assuming square geometry)",
    			      m_HoughParameterSetting.padSize,
    			      (float) 10.0 );

  registerProcessorParameter( "Hough::UseAnalogEnergy" ,
    			      "Set true to use 2D cluster deposited energy in active detectors (true for hgcal; false for sdhcal); default value=false",
    			      m_HoughParameterSetting.useAnalogEnergy,
    			      (bool) true );
  
  registerProcessorParameter( "Hough::MaxEnergy" ,
    			      "Maximum deposited energy in GeV in 2D cluster for mip candidate; default value = 1 MeV",
    			      m_HoughParameterSetting.maxEnergy,
    			      (float) 0.00029 );

  registerProcessorParameter( "Hough::PrintDebug" ,
    			      "If true, Hough algorithm will print some debug information",
    			      m_HoughParameterSetting.printDebug,
    			      (bool) false );

  /*------------algorithm::Cluster3D------------*/
  registerProcessorParameter( "Cluster3D::MaxLongitudinal" ,
 			      "Maximum difference between two hits cellID (2) to build a cluster",
 			      m_Cluster3DParameterSetting.maxLongitudinal,
 			      (int) 3 ); 

  registerProcessorParameter( "Cluster3D::MaxTransverseDistance" ,
 			      "Maximum transverse distance (in mm) between two hits to gathered them in one cluster",
 			      m_Cluster3DParameterSetting.maxTransverseDistance,
 			      (float) 30.0 ); 

}


void testClusteringProcessor::init()
{ 
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  /*--------------------Geomttry initialisation--------------------*/
  m_Cluster3DParameterSetting.geometry=m_CaloGeomSetting;
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

  algo_Cluster3D=new algorithm::Cluster3D();
  algo_Cluster3D->SetCluster3DParameterSetting(m_Cluster3DParameterSetting);

  /*-----------------------------------------------------------------*/
  gStyle->SetOptStat(0);
  int argc=0;
  char* argv=(char*)"";
  app = new TApplication("toto",&argc,&argv);  std::vector< std::string > canName;
  canName.push_back( std::string("space") );
  TCanvas *cc=NULL;
  for(unsigned int i=0; i<canName.size(); i++){
    cc=new TCanvas();
    cc->SetWindowSize(800,400);
    cc->Divide(2,1);
    std::pair< std::string,TCanvas* > p(canName.at(i),cc);
    canvasMap.insert(p);  
  }

  std::vector< std::string > histoName;
  histoName.push_back( std::string("spaceXZ_normal") );
  histoName.push_back( std::string("spaceYZ_normal") );
  histoName.push_back( std::string("spaceXZ_track") );
  histoName.push_back( std::string("spaceYZ_track") );
  histoName.push_back( std::string("spaceXZ_cluster") );
  histoName.push_back( std::string("spaceYZ_cluster") );
  std::vector<int> colorVec;
  colorVec.push_back(1);
  colorVec.push_back(1);
  colorVec.push_back(2);
  colorVec.push_back(2);
  colorVec.push_back(4);
  colorVec.push_back(4);

  TH2D *h=NULL;
  for(unsigned int i=0; i<histoName.size(); ++i){
    if( histoName.at(i).find("space")<histoName.at(i).size() ){
      h=new TH2D(histoName.at(i).c_str(),"",4*m_CaloGeomSetting.nLayers,0.0,m_CaloGeomSetting.nLayers,4*m_CaloGeomSetting.nPixelsPerLayer,0.0,m_CaloGeomSetting.nPixelsPerLayer);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(.5);
      h->GetXaxis()->SetTitle("layer");
      if( histoName.at(i).find("X")<histoName.at(i).size() )
	h->GetYaxis()->SetTitle("X [channel number]");
      else if( histoName.at(i).find("Y")<histoName.at(i).size() )
	h->GetYaxis()->SetTitle("Y [channel number]");
      else{
	std::cout << "problem in histo name booking --> std::abort()" << std::endl;
	std::abort();
      }
      h->SetMarkerColor(colorVec.at(i));
    }
    else{
      std::cout << "problem in histo name booking --> std::abort()" << std::endl;
      std::abort();
    }
    std::pair< std::string,TH2D* > p( histoName.at(i),h );
    histo2DMap.insert(p);
  }
  
}

/*---------------------------------------------------------------------------*/

void testClusteringProcessor::resetHisto()
{
  for(std::map< std::string,TH2D* >::iterator it=histo2DMap.begin(); it!=histo2DMap.end(); ++it)
    it->second->Reset();
}

void testClusteringProcessor::fillSpaceHisto( )
{
  for( std::vector<HitAndTag>::iterator it=hitAndTagVec.begin(); it!=hitAndTagVec.end(); ++it ){
    if( (*it).tag==normal ){
      histo2DMap[ std::string("spaceXZ_normal") ]->Fill( (*it).hit->getCellID()[2], (*it).hit->getCellID()[0] );
      histo2DMap[ std::string("spaceYZ_normal") ]->Fill( (*it).hit->getCellID()[2], (*it).hit->getCellID()[1] );
    }
    else if( (*it).tag==cluster ){
      histo2DMap[ std::string("spaceXZ_cluster") ]->Fill( (*it).hit->getCellID()[2], (*it).hit->getCellID()[0] );
      histo2DMap[ std::string("spaceYZ_cluster") ]->Fill( (*it).hit->getCellID()[2], (*it).hit->getCellID()[1] );
    }
    else if( (*it).tag==track ){
      histo2DMap[ std::string("spaceXZ_track") ]->Fill( (*it).hit->getCellID()[2], (*it).hit->getCellID()[0] );
      histo2DMap[ std::string("spaceYZ_track") ]->Fill( (*it).hit->getCellID()[2], (*it).hit->getCellID()[1] );
    }
  }
}

void testClusteringProcessor::drawHistos()
{
  canvasMap[ std::string("space") ]->cd();

  canvasMap[ std::string("space") ]->cd(1);
  histo2DMap[ std::string("spaceXZ_normal") ]->Draw();
  histo2DMap[ std::string("spaceXZ_cluster") ]->Draw("same");
  histo2DMap[ std::string("spaceXZ_track") ]->Draw("same");
  
  canvasMap[ std::string("space") ]->cd(2);
  histo2DMap[ std::string("spaceYZ_normal") ]->Draw();
  histo2DMap[ std::string("spaceYZ_cluster") ]->Draw("same");
  histo2DMap[ std::string("spaceYZ_track") ]->Draw("same");

  canvasMap[ std::string("space") ]->Update();

  if(_pauseAfterDraw)
    canvasMap.begin()->second->WaitPrimitive();

  std::ostringstream name (std::ostringstream::ate);
  name.str("");
  name << std::string("./displays/") << _prefixPlotName << std::string("_event") << _nEvt << std::string(".png");
  canvasMap[ std::string("space") ]->SaveAs(name.str().c_str());

}

/*---------------------------------------------------------------------------*/

void testClusteringProcessor::DoClustering()
{
  std::vector<caloobject::CaloCluster2D*> clusters;
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    algo_Cluster->Run(it->second,clusters);

  std::sort(clusters.begin(), clusters.end(), algorithm::ClusteringHelper::SortClusterByLayer);

  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    if(algo_ClusteringHelper->IsIsolatedCluster(*it,clusters)){
      delete *it; 
      clusters.erase(it); 
      it--;
    }
  
  std::vector< caloobject::CaloTrack* > tracks;
  algo_Hough->runHough(clusters,tracks,algo_Tracking);
  
  for( std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it){
    for( std::vector<caloobject::CaloCluster2D*>::const_iterator jt=(*it)->getClusters().begin(); jt!=(*it)->getClusters().end(); ++jt){
      for(std::vector<caloobject::CaloHit*>::iterator kt=(*jt)->getHits().begin(); kt!=(*jt)->getHits().end(); ++kt){
	for(std::vector<HitAndTag>::iterator lt=hitAndTagVec.begin(); lt!=hitAndTagVec.end(); ++lt){
	  if( (*lt).hit==(*kt) ){
	    (*lt).tag=track;
	    if( std::find(hitVec.begin(), hitVec.end(), (*lt).hit)!=hitVec.end() )
	      hitVec.erase( std::find(hitVec.begin(), hitVec.end(), (*lt).hit) );
	  }
	}
      }
    }
  }
  
  std::vector<caloobject::CaloCluster3D*> clusters3D;
  algo_Cluster3D->Run(hitVec,clusters3D,gunPosition,gunMomentum);

  if(clusters3D.size()!=0){
    std::vector<caloobject::CaloCluster3D*>::iterator clit=clusters3D.begin();
    for(std::vector<HitAndTag>::iterator it=hitAndTagVec.begin(); it!=hitAndTagVec.end(); ++it){
      for(std::vector<caloobject::CaloHit*>::iterator jt=(*clit)->getHits().begin(); jt!=(*clit)->getHits().end(); ++jt)
  	if( (*it).hit==(*jt) )
  	  (*it).tag=cluster;
    }
  }
  
  for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it)
    delete (*it);
  tracks.clear();
  for(std::vector<caloobject::CaloCluster3D*>::iterator it=clusters3D.begin(); it!=clusters3D.end(); ++it)
    delete (*it);
  clusters3D.clear();

  std::cout << "before sorting cluster \t";
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    std::cout << (*it)->getHits().size() << "";
  std::cout << std::endl;

  algorithm::SortClusterBySize<caloobject::CaloCluster2D,caloobject::CaloCluster2D> sorter;
  std::sort( clusters.begin(), clusters.end(), sorter.sort );
  std::cout << "after sorting cluster \t";
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    std::cout << (*it)->getHits().size() << "";
  std::cout << std::endl;
  
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete (*it);
  clusters.clear();
}

/*---------------------------------------------------------------------------*/

void testClusteringProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void testClusteringProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HGCAL Collections of CalorimeterHits* 
  //
  //std::string initString;
  CLHEP::Hep3Vector posShift(0.,0.,0.);
  UTIL::CellIDDecoder<EVENT::CalorimeterHit> IDdecoder("M:3,S-1:3,I:9,J:9,K-1:6");

  std::vector<float> vecP;
  evt->parameters().getFloatVals(std::string("GunPosition"),vecP);
  std::vector<float> vecM;
  evt->parameters().getFloatVals(std::string("particleMomentum_0"),vecM);
  gunPosition = CLHEP::Hep3Vector(vecP.at(0),
				  vecP.at(1),
				  vecP.at(2));
  gunMomentum = CLHEP::Hep3Vector(vecM.at(0),
				  vecM.at(1),
				  vecM.at(2));

  for (unsigned int i(0); i < _hgcalCollections.size(); ++i) {
    std::string colName =  _hgcalCollections[i] ;
    try{
      col = evt->getCollection( _hgcalCollections[i].c_str() ) ;
      numElements = col->getNumberOfElements();
      if(numElements<10)continue;
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
	//	std::cout << "vec.z() = " << vec.z() << " " << hit->getPosition()[2] << std::endl;
	int cellID[]={IDdecoder(hit)["I"],IDdecoder(hit)["J"],IDdecoder(hit)["K-1"]};
	caloobject::CaloHit *aHit=new caloobject::CaloHit(cellID,vec,hit->getEnergy(),hit->getTime(),posShift);
	hitMap[ aHit->getCellID()[2] ].push_back(aHit);
	hitVec.push_back(aHit);
	HitAndTag hitAndTag(aHit);
	hitAndTagVec.push_back(hitAndTag);
      }
      DoClustering();
      fillSpaceHisto();
      drawHistos();
      resetHisto();
      clearVec();
    }
    catch(DataNotAvailableException &e){ 
      std::cout << "Exeption " << std::endl;
    }
  }
  _nEvt ++ ;
  std::cout << "Event processed : " << _nEvt << std::endl;
}

void testClusteringProcessor::clearVec()
{
  for(std::map<int, std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    for( std::vector<caloobject::CaloHit*>::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt)
      delete (*jt);

  hitVec.clear();
  hitMap.clear();
  hitAndTagVec.clear();
}


void testClusteringProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void testClusteringProcessor::end(){ 
  delete algo_Cluster3D;
  delete algo_Cluster;
  delete algo_ClusteringHelper;
  delete algo_Tracking;
  delete algo_InteractionFinder;
  delete algo_Hough;

}
