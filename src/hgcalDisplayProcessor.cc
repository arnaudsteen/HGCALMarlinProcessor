#include "hgcalDisplayProcessor.h"
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

hgcalDisplayProcessor ahgcalDisplayProcessor ;

hgcalDisplayProcessor::hgcalDisplayProcessor() : Processor("hgcalDisplayProcessor") {

  // modify processor description
  _description = "hgcalDisplayProcessor displays events in HGCAL" ;
  
  
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

  registerProcessorParameter( "EfficiencyDistance" ,
			      "Maximum distance between track expected projection and gun position" ,
			      efficiencyDistance,
			      (float)15.0 );
  std::vector<float> vec;
  registerProcessorParameter( "LayerZPosition" ,
			      "z position (in mm ) of each hgcal layers" ,
			      layerZPosition,
			      vec );
  
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
 			      (float) -209.55);  
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


void hgcalDisplayProcessor::init()
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
  /*-----------------------------------------------------------------*/
  gStyle->SetOptStat(0);
  int argc=0;
  char* argv=(char*)"";
  app = new TApplication("toto",&argc,&argv);  std::vector< std::string > canName;
  canName.push_back( std::string("spaceXZ") );
  canName.push_back( std::string("houghXZ") );
  canName.push_back( std::string("spaceYZ") );
  canName.push_back( std::string("houghYZ") );
  TCanvas *cc=NULL;
  for(unsigned int i=0; i<canName.size(); i++){
    cc=new TCanvas();
    cc->SetWindowSize(400,400);
    std::pair< std::string,TCanvas* > p(canName.at(i),cc);
    canvasMap.insert(p);  
  }
  cc=new TCanvas();
  cc->SetWindowSize(800,800);
  std::pair< std::string,TCanvas* > p( "space3D" ,cc);
  canvasMap.insert(p);  

  fZX = new TF1("fZX","pol1",layerZPosition.at(0),layerZPosition.at(m_CaloGeomSetting.nLayers-1));
  fZX->SetLineWidth(2);
  fZX->SetLineColor(kBlue-6);
  fZX->SetLineStyle(7);
  fZY = new TF1("fZY","pol1",layerZPosition.at(0),layerZPosition.at(m_CaloGeomSetting.nLayers-1));
  fZY->SetLineWidth(2);
  fZY->SetLineColor(kBlue-6);
  fZY->SetLineStyle(7);
  
  std::vector< std::string > histoName;
  histoName.push_back( std::string("houghXZ") );
  histoName.push_back( std::string("houghYZ") );
  std::vector<int> colorVec;
  colorVec.push_back(1);
  colorVec.push_back(2);
  colorVec.push_back(4);
  TH2D *h=NULL;
  for(unsigned int i=0; i<histoName.size(); ++i){
    if( histoName.at(i).find("hough")<histoName.at(i).size() ){
      h=new TH2D(histoName.at(i).c_str(),"",
		 m_HoughParameterSetting.thetaSteps,0,M_PI,
		 m_CaloGeomSetting.nPixelsPerLayer/2,-m_CaloGeomSetting.nPixelsPerLayer/4*m_CaloGeomSetting.pixelSize,m_CaloGeomSetting.nPixelsPerLayer/4*m_CaloGeomSetting.pixelSize);
      if( histoName.at(i).find("X")<histoName.at(i).size() )
	h->SetTitle("z-x hough space");
      else
	h->SetTitle("z-y hough space");
      h->GetXaxis()->SetTitle("theta");
      h->GetYaxis()->SetTitle("rho");
    }
    else{
      std::cout << "problem in histo name booking --> std::abort()" << std::endl;
      std::abort();
    }
    std::pair< std::string,TH2D* > p( histoName.at(i),h );
    histo2DMap.insert(p);
  }

  hZX = new TH1D( "hZX","",1000,layerZPosition.at(0),layerZPosition.at(m_CaloGeomSetting.nLayers-1) );
  hZX->GetXaxis()->SetTitle("z [mm]");
  hZX->GetYaxis()->SetTitle("x [mm]");
  hZX->GetXaxis()->SetTitleOffset(1.2);
  hZX->GetYaxis()->SetTitleOffset(1.2);
  hZX->SetMinimum( -m_CaloGeomSetting.nPixelsPerLayer/4*m_CaloGeomSetting.pixelSize );
  hZX->SetMaximum( +m_CaloGeomSetting.nPixelsPerLayer/4*m_CaloGeomSetting.pixelSize );


  hZY = new TH1D( "hZY","",1000,layerZPosition.at(0),layerZPosition.at(m_CaloGeomSetting.nLayers-1) );
  hZY->GetXaxis()->SetTitle("z [mm]");
  hZY->GetYaxis()->SetTitle("y [mm]");
  hZY->GetXaxis()->SetTitleOffset(1.2);
  hZY->GetYaxis()->SetTitleOffset(1.2);
  hZY->SetMinimum( -m_CaloGeomSetting.nPixelsPerLayer/4*m_CaloGeomSetting.pixelSize );
  hZY->SetMaximum( +m_CaloGeomSetting.nPixelsPerLayer/4*m_CaloGeomSetting.pixelSize );
  
  std::vector< std::string > graphName;
  graphName.push_back( std::string("XZnormal") );
  graphName.push_back( std::string("YZnormal") );
  graphName.push_back( std::string("XZtrack") );
  graphName.push_back( std::string("YZtrack") );
  graphName.push_back( std::string("XZrecomuon") );
  graphName.push_back( std::string("YZrecomuon") );
  TGraph* gr=NULL;
  for(unsigned int i=0; i<graphName.size(); ++i){
    gr = new TGraph();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(.5);
    gr->GetXaxis()->SetTitle("z [mm]");
    if( graphName.at(i).find("X")<graphName.at(i).size() )
      gr->GetYaxis()->SetTitle("x [mm]");
    else if( graphName.at(i).find("Y")<graphName.at(i).size() )
      gr->GetYaxis()->SetTitle("y [mm]");
    else{
      std::cout << "problem in graph name booking --> std::abort()" << std::endl;
      std::abort();
    }
    gr->SetMarkerColor(colorVec.at(i/2));
    gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetYaxis()->SetTitleOffset(1.2);
    std::pair< std::string,TGraph* > p( graphName.at(i),gr );
    graphMap.insert(p);
  }
}

/*---------------------------------------------------------------------------*/

void hgcalDisplayProcessor::bookGraph2D()
{

  std::vector<int> colorVec;
  colorVec.push_back(1);
  colorVec.push_back(2);
  colorVec.push_back(4);
  std::vector< std::string > graph2DName;
  graph2DName.push_back( std::string("axis") );
  graph2DName.push_back( std::string("normal") );
  graph2DName.push_back( std::string("track") );
  graph2DName.push_back( std::string("recomuon") );
  TGraph2D* gr2D=NULL;
  for(unsigned int i=0; i<graph2DName.size(); ++i){
    if( graph2DName.at(i).find("axis")<graph2DName.at(i).size() ){
      gr2D=new TGraph2D();
      gr2D->GetXaxis()->SetTitleOffset(1.5);
      gr2D->GetYaxis()->SetTitleOffset(2.0);
      gr2D->GetYaxis()->SetNdivisions(510);
      gr2D->GetZaxis()->SetTitleOffset(1.5);
      gr2D->GetXaxis()->SetTitle("z [mm]");
      gr2D->GetYaxis()->SetTitle("x [mm]");
      gr2D->GetZaxis()->SetTitle("y [mm]");
      std::pair< std::string,TGraph2D* > p( graph2DName.at(i),gr2D );
      graph2DMap.insert(p);
    }
    else{
      if( graph2DName.at(i).find("normal")>graph2DName.at(i).size() &&
	  graph2DName.at(i).find("track")>graph2DName.at(i).size() &&
	  graph2DName.at(i).find("recomuon")>graph2DName.at(i).size() ){
	std::cout << "problem in graph name booking --> std::abort()" << std::endl;
	std::abort();
      }
      gr2D=new TGraph2D();
      gr2D->SetMarkerStyle(20);
      gr2D->SetMarkerSize(.5);
      gr2D->SetMarkerColor(colorVec.at(i-1));
      std::pair< std::string,TGraph2D* > p( graph2DName.at(i),gr2D );
      graph2DMap.insert(p);
    }
  }
}

void hgcalDisplayProcessor::resetHisto()
{
  for(std::map< std::string,TH2D* >::iterator it=histo2DMap.begin(); it!=histo2DMap.end(); ++it)
    it->second->Reset();
  for(std::map< std::string,TGraph* >::iterator it=graphMap.begin(); it!=graphMap.end(); ++it)
    it->second->Set(0);
  for(std::map< std::string,TGraph2D* >::iterator it=graph2DMap.begin(); it!=graph2DMap.end(); ++it)
    it->second->Delete();
  graph2DMap.clear();
}

void hgcalDisplayProcessor::fillSpaceHisto( )
{
  double xA[2];
  xA[0]=layerZPosition.at(0); xA[1]=layerZPosition.at(m_CaloGeomSetting.nLayers-1);
  double yA[2];
  yA[0]=-m_CaloGeomSetting.nPixelsPerLayer/2*m_CaloGeomSetting.pixelSize; yA[1]=m_CaloGeomSetting.nPixelsPerLayer/2*m_CaloGeomSetting.pixelSize;
  double zA[2];
  zA[0]=-m_CaloGeomSetting.nPixelsPerLayer/2*m_CaloGeomSetting.pixelSize; zA[1]=m_CaloGeomSetting.nPixelsPerLayer/2*m_CaloGeomSetting.pixelSize;
  graph2DMap[ std::string("axis") ]->SetPoint(0,xA[0],yA[0],zA[0]);
  graph2DMap[ std::string("axis") ]->SetPoint(1,xA[1],yA[1],zA[1]);
  
  int nn=0;
  int nt=0;
  int nr=0;
  for( std::vector<HitAndTag>::iterator it=hitAndTagVec.begin(); it!=hitAndTagVec.end(); ++it ){
    if( (*it).tag==normal ){
      graphMap[ std::string("XZnormal") ]->SetPoint( nn, (*it).hit->getPosition().z(), (*it).hit->getPosition().x() );
      graphMap[ std::string("YZnormal") ]->SetPoint( nn, (*it).hit->getPosition().z(), (*it).hit->getPosition().y() );
      graph2DMap[ std::string("normal") ]->SetPoint( nn, (*it).hit->getPosition().z(), (*it).hit->getPosition().x(), (*it).hit->getPosition().y() );
      nn++;
    }
    else if( (*it).tag==track ){
      graphMap[ std::string("XZtrack") ]->SetPoint( nt, (*it).hit->getPosition().z(), (*it).hit->getPosition().x() );
      graphMap[ std::string("YZtrack") ]->SetPoint( nt, (*it).hit->getPosition().z(), (*it).hit->getPosition().y() );
      graph2DMap[ std::string("track") ]->SetPoint( nt, (*it).hit->getPosition().z(), (*it).hit->getPosition().x(), (*it).hit->getPosition().y() );
      nt++;
    }
    else if( (*it).tag==reco_muon ){
      graphMap[ std::string("XZrecomuon") ]->SetPoint( nr, (*it).hit->getPosition().z(), (*it).hit->getPosition().x() );
      graphMap[ std::string("YZrecomuon") ]->SetPoint( nr, (*it).hit->getPosition().z(), (*it).hit->getPosition().y() );
      graph2DMap[ std::string("recomuon") ]->SetPoint( nr, (*it).hit->getPosition().z(), (*it).hit->getPosition().x(), (*it).hit->getPosition().y() );
      nr++;
    }
  }
}

void hgcalDisplayProcessor::fillHoughSpaceHisto( std::vector<caloobject::CaloCluster2D*> &clusters )
{
  algorithm::Distance<caloobject::CaloCluster2D,caloobject::CaloCluster2D> dist;
  for( std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it ){
    // test mip candidate
    if( (*it)->getEnergy() > m_HoughParameterSetting.maxEnergy )
      continue;
    int nNeighbours=0;
    int nCoreNeighbours=0;
    for( std::vector<caloobject::CaloCluster2D*>::iterator jt=clusters.begin(); jt!=clusters.end(); ++jt ){
      if( (*it)==(*jt) || fabs( (*it)->getPosition().z()-(*jt)->getPosition().z() ) > std::numeric_limits<float>::epsilon() ) continue;
      if( dist.getDistance( (*it),(*jt) ) < m_HoughParameterSetting.transversalDistance ){
	nNeighbours++;
	if( (*jt)->getEnergy() > m_HoughParameterSetting.maxEnergy )
	  nCoreNeighbours++;
      }
    }
    if( nNeighbours > m_HoughParameterSetting.maximumNumberOfNeighboursForMip &&
	nCoreNeighbours > m_HoughParameterSetting.maximumNumberOfCoreNeighboursForMip )
      continue;
    // if arrive here, (*it) is a mip candidate
    for( unsigned int itheta=0; itheta<m_HoughParameterSetting.thetaSteps; itheta++ ){
      float theta = itheta*M_PI/m_HoughParameterSetting.thetaSteps ;
      double rhoX = ( (*it)->getPosition().z()*std::cos(theta) + 
	       (*it)->getPosition().x()*std::sin(theta) ) ;
      histo2DMap[ std::string("houghXZ") ]->Fill( theta, rhoX);

      double rhoY = ( (*it)->getPosition().z()*std::cos(theta) + 
	       (*it)->getPosition().y()*std::sin(theta) );
      histo2DMap[ std::string("houghYZ") ]->Fill( theta, rhoY);
    }
  }
}

void hgcalDisplayProcessor::drawHistos()
{
  std::ostringstream name (std::ostringstream::ate);

  name.str("");
  name << std::string("./displays/") << _prefixPlotName << std::string("_event") << _nEvt << std::string("_zx_houghSpace.png");
  canvasMap[ std::string("houghXZ") ]->cd();
  histo2DMap[ std::string("houghXZ") ]->Draw("colz");
  canvasMap[ std::string("houghXZ") ]->Update();
  //canvasMap[ std::string("houghXZ") ]->SaveAs(name.str().c_str());

  name.str("");
  name << std::string("./displays/") << _prefixPlotName << std::string("_event") << _nEvt << std::string("_zy_houghSpace.png");
  canvasMap[ std::string("houghYZ") ]->cd();
  histo2DMap[ std::string("houghYZ") ]->Draw("colz");
  canvasMap[ std::string("houghYZ") ]->Update();
  //canvasMap[ std::string("houghYZ") ]->SaveAs(name.str().c_str());

  name.str("");
  name << std::string("./displays/") << _prefixPlotName << std::string("_event") << _nEvt << std::string("_zx.png");
  canvasMap[ std::string("spaceXZ") ]->cd();
  hZX->Draw("axis");
  fZX->SetParameters( gunProjection.x(), gunMomentum.x()/gunMomentum.z());
  fZX->Draw("same");
  if( graphMap[ std::string("XZnormal") ]->GetN()>0 )graphMap[ std::string("XZnormal") ]->Draw("psame");
  if( graphMap[ std::string("XZtrack") ]->GetN()>0 )graphMap[ std::string("XZtrack") ]->Draw("psame");
  if( graphMap[ std::string("XZrecomuon") ]->GetN()>0 )graphMap[ std::string("XZrecomuon") ]->Draw("psame");
  canvasMap[ std::string("spaceXZ") ]->Update();
  //canvasMap[ std::string("spaceXZ") ]->SaveAs(name.str().c_str());

  name.str("");
  name << std::string("./displays/") << _prefixPlotName << std::string("_event") << _nEvt << std::string("_zy.png");
  canvasMap[ std::string("spaceYZ") ]->cd();
  hZY->Draw("axis");
  fZY->SetParameters( gunProjection.y(), gunMomentum.y()/gunMomentum.z());
  fZY->Draw("same");
  if( graphMap[ std::string("YZnormal") ]->GetN()>0 )graphMap[ std::string("YZnormal") ]->Draw("psame");
  if( graphMap[ std::string("YZtrack") ]->GetN()>0 )graphMap[ std::string("YZtrack") ]->Draw("psame");
  if( graphMap[ std::string("YZrecomuon") ]->GetN()>0 )graphMap[ std::string("YZrecomuon") ]->Draw("psame");
  canvasMap[ std::string("spaceYZ") ]->Update();
  //canvasMap[ std::string("spaceYZ") ]->SaveAs(name.str().c_str());

  canvasMap[ std::string("space3D") ]->cd();
  canvasMap[ std::string("space3D") ]->SetTheta(17);  //<========
  canvasMap[ std::string("space3D") ]->SetPhi(15);     //<========
  TPolyLine3D line3D(2);
  line3D.SetLineWidth(2);
  line3D.SetLineColor(kBlue);
  line3D.SetLineStyle(2);
  line3D.SetPoint(0,
		  layerZPosition.at(0),
		  gunProjection.x(),
		  gunProjection.y());
  line3D.SetPoint(1,
		  layerZPosition.at(m_CaloGeomSetting.nLayers-1),
		  gunProjection.x()+layerZPosition.at(m_CaloGeomSetting.nLayers-1)*gunMomentum.x()/gunMomentum.z(),
		  gunProjection.y()+layerZPosition.at(m_CaloGeomSetting.nLayers-1)*gunMomentum.y()/gunMomentum.z());
  graph2DMap[ std::string("axis") ]->Draw("p");
  if( graph2DMap[ std::string("normal")   ]->GetN()>0 ) graph2DMap[ std::string("normal") ]->Draw("psame");
  if( graph2DMap[ std::string("track")    ]->GetN()>0 ) graph2DMap[ std::string("track") ]->Draw("psame");
  if( graph2DMap[ std::string("recomuon") ]->GetN()>0 ) graph2DMap[ std::string("recomuon") ]->Draw("psame");
  line3D.Draw("same");
  canvasMap[ std::string("space3D") ]->Update();

  if(_pauseAfterDraw)
    canvasMap.begin()->second->WaitPrimitive();
}

/*---------------------------------------------------------------------------*/

void hgcalDisplayProcessor::DoHough()
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
  
  std::vector< caloobject::CaloTrack* > tracks;
  clock_t start,end;
  start=clock();
  algo_Hough->runHough(clusters,tracks,algo_Tracking);
  end=clock();
  std::cout << "time to perform hough transform = " <<  ((float)end-(float)start)/CLOCKS_PER_SEC*1000 << std::endl;
  fillHoughSpaceHisto( clusters );
  
  std::cout << "number of created tracks = " << tracks.size() << std::endl;
  for(std::vector<HitAndTag>::iterator it=hitAndTagVec.begin(); it!=hitAndTagVec.end(); ++it)
    for( std::vector<caloobject::CaloTrack*>::iterator jt=tracks.begin(); jt!=tracks.end(); ++jt)
      for( std::vector<caloobject::CaloCluster2D*>::const_iterator kt=(*jt)->getClusters().begin(); kt!=(*jt)->getClusters().end(); ++kt){
	if( (*kt)->getLayerID()!=(*it).hit->getCellID()[2] ) continue;
	for(std::vector<caloobject::CaloHit*>::iterator lt=(*kt)->getHits().begin(); lt!=(*kt)->getHits().end(); ++lt)
	  if( (*it).hit==(*lt) ){
	    (*it).tag=track;
	  }
      }

  tryToFindMuon(tracks);
  
  for(std::vector<caloobject::CaloCluster2D*>::iterator it=clusters.begin(); it!=clusters.end(); ++it)
    delete (*it);
  clusters.clear();
  for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it)
    delete (*it);
  tracks.clear();
  
}

/*---------------------------------------------------------------------------*/

void hgcalDisplayProcessor::tryToFindMuon( std::vector<caloobject::CaloTrack*> &tracks )
{
  algorithm::Distance<CLHEP::Hep3Vector,float> dist;
  float minDist=std::numeric_limits<float>::max();
  std::vector<caloobject::CaloTrack*>::iterator bestIt;
  for(std::vector<caloobject::CaloTrack*>::iterator it=tracks.begin(); it!=tracks.end(); ++it){
    float dist = (float)(gunProjection-(*it)->expectedTrackProjection(gunProjection.z())).mag() ;
    if( dist<minDist ){
      minDist=dist;
      bestIt=it;
    }
  }

  if( minDist<efficiencyDistance ){
    std::cout << "minDist = " << minDist << std::endl;
    for(std::vector<HitAndTag>::iterator it=hitAndTagVec.begin(); it!=hitAndTagVec.end(); ++it)
      for( std::vector<caloobject::CaloCluster2D*>::const_iterator jt=(*bestIt)->getClusters().begin(); jt!=(*bestIt)->getClusters().end(); ++jt){
	if( (*jt)->getLayerID()!=(*it).hit->getCellID()[2] ) continue;
	if( std::find( (*jt)->getHits().begin(),(*jt)->getHits().end(), (*it).hit ) != (*jt)->getHits().end() )
	  (*it).tag=reco_muon;
      }
  }
}


/*---------------------------------------------------------------------------*/

void hgcalDisplayProcessor::processRunHeader( LCRunHeader* run)
{
  _nRun++ ;
  _nEvt = 0;
} 

void hgcalDisplayProcessor::processEvent( LCEvent * evt )
{   
  //
  // * Reading HGCAL Collections of CalorimeterHits* 
  //
  //std::string initString;
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
      if(numElements<10)continue;
      for (int j=0; j < numElements; ++j) {
	CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>( col->getElementAt( j ) ) ;
	CLHEP::Hep3Vector vec(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
	//	std::cout << "vec.z() = " << vec.z() << " " << hit->getPosition()[2] << std::endl;
	int cellID[]={IDdecoder(hit)["I"],IDdecoder(hit)["J"],IDdecoder(hit)["K-1"]};
	caloobject::CaloHit *aHit=new caloobject::CaloHit(cellID,vec,hit->getEnergy(),hit->getTime(),posShift);
	hitMap[cellID[2]].push_back(aHit);
	HitAndTag hitAndTag(aHit);
	hitAndTagVec.push_back(hitAndTag);
      }
      DoHough();
      bookGraph2D();
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

void hgcalDisplayProcessor::clearVec()
{
  for(std::map<int,std::vector<caloobject::CaloHit*> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
    for( std::vector<caloobject::CaloHit*>::iterator jt=(it->second).begin(); jt!=(it->second).end(); ++jt)
      delete *(jt);

  hitMap.clear();
  hitAndTagVec.clear();
}


void hgcalDisplayProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void hgcalDisplayProcessor::end(){ 
  delete algo_Cluster;
  delete algo_ClusteringHelper;
  delete algo_Tracking;
  delete algo_InteractionFinder;
}
