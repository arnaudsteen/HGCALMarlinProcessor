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
  vec.push_back(1.0);
  registerProcessorParameter( "ShowerAnalyser::EnergyCalibrationFactors" ,
 			      "Calibration factors to reconstruct the shower energy",
 			      m_ShowerAnalyserSetting.energyCalibrationFactors,
 			      vec );
  /*-------------------------------------------------*/
}


void hgcalShowerProcessor::init()
{ 
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  /*--------------------Algorithms initialisation--------------------*/
  algo_ShowerAnalyser=new algorithm::ShowerAnalyser();
  algo_ShowerAnalyser->SetShowerAnalyserParameterSetting(m_ShowerAnalyserSetting);
  
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

}

/*---------------------------------------------------------------------------*/

void hgcalShowerProcessor::DoShower()
{
  caloobject::Shower* shower=new caloobject::Shower(hitVec);
  algo_ShowerAnalyser->Run(shower);

  energy=shower->getEnergy();
  edep=shower->getEdep();
  nlayer=shower->getNlayer();
  reconstructedCosTheta=shower->getReconstructedCosTheta();
  transverseRatio=shower->getTransverseRatio();
  eta=shower->getEta();
  phi=shower->getPhi();
  f1=shower->getF1();
  showerMax=shower->getShowerMax();
  edepAtMax=shower->getEdepAtMax();
  edepPerCell=shower->getEdepPerCell();

  outTree->Fill();
  delete shower;
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
	hitVec.push_back(aHit);
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
  for(std::vector<caloobject::CaloHit*>::iterator it=hitVec.begin(); it!=hitVec.end(); ++it)
    delete (*it);

  hitVec.clear();
}


void hgcalShowerProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void hgcalShowerProcessor::end(){ 
  outFile->Write();
  outFile->Close();
  delete algo_ShowerAnalyser;
}
