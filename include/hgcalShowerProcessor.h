#ifndef hgcalShowerProcessor_h
#define hgcalShowerProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>

#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloGeom.h"
#include "CaloObject/Shower.h"
#include "Algorithm/ShowerAnalyser.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

using namespace lcio ;
using namespace marlin ;

class hgcalShowerProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new hgcalShowerProcessor ; }
  
  
  hgcalShowerProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  void DoShower(); 
  void clearVec();
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hgcalCollections;

 private:
  std::vector<caloobject::CaloHit*> hitVec;
  
  /*--------------------Global parameters--------------------*/
  int numElements;
  LCCollection * col;
  std::string _outName;
  /*------------------------------------------------------------------------------*/

  /*--------------------CaloObjects setting parameter structure--------------------*/
  caloobject::GeomParameterSetting m_CaloGeomSetting;
  /*------------------------------------------------------------------------------*/

  /*--------------------Algorithms list to initialise--------------------*/
  algorithm::ShowerAnalyser *algo_ShowerAnalyser;
  /*------------------------------------------------------------------------------*/
  
  /*--------------------Algorithms setting parameter structure--------------------*/
  algorithm::ShowerAnalyserParameterSetting m_ShowerAnalyserSetting; 
  /*------------------------------------------------------------------------------*/
    
  /*--------------------Root output object--------------------*/
  TFile *outFile; 
  TTree* outTree;

  float energy;
  float edep;
  int nlayer;
  float reconstructedCosTheta;
  float transverseRatio;
  float eta;
  float phi;
    
  float f1; //edep in 10 first layers/total edep
  float showerMax; //x0 unit
  float edepAtMax;
  std::vector<double> edepPerCell;   
} ;

#endif
