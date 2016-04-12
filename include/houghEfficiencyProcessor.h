#ifndef houghEfficiencyProcessor_h
#define houghEfficiencyProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <map>

#include "CaloObject/CaloHit.h"
#include "Algorithm/Cluster.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Hough.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

using namespace lcio ;
using namespace marlin ;

class houghEfficiencyProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new houghEfficiencyProcessor ; }
  
  
  houghEfficiencyProcessor() ;
  
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

  void DoHough(); 
  void fillHistograms(std::vector<caloobject::CaloTrack*> &tracks);
  void clearVec();
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hgcalCollections;

 private:
  std::map<int,std::vector<caloobject::CaloHit*> > hitMap;
  
  /*--------------------Global parameters--------------------*/
  int _nActiveLayers;
  int numElements;
  LCCollection * col;
  int _nPixelsPerLayer;
  std::string _outName;
  float _efficiencyDistance;
  /*------------------------------------------------------------------------------*/

  /*--------------------Algorithms list to initialise--------------------*/
  algorithm::Cluster *algo_Cluster;
  algorithm::ClusteringHelper *algo_ClusteringHelper;
  algorithm::Tracking *algo_Tracking;
  algorithm::InteractionFinder *algo_InteractionFinder;
  algorithm::Hough *algo_Hough;
  /*------------------------------------------------------------------------------*/
  
  /*--------------------Algorithms setting parameter structure--------------------*/
   algorithm::clusterParameterSetting m_ClusterParameterSetting; 
   algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting; 
   algorithm::TrackingParameterSetting m_TrackingParameterSetting; 
   algorithm::InteractionFinderParameterSetting m_InteractionFinderParameterSetting; 
   algorithm::HoughParameterSetting m_HoughParameterSetting; 
  /*------------------------------------------------------------------------------*/
    
  /*--------------------Root output object--------------------*/
   TFile *outFile; 
   TTree* outTree;

   TH1D* hDistance;
   TH1D* hEfficiency;
   TH1D* hCosThetaSim;
   TH1D* hCosThetaRec;
   TH1D* hNtrack;
   
   float noiseRate; 
   float distanceToVertex; 
   int ntrack; 
   float cosTheta; 
   float eta; 
   float theta; 

   std::vector<double> gunPosition; 
   std::vector<double> gunMomentum; 
} ;

#endif
