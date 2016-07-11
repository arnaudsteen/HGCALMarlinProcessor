#ifndef fsrStudyProcessor_h
#define fsrStudyProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <map>

#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloGeom.h"
#include "CaloObject/Shower.h"
#include "Algorithm/ShowerAnalyser.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Cluster.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/Hough.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/Distance.h"
#include "Algorithm/Cluster3D.h"
#include "Algorithm/CalohitHelper.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

using namespace lcio ;
using namespace marlin ;

class fsrStudyProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new fsrStudyProcessor ; }
  
  
  fsrStudyProcessor() ;
  
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
  void DoHough(std::vector<caloobject::CaloCluster2D*> &clusters, caloobject::CaloTrack* &track); 
  void clearVec();
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hgcalCollections;

 private:
  std::map<int, std::vector<caloobject::CaloHit*> > hitMap;
  
  /*--------------------Global parameters--------------------*/
  int numElements;
  LCCollection * col;
  std::string _outName;
  float _minDistanceToProjection;
  bool _do3DClustering;
  /*------------------------------------------------------------------------------*/

  /*--------------------CaloObjects setting parameter structure--------------------*/
  caloobject::GeomParameterSetting m_CaloGeomSetting;
  /*------------------------------------------------------------------------------*/

  /*--------------------Algorithms list to initialise--------------------*/
  algorithm::ShowerAnalyser *algo_ShowerAnalyser;
  algorithm::InteractionFinder *algo_InteractionFinder;
  algorithm::Cluster *algo_Cluster;
  algorithm::ClusteringHelper *algo_ClusteringHelper;
  algorithm::Tracking *algo_Tracking;
  algorithm::Hough *algo_Hough;
  algorithm::Cluster3D *algo_Cluster3D;
  /*------------------------------------------------------------------------------*/
  
  /*--------------------Algorithms setting parameter structure--------------------*/
  algorithm::ShowerAnalyserParameterSetting m_ShowerAnalyserSetting;
  algorithm::InteractionFinderParameterSetting m_InteractionFinderSetting;
  algorithm::clusterParameterSetting m_ClusterParameterSetting; 
  algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting;
  algorithm::HoughParameterSetting m_HoughParameterSetting;
  algorithm::TrackingParameterSetting m_TrackingParameterSetting; 
  algorithm::cluster3DParameterSetting m_Cluster3DParameterSetting; 
  /*------------------------------------------------------------------------------*/
    
  /*--------------------Root output object--------------------*/
  TFile *outFile; 
  TTree* outTree;

  float energy;
  float edep;
  float meanEdep;
  float rmsEdep;
  int nlayer;
  float reconstructedCosTheta;
  float transverseRatio;
  float eta;
  float phi;
    
  float f1; //edep in 10 first layers/total edep
  float f2; //edep in 20 first layers/total edep
  float f3; //edep in 30 first layers/total edep
  float f4; //edep in 40 first layers/total edep
  float showerMax; //x0 unit
  float edepAtMax;
  float beginX;
  float beginY;
  float beginZ;
  bool findInteraction;
  std::vector<double> edepPerCell;
  std::vector<double> longitudinal;
  std::vector<double> transverse;
  std::vector<double> distanceToAxis;
  std::vector<double> clustersEnergy;
  std::vector<double> hitTimes;   

  float distanceToProjection; 
  int ntrack; 
  float muonCosTheta; 

  std::vector<double> muonClusterEnergy;
  std::vector<double> muonGunPosition; 
  std::vector<double> muonGunMomentum; 

} ;

#endif
