#ifndef hgcalMuonFinder_h
#define hgcalMuonFinder_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <map>

#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloGeom.h"
#include "Algorithm/Cluster.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Hough.h"
#include "Algorithm/Distance.h"

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>

using namespace lcio ;
using namespace marlin ;


class  myVerbosityLevel{
 public :
  myVerbosityLevel(std::string aString){
    if( aString==std::string("DEBUG") )  theLevel=DEBUG;
    else if( aString==std::string("MESSAGE") )  theLevel=MESSAGE;
    else if( aString==std::string("WARNING") )  theLevel=WARNING;
    else if( aString==std::string("ERROR") )  theLevel=ERROR;
  }
  enum VERBOSITY_LEVEL{ DEBUG, MESSAGE, WARNING, ERROR };
  VERBOSITY_LEVEL theLevel;
};

class hgcalMuonFinder : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new hgcalMuonFinder ; }
  
  
  hgcalMuonFinder() ;
  
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

  void DoHough(std::vector<caloobject::CaloTrack*> &tracks);
  void tryToFindMuon();
  void clearVec();
  void fillCylinder();
  void eventProperties( LCEvent * evt ) ; 
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hgcalCollections;

  myVerbosityLevel::VERBOSITY_LEVEL processorVerbosityLevel;

 private:
  std::map<int,std::vector<caloobject::CaloHit*> > hitMap;
  
  /*--------------------Global parameters--------------------*/
  int numElements;
  LCCollection * col;
  std::string outputName;
  float efficiencyDistance;
  std::vector<float> layerZPosition;
  float cylinderRadius;
  std::string verbLevel;

  /*------------------------------------------------------------------------------*/

  /*--------------------Algorithms setting parameter structure--------------------*/
  caloobject::GeomParameterSetting m_CaloGeomSetting;
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
  TH2D* trackPosition;
  TH2D* particlesPosition;
  TH2D* recoTracksChi2vsCosTheta;

  TH1D* particlesEta;
  TH1D* recoTracksCosTheta;
  TH1D* recoTracksChi2;
  
  float distanceToProjection;
  int ntrack; 
  float eta;
  float phi;
  float simEta;
  float simPhi;
  float angleSimRec;
  
  std::vector<double> muonClusterEnergy;
  std::vector<int> nhitInCylinder;
  std::vector<double> energyInCylinder;

  /*-------------------Event parameters-------------------*/
  CLHEP::Hep3Vector gunPosition;
  CLHEP::Hep3Vector gunProjection; 
  CLHEP::Hep3Vector gunMomentum;
  /*------------------------------------------------------*/

} ;

#endif

