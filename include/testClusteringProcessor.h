#ifndef testClusteringProcessor_h
#define testClusteringProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <cstring>
#include <EVENT/CalorimeterHit.h>
#include <vector>
#include <map>

#include "CaloObject/CaloHit.h"
#include "CaloObject/CaloGeom.h"
#include "Algorithm/Cluster3D.h"
#include "Algorithm/Cluster.h"
#include "Algorithm/Tracking.h"
#include "Algorithm/ClusteringHelper.h"
#include "Algorithm/InteractionFinder.h"
#include "Algorithm/Hough.h"

#include <TCanvas.h>
#include <TH2D.h>
#include <TApplication.h>

using namespace lcio ;
using namespace marlin ;

enum HitTag{
  normal,
  track,
  cluster
};

struct HitAndTag{
  caloobject::CaloHit* hit;
  HitTag tag;
  HitAndTag(caloobject::CaloHit* aHit) : hit(aHit),
                                         tag(normal){;}
};

class testClusteringProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new testClusteringProcessor ; }
  
  
  testClusteringProcessor() ;
  
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

  void DoClustering();
  void fillSpaceHisto();
  void drawHistos();
  void resetHisto();
  void clearVec();
 protected:

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hgcalCollections;

 private:
  std::vector<caloobject::CaloHit*> hitVec;
  std::map<int,std::vector<caloobject::CaloHit*> > hitMap;
  std::vector<HitAndTag> hitAndTagVec;
  
  /*--------------------Global parameters--------------------*/
  int numElements;
  LCCollection * col;
  std::string _prefixPlotName;
  bool _pauseAfterDraw;
  /*---------------------------------------------------------*/


  /*--------------------Algorithms setting parameter structure--------------------*/
  caloobject::GeomParameterSetting m_CaloGeomSetting;
  /*--------------------------------------------------------------------------------*/

  /*--------------------Algorithms list to initialise--------------------*/
  algorithm::Cluster3D *algo_Cluster3D;
  algorithm::Cluster *algo_Cluster;
  algorithm::ClusteringHelper *algo_ClusteringHelper;
  algorithm::Tracking *algo_Tracking;
  algorithm::InteractionFinder *algo_InteractionFinder;
  algorithm::Hough *algo_Hough;
  /*---------------------------------------------------------------------*/
  
  /*--------------------Algorithms setting parameter structure--------------------*/
   algorithm::cluster3DParameterSetting m_Cluster3DParameterSetting; 
   algorithm::clusterParameterSetting m_ClusterParameterSetting; 
   algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting; 
   algorithm::TrackingParameterSetting m_TrackingParameterSetting; 
   algorithm::InteractionFinderParameterSetting m_InteractionFinderParameterSetting;
   algorithm::HoughParameterSetting m_HoughParameterSetting; 
  /*------------------------------------------------------------------------------*/
  
   std::map< std::string,TCanvas* > canvasMap;
   std::map< std::string,TH2D* > histo2DMap;
   TApplication* app;
   CLHEP::Hep3Vector gunPosition;
   CLHEP::Hep3Vector gunMomentum;
} ;

#endif
