#ifndef hgcalDisplayProcessor_h
#define hgcalDisplayProcessor_h 1

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

#include <TCanvas.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TApplication.h>
#include <TPolyLine3D.h>

using namespace lcio ;
using namespace marlin ;

enum HitTag{
  normal,
  track,
  reco_muon
};

struct HitAndTag{
  caloobject::CaloHit* hit;
  HitTag tag;
  HitAndTag(caloobject::CaloHit* aHit) : hit(aHit),
                                         tag(normal){;}
};

class hgcalDisplayProcessor : public Processor {
  
 public:

  virtual Processor*  newProcessor() { return new hgcalDisplayProcessor ; }
  
  
  hgcalDisplayProcessor() ;
  
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
  void bookGraph2D();
  void fillSpaceHisto();
  void fillHoughSpaceHisto( std::vector<caloobject::CaloCluster2D*> &clusters );
  void drawHistos();
  void resetHisto();
  void clearVec();
  void tryToFindMuon( std::vector<caloobject::CaloTrack*> &tracks );

 protected:

  static double fitFunc(double *x, double *par)
  {
    return par[0] + x[0]*par[1];
  }

  int _nRun ;
  int _nEvt ;
  /** Input collection name.
   */
  std::vector<std::string> _hgcalCollections;

 private:
  std::map<int,std::vector<caloobject::CaloHit*> > hitMap;
  std::vector<HitAndTag> hitAndTagVec;
  
  /*--------------------Global parameters--------------------*/
  int numElements;
  LCCollection * col;
  std::string _prefixPlotName;
  bool _pauseAfterDraw;
  float efficiencyDistance;
  std::vector<float> layerZPosition;
  /*---------------------------------------------------------*/


  /*--------------------Algorithms setting parameter structure--------------------*/
  caloobject::GeomParameterSetting m_CaloGeomSetting;
  /*--------------------------------------------------------------------------------*/

  /*--------------------Algorithms list to initialise--------------------*/
  algorithm::Cluster *algo_Cluster;
  algorithm::ClusteringHelper *algo_ClusteringHelper;
  algorithm::Tracking *algo_Tracking;
  algorithm::InteractionFinder *algo_InteractionFinder;
  algorithm::Hough *algo_Hough;
  /*---------------------------------------------------------------------*/
  
  /*--------------------Algorithms setting parameter structure--------------------*/
   algorithm::clusterParameterSetting m_ClusterParameterSetting; 
   algorithm::ClusteringHelperParameterSetting m_ClusteringHelperParameterSetting; 
   algorithm::TrackingParameterSetting m_TrackingParameterSetting; 
   algorithm::InteractionFinderParameterSetting m_InteractionFinderParameterSetting;
   algorithm::HoughParameterSetting m_HoughParameterSetting; 
  /*------------------------------------------------------------------------------*/
  
   std::map< std::string,TCanvas* > canvasMap;
   std::map< std::string,TH2D* > histo2DMap;
   std::map< std::string,TGraph* > graphMap;
   std::map< std::string,TGraph2D* > graph2DMap;
   TH1D *hZX;
   TH1D *hZY;
   TF1 *fZX;
   TF1 *fZY;
   TApplication* app;

   /*-------------------Event parameters-------------------*/
   CLHEP::Hep3Vector gunPosition;
   CLHEP::Hep3Vector gunProjection; 
   CLHEP::Hep3Vector gunMomentum;
   /*------------------------------------------------------*/
} ;

#endif
