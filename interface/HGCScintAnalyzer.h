#ifndef _HGCScintAnalyzer_h_
#define _HGCScintAnalyzer_h_

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

//#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

// system include files
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"

//
// class declaration
//
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class HGCScintAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
 public:
  explicit HGCScintAnalyzer(const edm::ParameterSet&);
  ~HGCScintAnalyzer() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   
 private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void analyzeDigis(edm::Handle<HGCalDigiCollection> &, const HGCalGeometry *);

  // ----------member data ---------------------------
  int eventsCount_{0};
  int totalADCHits_{0};

  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::EDGetTokenT<HGCalDigiCollection> digisCEE_,digisCEH_;

  int layerIdxOffset_{26}; // scint layers go from 8-21 instead of 34-47


  std::vector<int> layerTdcHits_, layerToaHits_,layerAdcHits_,layerValidDetIds_;
  std::vector<std::vector<int>> tileTdcHits_, tileToaHits_,tileAdcHits_, tileValidDetIds_;

  TTree *t_events_;
  TTree *t_info_;
  int b_npu_;  
  bool DEBUG=false;
  int Nlayers_ = 51;
  
  // tile boards info
  std::map<std::string, int> b_tboard_TdcHits_;
  std::map<std::string, int> b_tboard_ToaHits_;
  std::map<std::string, int> b_tboard_AdcHits_;
  std::map<std::string, int> b_tboard_ValidDetIds_;
  

  TH1F *h_cellCount_;
  TProfile *h_tdcCountProf_,*h_toaCountProf_,*h_adcCountProf_;
  TProfile2D *h2_tdcCount_, *h2_toaCount_, *h2_adcCount_;
  TH2F *h_adcHitsVsPU_;

};


#endif