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

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

//#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

// system include files
#include <string>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"

// HGCalMapping
#include "CondFormats/DataRecord/interface/HGCalElectronicsMappingRcd.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingModuleIndexer.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingCellIndexer.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingParameterHostCollection.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "Geometry/HGCalMapping/interface/HGCalMappingTools.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include <nlohmann/json.hpp>

// ECOND Info
// ------------------------------------------------
struct Hits {
  unsigned nhits10 = 0;
  unsigned nhits20 = 0;
  unsigned nhits30 = 0;
  unsigned nhitsBXM1 = 0;

  void Fill(bool timeValid, bool passBXM1) {
    if (timeValid) {
      nhits30 += 1;
    } else {
      if (passBXM1) {
        nhits20 += 1;
      } else {
        nhits10 += 1;
      }
    }
    if (passBXM1) {
      nhitsBXM1 += 1;
    }
  }

  void resetHits(){
    nhits10 = 0;
    nhits20 = 0;
    nhits30 = 0;
    nhitsBXM1 = 0;
  }
};



struct ECONDInfo {
  Hits hits;
  //int u_orig = 0;
  //int v_orig = 0;
  //int type = 0;
  uint32_t nHalfROCS = 0;
  unsigned nDAQWords;
  std::vector<unsigned> vecNDAQWords;

  void FillHits(bool timeValid, bool passBXM1) { 
    hits.Fill(timeValid, passBXM1); 
  }

  void setNumDAQPacketWords() {
    nDAQWords = daqRawDataPacketWordsDec20(nHalfROCS, hits.nhits10, hits.nhits20, hits.nhits30);
    vecNDAQWords.push_back(nDAQWords);
  }

  void printInfo(std::ostream& ostream) {
    ostream << "n10 = " << hits.nhits10 << ", n20 = " << hits.nhits20 << ", n30 = " << hits.nhits30
            << " nDAQWords = " << nDAQWords << ", nHalfROCS = " << (uint32_t)nHalfROCS << std::endl;
  }

  void resetHitsInfo(){
    nDAQWords = 0;
    hits.resetHits();
  }

  // Taken from https://gitlab.cern.ch/agilbert/HGCalBufferModel/-/blob/master/src/EcondAsicHitCounts.cc#L851-909
  unsigned daqRawDataPacketWordsDec20(uint32_t numberOfHalfHgcrocs,  //const EcondAsicDefinition &e,
                                      unsigned n10,
                                      unsigned n20,
                                      unsigned n30) {
    unsigned nWords(0);
    if (numberOfHalfHgcrocs == 0) {
      throw std::runtime_error("numberOfHalfHgcrocs must be non-zero");
    }
    //std::vector<unsigned> subPacketVector(e.numberOfHalfHgcrocs());
    std::vector<unsigned> subPacketVector(numberOfHalfHgcrocs);

    unsigned nEach10(n10 / subPacketVector.size());
    unsigned nHigh10(n10 % subPacketVector.size());
    unsigned nEach20(n20 / subPacketVector.size());
    unsigned nHigh20(n20 % subPacketVector.size());
    unsigned nEach30(n30 / subPacketVector.size());
    unsigned nHigh30(n30 % subPacketVector.size());

    unsigned nCheck10(0);
    unsigned nCheck20(0);
    unsigned nCheck30(0);

    for (unsigned i(0); i < subPacketVector.size(); i++) {
      unsigned nSp10(nEach10 + (i < nHigh10 ? 1 : 0));
      unsigned nSp20(nEach20 + (i < nHigh20 ? 1 : 0));
      unsigned nSp30(nEach30 + (i < nHigh30 ? 1 : 0));

      nCheck10 += nSp10;
      nCheck20 += nSp20;
      nCheck30 += nSp30;

      if ((nSp10 + nSp20 + nSp30) == 0) {  // Empty

        // Header only
        subPacketVector[i] = 1;
        //subPacketVector[i]=0; // No packet at all

        //} else if((nSp10+nSp20+nSp30)==1) { // 1 hit

        //subPacketVector[i]=2;

      } else {
        // Header+channels
        subPacketVector[i] = 2 + (16 * nSp10 + 24 * nSp20 + 32 * nSp30 + 31) / 32;  // Real
        //subPacketVector[i]=2+(12*nSp10+22*nSp20+34*nSp30+31)/32; // Squeezed
      }

      nWords += 2 * subPacketVector[i];
    }

    assert(nCheck10 == n10);
    assert(nCheck20 == n20);
    assert(nCheck30 == n30);

    // Overall header + CRC
    nWords += 2 * 3;

    return nWords;
  }
};

//ZSideType, LocalFEDIdType, CaptureBlockType, EcondIdxType
using EconDType = std::tuple<bool, uint32_t, uint32_t, uint32_t>;
using ValueType = int;
using HitsEleMap = std::map<EconDType, ValueType>;


// HGCScintAnalyzer Class
// ------------------------------------------------
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
  std::map<uint32_t, uint32_t> mapGeoToElectronics(const hgcal::HGCalMappingModuleParamHostCollection &modules,
                                                   const hgcal::HGCalMappingCellParamHostCollection &cells,
                                                   bool geo2ele,
                                                   bool sipm);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void analyzeDigis(edm::Handle<HGCalDigiCollection> &,
                    const HGCalGeometry *);
  uint32_t getRawEleIdFromMap(HGCScintillatorDetId detId);
  void printDetIdBitValues(HGCScintillatorDetId &detId);
  void printEleIdBitValues(HGCalElectronicsId &eleId);
  std::vector<float> minMaxRingPerLayer(const HGCalGeometry *geo, int layer, bool compareRing);
  std::vector<float> minMaxRingPerLayer(std::map<uint32_t, uint32_t> sipm_geo2ele_, int layer, bool compareRing);
  //int daqRawDataPacketWordsDec20(uint32_t numberOfHalfHgcrocs, unsigned n10, unsigned n20, unsigned n30);

  // ----------member data ---------------------------
  int eventsCount_{0};
  int totalADCHits_{0};

  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::EDGetTokenT<HGCalDigiCollection> digisCEE_,digisCEH_;
  std::map<uint32_t, uint32_t> sipm_geo2ele_;
  std::map<uint32_t, uint32_t> sipm_ele2geo_;

  hgcal::RecHitTools rhtools_;

  // HGCal mapping
  edm::ESWatcher<HGCalElectronicsMappingRcd> cfgWatcher_;
  edm::ESGetToken<HGCalMappingCellIndexer, HGCalElectronicsMappingRcd> cellIndexTkn_;
  edm::ESGetToken<hgcal::HGCalMappingCellParamHostCollection, HGCalElectronicsMappingRcd> cellTkn_;
  edm::ESGetToken<HGCalMappingModuleIndexer, HGCalElectronicsMappingRcd> moduleIndexTkn_;
  edm::ESGetToken<hgcal::HGCalMappingModuleParamHostCollection, HGCalElectronicsMappingRcd> moduleTkn_;


  int layerIdxOffset_{26}; // scint layers go from 8-21 instead of 34-47

  std::vector<int> layerTdcHits_, layerToaHits_,layerAdcHits_,layerValidDetIds_;
  std::vector<std::vector<int>> tileTdcHits_, tileToaHits_,tileAdcHits_, tileValidDetIds_;

  TTree *t_events_;
  TTree *t_info_;
  int b_npu_;  
  bool DEBUG=false;
  bool outputECONDFile=true;
  int Nlayers_ = 51;

  // Map to store the maximum rocChannel for each unique set
  std::map<EconDType, ECONDInfo> EcondInfoMap_;
  std::ofstream ofStreamECONDHitsFile_;

  // tile boards info
  std::map<std::string, int> b_tboard_TdcHits_;
  std::map<std::string, int> b_tboard_ToaHits_;
  std::map<std::string, int> b_tboard_AdcHits_;
  std::map<std::string, int> b_tboard_AdcOcc_;
  std::map<std::string, int> b_tboard_ValidDetIds_;
  

  TH1F *h_cellCount_;
  TProfile *h_tdcCountProf_,*h_toaCountProf_,*h_adcCountProf_;
  TProfile2D *h2_tdcCount_, *h2_toaCount_, *h2_adcCount_;
  TH2F *h_adcHitsVsPU_;

  // only histograms within given indices are stored
  TH2D* hlistXYhits_HECback[50];
  TH2D* hlistGlobalXYhits_HECback[50];
  TH2D* hlistXYhits_HECback_Sipm2mm[50];
  TH2D* hlistGlobalXYhits_HECback_Sipm2mm[50];
  TH2D* hlistXYhits_HECback_Sipm4mm[50];
  TH2D* hlistGlobalXYhits_HECback_Sipm4mm[50];
  TH2D* hlistXYhits_HECback_detIdMiss[50];
  TH2D* hlistGlobalXYhits_HECback_detIdMiss[50];
  TH2D* hlistXYhits_HECback_polar[50];
};


#endif