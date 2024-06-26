// -*- C++ -*-
//
// Package:    UserCode/HGCScintAnalyzer
// Class:      HGCScintAnalyzer
//
/**\class HGCScintAnalyzer HGCScintAnalyzer.cc UserCode/HGCScintAnalyzer/plugins/HGCScintAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maria Mazza
//         Created:  Wed, 03 Apr 2024 15:11:09 GMT
//
//

// user include files
#include "UserCode/HGCElectronicsValidation/interface/HGCScintAnalyzer.h"
#include "UserCode/HGCElectronicsValidation/interface/HGCTileBoard.h"

// HGCalMapping
#include "CondFormats/DataRecord/interface/HGCalElectronicsMappingRcd.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingModuleIndexer.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingCellIndexer.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingParameterHostCollection.h"
#include "Geometry/HGCalMapping/interface/HGCalMappingTools.h"

// system include files
#include <memory>




HGCScintAnalyzer::HGCScintAnalyzer(const edm::ParameterSet &iConfig)
    : puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      digisCEH_(consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis", "HEback"))),
      cellIndexTkn_(esConsumes()), 
      cellTkn_(esConsumes()),
      moduleIndexTkn_(esConsumes()),
      moduleTkn_(esConsumes())

{
  // do whatever initialization is needed

    edm::Service<TFileService> fs;

    // Open output stream file
    ofStreamECONDHitsFile_.open("outputECONDHits.txt");

    // TTrees
    t_events_ = fs->make<TTree>("Events","Events");
    t_info_ = fs->make<TTree>("Info","Info"); // information that does not change between events

    t_events_->Branch("npu", &b_npu_, "pileup/I");
    for ( const auto &mpair : HGCTileBoards::tb_map ) {
        std::string tb_name = mpair.first;
        TString rs_tb_name = tb_name;
        t_events_->Branch(rs_tb_name+"_tdcHits", &b_tboard_TdcHits_[tb_name], rs_tb_name+"_tdcHits/I");
        t_events_->Branch(rs_tb_name+"_toaHits", &b_tboard_ToaHits_[tb_name], rs_tb_name+"_toaHits/I");
        t_events_->Branch(rs_tb_name+"_adcHits", &b_tboard_AdcHits_[tb_name], rs_tb_name+"_adcHits/I");
        t_info_->Branch(rs_tb_name+"_validDetIds", &b_tboard_ValidDetIds_[tb_name], rs_tb_name+"_validDetIds/I");
        t_events_->Branch(rs_tb_name + "_adcOcc", &b_tboard_AdcOcc_[tb_name], rs_tb_name + "_adcOcc/I");
    }

    
    // Histograms
    //h_validDetIds_=fs->make<TH1F>("validDetIds",";Layer;Number of valid detIDs/layer",50,0.5,50.5);
    h_cellCount_=fs->make<TH1F>("cellcount",";Layer;Number of cells/layer",50,0.5,50.5);
    h_adcHitsVsPU_=fs->make<TH2F>("adchitsvspu",";Number of PU interactions; Number of ADC hits",100,100,300,1000,0,4000); 
    h_tdcCountProf_=fs->make<TProfile>("tdccount",";Layer;Number of hits/layer",50,0.5,50.5);
    h_toaCountProf_=fs->make<TProfile>("toacount",";Layer;Number of hits/layer",50,0.5,50.5);
    h_adcCountProf_=fs->make<TProfile>("adccount",";Layer;Number of hits/layer",50,0.5,50.5);

    h2_tdcCount_=fs->make<TProfile2D>("tdccount2D",";Layer;Ring",50,0.5,50.5,42,0.5,42.5);
    h2_toaCount_=fs->make<TProfile2D>("toacount2D",";Layer;Ring",50,0.5,50.5,42,0.5,42.5);
    h2_adcCount_=fs->make<TProfile2D>("adccount2D",";Layer;Ring",50,0.5,50.5,42,0.5,42.5);

    for (int i = 0; i < 50; ++i) {
      TString layer = std::to_string(i);
      if(i>6 && i <23){
        hlistXYhits_HECback[i] = fs->make<TH2D>("hXYhits_HECback_L" + layer, "Hits in XY", 100, -50., 50., 100, -50., 50.);
        hlistGlobalXYhits_HECback[i] = fs->make<TH2D>("hGlobXYhits_HECback_L" + layer, "Hits in global XY", 600, -300., 300., 600, -300., 300.);
        hlistXYhits_HECback_Sipm2mm[i] = fs->make<TH2D>("hXYhits_HECback_Sipm2mm_L" + layer, "Hits in XY", 100, -50., 50., 100, -50., 50.);
        hlistGlobalXYhits_HECback_Sipm2mm[i] = fs->make<TH2D>("hGlobXYhits_HECback_Sipm2mm_L" + layer, "Hits in XY",600, -300., 300., 600, -300., 300.);
        hlistXYhits_HECback_Sipm4mm[i] = fs->make<TH2D>("hXYhits_HECback_Sipm4mm_L" + layer, "Hits in XY", 100, -50., 50., 100, -50., 50.);
        hlistGlobalXYhits_HECback_Sipm4mm[i] = fs->make<TH2D>("hGlobXYhits_HECback_Sipm4mm_L" + layer, "Hits in XY",600, -300., 300., 600, -300., 300.);
        hlistXYhits_HECback_detIdMiss[i] = fs->make<TH2D>("hXYhits_HECback_detIdMissingFromMap_L" + layer, "Hits in XY", 100, -50., 50., 100, -50., 50.);
        hlistGlobalXYhits_HECback_detIdMiss[i] = fs->make<TH2D>("hGlobXYhits_HECback_detIdMissingFromMap_L" + layer, "Hits in XY", 600, -300., 300., 600, -300., 300.);
        hlistXYhits_HECback_polar[i] = fs->make<TH2D>("hRPhiHits_HECback_polar_L" + layer, "Hits in iring phi", 288, 0, 360, 50, 0., 50.);
      }
    }

    // Vectors
    layerTdcHits_.resize(Nlayers_,0);
    layerToaHits_.resize(Nlayers_,0);
    layerAdcHits_.resize(Nlayers_,0);
    layerValidDetIds_.resize(Nlayers_,0);

    tileTdcHits_.resize(Nlayers_,std::vector<int>(42));
    tileToaHits_.resize(Nlayers_,std::vector<int>(42));
    tileAdcHits_.resize(Nlayers_,std::vector<int>(42));
    tileValidDetIds_.resize(Nlayers_,std::vector<int>(42));
}

// ------------ destructor  ------------
HGCScintAnalyzer::~HGCScintAnalyzer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    //
    // please remove this method altogether if it would be left empty
}

// ------------ method called once each job just before starting event loop  ------------
void HGCScintAnalyzer::beginJob() {}

// ------------ method called for each event  ------------
void HGCScintAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  eventsCount_++;
  using namespace edm;
  
  // reset global counters
  std::fill(layerTdcHits_.begin(), layerTdcHits_.end(), 0);
  std::fill(layerToaHits_.begin(), layerToaHits_.end(), 0);
  std::fill(layerAdcHits_.begin(), layerAdcHits_.end(), 0);

  std::for_each(std::begin(tileTdcHits_), std::end(tileTdcHits_),
   [](auto & v){ std::fill(std::begin(v), std::end(v), 0); });
  std::for_each(std::begin(tileToaHits_), std::end(tileToaHits_),
   [](auto & v){ std::fill(std::begin(v), std::end(v), 0); });
  std::for_each(std::begin(tileAdcHits_), std::end(tileAdcHits_),
   [](auto & v){ std::fill(std::begin(v), std::end(v), 0); });

  for ( const auto &mpair : HGCTileBoards::tb_map ) {
    std::string tb_name = mpair.first;
    b_tboard_TdcHits_.at(tb_name) = 0;
    b_tboard_ToaHits_.at(tb_name) = 0;
    b_tboard_AdcHits_.at(tb_name) = 0;
    b_tboard_ValidDetIds_.at(tb_name) = 0;
  }

  // setup HGCal geom-electronics mapping
  // >>
  //get cell indexers and SoA
  auto cellIdx = iSetup.getData(cellIndexTkn_);
  for (size_t idx = 0; idx < cellIdx.di_.size(); idx++) {
    //check typecode exists
    auto typecode = cellIdx.getTypecodeFromEnum(idx);
    std::cout << "typecode.c_str() = " << typecode.c_str() << std::endl;
  }
  auto const &cells = iSetup.getData(cellTkn_);
  printf("[HGCScintAnalyzer][analyze] Cell dense indexers and associated SoA retrieved for HGCAL\n");
  int nmodtypes = cellIdx.typeCodeIndexer_.size();
  printf("[HGCScintAnalyzer][analyze] module cell indexer has %d module types\n", nmodtypes);
  //get the module mapper SoA
  auto const &modules = iSetup.getData(moduleTkn_);
  sipm_geo2ele_ = this->mapGeoToElectronics(modules, cells, true, true);
  sipm_ele2geo_ = this->mapGeoToElectronics(modules, cells, false, true);
  // << HGCal mapping

  // get detector information 
  edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(caloGeomToken_);
  const HGCalGeometry *geo = static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(DetId::HGCalHSc, ForwardSubdetector::ForwardEmpty));

  // for rechittools
  const CaloGeometry &geomCalo = iSetup.getData(caloGeomToken_);
  rhtools_.setGeometry(geomCalo);

  // find min max rings per layer from geometry and from map
  // std::vector<int> layervals {8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  // for(auto lay: layervals){
  //    std::cout << "Layer " << lay << std::endl;
  //   std::vector<float> minmaxvals_geo = minMaxRingPerLayer(geo, lay, true);
  //   std::cout << "Geo : min ring = " << minmaxvals_geo.at(0) << ", max ring = " << minmaxvals_geo.at(1) << std::endl;
  //   std::vector<float> minmaxvals_geo_iphi = minMaxRingPerLayer(geo, lay, false);
  //   std::cout << "      min iphi = " << minmaxvals_geo_iphi.at(0) << ", max iphi = " << minmaxvals_geo_iphi.at(1)
  //             << std::endl;

  //   std::vector<float> minmaxvals_map = minMaxRingPerLayer(sipm_geo2ele_, lay+1, true);
  //   std::cout << "Map (lay+1) : min ring = " << minmaxvals_map.at(0) << ", max ring = " << minmaxvals_map.at(1) << std::endl;
  //   std::vector<float> minmaxvals_map_iphi = minMaxRingPerLayer(sipm_geo2ele_, lay + 1, false);
  //   std::cout << "              min iphi = " << minmaxvals_map_iphi.at(0)
  //             << ", max iphi = " << minmaxvals_map_iphi.at(1) << std::endl;
  // }

  

  // Fill validDetId histogram for first event only
  if (eventsCount_==1){
    const auto &validDetIds = geo->getValidDetIds();
    for(const auto &didIt : validDetIds) {
      HGCScintillatorDetId detId(didIt.rawId());
      // Choose negative because current module locator set only for negative CE
      if (detId.zside() > 0)
        continue;  // do only one side since symmetrical?

      // tile info
      auto rawDetId = detId.rawId();
      int zside{detId.zside()};
      int type{detId.type()};
      int layer = detId.layer() + layerIdxOffset_;
      int layerDetId = detId.layer();
      int ring{detId.ring()};
      int iradius{detId.iradius()};
      int ieta{detId.ieta()};
      int iphi{detId.iphi()};
      int trigger{detId.trigger()};
      int sipm{detId.sipm()};

      h_cellCount_->Fill(layer);
      // find its tileboar
      for (const auto &[tb_name, tb] : HGCTileBoards::tb_map) {
        if (layer != tb.plane) {  // the detIds taken from geo are numbered correctly
          continue;
        }
        if ((iradius > tb.irmin) && (iradius < tb.irmax)) {
          b_tboard_ValidDetIds_.at(tb_name) += 1;
          break;
        }
      }
    

      //
      // Get number of HGCROC channels per ECON-D and also set the keys in the ECOND map

      uint32_t rawEleID = getRawEleIdFromMap(detId);
      // if rawEleID==0, it means the detId was not found in the map. Skip this detId for now. //FIXME
      if (rawEleID == 0) {
        continue;
      }
      HGCalElectronicsId eleID(rawEleID);
      EconDType key = std::make_tuple(eleID.zSide(), eleID.localFEDId(), eleID.captureBlock(), eleID.econdIdx());
      uint32_t econdeRxValue = eleID.econdeRx() + 1;  // eRx start from 0
      if (EcondInfoMap_.find(key) == EcondInfoMap_.end()) {
        EcondInfoMap_[key].nHalfROCS = econdeRxValue;
      } else {
        if (econdeRxValue > EcondInfoMap_[key].nHalfROCS) {
          EcondInfoMap_[key].nHalfROCS = econdeRxValue;
        }
      }
    }

    // Output unique ECOND indices
    if (outputECONDFile) {
      ofStreamECONDHitsFile_ << "zSide localFEDId captureBlock econdIdx nHalfROCS" << std::endl;
      for (const auto &entry : EcondInfoMap_) {
        ofStreamECONDHitsFile_ << std::get<0>(entry.first) << " " << std::get<1>(entry.first) << " "
                            << std::get<2>(entry.first) << " " << std::get<3>(entry.first) << " "
                            << entry.second.nHalfROCS << std::endl;
      }
      ofStreamECONDHitsFile_ << std::endl << std::endl << std::endl;
    }

  } else {
    // if (eventsCount_ % 100 == 0) {
    //   std::cout << "Event " << eventsCount_ << std::endl;
    // }
  }
  std::cout << "Event " << eventsCount_ << std::endl;
  if (outputECONDFile) { 
    ofStreamECONDHitsFile_ << "Event " << eventsCount_ << std::endl;
  }

  //////////////////////////////////////////////////////////
  // Do mapping stuff 
  // if (eventsCount_ == 1) {
  //   const auto &validDetIds = geo->getValidDetIds();

  //   for (const auto &didIt : validDetIds) {
  //     HGCScintillatorDetId detId(didIt.rawId());

  //     // Choose negative because current module locator set only for negative CE
  //     if (detId.zside() > 0) {
  //       continue;  // do only one side since symmetrical
  //     }
  //     uint32_t rawEleID = getRawEleIdFromMap(detId);
  //     if (rawEleID == 0) {
  //       continue;
  //     }
  //     HGCalElectronicsId eleID(rawEleID);
  //     printEleIdBitValues(eleID);

  //   }
  // }

  // HitsEleMap ACDhitsMap;
  // for (const auto &pair : sipm_ele2geo_) {
  //   HGCalElectronicsId eleID(pair.first);
  //   if (eleID.zSide() > 0)
  //     continue;

  //   printEleIdBitValues(eleID);
  //   EconDType key = std::make_tuple(eleID.zSide(), eleID.localFEDId(), eleID.captureBlock(), eleID.econdIdx());
  //   // Update the map: if the key exists, add 1 to the value, otherwise insert it with value 1
  //   ACDhitsMap[key] += 1;
  // }

  // std::cout << "From map:" << std::endl;
  // for (const auto &entry : ACDhitsMap) {
  //   const auto &key = entry.first;
  //   const auto &value = entry.second;
  //   std::cout << "Key: (" << std::get<0>(key) << ", " << std::get<1>(key) << ", " << std::get<2>(key) << ", "
  //             << std::get<3>(key) << "), ";
  //   std::cout << "Value: " << value << std::endl;
  // }

  // Get number of HGCROC channels per ECON-D
  // Map to store the maximum rocChannel for each unique set
  // std::map<std::tuple<int, int, int, int>, uint8_t> nHalfRocsEcondMap;
  // for (const auto &pair : sipm_ele2geo_) {
  //   HGCalElectronicsId eleID(pair.first);
  //   if (eleID.zSide() > 0)
  //     continue;
  //   auto key = std::make_tuple(eleID.zSide(), eleID.localFEDId(), eleID.captureBlock(), eleID.econdIdx());
  //   if (nHalfRocsEcondMap.find(key) == nHalfRocsEcondMap.end()) {
  //     nHalfRocsEcondMap[key] = eleID.econdeRx();
  //   } else {
  //     nHalfRocsEcondMap[key] = std::max(nHalfRocsEcondMap[key], eleID.econdeRx());
  //   }
  // }
  // // Output the results
  // for (const auto &entry : nHalfRocsEcondMap) {
  //   std::cout << "zSide: " << std::get<0>(entry.first) << ", localFEDId: " << std::get<1>(entry.first)
  //             << ", captureBlock: " << std::get<2>(entry.first) << ", econdIdx: " << std::get<3>(entry.first)
  //             << ", max rocChannel: " << static_cast<int> (entry.second) << std::endl;
  // }
  //////////////////////////////////////////////////////////
  // <------- Do mapping stuff


  // analyze digi collections
  edm::Handle<HGCalDigiCollection> digisHandle;
  iEvent.getByToken(digisCEH_, digisHandle);
  analyzeDigis(digisHandle, geo);

  // get pileup information 
  edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
  iEvent.getByToken(puToken_,PupInfo);
  b_npu_ = 0;
  if(PupInfo.isValid()){
    std::vector<PileupSummaryInfo>::const_iterator ipu;
    for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu){
      if (ipu->getBunchCrossing()!=0) continue; // storing PU info only for BX=0
      b_npu_=ipu->getPU_NumInteractions();
      break;
    }
  }

  // Divide hits per tileboard by 360/10=36 to remove repetition in phi
  for (const auto &[tb_name, tb] : HGCTileBoards::tb_map) {
    b_tboard_ToaHits_.at(tb_name) /= 36;
    b_tboard_TdcHits_.at(tb_name) /= 36;
    b_tboard_AdcHits_.at(tb_name) /= 36;
    b_tboard_ValidDetIds_.at(tb_name) /= 36;
  }

  // fill profile histograms per layer (profiles are averages, so no phi redundancy)
  totalADCHits_ = 0;
  for(size_t layIdx=0; layIdx<layerTdcHits_.size(); layIdx++){
    int layer = layIdx; 
    if (DEBUG) {
      std::cout << "Filling 1D profile layer histograms:" 
                << " layIdx = " << layIdx 
                << " layer = " << layer 
                << "layerTdcHits_ = " << layerTdcHits_.at(layIdx)
                << "layerToaHits_ = " << layerToaHits_.at(layIdx)
                << "layerAdcHits_ = " << layerAdcHits_.at(layIdx)
                << std::endl;
    }
    h_tdcCountProf_->Fill(layer,layerTdcHits_.at(layIdx));
    h_toaCountProf_->Fill(layer,layerToaHits_.at(layIdx));
    h_adcCountProf_->Fill(layer,layerAdcHits_.at(layIdx));
    totalADCHits_+=layerAdcHits_.at(layIdx);
  }  
  h_adcHitsVsPU_->Fill(b_npu_,totalADCHits_, 1); 


  // removing phi redundancy from tiles hits and fill 2D tiles histograms
  for(size_t layIdx=0; layIdx<tileTdcHits_.size(); layIdx++){
    for(size_t rinIdx=0; rinIdx<tileTdcHits_[0].size(); rinIdx++){
      int layer = layIdx;
      tileTdcHits_.at(layIdx).at(rinIdx) /= 36;
      tileToaHits_.at(layIdx).at(rinIdx) /= 36;
      tileAdcHits_.at(layIdx).at(rinIdx) /= 36;

      if (DEBUG) {
        std::cout << "Filling 2D tiles histograms:" 
                  << " layIdx = " << layIdx 
                  << ", rinIdx = " << rinIdx 
                  << ", tileTdcHits_" << tileTdcHits_.at(layIdx).at(rinIdx)
                  << ", tileToaHits_" << tileToaHits_.at(layIdx).at(rinIdx)
                  << ", tileAdcHits_" << tileAdcHits_.at(layIdx).at(rinIdx)
                  << std::endl;
      }
      h2_tdcCount_->Fill(layer,rinIdx,tileTdcHits_.at(layer).at(rinIdx));
      h2_toaCount_->Fill(layer,rinIdx,tileToaHits_.at(layIdx).at(rinIdx));
      h2_adcCount_->Fill(layer,rinIdx,tileAdcHits_.at(layIdx).at(rinIdx));
    }  
  }

  if (DEBUG) {
    for ( const auto &mpair : HGCTileBoards::tb_map ) {
      std::cout << mpair.first;
      std::string tb_name = mpair.first;
      std::cout << ", b_tboard_TdcHits_ " << b_tboard_TdcHits_.at(tb_name);
      std::cout << ", b_tboard_ToaHits_ " << b_tboard_ToaHits_.at(tb_name);
      std::cout << ", b_tboard_AdcHits_ " << b_tboard_AdcHits_.at(tb_name) << std::endl;
    }
  }
   
  // fill info tree with constant information only once
  if (eventsCount_==1){
    t_info_->Fill(); 
  }
  // fill events tree for each event
  t_events_->Fill();
}

void HGCScintAnalyzer::analyzeDigis(edm::Handle<HGCalDigiCollection> &digiColl,
                                    const HGCalGeometry *geom ) {
  // check inputs
  if(!digiColl.isValid() || geom==NULL){
 	  std::cout << "HGCScintAnalyzer::analyzeDigis: digicoll not valid or no geom. Returning." << std::endl;
	  return;	
  }

  const int inTimeSample{2}; // in-time BX

  for(auto &hit : *digiColl){
    if(hit.size()==0) continue;

    HGCScintillatorDetId detId{hit.id()};
    auto rawDetId = detId.rawId();

    // take negative side of detector bc ModuleLocator currently setup for z<0 only
    if(detId.zside()>0) continue; // make sure same as inside analysis()

    // detId info
    int layer{detId.layer()+layerIdxOffset_};
    int iradiusAbs{detId.iradiusAbs()};
    int iradius{detId.iradius()};
    int ietaAbs{detId.ietaAbs()};
    int ieta{detId.ieta()};
    int iphi{detId.iphi()};

    float phi = iphi * 1.25;         //phi in degrees
    float theta = phi * M_PI / 180;  // phi in radians
    float x = iradiusAbs * cos(theta);
    float y = iradiusAbs * sin(theta);
    GlobalPoint global = rhtools_.getPosition(detId);
    float xglob = global.x();
    float yglob = global.y();
    //hlistXYhits_HECback.at(detId.layer())->Fill(x, y);
    hlistXYhits_HECback[detId.layer()]->Fill(x, y);
    hlistGlobalXYhits_HECback[detId.layer()]->Fill(xglob, yglob);
    hlistXYhits_HECback_polar[detId.layer()]->Fill(phi, iradiusAbs);

    if (DEBUG) {
      std::cout << "layer: " << layer << "detId.layer(): " << detId.layer() << ", iradiusAbs: " << iradiusAbs
                << ", iradius: " << iradius << ", ietaAbs: " << ietaAbs << ", ieta: " << ieta << ", iphi: " << iphi
                << std::endl;
    }
    if (!detId.sipm()) {
      hlistXYhits_HECback_Sipm4mm[detId.layer()]->Fill(x, y);
      hlistGlobalXYhits_HECback_Sipm4mm[detId.layer()]->Fill(xglob, yglob);
    } else {
      hlistXYhits_HECback_Sipm2mm[detId.layer()]->Fill(x, y);
      hlistGlobalXYhits_HECback_Sipm2mm[detId.layer()]->Fill(xglob, yglob);
    }

    // in-time BX info
    uint32_t rawData{hit.sample(inTimeSample).data()};
    bool isTOA{hit.sample(inTimeSample).getToAValid()};
    bool isTDC{hit.sample(inTimeSample).mode()};
    bool isBusy{isTDC && rawData == 0};

    //BX-1 info
    uint32_t rawDataBXm1(hit.sample(inTimeSample - 1).data());
    bool isTDCBXm1(hit.sample(inTimeSample - 1).mode());
    bool isTOAbxm1(hit.sample(inTimeSample - 1).getToAValid());
    bool passBXM1 = false;
    // FIX ME: how to get mipADC from noise map?
    // >>
    // if (thisdigi[bx - 1].mode() || (thisdigi[bx - 1].data() > std::floor(mipADC * 2.5))) {
    //   passBXM1 = true;
    // }
    // <<
    float thrbxm1 = 10;  // FIXME!!!
    if (isTDCBXm1 || rawDataBXm1 > thrbxm1) {
      passBXM1 = true;
    }

    // Find corresponding tileboard and fill vectors
    for (const auto &[tb_name, tb] : HGCTileBoards::tb_map) {
      // check if tile is in the same layer
      if (layer != tb.plane) {
        continue;
      }
      if (DEBUG) {
        std::cout << "tb_name: " << tb_name << std::endl;
        std::cout << "layer = " << tb.plane << ", tb.irmin = " << tb.irmin << ", tb.irmax = " << tb.irmax << std::endl;
      }
      // check if tile is in the same radius range
      if ((iradius>tb.irmin)&&(iradius<tb.irmax)) {
        if (!isBusy){
          if (isTDC) { b_tboard_TdcHits_.at(tb_name) += 1; }
          if (isTOA) { b_tboard_ToaHits_.at(tb_name) += 1; }
          if (!isTDC) { b_tboard_AdcHits_.at(tb_name) += 1; }          
        }
        break;
      }
    }

    int layidx = layer; //-1;
    if (!isBusy) {
        if (isTDC) {
            layerTdcHits_[layidx]+=1;
            tileTdcHits_[layidx][iradiusAbs]+=1;
        }
        if (isTOA) {  
            layerToaHits_[layidx]+=1;
            tileToaHits_[layidx][iradiusAbs]+=1;
        }
        if (!isTDC) { 
            layerAdcHits_[layidx]+=1;
            tileAdcHits_[layidx][iradiusAbs]+=1;
        }
    }
    
    // Get electronics information if possible
    uint32_t rawEleID = getRawEleIdFromMap(detId);
    // if rawEleID==0, it means the detId was not found in the map. Skip this detId for now. //FIXME
    if (rawEleID == 0) {
      hlistXYhits_HECback_detIdMiss[detId.layer()]->Fill(x, y);
      hlistGlobalXYhits_HECback_detIdMiss[detId.layer()]->Fill(xglob, yglob);
      continue;
    }

    // Get hits per ECON-D
    HGCalElectronicsId eleID(rawEleID);
    //printEleIdBitValues(eleID);
    EconDType key = std::make_tuple(eleID.zSide(), eleID.localFEDId(), eleID.captureBlock(), eleID.econdIdx());
    // Update the map: if the key exists, add 1 to the value, otherwise insert it with value 1
    if (EcondInfoMap_.find(key) == EcondInfoMap_.end()) {
      std::cout << "key does not exist in map " << std::get<0>(key)
                                                            << " " << std::get<1>(key) << " "
                                                            << std::get<2>(key) << " "
                                                            << std::get<3>(key) << std::endl;
    }

    if (!isBusy) {
      EcondInfoMap_[key].FillHits(isTOA, passBXM1);
    }
  } // end digiColl

  // Get number of words per ECON-D

  for (const auto &mpair : EcondInfoMap_) {
    EconDType key = mpair.first;
    EcondInfoMap_.at(key).setNumDAQPacketWords();
    if (outputECONDFile) {
      EcondInfoMap_.at(mpair.first).printInfo(ofStreamECONDHitsFile_);
    } else {
      EcondInfoMap_.at(mpair.first).printInfo(std::cout);
    }
    EcondInfoMap_.at(mpair.first).resetHitsInfo();
    std::cout << "after clear " << std::endl;
    EcondInfoMap_.at(mpair.first).printInfo(std::cout);
  }

}

std::map<uint32_t, uint32_t> HGCScintAnalyzer::mapGeoToElectronics(
    const hgcal::HGCalMappingModuleParamHostCollection &modules,
    const hgcal::HGCalMappingCellParamHostCollection &cells,
    bool geo2ele,
    bool sipm) {
  //loop over different modules
  std::map<uint32_t, uint32_t> idmap;
  uint32_t ndups(0);
  printf("\n");
  for (int i = 0; i < modules.view().metadata().size(); i++) {
    auto imod = modules.view()[i];

    if (!imod.valid())
      continue;

    //require match to si or SiPM
    if (sipm != imod.isSiPM())
      continue;

    if (imod.plane() == 0) {
      printf("WARNING: found plane=0 for i1=%d i2=%d siPM=%d @ index=%i\n", imod.i1(), imod.i2(), imod.isSiPM(), i);
      continue;
    }

    //loop over cells in the module
    for (int j = 0; j < cells.view().metadata().size(); j++) {
      auto jcell = cells.view()[j];

      //use only the information for cells which match the module type index
      if (jcell.typeidx() != imod.typeidx())
        continue;

      //require that it's a valid cell
      if (!jcell.valid())
        continue;

      //assert type of sensor
      assert(imod.isSiPM() == jcell.isSiPM());
      // make sure the cell is part of the module and it's not a calibration cell
      if (jcell.t() != 1)
        continue;

      // uint32_t elecid = ::hgcal::mappingtools::getElectronicsId(imod.zside(),
      //                                                         imod.fedid(),
      //                                                         imod.captureblockidx(),
      //                                                         imod.econdidx(),
      //                                                         jcell.chip(),
      //                                                         jcell.half(),
      //                                                         jcell.seq());

      uint32_t elecid = imod.eleid() + jcell.eleid();
      uint32_t geoid(0);

      if (sipm) {
        geoid = ::hgcal::mappingtools::getSiPMDetId(
            imod.zside(), imod.plane(), imod.i2(), imod.celltype(), jcell.i1(), jcell.i2());
        //std::cout << "::hgcal::mappingtools::getSiPMDetId  geoid = " << geoid << std::endl;

      } else {
        // geoid = ::hgcal::mappingtools::getSiDetId(imod.zside(),
        //                                         imod.plane(),
        //                                         imod.i1(),
        //                                         imod.i2(),
        //                                         imod.celltype(),
        //                                         jcell.i1(),
        //                                         jcell.i2());

        geoid = imod.detid() + jcell.detid();
      }

      if (geo2ele) {
        auto it = idmap.find(geoid);
        ndups += (it != idmap.end());
        if (!sipm && it != idmap.end() && imod.plane() <= 26) {
          HGCSiliconDetId detid(geoid);
          printf("WARNING duplicate found for plane=%d u=%d v=%d cellU=%d cellV=%d valid=%d -> detid=0x%x\n",
                 imod.plane(),
                 imod.i1(),
                 imod.i2(),
                 jcell.i1(),
                 jcell.i2(),
                 jcell.valid(),
                 detid.rawId());
        }
      }
      if (!geo2ele) {
        auto it = idmap.find(elecid);
        ndups += (it != idmap.end());
      }

      //map
      idmap[geo2ele ? geoid : elecid] = geo2ele ? elecid : geoid;
    }
  }

  if (ndups > 0) {
    printf("[HGCalMappingESSourceTester][mapGeoToElectronics] found %d duplicates with geo2ele=%d for sipm=%d\n",
           ndups,
           geo2ele,
           sipm);
  }

  return idmap;
}

void HGCScintAnalyzer::printDetIdBitValues(HGCScintillatorDetId &detId) {
  std::cout << "rawId: " << detId.rawId() 
            << ", iphi: " << detId.iphi() 
            << ", ring: " << detId.ring() 
            << ", iradius: " << detId.iradius()
            << ", ieta: " << detId.ieta() 
            << ", layer: " << detId.layer()
            << ", trigger: " << detId.trigger()
            << ", sipm: " << detId.sipm() << std::endl
            << ", zside: " << detId.zside()
            << ", type: " << detId.type() 
            << std::endl;
}

void HGCScintAnalyzer::printEleIdBitValues(HGCalElectronicsId &eleId) {
  std::cout << "raw: " << eleId.raw() 
            << ", zSide: " << eleId.zSide()
            << ", localFEDId: " << (uint32_t)eleId.localFEDId() 
            << ", captureBlock: " << (uint32_t)eleId.captureBlock()
            << ", econdIdx: " << (uint32_t)eleId.econdIdx() 
            << ", econdeRx: " << (uint32_t)eleId.econdeRx()
            << ", halfrocChannel: " << (uint32_t)eleId.halfrocChannel() 
            << ", rocChannel: " << (uint32_t)eleId.rocChannel()
            << ", cmWord: " << (uint32_t)eleId.cmWord() 
            << ", isCM: " << eleId.isCM() << std::endl; // "common mode"? Four dedicated cells for measuring electronics common mode noise?
  //eleId.print();
}

std::vector<float> HGCScintAnalyzer::minMaxRingPerLayer(const HGCalGeometry *geo, int layer, bool compareRing) {
  const auto &validDetIds = geo->getValidDetIds();
  int min_ring = 999;
  int max_ring = -9;
  std::vector<float> ringvals(2);
  for (const auto &didIt : validDetIds) {
    HGCScintillatorDetId detId(didIt.rawId());

    // Choose negative because current module locator set only for negative CE
    if (detId.zside() > 0)
      continue;  // do only one side since symmetrical?

    if (detId.layer() != layer)
      continue;

    int ring;
    if(compareRing){ ring = detId.ring();}
    else{ ring = detId.iphi();}

    if (ring > max_ring) {
      max_ring = ring;
    } else if (ring < min_ring) {
      min_ring = ring;
    }
  
  }
  ringvals[0] = min_ring;
  ringvals[1] = max_ring;
  return ringvals;
}

std::vector<float> HGCScintAnalyzer::minMaxRingPerLayer(std::map<uint32_t, uint32_t> sipm_geo2ele_, int layer, bool compareRing){
  int min_ring = 999;
  int max_ring = -9;
  std::vector<float> ringvals(2);
  for (const auto &pair : sipm_geo2ele_) {
    HGCScintillatorDetId detId(pair.first);
    if (detId.layer() != layer)
      continue;
    if (detId.zside() > 0)
      continue;
    int ring;
    if(compareRing){ ring = detId.ring();}
    else{ ring = detId.iphi();}

    if (ring > max_ring) {
      max_ring = ring;
    } else if (ring < min_ring) {
      min_ring = ring;
    }
  }
  ringvals[0] = min_ring;
  ringvals[1] = max_ring;
  return ringvals;
}

// Fix bits in detId to agree with older geometry version stored in the map
uint32_t HGCScintAnalyzer::getRawEleIdFromMap(HGCScintillatorDetId didIt) {
  if (DEBUG) {
    printDetIdBitValues(didIt);
  }

  // ------------------------------------------

  if (DEBUG) { std::cout << "==== Setting type bit to 0 ====" << std::endl; }
  // [26:27] Tile granularity and type (0 fine divisions of scintillators;
  //                                    1 coarse divisions of type "c";
  //                                    2 coarse divisions of type "m")
  static const int kHGCalTypeOffset = 26;
  static const int kHGCalTypeMask0 = 0xF3FFFFFF;
  uint32_t rawDetId_typeMask0 = didIt.rawId();
  if (DEBUG) { std::cout << "rawDetId_typeMask0 before = " << rawDetId_typeMask0 << std::endl; }
  rawDetId_typeMask0 &= kHGCalTypeMask0;
  if (DEBUG) {
    std::cout << "rawDetId_typeMask0 after = " << rawDetId_typeMask0 << std::endl << std::endl;
  }
  HGCScintillatorDetId detId_fixTypeBit(rawDetId_typeMask0);
  if (DEBUG) {
    printDetIdBitValues(detId_fixTypeBit);
  }

  if (DEBUG) {
    std::cout << "==== Setting sipm bit to 0 ====" << std::endl;
  }
  // [23]    SiPM type (0 for 2mm: 1 for 4mm)
  static const int kHGCalSiPMOffset = 23;
  static const int kHGCalSiPMMask0 = 0xFF7FFFFF;
  uint32_t rawDetId_sipmMask0 = rawDetId_typeMask0;
  if (DEBUG) {
    std::cout << "rawDetId_sipmMask0 before = " << rawDetId_sipmMask0 << std::endl;
  }
  rawDetId_sipmMask0 &= kHGCalSiPMMask0;
  if (DEBUG) {
    std::cout << "rawDetId_sipmMask0 after = " << rawDetId_sipmMask0 << std::endl << std::endl;
  }
  HGCScintillatorDetId detId_fixSipmBit(rawDetId_sipmMask0);
  if (DEBUG) {
    printDetIdBitValues(detId_fixSipmBit);
  }

  if (DEBUG) {
    std::cout << "==== Increase layer bit by 1 ====" << std::endl;
  }
  //  detId: min_layer = 8, max_layer= 21
  //  map: min_layer = 9, max_layer= 22
  //  [17:21] Layer
  static const int kHGCalLayerOffset = 17;
  static const int kHGCalLayerMask = 0x1F;
  static const int kHGCalLayerMask0 = 0xFFC1FFFF;
  uint32_t rawDetId_layerP1 = rawDetId_sipmMask0;
  // Get layer bit: right shift by offset, select significant bits with mask
  int layerBits = (rawDetId_layerP1 >> kHGCalLayerOffset) & kHGCalLayerMask;
  if (DEBUG) {
    std::cout << "layerBits = " << layerBits << std::endl;
  }
  // Check it's correct
  assert(layerBits == detId_fixSipmBit.layer());
  // Increase by +1
  layerBits++;
  if (DEBUG) {
    std::cout << "layerBits+1 = " << layerBits << std::endl;
  }
  rawDetId_layerP1 &= kHGCalLayerMask0;  // set relevant bits to 0
  HGCScintillatorDetId detId_test(rawDetId_layerP1);
  if (DEBUG) {
    std::cout << " detId layer = " << detId_test.layer() << std::endl;
  }
  rawDetId_layerP1 |= ((layerBits & kHGCalLayerMask) << kHGCalLayerOffset);  // OR with new bits
  HGCScintillatorDetId detId_fixLayerBit(rawDetId_layerP1);
  if (DEBUG) {
    printDetIdBitValues(detId_fixLayerBit);
  }

  if (DEBUG) {
    std::cout << "==== Decrease radius bit by 1 ====" << std::endl;
  }
  //  [9:16]  |ring| index (starting from a minimum radius depending on type)
  //  detId: min_ring = 1, max_ring= 42
  //  map: min_ring = 0, max_ring= 41
  static const int kHGCalRadiusOffset = 9;
  static const int kHGCalRadiusMask = 0xFF;
  static const int kHGCalRadiusMask0 = 0xFFFE01FF;
  uint32_t rawDetId_radiusM1 = rawDetId_layerP1;
  // Get layer bit: right shift by offset, select significant bits with mask
  int radiusBits = (rawDetId_radiusM1 >> kHGCalRadiusOffset) & kHGCalRadiusMask;
  if (DEBUG) {
    std::cout << "radiusBits = " << radiusBits << std::endl;
  }
  // Check it's correct
  assert(radiusBits == detId_fixLayerBit.ring());
  // Decrease by -1
  radiusBits--;
  if (DEBUG) {
    std::cout << "radiusBits-1 = " << radiusBits << std::endl;
  }
  rawDetId_radiusM1 &= kHGCalRadiusMask0;  // set relevant bits to 0
  HGCScintillatorDetId detId_testradiusM1(rawDetId_radiusM1);
  if (DEBUG) {
    std::cout << " detId ring = " << detId_testradiusM1.ring() << std::endl;
  }
  rawDetId_radiusM1 |= ((radiusBits & kHGCalRadiusMask) << kHGCalRadiusOffset);  // OR with new bits
  HGCScintillatorDetId detId_fixRadiusBit(rawDetId_radiusM1);
  if (DEBUG){
    printDetIdBitValues(detId_fixRadiusBit);
  }

  if (DEBUG) { 
    std::cout << "==== Get electronics ID from map ====" << std::endl;
  }

  uint32_t rawEleID = 0;
  try {
    if (DEBUG)
      std::cout << " sipm_geo2ele_.at(detId_fixRadiusBit.rawId()) " << sipm_geo2ele_.at(detId_fixRadiusBit.rawId())
                << std::endl;
    rawEleID = sipm_geo2ele_.at(detId_fixRadiusBit.rawId());
    //HGCalElectronicsId eleID(rawEleID);
    
  } catch (...) {
    std::cout << "Could not find a match for detId: " << detId_fixRadiusBit.rawId() << std::endl;
    //continue;
  }

  return rawEleID;
}



// ------------ method called once each job just after ending the event loop  ------------
void HGCScintAnalyzer::endJob() {
  ofStreamECONDHitsFile_.close();
  // please remove this method if not needed
  // for(size_t idx=0; idx<layerTdcHits_.size(); idx++){
  //   int layer = idx+1+layerIdxOffset_;
  //     std::cout << "idx = " << idx
  //               << "layerTdcHits_ = " << h_tdcCountProf_->GetBinContent(idx)
  //               << "layerToaHits_ = " << h_toaCountProf_->GetBinContent(idx)
  //               << "layerAdcHits_ = " << h_adcCountProf_->GetBinContent(idx)
  //               << std::endl;
  // }
  }

  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void HGCScintAnalyzer::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addWithDefaultLabel(desc);
  }

  //define this as a plug-in
  DEFINE_FWK_MODULE(HGCScintAnalyzer);
