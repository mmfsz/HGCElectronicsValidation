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

// system include files
#include <memory>



HGCScintAnalyzer::HGCScintAnalyzer(const edm::ParameterSet& iConfig): 
    puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"))),
    caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
    digisCEH_( consumes<HGCalDigiCollection>(edm::InputTag("simHGCalUnsuppressedDigis","HEback")) )
{
  // do what ever initialization is needed

    edm::Service<TFileService> fs;

    // TTree
    t_events_ = fs->make<TTree>("Events","Events");
    t_info_ = fs->make<TTree>("Info","Info"); // information that does not change between events

    t_events_->Branch("npu", &b_npu_, "pileup/I");
    for ( const auto &mpair : HGCTileBoards::tb_map ) {
        std::string tb_name = mpair.first;
        TString rs_tb_name = tb_name;
        t_events_->Branch(rs_tb_name+"_tdcHits", &b_tboard_TdcHits_[tb_name], rs_tb_name+"_tdcHits/I");
        t_events_->Branch(rs_tb_name+"_toaHits", &b_tboard_ToaHits_[tb_name], rs_tb_name+"_toaHits/I");
        t_events_->Branch(rs_tb_name+"_adcHits", &b_tboard_AdcHits_[tb_name], rs_tb_name+"_adcHits/I");
        t_events_->Branch(rs_tb_name + "_adcRawData", &b_tboard_AdcRawData_[tb_name], rs_tb_name + "_adcRawData/I");
        t_info_->Branch(rs_tb_name+"_validDetIds", &b_tboard_ValidDetIds_[tb_name], rs_tb_name+"_validDetIds/I");
    }

    // histograms
    //h_validDetIds_=fs->make<TH1F>("validDetIds",";Layer;Number of valid detIDs/layer",50,0.5,50.5);
    h_cellCount_=fs->make<TH1F>("cellcount",";Layer;Number of cells/layer",50,0.5,50.5);
    h_tdcCountProf_=fs->make<TProfile>("tdccount",";Layer;Number of hits/layer",50,0.5,50.5);
    h_toaCountProf_=fs->make<TProfile>("toacount",";Layer;Number of hits/layer",50,0.5,50.5);
    h_adcCountProf_=fs->make<TProfile>("adccount",";Layer;Number of hits/layer",50,0.5,50.5);
    h_adcRawDataProf_ = fs->make<TProfile>("adcRawData", ";Layer;Number of hits/layer", 50, 0.5, 50.5);

    h_adcHitsVsPU_=fs->make<TH2F>("adchitsvspu",";Number of PU interactions; Number of ADC hits",100,100,300,1000,0,4000); 

    h2_tdcCount_=fs->make<TProfile2D>("tdccount2D",";Layer;Ring",50,0.5,50.5,44,-0.5,43.5);
    h2_toaCount_=fs->make<TProfile2D>("toacount2D",";Layer;Ring",50,0.5,50.5,44,-0.5,43.5);
    h2_adcCount_=fs->make<TProfile2D>("adccount2D",";Layer;Ring",50,0.5,50.5,44,-0.5,43.5);
    h2_adcRawData_ = fs->make<TProfile2D>("adcRawData2D", ";Layer;Ring", 50, 0.5, 50.5, 44, -0.5, 43.5);
    h2_adcRawData_entries_ = fs->make<TH2F>("adcRawData2DEntries", ";Layer;Ring", 50, 0.5, 50.5, 44, -0.5, 43.5);

    layerTdcHits_.resize(Nlayers_,0);
    layerToaHits_.resize(Nlayers_,0);
    layerAdcHits_.resize(Nlayers_,0);
    layerAdcRawData_.resize(Nlayers_, 0);
    layerValidDetIds_.resize(Nlayers_,0);

    tileTdcHits_.resize(Nlayers_,std::vector<int>(NiRings_));
    tileToaHits_.resize(Nlayers_,std::vector<int>(NiRings_));
    tileAdcHits_.resize(Nlayers_,std::vector<int>(NiRings_));
    tileAdcRawData_.resize(Nlayers_,std::vector<int>(NiRings_));
    tileValidDetIds_.resize(Nlayers_,std::vector<int>(NiRings_));

    
    for(size_t layIdx=0; layIdx<14; layIdx++){
      int layer_num = layIdx + 34;
      TString title = "Layer" + std::to_string(layer_num) + "_AdcRawData";
      hlist_layerAdcRawData_.push_back(fs->make<TH1F>(title, ";ADC Raw Data;", 100, 0., 100));
    }
    
}

// ------------ destructor  ------------
HGCScintAnalyzer::~HGCScintAnalyzer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    //
    // please remove this method altogether if it would be left empty
}


// ------------ method called once each job just before starting event loop  ------------
void HGCScintAnalyzer::beginJob() {


}

//
// member functions
//

// ------------ method called for each event  ------------
void HGCScintAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  eventsCount_++;
  using namespace edm;
  

  // reset global counters
  std::fill(layerTdcHits_.begin(), layerTdcHits_.end(), 0);
  std::fill(layerToaHits_.begin(), layerToaHits_.end(), 0);
  std::fill(layerAdcHits_.begin(), layerAdcHits_.end(), 0);
  std::fill(layerAdcRawData_.begin(), layerAdcRawData_.end(), 0);

  std::for_each(std::begin(tileTdcHits_), std::end(tileTdcHits_),
   [](auto & v){ std::fill(std::begin(v), std::end(v), 0); });
  std::for_each(std::begin(tileToaHits_), std::end(tileToaHits_),
   [](auto & v){ std::fill(std::begin(v), std::end(v), 0); });
  std::for_each(std::begin(tileAdcHits_), std::end(tileAdcHits_),
   [](auto & v){ std::fill(std::begin(v), std::end(v), 0); });
  std::for_each(std::begin(tileAdcRawData_), std::end(tileAdcRawData_),
   [](auto & v){ std::fill(std::begin(v), std::end(v), 0); });

  for ( const auto &mpair : HGCTileBoards::tb_map ) {
    std::string tb_name = mpair.first;
    b_tboard_TdcHits_.at(tb_name) = 0;
    b_tboard_ToaHits_.at(tb_name) = 0;
    b_tboard_AdcHits_.at(tb_name) = 0;
    b_tboard_AdcRawData_.at(tb_name) = 0;
    b_tboard_ValidDetIds_.at(tb_name) = 0;
  }


  // get detector information 
  edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(caloGeomToken_);
  const HGCalGeometry *geo = static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(DetId::HGCalHSc, ForwardSubdetector::ForwardEmpty));

  // fill validDetId histogram for first event only
  if (eventsCount_==1){
    //std::string tb_name_tmp = "L34_J8";
    const auto &validDetIds = geo->getValidDetIds();
    for(const auto &didIt : validDetIds) {
      HGCScintillatorDetId detId(didIt.rawId());
      if(detId.zside()<0) continue; // CHECK: do only one side since symmetrical? 
      
      // tile info
      int layer = detId.layer() + layerIdxOffset_;  
      int iradius {detId.iradius()};
      int ieta {detId.ieta()};
      int iphi {detId.iphi()};
      h_cellCount_->Fill(layer);

      // find its tileboar
      for ( const auto &[tb_name, tb] : HGCTileBoards::tb_map ) {
        if (layer!=tb.plane) { continue; } // the detIds taken from geo are numbered correctly
        if ((iradius>tb.irmin)&&(iradius<tb.irmax)) {          
          // std::cout << "entered for layer " << tb.plane << ", iradius= " << iradius 
          //       << ", ieta= " << ieta 
          //       << ", iphi= " << iphi << std::endl;
          b_tboard_ValidDetIds_.at(tb_name) += 1;
          break;
        }
      }
    }
    // store values 
    // for ( const auto &[tb_name, tb] : HGCTileBoards::tb_map ) {
    //   tb.nValidDetIds = b_tboard_ValidDetIds_.at(tb_name) / 36;
    // }
  }
  else{ 
    if (eventsCount_%100==0){
      std::cout << "Event " << eventsCount_ << std::endl;
      }
  }

  // analyze digi collections
  edm::Handle<HGCalDigiCollection> digisHandle;
  iEvent.getByToken(digisCEH_, digisHandle);
  analyzeDigis(digisHandle,geo);

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

  // fill profile histograms per layer
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
    h_adcRawDataProf_->Fill(layer, layerAdcRawData_.at(layIdx));
    
    totalADCHits_+=layerAdcHits_.at(layIdx);
  }  
  h_adcHitsVsPU_->Fill(b_npu_,totalADCHits_, 1); 


  // fill 2D tiles histograms
  for(size_t layIdx=0; layIdx<tileTdcHits_.size(); layIdx++){
    for(size_t rinIdx=0; rinIdx<tileTdcHits_.at(0).size(); rinIdx++){
      int layer = layIdx;
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
      h2_adcCount_->Fill(layer, rinIdx, tileAdcHits_.at(layIdx).at(rinIdx));
      //h2_adcRawData_->Fill(layer, rinIdx, tileAdcRawData_.at(layIdx).at(rinIdx));
    }  
  }

  if (DEBUG) {
    for ( const auto &mpair : HGCTileBoards::tb_map ) {
      std::cout << mpair.first;
      std::string tb_name = mpair.first;
      // std::cout << ", b_tboard_TdcHits_ " << b_tboard_TdcHits_.at(tb_name);
      // std::cout << ", b_tboard_ToaHits_ " << b_tboard_ToaHits_.at(tb_name);
      // std::cout << ", b_tboard_AdcHits_ " << b_tboard_AdcHits_.at(tb_name) << std::endl;
    }
  }
   
  // fill tree for each event
  t_events_->Fill();
  if (eventsCount_==1){
    t_info_->Fill(); 
  }

}

void HGCScintAnalyzer::analyzeDigis(edm::Handle<HGCalDigiCollection> &digiColl, const HGCalGeometry *geom)
{
  // check inputs
  if(!digiColl.isValid() || geom==NULL){
 	  std::cout << "HGCScintAnalyzer::analyzeDigis: digicoll not valid or no geom. Returning." << std::endl;
	   return;	
  }

  const int inTimeSample{2}; // in-time BX sample

  for(auto &hit : *digiColl){
    if(hit.size()==0) continue;

    HGCScintillatorDetId detId{hit.id()};

    if(detId.zside()<0) continue; // CHECK: do only one side since symmetrical? 

    // tile info
    int layer{detId.layer()+layerIdxOffset_};
    int iradiusAbs{detId.iradiusAbs()};
    int iradius{detId.iradius()};
    int ietaAbs{detId.ietaAbs()};
    int ieta{detId.ieta()};
    int iphi{detId.iphi()};

    if (DEBUG) {
      std::cout << "layer: " << layer 
                << ", iradiusAbs: " << iradiusAbs
                << ", iradius: " << iradius
                << ", ietaAbs: " << ietaAbs
                << ", ieta: " << ieta
                << ", iphi: " << iphi << std::endl;
    }

    // in-time BX info
    uint32_t rawData{hit.sample(inTimeSample).data()};
    bool isTOA{hit.sample(inTimeSample).getToAValid()};
    bool isTDC{hit.sample(inTimeSample).mode()};
    bool isBusy{isTDC && rawData==0};
    //uint32_t thr{std::floor(mipADC*adcThrMIP_)};
    //bool passThr{isTDC||rawData>thr};
   
    // Find corresponding tileboard and fill vectors
    for ( const auto &[tb_name, tb] : HGCTileBoards::tb_map ) {
      // check if tile is in the same layer
      if (layer!=tb.plane) { continue; }
      if (DEBUG) {
        std::cout << "tb_name: " << tb_name << std::endl;
        std::cout << "layer = " << tb.plane << ", tb.irmin = " << tb.irmin << ", tb.irmax = " << tb.irmax<< std::endl;
      }
      // check if tile is in the same radius range
      if ((iradius>tb.irmin)&&(iradius<tb.irmax)) {
        if (!isBusy){
          if (isTDC) { b_tboard_TdcHits_.at(tb_name) += 1; }
          if (isTOA) { b_tboard_ToaHits_.at(tb_name) += 1; }
          if (!isTDC) { 
            b_tboard_AdcHits_.at(tb_name) += 1;
            b_tboard_AdcRawData_.at(tb_name) += rawData;
          }          
        }
        break;
      }
    }


    int layidx = layer; //-1;
    if (!isBusy) {
        if (isTDC) {
            layerTdcHits_.at(layidx)+=1;
            tileTdcHits_.at(layidx).at(iradiusAbs)+=1;
        }
        if (isTOA) {  
            layerToaHits_.at(layidx)+=1;
            tileToaHits_.at(layidx).at(iradiusAbs)+=1;
        }
        if (!isTDC) { 
            layerAdcHits_.at(layidx)+=1;
            layerAdcRawData_.at(layidx) += rawData;
            tileAdcHits_.at(layidx).at(iradiusAbs)+=1;
            //std::cout << "layidx = " << layidx << ", iradiusAbs = " << iradiusAbs << ", rawData = " << rawData << std::endl;
            //tileAdcRawData_.at(layidx).at(iradiusAbs) += rawData;
            h2_adcRawData_->Fill(layidx, iradiusAbs, rawData);
            h2_adcRawData_entries_->Fill(layidx, iradiusAbs);
            hlist_layerAdcRawData_.at(layidx - 34)->Fill(rawData);
        }
    }   
  }

  // Divide hits per tileboard by 360/10=36 to remove repetition in phi
  //  for ( const auto &[tb_name, tb] : HGCTileBoards::tb_map ) {
  //     std::cout << "Before: " << b_tboard_ToaHits_.at(tb_name) << std::endl;
  //       b_tboard_ToaHits_.at(tb_name) /= 36; 
  //       b_tboard_TdcHits_.at(tb_name) /= 36; 
  //       b_tboard_AdcHits_.at(tb_name) /= 36; 
  //       std::cout << "After: " << b_tboard_ToaHits_.at(tb_name) << std::endl;   
  //   }

}



// ------------ method called once each job just after ending the event loop  ------------
void HGCScintAnalyzer::endJob() {
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
void HGCScintAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
