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
  tree_ = fs->make<TTree>("ntuple","ntuple");
  tree_->Branch("npu", &b_npu_, "pileup/I");
  
  for ( const auto &mpair : HGCTileBoards::tb_map ) {
    std::string tb_name = mpair.first;
    TString rs_tb_name = tb_name;

    tree_->Branch(rs_tb_name+"_tdcHits", &b_boardTdcHits_[tb_name], rs_tb_name+"_tdcHits/I");
    tree_->Branch(rs_tb_name+"_toaHits", &b_boardToaHits_[tb_name], rs_tb_name+"_toaHits/I");
    tree_->Branch(rs_tb_name+"_adcHits", &b_boardAdcHits_[tb_name], rs_tb_name+"_adcHits/I");
    
  }

  // histograms
  h_cellCount_=fs->make<TH1F>("cellcount",";Layer;Number of cells/layer",50,0.5,50.5);
  h_tdcCountProf_=fs->make<TProfile>("tdccount",";Layer;Number of hits/layer",50,0.5,50.5);
  h_toaCountProf_=fs->make<TProfile>("toacount",";Layer;Number of hits/layer",50,0.5,50.5);
  h_adcCountProf_=fs->make<TProfile>("adccount",";Layer;Number of hits/layer",50,0.5,50.5);
  h_adcHitsVsPU_=fs->make<TH2F>("adchitsvspu",";Number of PU interactions; Number of ADC hits",100,100,300,1000,0.5e5,5e5);

  h2_tdcCount_=fs->make<TH2F>("tdccount",";Layer;Ring",50,0.5,50.5,42,0.5,42.5);
  h2_toaCount_=fs->make<TH2F>("toacount",";Layer;Ring",50,0.5,50.5,42,0.5,42.5);
  h2_adcCount_=fs->make<TH2F>("adccount",";Layer;Ring",50,0.5,50.5,42,0.5,42.5);

  layerTdcHits_.resize(50,0);
  layerToaHits_.resize(50,0);
  layerAdcHits_.resize(50,0);

  tileTdcHits_.resize(50,std::vector<int>(42));
  tileToaHits_.resize(50,std::vector<int>(42));
  tileAdcHits_.resize(50,std::vector<int>(42));
  
}

HGCScintAnalyzer::~HGCScintAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void HGCScintAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
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
    // std::fill(b_boardTdcHits_[tb_name].begin(), b_boardTdcHits_[tb_name].end(), 0);
    // std::fill(b_boardToaHits_[tb_name].begin(), b_boardToaHits_[tb_name].end(), 0);
    // std::fill(b_boardAdcHits_[tb_name].begin(), b_boardAdcHits_[tb_name].end(), 0);
    std::string tb_name = mpair.first;
    b_boardTdcHits_[tb_name] = 0;
    b_boardToaHits_[tb_name] = 0;
    b_boardAdcHits_[tb_name] = 0;
    
  }

  // get detector information //FIXME: currently not used
  edm::ESHandle<CaloGeometry> geom = iSetup.getHandle(caloGeomToken_);
  const HGCalGeometry *geo = static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(DetId::HGCalHSc, ForwardSubdetector::ForwardEmpty));

  const auto &validDetIds = geo->getValidDetIds();
  for(const auto &didIt : validDetIds) {
    HGCScintillatorDetId detId(didIt.rawId());
    if(detId.zside()<0) continue; // CHECK: do only one side since symmetrical? 
    int layer{detId.layer()};
    h_cellCount_->Fill(layer+layerIdxOffset_);
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

  // fill profile histograms
  int totalADC{0};
  for(size_t idx=0; idx<layerTdcHits_.size(); idx++){
    int layer = idx+1+layerIdxOffset_;
    if (DEBUG) {
      std::cout << "layerTdcHits_ = " << layerTdcHits_[idx] 
                << "layerToaHits_ = " << layerToaHits_[idx] 
                << "layerAdcHits_ = " << layerAdcHits_[idx] 
                << std::endl;
    }
    h_tdcCountProf_->Fill(layer,layerTdcHits_[idx]);
    h_toaCountProf_->Fill(layer,layerToaHits_[idx]);
    h_adcCountProf_->Fill(layer,layerAdcHits_[idx]);
    totalADC+=layerAdcHits_[idx];
  }  
  h_adcHitsVsPU_->Fill(b_npu_,totalADC, 1);

  // fill 2D tiles histograms
  for(size_t layIdx=0; layIdx<tileTdcHits_.size(); layIdx++){
    for(size_t rinIdx=0; rinIdx<tileTdcHits_[0].size(); rinIdx++){
      int layer = layIdx+1+layerIdxOffset_;
      if (DEBUG) {
        std::cout << "layIdx = " << layIdx 
                  << ", layer = " << layer 
                  << ", rinIdx = " << rinIdx 
                  << ", tileTdcHits_" << tileTdcHits_.at(layIdx).at(rinIdx)
                  << ", tileToaHits_" << tileToaHits_.at(layIdx).at(rinIdx)
                  << std::endl;
      }

      h2_tdcCount_->Fill(layer,tileTdcHits_.at(layIdx).at(rinIdx),1);
      h2_toaCount_->Fill(layer,tileToaHits_.at(layIdx).at(rinIdx),1);
      h2_adcCount_->Fill(layer,tileAdcHits_.at(layIdx).at(rinIdx),1);
    }  
  }

  if (DEBUG) {
    for ( const auto &mpair : HGCTileBoards::tb_map ) {
      std::cout << mpair.first;
      std::string tb_name = mpair.first;
      std::cout << ", b_boardTdcHits_ " << b_boardTdcHits_[tb_name];
      std::cout << ", b_boardToaHits_ " << b_boardToaHits_[tb_name];
      std::cout << ", b_boardAdcHits_ " << b_boardAdcHits_[tb_name] << std::endl;
    }
  }
   
  // fill tree for each event
  tree_->Fill();

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
    int layer{detId.layer()};
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
   
    for ( const auto &[tb_name, tb] : HGCTileBoards::tb_map ) {
      if (layer+layerIdxOffset_!=tb.plane) { continue; }
      if (DEBUG) {
        std::cout << "tb_name: " << tb_name << std::endl;
        std::cout << "layer = " << tb.plane << ", tb.irmin = " << tb.irmin << ", tb.irmax = " << tb.irmax<< std::endl;
      }
      if ((iradius>tb.irmin)&&(iradius<tb.irmax)) {
        if (!isBusy){
          if (isTDC) { b_boardTdcHits_.at(tb_name) += 1; }
          if (isTOA) { b_boardToaHits_.at(tb_name) += 1; }
          if (!isTDC) { b_boardAdcHits_.at(tb_name) += 1; }
        }
        break;
      }

    }

    //global counts for the in-time bunch
    int layidx{layer-1};
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

  }

}


// ------------ method called once each job just before starting event loop  ------------
void HGCScintAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void HGCScintAnalyzer::endJob() {
  // please remove this method if not needed
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
