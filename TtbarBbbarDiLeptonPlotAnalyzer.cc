#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/GenTop.h"
//#include "CATTools/DataFormats/interface/GenWeights.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "CATTools/CatAnalyzer/interface/CSVHelper.h"
#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TTree.h"
#include "TRandom.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;
using namespace cat;

template<class T>
struct bigger_second
: std::binary_function<T,T,bool>
{
   inline bool operator()(const T& lhs, const T& rhs)
   {
      return lhs.second > rhs.second;
   }
};
typedef std::pair<int,float> data_t;
typedef std::vector<float> vfloat;

class TtbarBbbarDiLeptonPlotAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TtbarBbbarDiLeptonPlotAnalyzer(const edm::ParameterSet&);
  ~TtbarBbbarDiLeptonPlotAnalyzer();

  enum sys_e {sys_nom,  sys_jer_n,sys_jer_u, sys_jer_d,
     
//     sys_mueff_u, sys_mueff_d, sys_eleff_u, sys_eleff_d,
//    sys_btag_u, sys_btag_d,
    nsys_e
  };


private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  typedef std::vector<const cat::Lepton*> LeptonPtrs;

  void selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, TtbarBbbarDiLeptonPlotAnalyzer::sys_e sys) const;
  void selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, TtbarBbbarDiLeptonPlotAnalyzer::sys_e sys) const;
  //  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonCollection& recolep) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonPtrs& recolep, TtbarBbbarDiLeptonPlotAnalyzer::sys_e sys);
  cat::JetCollection selectBJets(const cat::JetCollection& jets, double workingpoint) const;
  //float MuonSF(float pt, float eta);
  void book(TTree* tree);

  void resetBr();
  void resetBrReco();
  void resetBrGEN();
  void resetBrJets();

  ScaleFactorEvaluator muonSF_, muonSFGH_, elecSF_;
  float getMuEffSF(const cat::Lepton& p, int sys) const
  {
    const int aid = abs(p.pdgId());
    if ( aid == 13 ) {
      const double pt = p.pt(), eta = p.eta();
      if      ( sys == +1 ) return muonSF_(pt, eta,  1);
      else if ( sys == -1 ) return muonSF_(pt, eta, -1);
      else return muonSF_(pt, eta, 0);
    }
    return 1;
  }
  float getMuEffSFGH(const cat::Lepton& p, int sys) const
  {
    const int aid = abs(p.pdgId());
    if ( aid == 13 ) {
      const double pt = p.pt(), eta = p.eta();
      if      ( sys == +1 ) return muonSFGH_(pt, eta,  1);
      else if ( sys == -1 ) return muonSFGH_(pt, eta, -1);
      else return muonSFGH_(pt, eta, 0);
    }
    return 1;
  }
  float getElEffSF(const cat::Lepton& p, int sys) const
  { 
    const int aid = abs(p.pdgId());
    if ( aid == 11 ) {
      const auto& el = dynamic_cast<const cat::Electron&>(p);
      const double pt = p.pt(), eta = el.scEta();
      if      ( sys == +1 ) return elecSF_(pt, eta,  1);
      else if ( sys == -1 ) return elecSF_(pt, eta, -1);
      else return elecSF_(pt, eta, 0);
    }
    return 1;
  }

/*
  float getSF(const cat::Lepton& p, int sys) const
  {
    const int aid = abs(p.pdgId());
    if ( aid == 13 ) {
      const double pt = p.pt(), aeta = std::abs(p.eta());
      if      ( sys == sys_mueff_u ) return muonSF_(pt, aeta,  1);
      else if ( sys == sys_mueff_d ) return muonSF_(pt, aeta, -1);
      else return muonSF_(pt, aeta, 0);
    }
    else {
      const auto& el = dynamic_cast<const cat::Electron&>(p);
      const double pt = p.pt(), aeta = std::abs(el.scEta());
      if      ( sys == sys_eleff_u ) return elecSF_(pt, aeta,  1);
      else if ( sys == sys_eleff_d ) return elecSF_(pt, aeta, -1);
      else return elecSF_(pt, aeta, 0);
    }
    return 1;
  }
*/
  edm::EDGetTokenT<int> recoFiltersToken_, recoFiltersMCToken_, nGoodVertexToken_, lumiSelectionToken_;
  edm::EDGetTokenT<float> genweightToken_;
  edm::EDGetTokenT<vfloat> scaleupweightsToken_, scaledownweightsToken_;
  edm::EDGetTokenT<float> puweightToken_, puweightUpToken_, puweightDownToken_, topPtWeight_;
  edm::EDGetTokenT<int> trigTokenMUEL_, trigTokenMUMU_, trigTokenELEL_, trigTokenMUELMC_;

  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;

  edm::EDGetTokenT<cat::GenTopCollection> GenTopToken_;


  std::string uncertaintyTag_, payloadName_;

  TTree * ttree_, * ttree2_;
  TTree * ttree3_, * ttree4_,* ttree5_;

  //TTree * ttree11_, * ttree12_,* ttree13_,* ttree14_;
  //TTree * ttree15_;
  
  int b_nvertex, b_step, b_channel;
  bool b_filtered;
  float b_tri, b_tri_up, b_tri_dn;
  float b_met, b_metphi;
  int b_njet30, b_nbjetL30, b_nbjetM30, b_nbjetT30;

  float b_lep1_pt, b_lep1_eta, b_lep1_phi, b_lep1_RelIso;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi, b_lep2_RelIso;
  float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;
  int   b_lep1_q, b_lep2_q;

  //for ttbb
  std::vector<float> b_jets_pt;
  std::vector<float> b_jets_eta;
  std::vector<float> b_jets_bDiscriminatorCSV;
  std::vector<int>    b_csvd_jetid;


  //mc
  std::vector<float> b_scaleWeightsUp, b_scaleWeightsDown;
  std::vector<float> b_csvweights2;
//  std::vector<float> b_mvaweights, b_btagweightsMVAL,  b_btagweightsMVAM,  b_btagweightsMVAT;

  float b_weight, b_puweight, b_puweightUp, b_puweightDown;


  //mc
  ///float  b_lepweight;
  float  b_mueffweight, b_mueffweight_up, b_mueffweight_dn;
  float  b_eleffweight, b_eleffweight_up, b_eleffweight_dn;

  /////////
  float  lepton1_pt       ;//  cms.string("lepton1().Pt()"),
  float  lepton1_eta      ;//  cms.string("lepton1().Eta()"),
  float  lepton2_pt       ;//  cms.string("lepton2().Pt()"),
  float  lepton2_eta      ;//  cms.string("lepton2().Eta()"),

  int    NbJets20         ;//  cms.string("NbJets20(0)"),
  int    NaddbJets20      ;//  cms.string("NaddbJets20(0)"),


  int    NJets20          ;//  cms.string("NJets20"),
  int    NaddJets20       ;//  cms.string("NaddJets20"),

  TRandom *r = new TRandom();
  bool runOnMC_;
  CSVHelper *myCsvWeight;
/*
  BTagWeightEvaluator mvaWeight;
  BTagWeightEvaluator bTagWeightMVAL;
  BTagWeightEvaluator bTagWeightMVAM;
  BTagWeightEvaluator bTagWeightMVAT;
*/
};
//
// constructors and destructor
//
TtbarBbbarDiLeptonPlotAnalyzer::TtbarBbbarDiLeptonPlotAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  recoFiltersMCToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFiltersMC"));
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));
  genweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  scaleupweightsToken_ = consumes<vfloat>(iConfig.getParameter<edm::InputTag>("scaleupweights"));
  scaledownweightsToken_ = consumes<vfloat>(iConfig.getParameter<edm::InputTag>("scaledownweights"));

  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  puweightUpToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweightUp"));
  puweightDownToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweightDown"));

  trigTokenMUEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUEL"));
  trigTokenMUMU_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMU"));
  trigTokenELEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigELEL"));
  trigTokenMUELMC_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELMC"));

//  genTtbarIdToken_, NgenJet30Token_, genTtbarLeptonDecayToken_;

  GenTopToken_     = consumes<cat::GenTopCollection>(iConfig.getParameter<edm::InputTag>("GenTop"));
  //NgenJet30Token_           = consumes<int>(iConfig.getParameter<edm::InputTag>("NgenJet30"));
  //genTtbarLeptonDecayToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("genTtbarLeptonDecay"));

  muonToken_ = consumes<cat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<cat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  const auto elecSFSet = iConfig.getParameter<edm::ParameterSet>("elecSF");
  elecSF_.set(elecSFSet.getParameter<std::vector<double>>("pt_bins"),
              elecSFSet.getParameter<std::vector<double>>("eta_bins"),
              elecSFSet.getParameter<std::vector<double>>("values"),
              elecSFSet.getParameter<std::vector<double>>("errors"));

  const auto muonSFSet = iConfig.getParameter<edm::ParameterSet>("muonSF");
  muonSF_.set(muonSFSet.getParameter<std::vector<double>>("pt_bins"),
              muonSFSet.getParameter<std::vector<double>>("eta_bins"),
              muonSFSet.getParameter<std::vector<double>>("values"),
              muonSFSet.getParameter<std::vector<double>>("errors"));

  const auto muonSFGHSet = iConfig.getParameter<edm::ParameterSet>("muonSFGH");
  muonSFGH_.set(muonSFSet.getParameter<std::vector<double>>("pt_bins"),
              muonSFSet.getParameter<std::vector<double>>("eta_bins"),
              muonSFSet.getParameter<std::vector<double>>("values"),
              muonSFSet.getParameter<std::vector<double>>("errors"));


  myCsvWeight = new CSVHelper();
  //csvWeight.initCSVWeight(true, "csvv2");
  //bTagWeightMVAL.init(3, "mva", BTagEntry::OP_LOOSE , 1);
  //bTagWeightMVAM.init(3, "mva", BTagEntry::OP_MEDIUM, 1);
  //bTagWeightMVAT.init(3, "mva", BTagEntry::OP_TIGHT , 1);

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("nom", "nom");
  ttree2_ = fs->make<TTree>("nom2", "nom2");

  ttree3_ = fs->make<TTree>("nomJES_up", "nom2");
  ttree4_ = fs->make<TTree>("nomJES_dw", "nom2");
  ttree5_ = fs->make<TTree>("nomJER_n", "nom2");

  //ttree12_ = fs->make<TTree>("nomMueff_up", "nom2");
  //ttree13_ = fs->make<TTree>("nomMueff_dw", "nom2");
  //ttree14_ = fs->make<TTree>("nomEleff_up", "nom2");
  //ttree15_ = fs->make<TTree>("nomEleff_dw", "nom2");

  book(ttree_);
  book(ttree2_);
  book(ttree3_);
  book(ttree4_);
  //book(ttree12_);
  //book(ttree13_);
  //book(ttree14_);
  //book(ttree15_);


}

void TtbarBbbarDiLeptonPlotAnalyzer::book(TTree* tree){
  tree->Branch("step", &b_step, "step/I");
  tree->Branch("channel", &b_channel, "channel/I");
  tree->Branch("tri", &b_tri, "tri/F");
  tree->Branch("tri_up", &b_tri_up, "tri_up/F");
  tree->Branch("tri_dn", &b_tri_dn, "tri_dn/F");
  tree->Branch("filtered", &b_filtered, "filtered/O");

  tree->Branch("weight", &b_weight, "weight/F");
  tree->Branch("scaleWeightsUp","std::vector<float>",&b_scaleWeightsUp);
  tree->Branch("scaleWeightsDown","std::vector<float>",&b_scaleWeightsDown);

  tree->Branch("csvweights2","std::vector<float>",&b_csvweights2);

  /*tree->Branch("mvaweights","std::vector<float>",&b_mvaweights);
  tree->Branch("btagweightsMVAL","std::vector<float>",&b_btagweightsMVAL);
  tree->Branch("btagweightsMVAM","std::vector<float>",&b_btagweightsMVAM);
  tree->Branch("btagweightsMVAT","std::vector<float>",&b_btagweightsMVAT);*/


  tree->Branch("puweight", &b_puweight, "puweight/F");
  tree->Branch("puweightUp", &b_puweightUp, "puweightUp/F");
  tree->Branch("puweightDown", &b_puweightDown, "puweightDown/F");

  tree->Branch("mueffweight",    &b_mueffweight,    "mueffweight/F");
  tree->Branch("mueffweight_up", &b_mueffweight_up, "mueffweight_up/F");
  tree->Branch("mueffweight_dn", &b_mueffweight_dn, "mueffweight_dn/F");
  tree->Branch("eleffweight",    &b_eleffweight,    "eleffweight/F");
  tree->Branch("eleffweight_up", &b_eleffweight_up, "eleffweight_up/F");
  tree->Branch("eleffweight_dn", &b_eleffweight_dn, "eleffweight_dn/F");
//  tree->Branch("lepweight", &b_lepweight, "lepweight/F");




  tree->Branch("jets_bDiscriminatorCSV","std::vector<float>",&b_jets_bDiscriminatorCSV);
  tree->Branch("csvd_jetid","std::vector<int>",&b_csvd_jetid);


/////////////////////////////
  tree->Branch("lepton1_pt",    &lepton1_pt   , "lepton1_pt/F");
  tree->Branch("lepton1_eta",   &lepton1_eta  , "lepton1_eta/F");
  tree->Branch("lepton2_pt",    &lepton2_pt   , "lepton2_pt/F");
  tree->Branch("lepton2_eta",   &lepton2_eta  , "lepton2_eta/F");

  tree->Branch("NJets20",         &NJets20        , "NJets20/I");
  tree->Branch("NbJets20",        &NbJets20       , "NbJets20/I");

  tree->Branch("NaddJets20",      &NaddJets20     , "NaddJets20/I");
  tree->Branch("NaddbJets20",     &NaddbJets20    , "NaddbJets20/I");


}

TtbarBbbarDiLeptonPlotAnalyzer::~TtbarBbbarDiLeptonPlotAnalyzer()
{
}

void TtbarBbbarDiLeptonPlotAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();
  resetBr();
  ////


  ////////////
  edm::Handle<cat::GenTopCollection> genTop;
  //std::cout << "!!reading genTop" << std::endl;
  if ( iEvent.getByToken(GenTopToken_, genTop)){
    //for (auto& genTop : genTops) {

    //std::cout << "reading genTop" << std::endl;
    lepton1_pt  =genTop->at(0).lepton1().Pt();
    lepton1_eta =genTop->at(0).lepton1().Eta();
    lepton2_pt  =genTop->at(0).lepton2().Pt();
    lepton2_eta =genTop->at(0).lepton2().Eta();

    NJets20          =genTop->at(0).NJets20();
    NbJets20         =genTop->at(0).NbJets20(0);

    NaddJets20       =genTop->at(0).NaddJets20();
    NaddbJets20      =genTop->at(0).NaddbJets20(0);




  }
  ////////////


  if (runOnMC_){
    edm::Handle<float> puweightHandle;
    iEvent.getByToken(puweightToken_, puweightHandle);
    b_puweight = *puweightHandle;
    edm::Handle<float> puweightHandleUp;
    iEvent.getByToken(puweightUpToken_, puweightHandleUp);
    b_puweightUp = *puweightHandleUp;
    edm::Handle<float> puweightHandleDown;
    iEvent.getByToken(puweightDownToken_, puweightHandleDown);
    b_puweightDown = *puweightHandleDown;

    edm::Handle<float> genweightHandle;
    iEvent.getByToken(genweightToken_, genweightHandle);
    b_weight = *genweightHandle;

    //cout << pdfweightsHandle->size() << endl;    
    edm::Handle<vfloat> scaleupweightsHandle;
    iEvent.getByToken(scaleupweightsToken_, scaleupweightsHandle);
    for ( auto& sup : *scaleupweightsHandle ) b_scaleWeightsUp.push_back(sup); 

    edm::Handle<vfloat> scaledownweightsHandle;
    iEvent.getByToken(scaledownweightsToken_, scaledownweightsHandle);
    for ( auto& sdown : *scaledownweightsHandle ) b_scaleWeightsDown.push_back(sdown);

    //////

  }

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ // skip the event if no PV found
    ttree_->Fill();
    return;
  }

  // const reco::Vertex &PV = vertices->front();
  edm::Handle<int> nGoodVertexHandle;
  iEvent.getByToken(nGoodVertexToken_, nGoodVertexHandle);
  b_nvertex = *nGoodVertexHandle;

  edm::Handle<int> lumiSelectionHandle;
  iEvent.getByToken(lumiSelectionToken_, lumiSelectionHandle);
  if (!runOnMC_ and *lumiSelectionHandle == 0){
    ttree_->Fill();
    return;
  }

  edm::Handle<int> recoFiltersHandle;
  if (!runOnMC_) iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
  else iEvent.getByToken(recoFiltersMCToken_, recoFiltersHandle);
  b_filtered = *recoFiltersHandle == 0 ? false : true;
  // if (!b_filtered){
  //   ttree_->Fill();
  //   return;
  // }

  edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
  edm::Handle<cat::METCollection> mets;            iEvent.getByToken(metToken_, mets);

////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////
  for (int sys = 0; sys < nsys_e; ++sys){
    double x = r->Uniform()*36.773;
    if (sys > 0 && !runOnMC_) break;
    resetBrReco();

    // Find leptons and sort by pT
    cat::MuonCollection selMuons;
    cat::ElectronCollection selElecs;
    selectMuons(*muons, selMuons, (sys_e) sys);
    selectElecs(*electrons, selElecs, (sys_e) sys);
    
    if (selMuons.size()+selElecs.size() < 2){
      if(sys==0) ttree_->Fill();
      continue;
    }

    LeptonPtrs recolep;
    for ( const auto& x : selMuons ) recolep.push_back(&x);
    for ( const auto& x : selElecs ) recolep.push_back(&x);
    
    sort(recolep.begin(), recolep.end(), [](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
    int lep1_idx=0, lep2_idx=1;
    //for emu
    if(std::abs(recolep[1]->pdgId())==11 && std::abs(recolep[0]->pdgId())==13){
      lep1_idx = 1; lep2_idx = 0;
    }
    
    const cat::Lepton& recolep1 = *recolep[lep1_idx];
    const cat::Lepton& recolep2 = *recolep[lep2_idx];
    
    // Determine channel
    const int pdgIdSum = std::abs(recolep1.pdgId()) + std::abs(recolep2.pdgId());
    if (pdgIdSum == 24) b_channel = CH_MUEL; // emu
    if (pdgIdSum == 22) b_channel = CH_ELEL; // ee
    if (pdgIdSum == 26) b_channel = CH_MUMU; // mumu
    
    // Trigger results
    // Scale factors are from AN16-025 (v4) http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_025_v4.pdf
    b_tri = b_tri_dn = b_tri_up = 0;
    edm::Handle<int> trigHandle;
    if      ( b_channel == CH_ELEL ) iEvent.getByToken(trigTokenELEL_, trigHandle);
    else if ( b_channel == CH_MUMU ) iEvent.getByToken(trigTokenMUMU_, trigHandle);
    else if ( b_channel == CH_MUEL) iEvent.getByToken(trigTokenMUEL_, trigHandle);
    //else if ( b_channel == CH_MUEL && !runOnMC_) iEvent.getByToken(trigTokenMUEL_, trigHandle);
    //else if ( b_channel == CH_MUEL && runOnMC_) iEvent.getByToken(trigTokenMUELMC_, trigHandle);

    if ( *trigHandle != 0 ) {
       b_tri = computeTrigSF(recolep1, recolep2);
       b_tri_up = computeTrigSF(recolep1, recolep2,  1);
       b_tri_dn = computeTrigSF(recolep1, recolep2, -1);
    }

    b_lep1_pt = recolep1.pt(); b_lep1_eta = recolep1.eta(); b_lep1_phi = recolep1.phi(); b_lep1_q = recolep1.charge();
    b_lep2_pt = recolep2.pt(); b_lep2_eta = recolep2.eta(); b_lep2_phi = recolep2.phi(); b_lep2_q = recolep2.charge();
    
    //LeptonWeight LepWeight;
    //double sf1 = 1.0;
    //double sf2 = 1.0;
    if (pdgIdSum == 24) {
      b_lep1_RelIso = recolep1.relIso();     b_lep2_RelIso = recolep2.relIso(0.4);
    } // emu
    if (pdgIdSum == 22) {
      b_lep1_RelIso = recolep1.relIso();     b_lep2_RelIso = recolep2.relIso();
    } // ee
    if (pdgIdSum == 26) {
      b_lep1_RelIso = recolep1.relIso(0.4);  b_lep2_RelIso = recolep2.relIso(0.4);
    } // mumu
    if(runOnMC_){ //b_lepweight = getSF(recolep1, sys)*getSF(recolep2, sys);
      if (x<20.235) {
        b_mueffweight    = getMuEffSF(recolep1,  0)*getMuEffSF(recolep2,  0);
        b_mueffweight_up = getMuEffSF(recolep1, +1)*getMuEffSF(recolep2, +1);
        b_mueffweight_dn = getMuEffSF(recolep1, -1)*getMuEffSF(recolep2, -1);
      }
      else {
        b_mueffweight    = getMuEffSFGH(recolep1,  0)*getMuEffSFGH(recolep2,  0);
        b_mueffweight_up = getMuEffSFGH(recolep1, +1)*getMuEffSFGH(recolep2, +1);
        b_mueffweight_dn = getMuEffSFGH(recolep1, -1)*getMuEffSFGH(recolep2, -1);
      }
      b_eleffweight    = getElEffSF(recolep1,  0)*getElEffSF(recolep2,  0);
      b_eleffweight_up = getElEffSF(recolep1, +1)*getElEffSF(recolep2, +1);
      b_eleffweight_dn = getElEffSF(recolep1, -1)*getElEffSF(recolep2, -1);
    }
    
    const auto tlv_ll = recolep1.p4()+recolep2.p4();
    b_ll_pt = tlv_ll.Pt(); b_ll_eta = tlv_ll.Eta(); b_ll_phi = tlv_ll.Phi(); b_ll_m = tlv_ll.M();
    
    //if (b_ll_m > 20. && recolep1.charge() * recolep2.charge() < 0){
    if (b_ll_m < 20. || recolep1.charge() * recolep2.charge() > 0){
      if(sys==0) ttree_->Fill();
      continue;
    }
    b_step = 1;
    
    if ( (b_channel == CH_MUEL) || ((b_ll_m < 76) || (b_ll_m > 106)) ){
      b_step = 2;
    }
    
    if(sys>0 && b_step!=2)  continue;
    
    //JetCollection&& selectedJets = selectJets(*jets, recolep);
    JetCollection&& selectedJets = selectJets(*jets, recolep, (sys_e)sys);
    JetCollection&& selectedBJetsM = selectBJets(selectedJets,WP_BTAG_CSVv2M);

    
    int idx=0;
    std::map<int,float> mapJetBDiscriminator;


    for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
      float bDisCSV= (float) jet1->bDiscriminator(BTAG_CSVv2);
      //float bDisCSV= (float) jet1->bDiscriminator(BTAG_CSVv2+"AK4PFPuppi");
      mapJetBDiscriminator[idx] = bDisCSV;
      idx++;
      b_jets_pt.push_back(jet1->p4().pt());
      b_jets_eta.push_back(jet1->p4().eta());
      b_jets_bDiscriminatorCSV.push_back(bDisCSV);

    }
       
    if (runOnMC_){
      b_csvweights2.push_back(myCsvWeight->getCSVWeight(selectedJets,0));
      //b_mvaweights.push_back(mvaWeight.eventWeight(selectedJets,0));
      for (unsigned int iu=0; iu<18; iu++)
      {
         b_csvweights2.push_back(myCsvWeight->getCSVWeight(selectedJets,iu+7));
         //b_mvaweights.push_back(mvaWeight.eventWeight(selectedJets,iu+1));
      }
 
    }
    //csvd order
    std::vector<data_t> vecJetBDisc(mapJetBDiscriminator.begin(), mapJetBDiscriminator.end());
    std::sort(vecJetBDisc.begin(), vecJetBDisc.end(), bigger_second<data_t>());
    for ( const auto& x : vecJetBDisc ) b_csvd_jetid.push_back(x.first);

    //mvad order
 
    const auto met = mets->front().p4();
    b_met = met.pt();
    b_njet30 = selectedJets.size();

    b_nbjetM30 = selectedBJetsM.size();

 
    
    if ((b_channel == CH_MUEL) || (b_met > 40.)){
      if (b_step == 2){
        ++b_step;
      }
    }
    
    if (sys>0 && b_step!=3) continue;    
    if (selectedJets.size() >3 ){
      if (b_step == 3){
        ++b_step;
      }
      auto addbjet1 = selectedJets.at(b_csvd_jetid.at(2)); 
      auto addbjet2 = selectedJets.at(b_csvd_jetid.at(3)); 
      for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
        if (!jet1->bDiscriminator(BTAG_CSVv2)) cout << "invalid CSV" << endl;;
        //float bDisCSV= (float) jet1->bDiscriminator(BTAG_CSVv2+"AK4PFPuppi");
      }
         
    }
    if (sys>0 && b_step!=4) continue;    
    
    
    if (selectedBJetsM.size() > 1){
      if (b_step == 4){
        ++b_step;
      }
    }
    
    if (sys>0 && b_step!=5) continue;    
    
    if(sys==0){
      ttree_->Fill();
      ttree2_->Fill();
    }else if (sys==1) ttree3_->Fill();
    else if (sys==2) ttree4_->Fill();
    else if (sys==3) ttree5_->Fill();
    //else if (sys==10) ttree12_->Fill();
    //else if (sys==11) ttree13_->Fill();
    //else if (sys==12) ttree14_->Fill();
    //else if (sys==13) ttree15_->Fill();

  }
}


void TtbarBbbarDiLeptonPlotAnalyzer::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons,sys_e sys) const
{
  for (auto& m : muons) {
    cat::Muon mu(m);

    if (mu.pt() < 20.) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    //if (!mu.isMediumMuon()) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.relIso(0.4) > 0.15) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }
}

void TtbarBbbarDiLeptonPlotAnalyzer::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const
{
  for (auto& e : elecs) {
    cat::Electron el(e);

    if (el.pt() < 20.) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > 2.4) continue;
    //if (!el.electronID("cutBasedElectronID-Spring15-50ns-V1-standalone-medium")) continue;
    //if (el.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium") == 0) continue;
    //if (!el.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90")) continue;
    if ( !el.electronID("cutBasedElectronID-Summer16-80X-V1-medium") ) continue;
    //if (!el.passConversionVeto()) continue;
    //if (!el.isPF()) continue;
    //if ( !el.isTrigMVAValid() or !el.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") ) continue;
    //if (el.relIso(0.3) > 0.12) continue;


    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
}


cat::JetCollection TtbarBbbarDiLeptonPlotAnalyzer::selectJets(const cat::JetCollection& jets, const TtbarBbbarDiLeptonPlotAnalyzer::LeptonPtrs& recolep, sys_e sys)
{
  cat::JetCollection seljets;
  for (auto& j : jets) {
    cat::Jet jet(j);
    double scale = 1.0;
    if (scale != 1) cout << scale << endl;

    if ( runOnMC_ ) {
      if      ( sys == sys_jer_u ) scale *= j.smearedResUp();
      else if ( sys == sys_jer_d ) scale *= j.smearedResDown();
      else if ( sys != sys_jer_n ) scale *= j.smearedRes();
    }

    jet.setP4(j.p4()*scale);

    if (jet.pt() < 30.) continue;
    if (std::abs(jet.eta()) > 2.4)  continue;
    if (!jet.LooseId()) continue;

    bool hasOverLap = false;
    for (auto lep : recolep){
      if (deltaR(jet.p4(),lep->p4()) < 0.4) hasOverLap = true;
    }
    if (hasOverLap) continue;
    seljets.push_back(jet);
  }
  return seljets;
}

cat::JetCollection TtbarBbbarDiLeptonPlotAnalyzer::selectBJets(const JetCollection& jets, double workingpoint) const
{
  //https://twiki.cern.ch/twiki/bin/view/CMS/TopBTV
  // 25ns pfCombinedInclusiveSecondaryVertexV2BJetTags :L,M,T 0.605, 0.890, 0.970
  // 76x : 0.460, 0.800, 0.935
  cat::JetCollection selBjets;
  for (auto& jet : jets) {
    if (jet.bDiscriminator(BTAG_CSVv2) < workingpoint) continue;
    selBjets.push_back(jet);
  }
  return selBjets;
}
void TtbarBbbarDiLeptonPlotAnalyzer::resetBr()
{
    resetBrGEN();
    resetBrReco();
}
void TtbarBbbarDiLeptonPlotAnalyzer::resetBrReco()
{

  b_step = -1; b_channel = 0;
  b_njet30 = 0;
  b_nbjetM30 = 0;  b_tri = 0;
  b_met = -9; 

  //b_lepweight = 1;
  b_mueffweight = 1;b_mueffweight_up = 1;b_mueffweight_dn = 1;
  b_eleffweight = 1;b_eleffweight_up = 1;b_eleffweight_dn = 1;

  ///////
  b_lep1_pt = -9;b_lep1_eta = -9;b_lep1_phi = -9; b_lep1_RelIso = -9; b_lep1_q=0;
  b_lep2_pt = -9;b_lep2_eta = -9;b_lep2_phi = -9; b_lep2_RelIso = -9; b_lep2_q=0;
  b_ll_pt = -9;b_ll_eta = -9;b_ll_phi = -9;b_ll_m = -9;
  resetBrJets();
  //resetBrGEN();
  //////
}
void TtbarBbbarDiLeptonPlotAnalyzer::resetBrJets()
{

  b_jets_pt.clear();
  b_jets_eta.clear();


  b_jets_bDiscriminatorCSV.clear();
  b_csvd_jetid.clear();
  b_csvweights2.clear();

  //b_mvaweights.clear(); b_btagweightsMVAL.clear();  b_btagweightsMVAM.clear();  b_btagweightsMVAT.clear();


}
void TtbarBbbarDiLeptonPlotAnalyzer::resetBrGEN()
{
  b_scaleWeightsUp.clear();
  b_scaleWeightsDown.clear();
  b_weight = 1; 
  b_puweight = 1; b_puweightUp = 1; b_puweightDown =1;
  b_nvertex = 0;
  b_filtered = 0;


  lepton1_pt  =0.0      ;//  cms.string("lepton1().Pt()"),
  lepton1_eta =-9.0     ;//  cms.string("lepton1().Eta()"),
  lepton2_pt  =0.0      ;//  cms.string("lepton2().Pt()"),
  lepton2_eta =-9.0     ;//  cms.string("lepton2().Eta()"),

  NJets20          =0;//  cms.string("NJets20"),
  NbJets20         =0;//  cms.string("NbJets20(0)"),

  NaddJets20       =0;//  cms.string("NaddJets20"),
  NaddbJets20      =0;//  cms.string("NaddbJets20(0)"),

}
//define this as a plug-in
DEFINE_FWK_MODULE(TtbarBbbarDiLeptonPlotAnalyzer);
