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

class TtbarBbbarDiLeptonDoubleRDAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TtbarBbbarDiLeptonDoubleRDAnalyzer(const edm::ParameterSet&);
  ~TtbarBbbarDiLeptonDoubleRDAnalyzer();

  enum sys_e {sys_nom,  sys_jes_u, sys_jes_d, sys_jer_n,sys_jer_u, sys_jer_d,
     sys_mu_u, sys_mu_d, sys_el_u, sys_el_d,
//     sys_mueff_u, sys_mueff_d, sys_eleff_u, sys_eleff_d,
//    sys_btag_u, sys_btag_d,
    nsys_e
  };


private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  typedef std::vector<const cat::Lepton*> LeptonPtrs;

  void selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, TtbarBbbarDiLeptonDoubleRDAnalyzer::sys_e sys) const;
  void selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, TtbarBbbarDiLeptonDoubleRDAnalyzer::sys_e sys) const;
  //  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonCollection& recolep) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonPtrs& recolep, TtbarBbbarDiLeptonDoubleRDAnalyzer::sys_e sys);
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
  edm::EDGetTokenT<vfloat> pdfweightsToken_, scaleupweightsToken_, scaledownweightsToken_;
  edm::EDGetTokenT<float> puweightToken_, puweightUpToken_, puweightDownToken_, topPtWeight_;
  bool isRunH_;
  edm::EDGetTokenT<int> trigTokenMUEL_, trigTokenMUELRunBtoG_, trigTokenMUELRunH_;
  edm::EDGetTokenT<int> trigTokenMUMU_, trigTokenMUMURunBtoG_, trigTokenMUMURunH_;
  edm::EDGetTokenT<int> trigTokenELEL_;

  edm::EDGetTokenT<int> trigTokenMUELDiLepton_, trigTokenMUELRunBtoGDiLepton_, trigTokenMUELRunHDiLepton_;
  edm::EDGetTokenT<int> trigTokenMUMUDiLepton_, trigTokenMUMURunBtoGDiLepton_, trigTokenMUMURunHDiLepton_;
  edm::EDGetTokenT<int> trigTokenELELDiLepton_;

  edm::EDGetTokenT<int> trigTokenMUELSingleLepton_, trigTokenMUELRunBtoGSingleLepton_, trigTokenMUELRunHSingleLepton_;
  edm::EDGetTokenT<int> trigTokenMUELSingleMuon_, trigTokenMUELRunBtoGSingleMuon_, trigTokenMUELRunHSingleMuon_;
  edm::EDGetTokenT<int> trigTokenMUELSingleElectron_, trigTokenMUELRunBtoGSingleElectron_, trigTokenMUELRunHSingleElectron_;
  edm::EDGetTokenT<int> trigTokenMUMUSingleLepton_, trigTokenMUMURunBtoGSingleLepton_, trigTokenMUMURunHSingleLepton_;
  edm::EDGetTokenT<int> trigTokenELELSingleLepton_;

  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;

  edm::EDGetTokenT<cat::GenTopCollection> GenTopToken_;

  edm::EDGetTokenT<int> genTtbarIdToken_;

  std::string uncertaintyTag_, payloadName_;

  TTree * ttree_, * ttree2_;
  TTree * ttree3_, * ttree4_,* ttree5_,* ttree6_;

  TTree * ttree7_, * ttree8_,* ttree9_,* ttree10_;
  TTree * ttree11_;
  //TTree * ttree11_, * ttree12_,* ttree13_,* ttree14_;
  //TTree * ttree15_;
  
  long long b_run, b_event;
  int b_nvertex, b_step, b_channel;
  bool b_step1, b_step2, b_step3, b_step4, b_step5, b_step6, b_filtered;
  float b_tri, b_tri_up, b_tri_dn;
  float b_tri_di, b_tri_single, b_tri_singleM, b_tri_singleE;
  float b_met, b_metphi;
  int b_njet30, b_nbjetL30, b_nbjetM30, b_nbjetT30;

  float b_lep1_pt, b_lep1_eta, b_lep1_phi, b_lep1_RelIso;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi, b_lep2_RelIso;
  float b_ll_pt, b_ll_eta, b_ll_phi, b_ll_m;
  int   b_lep1_q, b_lep2_q;

  //for ttbb
  std::vector<float> b_jets_pt;
  std::vector<float> b_jets_eta;
  std::vector<float> b_jets_phi;
  std::vector<float> b_jets_m;
  std::vector<int>    b_jets_flavor;
  std::vector<float> b_jets_bDiscriminatorCSV;
  std::vector<int>    b_csvd_jetid;


  float b_bbM, b_bbDR;
  //mc
  std::vector<float> b_pdfWeights, b_scaleWeightsUp, b_scaleWeightsDown;
  std::vector<float> b_csvweights, b_btagweightsCSVL,  b_btagweightsCSVM,  b_btagweightsCSVT;
  std::vector<float> b_csvweights2;
//  std::vector<float> b_mvaweights, b_btagweightsMVAL,  b_btagweightsMVAM,  b_btagweightsMVAT;

  float b_weight, b_puweight, b_puweightUp, b_puweightDown;

  int b_genTtbarId;

  //mc
  float  b_topPtWeight;
  ///float  b_lepweight;
  float  b_mueffweight, b_mueffweight_up, b_mueffweight_dn;
  float  b_eleffweight, b_eleffweight_up, b_eleffweight_dn;

  /////////
  float  lepton1_pt       ;//  cms.string("lepton1().Pt()"),
  float  lepton1_eta      ;//  cms.string("lepton1().Eta()"),
  float  lepton1_phi      ;//  cms.string("lepton1().Phi()"),
  float  lepton2_pt       ;//  cms.string("lepton2().Pt()"),
  float  lepton2_eta      ;//  cms.string("lepton2().Eta()"),
  float  lepton2_phi      ;//  cms.string("lepton2().Phi()"),

  bool   allHadronic      ;//  cms.string("allHadronic"),

  bool   semiLeptonicM1     ;//  cms.string("semiLeptonic"),
  bool   semiLeptonic0     ;//  cms.string("semiLeptonic"),
  bool   semiLeptonicP1     ;//  cms.string("semiLeptonic"),

  bool   diLeptonicM1;
  bool   diLeptonic0;
  bool   diLeptonicP1;

  bool   diLeptonicMuoMuo ;//  cms.string("diLeptonicMuoMuo"),
  bool   diLeptonicMuoEle ;//  cms.string("diLeptonicMuoEle"),
  bool   diLeptonicEleEle ;//  cms.string("diLeptonicEleEle"),
  bool   diLeptonicTauMuo ;//  cms.string("diLeptonicTauMuo"),
  bool   diLeptonicTauEle ;//  cms.string("diLeptonicTauEle"),
  bool   diLeptonicTauTau ;//  cms.string("diLeptonicTauTau"),

  int    NbJets20         ;//  cms.string("NbJets30(0)"),
  int    NaddbJets20      ;//  cms.string("NaddbJets30(0)"),

  int    NaddcJets20      ;//  cms.string("NaddcJets30(0)"),

  int    NcJets20         ;//  cms.string("NcJets30(0)"),
  int    NJets20          ;//  cms.string("NJets30"),
  int    NaddJets20       ;//  cms.string("NaddJets30"),

  int    NbJets30         ;//  cms.string("NbJets30(0)"),
  int    NaddbJets30      ;//  cms.string("NaddbJets30(0)"),

  int    NaddcJets30      ;//  cms.string("NaddcJets30(0)"),

  int    NcJets30         ;//  cms.string("NcJets30(0)"),
  int    NJets30          ;//  cms.string("NJets30"),
  int    NaddJets30       ;//  cms.string("NaddJets30"),

  TRandom *r = new TRandom();
  const static int NCutflow = 10;
  std::vector<std::vector<int> > cutflow_;
  bool runOnMC_;
  CSVHelper *myCsvWeight;
  BTagWeightEvaluator csvWeight;
  BTagWeightEvaluator bTagWeightCSVL;
  BTagWeightEvaluator bTagWeightCSVM;
  BTagWeightEvaluator bTagWeightCSVT;
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
TtbarBbbarDiLeptonDoubleRDAnalyzer::TtbarBbbarDiLeptonDoubleRDAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  recoFiltersMCToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFiltersMC"));
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));
  genweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  pdfweightsToken_ = consumes<vfloat>(iConfig.getParameter<edm::InputTag>("pdfweights"));
  scaleupweightsToken_ = consumes<vfloat>(iConfig.getParameter<edm::InputTag>("scaleupweights"));
  scaledownweightsToken_ = consumes<vfloat>(iConfig.getParameter<edm::InputTag>("scaledownweights"));
  topPtWeight_ = consumes<float>(iConfig.getParameter<edm::InputTag>("topPtWeight"));

  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  puweightUpToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweightUp"));
  puweightDownToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweightDown"));

  isRunH_ = iConfig.getParameter<bool>("isRunH");

  trigTokenMUELRunBtoG_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunBtoG"));
  trigTokenMUELRunH_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunH"));
  trigTokenMUEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUEL"));
  trigTokenMUMURunBtoG_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMURunBtoG"));
  trigTokenMUMURunH_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMURunH"));
  trigTokenMUMU_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMU"));
  trigTokenELEL_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigELEL"));

  trigTokenMUELRunBtoGDiLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunBtoGDiLepton"));
  trigTokenMUELRunHDiLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunHDiLepton"));
  trigTokenMUELDiLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELDiLepton"));
  trigTokenMUMURunBtoGDiLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMURunBtoGDiLepton"));
  trigTokenMUMURunHDiLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMURunHDiLepton"));
  trigTokenMUMUDiLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMUDiLepton"));
  trigTokenELELDiLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigELELDiLepton"));

  trigTokenMUELRunBtoGSingleLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunBtoGSingleLepton"));
  trigTokenMUELRunHSingleLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunHSingleLepton"));
  trigTokenMUELSingleLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELSingleLepton"));
  trigTokenMUELRunBtoGSingleElectron_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunBtoGSingleElectron"));
  trigTokenMUELRunHSingleElectron_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunHSingleElectron"));
  trigTokenMUELSingleElectron_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELSingleElectron"));
  trigTokenMUELRunBtoGSingleMuon_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunBtoGSingleMuon"));
  trigTokenMUELRunHSingleMuon_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELRunHSingleMuon"));
  trigTokenMUELSingleMuon_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUELSingleMuon"));

  trigTokenMUMURunBtoGSingleLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMURunBtoGSingleLepton"));
  trigTokenMUMURunHSingleLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMURunHSingleLepton"));
  trigTokenMUMUSingleLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMUMUSingleLepton"));

  trigTokenELELSingleLepton_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigELELSingleLepton"));

//  genTtbarIdToken_, NgenJet30Token_, genTtbarLeptonDecayToken_;
  genTtbarIdToken_  = consumes<int>(iConfig.getParameter<edm::InputTag>("genTtbarId"));

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
  csvWeight.initCSVWeight(false, "csvv2");
  //mvaWeight.initCSVWeight(false, "mva");
  bTagWeightCSVL.init(3, "csvv2", BTagEntry::OP_LOOSE , 1);
  bTagWeightCSVM.init(3, "csvv2", BTagEntry::OP_MEDIUM, 1);
  bTagWeightCSVT.init(3, "csvv2", BTagEntry::OP_TIGHT , 1);
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
  ttree6_ = fs->make<TTree>("nomJER_up", "nom2");
  ttree7_ = fs->make<TTree>("nomJER_dw", "nom2");

  ttree8_ = fs->make<TTree>("nomMu_up", "nom2");
  ttree9_ = fs->make<TTree>("nomMu_dw", "nom2");
  ttree10_ = fs->make<TTree>("nomEl_up", "nom2");
  ttree11_ = fs->make<TTree>("nomEl_dw", "nom2");

  //ttree12_ = fs->make<TTree>("nomMueff_up", "nom2");
  //ttree13_ = fs->make<TTree>("nomMueff_dw", "nom2");
  //ttree14_ = fs->make<TTree>("nomEleff_up", "nom2");
  //ttree15_ = fs->make<TTree>("nomEleff_dw", "nom2");

  book(ttree_);
  book(ttree2_);
  book(ttree3_);
  book(ttree4_);
  book(ttree5_);
  book(ttree6_);

  book(ttree7_);
  book(ttree8_);
  book(ttree9_);
  book(ttree10_);

  book(ttree11_);
  //book(ttree12_);
  //book(ttree13_);
  //book(ttree14_);
  //book(ttree15_);

  for (int i = 0; i < NCutflow; i++) cutflow_.push_back({0,0,0,0});

}

void TtbarBbbarDiLeptonDoubleRDAnalyzer::book(TTree* tree){
  tree->Branch("event", &b_event, "event/I");
  tree->Branch("run", &b_run, "run/I");
  tree->Branch("nvertex", &b_nvertex, "nvertex/I");
  tree->Branch("step", &b_step, "step/I");
  tree->Branch("channel", &b_channel, "channel/I");
  tree->Branch("njet30", &b_njet30, "njet30/I");

  tree->Branch("nbjetL30", &b_nbjetL30, "nbjetL30/I");
  tree->Branch("nbjetM30", &b_nbjetM30, "nbjetM30/I");
  tree->Branch("nbjetT30", &b_nbjetT30, "nbjetT30/I");

  tree->Branch("step1", &b_step1, "step1/O");
  tree->Branch("step2", &b_step2, "step2/O");
  tree->Branch("step3", &b_step3, "step3/O");
  tree->Branch("step4", &b_step4, "step4/O");
  tree->Branch("step5", &b_step5, "step5/O");
  tree->Branch("step6", &b_step6, "step6/O");
  tree->Branch("tri", &b_tri, "tri/F");
  tree->Branch("tri_up", &b_tri_up, "tri_up/F");
  tree->Branch("tri_dn", &b_tri_dn, "tri_dn/F");
  tree->Branch("tri_di", &b_tri_di, "tri_di/F");
  tree->Branch("tri_single", &b_tri_single, "tri_single/F");
  tree->Branch("tri_singleM", &b_tri_singleM, "tri_singleM/F");
  tree->Branch("tri_singleE", &b_tri_singleE, "tri_singleE/F");
  tree->Branch("filtered", &b_filtered, "filtered/O");
  tree->Branch("met", &b_met, "met/F");
  tree->Branch("metphi", &b_metphi, "metphi/F");

  tree->Branch("weight", &b_weight, "weight/F");
  tree->Branch("pdfWeights","std::vector<float>",&b_pdfWeights);
  tree->Branch("scaleWeightsUp","std::vector<float>",&b_scaleWeightsUp);
  tree->Branch("scaleWeightsDown","std::vector<float>",&b_scaleWeightsDown);

  tree->Branch("csvweights","std::vector<float>",&b_csvweights);
  tree->Branch("csvweights2","std::vector<float>",&b_csvweights2);
  tree->Branch("btagweightsCSVL","std::vector<float>",&b_btagweightsCSVL);
  tree->Branch("btagweightsCSVM","std::vector<float>",&b_btagweightsCSVM);
  tree->Branch("btagweightsCSVT","std::vector<float>",&b_btagweightsCSVT);

  /*tree->Branch("mvaweights","std::vector<float>",&b_mvaweights);
  tree->Branch("btagweightsMVAL","std::vector<float>",&b_btagweightsMVAL);
  tree->Branch("btagweightsMVAM","std::vector<float>",&b_btagweightsMVAM);
  tree->Branch("btagweightsMVAT","std::vector<float>",&b_btagweightsMVAT);*/

  tree->Branch("topPtWeight", &b_topPtWeight, "topPtWeight/F");

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

  tree->Branch("genTtbarId", &b_genTtbarId, "genTtbarId/I");



  tree->Branch("lep1_pt", &b_lep1_pt, "lep1_pt/F");
  tree->Branch("lep1_eta", &b_lep1_eta, "lep1_eta/F");
  tree->Branch("lep1_phi", &b_lep1_phi, "lep1_phi/F");
  tree->Branch("lep1_RelIso", &b_lep1_RelIso, "lep1_RelIso/F");
  tree->Branch("lep1_q", &b_lep1_q, "lep1_q/I");
  tree->Branch("lep2_pt", &b_lep2_pt, "lep2_pt/F");
  tree->Branch("lep2_eta", &b_lep2_eta, "lep2_eta/F");
  tree->Branch("lep2_phi", &b_lep2_phi, "lep2_phi/F");
  tree->Branch("lep2_RelIso", &b_lep2_RelIso, "lep2_RelIso/F");
  tree->Branch("lep2_q", &b_lep2_q, "lep2_q/I");

  tree->Branch("ll_pt", &b_ll_pt, "ll_pt/F");
  tree->Branch("ll_eta", &b_ll_eta, "ll_eta/F");
  tree->Branch("ll_phi", &b_ll_phi, "ll_phi/F");
  tree->Branch("ll_m", &b_ll_m, "ll_m/F");

  tree->Branch("jets_pt","std::vector<float>",&b_jets_pt);
  tree->Branch("jets_eta","std::vector<float>",&b_jets_eta);
  tree->Branch("jets_phi","std::vector<float>",&b_jets_phi);
  tree->Branch("jets_m","std::vector<float>",&b_jets_m);
  tree->Branch("jets_flavor","std::vector<int>",&b_jets_flavor);
  tree->Branch("jets_bDiscriminatorCSV","std::vector<float>",&b_jets_bDiscriminatorCSV);
  tree->Branch("csvd_jetid","std::vector<int>",&b_csvd_jetid);

  tree->Branch("bbM", &b_bbM, "bbM/F");
  tree->Branch("bbDR", &b_bbDR, "bbDR/F");

/////////////////////////////
  tree->Branch("lepton1_pt",    &lepton1_pt   , "lepton1_pt/F");
  tree->Branch("lepton1_eta",   &lepton1_eta  , "lepton1_eta/F");
  tree->Branch("lepton1_phi",   &lepton1_phi  , "lepton1_phi/F");
  tree->Branch("lepton2_pt",    &lepton2_pt   , "lepton2_pt/F");
  tree->Branch("lepton2_eta",   &lepton2_eta  , "lepton2_eta/F");
  tree->Branch("lepton2_phi",   &lepton2_phi  , "lepton2_phi/F");

  tree->Branch("allHadronic",      &allHadronic        , "allHadronic/O");
  tree->Branch("semiLeptonicM1",     &semiLeptonicM1       , "semiLeptonicM1/O");
  tree->Branch("semiLeptonic0",     &semiLeptonic0       , "semiLeptonic0/O");
  tree->Branch("semiLeptonicP1",     &semiLeptonicP1       , "semiLeptonicP1/O");

  tree->Branch("diLeptonicM1", &diLeptonicM1   , "diLeptonicM1/O");
  tree->Branch("diLeptonic0", &diLeptonic0   , "diLeptonic0/O");
  tree->Branch("diLeptonicP1", &diLeptonicP1   , "diLeptonicP1/O");

  tree->Branch("diLeptonicMuoMuo", &diLeptonicMuoMuo   , "diLeptonicMuoMuo/O");
  tree->Branch("diLeptonicMuoEle", &diLeptonicMuoEle   , "diLeptonicMuoEle/O");
  tree->Branch("diLeptonicEleEle", &diLeptonicEleEle   , "diLeptonicEleEle/O");
  tree->Branch("diLeptonicTauMuo", &diLeptonicTauMuo   , "diLeptonicTauMuo/O");
  tree->Branch("diLeptonicTauEle", &diLeptonicTauEle   , "diLeptonicTauEle/O");
  tree->Branch("diLeptonicTauTau", &diLeptonicTauTau   , "diLeptonicTauTau/O");

  tree->Branch("NJets20",         &NJets20        , "NJets20/I");
  tree->Branch("NbJets20",        &NbJets20       , "NbJets20/I");
  tree->Branch("NcJets20",        &NcJets20       , "NcJets20/I");

  tree->Branch("NaddJets20",      &NaddJets20     , "NaddJets20/I");
  tree->Branch("NaddbJets20",     &NaddbJets20    , "NaddbJets20/I");
  tree->Branch("NaddcJets20",     &NaddcJets20    , "NaddcJets20/I");

  tree->Branch("NJets30",         &NJets30        , "NJets30/I");
  tree->Branch("NbJets30",        &NbJets30       , "NbJets30/I");
  tree->Branch("NcJets30",        &NcJets30       , "NcJets30/I");

  tree->Branch("NaddJets30",      &NaddJets30     , "NaddJets30/I");
  tree->Branch("NaddbJets30",     &NaddbJets30    , "NaddbJets30/I");
  tree->Branch("NaddcJets30",     &NaddcJets30    , "NaddcJets30/I");


}

TtbarBbbarDiLeptonDoubleRDAnalyzer::~TtbarBbbarDiLeptonDoubleRDAnalyzer()
{
  cout <<"cut flow         emu         ee         mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<"step "<< i << " "<< cutflow_[i][0] <<  " "<< cutflow_[i][1] << " " << cutflow_[i][2] << " " << cutflow_[i][3]<< endl;
  }
}

void TtbarBbbarDiLeptonDoubleRDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runOnMC_ = !iEvent.isRealData();
  resetBr();
  ////
  cutflow_[0][b_channel]++;

  b_run = iEvent.id().run();
  b_event = iEvent.id().event();

  ////////////
  edm::Handle<cat::GenTopCollection> genTop;
  //std::cout << "!!reading genTop" << std::endl;
  if ( iEvent.getByToken(GenTopToken_, genTop)){
    //for (auto& genTop : genTops) {

    //std::cout << "reading genTop" << std::endl;
    lepton1_pt  =genTop->at(0).lepton1().Pt();
    lepton1_eta =genTop->at(0).lepton1().Eta();
    lepton1_phi =genTop->at(0).lepton1().Phi();
    lepton2_pt  =genTop->at(0).lepton2().Pt();
    lepton2_eta =genTop->at(0).lepton2().Eta();
    lepton2_phi =genTop->at(0).lepton2().Phi();

    allHadronic       =genTop->at(0).allHadronic();
    semiLeptonicM1    =genTop->at(0).semiLeptonic(-1);
    semiLeptonic0     =genTop->at(0).semiLeptonic(0);
    semiLeptonicP1    =genTop->at(0).semiLeptonic(1);

    diLeptonicM1 =genTop->at(0).diLeptonic(-1);
    diLeptonic0 =genTop->at(0).diLeptonic(0);
    diLeptonicP1 =genTop->at(0).diLeptonic(1);

    diLeptonicMuoMuo =genTop->at(0).diLeptonicMuoMuo();
    diLeptonicMuoEle =genTop->at(0).diLeptonicMuoEle();
    diLeptonicEleEle =genTop->at(0).diLeptonicEleEle();
    diLeptonicTauMuo =genTop->at(0).diLeptonicTauMuo();
    diLeptonicTauEle =genTop->at(0).diLeptonicTauEle();
    diLeptonicTauTau =genTop->at(0).diLeptonicTauTau();


    NJets30          =genTop->at(0).NJets3024();
    NbJets30         =genTop->at(0).NbJets3024();
    NcJets30         =genTop->at(0).NcJets3024();

    NaddJets30       =genTop->at(0).NaddJets3024();
    NaddbJets30      =genTop->at(0).NaddbJets3024();
    NaddcJets30      =genTop->at(0).NaddcJets3024();

    NJets20          =genTop->at(0).NJets20();
    NbJets20         =genTop->at(0).NbJets20();
    NcJets20         =genTop->at(0).NcJets20();

    NaddJets20       =genTop->at(0).NaddJets20();
    NaddbJets20      =genTop->at(0).NaddbJets20();
    NaddcJets20      =genTop->at(0).NaddcJets20();




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

    edm::Handle<vfloat> pdfweightsHandle;
    iEvent.getByToken(pdfweightsToken_, pdfweightsHandle);
    std::vector<float> pdfTemp;
    for ( auto& w : *pdfweightsHandle ) {
      pdfTemp.push_back(w);
    }
    float av = 0;
    float sd = 0;
    if (pdfTemp.size()>=102) {
      for (int i=0; i<=99; ++i) {
        float temppdf = pdfTemp.at(i);
        av += temppdf; sd += temppdf*temppdf;
      }
      av /= 100.; sd /= 99.;
      sd -= av*av; sd = sqrt(sd);
      b_pdfWeights.push_back((av+sd)/1.108);
      b_pdfWeights.push_back((av-sd)/0.8925);
      b_pdfWeights.push_back(pdfTemp.at(100));
      b_pdfWeights.push_back(pdfTemp.at(101));
    }
    else {
      b_pdfWeights.push_back(1.);
      b_pdfWeights.push_back(1.);
      b_pdfWeights.push_back(1.);
      b_pdfWeights.push_back(1.);
    }
    //cout << av << " " << sd <<endl;
    //cout << pdfweightsHandle->size() << endl;    
    edm::Handle<vfloat> scaleupweightsHandle;
    iEvent.getByToken(scaleupweightsToken_, scaleupweightsHandle);
    for ( auto& sup : *scaleupweightsHandle ) b_scaleWeightsUp.push_back(sup); 

    edm::Handle<vfloat> scaledownweightsHandle;
    iEvent.getByToken(scaledownweightsToken_, scaledownweightsHandle);
    for ( auto& sdown : *scaledownweightsHandle ) b_scaleWeightsDown.push_back(sdown);

    //////
    edm::Handle<int> genTtbarIdHandle;
    if(iEvent.getByToken(genTtbarIdToken_, genTtbarIdHandle))     b_genTtbarId = *genTtbarIdHandle;
    edm::Handle<float> topPtWeightHandle;
    if (iEvent.getByToken(topPtWeight_, topPtWeightHandle)) {
      b_topPtWeight = *topPtWeightHandle;
    }

  }
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()){ // skip the event if no PV found
    ttree_->Fill();
    return;
  }
  cutflow_[1][b_channel]++;

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
  cutflow_[2][b_channel]++;

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
    double x = r->Uniform()*35.823;
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
    cutflow_[3][b_channel]++;

    LeptonPtrs recolep;
    for ( const auto& x : selMuons ) recolep.push_back(&x);
    for ( const auto& x : selElecs ) recolep.push_back(&x);
    
    sort(recolep.begin(), recolep.end(), [](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
    if (recolep[0]->pt()<25.) {
      if(sys==0) ttree_->Fill();
      continue;
    }
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
    b_tri_di = b_tri_single = b_tri_singleM = b_tri_singleE = 0;
    edm::Handle<int> trigHandle;
    edm::Handle<int> trigHandleDi;
    edm::Handle<int> trigHandleSingle;
    edm::Handle<int> trigHandleSingleM;
    edm::Handle<int> trigHandleSingleE;
    if (!isRunH_) {
      if      ( b_channel == CH_ELEL ) {
        iEvent.getByToken(trigTokenELEL_, trigHandle);
        iEvent.getByToken(trigTokenELELDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenELELSingleLepton_, trigHandleSingle);
      }
      else if ( b_channel == CH_MUMU && !runOnMC_) {
        iEvent.getByToken(trigTokenMUMURunBtoG_, trigHandle);
        iEvent.getByToken(trigTokenMUMURunBtoGDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenMUMURunBtoGSingleLepton_, trigHandleSingle);
      }
      else if ( b_channel == CH_MUMU && runOnMC_) {
        iEvent.getByToken(trigTokenMUMU_, trigHandle);
        iEvent.getByToken(trigTokenMUMUDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenMUMUSingleLepton_, trigHandleSingle);
      }
      else if ( b_channel == CH_MUEL && !runOnMC_) {
        iEvent.getByToken(trigTokenMUELRunBtoG_, trigHandle);
        iEvent.getByToken(trigTokenMUELRunBtoGDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenMUELRunBtoGSingleLepton_, trigHandleSingle);
        iEvent.getByToken(trigTokenMUELRunBtoGSingleMuon_, trigHandleSingleM);
        iEvent.getByToken(trigTokenMUELRunBtoGSingleElectron_, trigHandleSingleE);
        b_tri_singleM = *trigHandleSingleM;
        b_tri_singleE = *trigHandleSingleE;
      }
      else if ( b_channel == CH_MUEL && runOnMC_) {
        iEvent.getByToken(trigTokenMUEL_, trigHandle);
        iEvent.getByToken(trigTokenMUELDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenMUELSingleLepton_, trigHandleSingle);
        iEvent.getByToken(trigTokenMUELSingleMuon_, trigHandleSingleM);
        iEvent.getByToken(trigTokenMUELSingleElectron_, trigHandleSingleE);
        b_tri_singleM = *trigHandleSingleM;
        b_tri_singleE = *trigHandleSingleE;
      }
    }
    else {
      if      ( b_channel == CH_ELEL ) {
        iEvent.getByToken(trigTokenELEL_, trigHandle);
        iEvent.getByToken(trigTokenELELDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenELELSingleLepton_, trigHandleSingle);
      }
      else if ( b_channel == CH_MUMU && !runOnMC_) {
        iEvent.getByToken(trigTokenMUMURunH_, trigHandle);
        iEvent.getByToken(trigTokenMUMURunHDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenMUMURunHSingleLepton_, trigHandleSingle);
      }
      else if ( b_channel == CH_MUMU && runOnMC_) {
        iEvent.getByToken(trigTokenMUMU_, trigHandle);
        iEvent.getByToken(trigTokenMUMUDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenMUMUSingleLepton_, trigHandleSingle);
      }
      else if ( b_channel == CH_MUEL && !runOnMC_) {
        iEvent.getByToken(trigTokenMUELRunH_, trigHandle);
        iEvent.getByToken(trigTokenMUELRunHDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenMUELRunHSingleLepton_, trigHandleSingle);
        iEvent.getByToken(trigTokenMUELRunHSingleMuon_, trigHandleSingleM);
        iEvent.getByToken(trigTokenMUELRunHSingleElectron_, trigHandleSingleE);
        b_tri_singleM = *trigHandleSingleM;
        b_tri_singleE = *trigHandleSingleE;
      }
      else if ( b_channel == CH_MUEL && runOnMC_) {
        iEvent.getByToken(trigTokenMUEL_, trigHandle);
        iEvent.getByToken(trigTokenMUELDiLepton_, trigHandleDi);
        iEvent.getByToken(trigTokenMUELSingleLepton_, trigHandleSingle);
        iEvent.getByToken(trigTokenMUELSingleMuon_, trigHandleSingleM);
        iEvent.getByToken(trigTokenMUELSingleElectron_, trigHandleSingleE);
        b_tri_singleM = *trigHandleSingleM;
        b_tri_singleE = *trigHandleSingleE;
      } 
    }
    b_tri_di = *trigHandleDi;
    b_tri_single = *trigHandleSingle;
    if (b_tri_di>0) b_tri=1;
    else b_tri=0;
    //else if ( b_channel == CH_MUEL && !runOnMC_) iEvent.getByToken(trigTokenMUEL_, trigHandle);
    //else if ( b_channel == CH_MUEL && runOnMC_) iEvent.getByToken(trigTokenMUELMC_, trigHandle);

    if ( *trigHandle != 0 ) {
      if (x<26.933) {
        b_tri_up = computeTrigSF2016RunBtoG(recolep1, recolep2,  1);
        b_tri_dn = computeTrigSF2016RunBtoG(recolep1, recolep2, -1);
      }
      else {
        b_tri_up = computeTrigSF2016RunH(recolep1, recolep2,  1);
        b_tri_dn = computeTrigSF2016RunH(recolep1, recolep2, -1);
      }
    }

    b_lep1_pt = recolep1.pt(); b_lep1_eta = recolep1.eta(); b_lep1_phi = recolep1.phi(); b_lep1_q = recolep1.charge();
    b_lep2_pt = recolep2.pt(); b_lep2_eta = recolep2.eta(); b_lep2_phi = recolep2.phi(); b_lep2_q = recolep2.charge();
    

    //LeptonWeight LepWeight;
    //double sf1 = 1.0;
    //double sf2 = 1.0;
    if (pdgIdSum == 24) {
      b_lep1_RelIso = recolep1.relIso();     b_lep2_RelIso = recolep2.relIso(0.4);
      if (b_lep2_pt > b_lep1_pt) {
        double tempPt = b_lep1_pt;
        double tempEta = b_lep1_eta;
        double tempPhi = b_lep1_phi;
        double tempQ = b_lep1_q;
        double tempRelIso = b_lep1_RelIso;
        b_lep1_pt = b_lep2_pt;
        b_lep1_eta = b_lep2_eta;
        b_lep1_phi = b_lep2_phi;
        b_lep1_q = b_lep2_q;
        b_lep1_RelIso = b_lep2_RelIso;
        b_lep2_pt = tempPt; 
        b_lep2_eta = tempEta;
        b_lep2_phi = tempPhi;
        b_lep2_q = tempQ;
        b_lep2_RelIso = tempRelIso;
      }
    } // emu
    if (pdgIdSum == 22) {
      b_lep1_RelIso = recolep1.relIso();     b_lep2_RelIso = recolep2.relIso();
    } // ee
    if (pdgIdSum == 26) {
      b_lep1_RelIso = recolep1.relIso(0.4);  b_lep2_RelIso = recolep2.relIso(0.4);
    } // mumu
    //if (b_lep2_pt<20) cout << b_lep2_pt << endl;
    if(runOnMC_){ //b_lepweight = getSF(recolep1, sys)*getSF(recolep2, sys);
      if (x<19.677) {
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
    b_step1 = true;
    b_step = 1;
    cutflow_[4][b_channel]++;
    
    if ( (b_channel == CH_MUEL) || ((b_ll_m < 76) || (b_ll_m > 106)) ){
      b_step2 = true;
      b_step = 2;
      cutflow_[5][b_channel]++;
    }
    
    if(sys>0 && b_step1==false)  continue;
    
    //JetCollection&& selectedJets = selectJets(*jets, recolep);
    JetCollection&& selectedJets = selectJets(*jets, recolep, (sys_e)sys);
    JetCollection&& selectedBJetsL = selectBJets(selectedJets,WP_BTAG_CSVv2L);
    JetCollection&& selectedBJetsM = selectBJets(selectedJets,WP_BTAG_CSVv2M);
    JetCollection&& selectedBJetsT = selectBJets(selectedJets,WP_BTAG_CSVv2T);

    
    int idx=0;
    std::map<int,float> mapJetBDiscriminator;


    for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
      float bDisCSV= (float) jet1->bDiscriminator(BTAG_CSVv2);
      //float bDisCSV= (float) jet1->bDiscriminator(BTAG_CSVv2+"AK4PFPuppi");
      int flavor = jet1->partonFlavour();
      mapJetBDiscriminator[idx] = bDisCSV;
      idx++;
      b_jets_pt.push_back(jet1->p4().pt());
      b_jets_eta.push_back(jet1->p4().eta());
      b_jets_phi.push_back(jet1->p4().phi());
      b_jets_m.push_back(jet1->p4().M());
      b_jets_flavor.push_back(flavor);
      b_jets_bDiscriminatorCSV.push_back(bDisCSV);

    }
       
    if (runOnMC_){
      b_csvweights2.push_back(myCsvWeight->getCSVWeight(selectedJets,0));
      b_csvweights.push_back(csvWeight.eventWeight(selectedJets,0));
      //b_mvaweights.push_back(mvaWeight.eventWeight(selectedJets,0));
      for (unsigned int iu=0; iu<18; iu++)
      {
         b_csvweights2.push_back(myCsvWeight->getCSVWeight(selectedJets,iu+7));
         b_csvweights.push_back(csvWeight.eventWeight(selectedJets,iu+1));
         //b_mvaweights.push_back(mvaWeight.eventWeight(selectedJets,iu+1));
      }
      for (unsigned int iu=0; iu<3; iu++)
      {
         b_btagweightsCSVL.push_back(bTagWeightCSVL.eventWeight(selectedJets, iu));
         b_btagweightsCSVM.push_back(bTagWeightCSVM.eventWeight(selectedJets, iu));
         b_btagweightsCSVT.push_back(bTagWeightCSVT.eventWeight(selectedJets, iu));
         //b_btagweightsMVAL.push_back(bTagWeightMVAL.eventWeight(selectedJets, iu));
         //b_btagweightsMVAM.push_back(bTagWeightMVAM.eventWeight(selectedJets, iu));
         //b_btagweightsMVAT.push_back(bTagWeightMVAT.eventWeight(selectedJets, iu));
         
      }
 
    }
    //csvd order
    std::vector<data_t> vecJetBDisc(mapJetBDiscriminator.begin(), mapJetBDiscriminator.end());
    std::sort(vecJetBDisc.begin(), vecJetBDisc.end(), bigger_second<data_t>());
    for ( const auto& x : vecJetBDisc ) b_csvd_jetid.push_back(x.first);

    //mvad order
 
    const auto met = mets->front().p4();
    b_met = met.pt();
    b_metphi = met.phi();
    b_njet30 = selectedJets.size();

    b_nbjetL30 = selectedBJetsL.size();
    b_nbjetM30 = selectedBJetsM.size();
    b_nbjetT30 = selectedBJetsT.size();

 
    
    if ((b_channel == CH_MUEL) || (b_met > 40.)){
      b_step3 = true;
      if (b_step == 2){
        ++b_step;
        cutflow_[6][b_channel]++;
      }
    }
    
    if (selectedJets.size() >3 ){
      b_step4 = true;
      if (b_step == 3){
        ++b_step;
        cutflow_[7][b_channel]++;
      }
      auto addbjet1 = selectedJets.at(b_csvd_jetid.at(2)); 
      auto addbjet2 = selectedJets.at(b_csvd_jetid.at(3)); 
      b_bbM = (addbjet1.p4() + addbjet2.p4()).M();
      b_bbDR = deltaR(addbjet1.p4(), addbjet2.p4());
      for (auto jet1 = selectedJets.begin(), end = selectedJets.end(); jet1 != end; ++jet1){
        if (!jet1->bDiscriminator(BTAG_CSVv2)) cout << "invalid CSV" << endl;;
        //float bDisCSV= (float) jet1->bDiscriminator(BTAG_CSVv2+"AK4PFPuppi");
      }
         
    }
    
    
    if (selectedBJetsM.size() > 1){
      b_step5 = true;
      if (b_step == 4){
        ++b_step;
        cutflow_[8][b_channel]++;
      }
    }
    
    if (selectedBJetsT.size() > 1){
      b_step6 = true;
      if (b_step == 5){
        ++b_step;
        cutflow_[9][b_channel]++;
      }
    }
    
    if(sys==0){
      ttree_->Fill();
      ttree2_->Fill();
    }else if (sys==1) ttree3_->Fill();
    else if (sys==2) ttree4_->Fill();
    else if (sys==3) ttree5_->Fill();
    else if (sys==4) ttree6_->Fill();

    else if (sys==5) ttree7_->Fill();
    else if (sys==6) ttree8_->Fill();
    else if (sys==7) ttree9_->Fill();
    else if (sys==8) ttree10_->Fill();

    else if (sys==9) ttree11_->Fill();
    //else if (sys==10) ttree12_->Fill();
    //else if (sys==11) ttree13_->Fill();
    //else if (sys==12) ttree14_->Fill();
    //else if (sys==13) ttree15_->Fill();

  }
}


void TtbarBbbarDiLeptonDoubleRDAnalyzer::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons,sys_e sys) const
{
  for (auto& m : muons) {
    cat::Muon mu(m);
    if (sys == sys_mu_u) mu.setP4(m.p4() * m.shiftedEnUp());
    if (sys == sys_mu_d) mu.setP4(m.p4() * m.shiftedEnDown());

    if (mu.pt() < 20.) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    //if (!mu.isMediumMuon()) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.relIso(0.4) > 0.15) continue;
    //printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selmuons.push_back(mu);
  }
}

void TtbarBbbarDiLeptonDoubleRDAnalyzer::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const
{
  for (auto& e : elecs) {
    cat::Electron el(e);
    if (sys == sys_el_u) el.setP4(e.p4() * e.shiftedEnUp());
    if (sys == sys_el_d) el.setP4(e.p4() * e.shiftedEnDown());

    if (el.pt() < 20.) continue;
    if ((std::abs(el.scEta()) > 1.4442) && (std::abs(el.scEta()) < 1.566)) continue;
    if (std::abs(el.eta()) > 2.4) continue;
    //if (!el.electronID("cutBasedElectronID-Spring15-50ns-V1-standalone-medium")) continue;
    //if (el.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium") == 0) continue;
    //if (!el.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90")) continue;
    if ( !el.electronID("cutBasedElectronID-Summer16-80X-V1-tight") ) continue;
    //if (!el.passConversionVeto()) continue;
    //if (!el.isPF()) continue;
    //if ( !el.isTrigMVAValid() or !el.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp90") ) continue;
    //if (el.relIso(0.3) > 0.12) continue;


    //printf("electron with pt %4.1f\n", el.pt());
    selelecs.push_back(el);
  }
}


cat::JetCollection TtbarBbbarDiLeptonDoubleRDAnalyzer::selectJets(const cat::JetCollection& jets, const TtbarBbbarDiLeptonDoubleRDAnalyzer::LeptonPtrs& recolep, sys_e sys)
{
  cat::JetCollection seljets;
  for (auto& j : jets) {
    cat::Jet jet(j);
    double scale = 1.0;
    if      ( sys == sys_jes_u ) scale *= j.shiftedEnUp();
    else if ( sys == sys_jes_d ) scale *= j.shiftedEnDown();
    //if (scale != 1) cout << scale << endl;

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

cat::JetCollection TtbarBbbarDiLeptonDoubleRDAnalyzer::selectBJets(const JetCollection& jets, double workingpoint) const
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
void TtbarBbbarDiLeptonDoubleRDAnalyzer::resetBr()
{
    resetBrGEN();
    resetBrReco();
}
void TtbarBbbarDiLeptonDoubleRDAnalyzer::resetBrReco()
{

  b_step = -1; b_channel = 0;
  b_njet30 = 0;
  b_nbjetL30=0, b_nbjetM30 = 0; b_nbjetT30 = 0;
  b_step1 = 0;b_step2 = 0;    b_step3 = 0;b_step4 = 0;b_step5 = 0;b_step6 = 0;
  b_tri = 0; b_tri_di = 0; b_tri_single=0; b_tri_singleM=0; b_tri_singleE=0;
  b_met = -9; b_metphi = -9;

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
void TtbarBbbarDiLeptonDoubleRDAnalyzer::resetBrJets()
{

  b_jets_pt.clear();
  b_jets_eta.clear();
  b_jets_phi.clear();
  b_jets_m.clear();
  b_jets_flavor.clear();


  b_jets_bDiscriminatorCSV.clear();
  b_csvd_jetid.clear();
  b_csvweights.clear();
  b_csvweights2.clear();
  b_btagweightsCSVL.clear();  b_btagweightsCSVM.clear();  b_btagweightsCSVT.clear();

  //b_mvaweights.clear(); b_btagweightsMVAL.clear();  b_btagweightsMVAM.clear();  b_btagweightsMVAT.clear();


}
void TtbarBbbarDiLeptonDoubleRDAnalyzer::resetBrGEN()
{
  b_genTtbarId=0;
  b_pdfWeights.clear();
  b_scaleWeightsUp.clear();
  b_scaleWeightsDown.clear();
  b_topPtWeight = 1.;
  b_weight = 1; 
  b_puweight = 1; b_puweightUp = 1; b_puweightDown =1;
  b_nvertex = 0;
  b_filtered = 0;


  lepton1_pt  =0.0      ;//  cms.string("lepton1().Pt()"),
  lepton1_eta =-9.0     ;//  cms.string("lepton1().Eta()"),
  lepton1_phi =-9.0     ;//  cms.string("lepton1().Phi()"),
  lepton2_pt  =0.0      ;//  cms.string("lepton2().Pt()"),
  lepton2_eta =-9.0     ;//  cms.string("lepton2().Eta()"),
  lepton2_phi =-9.0     ;//  cms.string("lepton2().Phi()"),

  allHadronic      =false;//  cms.string("allHadronic"),
  semiLeptonicM1     =false;
  semiLeptonic0     =false;
  semiLeptonicP1     =false;


  diLeptonicM1 =false;
  diLeptonic0 =false;
  diLeptonicP1 =false;

  diLeptonicMuoMuo =false;//  cms.string("diLeptonicMuoMuo"),
  diLeptonicMuoEle =false;//  cms.string("diLeptonicMuoEle"),
  diLeptonicEleEle =false;//  cms.string("diLeptonicEleEle"),
  diLeptonicTauMuo =false;//  cms.string("diLeptonicTauMuo"),
  diLeptonicTauEle =false;//  cms.string("diLeptonicTauEle"),
  diLeptonicTauTau =false;//  cms.string("diLeptonicTauTau"),

  NJets20          =0;//  cms.string("NJets30"),
  NbJets20         =0;//  cms.string("NbJets30(0)"),
  NcJets20         =0;//  cms.string("NcJets30(0)"),

  NaddJets20       =0;//  cms.string("NaddJets30"),
  NaddbJets20      =0;//  cms.string("NaddbJets30(0)"),
  NaddcJets20      =0;//  cms.string("NaddcJets30(0)"),

  NJets30          =0;//  cms.string("NJets30"),
  NbJets30         =0;//  cms.string("NbJets30(0)"),
  NcJets30         =0;//  cms.string("NcJets30(0)"),

  NaddJets30       =0;//  cms.string("NaddJets30"),
  NaddbJets30      =0;//  cms.string("NaddbJets30(0)"),
  NaddcJets30      =0;//  cms.string("NaddcJets30(0)"),

}
//define this as a plug-in
DEFINE_FWK_MODULE(TtbarBbbarDiLeptonDoubleRDAnalyzer);
