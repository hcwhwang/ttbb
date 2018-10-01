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

class TtbarBbbarDiLeptonJESAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TtbarBbbarDiLeptonJESAnalyzer(const edm::ParameterSet&);
  ~TtbarBbbarDiLeptonJESAnalyzer();

  enum sys_e {Nom, AbsoluteStat_u, AbsoluteStat_d, AbsoluteScale_u, AbsoluteScale_d, AbsoluteMPFBias_u, AbsoluteMPFBias_d, Fragmentation_u,  Fragmentation_d, SinglePionECAL_u, SinglePionECAL_d, SinglePionHCAL_u, SinglePionHCAL_d, FlavorQCD_u, FlavorQCD_d, TimePtEta_u, TimePtEta_d, RelativeJEREC1_u, RelativeJEREC1_d, RelativeJEREC2_u, RelativeJEREC2_d, RelativeJERHF_u, RelativeJERHF_d, RelativePtBB_u, RelativePtBB_d,  RelativePtEC1_u, RelativePtEC1_d, RelativePtEC2_u, RelativePtEC2_d, RelativePtHF_u, RelativePtHF_d, RelativeBal_u, RelativeBal_d, RelativeFSR_u, RelativeFSR_d, PileUpDataMC_u, PileUpDataMC_d, PileUpPtRef_u, PileUpPtRef_d, PileUpPtBB_u, PileUpPtBB_d, PileUpPtEC1_u, PileUpPtEC1_d, PileUpPtEC2_u, PileUpPtEC2_d, PileUpPtHF_u, PileUpPtHF_d, nsys_e};
//     sys_mueff_u, sys_mueff_d, sys_eleff_u, sys_eleff_d,
//    sys_btag_u, sys_btag_d,

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  typedef std::vector<const cat::Lepton*> LeptonPtrs;

  void selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons, TtbarBbbarDiLeptonJESAnalyzer::sys_e sys) const;
  void selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, TtbarBbbarDiLeptonJESAnalyzer::sys_e sys) const;
  //  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonCollection& recolep) const;
  cat::JetCollection selectJets(const cat::JetCollection& jets, const LeptonPtrs& recolep, TtbarBbbarDiLeptonJESAnalyzer::sys_e sys);
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
  edm::EDGetTokenT<float> puweightToken_;
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


  std::string uncertaintyTag_, payloadName_;

  TTree *ttree_;
  TTree *ttree3_, *ttree4_,*ttree5_,*ttree6_, *ttree7_, *ttree8_,*ttree9_;
  TTree *ttree10_, *ttree11_, *ttree12_, *ttree13_, *ttree14_, *ttree15_, *ttree16_, *ttree17_, *ttree18_, *ttree19_;
  TTree *ttree20_, *ttree21_, *ttree22_, *ttree23_, *ttree24_, *ttree25_, *ttree26_, *ttree27_, *ttree28_, *ttree29_;
  TTree *ttree30_, *ttree31_, *ttree32_, *ttree33_, *ttree34_, *ttree35_, *ttree36_, *ttree37_, *ttree38_, *ttree39_;
  TTree *ttree40_, *ttree41_, *ttree42_, *ttree43_, *ttree44_, *ttree45_, *ttree46_, *ttree47_, *ttree48_;
  //TTree * ttree11_, * ttree12_,* ttree13_,* ttree14_;
  //TTree * ttree15_;
  
  int b_nvertex, b_step, b_channel;
  bool b_filtered;
  float b_tri;
  float b_tri_di, b_tri_single, b_tri_singleM, b_tri_singleE;

  float b_met;
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
  std::vector<float> b_csvweights2;
//  std::vector<float> b_mvaweights, b_btagweightsMVAL,  b_btagweightsMVAM,  b_btagweightsMVAT;

  float b_weight, b_puweight;


  //mc
  ///float  b_lepweight;
  float  b_mueffweight;
  float  b_eleffweight;

  /////////
  float  lepton1_pt       ;//  cms.string("lepton1().Pt()"),
  float  lepton1_eta      ;//  cms.string("lepton1().Eta()"),
  float  lepton2_pt       ;//  cms.string("lepton2().Pt()"),
  float  lepton2_eta      ;//  cms.string("lepton2().Eta()"),

  int    NbJets30         ;//  cms.string("NbJets30(0)"),
  int    NaddbJets30      ;//  cms.string("NaddbJets30(0)"),
  int    NaddcJets30      ;//  cms.string("NaddbJets30(0)"),
  int    NcJets30          ;//  cms.string("NJets30"),


  int    NJets30          ;//  cms.string("NJets30"),
  int    NaddJets30       ;//  cms.string("NaddJets30"),

  TRandom *r = new TRandom();
  TRandom *rTr = new TRandom();
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
TtbarBbbarDiLeptonJESAnalyzer::TtbarBbbarDiLeptonJESAnalyzer(const edm::ParameterSet& iConfig)
{
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  recoFiltersMCToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFiltersMC"));
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  lumiSelectionToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("lumiSelection"));
  genweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));

  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));

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

  ttree3_ = fs->make<TTree>("JES_AbsoluteStat_up", "nom2");
  ttree4_ = fs->make<TTree>("JES_AbsoluteStat_dw", "nom2");
  ttree5_ = fs->make<TTree>("JES_AbsoluteScale_up", "nom2");
  ttree6_ = fs->make<TTree>("JES_AbsoluteScale_dw", "nom2");
  ttree7_ = fs->make<TTree>("JES_AbsoluteMPFBias_up", "nom2");
  ttree8_ = fs->make<TTree>("JES_AbsoluteMPFBias_dw", "nom2");
  ttree9_ = fs->make<TTree>("JES_Fragmentation_up", "nom2");
  ttree10_ = fs->make<TTree>("JES_Fragmentation_dw", "nom2");
  ttree11_ = fs->make<TTree>("JES_SinglePionECAL_up", "nom2");
  ttree12_ = fs->make<TTree>("JES_SinglePionECAL_dw", "nom2");
  ttree13_ = fs->make<TTree>("JES_SinglePionHCAL_up", "nom2");
  ttree14_ = fs->make<TTree>("JES_SinglePionHCAL_dw", "nom2");
  ttree15_ = fs->make<TTree>("JES_FlavorQCD_up", "nom2");
  ttree16_ = fs->make<TTree>("JES_FlavorQCD_dw", "nom2");
  ttree17_ = fs->make<TTree>("JES_TimePtEta_up", "nom2");
  ttree18_ = fs->make<TTree>("JES_TimePtEta_dw", "nom2");
  ttree19_ = fs->make<TTree>("JES_RelativeJEREC1_up", "nom2");
  ttree20_ = fs->make<TTree>("JES_RelativeJEREC1_dw", "nom2");
  ttree21_ = fs->make<TTree>("JES_RelativeJEREC2_up", "nom2");
  ttree22_ = fs->make<TTree>("JES_RelativeJEREC2_dw", "nom2");
  ttree23_ = fs->make<TTree>("JES_RelativeJERHF_up", "nom2");
  ttree24_ = fs->make<TTree>("JES_RelativeJERHF_dw", "nom2");
  ttree25_ = fs->make<TTree>("JES_RelativePtBB_up", "nom2");
  ttree26_ = fs->make<TTree>("JES_RelativePtBB_dw", "nom2");
  ttree27_ = fs->make<TTree>("JES_RelativePtEC1_up", "nom2");
  ttree28_ = fs->make<TTree>("JES_RelativePtEC1_dw", "nom2");
  ttree29_ = fs->make<TTree>("JES_RelativePtEC2_up", "nom2");
  ttree30_ = fs->make<TTree>("JES_RelativePtEC2_dw", "nom2");
  ttree31_ = fs->make<TTree>("JES_RelativePtHF_up", "nom2");
  ttree32_ = fs->make<TTree>("JES_RelativePtHF_dw", "nom2");
  ttree33_ = fs->make<TTree>("JES_RelativeBal_up", "nom2");
  ttree34_ = fs->make<TTree>("JES_RelativeBal_dw", "nom2");
  ttree35_ = fs->make<TTree>("JES_RelativeFSR_up", "nom2");
  ttree36_ = fs->make<TTree>("JES_RelativeFSR_dw", "nom2");
  ttree37_ = fs->make<TTree>("JES_PileUpDataMC_up", "nom2");
  ttree38_ = fs->make<TTree>("JES_PileUpDataMC_dw", "nom2");
  ttree39_ = fs->make<TTree>("JES_PileUpPtRef_up", "nom2");
  ttree40_ = fs->make<TTree>("JES_PileUpPtRef_dw", "nom2");
  ttree41_ = fs->make<TTree>("JES_PileUpPtBB_up", "nom2");
  ttree42_ = fs->make<TTree>("JES_PileUpPtBB_dw", "nom2");
  ttree43_ = fs->make<TTree>("JES_PileUpPtEC1_up", "nom2");
  ttree44_ = fs->make<TTree>("JES_PileUpPtEC1_dw", "nom2");
  ttree45_ = fs->make<TTree>("JES_PileUpPtEC2_up", "nom2");
  ttree46_ = fs->make<TTree>("JES_PileUpPtEC2_dw", "nom2");
  ttree47_ = fs->make<TTree>("JES_PileUpPtHF_up", "nom2");
  ttree48_ = fs->make<TTree>("JES_PileUpPtHF_dw", "nom2");
  //ttree12_ = fs->make<TTree>("nomMueff_up", "nom2");
  //ttree13_ = fs->make<TTree>("nomMueff_dw", "nom2");
  //ttree14_ = fs->make<TTree>("nomEleff_up", "nom2");
  //ttree15_ = fs->make<TTree>("nomEleff_dw", "nom2");

  book(ttree_);
  book(ttree3_);
  book(ttree4_);
  book(ttree5_);
  book(ttree6_);
  book(ttree7_);
  book(ttree8_);
  book(ttree9_);
  book(ttree10_);
  book(ttree11_);
  book(ttree12_);
  book(ttree13_);
  book(ttree14_);
  book(ttree15_);
  book(ttree16_);
  book(ttree17_);
  book(ttree18_);
  book(ttree19_);
  book(ttree20_);
  book(ttree21_);
  book(ttree22_);
  book(ttree23_);
  book(ttree24_);
  book(ttree25_);
  book(ttree26_);
  book(ttree27_);
  book(ttree28_);
  book(ttree29_);
  book(ttree30_);
  book(ttree31_);
  book(ttree32_);
  book(ttree33_);
  book(ttree34_);
  book(ttree35_);
  book(ttree36_);
  book(ttree37_);
  book(ttree38_);
  book(ttree39_);
  book(ttree40_);
  book(ttree41_);
  book(ttree42_);
  book(ttree43_);
  book(ttree44_);
  book(ttree45_);
  book(ttree46_);
  book(ttree47_);
  book(ttree48_);

}

void TtbarBbbarDiLeptonJESAnalyzer::book(TTree* tree){
//
  tree->Branch("step", &b_step, "step/I");
  tree->Branch("channel", &b_channel, "channel/I");


  tree->Branch("tri", &b_tri, "tri/F");
  tree->Branch("tri_di", &b_tri_di, "tri_di/F");
  tree->Branch("tri_single", &b_tri_single, "tri_single/F");
  tree->Branch("tri_singleM", &b_tri_singleM, "tri_singleM/F");
  tree->Branch("tri_singleE", &b_tri_singleE, "tri_singleE/F");

  tree->Branch("filtered", &b_filtered, "filtered/O");

  tree->Branch("weight", &b_weight, "weight/F");

  tree->Branch("csvweights2","std::vector<float>",&b_csvweights2);

  /*tree->Branch("mvaweights","std::vector<float>",&b_mvaweights);
  tree->Branch("btagweightsMVAL","std::vector<float>",&b_btagweightsMVAL);
  tree->Branch("btagweightsMVAM","std::vector<float>",&b_btagweightsMVAM);
  tree->Branch("btagweightsMVAT","std::vector<float>",&b_btagweightsMVAT);*/


  tree->Branch("puweight", &b_puweight, "puweight/F");

  tree->Branch("mueffweight",    &b_mueffweight,    "mueffweight/F");
  tree->Branch("eleffweight",    &b_eleffweight,    "eleffweight/F");
//  tree->Branch("lepweight", &b_lepweight, "lepweight/F");





  tree->Branch("jets_bDiscriminatorCSV","std::vector<float>",&b_jets_bDiscriminatorCSV);
  tree->Branch("csvd_jetid","std::vector<int>",&b_csvd_jetid);


/////////////////////////////
  tree->Branch("lepton1_pt",    &lepton1_pt   , "lepton1_pt/F");
  tree->Branch("lepton1_eta",   &lepton1_eta  , "lepton1_eta/F");
  tree->Branch("lepton2_pt",    &lepton2_pt   , "lepton2_pt/F");
  tree->Branch("lepton2_eta",   &lepton2_eta  , "lepton2_eta/F");

  tree->Branch("NJets30",         &NJets30        , "NJets30/I");
  tree->Branch("NcJets30",        &NbJets30       , "NbJets30/I");
  tree->Branch("NbJets30",        &NbJets30       , "NbJets30/I");

  tree->Branch("NaddJets30",      &NaddJets30     , "NaddJets30/I");
  tree->Branch("NaddbJets30",     &NaddbJets30    , "NaddbJets30/I");
  tree->Branch("NaddcJets30",     &NaddbJets30    , "NaddbJets30/I");


}

TtbarBbbarDiLeptonJESAnalyzer::~TtbarBbbarDiLeptonJESAnalyzer()
{
}
void TtbarBbbarDiLeptonJESAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    NJets30          =genTop->at(0).NJets3024();
    NbJets30         =genTop->at(0).NbJets3024();
    NcJets30         =genTop->at(0).NcJets3024();

    NaddJets30       =genTop->at(0).NaddJets3024();
    NaddbJets30      =genTop->at(0).NaddbJets3024();
    NaddcJets30      =genTop->at(0).NaddcJets3024();


  }
  ////////////


  if (runOnMC_){
    edm::Handle<float> puweightHandle;
    iEvent.getByToken(puweightToken_, puweightHandle);
    b_puweight = *puweightHandle;
    edm::Handle<float> genweightHandle;
    iEvent.getByToken(genweightToken_, genweightHandle);
    b_weight = *genweightHandle;

    //cout << pdfweightsHandle->size() << endl;    
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
    b_tri  = 0;
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

    if ( *trigHandle != 0 ) {
      if (x<26.933) {
        b_tri = computeTrigSF2016RunBtoG(recolep1, recolep2);
      }
      else {
        b_tri = computeTrigSF2016RunH(recolep1, recolep2);
      }

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
      if (x<19.677) {
        b_mueffweight    = getMuEffSF(recolep1,  0)*getMuEffSF(recolep2,  0);
      }
      else {
        b_mueffweight    = getMuEffSFGH(recolep1,  0)*getMuEffSFGH(recolep2,  0);
      }
      b_eleffweight    = getElEffSF(recolep1,  0)*getElEffSF(recolep2,  0);
    }
    
    const auto tlv_ll = recolep1.p4()+recolep2.p4();
    b_ll_m = tlv_ll.M();
    
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
    
    if(sys>0 && b_step!=3)  continue;
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
    if(sys>0 && b_step!=4)  continue;
    
    
    if (selectedBJetsM.size() > 1){
      if (b_step == 4){
        ++b_step;
      }
    }
    if(sys>0 && b_step!=5)  continue;
    
   
    if(sys==0){
      ttree_->Fill();
    }
    else if (sys==1) ttree3_->Fill();
    else if (sys==2) ttree4_->Fill();
    else if (sys==3) ttree5_->Fill();
    else if (sys==4) ttree6_->Fill();

    else if (sys==5) ttree7_->Fill();
    else if (sys==6) ttree8_->Fill();
    else if (sys==7) ttree9_->Fill();
    else if (sys==8) ttree10_->Fill();

    else if (sys==9) ttree11_->Fill();
    else if (sys==10) ttree12_->Fill();
    else if (sys==11) ttree13_->Fill();
    else if (sys==12) ttree14_->Fill();
    else if (sys==13) ttree15_->Fill();
    else if (sys==14) ttree16_->Fill();
    else if (sys==15) ttree17_->Fill();
    else if (sys==16) ttree18_->Fill();
    else if (sys==17) ttree19_->Fill();
    else if (sys==18) ttree20_->Fill();
    else if (sys==19) ttree21_->Fill();
    else if (sys==20) ttree22_->Fill();
    else if (sys==21) ttree23_->Fill();
    else if (sys==22) ttree24_->Fill();
    else if (sys==23) ttree25_->Fill();
    else if (sys==24) ttree26_->Fill();
    else if (sys==25) ttree27_->Fill();
    else if (sys==26) ttree28_->Fill();
    else if (sys==27) ttree29_->Fill();
    else if (sys==28) ttree30_->Fill();
    else if (sys==29) ttree31_->Fill();
    else if (sys==30) ttree32_->Fill();
    else if (sys==31) ttree33_->Fill();
    else if (sys==32) ttree34_->Fill();
    else if (sys==33) ttree35_->Fill();
    else if (sys==34) ttree36_->Fill();
    else if (sys==35) ttree37_->Fill();
    else if (sys==36) ttree38_->Fill();
    else if (sys==37) ttree39_->Fill();
    else if (sys==38) ttree40_->Fill();
    else if (sys==39) ttree41_->Fill();
    else if (sys==40) ttree42_->Fill();
    else if (sys==41) ttree43_->Fill();
    else if (sys==42) ttree44_->Fill();
    else if (sys==43) ttree45_->Fill();
    else if (sys==44) ttree46_->Fill();
    else if (sys==45) ttree47_->Fill();
    else if (sys==46) ttree48_->Fill();
  }
}


void TtbarBbbarDiLeptonJESAnalyzer::selectMuons(const cat::MuonCollection& muons, cat::MuonCollection& selmuons,sys_e sys) const
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

void TtbarBbbarDiLeptonJESAnalyzer::selectElecs(const cat::ElectronCollection& elecs, cat::ElectronCollection& selelecs, sys_e sys) const
{
  for (auto& e : elecs) {
    cat::Electron el(e);

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


cat::JetCollection TtbarBbbarDiLeptonJESAnalyzer::selectJets(const cat::JetCollection& jets, const TtbarBbbarDiLeptonJESAnalyzer::LeptonPtrs& recolep, sys_e sys)
{
  cat::JetCollection seljets;
  for (auto& j : jets) {
    cat::Jet jet(j);
    double scale = 1.0;

    if ( runOnMC_ ) {
      scale *= j.smearedRes();
      if (sys>0) {
        const int nsrc = 23;
        const char* srcnames[nsrc] = {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeFSR", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF"};
        std::vector<JetCorrectionUncertainty*> vsrc(nsrc);
        int jes_sys = (sys-1)/2;
        int jes_sys_ud = (sys-1)%2;
        const char *name = srcnames[jes_sys];
        std::string uncFile = std::string(getenv("CMSSW_BASE"))+"/src/CATTools/CatProducer/data/JEC/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PF.txt";
        JetCorrectorParameters *p = new JetCorrectorParameters(uncFile, name);
        JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
        unc->setJetPt(jet.pt());
        unc->setJetEta(jet.eta());
        double jes = 0.;
        if (jes_sys_ud==0) jes = unc->getUncertainty(true);
        else if(jes_sys_ud==1) jes = unc->getUncertainty(false)*(-1);
        scale *= (1 + jes);
        delete p;
        delete unc;
      }
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

cat::JetCollection TtbarBbbarDiLeptonJESAnalyzer::selectBJets(const JetCollection& jets, double workingpoint) const
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
void TtbarBbbarDiLeptonJESAnalyzer::resetBr()
{
    resetBrGEN();
    resetBrReco();
}
void TtbarBbbarDiLeptonJESAnalyzer::resetBrReco()
{

  b_step = -1; b_channel = 0;
  b_njet30 = 0;
  b_nbjetM30 = 0;
  b_tri = 0;
  b_met = -9;

  //b_lepweight = 1;
  b_mueffweight = 1;
  b_eleffweight = 1;

  ///////
  b_lep1_pt = -9;b_lep1_eta = -9;
  b_lep2_pt = -9;b_lep2_eta = -9;
  b_ll_m = -9;
  resetBrJets();
  //resetBrGEN();
  //////
}
void TtbarBbbarDiLeptonJESAnalyzer::resetBrJets()
{

  b_jets_pt.clear();
  b_jets_eta.clear();


  b_jets_bDiscriminatorCSV.clear();
  b_csvd_jetid.clear();
  b_csvweights2.clear();



}
void TtbarBbbarDiLeptonJESAnalyzer::resetBrGEN()
{


  b_weight = 1; 
  b_puweight = 1;
  b_nvertex = 0;
  b_filtered = 0;


  lepton1_pt  =0.0      ;//  cms.string("lepton1().Pt()"),
  lepton1_eta =-9.0     ;//  cms.string("lepton1().Eta()"),
  lepton2_pt  =0.0      ;//  cms.string("lepton2().Pt()"),
  lepton2_eta =-9.0     ;//  cms.string("lepton2().Eta()"),

  NJets30          =0;//  cms.string("NJets30"),
  NbJets30         =0;//  cms.string("NbJets30(0)"),
  NcJets30         =0;//  cms.string("NbJets30(0)"),

  NaddJets30       =0;//  cms.string("NaddJets30"),
  NaddbJets30      =0;//  cms.string("NaddbJets30(0)"),
  NaddcJets30      =0;//  cms.string("NaddbJets30(0)"),

}
//define this as a plug-in
DEFINE_FWK_MODULE(TtbarBbbarDiLeptonJESAnalyzer);
