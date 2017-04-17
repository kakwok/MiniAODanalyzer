// -*- C++ -*-
//
// Package:    bbtoDijet/MiniAODanalyzer
// Class:      MiniAODanalyzer
// 
/**\class MiniAODanalyzer MiniAODanalyzer.cc bbtoDijet/MiniAODanalyzer/plugins/MiniAODanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ka Hei Martin Kwok
//         Created:  Tue, 14 Mar 2017 21:57:33 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/JetID.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include <map>
#include <vector>
#include <string>

#include "TH1F.h"
#include "TTree.h"
//
// class declaration
//
using namespace edm;
using namespace std;

//class MiniAODanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
class MiniAODanalyzer : public EDAnalyzer  {
   public:
      explicit MiniAODanalyzer(const edm::ParameterSet&);
      ~MiniAODanalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      double  FindBosonPT(edm::View<pat::Jet> jets);
      virtual void clear();

      // ----------member data ---------------------------
      edm::Service<TFileService> fs_;
      void createHistograms();

      bool DEBUG;
      edm::EDGetTokenT<TriggerResults>                       triggerToken_;
      edm::EDGetTokenT<edm::View<pat::Jet> >                     jetToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> >      prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

      TTree* tree;
      std::vector< std::map<std::string,TH1*> > histos_;

      static const size_t kMaxPFJet_           = 25   ;
      int           pfJetCounter_              =  0                                   ;             
      double        pfJetPt_                   [kMaxPFJet_]              ;             
      float         pfJetEta_                  [kMaxPFJet_]               ;             
      float         pfJetPhi_                  [kMaxPFJet_]               ;             
      double        pfJetEt_                   [kMaxPFJet_]              ;             
      double        pfJetEnergy_               [kMaxPFJet_]              ;             
      double        pfJetPx_                   [kMaxPFJet_]              ;             
      double        pfJetPy_                   [kMaxPFJet_]              ;             
      double        pfJetPz_                   [kMaxPFJet_]              ;             
      double        pfJetCSV_                  [kMaxPFJet_]              ;
      
      int passed_DoublePFJetsC160_p026 = 0;
      int passed_DoublePFJetsC100_p014 = 0;
      int passed_DoublePFJetsC172_p026 = 0;
      int passed_DoublePFJetsC112_p014 = 0;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MiniAODanalyzer::MiniAODanalyzer(const edm::ParameterSet& iConfig):
    DEBUG(iConfig.getUntrackedParameter<bool>("DEBUG",false))
{
    triggerToken_             = consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerTag"));
    jetToken_                 = mayConsume<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetTag"));
    prunedGenToken_           = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles"));
    packedGenToken_           = consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedGenParticles"));
    cout <<" Constructin analyzer"<<endl;
}


MiniAODanalyzer::~MiniAODanalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //delete[]  pfJetPt_                      ;
   //delete[]  pfJetEta_                     ;
   //delete[]  pfJetPhi_                     ;
   //delete[]  pfJetEt_                      ;
   //delete[]  pfJetEnergy_                  ;
   //delete[]  pfJetPx_                      ;
   //delete[]  pfJetPy_                      ;
   //delete[]  pfJetPz_                      ;
   //delete[]  pfJetCSV_                     ;

}

void MiniAODanalyzer::clear(){
    std::fill(std::begin(pfJetPt_      ),  std::end( pfJetPt_      ),   0.    );
    std::fill(std::begin(pfJetEta_     ),  std::end( pfJetEta_     ),   0.    );
    std::fill(std::begin(pfJetPhi_     ),  std::end( pfJetPhi_     ),   0.    );
    std::fill(std::begin(pfJetEt_      ),  std::end( pfJetEt_      ),   0.    );
    std::fill(std::begin(pfJetEnergy_  ),  std::end( pfJetEnergy_  ),   0.    );
    std::fill(std::begin(pfJetPx_      ),  std::end( pfJetPx_      ),   0.    );
    std::fill(std::begin(pfJetPy_      ),  std::end( pfJetPy_      ),   0.    );
    std::fill(std::begin(pfJetPz_      ),  std::end( pfJetPz_      ),   0.    );
    std::fill(std::begin(pfJetCSV_     ),  std::end( pfJetCSV_     ),   -1.    );
}

//
// member functions
//

// ------------ method called for each event  ------------
void MiniAODanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    clear();
    using namespace edm;
    using namespace reco;
    using namespace pat;
	TriggerResults tr;
    Handle<TriggerResults> h_trigRes;
	iEvent.getByToken(triggerToken_, h_trigRes);
    tr = *h_trigRes;
	std::vector<string> triggerList;

    edm::Handle<edm::View<pat::Jet> > jetHandle;
    iEvent.getByToken(jetToken_,jetHandle);
    const edm::View<pat::Jet> & jets = *jetHandle;

    // Pruned particles are the one containing "important" stuff
    Handle<edm::View<reco::GenParticle> > prunedGensHandle;
    iEvent.getByToken(prunedGenToken_,prunedGensHandle);
    const edm::View<reco::GenParticle> & prunedGens = *prunedGensHandle;

    // Packed particles are all the status 1, so usable to remake jets
    // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    Handle<edm::View<pat::PackedGenParticle> > packed;
    iEvent.getByToken(packedGenToken_,packed);

    //Lorentz vectors for jets, egamma, muons, leptons, and all objects
    math::XYZTLorentzVectorF pJet;

    std::map<std::string,TH1*>& histo = histos_[0];
    double bosonPT=FindBosonPT(jets);
    histo["bosonPT"]->Fill(bosonPT);


    int gencnt=0;
    double genPt =0;
    //Look for the Higgs/Phi gen pT
    for(edm::View<reco::GenParticle>::const_iterator gen = prunedGens.begin(); gen!=prunedGens.end(); ++gen){
        gencnt++;
        if(DEBUG){
            cout <<" gen["<<gencnt<<"]   | pdgID = " << gen->pdgId() << "  |  pT = " << gen->pt() <<"  |   status = " << gen->status()<<endl;
        }
        // This will pick up the last status of the particle
        if(gen->pdgId()==9900032){
            genPt= gen->pt();
        }
    }
    if(DEBUG){
       cout <<" Final genPT of this event is "<< genPt <<endl;
    }
    histo["genPT"]->Fill(genPt);

    Service<service::TriggerNamesService> tns;
    bool foundNames = tns->getTrigPaths(*h_trigRes, triggerList);
    if (!foundNames) std::cout << "Could not get trigger names!\n";
    if (tr.size()!=triggerList.size()){
		 std::cout << "ERROR: length of names and paths not the same: "   << triggerList.size() << "," << tr.size() << endl;
	}
	// dump trigger list at first event
    for (unsigned int i=0; i< tr.size(); i++) {
	    //std::cout << "["<<i<<"] = " << triggerList[i]<<setw(40) << ": " << (tr[i].accept() ? "Event Passed" : "Event Failed") << endl;
	    if ( !tr[i].accept() == 1 ) continue;
            if( triggerList[i].find("HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160")!=std::string::npos          ){
              if(DEBUG){       std::cout<<triggerList[i]<< " = Passed"<<std::endl;}
                histo["bbPT_C160_p026_online"]->Fill(bosonPT);
                histo["OnlineBtag"]->Fill("DoublePFJetsC160_p026",1);
            }
            if( triggerList[i].find("HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6")!=std::string::npos){
              if(DEBUG){       std::cout<<triggerList[i]<< " = Passed"<<std::endl;}
                histo["bbPT_C100_p014_online"]->Fill(bosonPT);
                histo["OnlineBtag"]->Fill("DoublePFJetsC100_p014",1);
            }
            if( triggerList[i].find("HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172")!=std::string::npos          ){
              if(DEBUG){       std::cout<<triggerList[i]<< " = Passed"<<std::endl;}
                histo["bbPT_C172_p026_online"]->Fill(bosonPT);
                histo["OnlineBtag"]->Fill("DoublePFJetsC172_p026",1);
            }
            if( triggerList[i].find("HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6")!=std::string::npos){
              if(DEBUG){       std::cout<<triggerList[i]<< " = Passed"<<std::endl;}
                histo["bbPT_C112_p014_online"]->Fill(bosonPT);
                histo["OnlineBtag"]->Fill("DoublePFJetsC112_p014",1);
            }
	}
    

    int jetcnt=0;
    passed_DoublePFJetsC160_p026 = 0;
    passed_DoublePFJetsC100_p014 = 0;
    passed_DoublePFJetsC172_p026 = 0;
    passed_DoublePFJetsC112_p014 = 0;

    pfJetCounter_ = jets.size();
    histo["Njets"]->Fill(pfJetCounter_);
    for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){
        jetcnt++;
        if(DEBUG){
            cout <<"jet["<<jetcnt<<"]" <<" | pt = "<<jet->pt()<<" | eta = "<<jet->eta()<< " | phi = " <<jet->phi() <<" | energy = "<< jet->energy()
                 <<" | csv = "<<jet->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") <<endl ;
        }
        double csv =  jet->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        if( fabs(jet->eta()) < 2.4 ){
            if (jet->pt()>160.0 &&  csv >0.78) passed_DoublePFJetsC160_p026++;
            if (jet->pt()>100.0 &&  csv >0.84) passed_DoublePFJetsC100_p014++;
            if (jet->pt()>172.0 &&  csv >0.78) passed_DoublePFJetsC172_p026++;
            if (jet->pt()>112.0 &&  csv >0.84) passed_DoublePFJetsC112_p014++;
        }
         pfJetPt_[jetcnt]     = jet->pt();
         pfJetEta_[jetcnt]    = jet->eta(); 
         pfJetPhi_[jetcnt]    = jet->phi();
         pfJetEt_[jetcnt]     = jet->et();
         pfJetEnergy_[jetcnt] = jet->energy();
         pfJetPx_[jetcnt]     = jet->px();
         pfJetPy_[jetcnt]     = jet->py();
         pfJetPz_[jetcnt]     = jet->pz();
         pfJetCSV_[jetcnt]    = csv;   
    }
    if(DEBUG){
       cout <<" passed_DoublePFJetsC160_p026 ="<< passed_DoublePFJetsC160_p026<<endl; 
       cout <<" passed_DoublePFJetsC100_p014 ="<< passed_DoublePFJetsC100_p014<<endl;
       cout <<" passed_DoublePFJetsC172_p026 ="<< passed_DoublePFJetsC172_p026<<endl;
       cout <<" passed_DoublePFJetsC112_p014 ="<< passed_DoublePFJetsC112_p014<<endl;
    }

    if (passed_DoublePFJetsC160_p026>=2){
        histo["OfflineBtag"]->Fill("DoublePFJetsC160_p026",1);
        histo["bbPT_C160_p026_offline"]->Fill(bosonPT);
    }
    if (passed_DoublePFJetsC100_p014>=2){
        histo["OfflineBtag"]->Fill("DoublePFJetsC100_p014",1);
        histo["bbPT_C100_p014_offline"]->Fill(bosonPT);
    }
    if (passed_DoublePFJetsC172_p026>=2){
        histo["OfflineBtag"]->Fill("DoublePFJetsC172_p026",1);
        histo["bbPT_C172_p026_offline"]->Fill(bosonPT);
    }
    if (passed_DoublePFJetsC112_p014>=2){
        histo["OfflineBtag"]->Fill("DoublePFJetsC112_p014",1);
        histo["bbPT_C112_p014_offline"]->Fill(bosonPT);
    }
    tree->Fill();
    
}

// ------------ method to calculate pT sum of best two b-jets
double MiniAODanalyzer::FindBosonPT(edm::View<pat::Jet> jets){
    float MaxCSV     = 0;
    float SecMaxCSV  = 0;
    int   iMaxCSV    = 0;
    int   iSecMaxCSV = 0;
    int   iJet = 0;
    double bosonPT   =0;
    for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){
        float iCSV = jet->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        if(iCSV>0){
            if (iCSV>MaxCSV){
                // save maxCSV's value before updating
                if(MaxCSV!=0){
                    SecMaxCSV = MaxCSV;
                    iSecMaxCSV= iJet;
                }  
                MaxCSV = iCSV;  //update MaxCSV
                iMaxCSV= iJet;
            }
            else{   
                if(iCSV>SecMaxCSV){
                    SecMaxCSV  = iCSV;   
                    iSecMaxCSV = iJet;
                }
            }
        }
        iJet++;
    }
    if (!(MaxCSV==0 && SecMaxCSV==0)){
        bosonPT = (jets[iMaxCSV].p4() + jets[iSecMaxCSV].p4()).Pt();
    }
    if(DEBUG){
        cout << "Max CSV at ["<<iMaxCSV<<"] = "<< MaxCSV << "  2nd Max CSV at ["<<iSecMaxCSV<<"] = "<< SecMaxCSV<< "      bosonPT =   "<<bosonPT<<endl;
    }

    return bosonPT;
}

// ------------ method called once each job just before starting event loop  ------------
void MiniAODanalyzer::beginJob()
{
    tree = fs_->make<TTree>("t","t");
    tree->Branch("pfJetCounter",             &pfJetCounter_,              "pfJetCounter/I")                                     ;
    tree->Branch("pfJetPt",                  &pfJetPt_,                   "pfJetPT[pfJetCounter]/D")                            ;
    tree->Branch("pfJetEta",                 &pfJetEta_,                  "pfJetEta[pfJetCounter]/F")                           ;
    tree->Branch("pfJetPhi",                 &pfJetPhi_,                  "pfJetPhi[pfJetCounter]/F")                           ;
    tree->Branch("pfJetEt",                  &pfJetEt_,                   "pfJetEt[pfJetCounter]/D")                            ;
    tree->Branch("pfJetEnergy",              &pfJetEnergy_,               "pfJetEnergy[pfJetCounter]/D")                        ;
    tree->Branch("pfJetPx",                  &pfJetPx_,                   "pfJetPx[pfJetCounter]/D")                            ;
    tree->Branch("pfJetPy",                  &pfJetPy_,                   "pfJetPy[pfJetCounter]/D")                            ;
    tree->Branch("pfJetPz",                  &pfJetPz_,                   "pfJetPz[pfJetCounter]/D")                            ;
    tree->Branch("pfJetCSV",                 &pfJetCSV_,                     "pfJetCSV[pfJetCounter]/D")                            ;
    createHistograms();
}

void MiniAODanalyzer::createHistograms(){
    TFileDirectory subDir = fs_->mkdir("plots");
    std::map<std::string, TH1*> container;

    container["bbPT_C160_p026_online"] =  subDir.make<TH1F>("bbPT_C160_p026_online", "bbPT_C160_p026_online", 200, 0, 2000);
    container["bbPT_C100_p014_online"] =  subDir.make<TH1F>("bbPT_C100_p014_online", "bbPT_C100_p014_online", 200, 0, 2000);
    container["bbPT_C172_p026_online"] =  subDir.make<TH1F>("bbPT_C172_p026_online", "bbPT_C172_p026_online", 200, 0, 2000);
    container["bbPT_C112_p014_online"] =  subDir.make<TH1F>("bbPT_C112_p014_online", "bbPT_C112_p014_online", 200, 0, 2000);
    container["bbPT_C160_p026_offline"] = subDir.make<TH1F>("bbPT_C160_p026_offline", "bbPT_C160_p026_offline", 200, 0, 2000);
    container["bbPT_C100_p014_offline"] = subDir.make<TH1F>("bbPT_C100_p014_offline", "bbPT_C100_p014_offline", 200, 0, 2000);
    container["bbPT_C172_p026_offline"] = subDir.make<TH1F>("bbPT_C172_p026_offline", "bbPT_C172_p026_offline", 200, 0, 2000);
    container["bbPT_C112_p014_offline"] = subDir.make<TH1F>("bbPT_C112_p014_offline", "bbPT_C112_p014_offline", 200, 0, 2000);
    container["bosonPT"]                = subDir.make<TH1F>("bosonPT", "bosonPT", 200, 0, 2000);
    container["genPT"]                  = subDir.make<TH1F>("genPT"  , "genPT"  , 200, 0, 2000);
    container["OfflineBtag"]            = subDir.make<TH1F>("OfflineBtag", "OfflineBtag", 4, 0, 4);
    container["OnlineBtag"]             = subDir.make<TH1F>("OnlineBtag" , "OnlineBtag" , 4, 0, 4);
    
    container["Njets"]                  = subDir.make<TH1F>("Njets" , "Njets" , 10, 0, 10);

    histos_.push_back(container);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAODanalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAODanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODanalyzer);
