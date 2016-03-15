// -*- C++ -*-
//
// Package:    QCDAnalyzer
// Class:      QCDAnalyzer
// 
/**\class QCDAnalyzer QCDAnalyzer.cc Appeltel/QCDAnalyzer/src/QCDAnalyzer.cc

   Description: Gen-Level analyzer for QCD processes (i.e. jets, charged particles, parton flavor)

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Eric Appelt
//         Created:  Tue Feb 25 14:25:23 CST 2014
// Modified by Raghav Kunnawalkam Elayavalli, Monday Mar 14th 2016
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "HepPDT/ParticleID.hh"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"


#include <TH1.h>
#include <TH2.h>

//
// class declaration
//

class QCDAnalyzer : public edm::EDAnalyzer {
public:
  explicit QCDAnalyzer(const edm::ParameterSet&);
  ~QCDAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  bool isDSEvent( const edm::Event&, const edm::EventSetup& );
  bool isInEtaRange( const reco::Candidate&, double, double ); 
  bool isInFlavor( const reco::MatchedPartons& ); 
  // bool isInFlavor_LeadJet( const reco::MatchedPartons& ); 
  // bool isInFlavor_SubLeadJet( const reco::MatchedPartons& );
  int GetFlavor( const reco::MatchedPartons& ); 
  bool isInSpecies( const reco::GenParticle& ); 

  // ----------member data ---------------------------

  edm::InputTag genJetSrc_;
  edm::InputTag genParticleSrc_;
  bool doFlavor_;
  bool doSpecies_;
  bool doCombBKGStudy_;
  double combBKGRMS_;
  bool onlyDS_;
  edm::InputTag flavorSrc_;
  std::vector<int> flavorId_;
  // std::vector<int> flavorId_LeadJet_;
  // std::vector<int> flavorId_SubLeadJet_;
  std::vector<int> speciesId_;
  bool useRapidity_;
  double jetEtaMin_;
  double jetEtaMax_;
  double hEtaMin_;
  double hEtaMax_;
  double jetRadius_;
  double pthatMin_;
  double pthatMax_;
  std::vector<double> jetPtBins_;
  std::vector<double> jetPhiBins_;
  std::vector<double> hPtBins_;
  std::vector<double> qScalePtBins_;
  std::vector<double> etaBins_;

  // for combinatorial BKG study
  double delPhi;
  double LeadPhi;
  double SubLeadPhi;
  double LeadJetpT;
  double SubLeadJetpT;
  // bool isFlavorLeadJet;
  // bool isFlavorSubLeadJet;

  std::vector <double> pt;
  std::vector <double> phi;
  std::vector <double> eta;
  std::vector <int> flav;

  std::vector <double> par_pt;
  std::vector <double> par_eta;
  std::vector <double> par_phi;
  std::vector <int> par_flav;

  std::map<std::string,TH1F*> hist_;
  std::map<std::string,TH2F*> hist2D_;
};


//
// constructors and destructor
//

QCDAnalyzer::QCDAnalyzer(const edm::ParameterSet& iConfig):
  genJetSrc_(iConfig.getParameter<edm::InputTag>("genJetSrc")),
  genParticleSrc_(iConfig.getParameter<edm::InputTag>("genParticleSrc")),
  doFlavor_(iConfig.getParameter<bool>("doFlavor")),
  doSpecies_(iConfig.getParameter<bool>("doSpecies")),
  doCombBKGStudy_(iConfig.getParameter<bool>("doCombBKGStudy")),
  combBKGRMS_(iConfig.getParameter<double>("combBKG_RMS")),
  onlyDS_(iConfig.getParameter<bool>("onlyDS")),
  flavorSrc_(iConfig.getParameter<edm::InputTag>("flavorSrc")),
  flavorId_(iConfig.getParameter<std::vector<int> >("flavorId")),
  // flavorId_LeadJet_(iConfig.getParameter<std::vector<int> >("flavorId_LeadJet")),
  // flavorId_SubLeadJet_(iConfig.getParameter<std::vector<int> >("flavorId_SubLeadJet")),
  speciesId_(iConfig.getParameter<std::vector<int> >("speciesId")),
  useRapidity_(iConfig.getParameter<bool>("useRapidity")),
  jetEtaMin_(iConfig.getParameter<double>("jetEtaMin")),
  jetEtaMax_(iConfig.getParameter<double>("jetEtaMax")),
  hEtaMin_(iConfig.getParameter<double>("hEtaMin")),
  hEtaMax_(iConfig.getParameter<double>("hEtaMax")),
  jetRadius_(iConfig.getParameter<double>("jetRadius")),
  pthatMin_(iConfig.getParameter<double>("pthatMin")),
  pthatMax_(iConfig.getParameter<double>("pthatMax")),
  jetPtBins_(iConfig.getParameter<std::vector<double> >("jetPtBins")),
  jetPhiBins_(iConfig.getParameter<std::vector<double> >("jetPhiBins")),
  hPtBins_(iConfig.getParameter<std::vector<double> >("hPtBins")),
  qScalePtBins_(iConfig.getParameter<std::vector<double> >("qScalePtBins")),
  etaBins_(iConfig.getParameter<std::vector<double> >("etaBins"))
{
  edm::Service<TFileService> fs;

  // Spectra histograms for eta (or y) within specified range

  // jet spectrum
  hist_["jetspectrum"] = fs->make<TH1F>("jetspectrum",";p_{T};counts",
					jetPtBins_.size()-1, &jetPtBins_[0]);
  // hadron spectrum
  hist_["hspectrum"] = fs->make<TH1F>("hspectrum",";p_{T};counts",
				      hPtBins_.size()-1, &hPtBins_[0]);
  // charged hadron spectrum
  hist_["chspectrum"] = fs->make<TH1F>("chspectrum",";p_{T};counts",
				       hPtBins_.size()-1, &hPtBins_[0]);

  // Combinatorial BKG study
  hist_["hAllJets"] = fs->make<TH1F>("hAllJets",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSignal"] = fs->make<TH1F>("hSignal",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hBKG"] = fs->make<TH1F>("hBKG",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hDeltaPhi"] = fs->make<TH1F>("hDeltaPhi",";p_{T};counts",
				     jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist_["hLeading"] = fs->make<TH1F>("hLeading",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSubLeading"] = fs->make<TH1F>("hSubLeading",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  // make 2D histograms for leading vs subleading and delta phi vs leading and subleading 
  hist2D_["hLeading_vs_SubLeading"] = fs->make<TH2F>("hLeading_vs_SubLeading",";SubLeading jet p_{T}; Leading Jet p_{T}",
						     jetPtBins_.size()-1, &jetPtBins_[0],
						     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist2D_["hDeltaPhi_vs_Leading"] = fs->make<TH2F>("hDeltaPhi_vs_Leading",";Leading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						   jetPtBins_.size()-1, &jetPtBins_[0],
						   jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist2D_["hDeltaPhi_vs_SubLeading"] = fs->make<TH2F>("hDeltaPhi_vs_SubLeading",";SubLeading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						      jetPtBins_.size()-1, &jetPtBins_[0],
						      jetPhiBins_.size()-1, &jetPhiBins_[0]);

  // qq
  hist_["hAllJets_qq"] = fs->make<TH1F>("hAllJets_qq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSignal_qq"] = fs->make<TH1F>("hSignal_qq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hBKG_qq"] = fs->make<TH1F>("hBKG_qq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hDeltaPhi_qq"] = fs->make<TH1F>("hDeltaPhi_qq",";#Delta Phi;counts",
				     jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist_["hLeading_qq"] = fs->make<TH1F>("hLeading_qq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSubLeading_qq"] = fs->make<TH1F>("hSubLeading_qq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);  
  // make 2D histograms for leading vs subleading and delta phi vs leading and subleading 
  hist2D_["hLeading_vs_SubLeading_qq"] = fs->make<TH2F>("hLeading_vs_SubLeading_qq",";SubLeading jet p_{T}; Leading Jet p_{T}",
						     jetPtBins_.size()-1, &jetPtBins_[0],
						     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist2D_["hDeltaPhi_vs_Leading_qq"] = fs->make<TH2F>("hDeltaPhi_vs_Leading_qq",";Leading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						   jetPtBins_.size()-1, &jetPtBins_[0],
						   jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist2D_["hDeltaPhi_vs_SubLeading_qq"] = fs->make<TH2F>("hDeltaPhi_vs_SubLeading_qq",";SubLeading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						      jetPtBins_.size()-1, &jetPtBins_[0],
						      jetPhiBins_.size()-1, &jetPhiBins_[0]);

  // gq
  hist_["hAllJets_gq"] = fs->make<TH1F>("hAllJets_gq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSignal_gq"] = fs->make<TH1F>("hSignal_gq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hBKG_gq"] = fs->make<TH1F>("hBKG_gq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hDeltaPhi_gq"] = fs->make<TH1F>("hDeltaPhi_gq",";#Delta Phi;counts",
				     jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist_["hLeading_gq"] = fs->make<TH1F>("hLeading_gq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSubLeading_gq"] = fs->make<TH1F>("hSubLeading_gq",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);  
  // make 2D histograms for leading vs subleading and delta phi vs leading and subleading 
  hist2D_["hLeading_vs_SubLeading_gq"] = fs->make<TH2F>("hLeading_vs_SubLeading_gq",";SubLeading jet p_{T}; Leading Jet p_{T}",
						     jetPtBins_.size()-1, &jetPtBins_[0],
						     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist2D_["hDeltaPhi_vs_Leading_gq"] = fs->make<TH2F>("hDeltaPhi_vs_Leading_gq",";Leading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						   jetPtBins_.size()-1, &jetPtBins_[0],
						   jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist2D_["hDeltaPhi_vs_SubLeading_gq"] = fs->make<TH2F>("hDeltaPhi_vs_SubLeading_gq",";SubLeading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						      jetPtBins_.size()-1, &jetPtBins_[0],
						      jetPhiBins_.size()-1, &jetPhiBins_[0]);

  // qg
  hist_["hAllJets_qg"] = fs->make<TH1F>("hAllJets_qg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSignal_qg"] = fs->make<TH1F>("hSignal_qg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hBKG_qg"] = fs->make<TH1F>("hBKG_qg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hDeltaPhi_qg"] = fs->make<TH1F>("hDeltaPhi_qg",";#Delta Phi;counts",
				     jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist_["hLeading_qg"] = fs->make<TH1F>("hLeading_qg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSubLeading_qg"] = fs->make<TH1F>("hSubLeading_qg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);  
  // make 2D histograms for leading vs subleading and delta phi vs leading and subleading 
  hist2D_["hLeading_vs_SubLeading_qg"] = fs->make<TH2F>("hLeading_vs_SubLeading_qg",";SubLeading jet p_{T}; Leading Jet p_{T}",
						     jetPtBins_.size()-1, &jetPtBins_[0],
						     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist2D_["hDeltaPhi_vs_Leading_qg"] = fs->make<TH2F>("hDeltaPhi_vs_Leading_qg",";Leading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						   jetPtBins_.size()-1, &jetPtBins_[0],
						   jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist2D_["hDeltaPhi_vs_SubLeading_qg"] = fs->make<TH2F>("hDeltaPhi_vs_SubLeading_qg",";SubLeading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						      jetPtBins_.size()-1, &jetPtBins_[0],
						      jetPhiBins_.size()-1, &jetPhiBins_[0]);

  //gg
  hist_["hAllJets_gg"] = fs->make<TH1F>("hAllJets_gg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSignal_gg"] = fs->make<TH1F>("hSignal_gg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hBKG_gg"] = fs->make<TH1F>("hBKG_gg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hDeltaPhi_gg"] = fs->make<TH1F>("hDeltaPhi_gg",";#Delta Phi;counts",
				     jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist_["hLeading_gg"] = fs->make<TH1F>("hLeading_gg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist_["hSubLeading_gg"] = fs->make<TH1F>("hSubLeading_gg",";p_{T};counts",
				     jetPtBins_.size()-1, &jetPtBins_[0]);  
  // make 2D histograms for leading vs subleading and delta phi vs leading and subleading 
  hist2D_["hLeading_vs_SubLeading_gg"] = fs->make<TH2F>("hLeading_vs_SubLeading_gg",";SubLeading jet p_{T}; Leading Jet p_{T}",
						     jetPtBins_.size()-1, &jetPtBins_[0],
						     jetPtBins_.size()-1, &jetPtBins_[0]);
  hist2D_["hDeltaPhi_vs_Leading_gg"] = fs->make<TH2F>("hDeltaPhi_vs_Leading_gg",";Leading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						   jetPtBins_.size()-1, &jetPtBins_[0],
						   jetPhiBins_.size()-1, &jetPhiBins_[0]);
  hist2D_["hDeltaPhi_vs_SubLeading_gg"] = fs->make<TH2F>("hDeltaPhi_vs_SubLeading_gg",";SubLeading jet p_{T}; #Delta #phi_{Lead-SubLead}",
						      jetPtBins_.size()-1, &jetPtBins_[0],
						      jetPhiBins_.size()-1, &jetPhiBins_[0]);

    
  // positively charged hadron spectrum
  hist_["pchspectrum"] = fs->make<TH1F>("pchspectrum",";p_{T};counts",
					hPtBins_.size()-1, &hPtBins_[0]);


  // 2D histograms of (eta[y],pT) as above
  hist2D_["jetspectrum2D"] = fs->make<TH2F>("jetspectrum2D",";p_{T};counts",
					    etaBins_.size()-1, &etaBins_[0],
					    jetPtBins_.size()-1, &jetPtBins_[0]);
  hist2D_["hspectrum2D"] = fs->make<TH2F>("hspectrum2D",";p_{T};counts",
					  etaBins_.size()-1, &etaBins_[0],
					  hPtBins_.size()-1, &hPtBins_[0]);
  hist2D_["chspectrum2D"] = fs->make<TH2F>("chspectrum2D",";p_{T};counts",
					   etaBins_.size()-1, &etaBins_[0],
					   hPtBins_.size()-1, &hPtBins_[0]);
  hist2D_["pchspectrum2D"] = fs->make<TH2F>("pchspectrum2D",";p_{T};counts",
					    etaBins_.size()-1, &etaBins_[0],
					    hPtBins_.size()-1, &hPtBins_[0]);

  // pt-hat or momentum transfer scale of PYTHIA process
  // this will be 0 for diffractive events
  hist_["qscale"] = fs->make<TH1F>("qscale",";p_{T}-hat;counts",
				   qScalePtBins_.size()-1, &qScalePtBins_[0]);
  hist2D_["ch_qscale2D"] = fs->make<TH2F>("ch_qscale2D",";p_{T}-hat;p_{T} h^{#pm}",
					  qScalePtBins_.size()-1, &qScalePtBins_[0],
					  hPtBins_.size()-1, &hPtBins_[0]);
  hist2D_["jet_qscale2D"] = fs->make<TH2F>("jet_qscale2D",";p_{T}-hat;p_{T} h^{#pm}",
					   qScalePtBins_.size()-1, &qScalePtBins_[0],
					   jetPtBins_.size()-1, &jetPtBins_[0]);

  // leading charged hadron in approximate tracker acceptance (|eta|<2.5)
  hist_["lead_track"] = fs->make<TH1F>("lead_track",";p_{T};counts",
				       hPtBins_.size()-1, &hPtBins_[0]);

  // number of recorded events
  hist_["events"] = fs->make<TH1F>("events",";;events",1,0.,2.);

  // fragmentation function matrix - i.e. 2D histogram of hadrons vs jets
  // also for just charged and positively charged hadrons
  hist2D_["ffmatrix"] = fs->make<TH2F>("ffmatrix",";p_{T} Jet;p_{T} h",
				       hPtBins_.size()-1, &hPtBins_[0],
				       jetPtBins_.size()-1, &jetPtBins_[0]);
  hist2D_["cffmatrix"] = fs->make<TH2F>("cffmatrix",";p_{T} Jet;p_{T} h^{#pm}",
					hPtBins_.size()-1, &hPtBins_[0],
					jetPtBins_.size()-1, &jetPtBins_[0]);
  hist2D_["pcffmatrix"] = fs->make<TH2F>("pcffmatrix",";p_{T} Jet;p_{T} h^{+}",
					 hPtBins_.size()-1, &hPtBins_[0],
					 jetPtBins_.size()-1, &jetPtBins_[0]);

}


QCDAnalyzer::~QCDAnalyzer()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
QCDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // Pythia does not respect max pt-hat for MB process  
  // where the minimum is 0. In this case we need to 
  // remove events over the pt-hat maximum by hand.
  // We skip and do not count such events.
  // Here it is only coded for the MB_0_to_20 process
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByLabel("generator", genEvtInfo);
   
  if( pthatMin_ < 1. )
    {
      if( genEvtInfo->qScale() > pthatMax_ ) return;
    }

  // Check if the event is DS and skip if configured
  if( ! isDSEvent(iEvent,iSetup) && onlyDS_ ) return;
  
  Handle<reco::GenParticleCollection> pcol;
  iEvent.getByLabel(genParticleSrc_,pcol);

  Handle<reco::GenJetCollection> gcol;
  if( !doFlavor_ ) iEvent.getByLabel(genJetSrc_,gcol);

  Handle<reco::JetMatchedPartonsCollection> fcol;
  if( doFlavor_) iEvent.getByLabel(flavorSrc_,fcol);

  hist_["events"]->Fill(1);
  hist_["qscale"]->Fill(genEvtInfo->qScale());

  // std::cout<<"Starting new Event "<<std::endl;
  int counter = 0;
  pt.clear();
  phi.clear();
  flav.clear();
  eta.clear();

  par_pt.clear();
  par_eta.clear();
  par_phi.clear();
  par_flav.clear();
  
  // std::cout<<"Leading Jet pT = "<<pt[0] <<"; eta = "<<eta[0]<<"; phi = "<<phi[0]<<std::endl;
  // std::cout<<"SubLeading Jet pT = "<<pt[1] <<"; eta = "<<eta[1]<<"; phi = "<<phi[1]<<std::endl;
  // std::cout<<"third  Jet pT = "<<pt[2] <<"; eta = "<<eta[2]<<"; phi = "<<phi[2]<<std::endl;
  
  // charged particle spectra and ff matrix
  double lead_track_pt = 0.0;
  for( const auto & h : *pcol )
    {
      counter++;
      // Only look at the initial scalttered particles 
      // if( h.status() == 2 || h.status() == 3  ){
      // 	if(fabs(h.pdgId())<=2 || h.pdgId()==21) std::cout<<"gen particle type = "<<h.pdgId()<<std::endl;
      // }
      // number 7 and 8 are the initial scattered particles. 
      if(counter == 7 || counter == 8)
	{
	  // std::cout<<"Gen Particle: "<<counter<<"; id = "<<h.pdgId()<<"; pt = "<<h.pt()<<"; eta = "<<h.eta()<<"; phi = "<<h.phi()<<std::endl<<"                DeltaR with Leading = "<<(float)deltaR(h.eta(), h.phi(), eta[0], phi[0])<<std::endl<<"                DeltaR with subleading = "<<(float)deltaR(h.eta(), h.phi(), eta[1], phi[1])<<std::endl<<"                DeltaR with third = "<<(float)deltaR(h.eta(), h.phi(), eta[2], phi[2])<<std::endl;
	  par_pt.push_back(h.pt());
	  par_eta.push_back(h.eta());
	  par_phi.push_back(h.phi());
	  par_flav.push_back(h.pdgId());
	}
      
    }
  
  
  counter = 0;

  // genjet spectra
  if ( doFlavor_ )
    {
      // std::cout<<" Passed doFlavor selection "<<std::endl;
      for( const auto & mjp : *fcol )
	{
	  counter++;	  
	  // std::cout<<" Counter of Jets = "<<counter<<std::endl;
	  // const reco::MatchedPartons & aMatch = mjp.second;
	  
	  // // std::cout<<"Flavor of matched parton (second) heaviest = "<<aMatch.heaviest().get()->pdgId()<<std::endl;
	  
	  // // std::cout<<"Flavor of matched parton (second) nearest_status2  = "<<aMatch.nearest_status2().get()->pdgId()<<std::endl;
	  // // std::cout<<"Flavor of matched parton (second) nearest_status3  = "<<aMatch.nearest_status3().get()->pdgId()<<std::endl;
	  // // std::cout<<"Flavor of matched parton (second) physicsDefinitionParton  = "<<aMatch.physicsDefinitionParton().get()->pdgId()<<std::endl;
	  // // std::cout<<"Flavor of matched parton (second) algoDefinitionParton  = "<<aMatch.algoDefinitionParton().get()->pdgId()<<std::endl;
	  
	  // if( ! isInFlavor( aMatch ) ) continue;
	  // // std::cout<<" Passed matching to quarks or gluons "<<std::endl;
	  const reco::Jet *aJet = mjp.first.get();
	  if( isInEtaRange( *aJet, jetEtaMin_, jetEtaMax_ ) ) 
	    {
	      // // std::cout<<" Passed eta selection  "<<std::endl;
	      if ( doCombBKGStudy_ )
		{
		  // std::cout<<" Loading information for Performing combinatorial study  "<<std::endl;
		  // if(counter == 0)
		  //   {
		  //     // std::cout<<" Lead Jet information collected  "<<std::endl;
		  //     isFlavorLeadJet = (bool) isInFlavor_LeadJet( aMatch );
		  //     LeadPhi = aJet->phi();
		  //     LeadJetpT = aJet->pt();
		  //   }
		  // if(counter == 1)
		  //   {
		  //     // std::cout<<" Sub-Lead Jet information collected  "<<std::endl;
		  //     isFlavorSubLeadJet = (bool) isInFlavor_SubLeadJet( aMatch );
		  //     SubLeadPhi = aJet->phi();		      
		  //     SubLeadJetpT = aJet->pt();
		  //   }

		  int flavor = (int)GetFlavor( mjp.second );
		  pt.push_back(aJet->pt());
		  phi.push_back(aJet->phi());
		  eta.push_back(aJet->eta());
		  flav.push_back(flavor);
		}
	      
	      hist_["jetspectrum"]->Fill(aJet->pt());
	      // std::cout<<"Jet pT = "<<aJet->pt()<<std::endl;		
	      hist2D_["jet_qscale2D"]->Fill(genEvtInfo->qScale(),aJet->pt());
	    }
	  hist2D_["jetspectrum2D"]->Fill(aJet->eta(),aJet->pt());
	}
    }
  else
    {
      for( const auto & jet : *gcol )
	{
	  if( isInEtaRange( jet, jetEtaMin_, jetEtaMax_ ) ) 
	    {
	      hist_["jetspectrum"]->Fill(jet.pt());
	      hist2D_["jet_qscale2D"]->Fill(genEvtInfo->qScale(),jet.pt());
	    }
	  hist2D_["jetspectrum2D"]->Fill(jet.eta(),jet.pt());
	}
    }

  // std::cout<<"Total number of jets in event = "<<pt.size()<<std::endl;
  if(pt.size()<=1) return;
  if(pt[0]<10 || pt[1]<10) return;
  // std::cout<<" Now checking if event is dijet parton matched  "<<std::endl;
  // check which flavor event is in and fill histograms for qq, qg, gq and gg

  // Now to find the flavor of the dijet event. need to match the partons from the initial scattering to the leading/subleading jets.

  bool isqq = false;
  bool isqg = false;
  bool isgq = false;
  bool isgg = false;

  // std::vector <std::vector <double> > delR_par_jet;

  // for(unsigned ij = 0; ij<pt.size(); ++ij)
  //   {
  //     for(unsigned ip = 0; ip<par_pt.size(); ++ip)
  // 	{
  // 	  delR_par_jet[ij][ip] = deltaR2(par_eta[ip], par_phi[ip], eta[ij], phi[ij]);
  // 	}
  //   }
    
  double delR_par1_jet1 = deltaR2(par_eta[0], par_phi[0], eta[0], phi[0]);
  double delR_par1_jet2 = deltaR2(par_eta[0], par_phi[0], eta[1], phi[1]);
  double delR_par2_jet1 = deltaR2(par_eta[1], par_phi[1], eta[0], phi[0]);
  double delR_par2_jet2 = deltaR2(par_eta[1], par_phi[1], eta[1], phi[1]);

  if((delR_par1_jet1 > 0.3 && delR_par1_jet2 > 0.3) || (delR_par2_jet1 > 0.3 && delR_par2_jet2 > 0.3)) return;

  if(fabs(par_flav[0])<=6 && fabs(par_flav[1])<=6 )
    isqq = true;
  if(par_flav[0] == 21 && par_flav[1] == 21 )
    isgg = true;
  
  if(fabs(par_flav[0])<=6 && par_flav[1] == 21) 
    {
      if((delR_par1_jet1 < delR_par1_jet2)) isqg = true;
      else isgq = true;
    }
  if(par_flav[0]==21 && fabs(par_flav[1])<=6 )
    {
      if((delR_par1_jet1 < delR_par1_jet2)) isgq = true;
      else isqg = true;
    }
      
  bool isSignal = true;
  bool isSignal_qq = true;
  bool isSignal_qg = true;
  bool isSignal_gq = true;
  bool isSignal_gg = true;

  // std::cout<<"Dijet Flavor  = "<<par_flav[0]<<"; "<<par_flav[1]<<std::endl;
  // std::cout<<"Delta R, par1 = ("<<delR_par1_jet1<<","<<delR_par1_jet2<<")"<<std::endl;
  // std::cout<<"Delta R, par2 = ("<<delR_par2_jet1<<","<<delR_par2_jet2<<")"<<std::endl;
  // std::cout<<"Boolean values = "<<std::endl;
  // std::cout<<"    isqq = "<<isqq<<std::endl;
  // std::cout<<"    isgq = "<<isgq<<std::endl;
  // std::cout<<"    isqg = "<<isqg<<std::endl;
  // std::cout<<"    isgg = "<<isgg<<std::endl;

  delPhi = deltaPhi(phi[0], phi[1]);

  if( pt[0] < 3 * combBKGRMS_ )
    isSignal = false;
  if( pt[0] > 3 * combBKGRMS_ && delPhi < 2.0943 )
    isSignal = false;
  if( pt[0] > 3 * combBKGRMS_ && delPhi > 2.0943 && pt[1] < 1.6 * combBKGRMS_ )
    isSignal = false;

  // std::cout<<"Checking for signal"<<std::endl;
  
  // qq
  if(isqq)
    {
      if( pt[0] < 3 * combBKGRMS_ )
	isSignal_qq = false;
      if( pt[0] > 3 * combBKGRMS_ && delPhi < 2.0943 )
	isSignal_qq = false;
      if( pt[0] > 3 * combBKGRMS_ && delPhi > 2.0943 && pt[1] < 1.6 * combBKGRMS_ )
	isSignal_qq = false;      
    }
  // qg
  if(isqg)
    {
      if( pt[0] < 3 * combBKGRMS_ )
	isSignal_qg = false;
      if( pt[0] > 3 * combBKGRMS_ && delPhi < 2.0943 )
	isSignal_qg = false;
      if( pt[0] > 3 * combBKGRMS_ && delPhi > 2.0943 && pt[1] < 1.6 * combBKGRMS_ )
	isSignal_qg = false;      
    }
  // gq
  if(isgq)
    {
      if( pt[0] < 3 * combBKGRMS_ )
	isSignal_gq = false;
      if( pt[0] > 3 * combBKGRMS_ && delPhi < 2.0943 )
	isSignal_gq = false;
      if( pt[0] > 3 * combBKGRMS_ && delPhi > 2.0943 && pt[1] < 1.6 * combBKGRMS_ )
	isSignal_gq = false;      
    }
  // gg
  if(isgg)
    {
      if( pt[0] < 3 * combBKGRMS_ )
	isSignal_gg = false;
      if( pt[0] > 3 * combBKGRMS_ && delPhi < 2.0943 )
	isSignal_gg = false;
      if( pt[0] > 3 * combBKGRMS_ && delPhi > 2.0943 && pt[1] < 1.6 * combBKGRMS_ )
	isSignal_gg = false;      
    }


  // std::cout<<"Filling event histograms"<<std::endl;

  hist_["hDeltaPhi"]->Fill(delPhi);
  hist_["hLeading"]->Fill(pt[0]);
  hist_["hSubLeading"]->Fill(pt[1]);
  hist2D_["hLeading_vs_SubLeading"]->Fill(pt[1], pt[0]);
  hist2D_["hDeltaPhi_vs_Leading"]->Fill(pt[0], delPhi);  
  hist2D_["hDeltaPhi_vs_SubLeading"]->Fill(pt[1], delPhi);

  // std::cout<<"Filling event histograms qq"<<std::endl;

  if(isqq)
    {
      hist_["hDeltaPhi_qq"]->Fill(delPhi);
      hist_["hLeading_qq"]->Fill(pt[0]);
      hist_["hSubLeading_qq"]->Fill(pt[1]);
      hist2D_["hLeading_vs_SubLeading_qq"]->Fill(pt[1], pt[0]);
      hist2D_["hDeltaPhi_vs_Leading_qq"]->Fill(pt[0], delPhi);  
      hist2D_["hDeltaPhi_vs_SubLeading_qq"]->Fill(pt[1], delPhi);
    }
  // std::cout<<"Filling event histograms gg"<<std::endl;

  if(isgg)
    {
      hist_["hDeltaPhi_gg"]->Fill(delPhi);
      hist_["hLeading_gg"]->Fill(pt[0]);
      hist_["hSubLeading_gg"]->Fill(pt[1]);
      hist2D_["hLeading_vs_SubLeading_gg"]->Fill(pt[1], pt[0]);
      hist2D_["hDeltaPhi_vs_Leading_gg"]->Fill(pt[0], delPhi);  
      hist2D_["hDeltaPhi_vs_SubLeading_gg"]->Fill(pt[1], delPhi);
    }
  // std::cout<<"Filling event histograms qg"<<std::endl;

  if(isqg)
    {
      hist_["hDeltaPhi_qg"]->Fill(delPhi);
      hist_["hLeading_qg"]->Fill(pt[0]);
      hist_["hSubLeading_qg"]->Fill(pt[1]);
      hist2D_["hLeading_vs_SubLeading_qg"]->Fill(pt[1], pt[0]);
      hist2D_["hDeltaPhi_vs_Leading_qg"]->Fill(pt[0], delPhi);  
      hist2D_["hDeltaPhi_vs_SubLeading_qg"]->Fill(pt[1], delPhi);
    }
  // std::cout<<"Filling event histograms qg"<<std::endl;

  if(isgq)
    {
      hist_["hDeltaPhi_gq"]->Fill(delPhi);
      hist_["hLeading_gq"]->Fill(pt[0]);
      hist_["hSubLeading_gq"]->Fill(pt[1]);
      hist2D_["hLeading_vs_SubLeading_gq"]->Fill(pt[1], pt[0]);
      hist2D_["hDeltaPhi_vs_Leading_gq"]->Fill(pt[0], delPhi);  
      hist2D_["hDeltaPhi_vs_SubLeading_gq"]->Fill(pt[1], delPhi);
    }
  
  // std::cout<<"Filling Jet histograms"<<std::endl;
  
  for(unsigned ij = 0; ij<pt.size(); ++ij)
    {
      // std::cout<<"counter = "<<ij<<std::endl;
      hist_["hAllJets"]->Fill(pt[ij]);
      if(isSignal) hist_["hSignal"]->Fill(pt[ij]);
      if(!isSignal) hist_["hBKG"]->Fill(pt[ij]);

      // std::cout<<"Filling Jet histograms qq"<<std::endl;
      if(isqq)
	{
	  hist_["hAllJets_qq"]->Fill(pt[ij]);
	  if(isSignal_qq) hist_["hSignal_qq"]->Fill(pt[ij]);
	  if(!isSignal_qq) hist_["hBKG_qq"]->Fill(pt[ij]);
	}
      // std::cout<<"Filling Jet histograms gg"<<std::endl;
      if(isgg)
	{
	  hist_["hAllJets_gg"]->Fill(pt[ij]);
	  if(isSignal_gg) hist_["hSignal_gg"]->Fill(pt[ij]);
	  if(!isSignal_gg) hist_["hBKG_gg"]->Fill(pt[ij]);
	}
      // std::cout<<"Filling Jet histograms qg"<<std::endl;
      if(isgq)
	{
	  hist_["hAllJets_qg"]->Fill(pt[ij]);
	  if(isSignal_qg) hist_["hSignal_qg"]->Fill(pt[ij]);
	  if(!isSignal_qg) hist_["hBKG_qg"]->Fill(pt[ij]);
	}
      // std::cout<<"Filling Jet histograms qq"<<std::endl;
      if(isgq)
	{
	  hist_["hAllJets_gq"]->Fill(pt[ij]);
	  if(isSignal_gq) hist_["hSignal_gq"]->Fill(pt[ij]);
	  if(!isSignal_gq) hist_["hBKG_gq"]->Fill(pt[ij]);
	}
    }
  
  
  
  // // do the combinatorial BKG checks and histogram fills
  // if ( isFlavorLeadJet && isFlavorSubLeadJet )
  //   {
  //     // std::cout<<" Event is dijet parton matched  "<<std::endl;

  //     delPhi = deltaPhi(LeadPhi, SubLeadPhi);

  //     hist_["hDeltaPhi"]->Fill(delPhi);
  //     hist_["hLeading"]->Fill(LeadJetpT);
  //     hist_["hSubLeading"]->Fill(SubLeadJetpT);
  //     hist2D_["hLeading_vs_SubLeading"]->Fill(SubLeadJetpT, LeadJetpT);
  //     hist2D_["hDeltaPhi_vs_Leading"]->Fill(LeadJetpT, delPhi);
  //     hist2D_["hDeltaPhi_vs_SubLeading"]->Fill(SubLeadJetpT, delPhi);
      
  //     bool isSignal = true;
  //     if( LeadJetpT < 3 * combBKGRMS_ )
  // 	isSignal = false;
  //     if( LeadJetpT > 3 * combBKGRMS_ && delPhi < 2.0943 )
  // 	isSignal = false;
  //     if( LeadJetpT > 3 * combBKGRMS_ && delPhi > 2.0943 && SubLeadJetpT < 1.6 * combBKGRMS_ )
  // 	isSignal = false;

  //     if( isSignal )
  // 	{
  // 	  // std::cout<<"Event is Signal"<<std::endl;
  // 	}
  //     else
  // 	{
  // 	  // std::cout<<"Event is BKG"<<std::endl;
  // 	}
  //     for( const auto & mjp : *fcol )
  // 	{
  // 	  if( ! isInFlavor( mjp.second ) ) continue;
  // 	  const reco::Jet *aJet = mjp.first.get();
  // 	  if( isInEtaRange( *aJet, jetEtaMin_, jetEtaMax_ ) ) 
  // 	    {
  // 	      hist_["AllJets"]->Fill(aJet->pt());
  // 	      if ( isSignal )
  // 		hist_["hSignal"]->Fill(aJet->pt());
  // 	      if ( !isSignal )
  // 		hist_["hBKG"]->Fill(aJet->pt());
  // 	    }
  // 	}
  
  //   }

  
  // std::cout<<"Leading Jet pT = "<<pt[0] <<"; eta = "<<eta[0]<<"; phi = "<<phi[0]<<std::endl;
  // std::cout<<"SubLeading Jet pT = "<<pt[1] <<"; eta = "<<eta[1]<<"; phi = "<<phi[1]<<std::endl;
  // std::cout<<"third  Jet pT = "<<pt[2] <<"; eta = "<<eta[2]<<"; phi = "<<phi[2]<<std::endl;
  
  // charged particle spectra and ff matrix
  lead_track_pt = 0.0;
  for( const auto & h : *pcol )
    {

      if(h.status() != 1) continue;
      
      // counter++;
      // // Only look at the initial scalttered particles 
      // // if( h.status() == 2 || h.status() == 3  ){
      // // 	if(fabs(h.pdgId())<=2 || h.pdgId()==21) std::cout<<"gen particle type = "<<h.pdgId()<<std::endl;
      // // }
      // // number 7 and 8 are the initial scattered particles. 
      // if(counter == 7 || counter == 8)
      // 	{
      // 	  // std::cout<<"Gen Particle: "<<counter<<"; id = "<<h.pdgId()<<"; pt = "<<h.pt()<<"; eta = "<<h.eta()<<"; phi = "<<h.phi()<<std::endl<<"                DeltaR with Leading = "<<(float)deltaR(h.eta(), h.phi(), eta[0], phi[0])<<std::endl<<"                DeltaR with subleading = "<<(float)deltaR(h.eta(), h.phi(), eta[1], phi[1])<<std::endl<<"                DeltaR with third = "<<(float)deltaR(h.eta(), h.phi(), eta[2], phi[2])<<std::endl;
      // 	  par_pt.push_back(h.pt());
      // 	  par_eta.push_back(h.eta());
      // 	  par_phi.push_back(h.phi());
      // 	  par_flav.push_back(h.pdfId());
      // 	}
      
      // update leading track
      if( h.charge() != 0 && fabs(h.eta()) < 2.5 && lead_track_pt < h.pt() )
	lead_track_pt = h.pt(); 

      if( doSpecies_ && ! isInSpecies(h) ) continue;
       
      hist2D_["hspectrum2D"]->Fill(h.eta(),h.pt());
      if( h.charge() != 0 ) hist2D_["chspectrum2D"]->Fill(h.eta(),h.pt());
      if( h.charge() > 0 ) hist2D_["pchspectrum2D"]->Fill(h.eta(),h.pt());

      if( isInEtaRange(h, hEtaMin_, hEtaMax_) )
	{
	  hist_["hspectrum"]->Fill(h.pt());
	  if( h.charge() != 0 ) hist_["chspectrum"]->Fill(h.pt());
	  if( h.charge() != 0 ) hist2D_["ch_qscale2D"]->Fill(genEvtInfo->qScale(),h.pt());
	  if( h.charge() > 0 ) hist_["pchspectrum"]->Fill(h.pt());

	  // associate with the highest-pt jet for which 
	  // the track is found in the jet cone
     
	  double maxPtJet = 0.;
	  if( doFlavor_)
	    {
	      for( const auto & mjp : *fcol )
		{
		  const reco::Jet *aJet = mjp.first.get();
		  if ( !isInEtaRange(*aJet,jetEtaMin_,jetEtaMax_) ) continue;
		  if( !isInFlavor( mjp.second ) ) continue;
		  double dr = deltaR(*aJet,h);
		  if( dr < jetRadius_ && aJet->pt() > maxPtJet)
		    maxPtJet = aJet->pt();
		}
	    }
	  else
	    {
	      for( const auto & jet : *gcol )
		{
		  if( ! isInEtaRange( jet, jetEtaMin_, jetEtaMax_ ) ) continue; 
		  double dr = deltaR(jet,h);
		  if( dr < jetRadius_ && jet.pt() > maxPtJet)
		    maxPtJet = jet.pt();
		}
	    }
	  hist2D_["ffmatrix"]->Fill( maxPtJet, h.pt());
	  if( h.charge() != 0 ) hist2D_["cffmatrix"]->Fill( maxPtJet, h.pt());
	  if( h.charge() > 0 ) hist2D_["pcffmatrix"]->Fill( maxPtJet, h.pt());
	}
    }

  hist_["lead_track"]->Fill(lead_track_pt);

  // std::cout<<"End of Event"<<std::endl;
  // std::cout<<"********************************* "<<std::endl;

}


bool 
QCDAnalyzer::isInEtaRange( const reco::Candidate& c, double etaMin, double etaMax )
{
  if( useRapidity_ == false &&  c.eta() <= etaMax && c.eta() >= etaMin )
    return true;
  if( useRapidity_ == true &&  c.y() <= etaMax && c.y() >= etaMin )
    return true;
  return false;
}

bool
QCDAnalyzer::isInFlavor( const reco::MatchedPartons & aMatch )
{
  int flavor = 0;
  if( aMatch.heaviest().isNonnull() )
    {
      flavor = aMatch.heaviest().get()->pdgId();
    }
  for( const auto inFlavor : flavorId_ )
    { if( flavor == inFlavor ) return true; }
 
  return false;
}

int
QCDAnalyzer::GetFlavor( const reco::MatchedPartons & aMatch )
{
  int flavor = 0;
  if( aMatch.heaviest().isNonnull() )
    {
      flavor = aMatch.heaviest().get()->pdgId();
    }
  return flavor;  
}

// bool
// QCDAnalyzer::isInFlavor_LeadJet( const reco::MatchedPartons & aMatch )
// {
//   int flavor = 0;
//   if( aMatch.heaviest().isNonnull() )
//     {
//       flavor = aMatch.heaviest().get()->pdgId();
//     }
//   for( const auto inFlavor : flavorId_LeadJet_ )
//     { if( flavor == inFlavor ) return true; }
 
//   return false;
// }

// bool
// QCDAnalyzer::isInFlavor_SubLeadJet( const reco::MatchedPartons & aMatch )
// {
//   int flavor = 0;
//   if( aMatch.heaviest().isNonnull() )
//     {
//       flavor = aMatch.heaviest().get()->pdgId();
//     }
//   for( const auto inFlavor : flavorId_SubLeadJet_ )
//     { if( flavor == inFlavor ) return true; }
 
//   return false;
// }

bool 
QCDAnalyzer::isInSpecies( const reco::GenParticle & h )
{
  for( const auto inSpecies : speciesId_ )
    { if( h.pdgId() == inSpecies ) return true; }

  return false;  
}

bool 
QCDAnalyzer::isDSEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  bool posDS = false; bool negDS = false;

  edm::ESHandle<ParticleDataTable> particleDataTable_;
  iSetup.getData(particleDataTable_);

  Handle<reco::GenParticleCollection> gcol;
  iEvent.getByLabel(genParticleSrc_, gcol);
  for( const auto & gen : * gcol )
    {
      // see if genpartice counts for DS
      HepPDT::ParticleID particleID(gen.pdgId());
      if (particleID.isValid())
	{
	  ParticleData const * particleData = particleDataTable_->particle(particleID);
	  if (particleData)
	    { 
	      double tau =  particleDataTable_->particle(particleID)->lifetime();
	      if ( tau  > 1e-18 || tau == 0.0 )
		{
		  if( gen.energy() > 3.0 && gen.eta() > 3.0 && gen.eta() < 5.0 ) posDS = true;
		  if( gen.energy() > 3.0 && gen.eta() < -3.0 && gen.eta() > -5.0 ) negDS = true;
		}
	    }
	}
    }
  if( posDS && negDS ) return true;
  else return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
QCDAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QCDAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
QCDAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
QCDAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
QCDAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
QCDAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QCDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(QCDAnalyzer);
