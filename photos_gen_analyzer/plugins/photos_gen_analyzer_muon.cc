// -*- C++ -*-
//
// Package:    photos_gen_analyzer_muon/photos_gen_analyzer_muon
// Class:      photos_gen_analyzer_muon
//
/**\class photos_gen_analyzer_muon photos_gen_analyzer_muon.cc photos_gen_analyzer_muon/photos_gen_analyzer_muon/plugins/photos_gen_analyzer_muon.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  JunHo Choi
//         Created:  Thu, 04 Oct 2018 18:28:04 GMT
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
 #include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using namespace edm;
using namespace reco;
using namespace std;

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
//#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"

#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>



class photos_gen_analyzer_muon : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit photos_gen_analyzer_muon(const edm::ParameterSet&);
      ~photos_gen_analyzer_muon();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<GenParticleCollection> genParticles_Token;
  edm::EDGetTokenT<GenEventInfoProduct> genInfo_Token;
  edm::EDGetTokenT<LHEEventProduct> LHEInfo_Token;

  
  TTree *Tphoton;
  TTree *Tlep1;
  TTree *Tlep2;
  TTree *Tboson; 
  TTree *Tweight;

  vector<double> photon_px, photon_py, photon_pz, photon_ee; vector<int> photon_isfsr;
  double lep1_px, lep1_py, lep1_pz, lep1_ee,lep1_dRboson,lep1_dPhiboson;
  double lep2_px, lep2_py, lep2_pz, lep2_ee,lep2_dRboson,lep2_dPhiboson;
  double boson_px, boson_py, boson_pz, boson_ee; //int isgamma, isz;   

  double dR1, dR2, dPhi1, dPhi2, dRboson, dPhiboson;
  double dR12, dPhi12, PT12, M12;
  double weight;



  ///Histograms for photons/////
  TH1D *h_dR1, *h_dR2, *h_dRboson;
  TH1D *h_dPhi1, *h_dPhi2, *h_dPhiboson;
  TH1D *h_Ephoton, *h_Nphoton;
 
  TH1D *h_dR1_fsr, *h_dR2_fsr, *h_dRboson_fsr;
  TH1D *h_dPhi1_fsr, *h_dPhi2_fsr,  *h_dPhiboson_fsr;
  TH1D *h_Ephoton_fsr, *h_Nphoton_fsr;


  /////Histograms for lepton pair and DYboson
  TH1D *h_dR12, *h_dPhi12, *h_PT12, *h_M12;
  TH1D *h_M_boson, *h_PT_boson;
  TH1D *h_PT_lep1, *h_dRboson_lep1;
  TH1D *h_PT_lep2, *h_dRboson_lep2;
  //double lep1_px, lep1_py, lep1_pz, lep1_ee; 
  // double lep2_px, lep2_py, lep2_pz, lep2_ee;
  //double boson_px, boson_py, boson_pz, boson_ee; //int isgamma, isz;
  /*
  TTree *Tweight; double weight;
  TTree *Tphoton; double dR1,dR2,dRboson,dPhi1,dPhi2,dPhiboson; int isfsr;
  TTree *Tlep1; double lep1_dRboson,lep1_dPhiboson;
  TTree *Tlep2; double lep2_dRboson,lep2_dPhiboson;
  TTree *Tboson; double dR12, dPhi12;
  */
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
photos_gen_analyzer_muon::photos_gen_analyzer_muon(const edm::ParameterSet& iConfig)
// :
// tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
   //now do what ever initialization is needed
  usesResource("TFileService");//genParticles
  genParticles_Token = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
  genInfo_Token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  LHEInfo_Token = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

}


photos_gen_analyzer_muon::~photos_gen_analyzer_muon()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
photos_gen_analyzer_muon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  photon_px.clear(), photon_py.clear(), photon_pz.clear(), photon_ee.clear(); photon_isfsr.clear();
  lep1_px=0, lep1_py=0, lep1_pz=0, lep1_ee=0; 
   lep2_px=0, lep2_py=0, lep2_pz=0, lep2_ee=0;
   boson_px=0, boson_py=0, boson_pz=0, boson_ee=0; //int isgamma, isz;
   weight=1; 

  int leppid=13;
  using namespace edm;
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genParticles_Token, genParticles);//genParticle

   edm::Handle<GenEventInfoProduct> genInfo;
   iEvent.getByToken(genInfo_Token, genInfo);
   edm::Handle<LHEEventProduct> LHEInfo;
   iEvent.getByToken(LHEInfo_Token, LHEInfo);


   weight = genInfo->weight();
   // Tweight->Fill();
   
vector<int> hardpid;
   hardpid.push_back(23);//Z boson 
   hardpid.push_back(22);//Z boson
   hardpid.push_back(leppid);
   hardpid.push_back(-leppid);
   vector<int> hardindex;
   int lep1_i=-99;
   int lep2_i=-99;
   int gensize= genParticles->size();
   int hardpidsize=hardpid.size();
   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int status = p.status();
     for(int j = 0 ; j < hardpidsize; j++){
       
       if(id==hardpid[j] && status>20 && status<25) hardindex.push_back(i);

     }
     
   }//End of push_back all hard particle.

   //Now just find if a particle's anscestor is hadparticle.
int hardindexsize=   hardindex.size();

   //Before that, let's make mother indices
   vector<int> gen_motherindex;
   
   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];
     
     int mother = -1;
     
     for( reco::GenParticleCollection::const_iterator mit = genParticles->begin(); mit != genParticles->end(); ++mit ) {
       if( p.mother()==&(*mit) ) {
         mother = std::distance(genParticles->begin(),mit);
         break;
       }
     }
     gen_motherindex.push_back(mother);




   }//end of find motherindex
   vector<int> gen_fromhardDY;

   for(int i = 0; i < gensize; ++ i) {

     int fromhardDY=0;
     int mother = gen_motherindex[i];

     if(mother<0){
       gen_fromhardDY.push_back(fromhardDY);
       continue;
     }
     while(((*genParticles)[mother]).status() != 4){
  

       
       for(int i_hard = 0 ; i_hard<hardindexsize; i_hard++){
	 if(mother==hardindex[i_hard]){
	   fromhardDY=1;	   
	 }	 
       }//for all hardparticles
       mother=gen_motherindex[mother];
       if(mother==-1) break;
     }
     gen_fromhardDY.push_back(fromhardDY);
   }//for all genptls



   //Let's find the DY leptons

   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int fromhardDY=gen_fromhardDY[i];
     if(fabs(id)==15 && fromhardDY) return; //veto tautau event

   }


   vector<int> vlep1_i;//why use container? => there could be an event where conversion pair is.
   vector<int> vlep2_i;

   lep1_px=0, lep1_py=0, lep1_pz=0, lep1_ee=0;
   lep2_px=0, lep2_py=0, lep2_pz=0, lep2_ee=0;
   
   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int status = p.status();
     int fromhardDY=gen_fromhardDY[i];
     int id = p.pdgId();
     if(status==1 && fromhardDY==1){
       if(id==leppid){
	 vlep1_i.push_back(i);
       }
       else if(id==-leppid){
	 vlep2_i.push_back(i);
       }
     }//if status==1 from DY
   }//for all genparticles

  
   if(vlep1_i.size()==0) return; //If it's not the channel we want
   else if(vlep1_i.size()==1) {
     lep1_i=vlep1_i[0];
     const GenParticle & p = (*genParticles)[vlep1_i[0]];
     double px = p.px();
     double py = p.py();
     double pz = p.pz();
     double ee = p.energy();
     lep1_px=px;lep1_py=py;lep1_pz=pz;lep1_ee=ee;

   }
   else if(vlep1_i.size() >1){
     double leppt=0;    //prompt== the muon whose pt is largest among muons    
     int vlep1size=vlep1_i.size();
     for(int i =0 ; i < vlep1size ; i++){
       cout<<"vlep1_i[i]="<<vlep1_i[i]<<endl;
       const GenParticle & p = (*genParticles)[vlep1_i[i]];
       double px = p.px();
       double py = p.py();
       double pz = p.pz();
       double ee = p.energy();
       TLorentzVector vtemp;
       vtemp.SetPxPyPzE(px,py,pz,ee);
       if(vtemp.Perp()>leppt){
         leppt=  vtemp.Perp();
         lep1_i= vlep1_i[i];
         lep1_px=px;lep1_py=py;lep1_pz=pz;lep1_ee=ee;
       }//end of if pt>pt 
     }//end of for vlep1  
   }//end of if sizevlep1>1  

   if(vlep2_i.size()==0) return; //If it's not the channel we want 
   else if(vlep2_i.size()==1){
     lep2_i=vlep2_i[0];
     const GenParticle & p = (*genParticles)[vlep2_i[0]];
     double px = p.px();
     double py = p.py();
     double pz = p.pz();
     double ee = p.energy();
     lep2_px=px;lep2_py=py;lep2_pz=pz;lep2_ee=ee;

   }
   else if(vlep2_i.size() >1){
     double leppt=0;    //prompt== the muon whose pt is largest among muons   
     int vlep2size=vlep2_i.size();
     for(int i =0 ; i < vlep2size ; i++){
       const GenParticle & p = (*genParticles)[vlep2_i[i]];
       double px = p.px();
       double py = p.py();
       double pz = p.pz();
       double ee = p.energy();
       TLorentzVector vtemp;
       vtemp.SetPxPyPzE(px,py,pz,ee);
       if(vtemp.Perp()>leppt){
         leppt=  vtemp.Perp();
         lep2_i=vlep2_i[i];
         lep2_px=px;lep2_py=py;lep2_pz=pz;lep2_ee=ee;
       }//end of if pt>pt   
     }//end of for all lep2
   }//end of if size >1


   TLorentzVector v1,v2, vboson;
   v1.SetPxPyPzE(lep1_px,lep1_py,lep1_pz,lep1_ee);
   v2.SetPxPyPzE(lep2_px,lep2_py,lep2_pz,lep2_ee);
   
   
   


   
   double fsr_px = 0, fsr_py=0, fsr_pz=0, fsr_ee=0;
   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];
     //     int id = p.pdgId();
     int status = p.status();
     double px = p.px();
     double py = p.py();
     double pz = p.pz();
     double ee = p.energy();
     int mother = gen_motherindex[i];
     int fromhardDY = gen_fromhardDY[i];
     TLorentzVector vfsr;
     vfsr.SetPxPyPzE(px,py,pz,ee);
     if(status!=1) continue;
     if(mother==-1) continue;
     if(i==lep1_i) continue;
     if(i==lep2_i) continue;
     if(!fromhardDY) continue;
     fsr_px=fsr_px+px;
     fsr_py=fsr_py+py;
     fsr_pz=fsr_pz+pz;
     fsr_ee=fsr_ee+ee;
     
     

   }//for all genparticles
   boson_px=lep1_px+lep2_px+fsr_px;
   boson_py=lep1_py+lep2_py+fsr_py;
   boson_pz=lep1_pz+lep2_pz+fsr_pz;
   boson_ee=lep1_ee+lep2_ee+fsr_ee;

   vboson.SetPxPyPzE(boson_px,boson_py,boson_pz,boson_ee);
   
   lep1_dRboson=v1.DeltaR(vboson);
   lep2_dRboson=v2.DeltaR(vboson);

   lep1_dPhiboson=v1.DeltaPhi(vboson);
   lep2_dPhiboson=v2.DeltaPhi(vboson);

   dR12=v1.DeltaR(v2);

   dPhi12=v1.DeltaPhi(v2);

   PT12=(v1+v2).Perp();
   M12=(v1+v2).M();
   h_dR12->Fill(dR12,weight);   
   h_dPhi12->Fill(dPhi12,weight);
   h_PT12->Fill(PT12,weight);
   h_M12->Fill(M12,weight);

   h_M_boson->Fill(vboson.M(),weight);
   h_PT_boson->Fill(vboson.Perp(),weight);

   h_PT_lep1->Fill(v1.Perp(),weight);
   h_dRboson_lep1->Fill(lep1_dPhiboson,weight);

   h_PT_lep2->Fill(v1.Perp(),weight);
   h_dRboson_lep2->Fill(lep2_dPhiboson,weight);

   //   cout<<"####event info###"<<endl;
   //cout<<"lep1_i="<<lep1_i<<" lep2_i="<<lep2_i<<endl;
   //cout<<"vboson.M()="<<vboson.M()<<" vboson.Perp()="<<vboson.Perp()<<endl;
  
   //   double dR12 = v1.DeltaR(v2);
   // double dPhi12=v1.DeltaPhi(v2);
  
   //Do analysis//
   int nphoton=0;
   int nphoton_fsr=0;
   for(int i = 0; i < gensize; ++ i) {
     // double dR1=-99, dR2=-99, dRboson=-99;
     //  double dPhi1=-99, dPhi2=-99, dPhiboson=-99;
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int status = p.status();
     double px = p.px();
     double py = p.py();
     double pz = p.pz();
     double ee = p.energy();
     //      TLorentzVector v; v.SetPxPyPzE(px,py,pz,ee);
     //double mm=v.M();
     // int mother = gen_motherindex[i];
      int fromhardDY = gen_fromhardDY[i];
     //     if(fromhardDY)     cout<<"i="<<i<<" id="<<id<<" status="<<status<<" mother="<<mother<<" fromhardDY="<<fromhardDY<<" mm="<<mm<<" px="<<px <<" py="<<py<<" pz="<<pz<<" ee="<<ee<<endl;


     if((id==22)&&(status==1)){//photons
       nphoton=nphoton+1;
       
       TLorentzVector v; v.SetPxPyPzE(px,py,pz,ee);
       dR1=v1.DeltaR(v);
       dR2=v2.DeltaR(v);
       dRboson=vboson.DeltaR(v);
       
       dPhi1=v1.DeltaR(v);
       dPhi2=v2.DeltaR(v);
       dPhiboson=vboson.DeltaR(v);
       
       
       
       photon_px.push_back(px);
       photon_py.push_back(py);
       photon_pz.push_back(pz);
       photon_ee.push_back(ee);

    
       h_dR1->Fill(dR1,weight);
       h_dR2->Fill(dR2,weight);
       h_dRboson->Fill(dRboson,weight);

       h_dPhi1->Fill(dPhi1,weight);
       h_dPhi2->Fill(dPhi2,weight);
       h_dPhiboson->Fill(dPhiboson,weight);

       h_Ephoton->Fill(ee,weight);

       if(fromhardDY){
	 photon_isfsr.push_back(1);
	 nphoton_fsr=nphoton_fsr+1;

	 h_dR1_fsr->Fill(dR1,weight);
	 h_dR2_fsr->Fill(dR2,weight);
	 h_dRboson_fsr->Fill(dRboson,weight);

	 h_dPhi1_fsr->Fill(dPhi1,weight);
	 h_dPhi2_fsr->Fill(dPhi2,weight);
	 h_dPhiboson_fsr->Fill(dPhiboson,weight);

	 h_Ephoton_fsr->Fill(ee,weight);


       }
       else{
	 photon_isfsr.push_back(0); 
       }       
       


     }//end of if photon
   }//end of all genptls


   h_Nphoton->Fill(nphoton,weight);
   h_Nphoton_fsr->Fill(nphoton_fsr,weight);
   // cout<<"nphoton_fsr="<<nphoton_fsr<<endl;
   Tphoton->Fill();
   Tlep1->Fill();
   Tlep2->Fill();
   Tboson->Fill();
   Tweight->Fill();

   /*
   cout<<"######"<<endl;
   cout<<"photon_px.size()="<<photon_px.size()<<endl;
   cout<<"photon_py.size()="<<photon_py.size()<<endl;
   cout<<"photon_pz.size()="<<photon_pz.size()<<endl;
   cout<<"photon_ee.size()="<<photon_ee.size()<<endl;
   cout<<"photon_isfsr.size()="<<photon_isfsr.size()<<endl;
   */
}


// ------------ method called once each job just before starting event loop  ------------
void
photos_gen_analyzer_muon::beginJob()
{
  edm::Service<TFileService> fs;
  Tphoton = fs->make<TTree>("Tphoton","Tphoton");
  Tlep1 = fs->make<TTree>("Tlep1","Tlep1");
  Tlep2 = fs->make<TTree>("Tlep2","Tlep2");
  Tboson = fs->make<TTree>("Tboson","Tboson");
  Tweight = fs->make<TTree>("Tweight","Tweight");



  h_dR1 = fs->make<TH1D>("h_dR1","h_dR1",100,0,0.1);
  h_dR2 = fs->make<TH1D>("h_dR2","h_dR2",100,0,0.1);
  h_dRboson = fs->make<TH1D>("h_dRboson","h_dRboson",100,0,5);
 
  h_dPhi1 = fs->make<TH1D>("h_dPhi1","h_dPhi1",100,0,0.1);
  h_dPhi2 = fs->make<TH1D>("h_dPhi2","h_dPhi2",100,0,0.1);
  h_dPhiboson = fs->make<TH1D>("h_dPhiboson","h_dPhiboson",100,0,5);

  h_Ephoton = fs->make<TH1D>("h_Ephoton","h_Ephoton",1000,0,100);
  h_Nphoton = fs->make<TH1D>("h_Nphoton","h_Nphoton",500,0,500);

  h_dR1_fsr = fs->make<TH1D>("h_dR1_fsr","h_dR1_fsr",100,0,0.1);
  h_dR2_fsr = fs->make<TH1D>("h_dR2_fsr","h_dR2_fsr",100,0,0.1);
  h_dRboson_fsr = fs->make<TH1D>("h_dRboson_fsr","h_dRboson_fsr",100,0,5);

  h_dPhi1_fsr = fs->make<TH1D>("h_dPhi1_fsr","h_dPhi1_fsr",100,0,0.1);
  h_dPhi2_fsr = fs->make<TH1D>("h_dPhi2_fsr","h_dPhi2_fsr",100,0,0.1);
  h_dPhiboson_fsr = fs->make<TH1D>("h_dPhiboson_fsr","h_dPhiboson_fsr",100,0,5);

  h_Ephoton_fsr = fs->make<TH1D>("h_Ephoton_fsr","h_Ephoton_fsr",1000,0,100);
  h_Nphoton_fsr = fs->make<TH1D>("h_Nphoton_fsr","h_Nphoton_fsr",10,0,10);

  h_dR12 = fs->make<TH1D>("h_dR12","h_dR12",100,0,5);
  h_dPhi12 = fs->make<TH1D>("h_dPhi12","h_dPhi12",100,0,5);
  h_PT12 = fs->make<TH1D>("h_PT12","h_PT12",300,0,300);
  h_M12 = fs->make<TH1D>("h_M12","h_M12",500,0,500);
  
  h_M_boson = fs->make<TH1D>("h_M_boson","h_M_boson",500,0,500);
  h_PT_boson = fs->make<TH1D>("h_PT_boson","h_PT_boson",300,0,300);

  h_PT_lep1 = fs->make<TH1D>("h_PT_lep1","h_PT_lep1",300,0,300);
  h_dRboson_lep1 = fs->make<TH1D>("h_dRboson_lep1","h_dRboson_lep1",100,0,5);

  h_PT_lep2 = fs->make<TH1D>("h_PT_lep2","h_PT_lep2",300,0,300);
  h_dRboson_lep2 = fs->make<TH1D>("h_dRboson_lep2","h_dRboson_lep2",100,0,5);

  /*
  Tlep1->Branch("lep1_dRboson",&lep1_dRboson,"lep1_dRboson/D");
  Tlep1->Branch("lep1_dPhiboson",&lep1_dPhiboson,"lep1_dPhiboson/D");
  
  Tlep2->Branch("lep2_dRboson",&lep2_dRboson,"lep2_dRboson/D");
  Tlep2->Branch("lep2_dPhiboson",&lep2_dPhiboson,"lep2_dPhiboson/D");

  Tboson->Branch("dR12",&dR12,"dR12/D");
  Tboson->Branch("dPhi12",&dPhi12,"dPhi12/D");
  */
  
  Tlep1->Branch("px",&lep1_px,"px/D");
  Tlep1->Branch("py",&lep1_py,"py/D");
  Tlep1->Branch("pz",&lep1_pz,"pz/D");
  Tlep1->Branch("ee",&lep1_ee,"ee/D");
  
  Tlep2->Branch("px",&lep2_px,"px/D");
  Tlep2->Branch("py",&lep2_py,"py/D");
  Tlep2->Branch("pz",&lep2_pz,"pz/D");
  Tlep2->Branch("ee",&lep2_ee,"ee/D");

  Tboson->Branch("px",&boson_px,"px/D");
  Tboson->Branch("py",&boson_py,"py/D");
  Tboson->Branch("pz",&boson_pz,"pz/D");
  Tboson->Branch("ee",&boson_ee,"ee/D");
  
  // Tboson->Branch("isgamma", &boson_isgamma, "isgamma/I");
  // Tboson->Branch("isz", &boson_isz, "isz/I");

  Tweight->Branch("weight",&weight,"weight/D");
  /*
  Tphoton->Branch("dR1",&dR1, "dR1/D");
  Tphoton->Branch("dR2",&dR2, "dR2/D");
  Tphoton->Branch("dRboson",&dRboson, "dRboson/D");
 
  Tphoton->Branch("dPhi1",&dPhi1, "dPhi1/D");
  Tphoton->Branch("dPhi2",&dPhi2, "dPhi2/D");
  Tphoton->Branch("dPhiboson",&dPhiboson, "dPhiboson/D");

  Tphoton->Branch("isfsr",&isfsr, "isfsr/I");

  */


  
  Tphoton->Branch("px","vector<double>",&photon_px);
  Tphoton->Branch("py","vector<double>",&photon_py);
  Tphoton->Branch("pz","vector<double>",&photon_pz);
  Tphoton->Branch("ee","vector<double>",&photon_ee);
  Tphoton->Branch("isfsr","vector<int>",&photon_isfsr);
  
  }

// ------------ method called once each job just after ending the event loop  ------------
void
photos_gen_analyzer_muon::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
photos_gen_analyzer_muon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(photos_gen_analyzer_muon);
