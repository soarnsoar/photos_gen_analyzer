// -*- C++ -*-
//
// Package:    photos_gen_analyzer/photos_gen_analyzer
// Class:      photos_gen_analyzer
//
/**\class photos_gen_analyzer photos_gen_analyzer.cc photos_gen_analyzer/photos_gen_analyzer/plugins/photos_gen_analyzer.cc

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



class photos_gen_analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit photos_gen_analyzer(const edm::ParameterSet&);
      ~photos_gen_analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<GenParticleCollection> genParticles_Token;
  edm::EDGetTokenT<GenEventInfoProduct> genInfo_Token;
  edm::EDGetTokenT<LHEEventProduct> LHEInfo_Token;

  TTree *photoninfo;
  TTree *lep1;
  TTree *lep2;
  TTree *weights;

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
photos_gen_analyzer::photos_gen_analyzer(const edm::ParameterSet& iConfig)
// :
// tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{
   //now do what ever initialization is needed
  usesResource("TFileService");//genParticles
  genParticles_Token = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
  genInfo_Token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  LHEInfo_Token = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

}


photos_gen_analyzer::~photos_gen_analyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
photos_gen_analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  int leppid=13;
  using namespace edm;
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genParticles_Token, genParticles);//genParticle

   edm::Handle<GenEventInfoProduct> genInfo;
   iEvent.getByToken(genInfo_Token, genInfo);
   edm::Handle<LHEEventProduct> LHEInfo;
   iEvent.getByToken(LHEInfo_Token, LHEInfo);


   double   weight=genInfo->weight();
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

   double lep1_px=0, lep1_py=0, lep1_pz=0, lep1_ee=0;
   double lep2_px=0, lep2_py=0, lep2_pz=0, lep2_ee=0;
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

   int nphoton=0;
   
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
     double dR1 = vfsr.DeltaR(v1);
     double dR2 = vfsr.DeltaR(v2);
     if(!fromhardDY) continue;
     fsr_px=fsr_px+px;
     fsr_py=fsr_py+py;
     fsr_pz=fsr_pz+pz;
     fsr_ee=fsr_ee+ee;
     
     

   }//for all genparticles
   vboson.SetPxPyPzE(lep1_px+lep2_px+fsr_px, lep1_py+lep2_py+fsr_py,lep1_pz+lep2_pz+fsr_pz,lep1_ee+lep2_ee+fsr_ee);

   cout<<"####event info###"<<endl;
   cout<<"lep1_i="<<lep1_i<<" lep2_i="<<lep2_i<<endl;
   cout<<"vboson.M()="<<vboson.M()<<" vboson.Perp()="<<vboson.Perp()<<endl;
   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int status = p.status();
     double px = p.px();
     double py = p.py();
     double pz = p.pz();
     double ee = p.energy();
     int mother = gen_motherindex[i];
     int fromhardDY = gen_fromhardDY[i];
     if(fromhardDY)     cout<<"i="<<i<<" id="<<id<<" status="<<status<<" mother="<<mother<<" fromhardDY="<<fromhardDY<<" px="<<px <<" py="<<py<<" pz="<<pz<<" ee="<<ee<<endl;
   }


}


// ------------ method called once each job just before starting event loop  ------------
void
photos_gen_analyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
photos_gen_analyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
photos_gen_analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(photos_gen_analyzer);
