#include <Normalizer.cc>
class JUnfold{

  //X :GEN / Y :RECO

private :
  
  TH2D *_hGenRecoMC;
  //TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
  TH1D *_hGenMC;
  TH1D *_hRecoMC;
  

  int _nbin_Gen = 0, _nbin_Reco=0;
  double _xmin_Gen=0, _xmin_Reco;
  double _xmax_Gen=0, _xmax_Reco;

  TString _hGenRecoMC_name="";
  TString _hGenMC_name="";
  TString _hRecoMC_name="";



public :

  void Set_nbin_xmin_xmax_Gen(int nbin_Gen, double xmin_Gen, double xmax_Gen);
  void Set_nbin_xmin_xmax_Reco(int nbin_Reco, double xmin_Reco, double xmax_Reco);  
  void Init_histos();
  void SetName_hGenRecoMC(TString name);
  void SetName_hGenMC(TString name);
  void SetName_hRecoMC(TString name);

  void FillhMC(double xgen, double xreco, double weight);//Fill(Double_t x, Double_t y, Double_t w);
  void ImportHists(TString dirname);
  

};

void JUnfold::Set_nbin_xmin_xmax_Gen(int nbin_Gen, double xmin_Gen, double xmax_Gen){
  _nbin_Gen=nbin_Gen; _xmin_Gen=xmin_Gen; _xmax_Gen=xmax_Gen; 
}
void JUnfold::Set_nbin_xmin_xmax_Reco(int nbin_Reco, double xmin_Reco, double xmax_Reco){
  _nbin_Reco=nbin_Reco; _xmin_Reco=xmin_Reco; _xmax_Reco=xmax_Reco;
}
void JUnfold::Init_histos(){
  // _hGenRecoMC->SetName(_hGenRecoMC_name);
  //_hGenMC->SetName(_hGenMC_name);
  //_hRecoMC->SetName(_hRecoMC_name);
  _hGenRecoMC=new TH2D(_hGenRecoMC_name,_hGenRecoMC_name,_nbin_Gen,_xmin_Gen,_xmax_Gen,_nbin_Reco,_xmin_Reco,_xmax_Reco);
  _hGenMC=new TH1D(_hGenMC_name,_hGenMC_name,_nbin_Gen,_xmin_Gen,_xmax_Gen);
  _hRecoMC=new TH1D(_hRecoMC_name,_hRecoMC_name,_nbin_Reco,_xmin_Reco,_xmax_Reco);
}
void JUnfold::SetName_hGenRecoMC(TString name){
  _hGenRecoMC_name=name;
  //  _hGenRecoMC->SetName(name);
}

void JUnfold::SetName_hGenMC(TString name){
  _hGenMC_name=name;  
//  _hGenMC->SetName(name);
}

void JUnfold::SetName_hRecoMC(TString name){
  _hRecoMC_name=name;
  //  _hRecoMC->SetName(name);
}



void JUnfold::FillhMC(double xgen, double xreco, double weight=1){
  _hGenRecoMC->Fill(xgen,xreco,weight);
  _hGenMC->Fill(xgen,weight);
  _hRecoMC->Fill(xreco,weight);
}


void JUnfold::ImportHists(TString dirname){

  gSystem->Exec("mkdir -p "+dirname);
  unsigned int Nhisto=3;
  TList **hList = new TList*[Nhisto];

  TCanvas *c2 = new TCanvas(); 
  for(int i =0 ; i < Nhisto; i++){
   
    hList[i] = new TList();
    TString hname="";

    if(i==0) hname= _hGenRecoMC->GetName();
    else if(i==1) hname= _hGenMC->GetName();
    else if(i==2) hname= _hRecoMC->GetName();
    TFile ft(dirname+"/"+hname+".root","recreate");


    if(i==0){
      TCanvas *c = new TCanvas(); 
      hList[i]->Add(_hGenRecoMC);
      _hGenRecoMC->Draw("colz");
      c->SetLogz();
      c->SaveAs(dirname+"/"+hname+"_image.pdf");
    }
    else if(i==1){
      
      hList[i]->Add(_hGenMC);
     
      
    }
    else if(i==2){
      hList[i]->Add(_hRecoMC);
      TCanvas *c2 = new TCanvas();
      _hGenMC->Draw();
      c2->Update();
      TPaveStats *st = ((TPaveStats *)(gPad->GetPrimitive("stats"))); 
      st->SetY1NDC(st->GetY1NDC() - 0.1);
      st->SetY2NDC(st->GetY2NDC() - 0.1);

      _hRecoMC->Draw("same");
      _hRecoMC->SetLineColor(kRed);
      c2->SetLogy();
      c2->Update();
      TPaveStats *st2 = ((TPaveStats *)(gPad->GetPrimitive("stats")));
      st2->SetY1NDC(st2->GetY1NDC() - 0.2);
      st2->SetY2NDC(st2->GetY2NDC() - 0.2);
      c2->Update();
      st2->Draw();
      c2->SaveAs(dirname+"/"+_hGenMC->GetName()+"_"+_hRecoMC->GetName()+"_image.pdf");
    }    

    hList[i]->Write();

    ft.Write();

    ft.Close();

    

  } 
  

}



class ISR_preFSR : public JUnfold{
private:
  TFile *_file;
  TTree *_lep1;
  TTree *_lep2;
  TTree *_fsr_t;
  TTree *_weights;

  double _pt_pre=0, _mass_pre=0;
  double _pt_post=0, _mass_post=0;
  double _weight=1;
  
  double _lep1_px, _lep1_py, _lep1_pz, _lep1_ee;
  double _lep2_px, _lep2_py, _lep2_pz, _lep2_ee;
  double _fsr_t_px, _fsr_t_py, _fsr_t_pz, _fsr_t_ee;

  unsigned int _Nevent=0;

public:
  //  void Set_Pt(TString filename, TString treename);
  // void Set_Mass();
  //void Set_Weight();
  //  void SetPtMassWeight(TString filename);
  void LoadPW(TString filename);
  void GetEntrySetPtMass(unsigned int i);
  double GetNevent();

  
  double Get_pt_pre();
  double Get_pt_post();
  double Get_mass_pre();
  double Get_mass_post();
  double Get_weight();

  void runMass(vector<TString> filenames, TString tag);

};

void ISR_preFSR::LoadPW(TString filename){
  _file=TFile::Open(filename);
  _lep1=(TTree*)_file->Get("DY/lep1");
  _lep2=(TTree*)_file->Get("DY/lep2");
  _fsr_t=(TTree*)_file->Get("DY/fsr_t");
  _weights=(TTree*)_file->Get("DY/truweight");
  


 
  _lep1->SetBranchAddress("px",&_lep1_px);
  _lep1->SetBranchAddress("py",&_lep1_py);
  _lep1->SetBranchAddress("pz",&_lep1_pz);
  _lep1->SetBranchAddress("ee",&_lep1_ee);

  _lep2->SetBranchAddress("px",&_lep2_px);
  _lep2->SetBranchAddress("py",&_lep2_py);
  _lep2->SetBranchAddress("pz",&_lep2_pz);
  _lep2->SetBranchAddress("ee",&_lep2_ee);

  _fsr_t->SetBranchAddress("px",&_fsr_t_px);
  _fsr_t->SetBranchAddress("py",&_fsr_t_py);
  _fsr_t->SetBranchAddress("pz",&_fsr_t_pz);
  _fsr_t->SetBranchAddress("ee",&_fsr_t_ee);



  _weights->SetBranchAddress("weight",&_weight);
  _Nevent = _weights->GetEntries();
  
  

}

void ISR_preFSR::GetEntrySetPtMass(unsigned int i){
 

  _lep1->GetEntry(i);
  _lep2->GetEntry(i);
  _weights->GetEntry(i);
  _fsr_t->GetEntry(i);

  //  double _pt_pre=0, _mass_pre=0;
  //double _pt_post=0, _mass_post=0;
  TLorentzVector v1,v2,vfsr;
  v1.SetPxPyPzE(_lep1_px,_lep1_py,_lep1_pz,_lep1_ee);
  v2.SetPxPyPzE(_lep2_px,_lep2_py,_lep2_pz,_lep2_ee);
  vfsr.SetPxPyPzE(_fsr_t_px,_fsr_t_py,_fsr_t_pz,_fsr_t_ee);

  _pt_pre=(v1+v2+vfsr).Perp();
  _pt_post=(v1+v2).Perp();
  _mass_pre=(v1+v2+vfsr).M();
  _mass_post=(v1+v2).M();

  _Nevent = _weights->GetEntries();
  
}

double ISR_preFSR::GetNevent(){
  return _Nevent;
}

double ISR_preFSR::Get_pt_pre(){
  return _pt_pre;
}
double ISR_preFSR::Get_pt_post(){
  return _pt_post;
}
double ISR_preFSR::Get_mass_pre(){
  return _mass_pre;
}
double ISR_preFSR::Get_mass_post(){
  return _mass_post;
}
double ISR_preFSR::Get_weight(){
  return _weight;
}


void ISR_preFSR::runMass(vector<TString> files,  TString tag){
  cout<<"##############################"<<endl;
  cout<<tag<<endl;


  Set_nbin_xmin_xmax_Gen(155,40,350);

  Set_nbin_xmin_xmax_Reco(155,40,350);

  SetName_hGenRecoMC(tag+" x_preFSR_y_postFSR");
  SetName_hGenMC(tag+" preFSR");
  SetName_hRecoMC(tag+" postFSR");

  Init_histos();
  cout<<"files.size()="<<files.size()<<endl;
  for(int i_f=0; i_f < files.size(); i_f++){ 
    TString filename=files[i_f];
    LoadPW(filename);
    
    Normalizer Norm;
    double xsec=1;
    if(filename.Contains("10to50"))  xsec=18610.;
    else if(filename.Contains("50plus"))  xsec=5760.;
    else if(filename.Contains("200to400"))  xsec=7.77;
    Norm.run(filename,"DY/truweight","weight",xsec);
    double normfactor=1/Norm.Get_efflumi();
    
    unsigned int Nevent=GetNevent();
    cout<<"Nevent="<<Nevent<<endl;
    for( int i = 0; i < Nevent; i++){
      /*      
if(i%10000==1) {

	cout<<Get_weight()<<endl;
	cout<<normfactor<<endl;


      }
      */
      GetEntrySetPtMass(i);
      FillhMC(Get_mass_pre(),Get_mass_post(), Get_weight()*normfactor );   
      
    }
  }
  cout<<"End of runMass"<<endl;  
  
}





void test_norm(){
  Normalizer Norm_mu10to50_Z;
  Norm_mu10to50_Z.run("/cms/ldap_home/jhchoi/PHOTOS_Study/for2015dataset/PHOTOS_Analyzer/combine_files_muon/DYPHOTOS_add_Zgen10to50.root","DY/truweight","weight",18610);
  //void Normalizer::run(TString filename, TString treename,TString wname, double xsec, double lumi=1){
  cout<<Norm_mu10to50_Z.Get_efflumi()<<endl;

}


void test_JUnfold(){
  vector<TString> vlep;
  vlep.push_back("muon");
  vlep.push_back("electron");

  vector<TString> vboson;
  vboson.push_back("Z");
  vboson.push_back("gamma");
  for(int i=0;i<vlep.size();i++){

    for(int j=0;j<vboson.size();j++){
      vector<TString> files; 
      files.push_back("/cms/ldap_home/jhchoi/PHOTOS_Study/for2015dataset/PHOTOS_Analyzer/combine_files_"+vlep[i]+"/DYPHOTOS_add_"+vboson[j]+"gen10to50.root");
      files.push_back("/cms/ldap_home/jhchoi/PHOTOS_Study/for2015dataset/PHOTOS_Analyzer/combine_files_"+vlep[i]+"/DYPHOTOS_add_"+vboson[j]+"gen50plus.root");
      files.push_back("/cms/ldap_home/jhchoi/PHOTOS_Study/for2015dataset/PHOTOS_Analyzer/combine_files_"+vlep[i]+"/DYPHOTOS_add_"+vboson[j]+"gen200to400.root");
      
      ISR_preFSR ana;
      ana.runMass(files,"PHOTOS_"+vlep[i]+"_add"+vboson[j]);
      ana.ImportHists("PHOTOS_"+vlep[i]+"_add"+vboson[j]);
      
      files.clear();
      
    }
  }

}
void MyUnfold(){
  //  test_norm();
  test_JUnfold();
}
