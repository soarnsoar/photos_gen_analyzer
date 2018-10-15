class Normalizer{
  //10to50 18610
  //50plus 5765
  //200to400 7.77
private:
  TFile *_file;
  TTree *_tree;
  double _weight=1;
  double _weightsum=1;
  double _xsec=1;
  double _efflumi=1;
  //  double _norm=1;
public:
  double Get_weight();
  double Get_weightsum();
  double Get_xsec();
  double Get_efflumi();
  // double Get_norm();

  void Set_xsec(double xsec);
  //void Set_lumi(double lumi);
  

  void Load_weight(TString filename, TString treename,TString wname);//Get weight info
  void Cal_weightsum();
  //  void Cal_norm();
  void Cal_efflumi();

  void run(TString filename, TString treename,TString wname, double xsec);
};

void Normalizer::run(TString filename, TString treename,TString wname, double xsec){
  Load_weight(filename,treename,wname);
  Set_xsec(xsec);
  //Set_lumi(lumi);
  Cal_weightsum();
  //  Cal_norm();
  Cal_efflumi();
}

double Normalizer::Get_xsec(){
  return _xsec;
}
//double Normalizer::Get_lumi(){
// return _lumi;
//}
double Normalizer::Get_efflumi(){
  return _efflumi; 
}  

//double Normalizer::Get_norm(){
//return _norm;
//}

void Normalizer::Set_xsec(double xsec){
  _xsec=xsec;
}
//void Normalizer::Set_lumi(double lumi){
// _lumi=lumi;
//}

void Normalizer::Load_weight(TString filename, TString treename,TString wname){
  //  file=TFile::Open(targetdir+ProcessName+".root");
  _file=TFile::Open(filename);
  //lep1=(TTree*)file->Get("demo/lep1")
  _tree=(TTree*)_file->Get(treename);
  //lep1->SetBranchAddress("px",&lep1_px);
  _tree->SetBranchAddress(wname,&_weight);
}


void Normalizer::Cal_weightsum(){
  unsigned int Nevent=_tree->GetEntries();
  _weightsum=0;
  for(int i = 0 ; i < Nevent; i++){

    _tree->GetEntry(i);
    _weightsum+=_weight;


  }
}
//
//void Normalizer::Cal_norm(){
//Make the same lumi
//weightsum = xsec*lumi
//must be satisfied
//=> norm * weightsum = xsec*lumi
//_norm = _xsec*_lumi/_weightsum;

//}
void Normalizer::Cal_efflumi(){
  //lumi = weightsum/xsec
  
  _efflumi=_weightsum/_xsec;
}

double Normalizer::Get_weight(){
  return _weight;
}
double Normalizer::Get_weightsum(){
  return _weightsum;
}
