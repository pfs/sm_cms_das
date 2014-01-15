//================================================================================================
//
// Perform fit to extract W->enu signal
//
//  * outputs plots and fit results summary
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TStyle.h>                       // class to handle ROOT plotting styles
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <TVectorD.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include <fstream>                        // functions for file I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include "Math/LorentzVector.h"           // 4-vector class

#include "Utils/MyTools.hh"	          // various helper functions
#include "Utils/CPlot.hh"	          // helper class for plots
#include "Utils/MitStyleRemix.hh"         // style settings for drawing
#include "Utils/WModels.hh"               // definitions of PDFs for fitting
#include "Utils/RecoilCorrector.hh"       // class to handle recoil corrections for MET

// RooFit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#endif

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

//=== FUNCTION DECLARATIONS ======================================================================================

//Generate Toys of a fit model
TTree* toyHist(int iNToys, bool iDoMu,std::string iName,TH1D *iData,TH1D *iW,TH1D *iEWK,TH1D *iAntiData,TH1D *iAntiW,TH1D *iAntiEWK,double *iResults);

//Fit histogram with the default function and plot it
double* fitHist(TCanvas *iC, bool iDoMu,int iPlot, std::string iName,TH1D *iData,TH1D *iW,TH1D *iEWK,TH1D *iAntiData,TH1D *iAntiW,TH1D *iAntiEWK,const Double_t METMAX,const Int_t NBINS,const Double_t lumi,const Int_t Ecm,int iAltQCD=-1);

// make data-fit difference plots
TH1D* makeDiffHist(TH1D* hData, TH1D* hFit, const TString name);

// print correlations of fitted parameters
void printCorrelations(ostream& os, RooFitResult *res);

// print chi2 test and KS test results
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe);

// make webpage
void makeHTML(const TString outDir);

// Global Hists for generating best fits from fit to alternative toy
TH1D *fBestFit     = 0;
TH1D *fAntiBestFit = 0;

//=== MAIN MACRO ================================================================================================= 
void fitWe(const TString  outputDir="test",  // output directory
           const Double_t lumi=18.7,         // integrated luminosity (/fb)
	   const Int_t    Ecm=8,
	   const Bool_t   doMu=true,         // center-of-mass energy
	   const Int_t    plotOption=1,      // Set Options below
	   const Int_t    fitOption=-1
) {
  gBenchmark->Start("fitWe");
  //plotOptions
  //plotOption == 0 fit and return results
  //plotOption == 1 plot and return results
  //plotOption == 2 plot, dump text file, and return results
  //plotOption == 3 plot, set fBestFit and fAntiBestFit to current best fits, generate toys of it 

  //fitOption
  //-1 Use the default fit
  // 0 Use the electron fit
  // 1 Use the muon fit
  // 2 Use the electron stystematics fit 

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  // MET histogram binning and range
  const Int_t    NBINS  = 50;
  const Double_t METMAX = 100;

  const Double_t PT_CUT  = 25;
  Double_t ETA_CUT = 2.4;
  if(doMu) ETA_CUT = 2.1;

  // file name with recoil correction
  TString recoilfname("Recoil/ZeeData/fits.root");
  if(doMu) recoilfname = "Recoil/ZmmData/fits.root";
  
  //
  // input ntuple file names
  //
  enum { eData, eW, eEWK, eAntiData, eAntiW, eAntiEWK  };  // data type enum
  vector<TString> fnamev;
  vector<Int_t>   typev;
  
  if(!doMu) { fnamev.push_back("SingleElectron.root"); typev.push_back(eData);}
  if(!doMu) { fnamev.push_back("WminusToENu.root");    typev.push_back(eW);}
  if(!doMu) { fnamev.push_back("WplusToENu.root");     typev.push_back(eW);}
  if(!doMu) { fnamev.push_back("DYToEE.root");         typev.push_back(eEWK);}
  if(doMu)  { fnamev.push_back("SingleMu.root");       typev.push_back(eData);}
  if(doMu)  { fnamev.push_back("WminusToMuNu.root");   typev.push_back(eW);}
  if(doMu)  { fnamev.push_back("WplusToMuNu.root");    typev.push_back(eW);}
  if(doMu)  { fnamev.push_back("DYToMuMu.root");       typev.push_back(eEWK);}
  fnamev.push_back("WToTauNu.root");                   typev.push_back(eEWK);
  fnamev.push_back("TTJets.root");                     typev.push_back(eEWK);
  fnamev.push_back("WW.root");                         typev.push_back(eEWK);
  fnamev.push_back("WZ.root");                         typev.push_back(eEWK);
  fnamev.push_back("ZZ.root");                         typev.push_back(eEWK);

  if(!doMu) { fnamev.push_back("SingleElectron.root"); typev.push_back(eAntiData);}
  if(!doMu) { fnamev.push_back("WminusToENu.root");    typev.push_back(eAntiW); }
  if(!doMu) { fnamev.push_back("WplusToENu.root");     typev.push_back(eAntiW); }
  if(!doMu) { fnamev.push_back("DYToEE.root");         typev.push_back(eAntiEWK);   } 
  if(doMu)  { fnamev.push_back("SingleMu.root");       typev.push_back(eAntiData);}
  if(doMu)  { fnamev.push_back("WminusToMuNu.root");   typev.push_back(eAntiW); }
  if(doMu)  { fnamev.push_back("WplusToMuNu.root");    typev.push_back(eAntiW); }
  if(doMu)  { fnamev.push_back("DYToMuMu.root");       typev.push_back(eAntiEWK);}
  fnamev.push_back("WToTauNu.root");       typev.push_back(eAntiEWK); 
  fnamev.push_back("TTJets.root");         typev.push_back(eAntiEWK); 
  fnamev.push_back("WW.root");             typev.push_back(eAntiEWK); 
  fnamev.push_back("WZ.root");             typev.push_back(eAntiEWK); 
  fnamev.push_back("ZZ.root");             typev.push_back(eAntiEWK); 


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir;  
  
   // Access recoil corrections
  //RecoilCorrector recoilCorr(recoilfname);
  
  //
  // Declare MET histograms
  //
  TH1D *hDataMet  = new TH1D("hDataMet", "",NBINS,0,METMAX); hDataMet->Sumw2();
  TH1D *hDataMetm = new TH1D("hDataMetm","",NBINS,0,METMAX); hDataMetm->Sumw2();  
  TH1D *hDataMetp = new TH1D("hDataMetp","",NBINS,0,METMAX); hDataMetp->Sumw2();
  TH1D *hWMet     = new TH1D("hWMet",    "",NBINS,0,METMAX); hWMet    ->Sumw2();
  TH1D *hWMetp    = new TH1D("hWMetp",   "",NBINS,0,METMAX); hWMetp   ->Sumw2();
  TH1D *hWMetm    = new TH1D("hWMetm",   "",NBINS,0,METMAX); hWMetm   ->Sumw2();
  TH1D *hEWKMet   = new TH1D("hEWKMet",  "",NBINS,0,METMAX); hEWKMet  ->Sumw2();
  TH1D *hEWKMetp  = new TH1D("hEWKMetp", "",NBINS,0,METMAX); hEWKMetp ->Sumw2();
  TH1D *hEWKMetm  = new TH1D("hEWKMetm", "",NBINS,0,METMAX); hEWKMetm ->Sumw2();

  TH1D *hAntiDataMet   = new TH1D("hAntiDataMet", "",NBINS,0,METMAX); hAntiDataMet->Sumw2();
  TH1D *hAntiDataMetm  = new TH1D("hAntiDataMetm","",NBINS,0,METMAX); hAntiDataMetm->Sumw2();  
  TH1D *hAntiDataMetp  = new TH1D("hAntiDataMetp","",NBINS,0,METMAX); hAntiDataMetp->Sumw2();
  TH1D *hAntiWMet      = new TH1D("hAntiWMet",    "",NBINS,0,METMAX); hAntiWMet    ->Sumw2();
  TH1D *hAntiWMetp     = new TH1D("hAntiWMetp",   "",NBINS,0,METMAX); hAntiWMetp   ->Sumw2();
  TH1D *hAntiWMetm     = new TH1D("hAntiWMetm",   "",NBINS,0,METMAX); hAntiWMetm   ->Sumw2();
  TH1D *hAntiEWKMet    = new TH1D("hAntiEWKMet",  "",NBINS,0,METMAX); hAntiEWKMet  ->Sumw2();
  TH1D *hAntiEWKMetp   = new TH1D("hAntiEWKMetp", "",NBINS,0,METMAX); hAntiEWKMetp ->Sumw2();
  TH1D *hAntiEWKMetm   = new TH1D("hAntiEWKMetm", "",NBINS,0,METMAX); hAntiEWKMetm ->Sumw2();

  //
  // Declare variables to read in ntuple
  //
  Float_t npv, cat;
  Float_t genWPt;//, genWPhi;
  Float_t weight,puweight,pt,eta;
  Float_t met, metPhi, sumEt, mt;
  
  Float_t ngenerated,xsec,total,npass;
  
  TFile *infile=0;
  TTree *intree=0;

  //
  // Loop over files
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "...";
    infile = new TFile(fnamev[ifile]);	  assert(infile);
    intree = (TTree*)infile->Get("data/data"); assert(intree);
    ngenerated = 0;
    for(int i0 = 0; i0 < 10; i0++) { 
      std::stringstream pSS; pSS << "iniEvents;" << i0;
      TVectorD* pNGen = (TVectorD*) infile->Get(pSS.str().c_str());
      if(pNGen == 0) continue;
      ngenerated   += (*pNGen)[0];
    }
    TVectorD* pXSec = (TVectorD*) infile->Get("crossSection");
    npass        = 0;
    xsec         = (*pXSec)[0];
    weight       = xsec*lumi/ngenerated;  
    if(typev[ifile] == eData || typev[ifile] == eAntiData) weight = 1.;
    
    intree->SetBranchAddress("nvtx",     &npv);       // number of primary vertices
    intree->SetBranchAddress("weight",   &puweight);       // number of in-time PU events (MC)
    intree->SetBranchAddress("genv_pt",  &genWPt);    // GEN W boson pT (signal MC)
    intree->SetBranchAddress("leg2_pt",  &met);       // MET
    intree->SetBranchAddress("leg2_phi", &metPhi);    // phi(MET)
    intree->SetBranchAddress("sumEt",    &sumEt);     // Sum ET
    intree->SetBranchAddress("v_mt",     &mt);        // transverse mass
    intree->SetBranchAddress("cat",      &cat);        // transverse mass
    intree->SetBranchAddress("leg1_pt",  &pt);
    intree->SetBranchAddress("leg1_eta", &eta);
    //intree->SetBranchAddress("q",        &q);         // lepton charge
    //intree->SetBranchAddress("lep",      &lep);       // lepton 4-vector
    //intree->SetBranchAddress("sc",       &sc);        // electron Supercluster 4-vector
  
    //
    // loop over events
    //
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      if(!doMu && (fabs(cat) != 11 && fabs(cat) != 1100))           continue; 
      if( doMu && (fabs(cat) != 13 && fabs(cat) != 1300))           continue; 
      if(pt        < PT_CUT)  continue;	
      if(fabs(eta) > ETA_CUT) continue;
      if( (typev[ifile]==eData     || typev[ifile]==eW     || typev[ifile]==eEWK    ) && fabs(cat) > 1000)  continue;
      if( (typev[ifile]==eAntiData || typev[ifile]==eAntiW || typev[ifile]==eAntiEWK) && fabs(cat) < 1000)  continue;
      puweight *= weight;
      total += puweight;
      npass ++;
      if(typev[ifile]==eData) {
        hDataMet->Fill(met);
	if(cat>0) { hDataMetp->Fill(met); } 
	else    { hDataMetm->Fill(met); }
      }
      if(typev[ifile]==eAntiData) {
        hAntiDataMet->Fill(met);
	if(cat>0) { hAntiDataMetp->Fill(met); } 
	else    { hAntiDataMetm->Fill(met); }
      }
      if(typev[ifile]==eW) {
	Double_t corrMet=met;//, corrMetPhi=metPhi;
	// apply recoil corrections to W MC
	//recoilCorr.Correct(corrMet,corrMetPhi,genWPt,genWPhi,lep->Pt(),lep->Phi());
	hWMet->Fill(corrMet,weight);
	if(cat>0) { hWMetp->Fill(corrMet,puweight); } 
	else    { hWMetm->Fill(corrMet,puweight); }
      }
      if(typev[ifile]==eAntiW) {
	Double_t corrMet=met;//, corrMetPhi=metPhi;
	// apply recoil corrections to W MC
	//recoilCorr.Correct(corrMet,corrMetPhi,genWPt,genWPhi,lep->Pt(),lep->Phi());
	hAntiWMet->Fill(corrMet,weight);
	if(cat>0) { hAntiWMetp->Fill(corrMet,puweight); } 
	else      { hAntiWMetm->Fill(corrMet,puweight); }
      }
      if(typev[ifile]==eEWK) {
	hEWKMet->Fill(met,weight);
	if(cat>0) { hEWKMetp->Fill(met,puweight); }
	else      { hEWKMetm->Fill(met,puweight); }
      }
      if(typev[ifile]==eAntiEWK) {
	hAntiEWKMet->Fill(met,weight);
	if(cat>0) { hAntiEWKMetp->Fill(met,puweight); }
	else      { hAntiEWKMetm->Fill(met,puweight); }
      }
    }
    cout << " total/pb : " << total << endl;
    total = 0;
  }  
  delete infile;
  infile=0, intree=0;   

  TCanvas *c = MakeCanvas("c","c",800,800);
  c->Divide(1,2,0,0);
  c->cd(1)->SetPad(0,0.3,1.0,1.0);
  c->cd(1)->SetTopMargin(0.1);
  c->cd(1)->SetBottomMargin(0.01);
  c->cd(1)->SetLeftMargin(0.18);  
  c->cd(1)->SetRightMargin(0.07);  
  c->cd(1)->SetTickx(1);
  c->cd(1)->SetTicky(1);  
  c->cd(2)->SetPad(0,0,1.0,0.3);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetBottomMargin(0.45);
  c->cd(2)->SetLeftMargin(0.18);
  c->cd(2)->SetRightMargin(0.07);
  c->cd(2)->SetTickx(1);
  c->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.400,"Y");
  
  double *lResults  = fitHist(c,doMu,plotOption,"W" ,hDataMet, hWMet, hEWKMet,  hAntiDataMet, hAntiWMet, hAntiEWKMet,METMAX,NBINS,lumi,Ecm,fitOption);
  if(plotOption == 3) {hDataMet = fBestFit; hAntiDataMet = fAntiBestFit;}
  double *lResultsP = fitHist(c,doMu,plotOption,"WP",hDataMetp,hWMetp,hEWKMetp,hAntiDataMetp,hAntiWMetp,hAntiEWKMetp,METMAX,NBINS,lumi,Ecm,fitOption);
  if(plotOption == 3) {hDataMetp = fBestFit; hAntiDataMetp = fAntiBestFit;}
  double *lResultsM = fitHist(c,doMu,plotOption,"WM",hDataMetm,hWMetm,hEWKMetm,hAntiDataMetm,hAntiWMetm,hAntiEWKMetm,METMAX,NBINS,lumi,Ecm,fitOption);
  if(plotOption == 3) {hDataMetm = fBestFit; hAntiDataMetm = fAntiBestFit;}
  //Do Toys
  //Alterative Fit Model
  //Lepton EScale
  //Met Scale
  //EWK?
  if(plotOption != 3) return;
  TFile *file = new TFile("toys.root","RECREATE"); 
  TTree *biasResults  = toyHist(100,doMu,"W"  ,hDataMet,  hWMet,  hEWKMet,   hAntiDataMet,  hAntiWMet,  hAntiEWKMet ,lResults);
  TTree *biasResultsp = toyHist(100,doMu,"WP" ,hDataMetp, hWMetp, hEWKMetp,  hAntiDataMetp, hAntiWMetp, hAntiEWKMetp,lResultsP);
  TTree *biasResultsm = toyHist(100,doMu,"WM" ,hDataMetm, hWMetm, hEWKMetm,  hAntiDataMetm, hAntiWMetm, hAntiEWKMetm,lResultsM);

  biasResults ->Write();
  biasResultsp->Write();
  biasResultsm->Write();

  makeHTML(outputDir);
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;     
  //gBenchmark->Show("fitWe");
}
//=== FUNCTION DEFINITIONS ======================================================================================
TTree* toyHist(int iNToys, bool iDoMu,std::string iName,TH1D *iData,TH1D *iW,TH1D *iEWK,TH1D *iAntiData,TH1D *iAntiW,TH1D *iAntiEWK,double *iOriginalResults) { 
  TCanvas *c = 0;
  TTree *tree = new TTree(("Toys"+iName).c_str(),("Toys"+iName).c_str());
  float *vals = new float[16];
  for(int i0 = 0; i0 < 16; i0++) vals[i0] = 0.;
  tree->Branch("nsig"      ,&vals[0] ,"Val0/F");
  tree->Branch("newk"      ,&vals[1] ,"Val1/F");
  tree->Branch("nqcd"      ,&vals[2] ,"Val2/F");
  tree->Branch("nantisig"  ,&vals[3] ,"Val3/F");
  tree->Branch("nantiewk"  ,&vals[4] ,"Val4/F");
  tree->Branch("nantiqcd"  ,&vals[5] ,"Val5/F");
  tree->Branch("sigma"     ,&vals[6] ,"Val6/F");
  tree->Branch("a1"        ,&vals[7] ,"Val7/F");
  tree->Branch("nsig_e"    ,&vals[8] ,"Val8/F");
  tree->Branch("newk_e"    ,&vals[9] ,"Val9/F");
  tree->Branch("nqcd_e"    ,&vals[10],"Val10/F");
  tree->Branch("nantisig_e",&vals[11],"Val11/F");
  tree->Branch("nantiewk_e",&vals[12],"Val12/F");
  tree->Branch("nantiqcd_e",&vals[13],"Val13/F");
  tree->Branch("sigma_e"   ,&vals[14],"Val14/F");
  tree->Branch("a1_e"      ,&vals[15],"Val15/F");
  //Set The first entry of the tree to be the default values
  for(int i1 = 0; i1 < 16; i1++) vals[i1] = float(iOriginalResults[i1]);
  tree->Fill();

  TH1D *pData       = 0;
  TH1D *pAntiData   = 0;
  for(int i0 = 0; i0 < iNToys; i0++) { 
    std::stringstream pSS; pSS << "Hist" << i0;
    pData       = (TH1D*) iData    ->Clone((pSS.str()+"A").c_str());
    pAntiData   = (TH1D*) iAntiData->Clone((pSS.str()+"B").c_str());
    for(int i1 = 0; i1 < pData    ->GetNbinsX()+1; i1++) pData    ->SetBinContent(i1,0);
    for(int i1 = 0; i1 < pAntiData->GetNbinsX()+1; i1++) pAntiData->SetBinContent(i1,0);
    pData    ->FillRandom(iData,    iData    ->Integral());
    pAntiData->FillRandom(iAntiData,iAntiData->Integral());
    cout <<"==> " << pData->Integral() << " -- " << iData->Integral() << endl;
    double* pVals = fitHist(c,iDoMu,0,pSS.str(),pData,iW,iEWK,pAntiData,iAntiW,iAntiEWK,100,100,100,100);
    for(int i1 = 0; i1 < 16; i1++) vals[i1] = float(pVals[i1]);
    tree->Fill();
    delete pData; delete pAntiData;
  }
  return tree;
}
double* fitHist(TCanvas *iC, bool iDoMu,int iPlot, std::string iName,TH1D *iData,TH1D *iW,TH1D *iEWK,TH1D *iAntiData,TH1D *iAntiW,TH1D *iAntiEWK,const Double_t METMAX,const Int_t NBINS,const Double_t lumi,const Int_t Ecm,int iAltQCD) { 
  //
  // Declare fit parameters for signal and background yields
  // Note: W signal and EWK+top PDFs are constrained to the ratio described in MC
  //
  RooRealVar nSig(("nSig"+iName).c_str(),("nSig"+iName).c_str(),0.7*(iData->Integral()),0,iData->Integral());
  RooRealVar nQCD(("nQCD"+iName).c_str(),("nQCD"+iName).c_str(),0.3*(iData->Integral()),0,iData->Integral());
  RooRealVar cewk(("cewk"+iName).c_str(),("cewk"+iName).c_str(),0.1,0,5) ;
  cewk.setVal(iEWK->Integral()/iW->Integral());
  cewk.setConstant(kTRUE);
  RooFormulaVar nEWK(("nEWK"+iName).c_str(),("nEWK"+iName).c_str(),("cewk"+iName+"*nSig"+iName).c_str(),RooArgList(nSig,cewk));

  RooRealVar nAntiSig(("nAntiSig"+iName).c_str(),("nAntiSig"+iName).c_str(),0.05*(iAntiData->Integral()),0,iAntiData->Integral());
  RooRealVar nAntiQCD(("nAntiQCD"+iName).c_str(),("nAntiQCD"+iName).c_str(),0.9 *(iAntiData->Integral()),0,iAntiData->Integral());
  RooRealVar dewk    (("dewk"    +iName).c_str(),("dewk"    +iName).c_str(),0.1,0,5) ;
  dewk.setVal(iAntiEWK->Integral()/iAntiW->Integral());
  dewk.setConstant(kTRUE);
  RooFormulaVar nAntiEWK(("nAntiEWK"+iName).c_str(),("nAntiEWK"+iName).c_str(),("dewk"+iName+"*nAntiSig"+iName).c_str(),RooArgList(nAntiSig,dewk));

  //
  // Construct PDFs for fitting
  //
  RooRealVar pfmet(("pfmet"+iName).c_str(),("pfmet"+iName).c_str(),0,METMAX);
  pfmet.setBins(NBINS);
   
  // Signal PDFs
  RooDataHist wMet  (("wMET" +iName).c_str(), ("wMET" +iName).c_str(),     RooArgSet(pfmet),iW);     RooHistPdf  pdfW(( "w"+iName).c_str(), ( "w"+iName).c_str(), pfmet, wMet, 1);
  RooDataHist awMet (("awMET"+iName).c_str(), ("awMET"+iName).c_str(),     RooArgSet(pfmet),iAntiW); RooHistPdf apdfW(("aw"+iName).c_str(), ("aw"+iName).c_str(), pfmet,awMet, 1);
  
  // EWK+top PDFs
  RooDataHist ewkMet (("ewkMET" +iName).c_str(),( "ewkMET"+iName).c_str(), RooArgSet(pfmet),iEWK);     RooHistPdf  pdfEWK (( "ewk"+iName).c_str(),( "ewk"+iName).c_str(), pfmet,ewkMet, 1);
  RooDataHist aewkMet(("aewkMET"+iName).c_str(),("aewkMET"+iName).c_str(), RooArgSet(pfmet),iAntiEWK); RooHistPdf apdfEWK (("aewk"+iName).c_str(),("aewk"+iName).c_str(), pfmet,aewkMet, 1);
  
  // QCD Pdfs
  CPepeModel0 qcd0 (("qcd0" +iName).c_str(),pfmet);
  CPepeModel1 qcd1 (("qcd1" +iName).c_str(),pfmet);
  CPepeModel2 qcd2 (("qcd2" +iName).c_str(),pfmet);
  CPepeModel0 aqcd0(("aqcd0"+iName).c_str(),pfmet);
  CPepeModel1 aqcd1(("aqcd1"+iName).c_str(),pfmet,qcd1.sigma);
  CPepeModel2 aqcd2(("aqcd2"+iName).c_str(),pfmet);
  RooGenericPdf *lQCD  =  qcd0.model;
  RooGenericPdf *lAQCD = aqcd0.model;
  if(iDoMu) lQCD  = qcd1.model; 
  if(iDoMu) lAQCD = aqcd1.model; 
  if(iAltQCD == 0) lQCD = qcd0.model;
  if(iAltQCD == 1) lQCD = qcd1.model;
  if(iAltQCD == 2) lQCD = qcd2.model;
  if(iAltQCD == 0) lAQCD = aqcd0.model;
  if(iAltQCD == 1) lAQCD = aqcd1.model;
  if(iAltQCD == 2) lAQCD = aqcd2.model;
  
  // Signal + Background PDFs
  RooAddPdf pdfMet (("pdfMet"+iName).c_str(), ("pdfMet" +iName).c_str(), RooArgList(pdfW ,pdfEWK ,*lQCD), RooArgList(nSig,nEWK,nQCD));  
  RooAddPdf apdfMet(("apdfMet"+iName).c_str(),("apdfMet"+iName).c_str(), RooArgList(apdfW,apdfEWK,*lAQCD), RooArgList(nAntiSig,nAntiEWK,nAntiQCD));  
  
  // PDF for simultaneous fit
  RooCategory rooCat("rooCat","rooCat");
  rooCat.defineType("Select");
  rooCat.defineType("Anti");
  
  RooSimultaneous pdfTotal("pdfTotal","pdfTotal",rooCat);
  pdfTotal.addPdf(pdfMet, "Select");
  if(iDoMu) pdfTotal.addPdf(apdfMet,"Anti");
      
  // Perform fits
  RooDataHist dataMet  (("dataMet"+iName).c_str(),("dataMet"+iName).c_str(), RooArgSet(pfmet),iData);
  RooDataHist antiMet  (("antiMet"+iName).c_str(),("antiMet"+iName).c_str(), RooArgSet(pfmet),iAntiData);
  RooDataHist dataTotal(("data"   +iName).c_str(),("data"   +iName).c_str(), RooArgList(pfmet), Index(rooCat),
                        Import("Select", dataMet),
                        Import("Anti",   antiMet));
  
  RooFitResult *fitRes = 0;
  bool runMinos = kTRUE; 
  if(iPlot == 0 || iPlot == 3) runMinos = kFALSE; //Remove Minos when running toys (too damn slow)
  if(!iDoMu) fitRes = pdfMet  .fitTo(dataMet  ,Extended(),Minos(runMinos),Save(kTRUE));
  if( iDoMu) fitRes = pdfTotal.fitTo(dataTotal,Extended(),Minos(runMinos),Save(kTRUE));

  double *lResults = new double[16];
  lResults[0]  = nSig.getVal(); 
  lResults[1]  = nEWK.getVal(); 
  lResults[2]  = nQCD.getVal(); 
  lResults[3]  = nAntiSig.getVal(); 
  lResults[4]  = nAntiEWK.getVal(); 
  lResults[5]  = nAntiQCD.getVal(); 
  if(!iDoMu) lResults[6]  = double(qcd0.sigma->getVal()); 
  if( iDoMu) lResults[6]  = double(qcd1.sigma->getVal()); 

  lResults[7]  = 0.;
  if(!iDoMu) lResults[7]  = qcd1.a1->getVal(); 
  lResults[8]   = nSig    .getPropagatedError(*fitRes);
  lResults[9]   = nEWK    .getPropagatedError(*fitRes);
  lResults[10]  = nQCD    .getPropagatedError(*fitRes);
  lResults[11]  = nAntiSig.getPropagatedError(*fitRes);
  lResults[12]  = nAntiEWK.getPropagatedError(*fitRes);
  lResults[13]  = nAntiQCD.getPropagatedError(*fitRes);
  if( iDoMu) lResults[14]  = qcd0.sigma->getError();
  if(!iDoMu) lResults[14]  = qcd1.sigma->getError();
  if( iDoMu) lResults[15]  = 0;
  if(!iDoMu) lResults[15]  = qcd1.a1   ->getError();
  if(iPlot == 0 ) return lResults;
  //
  // Use histogram version of fitted PDFs to make ratio plots
  // (Will also use PDF histograms later for Chi^2 and KS tests)
  //
  TH1D *hPdfMet = (TH1D*)(pdfMet.createHistogram(("hPdfMet"+iName).c_str(), pfmet));
  hPdfMet->Scale((nSig.getVal()+nEWK.getVal()+nQCD.getVal())/hPdfMet->Integral());
  TH1D *hMetDiff = makeDiffHist(iData,hPdfMet,"hMetDiff"+iName);
  hMetDiff->SetMarkerStyle(kFullCircle);
  hMetDiff->SetMarkerSize(0.9);
  
  TH1D *hPdfAntiMet = (TH1D*)(apdfMet.createHistogram(("hPdfAntiMet"+iName).c_str(), pfmet));
  hPdfAntiMet->Scale((nAntiSig.getVal()+nAntiEWK.getVal()+nAntiQCD.getVal())/hPdfAntiMet->Integral());
  TH1D *hAntiMetDiff = makeDiffHist(iAntiData,hPdfAntiMet,"hAntiMetDiff"+iName);
  hAntiMetDiff->SetMarkerStyle(kFullCircle);
  hAntiMetDiff->SetMarkerSize(0.9);
  if(iPlot == 3 ) {
    //Build best fit QCD with default W and EWK Shap
    TH1D *hPdfMetQCD  = (TH1D*)(lQCD ->createHistogram(("hPdfMetQCD"    +iName).c_str(), pfmet));
    TH1D *hPdfAMetQCD = (TH1D*)(lAQCD->createHistogram(("hPdfAntiMetQCD"+iName).c_str(), pfmet));
    hPdfMetQCD ->Scale(nQCD    .getVal()/hPdfMetQCD ->Integral());
    hPdfAMetQCD->Scale(nAntiQCD.getVal()/hPdfAMetQCD->Integral());

    TH1D *pW    = (TH1D*) iW      ->Clone("WForToys");    pW   ->Scale(nSig    .getVal()/pW   ->Integral());
    TH1D *pEWK  = (TH1D*) iEWK    ->Clone("EWKForToys");  pEWK ->Scale(nEWK    .getVal()/pEWK ->Integral());
    TH1D *pAW   = (TH1D*) iAntiW  ->Clone("AWForToys");   pAW  ->Scale(nAntiSig.getVal()/pAW  ->Integral());
    TH1D *pAEWK = (TH1D*) iAntiEWK->Clone("AEWKForToys"); pAEWK->Scale(nAntiEWK.getVal()/pAEWK->Integral());
    hPdfMetQCD ->Add(pW);
    hPdfMetQCD ->Add(pEWK);
    hPdfAMetQCD->Add(pAW);
    hPdfAMetQCD->Add(pAEWK);
    fBestFit     = hPdfMetQCD; 
    fAntiBestFit = hPdfAMetQCD; 
    return lResults;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  char ylabel[100];  // string buffer for y-axis label

  // file format for output plots
  const TString format("png"); 
  
  // label for lumi
  char lumitext[100];
  if(lumi<0.1) sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = %i TeV",lumi*1000.,Ecm);
  else         sprintf(lumitext,"%.2f fb^{-1}  at  #sqrt{s} = %i TeV",lumi      ,Ecm);
  
  // plot colors
  Int_t linecolorW   = kOrange-3;
  Int_t fillcolorW   = kOrange-2;
  Int_t linecolorEWK = kOrange+10;
  Int_t fillcolorEWK = kOrange+7;
  Int_t linecolorQCD = kViolet+2;
  Int_t fillcolorQCD = kViolet-5;
  Int_t ratioColor   = kGray+2;
  
  //
  // Dummy histograms for TLegend
  // (Nobody can figure out how to properly pass RooFit objects...)
  //
  TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  hDummyData->SetMarkerStyle(kFullCircle);
  hDummyData->SetMarkerSize(0.9);
  
  TH1D *hDummyW = new TH1D("hDummyW","",0,0,10);
  hDummyW->SetLineColor(linecolorW);
  hDummyW->SetFillColor(fillcolorW);
  hDummyW->SetFillStyle(1001);
  
  TH1D *hDummyEWK = new TH1D("hDummyEWK","",0,0,10);
  hDummyEWK->SetLineColor(linecolorEWK);
  hDummyEWK->SetFillColor(fillcolorEWK);
  hDummyEWK->SetFillStyle(1001);
  
  TH1D *hDummyQCD = new TH1D("hDummyQCD","",0,0,10);
  hDummyQCD->SetLineColor(linecolorQCD);
  hDummyQCD->SetFillColor(fillcolorQCD);
  hDummyQCD->SetFillStyle(1001);
   
  //
  // W MET plot
  //
  RooPlot *wmframe = pfmet.frame(Bins(NBINS));    
  dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
  pdfMet.plotOn(wmframe,FillColor(fillcolorW),DrawOption("F"));
  pdfMet.plotOn(wmframe,LineColor(linecolorW));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK,*lQCD)),FillColor(fillcolorEWK),DrawOption("F"));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfEWK,*lQCD)),LineColor(linecolorEWK));
  pdfMet.plotOn(wmframe,Components(RooArgSet(*lQCD)),FillColor(fillcolorQCD),DrawOption("F"));
  pdfMet.plotOn(wmframe,Components(RooArgSet(*lQCD)),LineColor(linecolorQCD));
  pdfMet.plotOn(wmframe,Components(RooArgSet(pdfW)),LineColor(linecolorW),LineStyle(2));
  dataMet.plotOn(wmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
  
  sprintf(ylabel,"Events / %.1f GeV",iData->GetBinWidth(1));
  CPlot plotMet(("fitmet"+iName).c_str(),wmframe,"","",ylabel);
  plotMet.SetLegend(0.68,0.57,0.93,0.77);
  plotMet.GetLegend()->AddEntry(hDummyData,"data","PL");
  plotMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
  plotMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
  plotMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
  plotMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
  plotMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
  plotMet.SetYRange(0.1,1.1*(iData->GetMaximum()));
  plotMet.Draw(iC,kFALSE,format,1);

  CPlot plotMetDiff(("fitmet"+iName).c_str(),"","#slash{E}_{T} [GeV]","#chi");
  plotMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
  plotMetDiff.SetYRange(-8,8);
  plotMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
  plotMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
  plotMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
  plotMetDiff.Draw(iC,kTRUE,format,2);
  plotMet.Draw(iC,kTRUE,format,1);

  plotMet.SetName(("fitmetlog"+iName).c_str());
  plotMet.SetLogy();
  plotMet.SetYRange(1e-3*(iData->GetMaximum()),10*(iData->GetMaximum()));
  plotMet.Draw(iC,kTRUE,format,1);
  
  if(iDoMu) {
    RooPlot *awmframe = pfmet.frame(Bins(NBINS));    
    antiMet.plotOn(awmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));
    apdfMet.plotOn(awmframe,FillColor(fillcolorW),DrawOption("F"));
    apdfMet.plotOn(awmframe,LineColor(linecolorW));
    apdfMet.plotOn(awmframe,Components(RooArgSet(apdfEWK,*lQCD)),FillColor(fillcolorEWK),DrawOption("F"));
    apdfMet.plotOn(awmframe,Components(RooArgSet(apdfEWK,*lQCD)),LineColor(linecolorEWK));
    apdfMet.plotOn(awmframe,Components(RooArgSet(*lQCD)),FillColor(fillcolorQCD),DrawOption("F"));
    apdfMet.plotOn(awmframe,Components(RooArgSet(*lQCD)),LineColor(linecolorQCD));
    apdfMet.plotOn(awmframe,Components(RooArgSet(apdfW)),LineColor(linecolorW),LineStyle(2));
    antiMet.plotOn(awmframe,MarkerStyle(kFullCircle),MarkerSize(0.9),DrawOption("ZP"));  
    
    sprintf(ylabel,"Events / %.1f GeV",iAntiData->GetBinWidth(1));
    CPlot plotAntiMet(("fitantimet"+iName).c_str(),awmframe,"","",ylabel);
    plotAntiMet.SetLegend(0.68,0.57,0.93,0.77);
    plotAntiMet.GetLegend()->AddEntry(hDummyData,"data","PL");
    plotAntiMet.GetLegend()->AddEntry(hDummyW,"W#rightarrow#mu#nu","F");
    plotAntiMet.GetLegend()->AddEntry(hDummyEWK,"EWK+t#bar{t}","F");
    plotAntiMet.GetLegend()->AddEntry(hDummyQCD,"QCD","F");
    plotAntiMet.AddTextBox(lumitext,0.55,0.80,0.90,0.86,0);
    plotAntiMet.AddTextBox("CMS Preliminary",0.63,0.92,0.95,0.99,0);
    plotAntiMet.SetYRange(0.1,1.1*(iAntiData->GetMaximum()));
    plotAntiMet.Draw(iC,kFALSE,format,1);
    
    CPlot plotAntiMetDiff(("fitantimet"+iName).c_str(),"","#slash{E}_{T} [GeV]","#chi");
    plotAntiMetDiff.AddHist1D(hMetDiff,"EX0",ratioColor);
    plotAntiMetDiff.SetYRange(-8,8);
    plotAntiMetDiff.AddLine(0, 0,METMAX, 0,kBlack,1);
    plotAntiMetDiff.AddLine(0, 5,METMAX, 5,kBlack,3);
    plotAntiMetDiff.AddLine(0,-5,METMAX,-5,kBlack,3);
    plotAntiMetDiff.Draw(iC,kTRUE,format,2);
    
    plotAntiMet.SetName(("fitantimetlog"+iName).c_str());
    plotAntiMet.SetLogy();
    plotAntiMet.SetYRange(1e-3*(iAntiData->GetMaximum()),10*(iAntiData->GetMaximum()));
    plotAntiMet.Draw(iC,kTRUE,format,1);
  }
  if(iPlot == 1) return lResults;
  
  ofstream txtfile;
  std::string txtfName = "fitres"+iName;
  if( iDoMu) txtfName + "Mu.txt";
  if(!iDoMu) txtfName + "Mu.txt";
  ios_base::fmtflags flags;
  cout << " --- test " << iData->Integral() << " -- " << hPdfMet->Integral() << endl;
  Double_t chi2prob = iData->Chi2Test(hPdfMet,"PUW");
  Double_t chi2ndf  = iData->Chi2Test(hPdfMet,"CHI2/NDFUW");
  Double_t ksprob   = iData->KolmogorovTest(hPdfMet);
  Double_t ksprobpe = 1;//iData->KolmogorovTest(hPdfMet,"DX");
  txtfile.open(txtfName.c_str());
  assert(txtfile.is_open());
  
  flags = txtfile.flags();
  txtfile << setprecision(10);
  txtfile << " *** Yields *** " << endl;
  txtfile << "Selected: " << iData->Integral() << endl;
  txtfile << "  Signal: " << nSig.getVal() << " +/- " << nSig.getPropagatedError(*fitRes) << endl;
  txtfile << "     QCD: " << nQCD.getVal() << " +/- " << nQCD.getPropagatedError(*fitRes) << endl;
  txtfile << "   Other: " << nEWK.getVal() << " +/- " << nEWK.getPropagatedError(*fitRes) << endl;
  txtfile << endl;
  txtfile.flags(flags);
  
  fitRes->printStream(txtfile,RooPrintable::kValue,RooPrintable::kVerbose);
  txtfile << endl;
  printCorrelations(txtfile, fitRes);
  txtfile << endl;
  printChi2AndKSResults(txtfile, chi2prob, chi2ndf, ksprob, ksprobpe);
  txtfile.close();
  return lResults;
}

//--------------------------------------------------------------------------------------------------
TH1D *makeDiffHist(TH1D* hData, TH1D* hFit, const TString name)
{
  TH1D *hDiff = new TH1D(name,"",hData->GetNbinsX(),hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax());
  for(Int_t ibin=1; ibin<=hData->GetNbinsX(); ibin++) {
    
    Double_t diff = (hData->GetBinContent(ibin)-hFit->GetBinContent(ibin));
    
    Double_t err = sqrt(hData->GetBinContent(ibin));
    if(err==0) err= sqrt(hFit->GetBinContent(ibin));
    
    if(err>0) hDiff->SetBinContent(ibin,diff/err);
    else      hDiff->SetBinContent(ibin,0);
    hDiff->SetBinError(ibin,1);   
  }
  
  hDiff->GetYaxis()->SetTitleOffset(0.48);
  hDiff->GetYaxis()->SetTitleSize(0.13);
  hDiff->GetYaxis()->SetLabelSize(0.10);
  hDiff->GetYaxis()->SetNdivisions(104);
  hDiff->GetYaxis()->CenterTitle();
  hDiff->GetXaxis()->SetTitleOffset(1.2);
  hDiff->GetXaxis()->SetTitleSize(0.13);
  hDiff->GetXaxis()->SetLabelSize(0.12);
  hDiff->GetXaxis()->CenterTitle();
  
  return hDiff;
}

//--------------------------------------------------------------------------------------------------
void printCorrelations(ostream& os, RooFitResult *res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList parlist = res->floatParsFinal();
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist.getSize(); i++) {
    for(Int_t j=0; j<parlist.getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

//--------------------------------------------------------------------------------------------------
void printChi2AndKSResults(ostream& os, 
                           const Double_t chi2prob, const Double_t chi2ndf, 
			   const Double_t ksprob, const Double_t ksprobpe)
{
  ios_base::fmtflags flags = os.flags();
  
  os << "  Chi2 Test" << endl;
  os << " -----------" << endl;
  os << "       prob = " << chi2prob << endl;
  os << "   chi2/ndf = " << chi2ndf << endl;
  os << endl;
  os << "  KS Test" << endl;
  os << " ---------" << endl;
  os << "   prob = " << ksprob << endl;
  os << "   prob = " << ksprobpe << " with 1000 pseudo-experiments" << endl;
  os << endl;
 
  os.flags(flags);
}
			   
//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/WenuFitPlots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Wenu</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmet.png\"><img src=\"fitmet.png\" alt=\"fitmet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetp.png\"><img src=\"fitmetp.png\" alt=\"fitmetp.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetm.png\"><img src=\"fitmetm.png\" alt=\"fitmetm.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetlog.png\"><img src=\"fitmetlog.png\" alt=\"fitmetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetplog.png\"><img src=\"fitmetplog.png\" alt=\"fitmetplog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"fitmetmlog.png\"><img src=\"fitmetmlog.png\" alt=\"fitmetmlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  
}
