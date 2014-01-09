#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TSystem.h"
#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  int numberPDF(int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}

struct stPDFval{
  stPDFval() {}
  stPDFval(const stPDFval& arg) : 
    qscale(arg.qscale),
    x1(arg.x1),
    x2(arg.x2),
    id1(arg.id1),
    id2(arg.id2){
    wgts.clear();
  }
  std::vector<Float_t> wgts;
  Float_t qscale;
  Float_t x1;
  Float_t x2;
  int id1;
  int id2;
};

void printHelp()
{
  cout <<"-h   --help    --> print this text and exit" << endl
       <<"-i   --input   --> input root file" << endl
       <<"-o   --output  --> directory to store the output" << endl
       <<"-p   --pdf     --> CSV list of PDFs to compute" << endl
       <<"-t   --tree    --> tree name with the gen level information" <<endl;
}


//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  //check arguments
  TString inUrl(""),outUrl(""),treeName("");
  std::vector<TString> pdfSets;
  std::vector<Int_t>   nPdfVars;
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg=="-h" || arg=="--help")                 { printHelp();        return 0; }
    if((arg=="-i" || arg=="--input") && i+1<argc)  { inUrl=argv[i+1];    i++; }
    if((arg=="-o" || arg=="--output") && i+1<argc) { outUrl=argv[i+1];   i++; }
    if((arg=="-t" || arg=="--tree") && i+1<argc)   { treeName=argv[i+1]; i++; }
    if((arg=="-p" || arg=="--pdf") && i+1<argc)    { 
      TString pdfList(argv[i+1]);
      TObjArray *tkns=pdfList.Tokenize(",");
      for(Int_t itkn=0; itkn<tkns->GetEntriesFast(); itkn++) pdfSets.push_back( tkns->At(itkn)->GetName() );
      i++;
    }
  }
  if(inUrl=="" || pdfSets.size()==0) { printHelp(); return 0; }
  
  //
  // READ ALL EVENTS TO MEMORY
  //
  //open the INPUT file and get events tree
  gSystem->ExpandPathName(inUrl);
  TFile *file = TFile::Open(inUrl);
  if(file==0) return -1;
  if(file->IsZombie()) return -1;
  TTree *inTree=(TTree *) file->Get(treeName);
  if(inTree==0){
    file->Close();
    return -1;
  }
  stPDFval valForPDF;
  inTree->SetBranchAddress("qscale",      &valForPDF.qscale);
  inTree->SetBranchAddress("x1",          &valForPDF.x1);
  inTree->SetBranchAddress("x2",          &valForPDF.x2);
  inTree->SetBranchAddress("id1",         &valForPDF.id1);
  inTree->SetBranchAddress("id2",         &valForPDF.id2);
  std::vector<stPDFval> pdfVals;
  cout << "Progres bars        :0%%       20%%       40%%       60%%       80%%       100%%" << endl
       << "Scanning the ntuple :" << flush;
  int treeStep = inTree->GetEntriesFast()/50;
  if(treeStep==0)treeStep=1;
  for( int iev=0; iev<inTree->GetEntriesFast(); iev++){
    if(iev%treeStep==0){cout <<"."<<flush;}
    inTree->GetEntry(iev);
    pdfVals.push_back(valForPDF);
  }
  cout<< endl << endl;
  file->Close();

  //
  //NOW DUMP THE PDF WEIGHTS TO ANOTHER TREE
  //
  gSystem->ExpandPathName(outUrl);
  gSystem->Exec("mkdir -p " + outUrl);
  TString outFile(gSystem->BaseName(inUrl));
  outFile.ReplaceAll(".root","_pdf.root");
  outUrl = outUrl + "/" + outFile;
  printf("PDF weights will be saved in %s\n", outUrl.Data());
  TFile *ofile=TFile::Open(outUrl, "recreate");  

  //init PDF tool
  for(size_t ipdf=0; ipdf<pdfSets.size(); ipdf++)  
    {
      LHAPDF::initPDFSet(ipdf+1, pdfSets[ipdf].Data());
      nPdfVars.push_back( LHAPDF::numberPDF(ipdf+1) );
    }
  cout << "Loop on PDF sets and variations" << endl;
  for(size_t ipdf=0; ipdf<pdfSets.size(); ipdf++){

    //compute PDF weights (notice we divide by x because xfx returnx x\times PDF)
    cout << "   - starting with " << pdfSets[ipdf] << " which has " << nPdfVars[ipdf] << " variations..." << flush;
    for(int i=0; i <(nPdfVars[ipdf]+1); ++i){
      for(unsigned int v=0; v<pdfVals.size(); v++){ 
	LHAPDF::usePDFMember(ipdf+1,i);
	Float_t xpdf1 = LHAPDF::xfx(ipdf+1, pdfVals[v].x1, pdfVals[v].qscale, pdfVals[v].id1)/pdfVals[v].x1;
	Float_t xpdf2 = LHAPDF::xfx(ipdf+1, pdfVals[v].x2, pdfVals[v].qscale, pdfVals[v].id2)/pdfVals[v].x2;
	pdfVals[i].wgts.push_back( xpdf1 * xpdf2 );
      }
    }

    //create the output tree (simple array with n entries)
    ofile->cd();    
    char nameBuf[256];sprintf(nameBuf,"%s", pdfSets[ipdf].Data());
    cout << "...done - filling " << nameBuf << " tree..." << flush;
    TTree *pdfT = new TTree(nameBuf,"pdf");
    Int_t nVars(nPdfVars[ipdf]+1);
    Float_t pdfWgt[nVars];
    pdfT->Branch("n", &nVars, "n/I");
    pdfT->Branch("w", pdfWgt, "w[n]/F");
    pdfT->SetDirectory(ofile);
    for(unsigned int v=0; v<pdfVals.size(); v++){
      for(size_t iwgt=0; iwgt<pdfVals[v].wgts.size(); iwgt++) pdfWgt[iwgt]=pdfVals[v].wgts[iwgt];
      pdfT->Fill();
      pdfVals[v].wgts.clear();
    }
    cout << "...done" << endl;
    
    //write to file
    ofile->cd();
    pdfT->Write();
  }
  
  ofile->Close();
  cout << "All done here" << endl;
}  


