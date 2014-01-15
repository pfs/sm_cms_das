void plotPull(std::string iTree="ToysW",std::string iPar="newk") { 
  TFile *lFile = new TFile("toys.root");
  TTree *lTree = (TTree*) lFile->FindObjectAny(iTree.c_str());
  float lPar  = 0; //Parameter
  float lEPar = 0; //Parameter Uncertainty 
  lTree->SetBranchAddress(iPar       .c_str(),&lPar);
  lTree->SetBranchAddress((iPar+"_e").c_str(),&lEPar);
  //Get the Mean
  lTree->GetEntry(0); 
  float lMean = lPar;
  float lErr  = lEPar;
  TH1F* lH    = new TH1F("Par"    ,"Par"    ,40, -lErr*7.,lErr*7.);
  TH1F* lPull = new TH1F("ParPull","ParPull",100,-10.              , 10.);
  TH1F* lRel  = new TH1F("ParRel" ,"ParRell",40, -0.1.            , 0.1);
  for(int i0 = 1; i0 < lTree->GetEntries(); i0++) {
    lTree->GetEntry(i0);
    lH   ->Fill(lPar-lMean);
    lRel ->Fill((lPar-lMean)/lMean);
    lPull->Fill((lPar-lMean)/lEPar);
  }
  cout << "Relative Difference : " << lRel->GetMean() << endl;
  lPull->Draw();
  lPull->Fit("gaus");
  
}
