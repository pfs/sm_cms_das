{    
  TString rfitpath("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.34.02-cms/include");
  TString path = gSystem->GetIncludePath();
  path += "-I. -I$ROOTSYS/src -I";
  path += rfitpath;
  gSystem->SetIncludePath(path.Data());
  
  gROOT->Macro("Utils/RooVoigtianShape.cc+");
  gROOT->Macro("Utils/RooCMSShape.cc+");
  
  gROOT->Macro("Utils/CPlot.cc++");
  gROOT->Macro("Utils/MitStyleRemix.cc++");  
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
