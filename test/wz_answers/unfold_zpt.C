{

  // ##################################################################################
  // # THESE DEFINITIONS NEED TO MATCH THOSE IN macroTemplate_AllSamples_unfolding.py #
  // ##################################################################################
  // Number of bins in final, unfolded distribution. Also used for histograms containing MC truth information.
  const int nbinsunfolded = 18;
  // Number of bins in reconstructed distribution. Using twice as many bins helps the stability of method.
  const int nbinsreconstructed = nbinsunfolded * 2;  
  // Bin edges for the generated and reconstructed histograms. Always one edge more than there are bins in the histogram.
  // Choosing these in such a way that all bins in the histograms contain the same numbers of events will increase the stability of the method. Depending on your use case, flattening the truth spectrum _after_ selection might be better than _before_ selection. 
  // But if you don't want to bother with this, using equidistant binning should work well, too.
  const double bin_edges_gen[] = {0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 30.0, 40.0, 50.0, 70.0, 90.0, 110.0, 150.0, 190.0, 250.0, 600.0};
  const double bin_edges_rec[] = {0.0, 1.5, 2.6, 3.5, 4.9, 5.9, 7.5, 8.5, 10.5, 11.25, 12.6, 13.9, 15.25, 16.4, 17.5, 18.75, 20.5, 25.8, 30.1, 36.0, 41.0, 46.0, 50.0, 61.0, 71.0, 82.0, 90.0, 101.0, 110.0, 132.0, 152.0, 173.0, 191.0, 226.0, 255.0, 500.0, 1000.0};
  // Alternative (VERY LARGE) set of bins
  // const double bin_edges_gen[] = {0, 26, 55, 104, 150};
  // const double bin_edges_rec[] = {0, 18, 28, 38, 52, 71, 101, 154, 200};
  // ##################################################################################
  // ##################################################################################

  TFile *fin = new TFile("plotter_unfolding.root");
  
  // Retrieve the migration/response matrix containing reconstructed values over MC truth values.
  TH2F* migrationmatrix = (TH2F*)fin->Get("DYToMuMu/DYToMuMu_vptmigrationmatrix");

  // Retrieve histogram with reconstructed events to be unfolded
  TH1F* data = (TH1F*)fin->Get("SingleMu/SingleMu_vpt");

  // TUnfoldSys provides methods to do systematic error propagation and to do unfolding with background subtraction
  // Regularization of curvature (TUnfold::kRegModeCurvature) is the least intrusive, but you can choose one of the other modes or give your own regularization conditions if this doesn't fit your needs.
  // Class Index: http://root.cern.ch/root/html/TUnfold.html
  TUnfoldSys tunfold(migrationmatrix, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeCurvature);
  
  // define the data histogram to be unfolded
  tunfold.SetInput(data);

  // Subtract your background templates (for WW, ZZ, WZ, TTJets, DYTauTau)
  // They need to use the reconstructed number of bins (bin_edges_rec) and you will need to subtract them from data before the unfolding
  TH1F* bg_all = (TH1F*)fin->Get("WW/WW_vpt");
  TH1F* bg_ZZ = (TH1F*)fin->Get("ZZ/ZZ_vpt");
  bg_all->Add(bg_ZZ);
  TH1F* bg_WZ = (TH1F*)fin->Get("WZ/WZ_vpt");
  bg_all->Add(bg_WZ);
  TH1F* bg_TTJets = (TH1F*)fin->Get("TTJets/TTJets_vpt");
  bg_all->Add(bg_TTJets);
  TH1F* bg_DYToTauTau = (TH1F*)fin->Get("DYToTauTau/DYToTauTau_vpt");
  bg_all->Add(bg_DYToTauTau);

  // Define the background histogram in TUnfold
  const double factor = 1; // additional normalization factor to be applied to template
  const double rel_error = 0.3; // uncertainty on normalization, in this example 30%
  tunfold.SubtractBackground(bg_all, "EWK_TT", factor, rel_error);
  
  
  TH1F* histogram_with_mc_truth = (TH1F*)fin->Get("DYToMuMu/DYToMuMu_genvpt");
  // Set a "bias" distribution, i.e. your MC truth distribution, corresponding to what you expect to see after unfolding.
  tunfold.SetBias(histogram_with_mc_truth);

  // Scale factor for the "bias" distribution. Choose such that normalization will correspond to expected normalization of result. 
  // For example, a simple approximation could be using expectedNUnfolded = n_data / ttbar_selection_efficiency
  // const double scaleBias = expectedNUnfolded / histogram_with_mc_truth->Integral();
  const double scaleBias = 1; // WE USE 1 IN THIS EXAMPLE

  // Regularization parameter, giving the strength of regularization. Will be roughly on the order of 1e-4.
  // You can determine this by performing unfolding with many different values, 
  // for example between 1e-3 and 1e-7 or such, and choosing that value of tau that minimizes tunfold.GetRhoAvg()
  // A different method would be to simply use tunfold.ScanLcurve(..), but this has proven far less reliable.
  
  // If the regularisation is strong, i.e. large parameter tau,
  // then the distribution x or its derivatives will look like the "bias" (mc truth) distribution. 
  // If the parameter tau is small, the distribution x is independent of the bias but will exhibit
  // many fluctuations and the high correlations between the bins.
  
  // Another way of crosschecking the tau parameter tuning is to look at the uncertainty of the bins after unfolding:
  // In absence of regularization they will be equal to the statistical (poisson) uncertainty
  // (which is the minimum error "possible", without considering correlations)
  // while introducing the regularization they will become a bit higher, due to the correlations,
  // but "not too much".

  const double tau = 1e-4; // value just an example, see above
    
  // Perform unfolding
  tunfold.DoUnfold(tau, data, scaleBias);
  // Get unfolded distrbution. 
  // WARNING: this will not use your own binning, so you might want to fill the contents into a correctly binned histogram.
  TH1D* unfolded = tunfold.GetOutput("unfolded", "unfolded", 0, 1);

  
  // THIS BOOL PARAMETER AND THE LOOP BELOW CAN BE USED TO TUNE THE TAU PARAMETER
  // BASED ON THE SUGGESTED tunfold.GetRhoAvg() method.
  // The loop compares also the final uncertainty of the unfolded bin number 5 (high statistics, no other real reason)
  // with its statistical uncertainty
  bool tune_tau_parameter = false;
    
  if(tune_tau_parameter){
    const int ntau = 16;
    const double tau_val[ntau] = {1e-1,5e-2,1e-2,7.5e-3,6.125e-3,5e-3,3.75e-3,2.5e-3,1e-3,1e-4,1e-5,1e-6,1e-7};
    double min_rho = 1e10;
    double i_min_rho = 0;
    TGraph *rho_vs_tau = new TGraph();
    for(int itau=0; itau<ntau; itau++){
      tunfold.DoUnfold(tau_val[itau], data, scaleBias);
      TH1D* unfolded = tunfold.GetOutput("unfolded", "unfolded", 0, 1);
      cout << "tau= " << tau_val[itau] << " tunfold.GetRhoAvg() " << tunfold.GetRhoAvg();
      cout << " bin5= " << unfolded->GetBinContent(5)<< " stat err. " << TMath::Sqrt(unfolded->GetBinContent(5)) << " unfolded err. " << unfolded->GetBinError(5) << endl;
      if(tunfold.GetRhoAvg()<min_rho){
        min_rho=tunfold.GetRhoAvg();
        i_min_rho=itau;
      }
      rho_vs_tau->SetPoint(itau,tau_val[itau],tunfold.GetRhoAvg());
    }
    cout << "best tau: " << tau_val[i_min_rho] << " min rho= " << min_rho << endl;
    rho_vs_tau->SaveAs("rho_vs_tau.root");
    return;
  }


  // Get covariance matrix. For description of error types, look at documentation of TUnfoldSys
  TH2D* cov_matrix = tunfold.GetEmatrix("cov_matrix", "cov_matrix"); // only error types d and e
  tunfold.GetEmatrixSysUncorr(cov_matrix, 0, false); // add error type a
  
  // Put the unfolded distribution in an histogram with the final binning (bin_edges_gen)
  TH1F* dataunfolded = (TH1F*)histogram_with_mc_truth->Clone("dataunfolded");
  dataunfolded->SetName("dataunfolded");
  dataunfolded->SetTitle("dataunfolded");
  for(int i=1; i<dataunfolded->GetNbinsX()+1; i++){
    dataunfolded->SetBinContent(i,unfolded->GetBinContent(i));
  }
  
  // save the histograms into a file
  TFile *fout = new TFile("vpt_unfolded.root","RECREATE");

  // subtract background from data before to save
  data->Add(bg_all,-1);
  // normalize bin to their bin width (in data)
  for(int i=1; i<data->GetNbinsX()+1; i++){
    	data->SetBinContent(i,data->GetBinContent(i)/data->GetBinWidth(i));
  }
  // then normalize all to 1
  data->Scale(1/data->Integral());
  data->Write();
  
  // normalize bin to their bin width (in reconstructed signal mc)
  for(int i=1; i<dataunfolded->GetNbinsX()+1; i++){
    	dataunfolded->SetBinContent(i,dataunfolded->GetBinContent(i)/dataunfolded->GetBinWidth(i));
  }
  // then normalize all to 1
  dataunfolded->Scale(1/dataunfolded->Integral());
  dataunfolded->Write();
  
  // normalize bin to their bin width (in reconstructed signal mc)
  for(int i=1; i<histogram_with_mc_truth->GetNbinsX()+1; i++){
    	histogram_with_mc_truth->SetBinContent(i,histogram_with_mc_truth->GetBinContent(i)/histogram_with_mc_truth->GetBinWidth(i));
  }
  // then normalize all to 1
  histogram_with_mc_truth->Scale(1/histogram_with_mc_truth->Integral());
  histogram_with_mc_truth->Write();
  
  // save the response matrix
  migrationmatrix->Write();
  // save the covariance matrix
  cov_matrix->Write();
  
  fout->Write();
  

}