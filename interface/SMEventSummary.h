#ifndef smeventsummary_h
#define smeventsummary_h

#include "TTree.h"
#include "TMath.h"
#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#define MAXDATAOBJECTS 1000

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class SMEventSummary
{
 public:

  Int_t run, lumi, event, cat;
  
  //trigger information
  Int_t tn;
  Bool_t t_bits[MAXDATAOBJECTS];
  Int_t t_prescale[MAXDATAOBJECTS];

  //pileup information
  Int_t nvtx;
  Float_t instLumi;
  Float_t rho;

  //gen information
  Int_t ngenITpu, ngenOOTpu, ngenOOTpum1, ngenTruepu;
  Float_t pthat, genWeight, qscale, x1,x2;
  Int_t id1, id2, nup;
  Int_t mcn, mc_id[MAXDATAOBJECTS], mc_status[MAXDATAOBJECTS];
  Float_t mc_px[MAXDATAOBJECTS],mc_py[MAXDATAOBJECTS],mc_pz[MAXDATAOBJECTS],mc_en[MAXDATAOBJECTS], mc_lxy[MAXDATAOBJECTS]; 

  //super clusters
  Int_t   scn, scn_et[MAXDATAOBJECTS], scn_eta[MAXDATAOBJECTS], scn_phi[MAXDATAOBJECTS], scn_sihih[MAXDATAOBJECTS], scn_sipip[MAXDATAOBJECTS], scn_r9[MAXDATAOBJECTS];

  //leptons
  Int_t ln;
  Int_t ln_id[MAXDATAOBJECTS],          ln_idbits[MAXDATAOBJECTS],    ln_genid[MAXDATAOBJECTS],   ln_Tbits[MAXDATAOBJECTS];
  Float_t ln_px[MAXDATAOBJECTS],        ln_py[MAXDATAOBJECTS],        ln_pz[MAXDATAOBJECTS],      ln_en[MAXDATAOBJECTS];
  Float_t ln_ecalIso[MAXDATAOBJECTS],   ln_hcalIso[MAXDATAOBJECTS],   ln_trkIso[MAXDATAOBJECTS];
  Float_t ln_gIso[MAXDATAOBJECTS],      ln_chIso[MAXDATAOBJECTS],     ln_puchIso[MAXDATAOBJECTS], ln_nhIso[MAXDATAOBJECTS]; 
  Float_t ln_ip3d[MAXDATAOBJECTS],      ln_ip3dsig[MAXDATAOBJECTS];
  
  //jets
  Int_t jn, jn_idbits[MAXDATAOBJECTS];
  Float_t jn_px[MAXDATAOBJECTS],    jn_py[MAXDATAOBJECTS],      jn_pz[MAXDATAOBJECTS],          jn_en[MAXDATAOBJECTS], jn_torawsf[MAXDATAOBJECTS], jn_csv[MAXDATAOBJECTS], jn_area[MAXDATAOBJECTS];
  Int_t   jn_genflav[MAXDATAOBJECTS], jn_genid[MAXDATAOBJECTS];
  Float_t jn_genpx[MAXDATAOBJECTS], jn_genpy[MAXDATAOBJECTS], jn_genpz[MAXDATAOBJECTS], jn_genen[MAXDATAOBJECTS];
  
  //met 
  Int_t metn;
  Float_t met_pt[MAXDATAOBJECTS], met_phi[MAXDATAOBJECTS],met_sig[MAXDATAOBJECTS],met_sigx2[MAXDATAOBJECTS],met_sigxy[MAXDATAOBJECTS],met_sigy2[MAXDATAOBJECTS];

  SMEventSummary() { } 
  ~SMEventSummary() { }

  //
  inline void reset()
  {
    run=0;    
    lumi=0;   
    event=0;  
    tn=0;
    scn=0;
    mcn=0;
    ln=0; 
    jn=0;
    metn=0;
  }

  //
  inline void fill() { for(size_t i=0; i<m_trees.size(); i++) m_trees[i]->Fill(); }
  
  //
  inline void attach(TTree *t)
  {
  
  }

  //
  inline void create(TTree *t)
  {

    m_trees.push_back(t);
  }


 private:
  std::vector<TTree *> m_trees;

};

#endif
