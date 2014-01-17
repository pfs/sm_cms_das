// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
// #include "Rivet/RivetYODA.hh"
// #include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ZFinder.hh"
// #include "Rivet/Particle.fhh"

namespace Rivet {


  /// @brief ATLAS Z pT in Drell-Yan events at 7 TeV
  /// @author Elena Yatsenko, Judith Katzy
  class CMS_2012_I941555 : public Analysis {
  public:

    /// Constructor
    CMS_2012_I941555()
    : Analysis("CMS_2012_I941555")
    {
      _sumw_mu_dressed_pt  = 0;
      _sumwpeak_mu_dressed = 0;
      _sumw_el_dressed_rap = 0;
      _sumw_el_dressed_pt  = 0;
      _sumwpeak_el_dressed = 0;
    }


    /// @name Analysis methods
    //@{

    void init() {

      // Set up projections
      // FinalState fs;
      ZFinder zfinder_dressed_mu_pt(-2.1, 2.1, 20, MUON, 60.0*GeV, 120.0*GeV, 0.2, false, true ); 
      addProjection(zfinder_dressed_mu_pt, "ZFinder_dressed_mu_pt");
      ZFinder zfinder_dressed_el_pt(-2.1, 2.1, 20, ELECTRON, 60.0*GeV, 120.0*GeV, 0.1, true, false);
      addProjection(zfinder_dressed_el_pt, "ZFinder_dressed_el_pt");
      
      ZFinder zfinder_dressed_mu_rap(-MAXRAPIDITY, MAXRAPIDITY, 0, MUON, 60.0*GeV, 120.0*GeV, 0.1, true, false);
      addProjection(zfinder_dressed_mu_rap, "ZFinder_dressed_mu_rap");
      ZFinder zfinder_dressed_el_rap(-MAXRAPIDITY, MAXRAPIDITY, 0, ELECTRON, 60.0*GeV, 120.0*GeV, 0.1, true, false);
      addProjection(zfinder_dressed_el_rap, "ZFinder_dressed_el_rap");

      // Book histograms
      // _hist_zrap_mu_dressed     = bookHistogram1D(1, 1, 1);  // muon "dressed" rapidity
      // _hist_zrap_el_dressed     = bookHistogram1D(1, 1, 2);  // electron "dressed" rapidity
      // _hist_zrap_comb_dressed     = bookHistogram1D(1, 1, 3);  // electron "dressed" rapidity
      
      // _hist_zpt_mu_dressed      = bookHistogram1D(2, 1, 1);  // muon "dressed" pt
      // _hist_zpt_el_dressed      = bookHistogram1D(2, 1, 2);  // electron "dressed" pt
      _hist_zpt_comb_dressed      = bookHistogram1D(2, 1, 3);  // electron "dressed" pt
      
      // _hist_zptpeak_mu_dressed  = bookHistogram1D(3, 1, 1);  // muon "dressed" pt peak
      // _hist_zptpeak_el_dressed  = bookHistogram1D(3, 1, 2);  // electron "dressed" pt peak
      // _hist_zptpeak_comb_dressed  = bookHistogram1D(3, 1, 3);  // electron "dressed" pt peak
    }


    /// Do the analysis
    void analyze(const Event& evt) {
    
      const double weight = evt.weight();

      const ZFinder& zfinder_dressed_mu_rap = applyProjection<ZFinder>(evt, "ZFinder_dressed_mu_rap");
      
      if (!zfinder_dressed_mu_rap.bosons().empty()) {
      
        _sumw_mu_dressed_rap += weight;
        const FourMomentum pZ = zfinder_dressed_mu_rap.bosons()[0].momentum();
        // _hist_zrap_mu_dressed->fill(pZ.rapidity()/GeV, weight);
        // _hist_zrap_comb_dressed->fill(pZ.rapidity()/GeV, weight);
      
      }

      const ZFinder& zfinder_dressed_mu_pt = applyProjection<ZFinder>(evt, "ZFinder_dressed_mu_pt");
      
      if(!zfinder_dressed_mu_pt.bosons().empty()) {
      
        _sumw_mu_dressed_pt += weight;
        const FourMomentum pZ = zfinder_dressed_mu_pt.bosons()[0].momentum();
        // _hist_zpt_mu_dressed->fill(pZ.pT()/GeV, weight);
        _hist_zpt_comb_dressed->fill(pZ.pT()/GeV, weight);
        
        if(pZ.pT()/GeV<30){
        
          _sumwpeak_mu_dressed += weight;
          // _hist_zptpeak_mu_dressed->fill(pZ.pT()/GeV, weight);
          // _hist_zptpeak_comb_dressed->fill(pZ.pT()/GeV, weight);
        
        }
        
      }

      const ZFinder& zfinder_dressed_el_rap = applyProjection<ZFinder>(evt, "ZFinder_dressed_el_rap");
      
      if(!zfinder_dressed_el_rap.bosons().empty()) {
      
        _sumw_el_dressed_rap += weight;
        const FourMomentum pZ = zfinder_dressed_el_rap.bosons()[0].momentum();
        // _hist_zrap_el_dressed->fill(pZ.rapidity()/GeV, weight);
        // _hist_zrap_comb_dressed->fill(pZ.rapidity()/GeV, weight);
      
      }

      const ZFinder& zfinder_dressed_el_pt = applyProjection<ZFinder>(evt, "ZFinder_dressed_el_pt");
      
      if(!zfinder_dressed_el_pt.bosons().empty()) {
      
        _sumw_el_dressed_pt += weight;
        const FourMomentum pZ = zfinder_dressed_el_pt.bosons()[0].momentum();
        // _hist_zpt_el_dressed->fill(pZ.pT()/GeV, weight);
        _hist_zpt_comb_dressed->fill(pZ.pT()/GeV, weight);
        
        if(pZ.pT()/GeV<30){
        
          _sumwpeak_el_dressed += weight;
          // _hist_zptpeak_el_dressed->fill(pZ.pT()/GeV, weight);
          // _hist_zptpeak_comb_dressed->fill(pZ.pT()/GeV, weight);
        
        }
        
      }

      
    } // close class definition


    void finalize() {
    
      // scale(_hist_zrap_mu_dressed, 1/(_sumw_mu_dressed_rap>0?_sumw_mu_dressed_rap:1));
      // scale(_hist_zpt_mu_dressed, 1/(_sumw_mu_dressed_pt>0?_sumw_mu_dressed_pt:1));
      // scale(_hist_zptpeak_mu_dressed, 1/(_sumwpeak_mu_dressed>0?_sumwpeak_mu_dressed:1));
      
      // scale(_hist_zrap_el_dressed, 1/(_sumw_el_dressed_rap>0?_sumw_el_dressed_rap:1));
      // scale(_hist_zpt_el_dressed, 1/(_sumw_el_dressed_pt>0?_sumw_el_dressed_pt:1));
      // scale(_hist_zptpeak_el_dressed, 1/(_sumwpeak_el_dressed>0?_sumwpeak_el_dressed:1));

      // scale(_hist_zrap_comb_dressed, 1/(_sumw_el_dressed_rap+_sumw_mu_dressed_rap>0?_sumw_el_dressed_rap+_sumw_mu_dressed_rap:1));
      scale(_hist_zpt_comb_dressed, 1/(_sumw_el_dressed_pt+_sumw_mu_dressed_pt>0?_sumw_el_dressed_pt+_sumw_mu_dressed_pt:1));
      // scale(_hist_zptpeak_comb_dressed, 1/(_sumwpeak_el_dressed+_sumwpeak_mu_dressed>0?_sumwpeak_el_dressed+_sumwpeak_mu_dressed:1));
    }

    //@}


  private:

    double _sumw_mu_dressed_rap;
    double _sumw_mu_dressed_pt;
    double _sumwpeak_mu_dressed;
    
    double _sumw_el_dressed_rap;
    double _sumw_el_dressed_pt;
    double _sumwpeak_el_dressed;

    // AIDA::IHistogram1D * _hist_zrap_mu_dressed;
    // AIDA::IHistogram1D * _hist_zpt_mu_dressed;
    // AIDA::IHistogram1D * _hist_zptpeak_mu_dressed;
    
    // AIDA::IHistogram1D * _hist_zrap_el_dressed;
    // AIDA::IHistogram1D * _hist_zpt_el_dressed;
    // AIDA::IHistogram1D * _hist_zptpeak_el_dressed;
    
    // AIDA::IHistogram1D * _hist_zrap_comb_dressed;
    AIDA::IHistogram1D * _hist_zpt_comb_dressed;
    // AIDA::IHistogram1D * _hist_zptpeak_comb_dressed;
    
  };


  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CMS_2012_I941555> plugin_CMS_2012_I941555;

}