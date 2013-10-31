import os,sys

isMC=True
gtag="START53_V23::All"
fileList=cms.untracked.vstring('/store/relval/CMSSW_5_3_12_patch2/RelValWM/GEN-SIM-RECO/START53_LV2-v1/00000/3E5A329C-B02B-E311-A039-003048FF9AA6.root')

execfile( os.path.expandvars('${CMSSW_BASE}/src/UserCode/sm_cms_das/test/runSManalyzer_cfg.py'))

