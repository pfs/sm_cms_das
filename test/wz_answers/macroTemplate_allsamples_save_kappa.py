import math,ROOT
from UserCode.sm_cms_das.Tools import *

h=ROOT.TH1F('leg2pt',';Missing transverse energy [GeV];Events',50,0,100)
h.Sumw2()
hsbmc=h.Clone("mcsideband_leg2pt")
hsbdata=h.Clone("datasideband_leg2pt")
normBin=h.GetXaxis().FindBin(20)
nMCBelowNormBin=0
nDataBelowNormBin=0

fOut=ROOT.TFile('plotter.root','RECREATE')

processes=getProcesses(jsonUrl="test/wz/wz_samples.json",getMC=True,getData=True)

for idx,p in enumerate(processes):
    files=getFilesForProcess(jsonUrl="test/wz/wz_samples.json",tag=p,inDir="/store/cmst3/user/psilva/CMSDAS_v2/summary/")
    lumi=19
    xsec=0
    totalSelected=0
    totalGenerated=0

    hproc=h.Clone(p+'_leg2pt')
    hproc.SetTitle(p)
    hproc.SetFillStyle(1001)
    hproc.SetFillColor(idx+1)
    hproc.SetDirectory(0)

    for f in files:
        inF=ROOT.TFile.Open(f)
        tree=inF.Get("data/data")
        totalSelected=totalSelected+tree.GetEntries("abs(cat)==13")
        totalGenerated=totalGenerated+inF.Get("iniEvents")[0]
        xsec=inF.Get("crossSection")[0]
        h.SetDirectory(inF)
        tree.Draw('leg2_pt>>leg2pt','abs(cat)==13 && leg1_pt>25','goff')
        hproc.Add(h)
        tree.Draw('leg2_pt>>leg2pt','abs(cat)==1300 && leg1_pt>25','goff')
        if p in ['SingleMu','SingleElectron']:
            hsbdata.Add(h)
        else:
            hsbmc.Add(h)
        h.SetDirectory(0)
        inF.Close()

    if xsec>=0 :
        nExpected=xsec*lumi*totalSelected/totalGenerated
        nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
        hproc.Scale(xsec*lumi/totalGenerated)
        nMCBelowNormBin=nMCBelowNormBin+hproc.Integral(0,normBin)
    else :
        nExpected=totalSelected
        nExpectedErr=0
        nDataBelowNormBin=nDataBelowNormBin+hproc.Integral(0,normBin)
        
    print "N(%s)=%3.3f +/- %3.3f"%(p,nExpected,nExpectedErr)
    
    fOut.mkdir(p).cd()
    hproc.Write()

hsbdata.Add(hsbmc,-1)
for xbin in xrange(1,hsbdata.GetXaxis().GetNbins()):
    y=hsbdata.GetBinContent(xbin)
    if y<0 :
        hsbdata.SetBinContent(xbin,0)
        #hsbdata.SetBinError(xbin,math.sqrt(math.fabs(y)))
        hsbdata.SetBinError(xbin,math.fabs(y))
kMJ=(nDataBelowNormBin-nMCBelowNormBin)/hsbdata.Integral(0,normBin)
print 'k_{MJ}=%3.4f'%(kMJ)
hsbdata.Scale(kMJ)
fOut.mkdir('MJ').cd()
hsbdata.Write('MJ_leg2pt')
    
fOut.Close()
