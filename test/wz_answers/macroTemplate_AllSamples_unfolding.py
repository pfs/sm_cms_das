import math,ROOT
from UserCode.sm_cms_das.Tools import *
from array import array

#########################################################################
######### THESE DEFINITIONS NEED TO MATCH THOSE IN unfold_zpt.C #########
#########################################################################
# Number of bins in final, unfolded distribution. Also used for histograms containing MC truth information.
nbinsunfolded = 18;
bin_edges_gen= [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 30.0, 40.0, 50.0, 70.0, 90.0, 110.0, 150.0, 190.0, 250.0, 600.0];
# Number of bins in reconstructed distribution. Using twice as many bins helps the stability of method.
nbinsreconstructed = nbinsunfolded * 2;
bin_edges_rec= [0.0, 1.5, 2.6, 3.5, 4.9, 5.9, 7.5, 8.5, 10.5, 11.25, 12.6, 13.9, 15.25, 16.4, 17.5, 18.75, 20.5, 25.8, 30.1, 36.0, 41.0, 46.0, 50.0, 61.0, 71.0, 82.0, 90.0, 101.0, 110.0, 132.0, 152.0, 173.0, 191.0, 226.0, 255.0, 500.0, 1000.0];
#########################################################################
#########################################################################

#declare the histograms needed for the unfolding:
# the reconstructed boson pT
hvpt=ROOT.TH1F('vpt',';Boson transverse momentum [GeV];Events',nbinsreconstructed, array('d',bin_edges_rec))
hvpt.Sumw2()
# the generated boson pT
hgenvpt=ROOT.TH1F('genvpt',';Boson transverse momentum [GeV];Events',nbinsunfolded, array('d',bin_edges_gen))
hgenvpt.Sumw2()
# reconstructed vs generated boson pT
hvptmigrationmatrix=ROOT.TH2F('vptmigrationmatrix',';Generated boson transverse momentum [GeV];Reconstructed boson transverse momentum [GeV]', nbinsunfolded, array('d',bin_edges_gen), nbinsreconstructed, array('d',bin_edges_rec))
hvptmigrationmatrix.Sumw2()

#open a ROOT file to contain the histograms
fOut=ROOT.TFile('plotter_unfolding.root','RECREATE')

#get all samples of interest for the analysis: both data and MC
processes=getProcesses(jsonUrl="../../test/wz/wz_samples.json",getMC=True,getData=True)

#loop over the different processes
for idx,p in enumerate(processes):
    files=getFilesForProcess(jsonUrl="../../test/wz/wz_samples.json",tag=p,inDir="/store/cmst3/user/psilva/CMSDAS_v2/summary/")
    lumi=19  #update the luminosity to the one in data
    xsec=0
    totalSelected=0
    totalGenerated=0

    #replicate the base histogram for a new process
    hvptproc=hvpt.Clone(p+'_'+hvpt.GetName())
    hvptproc.SetDirectory(0)
    hgenvptproc=hgenvpt.Clone(p+'_'+hgenvpt.GetName())
    hgenvptproc.SetDirectory(0)
    hvptmigrationmatrixproc=hvptmigrationmatrix.Clone(p+'_'+hvptmigrationmatrix.GetName())
    hvptmigrationmatrixproc.SetDirectory(0)
    
    #loop over files
    for f in files:
        inF=ROOT.TFile.Open(f)
        tree=inF.Get("data/data")

        #increment the total selected and generated
        totalSelected=totalSelected+tree.GetEntries("cat==-169")
        totalGenerated=totalGenerated+inF.Get("iniEvents")[0]
        xsec=inF.Get("crossSection")[0]
       
        #Method 1: project the tree to the base histogram and add it to the total expectations for this process
        hvpt.SetDirectory(inF)
        tree.Draw('v_pt>>vpt','cat==-169','goff') # plots the reco boson pT (v_pt) variable in the histogram named "vpt" for events with cat==-169 (Z di-muon)
        hvptproc.Add(hvpt)
        hvpt.SetDirectory(0)
        
        hgenvpt.SetDirectory(inF)
        tree.Draw('genv_pt>>genvpt','weight','goff') # plots the gen boson pT (genv_pt) variable in the histogram named "genvpt" for ALL THE GEN events (!!!) reweighting the pileup
        hgenvptproc.Add(hgenvpt)
        hgenvpt.SetDirectory(0)
        
        hvptmigrationmatrix.SetDirectory(inF)
        tree.Draw('genv_pt:v_pt>>vptmigrationmatrix','weight*(cat==-169)','goff') # plots the 2D response matrix reco vs gen boson pT (genv_pt:v_pt) variable in the histogram named "vptmigrationmatrix" for events with cat==-169 (Z di-muon) also reweighting the pileup
        hvptmigrationmatrixproc.Add(hvptmigrationmatrix)
        hvptmigrationmatrix.SetDirectory(0)

        #Method 2: run over the tree
        #for n in xrange(0,tree.GetEntriesFast()) :
        #    tree.GetEntry(n)
        #    if tree.cat!=-169 : continue
        #    effWeight=1.0
        #    puWeight=tree.weight
        #    eventWeight=effWeight*puWeight
        #    hvptproc.Fill               (tree.v_pt,            eventWeight)
        #    hgenvptproc.Fill            (tree.genv_pt,         eventWeight)
        #    hvptmigrationmatrixproc.Fill(tree.genv_pt,tree.v_pt,eventWeight)
        
        inF.Close()

    if xsec>=0 :
        #re-scale the yields according to the cross section
        nExpected=xsec*lumi*totalSelected/totalGenerated
        nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
        hvptproc.Scale(xsec*lumi/totalGenerated)
        hgenvptproc.Scale(xsec*lumi/totalGenerated)
        hvptmigrationmatrixproc.Scale(xsec*lumi/totalGenerated)
    else:
        nExpected=totalSelected
        nExpectedErr=0
    print "N(%s)=%3.3f +/- %3.3f"%(p,nExpected,nExpectedErr)

    #create a directory with the name of the process and store the histogram
    fOut.mkdir(p).cd()
    hvptproc.Write()
    hgenvptproc.Write()
    hvptmigrationmatrixproc.Write()

#close the output file
fOut.Close()
