#!/usr/bin/env python

import optparse
import os
import sys
import getopt
import math
import array
import random
import commands
from ROOT import gSystem
from ROOT import TTree, TFile, TLorentzVector, TH1F, TH2F, TGraph, TObjArray, TNtuple

class Monitor:
    def __init__(self,outUrl,norm):
        self.outUrl=outUrl
        self.norm=norm
        self.allHistos={}
        self.extra=TObjArray()
    def addToMonitor(self,h,tag):
        h.SetDirectory(0)
        name=h.GetName()
        if tag != 'base' : h.SetName(tag+'_'+name) 
        else :             h.Sumw2()
        if name in self.allHistos :
            self.allHistos[name].update({tag:h})
        else:
            self.allHistos[name]={'base':h}
    def addObject(self,obj):
        self.extra.Add(obj)
    def addHisto(self,name,title,nx,xmin,xmax) :
        h=TH1F(name,title,nx,xmin,xmax)
        self.addToMonitor(h,'base')
    def add2DHisto(self,name,title,nx,xmin,xmax,ny,ymin,ymax):
        h=TH2F(name,title,nx,xmin,xmax,ny,ymin,ymax)
        self.addToMonitor(h,'base')
    def fill(self,name,tags,valx,valy,valz=None):
        if name in self.allHistos:
            for t in tags:
                if not (t in self.allHistos[name]) :
                    self.addToMonitor(self.allHistos[name]['base'].Clone(),t)
                h=self.allHistos[name][t]
                if valz is None : h.Fill(valx,valy)
                else :            h.Fill(valx,valy,valz)
    def close(self):
        fOut=TFile.Open(self.outUrl,'RECREATE')
        for key in self.allHistos:
            for tag in self.allHistos[key] :
                self.allHistos[key][tag].Scale(self.norm)
                self.allHistos[key][tag].Write()
        for i in xrange(0,self.extra.GetEntriesFast()):
            obj=self.extra.At(i)
            outDir=fOut.mkdir(obj.GetName())
            outDir.cd()
            obj.SetDirectory(outDir)
            obj.Write()
            fOut.cd()
        fOut.Close()
                
"""
Wrapper for a lepton candidate object
"""
class LeptonCand:
    def __init__(self,id,px,py,pz,en):
        self.id=id
        self.p4=TLorentzVector(px,py,pz,en)
        self.Tbits=0
        self.relIso=99999.
        self.passLoose=False
        self.passLooseIso=False
        self.passTight=False
        self.passTightIso=False
    def triggerInfo(self,Tbits):
        self.Tbits=Tbits
    def selectionInfo(self,relIso,passLoose,passLooseIso,passTight,passTightIso):
        self.relIso
        self.passLoose=passLoose
        self.passLooseIso=passLooseIso
        self.passTight=passTight
        self.passTightIso=passTightIso

"""
Wrapper for a vector boson candidate
"""
class VectorBosonCand:
    def __init__(self, id,tag):
        self.id=id
        self.tag=tag
        self.m_legs=[]
    def addLeg(self,LeptonCand):
        if len(self.m_legs) >2 : return
        self.m_legs.append(LeptonCand)
        if len(self.m_legs)<2: return
        p1=self.m_legs[0].p4
        p2=self.m_legs[1].p4
        self.p4=p1+p2
        dphi=p1.DeltaPhi(p2)
        self.mt=math.sqrt(2*p1.Pt()*p2.Pt()*(1-math.cos(dphi)))
    def getTag(self) :
        return self.tag

def decodeTriggerWord(trigBits) :
    eFire     = (((trigBits >> 0) & 0x1)>0) or (((trigBits >> 1) & 0x1)>0)
    mFire     = (((trigBits >> 2) & 0x1)>0) or (((trigBits >> 3) & 0x1)>0)
    emFire    = (((trigBits >> 4) & 0x3)>0)
    return eFire,mFire,emFire

"""
Performs the selection of a good lepton
"""
def selectLepton(id, idbits, gIso, chIso, nhIso, puchIso, pt) :
    isLoose=False
    isTight=False
    isLooseIso=False
    isTightIso=False

    relIso=9999.
    if abs(id)==11 :
        isLoose = ((idbits >> 3) & 0x1)
        isTight = ((idbits >> 4) & 0x1)
        relIso=(chIso+nhIso+gIso)/pt
        if relIso<0.15: isLooseIso=True
        if relIso<0.15: isTightIso=True
        
    if abs(id)==13:
        isLoose = ((idbits >> 8) & 0x1)
        isTight = ((idbits >> 9) & 0x1)
        relIso=(chIso+nhIso+gIso)/pt
        if relIso<0.20: isLooseIso=True
        if relIso<0.12: isTightIso=True

    return relIso, isLoose, isLooseIso, isTight, isTightIso

"""
Returns a vector boson candidate depending on the kinematics, id and isolation of its decay leg candidates
"""
def buildVcand(eFire,mFire,emFire,leptonCands,met) :

    vCand=None
    if len(leptonCands)==0 : return vCand

    tightLeptons=[]
    tightNonIsoLeptons=[]
    vetoLeptons=[]
    for l in leptonCands :
        if   l.passTight and l.passTightIso     : tightLeptons.append(l)
        elif l.passLoose and l.passLooseIso     : vetoLeptons.append(l)
        elif l.passTight and not l.passLooseIso : tightNonIsoLeptons.append(l)

    # test first hypotheses: muon channels -> require muon trigger
    # a) 2 tight muons = Z->mm
    # b) 1 tight muon and no loose lepton = W->mv
    if vCand is None and mFire :
        if len(tightLeptons)>=2 and tightLeptons[0].id*tightLeptons[1].id==-13*13 :
            vCand = VectorBosonCand(23,'mumu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])
        elif len(tightLeptons)==1 and abs(tightLeptons[0].id)==13 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(24,'mu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(met)
        elif len(tightNonIsoLeptons)==1 and abs(tightNonIsoLeptons[0].id)==13 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(24,'munoniso')
            vCand.addLeg(tightNonIsoLeptons[0])
            vCand.addLeg(met)

    # test second hypotheses: electron channels -> require electron trigger
    # a) 2 tight electrons = Z->ee
    # b) 1 tight electron and no loose lepton = W->ev
    if vCand is None and eFire:
        if len(tightLeptons)>=2 and tightLeptons[0].id*tightLeptons[1].id==-11*11 :
            vCand = VectorBosonCand(23,'ee')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])
        elif len(tightLeptons)==1 and abs(tightLeptons[0].id)==11 and len(vetoLeptons)==0:
            vCand = VectorBosonCand(24,'e')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(met)
        elif len(tightNonIsoLeptons)==1 and abs(tightNonIsoLeptons[0].id)==11 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(24,'enoniso')
            vCand.addLeg(tightNonIsoLeptons[0])
            vCand.addLeg(met)

    # test third hypothesis with the emu trigger
    # a) 1 tight electron, 1 tight muon = Z->tt
    if vCand is None and emFire:
        if len(tightLeptons)>=2 and tightLeptons[0].id*tightLeptons[1].id==-11*13 :
            vCand = VectorBosonCand(23,'emu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])

    return vCand

"""
Runs the event selection over the trees
"""
def selectEvents(fileName,saveProbes=False,saveSummary=False,outputDir='./',xsec=-1,correctionsMap={}):

    gSystem.ExpandPathName(fileName)
    file=TFile.Open(fileName)
    yieldsNorm=1.0
    puWeightsGr=None
    
    if xsec>0 :
        origEvents=file.Get('smDataAnalyzer/cutflow').GetBinContent(1)
        if origEvents==0 :
            print '[Warning] normalization is requested but got 0 for original number of events - ignoring'
        else :
            yieldsNorm=xsec/origEvents
            print 'Applying overall normalization for xsec=%f pb and # original events=%d'%(xsec,origEvents)

        #derive pileup weights
        origPileup=file.Get('smDataAnalyzer/pileup')
        try:
            dataPileupFile=TFile.Open(correctionsMap['pu'])
            dataPileup=dataPileupFile.Get('pileup')
            normF=origPileup.Integral()/dataPileup.Integral()
            if normF>0 :
                puWeightsGr=TGraph()
                for xbin in xrange(1,origPileup.GetXaxis().GetNbins()+1) :
                    iweight=1.0
                    if origPileup.GetBinContent(xbin)>0 :
                        iweight=normF*dataPileup.GetBinContent(xbin)/origPileup.GetBinContent(xbin)
                        puWeightsGr.SetPoint( puWeightsGr.GetN(), origPileup.GetXaxis().GetBinCenter(xbin), iweight )
            dataPileupFile.Close()
        except : 
            print 'No data pileup filed provided or other error occurred. If you wish add -w pu,pu_file.root'
    

    tree=file.Get("smDataAnalyzer/data")
    nev = tree.GetEntries()

    outUrl=outputDir+'/'+os.path.basename(fileName)
    monitor=Monitor(outUrl,yieldsNorm)
    monitor.addHisto('nvtx',    ';Vertices;Events',                       50,0,50)
    monitor.addHisto('nvtxraw', ';Vertices;Events',                       50,0,50)
    monitor.addHisto('vmass',   ';Mass [GeV];Events',                     50,0,250)
    monitor.addHisto('vmt',     ';Transverse mass [GeV];Events',          50,0,250)
    monitor.addHisto('vpt',     ';Boson transverse momentum [GeV];Events',50,0,250)
    monitor.addHisto('leg1pt',  ';Transverse momentum [GeV];Events',      50,0,250)
    monitor.addHisto('leg2pt',  ';Transverse momentum [GeV];Events',      50,0,250)
    monitor.addHisto('leg1iso', ';Relative isolation;Events',             50,0,0.5)
    monitor.addHisto('leg2iso', ';Relative isolation;Events',             50,0,0.5)

    summaryTuple=None
    if saveSummary :
        summaryTuple=TNtuple('data','summary','cat:weight:v_mass:v_mt:v_pt:leg1_pt:leg1_eta:leg1_phi:leg2_pt:leg2_eta:leg2_phi')
        summaryTuple.SetDirectory(0)
        monitor.addObject(summaryTuple)

    probesTuple=None
    probesId   = array.array( 'f', [ 0 ] )
    probesPt   = array.array( 'f', [ 0 ] )
    probesEta  = array.array( 'f', [ 0 ] )
    probesPhi  = array.array( 'f', [ 0 ] )
    probesNvtx = array.array( 'f', [ 0 ] )
    probesMass = array.array( 'f', [ 0 ] )
    probesPassLoose = array.array( 'i', [ 0 ] )
    probesPassTight = array.array( 'i', [ 0 ] )
    probesFireTrigger = array.array( 'i', [ 0 ] )
    if saveProbes :
        probesTuple=TTree('tandp','summary for tandp')
        probesTuple.Branch( 'id', probesId, 'id/F' )
        probesTuple.Branch( 'pt', probesPt, 'pt/F' )
        probesTuple.Branch( 'eta', probesEta, 'eta/F' )
        probesTuple.Branch( 'phi', probesPhi, 'phi/F' )
        probesTuple.Branch( 'nvtx', probesNvtx, 'nvtx/F' )
        probesTuple.Branch( 'mass', probesMass, 'mass/F' )
        probesTuple.Branch( 'passLoose', probesPassLoose, 'passLoose/I' )
        probesTuple.Branch( 'passTight', probesPassTight, 'passTight/I' )
        probesTuple.Branch( 'fireTrigger', probesFireTrigger, 'fireTrigger/I' )
        probesTuple.SetDirectory(0)
        monitor.addObject(probesTuple)

    #
    # LOOP OVER THE EVENTS 
    #
    for iev in range(0,nev):
        tree.GetEntry(iev)
        if iev%10000 == 0 :
            sys.stdout.write("\r[ %d/100 ] completed" %(100.*iev/nev))
            sys.stdout.flush()

        #get triggers that fired
        eFire,mFire,emFire=decodeTriggerWord(tree.tbits)

        #select the leptons
        leptonCands=[]
        validTags=[]
        for l in xrange(0,tree.ln) :
            lep=LeptonCand(tree.ln_id[l],tree.ln_px[l],tree.ln_py[l],tree.ln_pz[l],tree.ln_en[l])
            if lep.p4.Pt()<20 : continue
            if abs(tree.ln_id[l])==11 :
                if math.fabs(lep.p4.Eta())>2.5 : continue
                if math.fabs(lep.p4.Eta())>1.4442 and math.fabs(lep.p4.Eta())<1.566 : continue
            if abs(tree.ln_id[l])==13 :
                if math.fabs(lep.p4.Eta())>2.1 : continue
            relIso, isLoose, isLooseIso, isTight, isTightIso = selectLepton(tree.ln_id[l],tree.ln_idbits[l],tree.ln_gIso[l],tree.ln_chIso[l],tree.ln_nhIso[l],tree.ln_puchIso[l],lep.p4.Pt())
            lep.selectionInfo(relIso,isLoose, isLooseIso, isTight, isTightIso)
            lep.triggerInfo(tree.ln_Tbits[l])
            leptonCands.append(lep)

            if not saveProbes: continue
            if not isTight or not isTightIso or lep.Tbits==0 : continue
            if abs(lep.id)==11 and not eFire: continue
            if abs(lep.id)==13 and not mFire: continue
            validTags.append( len(leptonCands)-1 )  

        #check if probes tree should be saved 
        if saveProbes and len(validTags)>0:

            # choose a random tag 
            tagIdx=random.choice(validTags)
            tag=leptonCands[tagIdx]

            #find probe
            probe=None
            for l in xrange(0,len(leptonCands)) :
                if l==tagIdx: continue
                if abs(tag.id)!=abs(leptonCands[l].id) : continue
                probe=leptonCands[l]
                break

            #for electrons save superclusters if probe is not found
            matchToEle=1
            if abs(tag.id)==11 and probe is None :
                hasScMatch=0
                for sc in xrange(0,tree.scn) :
                    sc_en=tree.scn_e[sc]
                    sc_eta=tree.scn_eta[sc]
                    sc_phi=tree.scn_phi[sc]
                    sc_pt=sc_en/math.cosh(sc_eta)
                    sc_p4=TLorentzVector(0,0,0,0)
                    sc_p4.SetPtEtaPhiE(sc_pt,sc_eta,sc_phi,sc_en)
                    lscp4=tag.p4+sc_p4
                    if math.fabs(lscp4.M()-91)>30 : continue
                    scCand=LeptonCand(tag.id,sc_p4.Px(),sc_p4.Py(),sc_p4.Pz(),sc_p4.E())
                    scCand.selectionInfo(0,0,0,0,0)
                    scCand.triggerInfo(0)
                    probe=scCand
                    break
            if abs(tag.id)==13 : matchToEle=0
            
            #save info
            if probe is not None:
                tpp4=tag.p4+probe.p4
                if math.fabs(tpp4.M()-91)<30 :
                    probesId[0]=probe.id
                    probesPt[0]=probe.p4.Pt()
                    probesEta[0]=probe.p4.Eta()
                    probesPhi[0]=probe.p4.Phi()
                    probesNvtx[0]=tree.nvtx
                    probesMass[0]=tpp4.M()
                    probesPassLoose[0]=(probe.passLoose and probe.passLooseIso)
                    probesPassTight[0]=(probe.passTight and probe.passTightIso)
                    probesFireTrigger[0]=(probe.Tbits>0)
                    probesTuple.Fill()

        # met
        metCand=LeptonCand(0,tree.met_pt[0]*math.cos(tree.met_phi[0]),tree.met_pt[0]*math.sin(tree.met_phi[0]),0,tree.met_pt[0])
        
        #build the candidate
        vCand=buildVcand(eFire,mFire,emFire,leptonCands,metCand)
        if vCand is None : continue

        #prepare to save
        weight=1.0
        if puWeightsGr is not None:
            weight=puWeightsGr.Eval(tree.ngenITpu)

        #show isolations
        for ileg in [0,1]:
            hname='leg'+str(ileg+1)+'iso'
            lid=''
            if abs(vCand.m_legs[ileg].id)==11 :   lid='e'
            elif abs(vCand.m_legs[ileg].id)==13 : lid='mu'
            else : continue
            monitor.fill(hname,[lid],vCand.m_legs[ileg].relIso,weight)

        tags=[vCand.tag]
        monitor.fill('nvtxraw',tags, tree.nvtx,               1.0)
        monitor.fill('nvtx',   tags, tree.nvtx,               weight)
        monitor.fill('vmass',  tags, vCand.p4.M(),            weight)
        monitor.fill('vmt',    tags, vCand.mt,                weight)
        monitor.fill('vpt',    tags, vCand.p4.Pt(),           weight)
        monitor.fill('leg1pt', tags, vCand.m_legs[0].p4.Pt(), weight)
        monitor.fill('leg2pt', tags, vCand.m_legs[1].p4.Pt(), weight)

        if saveSummary :
            values=[vCand.id, weight*yieldsNorm,
                    vCand.p4.M(), vCand.mt, vCand.p4.Pt(),
                    vCand.m_legs[0].p4.Pt(),vCand.m_legs[0].p4.Eta(),vCand.m_legs[0].p4.Phi(),
                    vCand.m_legs[1].p4.Pt(),vCand.m_legs[1].p4.Eta(),vCand.m_legs[1].p4.Phi()
                    ]
            summaryTuple.Fill(array.array("f",values))
               
    file.Close()
    monitor.close()


def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inputFile'  ,    dest='input'  ,     help='Name of the local input file',                                    default=None)
    parser.add_option('-p', '--saveProbes' ,    dest='saveProbes',  help='Save probes tree for tag and probe',                              default=False, action="store_true")
    parser.add_option('-t', '--saveSummary' ,   dest='saveSummary', help='Save summary with selected events',                               default=False, action="store_true")
    parser.add_option('-o', '--outputDir'  ,    dest='outputDir'  , help='Name of the local output directory (default=./)',                 default='./')
    parser.add_option('-w', '--weights',        dest='weights',     help='CSV weights to apply (e.g. pu,pu_file.root,eff,eff_file.root)',   default='')
    parser.add_option('-x', '--xsec',           dest='xsec'       , help='If given is used to normalize the final histograms (default=-1)', default=-1.0, type='float')
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    #decode weights to apply
    correctionsMap={}
    if len(opt.weights)>0 :
        tkns=opt.weights.split(',')
        for itkn in xrange(0,len(tkns)+1):
            try :
                key=tkns[itkn]
                wgtFile=tkns[itkn+1]
                if wgtFile.find('/store/')== 0:
                    wgtFile=commands.getstatusoutput('cmsPfn ' + wgtFile)[1]
                correctionsMap[key] = wgtFile
                itkn=itkn+1
            except:
                continue

    print '[runEventSelection] analyzing %s'%opt.input 
    selectEvents(fileName=opt.input, saveProbes=opt.saveProbes, saveSummary=opt.saveSummary, outputDir=opt.outputDir, xsec=opt.xsec, correctionsMap=correctionsMap)


if __name__ == "__main__":
    main()
