#!/usr/bin/env python

import optparse
import os
import sys
import getopt
import math

from ROOT import gSystem
from ROOT import TTree, TFile, TLorentzVector, TH1F, TH2F

class Monitor:
    def __init__(self,outUrl,norm):
        self.outUrl=outUrl
        self.norm=norm
        self.allHistos={}
    def addToMonitor(self,h,tag):
        h.SetDirectory(0)
        name=h.GetName()
        if tag != 'base' : h.SetName(tag+'_'+name) 
        else :             h.Sumw2()
        if name in self.allHistos :
            self.allHistos[name].update({tag:h})
        else:
            self.allHistos[name]={'base':h}
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
        fOut.Close()
                
            
class LeptonCand:
    def __init__(self,id,px,py,pz,en):
        self.id=id
        self.p4=TLorentzVector(px,py,pz,en)
        self.Tbits=0
        self.passLoose=False
        self.passLooseIso=False
        self.passTight=False
        self.passTightIso=False
    def triggerInfo(self,Tbits):
        self.Tbits=Tbits
    def selectionInfo(self,passLoose,passLooseIso,passTight,passTightIso):
        self.passLoose=passLoose
        self.passLooseIso=passLooseIso
        self.passTight=passTight
        self.passTightIso=passTightIso

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
    eFire     = (((trigBits >> 0) & 0x1)>0)
    eCtrlFire = (((trigBits >> 1) & 0x1)>0)
    mFire     = (((trigBits >> 2) & 0x1)>0)
    mCtrlFire = (((trigBits >> 3) & 0x1)>0)
    emFire    = (((trigBits >> 4) & 0x3)>0)
    return eFire,eCtrlFire,mFire,mCtrlFire,emFire

def selectLepton(id, idbits, gIso, chIso, nhIso, puchIso, pt) :
    isLoose=False
    isTight=False
    isLooseIso=False
    isTightIso=False
    
    if math.fabs(id)==11 :
        isLoose = ((idbits >> 3) & 0x1)
        isTight = ((idbits >> 4) & 0x1)
        relIso=(chIso+nhIso+gIso)/pt
        if relIso<0.15: isLooseIso=True
        if relIso<0.15: isTightIso=True
        
    if math.fabs(id)==13:
        isLoose = ((idbits >> 8) & 0x1)
        isTight = ((idbits >> 9) & 0x1)
        relIso=(chIso+nhIso+gIso)/pt
        if relIso<0.20: isLooseIso=True
        if relIso<0.12: isTightIso=True

    return isLoose, isLooseIso, isTight, isTightIso


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
        if len(tightLeptons)==2 and tightLeptons[0].id*tightLeptons[1].id==-13*13 :
            vCand = VectorBosonCand(23,'mumu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])
        elif len(tightLeptons)==1 and math.fabs(tightLeptons[0].id)==13 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(24,'mu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(met)
        elif len(tightNonIsoLeptons)==1 and math.fabs(tightNonIsoLeptons[0].id)==13 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(24,'munoniso')
            vCand.addLeg(tightNonIsoLeptons[0])
            vCand.addLeg(met)

    # test second hypotheses: electron channels -> require electron trigger
    # a) 2 tight electrons = Z->ee
    # b) 1 tight electron and no loose lepton = W->ev
    if vCand is None and eFire:
        if len(tightLeptons)==2 and tightLeptons[0].id*tightLeptons[1].id==-11*11 :
            vCand = VectorBosonCand(23,'ee')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])
        elif len(tightLeptons)==1 and math.fabs(tightLeptons[0].id)==11 and len(vetoLeptons)==0:
            vCand = VectorBosonCand(24,'e')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(met)
        elif len(tightNonIsoLeptons)==1 and math.fabs(tightNonIsoLeptons[0].id)==11 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(24,'enoniso')
            vCand.addLeg(tightNonIsoLeptons[0])
            vCand.addLeg(met)

    # test third hypothesis with the emu trigger
    # a) 1 tight electron, 1 tight muon = Z->tt
    if vCand is None and emFire:
        if len(tightLeptons)==2 and tightLeptons[0].id*tightLeptons[1].id==-11*13 :
            vCand = VectorBosonCand(23,'emu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])

    return vCand

"""
Runs the event selection over the trees
"""
def selectEvents(fileName,saveProbes=False,saveSummary=False,outputDir='./',xsec=-1) :

    gSystem.ExpandPathName(fileName)
    file=TFile.Open(fileName)
    yieldsNorm=1.0
    if xsec>0 :
        origEvents=file.Get('smDataAnalyzer/cutflow').GetBinContent(1)
        if origEvents==0 :
            print '[Warning] normalization is requested but got 0 for original number of events - ignoring'
        else :
            yieldsNorm=xsec/origEvents
            print 'Applying overall normalization for xsec=%f pb and # original events=%d'%(xsec,origEvents)
    
    tree=file.Get("smDataAnalyzer/data")
    nev = tree.GetEntries()

    outUrl=outputDir+'/'+os.path.basename(fileName)
    monitor=Monitor(outUrl,yieldsNorm)
    monitor.addHisto('nvtx',  ';Vertices;Events', 50,0,50)
    monitor.addHisto('vmass', ';Mass [GeV];Events',50,0,250)
    monitor.addHisto('vmt',   ';Transverse mass [GeV];Events',50,0,250)
    monitor.addHisto('vpt',   ';Boson transverse momentum [GeV];Events',50,0,250)
    monitor.addHisto('leg1pt',';Transverse momentum [GeV];Events',50,0,250)
    monitor.addHisto('leg2pt',';Transverse momentum [GeV];Events',50,0,250)

    for iev in range(0,nev) :
        tree.GetEntry(iev)

        #get triggers that fired
        eFire,eCtrlFire,mFire,mCtrlFire,emFire=decodeTriggerWord(tree.tbits)

        #select the leptons
        leptonCands=[]
        for l in xrange(0,tree.ln) :
            lep=LeptonCand(tree.ln_id[l],tree.ln_px[l],tree.ln_py[l],tree.ln_pz[l],tree.ln_en[l])
            if lep.p4.Pt()<20 or math.fabs(lep.p4.Eta())>2.5 : continue
            isLoose, isLooseIso, isTight, isTightIso = selectLepton(tree.ln_id[l],tree.ln_idbits[l],tree.ln_gIso[l],tree.ln_chIso[l],tree.ln_nhIso[l],tree.ln_puchIso[l],lep.p4.Pt())
            lep.selectionInfo(isLoose, isLooseIso, isTight, isTightIso)
            lep.triggerInfo(tree.ln_Tbits[l])
            leptonCands.append(lep)
            
        # met
        metCand=LeptonCand(0,tree.met_pt[0]*math.cos(tree.met_phi[0]),tree.met_pt[0]*math.sin(tree.met_phi[0]),0,tree.met_pt[0])
        
        #build the candidate
        vCand=buildVcand(eFire,mFire,emFire,leptonCands,metCand)
        if vCand is None : continue

        weight=1.0
        tags=[vCand.tag]
        monitor.fill('nvtx',   tags, tree.nvtx,               weight)
        monitor.fill('vmass',  tags, vCand.p4.M(),            weight)
        monitor.fill('vmt',    tags, vCand.mt,                weight)
        monitor.fill('vpt',    tags, vCand.p4.Pt(),           weight)
        monitor.fill('leg1pt', tags, vCand.m_legs[0].p4.Pt(), weight)
        monitor.fill('leg2pt', tags, vCand.m_legs[1].p4.Pt(), weight)
               
    file.Close()
    monitor.close()


def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inputFile'  ,    dest='input'  ,     help='Name of the local input file',                    default=None)
    parser.add_option('-p', '--saveProbes' ,    dest='saveProbes',  help='Save probes tree for tag and probe',              default=False, action="store_true")
    parser.add_option('-t', '--saveSummary' ,   dest='saveSummary', help='Save summary with selected events',               default=False, action="store_true")
    parser.add_option('-o', '--outputDir'  ,    dest='outputDir'  , help='Name of the local output directory (default=./)', default='./')
    parser.add_option('-x', '--xsec',           dest='xsec'       , help='If given is used to normalize the final histograms (default=-1)', default=-1.0, type='float')
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    print '[runEventSelection] analyzing %s'%opt.input 
    selectEvents(fileName=opt.input, saveProbes=opt.saveProbes, saveSummary=opt.saveSummary, outputDir=opt.outputDir, xsec=opt.xsec)


if __name__ == "__main__":
    main()
