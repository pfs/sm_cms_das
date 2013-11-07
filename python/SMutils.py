import math
from ROOT import TFile, TLorentzVector, TH1F, TH2F, TGraph, TObjArray

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
                print self.norm
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
        self.relIso=relIso
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
