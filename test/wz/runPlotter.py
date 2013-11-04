#!/usr/bin/env python
import os,sys
import json
import optparse
import commands
import math
import ROOT
from ROOT import TFile, TH1F, TH2F, THStack, TCanvas, TPad, TPaveText, TLegend

"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal


class Plot:
    def __init__(self,name):
        self.name=name
        self.mc=[]
        self.data=None
    def add(self,h,title,color,isData):
        h.SetTitle(title)
        h.SetDirectory(0)
        if isData:
            h.SetName('%sdata'%self.name)
            h.SetMarkerStyle(20)
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetLineWidth(2)
            h.SetFillColor(0)
            h.SetFillStyle(0)
            self.data=h
        else:
            h.SetName('%s%d'%(self.name,len(self.mc)))
            h.SetMarkerStyle(1)
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetLineWidth(1)
            h.SetFillColor(color)
            h.SetFillStyle(1001)
            self.mc.append(h)
    def show(self,outDir):
        canvas=TCanvas('c','c',500,500)
        canvas.cd()
        t1=TPad("t1","t1", 0.0, 0.20, 1.0, 1.0)
        t1.Draw()
        t1.cd()

        frame=None
        leg=TLegend(0.15,0.9,0.9,0.95)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        nlegCols=0
                
        maxY=1.0
        if self.data is not None:
            leg.AddEntry( self.data, self.data.GetTitle(),'p')
            frame=self.data.Clone('frame')
            maxY=self.data.GetMaximum()*1.1
            frame.Reset('ICE')
            
        stack=THStack('mc','mc')
        for h in self.mc :
            stack.Add(h,'hist')
            leg.AddEntry(h,h.GetTitle(),'f')
            nlegCols=nlegCols+1
            
        totalMC=None
        if nlegCols>0:
            totalMC=stack.GetStack().At( stack.GetStack().GetEntriesFast()-1 ).Clone('totalmc')
            totalMC.SetDirectory(0)
            maxY=max(totalMC.GetMaximum(),maxY)
            if frame is None:
                frame=totalMC.Clone('frame')
                frame.Reset('ICE')

        if self.data is not None: nlegCols=nlegCols+1

        if nlegCols==0:
            print '%s is empty'%self.name
            return

        frame.GetYaxis().SetRangeUser(1e-2,maxY)
        frame.SetDirectory(0)
        frame.Draw()
        frame.GetYaxis().SetTitleOffset(1.6)
        stack.Draw('hist same')
        if self.data is not None: self.data.Draw('same')
        leg.SetNColumns(nlegCols)
        leg.Draw()
        pt=TPaveText(0.12,0.95,0.9,0.99,'brNDC')
        pt.SetBorderSize(0)
        pt.SetFillStyle(0)
        pt.SetTextAlign(12)
        pt.AddText('CMS preliminary, #sqrt{s}=8 TeV')
        pt.Draw()

        if totalMC is None or self.data is None:
            t1.SetPad(0,0,1,1)
        else :
            canvas.cd()
            t2=TPad("t2","t2", 0.0, 0.0, 1.0, 0.2)
            t2.SetTopMargin(0)
            t2.SetBottomMargin(0.2)
            t2.Draw()
            t2.cd()
            
            ratio=self.data.Clone('ratio')
            ratio.Divide(totalMC)
            ratio.SetDirectory(0)
            ratio.Draw('e1')
            ratio.GetYaxis().SetRangeUser(0.42,1.38)
            ratio.GetYaxis().SetTitle('Data/#SigmaBkg')
            ratio.GetXaxis().SetTitle('')
            ratio.GetYaxis().SetNdivisions(5)
            ratio.GetYaxis().SetTitleOffset(0.5)
            ratio.GetYaxis().SetLabelSize(0.12)
            ratio.GetYaxis().SetTitleSize(0.12)
            ratio.GetXaxis().SetLabelSize(0.12)
            ratio.GetXaxis().SetTitleSize(0.12)
            

        canvas.cd()
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs(outDir+'/'+self.name+'.png')


"""
Loop over file to find all histograms
"""
def getAllPlotsFrom(dir):
    toReturn=[]
    allKeys=dir.GetListOfKeys()
    for tk in allKeys:
        k=tk.GetName()
        obj=dir.Get(k)
        if obj.InheritsFrom('TDirectory') :
            allKeysInSubdir = getAllPlotsFrom(obj)
            for kk in allKeysInSubdir : toReturn.append( k+'/'+kk )
        elif obj.InheritsFrom('TH1') :
            toReturn.append( k )
    return toReturn

        

"""
Loop over the inputs and launch jobs
"""
def runPlotter(inDir, jsonUrl, lumi ):

    jsonFile = open(jsonUrl,'r')
    procList=json.load(jsonFile,encoding='utf-8').items()

    #make a survey of *all* existing plots
    plots=[]
    for proc in procList :
    
        for desc in proc[1] :

            isData=getByLabel(desc,'isdata',False)
            data = desc['data']
            for d in data :
                dtag = getByLabel(d,'dtag','')
                split=getByLabel(d,'split',1)
            
                for segment in range(0,split) :
                    eventsFile=dtag
                    if split>1:
                        eventsFile=dtag + '_' + str(segment)
                    rootFileUrl=inDir+'/'+eventsFile+'.root'
                    if(rootFileUrl.find('/store/')==0)  :
                        rootFileUrl = commands.getstatusoutput('cmsPfn ' + rootFileUrl)[1]
                    
                    rootFile=TFile.Open(rootFileUrl)
                    try:
                        if rootFile.IsZombie() : continue
                    except:
                        continue
                    iplots=getAllPlotsFrom(dir=rootFile)
                    rootFile.Close()
                    plots=list(set(plots+iplots))

    #now plot them
    for p in plots:

        pName=p.replace('/','')
        newPlot=Plot(pName)
        
        for proc in procList :
            for desc in proc[1] :
                title=getByLabel(desc,'tag','unknown')
                isData=getByLabel(desc,'isdata',False)
                color=int(getByLabel(desc,'color',1))
                data = desc['data']
                h=None
                for d in data :
                    dtag = getByLabel(d,'dtag','')
                    split=getByLabel(d,'split',1)

                    for segment in range(0,split) :
                        eventsFile=dtag
                        if split>1:
                            eventsFile=dtag + '_' + str(segment)
                            rootFileUrl=inDir+'/'+eventsFile+'.root'
                        if(rootFileUrl.find('/store/')==0)  :
                            rootFileUrl = commands.getstatusoutput('cmsPfn ' + rootFileUrl)[1]
                                
                    rootFile=TFile.Open(rootFileUrl)
                    try:
                        if rootFile.IsZombie() : continue
                    except:
                        continue
                    if h is None :
                        h=rootFile.Get(p)
                        try:
                            h.SetDirectory(0)
                        except:
                            h=None
                            continue
                    else:
                        h.Add(rootFile.Get(p))
                    rootFile.Close()
                if h is None: continue
                if not isData: h.Scale(lumi)
                newPlot.add(h,title,color,isData)
        newPlot.show(inDir+'/plots')
        
                    
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in'         ,    dest='inDir'              , help='Input directory'                        , default=None)
    parser.add_option('-j', '--json'       ,    dest='json'               , help='A json file with the samples to analyze', default=None)
    parser.add_option('-l', '--lumi'       ,    dest='lumi'               , help='Re-scale to integrated luminosity [pb]',  default=1.0, type='float')
    (opt, args) = parser.parse_args()

    if opt.inDir is None or opt.json is None:
        parser.print_help()
        sys.exit(1)

    #style options (mostly from tdrStyle.C
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gStyle.SetPadTopMargin(0.1);
    ROOT.gStyle.SetPadBottomMargin(0.13);
    ROOT.gStyle.SetPadLeftMargin(0.15);
    ROOT.gStyle.SetPadRightMargin(0.02);
    ROOT.gStyle.SetLabelColor(1, "XYZ");
    ROOT.gStyle.SetLabelFont(42, "XYZ");
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
    ROOT.gStyle.SetLabelSize(0.05, "XYZ");
    ROOT.gStyle.SetAxisColor(1, "XYZ");
    ROOT.gStyle.SetStripDecimals(True);
    ROOT.gStyle.SetTickLength(0.03, "XYZ");
    ROOT.gStyle.SetNdivisions(510, "XYZ");
    ROOT.gStyle.SetPadTickX(0);
    ROOT.gStyle.SetPadTickY(0);
    ROOT.gStyle.SetMarkerStyle(20);
    ROOT.gStyle.SetHistLineColor(1);
    ROOT.gStyle.SetHistLineStyle(0);
    ROOT.gStyle.SetHistLineWidth(1);
    ROOT.gStyle.SetFrameBorderMode(0);
    ROOT.gStyle.SetFrameBorderSize(1);
    ROOT.gStyle.SetFrameFillColor(0);
    ROOT.gStyle.SetFrameFillStyle(0);
    ROOT.gStyle.SetFrameLineColor(1);
    ROOT.gStyle.SetFrameLineStyle(1);
    ROOT.gStyle.SetFrameLineWidth(1);

  
    os.system('mkdir -p %s/plots'%opt.inDir)
    runPlotter(inDir=opt.inDir, jsonUrl=opt.json, lumi=opt.lumi)
    print 'Plots have been saved to %s/plots'%opt.inDir

if __name__ == "__main__":
    main()
