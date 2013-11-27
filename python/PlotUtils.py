#!/usr/bin/env python

import math
import ROOT
from ROOT import TFile, TH1F, TH2F, THStack, TCanvas, TPad, TPaveText, TLegend

"""
A wrapper to store data and MC histograms for comparison
"""
class Plot:

    def __init__(self,name):
        self.name=name
        self.mc=[]
        self.data=None
        self.garbageList=[]

    def add(self,h,title,color,isData):
        self.garbageList.append(h)
        h.SetTitle(title)
        
        if isData:
            h.SetMarkerStyle(20)
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetLineWidth(2)
            h.SetFillColor(0)
            h.SetFillStyle(0)
            self.data=h
        else:
            h.SetMarkerStyle(1)
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetLineWidth(1)
            h.SetFillColor(color)
            h.SetFillStyle(1001)
            self.mc.append(h)

    def reset(self):
        for o in self.garbageList: o.Delete()
            
    def showTable(self,firstBin=1,lastBin=-1):
        if firstBin<1: firsBin=1

        print '-------------------------------------'
        print "Process".ljust(10),
        print "Events".ljust(10)

        tot=0
        err=0
        for h in self.mc:
            maxBins=h.GetXaxis().GetNbins()
            if lastBin<1 : lastBin=maxBins
            if lastBin>maxBins : lastBin=maxBins
            if firstBin>maxBins: firstBin=maxBins
            ierr=ROOT.Double(0)
            itot=h.IntegralAndError(firstBin,lastBin,ierr)
            pname=h.GetTitle()
            print pname.ljust(10)[:n],
            print '%3.3f\\pm%3.3f'%(itot,ierr)
            tot=tot+itot
            err=err+ierr*ierr
        print '-------------------------------------'
        print 'Total'.ljust(10),
        print '%3.3f\\pm%3.3f'%(tot,math.sqrt(err))

        if self.data is None : return
        print '-------------------------------------'
        print 'Data'.ljust(10),
        print '%3.3f'%(self.data.IntegralAndError(firstBin,lastBin))
                       
            
    def show(self,outDir):
        canvas=TCanvas('c','c',500,500)
        canvas.cd()
        t1=TPad("t1","t1", 0.0, 0.20, 1.0, 1.0)
        t1.Draw()
        t1.cd()
        self.garbageList.append(t1)
        
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
            self.garbageList.append(frame)
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
            self.garbageList.append(totalMC)
            totalMC.SetDirectory(0)
            maxY=max(totalMC.GetMaximum(),maxY)
            if frame is None:
                frame=totalMC.Clone('frame')
                frame.Reset('ICE')
                self.garbageList.append(frame)

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
            self.garbageList.append(t2)
            t2.SetTopMargin(0)
            t2.SetBottomMargin(0.2)
            t2.Draw()
            t2.cd()
            ratio=self.data.Clone('ratio')
            self.garbageList.append(ratio)
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
Style options mostly from CMS's tdrStyle.C
"""
def customROOTstyle() :
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

