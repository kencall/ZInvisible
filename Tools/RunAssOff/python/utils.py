from ROOT import *
import sys
from array import array
import numpy as np
import random

tl = TLatex()
tl.SetNDC()
cmsTextFont = 61
extraTextFont = 52
lumiTextSize = 0.6
lumiTextOffset = 0.2
cmsTextSize = 0.75
cmsTextOffset = 0.1
regularfont = 42
originalfont = tl.GetTextFont()



def histoStyler(h,color):
    h.SetLineWidth(2)
    h.SetLineColor(color)
    size = 0.065
    font = 132
    h.GetXaxis().SetLabelFont(font)
    h.GetYaxis().SetLabelFont(font)
    h.GetXaxis().SetTitleFont(font)
    h.GetYaxis().SetTitleFont(font)
    h.GetYaxis().SetTitleSize(size)
    h.GetXaxis().SetTitleSize(size)
    h.GetXaxis().SetLabelSize(size)   
    h.GetYaxis().SetLabelSize(size)
    h.GetXaxis().SetTitleOffset(1.0)
    h.GetYaxis().SetTitleOffset(0.9)
    h.Sumw2()


def histoStyler2d(h): 
    h.GetYaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetTitleOffset(1.0)
    size = 0.065
    font = 132
    h.GetXaxis().SetLabelFont(font)
    h.GetYaxis().SetLabelFont(font)
    h.GetXaxis().SetTitleFont(font)
    h.GetYaxis().SetTitleFont(font)
    h.GetYaxis().SetTitleSize(size)
    h.GetXaxis().SetTitleSize(size)
    h.GetXaxis().SetLabelSize(size)   
    h.GetYaxis().SetLabelSize(size)
    h.GetXaxis().SetTitleOffset(1.0)
    h.GetYaxis().SetTitleOffset(1.12)
    h.Sumw2()
    return h







units = {}
units['Mht']='Gev'
units['Met']=units['Mht']
units['Ht']='GeV'
units['NJets']='bin'
units['BTags']='bin'
units['Jet1Pt']='GeV'
units['Jet1Eta']='bin'
units['Jet2Pt']='GeV'
units['Jet2Eta']='bin'
units['Jet3Pt']='GeV'
units['Jet3Eta']='bin'
units['Jet4Pt']='GeV'
units['Jet4Eta']='bin'
units['MhtPhi']='rad'
units['DPhiJet1Mht']='rad'
units['DPhiJet2Mht']='rad'
units['DPhiJet3Mht']='rad'
units['DPhiJet4Mht']='rad'




def redoMET(originalMet,originalJets,newJets):
    if not len(originalJets)==len(newJets):
        print 'mismatch'; exit()
    newMET = originalMet.Clone()
    for jet in originalJets:
        newMET+=jet
    for jet in newJets:
        newMET-=jet
    return [newMET.Pt(), newMET.Phi()]



def getHT(tlvVec,thresh):
    ht = 0
    for tlv in tlvVec:
        if not (abs(tlv.Eta())<2.4): continue
        if not (tlv.Pt()>thresh): continue
        ht+=tlv.Pt()
    return ht

def mklegend(x1=.1105, y1=.53, x2=.3805, y2=.8, color=kWhite):
    lg = TLegend(x1, y1, x2, y2)
    lg.SetFillColor(color)
    lg.SetTextFont(42)
    lg.SetBorderSize(0)
    lg.SetShadowColor(kWhite)
    lg.SetFillStyle(0)
    return lg



def pause(thingy='please push enter'):
    import sys
    print thingy
    sys.stdout.flush()
    raw_input('')

def mkmet(metPt, metPhi):
    met = TLorentzVector()
    met.SetPtEtaPhiE(metPt, 0, metPhi, metPt)
    return met
    

def mkCutsLabel(title):
    archive = ['ht200','dphi']
    
    str_ = ''
    if 'ht200' in title: str_+= 'H_{T}>200 GeV, '
    if 'dphi' in title: str_+= '#Delta#phi cuts, '
    if '0b' in title: str_+= 'b-tags 0, '
    if 'g1b' in title: str_+= 'b-tags #geq 1, '
    if len(str_)>1:
        if str_[-2:]==', ':
            str_=str_[:-2]
    return str_
