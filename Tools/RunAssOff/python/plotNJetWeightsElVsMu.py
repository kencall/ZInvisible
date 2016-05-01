from ROOT import *
from utils import *

datamc = 'Data'
gStyle.SetOptStat(0)
f = TFile('dataMCweights.root')


histopairs = []
keys = f.GetListOfKeys()
for key in keys:
    name = key.GetName()
    print 'pairing',name,'and',name.replace('mu','el')
    if not '_muZ' in name: continue
    hMu = f.Get(name).Clone(name+'2')
    hEl = f.Get(name.replace('mu','el')).Clone(name.replace('mu','el')+'2')
    histopairs.append([hMu,hEl])

    
cGold = TCanvas('cEnchilada','cEnchilada', 800, 800)
for pair in histopairs:
    cGold.cd()
    histoStyler(pair[0],kViolet)
    histoStyler(pair[1],kGreen)
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.0)
    pad1.SetLeftMargin(0.12)
    #pad1.SetTopMargin(0.03)
    pad1.Draw()
    pad1.cd()
    pair[0].GetXaxis().SetRangeUser(0,8)
    pair[0].SetTitle('')
    pair[0].GetYaxis().SetTitle('Data/MC')
    pair[0].GetYaxis().SetRangeUser(0.1,1.4)
    pair[1].GetXaxis().SetRangeUser(0,8)
    pair[0].Draw('E1')
    pair[1].Draw('same E1')
    print 'painting',pair[0],'and',pair[1]
    tl.SetTextFont(cmsTextFont)
    tl.SetTextSize(1.12*tl.GetTextSize())
    tl.DrawLatex(0.16,0.82, 'CMS')
    tl.SetTextFont(extraTextFont)
    tl.SetTextSize(1.0/1.12*tl.GetTextSize())
    tl.DrawLatex(0.245,0.82, ('MC' in datamc)*' simulation '+'preliminary')
    tl.SetTextFont(regularfont)
    cutlabel = mkCutsLabel(pair[0].GetName())
    tl.DrawLatex(0.16,0.04, cutlabel)
    l = mklegend(x1=.15, y1=.28, x2=.8, y2=.42)
    l.AddEntry(pair[0],'Z#rightarrow#mu #mu')
    l.AddEntry(pair[1],'Z#rightarrowe e')
    l.Draw()
    pad1.Update()
    cGold.Update()
    cGold.cd()
    
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.0)
    pad2.SetBottomMargin(0.37)
    pad2.SetLeftMargin(0.12)
    pad2.SetGridx()
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()
    hRatio = pair[1].Clone('hRatio')
    hDen = pair[0].Clone('hDen')
    hRatio.SetMarkerStyle(20)
    hRatio.SetMarkerColor(1)
    histoStyler(hRatio, 1)
    histoStyler(hDen, 1)
    hRatio.Divide(hDen)
    hRatio.GetYaxis().SetRangeUser(0.1,1.9)
    hRatio.SetTitle('')
    hRatio.GetXaxis().SetTitle('n_{jet}')
    hRatio.GetYaxis().SetTitle('ee/#mu#mu')
    hRatio.GetXaxis().SetTitleSize(0.28)
    hRatio.GetXaxis().SetLabelSize(0.165)
    hRatio.GetYaxis().SetTitleSize(0.185)
    hRatio.GetYaxis().SetLabelSize(0.165)
    hRatio.GetYaxis().SetNdivisions(7)
    hRatio.GetYaxis().SetTitleOffset(0.3)
    hRatio.GetXaxis().SetTitleOffset(0.51)
    hRatio.Draw()
    cGold.Update()  
    print 'pausing here'
    pause() 

    
exit(0)



pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.SetBottomMargin(0.0)
pad1.SetLeftMargin(0.12)
pad1.SetLogy()
pad1.Draw()
pad1.cd()

l = mklegend(x1=.5, y1=.64, x2=.97, y2=.79)

hTruth.SetFillStyle(1002)
hTruth.SetFillColor(hTruth.GetLineColor())
#hTruth.SetMarkerStyle(21)
hTruth.SetTitle('')
