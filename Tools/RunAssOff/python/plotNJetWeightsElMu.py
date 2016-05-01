from ROOT import *
from utils import *




lep = 'emu'
lep = 'combined'
datamc = 'Data'
gStyle.SetOptStat(0)

f = TFile('dataMCweightsCombined.root')


histopairs = []
keys = f.GetListOfKeys()
for key in keys:
    name = key.GetName()
    print 'pairing',name,'and',name.replace('mu','el')
    if not '_muZ' in name: continue
    hMu = f.Get(name).Clone(name+'2')
    histopairs.append([hMu])

    
cGold = TCanvas('cEnchilada','cEnchilada', 800, 800)
for pair in histopairs:
    cGold.cd()
    histoStyler(pair[0],kRed)
    cGold.SetLeftMargin(0.145)
    cGold.SetBottomMargin(0.145)
    pair[0].GetXaxis().SetRangeUser(0,8)
    pair[0].SetTitle('')
    pair[0].GetYaxis().SetTitle('Data/MC')
    pair[0].GetYaxis().SetRangeUser(0.1,1.4)
    pair[0].Draw('E1')
    print 'painting',pair[0],'and'
    tl.SetTextFont(cmsTextFont)
    tl.SetTextSize(1.12*tl.GetTextSize())
    tl.DrawLatex(0.16,0.82, 'CMS')
    tl.SetTextFont(extraTextFont)
    tl.SetTextSize(1.0/1.12*tl.GetTextSize())
    tl.DrawLatex(0.3,0.82, ('MC' in datamc)*' simulation '+'preliminary')
    tl.SetTextFont(regularfont)
    cutlabel = mkCutsLabel(pair[0].GetName())
    tl.DrawLatex(0.18,0.17, cutlabel)
    l = mklegend(x1=.15, y1=.4, x2=.1, y2=.4)
    l.AddEntry(pair[0],'Z#rightarrow#mu #mu')
    l.Draw()
    cGold.Update()
    pause()

    
exit(0)
