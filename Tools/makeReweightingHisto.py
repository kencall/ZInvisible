from optparse import OptionParser

import array, sys
# Check if we want the help, if not import ROOT (otherwise ROOT overwrites the help)
if '-h' not in sys.argv and '--help' not in sys.argv:
    from ROOT import TH1D, TMath, TFile, TCanvas, TF1

from math import sqrt

############################
##  Some utilities first  ##
############################

c = TCanvas("c1", "c1", 800, 800)

def rebin1D(h, bins):
    """Rebin histo h to bins and recompute the errors."""
    new_h = TH1D("%s_rebin"%(h.GetName()), "%s_rebin"%(h.GetName()),
                 len(bins)-1, array.array('d', bins))
    new_h.Sumw2()
    for i in xrange(h.GetNbinsX()):
        bin_to_fill = new_h.FindBin(h.GetBinCenter(i+1))
        new_h.SetBinContent(bin_to_fill, h.GetBinContent(i+1) + new_h.GetBinContent(bin_to_fill))
        new_h.SetBinError(bin_to_fill, TMath.Sqrt(h.GetBinError(i+1)*h.GetBinError(i+1) + new_h.GetBinError(bin_to_fill)*new_h.GetBinError(bin_to_fill) ) )
    return new_h

def add(hs):
    """Add all histograms in list hs and return output."""
    for h in hs[1:]:
        hs[0].Add(h)
    return hs[0]

def makeRatio(h1, h2s, bins=None, newname="test"):
    """Make the data/stack ratio for the new bins."""
    # Add stack together
    h2 = add(h2s)
    if bins != None:
        # Rebin the histograms
        h1 = rebin1D(h1, bins)
        h2 = rebin1D(h2, bins)
    # Make the ratio
    h1.Divide(h2)
    h1.SetName(newname)
    h1.SetTitle(newname)
    return h1

def reweight(h, hsf):
    """Reweight histogram h according to histogram hsf: h_i * hsf_i, for each bin i"""
    new_h = h.Clone()
    for i in xrange(new_h.GetNbinsX()):
        # Get the scale factor
        sfbin = hsf.FindBin(new_h.GetBinCenter(i+1))
        sf = hsf.GetBinContent(sfbin)
        sf_e = hsf.GetBinError(sfbin)
        if sf == 0: continue
        # Get the old bin info
        bc_old = new_h.GetBinContent(i+1)
        be_old = new_h.GetBinError(i+1)
        if bc_old == 0: continue
        # Get the new bin content
        bc = bc_old*sf
        # Get the new bin error, just use basic error propagation without correlations: [E(AB)/(AB)]^2 = [E(A)/A]^2 + [E(B)/B]^2
        be = bc * TMath.Sqrt( (be_old/bc_old)*(be_old/bc_old) + (sf_e/sf)*(sf_e/sf) )
        # Write the new info
        new_h.SetBinContent(i+1,bc)
        new_h.SetBinError(i+1,be)
    return new_h

def subtract(h, hlist):
    """Subtract all histograms in hlist from h"""
    new_h = h.Clone()
    new_h.Sumw2()
    for hl in hlist:
        new_h.Add(hl,-1.)
    return new_h

##############################
##  Main scale factor code  ##
##############################

##################
## Njet weights ##
##################
def njetWeights(filename):
    # Get the file
    f = TFile.Open(filename)

    # Prepare writing to a file
    fout = TFile.Open("dataMCweights.root","RECREATE")
    fout.cd()

    # How to rebin
    bins = [0,1,2,3,4,5,6,7,8,20]
    bins_TT = [0,1,2,3,4,5,6,7,20]

    # Run over the relevant histograms
    cuts_DY = ["muZinv", "muZinv_0b", "muZinv_g1b"]
    cuts_TT = ["elmuZinv", "elmuZinv_0b", "elmuZinv_g1b"]
    selection = "ht200_dphi"
    # histo names
    hname1 = "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDatadata"
    hnames2 = ["cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDYstack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDY HT<100stack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvt#bar{t}stack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvsingle topstack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24Zinvt#bar{t}Zstack",
               "cntNJetsPt30Eta24Zinv/DataMC_SingleMuon_nj_%(cut)s_%(selection)scntNJetsPt30Eta24ZinvcntNJetsPt30Eta24ZinvDibosonstack"
               ]

    # dictionary to keep track of all the scale factors
    SFs = {}

    # First get TTbar reweighting, region is very pure, just take ratio of data/full stack and apply that to ttbar later
    for cut in cuts_TT:
        hname1_TT = hname1 % {"cut":cut, "selection":selection}
        hnames2_TT = [elem % {"cut":cut, "selection":selection} for elem in hnames2]
        # Get all histos
        h1 = f.Get(hname1_TT)
        h2s = [f.Get(hname2_TT) for hname2_TT in hnames2_TT]
        newname = "DataMC_nj_%s_%s"%(cut,selection)
        # Make the ratio
        newh = makeRatio(h1, h2s, bins_TT, newname)
        SFs["TT_%s"%(cut)] = newh
        newh.Write()

    # Now get DY reweighting:
    # Procedure: 1. Apply TTbar reweighting to TTbar MC
    #            2. Subtract non-DY MC from data
    #            3. Make ratio of subtracted data and DY
    for cut in cuts_DY:
        hname1_DY = hname1 % {"cut":cut, "selection":selection}
        hnames2_DY = [elem % {"cut":cut, "selection":selection} for elem in hnames2]
        # Get all histos
        h1 = f.Get(hname1_DY)
        h2s = [f.Get(hname2_DY) for hname2_DY in hnames2_DY]

        # apply weights to ttbar
        h2s[2] = reweight(h2s[2], SFs["TT_%s"%(cut.replace("mu","elmu"))])

        # subtract relevant histograms from data
        data_subtracted = subtract(h1, h2s[2:])

        newname = "DataMC_nj_%s_%s"%(cut,selection)
        newh = makeRatio(data_subtracted, h2s[:1], bins, newname)

        newh.Write()

    # Close files
    fout.Close()
    f.Close()

##################
## norm weights ##
##################
def normWeight(filename):
    # Get the file
    f = TFile.Open(filename)
    # Run over the relevant histograms
    cuts_DY = ["muZinv_0b"]
    selection = "blnotag"
    # histo names
    hname1 = "cntCSVSZinv/DataMC_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvDatadata"
    hnames2 = ["cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvDYstack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvDY HT<100stack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvt#bar{t}stack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvsingle topstack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvt#bar{t}Zstack",
               "cntCSVSZinv/DataMCw_SingleMuon_nb_%(cut)s_%(selection)scntCSVSZinvcntCSVSZinvDibosonstack"
               ]

    # Procedure: 1. Grab the njet reweighted MC
    #            2. Subtract non-DY MC from data
    #            3. Make ratio of subtracted data and DY
    for cut in cuts_DY:
        hname1_DY = hname1 % {"cut":cut, "selection":selection}
        hnames2_DY = [elem % {"cut":cut, "selection":selection} for elem in hnames2]
        # Get all histos
        h1 = f.Get(hname1_DY)
        h2s = [f.Get(hname2_DY) for hname2_DY in hnames2_DY]

        # subtract relevant histograms from data
        data_subtracted = subtract(h1, h2s[2:])

        newname = "DataMC_nb_%s_%s"%(cut,selection)
        newh = makeRatio(data_subtracted, h2s[:1], newname=newname)
        #newh = makeRatio(h1, h2s, newname=newname)

        print "Data/MC normalization scale factor in region %s_%s: %.3f +- %.3f" % (cuts_DY[0], selection, newh.GetBinContent(1), newh.GetBinError(1))

    f.Close()

def shapeSyst(filename):
    # Get the file
    f = TFile.Open(filename)
    fout = TFile.Open("syst_shape.root", "RECREATE")
    # Run over the relevant histograms
    # histo names
    hnameData = "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sDatadata"
    hnames2 =  ["%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sDYstack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sDY HT<100stack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)st#bar{t}stack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)ssingle topstack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)st#bar{t}Zstack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sDibosonstack",
                "%(var)s/DataMCw_SingleMuon_%(name)s_muZinv_loose0_mt2%(var)s%(var)sRarestack"]

    varList = [["met", "cleanMetPt",             [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 1500], "MET" ],
               ["mt2", "best_had_brJet_MT2Zinv", [0, 50, 100, 150, 200, 250, 300, 350, 400, 1500],      "M_{T2}" ],
               ["nt",  "nTopCandSortedCntZinv",  [0, 1, 2, 8 ],                                         "N(t)" ],
               ["nb",  "cntCSVSZinv",            [0, 1, 2, 3, 8 ],                                      "N(b)" ]]

    for var in varList:
        # Procedure: 1. Grab the njet reweighted MC
        #            2. Subtract non-DY MC from data
        #            3. Make ratio of subtracted data and DY
        # Get all histos
        hData = f.Get(hnameData%{"name":var[0], "var":var[1]})
        h2s   = [f.Get(hname2%{"name":var[0], "var":var[1]}) for hname2 in hnames2]
        
        # subtract relevant histograms from data
        data_subtracted = subtract(hData, h2s[2:])

        newname = "ShapeRatio_%s"%var[0]
        newh = makeRatio(data_subtracted, h2s[:1], newname=newname, bins=var[2])

        for i in xrange(1, newh.GetNbinsX() + 1):
            print newh.GetBinContent(i)

        fit = TF1("fit_%s"%var[0], "pol1")
        newh.SetLineWidth(2)
        newh.GetXaxis().SetTitle(var[3])
        newh.SetStats(0)
        newh.SetTitle("")
        if not "nt" in var[0] and not "nb" in var[0]:
            newh.Fit(fit, "")
        else:
            newh.Draw()
        #
        #    fit.Draw("same")
        #newh.DrawCopy()
        c.Print("shapeSyst_%s.png"%var[0])
        c.Print("shapeSyst_%s.pdf"%var[0])

        fout.cd()
        newh.Write()

    f.Close()
    fout.Close()

def systHarvest():
    # Get the file
    #f = TFile.Open(filename)
    #fout = TFile.Open("syst_shape.root", "RECREATE")
    # Run over the relevant histograms
    # histo names

    # Get shape central value uncertainty
    f = TFile("systematics.root")
    hShape_MET_Nom = f.Get("nSearchBin/systWgtMET_cleanMetPtnSearchBinnSearchBinNominalsingle")
    hShape_MET_Var = f.Get("nSearchBin/systWgtMET_cleanMetPtnSearchBinnSearchBinvariedsingle")
    hShape_MT2_Nom = f.Get("nSearchBin/systWgtMT2_best_had_brJet_MT2ZinvnSearchBinnSearchBinNominalsingle")
    hShape_MT2_Var = f.Get("nSearchBin/systWgtMT2_best_had_brJet_MT2ZinvnSearchBinnSearchBinvariedsingle")
    hShape_NT_Nom  = f.Get("nSearchBin/systWgtNT_nTopCandSortedCntZinvnSearchBinnSearchBinNominalsingle")
    hShape_NT_Var  = f.Get("nSearchBin/systWgtNT_nTopCandSortedCntZinvnSearchBinnSearchBinvariedsingle")
    hShape_NB_Nom  = f.Get("nSearchBin/systWgtNB_cntCSVSZinvnSearchBinnSearchBinNominalsingle")
    hShape_NB_Var  = f.Get("nSearchBin/systWgtNB_cntCSVSZinvnSearchBinnSearchBinvariedsingle")

    hShape_MET_ratio = hShape_MET_Var.Clone(hShape_MET_Nom.GetName()+"_ratio")
    hShape_MET_ratio.Divide(hShape_MET_Nom)

    hShape_MT2_ratio = hShape_MT2_Var.Clone(hShape_MT2_Nom.GetName()+"_ratio")
    hShape_MT2_ratio.Divide(hShape_MT2_Nom)

    hShape_NT_ratio = hShape_NT_Var.Clone(hShape_NT_Nom.GetName()+"_ratio")
    hShape_NT_ratio.Divide(hShape_NT_Nom)

    hShape_NB_ratio = hShape_NB_Var.Clone(hShape_NB_Nom.GetName()+"_ratio")
    hShape_NB_ratio.Divide(hShape_NB_Nom)

    hShape_final = hShape_MET_ratio.Clone("shape_central")

    for i in xrange(1, 46):
        uncertMET = hShape_MET_ratio.GetBinContent(i) - 1 if hShape_MET_ratio.GetBinContent(i) > 0 else 0
        uncertMT2 = hShape_MT2_ratio.GetBinContent(i) - 1 if hShape_MT2_ratio.GetBinContent(i) > 0 else 0
        uncertNT = hShape_NT_ratio.GetBinContent(i) - 1 if hShape_NT_ratio.GetBinContent(i) > 0 else 0
        uncertNB = hShape_NB_ratio.GetBinContent(i) - 1 if hShape_NB_ratio.GetBinContent(i) > 0 else 0
        hShape_final.SetBinContent(i, sqrt(uncertMET**2 + uncertMT2**2 + uncertNT**2 + uncertNB**2))

    fout = TFile("syst_all.root", "RECREATE")
    hShape_final.Write()

    # Get shape stats uncertainty
    f = TFile("syst_nJetWgt.root")
    hShapeStat = f.Get("syst68Max").Clone("shape_stat")
    fout.cd()
    hShapeStat.Write()

    # Get MC stats uncertainty
    f = TFile("/uscms/home/nstrobbe/nobackup/HadronicStop/DataTest/CMSSW_7_4_8/src/ZInvisible/Tools/condor/dataplots_muon_Jan20.root")
    hMC = f.Get("nSearchBin/NJetWgt_nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nusingle")
    hMCstats = hMC.Clone("MC_stats")
    for i in xrange(1, 46):
        if hMC.GetBinContent(i) > 0.00000001:
            hMCstats.SetBinContent(i, hMC.GetBinError(i) / hMC.GetBinContent(i))
        else:
            hMCstats.SetBinContent(i, 0.0)
    fout.cd()
    hMCstats.Write()

    hClosureZ = f.Get("nSearchBin/nSearchBinnSearchBinnSearchBinZ#rightarrow#nu#nusingle")
    hClosureDY = f.Get("nSearchBin/nSearchBinnSearchBinnSearchBinDY#rightarrow#mu#mu no #mu, Z eff+accsingle")
    hClosureRatio = hClosureZ.Clone("MC_closure")
    hClosureRatio.Add(hClosureDY, -1)
    hClosureRatio.Divide(hClosureZ)
    for i in xrange(1, 46):
        hClosureRatio.SetBinContent(i, abs(hClosureRatio.GetBinContent(i)))
    fout.cd()
    hClosureRatio.Write()

    fout.Close()

if __name__ ==  "__main__":

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="Grab histogram from FILE", metavar="FILE")
    parser.add_option("--norm", dest="normweight", default=False,
                      help="Compute the normalization weight", action='store_true')
    parser.add_option("--njet", dest="njetweight", default=False,
                      help="Compute the njet weights", action='store_true')
    parser.add_option("--shape", dest="shape", default=False,
                      help="Compute the shape systematics", action='store_true')
    parser.add_option("--systHarvest", dest="systHarvest", default=False,
                      help="Harvest systematics", action='store_true')

    (options, args) = parser.parse_args()

    if options.njetweight:
        njetWeights(options.filename)
    if options.normweight:
        normWeight(options.filename)
    if options.shape:
        shapeSyst(options.filename)
    if options.systHarvest:
        systHarvest()

