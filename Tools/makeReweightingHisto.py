from optparse import OptionParser
from ROOT import TH1D, TMath, TFile
import array

############################
##  Some utilities first  ##
############################

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
    selections = "ht200_dphi"
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


if __name__ ==  "__main__":
    
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="Grab histogram from FILE", metavar="FILE")
    parser.add_option("--norm", dest="normweight", default=False,
                      help="Compute the normalization weight", action='store_true')
    parser.add_option("--njet", dest="njetweight", default=False,
                      help="Compute the njet weights", action='store_true')


    (options, args) = parser.parse_args()

    if options.njetweight:
        njetweights(options.filename)
    if options.normweight:
        normWeight(options.filename)

    
