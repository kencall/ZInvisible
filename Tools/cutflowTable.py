from optparse import OptionParser
import ROOT as rt

# Define some useful latex things here
def header():
    return """\\documentclass[10pt]{article}

\\usepackage[margin=2cm]{geometry}

\\begin{document}

"""

def footer():
    return """
\\end{document}
"""

def tableheader(ncols, alignment='c', divisions=False):
    head = " | ".join(['c']*ncols) if divisions else "".join(['c']*ncols)
    return """\\begin{tabular}{%s}
""" % (head)

def tablefooter():
    return """\\end{tabular}
"""

def tablebody(table, cutlevels):
    # titles
    coltitles = sorted(table.keys())
    c = ["Cut & %s" % (" & ".join(coltitles))]
    # cutlevels
    for i,cut in enumerate(cutlevels):
        l = [cut]
        for sample in coltitles:
            l.append("%.3g"%(table[sample][i]))
        c.append(" & ".join(l))

    c.append("")
    return " \\\\ \n".join(c)


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="Grab cutflow histogram from FILE", metavar="FILE")
    parser.add_option("-o", "--out", dest="outfilename",
                      help="Output filename to store the latex cutflow", metavar="FILE")


    (options, args) = parser.parse_args()

    # Open file
    f = rt.TFile.Open(options.filename)

    print f

    # Get the relevant histograms
    hnames = ["CutFlows/Data",
              "CutFlows/MC",
              "CutFlows/Rare",
              "CutFlows/QCD",
              "CutFlows/Single top",
             # "CutFlows/Signal Points",
              "CutFlows/t#bar{t}Z",
              "CutFlows/t#bar{t}",
              "CutFlows/WJets",
              "CutFlows/ZJets",
              ]

    samples = ["Data HTMHT",
               "Total MC",
               "Rare",
               "QCD",
               "Single top",
               "t#bar{t}Z",
               "t#bar{t}",
               "WJets",
               "ZJets",
               #"Z$\\rightarrow\\nu\\nu$",
               #"Z$\\rightarrow\\nu\\nu$",
              ] 

    cuts = [" ",
            "Noise filters",
            "Njets",
            "Muon Veto",
            "Electron Veto",
            "IsoTrkVeto",
            "dPhi",
            "B Jets",
            "$\\textrm{MET}>200$",
            "$\\textrm{HT}>500$",
            "$\\textrm{Nt}>0$",
            "$\\textrm{MT2}>200$",
            "Baseline"]


    hs = [f.Get(hname) for hname in hnames]

    print hs

    # extract the info from the histograms
    info = {}
    for i, h in enumerate(hs):
        cutinfo = [h.GetBinContent(bin+1) for bin in xrange(h.GetNbinsX())]
        info[samples[i]] = cutinfo

    # open outputfile
    out = open(options.outfilename, 'w')

    # Write the actual tables
    out.write(header())
    out.write(tableheader(7, divisions=True))
    out.write(tablebody(info, cuts).replace("_","\_"))
    out.write(tablefooter())
    out.write(footer())

    # Close files
    out.close()
    f.Close()
