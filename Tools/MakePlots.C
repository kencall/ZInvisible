#include "Plotter.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "SusyAnaTools/Tools/samples.h"

#include <getopt.h>
#include <iostream>

int main(int argc, char* argv[])
{
    using namespace std;

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"plot",             no_argument, 0, 'p'},
        {"savehist",         no_argument, 0, 's'},
        {"savetuple",        no_argument, 0, 't'},
        {"fromFile",         no_argument, 0, 'f'},
        {"condor",           no_argument, 0, 'c'},
        {"histFile",   required_argument, 0, 'H'},
        {"dataSets",   required_argument, 0, 'D'},
        {"numFiles",   required_argument, 0, 'N'},
        {"startFile",  required_argument, 0, 'M'},
        {"numEvts",    required_argument, 0, 'E'},
        {"plotDir",    required_argument, 0, 'P'},
	{"luminosity", required_argument, 0, 'L'},
	{"sbEra",      required_argument, 0, 'S'}
    };

    bool doPlots = true, doSave = true, doTuple = true, fromTuple = true, runOnCondor = false;
    string histFile = "", dataSets = "", sampleloc = AnaSamples::fileDir, plotDir = "plots";
    int nFiles = -1, startFile = 0, nEvts = -1;
    double lumi = AnaSamples::luminosity;
    std::string sbEra = "SB_69_2016";

    while((opt = getopt_long(argc, argv, "pstfcH:D:N:M:E:P:L:S:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'p':
            if(doPlots) doSave  = doTuple = false;
            else        doPlots = true;
            break;

        case 's':
            if(doSave) doPlots = doTuple = false;
            else       doSave  = true;
            break;

        case 't':
            if(doTuple) doPlots = doSave = false;
            else        doTuple  = true;
            break;

        case 'f':
            fromTuple = false;
            break;

        case 'c':
            runOnCondor = true;
            break;

        case 'H':
            histFile = optarg;
            break;

        case 'D':
            dataSets = optarg;
            break;

        case 'N':
            nFiles = int(atoi(optarg));
            break;

        case 'M':
            startFile = int(atoi(optarg));
            break;

        case 'E':
            nEvts = int(atoi(optarg));
            break;

        case 'P':
            plotDir = optarg;
            break;

	case 'L':
            lumi = atof(optarg);
            break;

        case 'S':
            sbEra = optarg;
            break;
        }
    }

    //if running on condor override all optional settings
    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "histoutput_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
        //doSave = true;
        //doPlots = false;
        fromTuple = true;
        sampleloc = "condor";
    }

    AnaSamples::SampleSet        ss(sampleloc, lumi);
    AnaSamples::SampleCollection sc(ss);

    const double zAcc = 1.0;
//    const double zAcc = 0.5954;
//    const double zAcc = 0.855;
    const double znunu_mumu_ratio = 5.942;
    const double znunu_ee_ratio   = 5.942;

    map<string, vector<AnaSamples::FileSummary>> fileMap;

    //Select approperiate datasets here
    if(dataSets.compare("TEST") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_2500toInf"]};
        fileMap["DYJetsToLL_HT_600toInf"] = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["ZJetsToNuNu_HT_2500toInf"] = {ss["ZJetsToNuNu_HT_2500toInf"]};
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
        fileMap["TTbarNoHad"] = {ss["TTbarDiLep"]};
        fileMap["Data_SingleMuon"] = {ss["Data_SingleMuon_2016"]};
    }
    else if(dataSets.compare("TEST2") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["DYJetsToLL_HT_600toInf"] = {ss["DYJetsToLL_HT_600toInf"]};
        fileMap["IncDY"] = {ss["DYJetsToLL_Inc"]}; 
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
        fileMap["TTbarNoHad"] = {ss["TTbarDiLep"]};
        fileMap["Data_SingleMuon"] = {ss["Data_SingleMuon_2015C"]};
    }
    else
    {
        if(ss[dataSets] != ss.null())
        {
            fileMap[dataSets] = {ss[dataSets]};
            for(const auto& colls : ss[dataSets].getCollections())
            {
                fileMap[colls] = {ss[dataSets]};
            }
        }
        else if(sc[dataSets] != sc.null())
        {
            fileMap[dataSets] = {sc[dataSets]};
            int i = 0;
            for(const auto& fs : sc[dataSets])
            {
                fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
            }
        }
    }

    // Number of searchbins
    SearchBins sb(sbEra);
    int NSB = sb.nSearchBins();//37; // 45

    // Shortcuts for axis labels
    //std::string label_met = "p_{T}^{miss} [GeV]";
    std::string label_met = "#slash{E}_{T} [GeV]";
    std::string label_ht  = "H_{T} [GeV]";
    std::string label_mht = "MH_{T} [GeV]";
    std::string label_nj  = "N_{jets}";
    std::string label_nb  = "N_{b}";
    std::string label_nt  = "N_{t}";
    std::string label_mt2 = "M_{T2} [GeV]";
    std::string label_mupt  = "#mu p_{T} [GeV]";
    std::string label_mu1pt = "#mu_{1} p_{T} [GeV]";
    std::string label_mu2pt = "#mu_{2} p_{T} [GeV]";
    std::string label_elpt  = "e p_{T} [GeV]";
    std::string label_el1pt = "e_{1} p_{T} [GeV]";
    std::string label_el2pt = "e_{2} p_{T} [GeV]";
    std::string label_jpt  = "j p_{T} [GeV]";
    std::string label_j1pt = "j_{1} p_{T} [GeV]";
    std::string label_j2pt = "j_{2} p_{T} [GeV]";
    std::string label_j3pt = "j_{3} p_{T} [GeV]";
    std::string label_mll  = "m_{ll} [GeV]";


    vector<Plotter::HistSummary> vh;

    // -----------------
    // - Data/MC plots -
    // -----------------

    // Datasetsummaries we are using
    // no weight (genWeight deals with negative weights)
    //                                        legend label  sample commection           cuts               weights
    Plotter::DatasetSummary dsData_SingleMuon("Data",       fileMap["Data_SingleMuon"], "passMuTrigger",   "");
    Plotter::DatasetSummary dsData_DoubleEG(  "Data",       fileMap["Data_DoubleEG"],   "passElecTrigger", "");
    Plotter::DatasetSummary dsDY(             "DY",         fileMap["DYJetsToLL"],      "",                "");
    Plotter::DatasetSummary dsDYInc(          "DY HT<100",  fileMap["IncDY"],           "genHT<100",       "");
    Plotter::DatasetSummary dstt2l(           "t#bar{t}",   fileMap["TTbarNoHad"],      "",                "");
    Plotter::DatasetSummary dstW(             "Single top", fileMap["tW"],              "",                "");
    Plotter::DatasetSummary dsttZ(            "t#bar{t}Z",  fileMap["TTZ"],             "",                "genWeight");
    Plotter::DatasetSummary dsVV(             "Diboson",    fileMap["Diboson"],         "",                "");
    Plotter::DatasetSummary dsRare(           "Rare",       fileMap["Rare"],            "",                "genWeight");
    std::vector<std::vector<Plotter::DatasetSummary>> stack_MC = {{dsDY, dsDYInc}, {dstt2l}, {dstW}, {dsttZ}, {dsVV}, {dsRare}};
    // Collections for all variables, no cuts applied yet
    // met

    auto PDCMaker = [&](std::string var) 
    {
      //                                                   plotType  var  vector of datasetsummary
	return  std::vector<PDC>({ Plotter::DataCollection("data",   var, {dsData_SingleMuon}),
		                   Plotter::DataCollection("stack",  var, stack_MC)             });
    };

    // pair of cutlevels
    std::vector<std::pair<std::string,std::string>> cutlevels_muon = {
        //cut name                    cut string
	{"nosel",                     "passNoiseEventFilterZinv"},
	{"2mu",                       "passNoiseEventFilterZinv;passDiMuSel"},
	{"muZinv",                    "passNoiseEventFilterZinv;passMuZinvSel"},
	{"muZinv_ht200",              "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200"},
	{"muZinv_ht200_dphi",         "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passdPhisZinv"},
	{"muZinv_ht50_met50_dphi",    "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
	{"muZinv_loose0",             "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
	{"muZinv_loose0_mt2",         "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>0"},
	{"muZinv_loose0_mt21b",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv1b>0"},
	{"muZinv_loose0_mt22b",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv2b>0"},
	{"muZinv_loose0_mt23b",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv3b>0"},
	{"muZinv_loose0_ntop",        "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv>0"},
	{"muZinv_loose0_nb2",         "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv>=2"},
	{"muZinv_loose0_ht300",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>300;passnJetsZinv;passdPhisZinv"},
	{"muZinv_loose0_ht400",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>400;passnJetsZinv;passdPhisZinv"},
	{"muZinv_loose0_ht500",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>500;passnJetsZinv;passdPhisZinv"},
	{"muZinv_blnotagmt2",         "passMuZinvSel;passBaselineNoTagMT2Zinv"},
	{"muZinv_blnotag",            "passMuZinvSel;passBaselineNoTagZinv"},
	{"muZinv_bl",                 "passMuZinvSel;passBaselineZinv"},
	{"muZinv_0b",                 "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0"},
	{"muZinv_0b_ht200",           "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>200"},
	{"muZinv_0b_ht200_dphi",      "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>200;passdPhisZinv"},
	{"muZinv_0b_ht50_met50_dphi", "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
	{"muZinv_0b_loose0",          "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv"},
	{"muZinv_0b_loose0_ht300",    "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>300;passnJetsZinv;passdPhisZinv"},
	{"muZinv_0b_loose0_ht400",    "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>400;passnJetsZinv;passdPhisZinv"},
	{"muZinv_0b_loose0_ht500",    "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>500;passnJetsZinv;passdPhisZinv"},
	{"muZinv_0b_blnotagmt2",      "passMuZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv"},
	{"muZinv_0b_blnotag",         "passMuZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv"},
	{"muZinv_g1b",                "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv"},
	{"muZinv_g1b_ht200",          "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>200"},
	{"muZinv_g1b_ht200_dphi",     "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>200;passdPhisZinv"},
	{"muZinv_g1b_ht50_met50_dphi","passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
	{"muZinv_g1b_loose0",         "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv"},
	{"muZinv_g1b_loose0_ht300",   "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>300;passnJetsZinv;passdPhisZinv"},
	{"muZinv_g1b_loose0_ht400",   "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>400;passnJetsZinv;passdPhisZinv"},
	{"muZinv_g1b_loose0_ht500",   "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>500;passnJetsZinv;passdPhisZinv"},
	{"elmu",                      "passNoiseEventFilterZinv;passElMuSel"},
	{"elmuZinv",                  "passNoiseEventFilterZinv;passElMuZinvSel"},
	{"elmuZinv_ht200",            "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200"},
	{"elmuZinv_ht200_dphi",       "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200;passdPhisZinv"},
	{"elmuZinv_loose0",           "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
	{"elmuZinv_blnotagmt2",       "passElMuZinvSel;passBaselineNoTagMT2Zinv"},
	{"elmuZinv_bl",               "passElMuZinvSel;passBaselineZinv"},
	{"elmuZinv_0b",               "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0"},
	{"elmuZinv_0b_ht200",         "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200"},
	{"elmuZinv_0b_ht200_dphi",    "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passdPhisZinv"},
	{"elmuZinv_0b_loose0",        "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv"},
	{"elmuZinv_0b_blnotagmt2",    "passElMuZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv"},
	{"elmuZinv_0b_blnotag",       "passElMuZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv"},
	{"elmuZinv_g1b",              "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv"},
	{"elmuZinv_g1b_ht200",        "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200"},
	{"elmuZinv_g1b_ht200_dphi",   "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200;passdPhisZinv"},
	{"elmuZinv_g1b_loose0",       "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv"}
    };


    //push the histograms in a loop, save some copy-paste time
    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
	// unweighted
        //               histogram name                         vector of PDC           ratio  cutstring    bins ll ul   log   norm    x-axis      y-axis     
	vh.push_back(PHS("DataMC_SingleMuon_met_" + cut.first,  PDCMaker("cleamMetPt"), {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met,  "Events"));
	vh.push_back(PHS("DataMC_SingleMuon_ht_"  + cut.first,  PDCMaker("HTZinv"),     {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,   "Events"));
    }

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsNTuple(runOnCondor, sbEra);

    Plotter plotter(vh, vvf, fromTuple, histFile, nFiles, startFile, nEvts);
    //plotter.setCutFlows(cutFlowSummaries);
    plotter.setLumi(lumi);
    plotter.setPlotDir(plotDir);
    plotter.setDoHists(doSave || doPlots);
    plotter.setDoTuple(doTuple);
    plotter.setRegisterFunction(rf);
    plotter.read();
    if(doSave && fromTuple)  plotter.saveHists();
    if(doPlots)              plotter.plot();
}
