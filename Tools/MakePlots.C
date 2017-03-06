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
    std::string sbEra = "SB_v1_2017";//"SB_v1_2017";

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
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_600to800"]};
        fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_2500toInf"]};
        fileMap["DYJetsToLL_HT_600to800"] = {ss["DYJetsToLL_HT_600to800"]};
        fileMap["ZJetsToNuNu_HT_2500toInf"] = {ss["ZJetsToNuNu_HT_2500toInf"]};
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
        fileMap["TTbarNoHad"] = {ss["TTbarDiLep"]};
        fileMap["Data_SingleMuon"] = {ss["Data_SingleMuon"]};
    }
    else if(dataSets.compare("TEST2") == 0)
    {
        fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_600to800"]};
        fileMap["DYJetsToLL_HT_600to800"] = {ss["DYJetsToLL_HT_600to800"]};
        fileMap["IncDY"] = {ss["DYJetsToLL_Inc"]}; 
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
        fileMap["TTbarNoHad"] = {ss["TTbarDiLep"]};
        fileMap["Data_SingleMuon"] = {ss["Data_SingleMuon"]};
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
    std::string label_met = "E_{T}^{miss} [GeV]";
    //std::string label_met = "#slash{E}_{T} [GeV]";
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
    std::string label_toppt = "top p_{T} [GeV]";

    vector<Plotter::HistSummary> vh;

    // -----------------
    // - Data/MC plots -
    // -----------------

    // Datasetsummaries we are using
    // no weight (genWeight deals with negative weights); also add btag weights here
    Plotter::DatasetSummary dsData("Data",       fileMap["Data_MET"], "passSearchTrigger",   "");

    // Collections for all variables, no cuts applied yet

    // lambda functions for shorthand
    auto PDCData    = [&](std::string var){ return vector<PDC>({PDC("data", var, {dsData})});   };

    // pair of cutlevels
    std::vector<std::pair<std::string,std::string>> cutlevels_muon = {
	{"nosel",                     "passNoiseEventFilter"},
        {"baseline",         "passNoiseEventFilter;passBaseline"}
    };


    //push the histograms in a loop, save some copy-paste time
    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
	// unweighted
	vh.push_back(PHS("Data_MET_met_"        +cut.first,  PDCData("met"),                {1, 1}, cut.second, 50, 0, 1500,   true, false,  label_met,        "Events"));
	vh.push_back(PHS("Data_MET_ht_"         +cut.first,  PDCData("HT"),                 {1, 1}, cut.second, 50, 0, 1500,   true, false,  label_ht,         "Events"));
	vh.push_back(PHS("Data_MET_nt_"         +cut.first,  PDCData("nTopCandSortedCnt"),  {1, 1}, cut.second, 5,  0, 5,      true, false,  label_nt,         "Events"));
	vh.push_back(PHS("Data_MET_mt2_"        +cut.first,  PDCData("best_had_brJet_MT2"), {1, 1}, cut.second, 50, 0, 1500,   true, false,  label_mt2,        "Events"));
	vh.push_back(PHS("Data_MET_nb_"         +cut.first,  PDCData("cntCSVS"),            {1, 1}, cut.second, 10, 0, 10,     true, false,  label_nb,         "Events"));
	vh.push_back(PHS("Data_MET_nj_"         +cut.first,  PDCData("cntNJetsPt30"),       {1, 1}, cut.second, 20, 0, 20,     true, false,  label_nj,         "Events"));
	vh.push_back(PHS("Data_MET_nSearchBin_" +cut.first,  PDCData("nSearchBin"),         {1, 1}, cut.second, NSB, 0, NSB,   true, false,  "Search Bin",     "Events"));
    }

    vector<string> cfsSync = {"",
                              "passNoiseEventFilter",
                              "passNoiseEventFilter;passnJets",
                              "passNoiseEventFilter;passnJets;passMuonVeto",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets;passMET",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets;passMET;passHT",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets;passMET;passHT;passTagger",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets;passMET;passHT;passTagger;passMT2"};


    vector<Plotter::CutFlowSummary> cutFlowSummaries;

    cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("DataCutFlow",         PDC("", "", {dsData}),            cfsSync));

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    RegisterFunctions* rf = new RegisterFunctionsNTuple(runOnCondor, sbEra);

    Plotter plotter(vh, vvf, fromTuple, histFile, nFiles, startFile, nEvts);
    plotter.setCutFlows(cutFlowSummaries);
    plotter.setLumi(lumi);
    plotter.setPlotDir(plotDir);
    plotter.setDoHists(doSave || doPlots);
    plotter.setDoTuple(doTuple);
    plotter.setRegisterFunction(rf);
    plotter.read();
    if(doSave && fromTuple)  plotter.saveHists();
    if(doPlots)              plotter.plot();
}
