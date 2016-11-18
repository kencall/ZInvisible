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
        //fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_2500toInf"]};
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
    std::string label_nw  = "N_{W}";
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
    // no weight (genWeight deals with negative weights)
    //                                        legend label  sample commection           cuts               weights
    Plotter::DatasetSummary dsData_SingleMuon("Data",         fileMap["Data_HTMHT"], "passSingleMuTrigger",   "");
    Plotter::DatasetSummary dsData_SingleMuonNotrig("Data",   fileMap["Data_HTMHT"], "",   "");
    Plotter::DatasetSummary dsData_HTMHT(  "Data",            fileMap["Data_HTMHT"],   "passhtmhtTrigger", "");
    Plotter::DatasetSummary dsDY(             "DY",           fileMap["DYJetsToLL"],      "",                "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsT1tttt1500( "T1tttt(1500,100)", fileMap["Signal_T1tttt_mGluino1500_mLSP100"], "", "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsT1tttt1200( "T1tttt(1200,800)", fileMap["Signal_T1tttt_mGluino1200_mLSP800"],      "",  "TriggerEffMC");
    Plotter::DatasetSummary dsT2tt500(    "T2tt(500,325)",    fileMap["Signal_T2tt_mStop500_mLSP325"],      "",        "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsT2tt800( "T2tt(850,100)",         fileMap["Signal_T2tt_mStop850_mLSP100"],      "",       "TriggerEffMC");
    //Plotter::DatasetSummary dsDYInc(          "DY HT<100",  fileMap["IncDY"],           "",       "");
    Plotter::DatasetSummary dstt2l(           "t#bar{t}",   fileMap["TTbarNoHad"],      "",                "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dstW(             "Single top", fileMap["tW"],              "",                "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsttZ(            "t#bar{t}Z",  fileMap["TTZ"],             "",                "genWeight;TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsVV(             "Diboson",    fileMap["Diboson"],         "",                "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsVVV(             "Other",    fileMap["Triboson"],         "",                "genWeight;TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsRare(           "Rare",       fileMap["Rare"],            "",                "genWeight;TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsQCD(             "QCD",    fileMap["QCD"],         "passQCD",                "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsZJets(             "Z(#nu#nu)+jets",    fileMap["ZJetsToNuNu"],         "",                "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsWJets(             "W(l#nu)+jets",    fileMap["WJetsToLNu"],         "",                "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary dsTTW(             "TTW",    fileMap["TTW"],         "",                "genWeight;TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary ds_TTBar(      "t#bar{t} (1l)",  fileMap["TTbarAll"], "",     "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary ds_ZINV(      "Z(#nu#nu)+jets",  fileMap["ZJetsToNuNu"],         "",                "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    Plotter::DatasetSummary ds_TTBarHT(   "TTbarHT", fileMap["TTbarHT"],    "",  "TriggerEffMC");//;bTagSF_EventWeightSimple_Central");
    std::vector<std::vector<Plotter::DatasetSummary>> stack_MC = {{dstt2l}, {dsWJets}, {dstW}, {dsZJets}, {dsttZ}, {dsQCD}, {dsRare,dsDY,dsVV}};
    std::vector<std::vector<Plotter::DatasetSummary>> stack_MC_hong = {{dstt2l}, {dsWJets}, {dstW}, {dsZJets}, {dsttZ}, {dsQCD}, {dsVVV,dsVV,dsTTW}};
    std::vector<std::vector<Plotter::DatasetSummary>> singal_points= {{dsT2tt800},{dsT2tt500},{dsT1tttt1200},{dsT1tttt1500}};
     // Collections for all variables, no cuts applied yet
    // met

    auto PDCMaker = [&](std::string var) 
    {
      //                                                   plotType  var  vector of datasetsummary
	return  std::vector<PDC>({ Plotter::DataCollection("data",   var, {dsData_SingleMuon}),
		                   Plotter::DataCollection("stack",  var, stack_MC),             
                                   Plotter::DataCollection("single",   var, singal_points)  });
    };

    auto PDCMaker1 = [&](std::string var)
    {
        return  std::vector<PDC>({ Plotter::DataCollection("data",   var, {dsData_HTMHT}),
                                   Plotter::DataCollection("stack",  var, stack_MC), 
                                   Plotter::DataCollection("single",   var, singal_points)  });
    };

    auto PDCMaker_nonesig = [&](std::string var)
    {
         return  std::vector<PDC>({ Plotter::DataCollection("data",   var, {dsData_SingleMuon}),
                                   Plotter::DataCollection("stack",  var, stack_MC)  });
    };

    auto PDCMaker1_nonesig = [&](std::string var)
    {
        return  std::vector<PDC>({ Plotter::DataCollection("data",   var, {dsData_HTMHT}),
                                   Plotter::DataCollection("stack",  var, stack_MC)  });
    };

    auto PDCMaker1_nonesig_notrig = [&](std::string var)
    {
         return  std::vector<PDC>({ Plotter::DataCollection("data",   var, {dsData_SingleMuonNotrig}),
                                   Plotter::DataCollection("stack",  var, stack_MC)  });
    };
    auto PDCMaker1_rare = [&](std::string var)
    {
        return  std::vector<PDC>({ Plotter::DataCollection("data",   var, {dsData_HTMHT}),
                                   Plotter::DataCollection("stack",  var, stack_MC_hong),
                                   Plotter::DataCollection("single",   var, singal_points)  });
    };

    auto PDCMaker12 = [&](std::string var)
    {
        return std::vector<PDC>({Plotter::DataCollection("single",   var, {ds_ZINV}),
                                 Plotter::DataCollection("single",   var, {ds_TTBar}),
                                 Plotter::DataCollection("single",   var, {dsT1tttt1200}),
                                 //Plotter::DataCollection("single",   var, {dsT2tt800}),
                                 Plotter::DataCollection("single",   var, {dsT1tttt1500}),
                                 Plotter::DataCollection("single",   var, {dsT2tt500}) });
                                 //Plotter::DataCollection("single",   var, {ds_TTBarHT})  });
    };

    std::vector<std::pair<std::string,std::string>> cutlevels_muon = {
        //cut name                    cut string
	//{"Koushik",                     "passnJets;passdPhis;passMET;passBJets;passTagger;passHT;passMT2;passNoiseEventFilter"},//;passKoushik"},
	{"Empty",                       ""},
	//{"muZinv",                    "passNoiseEventFilterZinv;passMuZinvSel"},
	//{"muZinv_ht200",              "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200"},
	//{"muZinv_ht200_dphi",         "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passdPhisZinv"},
	//{"muZinv_ht50_met50_dphi",    "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
	//{"muZinv_loose0",             "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
	//{"muZinv_loose0_mt2",         "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>0"},
	//{"muZinv_loose0_mt21b",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv1b>0"},
	//{"muZinv_loose0_mt22b",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv2b>0"},
	//{"muZinv_loose0_mt23b",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv3b>0"},
	//{"muZinv_loose0_ntop",        "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv>0"},
	//{"muZinv_loose0_nb2",         "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv>=2"},
	//{"muZinv_loose0_ht300",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>300;passnJetsZinv;passdPhisZinv"},
	//{"muZinv_loose0_ht400",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>400;passnJetsZinv;passdPhisZinv"},
	//{"muZinv_loose0_ht500",       "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>500;passnJetsZinv;passdPhisZinv"},
	//{"muZinv_blnotagmt2",         "passMuZinvSel;passBaselineNoTagMT2Zinv"},
	//{"dr_min_2",              "ak82dRMin<0.2"},
  //      {"dr_med_2",              "0.2<ak82dRMin<=0.4"},
   //     {"dr_large_2",            "0.4<ak82dRMin<0.8"},
    //    {"dr_min_1",              "ak81dRMin<=0.2"},
     //   {"dr_med_1",              "0.2<ak81dRMin<=0.4"},
      //  {"dr_large_1",            "0.4<ak81dRMin<0.8"},
        //{"baseline",            "passBaseline"},
	//{"baseline_Zinv",                 "passBaselineZinv;passNoiseEventFilter"},//;passKoushik"},
	//{"muZinv_0b",                 "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0"},
/*
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
	{"muZinv_g1b_blnotag",        "passMuZinvSel;passBJetsZinv;passBaselineNoTagZinv"},
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
*/ 
   };


    //push the histograms in a loop, save some copy-paste time
    /*
    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
	// unweighted
        //               histogram name                         vector of PDC           ratio  cutstring    bins ll ul   log   norm    x-axis      y-axis     
       vh.push_back(PHS("DataMC_SingleMuon_model_metZinv_" + cut.first,  PDCMaker("cleanMetPt"), {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_ht_Zinv"  + cut.first,  PDCMaker("HTZinv"), {1, 2}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_ht_"  + cut.first,  PDCMaker("HT"),     {1, 2}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_met_"  + cut.first,  PDCMaker("met"),     {1, 2}, cut.second, 50, 0, 1500, true, true,  label_met,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_NBJetsZInv_"  + cut.first,  PDCMaker("cntCSVSZinv"),   {1, 2}, cut.second, 6, 0, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_NBJEts_"  + cut.first,  PDCMaker("cntCSVS"),     {1, 2}, cut.second, 6, 0, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_NJEtsZInv_" + cut.first, PDCMaker("cntNJetsPt30Eta24Zinv"), {1, 2}, cut.second, 20, 0, 20, false, true,  label_nj, "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_NJEts_"  + cut.first,  PDCMaker("cntNJetsPt30Eta24"),{1, 2}, cut.second, 20, 0, 20, false, true,  label_nj,   "Events")); 
       vh.push_back(PHS("DataMC_SingleMuon_model_NTopsZinv_"  + cut.first,  PDCMaker("nTopCandSortedCntZinv"), {1, 2}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_NTops_"  + cut.first,  PDCMaker("nTopCandSortedCnt"),     {1, 2}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_MT2ZInv_"  + cut.first,  PDCMaker("best_had_brJet_MT2Zinv"), {1, 2}, cut.second, 100, 0, 1500, true, true,  label_mt2, "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_MT2_"  + cut.first,  PDCMaker("best_had_brJet_MT2"),     {1, 2}, cut.second, 100, 0, 1500, true, true,  label_mt2,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_metNoratio_" + cut.first,  PDCMaker("met"), {1, 1}, cut.second, 50, 0, 1500, true, true,  label_met,  "Events")); 
       vh.push_back(PHS("DataMC_SingleMuon_model_htNoratio_"  + cut.first,  PDCMaker("HT"),     {1, 1}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_NBJEtsNoratio_"  + cut.first,  PDCMaker("cntCSVS"),     {1, 1}, cut.second, 6, 0, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_NJEtsNoratio_"  + cut.first,  PDCMaker("cntNJetsPt30Eta24"),{1, 1}, cut.second, 20, 0, 20, false, true,  label_nj,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_NTopsNoratio_"  + cut.first, PDCMaker("nTopCandSortedCnt"), {1, 1}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_model_MT2Noratio_"  + cut.first,  PDCMaker("best_had_brJet_MT2"), {1, 1}, cut.second, 100, 0, 1500, true, true,  label_mt2,   "Events"));
    }
*/
    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
       vh.push_back(PHS("DataMC_HTMHT_model_metZinv_" + cut.first,  PDCMaker12("cleanMetPt"), {1, 2}, cut.second, 26, 200, 1500, true, true,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_HTMTH_model_ht_Zinv"  + cut.first,  PDCMaker12("HTZinv"), {1, 2}, cut.second, 40, 500, 2500, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("_ht_"      + cut.first,  PDCMaker12("ht"),     {1, 2}, cut.second, 40, 500, 2500, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_model_met_"     + cut.first,  PDCMaker12("met"),     {1, 2}, cut.second, 26, 200, 1500, true, true,  label_met,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_model_NBJetsZInv_" + cut.first,  PDCMaker12("cntCSVSZinv"),   {1, 2}, cut.second, 5, 1, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_model_NBJEts_"     + cut.first,  PDCMaker12("cntCSVS"),     {1, 2}, cut.second, 5, 1, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_model_NJEtsZInv_"  + cut.first, PDCMaker12("cntNJetsPt30Eta24Zinv"), {1, 2}, cut.second, 12, 4, 16, false, true,  label_nj, "Events"));
     //  vh.push_back(PHS("DataMC_HTMTH_model_NJEts_"      + cut.first,  PDCMaker12("cntNJetsPt30Eta24"),{1, 2}, cut.second, 12, 4, 16, false, true,  label_nj,   "Events"));
     //  vh.push_back(PHS("DataMC_HTMTH_model_NTops_"      + cut.first,  PDCMaker12("nTopCandSortedCntZinv"), {1, 2}, cut.second, 5, 1, 6, false, true,  label_nt,   "Events"));
     //  vh.push_back(PHS("DataMC_HTMTH_model_NTops_"      + cut.first,  PDCMaker12("nTopCandSortedCnt"),     {1, 2}, cut.second, 5, 1, 6, false, true,  label_nt,   "Events"));
     //vh.push_back(PHS("DataMC_HTMTH_model_MT2ZInv_"  + cut.first,  PDCMaker12("best_had_brJet_MT2Zinv"),  {1, 2}, cut.second, 20, 200, 1200, false, true,  label_mt2,   "Events"));
       //vh.push_back(PHS("DataMC_HTMTH_model_MT2_"  + cut.first,  PDCMaker12("best_had_brJet_MT2"), {1, 2}, cut.second, 20, 200, 1200, false, true,  label_mt2,   "Events"));
       //vh.push_back(PHS("DataMC_HTMTH_model_metNoratio_" + cut.first,  PDCMaker12("met"), {1, 1}, cut.second, 26, 200, 1500, true, true,  label_met,  "Events"));
      // vh.push_back(PHS("DataMC_HTMTH_model_htNoratio_"  + cut.first,  PDCMaker12("HT"),     {1, 1}, cut.second, 40, 500, 2500, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("slimmedJetsAK8"  + cut.first,  PDCMaker12("slimmedJetsAK8"),  {1, 1}, cut.second, 20, 200, 1200, false, false,  label_mt2,   "Events"));
       vh.push_back(PHS("Tau1"  + cut.first,  PDCMaker12("tau1"),  {1, 1}, cut.second, 100, 0, 1, false, true,  "tau1",   "Events"));
       vh.push_back(PHS("TauD"  + cut.first,  PDCMaker12("puppitau2Dtau1"),  {1, 1}, cut.second, 100, 0, 1, false, true,  "tauD",   "Events"));
       vh.push_back(PHS("puppiTauD32"  + cut.first,  PDCMaker12("puppitau3Dtau2"),  {1, 1}, cut.second, 100, 0, 1, false, true,  "tauD",   "Events"));
       vh.push_back(PHS("Tau2"  + cut.first,  PDCMaker12("tau2"),  {1, 1}, cut.second, 100, 0, 1, false, true,  "tau2",   "Events"));
       vh.push_back(PHS("Tau3"  + cut.first,  PDCMaker12("tau3"),  {1, 1}, cut.second, 100, 0, 1, false, true,  "tau3",   "Events"));
       vh.push_back(PHS("SubjetB"  + cut.first,  PDCMaker12("subjetBdisc"),  {1, 1}, cut.second, 9, -11, 2, true, false,  "SubjetBdisc",   "Events"));
       vh.push_back(PHS("ak8pt"  + cut.first,  PDCMaker12("ak8pt"),  {1, 1}, cut.second, 20, 0, 400, false, false,  "ak8PT",   "Events"));
       vh.push_back(PHS("ak8mass"  + cut.first,  PDCMaker12("ak8mass"),  {1, 1}, cut.second, 10, 0, 150, false, false,  "ak8mass",   "Events"));
       vh.push_back(PHS("ak8rapi"  + cut.first,  PDCMaker12("ak8rapidity"),  {1, 1}, cut.second, 9, -4, 4, false, false,  "ak8rapidity",   "Events"));
       vh.push_back(PHS("ak8Pmass"  + cut.first,  PDCMaker12("prunedMass"),  {1, 1}, cut.second, 10, 0, 150, false, false,  "ak8prunedMass",   "Events"));
       vh.push_back(PHS("ak8softDrop"  + cut.first,  PDCMaker12("softDropMass"),  {1, 1}, cut.second, 10, 0, 150, false, true,  "ak8Soft",   "Events"));
       vh.push_back(PHS("puppisoftDrop"  + cut.first,  PDCMaker12("puppisoftDropMass"),  {1, 1}, cut.second, 10, 0, 150, false, true,  "ak8Soft",   "Events"));
       vh.push_back(PHS("nJets"+ cut.first,  PDCMaker12("nJets"),  {1, 1}, cut.second, 30, 0, 30, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsAk8"+ cut.first,  PDCMaker12("nJetsAk8"),  {1, 1}, cut.second, 15, 0, 15, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("ak81dRMin"+ cut.first,  PDCMaker12("ak81dRMin"),  {1, 1}, cut.second, 70, 0, 7, false, true,  "dR",   "Events"));
       vh.push_back(PHS("ak82dRMin"+ cut.first,  PDCMaker12("ak82dRMin"),  {1, 1}, cut.second, 70, 0, 7, false, true,  "dR",   "Events"));
       vh.push_back(PHS("ak81dRMin_0"+ cut.first,  PDCMaker12("ak81dRMin[0]"),  {1, 1}, cut.second, 40, 0, 4, false, true,  "dR",   "Events"));
       vh.push_back(PHS("ak82dRMin_0"+ cut.first,  PDCMaker12("ak82dRMin[0]"),  {1, 1}, cut.second, 40, 0, 4, false, true,  "dR",   "Events"));
       vh.push_back(PHS("nJetsAK4_l0.2"+ cut.first,  PDCMaker12("nJetsAK41_min"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsAK4_b0.20.4"+ cut.first,  PDCMaker12("nJetsAK41_med"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsAK4_b0.40.8"+ cut.first,  PDCMaker12("nJetsAK41_lar"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsAK41_l0.2"+ cut.first,  PDCMaker12("nJetsAK42_min"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsAK41_b0.20.4"+ cut.first,  PDCMaker12("nJetsAK42_med"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsAK41_b0.40.8"+ cut.first,  PDCMaker12("nJetsAK42_lar"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("puppi_Top_T_PT"+ cut.first,  PDCMaker12("puppiLVectight_top(pt)"),  {1, 1}, cut.second, 20, 0, 1200, false, true,  "PT",   "Events"));
       vh.push_back(PHS("puppi_Top_T_Num"+ cut.first,  PDCMaker12("puppiLVectight_top(size)"),  {1, 1}, cut.second, 5, 0, 5, false, true,  label_nt,   "Events")); 
       vh.push_back(PHS("puppi_Top_L_PT"+ cut.first,  PDCMaker12("puppiLVecLoose_top(pt)"),  {1, 1}, cut.second, 20, 0, 1200, false, true,  "PT",   "Events"));
       vh.push_back(PHS("puppi_Top_L_Num"+ cut.first,  PDCMaker12("puppiLVecLoose_top(size)"),  {1, 1}, cut.second, 5, 0, 5, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("puppi_w_T_PT"+ cut.first,  PDCMaker12("puppiLVectight_w(pt)"),  {1, 1}, cut.second, 50, 0, 1200, false, true,  "PT",   "Events"));
       vh.push_back(PHS("puppi_w_T_Num"+ cut.first,  PDCMaker12("puppiLVectight_w(size)"),  {1, 1}, cut.second, 5, 0, 5, false, true,  label_nw,   "Events"));
       vh.push_back(PHS("puppi_w_L_Num"+ cut.first,  PDCMaker12("puppiLVecLoose_w(size)"),  {1, 1}, cut.second, 5, 0, 5, false, true,  label_nw,   "Events"));
       vh.push_back(PHS("puppi_w_L_PT"+ cut.first,  PDCMaker12("puppiLVecLoose_w(pt)"),  {1, 1}, cut.second, 50, 0, 1200, false, true,  "PT",   "Events"));
       vh.push_back(PHS("nJetsPuppi_L_l0.2"+ cut.first,  PDCMaker12("nJetsPuppi_L1_min"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsPuppi_L_.20.4"+ cut.first,  PDCMaker12("nJetsPuppi_L1_med"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsPuppi_L_b0.40.8"+ cut.first,  PDCMaker12("nJetsPuppi_L1_lar"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsPuppi_T_l0.2"+ cut.first,  PDCMaker12("nJetsPuppi_T1_min"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsPuppi_T_.20.4"+ cut.first,  PDCMaker12("nJetsPuppi_T1_med"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("nJetsPuppi_T_b0.40.8"+ cut.first,  PDCMaker12("nJetsPuppi_T1_lar"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("puppiJetsLVec"  + cut.first,  PDCMaker12("puppiJetsLVec(pt)"),  {1, 1}, cut.second, 20, 200, 1200, false, false,  "pt",   "Events"));
       vh.push_back(PHS("puppiJetsLVec_size"  + cut.first,  PDCMaker12("puppiJetsLVec(size)"),  {1, 1}, cut.second, 7, 0, 7, false, true,  "nJets",   "Events"));
       vh.push_back(PHS("puppitau3Dtau2_SDM"+ cut.first,  PDCMaker12("puppitau3Dtau2_SDM"),  {1, 1}, cut.second, 100, 0, 1, false, true,  "tauD",   "Events"));
       vh.push_back(PHS("puppitau2Dtau1_SDM"+ cut.first,  PDCMaker12("puppitau2Dtau1_SDM"),  {1, 1}, cut.second, 100, 0, 1, false, true,  "tauD",   "Events"));
       }
/*
    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
       vh.push_back(PHS("DataMC_SingleMuon_metZinv_" + cut.first,  PDCMaker_nonesig("cleanMetPt"), {1, 2}, cut.second, 50, 0, 1500, true, true,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_ht_Zinv"  + cut.first,  PDCMaker_nonesig("HTZinv"), {1, 2}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_ht_"  + cut.first,  PDCMaker_nonesig("HT"),     {1, 2}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_met_"  + cut.first,  PDCMaker_nonesig("met"),     {1, 2}, cut.second, 50, 0, 1500, true, true,  label_met,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NBJetsZInv_"  + cut.first,  PDCMaker_nonesig("cntCSVSZinv"),   {1, 2}, cut.second, 6, 0, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NBJEts_"  + cut.first,  PDCMaker_nonesig("cntCSVS"),     {1, 2}, cut.second, 6, 0, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NJEtsZInv_" + cut.first, PDCMaker_nonesig("cntNJetsPt30Eta24Zinv"), {1, 2}, cut.second, 20, 0, 20, false, true,  label_nj, "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NJEts_"  + cut.first,  PDCMaker_nonesig("cntNJetsPt30Eta24"),{1, 2}, cut.second, 20, 0, 20, false, true,  label_nj,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NTops_"  + cut.first,  PDCMaker_nonesig("nTopCandSortedCntZinv"), {1, 2}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NTops_"  + cut.first,  PDCMaker_nonesig("nTopCandSortedCnt"),     {1, 2}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_MT2ZInv_"  + cut.first,  PDCMaker_nonesig("best_had_brJet_MT2Zinv"),{1, 2}, cut.second,100, 0, 1500, true, true,  label_mt2,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_MT2_"  + cut.first,  PDCMaker_nonesig("best_had_brJet_MT2"),     {1, 2}, cut.second, 100, 0, 1500, true, true,  label_mt2,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_metNoratio_" + cut.first,  PDCMaker_nonesig("met"), {1, 1}, cut.second, 50, 0, 1500, true, true,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_htNoratio_"  + cut.first,  PDCMaker_nonesig("HT"),     {1, 1}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NBJEtsNoratio_"  + cut.first,  PDCMaker_nonesig("cntCSVS"),     {1, 1}, cut.second, 6, 0, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NJEtsNoratio_"  + cut.first,  PDCMaker_nonesig("cntNJetsPt30Eta24"),{1, 1}, cut.second, 20, 0, 20, false, true,  label_nj,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_NTopsNoratio_"  + cut.first,  PDCMaker_nonesig("nTopCandSortedCnt"), {1, 1}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuon_MT2Noratio_"  + cut.first,  PDCMaker_nonesig("best_had_brJet_MT2"), {1, 1}, cut.second,100, 0, 1500, true, true,  label_mt2,   "Events"));
    }
*/
/*
    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
       vh.push_back(PHS("DataMC_HTMHT_metZinv_" + cut.first,  PDCMaker1_nonesig("cleanMetPt"), {1, 2}, cut.second, 26, 200, 1500, true, true,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_HTMTH_ht_Zinv"  + cut.first,  PDCMaker1_nonesig("HTZinv"), {1, 2}, cut.second, 40, 500, 2500, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__ht_"  + cut.first,  PDCMaker1_nonesig("HT"),     {1, 2}, cut.second, 40, 500, 2500, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__met_"  + cut.first,  PDCMaker1_nonesig("met"),     {1, 2}, cut.second, 26, 200, 1500, true, true,  label_met,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NBJetsZInv_"  + cut.first,  PDCMaker1_nonesig("cntCSVSZinv"),   {1, 2}, cut.second, 5, 1, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NBJEts_"  + cut.first,  PDCMaker1_nonesig("cntCSVS"),     {1, 2}, cut.second, 5, 1, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NJEtsZInv_" + cut.first, PDCMaker1_nonesig("cntNJetsPt30Eta24Zinv"), {1, 2}, cut.second, 12, 4, 16, false, true,  label_nj, "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NJEts_"  + cut.first,  PDCMaker1_nonesig("cntNJetsPt30Eta24"),{1, 2}, cut.second, 12, 4, 16, false, true,  label_nj,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NTops_"  + cut.first,  PDCMaker1_nonesig("nTopCandSortedCntZinv"), {1, 2}, cut.second, 5, 1, 6, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NTops_"  + cut.first,  PDCMaker1_nonesig("nTopCandSortedCnt"),     {1, 2}, cut.second, 5, 1, 6, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__MT2ZInv_"  + cut.first,  PDCMaker1_nonesig("best_had_brJet_MT2Zinv"), {1, 2}, cut.second, 20, 200, 1200, false, true,  label_mt2,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__MT2_"  + cut.first,  PDCMaker1_nonesig("best_had_brJet_MT2"),{1, 2}, cut.second, 20, 200, 1200, false, true,  label_mt2,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__metNoratio_" + cut.first,  PDCMaker1_nonesig("met"), {1, 1}, cut.second, 26, 200, 1500, true, true,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_HTMTH__htNoratio_"  + cut.first,  PDCMaker1_nonesig("HT"),     {1, 1}, cut.second, 40, 500, 2500, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NBJEtsNoratio_"  + cut.first,  PDCMaker1_nonesig("cntCSVS"),     {1, 1}, cut.second, 5, 1, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NJEtsNoratio_"  + cut.first,  PDCMaker1_nonesig("cntNJetsPt30Eta24"),{1, 1}, cut.second, 12, 4, 16, false, true,  label_nj,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__NTopsNoratio_"  + cut.first,  PDCMaker1_nonesig("nTopCandSortedCnt"),  {1, 1}, cut.second, 5, 1, 6, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH__MT2Noratio_"  + cut.first,  PDCMaker1_nonesig("best_had_brJet_MT2"),  {1, 1}, cut.second, 20, 200, 1200, false, true,  label_mt2,   "Events"));
    }

    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
       vh.push_back(PHS("DataMC_SingleMuonT_metZinv_" + cut.first,  PDCMaker1_nonesig_notrig("cleanMetPt"), {1, 2}, cut.second, 26, 200, 1500, true, false,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_SingleMuonT_ht_Zinv"  + cut.first,  PDCMaker1_nonesig_notrig("HTZinv"), {1, 2}, cut.second, 40, 500, 2500, true, false,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuonT_ht_"  + cut.first,  PDCMaker1_nonesig_notrig("HT"),     {1, 2}, cut.second, 40, 500, 2500, true, false,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuonT_met_"  + cut.first,  PDCMaker1_nonesig_notrig("met"),     {1, 2}, cut.second, 26, 200, 1500, true, false,  label_met,   "Events"));
    vh.push_back(PHS("DataMC_SingleMuonT_NBJetsZInv_"  + cut.first,  PDCMaker1_nonesig_notrig("cntCSVSZinv"),   {1, 2}, cut.second, 5, 1, 6, false, false,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuonT_NBJEts_"  + cut.first,  PDCMaker1_nonesig_notrig("cntCSVS"),     {1, 2}, cut.second, 5, 1, 6, false, false,  label_nb,   "Events"));
  vh.push_back(PHS("DataMC_SingleMuonT_NJEtsZInv_" + cut.first, PDCMaker1_nonesig_notrig("cntNJetsPt30Eta24Zinv"), {1, 2},cut.second, 12, 4, 16, false, false,  label_nj, "Events"));
      vh.push_back(PHS("DataMC_SingleMuonT_NJEts_"  + cut.first,  PDCMaker1_nonesig_notrig("cntNJetsPt30Eta24"),{1, 2}, cut.second, 12, 4, 16, false, false,  label_nj,   "Events"));
    vh.push_back(PHS("DataMC_SingleMuonT_NTopsZinv_"  + cut.first,PDCMaker1_nonesig_notrig("nTopCandSortedCntZinv"),{1, 2}, cut.second, 5, 1, 6, false, false, label_nt,"Events"));
       vh.push_back(PHS("DataMC_SingleMuonT_NTops_"  + cut.first,  PDCMaker1_nonesig_notrig("nTopCandSortedCnt"),{1, 2}, cut.second, 5, 1, 6, false, false,  label_nt, "Events"));
vh.push_back(PHS("DataMC_SingleMuonT_MT2ZInv_"  + cut.first, PDCMaker1_nonesig_notrig("best_had_brJet_MT2Zinv"),{1, 2}, cut.second,20, 200, 1200,false, false,  label_mt2, "Events"));
 vh.push_back(PHS("DataMC_SingleMuonT_MT2_"  + cut.first,  PDCMaker1_nonesig_notrig("best_had_brJet_MT2"),{1, 2}, cut.second, 20, 200, 1200, false, false,  label_mt2,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuonT_metNoratio_" + cut.first,  PDCMaker1_nonesig_notrig("met"), {1, 1}, cut.second, 26, 200, 1500, true, false,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_SingleMuonT_htNoratio_"  + cut.first,  PDCMaker1_nonesig_notrig("HT"),  {1, 1}, cut.second, 40, 500, 2500, true, false,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_SingleMuonT_NBJEtsNoratio_"  + cut.first,  PDCMaker1_nonesig_notrig("cntCSVS"), {1, 1}, cut.second, 5, 1, 6, false, false,  label_nb,   "Events"));
 vh.push_back(PHS("DataMC_SingleMuonT_NJEtsNoratio_"  + cut.first,  PDCMaker1_nonesig_notrig("cntNJetsPt30Eta24"),{1, 1}, cut.second, 12, 4, 16, false, false,  label_nj,"Events"));
 vh.push_back(PHS("DataMC_SingleMuonT_NTopsNoratio_"  + cut.first,  PDCMaker1_nonesig_notrig("nTopCandSortedCnt"),{1, 1}, cut.second, 5, 1, 6, false, false,  label_nt,"Events"));
 vh.push_back(PHS("DataMC_SingleMuonT_MT2Noratio_"  + cut.first,  PDCMaker1_nonesig_notrig("best_had_brJet_MT2"),{1, 1},cut.second, 20,200, 1200, false, false,  label_mt2, "Events"));
    }
*/
/*
    for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
       vh.push_back(PHS("DataMC_HTMHT_rare_metZinv_" + cut.first,  PDCMaker1_rare("cleanMetPt"), {1, 2}, cut.second, 50, 0, 1500, true, true,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_ht_Zinv"  + cut.first,  PDCMaker1_rare("HTZinv"), {1, 2}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_ht_"  + cut.first,  PDCMaker1_rare("HT"),     {1, 2}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_met_"  + cut.first,  PDCMaker1_rare("met"),     {1, 2}, cut.second, 50, 0, 1500, true, true,  label_met,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NBJetsZInv_"  + cut.first,  PDCMaker1_rare("cntCSVSZinv"),   {1, 2}, cut.second, 6, 0, 10, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NBJEts_"  + cut.first,  PDCMaker1_rare("cntCSVS"),     {1, 2}, cut.second, 6, 0, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NJEtsZInv_" + cut.first, PDCMaker1_rare("cntNJetsPt30Eta24Zinv"), {1, 2}, cut.second, 20, 0, 20, false, true,  label_nj, "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NJEts_"  + cut.first,  PDCMaker1_rare("cntNJetsPt30Eta24"),{1, 2}, cut.second, 20, 0, 20, false, true,  label_nj,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NTops_"  + cut.first,  PDCMaker1_rare("nTopCandSortedCntZinv"), {1, 2}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NTops_"  + cut.first,  PDCMaker1_rare("nTopCandSortedCnt"),     {1, 2}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_MT2ZInv_"  + cut.first,  PDCMaker1_rare("best_had_brJet_MT2Zinv"),{1, 2}, cut.second, 100, 0, 1500, true, true,  label_mt2,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_MT2_"  + cut.first,  PDCMaker1_rare("best_had_brJet_MT2"),{1, 2}, cut.second, 100, 0, 1500, true, true,  label_mt2,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_metNoratio_" + cut.first,  PDCMaker1_rare("met"), {1, 1}, cut.second, 50, 0, 1500, true, true,  label_met,  "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_htNoratio_"  + cut.first,  PDCMaker1_rare("HT"),     {1, 1}, cut.second, 50, 0, 2000, true, true,  label_ht,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NBJEtsNoratio_"  + cut.first,  PDCMaker1_rare("cntCSVS"),     {1, 1}, cut.second, 6, 0, 6, false, true,  label_nb,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NJEtsNoratio_"  + cut.first,  PDCMaker1_rare("cntNJetsPt30Eta24"),{1, 1}, cut.second, 20, 0, 20, false, true,  label_nj,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_NTopsNoratio_"  + cut.first,  PDCMaker1_rare("nTopCandSortedCnt"),   {1, 1}, cut.second, 10, 0, 10, false, true,  label_nt,   "Events"));
       vh.push_back(PHS("DataMC_HTMTH_rare_MT2Noratio_"  + cut.first,  PDCMaker1_rare("best_had_brJet_MT2"),{1, 1}, cut.second, 100, 0, 1500, true, true,  label_mt2,   "Events"));
    }
*/

vector<string> cfsSync = {"",
                              "passNoiseEventFilter",
                              "passNoiseEventFilter;passMuonVeto",
                              "passNoiseEventFilter;passMuonVeto;passEleVeto",
                              "passNoiseEventFilter;passMuonVeto;passEleVeto;passIsoTrkVeto",
                              "passNoiseEventFilter;passMuonVeto;passEleVeto;passIsoTrkVeto;passnJets",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passBJets",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passBJets;passHT",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passBJets;passHT;passMET",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passBJets;passMET;passHT;passdPhis",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets;passMET;passHT;passTagger",
                              "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets;passMET;passHT;passTagger;passMT2",
                              "passBaseline"};

vector<Plotter::CutFlowSummary> cutFlowSummaries;

cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Data",          PDC("", "" , {dsData_HTMHT}),        cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("MC",            PDC(" ", " ", {stack_MC}),            cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("t#bar{t}",            PDC(" ", " ", {dstt2l}),            cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("W(l#nu)+jets",            PDC(" ", " ", {dsWJets}),            cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Single top",            PDC(" ", " ", {dstW}),            cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Z(#nu#nu)+jets",            PDC(" ", " ", {dsZJets}),            cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("t#bar{t}Z",            PDC(" ", " ", {dsttZ}),            cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("QCD",            PDC(" ", " ", {dsQCD}),            cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Rare",            PDC(" ", " ", {dsRare,dsDY,dsVV}),            cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Signal Points", PDC(" ", " ", {singal_points}),       cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("T1tttt(1500,100)", PDC(" ", " ", {dsT1tttt1500}),       cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("T1tttt(1200,800)", PDC(" ", " ", {dsT1tttt1200}),       cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("T2tt(500,325)", PDC(" ", " ", {dsT2tt500}),       cfsSync));
cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("T2tt(850,100)", PDC(" ", " ", {dsT2tt800}),       cfsSync));

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
