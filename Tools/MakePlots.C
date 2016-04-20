#include "Plotter.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/samples.h"

#include <getopt.h>
#include <iostream>

int main(int argc, char* argv[])
{
  using namespace std;

  int opt;
  int option_index = 0;
  static struct option long_options[] = {
    {"plot",  no_argument, 0, 'p'},
    {"savehist",  no_argument, 0, 's'},
    {"savetuple", no_argument, 0, 't'},
    {"fromFile",  no_argument, 0, 'f'},
    {"condor", no_argument, 0, 'c'},
    {"histFile", required_argument, 0, 'H'},
    {"dataSets", required_argument, 0, 'D'},
    {"numFiles", required_argument, 0, 'N'},
    {"startFile",  required_argument, 0, 'M'},
    {"numEvts", required_argument, 0, 'E'},
    {"plotDir", required_argument, 0, 'P'},
    {"luminosity", required_argument, 0, 'L'}
  };

  bool doPlots = true, doSave = true, doTuple = true, fromTuple = true, runOnCondor = false;
  string histFile = "", dataSets = "", sampleloc = AnaSamples::fileDir, plotDir = "plots";
  int nFiles = -1, startFile = 0, nEvts = -1;
  double lumi = AnaSamples::luminosity;

  while((opt = getopt_long(argc, argv, "pstfcH:D:N:M:E:P:L:", long_options, &option_index)) != -1)
    {
      switch(opt)
	{
	case 'p':
	  if(doPlots) doSave  = doTuple = false;
	  else doPlots = true;
	  break;

	case 's':
	  if(doSave) doPlots = doTuple = false;
	  else doSave  = true;
	  break;

	case 't':
	  if(doTuple) doPlots = doSave = false;
	  else doTuple  = true;
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

  AnaSamples::SampleSet ss(sampleloc, lumi);
  AnaSamples::SampleCollection sc(ss);

  const double zAcc = 1.0;
  // const double zAcc = 0.5954;
  // const double zAcc = 0.855;
  const double znunu_mumu_ratio = 5.942;
  const double znunu_ee_ratio = 5.942;

  map<string, vector<AnaSamples::FileSummary>> fileMap;

  //Select approperiate datasets here
  if(dataSets.compare("TEST") == 0)
    {
      fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_600toInf"]};
      fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_600toInf"]};
      fileMap["DYJetsToLL_HT_600toInf"] = {ss["DYJetsToLL_HT_600toInf"]};
      fileMap["ZJetsToNuNu_HT_600toInf"] = {ss["ZJetsToNuNu_HT_600toInf"]};
      fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
      fileMap["TTbarNoHad"] = {ss["TTbarDiLep"]};
    }
  else if(dataSets.compare("TEST2") == 0)
    {
      fileMap["DYJetsToLL"]  = {ss["DYJetsToLL_HT_400to600"]};
      fileMap["ZJetsToNuNu"] = {ss["ZJetsToNuNu_HT_400to600"]};
      fileMap["DYJetsToLL_HT_400to600"] = {ss["DYJetsToLL_HT_400to600"]};
      fileMap["ZJetsToNuNu_HT_400to600"] = {ss["ZJetsToNuNu_HT_400to600"]};
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
  int NSB = 37; // 45

  // Shortcuts for axis labels
  std::string label_met = "p_{T}^{miss} [GeV]";
  std::string label_ht  = "H_{T} [GeV]";
  std::string label_mht = "MH_{T} [GeV]";
  std::string label_nj  = "N_{jets}";
  std::string label_nb  = "N_{b-jets}";
  std::string label_nt  = "N_{tops}";
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


  Plotter::DatasetSummary dsDY_ll_inc( "DY#rightarrow#mu#mu Inc", fileMap["IncDY"], "pdgIdZDec=13;passMuZinvSel", "");

  Plotter::DatasetSummary dsDY_ll_noSel(  "DY#rightarrow#mu#mu no sel", fileMap["DYJetsToLL"], "pdgIdZDec=13", "");
  Plotter::DatasetSummary dsDY_ll_gen_noSel(  "DY#rightarrow#mu#mu Gen no Sel", fileMap["DYJetsToLL"], "pdgIdZDec=13", "");

  Plotter::DatasetSummary dsTT_mm("t#bar{t}#rightarrow#mu#mu", fileMap["TTbarDiLep"], "", "");
  Plotter::DatasetSummary dsTT_ee("t#bar{t}#rightarrowee", fileMap["TTbarDiLep"], "", "");
  Plotter::DatasetSummary dsTT_em("t#bar{t}#rightarrowe#mu", fileMap["TTbarDiLep"], "", "");

  Plotter::DatasetSummary dsTT_mm_lepSel("t#bar{t}#rightarrow#mu#mu", fileMap["TTbarDiLep"], "passMuZinvSel", "");
  Plotter::DatasetSummary dsTT_ee_lepSel("t#bar{t}#rightarrowee", fileMap["TTbarDiLep"], "passElecZinvSel", "");
  Plotter::DatasetSummary dsTT_em_lepSel("t#bar{t}#rightarrowe#mu", fileMap["TTbarDiLep"], "passElMuZinvSel", "");

  Plotter::DatasetSummary dsDY_ee( "DY#rightarrowee",  fileMap["DYJetsToLL"], "pdgIdZDec=11;passElecZinvSel", "");
  Plotter::DatasetSummary dsDY_ee_scaled( "DY#rightarrowee",  fileMap["DYJetsToLL"], "pdgIdZDec=11;passElecZinvSel", "",  znunu_ee_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ee_zAcc_scaled( "DY#rightarrowee no e, Z eff+acc", fileMap["DYJetsToLL"], "pdgIdZDec=11;passElecZinvSel", "zEffWgtElec;zAccWgtElec", znunu_ee_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ee_zWeight_scaled("DY#rightarrowee no e, Z eff", fileMap["DYJetsToLL"], "pdgIdZDec=11;passElecZinvSel", "zEffWgtElec",  znunu_ee_ratio / zAcc);

  Plotter::DatasetSummary dsDY_ll( "DY#rightarrow#mu#mu",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "");
  Plotter::DatasetSummary dsDY_ll_clean(  "DY#rightarrow#mu#mu #mu cleaned", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "");
  Plotter::DatasetSummary dsDY_ll_base( "DY#rightarrow#mu#mu baseline",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;passBaselineZinv", "");
  Plotter::DatasetSummary dsDY_ll_baseNT( "DY#rightarrow#mu#mu baselineNoTag", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;passBaselineNoTagZinv", "");
  Plotter::DatasetSummary dsDY_ll_zWeight("DY#rightarrow#mu#mu no #mu, Z eff", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "zEffWgt");
  Plotter::DatasetSummary dsDY_ll_zAcc( "DY#rightarrow#mu#mu no #mu, Z eff+acc", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_zAccAcc("DY#rightarrow#mu#mu no #mu, Z acc", fileMap["DYJetsToLL"], "pdgIdZDec=13;ngenMuInAcc>1.5", "zAccWgt");
  Plotter::DatasetSummary dsDY_ll_zAccAcc_nw("DY#rightarrow#mu#mu no #mu, Z acc NW", fileMap["DYJetsToLL"], "pdgIdZDec=13;ngenMuInAcc>1.5", "");
  Plotter::DatasetSummary dsDY_ll_gen( "DY#rightarrow#mu#mu Gen", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "");
  Plotter::DatasetSummary dsDY_llclean( "DY#rightarrow#mu#mu no #mu", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "");
  Plotter::DatasetSummary dsDY_ll_mu30_mu20( "DY#rightarrow#mu#mu #mu p_{T} > 30, 20 GeV", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;cutMuPt1>30;cutMuPt2>20", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_mu30_mu30( "DY#rightarrow#mu#mu #mu p_{T} > 30, 30 GeV", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;cutMuPt1>30;cutMuPt2>30", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_mu40_mu30( "DY#rightarrow#mu#mu #mu p_{T} > 40, 30 GeV", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;cutMuPt1>40;cutMuPt2>30", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_mu35_mu20( "DY#rightarrow#mu#mu #mu p_{T} > 35, 20 GeV", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;cutMuPt1>35;cutMuPt2>20", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_mu40_mu20( "DY#rightarrow#mu#mu #mu p_{T} > 40, 20 GeV", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;cutMuPt1>40;cutMuPt2>20", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_mu45_mu20( "DY#rightarrow#mu#mu #mu p_{T} > 45, 20 GeV", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;cutMuPt1>45;cutMuPt2>20", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_singleMu45("DY#rightarrow#mu#mu SingleMu45eta2p1", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;passSingleMu45",  "zEffWgt;zAccWgt");

  Plotter::DatasetSummary dsDY_ll_scaled(  "DY#rightarrow#mu#mu", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel",  "", znunu_mumu_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ll_scaled_baseNT( "DY#rightarrow#mu#mu baselineNoTag", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;passBaselineNoTagZinv", "", znunu_mumu_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ll_zWeight_scaled("DY#rightarrow#mu#mu no #mu, Z eff", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "zEffWgt",  znunu_mumu_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ll_zAcc_scaled( "DY#rightarrow#mu#mu no #mu, Z eff+acc", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ll_trig_scaled( "DY#rightarrow#mu#mu no #mu, Z eff+acc+trig", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel;passDiMuIsoTrig", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
  Plotter::DatasetSummary dsDY_llclean_scaled( "DY#rightarrow#mu#mu no #mu", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuZinvSel", "", znunu_mumu_ratio / zAcc);

  Plotter::DatasetSummary dsDY_ll_forAcc( "DY#rightarrow#mu#mu no Sel forAcc", fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMu");
  Plotter::DatasetSummary dsDY_ll_forEff( "DY#rightarrow#mu#mu no Sel forEff", fileMap["DYJetsToLL"],  "pdgIdZDec=13", "");
  Plotter::DatasetSummary dsDY_ll_forEff_num( "DY#rightarrow#mu#mu no Sel N(#mu) match wgt", fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMatchMuInAcc");
  Plotter::DatasetSummary dsDY_ll_forEff_den( "DY#rightarrow#mu#mu no Sel N(#mu) wgt", fileMap["DYJetsToLL"],  "pdgIdZDec=13", "ngenMuInAcc");

  Plotter::DatasetSummary dsDY_nunu( "Z#rightarrow#nu#nu", fileMap["ZJetsToNuNu"], "passLeptVeto", "");
  Plotter::DatasetSummary dsDY_nunu_SB0b( "Z#rightarrow#nu#nu, N(b) = 0",  fileMap["ZJetsToNuNu"], "passLeptVeto", "");
  Plotter::DatasetSummary dsDY_nunu_SBMb( "Z#rightarrow#nu#nu, Direct MC", fileMap["ZJetsToNuNu"], "passLeptVeto", "");
  Plotter::DatasetSummary dsDY_test( "test", fileMap["ZJetsToNuNu"], "passLeptVeto", "");

  Plotter::DatasetSummary dsDY_ll_0t( "DY#rightarrow#mu#mu N(t) = 0", fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=0", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_1t( "DY#rightarrow#mu#mu N(t) = 1", fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=1", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_2t( "DY#rightarrow#mu#mu N(t) = 2", fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv=2", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_3t( "DY#rightarrow#mu#mu N(t) > 2", fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;nTopCandSortedCntZinv>2", "zEffWgt;zAccWgt");

  Plotter::DatasetSummary dsDY_ll_0b( "DY#rightarrow#mu#mu N(b) = 0", fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=0", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_1b( "DY#rightarrow#mu#mu N(b) = 1", fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=1", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_2b( "DY#rightarrow#mu#mu N(b) = 2", fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv=2", "zEffWgt;zAccWgt");
  Plotter::DatasetSummary dsDY_ll_3b( "DY#rightarrow#mu#mu N(b) > 2", fileMap["DYJetsToLL"],  "pdgIdZDec=13;passMuZinvSel;cntCSVSZinv>2", "zEffWgt;zAccWgt");

  Plotter::DatasetSummary dsDY_nunu_0t(  "Z#rightarrow#nu#nu N(t) = 0", fileMap["ZJetsToNuNu"], "passLeptVeto;nTopCandSortedCntZinv=0", "");
  Plotter::DatasetSummary dsDY_nunu_1t(  "Z#rightarrow#nu#nu N(t) = 1", fileMap["ZJetsToNuNu"], "passLeptVeto;nTopCandSortedCntZinv=1", "");
  Plotter::DatasetSummary dsDY_nunu_2t(  "Z#rightarrow#nu#nu N(t) = 2", fileMap["ZJetsToNuNu"], "passLeptVeto;nTopCandSortedCntZinv=2", "");
  Plotter::DatasetSummary dsDY_nunu_3t(  "Z#rightarrow#nu#nu N(t) > 2", fileMap["ZJetsToNuNu"], "passLeptVeto;nTopCandSortedCntZinv>2", "");

  Plotter::DatasetSummary dsDY_nunu_0b(  "Z#rightarrow#nu#nu N(b) = 0", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0", "");
  Plotter::DatasetSummary dsDY_nunu_1b(  "Z#rightarrow#nu#nu N(b) = 1", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=1", "");
  Plotter::DatasetSummary dsDY_nunu_2b(  "Z#rightarrow#nu#nu N(b) = 2", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=2", "");
  Plotter::DatasetSummary dsDY_nunu_3b(  "Z#rightarrow#nu#nu N(b) > 2", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv>2", "");

  Plotter::DatasetSummary dsDY_nunu_0b_blNoTag(  "Z#rightarrow#nu#nu N(b) = 0", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0;passBaselineNoTagZinv", "");
  Plotter::DatasetSummary dsDY_nunu_1b_blNoTag(  "Z#rightarrow#nu#nu N(b) = 1", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=1;passBaselineNoTagZinv", "");
  Plotter::DatasetSummary dsDY_nunu_2b_blNoTag(  "Z#rightarrow#nu#nu N(b) = 2", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=2;passBaselineNoTagZinv", "");
  Plotter::DatasetSummary dsDY_nunu_3b_blNoTag(  "Z#rightarrow#nu#nu N(b) > 2", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv>2;passBaselineNoTagZinv", "");

  Plotter::DatasetSummary dsDY_nunu_0b_1fake(  "Z#rightarrow#nu#nu N(b) = 0, 1 fake b", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0", "weight1fakeb");
  Plotter::DatasetSummary dsDY_nunu_0b_2fake(  "Z#rightarrow#nu#nu N(b) = 0, 2 fake b", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0", "weight2fakeb");
  Plotter::DatasetSummary dsDY_nunu_0b_3fake(  "Z#rightarrow#nu#nu N(b) = 0, 3 fake b", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0", "weight3fakeb");

  Plotter::DatasetSummary dsDY_nunu_0b_1fake_blNoTag(  "Z#rightarrow#nu#nu N(b) = 0, 1 fake b", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0;passBaselineNoTagZinv1b", "weight1fakeb");
  Plotter::DatasetSummary dsDY_nunu_0b_2fake_blNoTag(  "Z#rightarrow#nu#nu N(b) = 0, 2 fake b", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0;passBaselineNoTagZinv2b", "weight2fakeb");
  Plotter::DatasetSummary dsDY_nunu_0b_3fake_blNoTag(  "Z#rightarrow#nu#nu N(b) = 0, 3 fake b", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=0;passBaselineNoTagZinv3b", "weight3fakeb");

  Plotter::DatasetSummary dsDY_nunu_highCSV(  "Z#rightarrow#nu#nu csvMax > 0.5", fileMap["ZJetsToNuNu"], "passLeptVeto;maxCSV>0.5", "");
  Plotter::DatasetSummary dsDY_nunu_veryhighCSV1b(  "Z#rightarrow#nu#nu csvMax > 0.8 1b", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=1;maxCSV>0.8", "");
  Plotter::DatasetSummary dsDY_nunu_highCSV1b(  "Z#rightarrow#nu#nu csvMax > 0.5 1b", fileMap["ZJetsToNuNu"], "passLeptVeto;cntCSVSZinv=1;maxCSV>0.5", "");
  Plotter::DatasetSummary dsDY_nunu_lowCSV(  "Z#rightarrow#nu#nu csvMax < 0.5", fileMap["ZJetsToNuNu"], "passLeptVeto;maxCSV<0.5", "");
  Plotter::DatasetSummary dsDY_nunu_verylowCSV(  "Z#rightarrow#nu#nu csvMax < 0.4", fileMap["ZJetsToNuNu"], "passLeptVeto;maxCSV<0.4", "");

  Plotter::DatasetSummary dsDY_ll_HT_100to200(  "DY#rightarrow#mu#mu 100 < H_{T} < 200 GeV", fileMap["DYJetsToLL_HT_100to200"], "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ll_HT_200to400(  "DY#rightarrow#mu#mu 200 < H_{T} < 400 GeV", fileMap["DYJetsToLL_HT_200to400"], "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ll_HT_400to600(  "DY#rightarrow#mu#mu 400 < H_{T} < 600 GeV", fileMap["DYJetsToLL_HT_400to600"], "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);
  Plotter::DatasetSummary dsDY_ll_HT_600toInf(  "DY#rightarrow#mu#mu H_{T} > 600 GeV",  fileMap["DYJetsToLL_HT_600toInf"], "pdgIdZDec=13;passMuZinvSel", "zEffWgt;zAccWgt", znunu_mumu_ratio / zAcc);

  Plotter::DatasetSummary dsDY_nunu_HT_100to200("DY#rightarrow#nu#nu 100 < H_{T} < 200 GeV", fileMap["ZJetsToNuNu_HT_100to200"],  "passLeptVeto", "");
  Plotter::DatasetSummary dsDY_nunu_HT_200to400("DY#rightarrow#nu#nu 200 < H_{T} < 400 GeV", fileMap["ZJetsToNuNu_HT_200to400"],  "passLeptVeto", "");
  Plotter::DatasetSummary dsDY_nunu_HT_400to600("DY#rightarrow#nu#nu 400 < H_{T} < 600 GeV", fileMap["ZJetsToNuNu_HT_400to600"],  "passLeptVeto", "");
  Plotter::DatasetSummary dsDY_nunu_HT_600toInf("DY#rightarrow#nu#nu H_{T} > 600 GeV",  fileMap["ZJetsToNuNu_HT_600toInf"],  "passLeptVeto", "");

  Plotter::DatasetSummary dsDY_ll_passMuonVeto( "passMuonVeto",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuonVeto", "");
  Plotter::DatasetSummary dsDY_ll_passEleVeto(  "passEleVeto",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passEleVeto", "");
  Plotter::DatasetSummary dsDY_ll_passIsoTrkVeto( "passIsoTrkVeto", fileMap["DYJetsToLL"], "pdgIdZDec=13;passIsoTrkVeto", "");
  Plotter::DatasetSummary dsDY_ll_passnJets( "passnJets", fileMap["DYJetsToLL"], "pdgIdZDec=13;passnJets", "");
  Plotter::DatasetSummary dsDY_ll_passdPhis( "passdPhis", fileMap["DYJetsToLL"], "pdgIdZDec=13;passdPhis", "");
  Plotter::DatasetSummary dsDY_ll_passBJets( "passBJets", fileMap["DYJetsToLL"], "pdgIdZDec=13;passBJets", "");
  Plotter::DatasetSummary dsDY_ll_passMET(  "passMET", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMET", "");
  Plotter::DatasetSummary dsDY_ll_passTagger( "passTagger", fileMap["DYJetsToLL"], "pdgIdZDec=13;passTagger", "");
  //Plotter::DatasetSummary dsDY_ll_passNewCuts(  "passNewCuts",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passNewCuts", "");
  Plotter::DatasetSummary dsDY_ll_passBaseline( "passBaseline",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passBaseline", "");
  Plotter::DatasetSummary dsDY_ll_passBaselineNoTag("passBaselineNoTag", fileMap["DYJetsToLL"], "pdgIdZDec=13;passBaselineNoTag", "");

  Plotter::DatasetSummary dsDY_ll_passMuonVetoZinv( "passMuonVetoZinv",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passMuonVetoZinv", "");
  Plotter::DatasetSummary dsDY_ll_passEleVetoZinv(  "passEleVetoZinv",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passEleVetoZinv", "");
  Plotter::DatasetSummary dsDY_ll_passIsoTrkVetoZinv( "passIsoTrkVetoZinv", fileMap["DYJetsToLL"], "pdgIdZDec=13;passIsoTrkVetoZinv", "");
  Plotter::DatasetSummary dsDY_ll_passnJetsZinv( "passnJetsZinv", fileMap["DYJetsToLL"], "pdgIdZDec=13;passnJetsZinv", "");
  Plotter::DatasetSummary dsDY_ll_passdPhisZinv( "passdPhisZinv", fileMap["DYJetsToLL"], "pdgIdZDec=13;passdPhisZinv", "");
  Plotter::DatasetSummary dsDY_ll_passBJetsZinv( "passBJetsZinv", fileMap["DYJetsToLL"], "pdgIdZDec=13;passBJetsZinv", "");
  Plotter::DatasetSummary dsDY_ll_passMETZinv(  "passMETZinv", fileMap["DYJetsToLL"], "pdgIdZDec=13;passMETZinv", "");
  Plotter::DatasetSummary dsDY_ll_passTaggerZinv( "passTaggerZinv", fileMap["DYJetsToLL"], "pdgIdZDec=13;passTaggerZinv", "");
  //Plotter::DatasetSummary dsDY_ll_passNewCutsZinv(  "passNewCutsZinv",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passNewCutsZinv", "");
  Plotter::DatasetSummary dsDY_ll_passBaselineZinv( "passBaselineZinv",  fileMap["DYJetsToLL"], "pdgIdZDec=13;passBaselineZinv", "");
  Plotter::DatasetSummary dsDY_ll_passBaselineNoTagZinv("passBaselineNoTagZinv", fileMap["DYJetsToLL"], "pdgIdZDec=13;passBaselineNoTagZinv", "");

  Plotter::DataCollection scaled_chargedEMFrac( "single", "cleanChargedHadEFrac", {dsDY_nunu, dsDY_ll_zAcc_scaled});
  Plotter::DataCollection scaled_neutralEMFrac( "single", "cleanNeutralEMEFrac",  {dsDY_nunu, dsDY_ll_zAcc_scaled});
  Plotter::DataCollection scaled_chargedHadFrac("single", "cleanChargedEMEFrac",  {dsDY_nunu, dsDY_ll_zAcc_scaled});

  Plotter::DataCollection scaled_chargedEMFrac_j1( "single", "cleanChargedHadEFrac[0]", {dsDY_nunu, dsDY_ll_zAcc_scaled});
  Plotter::DataCollection scaled_neutralEMFrac_j1( "single", "cleanNeutralEMEFrac[0]",  {dsDY_nunu, dsDY_ll_zAcc_scaled});
  Plotter::DataCollection scaled_chargedHadFrac_j1("single", "cleanChargedEMEFrac[0]",  {dsDY_nunu, dsDY_ll_zAcc_scaled});

  Plotter::DataCollection scaled_chargedEMFrac_j2( "single", "cleanChargedHadEFrac[1]", {dsDY_nunu, dsDY_ll_zAcc_scaled});
  Plotter::DataCollection scaled_neutralEMFrac_J2( "single", "cleanNeutralEMEFrac[1]",  {dsDY_nunu, dsDY_ll_zAcc_scaled});
  Plotter::DataCollection scaled_chargedHadFrac_j2("single", "cleanChargedEMEFrac[1]",  {dsDY_nunu, dsDY_ll_zAcc_scaled});

  //Basic Raw compairson plots
  vh.push_back(PHS("inc_mht", {PDC("single", "mht",  {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  label_mht,  "Events"));
  vh.push_back(PHS("inc_met", {PDC("single", "met",  {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  label_met, "Events"));
  vh.push_back(PHS("inc_ht",  {PDC("single", "ht", {dsDY_ll_inc, dsDY_ll})}, {1, 2}, "", 100, 0, 2000,  true,  false,  label_ht,  "Events"));
  vh.push_back(PHS("inc_genht", {PDC("single", "genHt",  {dsDY_ll_inc, dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 2000,  true,  true,  "gen H_{T} [GeV]",  "Norm Events"));

  vh.push_back(PHS("minljdR", {PDC("single", "minljdR",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "", 100, 0, 2.0,  true,  true,  "min #DeltaR(l,j)", "Norm Events"));
  vh.push_back(PHS("minljdR_baseline", {PDC("single", "minljdR",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",  100, 0, 2.0,  true,  true,  "min #DeltaR(l,j)", "Norm Events"));
  vh.push_back(PHS("minljdR_baselineNoTag",  {PDC("single", "minljdR",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 2.0,  true,  true,  "min #DeltaR(l,j)", "Norm Events"));

  vh.push_back(PHS("METdeltaPtSyst", {PDC("single", "dPtMet",  {dsDY_nunu})}, {1, 1}, "", 200, -2.0, 2.0,  false,  false,  "#DeltaP_{T} p_{T}^{miss} [GeV]", "Events"));
  vh.push_back(PHS("METdeltaPhiSyst", {PDC("single", "dPhiMet", {dsDY_nunu})}, {1, 1}, "", 100, -3.14, 3.14,  false,  false,  "#Delta#phi p_{T}^{miss}", "Events"));

  vh.push_back(PHS("dPhi1",  {PDC("single", "dPhiVecZinv[0]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {3, 2}, "", 100, 0, 3.1415,  false, true,  "#Delta#phi(p_{T}^{miss},j1)", "Norm Events"));
  vh.push_back(PHS("dPhi2",  {PDC("single", "dPhiVecZinv[1]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {3, 2}, "", 100, 0, 3.1415,  false, true,  "#Delta#phi(p_{T}^{miss},j2)", "Norm Events"));
  vh.push_back(PHS("dPhi3",  {PDC("single", "dPhiVecZinv[2]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {3, 2}, "", 100, 0, 3.1415,  false, true,  "#Delta#phi(p_{T}^{miss},j3)", "Norm Events"));
  vh.push_back(PHS("dPhi1_comp", {PDC("single", {{"dPhiVecZinv[0]", dsDY_ll_clean}, {"dPhiVec[0]", dsDY_ll}})}, {1, 1}, "", 100, 0, 3.1415,  false, true,  "#Delta#phi(p_{T}^{miss},j1)", "Norm Events"));
  vh.push_back(PHS("dPhi2_comp", {PDC("single", {{"dPhiVecZinv[1]", dsDY_ll_clean}, {"dPhiVec[1]", dsDY_ll}})}, {1, 1}, "", 100, 0, 3.1415,  false, true,  "#Delta#phi(p_{T}^{miss},j2)", "Norm Events"));
  vh.push_back(PHS("dPhi3_comp", {PDC("single", {{"dPhiVecZinv[2]", dsDY_ll_clean}, {"dPhiVec[2]", dsDY_ll}})}, {1, 1}, "", 100, 0, 3.1415,  false, true,  "#Delta#phi(p_{T}^{miss},j3)", "Norm Events"));
  vh.push_back(PHS("dPhi1_nJet", {PDC("single", "dPhiVecZinv[0]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j1)", "Norm Events"));
  vh.push_back(PHS("dPhi2_nJet", {PDC("single", "dPhiVecZinv[1]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j2)", "Norm Events"));
  vh.push_back(PHS("dPhi3_nJet", {PDC("single", "dPhiVecZinv[2]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j3)", "Norm Events"));
  vh.push_back(PHS("dPhi1_nJetMET", {PDC("single", "dPhiVecZinv[0]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {3, 2}, "passnJetsZinv;passMETZinv", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j1)", "Norm Events"));
  vh.push_back(PHS("dPhi2_nJetMET", {PDC("single", "dPhiVecZinv[1]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {3, 2}, "passnJetsZinv;passMETZinv", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j2)", "Norm Events"));
  vh.push_back(PHS("dPhi3_nJetMET", {PDC("single", "dPhiVecZinv[2]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {3, 2}, "passnJetsZinv;passMETZinv", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j3)", "Norm Events"));
  vh.push_back(PHS("dPhi1_nJet_met_gt_600", {PDC("single", "dPhiVecZinv[0]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;cleanMetPt>600", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j1)", "Norm Events"));
  vh.push_back(PHS("dPhi2_nJet_met_gt_600", {PDC("single", "dPhiVecZinv[1]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;cleanMetPt>600", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j2)", "Norm Events"));
  vh.push_back(PHS("dPhi3_nJet_met_gt_600", {PDC("single", "dPhiVecZinv[2]",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passnJetsZinv;cleanMetPt>600", 100, 0, 3.1415,  true,  true,  "#Delta#phi(p_{T}^{miss},j3)", "Norm Events"));
  vh.push_back(PHS("cleanJetPt",  {PDC("single", "jetsLVecLepCleaned(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJetPt_baseline", {PDC("single", "jetsLVecLepCleaned(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",  100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJetPt_baselineNoTag", {PDC("single", "jetsLVecLepCleaned(pt)",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, 0, 1500,  true,  true,  "Cleaned Jet p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("mindPhigenZJ", {PDC("single", "mindPhiMetJ",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "", 100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
  vh.push_back(PHS("mindPhigenZJ_baseline",  {PDC("single", "mindPhiMetJ",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv",  100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
  vh.push_back(PHS("mindPhigenZJ_baselineNoTag", {PDC("single", "mindPhiMetJ",  {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv", 100, -3.14, 3.14,  true,  true,  "min #DeltaPhi(genZ, jet)", "Norm Events"));
  vh.push_back(PHS("cleanJet1pt", {PDC("single", "jetsLVecLepCleaned[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",  100, 0, 1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJet1pt_baseline",  {PDC("single", "jetsLVecLepCleaned[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJet1pt_baselineNoTag", {PDC("single", "jetsLVecLepCleaned[0](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",  100, 0, 1500,  true,  true, "p_{T}(j1) [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJet2pt", {PDC("single", "jetsLVecLepCleaned[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",  100, 0, 1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJet2pt_baseline",  {PDC("single", "jetsLVecLepCleaned[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJet2pt_baselineNoTag", {PDC("single", "jetsLVecLepCleaned[1](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",  100, 0, 1500,  true,  true, "p_{T}(j2) [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJet3pt", {PDC("single", "jetsLVecLepCleaned[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "",  100, 0, 1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJet3pt_baseline",  {PDC("single", "jetsLVecLepCleaned[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));
  vh.push_back(PHS("cleanJet3pt_baselineNoTag", {PDC("single", "jetsLVecLepCleaned[2](pt)", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {1, 2}, "passBaselineNoTagZinv",  100, 0, 1500,  true,  true, "p_{T}(j3) [GeV]", "Norm Events"));

  vh.push_back(PHS("mht",  {PDC("single", "mht",  {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mht,  "Norm Events"));
  vh.push_back(PHS("mt2",  {PDC("single", "best_had_brJet_MT2Zinv", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true, label_mt2, "Norm Events"));
  vh.push_back(PHS("met",  {PDC("single", "met",  {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_met, "Norm Events"));
  vh.push_back(PHS("ht",  {PDC("single", "HT", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 2000,  true,  true,  label_ht,  "Norm Events"));
  vh.push_back(PHS("nMuons",  {PDC("single", "cutMuVec(size)", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 10,  0, 10, true,  false, "N(#mu)", "Events"));
  vh.push_back(PHS("jetPt", {PDC("single", "jetsLVec(pt)",  {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "j p_{T} [GeV]",  "Norm Events"));
  vh.push_back(PHS("nJets", {PDC("single", {{"cntNJetsPt30Eta24", dsDY_ll}, {"cntNJetsPt30Eta24", dsDY_nunu}, {"cntNJetsPt30Eta24Zinv", dsDY_ll}})}, {1, 2}, "", 20,  0, 20, true,  false, label_nj,  "Events"));
  vh.push_back(PHS("muPt", {PDC("single", "cutMuVec(pt)",  {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mupt,  "Norm Events"));
  vh.push_back(PHS("muPt_zpt_gt_600", {PDC("single", "cutMuVec(pt)",  {dsDY_ll, dsDY_nunu})}, {1, 2}, "genZPt>800", 100, 0, 1500,  true,  true,  label_mupt,  "Norm Events"));
  vh.push_back(PHS("muPt_res", {PDC("sinlgle", "genMatchMuInAccRes", {dsDY_ll, dsDY_nunu})}, {1, 1}, "", 100, 0, 1.0, true,  false, "p_{T}(#mu) (gen - reco)/gen",  "Events"));
  vh.push_back(PHS("muMtw", {PDC("single", "muonsMtw", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 600, true,  false, "#mu M_{T}(W) [GeV]",  "Events"));
  vh.push_back(PHS("nMuGenMatch", {PDC("single", "genMatchMuInAcc(size)", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 10,  0, 10, true,  false, "N(#mu) gen Matched",  "Events"));
  vh.push_back(PHS("genZPt",  {PDC("single", "genZPt", {dsDY_ll, dsDY_nunu, dsDY_ll_forEff, dsDY_ll_zAcc})}, {4, 3}, "", 100, 0, 1500,  true,  false, "gen Z p_{T} [GeV]", "Events"));
  vh.push_back(PHS("genZPt2", {PDC("single", "genZPt", {dsDY_ll, dsDY_ll_forEff, dsDY_ll_zAccAcc, dsDY_ll_zAccAcc_nw})}, {3, 2}, "", 100, 0, 1500,  true,  true,  "gen Z p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("genZPt3", {PDC("single", "genZPt", {dsDY_ll_forEff, dsDY_ll_zAccAcc_nw})}, {2, 1}, "", 100, 0, 1500,  true,  false, "gen Z p_{T} [GeV]", "Events"));
  vh.push_back(PHS("genZPt_gr", {PDC("single", "genZPt", {dsDY_ll, dsDY_nunu, dsDY_ll_forEff, dsDY_ll_zAcc})}, {3, 2}, "", 200, 0, 3000,  true,  true,  "gen Z p_{T} [GeV]", "Norm Events"));
  //vh.push_back(PHS("genMuPt", {PDC("single", "genMu(pt)",  {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "gen #mu p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("genMuEta", {PDC("single", "genMu(eta)", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, -3, 3, true,  true,  "gen #mu #eta", "Norm Events"));
  vh.push_back(PHS("genZEta", {PDC("single", "genZEta", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, -3, 3, false, true,  "gen Z #eta", "Norm Events"));
  vh.push_back(PHS("genZmass", {PDC("single", "genZmass", {dsDY_ll, dsDY_nunu})}, {1, 2}, "", 100, 0, 1000,  true,  true,  "gen M(Z) [GeV]",  "Norm Events"));

  vh.push_back(PHS("ttbar_mht",  {PDC("single", "mht", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mht,  "Norm Events"));
  vh.push_back(PHS("ttbar_mt2",  {PDC("single", "best_had_brJet_MT2Zinv", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mt2, "Norm Events"));
  vh.push_back(PHS("ttbar_met",  {PDC("single", "met", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_met, "Norm Events"));
  vh.push_back(PHS("ttbar_cleanMet", {PDC("single", "cleanMetPt",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_met, "Norm Events"));
  vh.push_back(PHS("ttbar_bestZPt", {PDC("single", "bestRecoZPt", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "Z p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("ttbar_ht",  {PDC("single", "HT",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 2000,  true,  true,  label_ht,  "Norm Events"));
  vh.push_back(PHS("ttbar_cleanht", {PDC("single", "HTZinv",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 2000,  true,  true,  label_ht,  "Norm Events"));
  vh.push_back(PHS("ttbar_nMuons",  {PDC("single", "cutMuVec(size)",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 10,  0, 10, true,  false, "N(#mu)", "Events"));
  vh.push_back(PHS("ttbar_jetPt", {PDC("single", "jetsLVec(pt)", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  "j p_{T} [GeV]",  "Norm Events"));
  vh.push_back(PHS("ttbar_cleannJets",  {PDC("single", "cntNJetsPt30Eta24Zinv",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 20,  0, 20, true,  false, label_nj,  "Events"));
  vh.push_back(PHS("ttbar_muPt", {PDC("single", "cutMuVec(pt)", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mupt,  "Norm Events"));
  vh.push_back(PHS("ttbar_mu1Pt", {PDC("single", "cutMuVec[0](pt)", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mupt,  "Norm Events"));
  vh.push_back(PHS("ttbar_mu2Pt", {PDC("single", "cutMuVec[1](pt)", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mupt,  "Norm Events"));
  vh.push_back(PHS("ttbar_elecPt",  {PDC("single", "cutElecVec(pt)",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mupt,  "Norm Events"));
  vh.push_back(PHS("ttbar_elec1Pt", {PDC("single", "cutElecVec[0](pt)",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mupt,  "Norm Events"));
  vh.push_back(PHS("ttbar_elec2Pt", {PDC("single", "cutElecVec[1](pt)",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 100, 0, 1500,  true,  true,  label_mupt,  "Norm Events"));
  vh.push_back(PHS("ttbar_nb",  {PDC("single", "cntCSVSZinv", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 10,  0, 10, true,  false, label_nb, "Events"));
  vh.push_back(PHS("ttbar_nt",  {PDC("single", "nTopCandSortedCnt",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "", 10,  0, 10, true,  false, label_nt, "Events"));

  vh.push_back(PHS("ttbar_baseline_mht",  {PDC("single", "mht", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_mht,  "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_mt2",  {PDC("single", "best_had_brJet_MT2Zinv", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_mt2, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_met",  {PDC("single", "met", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_met, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_cleanMet", {PDC("single", "cleanMetPt",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_met, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_bestZPt", {PDC("single", "bestRecoZPt", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "Z p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_ht",  {PDC("single", "HT",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 2000,  true,  true,  label_ht, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_cleanht", {PDC("single", "HTZinv",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 2000,  true,  true,  label_ht, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_nMuons",  {PDC("single", "cutMuVec(size)",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 10,  0, 10, true,  false, "N(#mu)",  "Events"));
  vh.push_back(PHS("ttbar_baseline_jetPt", {PDC("single", "jetsLVec(pt)", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  "j p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_cleannJets",  {PDC("single", "cntNJetsPt30Eta24Zinv",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 20,  0, 20, true,  false, label_nj,  "Events"));
  vh.push_back(PHS("ttbar_baseline_muPt", {PDC("single", "cutMuVec(pt)", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_mupt, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_mu1Pt", {PDC("single", "cutMuVec[0](pt)", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_mupt, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_mu2Pt", {PDC("single", "cutMuVec[1](pt)", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_mupt, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_elecPt",  {PDC("single", "cutElecVec(pt)",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_mupt, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_elec1Pt", {PDC("single", "cutElecVec[0](pt)",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_mupt, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_elec2Pt", {PDC("single", "cutElecVec[1](pt)",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 100, 0, 1500,  true,  true,  label_mupt, "Norm Events"));
  vh.push_back(PHS("ttbar_baseline_nb",  {PDC("single", "cntCSVSZinv", {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 10,  0, 10, true,  false, label_nb, "Events"));
  vh.push_back(PHS("ttbar_baseline_nt",  {PDC("single", "nTopCandSortedCnt",  {dsTT_mm_lepSel, dsTT_ee_lepSel, dsTT_em_lepSel})}, {1, 2}, "passBaselineZinv", 10,  0, 10, true,  false, label_nt, "Events"));

  // Simple Reco vs. Gen plots
  vh.push_back(PHS("bestRecoZPt"  ,  {PDC("single",  {{"bestRecoZPt", dsDY_ll}, {"genZPt", dsDY_ll_gen}})},  {1, 2}, "", 100, 0, 1500,  true,  true,  "Z p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("bestRecoZMass_nosel", {PDC("single",  {{"bestRecoZM", dsDY_ll_noSel},  {"genZmass", dsDY_ll_gen_noSel}})},  {1, 2}, "", 100, 0, 1000,  true,  true,  "M(Z) [GeV]", "Norm Events"));
  vh.push_back(PHS("bestRecoZMass",  {PDC("single",  {{"bestRecoZM", dsDY_ll}, {"genZmass", dsDY_ll_gen}})}, {1, 2}, "", 100, 0, 1000,  true,  true,  "M(Z) [GeV]", "Norm Events"));
  vh.push_back(PHS("bestRecoZMass_dPhiCut", {PDC("single",  {{"bestRecoZM", dsDY_ll}, {"genZmass", dsDY_ll_gen}})}, {1, 2}, "dPhi0_CUT<0.5", 100, 0, 1000,  true,  true,  "M(Z) [GeV]", "Norm Events"));
  vh.push_back(PHS("zPtRes",  {PDC("single", "ZPtRes", {dsDY_ll, dsDY_ll_base})}, {1, 1}, "", 200,  -0.5,  0.5,  false,  true,  "Z p_{T} (reco - gen)/gen", "Norm Events"));
  vh.push_back(PHS("zPtRes_Zpt_lt_200",  {PDC("single", "ZPtRes", {dsDY_ll, dsDY_ll_base})}, {1, 1}, "genZPt<200", 200,  -0.5,  0.5,  false,  true,  "Z p_{T} (reco - gen)/gen", "Norm Events"));
  vh.push_back(PHS("zPtRes_Zpt_gt_200",  {PDC("single", "ZPtRes", {dsDY_ll, dsDY_ll_base})}, {1, 1}, "genZPt>200", 200,  -0.5,  0.5,  false,  true,  "Z p_{T} (reco - gen)/gen", "Norm Events"));
  vh.push_back(PHS("zEtaRes", {PDC("single", "ZEtaRes",  {dsDY_ll, dsDY_ll_base})}, {1, 1}, "", 200,  -2,  2,  false,  true,  "Z #eta (reco - gen)", "Norm Events"));
  vh.push_back(PHS("zPhiRes", {PDC("single", "ZPhiRes",  {dsDY_ll, dsDY_ll_base})}, {1, 1}, "", 200,  -2,  2,  false,  true,  "Z #phi (reco - gen)", "Norm Events"));
  vh.push_back(PHS("zMRes", {PDC("single", "ZMRes", {dsDY_ll, dsDY_ll_base})}, {1, 1}, "", 200,  -0.5,  0.5,  false,  true,  "M(Z) (reco - gen)/gen",  "Norm Events"));

  vh.push_back(PHS("mu1dRMin", {PDC("single", "mu1dRMin", {dsDY_ll, dsDY_ll_base})}, {1, 1}, "", 100, 0,  5,  false,  true,  "min dR (#mu1,jet)", "Norm Events"));
  vh.push_back(PHS("mu2dRMin", {PDC("single", "mu2dRMin", {dsDY_ll, dsDY_ll_base})}, {1, 1}, "", 100, 0,  5,  false,  true,  "min dR (#mu2,jet)", "Norm Events"));
  vh.push_back(PHS("mudR", {PDC("single", "mudR", {dsDY_ll, dsDY_ll_base})}, {1, 1}, "", 100, 0,  5,  false,  true,  "min dR (#mu1,#mu2)", "Norm Events"));
  vh.push_back(PHS("genMudR", {PDC("single", "genMudR",  {dsDY_ll_forEff, dsDY_ll_base})}, {1, 1}, "", 100, 0,  5,  false,  true,  "gen min dR (#mu1,#mu2)", "Norm Events"));
  vh.push_back(PHS("genMudR_Zpt_gt_800", {PDC("single", "genMudR",  {dsDY_ll_forEff, dsDY_ll_base})}, {1, 1}, "genZPt>800", 100, 0,  5,  false,  true,  "gen min dR (#mu1,#mu2)", "Norm Events"));
  vh.push_back(PHS("genMudPhi",  {PDC("single", "genMudPhi", {dsDY_ll_forEff, dsDY_ll_base})}, {1, 1}, "", 100,  -4,  4,  false,  true,  "gen min dPhi (#mu1,#mu2)", "Norm Events"));
  vh.push_back(PHS("genMudEta",  {PDC("single", "genMudEta", {dsDY_ll_forEff, dsDY_ll_base})}, {1, 1}, "", 100,  -5,  5,  false,  true,  "gen min dEta (#mu1,#mu2)", "Norm Events"));

  // Here are the interesting plots for the closure test

  // DataCollections which are reused multiple times are defined "upfront" to save time
  Plotter::DataCollection cleanht( "single", {{"HTZinv", dsDY_nunu},  {"HTZinv", dsDY_ll_zAcc}, {"HTZinv", dsDY_ll_trig_scaled} });
  Plotter::DataCollection cleanmht( "single", {{"cleanMHt", dsDY_nunu}, {"cleanMHt", dsDY_ll_zAcc}, {"cleanMHt", dsDY_ll_trig_scaled} });
  Plotter::DataCollection cleanmhtphi( "single", {{"cleanMHtPhi", dsDY_nunu}, {"cleanMHtPhi", dsDY_ll_zAcc}, {"cleanMHtPhi", dsDY_ll_trig_scaled} });
  Plotter::DataCollection cleanMet( "single", {{"cleanMetPt", dsDY_nunu},  {"cleanMetPt", dsDY_ll_zAcc}, {"cleanMetPt", dsDY_ll_trig_scaled} });
  Plotter::DataCollection nCleanJet( "single", {{"cntNJetsPt30Eta24", dsDY_nunu},  {"cntNJetsPt30Eta24Zinv", dsDY_ll_zAcc}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_trig_scaled} });
  Plotter::DataCollection nTop( "single", {{"nTopCandSortedCnt", dsDY_nunu},  {"nTopCandSortedCntZinv", dsDY_ll_zAcc}, {"nTopCandSortedCntZinv", dsDY_ll_trig_scaled} });
  Plotter::DataCollection nBottom( "single", {{"cntCSVS", dsDY_nunu}, {"cntCSVSZinv", dsDY_ll_zAcc}, {"cntCSVSZinv", dsDY_ll_trig_scaled} });

  Plotter::DataCollection scaled_ht( "single", "ht",  {{dsDY_nunu}, {dsDY_ll_zWeight_scaled}, {dsDY_ll_zAcc_scaled} });
  Plotter::DataCollection scaled_mht(  "single", "mht", {{dsDY_nunu}, {dsDY_ll_zWeight_scaled}, {dsDY_ll_zAcc_scaled} });
  Plotter::DataCollection scaled_met(  "single", "met", {{dsDY_nunu}, {dsDY_ll_zWeight_scaled}, {dsDY_ll_zAcc_scaled} });

  Plotter::DataCollection scaled_cleanht( "single", {{"HTZinv", dsDY_nunu},  {"HTZinv", dsDY_ll_zAcc_scaled}  });
  Plotter::DataCollection scaled_cleanmt2( "single", {{"best_had_brJet_MT2Zinv", dsDY_nunu}, {"best_had_brJet_MT2Zinv", dsDY_ll_zAcc_scaled} });
  Plotter::DataCollection scaled_cleanmht( "single", {{"cleanMHt", dsDY_nunu}, {"cleanMHt", dsDY_ll_zAcc_scaled} });
  Plotter::DataCollection scaled_cleanmhtphi( "single", {{"cleanMHtPhi", dsDY_nunu}, {"cleanMHtPhi", dsDY_ll_zAcc_scaled} });
  Plotter::DataCollection scaled_cleanMet( "single", {{"cleanMetPt", dsDY_nunu},  {"cleanMetPt", dsDY_ll_zAcc_scaled}  });
  Plotter::DataCollection scaled_nCleanJet( "single", {{"cntNJetsPt30Eta24Zinv", dsDY_nunu},  {"cntNJetsPt30Eta24Zinv", dsDY_ll_zAcc_scaled}  });
  Plotter::DataCollection scaled_nTop( "single", {{"nTopCandSortedCntZinv", dsDY_nunu},  {"nTopCandSortedCntZinv", dsDY_ll_zAcc_scaled}  });
  Plotter::DataCollection scaled_nBottom( "single", {{"cntCSVSZinv", dsDY_nunu}, {"cntCSVSZinv", dsDY_ll_zAcc_scaled} });

  Plotter::DataCollection scaled_elec_ht( "single", "ht",  {{dsDY_nunu}, {dsDY_ee_zWeight_scaled}, {dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_mht(  "single", "mht", {{dsDY_nunu}, {dsDY_ee_zWeight_scaled}, {dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_met(  "single", "met", {{dsDY_nunu}, {dsDY_ee_zWeight_scaled}, {dsDY_ee_zAcc_scaled} });

  Plotter::DataCollection scaled_elec_cleanht( "single", {{"HTZinv",  dsDY_nunu},  {"HTZinv", dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_cleanmt2( "single", {{"best_had_brJet_MT2Zinv", dsDY_nunu}, {"best_had_brJet_MT2Zinv", dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_cleanmht( "single", {{"cleanMHt", dsDY_nunu}, {"cleanMHt", dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_cleanmhtphi( "single", {{"cleanMHtPhi", dsDY_nunu}, {"cleanMHtPhi", dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_cleanMet( "single", {{"cleanMetPt",  dsDY_nunu},  {"cleanMetPt", dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_nCleanJet( "single", {{"cntNJetsPt30Eta24Zinv",  dsDY_nunu},  {"cntNJetsPt30Eta24Zinv", dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_nTop( "single", {{"nTopCandSortedCntZinv",  dsDY_nunu},  {"nTopCandSortedCntZinv", dsDY_ee_zAcc_scaled} });
  Plotter::DataCollection scaled_elec_nBottom( "single", {{"cntCSVSZinv", dsDY_nunu}, {"cntCSVSZinv", dsDY_ee_zAcc_scaled} });

  Plotter::DataCollection scaled_nSearchBin( "single", {{"nSearchBin", dsDY_nunu},  {"nSearchBin", dsDY_ll_zAcc_scaled}  });
  Plotter::DataCollection scaled_nSearchBin_stat("single", {{"nSearchBin", dsDY_ll}, {"nb0BinsNW",  dsDY_ll},  {"nb0Bins",  dsDY_ll} });
  Plotter::DataCollection scaled_nSearchBinNb0(  "single", {{"nSearchBin", dsDY_nunu},  {"nb0Bins", dsDY_ll_zAcc_scaled}  });
  Plotter::DataCollection scaled_nSearchBinNb0ZZ("single", {{"nSearchBin", dsDY_nunu_SBMb}, {"nb0Bins", dsDY_nunu_SB0b} });

  Plotter::DataCollection scaled_Elec_nSearchBin( "single", {{"nSearchBin", dsDY_nunu},  {"nSearchBin", dsDY_ee_zAcc_scaled}  });
  Plotter::DataCollection scaled_Elec_nSearchBin_stat("single", {{"nSearchBin", dsDY_ee}, {"nb0BinsNW",  dsDY_ee} });
  Plotter::DataCollection scaled_Elec_nSearchBinNb0(  "single", {{"nSearchBin", dsDY_nunu},  {"nb0Bins", dsDY_ee_zAcc_scaled}  });
  Plotter::DataCollection scaled_Elec_nSearchBinNb0ZZ("single", {{"nSearchBin", dsDY_nunu_SBMb}, {"nb0Bins", dsDY_nunu_SB0b} });

  Plotter::DataCollection scaled_nSearchBin_Sig_mm("stack", {{"nSearchBin", dsDY_ll}, {"nSearchBin",  dsTT_mm_lepSel} });
  Plotter::DataCollection scaled_nSearchBin_Con_mm("stack", {{"nb0BinsNW", dsDY_ll},  {"nb0BinsNW", dsTT_mm_lepSel} });
  Plotter::DataCollection scaled_nSearchBin_Sig_ee("stack", {{"nSearchBin", dsDY_ee}, {"nSearchBin",  dsTT_ee_lepSel} });
  Plotter::DataCollection scaled_nSearchBin_Con_ee("stack", {{"nb0BinsNW", dsDY_ee},  {"nb0BinsNW", dsTT_ee_lepSel} });


  PDC trg_cleanht(  "single","HTZinv", {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
  PDC trg_cleanmht( "single","cleanMHt",  {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
  PDC trg_cleanMet( "single","cleanMetPt", {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
  PDC trg_nCleanJet("single","cntNJetsPt30Eta24Zinv",{dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
  PDC trg_nTop( "single","nTopCandSortedCntZinv",{dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});
  PDC trg_nBottom(  "single","cntCSVSZinv",  {dsDY_ll_zAcc, dsDY_ll_mu30_mu20, dsDY_ll_mu30_mu30, dsDY_ll_mu40_mu30, dsDY_ll_mu45_mu20, dsDY_ll_singleMu45});

  Plotter::DataCollection scaled_stacked_DYtoll_genht(  "stack", "genHt", {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});

  Plotter::DataCollection scaled_stacked_DYtoll_cleanht(  "stack", "HTZinv", {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
  Plotter::DataCollection scaled_stacked_DYtoll_ht( "stack", "ht",  {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
  Plotter::DataCollection scaled_stacked_DYtoll_cleanmht( "stack", "cleanMHt", {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});
  Plotter::DataCollection scaled_stacked_DYtoll_cleanMet( "stack", "cleanMetPt", {dsDY_ll_HT_100to200, dsDY_ll_HT_200to400, dsDY_ll_HT_400to600, dsDY_ll_HT_600toInf});

  Plotter::DataCollection cleanht_znunu( "single", {{"HTZinv", dsDY_nunu}});
  Plotter::DataCollection ht_znunu(  "single", {{"ht",  dsDY_nunu}});
  Plotter::DataCollection cleanmht_znunu( "single", {{"cleanMHt", dsDY_nunu}});
  Plotter::DataCollection mht_znunu(  "single", {{"mht", dsDY_nunu}});
  Plotter::DataCollection cleanMet_znunu( "single", {{"cleanMetPt", dsDY_nunu}});
  Plotter::DataCollection met_znunu(  "single", {{"met", dsDY_nunu}});

  Plotter::DataCollection scaled_stacked_DYtonunu_cleanht(  "stack", "HTZinv", {dsDY_nunu_HT_100to200, dsDY_nunu_HT_200to400, dsDY_nunu_HT_400to600, dsDY_nunu_HT_600toInf});
  Plotter::DataCollection scaled_stacked_DYtonunu_cleanmht( "stack", "cleanMHt", {dsDY_nunu_HT_100to200, dsDY_nunu_HT_200to400, dsDY_nunu_HT_400to600, dsDY_nunu_HT_600toInf});
  Plotter::DataCollection scaled_stacked_DYtonunu_cleanMet( "stack", "cleanMetPt", {dsDY_nunu_HT_100to200, dsDY_nunu_HT_200to400, dsDY_nunu_HT_400to600, dsDY_nunu_HT_600toInf});

  vh.push_back(PHS("trigStudy_mu1Pt", {PDC("single", "cutMuPt1",  {dsDY_ll_baseNT, dsDY_ll})}, {2, 1}, "", 100, 0, 1500,  true,  false,  label_mupt, "Events"));
  vh.push_back(PHS("trigStudy_mu2Pt", {PDC("single", "cutMuPt2",  {dsDY_ll_baseNT, dsDY_ll})}, {2, 1}, "", 100, 0, 1500,  true,  false,  label_mupt, "Events"));

  vh.push_back(PHS("trigStudy_cleanht", {trg_cleanht  }, {1, 1}, "", 50,  0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("trigStudy_cleanmht",  {trg_cleanmht }, {1, 1}, "", 50,  0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("trigStudy_cleanmet",  {trg_cleanMet }, {1, 1}, "", 50,  0, 1500, true,  false,  label_met,  "Events"));
  vh.push_back(PHS("trigStudy_nTop",  {trg_nTop }, {1, 1}, "", 10,  0, 10, true,  false,  label_nt, "Events"));
  vh.push_back(PHS("trigStudy_cleannJet", {trg_nCleanJet}, {1, 1}, "", 20,  0, 20, true,  false,  label_nj,  "Events"));
  vh.push_back(PHS("trigStudy_nBottom", {trg_nBottom  }, {1, 1}, "", 10,  0, 10, true,  false,  label_nb, "Events"));

  vh.push_back(PHS("trigStudy_baseline_cleanht", {trg_cleanht  }, {1, 1}, "passBaselineZinv", 50,  0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("trigStudy_baseline_cleanmht",  {trg_cleanmht }, {1, 1}, "passBaselineZinv", 50,  0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("trigStudy_baseline_cleanmet",  {trg_cleanMet }, {1, 1}, "passBaselineZinv", 50,  0, 1500, true,  false,  label_met,  "Events"));
  vh.push_back(PHS("trigStudy_baseline_nTop",  {trg_nTop }, {1, 1}, "passBaselineZinv", 10,  0, 10, true,  false,  label_nt, "Events"));
  vh.push_back(PHS("trigStudy_baseline_cleannJet", {trg_nCleanJet}, {1, 1}, "passBaselineZinv", 20,  0, 20, true,  false,  label_nj,  "Events"));
  vh.push_back(PHS("trigStudy_baseline_nBottom", {trg_nBottom  }, {1, 1}, "passBaselineZinv", 10,  0, 10, true,  false,  label_nb, "Events"));

  vh.push_back(PHS("trigStudy_baselineNoTag_cleanht", {trg_cleanht  }, {1, 1}, "passBaselineNoTagZinv", 50,  0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("trigStudy_baselineNoTag_cleanmht",  {trg_cleanmht }, {1, 1}, "passBaselineNoTagZinv", 50,  0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("trigStudy_baselineNoTag_cleanmet",  {trg_cleanMet }, {1, 1}, "passBaselineNoTagZinv", 50,  0, 1500, true,  false,  label_met,  "Events"));
  vh.push_back(PHS("trigStudy_baselineNoTag_nTop",  {trg_nTop }, {1, 1}, "passBaselineNoTagZinv", 10,  0, 10, true,  false,  label_nt, "Events"));
  vh.push_back(PHS("trigStudy_baselineNoTag_cleannJet", {trg_nCleanJet}, {1, 1}, "passBaselineNoTagZinv", 20,  0, 20, true,  false,  label_nj,  "Events"));
  vh.push_back(PHS("trigStudy_baselineNoTag_nBottom", {trg_nBottom  }, {1, 1}, "passBaselineNoTagZinv", 10,  0, 10, true,  false,  label_nb, "Events"));

  vh.push_back(PHS("trigRatio_singleMu45_cleanht", {PDC("ratio", "HTZinv", {dsDY_ll_singleMu45, dsDY_ll_zAcc})}, {1, 1}, "passBaselineNoTagZinv",  50,  0, 2000, false,  false,  label_ht, "Ratio"));
  vh.push_back(PHS("trigRatio_singleMu45_cleanmht",  {PDC("ratio", "cleanMHt", {dsDY_ll_singleMu45, dsDY_ll_zAcc})}, {1, 1}, "passBaselineNoTagZinv",  50,  0, 1500, false,  false,  label_mht, "Ratio"));
  vh.push_back(PHS("trigRatio_singleMu45_cleanmet",  {PDC("ratio", "cleanMetPt", {dsDY_ll_singleMu45, dsDY_ll_zAcc})}, {1, 1}, "passBaselineNoTagZinv",  50,  0, 1500, false,  false,  label_met,  "Ratio"));
  vh.push_back(PHS("trigRatio_singleMu45_nTop",  {PDC("ratio", "nTopCandSortedCntZinv", {dsDY_ll_singleMu45, dsDY_ll_zAcc})}, {1, 1}, "passBaselineNoTagZinv",  10,  0, 10, false,  false,  label_nt, "Ratio"));
  vh.push_back(PHS("trigRatio_singleMu45_cleannJet", {PDC("ratio", "cntNJetsPt30Eta24Zinv", {dsDY_ll_singleMu45, dsDY_ll_zAcc})}, {1, 1}, "passBaselineNoTagZinv",  20,  0, 20, false,  false,  label_nj,  "Ratio"));
  vh.push_back(PHS("trigRatio_singleMu45_nBottom", {PDC("ratio", "cntCSVSZinv", {dsDY_ll_singleMu45, dsDY_ll_zAcc})}, {1, 1}, "passBaselineNoTagZinv",  10,  0, 10, false,  false,  label_nb, "Ratio"));

  vh.push_back(PHS("genht_DY_mumu_stack", {scaled_stacked_DYtoll_genht}, {1, 1}, "",  100, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("cleanht_DY_mumu_stack_er",  {scaled_stacked_DYtoll_cleanht,  cleanht_znunu }, {1, 1}, "",  150, 0, 3000, true,  false,  label_ht, "Events"));

  vh.push_back(PHS("cleanht_DY_mumu_stack", {scaled_stacked_DYtoll_cleanht,  cleanht_znunu }, {1, 1}, "", 100, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("ht_DY_mumu_stack", {scaled_stacked_DYtoll_ht, ht_znunu  }, {1, 1}, "", 100, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("cleanmht_DY_mumu_stack",  {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu}, {1, 1}, "", 100, 0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("cleanmet_DY_mumu_stack",  {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu}, {1, 1}, "", 100, 0, 1500, true,  false,  label_met,  "Events"));

  vh.push_back(PHS("baseline_cleanht_DY_mumu_stack",  {scaled_stacked_DYtoll_cleanht,  cleanht_znunu }, {1, 1}, "passBaselineZinv",  100, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("baseline_cleanmht_DY_mumu_stack", {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu}, {1, 1}, "passBaselineZinv",  100, 0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("baseline_cleanmet_DY_mumu_stack", {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu}, {1, 1}, "passBaselineZinv",  100, 0, 1500, true,  false,  label_met,  "Events"));

  vh.push_back(PHS("nb_0_cleanht_DY_mumu_stack",  {scaled_stacked_DYtoll_cleanht,  cleanht_znunu }, {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("nb_0_cleanmht_DY_mumu_stack", {scaled_stacked_DYtoll_cleanmht, cleanmht_znunu}, {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("nb_0_cleanmet_DY_mumu_stack", {scaled_stacked_DYtoll_cleanMet, cleanMet_znunu}, {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0, 1500, true,  false,  label_met,  "Events"));

  vh.push_back(PHS("cleanht_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanht}, {1, 1}, "", 100, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("cleanmht_DY_nunu_stack",  {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "", 100, 0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("cleanmet_DY_nunu_stack",  {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "", 100, 0, 1500, true,  false,  label_met,  "Events"));

  vh.push_back(PHS("baseline_cleanht_DY_nunu_stack",  {scaled_stacked_DYtonunu_cleanht}, {1, 1}, "passBaselineZinv",  100, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("baseline_cleanmht_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "passBaselineZinv",  100, 0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("baseline_cleanmet_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "passBaselineZinv",  100, 0, 1500, true,  false,  label_met,  "Events"));

  vh.push_back(PHS("nb_0_cleanht_DY_nunu_stack",  {scaled_stacked_DYtonunu_cleanht}, {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("nb_0_cleanmht_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanmht},  {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0, 1500, true,  false,  label_mht, "Events"));
  vh.push_back(PHS("nb_0_cleanmet_DY_nunu_stack", {scaled_stacked_DYtonunu_cleanMet},  {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0",  100, 0, 1500, true,  false,  label_met,  "Events"));

  vh.push_back(PHS("nSearchBin", {scaled_nSearchBin},  {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_log", {scaled_nSearchBin},  {2, 1}, "passBaselineZinv", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_stat",  {scaled_nSearchBin_stat}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_stat_log",  {scaled_nSearchBin_stat}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBinnb0", {scaled_nSearchBinNb0}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBinnb0_log", {scaled_nSearchBinNb0}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBinnb0ZZ",  {scaled_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBinnb0ZZ_log",  {scaled_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_pull",  {scaled_nSearchBin},  {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", false));
  vh.push_back(PHS("nSearchBin_pull_log",  {scaled_nSearchBin},  {2, 1}, "passBaselineZinv", NSB,  0, NSB, true,  false,  "Search Bin", "Events", false));
  vh.push_back(PHS("nSearchBinnb0_pull", {scaled_nSearchBinNb0}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", false));
  vh.push_back(PHS("nSearchBinnb0_pull_log", {scaled_nSearchBinNb0}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", false));
  vh.push_back(PHS("nSearchBinnb0ZZ_pull", {scaled_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", false));
  vh.push_back(PHS("nSearchBinnb0ZZ_pull_log", {scaled_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", false));

  vh.push_back(PHS("elec_nSearchBin", {scaled_Elec_nSearchBin},  {2, 1}, "passBaselineZinv",  NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("elec_nSearchBin_log", {scaled_Elec_nSearchBin},  {2, 1}, "passBaselineZinv",  NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("elec_nSearchBin_stat",  {scaled_Elec_nSearchBin_stat}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("elec_nSearchBin_stat_log",  {scaled_Elec_nSearchBin_stat}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("elec_nSearchBinnb0", {scaled_Elec_nSearchBinNb0}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("elec_nSearchBinnb0_log", {scaled_Elec_nSearchBinNb0}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("elec_nSearchBinnb0ZZ",  {scaled_Elec_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("elec_nSearchBinnb0ZZ_log",  {scaled_Elec_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("elec_nSearchBin_pull",  {scaled_Elec_nSearchBin},  {2, 1}, "passBaselineZinv",  NSB,  0, NSB, false, false,  "Search Bin", "Events", false));
  vh.push_back(PHS("elec_nSearchBin_pull_log",  {scaled_Elec_nSearchBin},  {2, 1}, "passBaselineZinv",  NSB,  0, NSB, true,  false,  "Search Bin", "Events", false));
  vh.push_back(PHS("elec_nSearchBinnb0_pull", {scaled_Elec_nSearchBinNb0}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", false));
  vh.push_back(PHS("elec_nSearchBinnb0_pull_log", {scaled_Elec_nSearchBinNb0}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", false));
  vh.push_back(PHS("elec_nSearchBinnb0ZZ_pull", {scaled_Elec_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", false));
  vh.push_back(PHS("elec_nSearchBinnb0ZZ_pull_log", {scaled_Elec_nSearchBinNb0ZZ}, {2, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", false));

  vh.push_back(PHS("nSearchBin_Sig_mm",  {scaled_nSearchBin_Sig_mm}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_Sig_mm_log",  {scaled_nSearchBin_Sig_mm}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_Con_mm",  {scaled_nSearchBin_Con_mm}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_Con_mm_log",  {scaled_nSearchBin_Con_mm}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_Sig_ee",  {scaled_nSearchBin_Sig_ee}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_Sig_ee_log",  {scaled_nSearchBin_Sig_ee}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_Con_ee",  {scaled_nSearchBin_Con_ee}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("nSearchBin_Con_ee_log",  {scaled_nSearchBin_Con_ee}, {1, 1}, "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));

  std::vector<std::pair<std::pair<std::string,std::string>, int>> cutlevels_MCCLosure = {
    {{"", ""  },  100},
    {{"baseline", "passBaselineZinv"  },  50},
    {{"baselineNoTag", "passBaselineNoTagZinv" },  50},
    {{"nb_0", "passBaselineNoTagZinv;cntCSVSZinv=0" },  100},
    {{"nb_1", "passBaselineNoTagZinv;cntCSVSZinv=1" },  50},
    {{"nb_2", "passBaselineNoTagZinv;cntCSVSZinv=2" },  25},
    {{"nb_3", "passBaselineNoTagZinv;cntCSVSZinv>2" },  15},
    {{"nt_0", "passBaselineNoTagMT2Zinv;nTopCandSortedCntZinv=0"  },  100},
    {{"nt_1", "passBaselineNoTagZinv;nTopCandSortedCntZinv=1" },  50},
    {{"nt_2", "passBaselineNoTagZinv;nTopCandSortedCntZinv=2" },  25},
    {{"nt_3", "passBaselineNoTagZinv;nTopCandSortedCntZinv>2" },  15},
    {{"nb_0_nt_0", "passBaselineNoTagMT2Zinv;cntCSVSZinv=0;nTopCandSortedCntZinv=0" },  50},
    {{"nb_1_nt_0", "passBaselineNoTagMT2Zinv;cntCSVSZinv=1;nTopCandSortedCntZinv=0" },  25},
    {{"nb_2_nt_0", "passBaselineNoTagMT2Zinv;cntCSVSZinv=2;nTopCandSortedCntZinv=0" },  15},
    {{"nb_3_nt_0", "passBaselineNoTagMT2Zinv;cntCSVSZinv>2;nTopCandSortedCntZinv=0" },  10},
    {{"nb_0_nt_1", "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=1"  },  50},
    {{"nb_1_nt_1", "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=1"  },  25},
    {{"nb_2_nt_1", "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=1"  },  15},
    {{"nb_3_nt_1", "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=1"  },  10},
    {{"nb_0_nt_2", "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv=2"  },  50},
    {{"nb_1_nt_2", "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv=2"  },  25},
    {{"nb_2_nt_2", "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv=2"  },  15},
    {{"nb_3_nt_2", "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv=2"  },  10},
    {{"nb_0_nt_3", "passBaselineNoTagZinv;cntCSVSZinv=0;nTopCandSortedCntZinv>2"  },  50},
    {{"nb_1_nt_3", "passBaselineNoTagZinv;cntCSVSZinv=1;nTopCandSortedCntZinv>2"  },  25},
    {{"nb_2_nt_3", "passBaselineNoTagZinv;cntCSVSZinv=2;nTopCandSortedCntZinv>2"  },  15},
    {{"nb_3_nt_3", "passBaselineNoTagZinv;cntCSVSZinv>2;nTopCandSortedCntZinv>2"  },  10},
    {{"loose0",  "HTZinv>200;passnJetsZinv;passMETZinv;passdPhisZinv"  },  80},
    {{"loose0_MT2",  "HTZinv>200;passnJetsZinv;passMETZinv;passdPhisZinv;best_had_brJet_MT2Zinv>0"  },  80}
  };

  for(auto& cut : cutlevels_MCCLosure)
    {
      //muons
      vh.push_back(PHS("MCClosure_muon_" + cut.first.first + "_mT2Zinv", {scaled_cleanmt2}, {2, 1}, cut.first.second, cut.second, 0, 2000, true,  false,  label_mt2, "Events"));
      vh.push_back(PHS("MCClosure_muon_" + cut.first.first + "_cleanht", {scaled_cleanht}, {2, 1}, cut.first.second, cut.second, 0, 2000, true,  false,  label_ht, "Events"));
      vh.push_back(PHS("MCClosure_muon_" + cut.first.first + "_cleanmht", {scaled_cleanmht}, {2, 1}, cut.first.second, cut.second, 0, 2000, true,  false,  label_mht, "Events"));
      vh.push_back(PHS("MCClosure_muon_" + cut.first.first + "_cleanmet", {scaled_cleanMet}, {2, 1}, cut.first.second, cut.second, 0, 2000, true,  false,  label_met,  "Events"));
      vh.push_back(PHS("MCClosure_muon_" + cut.first.first + "_nTop", {scaled_nTop}, {2, 1}, cut.first.second, 10, 0, 10, true,  false,  label_nt, "Events"));
      vh.push_back(PHS("MCClosure_muon_" + cut.first.first + "_cleannJet",  {scaled_nCleanJet},  {2, 1}, cut.first.second, 20, 0, 20, true,  false,  label_nj,  "Events"));
      vh.push_back(PHS("MCClosure_muon_" + cut.first.first + "_nBottom", {scaled_nBottom}, {2, 1}, cut.first.second, 10, 0, 10, true,  false,  label_nb, "Events"));

      //electron
      vh.push_back(PHS("MCClosure_elec_" + cut.first.first + "_mT2Zinv", {scaled_elec_cleanmt2}, {2, 1}, cut.first.second, cut.second, 0, 2000, true,  false,  label_mt2, "Events"));
      vh.push_back(PHS("MCClosure_elec_" + cut.first.first + "_cleanht", {scaled_elec_cleanht}, {2, 1}, cut.first.second, cut.second, 0, 2000, true,  false,  label_ht, "Events"));
      vh.push_back(PHS("MCClosure_elec_" + cut.first.first + "_cleanmht", {scaled_elec_cleanmht}, {2, 1}, cut.first.second, cut.second, 0, 2000, true,  false,  label_mht, "Events"));
      vh.push_back(PHS("MCClosure_elec_" + cut.first.first + "_cleanmet", {scaled_elec_cleanMet}, {2, 1}, cut.first.second, cut.second, 0, 2000, true,  false,  label_met,  "Events"));
      vh.push_back(PHS("MCClosure_elec_" + cut.first.first + "_nTop", {scaled_elec_nTop}, {2, 1}, cut.first.second, 10, 0, 10, true,  false,  label_nt, "Events"));
      vh.push_back(PHS("MCClosure_elec_" + cut.first.first + "_cleannJet", {scaled_elec_nCleanJet},  {2, 1}, cut.first.second, 20, 0, 20, true,  false,  label_nj,  "Events"));
      vh.push_back(PHS("MCClosure_elec_" + cut.first.first + "_nBottom", {scaled_elec_nBottom}, {2, 1}, cut.first.second, 10, 0, 10, true,  false,  label_nb, "Events"));
    }


  vh.push_back(PHS("fake1b_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_1b},  {"nTopCandSortedCntZinv1b", dsDY_nunu_0b_1fake}})},  {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));
  vh.push_back(PHS("fake1b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_1b_blNoTag}, {"nTopCandSortedCntZinv1b", dsDY_nunu_0b_1fake_blNoTag}})}, {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));

  vh.push_back(PHS("fake2b_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_2b},  {"nTopCandSortedCntZinv2b", dsDY_nunu_0b_2fake}})},  {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));
  vh.push_back(PHS("fake2b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_2b_blNoTag}, {"nTopCandSortedCntZinv2b", dsDY_nunu_0b_2fake_blNoTag}})}, {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));

  vh.push_back(PHS("fake3b_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_3b},  {"nTopCandSortedCntZinv3b", dsDY_nunu_0b_3fake}})},  {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));
  vh.push_back(PHS("fake3b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_3b_blNoTag}, {"nTopCandSortedCntZinv3b", dsDY_nunu_0b_3fake_blNoTag}})}, {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));

  vh.push_back(PHS("fake2bvs0b_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_0b},{"nTopCandSortedCntZinv2b", dsDY_nunu_0b_2fake}})}, {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));
  vh.push_back(PHS("fake2bvs0b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_0b},{"nTopCandSortedCntZinv2b", dsDY_nunu_0b_2fake}})}, {2, 1}, "passBaselineNoTagZinv", 10, 0, 10, true, true, label_nt, "Norm Events"));

  vh.push_back(PHS("fake3bvs0b_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_0b},{"nTopCandSortedCntZinv3b", dsDY_nunu_0b_3fake}})}, {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));
  vh.push_back(PHS("fake3bvs0b_baselineNoTag_nTop", {PDC("single", {{"nTopCandSortedCntZinv", dsDY_nunu_0b},{"nTopCandSortedCntZinv3b", dsDY_nunu_0b_3fake}})}, {2, 1}, "passBaselineNoTagZinv", 10, 0, 10, true, true, label_nt, "Norm Events"));

  vh.push_back(PHS("fake1b_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_1b},  {"best_had_brJet_MT2Zinv1b", dsDY_nunu_0b_1fake}})},  {2, 1}, "", 50, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("fake1b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_1b_blNoTag}, {"best_had_brJet_MT2Zinv1b", dsDY_nunu_0b_1fake_blNoTag}})}, {2, 1}, "", 50, 0, 2000, true, true, label_mt2, "Norm Events"));

  vh.push_back(PHS("fake2b_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_2b},  {"best_had_brJet_MT2Zinv2b", dsDY_nunu_0b_2fake}})},  {2, 1}, "", 50, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("fake2b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_2b_blNoTag}, {"best_had_brJet_MT2Zinv2b", dsDY_nunu_0b_2fake_blNoTag}})}, {2, 1}, "", 50, 0, 2000, true, true, label_mt2, "Norm Events"));

  vh.push_back(PHS("fake3b_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_3b},  {"best_had_brJet_MT2Zinv3b", dsDY_nunu_0b_3fake}})},  {2, 1}, "", 50, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("fake3b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv", dsDY_nunu_3b_blNoTag}, {"best_had_brJet_MT2Zinv3b", dsDY_nunu_0b_3fake_blNoTag}})}, {2, 1}, "", 50, 0, 2000, true, true, label_mt2, "Norm Events"));

  vh.push_back(PHS("fake2bvs0b_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv",dsDY_nunu_0b},{"best_had_brJet_MT2Zinv2b",dsDY_nunu_0b_2fake}})}, {2, 1}, "", 50, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("fake2bvs0b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv",dsDY_nunu_0b},{"best_had_brJet_MT2Zinv2b",dsDY_nunu_0b_2fake}})}, {2, 1}, "passBaselineNoTagZinv", 50, 0, 2000, true, true, label_mt2, "Norm Events"));

  vh.push_back(PHS("fake3bvs0b_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv",dsDY_nunu_0b},{"best_had_brJet_MT2Zinv3b",dsDY_nunu_0b_3fake}})}, {2, 1}, "", 50, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("fake3bvs0b_baselineNoTag_MT2", {PDC("single", {{"best_had_brJet_MT2Zinv",dsDY_nunu_0b},{"best_had_brJet_MT2Zinv3b",dsDY_nunu_0b_3fake}})}, {2, 1}, "passBaselineNoTagZinv", 50, 0, 2000, true, true, label_mt2, "Norm Events"));

  vh.push_back(PHS("fakedCSVByOrder",  {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "", 100, 0, 1, false, true, "CSV value", "Norm Events"));
  vh.push_back(PHS("fakedCSVByOrder0b", {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "cntCSVSZinv=0",  100, 0, 1, false, true, "CSV value", "Norm Events"));
  vh.push_back(PHS("fakedCSVByOrder_baseline", {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "passBaselineZinv",  100, 0, 1, false, true, "CSV value", "Norm Events"));
  vh.push_back(PHS("fakedCSVByOrder_baselineNoTag", {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "passBaselineNoTagZinv", 100, 0, 1, false, true, "CSV value", "Norm Events"));
  vh.push_back(PHS("fakedCSVByOrder0b_baselineNoTag", {PDC("single", {{"fakedCSVValues[0]", dsDY_nunu}, {"fakedCSVValues[1]", dsDY_nunu}, {"fakedCSVValues[2]", dsDY_nunu}})}, {1, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0", 100, 0, 1, false, true, "CSV value", "Norm Events"));

  vh.push_back(PHS("fakedCSVAll",  {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "cntCSVSZinv=0",  100, 0, 1, false, true, "CSV value", "Norm Events"));
  vh.push_back(PHS("fakedCSVAll0b", {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "",  100, 0, 1, false, true, "CSV value", "Norm Events"));
  vh.push_back(PHS("fakedCSVAll_baseline", {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "passBaselineZinv", 100, 0, 1, false, true, "CSV value", "Norm Events"));
  vh.push_back(PHS("fakedCSVAll_baselineNoTag", {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1, false, true, "CSV value", "Norm Events"));
  vh.push_back(PHS("fakedCSVAll0b_baselineNoTag", {PDC("single", "fakedCSVValues", {dsDY_ll, dsDY_nunu, dsDY_ll_zAcc})}, {2, 1}, "passBaselineNoTagZinv;cntCSVSZinv=0", 100, 0, 1, false, true, "CSV value", "Norm Events"));

  vh.push_back(PHS("bybTag_nJets", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "", 30, 0, 30, true, true, label_nj, "Norm Events"));
  vh.push_back(PHS("bybTag_nJets_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, label_nj, "Norm Events"));

  vh.push_back(PHS("bybTag_nJets_1f", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_1fake, dsDY_nunu_1b})}, {2, 1}, "", 30, 0, 30, true, true, label_nj, "Norm Events"));
  vh.push_back(PHS("bybTag_nJets_1f_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_1fake, dsDY_nunu_1b})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, label_nj, "Norm Events"));

  vh.push_back(PHS("bybTag_nJets_2f", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_2fake, dsDY_nunu_2b})}, {2, 1}, "", 30, 0, 30, true, true, label_nj, "Norm Events"));
  vh.push_back(PHS("bybTag_nJets_2f_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_2fake, dsDY_nunu_2b})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, label_nj, "Norm Events"));

  vh.push_back(PHS("bybTag_nJets_3f", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_3fake, dsDY_nunu_3b})}, {2, 1}, "", 30, 0, 30, true, true, label_nj, "Norm Events"));
  vh.push_back(PHS("bybTag_nJets_3f_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_0b_3fake, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, label_nj, "Norm Events"));

  vh.push_back(PHS("bybTag_cleanMET", {PDC("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "", 100, 0, 1500, true, true, label_met, "Norm Events"));
  vh.push_back(PHS("bybTag_cleanMET_baselineNoTag", {PDC("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, label_met, "Norm Events"));

  vh.push_back(PHS("bybTag_cleanjPt", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("bybTag_cleanjPt_baselineNoTag", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));

  vh.push_back(PHS("bybTag_MT2", {PDC("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "", 100, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("bybTag_MT2_baselineNoTag", {PDC("single", "cleanMetPt", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 2000, true, true, label_mt2, "Norm Events"));

  vh.push_back(PHS("bybTag_Ntops", {PDC("single", "nTopCandSortedCntZinv", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));
  vh.push_back(PHS("bybTag_Ntops_baselineNoTag", {PDC("single", "nTopCandSortedCntZinv", {dsDY_nunu_0b, dsDY_nunu_1b, dsDY_nunu_2b, dsDY_nunu_3b})}, {2, 1}, "passBaselineNoTagZinv", 10, 0, 10, true, true, label_nt, "Norm Events"));


  vh.push_back(PHS("byCSV_nJets", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "", 30, 0, 30, true, true, label_nj, "Norm Events"));
  vh.push_back(PHS("byCSV_nJets_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 30, 0, 30, true, true, label_nj, "Norm Events"));

  vh.push_back(PHS("byCSV_cleanMET", {PDC("single", "cleanMetPt", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "", 100, 0, 1500, true, true, label_met, "Norm Events"));
  vh.push_back(PHS("byCSV_cleanMET_baselineNoTag", {PDC("single", "cleanMetPt", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, label_met, "Norm Events"));

  vh.push_back(PHS("byCSV_cleanjPt", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("byCSV_cleanjPt_baselineNoTag", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));

  vh.push_back(PHS("byCSV_MT2", {PDC("single", "cleanMetPt", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "", 100, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("byCSV_MT2_baselineNoTag", {PDC("single", "cleanMetPt", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 2000, true, true, label_mt2, "Norm Events"));

  vh.push_back(PHS("byCSV_Ntops", {PDC("single", "nTopCandSortedCntZinv", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "", 10, 0, 10, true, true, label_nt, "Norm Events"));
  vh.push_back(PHS("byCSV_Ntops_baselineNoTag", {PDC("single", "nTopCandSortedCntZinv", {dsDY_nunu_highCSV, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 10, 0, 10, true, true, label_nt, "Norm Events"));


  vh.push_back(PHS("byCSV1b_nJets_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv",  {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv",  30, 0, 30, true, true, label_nj, "Norm Events"));
  vh.push_back(PHS("byCSV1b_cleanMET_baselineNoTag", {PDC("single", "cleanMetPt",  {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, label_met, "Norm Events"));
  vh.push_back(PHS("byCSV1b_cleanjPt_baselineNoTag", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("byCSV1b_MT2_baselineNoTag",  {PDC("single", "cleanMetPt",  {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("byCSV1b_Ntops_baselineNoTag", {PDC("single", "nTopCandSortedCntZinv",  {dsDY_nunu_highCSV1b, dsDY_nunu_lowCSV})}, {2, 1}, "passBaselineNoTagZinv",  10, 0, 10, true, true, label_nt, "Norm Events"));

  vh.push_back(PHS("byCSVvh1b_nJets_baselineNoTag", {PDC("single", "cntNJetsPt30Eta24Zinv",  {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv",  30, 0, 30, true, true, label_nj, "Norm Events"));
  vh.push_back(PHS("byCSVvh1b_cleanMET_baselineNoTag", {PDC("single", "cleanMetPt",  {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, label_met, "Norm Events"));
  vh.push_back(PHS("byCSVvh1b_cleanjPt_baselineNoTag", {PDC("single", "cleanJetpt30ArrVec(pt)", {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 1500, true, true, "j p_{T} [GeV]", "Norm Events"));
  vh.push_back(PHS("byCSVvh1b_MT2_baselineNoTag",  {PDC("single", "cleanMetPt",  {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv", 100, 0, 2000, true, true, label_mt2, "Norm Events"));
  vh.push_back(PHS("byCSVvh1b_Ntops_baselineNoTag", {PDC("single", "nTopCandSortedCntZinv",  {dsDY_nunu_veryhighCSV1b, dsDY_nunu_verylowCSV})}, {2, 1}, "passBaselineNoTagZinv",  10, 0, 10, true, true, label_nt, "Norm Events"));


  vh.push_back(PHS("muPt_met_lt_150", {PDC("single", "cutMuVec(pt)", {dsDY_ll, dsDY_nunu})}, {1, 2}, "cleanMetPt<150", 100, 0, 1500,  true,  true,  "#mu p_{T} [GeV]",  "Norm Events"));
  vh.push_back(PHS("cleanht_met_lt_150", {scaled_cleanht}, {1, 1}, "cleanMetPt<150", 50, 0, 2000, true,  false,  label_ht, "Events"));
  vh.push_back(PHS("cleanmet_met_lt_150", {scaled_cleanMet}, {1, 1}, "cleanMetPt<150", 50, 0, 1500, true,  false,  label_met,  "Events"));
  vh.push_back(PHS("cleanmht_met_lt_150", {scaled_cleanmht}, {1, 1}, "cleanMetPt<150", 50, 0, 1500, true,  false,  label_mht, "Events"));


  // Plots the information needed to calculate acceptance corrections, plot itself is meaningless
  vh.push_back(PHS("accInfo",  {PDC("single", {{"genMuInAcc(size)", dsDY_ll_forEff}, {"genMu(size)", dsDY_ll_forEff}})}, {1, 2}, "",  20,  0, 20, true, false,  "gen N(#mu)", "Events"));

  vh.push_back(PHS("genMuPt",  {PDC("single", {{"genMuInAcc(pt)", dsDY_ll_forEff}, {"genMu(pt)", dsDY_ll_forEff}})},  {1, 1}, "",  150,  0,  1500, true,  false,  "gen #mu p_{T} [GeV]", "Events"));
  vh.push_back(PHS("accHt", {PDC("ratio",  {{"HTZinv", dsDY_ll_forEff_den}, {"HTZinv", dsDY_ll_forAcc}})}, {1, 1}, "",  300,  0,  3000, false, false,  label_ht,  "Acc"));
  vh.push_back(PHS("accMt", {PDC("ratio",  {{"nTopCandSortedCntZinv", dsDY_ll_forEff_den}, {"nTopCandSortedCntZinv", dsDY_ll_forAcc}})}, {1, 1}, "",  10, 0,  10, false, false,  label_nt, "Acc"));
  vh.push_back(PHS("accNb", {PDC("ratio",  {{"cntCSVSZinv", dsDY_ll_forEff_den}, {"cntCSVSZinv", dsDY_ll_forAcc}})}, {1, 1}, "",  10, 0,  10, false, false,  label_nb, "Acc"));
  vh.push_back(PHS("accNjet",  {PDC("ratio",  {{"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_den}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_forAcc}})}, {1, 1}, "",  20, 0,  20, false, false,  label_nj, "Acc"));

  // Plots for calculation of muon efficiency
  vh.push_back(PHS("effMuPt",  {PDC("ratio",  {{"genMatchMuInAcc(pt)", dsDY_ll_forEff}, {"genMuInAcc(pt)", dsDY_ll_forEff}})},  {1, 1}, "", 200, 0, 2000, false, false,  "#mu p_{T} [GeV]",  "Efficiency"));
  vh.push_back(PHS("effMuEta", {PDC("ratio",  {{"genMatchMuInAcc(eta)", dsDY_ll_forEff},  {"genMuInAcc(eta)", dsDY_ll_forEff}})}, {1, 1}, "", 100, -3.0,  3.0,  false, false,  "#mu #eta",  "Efficiency"));
  vh.push_back(PHS("effMucleanHt", {PDC("ratio",  {{"HTZinv", dsDY_ll_forEff_num}, {"HTZinv", dsDY_ll_forEff_den}})},  {1, 1}, "", 300, 0, 3000, false, false,  label_ht,  "Efficiency"));
  vh.push_back(PHS("effMuNb",  {PDC("ratio",  {{"cntCSVSZinv", dsDY_ll_forEff_num}, {"cntCSVSZinv", dsDY_ll_forEff_den}})}, {1, 1}, "", 10, 0,  10.0, false, false,  label_nb,  "Efficiency"));
  vh.push_back(PHS("effMuNt",  {PDC("ratio",  {{"nTopCandSortedCntZinv", dsDY_ll_forEff_num}, {"nTopCandSortedCntZinv", dsDY_ll_forEff_den}})},  {1, 1}, "", 10, 0,  10.0, false, false,  label_nt,  "Efficiency"));
  vh.push_back(PHS("effMuNjet", {PDC("ratio",  {{"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_num}, {"cntNJetsPt30Eta24Zinv", dsDY_ll_forEff_den}})},  {1, 1}, "", 20, 0,  20.0, false, false,  label_nj,  "Efficiency"));

  // Plots of weight and efficiencies used as sanity checks
  vh.push_back(PHS("zEffWgt", {PDC("single", "zEffWgt",  {dsDY_ll})},  {1, 1}, "", 100, 0.0,  10.0,  false, false,  "Z Eff Weight",  ""));
  vh.push_back(PHS("zEff", {PDC("single", "zEff", {dsDY_ll})},  {1, 1}, "", 100, 0.0,  1.0, false, false,  "Z Eff Weight",  ""));
  vh.push_back(PHS("zAccWgt", {PDC("single", "zAccWgt",  {dsDY_ll})},  {1, 1}, "", 100, 0.0,  10.0,  false, false,  "Z Acc Weight",  ""));
  vh.push_back(PHS("zAcc", {PDC("single", "zAcc", {dsDY_ll})},  {1, 1}, "", 100, 0.0,  1.0, false, false,  "Z Acc Weight",  ""));

  PDC dcDY_cutFlow("single", "cleanMetPt", {dsDY_ll_passMuonVeto, dsDY_ll_passEleVeto, dsDY_ll_passIsoTrkVeto, dsDY_ll_passnJets, dsDY_ll_passdPhis, dsDY_ll_passMET, dsDY_ll_passBJets, dsDY_ll_passTagger, dsDY_ll_passBaseline, dsDY_ll_passBaselineNoTag});
  vh.push_back(PHS("gcf_met", {dcDY_cutFlow},  {1, 1}, "", 100, 0.0,  3000.0, true, false, label_met,  "Events"));

  PDC dcDY_cutFlowZinv("single", "cleanMetPt", {dsDY_ll_passMuonVetoZinv, dsDY_ll_passEleVetoZinv, dsDY_ll_passIsoTrkVetoZinv, dsDY_ll_passnJetsZinv, dsDY_ll_passdPhisZinv, dsDY_ll_passMETZinv, dsDY_ll_passBJetsZinv, dsDY_ll_passTaggerZinv, dsDY_ll_passBaselineZinv, dsDY_ll_passBaselineNoTagZinv});
  vh.push_back(PHS("gcfZinv_met", {dcDY_cutFlowZinv},  {1, 1}, "", 100, 0.0,  3000.0, true, false, label_met,  "Events"));

  // Test plots, these are simply meaningless demonstrations of features
  Plotter::DataCollection dcDY_test( "single", {{"mht", dsDY_ll}, {"jetsLVec(pt)", dsDY_nunu}, {"genZmass", dsDY_test}, {"genZPt", dsDY_test}});
  Plotter::DataCollection dcDY_ratio(  "ratio",  {{"cleanJetVec(pt)", dsDY_ll}, {"jetsLVec(pt)", dsDY_ll}});
  Plotter::DataCollection dcDY_stack(  "stack",  "mht", {dsDY_ll, dsDY_nunu});
  Plotter::DataCollection dcDY_dataTest(  "data",  "met", {dsDY_nunu});

  vh.push_back(PHS("test",  {dcDY_stack, dcDY_test, dcDY_dataTest},  {3, 2}, "genZmass>80", 100, 0, 1000, true, true,  "???",  "Norm Events"));
  vh.push_back(PHS("test2", {dcDY_ratio},  {1, 1}, "", 100, 0, 500, true, false, "jpt???", "Ratio"));



  // -----------------
  // - Data/MC plots -
  // -----------------

  // Datasetsummaries we are using
  // no weight (genWeight deals with negative weights)
  Plotter::DatasetSummary dsData_SingleMuon(            "Data",       fileMap["Data_SingleMuon"],     "passMuTrigger",   "");
  Plotter::DatasetSummary dsData_SingleElectron(        "Data",       fileMap["Data_SingleElectron"], "passElecTrigger", "");//Sam
  Plotter::DatasetSummary dsData_DoubleEG(              "Data",       fileMap["Data_DoubleEG"],       "passElecTrigger", "");
  Plotter::DatasetSummary dsDY(                         "DY",         fileMap["DYJetsToLL"],          "",                "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dsDYInc(                      "DY HT<100",  fileMap["IncDY"],               "genHT<100",       "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dstt2l(                       "t#bar{t}",   fileMap["TTbarNoHad"],          "",                "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dstW(                         "single top", fileMap["tW"],                  "",                "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dsttZ(                        "t#bar{t}Z",  fileMap["TTZ"],                 "",                "bTagSF_EventWeightSimple_Central;genWeight");
  Plotter::DatasetSummary dsVV(                         "Diboson",    fileMap["Diboson"],             "",                "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dsRare(                       "Rare",       fileMap["Rare"],                "",                "bTagSF_EventWeightSimple_Central;genWeight");
  std::vector<Plotter::DatasetSummary> stack_MC = {dsDY, dsDYInc, dstt2l, dstW, dsttZ, dsVV, dsRare};
  // eff*Acc*B(Z)/B(DY)
  Plotter::DatasetSummary dsData_SingleMuon_effAcc(     "Data",       fileMap["Data_SingleMuon"],     "passMuTrigger",   "zEffWgt;zAccWgt",                                            znunu_mumu_ratio);
  Plotter::DatasetSummary dsData_SingleElectron_effAcc( "Data",       fileMap["Data_SingleElectron"], "passElecTrigger", "zEffWgtElec;zAccWgtElec",                                    znunu_mumu_ratio);//Sam
  Plotter::DatasetSummary dsData_DoubleEG_effAcc(       "Data",       fileMap["Data_DoubleEG"],       "passElecTrigger", "zEffWgt;zAccWgt",                                            znunu_mumu_ratio);
  Plotter::DatasetSummary dsDY_effAcc(                  "DY",         fileMap["DYJetsToLL"],          "",                "bTagSF_EventWeightSimple_Central;zEffWgt;zAccWgt",           znunu_mumu_ratio);
  Plotter::DatasetSummary dsDYInc_effAcc(               "DY HT<100",  fileMap["IncDY"],               "genHT<100",       "bTagSF_EventWeightSimple_Central;zEffWgt;zAccWgt",           znunu_mumu_ratio);
  Plotter::DatasetSummary dstt2l_effAcc(                "t#bar{t}",   fileMap["TTbarNoHad"],          "",                "bTagSF_EventWeightSimple_Central;zEffWgt;zAccWgt",           znunu_mumu_ratio);
  Plotter::DatasetSummary dstW_effAcc(                  "single top", fileMap["tW"],                  "",                "bTagSF_EventWeightSimple_Central;zEffWgt;zAccWgt",           znunu_mumu_ratio);
  Plotter::DatasetSummary dsttZ_effAcc(                 "t#bar{t}Z",  fileMap["TTZ"],                 "",                "bTagSF_EventWeightSimple_Central;zEffWgt;zAccWgt;genWeight", znunu_mumu_ratio);
  Plotter::DatasetSummary dsVV_effAcc(                  "Diboson",    fileMap["Diboson"],             "",                "bTagSF_EventWeightSimple_Central;zEffWgt;zAccWgt",           znunu_mumu_ratio);
  Plotter::DatasetSummary dsRare_effAcc(                "Rare",       fileMap["Rare"],                "",                "bTagSF_EventWeightSimple_Central;zEffWgt;zAccWgt;genWeight", znunu_mumu_ratio);
  std::vector<Plotter::DatasetSummary> stack_MC_effAcc = {dsDY_effAcc, dsDYInc_effAcc, dstt2l_effAcc, dstW_effAcc, dsttZ_effAcc, dsVV_effAcc, dsRare_effAcc};

  // Apply data/mc njet weight for DY and ttbar
  Plotter::DatasetSummary dswDY(  "DY",  fileMap["DYJetsToLL"],  "", "bTagSF_EventWeightSimple_Central;nJetWgtDYZ");
  Plotter::DatasetSummary dswDYInc(  "DY HT<100",  fileMap["IncDY"], "genHT<100", "bTagSF_EventWeightSimple_Central;nJetWgtDYZ");
  Plotter::DatasetSummary dswtt2l( "t#bar{t}", fileMap["TTbarNoHad"],  "", "bTagSF_EventWeightSimple_Central;nJetWgtTTbar");
  Plotter::DatasetSummary dswtW(  "single top", fileMap["tW"], "", "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dswttZ( "t#bar{t}Z",  fileMap["TTZ"],  "", "bTagSF_EventWeightSimple_Central;genWeight");
  Plotter::DatasetSummary dswVV(  "Diboson", fileMap["Diboson"],  "", "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dswRare( "Rare", fileMap["Rare"], "", "bTagSF_EventWeightSimple_Central;genWeight");
  std::vector<Plotter::DatasetSummary> stackw_MC = {dswDY, dswDYInc, dswtt2l, dswtW, dswttZ, dswVV, dswRare};

  // Extra MC stack that have njet weight applied for ttbar but not for DY (to make plots for AN)
  std::vector<Plotter::DatasetSummary> stackwtt_MC = {dsDY, dsDYInc, dswtt2l, dstW, dsttZ, dsVV, dsRare};

  // Apply data/mc njet weight for DY and ttbar & apply the normalization weight
  Plotter::DatasetSummary dswwDY(  "DY",  fileMap["DYJetsToLL"],  "", "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b");
  Plotter::DatasetSummary dswwDYInc(  "DY HT<100",  fileMap["IncDY"], "genHT<100", "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b");
  std::vector<Plotter::DatasetSummary> stackww_MC = {dswwDY, dswwDYInc, dswtt2l, dswtW, dswttZ, dswVV, dswRare};

  // Collections for all variables, no cuts applied yet
  //create lambda function to make PDCs
  auto dcData_SingleMuon = [&](std::string variable){ return Plotter::DataCollection("data",  variable, {dsData_SingleMuon}); };
  auto dcData_DoubleEG   = [&](std::string variable){ return Plotter::DataCollection("data",  variable, {dsData_DoubleEG}); };
  auto dcMC              = [&](std::string variable){ return Plotter::DataCollection("stack", variable, stack_MC); };
  auto dcwMC             = [&](std::string variable){ return Plotter::DataCollection("stack", variable, stackw_MC); };
  auto dcwttMC           = [&](std::string variable){ return Plotter::DataCollection("stack", variable, stackwtt_MC); };
  auto dcwwMC            = [&](std::string variable){ return Plotter::DataCollection("stack", variable, stackww_MC); };
  
  
  ////////////// met
  ////////////Plotter::DataCollection dcData_SingleMuon_met("data", "cleanMetPt", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_met("data", "cleanMetPt", {dsData_SingleElectron});/////copied //dsData_SingleElectron
  ////////////Plotter::DataCollection dcData_DoubleEG_met(  "data", "cleanMetPt", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_met(  "stack",  "cleanMetPt", stack_MC);
  ////////////Plotter::DataCollection dcwMC_met( "stack",  "cleanMetPt", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_met(  "stack",  "cleanMetPt", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_met( "stack",  "cleanMetPt", stackww_MC);
  ////////////// ntops
  ////////////Plotter::DataCollection dcData_SingleMuon_nt("data", "nTopCandSortedCntZinv", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_nt("data", "nTopCandSortedCntZinv", {dsData_SingleElectron});//copied
  ////////////Plotter::DataCollection dcData_DoubleEG_nt(  "data", "nTopCandSortedCntZinv", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_nt(  "stack",  "nTopCandSortedCntZinv", stack_MC);
  ////////////Plotter::DataCollection dcwMC_nt( "stack",  "nTopCandSortedCntZinv", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_nt(  "stack",  "nTopCandSortedCntZinv", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_nt( "stack",  "nTopCandSortedCntZinv", stackww_MC);
  ////////////// MT2
  ////////////Plotter::DataCollection dcData_SingleMuon_mt2("data", "best_had_brJet_MT2Zinv", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_mt2("data", "best_had_brJet_MT2Zinv", {dsData_SingleElectron});//copied
  ////////////Plotter::DataCollection dcData_DoubleEG_mt2(  "data", "best_had_brJet_MT2Zinv", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_mt2(  "stack",  "best_had_brJet_MT2Zinv", stack_MC);
  ////////////Plotter::DataCollection dcwMC_mt2( "stack",  "best_had_brJet_MT2Zinv", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_mt2(  "stack",  "best_had_brJet_MT2Zinv", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_mt2( "stack",  "best_had_brJet_MT2Zinv", stackww_MC);
  ////////////// nb
  ////////////Plotter::DataCollection dcData_SingleMuon_nb("data", "cntCSVSZinv", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_nb("data", "cntCSVSZinv", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_nb(  "data", "cntCSVSZinv", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_nb(  "stack",  "cntCSVSZinv", stack_MC);
  ////////////Plotter::DataCollection dcwMC_nb( "stack",  "cntCSVSZinv", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_nb(  "stack",  "cntCSVSZinv", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_nb( "stack",  "cntCSVSZinv", stackww_MC);
  ////////////// nj
  ////////////Plotter::DataCollection dcData_SingleMuon_nj("data", "cntNJetsPt30Eta24Zinv", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_nj("data", "cntNJetsPt30Eta24Zinv", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_nj(  "data", "cntNJetsPt30Eta24Zinv", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_nj(  "stack",  "cntNJetsPt30Eta24Zinv", stack_MC);
  ////////////Plotter::DataCollection dcwMC_nj( "stack",  "cntNJetsPt30Eta24Zinv", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_nj(  "stack",  "cntNJetsPt30Eta24Zinv", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_nj( "stack",  "cntNJetsPt30Eta24Zinv", stackww_MC);
  ////////////// ht
  ////////////Plotter::DataCollection dcData_SingleMuon_ht("data", "HTZinv", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_ht("data", "HTZinv", {dsData_SingleElectron});
  ////////////Plotter::DataCollection dcData_DoubleEG_ht(  "data", "HTZinv", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_ht(  "stack",  "HTZinv", stack_MC);
  ////////////Plotter::DataCollection dcwMC_ht( "stack",  "HTZinv", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_ht(  "stack",  "HTZinv", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_ht( "stack",  "HTZinv", stackww_MC);
  ////////////// mht
  ////////////Plotter::DataCollection dcData_SingleMuon_mht("data", "cleanMHt", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_mht("data", "cleanMHt", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_mht(  "data", "cleanMHt", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_mht(  "stack",  "cleanMHt", stack_MC);
  ////////////Plotter::DataCollection dcwMC_mht( "stack",  "cleanMHt", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_mht(  "stack",  "cleanMHt", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_mht( "stack",  "cleanMHt", stackww_MC);
  ////////////// jpt
  ////////////Plotter::DataCollection dcData_SingleMuon_jpt("data", "jetsLVecLepCleaned(pt)", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_jpt("data", "jetsLVecLepCleaned(pt)", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_jpt(  "data", "jetsLVecLepCleaned(pt)", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_jpt(  "stack",  "jetsLVecLepCleaned(pt)", stack_MC);
  ////////////Plotter::DataCollection dcwMC_jpt( "stack",  "jetsLVecLepCleaned(pt)", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_jpt(  "stack",  "jetsLVecLepCleaned(pt)", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_jpt( "stack",  "jetsLVecLepCleaned(pt)", stackww_MC);
  ////////////// j1pt
  ////////////Plotter::DataCollection dcData_SingleMuon_j1pt("data", "jetsLVecLepCleaned[0](pt)", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_j1pt("data", "jetsLVecLepCleaned[0](pt)", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_j1pt(  "data", "jetsLVecLepCleaned[0](pt)", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_j1pt(  "stack",  "jetsLVecLepCleaned[0](pt)", stack_MC);
  ////////////Plotter::DataCollection dcwMC_j1pt( "stack",  "jetsLVecLepCleaned[0](pt)", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_j1pt(  "stack",  "jetsLVecLepCleaned[0](pt)", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_j1pt( "stack",  "jetsLVecLepCleaned[0](pt)", stackww_MC);
  ////////////// j2pt
  ////////////Plotter::DataCollection dcData_SingleMuon_j2pt("data", "jetsLVecLepCleaned[1](pt)", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_j2pt("data", "jetsLVecLepCleaned[1](pt)", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_j2pt(  "data", "jetsLVecLepCleaned[1](pt)", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_j2pt(  "stack",  "jetsLVecLepCleaned[1](pt)", stack_MC);
  ////////////Plotter::DataCollection dcwMC_j2pt( "stack",  "jetsLVecLepCleaned[1](pt)", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_j2pt(  "stack",  "jetsLVecLepCleaned[1](pt)", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_j2pt( "stack",  "jetsLVecLepCleaned[1](pt)", stackww_MC);
  ////////////// j3pt
  ////////////Plotter::DataCollection dcData_SingleMuon_j3pt("data", "jetsLVecLepCleaned[2](pt)", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_j3pt("data", "jetsLVecLepCleaned[2](pt)", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_j3pt(  "data", "jetsLVecLepCleaned[2](pt)", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_j3pt(  "stack",  "jetsLVecLepCleaned[2](pt)", stack_MC);
  ////////////Plotter::DataCollection dcwMC_j3pt( "stack",  "jetsLVecLepCleaned[2](pt)", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_j3pt(  "stack",  "jetsLVecLepCleaned[2](pt)", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_j3pt( "stack",  "jetsLVecLepCleaned[2](pt)", stackww_MC);
  ////////////// mupt
  ////////////Plotter::DataCollection dcData_SingleMuon_mupt("data", "cutMuVec(pt)", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_mupt("data", "cutMuVec(pt)", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_mupt(  "data", "cutMuVec(pt)", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_mupt(  "stack",  "cutMuVec(pt)", stack_MC);
  ////////////Plotter::DataCollection dcwMC_mupt( "stack",  "cutMuVec(pt)", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_mupt(  "stack",  "cutMuVec(pt)", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_mupt( "stack",  "cutMuVec(pt)", stackww_MC);
  ////////////// mu1pt
  ////////////Plotter::DataCollection dcData_SingleMuon_mu1pt("data", "cutMuPt1", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_mu1pt("data", "cutMuPt1", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_mu1pt(  "data", "cutMuPt1", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_mu1pt(  "stack",  "cutMuPt1", stack_MC);
  ////////////Plotter::DataCollection dcwMC_mu1pt( "stack",  "cutMuPt1", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_mu1pt(  "stack",  "cutMuPt1", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_mu1pt( "stack",  "cutMuPt1", stackww_MC);
  ////////////// mu2pt
  ////////////Plotter::DataCollection dcData_SingleMuon_mu2pt("data", "cutMuPt2", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_mu2pt("data", "cutMuPt2", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_mu2pt(  "data", "cutMuPt2", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_mu2pt(  "stack",  "cutMuPt2", stack_MC);
  ////////////Plotter::DataCollection dcwMC_mu2pt( "stack",  "cutMuPt2", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_mu2pt(  "stack",  "cutMuPt2", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_mu2pt( "stack",  "cutMuPt2", stackww_MC);
  ////////////// elpt
  ////////////Plotter::DataCollection dcData_SingleMuon_elpt("data", "cutElecVec(pt)", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_elpt("data", "cutElecVec(pt)", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_elpt(  "data", "cutElecVec(pt)", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_elpt(  "stack",  "cutElecVec(pt)", stack_MC);
  ////////////Plotter::DataCollection dcwMC_elpt( "stack",  "cutElecVec(pt)", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_elpt(  "stack",  "cutElecVec(pt)", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_elpt( "stack",  "cutElecVec(pt)", stackww_MC);
  ////////////// el1pt
  ////////////Plotter::DataCollection dcData_SingleMuon_el1pt("data", "cutElecPt1", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_el1pt("data", "cutElecPt1", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_el1pt(  "data", "cutElecPt1", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_el1pt(  "stack",  "cutElecPt1", stack_MC);
  ////////////Plotter::DataCollection dcwMC_el1pt( "stack",  "cutElecPt1", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_el1pt(  "stack",  "cutElecPt1", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_el1pt( "stack",  "cutElecPt1", stackww_MC);
  ////////////// el2pt
  ////////////Plotter::DataCollection dcData_SingleMuon_el2pt("data", "cutElecPt2", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_el2pt("data", "cutElecPt2", {dsData_SingleElectron});
  ////////////Plotter::DataCollection dcData_DoubleEG_el2pt(  "data", "cutElecPt2", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_el2pt(  "stack",  "cutElecPt2", stack_MC);
  ////////////Plotter::DataCollection dcwMC_el2pt( "stack",  "cutElecPt2", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_el2pt(  "stack",  "cutElecPt2", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_el2pt( "stack",  "cutElecPt2", stackww_MC);
  ////////////// mll
  ////////////Plotter::DataCollection dcData_SingleMuon_mll("data", "bestRecoZM", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_mll("data", "bestRecoZM", {dsData_SingleElectron});
  ////////////Plotter::DataCollection dcData_DoubleEG_mll(  "data", "bestRecoZM", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_mll(  "stack",  "bestRecoZM", stack_MC);
  ////////////Plotter::DataCollection dcwMC_mll( "stack",  "bestRecoZM", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_mll(  "stack",  "bestRecoZM", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_mll( "stack",  "bestRecoZM", stackww_MC);
  ////////////// nsearchbins
  ////////////Plotter::DataCollection dcData_SingleMuon_nSearchBin("data", "nSearchBin", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_nSearchBin("data", "nSearchBin", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_nSearchBin(  "data", "nSearchBin", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_nSearchBin(  "stack",  "nSearchBin", stack_MC);
  ////////////Plotter::DataCollection dcwMC_nSearchBin( "stack",  "nSearchBin", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_nSearchBin(  "stack",  "nSearchBin", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_nSearchBin( "stack",  "nSearchBin", stackww_MC);
  ////////////// nb0Bins
  ////////////Plotter::DataCollection dcData_SingleMuon_nb0Bins("data", "nb0Bins", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_nb0Bins("data", "nb0Bins", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_nb0Bins(  "data", "nb0Bins", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_nb0Bins(  "stack",  "nb0Bins", stack_MC);
  ////////////Plotter::DataCollection dcwMC_nb0Bins( "stack",  "nb0Bins", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_nb0Bins(  "stack",  "nb0Bins", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_nb0Bins( "stack",  "nb0Bins", stackww_MC);
  ////////////// nb0BinsNW
  ////////////Plotter::DataCollection dcData_SingleMuon_nb0BinsNW("data", "nb0BinsNW", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_nb0BinsNW("data", "nb0BinsNW", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_nb0BinsNW(  "data", "nb0BinsNW", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_nb0BinsNW(  "stack",  "nb0BinsNW", stack_MC);
  ////////////Plotter::DataCollection dcwMC_nb0BinsNW( "stack",  "nb0BinsNW", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_nb0BinsNW(  "stack",  "nb0BinsNW", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_nb0BinsNW( "stack",  "nb0BinsNW", stackww_MC);
  ////////////// nb0NJwBins
  ////////////Plotter::DataCollection dcData_SingleMuon_nb0NJwBins("data", "nb0NJwBins", {dsData_SingleMuon});
  ////////////Plotter::DataCollection dcData_SingleElectron_nb0NJwBins("data", "nb0NJwBins", {dsData_SingleElectron});//
  ////////////Plotter::DataCollection dcData_DoubleEG_nb0NJwBins(  "data", "nb0NJwBins", {dsData_DoubleEG});
  ////////////Plotter::DataCollection dcMC_nb0NJwBins(  "stack",  "nb0NJwBins", stack_MC);
  ////////////Plotter::DataCollection dcwMC_nb0NJwBins( "stack",  "nb0NJwBins", stackw_MC);
  ////////////Plotter::DataCollection dcwttMC_nb0NJwBins( "stack",  "nb0NJwBins", stackwtt_MC);
  ////////////Plotter::DataCollection dcwwMC_nb0NJwBins( "stack",  "nb0NJwBins", stackww_MC);
  ////////////// nb0NJwBins DY->Z
  ////////////Plotter::DataCollection dcData_SingleMuon_nb0NJwBins_effAcc("data", "nb0NJwBins", {dsData_SingleMuon_effAcc});
  ////////////Plotter::DataCollection dcData_SingleElectron_nb0NJwBins_effAcc("data", "nb0NJwBins", {dsData_SingleElectron_effAcc});//dsData_DoubleEG
  ////////////Plotter::DataCollection dcData_DoubleEG_nb0NJwBins_effAcc(  "data", "nb0NJwBins", {dsData_DoubleEG_effAcc});
  ////////////Plotter::DataCollection dcMC_nb0NJwBins_effAcc(  "stack",  "nb0NJwBins", stack_MC_effAcc);


  // pair of cutlevels
  std::vector<std::pair<std::string,std::string>> cutlevels_muon = {
    {"nosel",  "passNoiseEventFilterZinv"},
    {"2mu",  "passNoiseEventFilterZinv;passDiMuSel"},
    {"muZinv", "passNoiseEventFilterZinv;passMuZinvSel"},
    //{"muZinv_ht200", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200"},
    //{"muZinv_ht200_dphi",  "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passdPhisZinv"},
    //{"muZinv_ht50_met50_dphi", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
    {"muZinv_loose0",  "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    {"muZinv_loose0_mt2",  "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>0"},
    //{"muZinv_loose0_mt21b", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv1b>0"},
    //{"muZinv_loose0_mt22b", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv2b>0"},
    //{"muZinv_loose0_mt23b", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv3b>0"},
    {"muZinv_loose0_ntop", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv>0"},
    {"muZinv_loose0_nb2",  "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv>=2"},
    //{"muZinv_loose0_ht300", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>300;passnJetsZinv;passdPhisZinv"},
    //{"muZinv_loose0_ht400", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>400;passnJetsZinv;passdPhisZinv"},
    //{"muZinv_loose0_ht500", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>500;passnJetsZinv;passdPhisZinv"},
    {"muZinv_blnotagmt2",  "passMuZinvSel;passBaselineNoTagMT2Zinv"},
    {"muZinv_blnotag", "passMuZinvSel;passBaselineNoTagZinv"},
    {"muZinv_bl",  "passMuZinvSel;passBaselineZinv"},
    {"muZinv_0b",  "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0"},
    //{"muZinv_0b_ht200", "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>200"},
    //{"muZinv_0b_ht200_dphi",  "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>200;passdPhisZinv"},
    //{"muZinv_0b_ht50_met50_dphi", "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
    {"muZinv_0b_loose0",  "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    //{"muZinv_0b_loose0_ht300", "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>300;passnJetsZinv;passdPhisZinv"},
    //{"muZinv_0b_loose0_ht400", "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>400;passnJetsZinv;passdPhisZinv"},
    //{"muZinv_0b_loose0_ht500", "passNoiseEventFilterZinv;passMuZinvSel;cntCSVSZinv=0;HTZinv>500;passnJetsZinv;passdPhisZinv"},
    {"muZinv_0b_blnotagmt2",  "passMuZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv"},
    {"muZinv_0b_blnotag",  "passMuZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv"},
    {"muZinv_g1b", "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv"},
    {"muZinv_g1b_ht200",  "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>200"},
    //{"muZinv_g1b_ht200_dphi", "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>200;passdPhisZinv"},
    //{"muZinv_g1b_ht50_met50_dphi","passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
    {"muZinv_g1b_loose0",  "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    //{"muZinv_g1b_loose0_ht300", "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>300;passnJetsZinv;passdPhisZinv"},
    //{"muZinv_g1b_loose0_ht400", "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>400;passnJetsZinv;passdPhisZinv"},
    //{"muZinv_g1b_loose0_ht500", "passNoiseEventFilterZinv;passMuZinvSel;passBJetsZinv;HTZinv>500;passnJetsZinv;passdPhisZinv"},
    {"elmu", "passNoiseEventFilterZinv;passElMuSel"},
    {"elmuZinv", "passNoiseEventFilterZinv;passElMuZinvSel"},
    //{"elmuZinv_ht200", "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200"},
    //{"elmuZinv_ht200_dphi", "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200;passdPhisZinv"},
    {"elmuZinv_loose0", "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    {"elmuZinv_blnotagmt2", "passElMuZinvSel;passBaselineNoTagMT2Zinv"},
    {"elmuZinv_bl", "passElMuZinvSel;passBaselineZinv"},
    {"elmuZinv_0b", "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0"},
    //{"elmuZinv_0b_ht200",  "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200"},
    //{"elmuZinv_0b_ht200_dphi", "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passdPhisZinv"},
    {"elmuZinv_0b_loose0", "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    {"elmuZinv_0b_blnotagmt2", "passElMuZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv"},
    {"elmuZinv_0b_blnotag", "passElMuZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv"},
    {"elmuZinv_g1b", "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv"},
    //{"elmuZinv_g1b_ht200", "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200"},
    //{"elmuZinv_g1b_ht200_dphi", "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200;passdPhisZinv"},
    {"elmuZinv_g1b_loose0", "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv"}
  };

  //std::vector<std::pair<std::string,std::string>> cutlevels_muon_syst = {
  //  {"muZinv_loose0",  "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
  //
  //  {"muZinv_0b_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv=0"},
  //  {"muZinv_1b_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv=1"},
  //  {"muZinv_2b_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv=2"},
  //  {"muZinv_3b_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv>2"},
  //
  //  {"muZinv_0t_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv=0"},
  //  {"muZinv_1t_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv=1"},
  //  {"muZinv_2t_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv>1"},
  //
  //  {"muZinv_met_0_100_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>0;cleanMetPt<100"},
  //  {"muZinv_met_100_200_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100;cleanMetPt<200"},
  //  {"muZinv_met_200_300_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>200;cleanMetPt<300"},
  //  {"muZinv_met_300_400_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>300;cleanMetPt<400"},
  //  {"muZinv_met_gt400_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>400"},
  //
  //  {"muZinv_mt2_0_100_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>0;best_had_brJet_MT2Zinv<100"},
  //  {"muZinv_mt2_100_200_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>100;best_had_brJet_MT2Zinv<200"},
  //  {"muZinv_mt2_200_300_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>200;best_had_brJet_MT2Zinv<300"},
  //  {"muZinv_mt2_300_400_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>300;best_had_brJet_MT2Zinv<400"},
  //  {"muZinv_mt2_gt400_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>400"},
  //
  //  {"muZinv_ht_200_300_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;HTZinv>200;HTZinv<300"},
  //  {"muZinv_ht_300_400_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;HTZinv>300;HTZinv<400"},
  //  {"muZinv_ht_400_500_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;HTZinv>400;HTZinv<500"},
  //  {"muZinv_ht_gt500_loose0", "passNoiseEventFilterZinv;passMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;HTZinv>500"},
  //};

  // std::vector<std::pair<std::string,std::string>> cutlevels_electron = {
  // 	{"nosel",  "passNoiseEventFilterZinv"},
  // 	{"2el",  "passNoiseEventFilterZinv;passDiElecSel"},
  // 	{"elZinv", "passNoiseEventFilterZinv;passElecZinvSel"},
  // 	{"elZinv_ht200",  "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200"},
  // 	{"elZinv_blnotagmt2", "passElecZinvSel;passBaselineNoTagMT2Zinv"},
  // 	{"elZinv_bl",  "passElecZinvSel;passBaselineZinv"},
  // 	{"elZinv_0b",  "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0"},
  // 	{"elZinv_0b_ht200", "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>200"},
  // 	{"s_elZinv_0b_loose0", "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv"},
  // 	{"elZinv_0b_blnotagmt2",  "passElecZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv"},
  // 	{"elZinv_0b_blnotag", "passElecZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv"},
  // };

  // pair of cutlevels
  std::vector<std::pair<std::string,std::string>> cutlevels_electron = {
    {"nosel",                      "passNoiseEventFilterZinv"},
    {"2el",                        "passNoiseEventFilterZinv;passDiElecSel"},
    {"elZinv",                     "passNoiseEventFilterZinv;passElecZinvSel"},
    //{"elZinv_ht200",               "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200"},
    //{"elZinv_ht200_dphi",          "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passdPhisZinv"},
    //{"elZinv_ht50_met50_dphi",     "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
    {"elZinv_loose0",              "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    {"elZinv_loose0_mt2",          "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>0"},
    //{"elZinv_loose0_mt21b",        "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv1b>0"},
    //{"elZinv_loose0_mt22b",        "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv2b>0"},
    //{"elZinv_loose0_mt23b",        "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv3b>0"},
    //{"elZinv_loose0_ntop",         "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv>0"},
    {"elZinv_loose0_nb2",          "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv>=2"},
    //{"elZinv_loose0_ht300",        "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>300;passnJetsZinv;passdPhisZinv"},
    //{"elZinv_loose0_ht400",        "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>400;passnJetsZinv;passdPhisZinv"},
    //{"elZinv_loose0_ht500",        "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>500;passnJetsZinv;passdPhisZinv"},
    {"elZinv_blnotagmt2",          "passElecZinvSel;passBaselineNoTagMT2Zinv"},
    {"elZinv_blnotag",             "passElecZinvSel;passBaselineNoTagZinv"},
    {"elZinv_bl",                  "passElecZinvSel;passBaselineZinv"},
    {"elZinv_0b",                  "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0"},
    //{"elZinv_0b_ht200",            "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>200"},
    //{"elZinv_0b_ht200_dphi",       "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>200;passdPhisZinv"},
    //{"elZinv_0b_ht50_met50_dphi",  "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
    {"elZinv_0b_loose0",           "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    //{"elZinv_0b_loose0_ht300",     "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>300;passnJetsZinv;passdPhisZinv"},
    //{"elZinv_0b_loose0_ht400",     "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>400;passnJetsZinv;passdPhisZinv"},
    //{"elZinv_0b_loose0_ht500",     "passNoiseEventFilterZinv;passElecZinvSel;cntCSVSZinv=0;HTZinv>500;passnJetsZinv;passdPhisZinv"},
    {"elZinv_0b_blnotagmt2",       "passElecZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv"},
    {"elZinv_0b_blnotag",          "passElecZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv"},
    //{"elZinv_g1b",                 "passNoiseEventFilterZinv;passElecZinvSel;passBJetsZinv"},
    //{"elZinv_g1b_ht200",           "passNoiseEventFilterZinv;passElecZinvSel;passBJetsZinv;HTZinv>200"},
    //{"elZinv_g1b_ht200_dphi",      "passNoiseEventFilterZinv;passElecZinvSel;passBJetsZinv;HTZinv>200;passdPhisZinv"},
    //{"elZinv_g1b_ht50_met50_dphi", "passNoiseEventFilterZinv;passElecZinvSel;passBJetsZinv;HTZinv>50;cleanMetPt>50;passdPhisZinv"},
    {"elZinv_g1b_loose0",          "passNoiseEventFilterZinv;passElecZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    //{"elZinv_g1b_loose0_ht300",    "passNoiseEventFilterZinv;passElecZinvSel;passBJetsZinv;HTZinv>300;passnJetsZinv;passdPhisZinv"},
    //{"elZinv_g1b_loose0_ht400",    "passNoiseEventFilterZinv;passElecZinvSel;passBJetsZinv;HTZinv>400;passnJetsZinv;passdPhisZinv"},
    //{"elZinv_g1b_loose0_ht500",    "passNoiseEventFilterZinv;passElecZinvSel;passBJetsZinv;HTZinv>500;passnJetsZinv;passdPhisZinv"},
    //{"elmu",                       "passNoiseEventFilterZinv;passElMuSel"},
    //{"elmuZinv",                   "passNoiseEventFilterZinv;passElMuZinvSel"},
    //{"elmuZinv_ht200",             "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200"},
    //{"elmuZinv_ht200_dphi",        "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200;passdPhisZinv"},
    //{"elmuZinv_loose0",            "passNoiseEventFilterZinv;passElMuZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    //{"elmuZinv_blnotagmt2",        "passElMuZinvSel;passBaselineNoTagMT2Zinv"},
    //{"elmuZinv_bl",                "passElMuZinvSel;passBaselineZinv"},
    //{"elmuZinv_0b",                "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0"},
    //{"elmuZinv_0b_ht200",          "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200"},
    ////{"elmuZinv_0b_ht200_dphi",     "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passdPhisZinv"},
    //{"elmuZinv_0b_loose0",         "passNoiseEventFilterZinv;passElMuZinvSel;cntCSVSZinv=0;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    //{"elmuZinv_0b_blnotagmt2",     "passElMuZinvSel;cntCSVSZinv=0;passBaselineNoTagMT2Zinv"},
    //{"elmuZinv_0b_blnotag",        "passElMuZinvSel;cntCSVSZinv=0;passBaselineNoTagZinv"},
    //{"elmuZinv_g1b",               "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv"},
    //{"elmuZinv_g1b_ht200",         "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200"},
    //{"elmuZinv_g1b_ht200_dphi",    "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200;passdPhisZinv"},
    //{"elmuZinv_g1b_loose0",        "passNoiseEventFilterZinv;passElMuZinvSel;passBJetsZinv;HTZinv>200;passnJetsZinv;passdPhisZinv"}
  };

  //std::vector<std::pair<std::string,std::string>> cutlevels_electron_syst = {
  //  {"elZinv_loose0",  "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv"},
  //
  //  {"elZinv_0b_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv=0"},
  //  {"elZinv_1b_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv=1"},
  //  {"elZinv_2b_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv=2"},
  //  {"elZinv_3b_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cntCSVSZinv>2"},
  //
  //  {"elZinv_0t_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv=0"},
  //  {"elZinv_1t_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv=1"},
  //  {"elZinv_2t_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;nTopCandSortedCntZinv>1"},
  //
  //  {"elZinv_met_0_100_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>0;cleanMetPt<100"},
  //  {"elZinv_met_100_200_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100;cleanMetPt<200"},
  //  {"elZinv_met_200_300_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>200;cleanMetPt<300"},
  //  {"elZinv_met_300_400_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>300;cleanMetPt<400"},
  //  {"elZinv_met_gt400_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>400"},
  //
  //  {"elZinv_mt2_0_100_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>0;best_had_brJet_MT2Zinv<100"},
  //  {"elZinv_mt2_100_200_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>100;best_had_brJet_MT2Zinv<200"},
  //  {"elZinv_mt2_200_300_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>200;best_had_brJet_MT2Zinv<300"},
  //  {"elZinv_mt2_300_400_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>300;best_had_brJet_MT2Zinv<400"},
  //  {"elZinv_mt2_gt400_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;best_had_brJet_MT2Zinv>400"},
  //
  //  {"elZinv_ht_200_300_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;HTZinv>200;HTZinv<300"},
  //  {"elZinv_ht_300_400_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;HTZinv>300;HTZinv<400"},
  //  {"elZinv_ht_400_500_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;HTZinv>400;HTZinv<500"},
  //  {"elZinv_ht_gt500_loose0", "passNoiseEventFilterZinv;passElecZinvSel;HTZinv>200;passnJetsZinv;passdPhisZinv;HTZinv>500"},
  //};

  //push the histograms in a loop, save some copy-paste time
  for(std::pair<std::string,std::string>& cut : cutlevels_muon)
    {
      // unweighted
      vh.push_back(PHS("DataMC_SingleMuon_met" +cut.first,  {dcData_SingleMuon("cleanMetPt"), dcMC("cleanMetPt")},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
      vh.push_back(PHS("DataMC_SingleMuon_ht_"  +cut.first,  {dcData_SingleMuon("HTZinv"), dcMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
      vh.push_back(PHS("DataMC_SingleMuon_mht_" +cut.first,  {dcData_SingleMuon("cleanMHt"), dcMC("cleanMHt")},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
      vh.push_back(PHS("DataMC_SingleMuon_nt_"  +cut.first,  {dcData_SingleMuon("nTopCandSortedCntZinv"), dcMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
      //vh.push_back(PHS("DataMC_SingleMuon_nt1b_" +cut.first,  {dcData_SingleMuon_nt1b,  dcMC_nt1b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
      //vh.push_back(PHS("DataMC_SingleMuon_nt2b_" +cut.first,  {dcData_SingleMuon_nt2b,  dcMC_nt2b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
      //vh.push_back(PHS("DataMC_SingleMuon_nt3b_" +cut.first,  {dcData_SingleMuon_nt3b,  dcMC_nt3b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
      vh.push_back(PHS("DataMC_SingleMuon_mt2_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv"), dcMC("best_had_brJet_MT2Zinv")},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
      //vh.push_back(PHS("DataMC_SingleMuon_mt21b_"  +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")1b, dcMC("best_had_brJet_MT2Zinv")1b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_mt22b_"  +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")2b, dcMC("best_had_brJet_MT2Zinv")2b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_mt23b_"  +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")3b, dcMC("best_had_brJet_MT2Zinv")3b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
      vh.push_back(PHS("DataMC_SingleMuon_nb_"  +cut.first,  {dcData_SingleMuon("cntCSVSZinv"), dcMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
      vh.push_back(PHS("DataMC_SingleMuon_nj_"  +cut.first,  {dcData_SingleMuon("cntNJetsPt30Eta24Zinv"), dcMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_jpt_" +cut.first,  {dcData_SingleMuon_jpt, dcMC_jpt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_j1pt_" +cut.first,  {dcData_SingleMuon_j1pt,  dcMC_j1pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
      //vh.push_back(PHS("DataMC_SingleMuon_j2pt_" +cut.first,  {dcData_SingleMuon_j2pt,  dcMC_j2pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
      //vh.push_back(PHS("DataMC_SingleMuon_j3pt_" +cut.first,  {dcData_SingleMuon_j3pt,  dcMC_j3pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
      //vh.push_back(PHS("DataMC_SingleMuon_mupt_" +cut.first,  {dcData_SingleMuon_mupt,  dcMC_mupt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mupt,  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_mu1pt_"  +cut.first,  {dcData_SingleMuon_mu1pt, dcMC_mu1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu1pt,  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_mu2pt_"  +cut.first,  {dcData_SingleMuon_mu2pt, dcMC_mu2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu2pt,  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_elpt_" +cut.first,  {dcData_SingleMuon_elpt,  dcMC_elpt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_el1pt_"  +cut.first,  {dcData_SingleMuon_el1pt, dcMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_el2pt_"  +cut.first,  {dcData_SingleMuon_el2pt, dcMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
      //vh.push_back(PHS("DataMC_SingleMuon_mll_" +cut.first,  {dcData_SingleMuon_mll, dcMC_mll},  {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
      vh.push_back(PHS("DataMC_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon("nSearchBin"), dcMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMC_SingleMuon_nb0Bins_" +cut.first,  {dcData_SingleMuon_nb0Bins, dcMC_nb0Bins}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMC_SingleMuon_nb0BinsNW_"  +cut.first,  {dcData_SingleMuon_nb0BinsNW,  dcMC_nb0BinsNW},  {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMC_SingleMuon_nb0NJwBins_" +cut.first,  {dcData_SingleMuon_nb0NJwBins, dcMC_nb0NJwBins}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMC_SingleMuon_nb0NJwBins_effAcc_" +cut.first,  {dcData_SingleMuon_nb0NJwBins_effAcc, dcMC_nb0NJwBins_effAcc}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));

      // DataMC weights applied
      vector<double> metBins = {0, 50, 100, 150, 200, 275, 300, 350, 400, 450, 2000};
      vector<double> mt2Bins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 2000};
      vh.push_back(PHS("DataMCw_SingleMuon_met_" +cut.first,  {dcData_SingleMuon("cleanMetPt"), dcwMC("cleanMetPt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
      vh.push_back(PHS("DataMCw_SingleMuon_rebin_met_" +cut.first,  {dcData_SingleMuon("cleanMetPt"), dcwMC("cleanMetPt")}, {1, 2}, cut.second, metBins, true, false,  label_met, ""));
      vh.push_back(PHS("DataMCw_SingleMuon_ht_" +cut.first,  {dcData_SingleMuon("HTZinv"), dcwMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
      vh.push_back(PHS("DataMCw_SingleMuon_mht_" +cut.first,  {dcData_SingleMuon("cleanMHt"), dcwMC("cleanMHt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
      vh.push_back(PHS("DataMCw_SingleMuon_nt_" +cut.first,  {dcData_SingleMuon("nTopCandSortedCntZinv"), dcwMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_nt1b_"  +cut.first,  {dcData_SingleMuon_nt1b,  dcwMC_nt1b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_nt2b_"  +cut.first,  {dcData_SingleMuon_nt2b,  dcwMC_nt2b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_nt3b_"  +cut.first,  {dcData_SingleMuon_nt3b,  dcwMC_nt3b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
      vh.push_back(PHS("DataMCw_SingleMuon_mt2_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv"), dcwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_mt21b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")1b, dcwMC("best_had_brJet_MT2Zinv")1b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_mt22b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")2b, dcwMC("best_had_brJet_MT2Zinv")2b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_mt23b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")3b, dcwMC("best_had_brJet_MT2Zinv")3b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
      vh.push_back(PHS("DataMCw_SingleMuon_rebin_mt2_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv"), dcwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, mt2Bins, true, false,  label_mt2, ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_rebin_mt21b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")1b, dcwMC("best_had_brJet_MT2Zinv")1b}, {1, 2}, cut.second, mt2Bins, true, false,  "mt2 (1b fake)",  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_rebin_mt22b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")2b, dcwMC("best_had_brJet_MT2Zinv")2b}, {1, 2}, cut.second, mt2Bins, true, false,  "mt2 (2b fake)",  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_rebin_mt23b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")3b, dcwMC("best_had_brJet_MT2Zinv")3b}, {1, 2}, cut.second, mt2Bins, true, false,  "mt2 (3b fake)",  ""));
      vh.push_back(PHS("DataMCw_SingleMuon_nb_" +cut.first,  {dcData_SingleMuon("cntCSVSZinv"), dcwMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
      vh.push_back(PHS("DataMCw_SingleMuon_nj_" +cut.first,  {dcData_SingleMuon("cntNJetsPt30Eta24Zinv"), dcwMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_jpt_" +cut.first,  {dcData_SingleMuon_jpt, dcwMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_j1pt_"  +cut.first,  {dcData_SingleMuon_j1pt,  dcwMC_j1pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_j2pt_"  +cut.first,  {dcData_SingleMuon_j2pt,  dcwMC_j2pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_j3pt_"  +cut.first,  {dcData_SingleMuon_j3pt,  dcwMC_j3pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_mupt_"  +cut.first,  {dcData_SingleMuon_mupt,  dcwMC_mupt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mupt,  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_mu1pt_" +cut.first,  {dcData_SingleMuon_mu1pt, dcwMC_mu1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu1pt,  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_mu2pt_" +cut.first,  {dcData_SingleMuon_mu2pt, dcwMC_mu2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu2pt,  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_elpt_"  +cut.first,  {dcData_SingleMuon_elpt,  dcwMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_el1pt_" +cut.first,  {dcData_SingleMuon_el1pt, dcwMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_el2pt_" +cut.first,  {dcData_SingleMuon_el2pt, dcwMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_mll_" +cut.first,  {dcData_SingleMuon_mll, dcwMC_mll}, {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
      vh.push_back(PHS("DataMCw_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon("nSearchBin"), dcwMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_nb0Bins_" +cut.first,  {dcData_SingleMuon_nb0Bins, dcwMC_nb0Bins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_nb0BinsNW_"  +cut.first,  {dcData_SingleMuon_nb0BinsNW,  dcwMC_nb0BinsNW},  {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMCw_SingleMuon_nb0NJwBins_" +cut.first,  {dcData_SingleMuon_nb0NJwBins, dcwMC_nb0NJwBins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));

      // Normalization weight applied for ttbar only, only for loose0 selections
      if(cut.first.rfind("loose0") == (cut.first.size()-6))
	{
	  // DataMC weights applied
	  vh.push_back(PHS("DataMCwtt_SingleMuon_met_" +cut.first,  {dcData_SingleMuon("cleanMetPt"), dcwttMC("cleanMetPt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
	  vh.push_back(PHS("DataMCwtt_SingleMuon_ht_" +cut.first,  {dcData_SingleMuon("HTZinv"), dcwttMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
	  vh.push_back(PHS("DataMCwtt_SingleMuon_mht_" +cut.first,  {dcData_SingleMuon("cleanMHt"), dcwttMC("cleanMHt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
	  vh.push_back(PHS("DataMCwtt_SingleMuon_nt_" +cut.first,  {dcData_SingleMuon("nTopCandSortedCntZinv"), dcwttMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_nt1b_"  +cut.first,  {dcData_SingleMuon_nt1b,  dcwttMC_nt1b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_nt2b_"  +cut.first,  {dcData_SingleMuon_nt2b,  dcwttMC_nt2b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_nt3b_"  +cut.first,  {dcData_SingleMuon_nt3b,  dcwttMC_nt3b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
	  vh.push_back(PHS("DataMCwtt_SingleMuon_mt2_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv"), dcwttMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_mt21b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")1b, dcwttMC_mt21b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_mt22b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")2b, dcwttMC_mt22b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_mt23b_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv")3b, dcwttMC_mt23b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
	  vh.push_back(PHS("DataMCwtt_SingleMuon_nb_" +cut.first,  {dcData_SingleMuon("cntCSVSZinv"), dcwttMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
	  vh.push_back(PHS("DataMCwtt_SingleMuon_nj_" +cut.first,  {dcData_SingleMuon("cntNJetsPt30Eta24Zinv"), dcwttMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_jpt_" +cut.first,  {dcData_SingleMuon_jpt, dcwttMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_j1pt_"  +cut.first,  {dcData_SingleMuon_j1pt,  dcwttMC_j1pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_j2pt_"  +cut.first,  {dcData_SingleMuon_j2pt,  dcwttMC_j2pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_j3pt_"  +cut.first,  {dcData_SingleMuon_j3pt,  dcwttMC_j3pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_mupt_"  +cut.first,  {dcData_SingleMuon_mupt,  dcwttMC_mupt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mupt,  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_mu1pt_" +cut.first,  {dcData_SingleMuon_mu1pt, dcwttMC_mu1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu1pt,  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_mu2pt_" +cut.first,  {dcData_SingleMuon_mu2pt, dcwttMC_mu2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu2pt,  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_elpt_"  +cut.first,  {dcData_SingleMuon_elpt,  dcwttMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_el1pt_" +cut.first,  {dcData_SingleMuon_el1pt, dcwttMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_el2pt_" +cut.first,  {dcData_SingleMuon_el2pt, dcwttMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_mll_" +cut.first,  {dcData_SingleMuon_mll, dcwttMC_mll}, {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
	  vh.push_back(PHS("DataMCwtt_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon("nSearchBin"), dcwttMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_nb0Bins_" +cut.first,  {dcData_SingleMuon_nb0Bins, dcwttMC_nb0Bins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_nb0BinsNW_"  +cut.first,  {dcData_SingleMuon_nb0BinsNW,  dcwttMC_nb0BinsNW},  {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCwtt_SingleMuon_nb0NJwBins_" +cut.first,  {dcData_SingleMuon_nb0NJwBins, dcwttMC_nb0NJwBins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	 
	}

      // Normalization weight applied, only for blnotag selections
      if(cut.first.rfind("blnotag") == (cut.first.size()-7))
	{
	  // DataMC weights applied
	  vh.push_back(PHS("DataMCww_SingleMuon_met_" +cut.first,  {dcData_SingleMuon("cleanMetPt"), dcwwMC("cleanMetPt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
	  vh.push_back(PHS("DataMCww_SingleMuon_ht_" +cut.first,  {dcData_SingleMuon("HTZinv"), dcwwMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
	  vh.push_back(PHS("DataMCww_SingleMuon_mht_" +cut.first,  {dcData_SingleMuon("cleanMHt"), dcwwMC("cleanMHt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
	  vh.push_back(PHS("DataMCww_SingleMuon_nt_" +cut.first,  {dcData_SingleMuon("nTopCandSortedCntZinv"), dcwwMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_nt1b_"  +cut.first,  {dcData_SingleMuon_nt1b,  dcwwMC_nt1b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_nt2b_"  +cut.first,  {dcData_SingleMuon_nt2b,  dcwwMC_nt2b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_nt3b_"  +cut.first,  {dcData_SingleMuon_nt3b,  dcwwMC_nt3b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
	  vh.push_back(PHS("DataMCww_SingleMuon_mt2_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv"), dcwwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_mt21b_" +cut.first,  {dcData_SingleMuon_mt21b, dcwwMC_mt21b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_mt22b_" +cut.first,  {dcData_SingleMuon_mt22b, dcwwMC_mt22b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_mt23b_" +cut.first,  {dcData_SingleMuon_mt23b, dcwwMC_mt23b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
	  vh.push_back(PHS("DataMCww_SingleMuon_nb_" +cut.first,  {dcData_SingleMuon("cntCSVSZinv"), dcwwMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
	  vh.push_back(PHS("DataMCww_SingleMuon_nj_" +cut.first,  {dcData_SingleMuon("cntNJetsPt30Eta24Zinv"), dcwwMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_jpt_" +cut.first,  {dcData_SingleMuon_jpt, dcwwMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_j1pt_"  +cut.first,  {dcData_SingleMuon_j1pt,  dcwwMC_j1pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_j2pt_"  +cut.first,  {dcData_SingleMuon_j2pt,  dcwwMC_j2pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_j3pt_"  +cut.first,  {dcData_SingleMuon_j3pt,  dcwwMC_j3pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_mupt_"  +cut.first,  {dcData_SingleMuon_mupt,  dcwwMC_mupt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mupt,  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_mu1pt_" +cut.first,  {dcData_SingleMuon_mu1pt, dcwwMC_mu1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu1pt,  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_mu2pt_" +cut.first,  {dcData_SingleMuon_mu2pt, dcwwMC_mu2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu2pt,  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_elpt_"  +cut.first,  {dcData_SingleMuon_elpt,  dcwwMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_el1pt_" +cut.first,  {dcData_SingleMuon_el1pt, dcwwMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_el2pt_" +cut.first,  {dcData_SingleMuon_el2pt, dcwwMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_mll_" +cut.first,  {dcData_SingleMuon_mll, dcwwMC_mll}, {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
	  vh.push_back(PHS("DataMCww_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon("nSearchBin"), dcwwMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_nb0Bins_" +cut.first,  {dcData_SingleMuon_nb0Bins, dcwwMC_nb0Bins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_nb0BinsNW_"  +cut.first,  {dcData_SingleMuon_nb0BinsNW,  dcwwMC_nb0BinsNW},  {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCww_SingleMuon_nb0NJwBins_" +cut.first,  {dcData_SingleMuon_nb0NJwBins, dcwwMC_nb0NJwBins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	 
	}
    }

  //push the histograms in a loop, save some copy-paste time
  //for(std::pair<std::string,std::string>& cut : cutlevels_muon_syst)
  //  {
  //    // DataMC weights applied
  //    vector<double> metBins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 2000};
  //    vector<double> mt2Bins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 2000};
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_met_" +cut.first,  {dcData_SingleMuon("cleanMetPt"), dcwMC("cleanMetPt")}, {1, 2}, cut.second, metBins, true, false,  label_met, ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_ht_" +cut.first,  {dcData_SingleMuon("HTZinv"), dcwMC("HTZinv")}, {1, 2}, cut.second, 150, 0, 1500, true, false,  label_ht,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_nt_" +cut.first,  {dcData_SingleMuon("nTopCandSortedCntZinv"), dcwMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_mt2_" +cut.first,  {dcData_SingleMuon("best_had_brJet_MT2Zinv"), dcwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, mt2Bins, true, false,  label_mt2, ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_nb_" +cut.first,  {dcData_SingleMuon("cntCSVSZinv"), dcwMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_nj_" +cut.first,  {dcData_SingleMuon("cntNJetsPt30Eta24Zinv"), dcwMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_jpt_" +cut.first,  {dcData_SingleMuon_jpt, dcwMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_mupt_"  +cut.first,  {dcData_SingleMuon_mupt,  dcwMC_mupt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mupt,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_SingleMuon_nSearchBin_" +cut.first,  {dcData_SingleMuon("nSearchBin"), dcwMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
  //  }


  //electronics:
  for(std::pair<std::string,std::string>& cut : cutlevels_electron)
    {
      // unweighted
      vh.push_back(PHS("DataMC_DoubleEG_met_" +cut.first,  {dcData_DoubleEG("cleanMetPt"), dcMC("cleanMetPt")},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
      vh.push_back(PHS("DataMC_DoubleEG_ht_"  +cut.first,  {dcData_DoubleEG("HTZinv"), dcMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
      vh.push_back(PHS("DataMC_DoubleEG_mht_" +cut.first,  {dcData_DoubleEG("cleanMHt"), dcMC("cleanMHt")},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
      vh.push_back(PHS("DataMC_DoubleEG_nt_"  +cut.first,  {dcData_DoubleEG("nTopCandSortedCntZinv"), dcMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
      //vh.push_back(PHS("DataMC_DoubleEG_nt1b_" +cut.first,  {dcData_DoubleEG_nt1b,  dcMC_nt1b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
      //vh.push_back(PHS("DataMC_DoubleEG_nt2b_" +cut.first,  {dcData_DoubleEG_nt2b,  dcMC_nt2b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
      //vh.push_back(PHS("DataMC_DoubleEG_nt3b_" +cut.first,  {dcData_DoubleEG_nt3b,  dcMC_nt3b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
      vh.push_back(PHS("DataMC_DoubleEG_mt2_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv"), dcMC("best_had_brJet_MT2Zinv")},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
      //vh.push_back(PHS("DataMC_DoubleEG_mt21b_"  +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")1b, dcMC_mt21b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_mt22b_"  +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")2b, dcMC_mt22b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_mt23b_"  +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")3b, dcMC_mt23b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
      vh.push_back(PHS("DataMC_DoubleEG_nb_"  +cut.first,  {dcData_DoubleEG("cntCSVSZinv"), dcMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
      vh.push_back(PHS("DataMC_DoubleEG_nj_"  +cut.first,  {dcData_DoubleEG("cntNJetsPt30Eta24Zinv"), dcMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_jpt_" +cut.first,  {dcData_DoubleEG_jpt, dcMC_jpt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_j1pt_" +cut.first,  {dcData_DoubleEG_j1pt,  dcMC_j1pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
      //vh.push_back(PHS("DataMC_DoubleEG_j2pt_" +cut.first,  {dcData_DoubleEG_j2pt,  dcMC_j2pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
      //vh.push_back(PHS("DataMC_DoubleEG_j3pt_" +cut.first,  {dcData_DoubleEG_j3pt,  dcMC_j3pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
      //vh.push_back(PHS("DataMC_DoubleEG_elpt_" +cut.first,  {dcData_DoubleEG_elpt,  dcMC_elpt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_el1pt_"  +cut.first,  {dcData_DoubleEG_el1pt, dcMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_el2pt_"  +cut.first,  {dcData_DoubleEG_el2pt, dcMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_elpt_" +cut.first,  {dcData_DoubleEG_elpt,  dcMC_elpt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_el1pt_"  +cut.first,  {dcData_DoubleEG_el1pt, dcMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_el2pt_"  +cut.first,  {dcData_DoubleEG_el2pt, dcMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
      //vh.push_back(PHS("DataMC_DoubleEG_mll_" +cut.first,  {dcData_DoubleEG_mll, dcMC_mll},  {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
      vh.push_back(PHS("DataMC_DoubleEG_nSearchBin_" +cut.first,  {dcData_DoubleEG("nSearchBin"), dcMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMC_DoubleEG_nb0Bins_" +cut.first,  {dcData_DoubleEG_nb0Bins, dcMC_nb0Bins}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMC_DoubleEG_nb0BinsNW_"  +cut.first,  {dcData_DoubleEG_nb0BinsNW,  dcMC_nb0BinsNW},  {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMC_DoubleEG_nb0NJwBins_" +cut.first,  {dcData_DoubleEG_nb0NJwBins, dcMC_nb0NJwBins}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMC_DoubleEG_nb0NJwBins_effAcc_" +cut.first,  {dcData_DoubleEG_nb0NJwBins_effAcc, dcMC_nb0NJwBins_effAcc}, {1, 2}, cut.second, NSB, 0, NSB, true, false,  "Search Bin", ""));

      // DataMC weights applied
      vector<double> metBins = {0, 50, 100, 150, 200, 275, 300, 350, 400, 450, 2000};
      vector<double> mt2Bins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 2000};
      vh.push_back(PHS("DataMCw_DoubleEG_met_" +cut.first,  {dcData_DoubleEG("cleanMetPt"), dcwMC("cleanMetPt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
      vh.push_back(PHS("DataMCw_DoubleEG_rebin_met_" +cut.first,  {dcData_DoubleEG("cleanMetPt"), dcwMC("cleanMetPt")}, {1, 2}, cut.second, metBins, true, false,  label_met, ""));
      vh.push_back(PHS("DataMCw_DoubleEG_ht_" +cut.first,  {dcData_DoubleEG("HTZinv"), dcwMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
      vh.push_back(PHS("DataMCw_DoubleEG_mht_" +cut.first,  {dcData_DoubleEG("cleanMHt"), dcwMC("cleanMHt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_nt_" +cut.first,  {dcData_DoubleEG_nt, dcwMC_nt}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_nt1b_"  +cut.first,  {dcData_DoubleEG_nt1b,  dcwMC_nt1b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_nt2b_"  +cut.first,  {dcData_DoubleEG_nt2b,  dcwMC_nt2b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_nt3b_"  +cut.first,  {dcData_DoubleEG_nt3b,  dcwMC_nt3b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
      vh.push_back(PHS("DataMCw_DoubleEG_mt2_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv"), dcwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_mt21b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")1b, dcwMC("best_had_brJet_MT2Zinv")1b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_mt22b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")2b, dcwMC("best_had_brJet_MT2Zinv")2b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_mt23b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")3b, dcwMC("best_had_brJet_MT2Zinv")3b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
      vh.push_back(PHS("DataMCw_DoubleEG_rebin_mt2_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv"), dcwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, mt2Bins, true, false,  label_mt2, ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_rebin_mt21b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")1b, dcwMC("best_had_brJet_MT2Zinv")1b}, {1, 2}, cut.second, mt2Bins, true, false,  "mt2 (1b fake)",  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_rebin_mt22b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")2b, dcwMC("best_had_brJet_MT2Zinv")2b}, {1, 2}, cut.second, mt2Bins, true, false,  "mt2 (2b fake)",  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_rebin_mt23b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")3b, dcwMC("best_had_brJet_MT2Zinv")3b}, {1, 2}, cut.second, mt2Bins, true, false,  "mt2 (3b fake)",  ""));
      vh.push_back(PHS("DataMCw_DoubleEG_nb_" +cut.first,  {dcData_DoubleEG("cntCSVSZinv"), dcwMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
      vh.push_back(PHS("DataMCw_DoubleEG_nj_" +cut.first,  {dcData_DoubleEG("cntNJetsPt30Eta24Zinv"), dcwMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_jpt_" +cut.first,  {dcData_DoubleEG_jpt, dcwMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_j1pt_"  +cut.first,  {dcData_DoubleEG_j1pt,  dcwMC_j1pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_j2pt_"  +cut.first,  {dcData_DoubleEG_j2pt,  dcwMC_j2pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_j3pt_"  +cut.first,  {dcData_DoubleEG_j3pt,  dcwMC_j3pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_elpt_"  +cut.first,  {dcData_DoubleEG_elpt,  dcwMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_el1pt_" +cut.first,  {dcData_DoubleEG_el1pt, dcwMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_el2pt_" +cut.first,  {dcData_DoubleEG_el2pt, dcwMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_elpt_"  +cut.first,  {dcData_DoubleEG_elpt,  dcwMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_el1pt_" +cut.first,  {dcData_DoubleEG_el1pt, dcwMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_el2pt_" +cut.first,  {dcData_DoubleEG_el2pt, dcwMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_mll_" +cut.first,  {dcData_DoubleEG_mll, dcwMC_mll}, {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
      vh.push_back(PHS("DataMCw_DoubleEG_nSearchBin_" +cut.first,  {dcData_DoubleEG("nSearchBin"), dcwMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_nb0Bins_" +cut.first,  {dcData_DoubleEG_nb0Bins, dcwMC_nb0Bins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_nb0BinsNW_"  +cut.first,  {dcData_DoubleEG_nb0BinsNW,  dcwMC_nb0BinsNW},  {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
      //vh.push_back(PHS("DataMCw_DoubleEG_nb0NJwBins_" +cut.first,  {dcData_DoubleEG_nb0NJwBins, dcwMC_nb0NJwBins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));

      // Normalization weight applied for ttbar only, only for loose0 selections
      if(cut.first.rfind("loose0") == (cut.first.size()-6))
	{
	  // DataMC weights applied
	  vh.push_back(PHS("DataMCwtt_DoubleEG_met_" +cut.first,  {dcData_DoubleEG("cleanMetPt"), dcwttMC("cleanMetPt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
	  vh.push_back(PHS("DataMCwtt_DoubleEG_ht_" +cut.first,  {dcData_DoubleEG("HTZinv"), dcwttMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
	  vh.push_back(PHS("DataMCwtt_DoubleEG_mht_" +cut.first,  {dcData_DoubleEG("cleanMHt"), dcwttMC("cleanMHt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
	  vh.push_back(PHS("DataMCwtt_DoubleEG_nt_" +cut.first,  {dcData_DoubleEG("nTopCandSortedCntZinv"), dcwttMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_nt1b_"  +cut.first,  {dcData_DoubleEG_nt1b,  dcwttMC_nt1b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_nt2b_"  +cut.first,  {dcData_DoubleEG_nt2b,  dcwttMC_nt2b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_nt3b_"  +cut.first,  {dcData_DoubleEG_nt3b,  dcwttMC_nt3b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
	  vh.push_back(PHS("DataMCwtt_DoubleEG_mt2_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv"), dcwttMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_mt21b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")1b, dcwttMC_mt21b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_mt22b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")2b, dcwttMC_mt22b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_mt23b_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")3b, dcwttMC_mt23b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
	  vh.push_back(PHS("DataMCwtt_DoubleEG_nb_" +cut.first,  {dcData_DoubleEG("cntCSVSZinv"), dcwttMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
	  vh.push_back(PHS("DataMCwtt_DoubleEG_nj_" +cut.first,  {dcData_DoubleEG("cntNJetsPt30Eta24Zinv"), dcwttMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_jpt_" +cut.first,  {dcData_DoubleEG_jpt, dcwttMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_j1pt_"  +cut.first,  {dcData_DoubleEG_j1pt,  dcwttMC_j1pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_j2pt_"  +cut.first,  {dcData_DoubleEG_j2pt,  dcwttMC_j2pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_j3pt_"  +cut.first,  {dcData_DoubleEG_j3pt,  dcwttMC_j3pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_elpt_"  +cut.first,  {dcData_DoubleEG_elpt,  dcwttMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_el1pt_" +cut.first,  {dcData_DoubleEG_el1pt, dcwttMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_el2pt_" +cut.first,  {dcData_DoubleEG_el2pt, dcwttMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_elpt_"  +cut.first,  {dcData_DoubleEG_elpt,  dcwttMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_el1pt_" +cut.first,  {dcData_DoubleEG_el1pt, dcwttMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_el2pt_" +cut.first,  {dcData_DoubleEG_el2pt, dcwttMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_mll_" +cut.first,  {dcData_DoubleEG_mll, dcwttMC_mll}, {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
	  vh.push_back(PHS("DataMCwtt_DoubleEG_nSearchBin_" +cut.first,  {dcData_DoubleEG("nSearchBin"), dcwttMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_nb0Bins_" +cut.first,  {dcData_DoubleEG_nb0Bins, dcwttMC_nb0Bins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_nb0BinsNW_"  +cut.first,  {dcData_DoubleEG_nb0BinsNW,  dcwttMC_nb0BinsNW},  {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCwtt_DoubleEG_nb0NJwBins_" +cut.first,  {dcData_DoubleEG_nb0NJwBins, dcwttMC_nb0NJwBins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	 
	}

      // Normalization weight applied, only for blnotag selections
      if(cut.first.rfind("blnotag") == (cut.first.size()-7))
	{
	  // DataMC weights applied
	  vh.push_back(PHS("DataMCww_DoubleEG_met_" +cut.first,  {dcData_DoubleEG("cleanMetPt"), dcwwMC("cleanMetPt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
	  vh.push_back(PHS("DataMCww_DoubleEG_ht_" +cut.first,  {dcData_DoubleEG("HTZinv"), dcwwMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
	  vh.push_back(PHS("DataMCww_DoubleEG_mht_" +cut.first,  {dcData_DoubleEG("cleanMHt"), dcwwMC("cleanMHt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
	  vh.push_back(PHS("DataMCww_DoubleEG_nt_" +cut.first,  {dcData_DoubleEG("nTopCandSortedCntZinv"), dcwwMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_nt1b_"  +cut.first,  {dcData_DoubleEG_nt1b,  dcwwMC_nt1b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_nt2b_"  +cut.first,  {dcData_DoubleEG_nt2b,  dcwwMC_nt2b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_nt3b_"  +cut.first,  {dcData_DoubleEG_nt3b,  dcwwMC_nt3b},  {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
	  vh.push_back(PHS("DataMCww_DoubleEG_mt2_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv"), dcwwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_mt21b_" +cut.first,  {dcData_DoubleEG_mt21b, dcwwMC_mt21b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_mt22b_" +cut.first,  {dcData_DoubleEG_mt22b, dcwwMC_mt22b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_mt23b_" +cut.first,  {dcData_DoubleEG_mt23b, dcwwMC_mt23b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
	  vh.push_back(PHS("DataMCww_DoubleEG_nb_" +cut.first,  {dcData_DoubleEG("cntCSVSZinv"), dcwwMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
	  vh.push_back(PHS("DataMCww_DoubleEG_nj_" +cut.first,  {dcData_DoubleEG("cntNJetsPt30Eta24Zinv"), dcwwMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_jpt_" +cut.first,  {dcData_DoubleEG_jpt, dcwwMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_j1pt_"  +cut.first,  {dcData_DoubleEG_j1pt,  dcwwMC_j1pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_j2pt_"  +cut.first,  {dcData_DoubleEG_j2pt,  dcwwMC_j2pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_j3pt_"  +cut.first,  {dcData_DoubleEG_j3pt,  dcwwMC_j3pt},  {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_elpt_"  +cut.first,  {dcData_DoubleEG_elpt,  dcwwMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_el1pt_" +cut.first,  {dcData_DoubleEG_el1pt, dcwwMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_el2pt_" +cut.first,  {dcData_DoubleEG_el2pt, dcwwMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_elpt_"  +cut.first,  {dcData_DoubleEG_elpt,  dcwwMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_el1pt_" +cut.first,  {dcData_DoubleEG_el1pt, dcwwMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_el2pt_" +cut.first,  {dcData_DoubleEG_el2pt, dcwwMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_mll_" +cut.first,  {dcData_DoubleEG_mll, dcwwMC_mll}, {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
	  vh.push_back(PHS("DataMCww_DoubleEG_nSearchBin_" +cut.first,  {dcData_DoubleEG("nSearchBin"), dcwwMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_nb0Bins_" +cut.first,  {dcData_DoubleEG_nb0Bins, dcwwMC_nb0Bins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_nb0BinsNW_"  +cut.first,  {dcData_DoubleEG_nb0BinsNW,  dcwwMC_nb0BinsNW},  {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	  //vh.push_back(PHS("DataMCww_DoubleEG_nb0NJwBins_" +cut.first,  {dcData_DoubleEG_nb0NJwBins, dcwwMC_nb0NJwBins}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
	 
	}
    }

  ////push the histograms in a loop, save some copy-paste time
  //for(std::pair<std::string,std::string>& cut : cutlevels_electron_syst)
  //  {
  //    // DataMC weights applied
  //    vector<double> metBins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 2000};
  //    vector<double> mt2Bins = {0, 50, 100, 150, 200, 250, 300, 350, 400, 2000};
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_met_" +cut.first,  {dcData_DoubleEG("cleanMetPt"), dcwMC("cleanMetPt")}, {1, 2}, cut.second, metBins, true, false,  label_met, ""));//Elec on all these
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_ht_" +cut.first,  {dcData_DoubleEG("HTZinv"), dcwMC("HTZinv")}, {1, 2}, cut.second, 150, 0, 1500, true, false,  label_ht,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_nt_" +cut.first,  {dcData_DoubleEG("nTopCandSortedCntZinv"), dcwMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_mt2_" +cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv"), dcwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, mt2Bins, true, false,  label_mt2, ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_nb_" +cut.first,  {dcData_DoubleEG("cntCSVSZinv"), dcwMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_nj_" +cut.first,  {dcData_DoubleEG("cntNJetsPt30Eta24Zinv"), dcwMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_jpt_" +cut.first,  {dcData_DoubleEG_jpt, dcwMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_elpt_"  +cut.first,  {dcData_DoubleEG_elpt,  dcwMC_elpt},  {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
  //    vh.push_back(PHS("SystPlots_DataMCw_DoubleEG_nSearchBin_" +cut.first,  {dcData_DoubleEG("nSearchBin"), dcwMC("nSearchBin")}, {1, 2}, cut.second, NSB, 0, NSB,  true, false,  "Search Bin", ""));
  //  }


  // Do no look at electron plots for now. Save some time and space
  /*
    for(std::pair<std::string,std::string>& cut : cutlevels_electron)
    {
    // unweighted
    vh.push_back(PHS("DataMC_DoubleEG_met_"+cut.first, {dcData_DoubleEG("cleanMetPt"), dcMC("cleanMetPt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_"+cut.first, {dcData_DoubleEG("HTZinv"), dcMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_"+cut.first, {dcData_DoubleEG("cleanMHt"), dcMC("cleanMHt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_"+cut.first, {dcData_DoubleEG("nTopCandSortedCntZinv"), dcMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_"+cut.first, {dcData_DoubleEG_nt1b, dcMC_nt1b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_"+cut.first, {dcData_DoubleEG_nt2b, dcMC_nt2b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_"+cut.first, {dcData_DoubleEG_nt3b, dcMC_nt3b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_"+cut.first, {dcData_DoubleEG("best_had_brJet_MT2Zinv"), dcMC_mt2}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_"+cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")1b, dcMC_mt21b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_"+cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")2b, dcMC_mt22b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_"+cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")3b, dcMC_mt23b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_"+cut.first, {dcData_DoubleEG("cntCSVSZinv"), dcMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_"+cut.first, {dcData_DoubleEG("cntNJetsPt30Eta24Zinv"), dcMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_"+cut.first, {dcData_DoubleEG_jpt, dcMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_"+cut.first, {dcData_DoubleEG_j1pt, dcMC_j1pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_"+cut.first, {dcData_DoubleEG_j2pt, dcMC_j2pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_"+cut.first, {dcData_DoubleEG_j3pt, dcMC_j3pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_"+cut.first, {dcData_DoubleEG_mupt, dcMC_mupt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mupt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_"+cut.first,  {dcData_DoubleEG_mu1pt, dcMC_mu1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu1pt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_"+cut.first,  {dcData_DoubleEG_mu2pt, dcMC_mu2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu2pt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_"+cut.first, {dcData_DoubleEG_elpt, dcMC_elpt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_"+cut.first,  {dcData_DoubleEG_el1pt, dcMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_"+cut.first,  {dcData_DoubleEG_el2pt, dcMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_"+cut.first, {dcData_DoubleEG_mll, dcMC_mll}, {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
    // DataMC weights applied
    vh.push_back(PHS("DataMC_DoubleEG_met_"+cut.first, {dcData_DoubleEG("cleanMetPt"), dcwMC("cleanMetPt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_met, ""));
    vh.push_back(PHS("DataMC_DoubleEG_ht_"+cut.first, {dcData_DoubleEG("HTZinv"), dcwMC("HTZinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_ht,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mht_"+cut.first, {dcData_DoubleEG("cleanMHt"), dcwMC("cleanMHt")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mht, ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt_"+cut.first, {dcData_DoubleEG("nTopCandSortedCntZinv"), dcwMC("nTopCandSortedCntZinv")}, {1, 2}, cut.second, 5,  0, 5, true, false,  label_nt, ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt1b_"+cut.first, {dcData_DoubleEG_nt1b, dcwMC_nt1b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (1b fake)", ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt2b_"+cut.first, {dcData_DoubleEG_nt2b, dcwMC_nt2b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (2b fake)", ""));
    vh.push_back(PHS("DataMC_DoubleEG_nt3b_"+cut.first, {dcData_DoubleEG_nt3b, dcwMC_nt3b}, {1, 2}, cut.second, 5,  0, 5, true, false,  "Ntop (3b fake)", ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt2_"+cut.first, {dcData_DoubleEG("best_had_brJet_MT2Zinv"), dcwMC("best_had_brJet_MT2Zinv")}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_mt2, ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt21b_"+cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")1b, dcwMC("best_had_brJet_MT2Zinv")1b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (1b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt22b_"+cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")2b, dcwMC("best_had_brJet_MT2Zinv")2b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (2b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mt23b_"+cut.first,  {dcData_DoubleEG("best_had_brJet_MT2Zinv")3b, dcwMC("best_had_brJet_MT2Zinv")3b}, {1, 2}, cut.second, 50, 0, 1500, true, false,  "mt2 (3b fake)",  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nb_"+cut.first, {dcData_DoubleEG("cntCSVSZinv"), dcwMC("cntCSVSZinv")}, {1, 2}, cut.second, 10, 0, 10, true, false,  label_nb,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_nj_"+cut.first, {dcData_DoubleEG("cntNJetsPt30Eta24Zinv"), dcwMC("cntNJetsPt30Eta24Zinv")}, {1, 2}, cut.second, 20, 0, 20, true, false,  label_nj,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_jpt_"+cut.first, {dcData_DoubleEG_jpt, dcwMC_jpt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_jpt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_j1pt_"+cut.first, {dcData_DoubleEG_j1pt, dcwMC_j1pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j1pt, ""));
    vh.push_back(PHS("DataMC_DoubleEG_j2pt_"+cut.first, {dcData_DoubleEG_j2pt, dcwMC_j2pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j2pt, ""));
    vh.push_back(PHS("DataMC_DoubleEG_j3pt_"+cut.first, {dcData_DoubleEG_j3pt, dcwMC_j3pt}, {1, 2}, cut.second, 50, 0, 1500, true, false,  label_j3pt, ""));
    vh.push_back(PHS("DataMC_DoubleEG_mupt_"+cut.first, {dcData_DoubleEG_mupt, dcwMC_mupt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mupt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu1pt_"+cut.first,  {dcData_DoubleEG_mu1pt, dcwMC_mu1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu1pt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mu2pt_"+cut.first,  {dcData_DoubleEG_mu2pt, dcwMC_mu2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_mu2pt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_elpt_"+cut.first, {dcData_DoubleEG_elpt, dcwMC_elpt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_elpt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el1pt_"+cut.first,  {dcData_DoubleEG_el1pt, dcwMC_el1pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el1pt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_el2pt_"+cut.first,  {dcData_DoubleEG_el2pt, dcwMC_el2pt}, {1, 2}, cut.second, 50, 0, 1000, true, false,  label_el2pt,  ""));
    vh.push_back(PHS("DataMC_DoubleEG_mll_"+cut.first, {dcData_DoubleEG_mll, dcwMC_mll}, {1, 2}, cut.second, 40, 0, 200,  true, false,  label_mll, ""));
    }
  */

  //interlude for MC checks

  std::vector<std::pair<std::string,std::string>> cutlevels_MC = {
    {"nosel", "passNoiseEventFilterZinv"},
    {"baseline",  "passNoiseEventFilterZinv;passBaselineNoTagZinv;!nTopCandSortedCntZinv=0"},
    {"loose0", "passNoiseEventFilterZinv;passLeptVeto;HTZinv>200;passnJetsZinv;passdPhisZinv"},
    {"loose50", "passNoiseEventFilterZinv;passLeptVeto;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>50"},
    {"loose100",  "passNoiseEventFilterZinv;passLeptVeto;HTZinv>200;passnJetsZinv;passdPhisZinv;cleanMetPt>100"},
    {"loose200",  "passNoiseEventFilterZinv;passLeptVeto;HTZinv>200;passnJetsZinv;passdPhisZinv;passMETZinv"}
  };

  std::string s_znunu_baseline  = "passBaselineNoTagZinv;!nTopCandSortedCntZinv=0";  //THIS IS NOT PERFECT, clearly nTopCandSortedCntZinv and MT2 cuts need to be specialized by sample

  Plotter::DatasetSummary dsDY_nunu_SB0b_Wgt1b( "Z#rightarrow#nu#nu, N(b) = 0, 1 fake", fileMap["ZJetsToNuNu"], "cntCSVSZinv=0", "bTagSF_EventWeightSimple_Central;nJet1bfakeWgt");
  Plotter::DatasetSummary dsDY_nunu_SB0b_Wgt2b( "Z#rightarrow#nu#nu, N(b) = 0, 2 fake", fileMap["ZJetsToNuNu"], "cntCSVSZinv=0", "bTagSF_EventWeightSimple_Central;nJet2bfakeWgt");
  Plotter::DatasetSummary dsDY_nunu_SB0b_Wgt3b( "Z#rightarrow#nu#nu, N(b) = 0, 3 fake", fileMap["ZJetsToNuNu"], "cntCSVSZinv=0", "bTagSF_EventWeightSimple_Central;nJet3bfakeWgt");
  Plotter::DatasetSummary dsDY_nunu_SB0b_2( "Z#rightarrow#nu#nu, N(b) = 0",  fileMap["ZJetsToNuNu"], "cntCSVSZinv=0", "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dsDY_nunu_SB1b( "Z#rightarrow#nu#nu, Direct MC, N(b) = 1",  fileMap["ZJetsToNuNu"], "cntCSVSZinv=1", "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dsDY_nunu_SB2b( "Z#rightarrow#nu#nu, Direct MC, N(b) = 2",  fileMap["ZJetsToNuNu"], "cntCSVSZinv=2", "bTagSF_EventWeightSimple_Central");
  Plotter::DatasetSummary dsDY_nunu_SB3b( "Z#rightarrow#nu#nu, Direct MC, N(b) >= 3", fileMap["ZJetsToNuNu"], "cntCSVSZinv>2", "bTagSF_EventWeightSimple_Central");

  Plotter::DatasetSummary dsDY_nunu_njet( "Z#rightarrow#nu#nu DataMC weight",  fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ");
  Plotter::DatasetSummary dsDY_nunu_njetnorm( "Z#rightarrow#nu#nu Njet+norm weight", fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_bl(  "Z#rightarrow#nu#nu Njet+norm weight", fileMap["ZJetsToNuNu"], "passLeptVeto;passBaselineZinv",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_scaleUp( "Z#rightarrow#nu#nu Scale weight Up", fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b;Scaled_Variations_Up");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_scaleDown(  "Z#rightarrow#nu#nu Scale weight Down",  fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b;Scaled_Variations_Down");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_PDFUp(  "Z#rightarrow#nu#nu PDF weight Up",  fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b;NNPDF_From_Median_Up");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_PDFDown( "Z#rightarrow#nu#nu PDF weight Down", fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b;NNPDF_From_Median_Down");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_JEUUp(  "Z#rightarrow#nu#nu JEC Up",  fileMap["ZJetsToNuNu"], "passLeptVeto;passBaselineZinvJEUUp", "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_JEUDown( "Z#rightarrow#nu#nu JEC Down", fileMap["ZJetsToNuNu"], "passLeptVeto;passBaselineZinvJEUDn", "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_MEUUp(  "Z#rightarrow#nu#nu MEC Up",  fileMap["ZJetsToNuNu"], "passLeptVeto;passBaselineZinvMEUUp", "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_MEUDown( "Z#rightarrow#nu#nu MEC Down", fileMap["ZJetsToNuNu"], "passLeptVeto;passBaselineZinvMEUDn", "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerCentral( "Z#rightarrow#nu#nu Trigger weight Central", fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b;TriggerEffMC");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerUp(  "Z#rightarrow#nu#nu Trigger weight Up",  fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b;TriggerEffUpMC");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_TriggerDown( "Z#rightarrow#nu#nu Trigger weight Down", fileMap["ZJetsToNuNu"], "passLeptVeto",  "bTagSF_EventWeightSimple_Central;nJetWgtDYZ;normWgt0b;TriggerEffDownMC");

  Plotter::DatasetSummary dsDY_nunu_njetnorm_btagUp( "Z#rightarrow#nu#nu Trigger weight Central", fileMap["ZJetsToNuNu"], "passLeptVeto", "bTagSF_EventWeightSimple_Up;nJetWgtDYZ;normWgt0b;TriggerEffMC");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_btagDn( "Z#rightarrow#nu#nu Trigger weight Central", fileMap["ZJetsToNuNu"], "passLeptVeto", "bTagSF_EventWeightSimple_Down;nJetWgtDYZ;normWgt0b;TriggerEffMC");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_bmistagUp( "Z#rightarrow#nu#nu Trigger weight Central", fileMap["ZJetsToNuNu"], "passLeptVeto", "mistagSF_EventWeightSimple_Up;nJetWgtDYZ;normWgt0b;TriggerEffMC");
  Plotter::DatasetSummary dsDY_nunu_njetnorm_bmistagDn( "Z#rightarrow#nu#nu Trigger weight Central", fileMap["ZJetsToNuNu"], "passLeptVeto", "mistagSF_EventWeightSimple_Down;nJetWgtDYZ;normWgt0b;TriggerEffMC");


  Plotter::DatasetSummary dsTT_inc( "t#bar{t} Inc", fileMap["TTbar"], "", "");

  //norm plots
  Plotter::DataCollection dcMC_nunu_nj_1b( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_2, dsDY_nunu_SB1b});
  Plotter::DataCollection dcMC_nunu_nj_2b( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_2, dsDY_nunu_SB2b});
  Plotter::DataCollection dcMC_nunu_nj_3b( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_2, dsDY_nunu_SB3b});

  //1 fake b
  Plotter::DataCollection dcMC_nunu_Wgt1b_met("single", "cleanMetPt", {dsDY_nunu_SB0b_Wgt1b, dsDY_nunu_SB1b});
  Plotter::DataCollection dcMC_nunu_Wgt1b_ht( "single", "HTZinv", {dsDY_nunu_SB0b_Wgt1b, dsDY_nunu_SB1b});
  Plotter::DataCollection dcMC_nunu_Wgt1b_nb( "single", "cntCSVSZinv", {dsDY_nunu_SB0b_Wgt1b, dsDY_nunu_SB1b});
  Plotter::DataCollection dcMC_nunu_Wgt1b_nj( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_Wgt1b, dsDY_nunu_SB1b});
  Plotter::DataCollection dcMC_nunu_Wgt1b_nt( "single", {{ "nTopCandSortedCntZinv1b", dsDY_nunu_SB0b_Wgt1b}, { "nTopCandSortedCntZinv", dsDY_nunu_SB1b}});
  Plotter::DataCollection dcMC_nunu_Wgt1b_mt2("single", {{"best_had_brJet_MT2Zinv1b", dsDY_nunu_SB0b_Wgt1b}, {"best_had_brJet_MT2Zinv", dsDY_nunu_SB1b}});
  //2 fake b
  Plotter::DataCollection dcMC_nunu_Wgt2b_met("single", "cleanMetPt", {dsDY_nunu_SB0b_Wgt2b, dsDY_nunu_SB2b});
  Plotter::DataCollection dcMC_nunu_Wgt2b_ht( "single", "HTZinv", {dsDY_nunu_SB0b_Wgt2b, dsDY_nunu_SB2b});
  Plotter::DataCollection dcMC_nunu_Wgt2b_nb( "single", "cntCSVSZinv", {dsDY_nunu_SB0b_Wgt2b, dsDY_nunu_SB2b});
  Plotter::DataCollection dcMC_nunu_Wgt2b_nj( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_Wgt2b, dsDY_nunu_SB2b});
  Plotter::DataCollection dcMC_nunu_Wgt2b_nt( "single", {{ "nTopCandSortedCntZinv2b", dsDY_nunu_SB0b_Wgt2b}, { "nTopCandSortedCntZinv", dsDY_nunu_SB2b}});
  Plotter::DataCollection dcMC_nunu_Wgt2b_mt2("single", {{"best_had_brJet_MT2Zinv2b", dsDY_nunu_SB0b_Wgt2b}, {"best_had_brJet_MT2Zinv", dsDY_nunu_SB2b}});
  //3 fake b
  Plotter::DataCollection dcMC_nunu_Wgt3b_met("single", "cleanMetPt", {dsDY_nunu_SB0b_Wgt3b, dsDY_nunu_SB3b});
  Plotter::DataCollection dcMC_nunu_Wgt3b_ht( "single", "HTZinv", {dsDY_nunu_SB0b_Wgt3b, dsDY_nunu_SB3b});
  Plotter::DataCollection dcMC_nunu_Wgt3b_nb( "single", "cntCSVSZinv", {dsDY_nunu_SB0b_Wgt3b, dsDY_nunu_SB3b});
  Plotter::DataCollection dcMC_nunu_Wgt3b_nj( "single", "cntNJetsPt30Eta24Zinv", {dsDY_nunu_SB0b_Wgt3b, dsDY_nunu_SB3b});
  Plotter::DataCollection dcMC_nunu_Wgt3b_nt( "single", {{ "nTopCandSortedCntZinv3b", dsDY_nunu_SB0b_Wgt3b}, { "nTopCandSortedCntZinv", dsDY_nunu_SB3b}});
  Plotter::DataCollection dcMC_nunu_Wgt3b_mt2("single", {{"best_had_brJet_MT2Zinv3b", dsDY_nunu_SB0b_Wgt3b}, {"best_had_brJet_MT2Zinv", dsDY_nunu_SB3b}});

  Plotter::DataCollection dcMC_nunu_nSearchBin( "single", {{"nSearchBin", dsDY_nunu},  {"nb0Bins", dsDY_nunu},  {"nb0BinsNW", dsDY_nunu} });
  Plotter::DataCollection dcMC_nunu_nSearchBin_njW("single", {{"nSearchBin", dsDY_nunu_njet}, {"nb0NJwBins", dsDY_nunu_njet}, {"nb0BinsNW", dsDY_nunu_njet} });

  Plotter::DataCollection njetw_nSearchBin( "single", {{"nSearchBin", dsDY_nunu_njet}, {"nSearchBin", dsDY_nunu}  });
  Plotter::DataCollection scalew_nSearchBin( "single", {{"nSearchBin", dsDY_nunu_njetnorm}, {"nSearchBin", dsDY_nunu_njetnorm_scaleUp}, {"nSearchBin", dsDY_nunu_njetnorm_scaleDown}  });
  Plotter::DataCollection scalew_nj( "single", {{"cntNJetsPt30Eta24Zinv", dsDY_nunu_njetnorm}, {"cntNJetsPt30Eta24Zinv", dsDY_nunu_njetnorm_scaleUp}, {"cntNJetsPt30Eta24Zinv", dsDY_nunu_njetnorm_scaleDown}  });
  Plotter::DataCollection pdfw_nSearchBin( "single", {{"nSearchBin", dsDY_nunu_njetnorm}, {"nSearchBin", dsDY_nunu_njetnorm_PDFUp}, {"nSearchBin", dsDY_nunu_njetnorm_PDFDown} });
  Plotter::DataCollection pdfw_nj( "single", {{"cntNJetsPt30Eta24Zinv", dsDY_nunu_njetnorm}, {"cntNJetsPt30Eta24Zinv", dsDY_nunu_njetnorm_PDFUp}, {"cntNJetsPt30Eta24Zinv", dsDY_nunu_njetnorm_PDFDown} });
  Plotter::DataCollection JEU_nSearchBin(  "single", {{"nSearchBin", dsDY_nunu_njetnorm_bl}, {"nSearchBinJEUUp", dsDY_nunu_njetnorm_JEUUp}, {"nSearchBinJEUDn", dsDY_nunu_njetnorm_JEUDown}  });
  Plotter::DataCollection MEU_nSearchBin(  "single", {{"nSearchBin", dsDY_nunu_njetnorm_bl}, {"nSearchBinMEUUp", dsDY_nunu_njetnorm_MEUUp}, {"nSearchBinMEUDn", dsDY_nunu_njetnorm_MEUDown}  });
  Plotter::DataCollection trigger_nSearchBin( "single", {{"nSearchBin", dsDY_nunu_njetnorm_TriggerCentral}, {"nSearchBin", dsDY_nunu_njetnorm_TriggerUp}, {"nSearchBin", dsDY_nunu_njetnorm_TriggerDown}, {"nSearchBin", dsDY_nunu_njetnorm}  });
  Plotter::DataCollection btag_nSearchBin( "single", {{"nSearchBin", dsDY_nunu_njetnorm_TriggerCentral}, {"nSearchBin", dsDY_nunu_njetnorm_btagUp}, {"nSearchBin", dsDY_nunu_njetnorm_btagDn}});
  Plotter::DataCollection bmistag_nSearchBin( "single", {{"nSearchBin", dsDY_nunu_njetnorm_TriggerCentral}, {"nSearchBin", dsDY_nunu_njetnorm_bmistagUp}, {"nSearchBin", dsDY_nunu_njetnorm_bmistagDn}});


  vh.push_back(PHS("NJetWgt_nSearchBin", {njetw_nSearchBin}, {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("NJetWgt_nSearchBin_log", {njetw_nSearchBin}, {2, 1}, "passBaselineZinv", NSB,  0, NSB, true,  false,  "Search Bin", "Events", true));
  vh.push_back(PHS("NJetWgt_nSearchBin_pull", {njetw_nSearchBin}, {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", false));
  vh.push_back(PHS("NJetWgt_nSearchBin_pull_log", {njetw_nSearchBin}, {2, 1}, "passBaselineZinv", NSB,  0, NSB, true,  false,  "Search Bin", "Events", false));
  vh.push_back(PHS("syst_ScaleWgt_nSearchBin",  {scalew_nSearchBin},{2, 1}, "passBaselineZinv", NSB,  0, NSB, false, true, "Search Bin", "Events", true));
  vh.push_back(PHS("syst_PDFWgt_nSearchBin", {pdfw_nSearchBin},  {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, true, "Search Bin", "Events", true));
  vh.push_back(PHS("syst_JESUncert_nSearchBin", {JEU_nSearchBin}, {2, 1}, "", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("syst_MESUncert_nSearchBin", {MEU_nSearchBin}, {2, 1}, "", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("TriggerWgt_nSearchBin",  {trigger_nSearchBin},  {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  //b-tagging uncertainty
  vh.push_back(PHS("BTagUncert_nSearchBin",  {btag_nSearchBin}, {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  vh.push_back(PHS("BMistagUncert_nSearchBin",  {bmistag_nSearchBin},  {2, 1}, "passBaselineZinv", NSB,  0, NSB, false, false,  "Search Bin", "Events", true));
  // PDF and scale
  vh.push_back(PHS("DataMCwwscale_SingleMuon_nj_0b_blnotag",  {scalew_nj}, {1, 2}, "cntCSVSZinv=0;passBaselineNoTagZinv", 20, 0, 20, true, false,  label_nj, ""));
  vh.push_back(PHS("DataMCwwpdf_SingleMuon_nj_0b_blnotag", {pdfw_nj},  {1, 2}, "cntCSVSZinv=0;passBaselineNoTagZinv", 20, 0, 20, true, false,  label_nj, ""));

  vh.push_back(PHS("btagWeights", {PDC("single", "bTagSF_EventWeightSimple_Central", {dsDY, dstt2l, dsDY_nunu_njetnorm_TriggerCentral})},{1, 1}, "",  100, 0, 5,  false, false,  "b-tag correction weight", ""));
  vh.push_back(PHS("btagWeights_bl", {PDC("single", "bTagSF_EventWeightSimple_Central", {dsDY, dstt2l, dsDY_nunu_njetnorm_TriggerCentral})},{1, 1}, "passBaselineZinv", 100, 0, 5,  false, false,  "b-tag correction weight", ""));

  //Baseline cuts are flawed here, fix!
  for(std::pair<std::string,std::string>& cut : cutlevels_MC)
    {
      vh.push_back(PHS(  "ClosureNb_nj_nw_1fakeb_" + cut.first, {dcMC_nunu_nj_1b},  {1, 2}, cut.second,  20, 0, 20, true,  false,  label_nj, ""));
      vh.push_back(PHS(  "ClosureNb_nj_nw_2fakeb_" + cut.first, {dcMC_nunu_nj_2b},  {1, 2}, cut.second,  20, 0, 20, true,  false,  label_nj, ""));
      vh.push_back(PHS(  "ClosureNb_nj_nw_3fakeb_" + cut.first, {dcMC_nunu_nj_3b},  {1, 2}, cut.second,  20, 0, 20, true,  false,  label_nj, ""));
      vh.push_back(PHS( "ClosureNb_met_1fakeb_" + cut.first, {dcMC_nunu_Wgt1b_met},  {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_met,  ""));
      vh.push_back(PHS( "ClosureNb_ht_1fakeb_" + cut.first, {dcMC_nunu_Wgt1b_ht}, {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_ht, ""));
      vh.push_back(PHS( "ClosureNb_nt_1fakeb_" + cut.first, {dcMC_nunu_Wgt1b_nt}, {1, 2}, cut.second, 5, 0, 5, true,  false,  label_nt, ""));
      vh.push_back(PHS( "ClosureNb_nb_1fakeb_" + cut.first, {dcMC_nunu_Wgt1b_nb}, {1, 2}, cut.second,  10, 0, 10, true,  false,  label_nb, ""));
      vh.push_back(PHS( "ClosureNb_mt2_1fakeb_" + cut.first, {dcMC_nunu_Wgt1b_mt2},  {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_mt2,  ""));
      vh.push_back(PHS( "ClosureNb_nj_1fakeb_" + cut.first, {dcMC_nunu_Wgt1b_nj}, {1, 2}, cut.second,  20, 0, 20, true,  false,  label_nj, ""));
      vh.push_back(PHS( "ClosureNb_met_2fakeb_" + cut.first, {dcMC_nunu_Wgt2b_met},  {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_met,  ""));
      vh.push_back(PHS( "ClosureNb_ht_2fakeb_" + cut.first, {dcMC_nunu_Wgt2b_ht}, {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_ht, ""));
      vh.push_back(PHS( "ClosureNb_nt_2fakeb_" + cut.first, {dcMC_nunu_Wgt2b_nt}, {1, 2}, cut.second, 5, 0, 5, true,  false,  label_nt, ""));
      vh.push_back(PHS( "ClosureNb_nb_2fakeb_" + cut.first, {dcMC_nunu_Wgt2b_nb}, {1, 2}, cut.second,  10, 0, 10, true,  false,  label_nb, ""));
      vh.push_back(PHS( "ClosureNb_mt2_2fakeb_" + cut.first, {dcMC_nunu_Wgt2b_mt2},  {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_mt2,  ""));
      vh.push_back(PHS( "ClosureNb_nj_2fakeb_" + cut.first, {dcMC_nunu_Wgt2b_nj}, {1, 2}, cut.second,  20, 0, 20, true,  false,  label_nj, ""));
      vh.push_back(PHS( "ClosureNb_met_3fakeb_" + cut.first, {dcMC_nunu_Wgt3b_met},  {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_met,  ""));
      vh.push_back(PHS( "ClosureNb_ht_3fakeb_" + cut.first, {dcMC_nunu_Wgt3b_ht}, {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_ht, ""));
      vh.push_back(PHS( "ClosureNb_nt_3fakeb_" + cut.first, {dcMC_nunu_Wgt3b_nt}, {1, 2}, cut.second, 5, 0, 5, true,  false,  label_nt, ""));
      vh.push_back(PHS( "ClosureNb_nb_3fakeb_" + cut.first, {dcMC_nunu_Wgt3b_nb}, {1, 2}, cut.second,  10, 0, 10, true,  false,  label_nb, ""));
      vh.push_back(PHS( "ClosureNb_mt2_3fakeb_" + cut.first, {dcMC_nunu_Wgt3b_mt2},  {1, 2}, cut.second,  50, 0, 1500, true,  false,  label_mt2,  ""));
      vh.push_back(PHS( "ClosureNb_nj_3fakeb_" + cut.first, {dcMC_nunu_Wgt3b_nj}, {1, 2}, cut.second,  20, 0, 20, true,  false,  label_nj, ""));
      vh.push_back(PHS( "ClosureNb_nSearchBin_" + cut.first, {dcMC_nunu_nSearchBin}, {1, 2}, cut.second,  NSB, 0, NSB, true,  false,  "Search Bin", "Events"));
      vh.push_back(PHS("ClosureNb_nSearchBin_njW_" + cut.first, {dcMC_nunu_nSearchBin_njW}, {1, 2}, cut.second,  NSB, 0, NSB, true,  false,  "Search Bin", "Events"));
    }
  //MC interlude over



  // other plots for Zhenbin
  Plotter::DataCollection njetww_nj30( "single", "cntNJetsPt30Eta24", {dsDY_nunu_njetnorm});
  Plotter::DataCollection njetww_nj50( "single", "cntNJetsPt50Eta24", {dsDY_nunu_njetnorm});
  Plotter::DataCollection njetww_nt( "single", "nTopCandSortedCnt", {dsDY_nunu_njetnorm});
  Plotter::DataCollection njetww_nb( "single", "cntCSVS", {dsDY_nunu_njetnorm});
  Plotter::DataCollection njetww_met( "single", "met", {dsDY_nunu_njetnorm});
  Plotter::DataCollection njetww_mt2( "single", "best_had_brJet_MT2", {dsDY_nunu_njetnorm});
  Plotter::DataCollection njetww_ht( "single", "HT",  {dsDY_nunu_njetnorm});
  vh.push_back(PHS("NJetWgt_nj30",  {njetww_nj30}, {1, 1}, "passBaseline", 10,  0,  10,  false, false,  "nj", "Events"));
  vh.push_back(PHS("NJetWgt_nj50",  {njetww_nj50}, {1, 1}, "passBaseline", 10,  0,  10,  false, false,  "nj50", "Events"));
  vh.push_back(PHS("NJetWgt_nt", {njetww_nt}, {1, 1}, "passBaseline", 5,  0,  5, false, false,  "nt", "Events"));
  vh.push_back(PHS("NJetWgt_nb", {njetww_nb}, {1, 1}, "passBaseline", 5,  0,  5, false, false,  "nb", "Events"));
  vh.push_back(PHS("NJetWgt_met", {njetww_met},  {1, 1}, "passBaseline", 24,  200,  800, false, false,  "met",  "Events"));
  vh.push_back(PHS("NJetWgt_mt2", {njetww_mt2},  {1, 1}, "passBaseline", 24,  200,  800, false, false,  "mt2",  "Events"));
  vh.push_back(PHS("NJetWgt_ht", {njetww_ht}, {1, 1}, "passBaseline", 20,  500,  1000,  false, false,  "ht", "Events"));

  Plotter::DataCollection njetw_nj30( "single", "cntNJetsPt30Eta24", {dsDY_nunu_njet});
  Plotter::DataCollection njetw_nj50( "single", "cntNJetsPt50Eta24", {dsDY_nunu_njet});
  Plotter::DataCollection njetw_nt( "single", "nTopCandSortedCnt", {dsDY_nunu_njet});
  Plotter::DataCollection njetw_nb( "single", "cntCSVS", {dsDY_nunu_njet});
  Plotter::DataCollection njetw_met( "single", "met", {dsDY_nunu_njet});
  Plotter::DataCollection njetw_mt2( "single", "best_had_brJet_MT2", {dsDY_nunu_njet});
  Plotter::DataCollection njetw_ht( "single", "HT",  {dsDY_nunu_njet});
  std::string cutlevel_Zhenbin = "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets;passMET;passHT";
  vh.push_back(PHS("NJetWgt_ht_nj30",  {njetw_nj30}, {1, 1}, cutlevel_Zhenbin, 10,  0,  10,  false, false,  "nj", "Events"));
  vh.push_back(PHS("NJetWgt_ht_nj50",  {njetw_nj50}, {1, 1}, cutlevel_Zhenbin, 10,  0,  10,  false, false,  "nj50", "Events"));
  vh.push_back(PHS("NJetWgt_ht_nt", {njetw_nt}, {1, 1}, cutlevel_Zhenbin, 5,  0,  5, false, false,  "nt", "Events"));
  vh.push_back(PHS("NJetWgt_ht_nb", {njetw_nb}, {1, 1}, cutlevel_Zhenbin, 5,  0,  5, false, false,  "nb", "Events"));
  vh.push_back(PHS("NJetWgt_ht_met", {njetw_met},  {1, 1}, cutlevel_Zhenbin, 24,  200,  800, false, false,  "met",  "Events"));
  vh.push_back(PHS("NJetWgt_ht_ht", {njetw_ht}, {1, 1}, cutlevel_Zhenbin, 20,  500,  1000,  false, false,  "ht", "Events"));

  // top kinematics
  Plotter::DataCollection njetw_toppt( "single", "vTops[0](pt)",  {dsDY_nunu_njet});
  Plotter::DataCollection njetw_topeta( "single", "vTops[0](eta)", {dsDY_nunu_njet});
  Plotter::DataCollection njetw_topphi( "single", "vTops[0](phi)", {dsDY_nunu_njet});
  Plotter::DataCollection njetw_topmass(  "single", "vTops[0](M)", {dsDY_nunu_njet});
  std::string cutlevel_Zhenbin2 = "passNoiseEventFilter;passnJets;passMuonVeto;passEleVeto;passIsoTrkVeto;passdPhis;passBJets;passMET;passHT;nTopCandSortedCnt=1";
  vh.push_back(PHS("NJetWgt_ht_toppt", {njetw_toppt}, {1, 1}, cutlevel_Zhenbin2, 50,  0,  1000, false, false,  "top pt", "Events"));
  vh.push_back(PHS("NJetWgt_ht_topeta", {njetw_topeta}, {1, 1}, cutlevel_Zhenbin2, 20,  -5,  5,  false, false,  "top eta", "Events"));
  vh.push_back(PHS("NJetWgt_ht_topphi", {njetw_topphi}, {1, 1}, cutlevel_Zhenbin2, 20,  -5,  5,  false, false,  "top phi", "Events"));
  vh.push_back(PHS("NJetWgt_ht_topmass",  {njetw_topmass}, {1, 1}, cutlevel_Zhenbin2, 40,  80,  280, false, false,  "top mass", "Events"));

  //Generate cutflows
  vector<string> cfsZ = {"",
			 "passNoiseEventFilterZinv",
			 "passNoiseEventFilterZinv;passLeptVeto",
			 "passNoiseEventFilterZinv;passLeptVeto",
			 "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv",
			 "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv",
			 "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv",
			 "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv",
			 "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv",
			 "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
			 "passNoiseEventFilterZinv;passLeptVeto;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
			 "passLeptVeto;passBaselineZinv"};

  vector<string> cfsDYmm = {"",
			    "passNoiseEventFilterZinv",
			    "passNoiseEventFilterZinv;passMuZinvSel",
			    "passNoiseEventFilterZinv;passMuZinvSel",
			    "passNoiseEventFilterZinv;passMuZinvSel;passnJetsZinv",
			    "passNoiseEventFilterZinv;passMuZinvSel;passnJetsZinv;passdPhisZinv",
			    "passNoiseEventFilterZinv;passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv",
			    "passNoiseEventFilterZinv;passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv",
			    "passNoiseEventFilterZinv;passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv",
			    "passNoiseEventFilterZinv;passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
			    "passNoiseEventFilterZinv;passMuZinvSel;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
			    "passMuZinvSel;passBaselineZinv"};

  vector<string> cfsDatamm = {"",
			      "passNoiseEventFilterZinv",
			      "passNoiseEventFilterZinv;passMuZinvSel",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
			      "passNoiseEventFilterZinv;passMuZinvSel;passMuTrigger;passBaselineZinv"};

  vector<string> cfsDataem = {"",
			      "passNoiseEventFilterZinv",
			      "passNoiseEventFilterZinv;passElMuZinvSel",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passnJetsZinv;passdPhisZinv;passHTZinv;passMETZinv;passBJetsZinv;passTaggerZinv;passMT2Zinv",
			      "passNoiseEventFilterZinv;passElMuZinvSel;passMuTrigger;passBaselineZinv"};

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

  cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("ZtoNuNu", PDC("", "", {dsDY_nunu}), cfsZ));
  cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("DYtoMuMu",  PDC("", "", {dsDY_ll_zAcc_scaled}), cfsDYmm));
  cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("DYtoMuMu_Unscaled", PDC("", "", {dsDY_ll}),  cfsDYmm));
  cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Data_MuMu",  PDC("", "", {dsData_SingleMuon}), cfsDatamm));
  cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("Data_ElMu",  PDC("", "", {dsData_SingleMuon}), cfsDataem));
  cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("TTBar_ElMu", PDC("", "", {dstt2l}), cfsDataem));
  cutFlowSummaries.emplace_back(Plotter::CutFlowSummary("TTBar_Inc",  PDC("", "", {dsTT_inc}), cfsSync));

  set<AnaSamples::FileSummary> vvf;
  for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

  RegisterFunctions* rf = new RegisterFunctionsNTuple(runOnCondor);

  Plotter plotter(vh, vvf, fromTuple, histFile, nFiles, startFile, nEvts);
  plotter.setCutFlows(cutFlowSummaries);
  plotter.setLumi(lumi);
  plotter.setPlotDir(plotDir);
  plotter.setDoHists(doSave || doPlots);
  plotter.setDoTuple(doTuple);
  plotter.setRegisterFunction(rf);
  plotter.read();
  if(doSave && fromTuple)  plotter.saveHists();
  if(doPlots) plotter.plot();
}
