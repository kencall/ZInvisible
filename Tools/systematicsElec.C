#include "Plotter.h"
#include "RegisterFunctions.h"
#include "SusyAnaTools/Tools/samples.h"
#include "Systematic.h"

#include <getopt.h>
#include <iostream>

int main(int argc, char* argv[])
{
    using namespace std;

    TH1::AddDirectory(false);

    int opt;
    int option_index = 0;
    static struct option long_options[] = {
        {"plot",             no_argument, 0, 'p'},
        {"dataSets",   required_argument, 0, 'D'},
	{"luminosity", required_argument, 0, 'L'}
    };

    string dataSets = "ZJetsToNuNu";
    double lumi = AnaSamples::luminosity;

    while((opt = getopt_long(argc, argv, "pD:L:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
        case 'p':
            break;

        case 'D':
            dataSets = optarg;
            break;

	case 'L':
            lumi = atof(optarg);
            break;
        }
    }

    AnaSamples::SampleSet        ss("/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/ZInvisible/Tools/condor/", lumi);
    AnaSamples::SampleCollection sc(ss);

    map<string, vector<AnaSamples::FileSummary>> fileMap;

    //Select approperiate datasets here
    if(dataSets.compare("TEST") == 0)
    {
        fileMap["DYJetsToLL"]  = {sc["DYJetsToLL"]};
        fileMap["ZJetsToNuNu"] = {sc["ZJetsToNuNu"]};
        fileMap["TTbarDiLep"] = {ss["TTbarDiLep"]};
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

    // for(auto& vfs : fileMap)
    // {
    //     for(auto& fs : vfs.second)
    //     {
    //         size_t start = fs.filePath.rfind('/');
    //         size_t stop  = fs.filePath.rfind('.');
    //         if(int(stop) - (int(start) + 1) > 0)
    //         {
    //             fs.treePath = fs.filePath.substr(start + 1, stop - start - 1);
    //             fs.filePath = "stupidFileList.txt";
    //             fs.readFileList();
    //         }
    //     }
    // }

    TH1 *shapeMET, *shapeMT2, *shapeNT, *shapeNB;

    TH1 *shapeMETElec, *shapeMT2Elec, *shapeNTElec, *shapeNBElec;
    
    TFile *f = new TFile("syst_shape.root");
    if(f)
    {
        shapeMET = static_cast<TH1*>(f->Get("ShapeRatio_met")->Clone());
        shapeMT2 = static_cast<TH1*>(f->Get("ShapeRatio_mt2")->Clone());
        shapeNT  = static_cast<TH1*>(f->Get("ShapeRatio_nt")->Clone());
        shapeNB  = static_cast<TH1*>(f->Get("ShapeRatio_nb")->Clone());

        shapeMETElec = static_cast<TH1*>(f->Get("ShapeRatioElec_met")->Clone());
        shapeMT2Elec = static_cast<TH1*>(f->Get("ShapeRatioElec_mt2")->Clone());
        shapeNTElec  = static_cast<TH1*>(f->Get("ShapeRatioElec_nt")->Clone());
        shapeNBElec  = static_cast<TH1*>(f->Get("ShapeRatioElec_nb")->Clone());

        f->Close();
        delete f;
    }
    else
    {
        std::cout << "Failed to open: syst_shape.root" << std::endl;
    }


    RegisterFunctions *rf = new RegisterFunctionsSyst;

    vector<Plotter::HistSummary> vh;

    Plotter::DatasetSummary dsDY_mm( "DY#rightarrow#mu#mu", fileMap["DYJetsToLL"],  "", "");
    Plotter::DatasetSummary dsDY_ee( "DY#rightarrowee",     fileMap["DYJetsToLL"],  "", "");
    Plotter::DatasetSummary dsZ_nunu("Z#rightarrow#nu#nu",  fileMap["ZJetsToNuNu"], "", "");

    Plotter::DataCollection dcDY_mm_cleamMET( "single", "cleanMetPt", {dsDY_mm});
    Plotter::DataCollection dcDY_ee_cleamMET( "single", "cleanMetPt", {dsDY_ee});
    Plotter::DataCollection dcZ_nunu_cleamMET("single", "cleanMetPt", {dsZ_nunu});

    vh.push_back(PHS("cleanMet",     {dcDY_mm_cleamMET, dcDY_ee_cleamMET, dcZ_nunu_cleamMET}, {1, 2}, "",            100, 0, 1500,  true,  true,  "MET [GeV]",  "Norm Events"));
    vh.push_back(PHS("cleanMet_cut", {dcDY_mm_cleamMET, dcDY_ee_cleamMET, dcZ_nunu_cleamMET}, {1, 2}, "passMETZinv", 100, 0, 1500,  true,  true,  "MET [GeV]",  "Norm Events"));

    TF1 *fit = new TF1("fit","pol1");
    fit->SetParameters(0.98, 0.001);
    Systematic test("systWgtTest", "cleanMetPt", fit);
    test.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(test, std::placeholders::_1));

    //Met shape syst 
    Systematic METSyst("systWgtMET", "cleanMetPt", shapeMET);
    METSyst.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(METSyst, std::placeholders::_1));

    //MT2 shape syst 
    Systematic MT2Syst("systWgtMT2", "best_had_brJet_MT2Zinv", shapeMT2);
    MT2Syst.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(MT2Syst, std::placeholders::_1));

    //NT shape uncertainty
    Systematic NTSyst("systWgtNT", "nTopCandSortedCntZinv", shapeNT);
    NTSyst.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(NTSyst, std::placeholders::_1));

    //NB shape uncertainty
    Systematic NBSyst("systWgtNB", "cntCSVSZinv", shapeNB);
    NBSyst.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(NBSyst, std::placeholders::_1));



    ////    //Met shape syst electrons!
    Systematic METSystElec("systWgtMETElec", "cleanMetPt", shapeMETElec);
    METSystElec.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(METSystElec, std::placeholders::_1));

    //MT2 shape syst 
    Systematic MT2SystElec("systWgtMT2Elec", "best_had_brJet_MT2Zinv", shapeMT2Elec);
    MT2SystElec.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(MT2SystElec, std::placeholders::_1));

    //NT shape uncertainty
    Systematic NTSystElec("systWgtNTElec", "nTopCandSortedCntZinv", shapeNTElec);
    NTSystElec.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(NTSystElec, std::placeholders::_1));

    //NB shape uncertainty
    Systematic NBSystElec("systWgtNBElec", "cntCSVSZinv", shapeNBElec);
    NBSystElec.bookHist(vh, fileMap["ZJetsToNuNu"]);
    static_cast<RegisterFunctionsSyst*>(rf)->addFunction(std::bind(NBSystElec, std::placeholders::_1));
    ////


    Plotter::DataCollection dcDY_nunu_nJetWgtSyst( "single", {{"njSystWeightedSB", dsZ_nunu}, {"njSystUnweightedSB", dsZ_nunu}});
    vh.emplace_back(PHS("systNJetWgtStat", {dcDY_nunu_nJetWgtSyst}, {2, 1}, "", 45, 0, 45, true, true, "Search Bin", "Events"));

    set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    Plotter plotter(vh, vvf, true, "systematics.root");
    plotter.setLumi(lumi);
    plotter.setPlotDir("Syst_plots");
    plotter.setDoHists(true);
    plotter.setDoTuple(false);
    plotter.setPrintInterval(10000);
    plotter.setRegisterFunction(rf);
    plotter.read();
    plotter.saveHists();
    plotter.plot();
}
