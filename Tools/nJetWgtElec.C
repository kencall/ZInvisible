#include "../../SusyAnaTools/Tools/NTupleReader.h"
#include "../../SusyAnaTools/Tools/samples.h"
#include "derivedTupleVariables.h"

#include <iostream>
#include <string>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"

int main()
{
    TH1::AddDirectory(false);

    TH1* njWTTbar_0b;
    TH1* njWDYZ_0b, *njWDYZ_0bElec;//Sam
    TH1* njWTTbar_g1b;
    TH1* njWDYZ_g1b, *njWDYZ_g1bElec;//Sam

    TH1* shapeMET;
    TH1* shapeMT2;
    TH1* shapeNT;
    TH1* shapeNB;
    
    TRandom3 tr3(153474);

    TFile *f = new TFile("dataMCweights.root");
    if(f)
    {
        njWTTbar_0b  = static_cast<TH1*>(f->Get("DataMC_nj_elmuZinv_0b_ht200_dphi")->Clone());
        njWDYZ_0b    = static_cast<TH1*>(f->Get("DataMC_nj_muZinv_0b_ht200_dphi")->Clone());
        njWDYZ_0bElec    = static_cast<TH1*>(f->Get("DataMC_nj_elZinv_0b_ht200_dphi")->Clone());//Sam

        njWTTbar_g1b = static_cast<TH1*>(f->Get("DataMC_nj_elmuZinv_g1b_ht200_dphi")->Clone());
        njWDYZ_g1b   = static_cast<TH1*>(f->Get("DataMC_nj_muZinv_g1b_ht200_dphi")->Clone());
        njWDYZ_g1bElec   = static_cast<TH1*>(f->Get("DataMC_nj_elZinv_g1b_ht200_dphi")->Clone());//Sam
        f->Close();
        delete f;
    }
    else
    {
        std::cout << "Failed to open: dataMCweights.root" << std::endl;
    }

    f = new TFile("syst_shape.root");
    if(f)
    {
        shapeMET = static_cast<TH1*>(f->Get("ShapeRatio_met")->Clone());
        shapeMT2 = static_cast<TH1*>(f->Get("ShapeRatio_mt2")->Clone());
        shapeNT  = static_cast<TH1*>(f->Get("ShapeRatio_nt")->Clone());
        shapeNB  = static_cast<TH1*>(f->Get("ShapeRatio_nb")->Clone());
        f->Close();
        delete f;
    }
    else
    {
        std::cout << "Failed to open: syst_shape.root" << std::endl;
    }

    AnaSamples::SampleSet        ss("/uscms/home/pastika/nobackup/zinv/dev/CMSSW_7_4_8/src/ZInvisible/Tools/condor/", AnaSamples::luminosity);
    AnaSamples::SampleCollection sc(ss);

    //f = TFile::Open("condor/minituple-Feb5.root");

    //Prep variables for nJet and shape systematic weights
    const int NTRIALS = 1000;
    const int NSEARCHBINS = 37;



    TH1 *h[5][NSEARCHBINS];
    std::vector<std::string> hnames = {"njet", "met", "mt2", "nt", "nb"};
    float N0[NSEARCHBINS];
    float N0square[NSEARCHBINS];
    float N[5][NSEARCHBINS][NTRIALS];

    TH1 *hElec[5][NSEARCHBINS];//Sam
    float N0Elec[NSEARCHBINS];//Sam
    float N0squareElec[NSEARCHBINS];//Sam
    float NElec[5][NSEARCHBINS][NTRIALS];//Sam




    float variations[5][20][NTRIALS];
    for(int iT = 0; iT < NTRIALS; ++iT)
    {
        for(int i = 0; i < 5; ++i) for(int j = 0; j < 20; ++j) variations[i][j][iT] = 1.0;
        for(int i = 0; i<=njWDYZ_g1b->GetNbinsX() + 1; ++i) if(njWDYZ_g1b->GetBinContent(i) > 1e-10) variations[0][i][iT] = (float)tr3.Gaus(1.0, njWDYZ_g1b->GetBinError(i)/njWDYZ_g1b->GetBinContent(i));

        for(int i = 0; i<=shapeMET  ->GetNbinsX() + 1; ++i) if(shapeMET->GetBinContent(i) > 1e-10)   variations[1][i][iT] = (float)tr3.Gaus(1.0, shapeMET->GetBinError(i)/shapeMET->GetBinContent(i));
        for(int i = 0; i<=shapeMT2  ->GetNbinsX() + 1; ++i) if(shapeMT2->GetBinContent(i) > 1e-10)   variations[2][i][iT] = (float)tr3.Gaus(1.0, shapeMT2->GetBinError(i)/shapeMT2->GetBinContent(i));
        for(int i = 0; i<=shapeNT   ->GetNbinsX() + 1; ++i) if(shapeNT->GetBinContent(i) > 1e-10)    variations[3][i][iT] = (float)tr3.Gaus(1.0, shapeNT->GetBinError(i)/shapeNT->GetBinContent(i));
        for(int i = 0; i<=shapeNB   ->GetNbinsX() + 1; ++i) if(shapeNB->GetBinContent(i) > 1e-10)    variations[4][i][iT] = (float)tr3.Gaus(1.0, shapeNB->GetBinError(i)/shapeNB->GetBinContent(i));

        for(int i = 0; i<=njWDYZ_g1bElec->GetNbinsX() + 1; ++i) if(njWDYZ_g1bElec->GetBinContent(i) > 1e-10) variations[0][i][iT] = (float)tr3.Gaus(1.0, njWDYZ_g1bElec->GetBinError(i)/njWDYZ_g1bElec->GetBinContent(i));//Sam
        for(int i = 0; i<=shapeMETElec  ->GetNbinsX() + 1; ++i) if(shapeMETElec->GetBinContent(i) > 1e-10)   variations[1][i][iT] = (float)tr3.Gaus(1.0, shapeMETElec->GetBinError(i)/shapeMETElec->GetBinContent(i));//Sam
        for(int i = 0; i<=shapeMT2Elec  ->GetNbinsX() + 1; ++i) if(shapeMT2Elec->GetBinContent(i) > 1e-10)   variations[2][i][iT] = (float)tr3.Gaus(1.0, shapeMT2Elec->GetBinError(i)/shapeMT2Elec->GetBinContent(i));//Sam
        for(int i = 0; i<=shapeNTElec   ->GetNbinsX() + 1; ++i) if(shapeNTElec->GetBinContent(i) > 1e-10)    variations[3][i][iT] = (float)tr3.Gaus(1.0, shapeNTElec->GetBinError(i)/shapeNTElec->GetBinContent(i));//Sam
        for(int i = 0; i<=shapeNBElec   ->GetNbinsX() + 1; ++i) if(shapeNBElec->GetBinContent(i) > 1e-10)    variations[4][i][iT] = (float)tr3.Gaus(1.0, shapeNBElec->GetBinError(i)/shapeNBElec->GetBinContent(i));//Sam
    }

    for(int ih = 0; ih < 5; ++ih)
    {
        for(int i = 0; i < NSEARCHBINS; ++i)
        {
            char name[128];
            sprintf(name, "hSB_%s_%d", hnames[ih].c_str(), i);
            h[ih][i] = new TH1D(name, name, 2000, -1, 3);

            for(int j = 0; j < NTRIALS; ++j)
            {
                N[ih][i][j] = 0.0;
                NElec[ih][i][j] = 0.0;
            }
        }
    }

    for(int i = 0; i < NSEARCHBINS; ++i)
    {
        N0square[i] = 0.0;
        N0[i] = 0.0;
        N0squareElec[i] = 0.0;
        N0Elec[i] = 0.0;
    }

    // //prep histograms for MET-MT2 corrolation study
    // TH1* hMET_nom  = new TH1D("hMET_nom",  "hMET_nom",  200, 0, 2000);
    // TH1* hMET_Gaus = new TH1D("hMET_Gaus", "hMET_Gaus", 200, 0, 2000);
    // TH1* hMET_Logi = new TH1D("hMET_Logi", "hMET_Logi", 200, 0, 2000);

    // TH1* hMT2_nom  = new TH1D("hMT2_nom",  "hMT2_nom",  200, 0, 2000);
    // TH1* hMT2_Gaus = new TH1D("hMT2_Gaus", "hMT2_Gaus", 200, 0, 2000);
    // TH1* hMT2_Logi = new TH1D("hMT2_Logi", "hMT2_Logi", 200, 0, 2000);

    // TH2* h2D_nom  = new TH2D("hMT2vMET_nom",  "h2D_nom ;MET;MT2", 200, 0, 2000, 200, 0, 2000);
    // TH2* h2D_Gaus = new TH2D("hMT2vMET_Gaus", "h2D_Gaus;MET;MT2", 200, 0, 2000, 200, 0, 2000);
    // TH2* h2D_Logi = new TH2D("hMT2vMET_Logi", "h2D_Logi;MET;MT2", 200, 0, 2000, 200, 0, 2000);

    AnaFunctions::prepareTopTagger();

    //for(auto& fs : sc["DYJetsToLL"])
    for(auto& fs : sc["ZJetsToNuNu"])
    {
        //size_t start = fs.filePath.rfind('/');
        //size_t stop  = fs.filePath.rfind('.');
        //std::string treeName = fs.filePath.substr(start + 1, stop - start - 1);

        TChain *t = new TChain(fs.treePath.c_str());
        fs.addFilesToChain(t);

        std::cout << "Tree: " << fs.treePath << std::endl;
        std::cout << "sigma*lumi: " << fs.getWeight() << std::endl;

        plotterFunctions::PrepareMiniTupleVars pmt(false);
        plotterFunctions::NJetWeight njWeight;
        plotterFunctions::TriggerInfo triggerInfo(true);
        plotterFunctions::MetSmear metSmear;

        NTupleReader tr(t);
        tr.registerFunction(pmt);
        tr.registerFunction(njWeight);
        tr.registerFunction(triggerInfo);
        //tr.registerFunction(metSmear);

        while(tr.getNextEvent())
        {
            if(tr.getEvtNum() % 10000 == 0) std::cout << "Event #: " << tr.getEvtNum() << std::endl;

            const int& nSearchBin = tr.getVar<int>("nSearchBin");
            const int& cntNJetsPt30Eta24Zinv = tr.getVar<int>("cntNJetsPt30Eta24Zinv");
            const double& cleanMetPt = tr.getVar<double>("cleanMetPt");
            const double& best_had_brJet_MT2Zinv = tr.getVar<double>("best_had_brJet_MT2Zinv");
            const int& cntCSVSZinv            = tr.getVar<int>("cntCSVSZinv");
            const int& nTopCandSortedCntZinv  = tr.getVar<int>("nTopCandSortedCntZinv");
            const bool& passBaseline = tr.getVar<bool>("passBaseline");
            const bool& passBaselineZinv = tr.getVar<bool>("passBaselineZinv");
            const bool& passLeptVeto = tr.getVar<bool>("passLeptVeto");
            const double& bTagSF_EventWeightSimple_Central = tr.getVar<double>("bTagSF_EventWeightSimple_Central");

            const double& triggerEffMC = tr.getVar<double>("TriggerEffMC");
            const double& nJetWgtDYZ   = tr.getVar<double>("nJetWgtDYZ");
            const double& normWgt0b    = tr.getVar<double>("normWgt0b");



            //fill stat uncertainty histograms here
            if(passBaselineZinv && passLeptVeto)
            {
                double weight = triggerEffMC * nJetWgtDYZ * normWgt0b * bTagSF_EventWeightSimple_Central * fs.getWeight();

                if(nSearchBin >= 0 && nSearchBin < NSEARCHBINS)
                {
                    N0[nSearchBin] += weight;
                    N0square[nSearchBin] += weight*weight;

                    N0Elec[nSearchBin] += weight;//Sam
                    N0squareElec[nSearchBin] += weight*weight;//Sam

                    for(int iTrial = 0; iTrial < NTRIALS; ++iTrial)
                    {
                        if(iTrial < NTRIALS)
                        {
                            N[0][nSearchBin][iTrial] += variations[0][njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv)][iTrial] * weight;
                            N[1][nSearchBin][iTrial] += variations[1][shapeMET->FindBin(cleanMetPt)][iTrial]              * weight;
                            N[2][nSearchBin][iTrial] += variations[2][shapeMT2->FindBin(best_had_brJet_MT2Zinv)][iTrial]  * weight;
                            N[3][nSearchBin][iTrial] += variations[3][shapeNT->FindBin(nTopCandSortedCntZinv)][iTrial]    * weight;
                            N[4][nSearchBin][iTrial] += variations[4][shapeNB->FindBin(cntCSVSZinv)][iTrial]              * weight;

                            NElec[0][nSearchBin][iTrial] += variations[0][njWDYZ_g1bElec->FindBin(cntNJetsPt30Eta24Zinv)][iTrial] * weight;//Sam
                            NElec[1][nSearchBin][iTrial] += variations[1][shapeMETElec->FindBin(cleanMetPt)][iTrial]              * weight;//Sam
                            NElec[2][nSearchBin][iTrial] += variations[2][shapeMT2Elec->FindBin(best_had_brJet_MT2Zinv)][iTrial]  * weight;//Sam
                            NElec[3][nSearchBin][iTrial] += variations[3][shapeNTElec->FindBin(nTopCandSortedCntZinv)][iTrial]    * weight;//Sam
                            NElec[4][nSearchBin][iTrial] += variations[4][shapeNBElec->FindBin(cntCSVSZinv)][iTrial]              * weight;//Sam
                        }
                    }
                }
            }
        }
    }
    
    for(int ih = 0; ih < 5; ++ih)
    {
        for(int i = 0; i < NSEARCHBINS; ++i)
        {
            for(int iTrial = 0; iTrial < NTRIALS; ++iTrial)
            {
                h[ih][i]->Fill(N[ih][i][iTrial]/N0[i]);
                hElec[ih][i]->Fill(NElec[ih][i][iTrial]/N0[i]);
            }
        }
    }

    TFile fout("syst_nJetWgt.root", "RECREATE");
    TH1 *syst68Max = new TH1D("syst68Max", "syst68Max", NSEARCHBINS, 0, NSEARCHBINS);
    TH1 *syst68MaxElec = new TH1D("syst68MaxElec", "syst68MaxElec", NSEARCHBINS, 0, NSEARCHBINS);
    for(int ih = 0; ih < 5; ++ih)
    {
        char name[128];
        sprintf(name, "syst68_%s", hnames[ih].c_str());
        TH1 *syst68 = new TH1D(name, name, NSEARCHBINS, 0, NSEARCHBINS);
        sprintf(name, "systRMS_%s", hnames[ih].c_str());
        TH1 *systRMS = new TH1D(name, name, NSEARCHBINS, 0, NSEARCHBINS);

        sprintf(name, "syst68Elec_%s", hnames[ih].c_str());//Sam
        TH1 *syst68Elec = new TH1D(name, name, NSEARCHBINS, 0, NSEARCHBINS);//Sam
        sprintf(name, "systRMSElec_%s", hnames[ih].c_str());//Sam
        TH1 *systRMSElec = new TH1D(name, name, NSEARCHBINS, 0, NSEARCHBINS);//Sam

        for(int i = 0; i < NSEARCHBINS; ++i) 
        {
            h[ih][i]->Write();
            h[ih][i]->Scale(1/h[ih][i]->Integral(0, h[ih][i]->GetNbinsX() + 1));
            hElec[ih][i]->Write();///////////////////////////////////////////////////Sam
            hElec[ih][i]->Scale(1/h[ih][i]->Integral(0, h[ih][i]->GetNbinsX() + 1));//Sam

            double ll = -999.9, ul = -999.9;
            double llElec = -999.9, ulElec = -999.9;//
            TH1* hint = (TH1*)h[ih][i]->Clone((std::string(h[ih][i]->GetName())+"_int").c_str());
            TH1* hintElec = (TH1*)hElec[ih][i]->Clone((std::string(hElec[ih][i]->GetName())+"_int").c_str());//Elec
            for(int iBin = 1; iBin <= h[ih][i]->GetNbinsX(); ++iBin)
            {
                hint->SetBinContent(iBin, h[ih][i]->Integral(0, iBin));
                hintElec->SetBinContent(iBin, hElec[ih][i]->Integral(0, iBin));//Sam
                if     (ll < 0 && h[ih][i]->Integral(0, iBin) > 0.16) ll = h[ih][i]->GetBinCenter(iBin);
                else if(ul < 0 && h[ih][i]->Integral(0, iBin) > 0.84) ul = h[ih][i]->GetBinCenter(iBin);
                if     (llElec < 0 && hElec[ih][i]->Integral(0, iBin) > 0.16) llElec = h[ih][i]->GetBinCenter(iBin);
                else if(ulElec < 0 && h[ih][i]->Integral(0, iBin) > 0.84) ulElec = h[ih][i]->GetBinCenter(iBin);
            }
            if(ll < 0) ll = 0;
            if(ul < 0) ul = 2.0;
            if(llElec < 0) llElec = 0;//Sam
            if(ulElec < 0) ulElec = 2.0;//Sam
            syst68->SetBinContent(i + 1, (ul - ll) / 2.0);
            syst68Max->SetBinContent(i + 1, std::max(syst68Max->GetBinContent(i + 1), (ul - ll) / 2.0));
            systRMS->SetBinContent(i + 1, h[ih][i]->GetRMS());

            syst68Elec->SetBinContent(i + 1, (ulElec - llElec) / 2.0);//Sam
            syst68MaxElec->SetBinContent(i + 1, std::max(syst68MaxElec->GetBinContent(i + 1), (ulElec - llElec) / 2.0));//Sam
            systRMSElec->SetBinContent(i + 1, hElec[ih][i]->GetRMS());//Sam
            //std::cout << "bin: " << i << "\tll: " << ll << "\tul: " << ul << "\tsym err: " << (ul - ll) / 2.0 << "\t:RMS: " << h[ih][i]->GetRMS() << std::endl;
            hint->Write();
            hintElec->Write();
        }
        syst68->Write();
        systRMS->Write();
        syst68Elec->Write();//Sam
        systRMSElec->Write();//Sam
    }
    syst68Max->Write();
    syst68MaxElec->Write();//Sam

    TH1 *centralvalue = new TH1D("centralvalue", "centralvalue", NSEARCHBINS, 0, NSEARCHBINS);
    TH1 *avgWgt = new TH1D("avgWgt", "avgWgt", NSEARCHBINS, 0, NSEARCHBINS);
    TH1 *neff = new TH1D("neff", "neff", NSEARCHBINS, 0, NSEARCHBINS);

    TH1 *centralvalueElec = new TH1D("centralvalue", "centralvalue", NSEARCHBINS, 0, NSEARCHBINS);//Sam
    TH1 *avgWgtElec = new TH1D("avgWgt", "avgWgt", NSEARCHBINS, 0, NSEARCHBINS);//Sam
    TH1 *neffElec = new TH1D("neff", "neff", NSEARCHBINS, 0, NSEARCHBINS);//Sam
    
    for(int i = 0; i < NSEARCHBINS; ++i)
    {
        centralvalue->SetBinContent(i + 1, N0[i]);
        if(N0[i] > 1e-10) avgWgt->SetBinContent(i + 1, N0square[i]/N0[i]);
        if(N0square[i] > 1e-10) neff->SetBinContent(i + 1, N0[i]*N0[i]/N0square[i]);

        centralvalueElec->SetBinContent(i + 1, N0Elec[i]);
        if(N0Elec[i] > 1e-10) avgWgt->SetBinContent(i + 1, N0squareElec[i]/N0Elec[i]);
        if(N0squareElec[i] > 1e-10) neffElec->SetBinContent(i + 1, N0Elec[i]*N0Elec[i]/N0squareElec[i]);
    }

    centralvalue->Write();
    avgWgt->Write();
    neff->Write();

    centralvalueElec->Write();
    avgWgtElec->Write();
    neffElec->Write();

    fout.Close();

    //Write out MET-MT2 histograms
    // TFile fout2("Corrolation_MET-MT2.root", "RECREATE");

    // hMET_nom ->Write();
    // hMET_Gaus->Write();
    // hMET_Logi->Write();

    // hMT2_nom ->Write();
    // hMT2_Gaus->Write();
    // hMT2_Logi->Write();

    // h2D_nom ->Write();
    // h2D_Gaus->Write();
    // h2D_Logi->Write();

    //fout2.Close();
}
