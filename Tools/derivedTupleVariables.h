#ifndef DERIVEDTUPLEVARIABLES_H
#define DERIVEDTUPLEVARIABLES_H

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "ScaleFactors.h"

#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TTModule.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include "TVector2.h"

#include <vector>
#include <iostream>
#include <string>
#include <set>

namespace plotterFunctions
{
    class GenerateWeight
    {
    private:
        TH1* muEff;
        TH1* muEffReco;
        TH2* muEffIso;
        TH2* elecEffReco;
        TH2* elecEffIso;
        TH2* muAcc;
        TH1* hZEff;
        TH1* hZAcc;
        TH1* hZAccElec;
        TH1* hZAccLepPt;
        TH1* hZAccLepPtElec;

        void generateWeight(NTupleReader& tr)
        {
            const std::vector<TLorentzVector>& jetsLVec         = tr.getVec<TLorentzVector>("jetsLVecLepCleaned");
            const std::vector<TLorentzVector>& cutMuVec         = tr.getVec<TLorentzVector>("cutMuVec");
            const std::vector<TLorentzVector>& cutElecVec       = tr.getVec<TLorentzVector>("cutElecVec");
            const std::vector<double>& cutMuActivity            = tr.getVec<double>("cutMuActivity");
            const std::vector<double>& cutElecActivity          = tr.getVec<double>("cutElecActivity");
            const std::vector<TLorentzVector*>& genMu           = tr.getVec<TLorentzVector*>("genMu");
            const std::vector<TLorentzVector*>& genMuInAcc      = tr.getVec<TLorentzVector*>("genMuInAcc");
            const std::vector<TLorentzVector*>& genMatchMuInAcc = tr.getVec<TLorentzVector*>("genMatchMuInAcc");

            const int& pdgIdZDec      = tr.getVar<int>("pdgIdZDec");
            const bool& passMuZinvSel = tr.getVar<bool>("passMuZinvSel");
            const bool& passElecZinvSel = tr.getVar<bool>("passElecZinvSel");
            const double& ht          = tr.getVar<double>("ht");
            const double& bestRecoZPt = tr.getVar<double>("bestRecoZPt");
            const double& genZPt      = tr.getVar<double>("genZPt");

	    const double& stored_weight = tr.getVar<double>("stored_weight");

            // Calculate PU weight

            // Calculate Z-eff weight

            const double zMassMin = 71.0;
            const double zMass    = 91.0;
            const double zMassMax = 111.0;

            double zMassCurrent = 1.0e300, zEff = 1.0e100, zAcc = 1.0e100;
            for(int i = 0; i < cutMuVec.size(); ++i)
            {
                if(cutMuVec[i].Pt() < 10) continue;
                for(int j = 0; j < i && j < cutMuVec.size(); ++j)
                {
                    if(cutMuVec[j].Pt() < 10) continue;
                    double zm = (cutMuVec[i] + cutMuVec[j]).M();
                    if(fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    {
                        zMassCurrent = zm;
                        double mu1pt = cutMuVec[i].Pt();
                        double mu2pt = cutMuVec[j].Pt();

                        double mu1Act = cutMuActivity[i];
                        double mu2Act = cutMuActivity[j];

                        //set to not overflow histograms
                        if(mu1pt >= 2000.0) mu1pt = 1999.9;
                        if(mu2pt >= 2000.0) mu2pt = 1999.9;

                        //Get mu efficiencies
                        double muEff1 = 0.0, muEff2 = 0.0;

                        if(muEff && muEffReco && muEffIso)
                        {
                            //Fit to reco eff (eff = p0 + p1*pt + p2*pt^2
                            //PHYS14
                            //const double fitStart = 200.0; // extended to 1400 GeV
                            //const double p0 =     0.955847; // +/- 0.461944
                            //const double p1 = -2.24431e-05; // +/- 0.00128305
                            //const double p2 = -5.68907e-08; // +/- 7.85913e-07
                            //Sprint15
                            //const double fitStart = 200.0; // extended to 1400 GeV
                            //const double p0 =  9.83467e-01; // +/- 1.54469e+00
                            //const double p1 = -7.81897e-06; // +/- 4.16487e-03
                            //const double p2 = -1.22092e-08; // +/- 2.18556e-06
                            //Spring15 low stats
                            //const double fitStart = 200.0;  // extended to 2000 GeV
                            //const double p0 =     0.979238; //   +/-   1.17313
                            //const double p1 = -6.47338e-06; //   +/-   0.00342715
                            //const double p2 = -1.16258e-08; //   +/-   1.87837e-06
                            //Spring15 extended samples
                            const double fitStart = 200.0;  // extended to 2000 GeV
                            const double p0 =      0.97431; //   +/-   1.17119     
                            const double p1 =  2.87484e-05; //   +/-   0.00341371  
                            const double p2 = -5.37058e-08; //   +/-   1.86309e-06 

                            double muRecoEff = 0.0;

                            int recoPtBin = muEffReco->GetXaxis()->FindBin(mu1pt);
                            if(mu1pt > fitStart) muRecoEff = p0 + p1*mu1pt + p2*mu1pt*mu1pt;
                            else                 muRecoEff = muEffReco->GetBinContent(recoPtBin);
                            int isoPtBin = muEffIso->GetXaxis()->FindBin(mu1pt);
                            if(isoPtBin >= muEffIso->GetNbinsX()) isoPtBin = muEffIso->GetNbinsX();
                            int isoActBin = muEffIso->GetYaxis()->FindBin(mu1Act);
                            if(isoActBin >= muEffIso->GetNbinsY()) isoActBin = muEffIso->GetNbinsY();
                            muEff1 = muRecoEff * muEffIso->GetBinContent(isoPtBin, isoActBin);

                            muRecoEff = 0.0;
                            recoPtBin = muEffReco->GetXaxis()->FindBin(mu2pt);
                            if(mu2pt > fitStart) muRecoEff = p0 + p1*mu2pt + p2*mu2pt*mu2pt;
                            else                 muRecoEff = muEffReco->GetBinContent(recoPtBin);
                            isoPtBin = muEffIso->GetXaxis()->FindBin(mu2pt);
                            if(isoPtBin >= muEffIso->GetNbinsX()) isoPtBin = muEffIso->GetNbinsX();
                            isoActBin = muEffIso->GetYaxis()->FindBin(mu2Act);
                            if(isoActBin >= muEffIso->GetNbinsY()) isoActBin = muEffIso->GetNbinsY();
                            muEff2 = muRecoEff * muEffIso->GetBinContent(isoPtBin, isoActBin);
                        }


                        if((mu1pt > 20 && muEff1 < 1.0e-5) || (mu2pt > 20 && muEff2 < 1.0e-5))
                        {
                            std::cout << "SMALL muEff!!! muEff1: " << muEff1 << "\tmuEff2: " << muEff2 << "\t" << mu1pt << "\t" << mu2pt << std::endl;
                            zEff = 1.0e-10;
                        }
                        else zEff = muEff1 * muEff2;
                    }
                }
            }

            zMassCurrent = 1.0e300;
            double zEffElec = 1.0e100, zAccElec = 1.0e100;
            for(int i = 0; i < cutElecVec.size(); ++i)
            {
                if(cutElecVec[i].Pt() < 10) continue;
                for(int j = 0; j < i && j < cutElecVec.size(); ++j)
                {
                    if(cutElecVec[j].Pt() < 10) continue;
                    double zm = (cutElecVec[i] + cutElecVec[j]).M();
                    if(fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    {
                        zMassCurrent = zm;
                        double elec1pt = cutElecVec[i].Pt();
                        double elec2pt = cutElecVec[j].Pt();

                        double elec1Act = cutElecActivity[i];
                        double elec2Act = cutElecActivity[j];

                        //set to not overflow histograms
                        if(elec1pt >= 2000.0) elec1pt = 1999.9;
                        if(elec2pt >= 2000.0) elec2pt = 1999.9;

                        //Get elec efficiencies
                        double elecEff1 = 0.0, elecEff2 = 0.0;
                        if(elecEffReco && elecEffIso)
                        {
                            //Fit to iso eff (eff = p0 + p1*pt + p2*pt^2
                            const double fitStart =200.0;  // extended to 2000 GeV
                            const double p0 =     0.989411; //   +/-   1.18467
                            const double p1 =  3.66321e-06; //   +/-   0.00346729
                            const double p2 = -5.68292e-09; //   +/-   1.90334e-06

                            double elecIsoEff = 0.0;

                            //int isoPtBin = elecEffIso->GetXaxis()->FindBin(elec1pt);
                            //if(elec1pt > fitStart) elecIsoEff = p0 + p1*elec1pt + p2*elec1pt*elec1pt;
                            //else                   elecIsoEff = elecEffIso->GetBinContent(isoPtBin);

                            int recoPtBin = elecEffReco->GetXaxis()->FindBin(elec1pt);
                            if(recoPtBin >= elecEffReco->GetNbinsX()) recoPtBin = elecEffReco->GetNbinsX();

                            int recoActBin = elecEffReco->GetYaxis()->FindBin(elec1Act);
                            if(recoActBin >= elecEffReco->GetNbinsY()) recoActBin = elecEffReco->GetNbinsY();

                            elecIsoEff = elecEffIso->GetBinContent(recoPtBin, recoActBin);
                            elecEff1 = elecIsoEff * elecEffReco->GetBinContent(recoPtBin, recoActBin);
                            //std::cout << elec1pt << "\t" << recoPtBin << "\t" <<  elec1Act << "\t" << recoActBin << "\t" << elecEff1 << "\t" << elecIsoEff << "\t" << elecEffReco->GetBinContent(recoPtBin, recoActBin) << std::endl;

                            elecIsoEff = 0.0;
                            //isoPtBin = elecEffIso->GetXaxis()->FindBin(elec2pt);
                            //if(elec2pt > fitStart) elecIsoEff = p0 + p1*elec2pt + p2*elec2pt*elec2pt;
                            //else                   elecIsoEff = elecEffIso->GetBinContent(isoPtBin);

                            recoPtBin = elecEffReco->GetXaxis()->FindBin(elec2pt);
                            if(recoPtBin >= elecEffReco->GetNbinsX()) recoPtBin = elecEffReco->GetNbinsX();

                            recoActBin = elecEffReco->GetYaxis()->FindBin(elec2Act);
                            if(recoActBin >= elecEffReco->GetNbinsY()) recoActBin = elecEffReco->GetNbinsY();

                            elecIsoEff = elecEffIso->GetBinContent(recoPtBin, recoActBin);
                            elecEff2 = elecIsoEff * elecEffReco->GetBinContent(recoPtBin, recoActBin);
                            //std::cout << elec2pt << "\t" << recoPtBin << "\t" << elec2Act << "\t" << recoActBin << "\t" << elecEff2 << "\t" << elecIsoEff << "\t" << elecEffReco->GetBinContent(recoPtBin, recoActBin) << std::endl;
                        }


                        if((elec1pt > 33 && elec2pt) && (elecEff1 < 1.0e-5 || elecEff2 < 1.0e-5))
                        {
                            std::cout << "SMALL elecEff!!! elecEff1: " << elecEff1 << "\telecEff2: " << elecEff2 << "\t" << elec1pt << "\t" << elec2pt << std::endl;
                            zEffElec = 1.0e-10;
                        }
                        else zEffElec = elecEff1 * elecEff2;
                    }
                }
            }

            double genCleanHt = ht;
            for(auto& tlvp : genMuInAcc) if(tlvp->Pt() > 50) genCleanHt -= tlvp->Pt();

            const double MHT_jetPtMin = 30.0;
            TLorentzVector MHT;
            double mu1dRMin = 99.9, mu2dRMin = 99.9;
            for(auto& jet : jetsLVec)
            {
                double mu1dR = 999.9, mu2dR = 999.9;
                if(cutMuVec.size() >= 1) mu1dR = ROOT::Math::VectorUtil::DeltaR(jet, cutMuVec[0]);
                if(cutMuVec.size() >= 2) mu2dR = ROOT::Math::VectorUtil::DeltaR(jet, cutMuVec[1]);
                mu1dRMin = std::min(mu1dRMin, mu1dR);
                mu2dRMin = std::min(mu2dRMin, mu2dR);

                if(jet.Pt() > MHT_jetPtMin) MHT += jet;
            }

            double mudR = 99.9, genMudR = 99.9, genMudPhi = 99.9, genMudEta = 99.9;

            if(cutMuVec.size() >= 2) mudR    = ROOT::Math::VectorUtil::DeltaR(cutMuVec[0], cutMuVec[1]);
            if(genMu.size() > 2) std::cout << "MORE THAN 2 GEN MUONS: " << genMu.size() << std::endl;
            if(genMu.size() >= 2)
            {
                genMudR = ROOT::Math::VectorUtil::DeltaR(*(genMu[0]), *(genMu[1]));
                genMudPhi = ROOT::Math::VectorUtil::DeltaPhi(*(genMu[0]), *(genMu[1]));
                genMudEta = genMu[0]->Eta() - genMu[1]->Eta();
            }

            // muon Z pt acceptance corrections
            //functional form [2] - exp([0] + [1]*x)
            //PHYS14
            //double acc_p0 = -2.91374e-01;
            //double acc_p1 = -4.42884e-03;
            //double acc_p2 =  9.51190e-01;
            //Spring15 low stats
            //double acc_p0 = -2.64921e-01;
            //double acc_p1 = -4.65305e-03;
            //double acc_p2 =  9.50493e-01;
            //Spring15 extended sample
            const double acc_p0 = -2.60663e-01;
            const double acc_p1 = -4.60623e-03;
            const double acc_p2 = 9.49434e-01;

            if(hZAcc && hZAccLepPt)
            {
                //Eta portion of acceptance
                if(bestRecoZPt < 150) zAcc = hZAcc->GetBinContent(hZAcc->GetXaxis()->FindBin(bestRecoZPt));
                else                  zAcc = acc_p2 - exp(acc_p0 + acc_p1 * bestRecoZPt);

                //pt portion of acceptance
                zAcc *= hZAccLepPt->GetBinContent(hZAccLepPt->GetXaxis()->FindBin(bestRecoZPt));
            }

            if(passMuZinvSel && zAcc < 0.05)
            {
                std::cout << "WARNING: Muon Z acceptance < 0.05, forcing weight to zero! Acc_Z: " << zAcc << std::endl;
                zAcc = 1.0e101;
            }

            if(passMuZinvSel && zEff < 0.04)
            {
                std::cout << "WARNING: Muon Z efficiency < 0.05, forcing weight to zero! Eff_Z: " << zEff << "\tZ(Pt): " << bestRecoZPt <<  std::endl;
                zEff = 1.0e101;
            }

            // electron Z pt acceptance corrections
            //functional form [2] - exp([0] + [1]*x)
            //Spring15
            //double elecAcc_p0 = -3.18955e-02;
            //double elecAcc_p1 = -5.16522e-03;
            //double elecAcc_p2 = 8.92280e-01;
            double elecAcc_p0 = -2.50679e-01;
            double elecAcc_p1 = -5.34976e-03;
            double elecAcc_p2 = 9.33562e-01;

            if(hZAccElec && hZAccLepPtElec)
            {
                //Eta portion of the acceptance
                if(bestRecoZPt < 100) zAccElec = hZAccElec->GetBinContent(hZAccElec->GetXaxis()->FindBin(bestRecoZPt));
                else                  zAccElec = elecAcc_p2 - exp(elecAcc_p0 + elecAcc_p1 * bestRecoZPt);

                //Pt portion of the acceptance
                zAccElec = hZAccLepPtElec->GetBinContent(hZAccLepPtElec->GetXaxis()->FindBin(bestRecoZPt));
            }

            if(passElecZinvSel && zAccElec < 0.05)
            {
                std::cout << "WARNING: Elec Z acceptance < 0.05, forcing weight to zero! Acc_Z: " << zAccElec << std::endl;
                zAccElec = 1.0e101;
            }

            if(passElecZinvSel && zEffElec < 0.04)
            {
                std::cout << "WARNING: Elec Z efficiency < 0.05, forcing weight to zero! Eff_Z: " << zEffElec << "\tZ(Pt): " << bestRecoZPt <<  std::endl;
                zEffElec = 1.0e101;
            }

	    // Process the generator weight
	    double genWeight = 1.;
	    // Never apply this weight for data! In the old ntuple version <=3 this is "-1", in the newer ones it is "0"
	    if(stored_weight < 0)
	      genWeight = -1.;

            tr.registerDerivedVar("mu1dRMin", mu1dRMin);
            tr.registerDerivedVar("mu2dRMin", mu2dRMin);
            tr.registerDerivedVar("mudR", mudR);
            tr.registerDerivedVar("genMudR", genMudR);
            tr.registerDerivedVar("genMudPhi", genMudPhi);
            tr.registerDerivedVar("genMudEta", genMudEta);

            tr.registerDerivedVar("zEff", zEff);
            tr.registerDerivedVar("zEffWgt", 1.0/zEff);
            tr.registerDerivedVar("zAcc", zAcc);
            tr.registerDerivedVar("zAccWgt", 1.0/zAcc);

            tr.registerDerivedVar("zEffElec", zEffElec);
            tr.registerDerivedVar("zEffWgtElec", 1.0/zEffElec);
            tr.registerDerivedVar("zAccElec", zAccElec);
            tr.registerDerivedVar("zAccWgtElec", 1.0/zAccElec);

            tr.registerDerivedVar("cleanMHt", MHT.Pt());
            tr.registerDerivedVar("cleanMHtPhi", MHT.Phi());

	    tr.registerDerivedVar("genWeight", genWeight);
        }

    public:
        GenerateWeight()
        {
            muEff          = nullptr;
            muEffReco      = nullptr;
            muEffIso       = nullptr;
            muAcc          = nullptr;
            hZEff          = nullptr;
            hZAcc          = nullptr;
            elecEffReco    = nullptr;
            elecEffIso     = nullptr;
            hZAccElec      = nullptr;
            hZAccLepPt     = nullptr;
            hZAccLepPtElec = nullptr;

            TH1::AddDirectory(false);
            TFile *f = new TFile("lepEffHists.root");
            if(f)
            {
                muEff          = static_cast<TH1*>(f->Get("hMuEffPt_ratio"));
                muEffReco      = static_cast<TH1*>(f->Get("hMuEffPtReco_ratio"));
                muEffIso       = static_cast<TH2*>(f->Get("hMuEffPtActIso_ratio"));
                muAcc          = static_cast<TH2*>(f->Get("hMuAcc_ratio"));
                hZEff          = static_cast<TH1*>(f->Get("hZEffPt_ratio"));
                hZAcc          = static_cast<TH1*>(f->Get("hZAccPtSmear_ratio"));
                elecEffReco    = static_cast<TH2*>(f->Get("hElecEffPtActReco_ratio"));
                elecEffIso     = static_cast<TH2*>(f->Get("hElecEffPtActIso_ratio"));
                hZAccElec      = static_cast<TH1*>(f->Get("hZElecAccPtSmear_ratio"));
                hZAccLepPt     = static_cast<TH1*>(f->Get("hZAccPtMuPtSmear_ratio"));
                hZAccLepPtElec = static_cast<TH1*>(f->Get("hZElecAccPtPtSmear_ratio"));
                f->Close();
                delete f;
            }
            else
            {
                std::cout << "Failed to open: lepEffHists.root" << std::endl;
            }
        }

        ~GenerateWeight()
        {
            //if(muEff)          delete muEff;
            //if(muEffReco)      delete muEffReco;
            //if(muEffIso)       delete muEffIso;
            //if(muAcc)          delete muAcc;
            //if(hZEff)          delete hZEff;
            //if(hZAcc)          delete hZAcc;
            //if(elecEffReco)    delete elecEffReco;
            //if(elecEffIso)     delete elecEffIso;
            //if(hZAccElec)      delete hZAccElec;
            //if(hZAccLepPt)     delete hZAccLepPt;
            //if(hZAccLepPtElec) delete hZAccLepPtElec;
        }

        void operator()(NTupleReader& tr)
        {
            generateWeight(tr);
        }
    };

    class NJetWeight
    {
    private:
        TH1* njWTTbar_0b;
        TH1* njWDYZ_0b;
        TH1* njWTTbar_g1b;
        TH1* njWDYZ_g1b;

        TH1* MCfake1b;
        TH1* MCfake2b;
        TH1* MCfake3b;

        void generateWeight(NTupleReader& tr)
        {
            const int& cntNJetsPt30Eta24Zinv = tr.getVar<int>("cntNJetsPt30Eta24Zinv");
            const int& cntCSVSZinv = tr.getVar<int>("cntCSVSZinv");

            double wTT = 1.0;
            double wDY = 1.0;

	    if(cntCSVSZinv == 0)
	    {
		if(njWTTbar_0b)  wTT = njWTTbar_0b->GetBinContent(njWTTbar_0b->FindBin(cntNJetsPt30Eta24Zinv));
		if(njWDYZ_0b)    wDY = njWDYZ_0b->GetBinContent(njWDYZ_0b->FindBin(cntNJetsPt30Eta24Zinv));
	    }
	    else
	    {
		if(njWTTbar_g1b) wTT = njWTTbar_g1b->GetBinContent(njWTTbar_g1b->FindBin(cntNJetsPt30Eta24Zinv));
		if(njWDYZ_g1b)   wDY = njWDYZ_g1b->GetBinContent(njWDYZ_g1b->FindBin(cntNJetsPt30Eta24Zinv));
	    }

            double nJet1bfakeWgt = 1.0;
            double nJet2bfakeWgt = 1.0;
            double nJet3bfakeWgt = 1.0;

            if(MCfake1b)   nJet1bfakeWgt = MCfake1b->GetBinContent(MCfake1b->FindBin(cntNJetsPt30Eta24Zinv));
            if(MCfake2b)   nJet2bfakeWgt = MCfake2b->GetBinContent(MCfake2b->FindBin(cntNJetsPt30Eta24Zinv));
            if(MCfake3b)   nJet3bfakeWgt = MCfake3b->GetBinContent(MCfake3b->FindBin(cntNJetsPt30Eta24Zinv));

	    double normWgt0b = ScaleFactors::sf_norm0b();

            tr.registerDerivedVar("nJetWgtTTbar", wTT);
            tr.registerDerivedVar("nJetWgtDYZ",   wDY);

            tr.registerDerivedVar("nJet1bfakeWgt", nJet1bfakeWgt);
            tr.registerDerivedVar("nJet2bfakeWgt", nJet2bfakeWgt);
            tr.registerDerivedVar("nJet3bfakeWgt", nJet3bfakeWgt);

            tr.registerDerivedVar("normWgt0b", normWgt0b);
        }

    public:
        NJetWeight()
        {
            TH1::AddDirectory(false);

            njWTTbar_0b  = nullptr;
            njWDYZ_0b    = nullptr;
            njWTTbar_g1b = nullptr;
            njWDYZ_g1b   = nullptr;
            MCfake1b = nullptr;
            MCfake2b = nullptr;
            MCfake3b = nullptr;

            TFile *f = new TFile("njetWgtHists.root");
            if(f)
            {
                MCfake1b = static_cast<TH1*>(f->Get("h_njRatio_1fake"));
                MCfake2b = static_cast<TH1*>(f->Get("h_njRatio_2fake"));
                MCfake3b = static_cast<TH1*>(f->Get("h_njRatio_3fake"));
                f->Close();
                delete f;
            }
            else
            {
                std::cout << "Failed to open: njetWgtHists.root" << std::endl;
            }

            f = new TFile("dataMCweights.root");
            if(f)
            {
                njWTTbar_0b  = static_cast<TH1*>(f->Get("DataMC_nj_elmuZinv_ht200_dphi"));
                njWDYZ_0b    = static_cast<TH1*>(f->Get("DataMC_nj_muZinv_0b_ht200_dphi"));
                njWTTbar_g1b = static_cast<TH1*>(f->Get("DataMC_nj_elmuZinv_ht200_dphi"));
                njWDYZ_g1b   = static_cast<TH1*>(f->Get("DataMC_nj_muZinv_g1b_ht200_dphi"));
                f->Close();
                delete f;
            }
            else
            {
                std::cout << "Failed to open: dataMCweights.root" << std::endl;
            }
        }

        ~NJetWeight()
        {
            //if(njWTTbar) delete njWTTbar;
            //if(njWDYZ)   delete njWDYZ;
        }

        void operator()(NTupleReader& tr)
        {
            generateWeight(tr);
        }
    };

    class LepInfo
    {
    private:
        TRandom3 *tr3;

        void lepInfo(NTupleReader& tr)
        {
            const std::vector<int>& genDecayPdgIdVec        = tr.getVec<int>("genDecayPdgIdVec");
            const std::vector<TLorentzVector>& genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
            const std::vector<TLorentzVector>& muonsLVec    = tr.getVec<TLorentzVector>("muonsLVec");
            const std::vector<double>& muonsRelIso          = tr.getVec<double>("muonsRelIso");
            const std::vector<double>& muonsMiniIso         = tr.getVec<double>("muonsMiniIso");
            const std::vector<int>& W_emuVec                = tr.getVec<int>("W_emuVec");
            const std::vector<double>& muonsCharge          = tr.getVec<double>("muonsCharge");
            const std::vector<TLorentzVector>& jetsLVec     = tr.getVec<TLorentzVector>("jetsLVec");
            const std::vector<double>& recoJetschargedEmEnergyFraction     = tr.getVec<double>("recoJetschargedEmEnergyFraction");
            const std::vector<double>& recoJetschargedHadronEnergyFraction = tr.getVec<double>("recoJetschargedHadronEnergyFraction");
            const std::vector<int> & muonsFlagIDVec = tr.getVec<int>("muonsFlagMedium");
            const std::vector<int>&  elesFlagIDVec  = tr.getVec<int>("elesFlagVeto");

            const std::vector<double>& muonspfActivity      = tr.getVec<double>("muonspfActivity");
            const std::vector<double>& elespfActivity       = tr.getVec<double>("elespfActivity");
            const std::vector<double>& W_emu_pfActivityVec  = tr.getVec<double>("W_emu_pfActivityVec");

            //const double& ht                             = tr.getVar<double>("ht");
            const double& met                            = tr.getVar<double>("met");
            const double& metphi                         = tr.getVar<double>("metphi");

            const double& calomet                        = tr.getVar<double>("calomet");

            const std::vector<TLorentzVector, std::allocator<TLorentzVector> > elesLVec = tr.getVec<TLorentzVector>("elesLVec");
            const std::vector<double>& elesMiniIso          = tr.getVec<double>("elesMiniIso");
            const std::vector<double>& elesCharge           = tr.getVec<double>("elesCharge");
            const std::vector<unsigned int>& elesisEB       = tr.getVec<unsigned int>("elesisEB");

            int nElectrons = tr.getVar<int>("nElectrons_CUT");
            int nMuons = tr.getVar<int>("nMuons_CUT");

            bool passMuonVeto = false;
            bool passEleVeto = false;

            try
            {
                const bool& passMuonVetoTmp  = tr.getVar<bool>("passMuonVeto");
                const bool& passEleVetoTmp   = tr.getVar<bool>("passEleVeto");
                if(&passMuonVetoTmp != nullptr) passMuonVeto = passMuonVetoTmp;
                if(&passEleVetoTmp != nullptr) passEleVeto = passEleVetoTmp;
            }
            catch(const std::string e)
            {
                //std::cout << "void muInfo(NTupleReader& tr): Caught exception, variable \"" << e << "\" not found" << std::endl;
            }

            std::vector<const TLorentzVector*>* genMatchIsoElecInAcc = new std::vector<const TLorentzVector*>();
            std::vector<const TLorentzVector*>* genMatchElecInAcc = new std::vector<const TLorentzVector*>();
            std::vector<double>* genMatchElecInAccRes = new std::vector<double>();
            std::vector<const TLorentzVector*>* genElecInAcc = new std::vector<const TLorentzVector*>();
            std::vector<const TLorentzVector*>* genElec = new std::vector<const TLorentzVector*>();
            std::vector<double>* genMatchIsoElecInAccAct = new std::vector<double>();
            std::vector<double>* genMatchElecInAccAct = new std::vector<double>();
            std::vector<double>* genElecInAccAct = new std::vector<double>();
            std::vector<double>* genElecAct = new std::vector<double>();

            std::vector<TLorentzVector>* cutElecVec = new std::vector<TLorentzVector>();
            std::vector<double>* cutElecCharge = new std::vector<double>();
            std::vector<double>* cutElecActivity = new std::vector<double>();

            std::vector<TLorentzVector> cutElecVecRecoOnly;

            std::vector<const TLorentzVector*>* genMatchIsoMuInAcc = new std::vector<const TLorentzVector*>();
            std::vector<const TLorentzVector*>* genMatchMuInAcc = new std::vector<const TLorentzVector*>();
            std::vector<double>* genMatchMuInAccRes = new std::vector<double>();
            std::vector<const TLorentzVector*>* genMuInAcc = new std::vector<const TLorentzVector*>();
            std::vector<const TLorentzVector*>* genMu = new std::vector<const TLorentzVector*>();
            std::vector<double>* genMatchIsoMuInAccAct = new std::vector<double>();
            std::vector<double>* genMatchMuInAccAct = new std::vector<double>();
            std::vector<double>* genMuInAccAct = new std::vector<double>();
            std::vector<double>* genMuAct = new std::vector<double>();
            std::vector<TLorentzVector>* cutMuVec = new std::vector<TLorentzVector>();
            std::vector<double>* cutMuCharge = new std::vector<double>();
            std::vector<double>* cutMuActivity = new std::vector<double>();

            std::vector<TLorentzVector> cutMuVecRecoOnly;

            //muon selections
            int sumMuCharge = 0;
            int nTriggerMuons = 0;

            for(int i = 0; i < muonsLVec.size(); ++i)
            {
                if(AnaFunctions::passMuon( muonsLVec[i], 0.0, 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr)) // emulates muons with pt but no iso requirements.
                {
                    cutMuVecRecoOnly.push_back(muonsLVec[i]);
                }
                if(AnaFunctions::passMuon( muonsLVec[i], muonsMiniIso[i], 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr))
                {
                    if(AnaFunctions::passMuon( muonsLVec[i], muonsRelIso[i], 0.0, muonsFlagIDVec[i], AnaConsts::muonsMiniIsoArr))
                    {
                        if(nTriggerMuons == 0 && muonsLVec[i].Pt() > 17)  nTriggerMuons++;
                        else if(muonsLVec[i].Pt() > 8)  nTriggerMuons++;
                    }
                    cutMuVec->push_back(muonsLVec[i]);
                    cutMuCharge->push_back(muonsCharge[i]);
                    cutMuActivity->push_back(muonspfActivity[i]);
                    if(muonsCharge[i] > 0) sumMuCharge++;
                    else                   sumMuCharge--;
                }
            }

            //electron selection
            int sumElecCharge = 0;
            for(int i = 0; i < elesLVec.size(); ++i)
            {
                if(AnaFunctions::passElectron(elesLVec[i], 0.0, -1, elesisEB[i], elesFlagIDVec[i], AnaConsts::elesMiniIsoArr)) // emulates muons with pt but no iso requirements.
                {
                    cutElecVecRecoOnly.push_back(elesLVec[i]);
                }

                if(AnaFunctions::passElectron(elesLVec[i], elesMiniIso[i], -1, elesisEB[i], elesFlagIDVec[i], AnaConsts::elesMiniIsoArr))
                {
                    cutElecVec->push_back(elesLVec[i]);
                    cutElecCharge->push_back(elesCharge[i]);
                    cutElecActivity->push_back(AnaFunctions::getElectronActivity(elesLVec[i], jetsLVec, recoJetschargedHadronEnergyFraction, AnaConsts::elesAct));
                    if(elesCharge[i] > 0) sumElecCharge++;
                    else                  sumElecCharge--;
                }
            }


            //mu45 non-iso trigger emulation
            const double effsnom2012ABC[] = {0.928,0.8302,0.8018};
            const double upedge2012ABC[] = { 0.9, 1.2, 2.1};
            bool muTrigMu45 = false;
            for(TLorentzVector& mu : *cutMuVec)
            {
                if(mu.Pt() > 45)
                {
                    for(int iBin = 0; iBin < sizeof(effsnom2012ABC)/sizeof(double); ++iBin)
                    {
                        if(mu.Eta() < upedge2012ABC[iBin] && tr3 && tr3->Uniform() < effsnom2012ABC[iBin])
                        {
                            muTrigMu45 = true;
                            break;
                        }
                    }
                    if(muTrigMu45) break;
                }
            }

            double genHt = 0.0;

            const double   minMuPt = 20.0,   highMuPt = 45.0;
            const double minElecPt = 33.0, highElecPt = 33.0;
            double nuPt1 = -999.9, nuPt2 = -999.9;

            //Gen info parsing
            if(&genDecayPdgIdVec != nullptr && &genDecayLVec != nullptr)
            {
                for(int i = 0; i < genDecayPdgIdVec.size() && i < genDecayLVec.size(); ++i)
                {
                    if((abs(genDecayPdgIdVec[i]) != 0 &&  abs(genDecayPdgIdVec[i]) < 6) || (abs(genDecayPdgIdVec[i]) > 100 && abs(genDecayPdgIdVec[i]) < 10000)) genHt += genDecayLVec[i].Pt();

                    if(genDecayPdgIdVec[i] ==  13) nuPt1 = genDecayLVec[i].Pt();
                    if(genDecayPdgIdVec[i] == -13) nuPt2 = genDecayLVec[i].Pt();
                }

                for(int index = 0; index < W_emuVec.size(); ++index)
                {
                    int i = W_emuVec[index];

                    //muon efficiency and acceptance
                    if(abs(genDecayPdgIdVec[i]) == 13)
                    {
                        genMu->push_back(&genDecayLVec[i]);
                        genMuAct->push_back(W_emu_pfActivityVec[index]);
                        if(AnaFunctions::passMuonAccOnly(genDecayLVec[i], AnaConsts::muonsMiniIsoArr) && genDecayLVec[i].Pt() > minMuPt)
                        {
                            genMuInAcc->push_back(&genDecayLVec[i]);
                            genMuInAccAct->push_back(genMuAct->back());
                            double dRMin = 999.9;
                            double matchPt = -999.9;
                            for(int j = 0; j < cutMuVecRecoOnly.size(); ++j)
                            {
                                double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], cutMuVecRecoOnly[j]);
                                if(dR < dRMin)
                                {
                                    dRMin = dR;
                                    matchPt = cutMuVecRecoOnly[j].Pt();
                                }
                            }
                            if(dRMin < 0.02)
                            {
                                genMatchMuInAcc->push_back(&genDecayLVec[i]);
                                genMatchMuInAccAct->push_back(genMuAct->back());
                                genMatchMuInAccRes->push_back((genDecayLVec[i].Pt() - matchPt)/genDecayLVec[i].Pt());
                            }

                            dRMin = 999.9;
                            for(int j = 0; j < cutMuVec->size(); ++j)
                            {
                                double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], (*cutMuVec)[j]);
                                if(dR < dRMin)
                                {
                                    dRMin = dR;
                                }
                            }
                            if(dRMin < 0.02)
                            {
                                genMatchIsoMuInAcc->push_back(&genDecayLVec[i]);
                                genMatchIsoMuInAccAct->push_back(genMuAct->back());
                            }
                        }
                    }

                    //Elec efficiency and acceptance
                    if(abs(genDecayPdgIdVec[i]) == 11)
                    {
                        genElec->push_back(&genDecayLVec[i]);
                        genElecAct->push_back(W_emu_pfActivityVec[index]);
                        if(AnaFunctions::passElectronAccOnly(genDecayLVec[i], AnaConsts::elesMiniIsoArr) && genDecayLVec[i].Pt() > 20)
                        {
                            genElecInAcc->push_back(&genDecayLVec[i]);
                            genElecInAccAct->push_back(genElecAct->back());
                            double dRMin = 999.9;
                            double matchPt = -999.9;
                            for(int j = 0; j < cutElecVecRecoOnly.size(); ++j)
                            {
                                double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], cutElecVecRecoOnly[j]);
                                if(dR < dRMin)
                                {
                                    dRMin = dR;
                                    matchPt = cutElecVecRecoOnly[j].Pt();
                                }
                            }
                            if(dRMin < 0.02)
                            {
                                genMatchElecInAcc->push_back(&genDecayLVec[i]);
                                genMatchElecInAccAct->push_back(genElecAct->back());
                                genMatchElecInAccRes->push_back((genDecayLVec[i].Pt() - matchPt)/genDecayLVec[i].Pt());
                            }

                            dRMin = 999.9;
                            for(int j = 0; j < cutElecVec->size(); ++j)
                            {
                                double dR = ROOT::Math::VectorUtil::DeltaR(genDecayLVec[i], (*cutElecVec)[j]);
                                if(dR < dRMin)
                                {
                                    dRMin = dR;
                                }
                            }
                            if(dRMin < 0.02)
                            {
                                genMatchIsoElecInAcc->push_back(&genDecayLVec[i]);
                                genMatchIsoElecInAccAct->push_back(genElecAct->back());
                            }
                        }
                    }
                }
            }

            double genZPt = -999.9, genZEta = -999.9, genZmass = -999.9, genZPhi;
            int nZ = 0;
            TLorentzVector genZ;
            int pdgIdZDec = 0;
            if(&genDecayPdgIdVec != nullptr)
            {
                for(int j = 0; j <  genDecayPdgIdVec.size(); ++j)
                {
                    if(abs(genDecayPdgIdVec[j]) == 23)
                    {
                        nZ++;
                        genZ = genDecayLVec[j];
                        genZPt = genDecayLVec[j].Pt();
                        genZEta = genDecayLVec[j].Eta();
                        genZPhi = genDecayLVec[j].Phi();
                        genZmass = genDecayLVec[j].M();
                    }
                }
                if(nZ > 1) std::cout << "!!!WARNING MORE THAN 1 Z FOUND!!!" << std::endl;

                if(&W_emuVec != nullptr)
                {
                    if(W_emuVec.size() == 0) pdgIdZDec = 15;
                    else if(W_emuVec.size() == 2)
                    {
                        if(abs(genDecayPdgIdVec[W_emuVec[0]]) == 11) pdgIdZDec = 11;
                        else if(abs(genDecayPdgIdVec[W_emuVec[0]]) == 13) pdgIdZDec = 13;
                    }
                }

            }

            bool passDiMuTrig  = nTriggerMuons >= 2;

            const double zMassMin = 81.0;
            const double zMass    = 91.0;
            const double zMassMax = 101.0;

            double zMuMassCurrent = 1.0e300, zEff = 1.0e100, zAcc = 1.0e100;
            TLorentzVector bestRecoMuZ;
            for(int i = 0; i < cutMuVec->size(); ++i)
            {
                if((*cutMuVec)[i].Pt() < minMuPt) continue;
                for(int j = 0; j < i && j < cutMuVec->size(); ++j)
                {
                    if((*cutMuVec)[j].Pt() < minMuPt) continue;
                    double zm = ((*cutMuVec)[i] + (*cutMuVec)[j]).M();
                    //if(zm > zMassMin && zm < zMassMax && fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    if(fabs(zm - zMass) < fabs(zMuMassCurrent - zMass))
                    {
                        bestRecoMuZ = (*cutMuVec)[i] + (*cutMuVec)[j];
                        zMuMassCurrent = zm;
                    }
                }
            }

            double zElecMassCurrent = 1.0e300;
            TLorentzVector bestRecoElecZ;
            for(int i = 0; i < cutElecVec->size(); ++i)
            {
                if((*cutElecVec)[i].Pt() < minElecPt) continue;
                for(int j = 0; j < i && j < cutElecVec->size(); ++j)
                {
                    if((*cutElecVec)[j].Pt() < minElecPt) continue;
                    double zm = ((*cutElecVec)[i] + (*cutElecVec)[j]).M();
                    //if(zm > zMassMin && zm < zMassMax && fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    if(fabs(zm - zMass) < fabs(zElecMassCurrent - zMass))
                    {
                        bestRecoElecZ = (*cutElecVec)[i] + (*cutElecVec)[j];
                        zElecMassCurrent = zm;
                    }
                }
            }

            double zElMuMassCurrent = 1.0e300;
            TLorentzVector bestRecoElMuZ;
            for(int i = 0; i < cutMuVec->size(); ++i)
            {
                if((*cutMuVec)[i].Pt() < minMuPt) continue;
                for(int j = 0; j < cutElecVec->size(); ++j)
                {
                    if((*cutElecVec)[j].Pt() < minMuPt) continue;
                    double zm = ((*cutMuVec)[i] + (*cutElecVec)[j]).M();
                    //if(zm > zMassMin && zm < zMassMax && fabs(zm - zMass) < fabs(zMassCurrent - zMass))
                    if(fabs(zm - zMass) < fabs(zElMuMassCurrent - zMass))
                    {
                        bestRecoElMuZ = (*cutMuVec)[i] + (*cutElecVec)[j];
                        zElMuMassCurrent = zm;
                    }
                }
            }

            TLorentzVector metV, metZ;
            metV.SetPtEtaPhiM(met, 0.0, metphi, 0.0);

            TLorentzVector bestRecoZ = (fabs(bestRecoElecZ.M() - zMass) > fabs(bestRecoMuZ.M() - zMass))?(bestRecoMuZ):(bestRecoElecZ);
            if(fabs(bestRecoZ.M() - zMass) > fabs(bestRecoElMuZ.M() - zMass)) bestRecoZ = bestRecoElMuZ;

            metZ.SetPtEtaPhiM(bestRecoZ.Pt(), 0.0, bestRecoZ.Phi(), 0.0);
            TLorentzVector cleanMet = metV + metZ;

            bool passDiMuSel   =  passEleVeto && (cutMuVec->size() == 2   && sumMuCharge == 0   && (*cutMuVec)[0].Pt() > highMuPt     && (*cutMuVec)[1].Pt() > minMuPt);
            bool passDiElecSel = passMuonVeto && (cutElecVec->size() == 2 && sumElecCharge == 0 && (*cutElecVec)[0].Pt() > highElecPt && (*cutElecVec)[1].Pt() > minElecPt);
            bool passElMuSel = (cutMuVec->size() == 1 && cutElecVec->size() == 1 && sumElecCharge == -sumMuCharge && (*cutMuVec)[0].Pt() > highMuPt && (*cutElecVec)[0].Pt() > minMuPt);

            bool passMuZinvSel   =  passEleVeto && (cutMuVec->size() == 2   && sumMuCharge == 0   && (*cutMuVec)[0].Pt() > highMuPt     && (*cutMuVec)[1].Pt() > minMuPt)     && (bestRecoMuZ.M() > zMassMin)   && (bestRecoMuZ.M() < zMassMax);
            bool passElecZinvSel = passMuonVeto && (cutElecVec->size() == 2 && sumElecCharge == 0 && (*cutElecVec)[0].Pt() > highElecPt && (*cutElecVec)[1].Pt() > minElecPt) && (bestRecoElecZ.M() > zMassMin) && (bestRecoElecZ.M() < zMassMax);
            bool passElMuZinvSel = (cutMuVec->size() == 1 && cutElecVec->size() == 1 && sumElecCharge == -sumMuCharge && (*cutMuVec)[0].Pt() > highMuPt && (*cutElecVec)[0].Pt() > minMuPt) && (bestRecoElMuZ.M() > zMassMin) && (bestRecoElMuZ.M() < zMassMax);

            double cutMuPt1 = -999.9;
            double cutMuPt2 = -999.9;
            if(cutMuVec->size() >= 1) cutMuPt1 = cutMuVec->at(0).Pt();
            if(cutMuVec->size() >= 2) cutMuPt2 = cutMuVec->at(1).Pt();

            double cutElecPt1 = -999.9;
            double cutElecPt2 = -999.9;
            if(cutElecVec->size() >= 1) cutElecPt1 = cutElecVec->at(0).Pt();
            if(cutElecVec->size() >= 2) cutElecPt2 = cutElecVec->at(1).Pt();

            double mindPhiMetJ = 999.9;
            int jc = 0;
            for(const TLorentzVector& jet : jetsLVec)
            {
                if(jc >= 3) break;
                jc++;
                mindPhiMetJ = std::min(mindPhiMetJ, fabs(ROOT::Math::VectorUtil::DeltaPhi(genZ, jet)));
            }
            //if(genZPt > 600) std::cout << "HELLO THERE!!!!" << std::endl;
            //if(genZPt > 600 && mindPhiMetJ < 0.5) std::cout << "BONJOUR!!! \t" << genZPt << "\t" << mindPhiMetJ << "\t" << run << "\t" << lumi << "\t" << event << std::endl;

            double bestRecoZPt = bestRecoZ.Pt();
            double cleanMetPt = cleanMet.Pt();
            //double cleanMet2Pt = cleanMet2.Pt();
            tr.registerDerivedVar("bestRecoZPt", bestRecoZPt);
            tr.registerDerivedVar("bestRecoZM", bestRecoZ.M());
            tr.registerDerivedVar("cleanMetPt", cleanMetPt);
            tr.registerDerivedVar("cleanMetPhi", cleanMet.Phi());
            //tr.registerDerivedVar("cleanMet2Pt", cleanMet2Pt);
            tr.registerDerivedVar("genHt", genHt);
            tr.registerDerivedVar("cutMuPt1", cutMuPt1);
            tr.registerDerivedVar("cutMuPt2", cutMuPt2);
            tr.registerDerivedVar("cutElecPt1", cutElecPt1);
            tr.registerDerivedVar("cutElecPt2", cutElecPt2);
            tr.registerDerivedVar("mindPhiMetJ", mindPhiMetJ);

            tr.registerDerivedVar("ZPtRes", (bestRecoZPt - genZPt)/genZPt);
            tr.registerDerivedVar("ZEtaRes", bestRecoZ.Eta() - genZEta);
            tr.registerDerivedVar("ZPhiRes", bestRecoZ.Phi() - genZPhi);
            tr.registerDerivedVar("ZMRes", (bestRecoZ.M() - genZmass)/genZmass);

            tr.registerDerivedVec("cutMuVec", cutMuVec);
            tr.registerDerivedVec("cutElecVec", cutElecVec);
            tr.registerDerivedVec("cutMuActivity", cutMuActivity);
            tr.registerDerivedVec("cutElecActivity", cutElecActivity);
            tr.registerDerivedVec("genMu", genMu);
            tr.registerDerivedVar("ngenMu", static_cast<double>(genMu->size()));
            tr.registerDerivedVec("genMuInAcc", genMuInAcc);
            tr.registerDerivedVec("genMuAct", genMuAct);
            tr.registerDerivedVar("ngenMuInAcc", static_cast<double>(genMuInAcc->size()));
            tr.registerDerivedVec("genMuInAccAct", genMuInAccAct);
            tr.registerDerivedVec("genMatchMuInAcc", genMatchMuInAcc);
            tr.registerDerivedVec("genMatchMuInAccRes", genMatchMuInAccRes);
            tr.registerDerivedVec("genMatchIsoMuInAcc", genMatchIsoMuInAcc);
            tr.registerDerivedVar("ngenMatchMuInAcc", static_cast<double>(genMatchMuInAcc->size()));
            tr.registerDerivedVec("genMatchMuInAccAct", genMatchMuInAccAct);
            tr.registerDerivedVec("genMatchIsoMuInAccAct", genMatchIsoMuInAccAct);

            tr.registerDerivedVec("genElec", genElec);
            tr.registerDerivedVar("ngenElec", static_cast<double>(genElec->size()));
            tr.registerDerivedVec("genElecInAcc", genElecInAcc);
            tr.registerDerivedVec("genElecAct", genElecAct);
            tr.registerDerivedVar("ngenElecInAcc", static_cast<double>(genElecInAcc->size()));
            tr.registerDerivedVec("genElecInAccAct", genElecInAccAct);
            tr.registerDerivedVec("genMatchElecInAcc", genMatchElecInAcc);
            tr.registerDerivedVec("genMatchElecInAccRes", genMatchElecInAccRes);
            tr.registerDerivedVec("genMatchIsoElecInAcc", genMatchIsoElecInAcc);
            tr.registerDerivedVar("ngenMatchElecInAcc", static_cast<double>(genMatchElecInAcc->size()));
            tr.registerDerivedVec("genMatchElecInAccAct", genMatchElecInAccAct);
            tr.registerDerivedVec("genMatchIsoElecInAccAct", genMatchIsoElecInAccAct);

            tr.registerDerivedVar("genZPt", genZPt);
            tr.registerDerivedVar("genZEta", genZEta);
            tr.registerDerivedVar("genZmass", genZmass);
            tr.registerDerivedVar("pdgIdZDec", pdgIdZDec);
            tr.registerDerivedVar("passDiMuIsoTrig", passDiMuTrig);
            tr.registerDerivedVar("passSingleMu45", muTrigMu45);

            tr.registerDerivedVar("passDiMuSel", passDiMuSel);
            tr.registerDerivedVar("passDiElecSel", passDiElecSel);
	    tr.registerDerivedVar("passElMuSel", passElMuSel);

            tr.registerDerivedVar("passMuZinvSel", passMuZinvSel);
            tr.registerDerivedVar("passElecZinvSel", passElecZinvSel);
            tr.registerDerivedVar("passElMuZinvSel", passElMuZinvSel);
            bool passKoushik = false;
            bool passQCD = false;
            if(  nMuons == 1 && nElectrons == AnaConsts::nElectronsSel ) {
                 passKoushik = true;
            }
           if(met/calomet<3.0 && calomet>80){
                passQCD = true;
           }
           tr.registerDerivedVar("passKoushik", passKoushik);
           tr.registerDerivedVar("passQCD", passQCD);
        }
/*
       void Singlemuon(NTupleReader& tr)
      {
            if(  nMuons_CUT == 1 && nElectrons == AnaConsts::nElectronsSel ) {

            }
       }
*/
    public:
        LepInfo()
        {
            tr3 = new TRandom3();
        }

        void operator()(NTupleReader& tr)
        {
            lepInfo(tr);
        }
    };

    class Fakebtagvectors
    {
    private:
        void fakebtagvectors(NTupleReader& tr)
        {
            const std::vector<TLorentzVector>& jetsLVecLepCleaned = tr.getVec<TLorentzVector>("jetsLVecLepCleaned");
            const std::vector<double>& cleanJetpt30ArrBTag = tr.getVec<double>("recoJetsBtag_0_LepCleaned");

            double maxCSV = 0.0;
            double secCSV = 0.0;
            double tenCSV = 0.0;
            int iMaxCSV = -1;
            int iSecCSV = -1;
            int iTenCSV = -1;

            if(jetsLVecLepCleaned.size() != cleanJetpt30ArrBTag.size()) std::cout << "fakebtagvectors(...): Vector size missmatch!!!!" << std::endl;

            int njet = 0;

            //find index of 3 highest CSV values
            for(int i = 0; i < cleanJetpt30ArrBTag.size(); ++i)
            {
                //Skip jets which cannot pass bTag Acceptance requirements
                if(!AnaFunctions::jetPassCuts(jetsLVecLepCleaned[i], AnaConsts::bTagArr)) continue;

                //count possible fake b-jets
                njet++;

                if(cleanJetpt30ArrBTag[i] > maxCSV)
                {
                    tenCSV = secCSV;
                    secCSV = maxCSV;
                    maxCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = iSecCSV;
                    iSecCSV = iMaxCSV;
                    iMaxCSV = i;
                }
                else if(cleanJetpt30ArrBTag[i] > secCSV)
                {
                    tenCSV = secCSV;
                    secCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = iSecCSV;
                    iSecCSV = i;
                }
                else if(cleanJetpt30ArrBTag[i] > tenCSV)
                {
                    tenCSV = cleanJetpt30ArrBTag[i];
                    iTenCSV = i;
                }
            }

            std::vector<double>* cleanJetpt30ArrBTag1fake = new std::vector<double>(cleanJetpt30ArrBTag);
            std::vector<double>* cleanJetpt30ArrBTag2fake = new std::vector<double>(cleanJetpt30ArrBTag);
            std::vector<double>* cleanJetpt30ArrBTag3fake = new std::vector<double>(cleanJetpt30ArrBTag);
            std::vector<double>* fakedCSVValues = new std::vector<double>();

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag1fake)[iMaxCSV] = 0.99;

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag2fake)[iMaxCSV] = 0.99;
            if(iSecCSV >= 0) (*cleanJetpt30ArrBTag2fake)[iSecCSV] = 0.99;

            if(iMaxCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iMaxCSV] = 0.99;
            if(iSecCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iSecCSV] = 0.99;
            if(iTenCSV >= 0) (*cleanJetpt30ArrBTag3fake)[iTenCSV] = 0.99;

            if(iMaxCSV >= 0) fakedCSVValues->push_back(maxCSV);
            if(iSecCSV >= 0) fakedCSVValues->push_back(secCSV);
            if(iTenCSV >= 0) fakedCSVValues->push_back(tenCSV);

            //Calculate the combinatoric weights for b-jet faking
            double weight1fakeb = TMath::Binomial(njet, 1);
            double weight2fakeb = TMath::Binomial(njet, 2);
            double weight3fakeb = TMath::Binomial(njet, 3);
            //check for nans
            if(weight1fakeb != weight1fakeb) weight1fakeb = 0.0;
            if(weight2fakeb != weight2fakeb) weight2fakeb = 0.0;
            if(weight3fakeb != weight3fakeb) weight3fakeb = 0.0;

            tr.registerDerivedVar("weight1fakeb", weight1fakeb);
            tr.registerDerivedVar("weight2fakeb", weight2fakeb);
            tr.registerDerivedVar("weight3fakeb", weight3fakeb);

            tr.registerDerivedVec("cleanJetpt30ArrBTag1fake", cleanJetpt30ArrBTag1fake);
            tr.registerDerivedVec("cleanJetpt30ArrBTag2fake", cleanJetpt30ArrBTag2fake);
            tr.registerDerivedVec("cleanJetpt30ArrBTag3fake", cleanJetpt30ArrBTag3fake);
            tr.registerDerivedVec("fakedCSVValues", fakedCSVValues);
            tr.registerDerivedVar("maxCSV", maxCSV);
        }
    public:

        void operator()(NTupleReader& tr)
        {
            fakebtagvectors(tr);
        }
    };

    class GetSearchBin
    {
    private:
        SearchBins sbins;

        void getSearchBin(NTupleReader& tr)
        {
            const int& cntCSVS = tr.getVar<int>("cntCSVSZinv");
            const int& nTopCandSortedCnt = tr.getVar<int>("nTopCandSortedCntZinv");
            const int& nTopCandSortedCnt1b = tr.getVar<int>("nTopCandSortedCntZinv1b");
            const int& nTopCandSortedCnt2b = tr.getVar<int>("nTopCandSortedCntZinv2b");
            const int& nTopCandSortedCnt3b = tr.getVar<int>("nTopCandSortedCntZinv3b");
            const double& cleanMet = tr.getVar<double>("cleanMetPt");
            const double& cleanMetPhi = tr.getVar<double>("cleanMetPhi");
            const double& MT2 = tr.getVar<double>("best_had_brJet_MT2Zinv");
            const double& MT2_1b = tr.getVar<double>("best_had_brJet_MT2Zinv1b");
            const double& MT2_2b = tr.getVar<double>("best_had_brJet_MT2Zinv2b");
            const double& MT2_3b = tr.getVar<double>("best_had_brJet_MT2Zinv3b");
            const double& weight1fakeb = tr.getVar<double>("weight1fakeb");
            const double& weight2fakeb = tr.getVar<double>("weight2fakeb");
            const double& weight3fakeb = tr.getVar<double>("weight3fakeb");

            const double& nJet1bfakeWgt = tr.getVar<double>("nJet1bfakeWgt");
            const double& nJet2bfakeWgt = tr.getVar<double>("nJet2bfakeWgt");
            const double& nJet3bfakeWgt = tr.getVar<double>("nJet3bfakeWgt");

            int nSearchBin = sbins.find_Binning_Index(cntCSVS, nTopCandSortedCnt, MT2, cleanMet);

            std::vector<std::pair<double, double> > * nb0Bins = new std::vector<std::pair<double, double> >();
            std::vector<std::pair<double, double> > * nb0NJwBins = new std::vector<std::pair<double, double> >();
            std::vector<double> * nb0BinsNW = new std::vector<double>();

            //weights based on total N(b) yields vs. N(b) = 0 control region
            //These weights are derived from the rato of events in the N(t) = 1, 2, 3 bins after all baseline cuts except b tag between the
            //N(b) = 0 control region and each N(b) signal region using Z->nunu MC.  They account for both the combinatoric reweighting factor
            //as well as the different event yields between the control region and each signal region.
            const double wnb01 = 3.820752e-02;//3.478840e-02;//6.26687e-2;
            const double wnb02 = 2.946461e-03;//2.586369e-03;//5.78052e-3;
            const double wnb03 = 1.474770e-04;//1.640077e-04;//7.08235e-4;

            // weights to apply when doing b-faking
            const double w1b = wnb01 * weight1fakeb;
            const double w2b = wnb02 * weight2fakeb;
            const double w3b = wnb03 * weight3fakeb;

            if(cntCSVS == 0)
            {
                //nb0Bins->push_back(std::make_pair(find_Binning_Index(0, nTopCandSortedCnt, MT2, cleanMet), 1.0));
                nb0Bins->emplace_back(std::pair<double, double>(sbins.find_Binning_Index(1, nTopCandSortedCnt1b, MT2_1b, cleanMet), wnb01 * weight1fakeb));
                nb0Bins->emplace_back(std::pair<double, double>(sbins.find_Binning_Index(2, nTopCandSortedCnt2b, MT2_2b, cleanMet), wnb02 * weight2fakeb));
                nb0Bins->emplace_back(std::pair<double, double>(sbins.find_Binning_Index(3, nTopCandSortedCnt3b, MT2_3b, cleanMet), wnb03 * weight3fakeb));

                nb0NJwBins->emplace_back(std::pair<double, double>(sbins.find_Binning_Index(1, nTopCandSortedCnt1b, MT2_1b, cleanMet), nJet1bfakeWgt));
                nb0NJwBins->emplace_back(std::pair<double, double>(sbins.find_Binning_Index(2, nTopCandSortedCnt2b, MT2_2b, cleanMet), nJet2bfakeWgt));
                nb0NJwBins->emplace_back(std::pair<double, double>(sbins.find_Binning_Index(3, nTopCandSortedCnt3b, MT2_3b, cleanMet), nJet3bfakeWgt));

                nb0BinsNW->emplace_back(sbins.find_Binning_Index(1, nTopCandSortedCnt1b, MT2_1b, cleanMet));
                nb0BinsNW->emplace_back(sbins.find_Binning_Index(2, nTopCandSortedCnt2b, MT2_2b, cleanMet));
                nb0BinsNW->emplace_back(sbins.find_Binning_Index(3, nTopCandSortedCnt3b, MT2_3b, cleanMet));
            }

            tr.registerDerivedVar("nSearchBin", nSearchBin);
            tr.registerDerivedVec("nb0BinsNW", nb0BinsNW);
            tr.registerDerivedVec("nb0Bins", nb0Bins);
            tr.registerDerivedVec("nb0NJwBins", nb0NJwBins);
            tr.registerDerivedVar("weight1fakebComb", w1b);
            tr.registerDerivedVar("weight2fakebComb", w2b);
            tr.registerDerivedVar("weight3fakebComb", w3b);
        }

    public:

        GetSearchBin(std::string sb_era) : sbins(sb_era) {}

        void operator()(NTupleReader& tr)
        {
            getSearchBin(tr);
        }
    };


    class TriggerInfo
    {
    private:
	int indexMuTrigger;
	int indexElecTrigger;
	int indexSingleMuTrigger;
        int indexHTMHTTrigger;
        bool miniTuple_;

	double GetMuonTriggerEff(const double& muEta) 
	{
                if     (-2.6 <= muEta && muEta < -2.2) return 0.016; //+/-0.013/0.008
                else if(-2.2 <= muEta && muEta < -1.8) return 0.680; //+/-0.017/0.017
                else if(-1.8 <= muEta && muEta < -1.4) return 0.787; //+/-0.012/0.012
                else if(-1.4 <= muEta && muEta < -1.0) return 0.866; //+/-0.008/0.009
                else if(-1.0 <= muEta && muEta < -0.6) return 0.902; //+/-0.006/0.007
                else if(-0.6 <= muEta && muEta < -0.2) return 0.892; //+/-0.006/0.007
                else if(-0.2 <= muEta && muEta <  0.2) return 0.927; //+/-0.005/0.006
                else if( 0.2 <= muEta && muEta <  0.6) return 0.892; //+/-0.006/0.007
                else if( 0.6 <= muEta && muEta <  1.0) return 0.911; //+/-0.006/0.007
                else if( 1.0 <= muEta && muEta <  1.4) return 0.864; //+/-0.008/0.009
                else if( 1.4 <= muEta && muEta <  1.8) return 0.808; //+/-0.011/0.012
                else if( 1.8 <= muEta && muEta <  2.2) return 0.652; //+/-0.017/0.017
                else if( 2.2 <= muEta && muEta <  2.6) return 0.026; //+/-0.015/0.01
                else                                   return 0.000;
	}

	double GetTriggerEffWeight(const double& met, const double& ht) 
	{
	    if (ht<1000)
	    {
		if (met<25) return 0.002;
		else if (met<50) return 0.003;
		else if (met<75) return 0.011;
		else if (met<100) return 0.053;
		else if (met<125) return 0.211;
		else if (met<150) return 0.484;
		else if (met<175) return 0.751;
		else if (met<200) return 0.897;
		else if (met<275) return 0.970;
		else if (met<400) return 0.986; 
		else return 0.97;
	    } 
	    else 
	    {
		if (met<25) return 0.014;
		else if (met<50) return 0.019;
		else if (met<75) return 0.037;
		else if (met<100) return 0.089;
		else if (met<125) return 0.204;
		else if (met<150) return 0.358;
		else if (met<175) return 0.568;
		else if (met<200) return 0.747;
		else if (met<275) return 0.907;
		else if (met<400) return 0.988; 
		else return 0.972;
	    }
	}
	double GetTriggerEffStatUncHi(const double& met, const double& ht) 
	{
	    if (ht<1000)
	    {
		if (met<75) return 0.0;
		else if (met<100) return 0.001;
		else if (met<125) return 0.004;
		else if (met<150) return 0.006;
		else if (met<175) return 0.006;
		else if (met<200) return 0.005;
		else if (met<275) return 0.002;
		else if (met<400) return 0.002; 
		else return 0.007;
	    } 
	    else 
	    {
		if (met<25) return 0.002;
		else if (met<50) return 0.002;
		else if (met<75) return 0.003;
		else if (met<100) return 0.005;
		else if (met<125) return 0.010;
		else if (met<150) return 0.015;
		else if (met<175) return 0.020;
		else if (met<200) return 0.021;
		else if (met<275) return 0.010;
		else if (met<400) return 0.005; 
		else return 0.009;
	    }
	}
	double GetTriggerEffStatUncLo(const double& met, const double& ht) 
	{
	    if (ht<1000)
	    {
		if (met<75) return 0.0;
		else if (met<100) return 0.001;
		else if (met<125) return 0.004;
		else if (met<150) return 0.006;
		else if (met<175) return 0.007;
		else if (met<200) return 0.006;
		else if (met<275) return 0.003;
		else if (met<400) return 0.003; 
		else return 0.009;
	    } 
	    else 
	    {
		if (met<25) return 0.002;
		else if (met<50) return 0.002;
		else if (met<75) return 0.002;
		else if (met<100) return 0.005;
		else if (met<125) return 0.010;
		else if (met<150) return 0.015;
		else if (met<175) return 0.020;
		else if (met<200) return 0.022;
		else if (met<275) return 0.011;
		else if (met<400) return 0.007; 
		else return 0.012;
	    }
	}
	double GetTriggerEffSystUncHi(const double& met, const double& ht) 
	{
	    return 0.0;
	    /* if (met<100) return 0.0272; */
	    /* else if (met<150) return 0.0872; */
	    /* else if (met<175) return 0.1505; */
	    /* else if (met<200) return 0.0423; */
	    /* else if (met<275) return 0.0112; */
	    /* else if (met<400) return 0.0001;  */
	    /* else return 0.0; */
	}
	double GetTriggerEffSystUncLo(const double& met, const double& ht) 
	{
	    return 0.0;
	    /* if (met<100) return 0.0120; */
	    /* else if (met<150) return 0.0872; */
	    /* else if (met<175) return 0.1505; */
	    /* else if (met<200) return 0.0792; */
	    /* else if (met<275) return 0.0112; */
	    /* else if (met<400) return 0.0001;  */
	    /* else return 0.0018; */
	}

        void triggerInfo(NTupleReader& tr)
        {
            const std::vector<std::string>& triggerNames = tr.getVec<std::string>("TriggerNames");
            const std::vector<int>& passTrigger          = tr.getVec<int>("PassTrigger");

	    bool passMuTrigger = false;
	    bool passElecTrigger = false;
	    bool passSingleMuTrigger = false;
            bool passhtmhtTrigger = false;

	    const std::string muTrigName = "HLT_Mu45_eta2p1_v";
	    const std::string elecTrigName = "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v";
	    const std::string singleMuTrigName = "HLT_Mu15_IsoVVVL_PFHT350_v";
            const std::string htmhtTrigName = "HLT_PFHT300_PFMET100_v";

	    // Find the index of our triggers if we don't know them already
	    if(indexMuTrigger == -1 || indexElecTrigger == -1)
	    {
		for(int i = 0; i < triggerNames.size(); ++i)
		{
		    if(triggerNames[i].find(muTrigName) != std::string::npos)
		    {
			indexMuTrigger = i;
		    }
		    else if(triggerNames[i].find(elecTrigName) != std::string::npos)
		    {
			indexElecTrigger = i;
		    }
		    else if(triggerNames[i].find(singleMuTrigName) != std::string::npos)
		    {
			indexSingleMuTrigger = i;
		    }
                    else if(triggerNames[i].find(htmhtTrigName) != std::string::npos)
                    {
                        indexHTMHTTrigger = i;
                    }
		}
	    }
	    if(indexMuTrigger != -1 && indexElecTrigger != -1 && indexSingleMuTrigger != -1 && indexHTMHTTrigger != -1)
	    {
		// Check if the event passes the trigger, and double check that we are looking at the right trigger
		if(triggerNames[indexMuTrigger].find(muTrigName) != std::string::npos && passTrigger[indexMuTrigger])
		    passMuTrigger = true;
		if(triggerNames[indexElecTrigger].find(elecTrigName) != std::string::npos && passTrigger[indexElecTrigger])
		    passElecTrigger = true;
		if(triggerNames[indexSingleMuTrigger].find(singleMuTrigName) != std::string::npos && passTrigger[indexSingleMuTrigger])
		    passSingleMuTrigger = true;
                if(triggerNames[indexHTMHTTrigger].find(htmhtTrigName) != std::string::npos && passTrigger[indexHTMHTTrigger])
                    passhtmhtTrigger = true;
	    }
	    else
	    {
		std::cout << "Could not find trigger in the list of trigger names" << std::endl;
	    }

	    tr.registerDerivedVar("passMuTrigger",passMuTrigger);
	    tr.registerDerivedVar("passElecTrigger",passElecTrigger);
	    tr.registerDerivedVar("passSingleMuTrigger", passSingleMuTrigger);
            tr.registerDerivedVar("passhtmhtTrigger", passhtmhtTrigger);
        }

        void triggerInfoMC(NTupleReader& tr)
        {
            const double& met                            = tr.getVar<double>("cleanMetPt");
            const double& ht                             = tr.getVar<double>("HTZinv");
            const std::vector<TLorentzVector>& cutMuVec  = tr.getVec<TLorentzVector>("cutMuVec");

	    // MC trigger efficiencies
	    double triggerEff = GetTriggerEffWeight(met,ht);
	    double triggerEffStatUncUp = GetTriggerEffStatUncHi(met,ht);
	    double triggerEffSystUncUp = GetTriggerEffSystUncHi(met,ht);
	    double triggerEffUncUp     = TMath::Sqrt(triggerEffStatUncUp*triggerEffStatUncUp + triggerEffSystUncUp*triggerEffSystUncUp);
	    double triggerEffStatUncDown = GetTriggerEffStatUncLo(met,ht);
	    double triggerEffSystUncDown = GetTriggerEffSystUncLo(met,ht);
	    double triggerEffUncDown     = TMath::Sqrt(triggerEffStatUncDown*triggerEffStatUncDown + triggerEffSystUncDown*triggerEffSystUncDown);

            //Calculate muon trigger weights
            double muTrigWgt = 0.0;
            if(cutMuVec.size() >= 2 && cutMuVec[0].Pt() > 45 && cutMuVec[1].Pt() > 45)
            {
                double muEff1 = GetMuonTriggerEff(cutMuVec[0].Eta());
                double muEff2 = GetMuonTriggerEff(cutMuVec[1].Eta());

                muTrigWgt = 1 - (1 - muEff1)*(1 - muEff2);
            }
            else if(cutMuVec.size() >= 1 && cutMuVec[0].Pt() > 45)
            {
                //For events with only 1 muon (emu events in particular or events with a subleading muon below 45 GeV) just use the single muon eff
                muTrigWgt = GetMuonTriggerEff(cutMuVec[0].Eta());
            }

	    tr.registerDerivedVar("TriggerEffMC",triggerEff);
	    tr.registerDerivedVar("TriggerEffUpMC",triggerEff+triggerEffUncUp);
	    tr.registerDerivedVar("TriggerEffDownMC",triggerEff-triggerEffUncDown);

            tr.registerDerivedVar("muTrigWgt", muTrigWgt);
        }

    public:
	TriggerInfo(bool miniTuple = false)
	{
	    indexMuTrigger = -1;
	    indexElecTrigger = -1;
	    indexSingleMuTrigger = -1;
            indexHTMHTTrigger = -1;
            miniTuple_ = miniTuple;
	}

	void operator()(NTupleReader& tr)
	{
	    if(!miniTuple_) triggerInfo(tr);
            triggerInfoMC(tr);
	}

    };

    class SystematicPrep
    {
    private:

        void systematicPrep(NTupleReader& tr)
        {
            const std::vector<TLorentzVector>& jetsLVec  = tr.getVec<TLorentzVector>("jetsLVecLepCleaned");
            const std::vector<double>& recoJetsJecUnc    = tr.getVec<double>("recoJetsJecUncLepCleaned");

            const std::vector<double>& metMagUp   = tr.getVec<double>("metMagUp");
            const std::vector<double>& metMagDown = tr.getVec<double>("metMagDown");
            const std::vector<double>& metPhiUp   = tr.getVec<double>("metPhiUp");
            const std::vector<double>& metPhiDown = tr.getVec<double>("metPhiDown");

            const double& met    = tr.getVar<double>("met");
            const double& metphi = tr.getVar<double>("metphi");

            std::vector<TLorentzVector> *jetLVecUp = new std::vector<TLorentzVector>;
            std::vector<TLorentzVector> *jetLVecDn = new std::vector<TLorentzVector>;

            std::vector<double> *dPtMet = new std::vector<double>;
            std::vector<double> *dPhiMet = new std::vector<double>;

            double metUp = 0.0, metDn = 99990.0;

            for(int iMet = 0; iMet < metMagUp.size(); ++iMet)
            {
                metUp = std::max(metUp, metMagUp[iMet]);
                metDn = std::min(metDn, metMagDown[iMet]);
                
                dPtMet->push_back((metMagUp[iMet] - met)/met);
                dPtMet->push_back((metMagDown[iMet] - met)/met);
                dPhiMet->push_back(TVector2::Phi_mpi_pi(metPhiUp[iMet] - metphi));
                dPhiMet->push_back(TVector2::Phi_mpi_pi(metPhiDown[iMet] - metphi));
            }

            for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
            {
                jetLVecUp->push_back(jetsLVec[iJet] * (1 + recoJetsJecUnc[iJet]));
                jetLVecDn->push_back(jetsLVec[iJet] * (1 - recoJetsJecUnc[iJet]));
            }

            tr.registerDerivedVar("metMEUUp", metUp);
            tr.registerDerivedVar("metMEUDn", metDn);

            tr.registerDerivedVec("dPtMet", dPtMet);
            tr.registerDerivedVec("dPhiMet", dPhiMet);

            tr.registerDerivedVec("jetLVecUp", jetLVecUp);
            tr.registerDerivedVec("jetLVecDn", jetLVecDn);
        }

    public:
	SystematicPrep()
	{

	}

	void operator()(NTupleReader& tr)
	{
	    systematicPrep(tr);
	}

    };

    class SystematicCalc
    {
    private:
        SearchBins sbins;

        void systematicCalc(NTupleReader& tr)
        {
            const int& cntCSVSJEUUp = tr.getVar<int>("cntCSVSZinvJEUUp");
            const int& nTopCandSortedCntJEUUp = tr.getVar<int>("nTopCandSortedCntZinvJEUUp");
            const double& MT2JEUUp = tr.getVar<double>("best_had_brJet_MT2ZinvJEUUp");

            const int& cntCSVSJEUDn = tr.getVar<int>("cntCSVSZinvJEUDn");
            const int& nTopCandSortedCntJEUDn = tr.getVar<int>("nTopCandSortedCntZinvJEUDn");
            const double& MT2JEUDn = tr.getVar<double>("best_had_brJet_MT2ZinvJEUDn");

            const double& cleanMet = tr.getVar<double>("cleanMetPt");

            const int& cntCSVSMEUUp = tr.getVar<int>("cntCSVSZinvMEUUp");
            const int& nTopCandSortedCntMEUUp = tr.getVar<int>("nTopCandSortedCntZinvMEUUp");
            const double& MT2MEUUp = tr.getVar<double>("best_had_brJet_MT2ZinvMEUUp");
            const double& cleanMetMEUUp = tr.getVar<double>("metMEUUp");

            const int& cntCSVSMEUDn = tr.getVar<int>("cntCSVSZinvMEUDn");
            const int& nTopCandSortedCntMEUDn = tr.getVar<int>("nTopCandSortedCntZinvMEUDn");
            const double& MT2MEUDn = tr.getVar<double>("best_had_brJet_MT2ZinvMEUDn");
            const double& cleanMetMEUDn = tr.getVar<double>("metMEUDn");

            int nSearchBinJEUUp = sbins.find_Binning_Index(cntCSVSJEUUp, nTopCandSortedCntJEUUp, MT2JEUUp, cleanMet);
            int nSearchBinJEUDn = sbins.find_Binning_Index(cntCSVSJEUDn, nTopCandSortedCntJEUDn, MT2JEUDn, cleanMet);

            int nSearchBinMEUUp = sbins.find_Binning_Index(cntCSVSMEUUp, nTopCandSortedCntMEUUp, MT2MEUUp, cleanMetMEUUp);
            int nSearchBinMEUDn = sbins.find_Binning_Index(cntCSVSMEUDn, nTopCandSortedCntMEUDn, MT2MEUDn, cleanMetMEUDn);
            
            tr.registerDerivedVar("nSearchBinJEUUp", nSearchBinJEUUp);
            tr.registerDerivedVar("nSearchBinJEUDn", nSearchBinJEUDn);

            tr.registerDerivedVar("nSearchBinMEUUp", nSearchBinMEUUp);
            tr.registerDerivedVar("nSearchBinMEUDn", nSearchBinMEUDn);
        }

    public:
        SystematicCalc(std::string sb_era) : sbins(sb_era)
	{
            
	}

	void operator()(NTupleReader& tr)
	{
	    systematicCalc(tr);
	}

    };

    class PrepareTopVars
    {
    private:

        topTagger::type3TopTagger t3tagger;
        TopTagger* tt;
        Mt2::ChengHanBisect_Mt2_332_Calculator mt2Calculator;

        void prepareTopVars(NTupleReader& tr)
        {
            const std::vector<TLorentzVector>& jetsLVec  = tr.getVec<TLorentzVector>("jetsLVec");
            const std::vector<double>& recoJetsBtag      = tr.getVec<double>("recoJetsBtag_0");

            const double& met    = tr.getVar<double>("met");
            const double& metphi = tr.getVar<double>("metphi");
            
            TLorentzVector metLVec;
            metLVec.SetPtEtaPhiM(met, 0, metphi, 0);

            std::vector<TLorentzVector> jetsLVec_forTagger;
            std::vector<double> recoJetsBtag_forTagger;
            AnaFunctions::prepareJetsForTagger(jetsLVec, recoJetsBtag, jetsLVec_forTagger, recoJetsBtag_forTagger);

            int cntCSVS = AnaFunctions::countCSVS(jetsLVec_forTagger, recoJetsBtag_forTagger, AnaConsts::cutCSVS, AnaConsts::bTagArr);

            t3tagger.processEvent(jetsLVec_forTagger, recoJetsBtag_forTagger, metLVec);
            int nTops = t3tagger.nTopCandSortedCnt;

            std::vector<TLorentzVector> *vTops = new std::vector<TLorentzVector>();

            std::vector<std::vector<TLorentzVector> >* vTopConstituents = new std::vector<std::vector<TLorentzVector> >();;
            
            //variables for tri-jet tops 
            std::vector<std::vector<TLorentzVector> >* v3j_diJets = new std::vector<std::vector<TLorentzVector> >();;

            std::vector<double>* v3j_dr12 = new std::vector<double>();
            std::vector<double>* v3j_dr23 = new std::vector<double>();
            std::vector<double>* v3j_dr13 = new std::vector<double>();
            std::vector<double>* v3j_cm_dr12 = new std::vector<double>();
            std::vector<double>* v3j_cm_dr13 = new std::vector<double>();
            std::vector<double>* v3j_cm_dr23 = new std::vector<double>();
            
            for(int it = 0; it < nTops; it++)
            {
                TLorentzVector topLVec = t3tagger.buildLVec(jetsLVec_forTagger,
                                                            t3tagger.finalCombfatJets[t3tagger.ori_pickedTopCandSortedVec[it]]);
                vTops->push_back(topLVec);

                std::vector<TLorentzVector> tmpVec;
                for(const int& jetIndex : t3tagger.finalCombfatJets[t3tagger.ori_pickedTopCandSortedVec[it]])
                {
                    tmpVec.emplace_back(jetsLVec_forTagger[jetIndex]);
                }
                vTopConstituents->emplace_back(tmpVec);

                //Make properties of triplet tops here
                if(tmpVec.size() == 3)
                {
                    double dR12 = ROOT::Math::VectorUtil::DeltaR(tmpVec[0], tmpVec[1]);
                    double dR23 = ROOT::Math::VectorUtil::DeltaR(tmpVec[1], tmpVec[2]);
                    double dR13 = ROOT::Math::VectorUtil::DeltaR(tmpVec[0], tmpVec[2]);

                    v3j_dr12->push_back(dR12);
                    v3j_dr23->push_back(dR23);
                    v3j_dr13->push_back(dR13);

                    //jet copies for playing with
                    TLorentzVector j1 = tmpVec[0];
                    TLorentzVector j2 = tmpVec[1];
                    TLorentzVector j3 = tmpVec[2];

                    //subtract out the top boost to get to top CM frame
                    j1 = j1 - topLVec;
                    j2 = j2 - topLVec;
                    j3 = j3 - topLVec;

                    double cm_dR12 = ROOT::Math::VectorUtil::DeltaR(j1, j2);
                    double cm_dR23 = ROOT::Math::VectorUtil::DeltaR(j2, j3);
                    double cm_dR13 = ROOT::Math::VectorUtil::DeltaR(j1, j3);
                    
                    v3j_cm_dr12->push_back(cm_dR12);
                    v3j_cm_dr13->push_back(cm_dR23);
                    v3j_cm_dr23->push_back(cm_dR13);
                    
                    //make di-jet systems 
                    v3j_diJets->push_back({tmpVec[0] + tmpVec[1], tmpVec[1] + tmpVec[2], tmpVec[0] + tmpVec[2]});
                }
            }

            //New Tagger starts here
            //prep input object (constituent) vector
            std::vector<Constituent> constituents = ttUtility::packageConstituents(jetsLVec_forTagger, recoJetsBtag_forTagger);
            //run tagger
            tt->runTagger(constituents);
            //retrieve results
            const TopTaggerResults& ttr = tt->getResults();

            std::vector<TLorentzVector> *vTopsNew = new std::vector<TLorentzVector>();
            for(auto& top : ttr.getTops())
            {
                vTopsNew->push_back(top->p());
            }

            //compare tops from each tagger here
            std::vector<TLorentzVector> *vdTops = new std::vector<TLorentzVector>();
            if(vTops->size() == vTopsNew->size())
            {
                for(int it = 0; it < vTops->size(); it++)
                {
                    vdTops->push_back(vTops->at(it) - vTopsNew->at(it));
                }
            }
            else if(cntCSVS >= 1 && jetsLVec.size() >= 4)
            {
                //std::cout << std::endl;
                //std::cout << "nJets: " << jetsLVec_forTagger.size() << "   nb: " << cntCSVS << std::endl;
                //std::cout << "NEW AND OLD TOP VECTOR ARE NOT SAME SIZE!!!!! old: " << vTops->size() << "   new: " << vTopsNew->size() << std::endl;
                //if(vTopsNew->size() >= 1) std::cout << "NEW TOP nConstituents: " << ttr.getTops()[0]->getConstituents().size() << "   pt: " << ttr.getTops()[0]->p().Pt() << "   eta: " << ttr.getTops()[0]->p().Eta() << "   phi: " << ttr.getTops()[0]->p().Phi() << "   M: " << ttr.getTops()[0]->p().M() << "   dRmax: " << ttr.getTops()[0]->getDRmax() << std::endl;
                //int nbInTop = 0;
                //for(auto& jet : ttr.getTops()[0]->getConstituents()) if(jet->getBTagDisc() > 0.89) ++nbInTop;
                //std::cout << "NEW nb in top: " << nbInTop << std::endl;
                //
                //if(t3tagger.topJetCand_idx123Vec.size())
                //{
                //    TLorentzVector topLVec = t3tagger.buildLVec(jetsLVec_forTagger, t3tagger.topJetCand_idx123Vec[0]);
                //    std::cout << "OLD TOP nConstituents: " << t3tagger.topJetCand_idx123Vec[0].size() << "   pt: " << topLVec.Pt() << "   eta: " << topLVec.Eta() << "   phi: " << topLVec.Phi() << "   M: " << topLVec.M() << std::endl;
                //}
                //std::cout << "NEW TOP nCandidates: " << ttr.getTopCandidates().size() << std::endl;
                //std::cout << "OLD TOP nCandidates: " << t3tagger.finalCombfatJets.size() << std::endl;
                //std::cout << std::endl;
            }

            //Rsys variables 
            TLorentzVector oldRsysVec;// = t3tagger.best_had_brJet;
            TLorentzVector newRsysVec = ttr.getRsys().p();

            double oldMT2 = t3tagger.best_had_brJet_MT2;
            double newMT2 = -999.9;
            double deltaR = -999;
            if(ttr.getTops().size() >= 1)
            {
                const double massOfSystemA = ttr.getTops()[0]->P().M(); // GeV
                const double pxOfSystemA   = ttr.getTops()[0]->P().Px(); // GeV
                const double pyOfSystemA   = ttr.getTops()[0]->P().Py(); // GeV
            
                double massOfSystemB =  newRsysVec.M(); // GeV
                double pxOfSystemB   =  newRsysVec.Px(); // GeV
                double pyOfSystemB   =  newRsysVec.Py(); // GeV

                deltaR = ttr.getTops()[0]->P().DeltaR(newRsysVec);

                if(ttr.getTops().size() >= 2)
                {
                    massOfSystemB = ttr.getTops()[1]->P().M(); // GeV 
                    pxOfSystemB   = ttr.getTops()[1]->P().Px(); // GeV
                    pyOfSystemB   = ttr.getTops()[1]->P().Py(); // GeV

                    deltaR = ttr.getTops()[0]->P().DeltaR(ttr.getTops()[1]->P());
                }
            
                // The missing transverse momentum:
                const double pxMiss        = metLVec.Px(); // GeV
                const double pyMiss        = metLVec.Py(); // GeV
            
                const double invis_mass    = metLVec.M(); // GeV
          
                Mt2::LorentzTransverseVector  vis_A(Mt2::TwoVector(pxOfSystemA, pyOfSystemA), massOfSystemA);
                Mt2::LorentzTransverseVector  vis_B(Mt2::TwoVector(pxOfSystemB, pyOfSystemB), massOfSystemB);
                Mt2::TwoVector                pT_Miss(pxMiss, pyMiss);
          
                newMT2 = mt2Calculator.mt2_332(vis_A, vis_B, pT_Miss, invis_mass);
            }

            double oldptt1 = (vTops->size())?((*vTops)[0].Pt()):(0);
            double newptt1 = (vTops->size())?((*vTopsNew)[0].Pt()):(0);

            tr.registerDerivedVar("deltaRt1t2New", deltaR);

            tr.registerDerivedVar("oldRsysVec_pt", oldRsysVec.Pt());
            tr.registerDerivedVar("newRsysVec_pt", newRsysVec.Pt());
            tr.registerDerivedVar("oldRsysVec_eta", oldRsysVec.Eta());
            tr.registerDerivedVar("newRsysVec_eta", newRsysVec.Eta());
            tr.registerDerivedVar("oldRsysVec_phi", oldRsysVec.Phi());
            tr.registerDerivedVar("newRsysVec_phi", newRsysVec.Phi());
            tr.registerDerivedVar("oldRsysVec_m", oldRsysVec.M());
            tr.registerDerivedVar("newRsysVec_m", newRsysVec.M());

            tr.registerDerivedVar("oldMT2", oldMT2);
            tr.registerDerivedVar("newMT2", newMT2);

            tr.registerDerivedVar("oldt1pt", oldptt1);
            tr.registerDerivedVar("newt1pt", newptt1);

            tr.registerDerivedVar("nSnbJetRsys", int(ttr.getRsys().getNConstituents()));

            tr.registerDerivedVar("nTaggerJets", int(jetsLVec_forTagger.size()));

            tr.registerDerivedVar("nTops", nTops);

            tr.registerDerivedVar("nTopsNew", int(ttr.getTops().size()));

            tr.registerDerivedVar("cntCSVS", cntCSVS);

            tr.registerDerivedVec("vTops", vTops);

            tr.registerDerivedVec("vTopsNew", vTopsNew);

            tr.registerDerivedVec("vdTops", vdTops);

            tr.registerDerivedVec("v3j_dr12", v3j_dr12);
            tr.registerDerivedVec("v3j_dr23", v3j_dr23);
            tr.registerDerivedVec("v3j_dr13", v3j_dr13);
            tr.registerDerivedVec("v3j_cm_dr12", v3j_cm_dr12);
            tr.registerDerivedVec("v3j_cm_dr13", v3j_cm_dr13);
            tr.registerDerivedVec("v3j_cm_dr23", v3j_cm_dr23);
            tr.registerDerivedVec("v3j_diJets", v3j_diJets);
        }

    public:
        PrepareTopVars() : tt(nullptr)
	{
            t3tagger.setnJetsSel(1);
            t3tagger.setCSVS(AnaConsts::cutCSVS);

            std::string cfgDocText =
                "TopTagger\n"
                "{\n"
                "    module[0] = \"TTMLazyClusterAlgo\""
                "    module[1] = \"TTMHEPRequirements\""
                "    module[2] = \"TTMOverlapResolution\""
                "}\n";

            tt = new TopTagger();
            tt->setCfgFile("TopTagger.cfg");
            //tt = ttf.makeTopTagger();
	}

        ~PrepareTopVars()
        {
            //if(tt) delete tt;
        }

	void operator()(NTupleReader& tr)
	{
	    prepareTopVars(tr);
	}

    };


    class MetSmear
    {
    private:

        TRandom3 *trand;

        double calcB(double x1, double x2, double e1, double e2)
        {
            return (x2*x2*(e1-1) - x1*x1*(e2-1)) / (x1*x2*(x2-x1));
        }

        double calcC(double x1, double x2, double e1, double e2)
        {
            return (x2*(e1-1) - x1*(e2-1)) / (x1*x2*(x1-x2));
        }

        double calcQuad(double x, double x1, double x2, double e1, double e2)
        {
            return 1 + x*calcB(x1, x2, e1, e2) + x*x*calcC(x1, x2, e1, e2);
        }

        double logistical(double met, double A, double B, double C)
        {
            return 1 - A/(1 + exp(-(B*(met - C))));
        }
        
        void metSmear(NTupleReader& tr)
        {
            const double& met     = tr.getVar<double>("cleanMetPt");
            
            // Logistical smearing 
            double met_logi_1 = met * logistical(met, 0.15, 0.01, 300);
            double met_logi_2 = met * logistical(met, 0.20, 0.01, 400);
            double met_logi_3 = met * logistical(met, 0.25, 0.01, 500);
            double met_logi_4 = met * logistical(met, 0.20, 0.01, 400);
            double met_logi_5 = met * logistical(met, 0.20, 0.02, 400);
            double met_logi_6 = met * logistical(met, 0.20, 0.03, 400);
            double met_logi_7 = met * logistical(met, 0.20, 0.02, 300);
            double met_logi_8 = met * logistical(met, 0.20, 0.02, 400);
            double met_logi_9 = met * logistical(met, 0.20, 0.02, 500);

            // gaussian smearing 
            //double met_gaus_5  = trand->Gaus(met, 5);
            //double met_gaus_10 = trand->Gaus(met, 10);
            //double met_gaus_15 = trand->Gaus(met, 15);
            double met_gaus_20 = trand->Gaus(met, 20);
            //double met_gaus_25 = trand->Gaus(met, 25);
            double met_gaus_30 = trand->Gaus(met, 30);
            double met_gaus_40 = trand->Gaus(met, 40);
            double met_gaus_50 = trand->Gaus(met, 50);

            tr.registerDerivedVar("met_gaus_20", met_gaus_20);
            //tr.registerDerivedVar("met_gaus_25", met_gaus_25);
            tr.registerDerivedVar("met_gaus_30", met_gaus_30);
            tr.registerDerivedVar("met_gaus_40", met_gaus_40);
            tr.registerDerivedVar("met_gaus_50", met_gaus_50);

            tr.registerDerivedVar("met_logi_1", met_logi_1);
            tr.registerDerivedVar("met_logi_2", met_logi_2);
            tr.registerDerivedVar("met_logi_3", met_logi_3);
            tr.registerDerivedVar("met_logi_4", met_logi_4);
            tr.registerDerivedVar("met_logi_5", met_logi_5);
            tr.registerDerivedVar("met_logi_6", met_logi_6);
            tr.registerDerivedVar("met_logi_7", met_logi_7);
            tr.registerDerivedVar("met_logi_8", met_logi_8);
            tr.registerDerivedVar("met_logi_9", met_logi_9);
        }

        void mt2Smear(NTupleReader& tr)
        {
            const double& metphi       = tr.getVar<double>("cleanMetPhi");
            const double& met_logi_1   = tr.getVar<double>("met_logi_1");
            const double& met_gaus_30  = tr.getVar<double>("met_gaus_30");
            
            const std::vector<TLorentzVector>& jetsLVec_forTagger  = tr.getVec<TLorentzVector>("jetsLVec_forTaggerZinv");
            const std::vector<double>&     recoJetsBtag_forTagger  = tr.getVec<double>("recoJetsBtag_forTaggerZinv");

            //We choose 30 GeV gaussian smearing and logi_1 for the study

            // Form TLorentzVector of MET
            TLorentzVector metLVec_Logi;
            metLVec_Logi.SetPtEtaPhiM(met_logi_1, 0, metphi, 0);
            
            type3Ptr->processEvent(jetsLVec_forTagger, recoJetsBtag_forTagger, metLVec_Logi);
            double MT2_Logi = type3Ptr->best_had_brJet_MT2;

            TLorentzVector metLVec_Gaus;
            metLVec_Gaus.SetPtEtaPhiM(met_gaus_30, 0, metphi, 0);
            
            type3Ptr->processEvent(jetsLVec_forTagger, recoJetsBtag_forTagger, metLVec_Gaus); 
            double MT2_Gaus = type3Ptr->best_had_brJet_MT2;

            tr.registerDerivedVar("mt2_logi_1",  MT2_Logi);
            tr.registerDerivedVar("mt2_gaus_30", MT2_Gaus);
        }

    public:
        MetSmear()
        {
            trand = new TRandom3(452147);
        }

        //~MetSmear()
        // {
        //    if(trand) delete trand;
        // }

        void operator()(NTupleReader& tr)
        {
            metSmear(tr);
            mt2Smear(tr);
        }
    };

    class PrepareMiniTupleVars
    {
    private:

        static const int BIT_PASSLEPTVETO             = 0x00000001;
        static const int BIT_PASSMUONVETO             = 0x00000002;
        static const int BIT_PASSELEVETO              = 0x00000004;
        static const int BIT_PASSISOTRKVETO           = 0x00000008;
        static const int BIT_PASSNJETS                = 0x00000010;
        static const int BIT_PASSDPHIS                = 0x00000020;
        static const int BIT_PASSBJETS                = 0x00000040;
        static const int BIT_PASSMET                  = 0x00000080;
        static const int BIT_PASSMT2                  = 0x00000100;
        static const int BIT_PASSHT                   = 0x00000200;
        static const int BIT_PASSTAGGER               = 0x00000400;
        static const int BIT_PASSNOISEEVENTFILTER     = 0x00000800;
        static const int BIT_PASSBASELINE             = 0x00001000;
        static const int BIT_PASSBASELINENOTAGMT2     = 0x00002000;
        static const int BIT_PASSBASELINENOTAG        = 0x00004000;
        static const int BIT_PASSLEPTVETOZINV         = 0x00008000;
        static const int BIT_PASSMUONVETOZINV         = 0x00010000;
        static const int BIT_PASSELEVETOZINV          = 0x00020000;
        static const int BIT_PASSISOTRKVETOZINV       = 0x00040000;
        static const int BIT_PASSNJETSZINV            = 0x00080000;
        static const int BIT_PASSDPHISZINV            = 0x00100000;
        static const int BIT_PASSBJETSZINV            = 0x00200000;
        static const int BIT_PASSMETZINV              = 0x00400000;
        static const int BIT_PASSMT2ZINV              = 0x00800000;
        static const int BIT_PASSHTZINV               = 0x01000000;
        static const int BIT_PASSTAGGERZINV           = 0x02000000;
        static const int BIT_PASSNOISEEVENTFILTERZINV = 0x04000000;
        static const int BIT_PASSBASELINEZINV         = 0x08000000;
        static const int BIT_PASSBASELINENOTAGMT2ZINV = 0x10000000;
        static const int BIT_PASSBASELINENOTAGZINV    = 0x20000000;
        static const int BIT_PASSMUZINVSEL            = 0x40000000;
        static const int BIT_PASSELMUZINVSEL          = 0x80000000;

        bool pack_;

        void pack(NTupleReader& tr)
        {
            const bool& passLeptVeto =         tr.getVar<bool>("passLeptVeto");
            const bool& passMuonVeto =         tr.getVar<bool>("passMuonVeto");
            const bool& passEleVeto =          tr.getVar<bool>("passEleVeto");
            const bool& passIsoTrkVeto =       tr.getVar<bool>("passIsoTrkVeto");
            const bool& passnJets =            tr.getVar<bool>("passnJets");
            const bool& passdPhis =            tr.getVar<bool>("passdPhis");
            const bool& passBJets =            tr.getVar<bool>("passBJets");
            const bool& passMET =              tr.getVar<bool>("passMET");
            const bool& passMT2 =              tr.getVar<bool>("passMT2");
            const bool& passHT =               tr.getVar<bool>("passHT");
            const bool& passTagger =           tr.getVar<bool>("passTagger");
            const bool& passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilter");
            const bool& passBaseline =         tr.getVar<bool>("passBaseline");
            const bool& passBaselineNoTagMT2 = tr.getVar<bool>("passBaselineNoTagMT2");
            const bool& passBaselineNoTag =    tr.getVar<bool>("passBaselineNoTag");

            const bool& passLeptVetoZinv =         tr.getVar<bool>("passLeptVetoZinv");
            const bool& passMuonVetoZinv =         tr.getVar<bool>("passMuonVetoZinv");
            const bool& passEleVetoZinv =          tr.getVar<bool>("passEleVetoZinv");
            const bool& passIsoTrkVetoZinv =       tr.getVar<bool>("passIsoTrkVetoZinv");
            const bool& passnJetsZinv =            tr.getVar<bool>("passnJetsZinv");
            const bool& passdPhisZinv =            tr.getVar<bool>("passdPhisZinv");
            const bool& passBJetsZinv =            tr.getVar<bool>("passBJetsZinv");
            const bool& passMETZinv =              tr.getVar<bool>("passMETZinv");
            const bool& passMT2Zinv =              tr.getVar<bool>("passMT2Zinv");
            const bool& passHTZinv =               tr.getVar<bool>("passHTZinv");
            const bool& passTaggerZinv =           tr.getVar<bool>("passTaggerZinv");
            const bool& passNoiseEventFilterZinv = tr.getVar<bool>("passNoiseEventFilterZinv");
            const bool& passBaselineZinv =         tr.getVar<bool>("passBaselineZinv");
            const bool& passBaselineNoTagMT2Zinv = tr.getVar<bool>("passBaselineNoTagMT2Zinv");
            const bool& passBaselineNoTagZinv =    tr.getVar<bool>("passBaselineNoTagZinv");

            const bool& passMuZinvSel =            tr.getVar<bool>("passMuZinvSel");
            const bool& passElMuZinvSel =          tr.getVar<bool>("passElMuZinvSel");

            int cuts = 0;

            if(passLeptVeto)         cuts |= BIT_PASSLEPTVETO;
            if(passMuonVeto)         cuts |= BIT_PASSMUONVETO;
            if(passEleVeto)          cuts |= BIT_PASSELEVETO;
            if(passIsoTrkVeto)       cuts |= BIT_PASSISOTRKVETO;
            if(passnJets)            cuts |= BIT_PASSNJETS;
            if(passdPhis)            cuts |= BIT_PASSDPHIS;
            if(passBJets)            cuts |= BIT_PASSBJETS;
            if(passMET)              cuts |= BIT_PASSMET;
            if(passMT2)              cuts |= BIT_PASSMT2;
            if(passHT)               cuts |= BIT_PASSHT;
            if(passTagger)           cuts |= BIT_PASSTAGGER;
            if(passNoiseEventFilter) cuts |= BIT_PASSNOISEEVENTFILTER;
            if(passBaseline)         cuts |= BIT_PASSBASELINE;
            if(passBaselineNoTagMT2) cuts |= BIT_PASSBASELINENOTAGMT2;
            if(passBaselineNoTag)    cuts |= BIT_PASSBASELINENOTAG;

            if(passLeptVetoZinv)         cuts |= BIT_PASSLEPTVETOZINV;
            if(passMuonVetoZinv)         cuts |= BIT_PASSMUONVETOZINV;
            if(passEleVetoZinv)          cuts |= BIT_PASSELEVETOZINV;
            if(passIsoTrkVetoZinv)       cuts |= BIT_PASSISOTRKVETOZINV;
            if(passnJetsZinv)            cuts |= BIT_PASSNJETSZINV;
            if(passdPhisZinv)            cuts |= BIT_PASSDPHISZINV;
            if(passBJetsZinv)            cuts |= BIT_PASSBJETSZINV;
            if(passMETZinv)              cuts |= BIT_PASSMETZINV;
            if(passMT2Zinv)              cuts |= BIT_PASSMT2ZINV;
            if(passHTZinv)               cuts |= BIT_PASSHTZINV;
            if(passTaggerZinv)           cuts |= BIT_PASSTAGGERZINV;
            if(passNoiseEventFilterZinv) cuts |= BIT_PASSNOISEEVENTFILTERZINV;
            if(passBaselineZinv)         cuts |= BIT_PASSBASELINEZINV;
            if(passBaselineNoTagMT2Zinv) cuts |= BIT_PASSBASELINENOTAGMT2ZINV;
            if(passBaselineNoTagZinv)    cuts |= BIT_PASSBASELINENOTAGZINV;

            if(passMuZinvSel)            cuts |= BIT_PASSMUZINVSEL;
            if(passElMuZinvSel)          cuts |= BIT_PASSELMUZINVSEL;

            tr.registerDerivedVar("cuts", cuts);
        }

        void unpack(NTupleReader& tr)
        {
            const int& cuts = tr.getVar<int>("cuts");

            tr.registerDerivedVar("passLeptVeto",         static_cast<bool>(cuts & BIT_PASSLEPTVETO));
            tr.registerDerivedVar("passMuonVeto",         static_cast<bool>(cuts & BIT_PASSMUONVETO));
            tr.registerDerivedVar("passEleVeto",          static_cast<bool>(cuts & BIT_PASSELEVETO));
            tr.registerDerivedVar("passIsoTrkVeto",       static_cast<bool>(cuts & BIT_PASSISOTRKVETO));
            tr.registerDerivedVar("passnJets",            static_cast<bool>(cuts & BIT_PASSNJETS));
            tr.registerDerivedVar("passdPhis",            static_cast<bool>(cuts & BIT_PASSDPHIS));
            tr.registerDerivedVar("passBJets",            static_cast<bool>(cuts & BIT_PASSBJETS));
            tr.registerDerivedVar("passMET",              static_cast<bool>(cuts & BIT_PASSMET));
            tr.registerDerivedVar("passMT2",              static_cast<bool>(cuts & BIT_PASSMT2));
            tr.registerDerivedVar("passHT",               static_cast<bool>(cuts & BIT_PASSHT));
            tr.registerDerivedVar("passTagger",           static_cast<bool>(cuts & BIT_PASSTAGGER));
            tr.registerDerivedVar("passNoiseEventFilter", static_cast<bool>(cuts & BIT_PASSNOISEEVENTFILTER));
            tr.registerDerivedVar("passBaseline",         static_cast<bool>(cuts & BIT_PASSBASELINE));
            tr.registerDerivedVar("passBaselineNoTagMT2", static_cast<bool>(cuts & BIT_PASSBASELINENOTAGMT2));
            tr.registerDerivedVar("passBaselineNoTag",    static_cast<bool>(cuts & BIT_PASSBASELINENOTAG));

            tr.registerDerivedVar("passLeptVetoZinv",         static_cast<bool>(cuts & BIT_PASSLEPTVETOZINV));
            tr.registerDerivedVar("passMuonVetoZinv",         static_cast<bool>(cuts & BIT_PASSMUONVETOZINV));
            tr.registerDerivedVar("passEleVetoZinv",          static_cast<bool>(cuts & BIT_PASSELEVETOZINV));
            tr.registerDerivedVar("passIsoTrkVetoZinv",       static_cast<bool>(cuts & BIT_PASSISOTRKVETOZINV));
            tr.registerDerivedVar("passnJetsZinv",            static_cast<bool>(cuts & BIT_PASSNJETSZINV));
            tr.registerDerivedVar("passdPhisZinv",            static_cast<bool>(cuts & BIT_PASSDPHISZINV));
            tr.registerDerivedVar("passBJetsZinv",            static_cast<bool>(cuts & BIT_PASSBJETSZINV));
            tr.registerDerivedVar("passMETZinv",              static_cast<bool>(cuts & BIT_PASSMETZINV));
            tr.registerDerivedVar("passMT2Zinv",              static_cast<bool>(cuts & BIT_PASSMT2ZINV));
            tr.registerDerivedVar("passHTZinv",               static_cast<bool>(cuts & BIT_PASSHTZINV));
            tr.registerDerivedVar("passTaggerZinv",           static_cast<bool>(cuts & BIT_PASSTAGGERZINV));
            tr.registerDerivedVar("passNoiseEventFilterZinv", static_cast<bool>(cuts & BIT_PASSNOISEEVENTFILTERZINV));
            tr.registerDerivedVar("passBaselineZinv",         static_cast<bool>(cuts & BIT_PASSBASELINEZINV));
            tr.registerDerivedVar("passBaselineNoTagMT2Zinv", static_cast<bool>(cuts & BIT_PASSBASELINENOTAGMT2ZINV));
            tr.registerDerivedVar("passBaselineNoTagZinv",    static_cast<bool>(cuts & BIT_PASSBASELINENOTAGZINV));

            tr.registerDerivedVar("passMuZinvSel",   static_cast<bool>(cuts & BIT_PASSMUZINVSEL));
            tr.registerDerivedVar("passElMuZinvSel", static_cast<bool>(cuts & BIT_PASSELMUZINVSEL));
        }

    public:
        PrepareMiniTupleVars(bool pack)
        {
            pack_ = pack;
        }

        void operator()(NTupleReader& tr)
        {
            if(pack_) pack(tr);
            else      unpack(tr);
        }
    };

    //void printInterestingEvents(NTupleReader& tr)
    //{
    //    const unsigned int& run   = tr.getVar<unsigned int>("run");
    //    const unsigned int& event = tr.getVar<unsigned int>("event");
    //
    //    const double& met                            = tr.getVar<double>("met");
    //    const double& metphi                         = tr.getVar<double>("metphi");
    //
    //    const int& nMuons_CUT        = tr.getVar<int>("nMuons_CUT");
    //    const int& nElectrons_CUT    = tr.getVar<int>("nElectrons_CUT");
    //    const int& cntNJetsPt50Eta24 = tr.getVar<int>("cntNJetsPt50Eta24");
    //    const int& cntNJetsPt30Eta24 = tr.getVar<int>("cntNJetsPt30Eta24");
    //    const int& cntNJetsPt30      = tr.getVar<int>("cntNJetsPt30");
    //
    //    const double& mht    = tr.getVar<double>("mht");
    //    const double& mhtphi = tr.getVar<double>("mhtphi");
    //    const double& ht     = tr.getVar<double>("ht");
    //
    //    //if(met > 1000) std::cout << "run: " << run << "\tevent: " << event << "\tmet: " << met << "\tmetphi: " << metphi << "\tnMuons_CUT: " << nMuons_CUT << "\t nElectrons_CUT: " << nElectrons_CUT << "\tcntNJetsPt30: " << cntNJetsPt30 << "\tcntNJetsPt30Eta24: " << cntNJetsPt30Eta24 << "\tcntNJetsPt50Eta24: " << cntNJetsPt50Eta24 << "\tmht: " << mht << "\tmhtphi: " << mhtphi << "\tht: " << ht << std::endl;
    //}

}

#endif
