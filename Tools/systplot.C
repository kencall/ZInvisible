#include "tdrstyle.h"
#include "ScaleFactors.h"
#include "searchBins.h"

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

void smartMax(const TH1* const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax)
{
    const bool isLog = p->GetLogy();
    double min = 9e99;
    double max = -9e99;
    double pThreshMax = -9e99;
    int threshold = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        double bin = h->GetBinContent(i);
        if(bin > max) max = bin;
        else if(bin > 1e-10 && bin < min) min = bin;
        if(i >= threshold && bin > pThreshMax) pThreshMax = bin;
    }

    gpThreshMax = std::max(gpThreshMax, pThreshMax);
    gmax = std::max(gmax, max);
    gmin = std::min(gmin, min);
}

void makeplot(TFile* f, std::string hname, double xlow, double xhigh, double ylow, double yhigh, std::string prefix, bool log = false)
{

    TH1D* h = (TH1D*)f->Get(hname.c_str());

    // Set the style
    setTDRStyle();
    gStyle->SetErrorX(0.5);
    // Prepare canvas
    TCanvas *c;
    double fontScale;
    c = new TCanvas("c1", "c1", 1200, 800);
    c->Divide(1, 1);
    c->cd(1);
    gPad->SetPad("p1", "p1", 0, 0, 1, 1, kWhite, 0, 0);
    gPad->SetBottomMargin(0.12);
    fontScale = 6.5 / 8;
    
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.07 * (8.0 / 6.5) * fontScale);
        
    // Use dummy histogram to set the axes and stuff 
    TH1 *dummy = new TH1F("dummy", "dummy", 1000, h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX()));
    dummy->SetStats(0);
    dummy->SetTitle(0);
    dummy->GetYaxis()->SetTitle("Systematic uncertainty");
    dummy->GetYaxis()->SetTitleOffset(.9*1.05 / (fontScale));
    dummy->GetXaxis()->SetTitleOffset(1.05);
    dummy->GetXaxis()->SetTitle("Search Bin");
    dummy->GetXaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetXaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    if(dummy->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);

    double max = 0.0, lmax = 0.0, min = 1.0e300, minAvgWgt = 1.0e300;
    int iSingle = 0, iRatio = 0;
    char legEntry[128];

    // Format the histogram
    h->SetLineColor(kBlue+2);
    h->SetLineWidth(3);
    iSingle++;
    
    dummy->GetXaxis()->SetRangeUser(xlow, xhigh);
    dummy->GetYaxis()->SetRangeUser(ylow, yhigh);
    dummy->Draw();
    
    h->Draw("hist same");

    fixOverlay();

    c->cd(1);
    char lumistamp[128];
    sprintf(lumistamp, "%.1f fb^{-1} at 13 TeV", 2156. / 1000.0);
    TLatex mark;
    mark.SetTextSize(0.042 * fontScale);
    mark.SetTextFont(42);
    mark.SetNDC(true);
    mark.DrawLatex(gPad->GetLeftMargin(), 0.95, "CMS Preliminary");
    mark.SetTextAlign(31);
    mark.DrawLatex(1 - gPad->GetRightMargin(), 0.95, lumistamp);
    
    fixOverlay();
    drawSBregionDef(dummy->GetMinimum(),dummy->GetMaximum(),false);

    char outname[128];
    sprintf(outname, (prefix + "_%s.pdf").c_str(), hname.c_str());
    c->Print(outname);
    sprintf(outname, (prefix + "_%s.png").c_str(), hname.c_str());
    c->Print(outname);

    c->Close();
}

void makeplotratio(TFile* f, std::vector<std::string> hnames, std::vector<std::string> hnames_ratio, double xlow, double xhigh, double ylow, double yhigh, std::string prefix, double ratiomax = 0.08)
{

    std::vector<TH1D*> histos;
    for(std::string hname : hnames)
    {
	histos.push_back((TH1D*)f->Get(hname.c_str()));
    }
    std::vector<TH1D*> histosratio;
    for(std::string hname : hnames_ratio)
    {
	histosratio.push_back((TH1D*)f->Get(hname.c_str()));
    }

    // Set the style
    setTDRStyle();
    gStyle->SetErrorX(0.5);
    // Prepare canvas
    TCanvas *c;
    double fontScale;
    c = new TCanvas("c1", "c1", 1200, 800);
    c->Divide(1, 2);
    c->cd(1);
    gPad->SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, kWhite, 0, 0);
    gPad->SetBottomMargin(0.01);
    fontScale = 1.0;
    
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.07 * (8.0 / 6.5) * fontScale);
        
    // Use dummy histogram to set the axes and stuff 
    TH1 *dummy = new TH1F("dummy", "dummy", 1000, histos[0]->GetBinLowEdge(1), histos[0]->GetBinLowEdge(histos[0]->GetNbinsX()) + histos[0]->GetBinWidth(histos[0]->GetNbinsX()));
    dummy->SetStats(0);
    dummy->SetTitle(0);
    dummy->GetYaxis()->SetTitle("A.U.");
    dummy->GetYaxis()->SetTitleOffset(.9*1.05 / (fontScale));
    dummy->GetXaxis()->SetTitleOffset(1.05);
    dummy->GetXaxis()->SetTitle("");
    dummy->GetXaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetXaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetTitleSize(0.20 * 2 / 6.5 * fontScale);
    dummy->GetYaxis()->SetLabelSize(0.20 * 2 / 6.5 * fontScale);
    if(dummy->GetNdivisions() % 100 > 5) dummy->GetXaxis()->SetNdivisions(6, 5, 0);

    double max = 0.0, lmax = 0.0, min = 1.0e300, minAvgWgt = 1.0e300;
    int iSingle = 0, iRatio = 0;
    char legEntry[128];

    TLegend *leg = new TLegend(0.60, 0.88 - 3 * 2 * 0.045, 0.89, 0.88);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);

    // Format the histograms
    histos[0]->SetLineColor(kBlue+2);
    histos[1]->SetLineColor(kBlack);
    histos[2]->SetLineColor(kRed+2);
    for(TH1D* h : histos)
	h->SetLineWidth(3);
    iSingle++;
    leg->AddEntry(histos[0],"Up variation","l");
    leg->AddEntry(histos[1],"Nominal variation","l");
    leg->AddEntry(histos[2],"Down variation","l");
    
    dummy->GetXaxis()->SetRangeUser(xlow, xhigh);
    dummy->GetYaxis()->SetRangeUser(ylow, yhigh);
    dummy->Draw();
    
    for(TH1D* h : histos)
	h->Draw("hist same");
    leg->Draw();
    gPad->SetLogy(true);

    fixOverlay();

    c->cd(1);
    char lumistamp[128];
    sprintf(lumistamp, "%.1f fb^{-1} at 13 TeV", 2156. / 1000.0);
    TLatex mark;
    mark.SetTextSize(0.042 * fontScale);
    mark.SetTextFont(42);
    mark.SetNDC(true);
    mark.DrawLatex(gPad->GetLeftMargin(), 0.95, "CMS Preliminary");
    mark.SetTextAlign(31);
    mark.DrawLatex(1 - gPad->GetRightMargin(), 0.95, lumistamp);
    
    fixOverlay();
    drawSBregionDef(dummy->GetMinimum(),dummy->GetMaximum(),true);

    c->cd(2);
    gPad->SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, kWhite, 0, 0);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.37);
    gPad->SetLogy(false);

    TH1* dummy2 = new TH1F("dummy2", "dummy2", 1000, histosratio[0]->GetBinLowEdge(1), histosratio[0]->GetBinLowEdge(histosratio[0]->GetNbinsX()) + histosratio[0]->GetBinWidth(histosratio[0]->GetNbinsX()));
    dummy2->GetXaxis()->SetTitle("Search Bin");
    dummy2->GetXaxis()->SetTitleOffset(1.05);
    dummy2->GetYaxis()->SetTitle("#frac{syst - nominal}{nominal}");
    dummy2->GetYaxis()->SetTitleOffset(0.42);
    dummy2->GetYaxis()->SetNdivisions(3, 5, 0, true);

    dummy2->GetYaxis()->SetTitleSize(0.16 * 2 / 2.5);
    dummy2->GetYaxis()->SetLabelSize(0.20 * 2 / 2.5);
    dummy2->GetXaxis()->SetTitleSize(0.20 * 2 / 2.5);
    dummy2->GetXaxis()->SetLabelSize(0.20 * 2 / 2.5);

    dummy2->SetStats(0);
    dummy2->SetTitle(0);
    if(dummy2->GetNdivisions() % 100 > 5) dummy2->GetXaxis()->SetNdivisions(6, 5, 0);

    histosratio[0]->SetMarkerColor(kBlue+2);
    histosratio[1]->SetMarkerColor(kRed+2);
    for(TH1D* h : histosratio)
    {
	h->SetMarkerStyle(8);
	h->SetMarkerSize(1);
    }
    dummy2->GetXaxis()->SetRangeUser(xlow, xhigh);
    dummy2->GetYaxis()->SetRangeUser(-0.07, ratiomax);
    dummy2->Draw();
    for(TH1D* h : histosratio)
	h->Draw("histPsame");

    c->cd();

    char outname[128];
    sprintf(outname, (prefix + "_%s.pdf").c_str(), hnames[0].c_str());
    c->Print(outname);
    sprintf(outname, (prefix + "_%s.png").c_str(), hnames[0].c_str());
    c->Print(outname);

    c->Close();
}

int main(int argc, char* argv[])
{
    // // Get the relevant information
    // TFile* f1 = TFile::Open("/uscms_data/d3/nstrobbe/HadronicStop/DataTest/CMSSW_7_4_8/src/ZInvisible/Tools/dataMCweights.root");
    // std::vector<std::string> hnames = {"DataMC_nj_elmuZinv_ht200_dphi",
    // 				       "DataMC_nj_elmuZinv_0b_ht200_dphi",
    // 				       "DataMC_nj_elmuZinv_g1b_ht200_dphi",
    // 				       "DataMC_nj_muZinv_ht200_dphi",
    // 				       "DataMC_nj_muZinv_0b_ht200_dphi",
    // 				       "DataMC_nj_muZinv_g1b_ht200_dphi"};
    // for(std::string const& hname : hnames)
    // {
    // 	makeplot(f1, hname, 4, 20, 0.0, 2.0, "SF_njet");
    // }

    // f1->Close();

    TFile* f1 = TFile::Open("/uscms_data/d3/pastika/zinv/dev/CMSSW_7_4_8/src/ZInvisible/Tools/syst_all.root");
    std::vector<std::string> hnames = {"shape_central",
				       "shape_stat",
				       "MC_stats"};

    for(std::string const& hname : hnames)
    {
        makeplot(f1, hname, 0, 45, 0.0, 1.18, "Syst");
    }

    f1->Close();


    f1 = TFile::Open("/uscms_data/d3/nstrobbe/HadronicStop/DataTest/CMSSW_7_4_8/src/ZInvisible/Tools/syst_scalePDF.root");
    std::vector<std::string> hnames1 = {"nSearchBin_scale_up",
					"nSearchBin_scale_nominal",
					"nSearchBin_scale_down"};
    std::vector<std::string> hnamesratio1 = {"nSearchBin_ratio_scale_up",
					     "nSearchBin_ratio_scale_down"};
    std::vector<std::string> hnames2 = {"nSearchBin_pdf_up",
					"nSearchBin_pdf_nominal",
					"nSearchBin_pdf_down"};
    std::vector<std::string> hnamesratio2 = {"nSearchBin_pdf_up",
					     "nSearchBin_ratio_pdf_down"};

    makeplotratio(f1, hnames1, hnamesratio1, 0, 45, 0.00002, 0.3, "Syst");
    makeplotratio(f1, hnames2, hnamesratio2, 0, 45, 0.00002, 0.3, "Syst",0.25);


    f1->Close();


    return 0;

}