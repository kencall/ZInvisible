#include "RegisterFunctions.h"
#include "NTupleReader.h"
#include "derivedTupleVariables.h"
#include "baselineDef.h"
#include "Systematic.h"
#include "PDFUncertainty.h"
#include "BTagCorrector.h"

const std::set<std::string> RegisterFunctions::getMiniTupleSet()
{
    return std::set<std::string>({"HTZinv","cleanMetPt","cleanMetPhi","best_had_brJet_MT2Zinv","cntCSVSZinv","nTopCandSortedCntZinv","cntNJetsPt30Eta24Zinv","nSearchBin","cutMuVec","cutElecVec","jetsLVec_forTaggerZinv", "recoJetsBtag_forTaggerZinv","zEffWgt","zAccWgt","cuts","passMuTrigger","genHT","genWeight","bTagSF_EventWeightSimple_Central"});
    //return std::set<std::string>({"HTZinv","cleanMetPt","cleanMetPhi","best_had_brJet_MT2Zinv","cntCSVSZinv","nTopCandSortedCntZinv","cntNJetsPt30Eta24Zinv","nSearchBin","cutMuVec","cutElecVec","jetsLVec_forTaggerZinv", "recoJetsBtag_forTaggerZinv","zEffWgt","zAccWgt","cuts","passMuTrigger","genHT","genWeight"});
}

const std::set<std::string> RegisterFunctions::getMiniTupleSetData()
{
    return std::set<std::string>({"HTZinv","cleanMetPt","cleanMetPhi","best_had_brJet_MT2Zinv","cntCSVSZinv","nTopCandSortedCntZinv","cntNJetsPt30Eta24Zinv","nSearchBin","cutMuVec","cutElecVec","jetsLVec_forTaggerZinv", "recoJetsBtag_forTaggerZinv","zEffWgt","zAccWgt","cuts","passMuTrigger"});
}

void activateBranches(std::set<std::string>& activeBranches)
{
    for(auto& bn : AnaConsts::activatedBranchNames) activeBranches.insert(bn);
    for(auto& bn : AnaConsts::activatedBranchNames_DataOnly) activeBranches.insert(bn);
    activeBranches.insert("ht");
    activeBranches.insert("run");
    activeBranches.insert("lumi");
    activeBranches.insert("event");
    activeBranches.insert("mht");
    activeBranches.insert("mhtphi");
    activeBranches.insert("genDecayPdgIdVec");
    activeBranches.insert("genDecayLVec");
    activeBranches.insert("muonsLVec");
    activeBranches.insert("muonsRelIso");
    activeBranches.insert("muonsMiniIso");
    activeBranches.insert("W_emuVec");
    activeBranches.insert("muonsCharge");
    activeBranches.insert("muonsMtw");
    activeBranches.insert("met");
    activeBranches.insert("metphi");
    activeBranches.insert("jetsLVec");
    activeBranches.insert("elesLVec");
    activeBranches.insert("elesRelIso");
    activeBranches.insert("recoJetsBtag_0");
    activeBranches.insert("loose_isoTrks_mtw");
    activeBranches.insert("elesMtw");
    activeBranches.insert("loose_isoTrks_iso");
    activeBranches.insert("loose_isoTrksLVec");
    activeBranches.insert("muMatchedJetIdx");
    activeBranches.insert("eleMatchedJetIdx");
    activeBranches.insert("recoJetschargedEmEnergyFraction");
    activeBranches.insert("recoJetsneutralEmEnergyFraction");
    activeBranches.insert("recoJetschargedHadronEnergyFraction");
    activeBranches.insert("prodJetsNoMu_recoJetschargedEmEnergyFraction");
    activeBranches.insert("prodJetsNoMu_recoJetsneutralEmEnergyFraction");
    activeBranches.insert("prodJetsNoMu_recoJetschargedHadronEnergyFraction");
    activeBranches.insert("elesisEB");
    activeBranches.insert("elesMiniIso");
    activeBranches.insert("elesCharge");
    activeBranches.insert("calomet");
}

////////////////////////////////

RegisterFunctionsNTuple::RegisterFunctionsNTuple(bool isCondor, std::string sbEra) : RegisterFunctions()
{            
    AnaFunctions::prepareTopTagger();

    myBLV     = new BaselineVessel("", "");
    blvZinv   = new BaselineVessel("Zinv");
    blvZinv1b = new BaselineVessel("Zinv1b");
    blvZinv2b = new BaselineVessel("Zinv2b");
    blvZinv3b = new BaselineVessel("Zinv3b");
    blvZinvJEUUp = new BaselineVessel("ZinvJEUUp");
    blvZinvJEUDn = new BaselineVessel("ZinvJEUDn");
    blvZinvMEUUp = new BaselineVessel("ZinvMEUUp");
    blvZinvMEUDn = new BaselineVessel("ZinvMEUDn");

    weights              = new plotterFunctions::GenerateWeight;
    njWeight             = new plotterFunctions::NJetWeight;
    lepInfo              = new plotterFunctions::LepInfo;
    fakebtagvectors      = new plotterFunctions::Fakebtagvectors;
    getSearchBin         = new plotterFunctions::GetSearchBin(sbEra);
    triggerInfo          = new plotterFunctions::TriggerInfo;
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars(true);
    systematicPrep       = new plotterFunctions::SystematicPrep;
    systematicCalc       = new plotterFunctions::SystematicCalc(sbEra);

    myPDFUnc = new PDFUncertainty();
    bTagCorrector = nullptr;
    if(isCondor)
    {
        bTagCorrector = new BTagCorrector("bTagEffHists.root", "", false);
    }
    else
    {
        bTagCorrector = new BTagCorrector("bTagEffHists.root", "/uscms/home/pastika/nobackup/zinv/dev/CMSSW_7_4_8/src/SusyAnaTools/Tools/", false);
        //bTagCorrector = new BTagCorrector("bTagEffHists.root", "/uscms_data/d3/nstrobbe/HadronicStop/DataTest/CMSSW_7_4_8/src/SusyAnaTools/Tools/", false);
    }
}

RegisterFunctionsNTuple::~RegisterFunctionsNTuple()
{
    if(myBLV) delete myBLV;
    if(blvZinv) delete blvZinv;
    if(blvZinv1b) delete blvZinv1b;
    if(blvZinv2b) delete blvZinv2b;
    if(blvZinv3b) delete blvZinv3b;
    if(blvZinvJEUUp) delete blvZinvJEUUp;
    if(blvZinvJEUDn) delete blvZinvJEUDn;
    if(blvZinvMEUUp) delete blvZinvMEUUp;
    if(blvZinvMEUDn) delete blvZinvMEUDn;
    if(weights) delete weights;
    if(njWeight) delete njWeight;
    if(lepInfo) delete lepInfo;
    if(fakebtagvectors) delete fakebtagvectors;
    if(getSearchBin) delete getSearchBin;
    if(triggerInfo) delete triggerInfo;
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
    if(myPDFUnc) delete myPDFUnc;
    if(systematicPrep) delete systematicPrep;
    if(systematicCalc) delete systematicCalc;
    if(bTagCorrector) delete bTagCorrector;
}
        
void RegisterFunctionsNTuple::registerFunctions(NTupleReader& tr)
{
    //Make some global "constants" here

    //register functions with NTupleReader
    tr.registerFunction(*myBLV);
    tr.registerFunction(*lepInfo);
    tr.registerFunction(*weights);
    tr.registerFunction(*blvZinv);
    tr.registerFunction(*njWeight);
    tr.registerFunction(*fakebtagvectors);
    tr.registerFunction(*blvZinv1b);
    tr.registerFunction(*blvZinv2b);
    tr.registerFunction(*blvZinv3b);
    tr.registerFunction(*systematicPrep);
    tr.registerFunction(*blvZinvJEUUp);
    tr.registerFunction(*blvZinvJEUDn);
    tr.registerFunction(*blvZinvMEUUp);
    tr.registerFunction(*blvZinvMEUDn);
    tr.registerFunction(*getSearchBin);
    tr.registerFunction(*systematicCalc);
    tr.registerFunction(*triggerInfo);
    tr.registerFunction(*prepareMiniTupleVars);
    //tr.registerFunction(&printInterestingEvents);
    tr.registerFunction(*myPDFUnc);
    tr.registerFunction(*bTagCorrector);
}

void RegisterFunctionsNTuple::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}

void RegisterFunctionsNTuple::remakeBTagCorrector(std::string sampleName)
{
  if(sampleName.find("Data") == std::string::npos){
    if(bTagCorrector) bTagCorrector->resetEffs(sampleName);
  }
}

////////////////////////////////

RegisterFunctionsMiniTuple::RegisterFunctionsMiniTuple() : RegisterFunctions()
{            
    AnaFunctions::prepareTopTagger();

    njWeight             = new plotterFunctions::NJetWeight;
    triggerInfo          = new plotterFunctions::TriggerInfo(true);
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars(false);
    metSmear             = new plotterFunctions::MetSmear;
}

RegisterFunctionsMiniTuple::~RegisterFunctionsMiniTuple()
{
    if(njWeight) delete njWeight;
    if(triggerInfo) delete triggerInfo;
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
    if(metSmear) delete metSmear;
}
        
void RegisterFunctionsMiniTuple::registerFunctions(NTupleReader& tr)
{
    //Make some global "constants" here

    //register functions with NTupleReader
    tr.registerFunction(*njWeight);
    tr.registerFunction(*triggerInfo);
    tr.registerFunction(*prepareMiniTupleVars);
    tr.registerFunction(*metSmear);
}

void RegisterFunctionsMiniTuple::activateBranches(std::set<std::string>& activeBranches)
{
}

////////////////////////////////

RegisterFunctionsCalcEff::RegisterFunctionsCalcEff() : RegisterFunctions()
{
    AnaFunctions::prepareTopTagger();

    myBLV     = new BaselineVessel;
    lepInfo   = new plotterFunctions::LepInfo;
}

RegisterFunctionsCalcEff::~RegisterFunctionsCalcEff()
{
    if(myBLV) delete myBLV;
    if(lepInfo) delete lepInfo;
}

void RegisterFunctionsCalcEff::registerFunctions(NTupleReader& tr)
{
    //Make some global "constants" here

    //register functions with NTupleReader
    tr.registerFunction(*myBLV);
    tr.registerFunction(*lepInfo);
}

void RegisterFunctionsCalcEff::activateBranches(std::set<std::string>& activeBranches)
{
    ::activateBranches(activeBranches);
}

////////////////////////////////

void RegisterFunctionsSyst::addFunction(std::function<void(NTupleReader&)> func)
{
    funcs_.emplace_back(func);
}

RegisterFunctionsSyst::RegisterFunctionsSyst() : RegisterFunctions()
{
    njWeight             = new plotterFunctions::NJetWeight;
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars(false);
    triggerInfo          = new plotterFunctions::TriggerInfo(true);
    systWeights          = new SystWeights;
}

RegisterFunctionsSyst::~RegisterFunctionsSyst()
{
    if(njWeight) delete njWeight;
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
    if(triggerInfo) delete triggerInfo;
    if(systWeights) delete systWeights;
}

void RegisterFunctionsSyst::registerFunctions(NTupleReader& tr)
{
    for(auto& f : funcs_) tr.registerFunction(f);

    tr.registerFunction(*njWeight);
    tr.registerFunction(*prepareMiniTupleVars);
    tr.registerFunction(*triggerInfo);
    //tr.registerFunction(*systWeights);
}

////////////////////////////////

RegisterFunctions2Dplot::RegisterFunctions2Dplot()
{
    prepareMiniTupleVars = new plotterFunctions::PrepareMiniTupleVars(false);
}

RegisterFunctions2Dplot::~RegisterFunctions2Dplot()
{
    if(prepareMiniTupleVars) delete prepareMiniTupleVars;
}

void RegisterFunctions2Dplot::registerFunctions(NTupleReader& tr)
{
    tr.registerFunction(*prepareMiniTupleVars);
}

////////////////////////////////

RegisterFunctionsTopStudy::RegisterFunctionsTopStudy()
{
    prepareTopVars = new plotterFunctions::PrepareTopVars();
}

RegisterFunctionsTopStudy::~RegisterFunctionsTopStudy()
{
    if(prepareTopVars) delete prepareTopVars;
}

void RegisterFunctionsTopStudy::registerFunctions(NTupleReader& tr)
{
    tr.registerFunction(*prepareTopVars);
}

/////////////////////////////////

void drawSBregionDefCopy(const double ymin_Yields, const double ymax_Yields, const bool logscale)
{
    SearchBins::drawSBregionDef(ymin_Yields, ymax_Yields, logscale);
}
