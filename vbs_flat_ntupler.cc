/*
 * Implementation of the skimming step of the analysis
 */

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "Math/Vector4D.h"
#include "TStopwatch.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace ROOT::VecOps;

const float MUON_MASS = 0.1056583745;
const float ELE_MASS  = 0.000511;
const float W_MASS = 80.385;
const float Z_MASS = 91.1876;

//lepton cuts
const float LEP_PT_VETO_CUT = 20;
const float EL_PT_CUT = 35;
const float EL_ETA_CUT = 2.5;
const float MU_PT_CUT = 35;
const float MU_ETA_CUT = 2.4;

//ak8 jet cuts
const float AK8_MIN_PT = 200;
const float AK8_MAX_ETA = 2.4;
const float AK8_MIN_SDM = 40;
const float AK8_MAX_SDM = 150;

//ak4 jet cuts
//const float AK4_PT_VETO_CUT = 20;
const float AK4_ETA_CUT = 2.4;
const float AK4_PT_CUT = 30;
const float AK4_JJ_MIN_M = 40.0;
const float AK4_JJ_MAX_M = 150.0;
const float VBF_MJJ_CUT= 500;

//cleaning cuts
const float AK8_LEP_DR_CUT = 1.0;
const float AK4_AK8_DR_CUT = 0.8;
const float AK4_DR_CUT = 0.3;


void vbs_flat_ntupler() {
    ROOT::EnableImplicitMT();
    const auto poolSize = ROOT::GetImplicitMTPoolSize();
    std::cout << "Pool size: " << poolSize << std::endl;

    const std::string samplesBasePath = "root://cmseos.fnal.gov/";
    const std::vector<std::string> sampleNames = {
        "WminusTo2JZTo2LJJ_EWK_LO_SM_MJJ100PTJ10_TuneCP5_13TeV-madgraph-pythia8",
        "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",
    };

    std::ifstream sampleListFile;
    std::string fileLine;
    std::vector<std::string> sampleList;

    for (const auto &sample : sampleNames) {
        std::cout << ">>> Process sample " << sample << ":" << std::endl;
        TStopwatch time;
        time.Start();

        sampleListFile.open(sample + ".txt", std::ifstream::in);
        sampleList.clear();
        while (std::getline(sampleListFile, fileLine))
        {
           std::cout<< fileLine << std::endl;
           sampleList.push_back(samplesBasePath + fileLine);
        }
        sampleListFile.close();

        ROOT::RDataFrame df("Events", sampleList);
        std::cout << "Number of events: " << *df.Count() << std::endl;

        auto h_count = df.Histo1D("event");
        h_count.OnPartialResult(50000, [](TH1D &h_) { h_.Print();});

        auto df2 = df.Filter("HLT_IsoMu24 == true ||\
                              HLT_IsoMu27 == true ||\
                              HLT_Ele27_WPTight_Gsf == true ||\
                              HLT_Ele32_WPTight_Gsf == true ||\
                              HLT_Ele35_WPTight_Gsf == true",
                             "Passes trigger");

        auto df3 = df2.Define("vetoMuons",
                              "Muon_pt > LEP_PT_VETO_CUT && Muon_looseId == true")
                      .Define("tightMuons",
                              "abs(Muon_eta) < MU_ETA_CUT && Muon_pt > MU_PT_CUT && Muon_tightId == true")
                      .Define("nVetoMuons", "Sum(vetoMuons)")
                      .Define("nTightMuons", "Sum(tightMuons)");

        auto df4 = df3.Define("vetoElectrons", 
                              "Electron_pt > LEP_PT_VETO_CUT && Electron_cutBased == 1")
                      .Define("tightElectrons",
                              "abs(Electron_eta) < EL_ETA_CUT && Electron_pt > EL_PT_CUT && Electron_cutBased == 4")
                      .Define("nVetoElectrons", "Sum(vetoElectrons)")
                      .Define("nTightElectrons", "Sum(tightElectrons)");

        auto df5 = df4.Filter("!((nTightMuons + nTightElectrons) == 0)", "no tight lepton")
                      .Filter("!((nVetoMuons + nVetoElectrons) > 2)", "more than two veto lepton")
                      .Filter("!(nTightMuons > 0 && nVetoElectrons > 0)", "mix flavor leptons 1")
                      .Filter("!(nTightElectrons > 0 && nVetoMuons > 0)", "mix flavor leptons 2")
                      .Filter("!(nTightMuons == 1 && nVetoMuons > 1)", "more than ??")
                      .Filter("!(nTightElectrons == 1 && nVetoElectrons > 1)", "more than ??");


        auto lepton_idx = [](RVec<int>& tightLeptons, RVec<float>& pt) 
        {
            auto pt_sorted_indices = Reverse(Argsort(pt));
            auto sorted_tightLeptons = Take(tightLeptons, pt_sorted_indices);
            auto sorted_pt = Take(pt, pt_sorted_indices);
            
            RVec<int> idx = {-1, -1};
            for(size_t i = 0; i < sorted_tightLeptons.size(); i++) {
                 if(sorted_tightLeptons[i] == 1) {
                     if(idx[0] == -1) { idx[0] = i;}
                     if(idx[1] == -1 && idx[0] != int(i)) { idx[1] = i;}
                 }
             }
             return idx;
        };
        auto lepton_ch_idx = [](RVec<int>& tightMuons_idx, RVec<int>& tightElectrons_idx)
        {
            if(tightMuons_idx[0] != -1){
                 return RVec<int>({0, tightMuons_idx[0], tightMuons_idx[1]});
            }
            else
                 return RVec<int>({1, tightElectrons_idx[0], tightElectrons_idx[1]});
        };

        auto df6 = df5.Define("tightMuons_idx", lepton_idx, {"tightMuons", "Muon_pt"})
                      .Define("tightElectrons_idx", lepton_idx, {"tightElectrons", "Electron_pt"})
                      .Define("lepton_idx", lepton_ch_idx, {"tightMuons_idx", "tightElectrons_idx"});

        auto cleanedFatJets = [](RVec<int>& goodFatJets, RVec<float>& FatJet_eta, RVec<float>& FatJet_phi,
                                 RVec<int>& tightMuons, RVec<float>& Muon_eta, RVec<float>& Muon_phi,
                                 RVec<int>& tightElectrons, RVec<float>& Electron_eta, RVec<float>& Electron_phi)
        {
            auto cleanedFatJets_ = goodFatJets;
            auto isClean = true;
            for(size_t i = 0; i < FatJet_eta.size(); i++){
                isClean = true;
                cleanedFatJets_[i] = 0;
                if(goodFatJets[i] == 0) continue;
                
                for(size_t j = 0; j < Muon_eta.size(); j++){
                    if(tightMuons[j] == 0) continue;
                    if(DeltaR(FatJet_eta[i], Muon_eta[j], FatJet_phi[i], Muon_phi[j]) < AK8_LEP_DR_CUT) isClean = false;
                }
                
                for(size_t j = 0; j < Electron_eta.size(); j++){
                    if(tightElectrons[j] == 0) continue;
                    if(DeltaR(FatJet_eta[i], Electron_eta[j], FatJet_phi[i], Electron_phi[j]) < AK8_LEP_DR_CUT) isClean = false;
                }
                
                if(isClean == true) cleanedFatJets_[i] = 1;
            }
            return cleanedFatJets_;
        };

        auto df7 = df6.Define("goodFatJets", "(FatJet_pt > AK8_MIN_PT || FatJet_pt_jesTotalUp > AK8_MIN_PT ||\
                               FatJet_pt_jesTotalDown > AK8_MIN_PT) &&\
                               (abs(FatJet_eta) < AK8_MAX_ETA) &&\
                               (FatJet_msoftdrop > AK8_MIN_SDM || FatJet_msoftdrop_jesTotalUp > AK8_MIN_SDM ||\
                               FatJet_msoftdrop_jesTotalDown > AK8_MIN_SDM) &&\
                               (FatJet_msoftdrop < AK8_MAX_SDM || FatJet_msoftdrop_jesTotalUp < AK8_MAX_SDM ||\
                               FatJet_msoftdrop_jesTotalDown < AK8_MAX_SDM)")
                      .Define("goodCleanedFatJets", cleanedFatJets, {"goodFatJets", "FatJet_eta", "FatJet_phi",
                                                            "tightMuons", "Muon_eta", "Muon_phi",
                                                            "tightElectrons", "Electron_eta", "Electron_phi"});
        auto dfFinal = df7;
        auto report = dfFinal.Report();

        //dfFinal.Snapshot("Events", sample + "Skim.root", finalVariables);
        dfFinal.Snapshot("Events", sample + "Skim.root", "nTightMuons");

        time.Stop();
        report->Print();
        h_count->Print();
        time.Print();
    }
}






/*
 * Find Good Cleaned FatJets Index
 */
template <typename T>
auto FindGoodFatJetsIndex(T &df) {
    using namespace ROOT::VecOps;
    auto selectedFatJet_idx = [](RVec<int>& goodCleanedFatJets, RVec<float>& FatJet_msoftdrop)
            {
                auto idx = -1;
                if (Sum(goodCleanedFatJets) > 0){
                    auto sorted_dM_idx = Argsort(abs(FatJet_msoftdrop[goodCleanedFatJets] - W_MASS));
                    idx = sorted_dM_idx[0];
                }
                 return idx;
            };
    return df.Define("selectedFatJet_idx", selectedFatJet_idx, {"goodCleanedFatJets", "FatJet_msoftdrop"});
}


/*
 * Find Good Jets
 */
template <typename T>
auto FindGoodJets(T &df) {
using namespace ROOT::VecOps;
    auto cleanedJets = [](RVec<int>& goodJets, RVec<float>& Jet_eta, RVec<float>& Jet_phi,
                             RVec<int>& goodFatJets, RVec<float>& FatJet_eta, RVec<float>& FatJet_phi,
                             RVec<int>& tightMuons, RVec<float>& Muon_eta, RVec<float>& Muon_phi,
                             RVec<int>& tightElectrons, RVec<float>& Electron_eta, RVec<float>& Electron_phi)
            {
                auto cleanedJets_ = goodJets;
                auto isClean = true;
                for(size_t i = 0; i < Jet_eta.size(); i++){
                    isClean = true;
                    cleanedJets_[i] = 0;
                    if(goodJets[i] == 0) continue;

                    for(size_t j = 0; j < FatJet_eta.size(); j++){
                        if(goodFatJets[j] == 0) continue;
                        if(DeltaR(Jet_eta[i], FatJet_eta[j], Jet_phi[i], FatJet_phi[j]) < AK4_AK8_DR_CUT) isClean = false;
                    }
                    
                    for(size_t j = 0; j < Jet_eta.size(); j++){
                        if(goodJets[j] == 0) continue;
                        if(DeltaR(Jet_eta[i], Jet_eta[j], Jet_phi[i], Jet_phi[j]) < AK4_DR_CUT) isClean = false;
                    }

                    for(size_t j = 0; j < Muon_eta.size(); j++){
                        if(tightMuons[j] == 0) continue;
                        if(DeltaR(Jet_eta[i], Muon_eta[j], Jet_phi[i], Muon_phi[j]) < AK4_DR_CUT) isClean = false;
                    }
                    
                    for(size_t j = 0; j < Electron_eta.size(); j++){
                        if(tightElectrons[j] == 0) continue;
                        if(DeltaR(Jet_eta[i], Electron_eta[j], Jet_phi[i], Electron_phi[j]) < AK4_DR_CUT) isClean = false;
                    }
                    
                    if(isClean == true) cleanedJets_[i] = 1;
                }
                return cleanedJets_;
            };
    return df.Define("goodJets", "Jet_pt > AK4_PT_CUT || Jet_pt_jesTotalUp > AK4_PT_CUT ||\
                                     Jet_pt_jesTotalDown > AK4_PT_CUT")
             .Define("nBTag_loose", "abs(Jet_eta[goodJets]) < 2.4 && Jet_pt[goodJets] > 30 && Jet_btagDeepB[goodJets] > 0.1241")
             .Define("nBtag_medium", "nBTag_loose && Jet_btagDeepB[goodJets] > 0.4184")
             .Define("nBtag_tight", "nBTag_loose && Jet_btagDeepB[goodJets] > 0.7527")
             .Define("goodCleanedJets", cleanedJets, {"goodJets", "Jet_eta", "Jet_phi",
                                                      "goodCleanedFatJets", "FatJet_eta", "FatJet_phi",
                                                      "tightMuons", "Muon_eta", "Muon_phi",
                                                      "tightElectrons", "Electron_eta", "Electron_phi"});
}


/*
 * Declare leptons variables
 */
template <typename T>
auto DeclareLeptonVariables(T &df) {
    using namespace ROOT::VecOps;
    auto lep = [](RVec<int>& lep_ch_idx, RVec<float>& mu_, RVec<float>& ele_, int idx)
    {
        auto lep_ = mu_;
        if(lep_ch_idx[0] == 1) lep_ = ele_;
        
        if(lep_ch_idx[idx] != -1) return lep_[lep_ch_idx[idx]];
        else return -999.f;
    };
    return df.Define("idx_1", "0")
             .Define("idx_2", "1")
             .Define("lep1_pt", lep, {"lepton_idx", "Muon_pt", "Electron_pt", "idx_1"})
             .Define("lep2_pt", lep, {"lepton_idx", "Muon_pt", "Electron_pt", "idx_2"})
             .Define("lep1_eta", lep, {"lepton_idx", "Muon_eta", "Electron_eta", "idx_1"})
             .Define("lep2_eta", lep, {"lepton_idx", "Muon_eta", "Electron_eta", "idx_2"});
}


/*
 * Add the event weight to the dataset as the column "weight"
 */
/*
template <typename T>
auto AddEventWeight(T &df, const std::string& sample) {
    const auto weight = eventWeights[sample];
    return df.Define("weight", [weight](){ return weight; });
}
*/

/*
 * Declare all variables which shall end up in the final reduced dataset
 */
const std::vector<std::string> finalVariables = {
    "lep1_pt", "lep1_eta", 
    "lep2_pt", "lep2_eta"
};


/*
 * Main function of the skimming step of the analysis
 *
 * The function loops over all required samples, reduces the content to the
 * interesting events and writes them to new files.
 */

