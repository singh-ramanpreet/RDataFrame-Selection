/*
 *
 *
 *
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
const float ELE_MASS = 0.000511;
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
const float VBF_MJJ_CUT = 500;

//cleaning cuts
const float AK8_LEP_DR_CUT = 1.0;
const float AK4_AK8_DR_CUT = 0.8;
const float AK4_DR_CUT = 0.3;

void vbs_flat_ntupler(std::string sample, int year) {
  gErrorIgnoreLevel = kError;

  ROOT::EnableImplicitMT();
  const auto poolSize = ROOT::GetThreadPoolSize();
  std::cout << "Pool size: " << poolSize << std::endl;
  ROOT::EnableImplicitMT(poolSize >= 8 ? 8 : poolSize);

  const std::string samplesBasePath = "root://cmseos.fnal.gov/";

  std::cout << ">>> Process sample " << sample << std::endl;
  TStopwatch time;
  time.Start();

  std::ifstream sampleListFile(sample, std::ifstream::in);
  std::vector<std::string> sampleList;
  std::string fileLine;

  while (std::getline(sampleListFile, fileLine)) {
    //std::cout<< fileLine << std::endl;
    sampleList.push_back(samplesBasePath + fileLine);
  }
  sampleListFile.close();

  auto strip_sampleName = [&sample]() {
    std::string i(sample);
    auto i0 = i.find_last_of("/");
    auto i1 = i.find_last_of(".");
    return i.substr(i0 + 1, i1 - i0 - 1);
  };
  auto sample_basename = strip_sampleName();

  std::cout << ">>> Number of root files in sample " << sampleList.size() << std::endl;

  ROOT::RDataFrame dfRuns("Runs", sampleList);

  // check if sample is MC or Data
  const std::string genEventSumw = dfRuns.HasColumn("genEventSumw_") ? "genEventSumw_" : dfRuns.HasColumn("genEventSumw") ? "genEventSumw" : "None";
  const int isMC = dfRuns.HasColumn(genEventSumw);

  TH1F* totalEvents = new TH1F("TotalEvents", "TotalEvents", 2, -1, 1);
  if (isMC == 1) {
    //genEventCount_ -> Number of Events,
    //genEventSumw_ -> Sum of genWeight,
    //genEventSumw2_ -> Sum of (genWeight *genWeight)
    auto genEventCount_ = dfRuns.Sum(genEventSumw).GetValue();
    totalEvents->SetBinContent(2, genEventCount_);
  }

  ROOT::RDataFrame df("Events", sampleList);

  // make alias for scale variations branches for "Data" samples
  if (isMC == 0) {
    std::cout << ">>> This sample is Data, adding alias for Up/Down branches." << std::endl;
    df.Alias("FatJet_pt_jesTotalUp", "FatJet_pt")
        .Alias("FatJet_pt_jesTotalDown", "FatJet_pt")
        .Alias("FatJet_msoftdrop_jesTotalUp", "FatJet_msoftdrop")
        .Alias("FatJet_msoftdrop_jesTotalDown", "FatJet_msoftdrop")
        .Alias("Jet_pt_jesTotalUp", "Jet_pt")
        .Alias("Jet_pt_jesTotalDown", "Jet_pt");
  }

  //Load stuff, different for different year
  //1. Triggers
  //2. Read histograms as const for lookup
  if (year == 2018) {
    auto f_mu_SF_ID = TFile::Open("./data/SF2018/Muon_Run2018ABCD_SF_ID.root", "read");
    const TH2D* h_mu_SF_ID = (TH2D*)f_mu_SF_ID->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");

    auto f_mu_SF_ISO = TFile::Open("./data/SF2018/Muon_Run2018ABCD_SF_ISO.root", "read");
    const TH2D* h_mu_SF_ISO = (TH2D*)f_mu_SF_ISO->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

    auto f_ele_SF_IDIso = TFile::Open("./data/SF2018/2018_ElectronMVA90.root", "read");
    const TH2D* h_ele_SF_IDIso = (TH2D*)f_ele_SF_IDIso->Get("EGamma_SF2D");

    auto f_mu_FR = TFile::Open("data/SF2018/MuonFR_jet30.root", "read");
    const TH2D* h_mu_FR = (TH2D*)f_mu_FR->Get("FR_pT_eta");

    auto f_ele_FR = TFile::Open("data/SF2018/EleFR_jet30.root", "read");
    const TH2D* h_ele_FR = (TH2D*)f_ele_FR->Get("FR_pT_eta");
  }

  auto count = df.Count();
  int total_entries = count.GetValue();

  std::cout << ">>> Number of events: " << total_entries << std::endl;

  int progress = 0;
  std::mutex bar_mutex;
  auto bar_print = [&progress, &bar_mutex](unsigned int, ULong64_t&) {
    std::lock_guard<std::mutex> lg(bar_mutex);
    progress += 1;
    std::cout << ">>> Processed: " << progress << "%" << std::endl;
  };
  count.OnPartialResultSlot(total_entries / 100, bar_print);

  // chain dataframes in sequential order
  auto chainedDf = ROOT::RDF::RNode(df);
  chainedDf = chainedDf.Filter(
      "HLT_IsoMu24 == true ||\
                          HLT_IsoMu27 == true ||\
                          HLT_Ele27_WPTight_Gsf == true ||\
                          HLT_Ele32_WPTight_Gsf == true ||\
                          HLT_Ele35_WPTight_Gsf == true",
      "Passes trigger");

  chainedDf = chainedDf.Define("vetoMuons", "Muon_pt > LEP_PT_VETO_CUT && Muon_looseId == true")
                  .Define("tightMuons", "abs(Muon_eta) < MU_ETA_CUT && Muon_pt > MU_PT_CUT && Muon_tightId == true")
                  .Define("nVetoMuons", "Sum(vetoMuons)")
                  .Define("nTightMuons", "Sum(tightMuons)");

  chainedDf = chainedDf.Define("vetoElectrons", "Electron_pt > LEP_PT_VETO_CUT && Electron_cutBased == 1")
                  .Define("tightElectrons", "abs(Electron_eta) < EL_ETA_CUT && Electron_pt > EL_PT_CUT && Electron_cutBased == 4")
                  .Define("nVetoElectrons", "Sum(vetoElectrons)")
                  .Define("nTightElectrons", "Sum(tightElectrons)");

  chainedDf = chainedDf.Filter("!((nTightMuons + nTightElectrons) == 0)", "no tight lepton")
                  .Filter("!((nVetoMuons + nVetoElectrons) > 2)", "more than two veto lepton")
                  .Filter("!(nTightMuons > 0 && nVetoElectrons > 0)", "mix flavor leptons 1")
                  .Filter("!(nTightElectrons > 0 && nVetoMuons > 0)", "mix flavor leptons 2")
                  .Filter("!(nTightMuons == 1 && nVetoMuons > 1)", "more than ??")
                  .Filter("!(nTightElectrons == 1 && nVetoElectrons > 1)", "more than ??");

  //function: determine index of selected lepton
  //... first sorted in pt
  //... returned two indexes
  //... -1 => not selected
  auto lepton_idx = [](RVec<int>& tightLeptons, RVec<float>& pt) {
    auto pt_sorted_indices = Reverse(Argsort(pt));

    RVec<int> idx = {-1, -1};
    for (auto& i : pt_sorted_indices) {
      if (tightLeptons[i] == 1) {
        if (idx[0] == -1)
          idx[0] = int(i);
        if (idx[1] == -1 && idx[0] != int(i))
          idx[1] = int(i);
      }
    }
    return idx;
  };

  chainedDf = chainedDf.Define("mu_idx", lepton_idx, {"tightMuons", "Muon_pt"}).Define("ele_idx", lepton_idx, {"tightElectrons", "Electron_pt"});

  //function: determine clean fatjets
  //... clean w.r.t muons and, electrons
  //... returned fatjets Id
  auto cleanedFatJets = [](RVec<int>& fj_id,
                           RVec<float>& fj_eta,
                           RVec<float>& fj_phi,
                           RVec<int>& mu_id,
                           RVec<float>& mu_eta,
                           RVec<float>& mu_phi,
                           RVec<int>& el_id,
                           RVec<float>& el_eta,
                           RVec<float>& el_phi) {
    RVec<int> cleanedFatJets_(fj_id.size(), 0);
    auto isClean = true;
    for (auto& i : Nonzero(fj_id)) {
      isClean = true;

      for (auto& j : Nonzero(mu_id)) {
        if (DeltaR(fj_eta[i], mu_eta[j], fj_phi[i], mu_phi[j]) < AK8_LEP_DR_CUT)
          isClean = false;
      }
      if (isClean == false)
        continue;

      for (auto& j : Nonzero(el_id)) {
        if (DeltaR(fj_eta[i], el_eta[j], fj_phi[i], el_phi[j]) < AK8_LEP_DR_CUT)
          isClean = false;
      }

      if (isClean == true)
        cleanedFatJets_[i] = 1;
    }
    return cleanedFatJets_;
  };

  //Define which fatjets are good and clean
  //... no selection
  chainedDf =
      chainedDf
          .Define("goodFatJets",
                  "(FatJet_pt_nom > AK8_MIN_PT || FatJet_pt_jesTotalUp > AK8_MIN_PT ||\
                           FatJet_pt_jesTotalDown > AK8_MIN_PT) &&\
                           (abs(FatJet_eta) < AK8_MAX_ETA) &&\
                           (FatJet_msoftdrop > AK8_MIN_SDM || FatJet_msoftdrop_jesTotalUp > AK8_MIN_SDM ||\
                           FatJet_msoftdrop_jesTotalDown > AK8_MIN_SDM) &&\
                           (FatJet_msoftdrop < AK8_MAX_SDM || FatJet_msoftdrop_jesTotalUp < AK8_MAX_SDM ||\
                           FatJet_msoftdrop_jesTotalDown < AK8_MAX_SDM)")
          .Define(
              "goodCleanedFatJets",
              cleanedFatJets,
              {"goodFatJets", "FatJet_eta", "FatJet_phi", "tightMuons", "Muon_eta", "Muon_phi", "tightElectrons", "Electron_eta", "Electron_phi"});

  //function: determine selected fatjet index
  //... returned index
  //... -1 => no fatjet selected
  auto selectedFatJet_idx = [](RVec<int>& fj_id, RVec<float>& fj_m) {
    int idx = -1;
    if (Sum(fj_id) > 0) {
      auto sorted_dM_idx = Argsort(abs(fj_m - W_MASS));
      auto sorted_fj_id = Take(fj_id, sorted_dM_idx);
      if (sorted_dM_idx[sorted_fj_id].size() > 0) {
        idx = int(sorted_dM_idx[sorted_fj_id][0]);
      }
    }
    return idx;
  };

  //Define selected fatjet index
  //... no selection
  chainedDf = chainedDf.Define("selectedFatJet_idx", selectedFatJet_idx, {"goodCleanedFatJets", "FatJet_msoftdrop"});

  //function: determine clean jets
  //... w.r.t. muons, electrons, fatjet
  //... and with other jets
  //... returned clean jets ID
  auto cleanedJets = [](RVec<int>& jt_id,
                        RVec<float>& jt_eta,
                        RVec<float>& jt_phi,
                        RVec<int>& fj_id,
                        RVec<float>& fj_eta,
                        RVec<float>& fj_phi,
                        RVec<int>& mu_id,
                        RVec<float>& mu_eta,
                        RVec<float>& mu_phi,
                        RVec<int>& el_id,
                        RVec<float>& el_eta,
                        RVec<float>& el_phi) {
    RVec<int> cleanedJets_(jt_id.size(), 0);
    auto isClean = true;
    for (auto& i : Nonzero(jt_id)) {
      isClean = true;

      for (auto& j : Nonzero(jt_id)) {
        if (i == j)
          continue;
        if (DeltaR(jt_eta[i], jt_eta[j], jt_phi[i], jt_phi[j]) < AK4_DR_CUT)
          isClean = false;
      }
      if (isClean == false)
        continue;

      for (auto& j : Nonzero(fj_id)) {
        if (DeltaR(jt_eta[i], fj_eta[j], jt_phi[i], fj_phi[j]) < AK4_AK8_DR_CUT)
          isClean = false;
      }
      if (isClean == false)
        continue;

      for (auto& j : Nonzero(mu_id)) {
        if (DeltaR(jt_eta[i], mu_eta[j], jt_phi[i], mu_phi[j]) < AK4_DR_CUT)
          isClean = false;
      }
      if (isClean == false)
        continue;

      for (auto& j : Nonzero(el_id)) {
        if (DeltaR(jt_eta[i], el_eta[j], jt_phi[i], el_phi[j]) < AK4_DR_CUT)
          isClean = false;
      }

      if (isClean == true)
        cleanedJets_[i] = 1;
    }
    return cleanedJets_;
  };

  //Define good jets
  //Define number of btags, various IDs
  //Define goog and clean jets IDs
  //Select Events with atleast two jets
  chainedDf = chainedDf
                  .Define("goodJets",
                          "Jet_pt > AK4_PT_CUT ||\
                                       Jet_pt_jesTotalUp > AK4_PT_CUT ||\
                                       Jet_pt_jesTotalDown > AK4_PT_CUT")
                  .Define("nBtag_loose",
                          "Sum(abs(Jet_eta[goodJets]) < 2.4 && Jet_pt[goodJets] > 30 &&\
                                          Jet_btagDeepB[goodJets] > 0.1241)")
                  .Define("nBtag_medium",
                          "Sum(nBtag_loose &&\
                                           Jet_btagDeepB[goodJets] > 0.4184)")
                  .Define("nBtag_tight",
                          "Sum(nBtag_loose &&\
                                          Jet_btagDeepB[goodJets] > 0.7527)")
                  .Define("goodCleanedJets",
                          cleanedJets,
                          {"goodJets",
                           "Jet_eta",
                           "Jet_phi",
                           "goodCleanedFatJets",
                           "FatJet_eta",
                           "FatJet_phi",
                           "tightMuons",
                           "Muon_eta",
                           "Muon_phi",
                           "tightElectrons",
                           "Electron_eta",
                           "Electron_phi"})
                  .Filter("Sum(goodCleanedJets) > 2", "at least 2 AK4 cleaned jets");

  //function: determine selected jets indexes
  //... first figure out event is boosted or resolved
  //... if boosted, selected highest mjj vbf jets
  //... if resolved first select hadronic jets
  //... ... then mjj vbj jets
  //... returned indexes {boson jet1,jet2,vbf jet1, jet2 }
  auto selectedJets_idx = [](int selected_fj, RVec<int> jt_id, RVec<float> jt_pt, RVec<float> jt_eta, RVec<float> jt_phi, RVec<float> jt_m) {
    auto idx = Combinations(jt_id, 2);
    auto jt_pt_1 = Take(jt_pt, idx[0]);
    auto jt_eta_1 = Take(jt_eta, idx[0]);
    auto jt_phi_1 = Take(jt_phi, idx[0]);
    auto jt_m_1 = Take(jt_m, idx[0]);

    auto jt_pt_2 = Take(jt_pt, idx[1]);
    auto jt_eta_2 = Take(jt_eta, idx[1]);
    auto jt_phi_2 = Take(jt_phi, idx[1]);
    auto jt_m_2 = Take(jt_m, idx[1]);

    auto mass_jj = InvariantMasses(jt_pt_1, jt_eta_1, jt_phi_1, jt_m_1, jt_pt_2, jt_eta_2, jt_phi_2, jt_m_2);

    int selected_bos_j1 = -1, selected_bos_j2 = -1;
    int selected_vbf_j1 = -1, selected_vbf_j2 = -1;

    if (selected_fj == -1) {
      auto sorted_mass_bos = Argsort(abs(mass_jj - W_MASS));

      for (auto& i : sorted_mass_bos) {
        const int j1_idx = idx[0][i];
        const int j2_idx = idx[1][i];
        if (fabs(jt_eta[j1_idx]) > AK4_ETA_CUT)
          continue;
        if (fabs(jt_eta[j2_idx]) > AK4_ETA_CUT)
          continue;
        if (mass_jj[i] < AK4_JJ_MIN_M || mass_jj[i] > AK4_JJ_MAX_M)
          continue;
        selected_bos_j1 = j1_idx;
        selected_bos_j2 = j2_idx;
        break;
      }
    }

    auto sorted_mass_vbf = Reverse(Argsort(mass_jj));
    for (auto& i : sorted_mass_vbf) {
      const int j1_idx = idx[0][i];
      const int j2_idx = idx[1][i];
      if (mass_jj[i] < VBF_MJJ_CUT)
        continue;
      if (jt_eta[j1_idx] * jt_eta[j2_idx] > 0)
        continue;
      if (j1_idx == selected_bos_j1 || j1_idx == selected_bos_j2)
        continue;
      if (j2_idx == selected_bos_j1 || j2_idx == selected_bos_j2)
        continue;
      selected_vbf_j1 = j1_idx;
      selected_vbf_j2 = j2_idx;
      break;
    }

    return RVec<int>{selected_bos_j1, selected_bos_j2, selected_vbf_j1, selected_vbf_j2};
  };

  //Define selected jets indexes
  //Apply boosted or resolved Filter
  chainedDf =
      chainedDf.Define("selectedJets_idx", selectedJets_idx, {"selectedFatJet_idx", "goodCleanedJets", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
          .Define("bos_j1", "selectedJets_idx[0]")
          .Define("bos_j2", "selectedJets_idx[1]")
          .Define("vbf_j1", "selectedJets_idx[2]")
          .Define("vbf_j2", "selectedJets_idx[3]")
          .Filter(
              "selectedFatJet_idx != -1 ||\
                            (bos_j1 != -1 && bos_j2 != -1)",
              "Boosted or Resolved")
          .Filter("vbf_j1 != -1", "has selected AK4 vbf 1")
          .Filter("vbf_j2 != -1", "has selected AK4 vbf 2");

  //== METz Calculator ==
  //====== START ========
  int metz_type = 0;
  auto METzCalculator = [&metz_type](float l_pt, float l_eta, float l_phi, float l_m, float met_pt, float met_phi) {
    auto isComplex_ = false;

    auto met_x = met_pt * (TMath::Cos(met_phi));
    auto met_y = met_pt * (TMath::Sin(met_phi));
    auto lep_ = ROOT::Math::PtEtaPhiMVector(l_pt, l_eta, l_phi, l_m);
    auto met_ = ROOT::Math::PxPyPzEVector(met_x, met_y, 0.0, met_pt);

    double M_W = W_MASS;
    double M_mu = l_m, emu = lep_.E(), pxmu = lep_.Px(), pymu = lep_.Py(), pzmu = lep_.Pz();
    double pxnu = met_.Px(), pynu = met_.Py(), ptnu = met_.Pt(), pznu = 0.;
    double otherSol_ = 0.;
    double newPtneutrino1_ = -1, newPtneutrino2_ = -1;
    double a = M_W * M_W / 2 + pxnu * pxmu + pynu * pymu;
    double A = pzmu * pzmu - emu * emu;
    double B = 2 * pzmu * a;
    double C = a * a - emu * emu * ptnu * ptnu;

    double tmproot = B * B - 4.0 * A * C;
    if (tmproot < 0) {
      isComplex_ = true;
      pznu = -B / (2 * A);  // take real part of complex roots
      otherSol_ = pznu;
      // recalculate the neutrino pT
      // solve quadratic eq. discriminator = 0 for pT of nu
      double pnu = met_.E();
      double Delta = (M_W * M_W - M_mu * M_mu);
      double alpha = (pxmu * pxnu / pnu + pymu * pynu / pnu);
      double ptnu = TMath::Sqrt(pxnu * pxnu + pynu * pynu);  // old
      double AA = 4. * pzmu * pzmu - 4 * emu * emu + 4 * alpha * alpha;
      double BB = 4. * alpha * Delta;
      double CC = Delta * Delta;
      double tmpdisc = BB * BB - 4.0 * AA * CC;
      double tmpsolpt1 = (-BB + TMath::Sqrt(tmpdisc)) / (2.0 * AA);
      double tmpsolpt2 = (-BB - TMath::Sqrt(tmpdisc)) / (2.0 * AA);
      if (fabs(tmpsolpt1 - ptnu) < fabs(tmpsolpt2 - ptnu)) {
        newPtneutrino1_ = tmpsolpt1;
        newPtneutrino2_ = tmpsolpt2;
      } else {
        newPtneutrino1_ = tmpsolpt2;
        newPtneutrino2_ = tmpsolpt1;
      }
    } else {
      isComplex_ = false;
      double tmpsol1 = (-B + TMath::Sqrt(tmproot)) / (2.0 * A);
      double tmpsol2 = (-B - TMath::Sqrt(tmproot)) / (2.0 * A);
      //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;
      if (metz_type == 0) {
        // two real roots, pick the one closest to pz of muon
        if (TMath::Abs(tmpsol2 - pzmu) < TMath::Abs(tmpsol1 - pzmu)) {
          pznu = tmpsol2;
          otherSol_ = tmpsol1;
        } else {
          pznu = tmpsol1;
          otherSol_ = tmpsol2;
        }
        // if pznu is > 300 pick the most central root
        if (pznu > 300.) {
          if (TMath::Abs(tmpsol1) < TMath::Abs(tmpsol2)) {
            pznu = tmpsol1;
            otherSol_ = tmpsol2;
          } else {
            pznu = tmpsol2;
            otherSol_ = tmpsol1;
          }
        }
      }
    }
    return float(pznu);
  };

  //== METz Calculator ==
  //======= END =========

  // function: Calculate p4 vector for di-leptons
  // ... W -> lepton + MET
  // ... Z -> lepton + lepton same flavor
  // ... returned RVec pt, eta, phi, M
  auto dilep_p4 = [](float l_pt1,
                     float l_eta1,
                     float l_phi1,
                     float l_m1,
                     float l_pt2,
                     float l_eta2,
                     float l_phi2,
                     float l_m2,
                     float met_x,
                     float met_y,
                     float met_z,
                     float met_ett) {
    ROOT::Math::PtEtaPhiMVector vector1, vector2;
    vector1.SetPt(l_pt1);
    vector1.SetEta(l_eta1);
    vector1.SetPhi(l_phi1);
    vector1.SetM(l_m1);
    if (l_pt2 < 0) {
      float met_e = TMath::Sqrt(met_x * met_x + met_y * met_y + met_z * met_z);
      vector2.SetPxPyPzE(met_x, met_y, met_z, met_e);
    } else {
      vector2.SetPt(l_pt2);
      vector2.SetEta(l_eta2);
      vector2.SetPhi(l_phi2);
      vector2.SetM(l_m2);
    }

    auto dilep_ = vector1 + vector2;
    return RVec<double>{dilep_.Pt(), dilep_.Eta(), dilep_.Phi(), dilep_.M()};
  };

  //Collect list of output variables
  std::vector<std::string> finalVariables;
  finalVariables.push_back("run");

  //Some Utils
  TString select_lepton1 = "mu_idx[0] != -1? %s[mu_idx[0]]: ele_idx[0] != -1? %s[ele_idx[0]]: -999.f";
  TString select_lepton2 = "mu_idx[1] != -1? %s[mu_idx[1]]: ele_idx[1] != -1? %s[ele_idx[1]]: -999.f";

  //Define output variables for ntuples
  //do not define already existing variables
  //and save in std vector of strings 'finalVariables'
  chainedDf = chainedDf.Define("evt", "event");
  finalVariables.push_back("evt");

  chainedDf = chainedDf.Define("lep1_m", "mu_idx[0] != -1? MUON_MASS: ele_idx[0] != -1? ELE_MASS: -999.f");
  finalVariables.push_back("lep1_m");
  chainedDf = chainedDf.Define("lep2_m", "mu_idx[1] != -1? MUON_MASS: ele_idx[1] != -1? ELE_MASS: -999.f");
  finalVariables.push_back("lep2_m");

  chainedDf = chainedDf.Define("lep1_pt", Form(select_lepton1, "Muon_pt", "Electron_pt"));
  finalVariables.push_back("lep1_pt");
  chainedDf = chainedDf.Define("lep1_eta", Form(select_lepton1, "Muon_eta", "Electron_eta"));
  finalVariables.push_back("lep1_eta");
  chainedDf = chainedDf.Define("lep1_phi", Form(select_lepton1, "Muon_phi", "Electron_phi"));
  finalVariables.push_back("lep1_phi");
  chainedDf = chainedDf.Define("lep1_q", Form(select_lepton1, "Muon_charge", "Electron_charge"));
  finalVariables.push_back("lep1_q");
  chainedDf = chainedDf.Define("lep1_iso", Form(select_lepton1, "Muon_pfRelIso04_all", "Electron_pfRelIso03_all"));
  finalVariables.push_back("lep1_iso");

  chainedDf = chainedDf.Define("lep2_pt", Form(select_lepton2, "Muon_pt", "Electron_pt"));
  finalVariables.push_back("lep2_pt");
  chainedDf = chainedDf.Define("lep2_eta", Form(select_lepton2, "Muon_eta", "Electron_eta"));
  finalVariables.push_back("lep2_eta");
  chainedDf = chainedDf.Define("lep2_phi", Form(select_lepton2, "Muon_phi", "Electron_phi"));
  finalVariables.push_back("lep2_phi");
  chainedDf = chainedDf.Define("lep2_q", Form(select_lepton2, "Muon_charge", "Electron_charge"));
  finalVariables.push_back("lep2_q");
  chainedDf = chainedDf.Define("lep2_iso", Form(select_lepton2, "Muon_pfRelIso04_all", "Electron_pfRelIso03_all"));
  finalVariables.push_back("lep2_iso");

  finalVariables.push_back("MET_pt");
  finalVariables.push_back("MET_phi");
  chainedDf = chainedDf.Define("MET_px", "MET_pt * cos(MET_phi)");
  chainedDf = chainedDf.Define("MET_py", "MET_pt * sin(MET_phi)");
  chainedDf = chainedDf.Define("neu_pz_type0", METzCalculator, {"lep1_pt", "lep1_eta", "lep1_phi", "lep1_m", "MET_pt", "MET_phi"});
  finalVariables.push_back("neu_pz_type0");

  chainedDf = chainedDf.Define(
      "dilep_p4",
      dilep_p4,
      {"lep1_pt", "lep1_eta", "lep1_phi", "lep1_m", "lep2_pt", "lep2_eta", "lep2_phi", "lep2_m", "MET_px", "MET_py", "neu_pz_type0", "MET_pt"});
  chainedDf = chainedDf.Define("dilep_pt", "dilep_p4[0]");
  finalVariables.push_back("dilep_pt");
  chainedDf = chainedDf.Define("dilep_eta", "dilep_p4[1]");
  finalVariables.push_back("dilep_eta");
  chainedDf = chainedDf.Define("dilep_phi", "dilep_p4[2]");
  finalVariables.push_back("dilep_phi");
  chainedDf = chainedDf.Define("dilep_m", "dilep_p4[3]");
  finalVariables.push_back("dilep_m");

  chainedDf = chainedDf.Define("bos_PuppiAK8_m_sd0_corr", "selectedFatJet_idx != -1? FatJet_msoftdrop[selectedFatJet_idx]: -999.f");
  finalVariables.push_back("bos_PuppiAK8_m_sd0_corr");
  chainedDf = chainedDf.Define("bos_PuppiAK8_pt", "selectedFatJet_idx != -1? FatJet_pt[selectedFatJet_idx]: -999.f");
  finalVariables.push_back("bos_PuppiAK8_pt");

  finalVariables.push_back("nBtag_loose");
  finalVariables.push_back("nBtag_medium");
  finalVariables.push_back("nBtag_tight");

  if (isMC == 1) {
    chainedDf = chainedDf.Define("nScaleWeight", "nLHEScaleWeight");
    chainedDf = chainedDf.Define("scaleWeight", "LHEScaleWeight");
    finalVariables.push_back("nScaleWeight");
    finalVariables.push_back("scaleWeight");

    chainedDf = chainedDf.Define("nPdfWeight", "nLHEPdfWeight");
    chainedDf = chainedDf.Define("pdfWeight", "LHEPdfWeight");
    finalVariables.push_back("nPdfWeight");
    finalVariables.push_back("pdfWeight");
  }

  auto outputFileName = sample_basename + ".root";
  std::cout << ">>> Output Filename: " << outputFileName << std::endl;

  auto dfFinal = chainedDf;
  auto report = dfFinal.Report();

  dfFinal.Snapshot("Events", outputFileName, finalVariables);

  report->Print();

  TFile* outFile = TFile::Open(outputFileName.c_str(), "update");
  outFile->cd();
  totalEvents->Write();
  outFile->Close();

  time.Stop();
  time.Print();
}
