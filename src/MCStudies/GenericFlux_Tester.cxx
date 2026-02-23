// Copyright 2016-2021 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
 *    This file is part of NUISANCE.
 *
 *    NUISANCE is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    NUISANCE is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#include "GenericFlux_Tester.h"

//********************************************************************
/// @brief Class to perform MC Studies on a custom measurement
GenericFlux_Tester::GenericFlux_Tester(std::string name, std::string inputfile,
                                       FitWeight *rw, std::string type,
                                       std::string fakeDataFile) {
  //********************************************************************

  // Measurement Details
  fName = name;
  eventVariables = NULL;

  // Define our energy range for flux calcs
  EnuMin = 0.;
  EnuMax = 1E10; // Arbritrarily high energy limit

  // Set default fitter flags
  fIsDiag = true;
  fIsShape = false;
  fIsRawEvents = false;

  nu_4mom = new TLorentzVector(0, 0, 0, 0);
  pmu = new TLorentzVector(0, 0, 0, 0);
  ppip = new TLorentzVector(0, 0, 0, 0);
  ppim = new TLorentzVector(0, 0, 0, 0);
  ppi0 = new TLorentzVector(0, 0, 0, 0);
  pprot = new TLorentzVector(0, 0, 0, 0);
  pneut = new TLorentzVector(0, 0, 0, 0);

  // This function will sort out the input files automatically and parse all the
  // inputs,flags,etc.
  // There may be complex cases where you have to do this by hand, but usually
  // this will do.
  Measurement1D::SetupMeasurement(inputfile, type, rw, fakeDataFile);

  eventVariables = NULL;
  liteMode = Config::Get().GetParB("isLiteMode");

  if (Config::HasPar("EnuMin")) {
    EnuMin = Config::GetParD("EnuMin");
  }

  if (Config::HasPar("EnuMax")) {
    EnuMax = Config::GetParD("EnuMax");
  }

  // Setup fDataHist as a placeholder
  this->fDataHist = new TH1D(("empty_data"), ("empty-data"), 1, 0, 1);
  this->SetupDefaultHist();
  fFullCovar = StatUtils::MakeDiagonalCovarMatrix(fDataHist);
  covar = StatUtils::GetInvert(fFullCovar);

  // 1. The generator is organised in SetupMeasurement so it gives the
  // cross-section in "per nucleon" units.
  //    So some extra scaling for a specific measurement may be required. For
  //    Example to get a "per neutron" measurement on carbon
  //    which we do here, we have to multiple by the number of nucleons 12 and
  //    divide by the number of neutrons 6.
  // N.B. MeasurementBase::PredictedEventRate includes the 1E-38 factor that is
  // often included here in other classes that directly integrate the event
  // histogram. This method is used here as it now respects EnuMin and EnuMax
  // correctly.
  this->fScaleFactor =
      (this->PredictedEventRate("width") / double(fNEvents)) /
      this->TotalIntegratedFlux();
  if (fScaleFactor <= 0.0) {
    NUIS_ABORT("SCALE FACTOR TOO LOW: " << fScaleFactor);
  }

  NUIS_LOG(SAM, " Generic Flux Scaling Factor = "
                << fScaleFactor
                << " [= " << (GetEventHistogram()->Integral("width") * 1E-38)
                << "/(" << (fNEvents + 0.) << "*" << this->TotalIntegratedFlux()
                << ")]");

  // Setup our TTrees
  this->AddEventVariablesToTree();
  this->AddSignalFlagsToTree();

  Fill_ICARUS_QELike_Variable = Config::Get().GetParB("Add_ICARUS_QELike");
  if( Fill_ICARUS_QELike_Variable ){
    NUIS_LOG(SAM, " Generic Flux Adding ICARUS QELike variables");
    this->AddICARUS1muNp0piVariablesToTree();
    //this->AddICARUS1mu2p0piVariablesToTree();
  }

  Fill_SBND_QELike_Variable = Config::Get().GetParB("Add_SBND_QELike");
  if( Fill_SBND_QELike_Variable ){
    NUIS_LOG(SAM, " Generic Flux Adding SBND QELike variables");
    this->AddSBND1mu1p0piVariablesToTree();
  }

  Fill_ICARUS_1mu1pi0_Variable = Config::Get().GetParB("Add_ICARUS_pi0");
  if(Fill_ICARUS_1mu1pi0_Variable){
    this->AddICARUS1mu1pi0VariablesToTree();
  }

}

void GenericFlux_Tester::AddEventVariablesToTree() {
  // Setup the TTree to save everything
  if (!eventVariables) {
    Config::Get().out->cd();
    eventVariables = new TTree((this->fName + "_VARS").c_str(),
                               (this->fName + "_VARS").c_str());
  }

  NUIS_LOG(SAM, "Adding Event Variables");
  eventVariables->Branch("Mode", &Mode, "Mode/I");
  eventVariables->Branch("ResCode", &ResCode, "ResCode/I");

  eventVariables->Branch("PDGnu", &PDGnu, "PDGnu/I");
  eventVariables->Branch("Enu_true", &Enu_true, "Enu_true/F");

  eventVariables->Branch("Nleptons", &Nleptons, "Nleptons/I");
  // all sensible
  eventVariables->Branch("MLep", &MLep, "MLep/F");
  eventVariables->Branch("ELep", &ELep, "ELep/F");
  // negative -999
  eventVariables->Branch("TLep", &TLep, "TLep/F");
  eventVariables->Branch("CosLep", &CosLep, "CosLep/F");
  eventVariables->Branch("CosPmuPpip", &CosPmuPpip, "CosPmuPpip/F");
  eventVariables->Branch("CosPmuPpim", &CosPmuPpim, "CosPmuPpim/F");
  eventVariables->Branch("CosPmuPpi0", &CosPmuPpi0, "CosPmuPpi0/F");
  eventVariables->Branch("CosPmuPprot", &CosPmuPprot, "CosPmuPprot/F");
  eventVariables->Branch("CosPmuPneut", &CosPmuPneut, "CosPmuPneut/F");

  eventVariables->Branch("Nprotons", &Nprotons, "Nprotons/I");
  eventVariables->Branch("MPr", &MPr, "MPr/F");
  eventVariables->Branch("EPr", &EPr, "EPr/F");
  eventVariables->Branch("TPr", &TPr, "TPr/F");
  eventVariables->Branch("CosPr", &CosPr, "CosPr/F");
  eventVariables->Branch("CosPprotPneut", &CosPprotPneut, "CosPprotPneut/F");

  eventVariables->Branch("Nneutrons", &Nneutrons, "Nneutrons/I");
  eventVariables->Branch("MNe", &MNe, "MNe/F");
  eventVariables->Branch("ENe", &ENe, "ENe/F");
  eventVariables->Branch("TNe", &TNe, "TNe/F");
  eventVariables->Branch("CosNe", &CosNe, "CosNe/F");

  eventVariables->Branch("Npiplus", &Npiplus, "Npiplus/I");
  eventVariables->Branch("MPiP", &MPiP, "MPiP/F");
  eventVariables->Branch("EPiP", &EPiP, "EPiP/F");
  eventVariables->Branch("TPiP", &TPiP, "TPiP/F");
  eventVariables->Branch("CosPiP", &CosPiP, "CosPiP/F");
  eventVariables->Branch("CosPpipPprot", &CosPpipPprot, "CosPpipProt/F");
  eventVariables->Branch("CosPpipPneut", &CosPpipPneut, "CosPpipPneut/F");
  eventVariables->Branch("CosPpipPpim", &CosPpipPpim, "CosPpipPpim/F");
  eventVariables->Branch("CosPpipPpi0", &CosPpipPpi0, "CosPpipPpi0/F");

  eventVariables->Branch("Npineg", &Npineg, "Npineg/I");
  eventVariables->Branch("MPiN", &MPiN, "MPiN/F");
  eventVariables->Branch("EPiN", &EPiN, "EPiN/F");
  eventVariables->Branch("TPiN", &TPiN, "TPiN/F");
  eventVariables->Branch("CosPiN", &CosPiN, "CosPiN/F");
  eventVariables->Branch("CosPpimPprot", &CosPpimPprot, "CosPpimPprot/F");
  eventVariables->Branch("CosPpimPneut", &CosPpimPneut, "CosPpimPneut/F");
  eventVariables->Branch("CosPpimPpi0", &CosPpimPpi0, "CosPpimPpi0/F");

  eventVariables->Branch("Npi0", &Npi0, "Npi0/I");
  eventVariables->Branch("MPi0", &MPi0, "MPi0/F");
  eventVariables->Branch("EPi0", &EPi0, "EPi0/F");
  eventVariables->Branch("TPi0", &TPi0, "TPi0/F");
  eventVariables->Branch("CosPi0", &CosPi0, "CosPi0/F");
  eventVariables->Branch("CosPi0Pprot", &CosPi0Pprot, "CosPi0Pprot/F");
  eventVariables->Branch("CosPi0Pneut", &CosPi0Pneut, "CosPi0Pneut/F");

  eventVariables->Branch("Nother", &Nother, "Nother/I");

  eventVariables->Branch("Q2_true", &Q2_true, "Q2_true/F");
  eventVariables->Branch("q0_true", &q0_true, "q0_true/F");
  eventVariables->Branch("q3_true", &q3_true, "q3_true/F");
  eventVariables->Branch("Emiss", &Emiss, "Emiss/F");
  eventVariables->Branch("pmiss", &pmiss);
  eventVariables->Branch("Emiss_preFSI", &Emiss_preFSI, "Emiss_preFSI/F");
  eventVariables->Branch("pmiss_preFSI", &pmiss_preFSI);

  eventVariables->Branch("Enu_QE", &Enu_QE, "Enu_QE/F");
  eventVariables->Branch("Q2_QE", &Q2_QE, "Q2_QE/F");

  eventVariables->Branch("W_nuc_rest", &W_nuc_rest, "W_nuc_rest/F");
  eventVariables->Branch("bjorken_x", &bjorken_x, "bjorken_x/F");
  eventVariables->Branch("bjorken_y", &bjorken_y, "bjorken_y/F");

  eventVariables->Branch("Erecoil_true", &Erecoil_true, "Erecoil_true/F");
  eventVariables->Branch("Erecoil_charged", &Erecoil_charged,
                         "Erecoil_charged/F");
  eventVariables->Branch("Erecoil_minerva", &Erecoil_minerva,
                         "Erecoil_minerva/F");

  if (!liteMode) {
    eventVariables->Branch("nu_4mom", &nu_4mom);
    eventVariables->Branch("pmu_4mom", &pmu);
    eventVariables->Branch("hm_ppip_4mom", &ppip);
    eventVariables->Branch("hm_ppim_4mom", &ppim);
    eventVariables->Branch("hm_ppi0_4mom", &ppi0);
    eventVariables->Branch("hm_pprot_4mom", &pprot);
    eventVariables->Branch("hm_pneut_4mom", &pneut);
  }

  // Event Scaling Information
  eventVariables->Branch("Weight", &Weight, "Weight/F");
  eventVariables->Branch("InputWeight", &InputWeight, "InputWeight/F");
  eventVariables->Branch("RWWeight", &RWWeight, "RWWeight/F");
  eventVariables->Branch("FluxWeight", &FluxWeight, "FluxWeight/F");
  eventVariables->Branch("fScaleFactor", &fScaleFactor, "fScaleFactor/D");

  return;
}

void GenericFlux_Tester::AddSignalFlagsToTree() {
  if (!eventVariables) {
    Config::Get().out->cd();
    eventVariables = new TTree((this->fName + "_VARS").c_str(),
                               (this->fName + "_VARS").c_str());
  }

  NUIS_LOG(SAM, "Adding signal flags");

  // Signal Definitions from SignalDef.cxx
  eventVariables->Branch("flagCCINC", &flagCCINC, "flagCCINC/O");
  eventVariables->Branch("flagNCINC", &flagNCINC, "flagNCINC/O");
  eventVariables->Branch("flagCCQE", &flagCCQE, "flagCCQE/O");
  eventVariables->Branch("flagCC0pi", &flagCC0pi, "flagCC0pi/O");
  eventVariables->Branch("flagCCQELike", &flagCCQELike, "flagCCQELike/O");
  eventVariables->Branch("flagNCEL", &flagNCEL, "flagNCEL/O");
  eventVariables->Branch("flagNC0pi", &flagNC0pi, "flagNC0pi/O");
  eventVariables->Branch("flagCCcoh", &flagCCcoh, "flagCCcoh/O");
  eventVariables->Branch("flagNCcoh", &flagNCcoh, "flagNCcoh/O");
  eventVariables->Branch("flagCC1pip", &flagCC1pip, "flagCC1pip/O");
  eventVariables->Branch("flagNC1pip", &flagNC1pip, "flagNC1pip/O");
  eventVariables->Branch("flagCC1pim", &flagCC1pim, "flagCC1pim/O");
  eventVariables->Branch("flagNC1pim", &flagNC1pim, "flagNC1pim/O");
  eventVariables->Branch("flagCC1pi0", &flagCC1pi0, "flagCC1pi0/O");
  eventVariables->Branch("flagNC1pi0", &flagNC1pi0, "flagNC1pi0/O");
};

//-------------------------------------------
// ICARUS, 1muNp0pi

void GenericFlux_Tester::AddICARUS1muNp0piVariablesToTree() {
  if (!eventVariables) {
    Config::Get().out->cd();
    eventVariables = new TTree((this->fName + "_VARS").c_str(),
                               (this->fName + "_VARS").c_str());
  }

  NUIS_LOG(SAM, "Adding ICARUS 1muNp0pi variables");

  eventVariables->Branch("ICARUS_1muNp0pi_IsSignal", &ICARUS_1muNp0pi_IsSignal, "ICARUS_1muNp0pi_IsSignal/O");
  eventVariables->Branch("ICARUS_1muNp0pi_deltaPT", &ICARUS_1muNp0pi_deltaPT, "ICARUS_1muNp0pi_deltaPT/F");
  eventVariables->Branch("ICARUS_1muNp0pi_deltaalphaT", &ICARUS_1muNp0pi_deltaalphaT, "ICARUS_1muNp0pi_deltaalphaT/F");
  eventVariables->Branch("ICARUS_1muNp0pi_MuonCos", &ICARUS_1muNp0pi_MuonCos, "ICARUS_1muNp0pi_MuonCos/F");
  eventVariables->Branch("ICARUS_1muNp0pi_MuonProtonCos", &ICARUS_1muNp0pi_MuonProtonCos, "ICARUS_1muNp0pi_MuonProtonCos/F");
  eventVariables->Branch("ICARUS_1muNp0pi_ProtonP", &ICARUS_1muNp0pi_ProtonP, "ICARUS_1muNp0pi_ProtonP/F");

}

void GenericFlux_Tester::FillICARUS1muNp0piVariablesToTree(FitEvent *event) {

  unsigned int nMu_1muNp0pi(0), nP_1muNp0pi(0), nPi_1muNp0pi(0);
  unsigned int nPhoton_1muNp0pi(0), nElectron_1muNp0pi(0), nMesons_1muNp0pi(0), nBaryonsAndPi0_1muNp0pi(0);
  double maxMomentumP_1muNp0pi = -999.;
  bool passProtonPCut_1muNp0pi = false;

  std::vector<FitParticle *> protons;

  // Start Particle Loop
  UInt_t npart = event->Npart();
  for (UInt_t i = 0; i < npart; i++) {
    // Skip particles that weren't in the final state
    bool part_alive = event->PartInfo(i)->fIsAlive and
                      event->PartInfo(i)->Status() == kFinalState;
    if (!part_alive)
      continue;

    // PDG Particle
    int pdgc = event->PartInfo(i)->fPID;
    TLorentzVector part_4mom = event->PartInfo(i)->fP;

    // ICARUS 1muNp0pi

    // All FS protons for generic purpose
    if(pdgc==2212){
      protons.push_back(event->PartInfo(i));
    }  

    double momentum = part_4mom.Vect().Mag()/1000.;

    bool PassMuonPCut = (momentum > 0.226);
    if ( abs(pdgc) == 13 ) {
      if (PassMuonPCut) nMu_1muNp0pi+=1;
    }

    if ( abs(pdgc) == 2212 ) {
      nP_1muNp0pi+=1;
      if ( momentum > maxMomentumP_1muNp0pi ) {
        maxMomentumP_1muNp0pi = momentum;
        passProtonPCut_1muNp0pi = (momentum > 0.31 && momentum < 1.);
      }
    }

    // Pion veto with momentum threshold
    if ( (abs(pdgc) == 211) && momentum > 0.087 ) nPi_1muNp0pi+=1;
    // Photon veto with momentum threshold
    if ( abs(pdgc) == 22 && part_4mom.E()/1000. > 0.025 ) nPhoton_1muNp0pi+=1;
    // Electron veto with momentum threshold
    if ( abs(pdgc) == 11 && momentum > 0.0255 ) nElectron_1muNp0pi+=1;
    else if ( abs(pdgc) == 321 || abs(pdgc) == 323 ||
              pdgc == 111 || pdgc == 130 || pdgc == 310 || pdgc == 311 ||
              pdgc == 313 || abs(pdgc) == 221 || abs(pdgc) == 331 ) nMesons_1muNp0pi+=1;
    else if ( pdgc == 3112 || pdgc == 3122 || pdgc == 3212 || pdgc == 3222 ||
              pdgc == 4112 || pdgc == 4122 || pdgc == 4212 || pdgc == 4222 ||
              pdgc == 411 || pdgc == 421 || pdgc == 111 ) nBaryonsAndPi0_1muNp0pi+=1;


  }

  ICARUS_1muNp0pi_IsSignal = nMu_1muNp0pi==1 &&
                             nP_1muNp0pi>0 && passProtonPCut_1muNp0pi &&
                             nPi_1muNp0pi==0 &&
                             nPhoton_1muNp0pi==0 &&
                             nElectron_1muNp0pi==0
                             nMesons_1muNp0pi==0 &&
                             nBaryonsAndPi0_1muNp0pi==0;

  bool IsAntiNu = event->GetNeutrinoIn()->fPID<0;

  if(! event->GetHMFSParticle(IsAntiNu ? -13 : +13) ) return;
  if(! event->GetNeutrinoIn() ) return;

  TLorentzVector Pmu = event->GetHMFSParticle(IsAntiNu ? -13 : +13)->fP;
  TLorentzVector Pnu = event->GetNeutrinoIn()->fP;

  ICARUS_1muNp0pi_deltaPT = -999.;
  ICARUS_1muNp0pi_deltaalphaT = -999.;
  ICARUS_1muNp0pi_MuonCos = -999.;
  ICARUS_1muNp0pi_MuonProtonCos = -999.;
  ICARUS_1muNp0pi_ProtonP = -999.;

  if(protons.size()>0){

    // - Sort protons in descending order of KE
    std::sort(protons.begin(), protons.end(),
              [](FitParticle* a, FitParticle* b) {
                  return a->KE() > b->KE();
              });

    ICARUS_1muNp0pi_deltaPT = FitUtils::CalcTKI_deltaPT(Pmu.Vect(), protons[0]->fP.Vect(), Pnu.Vect())/1000.;
    ICARUS_1muNp0pi_deltaalphaT = FitUtils::CalcTKI_deltaalphaT(Pmu.Vect(), protons[0]->fP.Vect(), Pnu.Vect());

    ICARUS_1muNp0pi_MuonCos = cos( Pmu.Vect().Angle( Pnu.Vect() ) );
    ICARUS_1muNp0pi_MuonProtonCos = cos( Pmu.Vect().Angle( protons[0]->fP.Vect() ) );
    ICARUS_1muNp0pi_ProtonP = protons[0]->fP.Vect().Mag()/1000.;;

  }

}

//-------------------------------------------
// ICARUS, 1mu2p0pi

void GenericFlux_Tester::AddICARUS1mu2p0piVariablesToTree() {
  if (!eventVariables) {
    Config::Get().out->cd();
    eventVariables = new TTree((this->fName + "_VARS").c_str(),
                               (this->fName + "_VARS").c_str());
  }
  NUIS_LOG(SAM, "Adding ICARUS 1mu2p0pi variables");
  eventVariables->Branch("ICARUS_1mu2p0pi_IsSignal", &ICARUS_1mu2p0pi_IsSignal, "ICARUS_1mu2p0pi_IsSignal/O");
  // - Lab-frame opening angles
  eventVariables->Branch("ICARUS_1mu2p0pi_HadronicOpeningAngle", &ICARUS_1mu2p0pi_HadronicOpeningAngle, "ICARUS_1mu2p0pi_HadronicOpeningAngle/F");
  eventVariables->Branch("ICARUS_1mu2p0pi_MuonHadronAngle", &ICARUS_1mu2p0pi_MuonHadronAngle, "ICARUS_1mu2p0pi_MuonHadronAngle/F");
  // - Single-transverse kinematic imbalance
  eventVariables->Branch("ICARUS_1mu2p0pi_DeltaPT", &ICARUS_1mu2p0pi_DeltaPT, "ICARUS_1mu2p0pi_DeltaPT/F");
  eventVariables->Branch("ICARUS_1mu2p0pi_DeltaAlphaT", &ICARUS_1mu2p0pi_DeltaAlphaT, "ICARUS_1mu2p0pi_DeltaAlphaT/F");
  eventVariables->Branch("ICARUS_1mu2p0pi_DeltaPhiT", &ICARUS_1mu2p0pi_DeltaPhiT, "ICARUS_1mu2p0pi_DeltaPhiT/F");
  // - Double-transverse kinematic imbalance
  eventVariables->Branch("ICARUS_1mu2p0pi_DeltaPTT", &ICARUS_1mu2p0pi_DeltaPTT, "ICARUS_1mu2p0pi_DeltaPTT/F");
}
void GenericFlux_Tester::FillICARUS1mu2p0piVariablesToTree(FitEvent *event) {
  unsigned int nMu_1mu2p0pi(0), nP_1mu2p0pi(0), nPi_1mu2p0pi(0);
  unsigned int nPhoton_1mu2p0pi(0), nMesons_1mu2p0pi(0), nBaryonsAndPi0_1mu2p0pi(0);
  double maxMomentumP_1mu2p0pi = -999.;
  bool passProtonMaxPCut_1mu2p0pi = false;
  std::vector<FitParticle *> protons;
  // Start Particle Loop
  UInt_t npart = event->Npart();
  for (UInt_t i = 0; i < npart; i++) {
    // Skip particles that weren't in the final state
    bool part_alive = event->PartInfo(i)->fIsAlive and
                      event->PartInfo(i)->Status() == kFinalState;
    if (!part_alive)
      continue;
    // PDG Particle
    int pdgc = event->PartInfo(i)->fPID;
    TLorentzVector part_4mom = event->PartInfo(i)->fP;
    // ICARUS 1mu2p0pi
    double momentum = part_4mom.Vect().Mag()/1000.;
    bool PassMuonPCut = (momentum > 0.226);
    if ( abs(pdgc) == 13 ) {
      if (PassMuonPCut) nMu_1mu2p0pi+=1;
    }
    if ( abs(pdgc) == 2212 && momentum > 0.35) {
      nP_1mu2p0pi+=1;
      protons.push_back(event->PartInfo(i));
      if ( momentum > maxMomentumP_1mu2p0pi ) {
        maxMomentumP_1mu2p0pi = momentum;
        passProtonMaxPCut_1mu2p0pi = (momentum < 2.);
      }
    }
    if ( abs(pdgc) == 111 || abs(pdgc) == 211 ) nPi_1mu2p0pi+=1;
    else if ( abs(pdgc) == 211 || abs(pdgc) == 321 || abs(pdgc) == 323 ||
              pdgc == 111 || pdgc == 130 || pdgc == 310 || pdgc == 311 ||
              pdgc == 313 || abs(pdgc) == 221 || abs(pdgc) == 331 ) nMesons_1mu2p0pi+=1;
    else if ( pdgc == 3112 || pdgc == 3122 || pdgc == 3212 || pdgc == 3222 ||
              pdgc == 4112 || pdgc == 4122 || pdgc == 4212 || pdgc == 4222 ||
              pdgc == 411 || pdgc == 421 || pdgc == 111 ) nBaryonsAndPi0_1mu2p0pi+=1;
  }
  ICARUS_1mu2p0pi_IsSignal = nMu_1mu2p0pi==1 &&
                             nP_1mu2p0pi>1 && passProtonMaxPCut_1mu2p0pi &&
                             nPi_1mu2p0pi==0 &&
                             nMesons_1mu2p0pi==0 &&
                             nBaryonsAndPi0_1mu2p0pi==0;
  bool IsAntiNu = event->GetNeutrinoIn()->fPID<0;
  if(! event->GetHMFSParticle(IsAntiNu ? -13 : +13) ) return;
  if(! event->GetNeutrinoIn() ) return;
  TLorentzVector Pmu = event->GetHMFSParticle(IsAntiNu ? -13 : +13)->fP;
  TLorentzVector Pnu = event->GetNeutrinoIn()->fP;
  ICARUS_1mu2p0pi_DeltaPT = -999.;
  ICARUS_1mu2p0pi_DeltaPTT = -999.;
  ICARUS_1mu2p0pi_DeltaPhiT = -999.;
  ICARUS_1mu2p0pi_DeltaAlphaT = -999.;
  ICARUS_1mu2p0pi_HadronicOpeningAngle = -999.;
  ICARUS_1mu2p0pi_MuonHadronAngle = -999.;
  if(protons.size()>1){
    // - Sort protons in descending order of KE
    std::sort(protons.begin(), protons.end(),
              [](FitParticle* a, FitParticle* b) {
                  return a->KE() > b->KE();
              });
    // - Calculate lab-frame opening angles
    TVector3 dirNu = Pnu.Vect().Unit();
    TVector3 pMu = Pmu.Vect();
    TVector3 pP1 = protons[0]->fP.Vect();
    TVector3 pP2 = protons[1]->fP.Vect();
    TVector3 pHad = pP1 + pP2;
    ICARUS_1mu2p0pi_HadronicOpeningAngle = pP1.Dot(pP2) / (pP1.Mag() * pP2.Mag());
    ICARUS_1mu2p0pi_MuonHadronAngle = pHad.Dot(pMu) / (pMu.Mag() * pHad.Mag());
    // - Calculate single-transeverse quantities
    TVector3 pTMu = pMu - (pMu.Dot(dirNu) * dirNu);
    TVector3 pTHad = pHad - (pHad.Dot(dirNu) * dirNu);
    TVector3 dpT = pTMu + pTHad;
    ICARUS_1mu2p0pi_DeltaPT = dpT.Mag()/1000.;
    ICARUS_1mu2p0pi_DeltaAlphaT = TMath::ACos( -1. * (pTMu.Dot(dpT)) / (pTMu.Mag() * dpT.Mag()) );
    ICARUS_1mu2p0pi_DeltaPhiT = TMath::ACos( -1. * (pTMu.Dot(pTHad)) / (pTMu.Mag() * pTHad.Mag()) );
    // - Calculate double-transeverse quantities
    TVector3 dirTT = ( dirNu.Cross(pMu) ).Unit();
    double pTTP1 = pP1.Dot(dirTT);
    double pTTP2 = pP2.Dot(dirTT);
    ICARUS_1mu2p0pi_DeltaPTT = (pTTP1 + pTTP2)/1000.;
  }
}

//-------------------------------------------
// SBND, 1mu1p0pi
void GenericFlux_Tester::AddSBND1mu1p0piVariablesToTree(){

  if (!eventVariables) {
    Config::Get().out->cd();
    eventVariables = new TTree((this->fName + "_VARS").c_str(),
                               (this->fName + "_VARS").c_str());
  }

  NUIS_LOG(SAM, "Adding SBND 1mu1p0pi variables");

  eventVariables->Branch("SBND_1mu1p0pi_IsSignal", &SBND_1mu1p0pi_IsSignal, "SBND_1mu1p0pi_IsSignal/O");
  eventVariables->Branch("SBND_1mu1p0pi_deltaPT", &SBND_1mu1p0pi_deltaPT, "SBND_1mu1p0pi_deltaPT/F");
  eventVariables->Branch("SBND_1mu1p0pi_deltaalphaT", &SBND_1mu1p0pi_deltaalphaT, "SBND_1mu1p0pi_deltaalphaT/F");
  eventVariables->Branch("SBND_1mu1p0pi_MuonCos", &SBND_1mu1p0pi_MuonCos, "SBND_1mu1p0pi_MuonCos/F");
  eventVariables->Branch("SBND_1mu1p0pi_MuonProtonCos", &SBND_1mu1p0pi_MuonProtonCos, "SBND_1mu1p0pi_MuonProtonCos/F");

}

void GenericFlux_Tester::FillSBND1mu1p0piVariablesToTree(FitEvent *event) {

/*
- one muon with KE > 27 MeV, energy < 1.2 GeV
- no p+- with KE 30 MeV
- one proton with KE > 50 MeV
- no pi0
*/

// (df.nmu_27MeV == 1) & (df.npi_30MeV == 0) & (df.np_50MeV == 1) & (df.npi0 == 0) & (df.mu.genE < 1.2) 

  unsigned int nmu_27MeV(0), npi_30MeV(0), np_50MeV(0), npi0(0);

  UInt_t idx_proton(0);

  // Start Particle Loop
  UInt_t npart = event->Npart();
  for (UInt_t i = 0; i < npart; i++) {
    // Skip particles that weren't in the final state
    bool part_alive = event->PartInfo(i)->fIsAlive and
                      event->PartInfo(i)->Status() == kFinalState;
    if (!part_alive)
      continue;

    // PDG Particle
    int pdgc = event->PartInfo(i)->fPID;
    TLorentzVector part_4mom = event->PartInfo(i)->fP;

    double KE = event->PartInfo(i)->KE(); // MeV

    if(abs(pdgc)==13){
      if(KE>=27.0){
        nmu_27MeV++;
      }
    }

    if(pdgc==2212){
      if(KE>=50.0){
        np_50MeV++;
        idx_proton = i;
      }
    }

    if(abs(pdgc)==211){
      if(KE>=30.0){
        npi_30MeV++;
      }
    }

    if(abs(pdgc)==111){
      npi0++;
    }

  }

  if(nmu_27MeV==0){
    SBND_1mu1p0pi_IsSignal = false;
    return;
  }
  if(np_50MeV==0){
    SBND_1mu1p0pi_IsSignal = false;
    return;
  }

  bool IsAntiNu = event->GetNeutrinoIn()->fPID<0;

  FitParticle *fp_muon = event->GetHMFSParticle(IsAntiNu ? -13 : +13);
  FitParticle *fp_proton = event->PartInfo(idx_proton);
  FitParticle *fp_nu = event->GetNeutrinoIn();

  SBND_1mu1p0pi_IsSignal = (nmu_27MeV==1) && (npi_30MeV==0) && (np_50MeV==1) && (npi0==0) && (fp_muon->E() < 1200);

  if(! fp_muon ) return;
  if(! event->GetNeutrinoIn() ) return;

  SBND_1mu1p0pi_deltaPT = FitUtils::CalcTKI_deltaPT(fp_muon->fP.Vect(), fp_proton->fP.Vect(), fp_nu->fP.Vect())/1000.;
  SBND_1mu1p0pi_deltaalphaT = FitUtils::CalcTKI_deltaalphaT(fp_muon->fP.Vect(), fp_proton->fP.Vect(), fp_nu->fP.Vect());
  SBND_1mu1p0pi_MuonCos = cos( fp_muon->fP.Vect().Angle( fp_nu->fP.Vect() ) );
  SBND_1mu1p0pi_MuonProtonCos = cos( fp_muon->fP.Vect().Angle( fp_proton->fP.Vect() ) );

}

//-------------------------------------------
// ICARUS, 1mu1pi0

void GenericFlux_Tester::AddICARUS1mu1pi0VariablesToTree() {
  if (!eventVariables) {
    Config::Get().out->cd();
    eventVariables = new TTree((this->fName + "_VARS").c_str(),
                               (this->fName + "_VARS").c_str());
  }

  NUIS_LOG(SAM, "Adding ICARUS 1mu1pi0 variables");

  eventVariables->Branch("ICARUS_1mu1pi0_IsSignal", &ICARUS_1mu1pi0_IsSignal, "ICARUS_1mu1pi0_IsSignal/O");
  eventVariables->Branch("ICARUS_1mu1pi0_MuonP", &ICARUS_1mu1pi0_MuonP, "ICARUS_1mu1pi0_MuonP/F");
  eventVariables->Branch("ICARUS_1mu1pi0_NeutralPionP", &ICARUS_1mu1pi0_NeutralPionP, "ICARUS_1mu1pi0_NeutralPionP/F");

}

void GenericFlux_Tester::FillICARUS1mu1pi0VariablesToTree(FitEvent *event) {

  // Start Particle Loop
  UInt_t npart = event->Npart();

  bool PassMuonReq{false};

  int nPi0{0};
  unsigned int Pi0Index;

  int nPipmAboveThrs{0};


  for (UInt_t i = 0; i < npart; i++) {
    // Skip particles that weren't in the final state
    bool part_alive = event->PartInfo(i)->fIsAlive and
                      event->PartInfo(i)->Status() == kFinalState;
    if (!part_alive)
      continue;

    // PDG Particle
    int pdgc = event->PartInfo(i)->fPID;
    TLorentzVector part_4mom = event->PartInfo(i)->fP;
    double part_momentum = part_4mom.Vect().Mag()/1000.; // GeV
    double part_energy = part_4mom.E()/1000.; // GeV
    double part_mass = part_4mom.M()/1000.; // GeV
    double part_ke = part_energy - part_mass; // GeV 

    // muon
    if( abs(pdgc) == 13 ){
      if( part_ke > 0.143425 ) PassMuonReq = true;
    }
    // neutral pion
    if( abs(pdgc) == 111 ){
      nPi0++;
      Pi0Index = i;
    }
    // charged pion
    if( abs(pdgc) == 211 ){
      if( part_ke>0.025 ) nPipmAboveThrs++;
    }

/*
    // CHECK A SIMILAR DEFINITION AS MINERVA FOR EXTRA REJECTION OF UNWANTED THINGS IN SIGNAL DEFN.
    // MINERvA style
    if ( abs(pdgc) == 22 && part_4mom.E()/1000. > 0.01 ) nPhoton_1muNp0pi+=1;
    else if ( abs(pdgc) == 211 || abs(pdgc) == 321 || abs(pdgc) == 323 ||
              pdgc == 111 || pdgc == 130 || pdgc == 310 || pdgc == 311 ||
              pdgc == 313 || abs(pdgc) == 221 || abs(pdgc) == 331 ) nMesons_1muNp0pi+=1;
    else if ( pdgc == 3112 || pdgc == 3122 || pdgc == 3212 || pdgc == 3222 ||
              pdgc == 4112 || pdgc == 4122 || pdgc == 4212 || pdgc == 4222 ||
              pdgc == 411 || pdgc == 421 || pdgc == 111 ) nBaryonsAndPi0_1muNp0pi+=1;
*/


  }

  ICARUS_1mu1pi0_IsSignal = PassMuonReq &&
                            (nPipmAboveThrs==0) &&
                            (nPi0==1);

  // Init
  ICARUS_1mu1pi0_MuonP = -999.;
  ICARUS_1mu1pi0_NeutralPionP = -999.;

  if(!ICARUS_1mu1pi0_IsSignal) return;

  bool IsAntiNu = event->GetNeutrinoIn()->fPID<0;
  if(! event->GetHMFSParticle(IsAntiNu ? -13 : +13) ) return;
  if(! event->GetNeutrinoIn() ) return;

  TLorentzVector Pmu = event->GetHMFSParticle(IsAntiNu ? -13 : +13)->fP;
  //TLorentzVector Pnu = event->GetNeutrinoIn()->fP;
  TLorentzVector Ppi0 = event->PartInfo(Pi0Index)->fP;

  ICARUS_1mu1pi0_MuonP = Pmu.Vect().Mag()/1000.;
  ICARUS_1mu1pi0_NeutralPionP = Ppi0.Vect().Mag()/1000.;

}



//********************************************************************
void GenericFlux_Tester::ResetVariables() {
  //********************************************************************
  // Reset neutrino PDG
  PDGnu = 0;
  // Reset energies
  Enu_true = Enu_QE = __BAD_FLOAT__;

  // Reset auxillaries
  Q2_true = Q2_QE = W_nuc_rest = bjorken_x = bjorken_y = q0_true = q3_true = Emiss = Emiss_preFSI = 
      Erecoil_true = Erecoil_charged = Erecoil_minerva = __BAD_FLOAT__;

  // Reset particle counters
  Nparticles = Nleptons = Nother = Nprotons = Nneutrons = Npiplus = Npineg =
      Npi0 = 0;

  // Reset Lepton PDG
  PDGLep = 0;
  // Reset Lepton variables
  TLep = CosLep = ELep = PLep = MLep = __BAD_FLOAT__;

  // Rset proton variables
  PPr = CosPr = EPr = TPr = MPr = __BAD_FLOAT__;

  // Reset neutron variables
  PNe = CosNe = ENe = TNe = MNe = __BAD_FLOAT__;

  // Reset pi+ variables
  PPiP = CosPiP = EPiP = TPiP = MPiP = __BAD_FLOAT__;

  // Reset pi- variables
  PPiN = CosPiN = EPiN = TPiN = MPiN = __BAD_FLOAT__;

  // Reset pi0 variables
  PPi0 = CosPi0 = EPi0 = TPi0 = MPi0 = __BAD_FLOAT__;

  // Reset the cos angles
  CosPmuPpip = CosPmuPpim = CosPmuPpi0 = CosPmuPprot = CosPmuPneut =
      CosPpipPprot = CosPpipPneut = CosPpipPpim = CosPpipPpi0 = CosPpimPprot =
          CosPpimPneut = CosPpimPpi0 = CosPi0Pprot = CosPi0Pneut =
              CosPprotPneut = __BAD_FLOAT__;
  // Reset pmiss
  pmiss.SetXYZ(-999., -999., -999.);
  pmiss_preFSI.SetXYZ(-999., -999., -999.);
}

//********************************************************************
void GenericFlux_Tester::FillEventVariables(FitEvent *event) {
  //********************************************************************

  // Fill Signal Variables
  FillSignalFlags(event);
  NUIS_LOG(DEB, "Filling signal");

  // Reset the private variables (see header)
  ResetVariables();

  // Function used to extract any variables of interest to the event
  Mode = event->Mode;
  ResCode = event->fResCode;

  // Reset the highest momentum variables
  float proton_highmom = __BAD_FLOAT__;
  float neutron_highmom = __BAD_FLOAT__;
  float piplus_highmom = __BAD_FLOAT__;
  float pineg_highmom = __BAD_FLOAT__;
  float pi0_highmom = __BAD_FLOAT__;

  (*nu_4mom) = event->PartInfo(0)->fP;

  if (!liteMode) {
    (*pmu) = TLorentzVector(0, 0, 0, 0);
    (*ppip) = TLorentzVector(0, 0, 0, 0);
    (*ppim) = TLorentzVector(0, 0, 0, 0);
    (*ppi0) = TLorentzVector(0, 0, 0, 0);
    (*pprot) = TLorentzVector(0, 0, 0, 0);
    (*pneut) = TLorentzVector(0, 0, 0, 0);
  }

  Enu_true = nu_4mom->E();
  PDGnu = event->PartInfo(0)->fPID;

  bool cc = (abs(event->Mode) < 30);
  (void)cc;

  // Add all pion distributions for the event.

  // Add classifier for CC0pi or CC1pi or CCOther
  // Save Modes Properly
  // Save low recoil measurements

  // Start Particle Loop
  UInt_t npart = event->Npart();
  for (UInt_t i = 0; i < npart; i++) {
    // Skip particles that weren't in the final state
    bool part_alive = event->PartInfo(i)->fIsAlive and
                      event->PartInfo(i)->Status() == kFinalState;
    if (!part_alive)
      continue;

    // PDG Particle
    int PDGpart = event->PartInfo(i)->fPID;
    TLorentzVector part_4mom = event->PartInfo(i)->fP;

    Nparticles++;

    // Get Charged Lepton
    if (abs(PDGpart) == abs(PDGnu) - 1) {
      Nleptons++;

      PDGLep = PDGpart;

      TLep = FitUtils::T(part_4mom) * 1000.0;
      PLep = (part_4mom.Vect().Mag());
      ELep = (part_4mom.E());
      MLep = (part_4mom.Mag());
      CosLep = cos(part_4mom.Vect().Angle(nu_4mom->Vect()));
      (*pmu) = part_4mom;

      Q2_true = -1 * (part_4mom - (*nu_4mom)).Mag2();

      float ThetaLep = (event->PartInfo(0))
                           ->fP.Vect()
                           .Angle((event->PartInfo(i))->fP.Vect());

      q0_true = (part_4mom - (*nu_4mom)).E();
      q3_true = (part_4mom - (*nu_4mom)).Vect().Mag();

      Emiss = FitUtils::GetEmiss(event);
      pmiss = FitUtils::GetPmiss(event);

      Emiss_preFSI = FitUtils::GetEmiss(event, 1);
      pmiss_preFSI = FitUtils::GetPmiss(event, 1);

      // Get W_true with assumption of initial state nucleon at rest
      float m_n = (float)PhysConst::mass_proton * 1000.;
      W_nuc_rest = sqrt(-Q2_true + 2 * m_n * (Enu_true - ELep) + m_n * m_n);

      // Get the Bjorken x and y variables
      // Assume that E_had = Enu - Emu as in MINERvA
      bjorken_x = Q2_true / (2 * m_n * (Enu_true - ELep));
      bjorken_y = 1 - ELep / Enu_true;

      // Quasi-elastic ----------------------
      // ------------------------------------

      // Q2 QE Assuming Carbon Input. Should change this to be dynamic soon.
      Q2_QE =
          FitUtils::Q2QErec(part_4mom, cos(ThetaLep), 34., true) * 1000000.0;
      Enu_QE = FitUtils::EnuQErec(part_4mom, cos(ThetaLep), 34., true) * 1000.0;

      // Pion Production ----------------------
      // --------------------------------------

    } else if (PDGpart == 2212) {
      Nprotons++;
      if (part_4mom.Vect().Mag() > proton_highmom) {
        proton_highmom = part_4mom.Vect().Mag();

        PPr = (part_4mom.Vect().Mag());
        EPr = (part_4mom.E());
        TPr = FitUtils::T(part_4mom) * 1000.;
        MPr = (part_4mom.Mag());
        CosPr = cos(part_4mom.Vect().Angle(nu_4mom->Vect()));

        (*pprot) = part_4mom;
      }
    } else if (PDGpart == 2112) {
      Nneutrons++;
      if (part_4mom.Vect().Mag() > neutron_highmom) {
        neutron_highmom = part_4mom.Vect().Mag();

        PNe = (part_4mom.Vect().Mag());
        ENe = (part_4mom.E());
        TNe = FitUtils::T(part_4mom) * 1000.;
        MNe = (part_4mom.Mag());
        CosNe = cos(part_4mom.Vect().Angle(nu_4mom->Vect()));

        (*pneut) = part_4mom;
      }
    } else if (PDGpart == 211) {
      Npiplus++;
      if (part_4mom.Vect().Mag() > piplus_highmom) {
        piplus_highmom = part_4mom.Vect().Mag();

        PPiP = (part_4mom.Vect().Mag());
        EPiP = (part_4mom.E());
        TPiP = FitUtils::T(part_4mom) * 1000.;
        MPiP = (part_4mom.Mag());
        CosPiP = cos(part_4mom.Vect().Angle(nu_4mom->Vect()));

        (*ppip) = part_4mom;
      }
    } else if (PDGpart == -211) {
      Npineg++;
      if (part_4mom.Vect().Mag() > pineg_highmom) {
        pineg_highmom = part_4mom.Vect().Mag();

        PPiN = (part_4mom.Vect().Mag());
        EPiN = (part_4mom.E());
        TPiN = FitUtils::T(part_4mom) * 1000.;
        MPiN = (part_4mom.Mag());
        CosPiN = cos(part_4mom.Vect().Angle(nu_4mom->Vect()));

        (*ppim) = part_4mom;
      }
    } else if (PDGpart == 111) {
      Npi0++;
      if (part_4mom.Vect().Mag() > pi0_highmom) {
        pi0_highmom = part_4mom.Vect().Mag();

        PPi0 = (part_4mom.Vect().Mag());
        EPi0 = (part_4mom.E());
        TPi0 = FitUtils::T(part_4mom) * 1000.;
        MPi0 = (part_4mom.Mag());
        CosPi0 = cos(part_4mom.Vect().Angle(nu_4mom->Vect()));

        (*ppi0) = part_4mom;
      }
    } else {
      Nother++;
    }
  }

  // Get Recoil Definitions ------
  // -----------------------------
  Erecoil_true = FitUtils::GetErecoil_TRUE(event);
  Erecoil_charged = FitUtils::GetErecoil_CHARGED(event);
  Erecoil_minerva = FitUtils::GetErecoil_MINERvA_LowRecoil(event);

  // Do the angles between final state particles
  if (Nleptons > 0 && Npiplus > 0)
    CosPmuPpip = cos(pmu->Vect().Angle(ppip->Vect()));
  if (Nleptons > 0 && Npineg > 0)
    CosPmuPpim = cos(pmu->Vect().Angle(ppim->Vect()));
  if (Nleptons > 0 && Npi0 > 0)
    CosPmuPpi0 = cos(pmu->Vect().Angle(ppi0->Vect()));
  if (Nleptons > 0 && Nprotons > 0)
    CosPmuPprot = cos(pmu->Vect().Angle(pprot->Vect()));
  if (Nleptons > 0 && Nneutrons > 0)
    CosPmuPneut = cos(pmu->Vect().Angle(pneut->Vect()));

  if (Npiplus > 0 && Nprotons > 0)
    CosPpipPprot = cos(ppip->Vect().Angle(pprot->Vect()));
  if (Npiplus > 0 && Nneutrons > 0)
    CosPpipPneut = cos(ppip->Vect().Angle(pneut->Vect()));
  if (Npiplus > 0 && Npineg > 0)
    CosPpipPpim = cos(ppip->Vect().Angle(ppim->Vect()));
  if (Npiplus > 0 && Npi0 > 0)
    CosPpipPpi0 = cos(ppip->Vect().Angle(ppi0->Vect()));

  if (Npineg > 0 && Nprotons > 0)
    CosPpimPprot = cos(ppim->Vect().Angle(pprot->Vect()));
  if (Npineg > 0 && Nneutrons > 0)
    CosPpimPneut = cos(ppim->Vect().Angle(pneut->Vect()));
  if (Npineg > 0 && Npi0 > 0)
    CosPpimPpi0 = cos(ppim->Vect().Angle(ppi0->Vect()));

  if (Npi0 > 0 && Nprotons > 0)
    CosPi0Pprot = cos(ppi0->Vect().Angle(pprot->Vect()));
  if (Npi0 > 0 && Nneutrons > 0)
    CosPi0Pneut = cos(ppi0->Vect().Angle(pneut->Vect()));

  if (Nprotons > 0 && Nneutrons > 0)
    CosPprotPneut = cos(pprot->Vect().Angle(pneut->Vect()));

  if(Fill_ICARUS_QELike_Variable){
    FillICARUS1muNp0piVariablesToTree(event);
    FillICARUS1mu2p0piVariablesToTree(event);
  }

  if(Fill_SBND_QELike_Variable){
    FillSBND1mu1p0piVariablesToTree(event);
  }

  if(Fill_ICARUS_1mu1pi0_Variable){
    FillICARUS1mu1pi0VariablesToTree(event);
  }

  // Event Weights ----
  // ------------------
  Weight = event->RWWeight * event->InputWeight;
  RWWeight = event->RWWeight;
  InputWeight = event->InputWeight;
  FluxWeight =
      GetFluxHistogram()->GetBinContent(GetFluxHistogram()->FindBin(Enu)) /
      GetFluxHistogram()->Integral();

  // Fill the eventVariables Tree
  eventVariables->Fill();
  return;
};

//********************************************************************
void GenericFlux_Tester::Write(std::string drawOpt) {
  //********************************************************************

  // First save the TTree
  eventVariables->Write();

  // Save Flux and Event Histograms too
  GetInput()->GetFluxHistogram()->Write();
  GetInput()->GetEventHistogram()->Write();

  return;
}

//********************************************************************
void GenericFlux_Tester::FillSignalFlags(FitEvent *event) {
  //********************************************************************

  // Some example flags are given from SignalDef.
  // See src/Utils/SignalDef.cxx for more.
  int nuPDG = event->PartInfo(0)->fPID;

  // Generic signal flags
  flagCCINC = SignalDef::isCCINC(event, nuPDG);
  flagNCINC = SignalDef::isNCINC(event, nuPDG);
  flagCCQE = SignalDef::isCCQE(event, nuPDG);
  flagCCQELike = SignalDef::isCCQELike(event, nuPDG);
  flagCC0pi = SignalDef::isCC0pi(event, nuPDG);
  flagNCEL = SignalDef::isNCEL(event, nuPDG);
  flagNC0pi = SignalDef::isNC0pi(event, nuPDG);
  flagCCcoh = SignalDef::isCCCOH(event, nuPDG, 211);
  flagNCcoh = SignalDef::isNCCOH(event, nuPDG, 111);
  flagCC1pip = SignalDef::isCC1pi(event, nuPDG, 211);
  flagNC1pip = SignalDef::isNC1pi(event, nuPDG, 211);
  flagCC1pim = SignalDef::isCC1pi(event, nuPDG, -211);
  flagNC1pim = SignalDef::isNC1pi(event, nuPDG, -211);
  flagCC1pi0 = SignalDef::isCC1pi(event, nuPDG, 111);
  flagNC1pi0 = SignalDef::isNC1pi(event, nuPDG, 111);
}

// -------------------------------------------------------------------
// Purely MC Plot
// Following functions are just overrides to handle this
// -------------------------------------------------------------------
//********************************************************************
/// Everything is classed as signal...
bool GenericFlux_Tester::isSignal(FitEvent *event) {
  //********************************************************************
  (void)event;
  return true;
};

//********************************************************************
void GenericFlux_Tester::ScaleEvents() {
  //********************************************************************
  // Saving everything to a TTree so no scaling required
  return;
}

//********************************************************************
void GenericFlux_Tester::ApplyNormScale(float norm) {
  //********************************************************************

  // Saving everything to a TTree so no scaling required
  this->fCurrentNorm = norm;
  return;
}

//********************************************************************
void GenericFlux_Tester::FillHistograms() {
  //********************************************************************
  // No Histograms need filling........
  return;
}

//********************************************************************
void GenericFlux_Tester::ResetAll() {
  //********************************************************************
  eventVariables->Reset();
  return;
}

//********************************************************************
float GenericFlux_Tester::GetChi2() {
  //********************************************************************
  // No Likelihood to test, purely MC
  return 0.0;
}
