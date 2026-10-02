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

#ifndef GenericFlux_Tester_H_SEEN
#define GenericFlux_Tester_H_SEEN
#include "Measurement1D.h"

#ifndef __BAD__FLOAT__
#define __BAD_FLOAT__ -999.99
#endif

//********************************************************************
class GenericFlux_Tester : public Measurement1D {
//********************************************************************

public:

  GenericFlux_Tester(std::string name, std::string inputfile, FitWeight *rw, std::string type, std::string fakeDataFile);
  virtual ~GenericFlux_Tester() {};

  //! Clear private variables
  inline void ResetVariables();

  //! Grab info from event
  void FillEventVariables(FitEvent *event);

  //! Fill Custom Histograms
  void FillHistograms();

  //! ResetAll
  void ResetAll();

  //! Scale
  void ScaleEvents();

  //! Norm
  void ApplyNormScale(float norm);

  //! Define this samples signal
  bool isSignal(FitEvent *nvect);

  //! Write Files
  void Write(std::string drawOpt);

  //! Get Chi2
  float GetChi2();

  //! Fill all signal flags we currently have
  void FillSignalFlags(FitEvent *event);

  void AddEventVariablesToTree();
  void AddSignalFlagsToTree();

  void AddICARUS1muNp0piVariablesToTree();
  void FillICARUS1muNp0piVariablesToTree(FitEvent *event);
  void AddICARUS1mu2p0piVariablesToTree();
  void FillICARUS1mu2p0piVariablesToTree(FitEvent *event);

  void AddSBND1mu1p0piVariablesToTree();
  void FillSBND1mu1p0piVariablesToTree(FitEvent *event);

  void AddICARUS1mu1pi0VariablesToTree();
  void FillICARUS1mu1pi0VariablesToTree(FitEvent *event);

 private:

  // Lighter flat trees that don't include vectors
  bool liteMode;

  TTree* eventVariables;

  TLorentzVector *nu_4mom;
  TLorentzVector *pmu;
  TLorentzVector *ppip;
  TLorentzVector *ppim;
  TLorentzVector *ppi0;
  TLorentzVector *pprot;
  TLorentzVector *pneut;

  // Saved Variables
  float Enu_true;
  float Enu_QE;
  int PDGnu;

  // Auxillairies
  float Q2_true;
  float Q2_QE;
  float W_nuc_rest;
  float bjorken_x;
  float bjorken_y;
  float q0_true;
  float q3_true;
  float Emiss;
  float Emiss_preFSI;
  TVector3 pmiss;
  TVector3 pmiss_preFSI;
  float Erecoil_true;
  float Erecoil_charged;
  float Erecoil_minerva;

  // Interaction mode
  int Mode;

  // Res ID (for GENIE only)
  int ResCode;

  // Particle counters
  int Nparticles;
  int Nleptons;
  int Nother;
  int Nprotons;
  int Nneutrons;
  int Npiplus;
  int Npineg;
  int Npi0;

  // Lepton variables
  int PDGLep;
  float TLep;
  float CosLep;
  float ELep;
  float PLep;
  float MLep;

  // Proton variables
  float PPr;  //!< Highest Mom Proton
  float CosPr; //!< Highest Mom Proton
  float EPr;
  float TPr;
  float MPr;

  // Neutron variables
  float PNe;
  float CosNe;
  float ENe;
  float TNe;
  float MNe;

  // Pi+ variables
  float PPiP;
  float CosPiP;
  float EPiP;
  float TPiP;
  float MPiP;

  // Pi- variables
  float PPiN;
  float CosPiN;
  float EPiN;
  float TPiN;
  float MPiN;

  float PPi0;
  float CosPi0;
  float EPi0;
  float TPi0;
  float MPi0;

  // Angular variables
  float CosPmuPpip;
  float CosPmuPpim;
  float CosPmuPpi0;
  float CosPmuPprot;
  float CosPmuPneut;
  float CosPpipPprot;
  float CosPpipPneut;
  float CosPpipPpim;
  float CosPpipPpi0;
  float CosPpimPprot;
  float CosPpimPneut;
  float CosPpimPpi0;
  float CosPi0Pprot;
  float CosPi0Pneut;
  float CosPprotPneut;

  // Weights
  float Weight;
  float RWWeight;
  float InputWeight;
  float FluxWeight;

  // Generic signal flags
  bool flagCCINC;
  bool flagNCINC;
  bool flagCCQE;
  bool flagCC0pi;
  bool flagCCQELike;
  bool flagNCEL;
  bool flagNC0pi;
  bool flagCCcoh;
  bool flagNCcoh;
  bool flagCC1pip;
  bool flagNC1pip;
  bool flagCC1pim;
  bool flagNC1pim;
  bool flagCC1pi0;
  bool flagNC1pi0;

  // ICARUS 1muNp0pi variables
  bool Fill_ICARUS_QELike_Variable;
  // - Tree variables
  // 1) 1muNp0pi
  bool ICARUS_1muNp0pi_IsSignal;
  bool ICARUS_1muNp0pi_IsSignal_Howard;
  bool ICARUS_1muNp0pi_IsSignal_muonChanged;
  bool ICARUS_1muNp0pi_IsSignal_protonChanged;
  float ICARUS_1muNp0pi_deltaPT;
  float ICARUS_1muNp0pi_deltaalphaT;
  float ICARUS_1muNp0pi_MuonCos;
  float ICARUS_1muNp0pi_MuonProtonCos;
  float ICARUS_1muNp0pi_ProtonP;
  // 2) 1mu(N<1)p0pi
  bool ICARUS_1mu2p0pi_IsSignal;
  float ICARUS_1mu2p0pi_HadronicOpeningAngle;
  float ICARUS_1mu2p0pi_MuonHadronAngle;
  float ICARUS_1mu2p0pi_DeltaPT;
  float ICARUS_1mu2p0pi_DeltaAlphaT;
  float ICARUS_1mu2p0pi_DeltaPhiT;
  float ICARUS_1mu2p0pi_DeltaPTT;

  // SBND 1mu1p0pi
  bool Fill_SBND_QELike_Variable;
  // - Tree varaibles
  // 1) 1mu1p0pi
  bool SBND_1mu1p0pi_IsSignal;
  float SBND_1mu1p0pi_deltaPT;
  float SBND_1mu1p0pi_deltaalphaT;
  float SBND_1mu1p0pi_MuonCos;
  float SBND_1mu1p0pi_MuonProtonCos;
  
  // ICARUS 1mu1pi0
  bool Fill_ICARUS_1mu1pi0_Variable;
  // - Tree variables
  // 1) Lane's
  bool ICARUS_1mu1pi0_IsSignal;
  float ICARUS_1mu1pi0_MuonP;
  float ICARUS_1mu1pi0_NeutralPionP;

};

#endif
