////////////////////////////////////////////////////////////////////////
// Class:       generator_level_analyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        generator_level_analyzer_module.cc
//
// Generated at Wed Feb 14 10:17:20 2018 by Nicolo Foppiani using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <memory>

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ElectronEventSelectionAlg.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "TTree.h"
#include "TH1F.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfo/DetectorProperties.h"


class generator_level_analyzer;

class generator_level_analyzer : public art::EDAnalyzer
{
public:
  explicit ElectronNeutrinoFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ElectronNeutrinoFilter(ElectronNeutrinoFilter const &) = delete;
  ElectronNeutrinoFilter(ElectronNeutrinoFilter &&) = delete;
  ElectronNeutrinoFilter & operator = (ElectronNeutrinoFilter const &) = delete;
  ElectronNeutrinoFilter & operator = (ElectronNeutrinoFilter &&) = delete;

  // Required functions.
  bool analyze(art::Event & evt) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToOpenInputFile(art::FileBlock const &fb) override;
  bool endSubRun(art::SubRun &sr) override;
  void clear();

  /// Additional functions
  void storeBNBWeight(art::Event const &evt);
  void opticalInformation(art::Event const &evt);
  void PNCandidatesInformation(art::Event const &evt);
  void trueNeutrinoInformation(art::Event const &evt);

private:
  std::string _mctruthLabel = "generator";
  std::string _mcparticleLabel = "largeant";

  TTree *myPOTTTree;
  std::ofstream _run_subrun_list_file;
  unsigned int _run_sr, _subrun_sr;
  double _pot, _sum_pot;

  TTree *myTTree;

  double _bnbweight;

  std::vector<double> _flash_PE, _flash_time;

  unsigned int _n_total_candidates;
  std::vector<double> _nu_candidate_vx, _nu_candidate_vy, _nu_candidate_vz;

  unsigned int _n_true_nu;
  int _nu_pdg;
  double _nu_E;
  double _nu_theta;
  double _nu_phi;
  double _nu_T;
  double _nu_pt;
  double _nu_qsqr;
  double _nu_w;
  int _ccnc;
  double _true_vx;
  double _true_vy;
  double _true_vz;
  unsigned int _interaction_type;

  std::vector<double> _nu_daughters_E;
  std::vector<int> _nu_daughters_pdg;
  unsigned int _n_true_pions;
  unsigned int _n_true_protons;
  unsigned int _n_true_protons_above40;
  unsigned int _n_true_protons_above21;
  std::vector < std::vector<double> > _nu_daughters_p;
  std::vector < std::vector<double> > _nu_daughters_start_v;
  std::vector < std::vector<double> > _nu_daughters_end_v;
  double _true_daughter_E;
  double _true_daughter_theta;
  double _true_daughter_phi;
  double _true_daughter_T;
};
