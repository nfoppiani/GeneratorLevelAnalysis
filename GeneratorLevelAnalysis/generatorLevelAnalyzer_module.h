////////////////////////////////////////////////////////////////////////
// Class:       generatorLevelAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        generatorLevelAnalyzer_module.cc
//
// Generated at Wed Feb 14 10:17:20 2018 by Nicolo Foppiani using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <memory>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "TTree.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfo/DetectorProperties.h"

#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "uboone/EventWeight/MCEventWeight.h"

#include "PandoraInterfaceHelper.h"

class generatorLevelAnalyzer;

class generatorLevelAnalyzer : public art::EDAnalyzer
{
public:
  explicit generatorLevelAnalyzer(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  generatorLevelAnalyzer(generatorLevelAnalyzer const &) = delete;
  generatorLevelAnalyzer(generatorLevelAnalyzer &&) = delete;
  generatorLevelAnalyzer & operator = (generatorLevelAnalyzer const &) = delete;
  generatorLevelAnalyzer & operator = (generatorLevelAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const &evt) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const &p) override;
  void respondToOpenInputFile(art::FileBlock const &fb) override;
  void endSubRun(art::SubRun const &sr) override;
  void clear();

  /// Additional functions
  void storeBNBWeight(art::Event const &evt);
  void opticalInformation(art::Event const &evt);
  void PNCandidatesInformation(art::Event const &evt);
  void mcFluxInformation(art::Event const &evt);
  void trueNeutrinoInformation(art::Event const &evt);

private:
  std::string _mctruthLabel = "generator";
  std::string _mcparticleLabel = "largeant";

  lee::PandoraInterfaceHelper pandoraHelper;

  bool m_isOverlaidSample;
  bool m_isData;
  double m_beamStart;
  double m_beamEnd;
  std::string m_opticalFlashFinderLabel;
  std::string m_pfp_producer;

  TTree *myPOTTTree;
  std::ofstream _run_subrun_list_file;
  unsigned int _run_sr, _subrun_sr;
  double _pot, _sum_pot;

  TTree *myTTree;

  double _bnbweight;

  unsigned int _n_flash_simple, _n_flash_simple_over50, _n_flash_simple_beam, _n_flash_simple_over50_beam;
  std::vector<double> _flash_PE_simple, _flash_time_simple, _flash_y_simple, _flash_z_simple;
  unsigned int _n_flash_op, _n_flash_op_over50, _n_flash_op_beam, _n_flash_op_over50_beam;
  std::vector<double> _flash_PE_op, _flash_time_op, _flash_y_op, _flash_z_op;

  unsigned int _n_total_candidates;
  std::vector<double> _nu_candidate_vx, _nu_candidate_vy, _nu_candidate_vz;
  std::vector<int> _n_daughters_candidate;

  int _n_mcfluxes; //number of MC Fluxes in the handle
  double _dk_x; //x-position of decay point [cm]
  double _dk_y; //y-position of decay point [cm]
  double _dk_z; //z-position of decay point [cm]
  double _dk_px; //neutrino parent x-momentum at decay point [GeV]
  double _dk_py; //neutrino parent y-momentum at decay point [GeV]
  double _dk_pz; //neutrino parent z-momentum at decay point [GeV]
  int _dk_pdg; //neutrino parent PDG ID
  double _gen_x; //x-position of origin of ray from flux generator [m]
  double _gen_y; //y-position of origin of ray from flux generator [m]
  double _gen_z; //z-position of origin of ray from flux generator [m]
  double _dk2gen; //distance between decay point and origin of ray from flux generator [m]
  double _gen2vtx; //distance between origin of ray from flux generator and neutrino interaction point [m]
  double _totaltime; //total time between target and neutrino interaction [microseconds]
  double _delaytime; //difference between total time and time of prompt neutrino, i.e. travelling at light speed from the target to MicroBooNE's position [microseconds]

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
  unsigned int _n_true_daughter_candidates;
  double _true_daughter_E;
  double _true_daughter_theta;
  double _true_daughter_phi;
  double _true_daughter_T;

  std::map<int, int> _pdg_daughter = {{12, 11}, {14, 13}, {-12, -11}, {-14, -13}};
};
