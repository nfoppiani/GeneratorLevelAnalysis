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

  unsigned int _run, _subrun, _event;
  double _bnbweight;

  unsigned int _n_flash_simple, _n_flash_simple_over50, _n_flash_simple_beam, _n_flash_simple_over50_beam;
  std::vector<double> _flash_PE_simple, _flash_time_simple, _flash_y_simple, _flash_z_simple;
  unsigned int _n_flash_op, _n_flash_op_over50, _n_flash_op_beam, _n_flash_op_over50_beam;
  std::vector<double> _flash_PE_op, _flash_time_op, _flash_y_op, _flash_z_op;

  std::map<int, int> _pdg_daughter = {{12, 11}, {14, 13}, {-12, -11}, {-14, -13}};
};
