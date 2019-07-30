////////////////////////////////////////////////////////////////////////
// Class:       generatorLevelAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        generatorLevelAnalyzer_module.cc
//
// Generated at Wed Feb 14 10:17:20 2018 by Nicolo Foppiani using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "generatorLevelAnalyzer_module.h"

generatorLevelAnalyzer::generatorLevelAnalyzer(fhicl::ParameterSet const &p)
    : EDAnalyzer(p)
{
  art::ServiceHandle<art::TFileService> tfs;
  myPOTTTree = tfs->make<TTree>("pot", "POT Tree");
  myPOTTTree->Branch("pot", &_pot, "pot/D");
  myPOTTTree->Branch("run", &_run_sr, "run/i");
  myPOTTTree->Branch("subrun", &_subrun_sr, "subrun/i");

  myTTree = tfs->make<TTree>("tree", "Tree");
  myTTree->Branch("event", &_event, "event/i");
  myTTree->Branch("run", &_run, "run/i");
  myTTree->Branch("subrun", &_subrun, "subrun/i");

  myTTree->Branch("n_flash_simple", &_n_flash_simple, "n_flash_simple/I");
  myTTree->Branch("n_flash_simple_over50", &_n_flash_simple_over50, "n_flash_simple_over50/I");
  myTTree->Branch("n_flash_simple_beam", &_n_flash_simple_beam, "n_flash_simple_beam/I");
  myTTree->Branch("n_flash_simple_over50_beam", &_n_flash_simple_over50_beam, "n_flash_simple_over50_beam/I");
  myTTree->Branch("flash_time_simple", "std::vector< double >", &_flash_time_simple);
  myTTree->Branch("flash_PE_simple", "std::vector< double >", &_flash_PE_simple);
  myTTree->Branch("flash_y_simple", "std::vector< double >", &_flash_y_simple);
  myTTree->Branch("flash_z_simple", "std::vector< double >", &_flash_z_simple);

  // myTTree->Branch("n_flash_op", &_n_flash_op, "n_flash_op/I");
  // myTTree->Branch("n_flash_op_over50", &_n_flash_op_over50, "n_flash_op_over50/I");
  // myTTree->Branch("n_flash_op_beam", &_n_flash_op_beam, "n_flash_op_beam/I");
  // myTTree->Branch("n_flash_op_over50_beam", &_n_flash_op_over50_beam, "n_flash_op_over50_beam/I");
  // myTTree->Branch("flash_time_op", "std::vector< double >", &_flash_time_op);
  // myTTree->Branch("flash_PE_op", "std::vector< double >", &_flash_PE_op);
  // myTTree->Branch("flash_y_op", "std::vector< double >", &_flash_y_op);
  // myTTree->Branch("flash_z_op", "std::vector< double >", &_flash_z_op);

  _run_subrun_list_file.open("run_subrun_list_filter.txt", std::ofstream::out | std::ofstream::trunc);

  this->reconfigure(p);
}

void generatorLevelAnalyzer::reconfigure(fhicl::ParameterSet const &p)
{
  m_opticalFlashFinderLabel = p.get<std::string>("OpticalFlashFinderLabel", "simpleFlashBeam");
  m_pfp_producer = p.get<std::string>("PFParticleProducer", "pandoraNu");
  m_isOverlaidSample = p.get<bool>("isOverlaidSample", false);
  m_isData = p.get<bool>("isData");
  m_beamStart = p.get<double>("beamStart");
  m_beamEnd = p.get<double>("beamEnd");
}

void generatorLevelAnalyzer::respondToOpenInputFile(art::FileBlock const &fb)
{
  _sum_pot = 0;
}

void generatorLevelAnalyzer::endSubRun(art::SubRun const &sr)
{
  _run_subrun_list_file << sr.run() << " " << sr.subRun() << std::endl;

  _run_sr = sr.run();
  _subrun_sr = sr.subRun();

  art::Handle<sumdata::POTSummary> potListHandle;

  if (!m_isData || m_isOverlaidSample)
  {
    if (sr.getByLabel(_mctruthLabel, potListHandle)) {
      _pot = potListHandle->totpot - _sum_pot;
      if (m_isOverlaidSample) {
          _sum_pot += _pot;
      }
      std::cout << "Subrun POT " << _pot << " " << potListHandle->totpot << " " << _sum_pot << std::endl;
    } else {
        _pot = 0.;
    }
  }
  else
  {
    if (sr.getByLabel("beamdata", "bnbETOR860", potListHandle)) {
      _pot = potListHandle->totpot;
      std::cout << "Subrun POT " << _pot << std::endl;
    } else {
      _pot = 0.;
    }
  }

  myPOTTTree->Fill();
  return;
}

void generatorLevelAnalyzer::clear()
{
  // pot information
  _pot = std::numeric_limits<double>::lowest();
  _run_sr = std::numeric_limits<unsigned int>::lowest();
  _subrun_sr = std::numeric_limits<unsigned int>::lowest();

  _bnbweight = 1;

  //optical optical information
  _n_flash_simple = 0;
  _n_flash_simple_over50 = 0;
  _n_flash_simple_beam = 0;
  _n_flash_simple_over50_beam = 0;
  _flash_PE_simple.clear();
  _flash_time_simple.clear();
  _flash_y_simple.clear();
  _flash_z_simple.clear();

  // _n_flash_op = 0;
  // _n_flash_op_over50 = 0;
  // _n_flash_op_beam = 0;
  // _n_flash_op_over50_beam = 0;
  // _flash_PE_op.clear();
  // _flash_time_op.clear();
  // _flash_y_op.clear();
  // _flash_z_op.clear();
}

void generatorLevelAnalyzer::opticalInformation(art::Event const &evt)
{
  art::InputTag optical_tag_simple{"simpleFlashBeam"};
  auto const &optical_handle_simple = evt.getValidHandle<std::vector<recob::OpFlash>>(optical_tag_simple);

  for (unsigned int ifl = 0; ifl < optical_handle_simple->size(); ++ifl)
  {
    recob::OpFlash const &flash = optical_handle_simple->at(ifl);
    _flash_PE_simple.push_back(flash.TotalPE());
    _flash_time_simple.push_back(flash.Time());
    _flash_y_simple.push_back(flash.YCenter());
    _flash_z_simple.push_back(flash.ZCenter());

    _n_flash_simple ++;
    if (flash.TotalPE() > 50) _n_flash_simple_over50 ++;
    if (flash.Time() > m_beamStart && flash.Time() < m_beamEnd)
    {
      _n_flash_simple_beam ++;
      if (flash.TotalPE() > 50) _n_flash_simple_over50_beam ++;
    }
  }

  // art::InputTag optical_tag_op{"opflashBeam"};
  // auto const &optical_handle_op = evt.getValidHandle<std::vector<recob::OpFlash>>(optical_tag_op);
  //
  // for (unsigned int ifl = 0; ifl < optical_handle_op->size(); ++ifl)
  // {
  //   recob::OpFlash const &flash = optical_handle_op->at(ifl);
  //   _flash_PE_op.push_back(flash.TotalPE());
  //   _flash_time_op.push_back(flash.Time());
  //   _flash_y_op.push_back(flash.YCenter());
  //   _flash_z_op.push_back(flash.ZCenter());
  //
  //   _n_flash_op ++;
  //   if (flash.TotalPE() > 50) _n_flash_op_over50 ++;
  //   if (flash.Time() > m_beamStart && flash.Time() < m_beamEnd)
  //   {
  //     _n_flash_op_beam ++;
  //     if (flash.TotalPE() > 50) _n_flash_op_over50_beam ++;
  //   }
  // }
}

void generatorLevelAnalyzer::analyze(art::Event const &evt)
{
  clear();
  std::cout << "[GENERATOR LEVEL ANALYZER] : RUN " << evt.run()
            << " SUBRUN " << evt.subRun()
            << " EVENT " << evt.id().event()
            << std::endl;

  _event = evt.event();
  _run = evt.run();
  _subrun = evt.subRun();

  // optical information
  opticalInformation(evt);

  myTTree->Fill();
}

DEFINE_ART_MODULE(generatorLevelAnalyzer)
