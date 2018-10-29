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

  myTTree = tfs->make<TTree>("filtertree", "Filter Tree");
  myTTree->Branch("bnbweight", &_bnbweight, "bnbweight/D");
  myTTree->Branch("flash_time", "std::vector< double >", &_flash_time);
  myTTree->Branch("flash_PE", "std::vector< double >", &_flash_PE);

  myTTree->Branch("n_total_candidates", &_n_total_candidates, "n_total_candidates/i");
    
  myTTree->Branch("dk_x", &_dk_x, "dk_x/D");
  myTTree->Branch("dk_y", &_dk_y, "dk_y/D");
  myTTree->Branch("dk_z", &_dk_z, "dk_z/D");
  myTTree->Branch("dk_px", &_dk_px, "dk_px/D");
  myTTree->Branch("dk_py", &_dk_py, "dk_py/D");
  myTTree->Branch("dk_pz", &_dk_pz, "dk_pz/D");
  myTTree->Branch("dk_pdg", &_dk_pdg, "dk_pdg/I");
  myTTree->Branch("gen_x", &_gen_x, "gen_x/D");
  myTTree->Branch("gen_y", &_gen_y, "gen_y/D");
  myTTree->Branch("gen_z", &_gen_z, "gen_z/D");
  myTTree->Branch("dk2gen", &_dk2gen, "dk2gen/D");
  myTTree->Branch("gen2vtx", &_gen2vtx, "gen2vtx/D");
  myTTree->Branch("totaltime", &_totaltime, "totaltime/D");
  myTTree->Branch("delaytime", &_delaytime, "delaytime/D");
  
  myTTree->Branch("nu_candidate_vx", "std::vector< double >", &_nu_candidate_vx);
  myTTree->Branch("nu_candidate_vy", "std::vector< double >", &_nu_candidate_vy);
  myTTree->Branch("nu_candidate_vz", "std::vector< double >", &_nu_candidate_vz);

  myTTree->Branch("n_true_nu", &_n_true_nu, "n_true_nu/i");
  myTTree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  myTTree->Branch("nu_E", &_nu_E, "nu_E/D");
  myTTree->Branch("nu_theta", &_nu_theta, "nu_theta/D");
  myTTree->Branch("nu_phi", &_nu_phi, "nu_phi/D");
  myTTree->Branch("nu_T", &_nu_T, "nu_T/D");
  myTTree->Branch("nu_pt", &_nu_pt, "nu_pt/D");
  myTTree->Branch("nu_qsqr", &_nu_qsqr, "nu_qsqr/D");
  myTTree->Branch("nu_w", &_nu_w, "nu_w/D");
  myTTree->Branch("ccnc", &_ccnc, "ccnc/I");
  myTTree->Branch("true_vx", &_true_vx, "true_vx/D");
  myTTree->Branch("true_vy", &_true_vy, "true_vy/D");
  myTTree->Branch("true_vz", &_true_vz, "true_vz/D");
  myTTree->Branch("interaction_type", &_interaction_type, "interaction_type/i");

  myTTree->Branch("nu_daughters_E", "std::vector< double >", &_nu_daughters_E);
  myTTree->Branch("nu_daughters_pdg", "std::vector< int >", &_nu_daughters_pdg);
  myTTree->Branch("n_true_pions", &_n_true_pions, "n_true_pions/i");
  myTTree->Branch("n_true_protons", &_n_true_protons, "n_true_protons/i");
  myTTree->Branch("n_true_protons_above40", &_n_true_protons_above40, "n_true_protons_above40/i");
  myTTree->Branch("n_true_protons_above21", &_n_true_protons_above21, "n_true_protons_above21/i");
  myTTree->Branch("true_daughter_E", &_true_daughter_E, "true_daughter_E/d");
  myTTree->Branch("true_daughter_theta", &_true_daughter_theta, "true_daughter_theta/d");
  myTTree->Branch("true_daughter_phi", &_true_daughter_phi, "true_daughter_phi/d");
  myTTree->Branch("true_daughter_T", &_true_daughter_T, "true_daughter_T/d");

  _run_subrun_list_file.open("run_subrun_list_filter.txt", std::ofstream::out | std::ofstream::trunc);

  this->reconfigure(p);
}

void generatorLevelAnalyzer::reconfigure(fhicl::ParameterSet const &p)
{
  m_opticalFlashFinderLabel = p.get<std::string>("OpticalFlashFinderLabel", "simpleFlashBeam");
  m_pfp_producer = p.get<std::string>("PFParticleProducer", "pandoraNu");
  m_isOverlaidSample = p.get<bool>("isOverlaidSample", false);
  m_isData = p.get<bool>("isData", false);
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
  _flash_PE.clear();
  _flash_time.clear();

  //pandora neutrino candidates information
  _n_total_candidates = 0;
  _nu_candidate_vx.clear();
  _nu_candidate_vy.clear();
  _nu_candidate_vz.clear();

  //flux information
  _dk_x = std::numeric_limits<double>::lowest();
  _dk_y = std::numeric_limits<double>::lowest();
  _dk_z = std::numeric_limits<double>::lowest();
  _dk_px = std::numeric_limits<double>::lowest();
  _dk_py = std::numeric_limits<double>::lowest();
  _dk_pz = std::numeric_limits<double>::lowest();
  _dk_pdg = std::numeric_limits<int>::lowest();
  _gen_x = std::numeric_limits<double>::lowest();
  _gen_y = std::numeric_limits<double>::lowest();
  _gen_z = std::numeric_limits<double>::lowest();
  _dk2gen = std::numeric_limits<double>::lowest();
  _gen2vtx = std::numeric_limits<double>::lowest();
  _totaltime = std::numeric_limits<double>::lowest();
  _delaytime = std::numeric_limits<double>::lowest();
  
  //true neutrino information
  _n_true_nu = std::numeric_limits<unsigned int>::lowest();
  _nu_pdg = std::numeric_limits<int>::lowest();
  _nu_E = std::numeric_limits<double>::lowest();
  _nu_theta = std::numeric_limits<double>::lowest();
  _nu_phi = std::numeric_limits<double>::lowest();
  _nu_T = std::numeric_limits<double>::lowest();
  _nu_pt = std::numeric_limits<double>::lowest();
  _nu_qsqr = std::numeric_limits<double>::lowest();
  _nu_w = std::numeric_limits<double>::lowest();
  _ccnc = std::numeric_limits<int>::lowest();
  _true_vx = std::numeric_limits<double>::lowest();
  _true_vy = std::numeric_limits<double>::lowest();
  _true_vz = std::numeric_limits<double>::lowest();
  _interaction_type = std::numeric_limits<int>::lowest();

  _nu_daughters_E.clear();
  _nu_daughters_pdg.clear();
  _n_true_pions = std::numeric_limits<unsigned int>::lowest();
  _n_true_protons = std::numeric_limits<unsigned int>::lowest();
  _n_true_protons_above40 = std::numeric_limits<unsigned int>::lowest();
  _n_true_protons_above21 = std::numeric_limits<unsigned int>::lowest();
  _nu_daughters_p.clear();
  _nu_daughters_start_v.clear();
  _nu_daughters_end_v.clear();
  _true_daughter_E = std::numeric_limits<double>::lowest();
  _true_daughter_theta = std::numeric_limits<double>::lowest();
  _true_daughter_phi = std::numeric_limits<double>::lowest();
  _true_daughter_T = std::numeric_limits<double>::lowest();

}

void generatorLevelAnalyzer::storeBNBWeight(art::Event const &evt)
{
  if (!evt.isRealData() || m_isOverlaidSample)
  {
    // nu_e flux must be corrected by event weight
    art::InputTag eventweight_tag("eventweight");
    auto const &eventweights_handle =
        evt.getValidHandle<std::vector<evwgh::MCEventWeight>>(eventweight_tag);
    if (!eventweights_handle.isValid())
    {
      _bnbweight = 1;
    }
    else
    {
      auto const &eventweights(*eventweights_handle);
      if (eventweights.size() > 0)
      {
        for (auto last : eventweights.at(0).fWeight)
        {
          if (last.first.find("bnbcorrection") != std::string::npos && std::isfinite(last.second.at(0)))
          {
            _bnbweight = last.second.at(0);
          }
          else
          {
            _bnbweight = 1;
          }
        }
      }
      else
      {
        _bnbweight = 1;
      }
    }
  }
}

void generatorLevelAnalyzer::opticalInformation(art::Event const &evt)
{
  art::InputTag optical_tag{m_opticalFlashFinderLabel};
  auto const &optical_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(optical_tag);

  for (unsigned int ifl = 0; ifl < optical_handle->size(); ++ifl)
  {
    recob::OpFlash const &flash = optical_handle->at(ifl);
    _flash_PE.push_back(flash.TotalPE());
    _flash_time.push_back(flash.Time());
  }
}

void generatorLevelAnalyzer::PNCandidatesInformation(art::Event const &evt)
{
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);

  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt,
                                                 m_pfp_producer);

  for (size_t _i_pfp = 0; _i_pfp < pfparticle_handle->size(); _i_pfp++)
  {
    if ((abs(pfparticle_handle->at(_i_pfp).PdgCode()) == 12 ||
         abs(pfparticle_handle->at(_i_pfp).PdgCode()) == 14 ||
         abs(pfparticle_handle->at(_i_pfp).PdgCode()) == 16) &&
        pfparticle_handle->at(_i_pfp).IsPrimary())
    {
      _n_total_candidates++;
      auto const &nu_candidate_vertex_obj = vertex_per_pfpart.at(_i_pfp);
      double nu_candidate_vertex[3];
      nu_candidate_vertex_obj->XYZ(nu_candidate_vertex);
      _nu_candidate_vx.push_back(nu_candidate_vertex[0]);
      _nu_candidate_vy.push_back(nu_candidate_vertex[1]);
      _nu_candidate_vz.push_back(nu_candidate_vertex[2]);
    }
  }
}

void generatorLevelAnalyzer::mcFluxInformation(art::Event const &evt)
{
  auto const &fluxes_handle = evt.getValidHandle<std::vector<simb::MCFlux>>(_mctruthLabel);
  auto const &fluxes(*fluxes_handle);
  
  //Assume there's only one MCFlux in the fluxes vector
  if(fluxes.size()>0){
    
    if(fluxes.size()>1) std::cout<<"WARNING: More than one MCFlux present"<<std::endl;
    
    _dk_x = fluxes[0].fvx;
    _dk_y = fluxes[0].fvy;
    _dk_z = fluxes[0].fvz;

    _dk_px = fluxes[0].fpdpx;
    _dk_py = fluxes[0].fpdpy;
    _dk_pz = fluxes[0].fpdpz;
    
    _dk_pdg = fluxes[0].fptype;
    
    _gen_x = fluxes[0].fgenx;
    _gen_y = fluxes[0].fgeny;
    _gen_z = fluxes[0].fgenz;
    
    _dk2gen = fluxes[0].fdk2gen;
    _gen2vtx = fluxes[0].fgen2vtx;
  
    
    
    //Calculate total propagation time
    double const C_SPEED = 299792458;
    double const UBDIST = sqrt(pow(54.499,2)+pow(74.461,2)+pow(677.611,2));
    double parentMass = 0.;
    double parentMomentum = 0;
    double parentSpeed = 0.;

    switch (abs(_dk_pdg)) {
      case 130:
        parentMass = 0.498; //K0L
        break;
      case 321:
        parentMass = 0.494; //K+
        break;
      case 13:
        parentMass = 0.106; //muon
        break;
      case 211:
        parentMass = 0.140; //pi+
        break;
      default:
        std::cout<<"WARNING: Neutrino parent not recognized"<<std::endl;
        break;
    }
    
    parentMomentum = sqrt(pow(_dk_px,2)+pow(_dk_py,2)+pow(_dk_pz,2));
    if(parentMomentum > 0){
      parentSpeed = C_SPEED*sqrt(pow(parentMomentum,2)/(pow(parentMass,2)+pow(parentMomentum,2)));

      _totaltime  = (sqrt(pow(_dk_x,2)+pow(_dk_y,2)+pow(_dk_z,2))/100./parentSpeed + (_dk2gen + _gen2vtx)/C_SPEED)*pow(10,6);

      _delaytime = _totaltime - UBDIST/C_SPEED*pow(10,6);

      std::cout<<sqrt(pow(_dk_x,2)+pow(_dk_y,2)+pow(_dk_z,2))/100./parentSpeed<<" "<<parentSpeed<<" "<<parentSpeed/C_SPEED<<" "<<_delaytime<<" "<<_dk_px<<" "<<_dk_py<<" "<<_dk_pz<<" "<<std::endl;
    }
  }
}

void generatorLevelAnalyzer::trueNeutrinoInformation(art::Event const &evt)
{
  auto const &generator_handle = evt.getValidHandle<std::vector<simb::MCTruth>>(_mctruthLabel);
  auto const &generator(*generator_handle);
  _n_true_nu = generator.size();

  for (auto &gen : generator)
  {
    if (gen.Origin() == simb::kBeamNeutrino)
    {
      _nu_pdg = gen.GetNeutrino().Nu().PdgCode();
      _nu_E = gen.GetNeutrino().Nu().E();
      _nu_theta = gen.GetNeutrino().Nu().Momentum().Theta();
      _nu_phi = gen.GetNeutrino().Nu().Momentum().Phi();
      _nu_T = gen.GetNeutrino().Nu().T();

      _nu_pt = gen.GetNeutrino().Pt();
      _nu_qsqr = gen.GetNeutrino().QSqr();
      _nu_w = gen.GetNeutrino().W();

      _ccnc = gen.GetNeutrino().CCNC();

      _true_vx = gen.GetNeutrino().Nu().Vx();
      _true_vy = gen.GetNeutrino().Nu().Vy();
      _true_vz = gen.GetNeutrino().Nu().Vz();

      _interaction_type = gen.GetNeutrino().Mode();

      break; // In case of events with more than one neutrino (2% of the total) we take for the moment only the first one
      // It should probably be fixed to get the neutrino in the active volume, and see how many events have two neutrinos in the TPC volume
    }
  }

  auto const &mcparticles_handle = evt.getValidHandle<std::vector<simb::MCParticle>>(_mcparticleLabel);
  auto const &mcparticles(*mcparticles_handle);

  for (auto &mcparticle : mcparticles)
  {
    if (!(mcparticle.Process() == "primary" &&
          mcparticle.T() != 0 &&
          mcparticle.StatusCode() == 1))
      continue;

    const auto mc_truth = pandoraHelper.TrackIDToMCTruth(evt, _mcparticleLabel, mcparticle.TrackId());
    if (mc_truth->Origin() == simb::kBeamNeutrino)
    {
      _nu_daughters_E.push_back(mcparticle.E());
      _nu_daughters_pdg.push_back(mcparticle.PdgCode());

      if (abs(mcparticle.PdgCode()) == 211)
      {
        _n_true_pions += 1;
      }
      else if (abs(mcparticle.PdgCode()) == 2212)
      {
        _n_true_protons += 1;
        if (mcparticle.E() > (0.938272 + 0.040))
        {
          _n_true_protons_above40 += 1;
        }
        if (mcparticle.E() > (0.938272 + 0.02108))
        {
          _n_true_protons_above21 += 1;
        }
      }

      std::vector<double> p;
      p.push_back(mcparticle.Px());
      p.push_back(mcparticle.Py());
      p.push_back(mcparticle.Pz());

      _nu_daughters_p.push_back(p);

      std::vector<double> start_v;
      start_v.push_back(mcparticle.Vx());
      start_v.push_back(mcparticle.Vy());
      start_v.push_back(mcparticle.Vz());

      _nu_daughters_start_v.push_back(start_v);

      std::vector<double> end_v;
      end_v.push_back(mcparticle.EndX());
      end_v.push_back(mcparticle.EndY());
      end_v.push_back(mcparticle.EndZ());

      _nu_daughters_end_v.push_back(end_v);

      if (_ccnc==0)
      {
        if (mcparticle.PdgCode() == _pdg_daughter[_nu_pdg])
        {
          _true_daughter_E = mcparticle.E();
          _true_daughter_theta = mcparticle.Momentum().Theta();
          _true_daughter_phi = mcparticle.Momentum().Phi();
          _true_daughter_T = mcparticle.T();
        }
      }

    }
  }
}

void generatorLevelAnalyzer::analyze(art::Event const &evt)
{
  clear();
  std::cout << "RUN " << evt.run() << " SUBRUN " << evt.subRun() << " EVENT " << evt.id().event()
            << std::endl;

  bool is_data = evt.isRealData();

  // BNB fWeight
  storeBNBWeight(evt);

  // optical information
  opticalInformation(evt);

  //PNC information
  PNCandidatesInformation(evt);

  // gen level information
  if (~is_data)
  {
    mcFluxInformation(evt);
    trueNeutrinoInformation(evt);
  }
        
  myTTree->Fill();
}

DEFINE_ART_MODULE(generatorLevelAnalyzer)
