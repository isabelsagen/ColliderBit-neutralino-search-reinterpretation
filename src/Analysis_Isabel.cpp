
#include "gambit/ColliderBit/analyses/Analysis.hpp" 
#include "gambit/ColliderBit/ATLASEfficiencies.hpp" 
#include <cfloat>


// tagging b-jets at b_Efficiency effeciency

inline std::vector<const HEPUtils::Jet*> applyBTaggingEfficiency(std::vector<const HEPUtils::Jet*>& jets, double b_Efficiency){
  std::vector<const HEPUtils::Jet*> updatedJets;
    for (const HEPUtils::Jet* jet : jets){ 
      HEPUtils::Jet* newjet = new HEPUtils::Jet(*jet); 
      if (jet->btag()){ 
        double rnd = HEPUtils::rand01();
          if (rnd >= b_Efficiency){ 
            newjet->set_btag(false);
          }
        } 
      updatedJets.push_back(newjet);
      } 
    return updatedJets; 
  }


//Figure 18a in https://arxiv.org/abs/1908.00005

inline void applyTightPhoton(std::vector<const HEPUtils::Particle*>& photons){
    
  const static std::vector<double> binedges_Et = {0., 15., 20., 25., 30., 35., 40., 45., 50., 60., 80., 100., 130., 160., 180., 250., DBL_MAX};
  const static std::vector<double> bineffs_tight = {0.69,  0.78,  0.83,  0.86,  0.87,  0.89,  0.89,  0.9,  0.91,  0.92,  0.92,  0.93,  0.93,  0.93,  0.93,  0.92};

  HEPUtils::BinnedFn1D<double> photon_eff(binedges_Et, bineffs_tight);
  Gambit::ColliderBit::filtereff_pt(photons, photon_eff);
  }


//Figure 17 https://arxiv.org/pdf/1908.00005

inline void applyLooseElectron(std::vector<const HEPUtils::Particle*>& electrons){

  const static std::vector<double> binedges_Et = { 0,  6,  10,  15,  20,  25,  30,  35,  40,  45,  50,  60,  80,  DBL_MAX };
  const static std::vector<double> bineffs_Et  = { 0.98,  0.93,  0.88,  0.83,  0.86,  0.88,  0.90,  0.92,  0.93,  0.94,  0.94,  0.945,  0.95 };
  const static HEPUtils::BinnedFn1D<double> _eff_et(binedges_Et, bineffs_Et);

  const static std::vector<double> binedges_eta = {-DBL_MAX,  -2.5, -2., -1.8,  -1.5,  -1.35,  -1.15,  -0.8,  -0.6,  -0.1,  0.,  0.1,  0.6,  0.8,  1.15,  1.4,  1.5,  1.8,  2.,  2.5, DBL_MAX };
  const static std::vector<double> bineffs_eta  = {0.82,  0.87,  0.9,  0.92,  0.86,  0.92,  0.93,  0.94,  0.94,  0.91,  0.9,  0.94,  0.94,  0.94,  0.92,  0.87,  0.92,  0.9,  0.87,  0.83 };
  const static HEPUtils::BinnedFn1D<double> _eff_eta(binedges_eta, bineffs_eta);

  auto keptElectronsEnd = std::remove_if(electrons.begin(), electrons.end(),[](const HEPUtils::Particle* electron) {
    const double e_pt = electron->pT();
    const double e_aeta = electron->eta();


    const double eff1 = _eff_eta.get_at(e_aeta), eff2 = _eff_et.get_at(e_pt);
    const double eff = std::min(eff1 * eff2 / 0.95, 1.0); //< norm factor as approximate double differential
    return Gambit::ColliderBit::random_bool(1-eff);
    } );
  electrons.erase(keptElectronsEnd, electrons.end());
}



namespace Gambit {
  namespace ColliderBit {

    //#define CUTFLOW

    using namespace std; 


    class Analysis_ATLAS_Isabel : public Analysis {
    private:

      double _numSR;

    public:

      static constexpr const char* detector = "ATLAS";


      std::map<string, EventCounter> _counters = {
        {"SR1h", EventCounter("SR1h")},
        {"SR1Z", EventCounter("SR1Z")},
        {"SR2", EventCounter("SR2")},

        {"common_selection", EventCounter("common_selection")},

        {"SR1h_ETmiss_100", EventCounter("SR1h_ETmiss_100")},
        {"SR1h_m_diphoton", EventCounter("SR1h_m_diphoton")},
        {"SR1h_m_bb", EventCounter("SR1h_m_bb")},
        {"SR1h_pm_04", EventCounter("SR1h_pm_04")},
        {"SR1h_p_diphoton_90", EventCounter("SR1h_p_diphoton_90")},

        {"SR1Z_ETmiss_100", EventCounter("SR1Z_ETmiss_100")},
        {"SR1Z_m_diphoton", EventCounter("SR1Z_m_diphoton")},
        {"SR1Z_m_bb", EventCounter("SR1Z_m_bb")},
        {"SR1Z_pm_04", EventCounter("SR1Z_pm_04")},
        {"SR1Z_p_diphoton_90", EventCounter("SR1Z_p_diphoton_90")},

        {"SR2_ETmiss_100", EventCounter("SR2_ETmiss_100")},
        {"SR2_m_diphoton", EventCounter("SR2_m_diphoton")},
        {"SR2_m_bb", EventCounter("SR2_m_bb")},




      };
      

      Analysis_ATLAS_Isabel() {

        // Set the analysis name
        set_analysis_name("ATLAS_Isabel");

        // Set the LHC luminosity
        set_luminosity(139); 

        // Set number of events passing cuts to zero upon initialisation
        _numSR = 0;




      }


      void run(const HEPUtils::Event* event){ 
        
        double met = event->met(); 


        //-------------------------------------------------------------------------------------
        //                       Initial reqirements
        //-------------------------------------------------------------------------------------



        // Electrons

        vector<const HEPUtils::Particle*> baselineElectrons;
        for (const HEPUtils::Particle* electron : event->electrons()) {
          if (electron->pT() > 4.5 && fabs(electron->eta()) < 2.47){
            baselineElectrons.push_back(electron);
          } 
        }


        // Muons

        vector<const HEPUtils::Particle*> baselineMuons;
        for (const HEPUtils::Particle* muon : event->muons()) {
          if (muon->pT() > 4. && fabs(muon->eta()) < 2.47){
            baselineMuons.push_back(muon);
          } 
        }


        // Jets

        vector<const HEPUtils::Jet*> baselineJets;
          for (const HEPUtils::Jet* jet : event->jets("antikt_R04")) {
            if (jet->pT() > 30. && fabs(jet->eta()) < 2.8){ 
              baselineJets.push_back(jet);
            } 
          }
        

        // Photons

        vector<const HEPUtils::Particle*> baselinePhotons;
        for (const HEPUtils::Particle* photon : event->photons()){
          if (fabs(photon->eta()) < 2.37 && photon->pT() > 25.  && !( (fabs(photon->eta()) > 1.37) && (fabs(photon->eta()) < 1.52)) ){
            baselinePhotons.push_back(photon);
          } 
        }



        applyLooseElectron(baselineElectrons); 
        ATLAS::applyMuonIDEfficiency2020(baselineMuons, "Medium"); 

        removeOverlap(baselineElectrons, baselinePhotons, 0.4);
        removeOverlap(baselineMuons, baselinePhotons, 0.4);
        removeOverlap(baselineJets, baselinePhotons, 0.4);

        vector<const HEPUtils::Jet*> baselineJets_ = applyBTaggingEfficiency(baselineJets, 0.77); 

        applyTightPhoton(baselinePhotons); 



        vector<const HEPUtils::Jet*> signalBJets;
        for (const HEPUtils::Jet* jet : baselineJets_){
          if (jet->btag()){
            signalBJets.push_back(jet);
          }
        }


        
        vector<const HEPUtils::Particle*> signalPhotons = baselinePhotons; 


        //-------------------------------------------------------------------------------------
        //                                PRESELECTION
        //-------------------------------------------------------------------------------------

        int nElectrons = baselineElectrons.size();
        int nMuons = baselineMuons.size();
        int nJets = baselineJets_.size();
        int nBJets  = signalBJets.size();
        int nPhotons = signalPhotons.size(); 

      
      

        bool canidate_event = false;
        
        if(nPhotons >= 2){ 

          if (nElectrons == 0 && nMuons == 0){
        
            if (nPhotons == 2){ 
              sortByPt(signalPhotons); 

              if( (signalPhotons.at(0)->pT() > 35.)){ 
                      
                double Diphoton_invm = (signalPhotons.at(0)->mom() + signalPhotons.at(1)->mom()).m();
                if ( (Diphoton_invm > 95.) && (Diphoton_invm < 160.) ){ 

                  if ( ( (signalPhotons.at(1)->pT())/Diphoton_invm) > 0.2){

                    if (nBJets  == 2){ 

                      if (nJets == 2){ 

                        _counters.at("common_selection").add_event(event); 

                        canidate_event = true;

                      }
                    }
                  }
                }              
              }
            }
          }
        }
        
        //-------------------------------------------------------------------------------------
        //                                SIGNAL EVENTS
        //-------------------------------------------------------------------------------------
            
        if (canidate_event){

          double Diphoton_invm = (signalPhotons.at(0)->mom() + signalPhotons.at(1)->mom()).m();
          double bb_invm = (signalBJets.at(0)->mom() + signalBJets.at(1)->mom()).m();
          double diphoton_pT = (signalPhotons.at(0)->mom() + signalPhotons.at(1)->mom()).pT();
          double ratio_m_pT = signalPhotons.at(1)->mom().pT()/Diphoton_invm; 


          if (met <= 100.){

            _counters.at("SR1h_ETmiss_100").add_event(event);
            _counters.at("SR1Z_ETmiss_100").add_event(event);

            if (Diphoton_invm > 120. && Diphoton_invm < 130.){

              _counters.at("SR1h_m_diphoton").add_event(event);
              _counters.at("SR1Z_m_diphoton").add_event(event);

              //SR1h
              if (bb_invm > 100. && bb_invm < 140.){

                _counters.at("SR1h_m_bb").add_event(event);

                if (ratio_m_pT > 0.4){

                  _counters.at("SR1h_pm_04").add_event(event);

                  if (diphoton_pT > 90.){

                    _counters.at("SR1h_p_diphoton_90").add_event(event);
                    _counters.at("SR1h").add_event(event);

                  }
                }
              }

              //SR1Z
              if (bb_invm > 60. && bb_invm < 100.){

                _counters.at("SR1Z_m_bb").add_event(event);

                if (ratio_m_pT > 0.4){

                  _counters.at("SR1Z_pm_04").add_event(event);

                  if (diphoton_pT > 90.){

                    _counters.at("SR1Z_p_diphoton_90").add_event(event);
                    _counters.at("SR1Z").add_event(event);

                  }
                }
              }
            }
          }

          // SR2
          if (met > 100.){ 

            _counters.at("SR2_ETmiss_100").add_event(event);

            if (Diphoton_invm > 120. && Diphoton_invm < 130.){

              _counters.at("SR2_m_diphoton").add_event(event);

              if (ratio_m_pT > 0.2){ 

                if (bb_invm > 35. && bb_invm <= 145.){

                  _counters.at("SR2_m_bb").add_event(event);
                  _counters.at("SR2").add_event(event);

                }
              } 
            }
          }
        }
        
        // Delete new'd pointers
        for (const HEPUtils::Jet* jet : baselineJets_){ 
          delete jet;
        }
      

      }


      /// Combine the variables of another copy of this analysis (typically on another thread) into this one.
      void combine(const Analysis* other)
      {
        const Analysis_ATLAS_Isabel* specificOther = dynamic_cast<const Analysis_ATLAS_Isabel*>(other);
        _numSR += specificOther->_numSR;
      }


      void collect_results() {


        std::cout  << "SR1h=" << _counters.at("SR1h").sum() << ", " << "SR1Z=" << _counters.at("SR1Z").sum() << ", "  << "SR2=" << _counters.at("SR2").sum() << std::endl;
        add_result(SignalRegionData("SR1h", 3., {_counters.at("SR1h").sum(), 0.}, {3.9, 0.6}));
        add_result(SignalRegionData("SR1Z", 5., {_counters.at("SR1Z").sum(), 0.}, {6.4, 1.0}));
        add_result(SignalRegionData("SR2",  2., {_counters.at("SR2").sum(),  0.},  {1.7, 0.7}));



        #ifdef CUTFLOW 
        
        std::cout << "Candidate events " << _counters.at("common_selection").sum() << std::endl;
        std::cout << "=====================================" << std::endl;
         
        std::cout << "=== SR1h ===" << std::endl;
        std::cout << "ETmiss>100: " << _counters.at("SR1h_ETmiss_100").sum() << std::endl;
        std::cout << "Diphoton mass: " << _counters.at("SR1h_m_diphoton").sum() << std::endl;
        std::cout << "bb invariant mass: " << _counters.at("SR1h_m_bb").sum() << std::endl;
        std::cout << "Photon pT ratio>0.4: " << _counters.at("SR1h_pm_04").sum() << std::endl;
        std::cout << "Diphoton pT>90: " << _counters.at("SR1h_p_diphoton_90").sum() << std::endl;
        std::cout << "=====================================" << std::endl;

        std::cout << "=== SR1Z ===" << std::endl;
        std::cout << "ETmiss>100: " << _counters.at("SR1Z_ETmiss_100").sum() << std::endl;
        std::cout << "Diphoton mass: " << _counters.at("SR1Z_m_diphoton").sum() << std::endl;
        std::cout << "bb invariant mass: " << _counters.at("SR1Z_m_bb").sum() << std::endl;
        std::cout << "Photon pT ratio>0.4: " << _counters.at("SR1Z_pm_04").sum() << std::endl;
        std::cout << "Diphoton pT>90: " << _counters.at("SR1Z_p_diphoton_90").sum() << std::endl;
        std::cout << "=====================================" << std::endl;

        std::cout << "=== SR2 ===" << std::endl;
        std::cout << "ETmiss>100: " << _counters.at("SR2_ETmiss_100").sum() << std::endl;
        std::cout << "Diphoton mass: " << _counters.at("SR2_m_diphoton").sum() << std::endl;
        std::cout << "bb invariant mass: " << _counters.at("SR2_m_bb").sum() << std::endl;
        std::cout << "=====================================" << std::endl;
        
        std::cout << "=== Final ===" << std::endl;
        std::cout << "SR1h: " << _counters.at("SR1h").sum() << std::endl;
        std::cout << "SR1Z: " << _counters.at("SR1Z").sum() << std::endl;
        std::cout << "SR2 : " << _counters.at("SR2").sum() << std::endl;
        std::cout << "=====================================" << std::endl;

        #endif

        

      }



    protected:
      void analysis_specific_reset() {
        _numSR = 0;
        for (auto& pair : _counters) { pair.second.reset(); } 

      }

      ///////////////////

    };

    DEFINE_ANALYSIS_FACTORY(ATLAS_Isabel)

  }
}

