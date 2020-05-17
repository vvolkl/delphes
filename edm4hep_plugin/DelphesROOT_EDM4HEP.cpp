/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <map>
#include <vector>

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

#include "TApplication.h"
#include "TROOT.h"

#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"


using std::cout;
using std::cerr;
using std::endl;

static bool interrupted = false;

void SignalHandler(int sig) {
  interrupted = true;
}


int main(int argc, char *argv[]) {
  std::string appName = "DelphesROOT";
  std::unique_ptr<TFile> outputFile = nullptr;
  std::unique_ptr<TFile> inputFile = nullptr;
  ExRootTreeWriter *treeWriter = 0; // todo: edm4hep
  ExRootTreeBranch *branchEvent = 0; // todo: edm4hep


  Int_t i;
  Long64_t eventCounter, numberOfEvents;

  if(argc < 4)
  {
    cout << " Usage: " << appName << " config_file"
         << " output_file"
         << " input_file(s)" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " input_file(s) - input file(s) in ROOT format." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);


  try {
     outputFile = std::unique_ptr<TFile>(TFile::Open(argv[2], "RECREATE"));

    treeWriter = new ExRootTreeWriter(outputFile.get(), "Delphes");
    branchEvent = treeWriter->NewBranch("Event", HepMCEvent::Class());

    auto confReader = std::make_unique<ExRootConfReader>();
    confReader->ReadFile(argv[1]);

    auto modularDelphes = std::make_unique<Delphes>("Delphes");
    modularDelphes->SetConfReader(confReader.get());
    //modularDelphes->SetTreeWriter(treeWriter);

    auto chain = std::make_unique<TChain>("Delphes");

    auto factory = modularDelphes->GetFactory(); // memory managed in Delphes class

    // has to happen before InitTask
    //TObjArray *allParticleOutputArray = 0, *stableParticleOutputArray = 0, *partonOutputArray = 0;
    TObjArray* allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    TObjArray* stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    TObjArray* partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();


    for(i = 3; i < argc && !interrupted; ++i) {
      cout << "** Reading " << argv[i] << endl;
      chain->Add(argv[i]);
    }

      ExRootTreeReader *treeReader = new ExRootTreeReader(chain.get());

      numberOfEvents = treeReader->GetEntries();

      //if(numberOfEvents <= 0) continue;
      //ExRootProgressBar progressBar(numberOfEvents - 1);
      ExRootProgressBar progressBar(-1);
      // Loop over all objects
      eventCounter = 0;
      modularDelphes->Clear();
      //treeWriter->Clear();

      ////****************************************************** INPUT ^^^^^^^^^^^^^^^^^^^^^
      constexpr double c_light = 2.99792458E8;
      GenParticle *gen;
      Candidate *candidate;
      Int_t pdgCode;
      TClonesArray *branchParticle = treeReader->UseBranch("Particle");
      TClonesArray *branchHepMCEvent = treeReader->UseBranch("Event");
      for (Int_t entry = 0; entry < numberOfEvents && !interrupted; ++entry) {
        treeReader->ReadEntry(entry);
        for(Int_t j = 0; j < branchParticle->GetEntriesFast(); j++) {
          gen = (GenParticle *)branchParticle->At(j);
          candidate = modularDelphes->GetFactory()->NewCandidate();
          candidate->Momentum = gen->P4();
          candidate->Position.SetXYZT(gen->X, gen->Y, gen->Z, gen->T * 1.0E3 * c_light);
          candidate->PID = gen->PID;
          candidate->Status = gen->Status;
          candidate->M1 = gen->M1;
          candidate->M2 = gen->M2;
          candidate->D1 = gen->D1;
          candidate->D2 = gen->D2;
          candidate->Charge = gen->Charge;
          candidate->Mass = gen->Mass;
          allParticleOutputArray->Add(candidate);
          pdgCode = TMath::Abs(gen->PID);
          if(gen->Status == 1) {
            stableParticleOutputArray->Add(candidate);
          } else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15) {
            partonOutputArray->Add(candidate);
          }
        }
        ////****************************************************** INPUT ^^^^^^^^^^^^^^^^^^^^^

        modularDelphes->ProcessTask();

        //******************************************************************************** todo: edm4hep
        const TObjArray* delphesColl = modularDelphes->ImportArray("FatJetFinder/jets");
        for (int j = 0; j < delphesColl->GetEntries(); j++) {
          auto cand = static_cast<Candidate*>(delphesColl->At(j));
          std::cout << cand->Momentum.Pt() << " ";
        }
        std::cout << std::endl;
        //******************************************************************************** todo: edm4hep

        //treeWriter->Fill();
        modularDelphes->Clear();
        //treeWriter->Clear();
        progressBar.Update(eventCounter, eventCounter);
        ++eventCounter;
      }
      progressBar.Update(eventCounter, eventCounter, kTRUE);
      progressBar.Finish();
      delete treeReader; //todo: edm4hep

      modularDelphes->FinishTask();
      //treeWriter->Write(); 
      //todo: edm4hep
      cout << "** Exiting..." << endl;
      return 0;
    } catch(std::runtime_error &e) {
      cerr << "** ERROR: " << e.what() << endl;
      return 1;
    }
}
