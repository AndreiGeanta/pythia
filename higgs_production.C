#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TPad.h"

void higgs_production(Int_t nev  = 10000 /*<---- no. of events*/, Int_t ndeb = 1) 
/*!!!! --- for a quick test I recommend using nev = 10000 --- !!!!!*/

{
// Load libraries
   gSystem->Load("libEG");
   gSystem->Load("libEGPythia8");
   gSystem->Load("$PYTHIA8/lib/libpythia8");   
   gROOT->SetBatch(kTRUE);   
   
   TFile *file = TFile::Open("invariant_mass_bbar.root"/*<--- filename*/, "RECREATE"); /*produce a root file*/
   if (!file || !file->IsOpen()){
   		cout << "Couldn't open file...";
   }
   
   
   TCanvas *c1 = new TCanvas("c1", "My canvas", 200, 10, 600, 400);
   
// Histograms

   /*TH1F(histo name - this is the name of the histogram that you will find in the root file, 
   		  histo title, 
   		  no. of bins, 
   		  histo lower_limit on x-axis,
   		  histo upper_limit on x-axis)	
   */
   
   TH1F* invariant_H = new TH1F("invariant_H", "", 50, 124.9, 125.1);    
   

// Array of particles
   TClonesArray* particles = new TClonesArray("TParticle", 10000);
// Create pythia8 object
   TPythia8* pythia8 = new TPythia8();

#if PYTHIA_VERSION_INTEGER == 8235
   // Pythia 8.235 is known to cause crashes:
   printf("ABORTING PYTHIA8 TUTORIAL!\n");
   printf("The version of Pythia you use is known to case crashes due to memory errors.\n");
   printf("They have been reported to the authors; the Pythia versions 8.1... are known to work.\n");
   return;
#endif

// Configure
    pythia8->ReadString("HiggsSM:gg2H=on");   //the physics process: gluon fusion
    pythia8->ReadString("PartonLevel:MPI = off");  
	pythia8->ReadString("HadronLevel:Hadronize=off");    //no hadronization
	pythia8->ReadString("PartonLevel:ISR=off");         // no initial-state showers
	pythia8->ReadString("PartonLevel:FSR=off");          // no final-state showers

//	pythia8->ReadString("25:onMode =off");          //cancel all decays channels   
	pythia8->ReadString("25:OnIfMatch = 5 -5");    
    pythia8->ReadString("Random:setSeed = on");
   // use a reproducible seed: always the same results for the tutorial.
//    pythia8->ReadString("Random:seed = 42");


// Initialize

   pythia8->Initialize(2212 /* p */, 2212 /* p */, 13000. /* TeV */); /*we collide protons at 13 TeV CM*/

	
// Event loop
   for (Int_t iev = 0; iev < nev; iev++) {
      pythia8->GenerateEvent();
      if (iev < ndeb) pythia8->EventListing();
      pythia8->ImportParticles(particles,"All");
      Int_t np = particles->GetEntriesFast();
      TLorentzVector b, bbar, higgs; /*define the 4-vectors*/

   	  Double_t e_pT;      
      
// Particle loop
      for (Int_t ip = 0; ip < np; ip++) {// begin for loop
         
      
         TParticle* part = (TParticle*) particles->At(ip);
         
         Int_t ist = part->GetStatusCode();

         // Positive codes are final particles.
         if (ist <= 0) continue;
         Int_t pdg = part->GetPdgCode();
         

         
  		 if (pdg==5 && ist > 0){ //if it is an b in the final state
 			Double_t px1 = part->Px(); //take the x component of momentum
 			Double_t py1 = part->Py(); //take the y component of momentum
 			Double_t pz1 = part->Pz(); //take the z component of momentum
 			Double_t e1 = part->Energy(); //take the energy - time component of the 4-momentum 			 			 			
 			b.SetPxPyPzE(px1, py1, pz1, e1);
 			
 		}
 		
 		 if (pdg==-5 && ist > 0){//if it is a bbar in the final state
 			Double_t px2 = part->Px();
 			Double_t py2 = part->Py();
 			Double_t pz2 = part->Pz();
 			Double_t e2 = part->Energy(); 			 			 			
 			bbar.SetPxPyPzE(px2, py2, pz2, e2);
 		} 
      
} // end for loop
       

      TLorentzVector sum;
      sum = b + bbar;
      Double_t hpx = sum.Px();
      Double_t hpy = sum.Py();
      Double_t hpz = sum.Pz();
      Double_t he = sum.Energy();     
      higgs.SetPxPyPzE(hpx, hpy, hpz, he);     

	  invariant_H->Fill(higgs.Mag()); /*fill the histogram with the magnitude of the 4-vector*/	 
  	  invariant_H->SetXTitle("m_{bbar} [GeV]");
	  invariant_H->SetYTitle("Entries");
	
    
	  
}

	  invariant_H->Draw();
	  c1->Modified(); 
	  c1->Update();	  
	  c1->SaveAs("invariant_mass_bbar.pdf");
	  c1->Update();	  
	  
file->Write();
file->Close();
   
return;
}
