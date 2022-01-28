#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TPad.h"

void pythia8_Z2ee(Int_t nev  = 100000 /*<---- no. of events*/, Int_t ndeb = 1) 
/*!!!! --- for a quick test I recommend using nev = 10000 --- !!!!!*/

{
// Load libraries
   gSystem->Load("libEG");
   gSystem->Load("libEGPythia8");
   gSystem->Load("$PYTHIA8/lib/libpythia8");   
   gROOT->SetBatch(kTRUE);   
   
   TFile *file = TFile::Open("invariant_mass_ee.root"/*<--- filename*/, "RECREATE"); /*produce a root file*/
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
   
   TH1F* Zhisto = new TH1F("invariant_mass_Zee", "", 100, 80, 100);    
   TH1F* electron_pT = new TH1F("electron_pT", "", 100, 0, 60);    
   

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
    pythia8->ReadString("WeakSingleBoson:ffbar2gmZ = on");   //the physics process: ffbar -> gamma/Z 
    pythia8->ReadString("PartonLevel:MPI = off");  
	pythia8->ReadString("HadronLevel:Hadronize=off");    //no hadronization
	pythia8->ReadString("PartonLevel:ISR=off");         // no initial-state showers
	pythia8->ReadString("PartonLevel:FSR=off");          // no final-state showers

//	pythia8->ReadString("23:onMode =on");          //cancel all decays channels   
//	pythia8->ReadString("23:OnIfMatch = 11 -11");    
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
      TLorentzVector electron, positron, Z; /*define the 4-vectors*/

   	  Double_t e_pT;      
      
// Particle loop
      for (Int_t ip = 0; ip < np; ip++) {// begin for loop
         
      
         TParticle* part = (TParticle*) particles->At(ip);
         
         Int_t ist = part->GetStatusCode();

         // Positive codes are final particles.
         if (ist <= 0) continue;
         Int_t pdg = part->GetPdgCode();
         

         
  		 if (pdg==11 && ist > 0){ //if it is an electron in the final state
 			Double_t px1 = part->Px(); //take the x component of momentum
 			Double_t py1 = part->Py(); //take the y component of momentum
 			Double_t pz1 = part->Pz(); //take the z component of momentum
 			Double_t e1 = part->Energy(); //take the energy - time component of the 4-momentum 			 			 			
 			electron.SetPxPyPzE(px1, py1, pz1, e1);
 			
 			e_pT = part->Pt();
 		}
 		
 		 if (pdg==-11 && ist > 0){//if it is a positron in the final state
 			Double_t px2 = part->Px();
 			Double_t py2 = part->Py();
 			Double_t pz2 = part->Pz();
 			Double_t e2 = part->Energy(); 			 			 			
 			positron.SetPxPyPzE(px2, py2, pz2, e2);
 		} 
      
} // end for loop
       

      TLorentzVector sum;
      sum = electron + positron;
      Double_t hpx = sum.Px();
      Double_t hpy = sum.Py();
      Double_t hpz = sum.Pz();
      Double_t he = sum.Energy();     
      Z.SetPxPyPzE(hpx, hpy, hpz, he);     

	  Zhisto->Fill(Z.Mag()); /*fill the histogram with the magnitude of the 4-vector*/	 
  	  Zhisto->SetXTitle("m_{ee} [GeV]");
	  Zhisto->SetYTitle("Entries");
	  
	  electron_pT->Fill(e_pT); /*fill the histogram with electron transverse momentum*/	 
  	  electron_pT->SetXTitle("p_{T} [GeV]");
	  electron_pT->SetYTitle("Entries");	  
	  
}

	  Zhisto->Draw();
	  c1->Modified(); 
	  c1->Update();	  
	  c1->SaveAs("invariant_mass_ee.pdf");
	  c1->Update();	  
	  electron_pT->Draw();
	  c1->SaveAs("e_pT.pdf");
	  
file->Write();
file->Close();
   
return;
}
