#define EtaMass_cxx
#include "EtaMass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void EtaMass::Loop(){
  
  TCanvas *canvas = new TCanvas("canvas", "", 1000, 750);
  TH1D *Eta_Mass = new TH1D("Eta_Mass", "", 200, 0, 1);
  TH1D *MinPhi = new TH1D("MinPhi", "", 120, 0, 1.2);
  TH1D *MinNu = new TH1D("MinNu", "", 500, 0, 1);
  TH2D *MinDelta = new TH2D("MinDelta", "#Delta #phi vs #Delta #eta", 500, 0, 1, 50, 0, 0.01);

  Eta_Mass->Sumw2();
  Eta_Mass->SetMarkerStyle(20);
  Eta_Mass->SetMarkerColor(EColor::kBlack);
  Eta_Mass->SetLineColor(EColor::kBlack);
  Eta_Mass->SetXTitle("Mass (GeV)");

  TFile* fout = new TFile("fout810TA.root","RECREATE");
  float pi_m;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  for(Long64_t i = 0; i < nentries ; i++){
    cout << setw(7) <<i << '\r' << flush;
    fChain->GetEntry(i);     
    if (ncluster < 2) continue;
    map <Float_t, int> clusterpts;
 
   for (ULong64_t n = 0; n < ncluster; n++) {
     if (cluster_pt[n]<8) continue;
      if (cluster_tof[n] > 5e-8) continue;
      Float_t delta_phi,delta_eta = 1;
      Float_t temp_phi,temp_eta;
      bool trackveto = false;
      int t = 0;
      while ( !trackveto && (t<ntrack)){ //FIXME: Try pT dependant track veto
	 if ((abs((cluster_eta[n] - track_eta_emcal[t])) <= 0.03) && 
	     (abs((cluster_phi[n] - track_phi_emcal[t])) <= 0.035)&&
	     (track_pt[t]>0.05) && (abs(track_eta[t])<0.9)) trackveto = true;
	 
	 temp_phi = abs(cluster_phi[n] - track_phi_emcal[t]);
	 temp_eta = abs(cluster_eta[n] - track_eta_emcal[t]);
	 if (temp_phi < delta_phi) delta_phi = temp_phi;
	 if (temp_eta < delta_eta) delta_eta = temp_eta;
	 t++;  }

       MinDelta->Fill(delta_eta,delta_phi);
       MinNu->Fill(delta_eta);
       MinPhi->Fill(delta_phi);

       if (trackveto) continue;
       clusterpts.insert(make_pair(cluster_pt[n],n));
    } //cluster loop

    if (clusterpts.size() < 2) continue;    
    map<Float_t,int>::reverse_iterator rit;
    int c = 0;
    int n1,n2;

    for (rit=clusterpts.rbegin(); rit!=clusterpts.rend(); rit++, c++){
      if (c==0) n1 = rit->second;      
      if (c==1) n2 = rit->second;     
      if(c==2) break;  }
    if (cluster_pt[n1] < 10 || cluster_pt[n2] < 8) continue;

    TLorentzVector v1;
    v1.SetPtEtaPhiM(cluster_pt[n1], cluster_eta[n1], cluster_phi[n1], 0.0);    
    TLorentzVector v2;    v2.SetPtEtaPhiM(cluster_pt[n2], cluster_eta[n2], cluster_phi[n2], 0.0);
    TLorentzVector pi0 = v1+v2;
    if (pi0.Pt() < 6) continue;    
    if(v1.Vect().Angle(v2.Vect()) < 0.022) continue; //22mrad
    Eta_Mass->Fill(pi0.M()); 
   
    if (i % 10000 == 0) {
      Eta_Mass->Draw();
      canvas->Update();}    

    clusterpts.clear(); 

  }//end event loop
  
  Eta_Mass->Draw();
  Eta_Mass->Write();
  
  TCanvas *Del = new TCanvas("Del","",1000,750);
  MinDelta->GetXaxis()->SetTitle("#Delta #eta");
  MinDelta->GetYaxis()->SetTitle("#Delta #phi");
  MinDelta->Draw();
  MinDelta->Write();
  MinNu->Write();
  MinPhi->Write();
  fout->Close();
}
