//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 14 14:18:12 2017 by ROOT version 6.08/00
// from TTree _tree_event/
// found on file: 13d_v4_small.root
//////////////////////////////////////////////////////////

#ifndef EtaMass_h
#define EtaMass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class EtaMass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Char_t          id_git;
   Char_t          version_aliroot;
   Char_t          version_aliphysics;
   Char_t          version_jec;
   Char_t          grid_data_dir;
   Char_t          grid_data_pattern;
   Int_t           beam_particle[2];
   UChar_t         ntrigger_class;
   Char_t          trigger_class[1];   //[ntrigger_class]
   Int_t           run_number;
   UInt_t          period_number;
   UInt_t          orbit_number;
   UShort_t        bunch_crossing_number;
   UInt_t          time_stamp;
   ULong64_t       trigger_mask[2];
   Float_t         multiplicity_v0[64];
   Float_t         centrality_v0m;
   Float_t         centrality[9];
   Float_t         event_plane_psi_v0[3];
   Double_t        event_plane_q_v0[3][2];
   Bool_t          has_misalignment_matrix;
   Float_t         cell_eta[17664];
   Float_t         cell_phi[17664];
   Float_t         cell_voronoi_area[17664];
   Double_t        primary_vertex[3];
   Double_t        primary_vertex_sigma[3];
   Int_t           primary_vertex_ncontributor;
   Double_t        primary_vertex_spd[3];
   Double_t        primary_vertex_spd_sigma[3];
   Int_t           primary_vertex_spd_ncontributor;
   Int_t           npileup_vertex_spd;
   Int_t           pileup_vertex_spd_ncontributor;
   Double_t        pileup_vertex_spd_min_z_distance;
   Int_t           eg_signal_process_id;
   Int_t           eg_mpi;
   Float_t         eg_pt_hat;
   Float_t         eg_cross_section;
   Float_t         eg_weight;
   Float_t         eg_primary_vertex[3];
   Int_t           eg_ntrial;
   Float_t         eg_scale_pdf;
   Float_t         eg_alpha_qcd;
   Float_t         eg_alpha_qed;
   Int_t           eg_pdf_id[2];
   Float_t         eg_pdf_x[2];
   Float_t         eg_pdf_x_pdf[2];
   Bool_t          debug_libmklml_gnu_loaded;
   Char_t          debug_libmklml_gnu_error[8197];
   UInt_t          ncluster;
   Float_t         cluster_e[128];   //[ncluster]
   Float_t         cluster_pt[128];   //[ncluster]
   Float_t         cluster_eta[128];   //[ncluster]
   Float_t         cluster_phi[128];   //[ncluster]
   Float_t         cluster_lambda_square[128][2];   //[ncluster]
   Float_t         cluster_tof[128];   //[ncluster]
   Int_t           cluster_ncell[128];   //[ncluster]
   UShort_t        cluster_cell_id_max[128];   //[ncluster]
   Float_t         cluster_e_max[128];   //[ncluster]
   Float_t         cluster_e_cross[128];   //[ncluster]
   UChar_t         cluster_nlocal_maxima[128];   //[ncluster]
   UInt_t          cluster_nmc_truth[128];   //[ncluster]
   UShort_t        cluster_mc_truth_index[128][32];   //[ncluster]
   Float_t         cluster_iso_tpc_01[128];   //[ncluster]
   Float_t         cluster_iso_tpc_01_ue[128];   //[ncluster]
   Float_t         cluster_iso_tpc_02[128];   //[ncluster]
   Float_t         cluster_iso_tpc_02_ue[128];   //[ncluster]
   Float_t         cluster_iso_tpc_03[128];   //[ncluster]
   Float_t         cluster_iso_tpc_03_ue[128];   //[ncluster]
   Float_t         cluster_iso_tpc_04[128];   //[ncluster]
   Float_t         cluster_iso_tpc_04_ue[128];   //[ncluster]
   Float_t         cluster_iso_its_01[128];   //[ncluster]
   Float_t         cluster_iso_its_01_ue[128];   //[ncluster]
   Float_t         cluster_iso_its_02[128];   //[ncluster]
   Float_t         cluster_iso_its_02_ue[128];   //[ncluster]
   Float_t         cluster_iso_its_03[128];   //[ncluster]
   Float_t         cluster_iso_its_03_ue[128];   //[ncluster]
   Float_t         cluster_iso_its_04[128];   //[ncluster]
   Float_t         cluster_iso_its_04_ue[128];   //[ncluster]
   Float_t         cluster_frixione_tpc_04_02[128];   //[ncluster]
   Float_t         cluster_frixione_tpc_04_05[128];   //[ncluster]
   Float_t         cluster_frixione_tpc_04_10[128];   //[ncluster]
   Float_t         cluster_frixione_its_04_02[128];   //[ncluster]
   Float_t         cluster_frixione_its_04_05[128];   //[ncluster]
   Float_t         cluster_frixione_its_04_10[128];   //[ncluster]
   Float_t         cluster_iso_01_truth[128];   //[ncluster]
   Float_t         cluster_iso_02_truth[128];   //[ncluster]
   Float_t         cluster_iso_03_truth[128];   //[ncluster]
   Float_t         cluster_iso_04_truth[128];   //[ncluster]
   Float_t         cluster_frixione_04_02_truth[128];   //[ncluster]
   Float_t         cluster_frixione_04_05_truth[128];   //[ncluster]
   Float_t         cluster_frixione_04_10_truth[128];   //[ncluster]
   Float_t         cluster_s_nphoton[128][4];   //[ncluster]
   Float_t         cluster_s_ncharged_hadron[128][4];   //[ncluster]
   Float_t         cell_e[17664];
   Float_t         cell_tof[17664];
   UShort_t        cell_cluster_index[17664];
   UShort_t        cell_mc_truth_index[17664];
   UInt_t          ntrack;
   Float_t         track_e[505];   //[ntrack]
   Float_t         track_pt[505];   //[ntrack]
   Float_t         track_eta[505];   //[ntrack]
   Float_t         track_phi[505];   //[ntrack]
   Float_t         track_eta_emcal[505];   //[ntrack]
   Float_t         track_phi_emcal[505];   //[ntrack]
   Char_t          track_charge[505];   //[ntrack]
   UChar_t         track_quality[505];   //[ntrack]
   Float_t         track_tpc_dedx[505];   //[ntrack]
   Float_t         track_tpc_length_active_zone[505];   //[ntrack]
   UChar_t         track_tpc_xrow[505];   //[ntrack]
   UChar_t         track_tpc_ncluster[505];   //[ntrack]
   UChar_t         track_tpc_ncluster_dedx[505];   //[ntrack]
   UChar_t         track_tpc_ncluster_findable[505];   //[ntrack]
   UChar_t         track_its_ncluster[505];   //[ntrack]
   Float_t         track_its_chi_square[505];   //[ntrack]
   Float_t         track_dca_xy[505];   //[ntrack]
   Float_t         track_dca_z[505];   //[ntrack]
   UShort_t        track_mc_truth_index[505];   //[ntrack]
   Float_t         track_voronoi_area[505];   //[ntrack]
   UInt_t          nmuon_track;
   Float_t         muon_track_e[8];   //[nmuon_track]
   Float_t         muon_track_pt[8];   //[nmuon_track]
   Float_t         muon_track_eta[8];   //[nmuon_track]
   Float_t         muon_track_phi[8];   //[nmuon_track]
   Float_t         muon_track_r_abs[8];   //[nmuon_track]
   Float_t         muon_track_p_dca[8];   //[nmuon_track]
   Float_t         muon_track_sigma_p_dca[8];   //[nmuon_track]
   Float_t         muon_track_delta_sagitta_p[8];   //[nmuon_track]
   Float_t         muon_track_distance_sigma_slope_p[8];   //[nmuon_track]
   UShort_t        muon_track_mc_truth_index[8];   //[nmuon_track]
   UInt_t          nmc_truth;
   Float_t         mc_truth_e[1];   //[nmc_truth]
   Float_t         mc_truth_pt[1];   //[nmc_truth]
   Float_t         mc_truth_eta[1];   //[nmc_truth]
   Float_t         mc_truth_phi[1];   //[nmc_truth]
   Char_t          mc_truth_charge[1];   //[nmc_truth]
   Short_t         mc_truth_pdg_code[1];   //[nmc_truth]
   UChar_t         mc_truth_status[1];   //[nmc_truth]
   UChar_t         mc_truth_generator_index[1];   //[nmc_truth]
   Short_t         mc_truth_first_parent_pdg_code[1];   //[nmc_truth]
   Float_t         mc_truth_first_parent_e[1];   //[nmc_truth]
   Float_t         mc_truth_first_parent_pt[1];   //[nmc_truth]
   Float_t         mc_truth_first_parent_eta[1];   //[nmc_truth]
   Float_t         mc_truth_first_parent_phi[1];   //[nmc_truth]
   UShort_t        mc_truth_sibling_index[1];   //[nmc_truth]
   UInt_t          debug_njet_ue_estimation;
   Float_t         debug_jet_ue_estimation_pt_raw[56];   //[debug_njet_ue_estimation]
   Float_t         debug_jet_ue_estimation_eta_raw[56];   //[debug_njet_ue_estimation]
   Float_t         debug_jet_ue_estimation_phi_raw[56];   //[debug_njet_ue_estimation]
   Float_t         debug_jet_ue_estimation_area_raw[56];   //[debug_njet_ue_estimation]
   UInt_t          njet_ak04tpc;
   Float_t         debug_jet_ak04tpc_tag_dr_square[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_e_raw[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_e[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_e_charged[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_pt_raw_ue[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_pt_raw[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_pt[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_pt_charged[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_eta_raw[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_eta[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_phi[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_area_raw[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_area[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_emf_raw[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_emf[41];   //[njet_ak04tpc]
   UShort_t        jet_ak04tpc_multiplicity_raw[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_multiplicity[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_width_sigma_raw[41][2];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_width_sigma[41][2];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_ptd_raw[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_ptd[41];   //[njet_ak04tpc]
   Int_t           jet_ak04tpc_truth_index_z_truth[41][2];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_truth_z_truth[41][2];   //[njet_ak04tpc]
   Int_t           jet_ak04tpc_truth_index_z_reco[41][2];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_truth_z_reco[41][2];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_e_truth[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_pt_truth[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_eta_truth[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_phi_truth[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_area_truth[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_emf_truth[41];   //[njet_ak04tpc]
   UShort_t        jet_ak04tpc_multiplicity_truth[41];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_width_sigma_truth[41][2];   //[njet_ak04tpc]
   Float_t         jet_ak04tpc_ptd_truth[41];   //[njet_ak04tpc]
   UInt_t          njet_ak04its;
   Float_t         debug_jet_ak04its_tag_dr_square[47];   //[njet_ak04its]
   Float_t         jet_ak04its_e_raw[47];   //[njet_ak04its]
   Float_t         jet_ak04its_e[47];   //[njet_ak04its]
   Float_t         jet_ak04its_e_charged[47];   //[njet_ak04its]
   Float_t         jet_ak04its_pt_raw_ue[47];   //[njet_ak04its]
   Float_t         jet_ak04its_pt_raw[47];   //[njet_ak04its]
   Float_t         jet_ak04its_pt[47];   //[njet_ak04its]
   Float_t         jet_ak04its_pt_charged[47];   //[njet_ak04its]
   Float_t         jet_ak04its_eta_raw[47];   //[njet_ak04its]
   Float_t         jet_ak04its_eta[47];   //[njet_ak04its]
   Float_t         jet_ak04its_phi[47];   //[njet_ak04its]
   Float_t         jet_ak04its_area_raw[47];   //[njet_ak04its]
   Float_t         jet_ak04its_area[47];   //[njet_ak04its]
   Float_t         jet_ak04its_emf_raw[47];   //[njet_ak04its]
   Float_t         jet_ak04its_emf[47];   //[njet_ak04its]
   UShort_t        jet_ak04its_multiplicity_raw[47];   //[njet_ak04its]
   Float_t         jet_ak04its_multiplicity[47];   //[njet_ak04its]
   Float_t         jet_ak04its_width_sigma_raw[47][2];   //[njet_ak04its]
   Float_t         jet_ak04its_width_sigma[47][2];   //[njet_ak04its]
   Float_t         jet_ak04its_ptd_raw[47];   //[njet_ak04its]
   Float_t         jet_ak04its_ptd[47];   //[njet_ak04its]
   Int_t           jet_ak04its_truth_index_z_truth[47][2];   //[njet_ak04its]
   Float_t         jet_ak04its_truth_z_truth[47][2];   //[njet_ak04its]
   Int_t           jet_ak04its_truth_index_z_reco[47][2];   //[njet_ak04its]
   Float_t         jet_ak04its_truth_z_reco[47][2];   //[njet_ak04its]
   Float_t         jet_ak04its_e_truth[47];   //[njet_ak04its]
   Float_t         jet_ak04its_pt_truth[47];   //[njet_ak04its]
   Float_t         jet_ak04its_eta_truth[47];   //[njet_ak04its]
   Float_t         jet_ak04its_phi_truth[47];   //[njet_ak04its]
   Float_t         jet_ak04its_area_truth[47];   //[njet_ak04its]
   Float_t         jet_ak04its_emf_truth[47];   //[njet_ak04its]
   UShort_t        jet_ak04its_multiplicity_truth[47];   //[njet_ak04its]
   Float_t         jet_ak04its_width_sigma_truth[47][2];   //[njet_ak04its]
   Float_t         jet_ak04its_ptd_truth[47];   //[njet_ak04its]
   UInt_t          njet_truth_ak04;
   Float_t         jet_truth_ak04_e[1];   //[njet_truth_ak04]
   Float_t         jet_truth_ak04_pt[1];   //[njet_truth_ak04]
   Float_t         jet_truth_ak04_eta[1];   //[njet_truth_ak04]
   Float_t         jet_truth_ak04_phi[1];   //[njet_truth_ak04]
   Float_t         jet_truth_ak04_area[1];   //[njet_truth_ak04]
   Float_t         jet_truth_ak04_emf[1];   //[njet_truth_ak04]
   UShort_t        jet_truth_ak04_multiplicity[1];   //[njet_truth_ak04]
   Float_t         jet_truth_ak04_width_sigma[1][2];   //[njet_truth_ak04]
   Float_t         jet_truth_ak04_ptd[1];   //[njet_truth_ak04]
   Double_t        met_tpc[2];
   Double_t        met_its[2];
   Double_t        met_truth[2];

   // List of branches
   TBranch        *b_id_git;   //!
   TBranch        *b_version_aliroot;   //!
   TBranch        *b_version_aliphysics;   //!
   TBranch        *b_version_jec;   //!
   TBranch        *b_grid_data_dir;   //!
   TBranch        *b_grid_data_pattern;   //!
   TBranch        *b_beam_particle;   //!
   TBranch        *b_ntrigger_class;   //!
   TBranch        *b_trigger_class;   //!
   TBranch        *b_run_number;   //!
   TBranch        *b_period_number;   //!
   TBranch        *b_orbit_number;   //!
   TBranch        *b_bunch_crossing_number;   //!
   TBranch        *b_time_stamp;   //!
   TBranch        *b_trigger_mask;   //!
   TBranch        *b_multiplicity_v0;   //!
   TBranch        *b_centrality_v0m;   //!
   TBranch        *b_centrality;   //!
   TBranch        *b_event_plane_psi_v0;   //!
   TBranch        *b_event_plane_q_v0;   //!
   TBranch        *b_has_misalignment_matrix;   //!
   TBranch        *b_cell_eta;   //!
   TBranch        *b_cell_phi;   //!
   TBranch        *b_cell_voronoi_area;   //!
   TBranch        *b_primary_vertex;   //!
   TBranch        *b_primary_vertex_sigma;   //!
   TBranch        *b_primary_vertex_ncontributor;   //!
   TBranch        *b_primary_vertex_spd;   //!
   TBranch        *b_primary_vertex_spd_sigma;   //!
   TBranch        *b_primary_vertex_spd_ncontributor;   //!
   TBranch        *b_npileup_vertex_spd;   //!
   TBranch        *b_pileup_vertex_spd_ncontributor;   //!
   TBranch        *b_pileup_vertex_spd_min_z_distance;   //!
   TBranch        *b_eg_signal_process_id;   //!
   TBranch        *b_eg_mpi;   //!
   TBranch        *b_eg_pt_hat;   //!
   TBranch        *b_eg_cross_section;   //!
   TBranch        *b_eg_weight;   //!
   TBranch        *b_eg_primary_vertex;   //!
   TBranch        *b_eg_ntrial;   //!
   TBranch        *b_eg_scale_pdf;   //!
   TBranch        *b_eg_alpha_qcd;   //!
   TBranch        *b_eg_alpha_qed;   //!
   TBranch        *b_eg_pdf_id;   //!
   TBranch        *b_eg_pdf_x;   //!
   TBranch        *b_eg_pdf_x_pdf;   //!
   TBranch        *b_debug_libmklml_gnu_loaded;   //!
   TBranch        *b_debug_libmklml_gnu_error;   //!
   TBranch        *b_ncluster;   //!
   TBranch        *b_cluster_e;   //!
   TBranch        *b_cluster_pt;   //!
   TBranch        *b_cluster_eta;   //!
   TBranch        *b_cluster_phi;   //!
   TBranch        *b_cluster_lambda_square;   //!
   TBranch        *b_cluster_tof;   //!
   TBranch        *b_cluster_ncell;   //!
   TBranch        *b_cluster_cell_id_max;   //!
   TBranch        *b_cluster_e_max;   //!
   TBranch        *b_cluster_e_cross;   //!
   TBranch        *b_cluster_nlocal_maxima;   //!
   TBranch        *b_cluster_nmc_truth;   //!
   TBranch        *b_cluster_mc_truth_index;   //!
   TBranch        *b_cluster_iso_tpc_01;   //!
   TBranch        *b_cluster_iso_tpc_01_ue;   //!
   TBranch        *b_cluster_iso_tpc_02;   //!
   TBranch        *b_cluster_iso_tpc_02_ue;   //!
   TBranch        *b_cluster_iso_tpc_03;   //!
   TBranch        *b_cluster_iso_tpc_03_ue;   //!
   TBranch        *b_cluster_iso_tpc_04;   //!
   TBranch        *b_cluster_iso_tpc_04_ue;   //!
   TBranch        *b_cluster_iso_its_01;   //!
   TBranch        *b_cluster_iso_its_01_ue;   //!
   TBranch        *b_cluster_iso_its_02;   //!
   TBranch        *b_cluster_iso_its_02_ue;   //!
   TBranch        *b_cluster_iso_its_03;   //!
   TBranch        *b_cluster_iso_its_03_ue;   //!
   TBranch        *b_cluster_iso_its_04;   //!
   TBranch        *b_cluster_iso_its_04_ue;   //!
   TBranch        *b_cluster_frixione_tpc_04_02;   //!
   TBranch        *b_cluster_frixione_tpc_04_05;   //!
   TBranch        *b_cluster_frixione_tpc_04_10;   //!
   TBranch        *b_cluster_frixione_its_04_02;   //!
   TBranch        *b_cluster_frixione_its_04_05;   //!
   TBranch        *b_cluster_frixione_its_04_10;   //!
   TBranch        *b_cluster_iso_01_truth;   //!
   TBranch        *b_cluster_iso_02_truth;   //!
   TBranch        *b_cluster_iso_03_truth;   //!
   TBranch        *b_cluster_iso_04_truth;   //!
   TBranch        *b_cluster_frixione_04_02_truth;   //!
   TBranch        *b_cluster_frixione_04_05_truth;   //!
   TBranch        *b_cluster_frixione_04_10_truth;   //!
   TBranch        *b_cluster_s_nphoton;   //!
   TBranch        *b_cluster_s_ncharged_hadron;   //!
   TBranch        *b_cell_e;   //!
   TBranch        *b_cell_tof;   //!
   TBranch        *b_cell_cluster_index;   //!
   TBranch        *b_cell_mc_truth_index;   //!
   TBranch        *b_ntrack;   //!
   TBranch        *b_track_e;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_eta_emcal;   //!
   TBranch        *b_track_phi_emcal;   //!
   TBranch        *b_track_charge;   //!
   TBranch        *b_track_quality;   //!
   TBranch        *b_track_tpc_dedx;   //!
   TBranch        *b_track_tpc_length_active_zone;   //!
   TBranch        *b_track_tpc_xrow;   //!
   TBranch        *b_track_tpc_ncluster;   //!
   TBranch        *b_track_tpc_ncluster_dedx;   //!
   TBranch        *b_track_tpc_ncluster_findable;   //!
   TBranch        *b_track_its_ncluster;   //!
   TBranch        *b_track_its_chi_square;   //!
   TBranch        *b_track_dca_xy;   //!
   TBranch        *b_track_dca_z;   //!
   TBranch        *b_track_mc_truth_index;   //!
   TBranch        *b_track_voronoi_area;   //!
   TBranch        *b_nmuon_track;   //!
   TBranch        *b_muon_track_e;   //!
   TBranch        *b_muon_track_pt;   //!
   TBranch        *b_muon_track_eta;   //!
   TBranch        *b_muon_track_phi;   //!
   TBranch        *b_muon_track_r_abs;   //!
   TBranch        *b_muon_track_p_dca;   //!
   TBranch        *b_muon_track_sigma_p_dca;   //!
   TBranch        *b_muon_track_delta_sagitta_p;   //!
   TBranch        *b_muon_track_distance_sigma_slope_p;   //!
   TBranch        *b_muon_track_mc_truth_index;   //!
   TBranch        *b_nmc_truth;   //!
   TBranch        *b_mc_truth_e;   //!
   TBranch        *b_mc_truth_pt;   //!
   TBranch        *b_mc_truth_eta;   //!
   TBranch        *b_mc_truth_phi;   //!
   TBranch        *b_mc_truth_charge;   //!
   TBranch        *b_mc_truth_pdg_code;   //!
   TBranch        *b_mc_truth_status;   //!
   TBranch        *b_mc_truth_generator_index;   //!
   TBranch        *b_mc_truth_first_parent_pdg_code;   //!
   TBranch        *b_mc_truth_first_parent_e;   //!
   TBranch        *b_mc_truth_first_parent_pt;   //!
   TBranch        *b_mc_truth_first_parent_eta;   //!
   TBranch        *b_mc_truth_first_parent_phi;   //!
   TBranch        *b_mc_truth_sibling_index;   //!
   TBranch        *b_debug_njet_ue_estimation;   //!
   TBranch        *b_debug_jet_ue_estimation_pt_raw;   //!
   TBranch        *b_debug_jet_ue_estimation_eta_raw;   //!
   TBranch        *b_debug_jet_ue_estimation_phi_raw;   //!
   TBranch        *b_debug_jet_ue_estimation_area_raw;   //!
   TBranch        *b_njet_ak04tpc;   //!
   TBranch        *b_debug_jet_ak04tpc_tag_dr_square;   //!
   TBranch        *b_jet_ak04tpc_e_raw;   //!
   TBranch        *b_jet_ak04tpc_e;   //!
   TBranch        *b_jet_ak04tpc_e_charged;   //!
   TBranch        *b_jet_ak04tpc_pt_raw_ue;   //!
   TBranch        *b_jet_ak04tpc_pt_raw;   //!
   TBranch        *b_jet_ak04tpc_pt;   //!
   TBranch        *b_jet_ak04tpc_pt_charged;   //!
   TBranch        *b_jet_ak04tpc_eta_raw;   //!
   TBranch        *b_jet_ak04tpc_eta;   //!
   TBranch        *b_jet_ak04tpc_phi;   //!
   TBranch        *b_jet_ak04tpc_area_raw;   //!
   TBranch        *b_jet_ak04tpc_area;   //!
   TBranch        *b_jet_ak04tpc_emf_raw;   //!
   TBranch        *b_jet_ak04tpc_emf;   //!
   TBranch        *b_jet_ak04tpc_multiplicity_raw;   //!
   TBranch        *b_jet_ak04tpc_multiplicity;   //!
   TBranch        *b_jet_ak04tpc_width_sigma_raw;   //!
   TBranch        *b_jet_ak04tpc_width_sigma;   //!
   TBranch        *b_jet_ak04tpc_ptd_raw;   //!
   TBranch        *b_jet_ak04tpc_ptd;   //!
   TBranch        *b_jet_ak04tpc_truth_index_z_truth;   //!
   TBranch        *b_jet_ak04tpc_truth_z_truth;   //!
   TBranch        *b_jet_ak04tpc_truth_index_z_reco;   //!
   TBranch        *b_jet_ak04tpc_truth_z_reco;   //!
   TBranch        *b_jet_ak04tpc_e_truth;   //!
   TBranch        *b_jet_ak04tpc_pt_truth;   //!
   TBranch        *b_jet_ak04tpc_eta_truth;   //!
   TBranch        *b_jet_ak04tpc_phi_truth;   //!
   TBranch        *b_jet_ak04tpc_area_truth;   //!
   TBranch        *b_jet_ak04tpc_emf_truth;   //!
   TBranch        *b_jet_ak04tpc_multiplicity_truth;   //!
   TBranch        *b_jet_ak04tpc_width_sigma_truth;   //!
   TBranch        *b_jet_ak04tpc_ptd_truth;   //!
   TBranch        *b_njet_ak04its;   //!
   TBranch        *b_debug_jet_ak04its_tag_dr_square;   //!
   TBranch        *b_jet_ak04its_e_raw;   //!
   TBranch        *b_jet_ak04its_e;   //!
   TBranch        *b_jet_ak04its_e_charged;   //!
   TBranch        *b_jet_ak04its_pt_raw_ue;   //!
   TBranch        *b_jet_ak04its_pt_raw;   //!
   TBranch        *b_jet_ak04its_pt;   //!
   TBranch        *b_jet_ak04its_pt_charged;   //!
   TBranch        *b_jet_ak04its_eta_raw;   //!
   TBranch        *b_jet_ak04its_eta;   //!
   TBranch        *b_jet_ak04its_phi;   //!
   TBranch        *b_jet_ak04its_area_raw;   //!
   TBranch        *b_jet_ak04its_area;   //!
   TBranch        *b_jet_ak04its_emf_raw;   //!
   TBranch        *b_jet_ak04its_emf;   //!
   TBranch        *b_jet_ak04its_multiplicity_raw;   //!
   TBranch        *b_jet_ak04its_multiplicity;   //!
   TBranch        *b_jet_ak04its_width_sigma_raw;   //!
   TBranch        *b_jet_ak04its_width_sigma;   //!
   TBranch        *b_jet_ak04its_ptd_raw;   //!
   TBranch        *b_jet_ak04its_ptd;   //!
   TBranch        *b_jet_ak04its_truth_index_z_truth;   //!
   TBranch        *b_jet_ak04its_truth_z_truth;   //!
   TBranch        *b_jet_ak04its_truth_index_z_reco;   //!
   TBranch        *b_jet_ak04its_truth_z_reco;   //!
   TBranch        *b_jet_ak04its_e_truth;   //!
   TBranch        *b_jet_ak04its_pt_truth;   //!
   TBranch        *b_jet_ak04its_eta_truth;   //!
   TBranch        *b_jet_ak04its_phi_truth;   //!
   TBranch        *b_jet_ak04its_area_truth;   //!
   TBranch        *b_jet_ak04its_emf_truth;   //!
   TBranch        *b_jet_ak04its_multiplicity_truth;   //!
   TBranch        *b_jet_ak04its_width_sigma_truth;   //!
   TBranch        *b_jet_ak04its_ptd_truth;   //!
   TBranch        *b_njet_truth_ak04;   //!
   TBranch        *b_jet_truth_ak04_e;   //!
   TBranch        *b_jet_truth_ak04_pt;   //!
   TBranch        *b_jet_truth_ak04_eta;   //!
   TBranch        *b_jet_truth_ak04_phi;   //!
   TBranch        *b_jet_truth_ak04_area;   //!
   TBranch        *b_jet_truth_ak04_emf;   //!
   TBranch        *b_jet_truth_ak04_multiplicity;   //!
   TBranch        *b_jet_truth_ak04_width_sigma;   //!
   TBranch        *b_jet_truth_ak04_ptd;   //!
   TBranch        *b_met_tpc;   //!
   TBranch        *b_met_its;   //!
   TBranch        *b_met_truth;   //!

   EtaMass(TTree *tree=0);
   virtual ~EtaMass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EtaMass_cxx
EtaMass::EtaMass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("13d_v4_small.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("13d_v4_small.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("13d_v4_small.root:/AliAnalysisTaskNTGJ");
      dir->GetObject("_tree_event",tree);

   }
   Init(tree);
}

EtaMass::~EtaMass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EtaMass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EtaMass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EtaMass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("id_git", &id_git, &b_id_git);
   fChain->SetBranchAddress("version_aliroot", &version_aliroot, &b_version_aliroot);
   fChain->SetBranchAddress("version_aliphysics", &version_aliphysics, &b_version_aliphysics);
   fChain->SetBranchAddress("version_jec", &version_jec, &b_version_jec);
   fChain->SetBranchAddress("grid_data_dir", &grid_data_dir, &b_grid_data_dir);
   fChain->SetBranchAddress("grid_data_pattern", &grid_data_pattern, &b_grid_data_pattern);
   fChain->SetBranchAddress("beam_particle", beam_particle, &b_beam_particle);
   fChain->SetBranchAddress("ntrigger_class", &ntrigger_class, &b_ntrigger_class);
   fChain->SetBranchAddress("trigger_class", &trigger_class, &b_trigger_class);
   fChain->SetBranchAddress("run_number", &run_number, &b_run_number);
   fChain->SetBranchAddress("period_number", &period_number, &b_period_number);
   fChain->SetBranchAddress("orbit_number", &orbit_number, &b_orbit_number);
   fChain->SetBranchAddress("bunch_crossing_number", &bunch_crossing_number, &b_bunch_crossing_number);
   fChain->SetBranchAddress("time_stamp", &time_stamp, &b_time_stamp);
   fChain->SetBranchAddress("trigger_mask", trigger_mask, &b_trigger_mask);
   fChain->SetBranchAddress("multiplicity_v0", multiplicity_v0, &b_multiplicity_v0);
   fChain->SetBranchAddress("centrality_v0m", &centrality_v0m, &b_centrality_v0m);
   fChain->SetBranchAddress("centrality", centrality, &b_centrality);
   fChain->SetBranchAddress("event_plane_psi_v0", event_plane_psi_v0, &b_event_plane_psi_v0);
   fChain->SetBranchAddress("event_plane_q_v0", event_plane_q_v0, &b_event_plane_q_v0);
   fChain->SetBranchAddress("has_misalignment_matrix", &has_misalignment_matrix, &b_has_misalignment_matrix);
   fChain->SetBranchAddress("cell_eta", cell_eta, &b_cell_eta);
   fChain->SetBranchAddress("cell_phi", cell_phi, &b_cell_phi);
   fChain->SetBranchAddress("cell_voronoi_area", cell_voronoi_area, &b_cell_voronoi_area);
   fChain->SetBranchAddress("primary_vertex", primary_vertex, &b_primary_vertex);
   fChain->SetBranchAddress("primary_vertex_sigma", primary_vertex_sigma, &b_primary_vertex_sigma);
   fChain->SetBranchAddress("primary_vertex_ncontributor", &primary_vertex_ncontributor, &b_primary_vertex_ncontributor);
   fChain->SetBranchAddress("primary_vertex_spd", primary_vertex_spd, &b_primary_vertex_spd);
   fChain->SetBranchAddress("primary_vertex_spd_sigma", primary_vertex_spd_sigma, &b_primary_vertex_spd_sigma);
   fChain->SetBranchAddress("primary_vertex_spd_ncontributor", &primary_vertex_spd_ncontributor, &b_primary_vertex_spd_ncontributor);
   fChain->SetBranchAddress("npileup_vertex_spd", &npileup_vertex_spd, &b_npileup_vertex_spd);
   fChain->SetBranchAddress("pileup_vertex_spd_ncontributor", &pileup_vertex_spd_ncontributor, &b_pileup_vertex_spd_ncontributor);
   fChain->SetBranchAddress("pileup_vertex_spd_min_z_distance", &pileup_vertex_spd_min_z_distance, &b_pileup_vertex_spd_min_z_distance);
   fChain->SetBranchAddress("eg_signal_process_id", &eg_signal_process_id, &b_eg_signal_process_id);
   fChain->SetBranchAddress("eg_mpi", &eg_mpi, &b_eg_mpi);
   fChain->SetBranchAddress("eg_pt_hat", &eg_pt_hat, &b_eg_pt_hat);
   fChain->SetBranchAddress("eg_cross_section", &eg_cross_section, &b_eg_cross_section);
   fChain->SetBranchAddress("eg_weight", &eg_weight, &b_eg_weight);
   fChain->SetBranchAddress("eg_primary_vertex", eg_primary_vertex, &b_eg_primary_vertex);
   fChain->SetBranchAddress("eg_ntrial", &eg_ntrial, &b_eg_ntrial);
   fChain->SetBranchAddress("eg_scale_pdf", &eg_scale_pdf, &b_eg_scale_pdf);
   fChain->SetBranchAddress("eg_alpha_qcd", &eg_alpha_qcd, &b_eg_alpha_qcd);
   fChain->SetBranchAddress("eg_alpha_qed", &eg_alpha_qed, &b_eg_alpha_qed);
   fChain->SetBranchAddress("eg_pdf_id", eg_pdf_id, &b_eg_pdf_id);
   fChain->SetBranchAddress("eg_pdf_x", eg_pdf_x, &b_eg_pdf_x);
   fChain->SetBranchAddress("eg_pdf_x_pdf", eg_pdf_x_pdf, &b_eg_pdf_x_pdf);
   fChain->SetBranchAddress("debug_libmklml_gnu_loaded", &debug_libmklml_gnu_loaded, &b_debug_libmklml_gnu_loaded);
   fChain->SetBranchAddress("debug_libmklml_gnu_error", debug_libmklml_gnu_error, &b_debug_libmklml_gnu_error);
   fChain->SetBranchAddress("ncluster", &ncluster, &b_ncluster);
   fChain->SetBranchAddress("cluster_e", cluster_e, &b_cluster_e);
   fChain->SetBranchAddress("cluster_pt", cluster_pt, &b_cluster_pt);
   fChain->SetBranchAddress("cluster_eta", cluster_eta, &b_cluster_eta);
   fChain->SetBranchAddress("cluster_phi", cluster_phi, &b_cluster_phi);
   fChain->SetBranchAddress("cluster_lambda_square", cluster_lambda_square, &b_cluster_lambda_square);
   fChain->SetBranchAddress("cluster_tof", cluster_tof, &b_cluster_tof);
   fChain->SetBranchAddress("cluster_ncell", cluster_ncell, &b_cluster_ncell);
   fChain->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max, &b_cluster_cell_id_max);
   fChain->SetBranchAddress("cluster_e_max", cluster_e_max, &b_cluster_e_max);
   fChain->SetBranchAddress("cluster_e_cross", cluster_e_cross, &b_cluster_e_cross);
   fChain->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima, &b_cluster_nlocal_maxima);
   fChain->SetBranchAddress("cluster_nmc_truth", cluster_nmc_truth, &b_cluster_nmc_truth);
   fChain->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index, &b_cluster_mc_truth_index);
   fChain->SetBranchAddress("cluster_iso_tpc_01", cluster_iso_tpc_01, &b_cluster_iso_tpc_01);
   fChain->SetBranchAddress("cluster_iso_tpc_01_ue", cluster_iso_tpc_01_ue, &b_cluster_iso_tpc_01_ue);
   fChain->SetBranchAddress("cluster_iso_tpc_02", cluster_iso_tpc_02, &b_cluster_iso_tpc_02);
   fChain->SetBranchAddress("cluster_iso_tpc_02_ue", cluster_iso_tpc_02_ue, &b_cluster_iso_tpc_02_ue);
   fChain->SetBranchAddress("cluster_iso_tpc_03", cluster_iso_tpc_03, &b_cluster_iso_tpc_03);
   fChain->SetBranchAddress("cluster_iso_tpc_03_ue", cluster_iso_tpc_03_ue, &b_cluster_iso_tpc_03_ue);
   fChain->SetBranchAddress("cluster_iso_tpc_04", cluster_iso_tpc_04, &b_cluster_iso_tpc_04);
   fChain->SetBranchAddress("cluster_iso_tpc_04_ue", cluster_iso_tpc_04_ue, &b_cluster_iso_tpc_04_ue);
   fChain->SetBranchAddress("cluster_iso_its_01", cluster_iso_its_01, &b_cluster_iso_its_01);
   fChain->SetBranchAddress("cluster_iso_its_01_ue", cluster_iso_its_01_ue, &b_cluster_iso_its_01_ue);
   fChain->SetBranchAddress("cluster_iso_its_02", cluster_iso_its_02, &b_cluster_iso_its_02);
   fChain->SetBranchAddress("cluster_iso_its_02_ue", cluster_iso_its_02_ue, &b_cluster_iso_its_02_ue);
   fChain->SetBranchAddress("cluster_iso_its_03", cluster_iso_its_03, &b_cluster_iso_its_03);
   fChain->SetBranchAddress("cluster_iso_its_03_ue", cluster_iso_its_03_ue, &b_cluster_iso_its_03_ue);
   fChain->SetBranchAddress("cluster_iso_its_04", cluster_iso_its_04, &b_cluster_iso_its_04);
   fChain->SetBranchAddress("cluster_iso_its_04_ue", cluster_iso_its_04_ue, &b_cluster_iso_its_04_ue);
   fChain->SetBranchAddress("cluster_frixione_tpc_04_02", cluster_frixione_tpc_04_02, &b_cluster_frixione_tpc_04_02);
   fChain->SetBranchAddress("cluster_frixione_tpc_04_05", cluster_frixione_tpc_04_05, &b_cluster_frixione_tpc_04_05);
   fChain->SetBranchAddress("cluster_frixione_tpc_04_10", cluster_frixione_tpc_04_10, &b_cluster_frixione_tpc_04_10);
   fChain->SetBranchAddress("cluster_frixione_its_04_02", cluster_frixione_its_04_02, &b_cluster_frixione_its_04_02);
   fChain->SetBranchAddress("cluster_frixione_its_04_05", cluster_frixione_its_04_05, &b_cluster_frixione_its_04_05);
   fChain->SetBranchAddress("cluster_frixione_its_04_10", cluster_frixione_its_04_10, &b_cluster_frixione_its_04_10);
   fChain->SetBranchAddress("cluster_iso_01_truth", cluster_iso_01_truth, &b_cluster_iso_01_truth);
   fChain->SetBranchAddress("cluster_iso_02_truth", cluster_iso_02_truth, &b_cluster_iso_02_truth);
   fChain->SetBranchAddress("cluster_iso_03_truth", cluster_iso_03_truth, &b_cluster_iso_03_truth);
   fChain->SetBranchAddress("cluster_iso_04_truth", cluster_iso_04_truth, &b_cluster_iso_04_truth);
   fChain->SetBranchAddress("cluster_frixione_04_02_truth", cluster_frixione_04_02_truth, &b_cluster_frixione_04_02_truth);
   fChain->SetBranchAddress("cluster_frixione_04_05_truth", cluster_frixione_04_05_truth, &b_cluster_frixione_04_05_truth);
   fChain->SetBranchAddress("cluster_frixione_04_10_truth", cluster_frixione_04_10_truth, &b_cluster_frixione_04_10_truth);
   fChain->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton, &b_cluster_s_nphoton);
   fChain->SetBranchAddress("cluster_s_ncharged_hadron", cluster_s_ncharged_hadron, &b_cluster_s_ncharged_hadron);
   fChain->SetBranchAddress("cell_e", cell_e, &b_cell_e);
   fChain->SetBranchAddress("cell_tof", cell_tof, &b_cell_tof);
   fChain->SetBranchAddress("cell_cluster_index", cell_cluster_index, &b_cell_cluster_index);
   fChain->SetBranchAddress("cell_mc_truth_index", cell_mc_truth_index, &b_cell_mc_truth_index);
   fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
   fChain->SetBranchAddress("track_e", track_e, &b_track_e);
   fChain->SetBranchAddress("track_pt", track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_eta", track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", track_phi, &b_track_phi);
   fChain->SetBranchAddress("track_eta_emcal", track_eta_emcal, &b_track_eta_emcal);
   fChain->SetBranchAddress("track_phi_emcal", track_phi_emcal, &b_track_phi_emcal);
   fChain->SetBranchAddress("track_charge", track_charge, &b_track_charge);
   fChain->SetBranchAddress("track_quality", track_quality, &b_track_quality);
   fChain->SetBranchAddress("track_tpc_dedx", track_tpc_dedx, &b_track_tpc_dedx);
   fChain->SetBranchAddress("track_tpc_length_active_zone", track_tpc_length_active_zone, &b_track_tpc_length_active_zone);
   fChain->SetBranchAddress("track_tpc_xrow", track_tpc_xrow, &b_track_tpc_xrow);
   fChain->SetBranchAddress("track_tpc_ncluster", track_tpc_ncluster, &b_track_tpc_ncluster);
   fChain->SetBranchAddress("track_tpc_ncluster_dedx", track_tpc_ncluster_dedx, &b_track_tpc_ncluster_dedx);
   fChain->SetBranchAddress("track_tpc_ncluster_findable", track_tpc_ncluster_findable, &b_track_tpc_ncluster_findable);
   fChain->SetBranchAddress("track_its_ncluster", track_its_ncluster, &b_track_its_ncluster);
   fChain->SetBranchAddress("track_its_chi_square", track_its_chi_square, &b_track_its_chi_square);
   fChain->SetBranchAddress("track_dca_xy", track_dca_xy, &b_track_dca_xy);
   fChain->SetBranchAddress("track_dca_z", track_dca_z, &b_track_dca_z);
   fChain->SetBranchAddress("track_mc_truth_index", track_mc_truth_index, &b_track_mc_truth_index);
   fChain->SetBranchAddress("track_voronoi_area", track_voronoi_area, &b_track_voronoi_area);
   fChain->SetBranchAddress("nmuon_track", &nmuon_track, &b_nmuon_track);
   fChain->SetBranchAddress("muon_track_e", muon_track_e, &b_muon_track_e);
   fChain->SetBranchAddress("muon_track_pt", muon_track_pt, &b_muon_track_pt);
   fChain->SetBranchAddress("muon_track_eta", muon_track_eta, &b_muon_track_eta);
   fChain->SetBranchAddress("muon_track_phi", muon_track_phi, &b_muon_track_phi);
   fChain->SetBranchAddress("muon_track_r_abs", muon_track_r_abs, &b_muon_track_r_abs);
   fChain->SetBranchAddress("muon_track_p_dca", muon_track_p_dca, &b_muon_track_p_dca);
   fChain->SetBranchAddress("muon_track_sigma_p_dca", muon_track_sigma_p_dca, &b_muon_track_sigma_p_dca);
   fChain->SetBranchAddress("muon_track_delta_sagitta_p", muon_track_delta_sagitta_p, &b_muon_track_delta_sagitta_p);
   fChain->SetBranchAddress("muon_track_distance_sigma_slope_p", muon_track_distance_sigma_slope_p, &b_muon_track_distance_sigma_slope_p);
   fChain->SetBranchAddress("muon_track_mc_truth_index", muon_track_mc_truth_index, &b_muon_track_mc_truth_index);
   fChain->SetBranchAddress("nmc_truth", &nmc_truth, &b_nmc_truth);
   fChain->SetBranchAddress("mc_truth_e", &mc_truth_e, &b_mc_truth_e);
   fChain->SetBranchAddress("mc_truth_pt", &mc_truth_pt, &b_mc_truth_pt);
   fChain->SetBranchAddress("mc_truth_eta", &mc_truth_eta, &b_mc_truth_eta);
   fChain->SetBranchAddress("mc_truth_phi", &mc_truth_phi, &b_mc_truth_phi);
   fChain->SetBranchAddress("mc_truth_charge", &mc_truth_charge, &b_mc_truth_charge);
   fChain->SetBranchAddress("mc_truth_pdg_code", &mc_truth_pdg_code, &b_mc_truth_pdg_code);
   fChain->SetBranchAddress("mc_truth_status", &mc_truth_status, &b_mc_truth_status);
   fChain->SetBranchAddress("mc_truth_generator_index", &mc_truth_generator_index, &b_mc_truth_generator_index);
   fChain->SetBranchAddress("mc_truth_first_parent_pdg_code", &mc_truth_first_parent_pdg_code, &b_mc_truth_first_parent_pdg_code);
   fChain->SetBranchAddress("mc_truth_first_parent_e", &mc_truth_first_parent_e, &b_mc_truth_first_parent_e);
   fChain->SetBranchAddress("mc_truth_first_parent_pt", &mc_truth_first_parent_pt, &b_mc_truth_first_parent_pt);
   fChain->SetBranchAddress("mc_truth_first_parent_eta", &mc_truth_first_parent_eta, &b_mc_truth_first_parent_eta);
   fChain->SetBranchAddress("mc_truth_first_parent_phi", &mc_truth_first_parent_phi, &b_mc_truth_first_parent_phi);
   fChain->SetBranchAddress("mc_truth_sibling_index", &mc_truth_sibling_index, &b_mc_truth_sibling_index);
   fChain->SetBranchAddress("debug_njet_ue_estimation", &debug_njet_ue_estimation, &b_debug_njet_ue_estimation);
   fChain->SetBranchAddress("debug_jet_ue_estimation_pt_raw", debug_jet_ue_estimation_pt_raw, &b_debug_jet_ue_estimation_pt_raw);
   fChain->SetBranchAddress("debug_jet_ue_estimation_eta_raw", debug_jet_ue_estimation_eta_raw, &b_debug_jet_ue_estimation_eta_raw);
   fChain->SetBranchAddress("debug_jet_ue_estimation_phi_raw", debug_jet_ue_estimation_phi_raw, &b_debug_jet_ue_estimation_phi_raw);
   fChain->SetBranchAddress("debug_jet_ue_estimation_area_raw", debug_jet_ue_estimation_area_raw, &b_debug_jet_ue_estimation_area_raw);
   fChain->SetBranchAddress("njet_ak04tpc", &njet_ak04tpc, &b_njet_ak04tpc);
   fChain->SetBranchAddress("debug_jet_ak04tpc_tag_dr_square", debug_jet_ak04tpc_tag_dr_square, &b_debug_jet_ak04tpc_tag_dr_square);
   fChain->SetBranchAddress("jet_ak04tpc_e_raw", jet_ak04tpc_e_raw, &b_jet_ak04tpc_e_raw);
   fChain->SetBranchAddress("jet_ak04tpc_e", jet_ak04tpc_e, &b_jet_ak04tpc_e);
   fChain->SetBranchAddress("jet_ak04tpc_e_charged", jet_ak04tpc_e_charged, &b_jet_ak04tpc_e_charged);
   fChain->SetBranchAddress("jet_ak04tpc_pt_raw_ue", jet_ak04tpc_pt_raw_ue, &b_jet_ak04tpc_pt_raw_ue);
   fChain->SetBranchAddress("jet_ak04tpc_pt_raw", jet_ak04tpc_pt_raw, &b_jet_ak04tpc_pt_raw);
   fChain->SetBranchAddress("jet_ak04tpc_pt", jet_ak04tpc_pt, &b_jet_ak04tpc_pt);
   fChain->SetBranchAddress("jet_ak04tpc_pt_charged", jet_ak04tpc_pt_charged, &b_jet_ak04tpc_pt_charged);
   fChain->SetBranchAddress("jet_ak04tpc_eta_raw", jet_ak04tpc_eta_raw, &b_jet_ak04tpc_eta_raw);
   fChain->SetBranchAddress("jet_ak04tpc_eta", jet_ak04tpc_eta, &b_jet_ak04tpc_eta);
   fChain->SetBranchAddress("jet_ak04tpc_phi", jet_ak04tpc_phi, &b_jet_ak04tpc_phi);
   fChain->SetBranchAddress("jet_ak04tpc_area_raw", jet_ak04tpc_area_raw, &b_jet_ak04tpc_area_raw);
   fChain->SetBranchAddress("jet_ak04tpc_area", jet_ak04tpc_area, &b_jet_ak04tpc_area);
   fChain->SetBranchAddress("jet_ak04tpc_emf_raw", jet_ak04tpc_emf_raw, &b_jet_ak04tpc_emf_raw);
   fChain->SetBranchAddress("jet_ak04tpc_emf", jet_ak04tpc_emf, &b_jet_ak04tpc_emf);
   fChain->SetBranchAddress("jet_ak04tpc_multiplicity_raw", jet_ak04tpc_multiplicity_raw, &b_jet_ak04tpc_multiplicity_raw);
   fChain->SetBranchAddress("jet_ak04tpc_multiplicity", jet_ak04tpc_multiplicity, &b_jet_ak04tpc_multiplicity);
   fChain->SetBranchAddress("jet_ak04tpc_width_sigma_raw", jet_ak04tpc_width_sigma_raw, &b_jet_ak04tpc_width_sigma_raw);
   fChain->SetBranchAddress("jet_ak04tpc_width_sigma", jet_ak04tpc_width_sigma, &b_jet_ak04tpc_width_sigma);
   fChain->SetBranchAddress("jet_ak04tpc_ptd_raw", jet_ak04tpc_ptd_raw, &b_jet_ak04tpc_ptd_raw);
   fChain->SetBranchAddress("jet_ak04tpc_ptd", jet_ak04tpc_ptd, &b_jet_ak04tpc_ptd);
   fChain->SetBranchAddress("jet_ak04tpc_truth_index_z_truth", jet_ak04tpc_truth_index_z_truth, &b_jet_ak04tpc_truth_index_z_truth);
   fChain->SetBranchAddress("jet_ak04tpc_truth_z_truth", jet_ak04tpc_truth_z_truth, &b_jet_ak04tpc_truth_z_truth);
   fChain->SetBranchAddress("jet_ak04tpc_truth_index_z_reco", jet_ak04tpc_truth_index_z_reco, &b_jet_ak04tpc_truth_index_z_reco);
   fChain->SetBranchAddress("jet_ak04tpc_truth_z_reco", jet_ak04tpc_truth_z_reco, &b_jet_ak04tpc_truth_z_reco);
   fChain->SetBranchAddress("jet_ak04tpc_e_truth", jet_ak04tpc_e_truth, &b_jet_ak04tpc_e_truth);
   fChain->SetBranchAddress("jet_ak04tpc_pt_truth", jet_ak04tpc_pt_truth, &b_jet_ak04tpc_pt_truth);
   fChain->SetBranchAddress("jet_ak04tpc_eta_truth", jet_ak04tpc_eta_truth, &b_jet_ak04tpc_eta_truth);
   fChain->SetBranchAddress("jet_ak04tpc_phi_truth", jet_ak04tpc_phi_truth, &b_jet_ak04tpc_phi_truth);
   fChain->SetBranchAddress("jet_ak04tpc_area_truth", jet_ak04tpc_area_truth, &b_jet_ak04tpc_area_truth);
   fChain->SetBranchAddress("jet_ak04tpc_emf_truth", jet_ak04tpc_emf_truth, &b_jet_ak04tpc_emf_truth);
   fChain->SetBranchAddress("jet_ak04tpc_multiplicity_truth", jet_ak04tpc_multiplicity_truth, &b_jet_ak04tpc_multiplicity_truth);
   fChain->SetBranchAddress("jet_ak04tpc_width_sigma_truth", jet_ak04tpc_width_sigma_truth, &b_jet_ak04tpc_width_sigma_truth);
   fChain->SetBranchAddress("jet_ak04tpc_ptd_truth", jet_ak04tpc_ptd_truth, &b_jet_ak04tpc_ptd_truth);
   fChain->SetBranchAddress("njet_ak04its", &njet_ak04its, &b_njet_ak04its);
   fChain->SetBranchAddress("debug_jet_ak04its_tag_dr_square", debug_jet_ak04its_tag_dr_square, &b_debug_jet_ak04its_tag_dr_square);
   fChain->SetBranchAddress("jet_ak04its_e_raw", jet_ak04its_e_raw, &b_jet_ak04its_e_raw);
   fChain->SetBranchAddress("jet_ak04its_e", jet_ak04its_e, &b_jet_ak04its_e);
   fChain->SetBranchAddress("jet_ak04its_e_charged", jet_ak04its_e_charged, &b_jet_ak04its_e_charged);
   fChain->SetBranchAddress("jet_ak04its_pt_raw_ue", jet_ak04its_pt_raw_ue, &b_jet_ak04its_pt_raw_ue);
   fChain->SetBranchAddress("jet_ak04its_pt_raw", jet_ak04its_pt_raw, &b_jet_ak04its_pt_raw);
   fChain->SetBranchAddress("jet_ak04its_pt", jet_ak04its_pt, &b_jet_ak04its_pt);
   fChain->SetBranchAddress("jet_ak04its_pt_charged", jet_ak04its_pt_charged, &b_jet_ak04its_pt_charged);
   fChain->SetBranchAddress("jet_ak04its_eta_raw", jet_ak04its_eta_raw, &b_jet_ak04its_eta_raw);
   fChain->SetBranchAddress("jet_ak04its_eta", jet_ak04its_eta, &b_jet_ak04its_eta);
   fChain->SetBranchAddress("jet_ak04its_phi", jet_ak04its_phi, &b_jet_ak04its_phi);
   fChain->SetBranchAddress("jet_ak04its_area_raw", jet_ak04its_area_raw, &b_jet_ak04its_area_raw);
   fChain->SetBranchAddress("jet_ak04its_area", jet_ak04its_area, &b_jet_ak04its_area);
   fChain->SetBranchAddress("jet_ak04its_emf_raw", jet_ak04its_emf_raw, &b_jet_ak04its_emf_raw);
   fChain->SetBranchAddress("jet_ak04its_emf", jet_ak04its_emf, &b_jet_ak04its_emf);
   fChain->SetBranchAddress("jet_ak04its_multiplicity_raw", jet_ak04its_multiplicity_raw, &b_jet_ak04its_multiplicity_raw);
   fChain->SetBranchAddress("jet_ak04its_multiplicity", jet_ak04its_multiplicity, &b_jet_ak04its_multiplicity);
   fChain->SetBranchAddress("jet_ak04its_width_sigma_raw", jet_ak04its_width_sigma_raw, &b_jet_ak04its_width_sigma_raw);
   fChain->SetBranchAddress("jet_ak04its_width_sigma", jet_ak04its_width_sigma, &b_jet_ak04its_width_sigma);
   fChain->SetBranchAddress("jet_ak04its_ptd_raw", jet_ak04its_ptd_raw, &b_jet_ak04its_ptd_raw);
   fChain->SetBranchAddress("jet_ak04its_ptd", jet_ak04its_ptd, &b_jet_ak04its_ptd);
   fChain->SetBranchAddress("jet_ak04its_truth_index_z_truth", jet_ak04its_truth_index_z_truth, &b_jet_ak04its_truth_index_z_truth);
   fChain->SetBranchAddress("jet_ak04its_truth_z_truth", jet_ak04its_truth_z_truth, &b_jet_ak04its_truth_z_truth);
   fChain->SetBranchAddress("jet_ak04its_truth_index_z_reco", jet_ak04its_truth_index_z_reco, &b_jet_ak04its_truth_index_z_reco);
   fChain->SetBranchAddress("jet_ak04its_truth_z_reco", jet_ak04its_truth_z_reco, &b_jet_ak04its_truth_z_reco);
   fChain->SetBranchAddress("jet_ak04its_e_truth", jet_ak04its_e_truth, &b_jet_ak04its_e_truth);
   fChain->SetBranchAddress("jet_ak04its_pt_truth", jet_ak04its_pt_truth, &b_jet_ak04its_pt_truth);
   fChain->SetBranchAddress("jet_ak04its_eta_truth", jet_ak04its_eta_truth, &b_jet_ak04its_eta_truth);
   fChain->SetBranchAddress("jet_ak04its_phi_truth", jet_ak04its_phi_truth, &b_jet_ak04its_phi_truth);
   fChain->SetBranchAddress("jet_ak04its_area_truth", jet_ak04its_area_truth, &b_jet_ak04its_area_truth);
   fChain->SetBranchAddress("jet_ak04its_emf_truth", jet_ak04its_emf_truth, &b_jet_ak04its_emf_truth);
   fChain->SetBranchAddress("jet_ak04its_multiplicity_truth", jet_ak04its_multiplicity_truth, &b_jet_ak04its_multiplicity_truth);
   fChain->SetBranchAddress("jet_ak04its_width_sigma_truth", jet_ak04its_width_sigma_truth, &b_jet_ak04its_width_sigma_truth);
   fChain->SetBranchAddress("jet_ak04its_ptd_truth", jet_ak04its_ptd_truth, &b_jet_ak04its_ptd_truth);
   fChain->SetBranchAddress("njet_truth_ak04", &njet_truth_ak04, &b_njet_truth_ak04);
   fChain->SetBranchAddress("jet_truth_ak04_e", &jet_truth_ak04_e, &b_jet_truth_ak04_e);
   fChain->SetBranchAddress("jet_truth_ak04_pt", &jet_truth_ak04_pt, &b_jet_truth_ak04_pt);
   fChain->SetBranchAddress("jet_truth_ak04_eta", &jet_truth_ak04_eta, &b_jet_truth_ak04_eta);
   fChain->SetBranchAddress("jet_truth_ak04_phi", &jet_truth_ak04_phi, &b_jet_truth_ak04_phi);
   fChain->SetBranchAddress("jet_truth_ak04_area", &jet_truth_ak04_area, &b_jet_truth_ak04_area);
   fChain->SetBranchAddress("jet_truth_ak04_emf", &jet_truth_ak04_emf, &b_jet_truth_ak04_emf);
   fChain->SetBranchAddress("jet_truth_ak04_multiplicity", &jet_truth_ak04_multiplicity, &b_jet_truth_ak04_multiplicity);
   fChain->SetBranchAddress("jet_truth_ak04_width_sigma", &jet_truth_ak04_width_sigma, &b_jet_truth_ak04_width_sigma);
   fChain->SetBranchAddress("jet_truth_ak04_ptd", &jet_truth_ak04_ptd, &b_jet_truth_ak04_ptd);
   fChain->SetBranchAddress("met_tpc", met_tpc, &b_met_tpc);
   fChain->SetBranchAddress("met_its", met_its, &b_met_its);
   fChain->SetBranchAddress("met_truth", met_truth, &b_met_truth);
   Notify();
}

Bool_t EtaMass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EtaMass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EtaMass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EtaMass_cxx
