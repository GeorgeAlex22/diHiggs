#define MASSBQUARK      4.8

double Energy(double px, double py, double pz, double m)
{
    double E = pow(px*px + py*py +pz*pz + m*m, 0.5);
    return E;
}

double Energy2(double px, double py, double pz, double m)
{
    double E2 = px*px + py*py +pz*pz + m*m;
    return E2;
}

//--------------------------------------------------------//


void higgsReco()
{
    TFile* f = new TFile("hhDilep_pseudoanasol/tree_1970.root");

    TTree* t1 = (TTree *) f->Get("t1");

    Double_t MET_x0, MET_y0;

    Double_t jet1_px0, jet1_py0, jet1_pz0, jet1_E0;                 
    Double_t jet2_px0, jet2_py0, jet2_pz0, jet2_E0;           

    Double_t lepton1_px0, lepton1_py0, lepton1_pz0, lepton1_E0;
    Int_t lepton1_charge0;

    Double_t lepton2_px0, lepton2_py0, lepton2_pz0, lepton2_E0;
    Int_t lepton2_charge0;

    Double_t neutrino1_px0, neutrino1_py0, neutrino1_pz0, neutrino1_E0;
    Double_t neutrino2_px0, neutrino2_py0, neutrino2_pz0, neutrino2_E0;
    Double_t neutrino1_px_sol, neutrino1_py_sol, neutrino1_pz_sol, neutrino1_E_sol;
    Double_t neutrino2_px_sol, neutrino2_py_sol, neutrino2_pz_sol, neutrino2_E_sol;

    Double_t result_WStarmass;    

    t1->SetBranchAddress("MET_x0", &MET_x0); 
    t1->SetBranchAddress("MET_y0", &MET_y0);
    t1->SetBranchAddress("jet1_px0", &jet1_px0);
    t1->SetBranchAddress("jet1_px0", &jet1_px0);
    t1->SetBranchAddress("jet1_py0", &jet1_py0);
    t1->SetBranchAddress("jet1_py0", &jet1_py0);
    t1->SetBranchAddress("jet1_pz0", &jet1_pz0);
    t1->SetBranchAddress("jet1_E0" , &jet1_E0 );
    t1->SetBranchAddress("jet2_px0", &jet2_px0);
    t1->SetBranchAddress("jet2_py0", &jet2_py0);
    t1->SetBranchAddress("jet2_pz0", &jet2_pz0);
    t1->SetBranchAddress("jet2_E0" , &jet2_E0 );
    t1->SetBranchAddress("lepton1_px0", &lepton1_px0);
    t1->SetBranchAddress("lepton1_py0", &lepton1_py0);
    t1->SetBranchAddress("lepton1_pz0", &lepton1_pz0);
    t1->SetBranchAddress("lepton1_E0" , &lepton1_E0 );
    t1->SetBranchAddress("lepton2_px0", &lepton2_px0);
    t1->SetBranchAddress("lepton2_py0", &lepton2_py0);
    t1->SetBranchAddress("lepton2_pz0", &lepton2_pz0);
    t1->SetBranchAddress("lepton2_E0" , &lepton2_E0 );
    t1->SetBranchAddress("neutrino1_px0", &neutrino1_px0);
    t1->SetBranchAddress("neutrino1_py0", &neutrino1_py0);
    t1->SetBranchAddress("neutrino1_pz0", &neutrino1_pz0);
    t1->SetBranchAddress("neutrino1_E0" , &neutrino1_E0 );
    t1->SetBranchAddress("neutrino2_px0", &neutrino2_px0);
    t1->SetBranchAddress("neutrino2_py0", &neutrino2_py0);
    t1->SetBranchAddress("neutrino2_pz0", &neutrino2_pz0);
    t1->SetBranchAddress("neutrino2_E0" , &neutrino2_E0 );
    t1->SetBranchAddress("neutrino1_px_sol", &neutrino1_px_sol);
    t1->SetBranchAddress("neutrino1_py_sol", &neutrino1_py_sol);
    t1->SetBranchAddress("neutrino1_pz_sol", &neutrino1_pz_sol);
    t1->SetBranchAddress("neutrino1_E_sol" , &neutrino1_E_sol );
    t1->SetBranchAddress("neutrino2_px_sol", &neutrino2_px_sol);
    t1->SetBranchAddress("neutrino2_py_sol", &neutrino2_py_sol);
    t1->SetBranchAddress("neutrino2_pz_sol", &neutrino2_pz_sol);
    t1->SetBranchAddress("neutrino2_E_sol" , &neutrino2_E_sol );

    Double_t allEntries = t1->GetEntries();
    int entry;

    TH1D* higgs_mass_true_hist = new TH1D("higgs_mass_true_hist", "", 100, 0., 0.);
    TH1D* higgs_mass_sol_hist = new TH1D("higgs_mass_sol_hist", "", 100, 0., 0.);

    TLorentzVector jet1, jet2, lepton1, lepton2, neutrino1_true , neutrino2_true, neutrino1_sol, neutrino2_sol;
    /* Event Loop */
    for(entry = 0; entry < allEntries; entry++)
    {
        t1->GetEntry(entry);

        jet1.SetPxPyPzE(jet1_px0, jet1_py0, jet1_pz0, Energy(jet1_px0, jet1_py0, jet1_pz0, MASSBQUARK));
        jet2.SetPxPyPzE(jet2_px0, jet2_py0, jet2_pz0, Energy(jet2_px0, jet2_py0, jet2_pz0, MASSBQUARK));
        lepton1.SetPxPyPzE(lepton1_px0, lepton1_py0, lepton1_pz0, Energy(lepton1_px0, lepton1_py0, lepton1_pz0, 0.0));
        lepton2.SetPxPyPzE(lepton2_px0, lepton2_py0, lepton2_pz0, Energy(lepton2_px0, lepton2_py0, lepton2_pz0, 0.0));
        neutrino1_true.SetPxPyPzE(neutrino1_px0, neutrino1_py0, neutrino1_pz0, Energy(neutrino1_px0, neutrino1_py0, neutrino1_pz0, 0.0));
        neutrino2_true.SetPxPyPzE(neutrino2_px0, neutrino2_py0, neutrino2_pz0, Energy(neutrino2_px0, neutrino2_py0, neutrino2_pz0, 0.0));
        neutrino1_sol.SetPxPyPzE(neutrino1_px_sol, neutrino1_py_sol, neutrino1_pz_sol, Energy(neutrino1_px_sol, neutrino1_py_sol, neutrino1_pz_sol, 0.0));
        neutrino2_sol.SetPxPyPzE(neutrino2_px_sol, neutrino2_py_sol, neutrino2_pz_sol, Energy(neutrino2_px_sol, neutrino2_py_sol, neutrino2_pz_sol, 0.0));

        Double_t higgs_mass_true = (lepton1 + lepton2 + neutrino1_true + neutrino2_true).M();
        Double_t higgs_mass_sol = (lepton1 + lepton2 + neutrino1_sol + neutrino2_sol).M();
        higgs_mass_true_hist->Fill(higgs_mass_true);
        higgs_mass_sol_hist->Fill(higgs_mass_sol);

    }/* Event Loop */

    char buffer1[20];
    Int_t true_bin_width = higgs_mass_true_hist->GetBinWidth(1);
    sprintf(buffer1, "Events/%dGeV", true_bin_width);

    char buffer2[20];
    Int_t sol_bin_width = higgs_mass_sol_hist->GetBinWidth(1);
    sprintf(buffer2, "Events/%dGeV", sol_bin_width);

    higgs_mass_true_hist->GetXaxis()->SetTitle("M(ll#nu#nu)");
    higgs_mass_true_hist->GetYaxis()->SetTitle(buffer1);
    higgs_mass_sol_hist->GetXaxis()->SetTitle("M(ll#nu#nu)");
    higgs_mass_sol_hist->GetYaxis()->SetTitle(buffer2);

    TFile* f_out = new TFile("higgs_histos.root", "recreate");
    higgs_mass_true_hist->Write();
    higgs_mass_sol_hist->Write();

    gROOT->ProcessLine(".q");
}