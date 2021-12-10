#include <dirent.h>

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


const bool generatedInfo = !false;

void createHistos(const char* physics_process)
{
    TChain* chain = new TChain("t1");
    // TFile* f = new TFile("hhDilep_pseudoanasol/tree_1970.root");

    DIR* dir;
    struct dirent* ent;
    char eos_path[200];
    sprintf(eos_path, "/eos/user/a/alexandg/public/EOS.diHiggs/rootFiles/%s/trees/", physics_process);
    char file[100];
    char filename[200];
    TFile* tmp_file;
    TTree* tmp_tree;
    Long64_t allEntries;
    int i=0;
    if ((dir = opendir (eos_path)) != NULL) 
    {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL)
        {
            strcpy(file, ent->d_name);
            if( ((string)file != ".") && ((string)file != ".."))
            {
                strcpy(filename, eos_path);
                strcat(filename, file);
                tmp_file = new TFile(filename);
                tmp_tree = (TTree*) tmp_file->Get("t1");
                allEntries = tmp_tree->GetEntries();
                if(allEntries > 0) 
                {
                    cout << "**Adding " << filename << " to chain.." << endl;
                    chain->Add(filename);
                    i++;
                    if(i>900){break;}
                }   
            }
        }
        closedir (dir);
	} 

    else 
    {
        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
	}	


    // TTree* t1 = (TTree *) f->Get("t1");

    // Double_t MET_x0, MET_y0;

    // Double_t jet1_px0, jet1_py0, jet1_pz0, jet1_E0;                 
    // Double_t jet2_px0, jet2_py0, jet2_pz0, jet2_E0;           

    // Double_t lepton1_px0, lepton1_py0, lepton1_pz0, lepton1_E0;
    // Int_t lepton1_charge0;

    // Double_t lepton2_px0, lepton2_py0, lepton2_pz0, lepton2_E0;
    // Int_t lepton2_charge0;

    // Double_t neutrino1_px0, neutrino1_py0, neutrino1_pz0, neutrino1_E0;
    // Double_t neutrino2_px0, neutrino2_py0, neutrino2_pz0, neutrino2_E0;

    Double_t MET_x, MET_y;

    Double_t jet1_px, jet1_py, jet1_pz, jet1_E;                 
    Double_t jet2_px, jet2_py, jet2_pz, jet2_E;                 

    Double_t lepton1_px, lepton1_py, lepton1_pz, lepton1_E;
    Int_t lepton1_charge;

    Double_t lepton2_px, lepton2_py, lepton2_pz, lepton2_E;
    Int_t lepton2_charge;

    Double_t neutrino1_px, neutrino1_py, neutrino1_pz, neutrino1_E;
    Double_t neutrino2_px, neutrino2_py, neutrino2_pz, neutrino2_E;
    Double_t neutrino1_px_sol, neutrino1_py_sol, neutrino1_pz_sol, neutrino1_E_sol;
    Double_t neutrino2_px_sol, neutrino2_py_sol, neutrino2_pz_sol, neutrino2_E_sol;

    Double_t WStarmass_sol;

    TTreeReader* treeReader = new TTreeReader(chain);
    TTree* t1 = treeReader->GetTree();

    t1->SetBranchAddress("MET_x", &MET_x); 
    t1->SetBranchAddress("MET_y", &MET_y);
    t1->SetBranchAddress("jet1_px", &jet1_px);
    t1->SetBranchAddress("jet1_py", &jet1_py);
    t1->SetBranchAddress("jet1_pz", &jet1_pz);
    t1->SetBranchAddress("jet1_E" , &jet1_E );
    t1->SetBranchAddress("jet2_px", &jet2_px);
    t1->SetBranchAddress("jet2_py", &jet2_py);
    t1->SetBranchAddress("jet2_pz", &jet2_pz);
    t1->SetBranchAddress("jet2_E" , &jet2_E );
    t1->SetBranchAddress("lepton1_px", &lepton1_px);
    t1->SetBranchAddress("lepton1_py", &lepton1_py);
    t1->SetBranchAddress("lepton1_pz", &lepton1_pz);
    t1->SetBranchAddress("lepton1_E" , &lepton1_E );
    t1->SetBranchAddress("lepton2_px", &lepton2_px);
    t1->SetBranchAddress("lepton2_py", &lepton2_py);
    t1->SetBranchAddress("lepton2_pz", &lepton2_pz);
    t1->SetBranchAddress("lepton2_E" , &lepton2_E );

    if(generatedInfo)
    {
        t1->SetBranchAddress("neutrino1_px", &neutrino1_px);
        t1->SetBranchAddress("neutrino1_py", &neutrino1_py);
        t1->SetBranchAddress("neutrino1_pz", &neutrino1_pz);
        t1->SetBranchAddress("neutrino1_E" , &neutrino1_E );
        t1->SetBranchAddress("neutrino2_px", &neutrino2_px);
        t1->SetBranchAddress("neutrino2_py", &neutrino2_py);
        t1->SetBranchAddress("neutrino2_pz", &neutrino2_pz);
        t1->SetBranchAddress("neutrino2_E" , &neutrino2_E );
    }

    t1->SetBranchAddress("neutrino1_px_sol", &neutrino1_px_sol);
    t1->SetBranchAddress("neutrino1_py_sol", &neutrino1_py_sol);
    t1->SetBranchAddress("neutrino1_pz_sol", &neutrino1_pz_sol);
    t1->SetBranchAddress("neutrino1_E_sol" , &neutrino1_E_sol );
    t1->SetBranchAddress("neutrino2_px_sol", &neutrino2_px_sol);
    t1->SetBranchAddress("neutrino2_py_sol", &neutrino2_py_sol);
    t1->SetBranchAddress("neutrino2_pz_sol", &neutrino2_pz_sol);
    t1->SetBranchAddress("neutrino2_E_sol" , &neutrino2_E_sol );
    t1->SetBranchAddress("WStarmass" , &WStarmass_sol);

    allEntries = t1->GetEntries();
    Long64_t entry;

    TH1D* mWWs_true_hist;
    if(generatedInfo) {mWWs_true_hist = new TH1D("mWWs_true_hist", "", 100, 0., 0.);}
    TH1D* mHH_hist = new TH1D("mHH_hist", "", 500, 0., 0.);
    TH1D* mJJ_hist             = new TH1D("mJJ_hist", "", 1000, 0., 0.);
    TH3D* mWWs_p1x_p1y_hist      = new TH3D("mWWs_p1x_p1y_hist", "", 500, -250., 250.,  500, -250., 250., 400, 0., 400.);
    TH2D* mWWs_mWs_hist = new TH2D("mWWs_mWs_hist", "", 400, 0., 400., 80, 0., 80.);

    TLorentzVector jet1, jet2, lepton1, lepton2, neutrino1_true , neutrino2_true, neutrino1_sol, neutrino2_sol;
    /* Event Loop */
    for(entry = 0; entry < allEntries; entry++)
    {
        t1->GetEntry(entry);

        jet1.SetPxPyPzE(jet1_px, jet1_py, jet1_pz, Energy(jet1_px, jet1_py, jet1_pz, MASSBQUARK));
        jet2.SetPxPyPzE(jet2_px, jet2_py, jet2_pz, Energy(jet2_px, jet2_py, jet2_pz, MASSBQUARK));
        lepton1.SetPxPyPzE(lepton1_px, lepton1_py, lepton1_pz, Energy(lepton1_px, lepton1_py, lepton1_pz, 0.0));
        lepton2.SetPxPyPzE(lepton2_px, lepton2_py, lepton2_pz, Energy(lepton2_px, lepton2_py, lepton2_pz, 0.0));

        if(generatedInfo)
        {
            neutrino1_true.SetPxPyPzE(neutrino1_px, neutrino1_py, neutrino1_pz, Energy(neutrino1_px, neutrino1_py, neutrino1_pz, 0.0));
            neutrino2_true.SetPxPyPzE(neutrino2_px, neutrino2_py, neutrino2_pz, Energy(neutrino2_px, neutrino2_py, neutrino2_pz, 0.0));
        }

        neutrino1_sol.SetPxPyPzE(neutrino1_px_sol, neutrino1_py_sol, neutrino1_pz_sol, Energy(neutrino1_px_sol, neutrino1_py_sol, neutrino1_pz_sol, 0.0));
        neutrino2_sol.SetPxPyPzE(neutrino2_px_sol, neutrino2_py_sol, neutrino2_pz_sol, Energy(neutrino2_px_sol, neutrino2_py_sol, neutrino2_pz_sol, 0.0));

        Double_t mWWs_mass_true;
        if(generatedInfo)
        {
            mWWs_mass_true = (lepton1 + lepton2 + neutrino1_true + neutrino2_true).M();
            mWWs_true_hist->Fill(mWWs_mass_true);
        }
        Double_t mWWs_mass_sol = (lepton1 + lepton2 + neutrino1_sol + neutrino2_sol).M();
        Double_t mJJ = (jet1 + jet2).M();
        Double_t mHH = (jet1 + jet2 + lepton1 + lepton2 + neutrino1_sol + neutrino2_sol).M();
        
        mWWs_mWs_hist->Fill(mWWs_mass_sol, WStarmass_sol);
        mWWs_p1x_p1y_hist->Fill(neutrino1_sol.Px(), neutrino1_sol.Py(), mWWs_mass_sol);
        mHH_hist->Fill(mHH);
        mJJ_hist->Fill(mJJ);

    }/* Event Loop */

    char buffer1[20];
    Int_t true_bin_width;
    if(generatedInfo)
    {
        true_bin_width = mWWs_true_hist->GetBinWidth(1);
        sprintf(buffer1, "Events/%dGeV", true_bin_width);
        mWWs_true_hist->GetXaxis()->SetTitle("M(ll#nu#nu)");
        mWWs_true_hist->GetYaxis()->SetTitle(buffer1);
    }

    char buffer2[20];
    TH1D* mWWs_sol_hist = mWWs_mWs_hist->ProjectionX();
    Int_t sol_bin_width = mWWs_sol_hist->GetBinWidth(1);
    sprintf(buffer2, "Events/%dGeV", sol_bin_width);
    mWWs_sol_hist->GetXaxis()->SetTitle("M(ll#nu#nu)");
    mWWs_sol_hist->GetYaxis()->SetTitle(buffer2);

    mWWs_p1x_p1y_hist->GetXaxis()->SetTitle("p1x");
    mWWs_p1x_p1y_hist->GetYaxis()->SetTitle("p1y");
    mWWs_p1x_p1y_hist->GetZaxis()->SetTitle("M(ll#nu#nu)");

    TH1D* p1x_hist = mWWs_p1x_p1y_hist->ProjectionX();
    TH1D* p1y_hist = mWWs_p1x_p1y_hist->ProjectionY();
    TH1D* mWWs_hist  = mWWs_p1x_p1y_hist->ProjectionZ();

    TH1* p1x_p1y_hist = mWWs_p1x_p1y_hist->Project3D("xy");
    TH1* p1y_mWWs_hist = mWWs_p1x_p1y_hist->Project3D("yz");
    TH1* p1x_mWWs_hist  = mWWs_p1x_p1y_hist->Project3D("xz");
    
    TH1D* mWs_sol_hist = mWWs_mWs_hist->ProjectionY();
    mWs_sol_hist->GetXaxis()->SetTitle("M(W^{*})");

    mJJ_hist->GetXaxis()->SetTitle("M(jj)");

    // char eos_path_hist[200];
    // sprintf(eos_path_hist, "/eos/user/a/alexandg/public/EOS.diHiggs/rootFiles/%s/histos/", physics_process);
    // char histo_file[200];
    // sprintf(histo_file, "histos_%s.root", physics_process);
    // strcat(eos_path_hist, histo_file);

    char eos_path_hist[200];
    sprintf(eos_path_hist, "/eos/user/a/alexandg/public/EOS.diHiggs/rootFiles/%s/histos/histos_%s.root", physics_process, physics_process);

    TFile* f_out = new TFile(eos_path_hist, "recreate");
    if(generatedInfo) {mWWs_true_hist->Write();}
    mWWs_sol_hist->Write();
    mWWs_p1x_p1y_hist->Write();
    p1x_hist->Write();
    p1y_hist->Write();
    mWWs_hist->Write();
    p1x_p1y_hist->Write();
    p1y_mWWs_hist->Write();
    p1x_mWWs_hist->Write();
    mWs_sol_hist->Write();
    mJJ_hist->Write();
    mHH_hist->Write();
    mWWs_mWs_hist->Write();

    gROOT->ProcessLine(".q");
}
