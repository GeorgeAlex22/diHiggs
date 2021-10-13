#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "TH2D.h"
#include "TH2.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h" 
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"


#include "headerFor2DSearch.h"
#include "Particle.h"
#include "LHAPDF/LHAPDF.h"

using namespace std;
using namespace TMath;

/* Input Flags */

const bool input_ToyMC = false;
const bool input_Delphes = !input_ToyMC;

/* ***************************************************** */
LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT10", 0);
/* ***************************************************** */

const bool printInfo = false;
const bool randomSolution = false;
const bool printLoopInfo = false;
const double Q = 100.; //GeV -- Electroweak Scale (needed to calculate the PDFs)
// const double ECM = 13000.; // GeV 




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

double sqrt_sum_sq(double* components)
{
    int size = sizeof(components)/sizeof(components[0]);
    double sum2 = 0.;

    for(int i = 0; i < size; i++)
        sum2 += (components[i],2);
    
    return pow(sum2, 0.5);
}
//--------------------------------------------------------//

double PDFWeight(LHAPDF::PDF* pdf, double x1, double x2, double Q)
{

    double g1_g, g1_d, g1_db, g1_u, g1_ub;
    double g2_g, g2_d, g2_db, g2_u, g2_ub;
    double pdfWeight;

    if (x1!=-1 && x2!=-1)
    {
        g1_g   = pdf->xfxQ( 0, x1, Q)/x1;
        g1_d   = pdf->xfxQ( 1, x1, Q)/x1;
        g1_db  = pdf->xfxQ(-1, x1, Q)/x1;
        g1_u   = pdf->xfxQ( 2, x1, Q)/x1;
        g1_ub  = pdf->xfxQ(-2, x1, Q)/x1;

        g2_g   = pdf->xfxQ( 0, x2, Q)/x2;
        g2_d   = pdf->xfxQ( 1, x2, Q)/x2;
        g2_db  = pdf->xfxQ(-1, x2, Q)/x2;
        g2_u   = pdf->xfxQ( 2, x2, Q)/x2;
        g2_ub  = pdf->xfxQ(-2, x2, Q)/x2;

        //int Solvability = 1;

        pdfWeight= g1_g * g2_g + g1_d * g2_db + g1_db * g2_d + g1_u * g2_ub + g1_ub * g2_u;

        return pdfWeight;
    }
}

//--------------------------------------------------------//

void solveW(double W_m, double muon1_px, double muon1_py, double muon1_pz, double MET_x, double MET_y, std::vector<double>* MET_z) 
{
    (*MET_z).clear();

    // double muon1_m=0.106;
    double muon1_m=0.;
    double muon1_e   = sqrt( muon1_m*muon1_m+ muon1_px*muon1_px+muon1_py*muon1_py+muon1_pz*muon1_pz);


    double alpha = W_m*W_m - muon1_m * muon1_m + 2 * muon1_px * MET_x + 2 * muon1_py * MET_y ;
    double beta  = 4 * muon1_e * muon1_e * (MET_x * MET_x + MET_y * MET_y);

    double   A = 4* (muon1_e * muon1_e - muon1_pz * muon1_pz);
    double   B = (-4) * alpha * muon1_pz;
    double   C = beta - alpha* alpha;

    double   D = B*B-4*A*C;

    if (D>0)
    {
        double   pnz1 = ((-1)*B+sqrt(D)) / (2*A);
        double   pnz2 = ((-1)*B-sqrt(D)) / (2*A);

        (*MET_z).push_back(pnz1);
        (*MET_z).push_back(pnz2);
        if(printInfo)
        {
        printf("\n W boson mass %lf ",W_m);
        printf("\n pnz1 = %f ", pnz1);
        printf("\n pnz2 = %f \n\n", pnz2);      
        }
    }
    else if(printInfo) printf("\n no solutions %lf \n",D);
}

//--------------------------------------------------------//

void resolve_combinatorics_neutrinos(std::vector<TLorentzVector>* neutrinos_p4, TLorentzVector jj_p4, TLorentzVector ll_p4, TLorentzVector neutrino1_p4_1, TLorentzVector neutrino1_p4_2, TLorentzVector neutrino2_p4_1, TLorentzVector neutrino2_p4_2)
{
    // std::vector<TLorentzVector> neutrinos_p4;
    (*neutrinos_p4).clear();

    TLorentzVector nn_solutions[4][2] = {{neutrino1_p4_1, neutrino2_p4_1}, {neutrino1_p4_1, neutrino2_p4_2}, {neutrino1_p4_2, neutrino2_p4_1}, {neutrino1_p4_2, neutrino2_p4_2}};

    TLorentzVector arr_WWs_p4[4]; 
    TLorentzVector jjWWs_p4;

    int result_index = -1;

    double x1 = -1, x2 = -1;
    double pdfWeight = -1;
    double tmppdfWeight = -1;

    bool isUnphysical = true;

    for(int i = 0; i < 4; i++)
    {
        arr_WWs_p4[i] = ll_p4 + nn_solutions[i][0] + nn_solutions[i][1];

        jjWWs_p4 = jj_p4 + arr_WWs_p4[i];

        isUnphysical = ( TMath::IsNaN(jjWWs_p4.Px()) || !TMath::Finite(jjWWs_p4.Px()))
                    && ( TMath::IsNaN(jjWWs_p4.Py()) || !TMath::Finite(jjWWs_p4.Py())) 
                    && ( TMath::IsNaN(jjWWs_p4.Pz()) || !TMath::Finite(jjWWs_p4.Pz())) 
                    && ( TMath::IsNaN(jjWWs_p4.E() ) || !TMath::Finite(jjWWs_p4.E()) )
                    && (jjWWs_p4.M() > ECM);

        if(!isUnphysical)
        {
            x1 = (jjWWs_p4.E() + jjWWs_p4.Pz())/ECM;
            if(x1 > 0.999) x1 = 0.999;
            x2 = (jjWWs_p4.E() - jjWWs_p4.Pz())/ECM;
            if(x2 > 0.999) x2 = 0.999;


            // Select combination with highest pdf score.
            tmppdfWeight = PDFWeight(pdf, x1, x2, Q);
            if(pdfWeight < tmppdfWeight)
            {
                pdfWeight = tmppdfWeight;
                result_index = i;
            }


        }// if(!isUnphysical)
    }// Solution Loop
    if(result_index != -1)
    {
        (*neutrinos_p4).push_back(nn_solutions[result_index][0]);
        (*neutrinos_p4).push_back(nn_solutions[result_index][1]);
    }

    // return neutrinos_p4;
}


//--------------------------------------------------------//

void solveWW(std::vector<TLorentzVector>* neutrinos_p4_sol,double Wmass, double WStarmass, double* l1, double* l2, double n1x, double n1y, double n2x, double n2y, double* j1, double* j2)
{
    std::vector<double> n1z;
    std::vector<double> n2z;

    // Jet 4-Momenta
    TLorentzVector jet1_p4;
    jet1_p4.SetPxPyPzE(j1[1], j1[2], j1[3], j1[0]);
    TLorentzVector jet2_p4;
    jet2_p4.SetPxPyPzE(j2[1], j2[2], j2[3], j2[0]);
    TLorentzVector jj_p4 = jet1_p4 + jet2_p4;
    double mJJ = jj_p4.M();

    // Lepton 4-Momenta
    TLorentzVector lepton1_p4;
    lepton1_p4.SetPxPyPzE(l1[1], l1[2], l1[3], l1[0]);
    TLorentzVector lepton2_p4;
    lepton2_p4.SetPxPyPzE(l2[1], l2[2], l2[3], l2[0]); 
    TLorentzVector ll_p4 = lepton1_p4 + lepton2_p4;
    double mll = ll_p4.M(); // in case we need it    

    solveW(Wmass    , l1[1], l1[2], l1[3], n1x, n1y, &n1z);
    solveW(WStarmass, l2[1], l2[2], l2[3], n2x, n2y, &n2z);

    if(n1z.size()*n2z.size() > 0)
    {

        // if(randomSolution)
        // {
        //     TRandom* rndm = new TRandom(); // For now choose randomly a solution combination
        //     int index1 = rndm->Integer(n1z.size());
        //     TLorentzVector neutrino1_p4_solRndm; 
        //     neutrino1_p4_solRndm.SetPxPyPzE(n1x    , n1y    , n1z[index1], Energy(n1x    , n1y    , n1z[index1], 0.0));
        //     int index2 = rndm->Integer(n2z.size());
        //     TLorentzVector neutrino2_p4_solRndm;
        //     neutrino2_p4_solRndm.SetPxPyPzE(n2x, n2y, n2z[index2], Energy(n2x, n2y, n2z[index2], 0.0));

        //     double mWWsRndm = (ll_p4 + neutrino1_p4_solRndm + neutrino2_p4_solRndm).M();

        // }// if(randomSolution)

        // else
        // {
            TLorentzVector neutrino1_p4_1;
            neutrino1_p4_1.SetPxPyPzE(n1x    , n1y    , n1z[0], Energy(n1x    , n1y    , n1z[0], 0.0));
            TLorentzVector neutrino1_p4_2;
            neutrino1_p4_2.SetPxPyPzE(n1x    , n1y    , n1z[1], Energy(n1x    , n1y    , n1z[1], 0.0));
            TLorentzVector neutrino2_p4_1;
            neutrino2_p4_1.SetPxPyPzE(n2x    , n2y    , n2z[0], Energy(n2x    , n2y    , n2z[0], 0.0));
            TLorentzVector neutrino2_p4_2;
            neutrino2_p4_2.SetPxPyPzE(n2x    , n2y    , n2z[1], Energy(n2x    , n2y    , n2z[1], 0.0));

            std::vector<TLorentzVector> neutrinos_p4 ;
            resolve_combinatorics_neutrinos(&neutrinos_p4, jj_p4, ll_p4, neutrino1_p4_1, neutrino1_p4_2, neutrino2_p4_1, neutrino2_p4_2);
            
            // Check if there are any (physical) solutions for these (WStarmass, p1x, p1y)
            // Physicallity is already checked in resolve_combinatorics_neutrinos(...)
            if(printLoopInfo) {cout << "Size of neutrinos_p4 in solveWW: " << neutrinos_p4.size() << endl;}
            if(neutrinos_p4.size() > 0)
            {
                (*neutrinos_p4_sol) = neutrinos_p4;                
            }                                                      
        // }// if(randomSolution)                                                                                                                                                                                                                                                                                                   
    } // if(n1z.size()*n2z.size() > 0)    
    else if(printLoopInfo) {cout << "... solveW or solveW* did not find solution." << endl;}
}

//--------------------------------------------------------//

void solve_system(std::vector<TLorentzVector>* neutrinos_p4_sol, double Wmass, double WStarmass, double* l1, double* l2, double n1x, double n1y, double n2x, double n2y, double* j1, double* j2)
{
    double x1 = -1, x2 = -1;
    double tmppdfWeight = -1;
    double pdfWeight = -1;
    int result_index = -1;

    std::vector<TLorentzVector> neutrinos_p4_1122;
    solveWW(&neutrinos_p4_1122, Wmass, WStarmass, l1, l2, n1x, n1y, n2x, n2y, j1, j2);

    std::vector<TLorentzVector> neutrinos_p4_1221;
    solveWW(&neutrinos_p4_1221, Wmass, WStarmass, l2, l1, n1x, n1y, n2x, n2y, j1, j2);

    // Jet 4-Momenta
    TLorentzVector jet1_p4;
    jet1_p4.SetPxPyPzE(j1[1], j1[2], j1[3], j1[0]);
    TLorentzVector jet2_p4;
    jet2_p4.SetPxPyPzE(j2[1], j2[2], j2[3], j2[0]);
    TLorentzVector jj_p4 = jet1_p4 + jet2_p4;
    double mJJ = jj_p4.M();

    // Lepton 4-Momenta
    TLorentzVector lepton1_p4;
    lepton1_p4.SetPxPyPzE(l1[1], l1[2], l1[3], l1[0]);
    TLorentzVector lepton2_p4;
    lepton2_p4.SetPxPyPzE(l2[1], l2[2], l2[3], l2[0]); 
    TLorentzVector ll_p4 = lepton1_p4 + lepton2_p4;
    double mll = ll_p4.M(); // in case we need it        

    std::vector<TLorentzVector> neutrinos_p4[2] =  {neutrinos_p4_1122, neutrinos_p4_1221};
    for(int i = 0; i < 2; i++)
    {   
        if(printLoopInfo) {cout << "Size of neutrinos_p4[" << i << "] in solve_system: " << neutrinos_p4[i].size() << endl;}
        if(neutrinos_p4[i].size() > 0)
        {
            TLorentzVector nn_p4 = neutrinos_p4[i][0] + neutrinos_p4[i][1];

            TLorentzVector WWs_p4 = ll_p4 + nn_p4;
            double mWWs = WWs_p4.M();

            TLorentzVector event_p4 = jj_p4 + ll_p4 + nn_p4;

            x1 = (event_p4.E() + event_p4.Pz())/ECM;
            if(x1 > 0.999) x1 = 0.999;
            x2 = (event_p4.E() - event_p4.Pz())/ECM;
            if(x2 > 0.999) x2 = 0.999;

            tmppdfWeight = PDFWeight(pdf, x1, x2, Q);
            if(pdfWeight < tmppdfWeight)
            {
                pdfWeight = tmppdfWeight;
                result_index = i;

            }           
        }                        
    }
    // Set solution
    if(result_index != -1)
    {
        (*neutrinos_p4_sol) = neutrinos_p4[result_index];
    }    

}


//--------------------------------------------------------//



int main() {

    double l1[4], l2[4], j1[4], j2[4];
    double MET[2], n1[4], n2[4];

    double  wp_mm,tp_mm;
    double  wm_mm,tm_mm;

    TFile *f;
    f = new TFile("data/pdf_input.root");

    /* Read Input Tree */

    TTree* t1 = (TTree *) f->Get("t1");

    Double_t MET_x, MET_y;

    Double_t jet1_px, jet1_py, jet1_pz;                 
    Double_t jet2_px, jet2_py, jet2_pz;                 

    Double_t lepton1_px, lepton1_py, lepton1_pz;
    Int_t lepton1_charge;

    Double_t lepton2_px, lepton2_py, lepton2_pz;
    Int_t lepton2_charge;

    Double_t neutrino1_px, neutrino1_py, neutrino1_pz;
    Double_t neutrino2_px, neutrino2_py, neutrino2_pz;

    t1->SetBranchAddress("MET_x", &MET_x);
    t1->SetBranchAddress("MET_y", &MET_y); 
    t1->SetBranchAddress("jet1_px", &jet1_px);
    t1->SetBranchAddress("jet1_py", &jet1_py);
    t1->SetBranchAddress("jet1_pz", &jet1_pz);
    t1->SetBranchAddress("jet2_px", &jet2_px);
    t1->SetBranchAddress("jet2_py", &jet2_py);
    t1->SetBranchAddress("jet2_pz", &jet2_pz);
    t1->SetBranchAddress("lepton1_px", &lepton1_px);
    t1->SetBranchAddress("lepton1_py", &lepton1_py);
    t1->SetBranchAddress("lepton1_pz", &lepton1_pz);
    t1->SetBranchAddress("lepton2_px", &lepton2_px);
    t1->SetBranchAddress("lepton2_py", &lepton2_py);
    t1->SetBranchAddress("lepton2_pz", &lepton2_pz);
    t1->SetBranchAddress("neutrino1_px", &neutrino1_px);
    t1->SetBranchAddress("neutrino1_py", &neutrino1_py);
    t1->SetBranchAddress("neutrino1_pz", &neutrino1_pz);
    t1->SetBranchAddress("neutrino2_py", &neutrino2_py);
    t1->SetBranchAddress("neutrino2_px", &neutrino2_px);
    t1->SetBranchAddress("neutrino2_pz", &neutrino2_pz);

    Long64_t allEntries = t1->GetEntries();

    /* If non-empty input, create output file */

    TFile* treeResults;
    TTree *t2 ;
    /* suffix "0" refers to the variables of the output tree that would otherwise have the same name as the ones on the input tree. */
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

    if(allEntries != 0)
    {    
        treeResults = new TFile("data/pdf_output.root", "RECREATE");       
        t2  = new TTree("t1","RECREATE");

        t2->Branch("MET_x", &MET_x0, "MET_x/D");
        t2->Branch("MET_y", &MET_y0, "MET_y/D");
        t2->Branch("jet1_px", &jet1_px0, "jet1_px/D");
        t2->Branch("jet1_py", &jet1_py0, "jet1_py/D");
        t2->Branch("jet1_pz", &jet1_pz0, "jet1_pz/D");
        t2->Branch("jet1_E" , &jet1_E0 , "jet1_E/D" );
        t2->Branch("jet2_px", &jet2_px0, "jet2_px/D");
        t2->Branch("jet2_py", &jet2_py0, "jet2_py/D");
        t2->Branch("jet2_pz", &jet2_pz0, "jet2_pz/D");
        t2->Branch("jet2_E" , &jet2_E0 , "jet2_E/D" );
        t2->Branch("lepton1_px", &lepton1_px0, "lepton1_px/D");
        t2->Branch("lepton1_py", &lepton1_py0, "lepton1_py/D");
        t2->Branch("lepton1_pz", &lepton1_pz0, "lepton1_pz/D");
        t2->Branch("lepton1_E" , &lepton1_E0 , "lepton1_E/D" );
        t2->Branch("lepton2_px", &lepton2_px0, "lepton2_px/D");
        t2->Branch("lepton2_py", &lepton2_py0, "lepton2_py/D");
        t2->Branch("lepton2_pz", &lepton2_pz0, "lepton2_pz/D");
        t2->Branch("lepton2_E" , &lepton2_E0 , "lepton2_E/D" );
        t2->Branch("neutrino1_px", &neutrino1_px0, "neutrino1_px/D");
        t2->Branch("neutrino1_py", &neutrino1_py0, "neutrino1_py/D");
        t2->Branch("neutrino1_pz", &neutrino1_pz0, "neutrino1_pz/D");
        t2->Branch("neutrino1_E" , &neutrino1_E0 , "neutrino1_E/D" );
        t2->Branch("neutrino2_px", &neutrino2_px0, "neutrino2_px/D");
        t2->Branch("neutrino2_py", &neutrino2_py0, "neutrino2_py/D");
        t2->Branch("neutrino2_pz", &neutrino2_pz0, "neutrino2_pz/D");
        t2->Branch("neutrino2_E" , &neutrino2_E0 , "neutrino2_E/D" );        
        t2->Branch("neutrino1_px_sol", &neutrino1_px_sol, "neutrino1_px_sol/D");
        t2->Branch("neutrino1_py_sol", &neutrino1_py_sol, "neutrino1_py_sol/D");
        t2->Branch("neutrino1_pz_sol", &neutrino1_pz_sol, "neutrino1_pz_sol/D");
        t2->Branch("neutrino1_E_sol" , &neutrino1_E_sol , "neutrino1_E_sol/D" );
        t2->Branch("neutrino2_px_sol", &neutrino2_px_sol, "neutrino2_px_sol/D");
        t2->Branch("neutrino2_py_sol", &neutrino2_py_sol, "neutrino2_py_sol/D");
        t2->Branch("neutrino2_pz_sol", &neutrino2_pz_sol, "neutrino2_pz_sol/D");
        t2->Branch("neutrino2_E_sol" , &neutrino2_E_sol , "neutrino2_E_sol/D" );            
        t2->Branch("WStarmass" , &result_WStarmass , "WStarmass/D" );
    }

    int entry;
    int solvable_events = 0;

    bool event_isTotallySolvable = false;

    // allEntries = 3; // for debugging
    // Event Loop
    for(entry = 0; entry < allEntries; entry++)
    {

        cout << "------------ Event " << entry << "------------" << endl;

        t1->GetEntry(entry);

        MET[0] = MET_x;
        MET[1] = MET_y;

        j1[1] = jet1_px;
        j1[2] = jet1_py;
        j1[3] = jet1_pz;
        j1[0] = Energy(jet1_px, jet1_py, jet1_pz, MASSBQUARK);

        j2[1] = jet2_px;
        j2[2] = jet2_py;
        j2[3] = jet2_pz;
        j2[0] = Energy(jet2_px, jet2_py, jet2_pz, MASSBQUARK);

        l1[1] = lepton1_px;
        l1[2] = lepton1_py;
        l1[3] = lepton1_pz;
        l1[0] = Energy(lepton1_px, lepton1_py, lepton1_pz, 0.);

        l2[1] = lepton2_px;
        l2[2] = lepton2_py;
        l2[3] = lepton2_pz;
        l2[0] = Energy(lepton2_px, lepton2_py, lepton2_pz, 0.);

        n1[1] = neutrino1_px;
        n1[2] = neutrino1_py;
        n1[3] = neutrino1_pz;
        n1[0] = Energy(neutrino1_px, neutrino1_py, neutrino1_pz, 0.);

        n2[1] = neutrino2_px;
        n2[2] = neutrino2_py;
        n2[3] = neutrino2_pz;
        n2[0] = Energy(neutrino2_px, neutrino2_py, neutrino2_pz, 0.);

        double Wmass     = 80.; // GeV
        double WStarmass = 45.; // GeV

        double x1 = -1, x2 = -1;
        double pdfWeight = -1;
        double tmppdfWeight = -1;

        // Jet 4-Momenta
        TLorentzVector jet1_p4;
        jet1_p4.SetPxPyPzE(j1[1], j1[2], j1[3], j1[0]);
        TLorentzVector jet2_p4;
        jet2_p4.SetPxPyPzE(j2[1], j2[2], j2[3], j2[0]);
        TLorentzVector jj_p4 = jet1_p4 + jet2_p4;
        double mJJ = jj_p4.M();

        // Lepton 4-Momenta
        TLorentzVector lepton1_p4;
        lepton1_p4.SetPxPyPzE(l1[1], l1[2], l1[3], l1[0]);
        TLorentzVector lepton2_p4;
        lepton2_p4.SetPxPyPzE(l2[1], l2[2], l2[3], l2[0]); 
        TLorentzVector ll_p4 = lepton1_p4 + lepton2_p4;
        double mll = ll_p4.M(); // in case we need it

        // Neutrino 4-Momenta Solution for this event (initialization)
        TLorentzVector neutrino1_p4;
        neutrino1_p4.SetPxPyPzE(0., 0., 0., 777.); // for  debugging
        TLorentzVector neutrino2_p4;
        neutrino2_p4.SetPxPyPzE(0., 0., 0., 777.); // for  debugging

        // double result_WStarmass = -1;

        double WStarmass0 = 0.;
        double WStarmassMAX = 80.;
        double WStarmassSTEP = 1.;
        double momentumSTEP = 1.;
        double momentumMAX = 250.;
        for(WStarmass = WStarmass0; WStarmass < WStarmassMAX; WStarmass += WStarmassSTEP)
        {
            // cout << WStarmass << endl;
            double n1x, n1y;
            // initialize
            n1x = n1[0];
            n1y = n1[1];
            for(n1x = -momentumMAX; n1x < momentumMAX+1; n1x++)
            {
                for(n1y = -momentumMAX; n1y < momentumMAX+1; n1y++)
                {
                    if(printLoopInfo) {cout << "M(W*) = " << WStarmass << "\tp1x = " << n1x << "\tp1y = " << n1y << endl;}
                    double n2x = MET_x - n1x;
                    double n2y = MET_y - n1y;

                    std::vector<TLorentzVector> neutrinos_p4;
                    neutrinos_p4.clear();                    
                    solve_system(&neutrinos_p4, Wmass, WStarmass, l1, l2, n1x, n1y, n2x, n2y, j1, j2);

                    if(printLoopInfo) {cout << "\tSolution " << neutrinos_p4.size() << endl;}

                    if(neutrinos_p4.size() > 0)
                    {
                        TLorentzVector tmp_neutrino1_p4 = neutrinos_p4[0];
                        TLorentzVector tmp_neutrino2_p4 = neutrinos_p4[1];

                        TLorentzVector tmp_nn_p4 = tmp_neutrino1_p4 + tmp_neutrino2_p4;
                        TLorentzVector tmp_jjWWs = jj_p4 + ll_p4 + tmp_nn_p4;                    

                        x1 = (tmp_jjWWs.E() + tmp_jjWWs.Pz())/ECM;
                        if(x1 > 0.999) x1 = 0.999;
                        x2 = (tmp_jjWWs.E() - tmp_jjWWs.Pz())/ECM;
                        if(x2 > 0.999) x2 = 0.999;

                        tmppdfWeight = PDFWeight(pdf, x1, x2, Q);

                        // cout << "\t\tpdfWeight: " << pdfWeight << "\ttmppdfWeight: " << tmppdfWeight << endl;

                        if(pdfWeight < tmppdfWeight)
                        {
                            // cout << "\t\tpdfWeight: " << pdfWeight << "\ttmppdfWeight: " << tmppdfWeight << endl;
                            pdfWeight = tmppdfWeight;
                            result_WStarmass = WStarmass;
                            neutrino1_p4 = tmp_neutrino1_p4;
                            neutrino2_p4 = tmp_neutrino2_p4;
                            event_isTotallySolvable = true;
                        }

                    }

                } // Momentum Y Loop                
            } // Momentum X Loop
        } // W* Mass Loop

        if(event_isTotallySolvable)
        {
            solvable_events++ ;
            // cout << "pdfWeight: " << pdfWeight << endl;

            TLorentzVector nn_p4 = neutrino1_p4 + neutrino2_p4;
            TLorentzVector WWs_p4 = ll_p4 + nn_p4;
            TLorentzVector jjWWs_p4 = jj_p4 + WWs_p4;
            double mWWs = WWs_p4.M();
            double mjjWWs = jjWWs_p4.M();

            // cout << " True --  M(W*): N/A" <<        ""        << "\t p1x: " <<          n1[1]          << "\t p1y: " <<          n1[2]          << endl;
            // cout << "Solved--  M(W*):    " << result_WStarmass << "\t p1x: " <<    neutrino1_p4.Px()    << "\t p1y: " <<    neutrino1_p4.Py()    << endl;
            // cout << " Diff -- DM(W*): N/A" <<        ""        << "\tDp1x: " << n1[1]-neutrino1_p4.Px() << "\tDp1y: " << n1[2]-neutrino1_p4.Py() << endl;
            // cout << "mWWs = " << mWWs << endl;

            MET_x0 = MET_x, MET_y0 = MET_y;

            jet1_px0 = jet1_p4.Px(), jet1_py0 = jet1_p4.Py(), jet1_pz0 = jet1_p4.Pz(), jet1_E0 = jet1_p4.E();
            jet2_px0 = jet2_p4.Px(), jet2_py0 = jet2_p4.Py(), jet2_pz0 = jet2_p4.Pz(), jet2_E0 = jet2_p4.E();   

            lepton1_px0 = lepton1_p4.Px(), lepton1_py0 = lepton1_p4.Py(), lepton1_pz0 = lepton1_p4.Pz(), lepton1_E0 = lepton1_p4.E();
            //lepton1_charge0 = ;

            lepton2_px0 = lepton2_p4.Px(), lepton2_py0 = lepton2_p4.Py(), lepton2_pz0 = lepton2_p4.Pz(), lepton2_E0 = lepton2_p4.E();
            //lepton2_charge0 = ;

            neutrino1_px0 = neutrino1_px, neutrino1_py0 = neutrino1_py, neutrino1_pz0 = neutrino1_pz, neutrino1_E0 = Energy(neutrino1_px, neutrino1_py, neutrino1_pz, 0.);
            neutrino2_px0 = neutrino2_px, neutrino2_py0 = neutrino2_py, neutrino2_pz0 = neutrino2_pz, neutrino2_E0 = Energy(neutrino2_px, neutrino2_py, neutrino2_pz, 0.);
            neutrino1_px_sol = neutrino1_p4.Px(), neutrino1_py_sol = neutrino1_p4.Py(), neutrino1_pz_sol = neutrino1_p4.Pz(), neutrino1_E_sol = neutrino1_p4.E();
            neutrino2_px_sol = neutrino2_p4.Px(), neutrino2_py_sol = neutrino2_p4.Py(), neutrino2_pz_sol = neutrino2_p4.Pz(), neutrino2_E_sol = neutrino2_p4.E();

            t2->Fill();
            
        }

       
        


    } //Event Loop

    if(allEntries != 0)
    {
        t2->Write();
        treeResults->Close();
        cout << "** Input contains " << allEntries << " events from which " << solvable_events << " were solvable" <<  endl;

    }

    // f->Close();


    return 0;

}


