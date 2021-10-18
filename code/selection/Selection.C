/* 
This macro takes a delphes_output.root file and performs a selection based
on the diHiggs signal. It then outputs a selection_output.root file containing a TTree
with no Delphes dependence.

The selection_input.root files I am using are inside 
/eos/cms/store/group/phys_b2g/giorgos/pheno/rootFiles/$(PHYSICS_PROCCESS)/trees/
where $(PHYSICS_PROCCESS) can be pptohh, WW, WZ, ZZ, ttbar1 etc.. (just ls that directory up to /eos/cms/......./rootFiles/ and pich a subdirectory)

and have the following format: run_%c.lhe.hep.root
where %c is 01, 02, 03, ..., 09, 10, 11, 12, ..., 500

version: Delphes-3.4.0

cp /eos/cms/store/group/phys_b2g/giorgos/pheno/rootFiles/pptohh/trees/run_01.lhe.hep.root ./selection_input.root
root -l -b Selection.C
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif
#include "TTree.h"
//#include <vector>

//------------------------------------------------------------------------------

class Lepton // public SortableObject
{
public:
  Float_t PT; // lepton transverse momentum
  Float_t Eta; // lepton pseudorapidity
  Float_t Phi; // lepton azimuthal angle

  Float_t T; // particle arrival time of flight

  Int_t Charge; // lepton charge

  Float_t EhadOverEem; // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter

  TRef Particle; // reference to generated particle

  Float_t IsolationVar; // isolation variable
  Float_t IsolationVarRhoCorr; // isolation variable
  Float_t SumPtCharged; // isolation variable
  Float_t SumPtNeutral; // isolation variable
  Float_t SumPtChargedPU; // isolation variable
  Float_t SumPt; // isolation variable

  Float_t D0; // track transverse impact parameter
  Float_t DZ; // track longitudinal impact parameter
  Float_t ErrorD0; // track transverse impact parameter error
  Float_t ErrorDZ; // track longitudinal impact parameter error

  Int_t LeptonID; // lepton PID

  //static CompBase *fgCompare; //!
  //const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const
  {
    TLorentzVector vec;
    vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
    return vec;
  }

  ClassDef(Lepton, 4)
};

//------------------------------------------------------------------------------

bool sortLeptonsOperator(const Lepton& l1, const Lepton& l2)                          { return l1.PT > l2.PT;}
bool sortGenParticleOperator(const GenParticle& p1, const GenParticle& p2)            { return p1.PT > p2.PT;}
bool sortJetsOperator(const Jet& j1, const Jet& j2)                                   { return j1.PT > j2.PT;}
bool sortLorentzOperator(const TLorentzVector& v1, const TLorentzVector& v2)          { return v1.Pt() > v2.Pt();}
bool sortElectrons(const Electron& e1, const Electron& e2)                            { return e1.PT > e2.PT;}
bool sortMuons(const Muon&     m1, const Muon&     m2)                                { return m1.PT > m2.PT;}


bool generatedInfo  = false;
bool printInfo      = !true;
bool debug          = !true;
bool signal_flag    = true;
int  selectedEvents = 0;
//------------------------------------------------------------------------------



void Selection()
{

  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");

  chain->Add("selection_input.root");

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

//////////////////////////////////////////////////////////////////  
// root tree variables, branches 
//////////////////////////////////////////////////////////////////

  // variables for filling /
  Double_t MET_x, MET_y;

  Double_t jet1_px, jet1_py, jet1_pz, jet1_btagging, jet1_HoE;                 Int_t jet1_charge, jet1_ncharged, jet1_nneutral;
  Double_t jet2_px, jet2_py, jet2_pz, jet2_btagging, jet2_HoE;                 Int_t jet2_charge, jet2_ncharged, jet2_nneutral;

  Double_t jet3_px, jet3_py, jet3_pz, jet3_btagging, jet3_HoE;                 Int_t jet3_charge, jet3_ncharged, jet3_nneutral;
  Double_t jet4_px, jet4_py, jet4_pz, jet4_btagging, jet4_HoE;                 Int_t jet4_charge, jet4_ncharged, jet4_nneutral;


  Double_t lepton1_px, lepton1_py, lepton1_pz, lepton1_HoE, lepton1_IsoVar, lepton1_IsoRho, lepton1_SumPtCharged, lepton1_SumPtNeutral, lepton1_SumPtChargedPU, lepton1_SumPt;
  Int_t lepton1_charge;

  Double_t lepton2_px, lepton2_py, lepton2_pz, lepton2_HoE, lepton2_IsoVar, lepton2_IsoRho, lepton2_SumPtCharged, lepton2_SumPtNeutral, lepton2_SumPtChargedPU, lepton2_SumPt;
  Int_t lepton2_charge;

  Double_t neutrino1_px, neutrino1_py, neutrino1_pz;
  Double_t neutrino2_px, neutrino2_py, neutrino2_pz;


  // tree declaration, branches /
  TFile* treeResults;
  treeResults = new TFile("selection_output.root", "RECREATE");

  TTree *t1 ;
  t1  = new TTree("t1","RECREATE");

  t1->Branch("MET_x",&MET_x,"MET_x/D");
  t1->Branch("MET_y",&MET_y,"MET_y/D");
  t1->Branch("jet1_px",&jet1_px,"jet1_px/D");
  t1->Branch("jet1_py",&jet1_py,"jet1_py/D");
  t1->Branch("jet1_pz",&jet1_pz,"jet1_pz/D");
  t1->Branch("jet2_px",&jet2_px,"jet2_px/D");
  t1->Branch("jet2_py",&jet2_py,"jet2_py/D");
  t1->Branch("jet2_pz",&jet2_pz,"jet2_pz/D");
  t1->Branch("lepton1_px", &lepton1_px, "lepton1_px/D");
  t1->Branch("lepton1_py", &lepton1_py, "lepton1_py/D");
  t1->Branch("lepton1_pz", &lepton1_pz, "lepton1_pz/D");
  t1->Branch("lepton2_px", &lepton2_px, "lepton2_px/D");
  t1->Branch("lepton2_py", &lepton2_py, "lepton2_py/D");
  t1->Branch("lepton2_pz", &lepton2_pz, "lepton2_pz/D");
  t1->Branch("neutrino1_px", &neutrino1_px, "neutrino1_px/D");
  t1->Branch("neutrino1_py", &neutrino1_py, "neutrino1_py/D");
  t1->Branch("neutrino1_pz", &neutrino1_pz, "neutrino1_pz/D");
  t1->Branch("neutrino2_px", &neutrino2_px, "neutrino2_px/D");
  t1->Branch("neutrino2_py", &neutrino2_py, "neutrino2_py/D");
  t1->Branch("neutrino2_pz", &neutrino2_pz, "neutrino2_pz/D");

 if (!generatedInfo)
 {
  t1->Branch("jet1_btagging",&jet1_btagging,"jet1_btagging/D");
  t1->Branch("jet1_charge",&jet1_charge,"jet1_charge/I");
  t1->Branch("jet1_HoE",&jet1_HoE,"jet1_HoE/D");
  t1->Branch("jet1_ncharged",&jet1_ncharged,"jet1_ncharged/I");
  t1->Branch("jet1_nneutral",&jet1_nneutral,"jet1_nneutral/I");
  
  t1->Branch("jet2_btagging",&jet2_btagging,"jet2_btagging/D");
  t1->Branch("jet2_charge",&jet2_charge,"jet2_charge/I");
  t1->Branch("jet2_HoE",&jet2_HoE,"jet2_HoE/D");
  t1->Branch("jet2_ncharged",&jet2_ncharged,"jet2_ncharged/I");
  t1->Branch("jet2_nneutral",&jet2_nneutral,"jet2_nneutral/I");
 
  t1->Branch("lepton1_charge", &lepton1_charge, "lepton1_charge/I");
  t1->Branch("lepton1_HoE", &lepton1_HoE, "lepton1_HoE/D"); 
  t1->Branch("lepton1_IsoVar", &lepton1_IsoVar, "lepton1_IsoVar/D");
  t1->Branch("lepton1_IsoRho", &lepton1_IsoRho, "lepton1_IsoRho/D");
  t1->Branch("lepton1_SumPtCharged", &lepton1_SumPtCharged, "lepton1_SumPtCharged/D");
  t1->Branch("lepton1_SumPtNeutral", &lepton1_SumPtNeutral, "lepton1_SumPtNeutral/D");
  t1->Branch("lepton1_SumPtChargedPU", &lepton1_SumPtChargedPU, "lepton1_SumPtChargedPU/D");
  t1->Branch("lepton1_SumPt", &lepton1_SumPt, "lepton1_SumPt/D");
  
  t1->Branch("lepton2_charge", &lepton2_charge, "lepton2_charge/I");
  t1->Branch("lepton2_HoE", &lepton2_HoE, "lepton2_HoE/D");
  t1->Branch("lepton2_IsoVar", &lepton2_IsoVar, "lepton2_IsoVar/D");
  t1->Branch("lepton2_IsoRho", &lepton2_IsoRho, "lepton2_IsoRho/D");
  t1->Branch("lepton2_SumPtCharged", &lepton2_SumPtCharged, "lepton2_SumPtCharged/D");
  t1->Branch("lepton2_SumPtNeutral", &lepton2_SumPtNeutral, "lepton2_SumPtNeutral/D");
  t1->Branch("lepton2_SumPtChargedPU", &lepton2_SumPtChargedPU, "lepton2_SumPtChargedPU/D");
  t1->Branch("lepton2_SumPt", &lepton2_SumPt, "lepton2_SumPt/D");

 }
/////////////////////////////////////////////////////////////////////  

  /* for reading input tree */
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  Long64_t allEntries = treeReader->GetEntries();

  //cout << "** Chain contains " << allEntries << " events" << endl;
  if (printInfo) cout << "** Chain contains " << allEntries << " events" << endl;

  TObject *object;
  Long64_t entry;

// ********************************* /
// Loop over all events              /
// ********************************* /

  // allEntries = 10; // only proccess a small number of events for debugging


  for(entry = 0; entry < allEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    int njets=0;
    int nelectrons=0;
    int nmuons=0;
    int nleptons=0;

    double   missingEnergyX, missingEnergyY;
    Jet      mostEnergeticJet,      secEnergeticJet; 
    Electron mostEnergeticElectron, secEnergeticElectron;
    Muon     mostEnergeticMuon,     secEnergeticMuon; 
    Lepton   mostEnergeticLepton,   secEnergeticLepton;
 
    vector<Lepton> Leptons; Leptons.clear();
    vector<Jet>    Jets;    Jets.clear();
 
    // generated values /
    double genMETx=0; 
    double genMETy=0;
    double genMEz =0;
    vector<GenParticle> genBquarks;   genBquarks.clear();
    vector<GenParticle> genLeptons;   genLeptons.clear();
    vector<GenParticle> genNeutrinos;   genNeutrinos.clear();

    GenParticle   mostEnergeticGenLepton;
    GenParticle   secEnergeticGenLepton; 
    GenParticle   mostEnergeticBquark;
    GenParticle   secEnergeticBquark;
    GenParticle   mostEnergeticGenNeutrino;
    GenParticle   secEnergeticGenNeutrino; 

    GenParticle *particle;
    Int_t i, j, pdgCode;
    
    // if (printInfo){ cout << endl; cout << "------------  event " << entry << " -------------- " << endl; cout << endl;}


//*************************** /
//     generated event        /
//*************************** /

    // if (printInfo) { cout << " genereated event " << endl; cout << endl; }



    Int_t branchParticleEntries = branchParticle->GetEntriesFast();
    // if(debug) {branchParticleEntries = 100;}// print only the first 100 particles of an event for debugging

    for(i=0; i < branchParticleEntries; ++i)
    {
       GenParticle *gen = (GenParticle*) branchParticle->At(i);
       // if (printInfo)  cout<< " particle PID "  <<  gen->PID << " Status " << gen->Status << " Mass " <<   gen->Mass << endl;

       // **************************** / 
       // bquarks from hard scattering /
       // **************************** /  

       if (fabs(gen->PID)==5 && (gen->Status>=21 && gen->Status<=29)  ) 
       { 
         genBquarks.push_back(*gen);
         // if (printInfo && (gen->PT>10) ) cout<< " bquark "  <<  gen->PID << " " <<   gen->PT <<  " " << gen->Eta  << " " << gen->Phi << endl;
       }
       
       // ******************************************* /   
       // lepton is 11,13,15 PDG and final (status=1) /
       // ******************************************* /

       //if (((fabs(gen->PID)==11)||(fabs(gen->PID)==13)||(fabs(gen->PID)==15)) && (gen->Status==1) )
       if (((fabs(gen->PID)==11)||(fabs(gen->PID)==13)) && (gen->Status==1) )
       { 
         genLeptons.push_back(*gen);
         // if (printInfo && (gen->PT>10) ) cout<< " lepton "  <<  gen->PID << " " <<   gen->PT <<  " " << gen->Eta  << " " << gen->Phi << endl;
       }

       // ********************************************* /
       // neutrino is 12,14,16 PDG and final (status=1) /
       // ********************************************* /

       if (((fabs(gen->PID)==12)||(fabs(gen->PID)==14)||(fabs(gen->PID)==16) ) && (gen->Status==1)  )
       { 
         genNeutrinos.push_back(*gen);
         genMETx += gen->Px; genMETy += gen->Py; genMEz += gen->Pz;
         // if (printInfo) cout<< " neutrino "  <<  gen->PID << " " <<   gen.Px() <<  " " << gen.Py()  << " " << gen.Pz() << endl;
       }
    }

    // if (printInfo) cout << " genMET " << genMETx << " " << genMETy << endl;
   
// ***************************** /  
//       Electrons               /
// ***************************** / 
  
    // if (printInfo) {  cout << endl;  cout << " reconstructed event " << endl; cout << endl; }


    nelectrons=branchElectron->GetEntriesFast();   
    
    /* ******************************** */ 
    /* Loop over all electrons in event */
    /* ******************************** */

    // if (printInfo) printf("\n electrons=%d\n",nelectrons);

    for(i = 0; i < nelectrons; ++i)
    { 
      Electron *myelectron = (Electron*) branchElectron->At(i);
      //particle = (GenParticle*) electron->Particle.GetObject(); /*  particle->Eta - electron->Eta */
      // if (printInfo) cout<< " electron " << i << " " <<   myelectron->PT <<  " " << myelectron->Eta  << " " << myelectron->Phi << endl;

      Lepton mylepton;
      mylepton.PT = myelectron->PT;
      mylepton.Eta =(myelectron->Eta);
      mylepton.Phi = (myelectron->Phi);
      mylepton.T = (myelectron->T);
      mylepton.Charge = (myelectron->Charge);
      mylepton.LeptonID = (11);
      mylepton.Particle = (myelectron->Particle);
      // mylepton.EhadOverEem = (myelectron->EhadOverEem);
      // mylepton.IsolationVar = (myelectron->IsolationVar);
      // mylepton.IsolationVarRhoCorr = (myelectron->IsolationVarRhoCorr);
      // mylepton.SumPtCharged = (myelectron->SumPtCharged);
      // mylepton.SumPtNeutral = (myelectron->SumPtNeutral);
      // mylepton.SumPtChargedPU = (myelectron->SumPtChargedPU);
      // mylepton.SumPt = (myelectron->SumPt);

      //if (printInfo) cout<< " lepton " << i << " " <<   mylepton.PT <<  " " << mylepton.Eta  << " " << mylepton.Phi << endl;

      Leptons.push_back(mylepton); 
     
    }

// ***************************** /  
//       Muons                   /
// ***************************** / 

    nmuons = branchMuon->GetEntriesFast(); 

    // **************************** /
    // Loop over all muons in event /
    // **************************** /

    // if (printInfo) printf("\n muons=%d\n",nmuons);

    for(i = 0; i < nmuons; ++i)
    {
      Muon *mymuon = (Muon*) branchMuon->At(i);
      //particle = (GenParticle*) muon->Particle.GetObject();  // particle->Eta - muon->Eta 
      // if (printInfo) cout<< " muon " <<  i << " " <<  mymuon->PT <<  " " << mymuon->Eta  << " " << mymuon->Phi << endl;

      Lepton yourlepton;
      yourlepton.PT=(mymuon->PT);
      yourlepton.Eta=(mymuon->Eta);
      yourlepton.Phi=(mymuon->Phi);
      yourlepton.T=(mymuon->T);
      yourlepton.Charge=(mymuon->Charge);
      yourlepton.LeptonID=(13);
      yourlepton.Particle=(mymuon->Particle);
      // yourlepton.EhadOverEem=(0.0);
      // yourlepton.IsolationVar=(mymuon->IsolationVar);
      // yourlepton.IsolationVarRhoCorr=(mymuon->IsolationVarRhoCorr);
      // yourlepton.SumPtCharged=(mymuon->SumPtCharged);
      // yourlepton.SumPtNeutral=(mymuon->SumPtNeutral);
      // yourlepton.SumPtChargedPU=(mymuon->SumPtChargedPU);
      // yourlepton.SumPt=(mymuon->SumPt);

      //if (printInfo) cout<< " lepton " << i << " " <<   yourlepton.PT <<  " " << yourlepton.Eta  << " " << yourlepton.Phi << endl;

      Leptons.push_back(yourlepton);

    }


// ***************************** /
//       Leptons                 /
// ***************************** /

   nleptons = Leptons.size();

 
   // sort them 
   std::sort(Leptons.begin(), Leptons.end(), sortLeptonsOperator);

   // // loop leptons 
   // for (unsigned i=0; i< Leptons.size(); i++)  
   // { 
   //   //if (printInfo) cout << " lepton " <<  Leptons[i].PT  << " "  << Leptons[i].Eta << " " <<  Leptons[i].Phi << endl;
   //   if (printInfo) printf(" lepton %f %f %f \n", Leptons[i].PT,Leptons[i].Eta, Leptons[i].Phi);
   // }


   // analyze leading leptons 
   if(nleptons >= 2)
   {
      mostEnergeticLepton = Leptons[0];
      secEnergeticLepton  = Leptons[1];
   }


// ***************************** /  
//       Jets                    /
// ***************************** / 

    njets=branchJet->GetEntriesFast();

    // Loop over all jets in event 
    for(i = 0; i < njets; ++i) 
    {  
       Jet *myjet = (Jet*) branchJet->At(i);  
       //TLorentzVector myjet; myjet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
       // if (printInfo) cout<< " jet " <<   myjet->PT <<  " " << myjet->Eta  << " " << myjet->Phi << endl;
       Jets.push_back(*myjet);
    }

    // sort them 
    std::sort(Jets.begin(), Jets.end(), sortJetsOperator);

    // Analyse two leading jets 
    if(njets >= 2)
    {
     mostEnergeticJet = Jets[0];
     secEnergeticJet  = Jets[1];
    }
   
// *********************************  / 
// met in event.  Carefull: MissingET /
// has 3 components as it refers to   /
// missing energy in 3D not MET!      /
// *********************************  /

    MissingET *missingEnergy;

    if ( branchMissingET->GetEntriesFast()> 0)
    {
       missingEnergy = (MissingET*)  branchMissingET->At(0);
       TVector3 temp = missingEnergy->P4().Vect();
       missingEnergyX = temp.Px(); missingEnergyY = temp.Py();
       // if (printInfo) cout << " missing energy "<< temp.Px()<< " "<< temp.Py() << " "<<temp.Pz() << endl;
    }




// *********************************** /
//                                     /
// event selection - generated         /
//                                     /
// *********************************** /
  
  int ngeneratedLeptons = Leptons.size();
  int ngeneratedBquarks = genBquarks.size();

  std::sort(genBquarks.begin(), genBquarks.end(), sortGenParticleOperator);
  if (ngeneratedBquarks>=2)
  {
    mostEnergeticBquark  = genBquarks[0];
    secEnergeticBquark   = genBquarks[1];  
  }

  std::sort(genLeptons.begin(), genLeptons.end(), sortGenParticleOperator);
  if (ngeneratedLeptons>=2)
  {
    mostEnergeticGenLepton = genLeptons[0];     
    secEnergeticGenLepton  = genLeptons[1];

  }

  Double_t Mbb_gen = (mostEnergeticBquark.P4()+secEnergeticBquark.P4()).M();
  Double_t Mll_gen = (mostEnergeticGenLepton.P4()+secEnergeticGenLepton.P4()).M();

  int nneutrinos = genNeutrinos.size();

  //  if (generatedInfo && (ngeneratedBquarks>=2) && (ngeneratedLeptons ==2) && (nneutrinos == 2) && Mbb_gen > 95 && Mbb_gen < 140 && Mll_gen < 65 && pow(genMETx*genMETx + genMETy*genMETy, 0.5) > 30)
  if (generatedInfo && (ngeneratedBquarks>=2) && (ngeneratedLeptons ==2) && Mbb_gen > 95 && Mbb_gen < 140 && Mll_gen < 65 && pow(genMETx*genMETx + genMETy*genMETy, 0.5) > 30)
   {
      if (printInfo){ cout << endl; cout << "------------  event " << entry << " -------------- " << endl; cout << endl;}

      if (printInfo) { cout << " genereated event " << endl; cout << endl; }


      Int_t branchParticleEntries = branchParticle->GetEntriesFast();
      if(debug) {branchParticleEntries = 100;}// print only the first 20 particles of an event for debugging

      for(i=0; i < branchParticleEntries; ++i)
      {
         GenParticle *gen = (GenParticle*) branchParticle->At(i);
         // if (printInfo)  cout<< " particle PID "  <<  gen->PID << " Status " << gen->Status << " Mass " <<   gen->Mass << endl;

         // ************************** / 
         // Higgs from hard scattering /
         // ************************** /  

         if (fabs(gen->PID)==25 && (gen->Status>=21 && gen->Status<=29)  ) 
         { 
           TLorentzVector myparticle; myparticle.SetPtEtaPhiM(gen->PT, gen->Eta, gen->Phi, gen->Mass); 
           if (printInfo && (gen->PT>10) ) cout<< " Higgs "  <<  gen->PID << " " <<   myparticle.Px() <<  " " << myparticle.Py()  << " " << myparticle.Pz() << endl;
         }

         // **************************** / 
         // bquarks from hard scattering /
         // **************************** /  

         if (fabs(gen->PID)==5 && (gen->Status>=21 && gen->Status<=29)  ) 
         { 
           TLorentzVector myparticle; myparticle.SetPtEtaPhiM(gen->PT, gen->Eta, gen->Phi, gen->Mass); 
           if (printInfo && (gen->PT>10) ) cout<< " bquark "  <<  gen->PID << " " <<   myparticle.Px() <<  " " << myparticle.Py()  << " " << myparticle.Pz() << endl;
         }
         
         // ******************************************* /   
         // lepton is 11,13,15 PDG and final (status=1) /
         // ******************************************* /

         //if (((fabs(gen->PID)==11)||(fabs(gen->PID)==13)||(fabs(gen->PID)==15)) && (gen->Status==1) )
         if (((fabs(gen->PID)==11)||(fabs(gen->PID)==13)) && (gen->Status==1) )
         { 
           TLorentzVector myparticle; myparticle.SetPtEtaPhiM(gen->PT, gen->Eta, gen->Phi, gen->Mass); 
           if (printInfo && (gen->PT>10) ) cout<< " lepton "  <<  gen->PID << " " <<   myparticle.Px() <<  " " << myparticle.Py()  << " " << myparticle.Pz() << endl;
         }

         // ********************************************* /
         // neutrino is 12,14,16 PDG and final (status=1) /
         // ********************************************* /

         if (((fabs(gen->PID)==12)||(fabs(gen->PID)==14)||(fabs(gen->PID)==16))&& (gen->Status==1)  )
         { 
           TLorentzVector myparticle; myparticle.SetPtEtaPhiM(gen->PT, gen->Eta, gen->Phi, gen->Mass);
           // genMETx += myparticle.Px(); genMETy += myparticle.Py();
           if (printInfo) cout<< " neutrino "  <<  gen->PID << " " <<   myparticle.Px() <<  " " << myparticle.Py()  << " " << myparticle.Pz() << endl;
         }
      }
      if (printInfo) cout << " genMET " << genMETx << " " << genMETy << endl;

      // //  sort generated leptons - get first 2
      // std::sort(genLeptons.begin(), genLeptons.end(), sortGenParticleOperator);
      // mostEnergeticGenLepton = genLeptons[0];     
      // secEnergeticGenLepton  = genLeptons[1];
 
      // //  sort generated bquarks - get first 2 
      // std::sort(genBquarks.begin(), genBquarks.end(), sortGenParticleOperator);
      // mostEnergeticBquark  = genBquarks[0];
      // secEnergeticBquark   = genBquarks[1];  

      //  sort generated bquarks - get first 2 
      std::sort(genNeutrinos.begin(), genNeutrinos.end(), sortGenParticleOperator);
      mostEnergeticGenNeutrino  = genNeutrinos[0];
      secEnergeticGenNeutrino   = genNeutrinos[1];  

      // now fill in vars 
      MET_x      = genMETx;
      MET_y      = genMETy;
      jet1_px    = mostEnergeticBquark.Px;    jet1_py    = mostEnergeticBquark.Py;    jet1_pz    = mostEnergeticBquark.Pz;
      jet2_px    = secEnergeticBquark.Px;     jet2_py    = secEnergeticBquark.Py;     jet2_pz    = secEnergeticBquark.Pz;
      lepton1_px = mostEnergeticGenLepton.Px; lepton1_py = mostEnergeticGenLepton.Py; lepton1_pz = mostEnergeticGenLepton.Pz;
      lepton2_px = secEnergeticGenLepton.Px;  lepton2_py = secEnergeticGenLepton.Py;  lepton2_pz = secEnergeticGenLepton.Pz;
      neutrino1_px = mostEnergeticGenNeutrino.Px; neutrino1_py = mostEnergeticGenNeutrino.Py; neutrino1_pz = mostEnergeticGenNeutrino.Pz;
      neutrino2_px = secEnergeticGenNeutrino.Px; neutrino2_py = secEnergeticGenNeutrino.Py; neutrino2_pz = secEnergeticGenNeutrino.Pz;

      // print generated event numbers 
      if (printInfo) 
      { 
        printf("\n\n generated event selected\n");
        printf("\n MET %lf %lf", MET_x,MET_y); 
        printf("\n Jet1 %lf %lf %lf", jet1_px, jet1_py, jet1_pz);
        printf("\n Jet2 %lf %lf %lf", jet2_px, jet2_py, jet2_pz);
        printf("\n Lep1 %lf %lf %lf", lepton1_px, lepton1_py, lepton1_pz);
        printf("\n Lep2 %lf %lf %lf", lepton2_px, lepton2_py, lepton2_pz); 
        printf("\n Neutr1 %lf %lf %lf", neutrino1_px, neutrino1_py, neutrino1_pz);
        printf("\n Neutr2 %lf %lf %lf \n", neutrino2_px, neutrino2_py, neutrino2_pz); 

      }

      // fill tree 
      selectedEvents++;
      t1->Fill();
      // if(selectedEvents > 10) {break;}

   }


// *********************************** /
//                                     /
// event selection - reconstructed     /
//                                     /
// *********************************** /
Double_t Mjj = (mostEnergeticJet.P4()+secEnergeticJet.P4()).M();
Double_t Mll = (mostEnergeticLepton.P4()+secEnergeticLepton.P4()).M();

    // if (!(generatedInfo) && (njets >= 2) && (nleptons >=2) && (mostEnergeticJet.P4().Pt()>30) && (secEnergeticJet.P4().Pt()>30) && (mostEnergeticLepton.P4().Pt()>30)  && (secEnergeticLepton.P4().Pt()>30)&& missingEnergyX*missingEnergyX + missingEnergyY*missingEnergyY > 1)
    if (!(generatedInfo) && (njets >= 2) && (nleptons ==2) && Mjj > 95 && Mjj < 140 && Mll < 65 && pow(missingEnergyX*missingEnergyX + missingEnergyY*missingEnergyY, 0.5) > 30)
    {

      if (printInfo){ cout << endl; cout << "------------  event " << entry << " -------------- " << endl; cout << endl;}

      if (printInfo) { cout << " genereated event " << endl; cout << endl; }


      Int_t branchParticleEntries = branchParticle->GetEntriesFast();
      if(debug) {branchParticleEntries = 100;}// print only the first 20 particles of an event for debugging

      for(i=0; i < branchParticleEntries; ++i)
      {
         GenParticle *gen = (GenParticle*) branchParticle->At(i);
         // if (printInfo)  cout<< " particle PID "  <<  gen->PID << " Status " << gen->Status << " Mass " <<   gen->Mass << endl;

         // ************************** / 
         // Higgs from hard scattering /
         // ************************** /  

         if (fabs(gen->PID)==25 && (gen->Status>=21 && gen->Status<=29)  ) 
         { 
           TLorentzVector myparticle; myparticle.SetPtEtaPhiM(gen->PT, gen->Eta, gen->Phi, gen->Mass); 
           if (printInfo && (gen->PT>10) ) cout<< " Higgs "  <<  gen->PID << " " <<   myparticle.Pt() <<  " " << myparticle.Eta()  << " " << myparticle.Phi() << endl;
         }

         // **************************** / 
         // bquarks from hard scattering /
         // **************************** /  

         if (fabs(gen->PID)==5 && (gen->Status>=21 && gen->Status<=29)  ) 
         { 
           TLorentzVector myparticle; myparticle.SetPtEtaPhiM(gen->PT, gen->Eta, gen->Phi, gen->Mass); 
           if (printInfo && (gen->PT>10) ) cout<< " bquark "  <<  gen->PID << " " <<   myparticle.Pt() <<  " " << myparticle.Eta()  << " " << myparticle.Phi() << endl;
         }
         
         // ******************************************* /   
         // lepton is 11,13,15 PDG and final (status=1) /
         // ******************************************* /

         //if (((fabs(gen->PID)==11)||(fabs(gen->PID)==13)||(fabs(gen->PID)==15)) && (gen->Status==1) )
         if (((fabs(gen->PID)==11)||(fabs(gen->PID)==13)) && (gen->Status==1) )
         { 
           TLorentzVector myparticle; myparticle.SetPtEtaPhiM(gen->PT, gen->Eta, gen->Phi, gen->Mass); 
           if (printInfo && (gen->PT>10) ) cout<< " lepton "  <<  gen->PID << " " <<   myparticle.Pt() <<  " " << myparticle.Eta()  << " " << myparticle.Phi() << endl;
         }

         // ********************************************* /
         // neutrino is 12,14,16 PDG and final (status=1) /
         // ********************************************* /

         if (((fabs(gen->PID)==12)||(fabs(gen->PID)==14)||(fabs(gen->PID)==16) ) && (gen->Status==1)  )
         { 
           TLorentzVector myparticle; myparticle.SetPtEtaPhiM(gen->PT, gen->Eta, gen->Phi, gen->Mass);
           // genMETx += myparticle.Px(); genMETy += myparticle.Py();
           if (printInfo) cout<< " neutrino "  <<  gen->PID << " " <<   myparticle.Px() <<  " " << myparticle.Py()  << " " << myparticle.Pz() << endl;
         }
      }
      if (printInfo) cout << " genMET " << genMETx << " " << genMETy << endl;

         // now fill in vars 
         MET_x      = missingEnergyX;
         MET_y      = missingEnergyY;

         jet1_px    = mostEnergeticJet.P4().Px();    jet1_py    = mostEnergeticJet.P4().Py();    jet1_pz    = mostEnergeticJet.P4().Pz(); 
         jet1_btagging = mostEnergeticJet.BTag;
         jet1_charge   = mostEnergeticJet.Charge;
         jet1_HoE      = mostEnergeticJet.EhadOverEem;
         jet1_ncharged = mostEnergeticJet.NCharged;
         jet1_nneutral = mostEnergeticJet.NNeutrals;

         jet2_px    = secEnergeticJet.P4().Px();     jet2_py    = secEnergeticJet.P4().Py();     jet2_pz    = secEnergeticJet.P4().Pz();
         jet2_btagging = secEnergeticJet.BTag;
         jet2_charge   = secEnergeticJet.Charge;
         jet2_HoE      = secEnergeticJet.EhadOverEem;
         jet2_ncharged = secEnergeticJet.NCharged;
         jet2_nneutral = secEnergeticJet.NNeutrals;


         lepton1_px = mostEnergeticLepton.P4().Px(); lepton1_py = mostEnergeticLepton.P4().Py(); lepton1_pz = mostEnergeticLepton.P4().Pz();
         lepton1_charge         =  mostEnergeticLepton.Charge;
         lepton1_HoE            =  mostEnergeticLepton.EhadOverEem;
         lepton1_IsoVar         =  mostEnergeticLepton.IsolationVar;
         lepton1_IsoRho         =  mostEnergeticLepton.IsolationVarRhoCorr;
         lepton1_SumPtCharged   =  mostEnergeticLepton.SumPtCharged;
         lepton1_SumPtNeutral   =  mostEnergeticLepton.SumPtNeutral;
         lepton1_SumPtChargedPU =  mostEnergeticLepton.SumPtChargedPU;
         lepton1_SumPt          =  mostEnergeticLepton.SumPt;

         lepton2_px = secEnergeticLepton.P4().Px();  lepton2_py = secEnergeticLepton.P4().Py();  lepton2_pz = secEnergeticLepton.P4().Pz();
         lepton2_charge         =  secEnergeticLepton.Charge;
         lepton2_HoE            =  secEnergeticLepton.EhadOverEem;
         lepton2_IsoVar         =  secEnergeticLepton.IsolationVar;
         lepton2_IsoRho         =  secEnergeticLepton.IsolationVarRhoCorr;
         lepton2_SumPtCharged   =  secEnergeticLepton.SumPtCharged;
         lepton2_SumPtNeutral   =  secEnergeticLepton.SumPtNeutral;
         lepton2_SumPtChargedPU =  secEnergeticLepton.SumPtChargedPU;
         lepton2_SumPt          =  secEnergeticLepton.SumPt;

         neutrino1_px = 0; neutrino1_py = 0; neutrino1_pz = 0;
         neutrino2_px = 0; neutrino2_py = 0; neutrino2_pz = 0;

         /// print reco event 
         if (printInfo) {
         printf("\n\n reco event selected\n");
         printf("\n MET %lf %lf", MET_x,MET_y);
         //printf("\n %lf %lf %lf", jet1_px, jet1_py, jet1_pz);
         //printf("\n %lf %lf %lf", jet2_px, jet2_py, jet2_pz);
         //printf("\n %lf %lf %lf", lepton1_px, lepton1_py, lepton1_pz);
         //printf("\n %lf %lf %lf", lepton2_px, lepton2_py, lepton2_pz);
         printf("\n Jet1 %lf %lf %lf", mostEnergeticJet.P4().Pt(), mostEnergeticJet.P4().Eta(), mostEnergeticJet.P4().Phi());
         printf("\n Jet2 %lf %lf %lf", secEnergeticJet.P4().Pt(),  secEnergeticJet.P4().Eta(),  secEnergeticJet.P4().Phi());
         printf("\n Lep1 %lf %lf %lf", mostEnergeticLepton.P4().Pt(), mostEnergeticLepton.P4().Eta(), mostEnergeticLepton.P4().Phi());
         printf("\n Lep2 %lf %lf %lf \n ", secEnergeticLepton.P4().Pt(),  secEnergeticLepton.P4().Eta(),  secEnergeticLepton.P4().Phi());
  
         }

         // fill tree 
         selectedEvents++;
         t1->Fill(); 

     } 

   

 
 

  } /* event loop */


 /* ****************** */
 /* store tree results */
 /* ****************** */

 t1->Write();
 treeResults->Close();

 cout << "** Chain contains " << allEntries << " events from which " << selectedEvents << " were selected" <<  endl;
    
 if (printInfo) {
 cout << endl;
 cout << "** Exiting..." << endl; }

 delete treeReader;
 delete chain;
 gROOT->ProcessLine(".q");

}

//------------------------------------------------------------------------------

