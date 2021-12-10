#include <dirent.h>

void renameTree()
{
    DIR* dir;
    struct dirent* ent;
    char eos_path[200];
    strcpy(eos_path,"/eos/user/a/alexandg/public/DelphesTrees/pptohhsm/trees/");
    char file[100];
    char filename[200];
    TFile* tmp_file;
    TTree* tmp_tree;
    Long64_t allEntries;
    int i =0;

    if ((dir = opendir (eos_path)) != NULL) 
    {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL)
        {
            strcpy(filename, " ");
            strcpy(file, ent->d_name);
            cout << "ent->d_name: " << ent->d_name << endl;
            if( ((string)file != ".") && ((string)file != ".."))
            {
                
                strcpy(filename, eos_path);
                strcat(filename, file);
                cout << "Processing file: " << filename << endl; 
                TFile *f = new TFile(filename,"update");
                TTree *T = (TTree*)f->Get("Delphes");
                T->Write("t1", TObject::kOverwrite);
                gDirectory->Delete("Delphes;1");
                f->Close();
            }

            // i++;
            // if(i > 5) {break;}
        }
        closedir (dir);
	} 

    else 
    {
        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
	}	

    gROOT->ProcessLine(".q");
}
