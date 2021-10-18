## Dependencies
* ### GCC
    We will be using [GCC][gcc]-4.9.1 for this project[^1].  
    If you are on _lxplus_, you can get it by sourcing
    ```
    source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6-gcc48-opt/setup.csh
    ```

* ### ROOT
    We will be using [ROOT][root]-6.06/06 for this project.  
    If you are on _lxplus_, you can get it by sourcing
    ```
    source source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.06.06-5e975/x86_64-slc6-gcc49-opt/bin/thisroot.csh
    ```
* ### LHAPDF
    This analysis uses [LHAPDF][lhapdf]-6.1.2 which has been installed in the $HOME directory.   
    It is needed for compiling and running main.x in ./hhDilep_pseudoanasol/ 

    ```
    gmake clean && gmake main.x && ./main.x 
    ```

    You can install LHAPDF-6.1.2 by running

    ```
    source /bin/installLHAPDF-6.1.2.csh     
    ```

* ### Delphes
    In order to run the Selection.C root script and installation of [Delphes][delphes]-3.4.0 in the $HOME directory is needed.
    You can install it by running 

    ```
    source /bin/installDelphes-3.4.0.csh 
    ```

    For me, it works with ROOT 6.06/06 and gcc 4.9.1 which I load form lxplus.  
    If you are not on lxplus, it should work by just downloading these versions on your local machine  
    and editing the appropriate lines in DelphesEnv.csh and installDelphes.csh 

## Analysis Workflow
Due to the relatively large phase-space scan of _main.x_, the analysis takes quite a bit of time to run on a single machine. That's why we rely on [HTCondor][condor]'s resources to allow us to use parallel processing by submitting multiple jobs at a time.  
If you want to study a certain _PHYSICS\_PROCESS_ = (_pptohh_, _ttbar1_, _pptoohhsm_, etc...) run 
```
source create_jobs.csh pptohhsm
source submit_jobs.csh pptohhsm
```
if _PHYSICS\_PROCESS_ = _pptohhsm_. Else, change it accordingly.  
A full list of the available _PHYSICS\_PROCESS_ is in [physics_processes.txt][physics_processes.txt]

After the jobs finish running, use 
```
root -l -b code/createHistos.C'("pptohhsm")'
```
to create histograms for the process _pptohhsm_.






[^1]:Except for when we are compiling LHAPDF, where we will use 4.7, but this process is essentially outside the project.  

[gcc]:https://gcc.gnu.org/
[root]:https://root.cern/
[lhapdf]:https://lhapdf.hepforge.org/
[delphes]:https://cp3.irmp.ucl.ac.be/projects/delphes
[condor]:https://htcondor.readthedocs.io/en/latest/
[physics_processes.txt]:https://github.com/GeorgeAlex22/diHiggs/blob/main/create_jobs.csh

