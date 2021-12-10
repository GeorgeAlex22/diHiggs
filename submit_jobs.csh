#! usr/bin/csh

setenv STARTING_DIR $PWD

# This could later become a foreach loop in the various physics processes
# setenv PHYSICS_PROCESS pptohh
# setenv PHYSICS_PROCESS ttbar1
setenv PHYSICS_PROCESS $1

setenv SUBMIT_DIR ${STARTING_DIR}/cfgFiles/submit
setenv SUB_INPUT_LIST `ls -l ${SUBMIT_DIR}/${PHYSICS_PROCESS} | awk '{print $9}' | grep .sub`
# | awk 'FNR<=3'`

setenv X509_USER_PROXY /afs/cern.ch/user/a/alexandg/x509up_u17268
voms-proxy-init

@ index=0
foreach INPUT (${SUB_INPUT_LIST})
    # echo "SUBMITTING ${INPUT} in my imagination..."
    condor_submit ${STARTING_DIR}/cfgFiles/submit/${PHYSICS_PROCESS}/${INPUT}
    @ index += 1
end

#echo "*Returning to ${STARTING_DIR}\n"
echo "Submited $index jobs."
cd $STARTING_DIR
