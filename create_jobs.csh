#!/usr/bin/csh

setenv STARTING_DIR $PWD

# This could later become a foreach loop in the various physics processes
# setenv PHYSICS_PROCESS pptohh
# setenv PHYSICS_PROCESS ttbar1
setenv PHYSICS_PROCESS $1

setenv INPUT_PATH /eos/cms/store/group/phys_b2g/giorgos/pheno/rootFiles/$PHYSICS_PROCESS/trees 
# setenv INPUT_PATH  /eos/user/a/alexandg/public/DelphesTrees/$PHYSICS_PROCESS/trees
setenv EOS_PATH /eos/user/a/alexandg/public/EOS.diHiggs

cd $EOS_PATH
echo "**Entering ${EOS_PATH}"
rm -rf rootFiles/$PHYSICS_PROCESS/trees 
mkdir -p rootFiles/$PHYSICS_PROCESS/trees rootFiles/$PHYSICS_PROCESS/histos

cd ${STARTING_DIR}
echo "**Entering ${STARTING_DIR}"
mkdir -p cfgFiles/error cfgFiles/jobs cfgFiles/log cfgFiles/output cfgFiles/submit

foreach subDIR (`ls ${STARTING_DIR}/cfgFiles`)
    cd ${STARTING_DIR}/cfgFiles/${subDIR}
    rm -rf ${PHYSICS_PROCESS}
    mkdir -p ${PHYSICS_PROCESS}
    echo "**Created ${STARTING_DIR}/cfgFiles/${subDIR}/${PHYSICS_PROCESS}"
    cd ${STARTING_DIR}
end

setenv JOB_INPUT_PATH ${STARTING_DIR}/cfgFiles/jobs/${PHYSICS_PROCESS}
setenv SUBMIT_INPUT_PATH ${STARTING_DIR}/cfgFiles/submit/${PHYSICS_PROCESS}

cd $STARTING_DIR
echo "**Entering ${STARTING_DIR}"

# Clear any files in the job folder
foreach file (`ls ${JOB_INPUT_PATH}`)
    echo "**Deleting ${JOB_INPUT_PATH}/${file}\n"
end

# rm -f $JOB_INPUT_PATH/*
rm -r $JOB_INPUT_PATH
mkdir $JOB_INPUT_PATH


### REMOVE ME
# Later remove the  | awk 'FNR<=X' at the end of next line.
if ($1 == "ttbar1") then
    setenv INPUT_LIST `ls -l $INPUT_PATH | awk '{print $9}' | grep .hep.root | awk 'FNR<=300' `    
else
    setenv INPUT_LIST `ls -l $INPUT_PATH | awk '{print $9}' | grep .hep.root`
endif

# | awk 'FNR<=3'`

foreach INPUT ( ${INPUT_LIST} )
    setenv jobfile ${JOB_INPUT_PATH}/${INPUT}.job
    setenv inputfile root://eoscms.cern.ch/${INPUT_PATH}/${INPUT}
    setenv outputfile "root://eosuser.cern.ch/${EOS_PATH}/rootFiles/${PHYSICS_PROCESS}/trees/solved_${INPUT}"
    sed 's@INPUTFILE@'"${inputfile}"'@g' templates/job_template >  tmp_job_template
    sed 's@OUTPUTFILE@'"${outputfile}"'@g' tmp_job_template > ${jobfile}
    echo "CREATED ${jobfile} \n\tFROM ${INPUT_PATH}/${INPUT}\n"

    setenv subfile ${SUBMIT_INPUT_PATH}/${INPUT}.sub

    sed 's@INPUTFILE@'"${INPUT}"'@g' templates/sub_template > tmp_sub_template
    sed 's@PHYSICSPROCESS@'"${PHYSICS_PROCESS}"'@g' tmp_sub_template > ${subfile}
    echo "CREATED ${subfile} \n\tFOR ${jobfile}\n"
    
end
rm -f tmp_job_template

rm -f tmp_sub_template


echo "*Returning to ${STARTING_DIR}\n"
cd $STARTING_DIR

