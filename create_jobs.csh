#!/usr/bin/csh

setenv STARTING_DIR $PWD

# This could later become a foreach loop in the various physics processes
# setenv PHYSICS_PROCESS pptohh
# setenv PHYSICS_PROCESS ttbar1
setenv PHYSICS_PROCESS $1
setenv INPUT_PATH  /eos/cms/store/group/phys_b2g/giorgos/pheno/rootFiles/$PHYSICS_PROCESS/trees 
setenv EOS_PATH /eos/user/a/alexandg/public/EOS.diHiggs

cd $EOS_PATH
echo "**Entering ${EOS_PATH}"
mkdir -p rootFiles/$PHYSICS_PROCESS/trees rootFiles/$PHYSICS_PROCESS/histos

cd ${STARTING_DIR}
echo "**Entering ${STARTING_DIR}"
mkdir -p cfgFiles/error cfgFiles/jobs cfgFiles/log cfgFiles/output cfgFiles/submit

foreach subDIR (`ls ${STARTING_DIR}/cfgFiles`)
    cd ${STARTING_DIR}/cfgFiles/${subDIR}
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
rm -f $JOB_INPUT_PATH/*


### REMOVE ME
# Later remove the  | awk 'FNR<=X' at the end of next line.
setenv INPUT_LIST `ls -l $INPUT_PATH | awk '{print $9}' | grep .hep.root`

# | awk 'FNR<=3'`

foreach INPUT ( ${INPUT_LIST} )
    setenv jobfile ${JOB_INPUT_PATH}/${INPUT}.job
    setenv inputfile ${INPUT_PATH}/${INPUT}
    setenv outputfile "${EOS_PATH}/rootFiles/${PHYSICS_PROCESS}/trees/solved_${INPUT}"
    sed 's@INPUTFILE@'"${inputfile}"'@g' job_template >  tmp_job_template
    sed 's@OUTPUTFILE@'"${outputfile}"'@g' tmp_job_template > ${jobfile}
    echo "CREATED ${jobfile} \n\tFROM ${INPUT_PATH}/${INPUT}\n"

    setenv subfile ${SUBMIT_INPUT_PATH}/${INPUT}.sub

    sed 's@INPUTFILE@'"${INPUT}"'@g' sub_template > tmp_sub_template
    sed 's@PHYSICSPROCESS@'"${PHYSICS_PROCESS}"'@g' tmp_sub_template > ${subfile}
    echo "CREATED ${subfile} \n\tFOR ${jobfile}\n"
    
end
rm -f tmp_job_template
rm -f tmp_sub_template

echo "*Returning to ${STARTING_DIR}\n"
cd $STARTING_DIR

