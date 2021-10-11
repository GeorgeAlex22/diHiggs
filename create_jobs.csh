#!/usr/bin/csh

setenv STARTING_DIR $PWD

# This could later become a foreach loop in the various physics processes
setenv PHYSICS_PROCESS pptohh
setenv INPUT_PATH  /eos/cms/store/group/phys_b2g/giorgos/pheno/rootFiles/$PHYSICS_PROCESS/trees
setenv EOS_PATH /eos/user/a/alexandg/public/EOS.diHiggs

cd $EOS_PATH
echo "**Entering ${EOS_PATH}"
mkdir -p rootFiles/$PHYSICS_PROCESS/trees rootFiles/$PHYSICS_PROCESS/histos
mkdir -p cfgFiles/error cfgFiles/jobs cfgFiles/log cfgFiles/output

foreach subDIR (`ls ${EOS_PATH}/cfgFiles`)
    cd ${EOS_PATH}/cfgFiles/${subDIR}
    mkdir -p ${PHYSICS_PROCESS}
    echo "**Created ${EOS_PATH}/cfgFiles/${subDIR}/${PHYSICS_PROCESS}"
    cd ${EOS_PATH}
end

setenv JOB_INPUT_PATH ${EOS_PATH}/cfgFiles/jobs/${PHYSICS_PROCESS}

echo "**Leaving ${EOS_PATH}"
cd $STARTING_DIR
echo "**Entering ${STARTING_DIR}"

# Clear any files in the job folder
foreach file (`ls ${JOB_INPUT_PATH}`)
    echo "**Deleting ${JOB_INPUT_PATH}/${file}\n"
end
rm $JOB_INPUT_PATH/*


### REMOVE ME
# Later remove the  | awk 'FNR<=X' at the end of next line.
setenv INPUT_LIST `ls -l $INPUT_PATH | awk '{print $9}' | grep lhe.hep.root | awk 'FNR<=3'`

foreach INPUT ( ${INPUT_LIST} )
    setenv jobfile ${JOB_INPUT_PATH}/${INPUT}.job
    setenv inputfile ${INPUT_PATH}/${INPUT}
    setenv outputfile "${EOS_PATH}/rootFiles/${PHYSICS_PROCESS}/trees/solved_${INPUT}"
    sed 's@INPUTFILE@'"${inputfile}"'@g' job_template >  tmp_job_template
    sed 's@OUTPUTFILE@'"${outputfile}"'@g' tmp_job_template > ${jobfile}
    echo "CREATED ${jobfile} \n\tFROM ${INPUT_PATH}/${INPUT}\n"
end
rm -f tmp_job_template

echo "*Returning to ${STARTING_DIR}\n"
cd $STARTING_DIR

