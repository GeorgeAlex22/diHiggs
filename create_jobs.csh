#!/usr/bin/csh

setenv STARTING_DIR $PWD

setenv PHYSICS_PROCCESS pptohh
setenv INPUT_PATH  /eos/cms/store/group/phys_b2g/giorgos/pheno/rootFiles/$PHYSICS_PROCCESS/trees
setenv EOS_PATH /eos/user/a/alexandg/public/DemokritosInternship2021

cd $EOS_PATH
mkdir -p diHiggs/rootFiles/$PHYSICS_PROCCESS
setenv PHYSICS_PROCCESS_PATH $EOS_PATH/diHiggs/rootFiles/$PHYSICS_PROCCESS

cd $PHYSICS_PROCCESS_PATH
mkdir -p jobs output
setenv OUTPUT_PATH $PHYSICS_PROCCESS_PATH/output
setenv JOB_INPUT_PATH  $PHYSICS_PROCCESS_PATH/jobs


cd $STARTING_DIR

# Clear any files in the job folder
rm $JOB_INPUT_PATH/*

### REMOVE ME
# Later remove the  | awk 'FNR<=X' at the end of next line.
setenv INPUT_LIST `ls -l $INPUT_PATH | awk '{print $9}' | grep lhe.hep.root | awk 'FNR<=3'`
setenv SELECTION_PREFIX "selection"

foreach INPUT ( ${INPUT_LIST} )
    setenv jobfile ${JOB_INPUT_PATH}/${INPUT}.job

    echo '#\!/usr/bin/csh' > ${jobfile}
    echo 'cp -r selection ${_CONDOR_SCRATCH_DIR}' >> ${jobfile}
    echo 'cd ${_CONDOR_SCRATCH_DIR}' >> ${jobfile}
    echo 'echo "entering directory ./${_CONDOR_SCRATCH_DIR}"' >> ${jobfile}

cat >> ${jobfile} << @EOI    

cp ${INPUT_PATH}/${INPUT} selection_input.root

echo "..Copying ${INPUT_PATH}/${INPUT} -> ./selection_input.root"

source bin/setup.csh

root -l -b Selection.C

cd ..
@EOI

    echo '"exiting directory ./${_CONDOR_SCRATCH_DIR}"' >> ${jobfile}

cat >> $jobfile << @EOI

cp ./selection/selection_output.root ${OUTPUT_PATH}/${SELECTION_PREFIX}.${INPUT}'
echo "..Copying ./selection/selection_output.root -> ${OUTPUT_PATH}/${SELECTION_PREFIX}.${INPUT}"

@EOI
    echo "CREATED ${jobfile} \n\tFROM ${INPUT_PATH}/${INPUT}\n"
end

cd $STARTING_DIR