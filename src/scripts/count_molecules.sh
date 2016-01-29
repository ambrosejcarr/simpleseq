#!/bin/bash

# script to create count matrices

# Use > 1 to consume two arguments per pass in the loop (e.g. each
# argument has a corresponding value to go with it).
# Use > 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
# note: if this is set to > 0 the /etc/hosts part is not recognized ( may be a bug )
while [[ $# > 1 ]]
do
key="$1"

case $key in
    -s|--samfile)
    SAMFILE="$2"
    shift # past argument
    ;;
    -g|--gtf)
    GTF="$2"
    shift # past argument
    ;;
#    -a|--alignment-summary)
#    ALIGNMENT_SUMMARY="$2"
#    shift # past argument
#    ;;
#    --h5-output)
#    H5_OUTPUT="$2"
#    shift # past argument
#    ;;
    -o|--output)
    COUNT_OUTPUT="$2"
    shift # past argument
    ;;
    --sample)
    SAMPLE="$2"
    shift # past argument
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done
echo SAMFILE = "${SAMFILE}"
echo GTF = "${GTF}"
#echo ALIGNMENT_SUMMARY = "${ALIGNMENT_SUMMARY}"
#echo H5_OUTPUT = "${H5_OUTPUT}"
echo COUNT_OUTPUT = "${COUNT_OUTPUT}"
echo SAMPLE = "${SAMPLE}"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

# get seqc-0.1.2_cm_mod to create readarrays

echo "removing existing installation of seqc and simpleseq"
#pip3 uninstall -y seqc
pip3 uninstall -y simpleseq

#echo "downloading and installing seqc_cm_mod"
#aws s3 cp s3://ajc-data/software/seqc-0.1.2_cm_mod.tar.gz ./seqc-0.1.2_cm_mod.tar.gz
#tar -xzf seqc-0.1.2_cm_mod.tar.gz
#cd seqc-0.1.2_cm_mod
#pip3 install -e ./
#cd ../

echo "downloading and installing simpleseq"
git clone -q https://github.com/ambrosejcarr/simpleseq.git simpleseq/
cd simpleseq
pip3 install -e ./ 
cd ../
mkdir -r ${FOLDER}

echo "creating molecule counts object"
count_molecules -s ${SAMFILE} -g ${GTF} -o ${COUNT_OUTPUT} # -a ${ALIGNMENT_SUMMARY}

#echo "creating readarray"
#sam2ra.py -s ${SAMFILE} -g ${GTF} -o ${H5_OUTPUT}

#aws s3 cp ${H5_OUTPUT} s3://ajc-data/${SAMPLE}/${H5_OUTPUT}
aws s3 cp ${COUNT_OUTPUT} s3://ajc-data/${SAMPLE}/${COUNT_OUTPUT}

echo "run completed"
