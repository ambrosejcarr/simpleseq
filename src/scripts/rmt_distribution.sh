#!/bin/bash


# parse command line arguments
while [[ $# > 1 ]]
do
key="$1"

case $key in
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
echo SAMPLE = "${SAMPLE}"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

mkdir -r ${SAMPLE}

if [ ! -f ${SAMPLE}/Aligned.out.sam ]; then
    echo "downloading aligment file from s3"
    aws s3 cp s3://ajc-data/${SAMPLE}/Aligned.out.sam ${SAMPLE}/ 2>&1 > /dev/null
fi

# make sure simpleseq exists
if hash rmt_distribution 2>/dev/null; then
    :
else
    echo "installing simpleseq"
    pip3 uninstall -y simpleseq
    git clone -q https://github.com/ambrosejcarr/simpleseq.git simpleseq/
    cd simpleseq
    pip3 install -e ./
    cd ../
fi

echo "creating molecule counts object"
rmt_distribution ${SAMPLE}/Aligned.out.sam

echo "uploading data"
aws s3 sync ${SAMPLE}/ s3://ajc-data/${SAMPLE}/ 2>&1 > /dev/null

echo "run completed"

