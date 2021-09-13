#!/usr/bin/env bash
set -x
while getopts b:v:r:c:p:o:t: option
do
    case "${option}"
        in
        b) BAM_PATH=("`readlink -f ${OPTARG}`");;
        v) VCF_PATH=("`readlink -f ${OPTARG}`");;
        r) REF_PATH=`readlink -f ${OPTARG}`;;
        c) CONTIGS=${OPTARG};;
        p) PLATFORM=${OPTARG};;
		o) OUT_FOLDER=${OPTARG};;
        t) THREADS=${OPTARG};;
    esac
done

PYPY=pypy3
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
SCRIPT_DIR=$(dirname ${SCRIPT_PATH})
DATE_TIME=`date "+%Y%m%d_%H%M%S"`
mkdir -p ${OUT_FOLDER}
mkdir -p ${OUT_FOLDER}/tmp
cd ${OUT_FOLDER}

if [ "${PLATFORM}" = "illumina" ]
then
    echo "[INFO] Realign Illumina BAM"
    ILLUMINA_REALIGN_BAM_FOLDER=${OUT_FOLDER}/illumina_realignment
    IF_LOCAL_REALIGNMENT_ILLUMINA_REALIGN_FOLDER_FLAG="--ilmn_realign_folder"
    mkdir -p ${ILLUMINA_REALIGN_BAM_FOLDER}/
    time parallel --joblog ./illumina_reads_realignment.log -j${THREADS} \
    "${PYPY} ${SCRIPT_PATH}/realign_illumina_reads.py \
        --bam_fn {1} \
        --ref_fn ${REF_PATH} \
        --read_fn ${ILLUMINA_REALIGN_BAM_FOLDER}/{1/} \
        --samtools samtools \
        --ctgName ${CONTIGS}" ::: ${BAM_PATH[@]}

    ls ${ILLUMINA_REALIGN_BAM_FOLDER}/*.bam | parallel -j20 samtools index {}
fi

for i in ${!VCF_PATH[@]}
do
BAM_FILE=${BAM_PATH[i]}
gzip -fdc ${VCF_PATH[i]} | grep -v '#' | cut -f2 > ${OUT_FOLDER}/tmp/${i}_all_pos
readarray ALL_POS < ${OUT_FOLDER}/tmp/${i}_all_pos
echo "[INFO] Extract Candidates"
time parallel --joblog  ./${i}_create_tensor_pileup.log -j${THREADS} \
"${PYPY} ${SCRIPT_PATH}/local_realignment.py \
--bam_fn ${BAM_FILE} \
--ref_fn ${REF_PATH} \
--pos {1} \
--platform ${PLATFORM} \
--samtools ${SCRIPT_DIR}/samtools-1.13/samtools \
${IF_LOCAL_REALIGNMENT_ILLUMINA_REALIGN_FOLDER_FLAG} ${ILLUMINA_REALIGN_BAM_FOLDER} \
--ctgName ${CONTIGS} \
--output_fn ${OUT_FOLDER}/tmp/${i}_output_{#}" ::: ${ALL_POS[@]} |& tee  ./${i}_CTP.log

cat ${OUT_FOLDER}/tmp/${i}_output_* > ${OUT_FOLDER}/tmp/${i}_all_alt_info

${PYPY} ${SCRIPT_PATH}/extract_vcf_position.py \
--vcf_fn ${VCF_PATH[i]} \
--alt_fn ${OUT_FOLDER}/tmp/${i}_all_alt_info \
--output_fn ${OUT_FOLDER}/${i}_realign.vcf |& tee  ./${i}_EVP.log
done
