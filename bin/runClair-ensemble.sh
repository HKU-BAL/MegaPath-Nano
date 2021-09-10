# main purpose of this script:
# handle 3 dimensions for flexibility
# - multi-sample (does not like MES, just one sample) (SAMPLE_INDEX)
# - multi-depths (downsampled by me, or downsampled by variant bam) (DEPTH_INDEX / BAM_INDEX)
# - multi-models (different clair models) (MODEL_INDEX)
while getopts b:d:r:m:t: option
do
    case "${option}"
        in
        # assumed given bam file is sorted and indexed
        b) BAM_FILE_PATH=`readlink -f ${OPTARG}`;;
        #d) BED_FILE_PATH=`readlink -f ${OPTARG}`;;
        r) REFERENCE_FASTA_FILE_PATH=`readlink -f ${OPTARG}`;;
        m) old_IFS=$IFS
           echo $IFS
           IFS=','
           CLAIR_MODELS=($OPTARG)
           IFS=${old_IFS};;
        t) PARALLEL_THREADS=${OPTARG};;
    esac
done

SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
NANO_DIR=$(dirname ${SCRIPT_PATH})
DB=${SCRIPT_PATH}/db
#REFERENCE_FASTA_FILE_PATH=$DB/amplicon/Mycobacterium_tuberculosis_H37Rv_genome_v3.fasta
#CLAIR_MODELS=("$SCRIPT_PATH/Clair-ensemble/model/model-000016")
ROOT_FOLDER_PATH="$BAM_FILE_PATH.clair-ensemble"
CLAIR="$SCRIPT_PATH/Clair-ensemble/Clair.beta.ensemble.cpu/clair.py"
CLAIR_HOME_DIR=`dirname ${CLAIR}`
OVERLAP_VARIANT="${CLAIR_HOME_DIR}/clair/post_processing/overlap_variant.py"
ORIGINAL_BAM_PARENT_FOLDER_PATH=`dirname $BAM_FILE_PATH`
ORIGINAL_BAM_FILE_NAMES=(`basename $BAM_FILE_PATH`)
NO_OF_SAMPLES=${#ORIGINAL_BAM_FILE_NAMES[@]}


NO_OF_CLAIR_MODELS=${#CLAIR_MODELS[@]}

PARALLEL_THREADS="32"
SAMPLE_NAME="Amplicon"
AF_THRESHOLD="0.125"
NAME=16

WORKING_DIRECTORY="${ROOT_FOLDER_PATH}/"
mkdir -p ${ROOT_FOLDER_PATH}
mkdir ${WORKING_DIRECTORY}
PYTHON="/autofs/bal33/zxzheng/miniconda2/envs/clair3/bin/python"


# STEP 0: downsample bams (and store the downsample bams into one folder per bam)
TARGET_BAMS_PARENT_FOLDER_PATH="${WORKING_DIRECTORY}/input_bams"
mkdir ${TARGET_BAMS_PARENT_FOLDER_PATH}

#NO_OF_BAMS_PER_SAMPLE=${#DEPTHS[@]}
####
NO_OF_BAMS_PER_SAMPLE=1

for k in `seq 0 1 $((NO_OF_SAMPLES - 1))`
do
  SAMPLE_INDEX=`printf "%03d" ${k}`
  ORIGINAL_BAM_FILE_NAME=${ORIGINAL_BAM_FILE_NAMES[k]}
  ORIGINAL_BAM_FILE_PATH="${ORIGINAL_BAM_PARENT_FOLDER_PATH}/${ORIGINAL_BAM_FILE_NAME}"


  TARGET_BAM_FOLDER_PATH="${TARGET_BAMS_PARENT_FOLDER_PATH}/${SAMPLE_INDEX}"
  mkdir ${TARGET_BAM_FOLDER_PATH}

  for i in `seq 0 1 $((NO_OF_BAMS_PER_SAMPLE - 1))`
  do
    BAM_INDEX=`printf "%03d" ${i}`

    TARGET_BAM_FILE_PATH="${TARGET_BAM_FOLDER_PATH}/${BAM_INDEX}.bam"

   ln -s ${ORIGINAL_BAM_FILE_PATH} ${TARGET_BAM_FILE_PATH}
   ln -s ${ORIGINAL_BAM_FILE_PATH}.bai ${TARGET_BAM_FILE_PATH}.bai

  done

  # TODO: add symbolic link for 1.000 bam and 1.000.bam.bai
done

#

# STEP 1: output one call.sh for one bam, and call all call.sh with parallel
INTERMEDIATE_OUTPUT_FOLDER="${WORKING_DIRECTORY}/tmp_output"
mkdir ${INTERMEDIATE_OUTPUT_FOLDER}
cd ${INTERMEDIATE_OUTPUT_FOLDER}

for k in `seq 0 1 $((NO_OF_SAMPLES - 1))`
do
  SAMPLE_INDEX=`printf "%03d" ${k}`
  TARGET_BAM_FOLDER_PATH="${TARGET_BAMS_PARENT_FOLDER_PATH}/${SAMPLE_INDEX}"

  for i in `seq 0 1 $((NO_OF_BAMS_PER_SAMPLE - 1))`
  do
    BAM_INDEX=`printf "%03d" ${i}`
    INPUT_BAM_FILE_PATH="${TARGET_BAM_FOLDER_PATH}/${BAM_INDEX}.bam"

    for j in `seq 0 1 $((NO_OF_CLAIR_MODELS - 1))`
    do
      MODEL_INDEX=`printf "%03d" ${j}`

      CLAIR_MODEL="${CLAIR_MODELS[j]}"
      SCRIPT_OUTPUT_FOLDER="s${SAMPLE_INDEX}_m${MODEL_INDEX}_b${BAM_INDEX}"

      mkdir ${SCRIPT_OUTPUT_FOLDER}
      OUTPUT_PREFIX="${SCRIPT_OUTPUT_FOLDER}/tmp"

      ${PYTHON} ${CLAIR} callVarBamParallel \
      --chkpnt_fn "${CLAIR_MODEL}" \
      --ref_fn "${REFERENCE_FASTA_FILE_PATH}" \
      --bam_fn "${INPUT_BAM_FILE_PATH}" \
      --pysam_for_all_indel_bases \
      --output_for_ensemble \
      --includingAllContigs \
      --refChunkSize "10000000" \
      --threshold ${AF_THRESHOLD} \
      --sampleName "${SAMPLE_NAME}" \
      --output_prefix "${OUTPUT_PREFIX}" > ${SCRIPT_OUTPUT_FOLDER}/call.sh
    done
  done
done
( time cat */call.sh | parallel -j${PARALLEL_THREADS} ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.call.txt



# STEP 2: ensemble the result per sample using ensemble.cpp
ENSEMBLE_OUTPUT_FOLDER="${INTERMEDIATE_OUTPUT_FOLDER}/ensemble"
mkdir ${ENSEMBLE_OUTPUT_FOLDER}
MININUM_NO_OF_VOTE_FOR_VARIANT="$(((NO_OF_BAMS_PER_SAMPLE * NO_OF_CLAIR_MODELS + 2) / 2))"
rm ensemble_command.sh

for s in `seq 0 1 $((NO_OF_SAMPLES - 1))`
do
  SAMPLE_INDEX=`printf "%03d" ${s}`
  FILES=(`ls s${SAMPLE_INDEX}_m000_b000/*.vcf`)
  ORIGINAL_BAM_FILE_NAME=${ORIGINAL_BAM_FILE_NAMES[s]}

  ENSEMBLE_FOLDER_FOR_SAMPLE="${ENSEMBLE_OUTPUT_FOLDER}/${ORIGINAL_BAM_FILE_NAME}"
  mkdir "${ENSEMBLE_FOLDER_FOR_SAMPLE}"

  for i in "${!FILES[@]}"
  do
    TARGET_FILE_NAME=`basename ${FILES[i]}`
    CAT_COMMAND=""

    for j in `seq 0 1 $((NO_OF_BAMS_PER_SAMPLE - 1))`
    do
      BAM_INDEX=`printf "%03d" ${j}`

      for k in `seq 0 1 $((NO_OF_CLAIR_MODELS - 1))`
      do
        MODEL_INDEX=`printf "%03d" ${k}`
        FOLDER_NAME="s${SAMPLE_INDEX}_m${MODEL_INDEX}_b${BAM_INDEX}"
        CAT_COMMAND="${CAT_COMMAND} ${FOLDER_NAME}/${TARGET_FILE_NAME}"
      done
    done

    echo "cat ${CAT_COMMAND:1} | pypy3 ${CLAIR} ensemble --minimum_count_to_output ${MININUM_NO_OF_VOTE_FOR_VARIANT} > ${ENSEMBLE_FOLDER_FOR_SAMPLE}/${TARGET_FILE_NAME}" >> ensemble_command.sh
  done
done

( time cat ensemble_command.sh | parallel -j${PARALLEL_THREADS} ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.ensemble.txt
# rm ensemble_command.sh



# STEP 3: call_var the ensemble-d result into one vcf per one output file from previous step
VCF_OUTPUT_FOLDER="${WORKING_DIRECTORY}/output"
mkdir ${VCF_OUTPUT_FOLDER}
cd ${WORKING_DIRECTORY}

for s in `seq 0 1 $((NO_OF_SAMPLES - 1))`
do
  ORIGINAL_BAM_FILE_NAME=${ORIGINAL_BAM_FILE_NAMES[s]}
  ENSEMBLE_FOLDER_FOR_SAMPLE="${ENSEMBLE_OUTPUT_FOLDER}/${ORIGINAL_BAM_FILE_NAME}"

  TARGET_OUTPUT_FOLDER_PATH="${VCF_OUTPUT_FOLDER}/${ORIGINAL_BAM_FILE_NAME}"
  mkdir ${TARGET_OUTPUT_FOLDER_PATH}

  INPUT_FILES=(`ls ${ENSEMBLE_FOLDER_FOR_SAMPLE}/*.vcf`)

  for i in "${!INPUT_FILES[@]}"
  do
    FILE_NAME=`basename ${INPUT_FILES[i]}`
    ORIGINAL_BAM_FILE_PATH="${ORIGINAL_BAM_PARENT_FOLDER_PATH}/${ORIGINAL_BAM_FILE_NAME}"

    echo "cat ${ENSEMBLE_FOLDER_FOR_SAMPLE}/${FILE_NAME} | \
    python ${CLAIR} call_var \
    --chkpnt_fn "${CLAIR_MODEL}" \
    --ref_fn "${REFERENCE_FASTA_FILE_PATH}" \
    --bam_fn "${ORIGINAL_BAM_FILE_PATH}" \
    --call_fn "${TARGET_OUTPUT_FOLDER_PATH}/${FILE_NAME}" \
    --sampleName ${SAMPLE_NAME} \
    --pysam_for_all_indel_bases \
    --input_probabilities" >> output.sh
  done

  ( time cat output.sh | parallel -j${PARALLEL_THREADS} ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.output.txt
  rm output.sh
done



# STEP 4: merge vcfs into one vcf per sample (with overlap variant handling)
for s in `seq 0 1 $((NO_OF_SAMPLES - 1))`
do
  ORIGINAL_BAM_FILE_NAME=${ORIGINAL_BAM_FILE_NAMES[s]}
  ENSEMBLE_FOLDER_FOR_SAMPLE="${ENSEMBLE_OUTPUT_FOLDER}/${ORIGINAL_BAM_FILE_NAME}"

  TARGET_OUTPUT_FOLDER_PATH="${VCF_OUTPUT_FOLDER}/${ORIGINAL_BAM_FILE_NAME}"

  cd ${TARGET_OUTPUT_FOLDER_PATH}
  vcfcat ${TARGET_OUTPUT_FOLDER_PATH}/*.vcf | sort -k1,1V -k2,2n | bgziptabix snp_and_indel.vcf.gz

  zcat snp_and_indel.vcf.gz | \
  python ${OVERLAP_VARIANT} | \
  bgziptabix snp_and_indel.filtered.vcf.gz
done

