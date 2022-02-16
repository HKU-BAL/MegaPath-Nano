#!/usr/bin/env bash
set -e
set -m
set -o pipefail

SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
NANO_DIR=$(dirname ${SCRIPT_PATH})

# configurations
DB=${NANO_DIR}/db

# default parameters
THREADS=24
TARGET_SEQ_ID=NC_000962.3
CLAIR_ENSEMBLE_MODEL=$SCRIPT_PATH/Clair-ensemble/model/model-000016

while getopts "1:r:t:d:" option; do
	case "${option}" in
		r) READ=${OPTARG};;
		t) THREADS=${OPTARG};; 
		d) DB=${OPTARG};;
		*) exit 1;;
	esac
done


if [ -z "${READ}" ] ; then
   echo "Usage: $0 -r <read.fq> [options]"
   echo "    -t  number of threads [24]"
   echo "    -d  database directory [${SCRIPT_PATH}/db]"
   exit 1
fi

PREFIX=`basename $READ`
TARGET_IDX=$DB/amplicon/Mycobacterium_tuberculosis_H37Rv_genome_v3.fasta

run_minimap2(){
MINIMAP2_PREFIX=`basename $2`
minimap2 -ax map-ont $1 $2 | samtools sort -o $MINIMAP2_PREFIX.bam 
samtools index $MINIMAP2_PREFIX.bam
}

# 0. MegaPath-Nano filtering
if [ -e ${PREFIX}.mpn.done ]; then
	echo "Skipping MegaPath-Nano filtering";
else
	echo "[TIMESTAMP] $(date) Running MegaPath-Nano filtering..."
	STARTTIME=$(date +%s)
    
    $SCRIPT_PATH/megapath_nano.py --query ${READ} --amplicon_filter_module --decoy_filter_alignment_score_percent_threshold 150 || true
    python $SCRIPT_PATH/lib/get_highestAS_read_match_target.py `echo $PREFIX|sed 's/\..*//'`.species.bam  
    seqtk subseq $PREFIX.adaptor_trimmed.trimmed_and_filtered.human_and_decoy_filtered `echo $PREFIX|sed 's/\..*//'`.species.bam.highestAS_read_match_target.list > $PREFIX.adaptor_trimmed.trimmed_and_filtered.human_and_decoy_filtered.taxonfiltered
    run_minimap2 $TARGET_IDX  $PREFIX.adaptor_trimmed.trimmed_and_filtered.human_and_decoy_filtered.taxonfiltered
	echo "[TIMESTAMP] $(date) Running MegaPath-Nano filtering... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] MegaPath-Nano filtering took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.mpn.done
fi

# 1 Variant calling with clair-ensemble
if [ -e ${PREFIX}.clair-ensemble.done ]; then
	echo "Skipping variant calling with clair-ensemble";
else
	STARTTIME=$(date +%s)
	CLAIR_ENS_INPUT=${PREFIX}.adaptor_trimmed.trimmed_and_filtered.human_and_decoy_filtered.taxonfiltered.bam

	echo "[TIMESTAMP] $(date) Variant calling with clair-ensemble..."
    bash $SCRIPT_PATH/runClair-ensemble.sh -b ${PREFIX}.adaptor_trimmed.trimmed_and_filtered.human_and_decoy_filtered.taxonfiltered.bam -r $TARGET_IDX -m $CLAIR_ENSEMBLE_MODEL -t $THREADS
	echo "[TIMESTAMP] $(date) Variant calling with clair-ensemble... Done"

	ENDTIME=$(date +%s)
	echo "[TIMER] Variant calling with clair-ensemble took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.clair-ensemble.done
fi

# 2 Realignment
if [ -e ${PREFIX}.realignment.done ]; then
	echo "Skipping realignment";
else
	STARTTIME=$(date +%s)
	REALIGNMENT_INPUT_BAM=${PREFIX}.adaptor_trimmed.trimmed_and_filtered.human_and_decoy_filtered.taxonfiltered.bam
	REALIGNMENT_INPUT_VCF=$REALIGNMENT_INPUT_BAM.clair-ensemble/output/$REALIGNMENT_INPUT_BAM/snp_and_indel.filtered.vcf.gz

	echo "[TIMESTAMP] $(date) Realignment..."
    bash $SCRIPT_PATH/realignment/realignment.sh \
        -b $REALIGNMENT_INPUT_BAM \
        -r $TARGET_IDX \
        -v $REALIGNMENT_INPUT_VCF \
        -c $TARGET_SEQ_ID \
        -p "ont" \
        -t $THREADS \
        -o `readlink -f ${REALIGNMENT_INPUT_BAM}`.clair-ensemble_realigned

	echo "[TIMESTAMP] $(date) Realignment... Done"

	ENDTIME=$(date +%s)
	echo "[TIMER] Realignment took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.realignment.done
fi
