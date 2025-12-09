#!/bin/bash

#set -x #for complete debugging

set -euo pipefail

# STAGING_DIR is used to find the Singularity image (and databases)
# It is not used if this script is run inside a container
# Also, since this script does inference only, it doesn't need the databases
# this varible doesn't need to be changed

echo "$myjobdir"
echo "$myjobdir"
echo "$myjobdir"


readonly STAGING_DIR=/staging/groups/glbrc_alphafold/af3

# SINGIMG is used if we want to find a container. The script will 
# first look in a local directory and then the staging directory for it.
SINGIMG=""

# This file should be protected by your staging dir unix permissions
# can also set it with --model_param_file script option 
# or change it in the line below
# MODEL_PARAM_FILE=/staging/USERNAME/af3/weights/af3.bin.zst

# By default we copy containers and params to a working directory 
# on the local execute node. If working off a single file system 
# (scarcity) we can turn this off with --no_copy
COPY_BINARIES=1

VERBOSE_LEVEL=1 # 0 = silent, 1 = info, 2 = verbose

# By default we create a working directory called work.random
# Ideally this should be set to 
# --work_dir_ext $(ClusterId)_$(ProcID)  in the submit file
# but not needed if we are sure that multiple copies of this script
# will not overwrite each other
WORK_DIR_EXT="random"


function printstd() { echo "$@"; }
function printerr() { echo "ERROR: $@" 1>&2; }

function printinfo() {
  if [[ $VERBOSE_LEVEL -ge 1 ]]; then
    printstd "INFO: $@"
  fi
}
function printverbose() {
  if [[ $VERBOSE_LEVEL -ge 2 ]]; then
    printstd "DEBUG: $@"
  fi
}

ARGS="$@"

while [[ $# -gt 0 ]]; do
  case $1 in
     -w|--work_dir_ext)
      WORK_DIR_EXT="$2"
      printinfo "Setting WORK_DIR_EXT     : ${WORK_DIR_EXT}"
      shift # past argument
      shift # past value
      ;;
     -v|--verbose)
      VERBOSE=2
      printinfo "Setting PRINT_SUMMARY and VERBOSE on"
      shift # past argument
      ;;
     -s|--silent)
      VERBOSE=0
      printinfo "Setting PRINT_SUMMARY and VERBOSE off" # this will not print
      shift # past argument
      ;;
     -n|--no_copy)
      COPY_BINARIES=0
      printinfo "Not copying singularity container or model parameters"
      shift # past argument
      ;;
    -c|--container)
      SINGIMG="$2"
      printinfo "Will run inside container: $SINGIMG"
      shift # past argument
      shift # past value
      ;;
    -m|--model_param_file)
      MODEL_PARAM_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    -x|--user_specified_alphafold_options)
      USER_SPECIFIED_AF3_OPTIONS="$2"
      shift # past argument
      shift # past value
      ;;
    -u|--enable_unified_memory)
      USE_UNIFIED_MEMORY="true"
      shift # past argument
      ;;
    -*|--*)
      printerr "Unknown option $1"
      exit 1
      ;;
  esac
done

printinfo "Script         : $0"
printinfo "Running on     : `whoami`@`hostname`"
printinfo "Arguments      : $ARGS"
printinfo "Script dir     : $(dirname $0)"


readonly WORK_DIR="work.${WORK_DIR_EXT}"

printinfo "WORK_DIR       : `realpath $WORK_DIR`"
printverbose "Creating workdir and subdirectories : ${WORK_DIR}"
mkdir -p "${WORK_DIR}"
pushd "${WORK_DIR}" > /dev/null
mkdir -p af_input af_output models public_databases
popd

readonly WORK_INPUT_DIR="${WORK_DIR}/af_input"
# check that input files exist in the INPUT_DIR
printverbose "Extracting input files to : ${WORK_INPUT_DIR}" 
if compgen -G "*.data_pipeline.tar.gz"  > /dev/null; then
  # prepare input directory by expanding all the outputs from the data pipeline
  for filename in *.data_pipeline.tar.gz ;
  do
   printverbose "Extracting to WORK_INPUT_DIR : ${filename}"
   tar zxf "${filename}" -C ${WORK_INPUT_DIR}/
   rm "${filename}" # cleanup so that this file isn't accidentally sent back
  done
else
  printerr "Cannot find any input files matching " \
           "*.data_pipeline.tar.gz in directory : $(dirname $0)"
  exit 1
fi

## copy the container if we are going to run commands inside it
IMG_EXE_CMD="" # default is to not pipe commands through container
SINGIMG_PATH=""
if [[ -n "$SINGIMG" ]] ; then
  printverbose "Calling apptainer externally : ${SINGIMG}"
  if [[ "$COPY_BINARIES" -ne 0 ]] ; then
    printverbose "Copying container to WORK_DIR"
    if [ -f "$SINGIMG" ]; then 
      printverbose "Copying container from local directory"
      cp "${SINGIMG}" "${WORK_DIR}"/
      SINGIMG_PATH="${WORK_DIR}/${SINGIMG}"
    else # container is not in the local directory, check if it is in staging
      if [ -f ${STAGING_DIR}/${SINGIMG} ]; then
        printverbose "Copying container from staging directory"
        cp "${STAGING_DIR}/${SINGIMG}" "${WORK_DIR}"/ 
        SINGIMG_PATH="${WORK_DIR}/${SINGIMG}"
      else #not in staging
        printerr "Cannot find container to copy : $SINGIMG"
        exit 1
      fi # SINGIMG is not available to copy
    fi
  else # Do not copy binaries
    if [ -f "$SINGIMG" ]; then 
      printverbose "Container found. Not copying to workdir : $SINGIMG"
      SINGIMG_PATH="${SINGIMG}"
    else # container not found
      printerr "Trying to run in container (not found) : $SINGIMG"
      exit 1
    fi
  fi
  IMG_EXE_CMD="apptainer exec --nv ${SINGIMG_PATH}"
else
  printverbose "Not calling apptainer as we are inside the container"
fi

printinfo "SINGIMG_PATH   : $SINGIMG_PATH"
printinfo "IMG_EXE_CMD    : $IMG_EXE_CMD"

MODEL_PARAM_DIR=`realpath ${WORK_DIR}/models`
if [[ ${MODEL_PARAM_FILE} == *.zst ]]; then
  cat "${MODEL_PARAM_FILE}"  | \
        ${IMG_EXE_CMD} zstd  --decompress > ${WORK_DIR}/models/af3.bin
  printverbose "Decompressing model weights : ${MODEL_PARAM_FILE}"
else
  if [[ "$COPY_BINARIES" -ne 0 ]] ; then
    cp "${MODEL_PARAM_FILE}" ${WORK_DIR}/models/af3.bin
    printverbose "Copying model weights to workdir: ${MODEL_PARAM_FILE}"
  else
    printverbose "Not copying model weights"
    MODEL_PARAM_DIR=$(dirname ${MODEL_PARAM_FILE})
  fi
fi
printinfo "MODEL_PARAM_DIR: ${MODEL_PARAM_DIR}"

# Extra hand holding for CUDA_CAPABILITY 7.x devices
PYTHON_CC_EXEC_STR="import jax;
print(jax.local_devices(backend='gpu')[0].compute_capability)"
CUDA_CAPABILITY=`${IMG_EXE_CMD} python -c "${PYTHON_CC_EXEC_STR}"`
printinfo "CUDA_CAPABILITY: $CUDA_CAPABILITY"
EXTRA_RUN_ALPHAFOLD_FLAGS=""
EXTRA_APPTAINER_ENV=" " # add a space here so it isn't an empty newline below
if [  ${CUDA_CAPABILITY%.*} -eq 7 ] ; then
  printverbose "Setting extra flags for CUDA_CAPABILITY 7.x devices"
  EXTRA_RUN_ALPHAFOLD_FLAGS="--flash_attention_implementation=xla"
  export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
  EXTRA_APPTAINER_ENV="--env XLA_FLAGS=${XLA_FLAGS}"
fi

# Enable unified memory if requested
if [[ "${USE_UNIFIED_MEMORY:-}" == "true" ]] ; then
  printverbose "Enabling unified memory for Alphafold3 inference"
  export XLA_PYTHON_CLIENT_PREALLOCATE=false
  export TF_FORCE_UNIFIED_MEMORY=true
  export XLA_CLIENT_MEM_FRACTION=3.2
fi

printinfo "EXTRA_RUN_ALPHAFOLD_FLAGS: $EXTRA_RUN_ALPHAFOLD_FLAGS"
printinfo "EXTRA_APPTAINER_ENV      : $EXTRA_APPTAINER_ENV"
printinfo "USER_SPECIFIED_AF3_OPTIONS: $USER_SPECIFIED_AF3_OPTIONS"

exitcode=0

if [[ -n "$SINGIMG" ]] ; then # use apptainer to run the container
  apptainer exec \
     --bind "${WORK_DIR}/af_input":/root/af_input \
     --bind "${WORK_DIR}/af_output":/root/af_output \
     --bind "${MODEL_PARAM_DIR}":/root/models \
     --bind "${WORK_DIR}/public_databases":/root/public_databases \
     ${EXTRA_APPTAINER_ENV} \
     --cwd /app/alphafold \
     --nv \
     ${SINGIMG_PATH} \
     python run_alphafold.py \
     --db_dir=/root/public_databases \
     --run_data_pipeline=false \
     --run_inference=true \
     --input_dir=/root/af_input \
     --model_dir=/root/models \
     --output_dir=/root/af_output \
     $EXTRA_RUN_ALPHAFOLD_FLAGS $USER_SPECIFIED_AF3_OPTIONS \
    || exitcode=$?
else # we must already be in the container
  WORK_DIR_FULL_PATH=`realpath ${WORK_DIR}` # full path to working directory
  pushd /app/alphafold
  python run_alphafold.py \
       --db_dir="${WORK_DIR_FULL_PATH}/public_databases" \
       --model_dir="${MODEL_PARAM_DIR}" \
       --run_data_pipeline=false \
       --run_inference=true \
       --input_dir="${WORK_DIR_FULL_PATH}/af_input" \
       --output_dir="${WORK_DIR_FULL_PATH}/af_output" \
       $EXTRA_RUN_ALPHAFOLD_FLAGS $USER_SPECIFIED_AF3_OPTIONS \
      || exitcode=$?
  popd # back to execution directory
fi

# PROPAGATE ERROR BACK TO HTCONDOR
if (( exitcode != 0 )); then
    printerr "Alphafold3 FAILED with exit code ${exitcode}"
    exit $exitcode
fi

printverbose "Finished running Alphafold3 inference pipeline. Packing up output dir"
shopt -s nullglob # we do not want an empty match below
for output_dir in "${WORK_DIR}/af_output"/*/ ;
do
  output_name_base="$(basename ${output_dir})"
  printverbose "Compressing : $output_name_base"
  tar zcf "${output_name_base}".inference_pipeline.tar.gz -C "${output_dir}" .
done

# clean up
printverbose "Cleaning up working directory"
rm -rf "${WORK_DIR}"
rm -rf .bash_history .bashrc .lesshst .viminfo

printverbose "Done"
