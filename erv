#!/usr/bin/env bash
set -eo pipefail
module purge

# Author: Vishal Koparde, PhD
# CCBR, NCI
# (c) 2023
#
# wrapper script to run the snakemake pipeline
# a) on an interactive node (runlocal) OR
# b) submit to the slurm load scheduler (run)
#

# setting python and snakemake versions
PYTHONVERSION="3.10"
SNAKEMAKEVERSION="7.32.3"

##########################################################################################
# functions
##########################################################################################

function get_git_commitid_tag() {
  # This function gets the latest git commit id and tag
  # @Input:
  #   @param1: PIPELINE_HOME folder which is a initialized git repo folder path
  # @Output:
  #   @param1: tab-delimited commit id and tag

    cd $1
    gid=$(git rev-parse HEAD)
    tag=$(git describe --tags $gid 2>/dev/null)
    echo -ne "$gid\t$tag"
}

##########################################################################################

function print_banner() {
  versionnumber=$1
cat << EOF

╭━━━╮╱╱╱╱╭╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╭━━━╮╱╱╭╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╭━━━╮╱╱╱╱╱╱╭╮
┃╭━━╯╱╱╱╱┃┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃╭━╮┃╱╭╯╰╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃╭━╮┃╱╱╱╱╱╱┃┃
┃╰━━┳━╮╭━╯┣━━┳━━┳━━┳━╮╭━━┳━━┳╮╭┳━━╮┃╰━╯┣━┻╮╭╋━┳━━┳╮╭┳┳━┳╮╭┳━━╮┃╰━╯┣┳━━┳━━┫┃╭┳━╮╭━━╮
┃╭━━┫╭╮┫╭╮┃╭╮┃╭╮┃┃━┫╭╮┫┃━┫╭╮┃┃┃┃━━┫┃╭╮╭┫┃━┫┃┃╭┫╭╮┃╰╯┣┫╭┫┃┃┃━━┫┃╭━━╋┫╭╮┃┃━┫┃┣┫╭╮┫┃━┫
┃╰━━┫┃┃┃╰╯┃╰╯┃╰╯┃┃━┫┃┃┃┃━┫╰╯┃╰╯┣━━┃┃┃┃╰┫┃━┫╰┫┃┃╰╯┣╮╭┫┃┃┃╰╯┣━━┃┃┃╱╱┃┃╰╯┃┃━┫╰┫┃┃┃┃┃━┫
╰━━━┻╯╰┻━━┻━━┻━╮┣━━┻╯╰┻━━┻━━┻━━┻━━╯╰╯╰━┻━━┻━┻╯╰━━╯╰╯╰┻╯╰━━┻━━╯╰╯╱╱╰┫╭━┻━━┻━┻┻╯╰┻━━╯
╱╱╱╱╱╱╱╱╱╱╱╱╱╭━╯┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃┃
╱╱╱╱╱╱╱╱╱╱╱╱╱╰━━╯╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╰╯
... v${versionnumber}
EOF
# generated using https://fsymbols.com/generators/carty/
}

##########################################################################################

function _set_rand_str() {
  x=$(mktemp)
  rm -rf $x
  RAND_STR=$(echo $x|awk -F"." '{print $NF}')
}

##########################################################################################

function _set_config() {
  f="${PIPELINE_HOME}/config/config.yaml"
  echo "Copying essential file: $f"
  fbn=$(basename $f)
sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" \
  -e "s/WORKDIR/${WORKDIR//\//\\/}/g" \
  -e "s/GENOME/${GENOME}/g" $f > $WORKDIR/$fbn
}

##########################################################################################

function _cp_other_resources() {

for f in ${PIPELINE_HOME}/resources/cluster.yaml ${PIPELINE_HOME}/resources/tools.yaml
do
  echo "Copying essential file: $f"
  fbn=$(basename $f)
  cp $f $WORKDIR/$fbn
done

# copy essential folders
for f in $ESSENTIAL_FOLDERS;do
  # rsync -az --progress ${PIPELINE_HOME}/$f $WORKDIR/
  cp -rv ${PIPELINE_HOME}/$f ${WORKDIR}/
done

}


##########################################################################################

function _exe_in_path() {

  name_of_exe=$1
  path_to_exe=$(which $name_of_exe 2>/dev/null)
 if [ ! -x "$path_to_exe" ] ; then
    err $path_to_exe NOT FOUND!
 fi

}

##########################################################################################

function _get_file_modtime() {

# get the modification time for a file

  filename=$1
  modtime=$(stat $filename|grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
  echo $modtime

}

##########################################################################################
# initial setup
##########################################################################################

SCRIPTNAME="$0"
SCRIPTDIRNAME=$(readlink -f $(dirname $0))
SCRIPTBASENAME=$(readlink -f $(basename $0))

GENOME="hg38"
SUPPORTED_GENOMES="hg38 mm10"

# essential files
# these are relative to the workflows' base folder
# these are copied into the WORKDIR
# dealt with in init function
# ESSENTIAL_FILES="config/config.yaml config/samples.tsv config/fqscreen_config.conf resources/cluster.yaml resources/tools.json"
ESSENTIAL_FOLDERS="workflow/scripts"

## setting PIPELINE_HOME
PIPELINE_HOME=$(readlink -f $(dirname "$0"))

# set snakefile
SNAKEFILE="${PIPELINE_HOME}/workflow/Snakefile"

# get github commit tag
GIT_COMMIT_TAG=$(get_git_commitid_tag $PIPELINE_HOME)

VERSIONFILE="${PIPELINE_HOME}/VERSION"
VERSION=$(head -n1 $VERSIONFILE|awk '{print $1}')

##########################################################################################
# USAGE
##########################################################################################

function usage() { cat << EOF

##########################################################################################

Welcome to
EOF
print_banner $VERSION
cat << EOF

##########################################################################################

This pipeline was built by CCBR (https://bioinformatics.ccr.cancer.gov/ccbr)
Please contact Vishal Koparde for comments/questions (vishal.koparde@nih.gov)

##########################################################################################

Here is a list of genome supported by this pipeline:

  * hg38          [Human]
  * mm10          [Mouse]

USAGE:
  /path/to/erv -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>

Required Arguments:
1.  WORKDIR     : [Type: String]: Absolute or relative path to the output folder with write permissions.

2.  RUNMODE     : [Type: String] Valid options:
    * init      : initialize workdir
    * dryrun    : dry run snakemake to generate DAG
    * run       : run with slurm
    * runlocal  : run without submitting to sbatch
    ADVANCED RUNMODES (use with caution!!)
    * unlock    : unlock WORKDIR if locked by snakemake NEVER UNLOCK WORKDIR WHERE PIPELINE IS CURRENTLY RUNNING!
    * reconfig  : recreate config file in WORKDIR (debugging option) EDITS TO config.yaml WILL BE LOST!
    * recopy    : recreate tools.yaml, cluster.yaml and scriptsdir in WORKDIR (debugging option) EDITS TO these files WILL BE LOST!
    * reset     : DELETE workdir dir and re-init it (debugging option) EDITS TO ALL FILES IN WORKDIR WILL BE LOST!
    * local     : same as runlocal

Optional Arguments:

--genome|-g     : genome eg. hg38(default) or mm10
--manifest|-s   : absolute path to samples.tsv. This will be copied to output folder  (--runmode=init only)
--help|-h       : print this help

Example commands:
  ${SCRIPTNAME} -w=/my/output/folder -m=init [ -g="mm10" -s="/path/to/sample.tsv" ]
  ${SCRIPTNAME} -w=/my/output/folder -m=dryrun
  ${SCRIPTNAME} -w=/my/output/folder -m=run

##########################################################################################

VersionInfo:
  python          : $PYTHONVERSION
  snakemake       : $SNAKEMAKEVERSION
  pipeline_home   : $PIPELINE_HOME
  git commit/tag  : $GIT_COMMIT_TAG
  pipeline_version: v${VERSION}

##########################################################################################

EOF
}

##########################################################################################
# ERR
##########################################################################################

function err() { usage && cat <<< "
#
# ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR
#
  $@
#
# ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR
#
" && exit 1 1>&2;
}

##########################################################################################
# INIT
##########################################################################################

function init() {

# This function initializes the workdir by:
# 1. creating the working dir
# 2. copying essential files like config.yaml and samples.tsv into the workdir
# 3. setting up logs and stats folders

print_banner $VERSION

# create output folder
if [ -d $WORKDIR ];then err "Folder $WORKDIR already exists!"; fi
mkdir -p $WORKDIR

# copy essential files
# for f in $ESSENTIAL_FILES;do
_set_config
echo "Set config.yaml"

cp $MANIFEST $WORKDIR/samples.tsv

_cp_other_resources

cd ${WORKDIR}
curl https://raw.githubusercontent.com/CCBR/Tools/master/Biowulf/jobby -O
curl https://github.com/CCBR/Tools/blob/master/Biowulf/run_jobby_on_snakemake_log -O

chmod a+x ${WORKDIR}/jobby
chmod a+x ${WORKDIR}/run_jobby_on_snakemake_log

#create log and stats folders
if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi

cat << EOF
Done Initializing   : $WORKDIR
You can now edit    : $WORKDIR/config.yaml and
                      $WORKDIR/samples.tsv
EOF

}

##########################################################################################
# CHECK ESSENTIAL FILES
##########################################################################################

function check_essential_files() {

# Checks if files essential to start running the pipeline exist in the workdir

  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!"; fi
  for f in config.yaml samples.tsv cluster.yaml tools.yaml; do
    if [ ! -f $WORKDIR/$f ]; then err "Error: '${f}' file not found in workdir ... initialize first!";fi
  done

}

##########################################################################################
# RECONFIG ... recreate config.yaml and overwrite old version
##########################################################################################

function reconfig(){

  # Rebuild config file and replace the config.yaml in the WORKDIR
  # this is only for dev purposes when new key-value pairs are being
  # updated in the config file in PIPELINE_HOME

  _set_config
  echo "$WORKDIR/config.yaml has been updated!"

}

##########################################################################################
# RECOPY ... recopy tools.yaml, cluster.yaml and scripts folder
##########################################################################################

function recopy() {
  _cp_other_resources
  echo "Done!"
}

##########################################################################################
# RUNCHECK ... check essential files and load required packages
##########################################################################################

function runcheck(){
  # Check "job-essential" files and load required modules
  check_essential_files
  MODULE_STR="module purge && module load python/$PYTHONVERSION snakemake/$SNAKEMAKEVERSION"

}

##########################################################################################
# DRYRUN ... also run automatically before actual run
##########################################################################################

function dryrun() {
  # check essential files, load modules and do Dry-run
  runcheck
  # can add controlcheck if needed

  if [ ! -d ${WORKDIR}/logs/dryrun/ ]; then mkdir ${WORKDIR}/logs/dryrun/; fi

  if [ -f ${WORKDIR}/dryrun.log ]; then
    modtime=$(stat ${WORKDIR}/dryrun.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/dryrun.log ${WORKDIR}/logs/dryrun/dryrun.${modtime}.log
  fi

  run "--dry-run"

}

##########################################################################################
# UNLOCK
##########################################################################################

function unlock() {
  # check essential files, load modules and
  # unlock the workdir if previous snakemake run ended abruptly

  runcheck
  run "--unlock"

}

##########################################################################################
# DAG
##########################################################################################

function dag() {
  runcheck
  snakemake -s $SNAKEFILE --configfile ${WORKDIR}/config.yaml --forceall --dag |dot -Teps > ${WORKDIR}/dag.eps
}

##########################################################################################
# RUNLOCAL ... run directly on local interactive node ... no submission to SLURM
##########################################################################################

function runlocal() {
# If the pipeline is fired up on an interactive node (with sinteractive), this function runs the pipeline

  runcheck
  if [ "$SLURM_JOB_ID" == "" ];then err "runlocal can only be done on an interactive node"; exit 1; fi
  run "--dry-run" && echo "Dry-run was successful .... starting local execution" && \
  run "local"
}

function runslurm() {

  # Submit the execution of the pipeline to the biowulf job scheduler (slurm)
  runcheck
  run "slurm"

}

##########################################################################################
# CREATE RUNINFO ... create runinfo.yaml in workdir
##########################################################################################

function create_runinfo() {

# Create a runinfo.yaml file in the WORKDIR

  if [ -f ${WORKDIR}/runinfo.yaml ];then
    modtime=$(_get_file_modtime ${WORKDIR}/runinfo.yaml)
    mv ${WORKDIR}/runinfo.yaml ${WORKDIR}/runinfo.yaml.${modtime}
  fi

  echo "Pipeline: ERVPipeline" > ${WORKDIR}/runinfo.yaml
  echo "Pipeline Version: v${VERSION}" >> ${WORKDIR}/runinfo.yaml
  echo "Pipeline Dir: $PIPELINE_HOME" >> ${WORKDIR}/runinfo.yaml
  echo "Git Commit/Tag: $GIT_COMMIT_TAG" >> ${WORKDIR}/runinfo.yaml
  echo "Work Dir: $WORKDIR" >> ${WORKDIR}/runinfo.yaml
  userlogin=$(whoami)
  IFS=: read user x uid gid gecos hm sh < <( getent passwd $userlogin )
  # fullname=${gecos%%,*}
  fullname=${gecos}
  # username=$(finger $userlogin|grep ^Login|awk -F"Name: " '{print $2}')
  echo "Username: $userlogin" >> ${WORKDIR}/runinfo.yaml
  echo "Full Name: $fullname" >> ${WORKDIR}/runinfo.yaml
  echo "UID: $uid" >> ${WORKDIR}/runinfo.yaml
  echo "GID: $gid" >> ${WORKDIR}/runinfo.yaml
  group=$(getent group $gid | cut -d: -f1)
  echo "Group: $group" >> ${WORKDIR}/runinfo.yaml
  d=$(date)
  echo "Date/Time: $d" >> ${WORKDIR}/runinfo.yaml

}

##########################################################################################
# PRERUN CLEANUP ... get ready to run .. park old logs/stats etc.
##########################################################################################

function preruncleanup() {

  # Cleanup function to rename/move files related to older runs to prevent overwriting them.
  echo "Running..."

  # check initialization
  check_essential_files

  cd $WORKDIR
  ## Archive previous run files
  if [ -f ${WORKDIR}/snakemake.log ];then
    modtime=$(_get_file_modtime ${WORKDIR}/snakemake.log)
    mv ${WORKDIR}/snakemake.log ${WORKDIR}/stats/snakemake.${modtime}.log
    if [ -f ${WORKDIR}/snakemake.log.jobby ];then
      mv ${WORKDIR}/snakemake.log.jobby ${WORKDIR}/stats/snakemake.${modtime}.log.jobby
    fi
    if [ -f ${WORKDIR}/snakemake.log.jobby ];then
      mv ${WORKDIR}/snakemake.log.jobby_short ${WORKDIR}/stats/snakemake.${modtime}.log.jobby_short
    fi
    if [ -f ${WORKDIR}/snakemake.stats ];then
      mv ${WORKDIR}/snakemake.stats ${WORKDIR}/stats/snakemake.${modtime}.stats
    fi
  fi
  nslurmouts=$(find ${WORKDIR} -maxdepth 1 -name "slurm-*.out" |wc -l)
  if [ "$nslurmouts" != "0" ];then
    for f in $(ls ${WORKDIR}/slurm-*.out);do mv ${f} ${WORKDIR}/logs/;done
  fi

  create_runinfo

}

function run() {
  # RUN function
  # argument1 can be:
  # 1. local or
  # 2. dryrun or
  # 3. unlock or
  # 4. slurm

  print_banner $VERSION

  if [ "$1" == "local" ];then

  preruncleanup
  _set_rand_str

  cat > ${WORKDIR}/.${RAND_STR} << EOF
#/bin/bash
set -exo pipefail

$MODULE_STR

snakemake -s $SNAKEFILE \
--directory $WORKDIR \
--printshellcmds \
--use-envmodules \
--latency-wait 120 \
--configfile ${CONFIGFILE} \
--cores all \
--restart-times ${RETRIES} \
${RERUNTRIGGERS} \
--keep-going \
--stats ${WORKDIR}/snakemake.stats \
2>&1|tee ${WORKDIR}/snakemake.log

# if [ "$?" -eq "0" ];then
#   snakemake -s $SNAKEFILE \
#   --report ${WORKDIR}/runlocal_snakemake_report.html \
#   --directory $WORKDIR \
#   --configfile ${WORKDIR}/config.yaml
# fi
EOF

bash ${WORKDIR}/.${RAND_STR}

  elif [ "$1" == "slurm" ];then

    preruncleanup
    # if QOS is other than "global" and is supplied in the cluster.yaml file then add " --qos={cluster.qos}" to the
    # snakemake command below
  cat > ${WORKDIR}/submit_script.sbatch << EOF
#!/bin/bash
#SBATCH --job-name="ERV"
#SBATCH --mem=40g
#SBATCH --partition="norm"
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2
#SBATCH --gres=lscratch:48

$MODULE_STR

cd \$SLURM_SUBMIT_DIR

snakemake -s $SNAKEFILE \
--directory $WORKDIR \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--configfile $CONFIGFILE \
--cluster-config $CLUSTERFILE \
--cluster "$CLUSTER_SBATCH_CMD" \
-j 500 \
--rerun-incomplete \
--restart-times ${RETRIES} \
${RERUNTRIGGERS} \
--keep-going \
--stats ${WORKDIR}/snakemake.stats \
2>&1|tee ${WORKDIR}/snakemake.log

if [ "\$?" -eq "0" ];then
  snakemake -s $SNAKEFILE \
  --directory $WORKDIR \
  --report ${WORKDIR}/runslurm_snakemake_report.html \
  --configfile $CONFIGFILE
fi

EOF

  sbatch ${WORKDIR}/submit_script.sbatch

  else # for unlock and dryrun
  _set_rand_str

  cat > ${WORKDIR}/.${RAND_STR} << EOF
#/bin/bash
set -exo pipefail
$MODULE_STR

snakemake $1 -s $SNAKEFILE \
--directory $WORKDIR \
--printshellcmds \
--latency-wait 120 \
--configfile $CONFIGFILE \
--cluster-config $CLUSTERFILE \
--cluster "$CLUSTER_SBATCH_CMD" \
${RERUNTRIGGERS} \
-j 500 \
--rerun-incomplete \
--keep-going \
--reason \
--stats ${WORKDIR}/snakemake.stats
EOF

if [ "$1" == "--dry-run" ];then
bash ${WORKDIR}/.${RAND_STR} > ${WORKDIR}/dryrun.log && rm -f ${WORKDIR}/.${RAND_STR}
less ${WORKDIR}/dryrun.log
else # unlock
bash ${WORKDIR}/.${RAND_STR} && rm -f ${WORKDIR}/.${RAND_STR}
fi
  fi

}

##########################################################################################
# RESET ... delete workdir and then initialize
##########################################################################################

function reset() {
  # Delete the workdir and re-initialize it
  print_banner $VERSION
  echo "Working Dir: $WORKDIR"
  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!";fi
  echo "Deleting $WORKDIR"
  rm -rf $WORKDIR
  echo "Re-Initializing $WORKDIR"
  init

}

function main(){
  # Main function which parses all arguments

  if [ $# -eq 0 ]; then usage && exit 1; fi

  for i in "$@"; do
  case $i in
      -m=*|--runmode=*)
        RUNMODE="${i#*=}"
      ;;
      -w=*|--workdir=*)
        WORKDIR="${i#*=}"
      ;;
      -s=*|--manifest=*)
        MANIFEST="${i#*=}"
        if [ ! -f $MANIFEST ];then err "File $MANIFEST does NOT exist!";fi
      ;;
      -g=*|--genome=*)
        GENOME="${i#*=}"
        found=0
        for g in $SUPPORTED_GENOMES;do
          if [[ "$GENOME" == "$g" ]];then
            found=1
            break
          fi
        done
        if [[ "$found" == "0" ]];then
          err "$GENOME is not supported by this pipeline; Supported genomes are: $SUPPORTED_GENOMES"
          exit 1
        fi
      ;;
      --version)
	cat ${SCRIPTDIRNAME}/VERSION && exit 0;
      ;;
      -h|--help)
        usage && exit 0;;
      *)
        err "Unknown argument:     $i!"    # unknown option
      ;;
  esac
  done
  WORKDIR=$(readlink -f "$WORKDIR")
  if [[ -z $MANIFEST ]];then
    MANIFEST=${PIPELINE_HOME}/config/samples.tsv
  fi
  echo "Working Dir       : $WORKDIR"
  echo "Samples Manifest  : $MANIFEST"
  CLUSTER_SBATCH_CMD="sbatch --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} --gres {cluster.gres}"
  RETRIES="2"
  RERUNTRIGGERS="--rerun-triggers mtime"
  # RERUNTRIGGERS="--rerun-triggers input"
  CONFIGFILE="${WORKDIR}/config.yaml"
  CLUSTERFILE="${WORKDIR}/cluster.yaml"
  JOBBY="${WORKDIR}/jobby"
  JOBBY2="${WORKDIR}/run_jobby_on_snakemake_log"

  case $RUNMODE in
    init) init && exit 0;;
    dryrun) dryrun && exit 0;;
    unlock) unlock && exit 0;;
    run) runslurm && exit 0;;
    runlocal) runlocal && exit 0;;
    reset) reset && exit 0;;
    dry) dryrun && exit 0;;                      # hidden option
    local) runlocal && exit 0;;                  # hidden option
    reconfig) reconfig && exit 0;;               # hidden option for debugging
    recopy) recopy && exit 0;;                   # hidden option
    *) err "Unknown RUNMODE \"$RUNMODE\"";;
  esac


}

# call the main function

main "$@"
