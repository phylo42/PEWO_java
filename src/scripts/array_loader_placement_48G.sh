#!/bin/sh
# Tell the SGE that this is an array job where tasks are run 20 per 20
# Also that it needs a 20G of free ram
# and finally give error/output log destinations
#$ -tc 10
#$ -l "mem_free=50G"
COMMANDFILE=$1"/"$2
CURRENT_COMMAND=$(awk "NR==$SGE_TASK_ID" $COMMANDFILE)
$CURRENT_COMMAND
