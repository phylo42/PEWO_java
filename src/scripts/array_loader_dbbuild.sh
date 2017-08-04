#!/bin/sh
# Tell the SGE that this is an array job where tasks are run 20 per 20
# Also that it needs a 20G of free ram
# and finally give error/output log destinations
# removed: -tc 20  (limit to 20 simultaneous runs)
#$ -l "mem_total=16G"
COMMANDFILE=$1"/dbbuild_commands.list"
CURRENT_COMMAND=$(awk "NR==$SGE_TASK_ID" $COMMANDFILE)
$CURRENT_COMMAND
