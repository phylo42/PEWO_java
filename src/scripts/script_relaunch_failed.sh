#!/usr/bin/bash

echo "must be launched from Dx directory"
echo "arguments: job_id associated_k"

job_id=$1
files="ls -lh qsub_logs/*e`echo $job_id`* | grep -v ' 0 ' | cut -d '.' -f4 | tr -d 'e'"
ids=`eval "$files"`
for i in $ids
do
	echo "array_id to relaunch:"$i
	awk "NR==$i" "k"$2"_placement_commands.list" >> "k"$2"_to_relaunch"
done

if [ ! -f "k"$2"_placement_commands.list_backup" ]
then
    echo "backup file created"
    cp "k"$2"_placement_commands.list" "k"$2"_placement_commands.list_backup"
fi

cp "k"$2"_to_relaunch" "k"$2"_placement_commands.list"

exit 0
