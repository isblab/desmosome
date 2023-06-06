#!/bin/bash
SECONDS=0
echo Doing Progress check for $1 \($4 output logs\) at intervals of $2-seconds \($3 iterations\)
echo At $SECONDS seconds since starting:
for ((j = 0; j < $4; j++))
do
	if test -f ./$1/output_log_$j.log; then
		val=$(cat ./$1/output_log_$j.log | grep -E "frame [0-9]+ score" | wc -l)
		echo For output log $j frames: $val
	fi
done
echo ------------------------------------------------
for ((i = 0; i < $3; i++))
do
	sleep $2s
	echo At $SECONDS second since starting:
	for ((j = 0; j < $4; j++))
	do
		if test -f ./$1/output_log_$j.log; then
			val=$(cat ./$1/output_log_$j.log | grep -E "frame [0-9]+ score" | wc -l)
			echo For output log $j frames: $val
		fi
	done
	echo ------------------------------------------------
done

