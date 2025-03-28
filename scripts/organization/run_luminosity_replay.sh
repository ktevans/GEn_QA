#!/bin/bash

runnum=$1
first_event=$2
last_event=$3


cd ../../replay

launch_GEN_replay_swif2_test.sh $runnum $first_event $last_event
#echo "launch_GEN_replay_sbatch_test.sh "$runnum" "$first_event" "$last_event

