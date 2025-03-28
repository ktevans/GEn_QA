#!/bin/bash

kinematic=$1
type=$2
runnum=$3

outfile=/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/$kinematic/$type
workflow=jeffas_${kinematic}_$type

cd ../../replay

swif2 create $workflow

launch_GEN_replay_swif2.sh $runnum 30 1 $outfile $workflow
#echo "launch_GEN_replay_swif2.sh "$kinematic" "$runnum" 0 30 1 "$outfile" "$workflow

swif2 run $workflow
