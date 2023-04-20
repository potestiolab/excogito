#!/bin/bash
# test suite for build/excogito
BUILDDIR=../../build

declare -a HashTable=( "test0:$BUILDDIR/excogito optimize -p ../files/parameters/parameters_optimize_6d93_N31_small.ini -t ../files/trajectories/6d93_100frames.xyz -e ../files/energies/6d93_energies_100frames.txt -c 6d93"
"test1:$BUILDDIR/excogito random -p ../files/parameters/parameters_random_6d93_N31_small.ini -t ../files/trajectories/6d93_100frames.xyz -e ../files/energies/6d93_energies_100frames.txt -c 6d93"
"test2:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt"
"test3:$BUILDDIR/excogito optimize -p ../files/parameters/parameters_optimize_6d93_N31_notemp.ini -t ../files/trajectories/6d93_100frames.xyz -e ../files/energies/6d93_energies_100frames.txt -c 6d93"
"test4:$BUILDDIR/excogito optimizer -p ../files/parameters/parameters_test_4.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93"
"test5:$BUILDDIR/excogito optimize -p ../files/parameters/parameters_test_5.ini -t ../files/trajectories/6d93_100frames.xyz -e ../files/energies/6d93_energies_100frames.txt -c 6d93"
"test6:$BUILDDIR/excogito optimize -p ../files/parameters/parameters_test_5.ini -t ../files/trajectories/6d93_100frames_missing.xyz -e ../files/energies/6d93_energies_100frames.txt -c 6d93"
"test7:$BUILDDIR/excogito optimize -p ../files/parameters/parameters_test_5.ini -t ../files/trajectories/6d93_100frames_missing_last.xyz -e ../files/energies/6d93_energies_100frames.txt -c 6d93"
"test8:$BUILDDIR/excogito optimize -p ../files/parameters/parameters_test_5.ini -t ../files/trajectories/6d93_100frames.xyz -e ../files/energies/6d93_energies_100frames_missing.txt -c 6d93"
"test9:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e  ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_mapping_shorter.txt"
"test10:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_mapping_longer.txt"
"test11:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_mapping_wrong.txt"
"test12:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_mapping_wrong_2.txt"
"test13:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_mapping_duplicate.txt"
"test14:$BUILDDIR/excogito optimize -p ../files/parameters/parameters_optimize_6d93_N31_small_stringError.ini -t ../files/trajectories/6d93_100frames.xyz -e ../files/energies/6d93_energies_100frames.txt -c 6d93"
"test15:$BUILDDIR/excogito measure -p ../files//parameters/parameters_measure_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames_IntValue.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt"
"test16:$BUILDDIR/excogito measure -p ../files//parameters/parameters_measure_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames_MoreColumns.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt"
"test17:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31_10frames.ini -t ../files/trajectories/6d93_10frames_Error_nAT.xyz -e ../files/energies/6d93_energies_10frames.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt"
"test18:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31_10frames.ini -t ../files/trajectories/6d93_10frames_StringSecondCol.xyz -e ../files/energies/6d93_energies_10frames.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt"
"test19:$BUILDDIR/excogito norm -p ../files/parameters/parameters_norm_4ake_N214.ini -t ../files/trajectories/4ake_coord_xyz.xyz -e ../files/energies/4ake_energies_1frames.txt -c 4ake -m ../files/mappings/adenylate_ca_214.txt"
"test20:$BUILDDIR/excogito cosine -p ../files/parameters/parameters_cosine_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt --m2 ../files/mappings/tamapin_ca_mapping.txt"
"test21:$BUILDDIR/excogito norm -p ../files/parameters/parameters_norm_6d93_N31.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt"
"test22:$BUILDDIR/excogito distance -p ../files/parameters/parameters_distance_4ake_N214.ini -t ../files/trajectories/4ake_coord_xyz.xyz -e ../files/energies/4ake_energies_1frames.txt -c 4ake -x ../files/mappings/4ake_pdb_mappings.txt"
"test23:$BUILDDIR/excogito measure -p ../files/parameters/parameters_measure_6d93_N31_fastclust.ini -t ../files/trajectories/6d93_1000frames.xyz -e ../files/energies/6d93_energies_1000frames.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt"
"test24:$BUILDDIR/excogito optimize_kl -p ../files/parameters/parameters_optimizekl_6d93_N31_notemp.ini -t ../files/trajectories/6d93_100frames.xyz -r ../files/energies/6d93_probs_100frames.txt -c 6d93"
"test25:$BUILDDIR/excogito measure_kl -p ../files/parameters/parameters_measurekl_6d93_N31.ini -t ../files/trajectories/6d93_100frames.xyz -r ../files/energies/6d93_probs_100frames.txt -c 6d93 -m ../files/mappings/tamapin_ca_mapping.txt"
"test26:$BUILDDIR/excogito random_kl -p ../files/parameters/parameters_randomkl_6d93_N31.ini -t ../files/trajectories/6d93_100frames.xyz -r ../files/energies/6d93_probs_100frames.txt -c 6d93")

for test_line in "${HashTable[@]}"; do
test_id="${test_line%%:*}"
echo $test_id
if [ -d "$test_id" ];
then
echo "directory ${test_id} exists"
rm -r $test_id/*
fi
done

for test_line in "${HashTable[@]}"; do
test_id="${test_line%%:*}"
command="${test_line##*:}"
# check if directory exists
dir=${test_id}
[ ! -d "$dir" ] && mkdir $dir
cd ${test_id}
$command > log_${test_id}.log
cd ..
done
