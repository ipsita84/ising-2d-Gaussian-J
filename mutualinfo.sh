#!/bin/bash
test ${#} -ne 3 && (echo "Insufficient number of arguments: 3 expected, got ${#}"; exit 1)

xmin=${1}
xmax=${2}
xdel=${3}
xmaxby2=$(bc <<< "scale=3; ${xmax}/2")

bindis=../disorder
binnor=../normal
binrepA=../replicaA
binrepB=../replicaB
binmut=../mutualinfo

dirname=`date +%F-%H-%M-%S`
echo
echo "Mutual-Info calculator"
echo "============================================================"
dirtry=1
while [ -d ${dirname} ]; do
  dirname=${dirname}-${dirtry}
  echo "Directory with same name exists, using different dir name ${dirname}"
done

echo ""
echo "Generating executables..."
make all
echo "All executables generated."

mkdir ${dirname}
pushd ${dirname}

echo
echo "Generating disorder realisation: ${bindis}..."
${bindis}
echo "${bindis} done."

echo
echo "Running ${binnor} from ${xmin}--${xmax} with ${xdel} increments..."
${binnor} ${xmin} ${xmax} ${xdel}
echo "${binnor} done."

echo
echo "Running ${binrepA} from ${xmin}--${xmaxby2} with ${xdel} increments..."
${binrep} ${xmin} ${xmaxby2} ${xdel}
echo "${binrep} done."

echo
echo "Running ${binrepB} from ${xmin}--${xmaxby2} with ${xdel} increments..."
${binrep} ${xmin} ${xmaxby2} ${xdel}
echo "${binrep} done."

echo
echo "Computing mutual-info..."
${binmut} ${xmin} ${xmaxby2} ${xdel}
echo "${binmut} done."

sleep 2

echo
echo "All done. Bye!"
echo "============================================================"

popd
