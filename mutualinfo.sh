#!/bin/bash
test ${#} -ne 3 && (echo "Insufficient number of arguments: 3 expected, got ${#}"; exit 1)

xmin=${1}
xmax=${2}
xdel=${3}
xmaxby2=$((${xmax}/2))

bindis=../disorder
binnor=../normal
binrep=../replica
binmut=../mutualinfo

dirname=`date +%F-%H-%M-%S`
echo
echo "Mutual-Info calculator"
echo "============================================================"
if [ -d ${dirname} ]
then
  echo "Directory with same exists, moving existing dir to ${dirname}.old"
  mv ${dirname} ${dirname}.old
fi

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
echo "Running ${binrep} from ${xmin}--${xmaxby2} with ${xdel} increments..."
${binrep} ${xmin} ${xmaxby2} ${xdel}
echo "${binrep} done."

echo
echo "Computing mutual-info..."
${binmut} ${xmin} ${xmaxby2} ${xdel}
echo "${binmut} done."

# TO GUARANTEE SAFETY WHEN THE SCRIPT IS CALLED SEVERAL TIMES, E.G., IN A LOOP, MAKE SURE THE DIR NAMES DIFFER BY AT LEAST 2s
sleep 2

echo
echo "All done. Bye!"
echo "============================================================"

popd
