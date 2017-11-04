#!/bin/bash
cd /usr/local/geant4/applications/share

echo "Running eff10_mod version: "`eff10_mod --version`" in: "`pwd`
echo "Geant4 Version: "`geant4-config --version`
echo "Macro file loaded: "${1}
echo "Flags: " ${2}
date

echo "Starting, output in $PWD/cout.log and $PWD/cerr.log"
eff10_mod /macros/${1} ${2} > >(tee $PWD/.cerr.log) 1> >(tee $PWD/cout.log)
set -e

echo "Done, copying output to /output/out.tgz"
mkdir -p /output

set +e
files=`ls *.csv *.root 2> /dev/null` 
set -e
tar czf /output/out.tgz ${files}
date
echo "All done, bye"
