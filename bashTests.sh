#!/bin/bash

# Run with configuration file and output log.
#python3 $1 >& $2
python3 $1

if [ -s ZZ4lAnalysis.root ]; then
    root -q -b 'AnalysisStep/test/prod/rootFileIntegrity.r("ZZ4lAnalysis.root")'
else
    if [ -f ZZ4lAnalysis.root ]; then
        echo moving empty file
        mv ZZ4lAnalysis.root ZZ4lAnalysis.root.empty
    else
        echo "ZZ4lAnalysis.root does not exist. Writing fake empty file."
        touch ZZ4lAnalysis.root.empty
        exit 1
    fi
fi