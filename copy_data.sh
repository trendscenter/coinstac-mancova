#!/usr/bin/env bash

mkdir -p test/input/remote/simulatorRun
echo "Created test/remote/simulatorRun"

cp local_data/* test/input/remote/simulatorRun
cp local_data/* test/input/local0/simulatorRun
cp local_data/* test/input/local1/simulatorRun

echo "Copied data to remote, local0, and local1"
