#!/bin/bash

export PATH=/e/br07/build/bin/Release:/e/br07/third_party/Pre-built.2/lib:$PATH

ls *.ply | preprocess -bbox bbox.dat