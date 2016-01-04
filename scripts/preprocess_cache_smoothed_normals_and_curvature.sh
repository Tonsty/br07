#!/bin/bash

export PATH=/e/br07/build/bin/Release:/e/br07/third_party/Pre-built.2/lib:$PATH

(for i in *.ply; do echo $i; echo ${i%ply}pre; done) | \
	preprocess -curve bbox.dat \
			   -dn 3 \
			   -dc 3
