#!/bin/bash

iqtree2 -s concatenated.phy -p concatenated.phy.nex -B 1000 -T AUTO --prefix rev_dna
