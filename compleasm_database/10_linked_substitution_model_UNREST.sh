#!/bin/bash

iqtree2 -s concatenated.phy -p rev_dna.best_scheme.nex --model-joint UNREST -B 1000 -T AUTO --prefix nonrev_dna
