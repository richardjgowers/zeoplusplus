#!/bin/bash

for i in Molecules/*.xyz
do
	./single_mol_analysis $i
done



