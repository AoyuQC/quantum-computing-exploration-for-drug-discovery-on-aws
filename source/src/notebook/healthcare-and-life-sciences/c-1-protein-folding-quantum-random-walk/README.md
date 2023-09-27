# Introduction

The qfold here is mainly contributed by [Roberto Campos](https://github.com/roberCO) 
based on his work and implementation, 
[QFold: quantum walk and deep learning to solve protein folding](https://iopscience.iop.org/article/10.1088/2058-9565/ac4f2f) 
and 
[QFold Github Repo](https://github.com/roberCO/QFold)

This project is also based on the work MiniFold of Eric Alcaide (https://github.com/hypnopump/MiniFold)

# Instructions

Data file for training the proteins model should be place in /protein-folding-data/training_data/training_70.txt. It can be downloaded from [proteinnet](https://github.com/aqlaboratory/proteinnet) (select CASP7 text-based format). Then extract, decompress and change extension to .txt.

PSI4 library (https://psicode.org/) is executed using a binary in the folder psi4/psi4conda/bin/psi4 and the executable is for python3.10. If you want to execute with a different version, you can download a different executable from https://psicode.org/installs/v18/. PSI4 is the library to extract the atom structure given a sequence of aminoacids.