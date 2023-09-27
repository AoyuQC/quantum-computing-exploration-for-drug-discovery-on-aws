$Env:CONDA_EXE = "/home/rco/Desktop/projects/qfold/quantum-computing-exploration-for-drug-discovery-on-aws/source/src/notebook/healthcare-and-life-sciences/c-1-protein-folding-quantum-random-walk/quantum-random-walk/psi4/psi4conda/bin/conda"
$Env:_CE_M = ""
$Env:_CE_CONDA = ""
$Env:_CONDA_ROOT = "/home/rco/Desktop/projects/qfold/quantum-computing-exploration-for-drug-discovery-on-aws/source/src/notebook/healthcare-and-life-sciences/c-1-protein-folding-quantum-random-walk/quantum-random-walk/psi4/psi4conda"
$Env:_CONDA_EXE = "/home/rco/Desktop/projects/qfold/quantum-computing-exploration-for-drug-discovery-on-aws/source/src/notebook/healthcare-and-life-sciences/c-1-protein-folding-quantum-random-walk/quantum-random-walk/psi4/psi4conda/bin/conda"
$CondaModuleArgs = @{ChangePs1 = $True}
Import-Module "$Env:_CONDA_ROOT\shell\condabin\Conda.psm1" -ArgumentList $CondaModuleArgs

Remove-Variable CondaModuleArgs