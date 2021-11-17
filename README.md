# Metalloproteome


## Compiling the code
Assuming a GNU fortran compiler 
```
gfortran src/analyze_alphafold2_metal_clusters_for_release.f -o analyze_alphafold2_metal_clusters_for_release.exe
```

## Running the code

All geometric values are in angstroms
```
analyze_alphafold2_metal_clusters_for_release.exe in.pdb ligand_list_file regionDistanceMinium regionDistanceMaximum stericClashCutoff disulfideCutoff addBbnAtomCutoff ligandLigandClashCutoff jobID
```

In order the arguements are:
```
in.pdb                     = A file in PDB format with ATOM entries
ligand_list_file           = A file listing all of the ligand PDB files to try fitting with the RMSD cutoff
regionDistanceMinium       = Minimum distance that residues are considered in a "region"
regionDistanceMaximum      = Maximum distance that residues are considered in a "region"
stericClashCutoff          = Distance between two atoms in the template and structure before it is considered a clash
disulfideCutoff            = Distance between two SG atoms in cysteines before being considered a disulfide 
addBbnAtomCutoff           = Cutoff distance for adding backbone atoms, i.e. carboxyl oxygens, to a "region" 
ligandLigandClashCutoff    = Distance between two ligands before they are considered a steric clash
jobID                      = Integer for job identification 
```

The format of the `ligand_list_file` should be as follows
```
PDB-FILE-NAME NUMBER-OF-COORDINATING-RESIDUES RMSD-CUTOFF
```

An example run with a AlphaFold2 PDB file would be as follows
```
mkdir TMP
cd TMP
wget https://alphafold.ebi.ac.uk/files/AF-P08201-F1-model_v1.pdb
cp ../TEMPLATES/* . # directory containing template PDB structures
analyze_alphafold2_metal_clusters_for_release.exe AF-P08201-F1-model_v1.pdb ligand_list_FES_ZINC_RMSD_0.5_12_LIGANDS 0.0 8.0 2.0 2.5 0.0 2.5 1
```

## Data
`analyze_alphafold2_metal_clusters_for_release.exe` was run on the 21 proteomes that were predicted with AlphaFold2. The results of this effort can be found in the `DATA` folder, which contains an excel file (`supplementary_file_ligand_binding_sites.xlsx`) with summary information about each ligand binding site that we identified and folders for each organism containing a `.tar` file of the pdb files containing the ligand coordinates.

## Citing us
Please cite the following work if you use our code or our data:
