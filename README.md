# Metalloproteome


## Compiling the code
A precompiled executable can be found in the /bin folder. If you wish to recompile the code, e.g. after making some modifications, you can recompile it using the GNU fortran compiler thus: 
```
gfortran src/analyze_alphafold2_metal_clusters_for_release.f -o analyze_alphafold2_metal_clusters_for_release.exe
```

## Running the code

The code requires that a total of nine(!) command line arguments be provided by the user. They are described below in the order that they are required; the values used in our paper are provided in parentheses. All distances listed are in Angstroms. 
A variety of information is written to the screen (and can be piped to an output file if desired). 
```
analyze_alphafold2_metal_clusters_for_release.exe input.pdb ligand_list_file regionDistanceMinimum regionDistanceMaximum stericClashCutoff disulfideCutoff addBbnAtomCutoff ligandLigandClashCutoff jobID
```

In order, the arguments are:
```
input.pdb                  = A file containing the protein of interest in PDB format
ligand_list_file           = A file listing all of the ligand types to search for, together with their RMSD cutoffs
regionDistanceMinimum      = (0.0) "regions" are merged if any pair of atoms are further than this distance and...
regionDistanceMaximum      = (8.0) ...no greater than this distance...
stericClashCutoff          = (2.0) reject ligand placements if any ligand atom is within this distance of a protein atom
disulfideCutoff            = (2.5) SG atoms in cysteines are considered disulfides if within this distance 
addBbnAtomCutoff           = (0.0) backbone atoms, e.g. carboxyl oxygens, are added to a "region" if within this distance
ligandLigandClashCutoff    = (2.5) reject ligand placements if any ligand atom is within this distance of a previously placed ligand atom
jobID                      = (999) an integer to identify this job if multiple jobs are run in the same folder simultaneously  
```

The format of the `ligand_list_file` should be as follows
```
PDB-FILE-NAME NUMBER-OF-ATOMS-FOR-SUPERPOSITION RMSD-CUTOFF
```
For example:
```
4Fe4S_4CYS.pdb        4  0.5
```

An example run, including making a temporary directory, downloading an AlphaFold2 PDB file, and copying the contents of the `TEMPLATES` folder, would be performed as follows:
```
mkdir TMP
cd TMP
wget https://alphafold.ebi.ac.uk/files/AF-P08201-F1-model_v1.pdb
cp ../TEMPLATES/* .
analyze_alphafold2_metal_clusters_for_release.exe AF-P08201-F1-model_v1.pdb ligand_list_FES_ZINC_RMSD_0.5_12_LIGANDS 0.0 8.0 2.0 2.5 0.0 2.5 999
```

The code generates several output files, most of which were used in debugging the code. The most critical files are as follows:
```
ligand_summary_info_XXXXXX.txt             = Contains information on each ligand type placed, including the coordinating residues, pLDDT scores, RMSDs, etc
ligand_sites_XXXXXX.pdb                    = A PDB file containing the superimposed coordinates of all placed ligands - this can be combined with the ATOM entries from the original AlphaFold2 pdb file
```
The `XXXXXX` in the files above corresponds to the `jobID` input argument.

## Data
The code was run on the 21 proteomes that were predicted with AlphaFold2. The results of this effort can be found in the `DATA` folder, which contains both a single excel file (`ligand_binding_sites_summary.xlsx`) with summary information about each ligand binding site that we identified, and separate folders for each organism containing a `.tar` file combining all of the ligand_sites_XXXXXX.pdb files.

## Citing us
Please cite the following work if you use our code or our data:
