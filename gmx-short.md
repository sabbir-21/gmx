# Gromacs tutorial

## Requirements
- Ubuntu 22.04.3 LTS (Download from microsoft store or dual boot) [☞Tutorial](https://youtu.be/RQKp_RA_y2k)
- [Gromacs](https://manual.gromacs.org/documentation/current/download.html) 2023.3 [☞Tutorial](https://youtu.be/JzavO2jt7Pk)
- [Chimera](https://www.cgl.ucsf.edu/chimera/download.html)
- [Discovery Studio](https://discover.3ds.com/discovery-studio-visualizer-download)
- [Swiss Pdb viewer](https://spdbv.unil.ch/download/binaries/SPDBV_4.10_PC.zip)
- [Notepad++](https://notepad-plus-plus.org/downloads/)

## Protein preparation
- Open pdb file in discovery studio
- Delete Ligand groups, Water, Hetatm, Unnecessary protein chain
- Save as ```REC.pdb```
- Save the file to (Working directory)

- Go to swiss pdb viewer
Open spdbv.exe file > File > OPen pdb file > select REC.pdb > OPen > Tools > Fix selected side chains > Quick and dirty > File > save > Current layer > select REC.pdb > save and replace
## Ligand preparation
- Open chimera
- Open ```ligand.sdf```

- Tools>Structure Editing>Add H>OK

- File>Save Mol2>"File name" LIG> Save ```LIG.mol2```
- File should be saved at (Working directory)

### Parameterization
- Go to swiss param http://old.swissparam.ch/
- Select ```LIG.mol2```
- Submit
- Download ```LIG.zip```
- Unzip 7 files and paste it in the working directory
### Error Handling
If error occurs in parametizing in swiss param, 
- Go to http://www.swissparam.ch/ > Draw a molecule using the sketcher > Paste the SMILES copied from pubchem > Add H atom > Parameterize

- Download the zip file > Extract > MATCH/smiles.mol2 > rename LIG.mol2 > Paste in working directory
(if the zip file is not downloading, copy the link and paste in firefox)

- Then again follow the Parameterization step again

## MD Simulation
- Paste the 5 common files ions.mdp, em.mdp, npt.mdp, nvt.mdp, md.mdp. Download from [Here](https://github.com/sabbir-21/gmx/releases/download/v1.0/main.zip)


### Receptor Topology Generation
- Open ubuntu > Go to working directory 
eg: ```cd /mnt/e/working folder name```
to check which files are in a folder give ls command

### Gromacs
```
source /usr/local/gromacs/bin/GMXRC
```
```
printf "8\n1\n" | gmx pdb2gmx -f REC.pdb -o REC_processed.gro
```

### Ligand topology Generation
```
gmx editconf -f LIG.pdb -o LIG.gro
```
- Open ```LIG.gro``` and ```REC_processed.gro``` in npp++

- Copy the coordinates from ```LIG.gro``` (all 1LIG sections) to ```REC_processed.gro``` before the last line and fix the allignment of the first line of 1LIG section

- Observe the line number before the last line and substruct 2 from it and place it to the second line eg: the line before the last line is 2314 ,the line becomes 2312. replace it in the second line.

- Open ```topol.top``` > go to last line and add
```
LIG		    1
```
- Add ```#include "LIG.itp"``` just below the ```#include "charmm27.ff/forcefield.itp"```

eg: it might be in 23 line
```
; Include forcefield parameters
#include "charmm27.ff/forcefield.itp"
#include "LIG.itp"
```
### Define box and solvate 
```
gmx editconf -f REC_processed.gro -o newbox.gro -bt dodecahedron -d 1.0
```
```
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```
```
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
```
Check the topol.top file whether the ```SOL     1``` is in the last line or not. if not ,correct it.
```
printf "SOL\n" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

- Check the ```topol.top``` file whether the CL  is in the last line or not. If not, correct it.

### Energy minimization
```
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
```
```
gmx mdrun -v -deffnm em
```
(about 700-800 steps)
```
printf "0 & ! a H*\nq\n" | gmx make_ndx -f LIG.gro -o index_LIG.ndx
```

- Go to ```topol.top``` add some line before ```; Include water topology```

eg:
```
; Ligand position restraints
#ifdef POSRES_LIG
#include "posre_LIG.itp"
#endif
```
```
printf "System_&_!H*\n" | gmx genrestr -f LIG.gro -n index_LIG.ndx -o posre_LIG.itp -fc 1000 1000 1000
```


### Equilibration step
```
printf "1 | 13\nq\n" | gmx make_ndx -f em.gro -o index.ndx
```

```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
```
```
gmx mdrun -v -deffnm nvt
```
(50,000 steps)

~~gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr~~

```
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 1
```
```
gmx mdrun -v -deffnm npt
```
(50,000 steps)

```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
```
```
gmx mdrun -v -deffnm md_0_10
```
(5,000,000 steps = 10,000 ps = 10ns)

## MD Analysis
```
printf "Protein\nSystem\n" | gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
```

### RMSD
```
printf "Backbone\nBackbone\n" | gmx rms -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd.xvg
```

### RMSF
```
printf "Backbone\n" | gmx rmsf -s em.tpr -f md_0_10_center.xtc -n index.ndx -o rmsf.xvg
```

### Radius of gyration
```
printf "Backbone\n" | gmx gyrate -s em.tpr -f md_0_10_center.xtc -n index.ndx -o RG.xvg
```

### H bond
```
printf "Protein\nLIG\n" | gmx hbond -s em.tpr -f md_0_10_center.xtc -n index.ndx -num Hbond.xvg
```
