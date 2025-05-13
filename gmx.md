# Gromacs tutorial

## Requirements
- Ubuntu 22.04.3 LTS (Download from microsoft store or dual boot) [☞Tutorial](https://youtu.be/RQKp_RA_y2k)
- [Gromacs](https://manual.gromacs.org/documentation/current/download.html) 2023.3 [☞Tutorial](https://youtu.be/JzavO2jt7Pk)
- [Chimera](https://www.cgl.ucsf.edu/chimera/download.html)
- [Discovery Studio](https://discover.3ds.com/discovery-studio-visualizer-download)
- [Swiss Pdb viewer](https://spdbv.unil.ch/download/binaries/SPDBV_4.10_PC.zip)
- [Notepad++](https://notepad-plus-plus.org/downloads/)

## Protein preparation

### Fix Side Chain / Add missing chain
- Open pdb file in discovery studio
- Delete Ligand groups, Water, Hetatm, Unnecessary protein chain
- Save as ```REC.pdb```
- Save the file to (Working directory)

- Go to swiss pdb viewer
Open spdbv.exe file > File > OPen pdb file > select REC.pdb > OPen > Tools > Fix selected side chains > Quick and dirty > File > save > Current layer > select REC.pdb > save and replace
## Ligand preparation

### Add Hydrogen
- Open chimera
- Open ```ligand.sdf```

- Tools>Structure Editing>Add H>OK

- File>Save Mol2>"File name" LIG> Save ```LIG.mol2```
- File should be saved at (Working directory)

### Parameterization
(a)
- Go to swiss param http://old.swissparam.ch/
- Select ```LIG.mol2```
- Submit
- Download ```LIG.zip```
- Unzip 7 files and paste it in the working directory

(b)

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

(optional) if complex file provided, separate ligand using, 
```
grep LIG REC.pdb > LIG.pdb
```
# 1 pdb2gmx
```
gmx pdb2gmx -f REC.pdb -o REC_processed.gro
```
- Select group: 8
- Select group: 1


# 2 Ligand topology
```
gmx editconf -f LIG.pdb -o LIG.gro
```
# 2(a)
- Open ```LIG.gro``` and ```REC_processed.gro``` in npp++

- Copy the coordinates from ```LIG.gro``` (all 1LIG sections) to ```REC_processed.gro``` before the last line and fix the allignment of the first line of 1LIG section

- Observe the line number before the last line and substruct 2 from it and place it to the second line eg: the line before the last line is 2314 ,the line becomes 2312. replace it in the second line.

# 2(b)
- Open ```topol.top``` > go to last line and add
```
LIG		    1
```
- Add ```#include "LIG.itp"``` just below the ```#include "charmm27.ff/forcefield.itp"```

eg: it might be in 22 line
```
; Include forcefield parameters
#include "charmm27.ff/forcefield.itp"
#include "LIG.itp"
```

# 3 Define box and solvate 
```
gmx editconf -f REC_processed.gro -o newbox.gro -bt dodecahedron -d 1.0
```
# 4
```
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```

Check the topol.top file whether the ```SOL     1``` is in the last line or not. if not ,correct it.
 Must check LIG.itp, search for 
```
[ moleculetype ]
; Name nrexcl 
LIG 3
```
it should show `LIG 3`

# 5
```
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
```
# 6
```
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```
- Select group: 15 (SOL)

- Check the ```topol.top``` file whether the CL / Na  is in the last line or not. If not, correct it.

# 7 Energy minimization

```
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
```
# 8
```
gmx mdrun -v -deffnm em
```
(about 700-800 steps)

# 9
```
gmx make_ndx -f LIG.gro -o index_LIG.ndx
```
```
0 & ! a H*
```
```
q
```
# 9(a)

- Go to ```topol.top``` add some line before ```; Include water topology```

eg:
```
; Ligand position restraints
#ifdef POSRES_LIG
#include "posre_LIG.itp"
#endif
```
# 10
```
gmx genrestr -f LIG.gro -n index_LIG.ndx -o posre_LIG.itp -fc 1000 1000 1000
```
- Select group: 3 (System_&_!H*)


# 11 Equilibration step


```
gmx make_ndx -f em.gro -o index.ndx
```
```
1 | 13
```
```
q
```
# 12
```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
```
# 13
```
gmx mdrun -v -deffnm nvt
```
(50,000 steps)

~~gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr~~

# 14

```
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 1
```
# 15
```
gmx mdrun -v -deffnm npt
```
(50,000 steps)

# 16
```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
```
# 17
```
gmx mdrun -v -deffnm md_0_10 -nsteps <steps>
```
(5,000,000 steps = 10,000 ps = 10ns)

# MD Analysis
```
gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
```

- Select group: 1 (Protein)
- Select group: 0 (System)

## RMSD
```
gmx rms -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd.xvg
```
- Select a group: 4 (Backbone)
- Select a group: 4 (Backbone)

## RMSF
```
gmx rmsf -s em.tpr -f md_0_10_center.xtc -n index.ndx -o rmsf.xvg
```
- Select a group: 4 (Backbone)

## Radius of gyration
```
gmx gyrate -s em.tpr -f md_0_10_center.xtc -n index.ndx -o RG.xvg
```
- Select a group: 4 (Backbone)

## H bond
```
gmx hbond -s em.tpr -f md_0_10_center.xtc -n index.ndx -num Hbond.xvg
```
- Select group: 1 (Protein)
- Select group: 13 (LIG)
