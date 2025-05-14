# GMX_MMPBSA tutorial

## Requirements

md_0_100.tpr

md_0_10_center.xtc

index.ndx

topol.top

LIG.itp

topol_Protein_chain_A.itp

## MMPBSA input creation

For All calculation
```
gmx_MMPBSA --create_input
```

For GB,
```
gmx_MMPBSA --create_input gb
```

For PB,
```
gmx_MMPBSA --create_input pb
```

For GB, PB, decomposition
```
gmx_MMPBSA --create_input gb pb decomp
```

modify `mmpbsa.in` 

## Run Calculation
```
gmx_MMPBSA -O -i mmpbsa.in -cs md_0_100.tpr -ct md_0_100_center.xtc -ci index.ndx -cg 1 13 -cp topol.top -o FINAL_MMPBSA.dat -eo FINAL_MMPBSA.csv
```