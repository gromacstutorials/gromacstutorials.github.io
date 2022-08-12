# Install Instruction

1. Untar or copy the downloaded force field directory where you want to store it.
The best place is by far be in the standard GROMACS force field directory (`{$GMX_DIR}/share/top/`), next to the original versions of the GROMOS force field.
Alternatively, if you are on a shared machine (e.g cluster) where you don't have right access to the original topology directory,
you may want to store them in your local user partition. In that case, you should always provided absolute paths in your `#include` statements (like in the example provided below).

2. Modify your GROMACS topology ('.top') file to include the ATB modified topology
instead of the standard provided with GROMACS.

An example is provided below.

# Example Top File


```
; Include forcefield parameters
#include "/marksw/gromacs-5.0.4/share/gromacs/top/top_atb/ffG54a7.itp"

; Include water topology
#include "/marksw/gromacs-5.0.4/share/gromacs/top/top_atb/spc.itp"

; Include generic topology for ions
#include "/marksw/gromacs-5.0.4/share/gromacs/top/top_atb/ions.itp"

[ system ]
; Name
Protein in water

[ molecules ]
; Compound        #mols
Protein         1
SOL         22844
```


