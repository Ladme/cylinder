# cylinder: Calculating Density Inside A Cylinder 

Calculates density of atoms inside a cylinder as a function of the cylinder height.

## Dependencies

`cylinder` requires you to have groan library installed. You can get groan from [here](https://github.com/Ladme/groan). See also the [installation instructions](https://github.com/Ladme/groan#installing) for groan.

## Installation

1) Run `make groan=PATH_TO_GROAN` to create a binary file `cylinder` that you can place wherever you want. `PATH_TO_GROAN` is a path to the directory containing groan library (containing `groan.h` and `libgroan.a`).
2) (Optional) Run `make install` to copy the the binary file `cylinder` into `${HOME}/.local/bin`.

## Options

```
Usage: cylinder -c GRO_FILE -s SELECTION -p CYLINDER_POSITION [OPTION]...

OPTIONS
-h               print this message and exit
-c STRING        gro file to read
-f STRING        xtc file to read (optional)
-n STRING        ndx file to read (optional, default: index.ndx)
-o STRING        output file (default: cylinder.xvg)
-s STRING        selection of atoms to analyze
-p STRING        selection of atoms defining cylinder position
-x/y/z           direction of the main axis of the cylinder (default: z)
-r FLOAT         radius of the cylinder in nm (default: 2.5)
-e FLOAT         height of the cylinder in nm (default: height of the box in gro file)
-d INTEGER       grid spacing along the specified axis in points per nm (default: 10)
```

Use the [groan selection language](https://github.com/Ladme/groan#groan-selection-language) specify selections of atoms (flags `-s` and `-p`).

## Examples

```
cylinder -c md.gro -f md.xtc -s "name OH" -p Protein -z -r 3.1 -d 20
```

The program will calculate the average density of atoms named `OH` (flag `-s`) from file `md.xtc` (flag `-f`) along the main axis of a cylinder that is oriented along the z-axis (flag `-z`) and has its center positioned in the geometric center of selection `Protein` (flag `-p`). The cylinder has a radius of 3.1 nm (flag `-r`) and a height corresponding to the z-dimension of the simulation box in `md.gro`. Information about the ndx groups (such as `Protein`) will be read from `index.ndx` (default option). Density will be calculated for z-axis slices separated by 0.05 (1/20) nm (flag `-d`). 

The output will be written into `cylinder.xvg` (default option) and can be visualized using `xmgrace` (`xmgrace cylinder.xvg`). `cylinder` also calculates the average number of atoms inside the entire cylinder and writes this information to the end of the output file.

```
cylinder -c md.gro -f md.xtc -s "name OH" -p "Protein; Protein; resname POPC" -z -r 3.1 -d 20
```

Same as above except the cylinder is positioned in the geometric center of `Protein` in the xy-plane, while on the z-axis, it will be placed in the geometric center of the selection `resname POPC`.

```
cylinder -c md.gro -f md.xtc -s "name OH" -p "Protein; Protein; resname POPC" -x -r 3.1 -d 20 -e 5
```

Same as above except the cylinder is oriented along the x-axis (flag `-x`) and has a height of 5 nm (flag `-e`). That means that the density will only be calculated for x-dimension range of -2.5 to 2.5 nm relative to the cylinder center.

## Limitations

Assumes that the simulation box is rectangular and that periodic boundary conditions are applied in all three dimensions.

Always uses center of geometry instead of center of mass.

Only tested on Linux. Probably will not work on anything that is not UNIX-like.

