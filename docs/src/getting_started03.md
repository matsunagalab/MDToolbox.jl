# Atom selection

Using `TrjArray`-type for storing MD trajectory data, users can flexibly select specific atoms and their trajectory data. 
Users can select atoms as like slicing arrays from the buit-in `Array` type in Julia. 

```julia
using MDToolbox
ENV["COLUMNS"] = 130; # if you are using Jupyter, the column width should be set for the message width from MDToolbox.jl

t = mdload("examples/data/3gb1.pdb")
    32x855 TrjArray{Float64, Int64}
    | A                          | A                          |  …   A                          |
    | 1MET                       | 1MET                       |  …   56GLU                      |
    | 1N                         | 2CA                        |  …   855HG3                     |
    |   -13.15    -1.71     5.51 |   -12.20    -2.85     5.70 |  …      10.24     2.24    -4.53 |
    |   -13.16    -2.75     4.19 |   -12.38    -3.43     5.26 |         10.29     1.79    -3.20 |
    |   -13.26    -2.20     5.12 |   -12.23    -3.26     5.29 |         11.09     3.24    -6.13 |
    |   -12.45    -4.30     6.28 |   -12.25    -3.22     5.28 |         11.75     2.57    -3.31 |
    |   -13.05    -1.52     5.31 |   -12.15    -2.68     5.56 |         10.46     1.51    -3.50 |
    |   -12.98    -3.74     6.65 |   -12.43    -3.31     5.34 |  …      11.26     3.06    -3.59 |
    |   -13.03    -2.10     3.64 |   -12.47    -2.79     4.83 |         10.38     1.53    -3.30 |
    |   -13.00    -2.17     5.21 |   -12.09    -3.19     5.78 |         11.32     3.37    -5.78 |
    |   -12.57    -2.96     3.87 |   -12.12    -3.42     5.21 |         10.36     1.38    -3.53 |
    |   -12.60    -2.85     4.27 |   -12.16    -3.02     5.68 |         10.44     1.40    -3.37 |
    |   -13.12    -1.61     4.72 |   -12.34    -2.60     5.52 |  …      10.96     3.58    -5.68 |
    |             ⋮              |             ⋮              |  ⋱               ⋮              |
    |   -13.46    -2.54     6.05 |   -12.41    -3.14     5.18 |         10.33     1.63    -3.72 |
    |   -12.61    -4.31     6.34 |   -12.22    -3.51     5.15 |         10.26     1.56    -3.59 |
    |   -13.24    -2.72     6.28 |   -12.16    -3.34     5.45 |         11.30     2.91    -3.41 |
    |   -13.21    -2.29     5.61 |   -12.13    -3.31     5.53 |         11.76     2.84    -5.65 |
    |   -12.75    -2.42     4.33 |   -12.22    -2.96     5.61 |  …      11.53     3.05    -5.67 |
    |   -13.32    -2.12     5.10 |   -12.27    -3.17     5.25 |         10.38     1.88    -3.32 |
    |   -13.07    -2.31     4.32 |   -12.27    -3.22     5.19 |         11.65     3.11    -5.90 |
    |   -12.76    -3.69     6.68 |   -12.11    -3.33     5.38 |         11.51     2.97    -5.63 |
    |   -13.21    -2.55     5.53 |   -12.06    -3.41     5.11 |         10.26     1.60    -3.20 |
    |   -13.13    -2.12     4.11 |   -12.34    -3.00     5.02 |  …      12.68     3.28    -4.88 |
    |   -12.81    -2.86     4.32 |   -12.23    -3.33     5.62 |          9.98     2.33    -5.34 |
```

First, specific frames of `TrjArray` variable can be selected by row indices. For example, the 3rd to 5th frames can be sliced by,

```julia
t[3:5, :]
    3x855 TrjArray{Float64,Int64}
    | A                          | A                          |  …   A                          |
    | 1MET                       | 1MET                       |  …   56GLU                      |
    | 1N                         | 2CA                        |  …   855HG3                     |
    |   -13.26    -2.20     5.12 |   -12.23    -3.26     5.29 |  …      11.09     3.24    -6.13 |
    |   -12.45    -4.30     6.28 |   -12.25    -3.22     5.28 |         11.75     2.57    -3.31 |
    |   -13.05    -1.52     5.31 |   -12.15    -2.68     5.56 |         10.46     1.51    -3.50 |
```

Next, specific atoms of `TrjArray` variable can be selected by column indices. For example, the 1st to 3rd atoms can be sliced by,

```julia
t[:, 1:3]
    32x3 TrjArray{Float64,Int64}
    | A                          | A                          | A                          |
    | 1MET                       | 1MET                       | 1MET                       |
    | 1N                         | 2CA                        | 3C                         |
    |   -13.15    -1.71     5.51 |   -12.20    -2.85     5.70 |   -10.79    -2.30     5.87 |
    |   -13.16    -2.75     4.19 |   -12.38    -3.43     5.26 |   -11.06    -2.70     5.50 |
    |   -13.26    -2.20     5.12 |   -12.23    -3.26     5.29 |   -10.88    -2.62     5.57 |
    |             ⋮              |             ⋮              |             ⋮              |
    |   -13.21    -2.55     5.53 |   -12.06    -3.41     5.11 |   -10.74    -2.72     5.46 |
    |   -13.13    -2.12     4.11 |   -12.34    -3.00     5.02 |   -11.00    -2.35     5.33 |
    |   -12.81    -2.86     4.32 |   -12.23    -3.33     5.62 |   -10.89    -2.63     5.85 |
```

Of couse, both frame and atom selections can be specified together, 

```julia
t[3:5, 1:3]
    3x3 TrjArray{Float64,Int64}
    | A                          | A                          | A                          |
    | 1MET                       | 1MET                       | 1MET                       |
    | 1N                         | 2CA                        | 3C                         |
    |   -13.26    -2.20     5.12 |   -12.23    -3.26     5.29 |   -10.88    -2.62     5.57 |
    |   -12.45    -4.30     6.28 |   -12.25    -3.22     5.28 |   -10.90    -2.54     5.53 |
    |   -13.05    -1.52     5.31 |   -12.15    -2.68     5.56 |   -10.73    -2.16     5.80 |
```

Also, atoms can be selected by strings. For example, residue ID == 10 can be selected by

```julia
t[:, "resid 10"]
    32x22 TrjArray{Float64,Int64}
    | A                          | A                          |  …   A                          |
    | 10LYS                      | 10LYS                      |  …   10LYS                      |
    | 1N                         | 2CA                        |  …   22HZ3                      |
    |    14.51     1.47    -1.06 |    15.96     1.11    -0.95 |  …      15.72     3.47    -6.83 |
    |    14.34     0.95    -1.15 |    15.79     0.64    -1.11 |         18.17     4.93    -4.94 |
    |    14.53     1.21    -1.25 |    15.98     0.88    -1.17 |         13.99     4.58    -3.58 |
    |             ⋮              |             ⋮              |  ⋱               ⋮              |
    |    14.35     1.24    -1.20 |    15.81     0.95    -1.20 |         14.40     4.35    -5.36 |
    |    14.48     1.14    -1.19 |    15.92     0.79    -1.15 |  …      18.17     5.89    -5.06 |
    |    14.56     1.13    -1.20 |    16.01     0.83    -1.11 |         14.94     2.44    -7.45 |
```

The colon in the row can be omitted,

```julia
t["resid 10"]
    32x22 TrjArray{Float64,Int64}
    | A                          | A                          |  …   A                          |
    | 10LYS                      | 10LYS                      |  …   10LYS                      |
    | 1N                         | 2CA                        |  …   22HZ3                      |
    |    14.51     1.47    -1.06 |    15.96     1.11    -0.95 |  …      15.72     3.47    -6.83 |
    |    14.34     0.95    -1.15 |    15.79     0.64    -1.11 |         18.17     4.93    -4.94 |
    |    14.53     1.21    -1.25 |    15.98     0.88    -1.17 |         13.99     4.58    -3.58 |
    |             ⋮              |             ⋮              |  ⋱               ⋮              |
    |    14.35     1.24    -1.20 |    15.81     0.95    -1.20 |         14.40     4.35    -5.36 |
    |    14.48     1.14    -1.19 |    15.92     0.79    -1.15 |  …      18.17     5.89    -5.06 |
    |    14.56     1.13    -1.20 |    16.01     0.83    -1.11 |         14.94     2.44    -7.45 |
```

Ranges can be specified using colons. For example, residue IDs rom 10 to 13 can be selected by

```julia
t["resid 10:13"]
32x77 TrjArray{Float64,Int64}
| A                          | A                          |  …   A                          |
| 10LYS                      | 10LYS                      |  …   13LYS                      |
| 1N                         | 2CA                        |  …   77HZ3                      |
|    14.51     1.47    -1.06 |    15.96     1.11    -0.95 |  …      15.42    -3.79     2.49 |
|    14.34     0.95    -1.15 |    15.79     0.64    -1.11 |         16.92    -1.15     2.36 |
|    14.53     1.21    -1.25 |    15.98     0.88    -1.17 |         14.80    -0.42     7.32 |
|             ⋮              |             ⋮              |  ⋱               ⋮              |
|    14.35     1.24    -1.20 |    15.81     0.95    -1.20 |         16.18    -2.37     2.19 |
|    14.48     1.14    -1.19 |    15.92     0.79    -1.15 |  …      12.10    -4.45     4.27 |
|    14.56     1.13    -1.20 |    16.01     0.83    -1.11 |         16.27    -3.58     4.58 |
```

Not only IDs, names can be specified. For example, residue name of GLU ASP can be selected by

```julia
t["resname GLU ASP"]
32x136 TrjArray{Float64,Int64}
| A                          | A                          |  …   A                          |
| 15GLU                      | 15GLU                      |  …   56GLU                      |
| 1N                         | 2CA                        |  …   136HG3                     |
|     5.87    -0.14     7.35 |     4.52    -0.66     6.99 |  …      10.24     2.24    -4.53 |
|     5.97    -0.41     7.47 |     4.59    -0.80     7.07 |         10.29     1.79    -3.20 |
|     5.82    -0.05     7.35 |     4.47    -0.58     6.99 |         11.09     3.24    -6.13 |
|             ⋮              |             ⋮              |  ⋱               ⋮              |
|     5.95    -0.37     7.62 |     4.60    -0.89     7.26 |         10.26     1.60    -3.20 |
|     5.96    -0.12     7.34 |     4.61    -0.66     7.02 |  …      12.68     3.28    -4.88 |
|     5.91    -0.11     7.29 |     4.54    -0.56     6.92 |          9.98     2.33    -5.34 |
```

These can be combined by using `and`, `or`, and brackets `()`. For example, C_alpha atoms of residue IDs from 1 to 20 can be selected by,

```julia
t["atomname CA and resid 1:20"]
    32x20 TrjArray{Float64,Int64}
    | A                          | A                          |  …   A                          |
    | 1MET                       | 2THR                       |  …   20ALA                      |
    | 1CA                        | 2CA                        |  …   20CA                       |
    |   -12.20    -2.85     5.70 |    -8.56    -2.53     6.83 |  …     -11.12     2.18     3.06 |
    |   -12.38    -3.43     5.26 |    -8.87    -2.64     6.54 |        -11.10     2.14     3.34 |
    |   -12.23    -3.26     5.29 |    -8.72    -2.69     6.68 |        -11.29     2.23     3.10 |
    |             ⋮              |             ⋮              |  ⋱               ⋮              |
    |   -12.06    -3.41     5.11 |    -8.63    -2.73     6.65 |        -11.07     2.11     3.15 |
    |   -12.34    -3.00     5.02 |    -8.86    -2.42     6.49 |  …     -11.11     2.25     3.12 |
    |   -12.23    -3.33     5.62 |    -8.71    -2.59     6.92 |        -11.02     2.10     3.14 |
```

The following keywords are avaiable for atom selections by strings

|  keywords       | description                        |  examples                                  |
| :-------------- | :--------------------------------- | :----------------------------------------- |
| chainname       | specify chain names                | `chainname A` `chainname A B`              |
| chainid         | specify chain IDs                  | `chainid 1`  `chainid 1:3` `chainid 1:3 5` |
| resname         | specify residue names              | `resname ARG` `resname GLU ASP`            |
| resid           | specify residue IDs                | `resid 1` `resid 1:3` `resid 1:3 5`        |
| atomname        | specify atom names                 | `atomname CA` `atomname C N O CA`          |
| atomid          | specify atom IDs                   | `atomid 1` `atomid 1:3` `atomid 1:3 5`     |
| protein         | select proteins                    | `protein`                                  |
| solvent         | select solvents                    | `solvent`                                  |
| water           | select water molecues              | `water`                                    |
| hydrogen        | select hydrogens                   | `hydrogen`                                 |
