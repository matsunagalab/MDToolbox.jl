# Workflow for developers

This note explains general workflows for developers who are interested in adding their own functions in MDToolbox.jl

## Setup of MDToolbox.jl source codes

1. Fork the MDToolbox.jl repository to your space

Just press the fork button on GitHub

2. Download (clone) the forked repository

In your shell terminal, please clone the forked repository with a git command:

```
$ git clone https://github.com/your_account_name/MDToolbox.jl.git
```

3. Register the downloaded source codes as developmental codes in Julia

```
$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.1.0 (2019-01-21)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> 
# enter the package mode by pressing ]
pkg> develop /path/to/MDToolbox.jl
# return to the REPL mode by pressing BACKSPACE
# now you can use MDToolbox
julia> using MDToolbox
```

## Workflow for adding your original function

1. Edit `src/MDToolbox.jl`

Suppose that your new function name is `wham()`, and the definition of `wham()` is now being coded in `wham.jl`.

Then, please write your original function name and filename in `src/MDToolbox.jl` file as follows:

```
###

export wham

###
include("wham.jl")
```

2. Write your codes in `wham.jl`

During the coding, it is recommended to execute test codes in your REPL or jupyter. Also, it is convenient to use `Reivse` package which enables you to modify codes, making changes active without restarting Julia.

```
julia> using Revise, MDToolbox
julia> # run your test codes
```
