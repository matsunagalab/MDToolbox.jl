# Workflow for developers

This note explains typical workflows for developers who are interested in adding their own functions in MDToolbox.jl

## Setup of MDToolbox.jl source codes

### 1. Fork the MDToolbox.jl repository https://github.com/matsunagalab/MDToolbox.jl to your in GitHub. Just press the fork button to fork the repo. 

### 2. Clone (download) the forked repository. In your terminal, run the following command: 

```bash
$ git clone https://github.com/your_account_name/MDToolbox.jl.git
```

### 3. Add the downloaded source codes as developmental codes in Julia

```julia
$ julia
julia> 
# enter the package mode by pressing ]
pkg> develop /path/to/MDToolbox.jl
# return to the REPL mode by pressing BACKSPACE or DELETE
julia> using MDToolbox
```

## Workflow for adding your original function

### 1. Edit `src/MDToolbox.jl`

Suppose that your new function name is `wham()`, and the definition of `wham()` is now being coded in `wham.jl`.

Then, let's write your original function name and filename in `src/MDToolbox.jl` file as follows:

```julia
###

export wham

###
include("wham.jl")
```

### 2. Write your codes in `wham.jl`

During the coding, it is recommended to execute test codes in your REPL or jupyter. Also, it is convenient to use `Reivse` package which enables you to modify codes, making changes active without restarting Julia.

```julia
julia> using Revise, MDToolbox
julia> # run your test codes
```
