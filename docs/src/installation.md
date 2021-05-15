# Installation

## Installation of the package

Julia version 1.6 or later is required.

```julia
]add https://github.com/matsunagalab/MDToolbox.jl.git
using MDToolbox
```

## Docker image

A docker image for MDToolbox.jl is available,

```bash
## REPL
$ docker run -it --rm -v "$PWD":/home/jovyan/work matsunagalab/mdtoolbox julia

## JupyterLab
$ docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan matsunagalab/mdtoolbox
```
