# Installation

## Installation of the package

Julia version 1.0 or later is required.

```julia
]add https://github.com/matsunagalab/MDToolbox.jl.git
using MDToolbox
```

## Docker image

Using the docker image, MDToolbox.jl can be readily used within REPL or Jupyter(Lab).

```bash
## REPL
$ docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work matsunagalab/mdtoolbox.jl julia

## JupyterLab
$ docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work matsunagalab/mdtoolbox.jl
```
