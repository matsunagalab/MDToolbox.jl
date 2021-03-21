<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.bio.ics.saitama-u.ac.jp/MDToolbox.jl/stable/) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://www.bio.ics.saitama-u.ac.jp/MDToolbox.jl/dev/)
<!-- [![Build Status](https://travis-ci.org/ymatsunaga/MDToolbox.jl.svg?branch=master)](https://travis-ci.org/ymatsunaga/MDToolbox.jl) -->

# MDToolbox.jl

A Julia package for molecular dynamics trajectories analysis and modeling of biomolecules

Some function are ported into Julia from an old MATLAB package https://github.com/ymatsunaga/mdtoolbox

Julia version 1.0 or later is required. 
```
julia> ]add https://github.com/matsunagalab/MDToolbox.jl.git
julia> using MDToolbox
```

Docker image is available. MDToolbox can be used in Jupyter Lab.
```bash
# Terminal
$ docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work matsunagalab/mdtoolbox.jl
# in Jupyter Lab
using MDToolbox

```
