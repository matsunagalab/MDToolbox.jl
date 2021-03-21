# MDToolbox.jl
MDToolbox: A Julia toolbox for analysis and modeling of molecular dynamics trajectories

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.bio.ics.saitama-u.ac.jp/MDToolbox.jl/stable/) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://www.bio.ics.saitama-u.ac.jp/MDToolbox.jl/dev/)

[![Build Status](https://travis-ci.org/ymatsunaga/MDToolbox.jl.svg?branch=master)](https://travis-ci.org/ymatsunaga/MDToolbox.jl)

Some function are ported into Julia from the original MATLAB version https://github.com/ymatsunaga/mdtoolbox

Julia version 1.0 or later is required. 
```
julia> ]add https://github.com/matsunagalab/MDToolbox.jl.git
julia> using MDToolbox
```

Docker image is available. MDToolbox can be used in Jupyter Lab.
```bash
# Terminal
$ docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work ymatsunaga/mdtoolbox.jl
# in Jupyter Lab
using MDToolbox

```
