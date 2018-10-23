FROM jupyter/datascience-notebook
MAINTAINER Yasuhiro Matsunaga <ymatsunaga@riken.jp>

USER $NB_UID
ENV DEBIAN_FRONTEND noninteractive
ENV OMP_NUM_THREADS 8
ENV JULIA_NUM_THREADS 8

RUN julia -e 'import Pkg; Pkg.update()' && \
    julia -e 'import Pkg; Pkg.add("PyPlot")' && \
    julia -e 'import Pkg; Pkg.add("Flux")' && \
    julia -e 'import Pkg; Pkg.add("Optim")' && \
    julia -e 'import Pkg; Pkg.add("NLsolve")' && \
    julia -e 'import Pkg; Pkg.add("FFTW")' && \
    julia -e 'import Pkg; Pkg.add("Revise")' && \
    julia -e 'import Pkg; Pkg.add("ProfileView")' && \
    julia -e 'import Pkg; Pkg.add("BenchmarkTools")' && \
    julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(url="https://github.com/ymatsunaga/MDToolbox.jl.git", rev="master"))' && \
    julia -e 'using IJulia' && \
    julia -e 'using PyPlot' && \
    julia -e 'using Flux' && \
    julia -e 'using Optim' && \
    julia -e 'using NLsolve' && \
    julia -e 'using FFTW' && \
    julia -e 'using Revise' && \
    julia -e 'using ProfileView' && \
    julia -e 'using BenchmarkTools' && \
    julia -e 'using HDF5' && \
    julia -e 'using JLD' && \
    julia -e 'using MDToolbox' && \
    fix-permissions $JULIA_PKGDIR $CONDA_DIR/share/jupyter
