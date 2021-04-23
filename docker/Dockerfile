FROM jupyter/datascience-notebook
MAINTAINER Yasuhiro Matsunaga <ymatsunaga@mail.saitama-u.ac.jp>

USER $NB_UID
ENV DEBIAN_FRONTEND noninteractive
ENV OMP_NUM_THREADS 8
ENV JULIA_NUM_THREADS 8

USER $NB_UID
RUN conda install --quiet --yes \
    'qt' \
    'pyqt' && \
    conda clean -tipsy && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf /home/$NB_USER/.cache/yarn && \
    rm -rf /home/$NB_USER/.node-gyp && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

RUN julia -e 'import Pkg; Pkg.update()' && \
    julia -e 'import Pkg; Pkg.add("PyPlot")' && \
    julia -e 'import Pkg; Pkg.add("Plots")' && \
    julia -e 'import Pkg; Pkg.add("Flux")' && \
    julia -e 'import Pkg; Pkg.add("DifferentialEquations.jl")' && \
    julia -e 'import Pkg; Pkg.add("Optim")' && \
    julia -e 'import Pkg; Pkg.add("NLsolve")' && \
    julia -e 'import Pkg; Pkg.add("StatsBase")' && \
    julia -e 'import Pkg; Pkg.add("FFTW")' && \
    julia -e 'import Pkg; Pkg.add("Revise")' && \
    julia -e 'import Pkg; Pkg.add("ProfileView")' && \
    julia -e 'import Pkg; Pkg.add("BenchmarkTools")' && \
    julia -e 'import Pkg; Pkg.add("HDF5")' && \
    julia -e 'import Pkg; Pkg.add("JLD2")' && \
    julia -e 'import Pkg; Pkg.add("MAT")' && \
    #julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(name="GR", version="0.35.0"))' && \
    #julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(url="https://github.com/halleysfifthinc/MAT.jl", rev="v0.7-update"))' && \
    julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(url="https://github.com/matsunagalab/MDToolbox.jl.git", rev="master"))' && \
    julia -e 'import Pkg; Pkg.add("FileIO")'

RUN julia -e 'using IJulia' && \
    #julia -e 'pygui(false)' && \
    julia -e 'using PyPlot' && \
    julia -e 'using Plots' && \
    julia -e 'using Flux' && \
    julia -e 'using DifferentialEquations.jl' && \
    julia -e 'using Optim' && \
    julia -e 'using NLsolve' && \
    julia -e 'using StatsBase' && \
    julia -e 'using FFTW' && \
    julia -e 'using Revise' && \
    julia -e 'using ProfileView' && \
    julia -e 'using BenchmarkTools' && \
    julia -e 'using HDF5' && \
    julia -e 'using JLD2' && \
    #julia -e 'using MAT' && \
    julia -e 'using MDToolbox' && \
    julia -e 'using FileIO'

RUN fix-permissions $JULIA_PKGDIR $CONDA_DIR/share/jupyter
