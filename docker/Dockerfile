FROM jupyter/datascience-notebook
MAINTAINER Yasuhiro Matsunaga <ymatsunaga@riken.jp>

USER $NB_UID
ENV DEBIAN_FRONTEND noninteractive
ENV OMP_NUM_THREADS 8

RUN julia --version

RUN julia -e 'import Pkg; Pkg.update()' && \
    julia -e 'import Pkg; Pkg.add("PyPlot")' && \
    julia -e 'import Pkg; Pkg.add("Revise")' && \
    julia -e 'import Pkg; Pkg.add("BenchmarkTools")' && \
    julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(url="https://github.com/ymatsunaga/MDToolbox.jl.git", rev="master"))' && \
    julia -e 'using IJulia' && \
    fix-permissions $JULIA_PKGDIR $CONDA_DIR/share/jupyter
