FROM jupyter/datascience-notebook

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
    julia -e 'import Pkg; Pkg.add(["PyPlot", "Plots", "Flux", "Optim", "NLsolve", "StatsBase", "FFTW", "Revise", "BenchmarkTools", "MAT", "BSON", "JLD2", "FileIO"])' && \
    julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(url="https://github.com/matsunagalab/MDToolbox.jl.git", rev="master"))'

RUN julia -e 'using IJulia, Plots, Flux, Optim, NLsolve, StatsBase, FFTW, BenchmarkTools, MAT, BSON, JLD2, FileIO' && \
    julia -e 'using MDToolbox'

RUN fix-permissions $JULIA_PKGDIR $CONDA_DIR/share/jupyter
