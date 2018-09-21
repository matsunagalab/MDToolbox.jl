run with jupyter
```bash
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work ymatsunaga/mdtoolbox.jl
```

run with console
```bash
docker run -it --rm -v "$PWD":/home/jovyan/work ymatsunaga/mdtoolbox.jl julia
```

