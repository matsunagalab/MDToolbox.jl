# Makorv state models

```@docs
msmgenerate(nframe::Int, T, pi_i)
```

```@docs
msmgenerate(nframe::Int, T, pi_i, emission)
```

```@docs
msmcountmatrix(indexOfCluster; tau=1)
```

```@docs
msmtransitionmatrix(C; TOLERANCE=10^(-4), verbose=true)
```

```@docs
msmviterbi(observation, T, pi_i, emission)
```

```@docs
msmbaumwelch(data_list, T0, pi_i0, emission0; TOLERANCE=10.0^(-4), MAXITERATION=Inf64)
```

```@docs
msmplot(T; pi_i=nothing, x=nothing, y=nothing, filename=nothing, edgewidth_scale=3.0, arrow_scale=0.0001, nodesize=0.5, fontsize=10, names=[], dpi=100)
```
