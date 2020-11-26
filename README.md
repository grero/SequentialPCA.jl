# Sequential PCA

## Introduction

This implements the sequential PCA algorithm mentioned in [[1]](#aoietal).

## Usage 

The algorithm expects an input matrix of dimensions *Dxt* where *D* represents the dimension of the data and *t* represents time. The output is a matrix of projections that capture successive temporal components of the signal. Have a look at this example from the tests.

```julia
using SequentialPCA
using StatsBase
using JLD2
fname = joinpath(dirname(pathof(SequentialPCA)), "..","test","seqpca_testdata.jd2")

@load fname Y t
spca = fit(SeqPCA, permutedims(Y,[2,1]);ncomps=2)
```


## References
<a name="aoietal">[1]</a> Aoi, M. C., Mante, V., & Pillow, J. W. (2020). Prefrontal cortex exhibits multidimensional dynamic encoding during decision-making. Nature Neuroscience, 23(11), 1410–1420. http://doi.org/10.1038/s41593-020-0696-5

