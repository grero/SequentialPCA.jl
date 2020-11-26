using LinearAlgebra
using StatsBase

struct SeqPCA
    proj::Matrix{Float64}
    idx::Vector{Int64}
end


function StatsBase.fit(::Type{SeqPCA}, Y::Matrix{T};ncomps=size(Y,1), kvs...) where T <: Real
    U,idx = seqPCA(Y,ncomps)
    SeqPCA(U,idx)
end

"""
Find the axis that accounts for the larget amount of variance up to a specifc time in the signal `Y`.
"""
function seqPCA(Y::Matrix{TT}) where TT <: Real
	D,T = size(Y)
	vmax = -Inf
	idx = 0
	um = fill(0.0, D)
	for i in D:T
		u,s,vt = svd(Y[:,1:i],full=true)
		v = s[1]*s[1]/sum(s.*s)
		if v > vmax
			vmax = v
			idx = D+i-1
			um .= u[:,1]
		end
	end
	um, idx
end

function seqPCA(Y, ncomp::Int64)
	D,T = size(Y)
    U = fill(0.0, D, ncomp)
    idx = fill(0, ncomp)
    U[:,1], idx[1] = seqPCA(Y)
    Yp = copy(Y)
    for i in 2:ncomp
        Yp .= (I - U[:,i-1]*U[:,i-1]')*Yp
        U[:,i], idx[i] = seqPCA(Yp)
    end
    U, idx
end

"""
Orthogonalizes the columsn on `U` usnig the Gramm-Schmitt method
"""
function orthogonalize(U::Matrix{T}) where T <: Real
    d,n = size(U)
    V = similar(U)
    V[:,1] = U[:,1]
    for i in 2:n
        V[:,i] = U[:,i] - (U[:,i]'*V[:,1:(i-1)]).*V[:,1:(i-1)]
    end
    V
end
