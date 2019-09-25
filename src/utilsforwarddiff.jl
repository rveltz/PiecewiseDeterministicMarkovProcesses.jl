# Code to get PiecewiseDeterministicMarkovProcesses working with ForwardDiff
# Makes use of the trick in http://docs.juliadiffeq.org/latest/basics/faq.html#Are-the-native-Julia-solvers-compatible-with-autodifferentiation?-1
using ForwardDiff

struct DiffCache{T<:AbstractArray, S<:AbstractArray}
	rate::T
	dual_rate::S
end

function DiffCache(u::AbstractArray{T}, siz, ::Type{Val{chunk_size}}) where {T, chunk_size}
	DiffCache(u, zeros(ForwardDiff.Dual{nothing,T,chunk_size}, siz...))
end

dualcache(u::AbstractArray, N=Val{ForwardDiff.pickchunksize(length(u))}) = DiffCache(u, size(u), N)

get_rate(dc::DiffCache, u::AbstractArray{T}) where T<:ForwardDiff.Dual = reinterpret(T, dc.dual_rate)
get_rate(dc::DiffCache, u::AbstractArray) = dc.rate
