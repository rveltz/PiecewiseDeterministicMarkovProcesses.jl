# Code to get PiecewiseDeterministicMarkovProcesses working with ForwardDiff
# Makes use of the trick in http://docs.juliadiffeq.org/latest/basics/faq.html#Are-the-native-Julia-solvers-compatible-with-autodifferentiation?-1
using ForwardDiff

struct DiffCache{T<:AbstractArray, S<:AbstractArray}
	rate::T
	dual_rate::S
end

function DiffCache(T, size, ::Type{Val{chunk_size}}) where chunk_size
	DiffCache(zeros(T, size...), zeros(ForwardDiff.Dual{nothing,T,chunk_size}, size...))
end

DiffCache(u::AbstractArray) = DiffCache(eltype(u),size(u),Val{ForwardDiff.pickchunksize(length(u))})

get_rate(dc::DiffCache, ::Type{T}) where {T<:ForwardDiff.Dual} = dc.dual_rate
get_rate(dc::DiffCache, T) = dc.rate
