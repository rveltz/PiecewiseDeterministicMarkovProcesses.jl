abstract type AbstractJump end

# Dummy Jump function
function Delta_dummy(xc, xd, t, parms, ind_reaction)
	return nothing
end

struct RateJump{Td, Tnu <: AbstractArray{Td}, TD} <: AbstractJump
	nu::Tnu						# implements jumps on the discrete variable with a matrix
	Delta::TD		    		# function to implement the jumps (optional)

	function RateJump(nu::Tnu, DX::TD) where {Td, Tnu <: AbstractArray{Td}, TD}
		return new{Td, Tnu, TD}(nu, DX)
	end

	function RateJump(DX::TD) where {TD}
		return new{Int64, Array{Int64,2}, TD}(zeros(Int64, 0, 0), DX)
	end

	function RateJump(nu::Tnu) where {Td, Tnu <: AbstractArray{Td}}
		return new{Td, Tnu, typeof(Delta_dummy)}(nu, Delta_dummy)
	end

end

# perform the jump on the discrete variable
function affect!(ratejump::RateJump, ev, xd)
	deltaxd = view(ratejump.nu, ev, :)
	@inbounds for ii in eachindex(xd)
		xd[ii] += deltaxd[ii]
	end
end

# perform the jump on the continuous variable
function affect!(ratejump::RateJump, ev, u, xd, t, parms)
	ratejump.Delta(u, xd, t, parms, ev)
end
