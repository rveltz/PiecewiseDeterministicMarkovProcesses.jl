abstract type AbstractJump end

# Dummy Jump function
function Delta_dummy(xc, xd, parms, t, ind_reaction)
	return nothing
end

struct Jump{Td, Tnu <: AbstractArray{Td}, TD} <: AbstractJump
	nu::Tnu						# implements jumps on the discrete variable with a matrix
	Delta::TD					# function to implement the jumps (optional)

	function Jump(nu::Tnu, DX::TD) where {Td, Tnu <: AbstractArray{Td}, TD}
		return new{Td, Tnu, TD}(nu, DX)
	end

	function Jump(DX::TD) where {TD}
		return new{Int64, Array{Int64,2}, TD}(zeros(Int64, 0, 0), DX)
	end

	function Jump(nu::Tnu) where {Td, Tnu <: AbstractArray{Td}}
		return new{Td, Tnu, typeof(Delta_dummy)}(nu, Delta_dummy)
	end

end

get_rate_prototype(jp::Jump, Tc) = zeros(Tc, size(jp.nu, 1))

function affect!(ratejump::Jump, ev, xc, xd, parms, t)
	# perform the jump on the discrete variable
	deltaxd = view(ratejump.nu, ev, :)
	@inbounds for ii in eachindex(xd)
		xd[ii] += deltaxd[ii]
	end
	# perform the jump on the continuous variable
	ratejump.Delta(xc, xd, parms, t, ev)
end
