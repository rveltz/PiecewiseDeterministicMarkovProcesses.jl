abstract type AbstractRate end

# These are different structs introduced to deal with the case where the rate function is constant in between jumps: this are defined through the following type `ConstantRate`. This leads us to treat the case where the user can provide two rate functions, the first one being a `ConstantRate` and a second one `VariableRate` where no a-priori knowledge in infused by the user. These two functions are encapsulated in a `CompositeRate` structure. A composite rate `r::CompositeRate` is called like `r(rate, xc, xd, p, t, issum)`. In this case, the two components of `r` act on the same rate vector `rate` so the indexing of `rate` inside `r.Rcst` and `r.Rvar` should be considered global by the user and not local.

# this is used to initialise the rate component of the structure `PDMPCaracteristics`. This is useful so that a call to `solve` with identical seed always return the same trajectory
init!(R) = nothing

struct VariableRate{TR} <: AbstractRate
	R::TR
end

function (vr::VariableRate)(rate, xc, xd, p, t, issum)
	return vr.R(rate, xc, xd, p, t, issum)
end

# Structure meant to deal with rate functions which are constant in between jumps. The only way the rate function can change is when `xd` changes. Hence, while c::ConstantRate is called like `c(rate, xc, xd, p, t, true)`, it returns `c.totalrate`. In the codes CHV and Rejection, a call to `c(rate, xc, xd, p, t, false)` that a jump has occurred and one wants to (possibly) change `xd`. We use this call to trigger the update of `c.totalrate`. This update is also triggered whenever `c.totalrate < 0` like for initialisation purposes.
mutable struct ConstantRate{TR} <: AbstractRate
	R::TR
	totalrate::Float64
	function ConstantRate(R)
		return new{typeof(R)}(R, -1.0)
	end
end

init!(r::ConstantRate) = r.totalrate = -1.0

function (cr::ConstantRate)(rate, xc, xd, p, t, issum)
	if issum == true
		if cr.totalrate < 0
			# update the catched value
			cr.totalrate = cr.R(rate, xc, xd, p, t, issum)[1]
		end
		return cr.totalrate, cr.totalrate
	else
		# the following call will be amortized if we call the method twice
		cr.totalrate = -1
		cr.R(rate, xc, xd, p, t, issum)
	end
end

struct CompositeRate{TRc, TRv} <: AbstractRate
	Rcst::TRc
	Rvar::TRv
	function CompositeRate(Rc, Rv)
		rc = ConstantRate(Rc)
		rv = VariableRate(Rv)
		return new{typeof(rc), typeof(rv)}(rc, rv)
	end
end

function (vr::CompositeRate)(rate, xc, xd, p, t, issum)
	if issum == false
		vr.Rcst(rate, xc, xd, p, t, issum)
		vr.Rvar(rate, xc, xd, p, t, issum)
	else
		out_cst = vr.Rcst(rate, xc, xd, p, t, issum)
		out_var = vr.Rvar(rate, xc, xd, p, t, issum)
		return out_cst .+ out_var
	end
end
