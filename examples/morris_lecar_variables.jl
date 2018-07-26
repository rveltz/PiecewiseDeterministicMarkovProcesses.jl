struct mlParams
  v_Na::Float64
  g_Na::Float64
  v_K::Float64
  g_K::Float64
  v_L::Float64
  g_L::Float64
  I_app::Float64
  gamma_na::Float64
  k_na::Float64
  beta_na::Float64
  gamma_k::Float64
  k_k::Float64
  beta_k::Float64
  M::Float64
  N::Float64
end

function ml(p)
  return mlParams(p["v_Na"] , p["g_Na"] , p["v_K"],p["g_K"] , p["v_L"] , p["g_L"] , p["I_app"] , p["gamma_na"] , p["k_na"] , p["beta_na"] , p["gamma_k"] , p["k_k"] , p["beta_k"] , p["M"] , p["N"])
end
