@ilmproblem NeumannPoisson scalar

struct NeumannPoissonCache{SMT,DVT,VNT,FT,ST} <: AbstractExtraILMCache
 S :: SMT
 dvn :: DVT
 vn :: VNT
 fstar :: FT
 sstar :: ST
end

function ImmersedLayers.prob_cache(prob::NeumannPoissonProblem,base_cache::BasicILMCache)
  S = create_CLinvCT(base_cache)
  dvn = zeros_surface(base_cache)
  vn = zeros_surface(base_cache)
  fstar = zeros_grid(base_cache)
  sstar = zeros_gridcurl(base_cache)
  NeumannPoissonCache(S,dvn,vn,fstar,sstar)
end

function ImmersedLayers.solve(prob::NeumannPoissonProblem,sys::ILMSystem)
  @unpack extra_cache, base_cache, bc, phys_params = sys
  @unpack S, dvn, vn, fstar, sstar = extra_cache

  fill!(fstar,0.0)
  fill!(sstar,0.0)

  f = zeros_grid(base_cache)
  s = zeros_gridcurl(base_cache)
  df = zeros_surface(base_cache)
  ds = zeros_surface(base_cache)

  # Get the precribed jump and average of the surface normal derivatives
  prescribed_surface_jump!(dvn,sys)
  prescribed_surface_average!(vn,sys)

  # Find the potential
  regularize!(fstar,dvn,base_cache)
  inverse_laplacian!(fstar,base_cache)

  surface_grad!(df,fstar,base_cache)
  df .= vn - df
  df .= -(S\df);

  surface_divergence!(f,df,base_cache)
  inverse_laplacian!(f,base_cache)
  f .+= fstar

  # Find the streamfunction
  surface_curl!(sstar,df,base_cache)

  surface_grad_cross!(ds,fstar,base_cache)
  ds .= S\ds

  surface_curl_cross!(s,ds,base_cache)
  s .-= sstar
  s .*= -1.0

  inverse_laplacian!(s,base_cache)

  return f, df, s, ds
end

# Boundary condition functions for x, y, and rotation

function dfdn_nx(base_cache,phys_params)
    nrm = normals(base_cache)
    dfdn = zeros_surface(base_cache)
    copyto!(dfdn,nrm.u,base_cache,1)
    return dfdn
end

function dfdn_ny(base_cache,phys_params)
    nrm = normals(base_cache)
    dfdn = zeros_surface(base_cache)
    copyto!(dfdn,nrm.v,base_cache,1)
    return dfdn
end

function dfdn_xcrossn(base_cache,phys_params)
    bl = base_cache.bl
    nrm = normals(base_cache)
    pts = points(base_cache)

    # Center of body
    Xc, Yc = bl[1].cent
    pts.u .-= Xc
    pts.v .-= Yc

    dfdn = zeros_surface(base_cache)
    copyto!(dfdn,pointwise_cross(pts,nrm),base_cache,1)
    return dfdn
end
