function source(xc::Array{Float64}, par::Array{Float64}, sradius::Float64)
	#=
	    This function computes the source term of the matrix system (K,G)(v,p)=(f).

	    	INPUT:
	    		xc:			Coordinates of the particles in cartesian base
	    		sradius:    Radius to choose the mean of density and viscosity
	    		par:		Set of particles

	    	OUTPUT:
	    		f: 			Array of the source term
	=#

  g::Float64 = -9.8
  nPar::Int64 = size(par, 1)
  ix::Array{Float64} = []

  ix = find(sum( [par[:,1:2] - repmat(xc,nPar,1)]'.^2,1 ) .< sradius^2 )

  meanrho::Float64 = mean(par[ix,4])
  f::Array{Float64} = [0; meanrho*g]

  return f
  
end
