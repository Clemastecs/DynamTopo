function setparticles(x1::Int64, x2::Int64, y2::Int64, nx::Int64, ny::Int64, ppe::Int64, radius::Float64, visc::Array{Float64}, rho::Array{Float64}, pto::Array{Float64}, air::Bool, gap::Float64)
	#=
		This function creates the set of particles.

			INPUT:
	    		x1,x2,y2: 	Domain [x1,x2]x[y1,y2]
				ppe: 			Number of particler per element
	    		nx,ny: 		Mesh dimensions (nx)x(ny)-elements
	    		radius:		The radius of the sphere
	    		rho:			Array of densities
	    		visc:			Array of viscosities
	    		pto:			Center of the sphere
	    		gap			Layer dimensions (nx)x(gap)
	    		air:			Bool index for air layer

	    	OUTPUT:
 				nPar: 		Dimension of the array of particles
				par: 			Array of particles with own properties
	=#

	# Variable declaration.
	ix::Array{Int64} = []
	par::Array{Float64} = []
	nPar::Int64 = 0

	xx::Array{Float64} = linspace(x1, x2, itrunc(sqrt(ppe)*nx)) # position grid of the particles

	# Build the set of particles
	par = [xx repmat([xx[1]], length(xx), 1) ones(length(xx)) rho[1]*ones(length(xx)) visc[1]*ones(length(xx))]
	for i = 2:length(xx)
		par = [ par; xx repmat([xx[i]], length(xx), 1) ones(length(xx)) rho[1]*ones(length(xx)) visc[1]*ones(length(xx))]
	end

	nPar = size(par, 1) #number of particles

	# Properties of the sphere particles
	ix = find(sum((par[:,1:2] - repmat(pto,nPar,1))'.^2,1) .< radius^2) # sphere
	par[ix,3] = 2
	# density/viscosity of the sphere particles
	par[ix,4] = rho[2]
	par[ix,5] = visc[2]

	# Air layer
	if air == true
		# Properties of the particles
		#ix = find(par[:,2] .> (y2-gap/2)) # remeshed air layer
		ix = find(par[:,2] .> 35 )# no remeshed air layer case 40x40 only
		par[ix,3] = 3
		# density/viscosity of the air layer particles
		par[ix,4] = rho[3]
		par[ix,5] = visc[3]
	end

	return nPar, par
end
