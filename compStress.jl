function compStress(velo::Array{Float64},pres::Array{Float64},X::Array{Float64},T::Array{Float64},TP::Array{Float64},par::Array{Float64},nx::Int64,sradius::Float64)
	"""
	    Compute the sigmazz on the top boundary

	    	INPUT:
	    		velo:    Array of velocities in x and y
	    		pres:    Array of pressures
	    		X:			Nodal Coordenates of velocity
	    		T:       Connectivity Matrix of velocity
	    		TP:      Connectivity Matrix of pressure

	    	OUTPUT:
	    		sigmazz: The solution of sigma on the top boundary

	    	ANNOTATIONS:
				sigma(z,z) = -(pressure)*NP_boundary+2*(visc)*N_boundary*(diff(velo))
	"""

	# Variable declaration
	sigmazz::Array{Float64} = []
	spres::Array{Float64} = []
	svelo::Array{Float64} = []
	numel::Int64 = size(T,1)
	nen::Int64 = size(T,2)
	nenP::Int64 = size(TP,2)
	ngeom::Int64 = 9
	ssol::Array{Float64}= [0,0]

	pospg::Array{Float64} = []
	N::Array{Float64} = []
	Nxi::Array{Float64} = []
	Neta::Array{Float64} = []
	NP::Array{Float64} = []
	TeP::Array{Float64} = []
	Pe::Array{Float64} = []
	res::Array{Float64} = []
	Vze::Array{Float64} = []
	Xe::Array{Float64} = []
	Te::Array{Float64} = []
	Jacob::Array{Float64} = []
	nPar::Int64 = size(par, 1)
	j::Int64 = 0
	i::Int64 = 0
	ix::Array{Float64} = {}
	
	N_igaus::Array{Float64} = []
	Nxi_igaus::Array{Float64} = []
	Neta_igaus::Array{Float64} = []
	
	
	pospg=[-1/2 1; 1/2 1] # quadrature for gauss points on the top boundary
	ngaus::Int64 = size(pospg,1)

	(N,Nxi,Neta) = shapeFunc(nen,pospg)
	NP = shapeFunc(nenP,pospg)[1]

	for i = (numel-nx+1):numel # Only the top elements

		# Velocity solution on the top boundary
		Te = vec(T[i,:])		
    	Xe = X[Te,:]
		Vze = velo[Te,2] # Vertical velocities

		# Pressure solution on the top boundary
		TeP = vec(TP[i,:])
		Pe = pres[TeP]
		spres = NP*Pe


		for igaus = 1:ngaus

			# Function form values & Local Gauss points derivatives
			N_igaus = N[igaus,:]
			Nxi_igaus = Nxi[igaus,:]
      	Neta_igaus = Neta[igaus,:]

     		# Jacobian Matrix on Gauss point
     		Jacob = [ Nxi_igaus[1:ngeom]'*Xe[1:ngeom,1]	Nxi_igaus[1:ngeom]'*Xe[1:ngeom,2]; Neta_igaus[1:ngeom]'*Xe[1:ngeom,1]	Neta_igaus[1:ngeom]'*Xe[1:ngeom,2] ]

     		# Functions form derivatives on global coords
     		res = Jacob\[Nxi_igaus; Neta_igaus]
     		
     		# Choose a viscosity
     		#j = 1
     		#while ix == []
		  		ix = find(sum( [par[:,1:2] - repmat(N_igaus*Xe,nPar,1)]'.^2,1 ) .< sradius^2 )
		  	#	sradius = 2*sradius
		  	#	j = j + 1
		  	#	if j == 100
		  	#		continue
		  	#	end
     		#end
     		
  	  		meanvisc = mean(par[ix,5])
  	  		      	
  			svelo = 2*meanvisc*res[2,:]*Vze

  			# Solution of sigmazz
			ssol[igaus] =  svelo[1] - spres[igaus]

		end

		sigmazz = [sigmazz; ssol]

	end

	return sigmazz

end
