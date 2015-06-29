function localmat(Xe::Array{Float64},nen::Int64,ndofn::Int64,pospg::Array{Float64},weipg::Array{Float64},N::Array{Float64},Nxi::Array{Float64},Neta::Array{Float64},nenP::Int64,NP::Array{Float64},par::Array{Float64},sradius::Float64)
	"""
	    This function computes the local components of the system matrix fo FEM: (K,G)(v,p)=(f) 

	    	INPUT:
	    		Xe:         Nodal Coordenates of velocity
	    		nen:        Number of the volocity nodes in each element
	    		ndofn:      Number of dofs in velocity element
	    		pospg:	   Position of the Gauss points on reference element
	    		weipg:      Weigth of the Gauss points on reference element
	    		Nxi:  		Array of the derivatives shape functions of the velocity
	    		Neta:			Array of the derivatives shape functions of the velocity
	    		nenP:			Number of the pressure nodes in each element
	    		NP:			Array of the shape functions of the pressure
	    		N:				Array of the shape functions of the velocity
	    		sradius:    Radius to choose the mean of density and viscosity
	    		par:			Set of particles

	    	OUTPUT:
	    		Ke:	  		 Element diffussion matrix
	    		Ge:   		 Element gradient matrix
	    		fe:   		 Element font vector
	"""

  #Declaration of Variables
  Ke::Array{Float64} = zeros(ndofn,ndofn)
  Ge::Array{Float64} = zeros(nenP,ndofn)
  fe::Array{Float64} = zeros(ndofn,1)

  igaus::Int64 = 0
  ngaus::Int64 = size(pospg,1)
  N_igaus::Array{Float64} = []
  Nxi_igaus::Array{Float64} = []
  Neta_igaus::Array{Float64} = []
  NP_igaus::Array{Float64} = []
  f_igaus::Array{Float64} = []
  Jacob::Array{Float64} = []
  dvolu::Float64 = 0
  res::Array{Float64} = []
  Nx::Array{Float64} = []
  Ny::Array{Float64} = []
  dN::Array{Float64} = []
  Ngp::Array{Float64} = []
  meanvisc::Float64 = 1
  ix::Array{Float64} = []
  nPar::Int64 = size(par, 1)
  j::Int64 = 0

  # Loop on Quadrature Gauss points
  for igaus = 1:ngaus

     # Function shape values & Local Gauss points derivatives
     N_igaus = N[igaus,:]
     Nxi_igaus = Nxi[igaus,:]
     Neta_igaus = Neta[igaus,:]
     NP_igaus = NP[igaus,:]

     # Jacobian Matrix on Gauss point
     Jacob = [ Nxi_igaus[1:nen]'*Xe[1:nen,1]	Nxi_igaus[1:nen]'*Xe[1:nen,2]; Neta_igaus[1:nen]'*Xe[1:nen,1]	Neta_igaus[1:nen]'*Xe[1:nen,2] ]

     # Weight per Jacobian determinant
     dvolu=weipg[igaus].*det(Jacob)

     # Functions form derivatives on global coords
     res = Jacob\[Nxi_igaus;Neta_igaus]

     #Functions form gradient
     Nx = [reshape([1;0]*res[1,:],1,ndofn); reshape([0;1]*res[1,:],1,ndofn)]
     Ny = [reshape([1;0]*res[2,:],1,ndofn); reshape([0;1]*res[2,:],1,ndofn)]

     # Local Matrix
	 ix = find(sum( [par[:,1:2] - repmat(N_igaus*Xe,nPar,1)]'.^2,1 ) .< sradius^2 ) # choose a viscosity
  	 meanvisc = mean(par[ix,5])

     Ke = Ke + meanvisc*(Nx'*Nx+Ny'*Ny)*dvolu

     # Function form divergence
     dN = reshape(res,1,ndofn)
     Ge = Ge - NP_igaus'*dN*dvolu

     f_igaus = source(N_igaus*Xe, par, sradius)

     Ngp = [reshape([1;0].*N_igaus,1,ndofn); reshape([0;1].*N_igaus,1,ndofn)]

     fe = fe + Ngp'*f_igaus*dvolu

  end

  return Ke, Ge, fe;

end
