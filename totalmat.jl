function totalmat(X::Array{Float64},T::Array{Float64},XP::Array{Float64},TP::Array{Float64},par::Array{Float64},sradius::Float64)
	#=
	    This function computes the total components of the system matrix (K,G)(v,p)=(f).

			INPUT:
	    		X:        Nodal Coordenates of velocity
	    		T:        Nodal Connectivities of velocity
	    		XP:       Nodal Coordenates of pressure
	    		TP:       Nodal Connectivities of pressure
	    		par:      Set of particles
	    		sradius:  Radius to choose the mean of density and viscosity

	    	OUTPUT:
	    		K:	  		 Diffussion matrix
	    		G:   		 Gradient matrix
	    		f:   		 Font vector

	=#

  # Variables Declaration
  numel::Int64 = size(T,1)
  nen::Int64 = size(T,2)
  nenP::Int64 = size(TP,2)
  numnp::Int64 = size(X,1)
  nunk::Int64  = 2*numnp
  ndofn::Int64 = 2*nen
  nunkP::Int64 = size(XP,1)
  i::Int64 = 0

  (pospg,weipg) = quadrature()
  (N,Nxi,Neta) = shapefunc(nen,pospg)
  NP = shapefunc(nenP,pospg)[1]

  Te::Vector{Int64} = []
  TeP::Vector{Int64} = []
  Xe::Array{Float64} = []
  Ke::Array{Float64} = []
  Ge::Array{Float64} = []
  fe::Array{Float64} = []

  K::SparseMatrixCSC{Float64,Int64} = spzeros(nunk,nunk)
  G::SparseMatrixCSC{Float64,Int64} = spzeros(nunkP,nunk)
  f::Array{Float64} = zeros(nunk,1)

  # Elements Loop
  for i = 1:numel
      # Reorder
      Te = vec( reshape([2*T[i,:]-1; 2*T[i,:]],1,ndofn) )
      TeP = vec( TP[i,:] )

      # Element Nodes Coords
      T2 = vec(T[i,1:nen])
      Xe = X[T2,:]

      # Local Matrix
      (Ke,Ge,fe) = localmat(Xe,nen,ndofn,pospg,weipg,N,Nxi,Neta,nenP,NP,par,sradius)

      # Local Matrix Assembly
      K[Te,Te]  = K[Te,Te] + Ke
      G[TeP,Te] = G[TeP,Te]+ Ge
      f[Te] = f[Te] + fe
  end

  return K, G, f
  
end
