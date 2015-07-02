function shapefunc(geom::Int64,pospg::Array{Float64})
	#=
	    This function computes the shape functions for Q2Q1-elements.

	    	INPUT:
	    		geom:      The nodal geometry of the element (nen, nenP)
	    		pospg:	   Position of the Gauss points on reference element

	    	OUTPUT:
	    		Nxi:  		Array of the derivatives shape functions of the velocity
	    		Neta:		Array of the derivatives shape functions of the velocity
	    		N:			Array of the shape functions of the velocity/pressure
	=#

  # Variable Declaration
  xi::Array{Float64} = pospg[:,1]
  eta::Array{Float64} = pospg[:,2]

  N::Array{Float64} = []
  Nxi::Array{Float64} = []
  Neta::Array{Float64} = []

  # Compute the funcform
  if geom == 4  #element Q1
     N    = [ (1-xi).*(1-eta)/4 (1+xi).*(1-eta)/4 (1+xi).*(1+eta)/4 (1-xi).*(1+eta)/4 ]
     Nxi  = [ (eta-1)/4 (1-eta)/4 (1+eta)/4 -(1+eta)/4 ]
     Neta = [ (xi-1)/4 -(1+xi)/4 (1+xi)/4 (1-xi)/4 ]
  elseif geom == 9  #element Q2
     N    = [ xi.*(xi-1).*eta.*(eta-1)/4 xi.*(xi+1).*eta.*(eta-1)/4 xi.*(xi+1).*eta.*(eta+1)/4 xi.*(xi-1).*eta.*(eta+1)/4 (1-xi.^2).*eta.*(eta-1)/2 xi.*(xi+1).*(1-eta.^2)/2 (1-xi.^2).*eta.*(eta+1)/2 xi.*(xi-1).*(1-eta.^2)/2 (1-xi.^2).*(1-eta.^2) ]
     Nxi  = [ (xi-1/2).*eta.*(eta-1)/2 (xi+1/2).*eta.*(eta-1)/2 (xi+1/2).*eta.*(eta+1)/2 (xi-1/2).*eta.*(eta+1)/2 -xi.*eta.*(eta-1) (xi+1/2).*(1-eta.^2) -xi.*eta.*(eta+1) (xi-1/2).*(1-eta.^2) -2*xi.*(1-eta.^2) ]
     Neta = [ xi.*(xi-1).*(eta-1/2)/2 xi.*(xi+1).*(eta-1/2)/2 xi.*(xi+1).*(eta+1/2)/2 xi.*(xi-1).*(eta+1/2)/2 (1-xi.^2).*(eta-1/2) xi.*(xi+1).*(-eta) (1-xi.^2).*(eta+1/2) xi.*(xi-1).*(-eta) (1-xi.^2).*(-2*eta) ]
  end

  return N, Nxi, Neta
  
end
