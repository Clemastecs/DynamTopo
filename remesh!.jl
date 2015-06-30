function remesh!(XP::Array{Float64}, X::Array{Float64}, nx::Int64, ngap::Int64, gap::Float64)
	#=
		This function Re-meshes the mesh on the top subdomain depending on ngap value.

			INPUT:
				XP: 		Nodal Coordenates of pressure
				X: 			Nodal Coordenates of velocity
	    		nx,ny: 		Mesh dimensions (nx)x(ny)-elements
	    		ngap:		Layer mesh dimensions (nx)x(ngap)-elements

	    	OUTPUT:
 				This function modify XP and X.

 			ANNOTATIONS:
 				The remesh is based on the intervals projection rule py=(px-x1)/(x2-x1)*(y2-y1)+y1
	=#

	# Variable declaration
  	i::Int64 = 0
  	x1::Float64 = 0
  	x2::Float64 = 0
  	y1::Float64 = 0
  	y2::Float64 = 0

  	# Remeshing the pressure nodes
  	x1 = XP[end-(nx+1)*(2*ngap+1)+1,2]
  	x2 = XP[end,2]
  	y1 = XP[end-(nx+1)*(2*ngap+1)+1,2]
  	y2 = XP[end-(nx+1)*(ngap+1)+1,2]
  	y2=divrem(itrunc(y2),10)[1]*10  # Mantain the domain at (x1,x2)x(y1,y2+2/gap)

  	for i = (size(XP,1)-(nx+1)*(2*ngap+1)+1):size(XP,1)
  		XP[i,2] = (XP[i,2]-x1)/(x2-x1) * (y2-y1) + y1
  	end

  	# Remeshing the velocity nodes
  	x1 = X[end-(2*nx+1)*(4*ngap+1)+1,2]
  	x2 = X[end,2]
  	y1 = X[end-(2*nx+1)*(4*ngap+1)+1,2]
  	y2 = X[end-(2*nx+1)*(2*ngap+1)+1,2]
  	y2=divrem(itrunc(y2),10)[1]*10  # Mantain the domain at (x1,x2)x(y1,y2+2/gap)

  	for i = (size(X,1)-(2*nx+1)*(4*ngap+1)+1):size(X,1)
  		X[i,2] = (X[i,2]-x1)/(x2-x1) * (y2-y1) + y1
  	end

end
