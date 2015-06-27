function  mpresnod(X::Array{Float64},T::Array{Float64},nx::Int64,ny::Int64,nenP::Int64)
  """
		This function creates the topology of the mesh of pressure in [x1,x2]x[y1,y2] domain with (nx)*(ny) Q2Q1-elements

			INPUT:
	    		X: 		Nodal Coordenates of velocity
				T: 		Nodal Connectivities of velocity
	    		nx,ny: 	Mesh dimensions (nx)x(ny)-elements
	    		nenP:		Number of the pressure nodes in each element

	    	OUTPUT:
 				XP: 		Nodal Coordenates of pressure
				TP: 		Nodal Connectivities of pressure
	"""

  # Variables Declaration
  numel::Int64 = size(T,1)
  npx::Int64 = nx+1
  npy::Int64 = ny+1
  i::Int64 = 0
  j::Int64 = 0
  ielem::Int64 = 0
  inode::Int64 = 0

  XP::Array{Float64} = zeros(npx*npy,2)
  TP::Array{Float64} = zeros(numel, nenP)

  # Compute of XP
  for i=1:npy
    iodd = 2*(i-1)+1
    XP[(i-1)*npx+1:i*npx,:] = X[(iodd-1)*(2*nx+1)+1 : 2 : iodd*(2*nx+1),:]
  end

  # Compute of TP
  for i=1:ny, j=1:nx
      ielem = (i-1)*nx+j;
      inode = (i-1)*(npx)+j;
      TP[ielem,:] = [inode inode+1 inode+(npx+1) inode+(npx)]
  end

  return XP, TP

end
