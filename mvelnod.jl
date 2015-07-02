function mvelnod(x1::Int64,x2::Int64,y1::Int64,y2::Int64,nx::Int64,ny::Int64,nen::Int64)
	#=
		This function creates the topology of the mesh of velocity in [x1,x2]x[y1,y2] domain with (nx)*(ny) Q2Q1-elements

			INPUT:
	    		nen:		  Number of the velocity nodes in each element
	    		nx,ny: 		  Mesh dimensions (nx)x(ny)-elements
	    		x1,x2,y1,y2:  Domain [x1,x2]x[y1,y2]

	    	OUTPUT:
 				X: 			  Nodal Coordenates of velocity
				T: 			  Nodal Connectivities of velocity
	=#

  # Variables Declaration
  npx::Int64 = 2*nx+1 # Auxiliar dimensions to create the mesh with 9-nodes
  npy::Int64 = 2*ny+1
  i::Int64 = 0
  j::Int64 = 0
  ielem::Int64 = 0
  inode::Int64 = 0
  k::Int16 = 0

  X::Array{Float64} = zeros((npx)*(npy),2)
  xs::Array{Float64} = linspace(x1,x2,npx)'
  ys::Array{Float64} = zeros(npy,1)
  uns::Array{Float64} = ones(npx,1)
  T::Array{Float64} = zeros(nx*ny,nen)
  yys::Array{Float64} =[]

  yys = linspace(y1,y2,npy)

  # Compute of X
  for i=1:npy
    ys = yys[i]*uns
    k = 1
    for j=((i-1)*npx+1):(i*npx)
        X[j,:] = [xs[k], ys[k]]
        k = k+1
    end
  end

  # Compute of T
  for i = 1:ny, j = 1:nx
      ielem = (i-1)*nx+j
      inode = (i-1)*2*npx+2*(j-1)+1
      T[ielem,:] =[ inode inode+2 inode+2*(npx)+2 inode+2*(npx) inode+1 inode+2+(npx) inode+2*(npx)+1 inode+(npx) inode+(npx)+1 ]
  end

  return X, T

end
