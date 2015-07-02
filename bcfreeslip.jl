function bcfreeslip(X::Array{Float64},nunk::Int64)
	#=
	    This function imposes the boundary conditions using lagrange multipliers.

	    	INPUT:
	    		X:      Nodal Coordenates of velocity
	    		nunk:   Number of unknows of the velocity

	    	OUTPUT:
	    		Accd:	  BC matrix
	    		bccd:     BC vector

	=#

  # Variables declaration
  xmin::Float64 = minimum(X[:,1])
  xmax::Float64 = maximum(X[:,1])
  ymin::Float64 = minimum(X[:,2])
  ymax::Float64 = maximum(X[:,2])

  nodesDx0::Array{Float64} = find(X[:,1].==xmin)
  nodesDx1::Array{Float64} = find(X[:,1].==xmax)
  nodesDy0::Array{Float64} = find(X[:,2].==ymin)
  nodesDy1::Array{Float64} = find(X[:,2].==ymax)

  Cy0::Array{Float64} = [ reshape([2*nodesDy0'],size(nodesDy0,1),1) reshape([zeros(1,size(nodesDy0,1))],size(nodesDy0,1),1) ]
  Cx1::Array{Float64} = [ reshape([2*nodesDx1'-1],size(nodesDx1,1),1) reshape([zeros(1,size(nodesDx1,1))],size(nodesDx1,1),1) ]
  Cx0::Array{Float64} = [ reshape([2*nodesDx0'-1],size(nodesDx0,1),1) reshape([zeros(1,size(nodesDx0,1))],size(nodesDx0,1),1) ]
  Cy1::Array{Float64} = [ reshape([2*nodesDy1'],size(nodesDy1,1),1) reshape([zeros(1,size(nodesDy1,1))],size(nodesDy1,1),1) ]

  C::Array{Float64} = [Cy0; Cx1; Cx0; Cy1]

  # BC Matrix
  nDir::Int64 = size(C,1)
  Accd::Array{Float64} = zeros(nDir,nunk)
  Accd[:,C[:,1]] = eye(nDir)

  # BC Vector
  bccd::Array{Float64}=[]
  bccd = C[:,2]

  return Accd, bccd
end
