function updateparticles!( par::Array{Float64}, velo::Array{Float64}, nx::Int64, ny::Int64, x1::Int64, x2::Int64, y1::Int64, y2::Int64, rad::Float64)
	"""
	    This function updates the set of particles using the velocity field.

	    	INPUT:
	    		par:     	  Array of particles with own properties
	    		velo: 		  Array of velocities
	    		nx,ny: 		  Mesh dimensions (nx)x(ny)-elements
	    		x1,x2,y1,y2:  Domain [x1,x2]x[y1,y2]
	    		rad:		  Dimension of the smallest node's side

	    	OUTPUT:
	    		This function modify the par

	    	ANNOTATIONS:
				Needs the Grid.jl pkg to uses CoordInterpGrid.
				Remember the construction of velocity mesh to reshape and interpolate
	"""

  alpha::Float64 = 0.5
  xrange::FloatRange{Float64} = x1:x2/(2*nx):x2
  yrange::FloatRange{Float64} = y1:y2/(2*ny):y2
  dt::Float64 = rad / sqrt( maximum( sum( velo'.^2 ,1) ) ) * alpha

  vx::Array{Float64} = reshape(velo[:,1], 2*nx+1, 2*ny+1 )
  vy::Array{Float64} = reshape(velo[:,2], 2*nx+1, 2*ny+1 )

  #Interpolation
  fvx_i = CoordInterpGrid((xrange,yrange), vx, BCnearest, InterpLinear)
  fvy_i = CoordInterpGrid((xrange,yrange), vy, BCnearest, InterpLinear)

  vx_i::Array{Float64} = {fvx_i[par[I,1:2]...] for I = 1:size(par,1)}
  vy_i::Array{Float64} = {fvy_i[par[I,1:2]...] for I = 1:size(par,1)}

  #Euler explicit (y_n=y_n+dt*f(y_n))
  par[:,1:2] = par[:,1:2] + [vx_i vy_i] * dt

end
