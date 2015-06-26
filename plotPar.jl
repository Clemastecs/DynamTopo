function plotPar(par::Array{Float64})
	"""
		Plot the set of particles
		
			INPUT:
	    		X:      Set of particles
	    		
	    	OUTPUT:
	    		
	    	ANNOTATIONS:
	    		Needs a plotting pkg to draw
	"""

	#Pyplot
	p = par[par[:,3].==1,:];
	PyPlot.plot( p[:,1], p[:,2],"b.", alpha= 0.6, antialiased=true)

	p = par[par[:,3].==2,:];
	PyPlot.plot( p[:,1], p[:,2],"r.", alpha= 0.6, antialiased=true)
	
	p = par[par[:,3].==3,:];
	PyPlot.plot( p[:,1], p[:,2],"c.", alpha= 0.6, antialiased=true)

end