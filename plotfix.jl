function plotfix(sigmazz::Array{Float64}, nx::Int64, nSteps::Int64, rho::Array{Float64}, x2::Int64)
	"""
		This function computes the height and plots the results of the fixed case.

			INPUT:
	    		sigmazz:  The stress on the surface for each step
	    		nx:       Mesh dimensions (nx)x(ny)-elements
	    		nSteps:	 Number of the time steps
	    		rho:      Array of the densities
	    		x2:       Domain [x1,x2]x[y1,y2]

	    	OUTPUT:


	    	ANNOTATIONS:
	    		Needs a plotting pkg to draw
	"""

	  # geological height and aux variables
	  h::Array{Float64} = zeros(2*nx) # Geological height (2*nx: number of gauss points in compStress)
	  maxh::Array{Float64} = zeros(nSteps)
	  i::Int64 = 0

	  h =  sigmazz/(-9.8*rho[1]) # Expression to compute the heigth

	  for i=2:nSteps
	  		maxh[i] = h[nx,i] # heigth value in the middle in each step
	  end

	  sigmadom::Array{Float64}= linspace((x2/nx)/4,x2-(x2/nx)/4,size(sigmazz,1))
	  alpha::Array{Float64} = linspace(0,1,nSteps)

	  # draw stresses in each step
	  PyPlot.figure(2)
	  PyPlot.clf()
	  PyPlot.subplot(2,2,1)
	  PyPlot.ylabel("Stress \$\\sigma_{zz}\$")
	  PyPlot.xlabel("Surface")
	  PyPlot.title("Surface Stress")
	  for i=2:nSteps-1
	  		PyPlot.plot(sigmadom,sigmazz[:,i], linewidth=1.2, alpha= alpha[i], antialiased=true, color="cyan")
	  end
	  PyPlot.plot(sigmadom,sigmazz[:,i+1], linewidth=2, alpha= 1, antialiased=true, color="black")
	  PyPlot.grid()

	  #PyPlot.subplot(2,2,2,)

	  # draw height in each step
	  PyPlot.subplot(2,2,3)
	  PyPlot.ylabel("Height h")
	  PyPlot.xlabel("Surface")
	  PyPlot.title("Geological Height")
	  for i=2:nSteps-1
	  		PyPlot.plot(sigmadom,h[:,i], linewidth=1.2, alpha= alpha[i], antialiased=true, color="red")
	  end
	  PyPlot.plot(sigmadom,h[:,i+1], linewidth=2, alpha= 1, antialiased=true,color="black")
	  PyPlot.grid()

	  PyPlot.subplot(2,2,4)
	  PyPlot.ylabel("Max h")
	  PyPlot.xlabel("Steps")
	  PyPlot.title("Maximum Height")
	  PyPlot.plot([2:nSteps], maxh[2:end], linewidth=2,color="red")
	  PyPlot.grid()
	  print_with_color(:red,"h max: "*string(maxh[end])*"\n")

	  #PyPlot.savefig("Geo_High.jpg")

end
