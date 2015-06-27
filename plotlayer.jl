function plotlayer(x1::Int64, x2::Int64, surface::Array{Float64}, ppe::Int64, nx::Int64, nSteps::Int64)
	"""
		This function computes the height surface using plotfit() and plots the results of the air layer case.

			INPUT:
	    		surface:  The height on the surface for each step
	    		nx:       Mesh dimensions (nx)x(ny)-elements
	    		nSteps:	 	 Number of the time steps
	    		rho:      	 Array of the densities
	    		x1,x2:       Domain [x1,x2]x[y1,y2]

	    	OUTPUT:


	    	ANNOTATIONS:
	    		Needs a plotting pkg to draw
	"""
	i::Int64 = 0
	sdom::Array{Float64}= linspace(x1,x2,size(surface,1))
	alpha::Array{Float64} = linspace(0,1,nSteps)
	x=linspace(x1+1,x2-1)
	coef::Array{Float64} = []
	maxh::Array{Float64} = zeros(size(surface,2))

	PyPlot.figure(3)
	PyPlot.clf()
	PyPlot.subplot(1,2,1)
	PyPlot.ylabel("Height h")
	PyPlot.xlabel("Surface")
	PyPlot.title("Geological Height with layer")
	for i = 2:nSteps-1
		coef = polyfit(sdom,surface[:,i]-surface[:,2],4)
		p(x) = coef[5].*x.^4+coef[4].*x.^3+coef[3].*x.^2+coef[2].*x+coef[1]
		maxh[i] = p(20)

		#PyPlot.plot(x,p(x),linewidth=1.2, alpha= alpha[i], antialiased=true, color="red")
		PyPlot.plot(sdom,surface[:,i]-surface[:,2], linewidth=1.2, alpha= alpha[i], antialiased=true, color="red")
	end
	coef = polyfit(sdom,surface[:,i+1]-surface[:,2],4)
	p(x) = coef[5].*x.^4+coef[4].*x.^3+coef[3].*x.^2+coef[2].*x+coef[1]
	maxh[i+1] = p(20)

	#PyPlot.plot(x,p(x),linewidth=2, alpha= 1, antialiased=true,color="black")
	PyPlot.plot(sdom,surface[:,i+1]-surface[:,2], linewidth=2, alpha= 1, antialiased=true,color="black")
	PyPlot.grid()


	PyPlot.subplot(1,2,2)
	PyPlot.ylabel("Max h")
	PyPlot.xlabel("Steps")
	PyPlot.title("Maximum Height")
	PyPlot.plot([2:size(surface,2)], maxh[2:end], linewidth=2, color="red")
	PyPlot.grid()

end

function polyfit(x::Array{Float64}, y::Array{Float64}, n::Int64)
	"""
		Compute the coeficients of n-polynomial associated to (x,y) points

			INPUT:
	    		x, y: 	Points to interpolate
	    		n: 		Dregree of the interpolation

	    	OUTPUT:
	    		Array with polynomial's coeficients

	    	ANNOTATIONS:
	    		Based on least-square method
	"""
   A::Array{Float64} = []
	i::Int64 = 0

   A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
   return A \ y

end
