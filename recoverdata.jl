function recoverdata(nx::Int64, case::Int64,nSteps::Int64=1)
"""
	This function recovers the data of the computations.

			INPUT:
	    		case:  1-fixed, 2-airlayer
	    		nx:    Mesh dimensions (nx)x(ny)-elements

	    	OUTPUT:

"""

	if case == 1

		sigmazz = recoversigma("sigmazz_elem"*string(nx)*".dat")
		plotfix(sigmazz,nx,nSteps,[1420.],40)

	elseif case == 2

		surface = recoversurface("surface_elem"*string(nx)*".dat")
		plotlayer(0,40,surface,20,nx,nSteps)

	else

		println("The case"*string(case)*"do not exists. Try again.")

	end

end

function recoversigma(f::String)
"""
	Recover the height-plots. Fixed case.

			INPUT:
	    		f:				String with the path of file

	    	OUTPUT:
	    		sigmazz:		Array with sigmazz values

"""

	file = open("./results/"*f)# openfile
	simgazz = readdlm( file, '\t' )
	close(file)

return simgazz'

end

function recoversurface(f::String)
"""
	Recover the height-plots. Fixed case.

			INPUT:
	    		f:				String with the path of file

	    	OUTPUT:
	    		surface:		Array with surface values

"""
	file = open("./results/"*f)# openfile
	surface = readdlm( file, '\t' )
	close(file)

return surface'

end
