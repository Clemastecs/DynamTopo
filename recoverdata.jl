function recoverdata(nx::Int64, case::Int64,nSteps::Int64 = 1)
	#=
		This function recovers the data of the computations and plot.

			INPUT:
	    		case:  1-fixed, 2-airlayer
	    		nx:    Mesh dimensions (nx)x(ny)-elements

	    	OUTPUT:

	=#

	if case == 1
		sigmazz = recoverfile("sigmazz_elem"*string(nx)*".dat")
		plotfix(sigmazz,nx,nSteps,[1420.],40)
	elseif case == 2
		surface = recoverfile("surface_elem"*string(nx)*".dat")
		plotlayer(0,40,surface,20,nx,nSteps)
	else
		println("The case"*string(case)*"do not exists. Try again.")
	end

end

function recoverfile(f::String)
	#=
		Recover the height-plots. Fixed case.

			INPUT:
	    		f:				String with the path of file

	    	OUTPUT:
	    		values:		Array with sigmazz/surface values

    =#

	file = open("./results/"*f)# openfile
	values = readdlm( file, '\t' )
	close(file)

	return values'
end
