function quadrature()
	"""
	    This function constructs the quadrature on reference element.

	    	INPUT:

	    	OUTPUT:
	    		pospg:	Position of the Gauss points on reference element
	    		weipg:   Weigth of the Gauss points on reference element

	"""
  pos1::Float64 = sqrt(3/5)
  pospg::Array{Float64}=[-pos1   -pos1;
                             0   -pos1;
                          pos1   -pos1;
                         -pos1       0;
                             0       0;
                          pos1       0;
                         -pos1    pos1;
                             0    pos1;
                          pos1    pos1 ]
  pg1::Float64 = 5/9
  pg2::Float64 = 8/9
  pg3::Float64 = pg1
  weipg::Array{Float64} = [pg1*pg1 pg2*pg1 pg3*pg1 pg1*pg2 pg2*pg2 pg3*pg2 pg1*pg3 pg2*pg3 pg3*pg3]

  return pospg, weipg

end
