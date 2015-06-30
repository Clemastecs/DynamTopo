function Stokes(nx::Int64 = 10, nSteps::Int64 = 20, air::Bool = false)
  #=
           This function solves the 2D-Stokes problem in a square 40x40 of Q2Q1-elements with FEM/Mark-in-Cell.

               INPUT:
                air:     Boolean to apply the layer
                nx:      Mesh dimensions (nx)x(ny)-elements
                nSteps:  Number of the time steps

            OUTPUT:
                One data files on ./results/file.data
                One figure of the last map of the particles on ./results/figs/fig.png

            ANNOTATIONS:
                Sphere density:  1150
                Mantle density: 1420
                Air density: 1
                Sphere viscosity: 69000
                Mantle viscosity: 50
                Air viscosity: Sticky air
  =#

  ## Initial Data ##
  ##################
  i::Int64 = 0

  visc::Array{Float64} = [50 69000 50] # [mantle viscosity, sphere viscosity, air viscosity]
  rho::Array{Float64} = [1420 1150 0] # [mantle density, sphere desnity, air density]
  pto::Array{Float64} = [ 20 5 ] # initial center position of the body
  radius::Float64 = 3 # initial radius of the body
  ppe::Int64 = 30 # particles per element

  # domain [x1,x2]x[y1,y2]
  x1::Int64 = 0
  x2::Int64 = 40
  y1::Int64 = x1
  y2::Int64 = x2

  # mesh dimensions (nx)x(ny)-elements
  ny::Int64 = nx

  nxx::Int64 = 2*nx # data for creation/ploting the mesh
  nyy::Int64 = nxx

  # element Q2Q1
  nen::Int64 = 9 # num nodes of velocity element
  nenP::Int64 = 4 # num nodes of pressure element

  # adding air layer
  gap::Float64 = 0.

   if air == true
          pto[2] = pto[2] - 2.5
          ngap::Int64 = itrunc((1/4)*ny) # number of the elements of the layer in y direction
          gap = (1/4)*y2 # y dimensions of the layer
          ny = ny + ngap
          y2 = y2 + gap
  end

  ## Velocity/Pressure coords and connectivity matrixs ##
  #######################################################

  X::Array{Float64} = []
  T::Array{Float64} = []
  (X,T) = mvelnod(x1,x2,y1,y2,nx,ny,nen)

  XP::Array{Float64} = []
  TP::Array{Float64} = []
  (XP,TP) = mpresnod(X,T,nx,ny,nenP)

  # Remesh in case of air layer
  if air == true
          remesh!(XP,X,nx,ngap,gap)
          y2 = x2 # recover the dimensions due to the remesh
  end

  # radius to choose the mean of density and viscosity
  Te::Array{Float64} = vec(T[end,:])
  Ye::Array{Float64} = X[Te,2]
  rad::Float64 = Ye[3] - Ye[2]
  sradius::Float64 = rad / 2

  ## Initial set of particles ##
  ##############################

  nPar::Int64 = 0 # num of particles
  par::Array{Float64} = [] # grid of particles
  (nPar,par) = setparticles(x1,x2,y2,nx,ny,ppe,radius,visc,rho,pto,air,gap)

  if air == true
       isurface::Int64 = findfirst(par[:,3] .== 3) # mark the surface for the par
       surface::Array{Float64} = zeros(itrunc(sqrt(ppe)*nx)) # height of the set of particles surface
       file2 = open("./results/surface_elem"*string(nx)*".dat","w")# createfile
       close(file2)
  else
          file1 = open("./results/sigmazz_elem"*string(nx)*".dat","w")# openfile
          close(file1)
  end

  ## Number of Nodes and Unknowns ##

  nNodesVel::Int64 = size(X,1)
  nUnkVel::Int64 = 2*nNodesVel
  nUnkPre::Int64 = size(XP,1)

  ## Boundary conditions: impose BC via Lagrange ##

  Accd::Array{Float64} = []
  bccd::Array{Float64} = []
  (Accd, bccd) = bcfreeslip(X,nUnkVel)
  nDirichletBC::Int64 = size(Accd,1)

  ## Computation of (K,G)(v,p)=(f) system ##
  ##########################################

  nUnkPre = nUnkPre-1 # necessary to impose a 0-row
  velo::Array{Float64} = [] # velocity vector
  pres::Array{Float64} = [] # pressure vector
  K::SparseMatrixCSC{Float64,Int64} = spzeros(1,1)
  G::SparseMatrixCSC{Float64,Int64} = spzeros(1,1)
  f::Array{Float64} = []
  M::Array{Float64} = zeros(nUnkPre, nUnkPre)
  Atot::Array{Float64} = []
  btot::Array{Float64} = []
  aux::Array{Float64} = []
  sigmazz::Array{Float64} = zeros(2*nx) # stress_zz per columns (2*nx: number of gauss points in compStress)

  timestep::Float64 = 0 # auxiliar variable to show the completed percentantge of the computation
  print_with_color(:red,"Elements = "*string(nx*ny)*"\n")
  print_with_color(:red,"Particles = "*string(nPar)"\n")
  time::Float64 = 0

  PyPlot.figure(1) # start plotting

  for theStep = 1:nSteps
         tic()
         # draw material map
         PyPlot.clf()
         PyPlot.subplot(1,2,1,aspect=1)
         plotpar(par); # plot the grid particles with materials
         PyPlot.plot(XP[:,1],XP[:,2],"ko", alpha=1) # plot nodes
         PyPlot.axis([x1, x2, y1, y2])
         PyPlot.title("Materials")

         # compute the matrix components
         (K, G, f) = totalmat(X,T,XP,TP,par,sradius)

         # boundary conditions for pressure
         G = G[ [1:end-1], :] # Erase last row

         # global matrix generation with Dirichlet BC
         Atot = [K Accd' G'; Accd zeros(nDirichletBC,nDirichletBC) zeros(nDirichletBC,nUnkPre); G zeros(nUnkPre,nDirichletBC) M]
         btot = [f; bccd ; zeros(nUnkPre,1)]

         aux = Atot\btot # solve the system

         velo = reshape(aux[1:nUnkVel],2,nNodesVel)' # velocity solution

         pres = [aux[(nUnkVel+nDirichletBC+1):end]; 0] # pressure solution

         # draw Velocity Field
         PyPlot.subplot(1,2,2,aspect=1)
         PyPlot.quiver(X[:,1],X[:,2],velo[:,1],velo[:,2],1)
         PyPlot.title("Velocity Field")

         PyPlot.savefig("./results/figs/map"*string(nx)*".png",dpi=250)

         if air == true
             rank2m = (isurface-2*itrunc(sqrt(ppe)*nx)):(isurface-itrunc(sqrt(ppe)*nx)-1)
             rank1m = (isurface-itrunc(sqrt(ppe)*nx)):(isurface-1)
             rank1a = (isurface):(isurface+itrunc(sqrt(ppe)*nx)-1)
             rank2a = (isurface+itrunc(sqrt(ppe)*nx)):(isurface+2*itrunc(sqrt(ppe)*nx)-1)
             surfacemean = ( par[rank2m,2] .+ par[rank1m,2] .+ par[rank2a,2] .+ par[rank1a,2]  )./4
             surface = [surface surfacemean ] # mean of the 2 first air/mantle rows
             file2 = open("./results/surface_elem"*string(nx)*".dat","a")# openfile
             writedlm(file2, surfacemean',"\t") # write in file
             close(file2)
         else
             # compute stresses on the top boundary
             sigmazze = compstress(velo,pres,X,T,TP,par,nx,sradius);
             sigmazz = [sigmazz sigmazze] # the 1st row are zeros
             file1 = open("./results/sigmazz_elem"*string(nx)*".dat","a")
             writedlm(file1, sigmazze',"\t") # write in file
             close(file1)
         end

         # move particles
         updateparticles!(par,velo,nx,ny,x1,x2,y1,y2,rad)

         # show the computational time
         timestep = toq()
         time = time + timestep
         timestep = round(nSteps*timestep-theStep*timestep,2)
         print_with_color(:black,"Time: "*string(timestep)*" s (Total:"*string(time)*")\n")


  end


  ## Post-process ##
  ##################

  if air == true
        plotlayer(x1,x2,surface,ppe,nx,nSteps)
  else
        plotfix(sigmazz,nx,nSteps,rho,x2)
  end

  ## End of the the function ##
  #############################

  print_with_color(:green,"Successful computation.\n")

end
