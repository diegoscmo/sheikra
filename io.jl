
################################################################################
#                         Data Inputs and Handling                             #
################################################################################

# Read Resource File into several arrays downsizes x,y
 function Load_RSF(file)

   # Opens RSF file handler
   fhand=open(readdlm,file)

   # Reduces the x,y dimensions
   xmin=minimum(fhand[:,1])
   ymin=minimum(fhand[:,2])
   fhand[:,1]=fhand[:,1]-xmin
   fhand[:,2]=fhand[:,2]-ymin

   # Generates the gridsize
   gridsize=abs(fhand[2,1]-fhand[1,1])

   # Loads number of sectors
   numsec=convert(Int64,fhand[1,8])

   # Defines size of x and y uniques and a x y grid list
   stx = size(unique(fhand[:,1]),1)
   sty = size(unique(fhand[:,2]),1)
   xylist=[fhand[:,1] fhand[:,2]]/gridsize

   # Grids for power, height, weibulls and frequency
   p_grid = zeros(stx,sty)
   z_grid = zeros(stx,sty)
   A_grid = zeros(stx,sty,numsec)
   k_grid = zeros(stx,sty,numsec)
   f_grid = zeros(stx,sty,numsec)

   # Gets all of them for each x and y in the list
   for i=1:size(fhand,1)
     a=convert(Int64,xylist[i,1]+1)
     b=convert(Int64,xylist[i,2]+1)
     p_grid[a,b] = fhand[i,7]
     z_grid[a,b] = fhand[i,3]

     # For each sector too
     for s=1:numsec
       f_grid[a,b,s] = fhand[i,(3*s+6)]
       A_grid[a,b,s] = fhand[i,(3*s+7)]/10.
       k_grid[a,b,s] = fhand[i,(3*s+8)]/100.
     end # for s

     f_sum = sum(f_grid[a,b,:])

     for s=1:numsec
       f_grid[a,b,s]=f_grid[a,b,s]/f_sum
     end

   end # for i

   return p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,xmin,ymin
 end # function Load_RSF

 # Read Park layout into a matrix and downsizes x,y
  function Load_Layout(Lay,xmin,ymin)
    L=open(readdlm,Lay)

    L[:,1]=L[:,1]-xmin
    L[:,2]=L[:,2]-ymin

    return  L
  end #Load_Map



# Creates a plot from grid
function Plot_Grid(any_grid)

  # New window and plot
  fig = plot()
  any_grid = any_grid'
  p = imshow(any_grid)

  # clamp the color limits and colorbar
  clim(minimum(any_grid),maximum(any_grid))
  colorbar()
  title("")
  ax = gca()
  ax[:set_ylim]([0.0,size(any_grid,1)])
  ax[:set_xlim]([0.0,size(any_grid,2)])

end

# Plots WT in the already plotted map
 function Plot_WT(apos,gridsize)

   apos = apos/gridsize
   plot(apos[:,1],apos[:,2],"w*",markersize=15)

 end #Plot_WT

#Writes to file and outputs
function Disp_Res(layf,rsf,AEP)

  # open file to write
  f = open(string(layf,"_",rsf,"_res.txt"),"w")

  println("\n Computation done:")
  println("\n      X[m]        Y[m]       AEP[MWh]     c/wake       Eff[%]")

  write(f," X[m]   Y[m]   AEP[MWh] c/wake Eff[%]\n")

  for i=1:size(AEP,1)
    for j=1:size(AEP,2)

      @printf("%12.2f", AEP[i,j])
      s = @sprintf("%.2f  ", AEP[i,j])
      write(f,s)

    end #for j

    @printf("\n")
    write(f,"\n")

  end #for i

  close(f)
end #Disp_Res


@everywhere function Interpola1D(dados,ponto)


    xl = convert(Int64,floor(ponto))
    xu = convert(Int64,ceil(ponto))

    # Testa por um valor nulo em xl
    if xl==0
        return 0.0
    end

    # duvido !
    if xl==xu
        return dados[xl]
    end

    if xu>length(dados)
        return dados[end]
    end

    if xl<0
        return dados[1]
    end


    # Mais provável...
    fl = dados[xl]
    fu = dados[xu]

    saida = fl + ((fu-fl)/(xu-xl))*(ponto-xl)

    return saida

end
