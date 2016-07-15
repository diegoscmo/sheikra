
################################################################################
#                                Park Calculation                              #
################################################################################

# Park AEP calculation
@everywhere function Calc_Park(layout,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

  ### Future turbine inputs ####
  # assuming here same turbine for all layout
  const rot_rad = 55.0                   # Rotor radius
  const hub_hei = 120.0                  # Hub Height
  const wake_dk = 0.11                   # Wake decay (wake spread angle)
                                         # Air density alrays the same

  ### Future calculation inputs ####
  const vini    = 0.5                    # Initial velocity #0.5
  const vstep   = 1.0                    # Step for loop through velocities #1.0
  const vfin    = 20.5                   # Final velocity #30.5 ##MODIFICADO
  const nangles = 31                     # Number of angles/sector (for wake) #31, can't be 1  ##MODIFICADO


  # Velocities
  const vels    = collect(vini:vstep:vfin)  # Velocities to be calculated
  const numvel  = size(vels,1)

  # Inicializes AEP
  const numturb = size(layout,1)
  const AEP     = zeros(Float64,numturb+1,5)

  # Calculates wake for the current layout
  const wakemtx = Jensen_Wake(layout,numturb,A_grid,z_grid,gridsize,numsec,nangles,vels,rot_rad,hub_hei,wake_dk,ctcurve)

  # Loops through each turbine and sector
  @inbounds for t=1:numturb
    @inbounds for sec=1:numsec

      # Acquires data from the grid
      const f_dir = Get_Grid(layout[t,:],f_grid[:,:,sec],gridsize)
      const A_dir = Get_Grid(layout[t,:],A_grid[:,:,sec],gridsize)
      const k_dir = Get_Grid(layout[t,:],k_grid[:,:,sec],gridsize)

      # Calculates Weibull distribution
      const dist_vel = Gen_Weibull(A_dir,k_dir,vini,vstep,vfin)

      @inbounds for v=1:numvel
        # Computes the frequency of wind for each velocity
        const veloc = dist_vel[v,1]
        const f_vel = dist_vel[v,2]
        const freq  = f_dir*f_vel*365.0*24.0/1E3

        # Acumulates AEP without wake
        AEP[t,3] += freq*Interpola1D(pcurve,veloc)

        # For wake calculation:
        @inbounds for ang = 1:nangles
          # Computes wake_velocity for each angle
          const wake_vel = veloc*(1.0 - wakemtx[t,v,sec,ang])

          # Acumulates AEP with wake (divided by number of angles)
          AEP[t,4] += freq*Interpola1D(pcurve,wake_vel)/nangles

        end #for ang

      end # for v
    end # for sec

    # Includes X, Y and Eff % to the final display
    AEP[t,1:2]=layout[t,:]
    AEP[t,5]=100*AEP[t,4]/AEP[t,3]

  end # for t
    AEP[numturb+1,3]=sum(AEP[:,3])
    AEP[numturb+1,4]=sum(AEP[:,4])
    AEP[numturb+1,5]=100*AEP[numturb+1,4]/AEP[numturb+1,3]

  return AEP
end # function Calc_Park


# Calculates velocity deficits for a layout using Jensen Wake model
@everywhere function Jensen_Wake(layout,numturb,A_grid,z_grid,gridsize,numsec,nangles,vels,rot_rad,hub_hei,wake_dk,ctcurve)

  # Rotor area calculation
  const rot_area = pi*rot_rad^2

  # Loads Linearly Interpolated CT curve (here for standard density)
  #CT = Load_Curve("others\\CT1.txt")

  # Includes turbine position to the layout
  const ord_layout = [layout 1:1:numturb]

  # Initializes wake matrix
  const numvel   = size(vels,1)
  const wakemtx  = zeros(Float64,numturb,numvel,numsec,nangles)

  # Loads height of each WT
  const z = Array(Float64,numturb)
  @inbounds for i=1:numturb
    z[i] = Get_Grid(ord_layout[i,1:2],z_grid,gridsize)
  end

  # Step angles
  const angstep = 360/(numsec*(nangles-1))

  # Iteration through sectors and angles (eg 1 (north), -15° to +15°)
  @inbounds for sec=1:numsec
    @inbounds for ang=1:nangles

      # Calculates angle tetha for each angle to be checked on the sector (north = -15 to +15)
      const theta = ( ((sec-1)*360/numsec) - 360.0/numsec/2.0 + (ang-1)*angstep )*pi/180.0

      # Includes turbine position to the layout
      const ord_layout = [layout 1:1:numturb]
      # Organizes the turbines according to the wind direction, upwinds first
      ord_layout = Order_Layout(ord_layout,numturb,theta)

      # Loop through each ith turbine
      @inbounds for i = 1:numturb

        # Only jth turbines downind of ith will be looped through
        @inbounds for j = (i+1):numturb
          # Checks the downwind distance from ith to the jth turbine
          const vector = ord_layout[i,1:2] - ord_layout[j,1:2]
          const ydist  = abs(vector[1]*sin(theta) + vector[2]*cos(theta))

          # Checks distance in x
          const xdist = abs(vector[1]*cos(theta) + vector[2]*sin(theta))

          # Identifies the turbine numbers for localization
          const tposi = convert(Int64,ord_layout[i,3])
          const tposj = convert(Int64,ord_layout[j,3])

          # Computes new xdist with x and z pitagoras
          const xdist = sqrt( xdist^2 + abs(z[tposi] - z[tposj])^2 )

          # Computes the Wake circle radius
          const r_wake = rot_rad + wake_dk*ydist
          # Inicializes wake multiplier
          const wake_mult = 0.0
          # If the turbine is 100% inside the wake, normal wake_mult calculation
          if xdist <= r_wake - rot_rad
            wake_mult = 1.0 / (1.0 + wake_dk*ydist/(rot_rad) )^2

            # If the turbine is partially inside the wake cone
          elseif xdist < r_wake + rot_rad

          # Calculates the intersection shadow area // interference of the circles
          area_mult = 0.5*(r_wake^2*(2.0*acos((r_wake^2+xdist^2-rot_rad^2)/(2.0*r_wake*xdist))  -
                           sin(2.0*acos((r_wake^2+xdist^2-rot_rad^2)/(2.0*r_wake*xdist)))) +
                           rot_rad^2*(2.0*acos((rot_rad^2+xdist^2-r_wake^2)/(2.0*rot_rad*xdist)) -
                            sin(2.0*acos((rot_rad^2+xdist^2-r_wake^2)/(2.0*r_wake*xdist)))))/rot_area

          # Calculates the intersection shadow area // interference of the squares
          #area_mult = ((r_wake + rot_rad)-xdist)  / (2*rot_rad)

            # And the area weighted wake multiplier
          wake_mult = area_mult / (1.0 + wake_dk*ydist/(rot_rad) )^2

          end #if xdist

          # If a wake multiplier was sucessfully computed
          if wake_mult > 0.0
            @inbounds for v=1:numvel
              # Velocity of the upwind turbine considering upwind wake
              vup = vels[v]  * (1.0 - sqrt(wakemtx[tposi,v,sec,ang]) )
              #
              # Velocity deficit relating to the upwind turbine speed
              vel_def = (1.0 - sqrt(1.0 - Interpola1D(ctcurve,vup))*(vup/vels[v]) ) * wake_mult

              # Acumulates the squared sum for each velocity
              wakemtx[tposj,v,sec,ang] += vel_def^2

            end # for v
          end # if wake_mult

        end # for j
      end # for t

    end # for ang
  end # for sec

  # Only the squares were computed, now everything needs to be rooted
  wakemtx = sqrt(wakemtx)
  return wakemtx
end #function

# Organizes the turbines according to the sector angle, auxiliary function for wake calculation
@everywhere function Order_Layout(ord_layout,numturb,theta)

  # Finds the extremities of the layout
  xmin = minimum(ord_layout[:,1])
  ymin = minimum(ord_layout[:,2])
  xmax = maximum(ord_layout[:,1])
  ymax = maximum(ord_layout[:,2])

  # and computes the center
  xcent = (xmax+xmin)/2.0
  ycent = (ymax+ymin)/2.0

  # Then create a point outside these extremities according to the angle
  allmax = max((ymax-ymin),(xmax-xmin))
  farpoint = [allmax*sin(theta)+xcent allmax*cos(theta)+ycent]

  # A new column is added to compute distances
  ord_layout = [ord_layout zeros(numturb)]

  # Then the distance of the farpoint to each turbine in theta is found
  for t = 1:numturb
    dist = farpoint - ord_layout[t,1:2]
    ord_layout[t,4] = norm(dist[1]*sin(theta) + dist[2]*cos(theta))

  end # for t

  # And the list gets sorted by distance, and the distances column is excluded
  ord_layout = sortrows(ord_layout, lt=(x,y)->isless(x[4],y[4]))[:,1:3]

  # Ele é passado por referência para esta rotina, portanto não precisamos dar um return
  return ord_layout
end #function Order_Layout


# Returns Weibull distribution from 3 to 20 m/s wind speed 1m/s bins
@everywhere function Gen_Weibull(A,k,vini,vstep,vfin)

  # Number of velocities
  const si=size(collect(vini:vstep:vfin),1)

  # Initializes distribution array
  const W = Array(Float64,convert(Int64,(vfin-vini))+1,2)

  # Loops through speeds
  const cont = 1
  @inbounds for x = vini:vstep:vfin

    const f =  ( (k/A) * ((x/A)^(k-1)) * e^-(x/A)^k )*vstep
    W[cont,:] = [x   f]

    # In case of NaN (i.e.: A = 0 inputs), W = 0
    if isnan(W[cont,2])
      W[cont,2] = 0.0
    end

    cont += 1
  end

  return W
end #function Gen_Weibull

## Searches for a interpolated property in the Grid
@everywhere function Get_Grid(pos,Mg,grd)

  # transforms the x,y positions into the arrays
  x,y=pos/grd + 1

  # Caso apareça uma variavel negativa.
  if x<1
    x=1
  end
  if y<1
    y=1
  end

  # to avoid bugs, negative coordinates are zero
  if x<0.0 || y<0.0 || x>size(Mg,1) || y>size(Mg,2)
    int_prop = 0.0
  else
    # transforms into array integers
    xl = convert(Int64,floor(x))
    xu = convert(Int64,ceil(x))
    yl = convert(Int64,floor(y))
    yu = convert(Int64,ceil(y))

    # If it is exactly in a node
    if xl==xu && yl==yu
      int_prop = Mg[xl,yl]

    # If it is in the middle of the four nodes
    elseif xl!=xu && yl!=yu

      # Bilinear interpolation
      A = [1 xl yl xl*yl
           1 xl yu xl*yl
           1 xu yl xu*yl
           1 xu yu xu*yu]

      Pts = convert(Array{Int64,2},A[:,2:3])
      B   = zeros(4,1)

      for i=1:4
        B[i,:]=Mg[Pts[i,1],Pts[i,2]]
      end

      if minimum(B) == 0.0
        int_prop = 0.0
      else
        C = A\B
        int_prop =  C[1]+C[2]*x+C[3]*y+C[4]*x*y
      end

    # Interpola em x
    elseif xl!=xu && yl==yu
            int_prop = Mg[xl,yl] + ((Mg[xu,yl]-Mg[xl,yl])/(xu-xl))*(x-xl)
    else
            # Interpola em y
            int_prop = Mg[xl,yl] + ((Mg[xu,yu]-Mg[xl,yl])/(yu-yl))*(y-yl)
    end

  end
  return int_prop
end
