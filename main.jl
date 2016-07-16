
################################################################################
#                                   Sheikra                                    #
################################################################################

#  cd(ENV["SHEIKRA_DIR"])  #Working directory for library and files

if nprocs()==1
   addprocs(4)
end

# Initializes libraries
using PyPlot
using Distributions
using ProgressMeter

include("io.jl")
include("fobj.jl")
include("park.jl")
include("pso.jl")
include("ch.jl")
include("powell2.jl")

########## FATOR DE PENALIZACAO #########
# Comentário no cálculo da fobj
#########################################
@everywhere  PENAL = 100.0
#########################################

########################## Inputs #############################
# RSF and Layout files input
rsf  = "Chicoloma"       # ReStource File, under maps directory
layfile = "Chico_ori" # Layout file, under maps directory

######################## Single Run ###########################

# Loads maps
p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,minx,miny=Load_RSF(string("maps/",rsf,".rsf")) ;
layout=Load_Layout(string("maps/",layfile,".txt"),minx,miny)

# Loads Linear Interpolation of the power curve (here for standard density)
pcurve  = open(readdlm,"others/power_curve.txt")[:,2]
ctcurve = open(readdlm,"others/CT1.txt")[:,2]

# Layout optimization
# layout = Opt_PSO(43,rsf,true)
  layout = Powell(12,1000,1E-5,rsf,true)

# Single park calculation
AEP = Calc_Park(layout,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve);

# Grava no disco
writedlm("resultado.dat",layout)

# Result output to console and file
Disp_Res(layfile,rsf,AEP)
