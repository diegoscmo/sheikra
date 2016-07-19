
################################################################################
#                                   Sheikra                                    #
################################################################################

# Variável de sistema com o diretório de trabalho
#  cd(ENV["SHEIKRA_DIR"])

# Adiciona processos
if nprocs()==1
   addprocs(7)
end

# Inicializa bibliotecas
using PyPlot
using Distributions
using ProgressMeter

# Inicializa rotinas
include("io.jl")
include("fobj.jl")
include("park.jl")
include("pso.jl")
include("ch.jl")
include("powell2.jl")
include("cmaes.jl")

########## FATOR DE PENALIZACAO #########
# Comentário no cálculo da fobj
#########################################
@everywhere  PENAL = 100.0
@everywhere  PENAL2 = 1000.0
#########################################


# Entrada do RSF e Layout
runname = "teste_cmaes43"    # Para nomear os arquivos de saída
rsf     = "Chicoloma"      # Resource File, deve estar no diretório maps/
layfile = "Chico_ori"      # Layout file, deve estar no diretório maps/


# Carrega RSF e layout
p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,minx,miny = Load_RSF(string("maps/",rsf,".rsf"))
layout = Load_Layout(string("maps/",layfile,".txt"),minx,miny)

# Carrega curvas de potência e empuxo (densidade padrão)
pcurve  = open(readdlm,"others/power_curve.txt")[:,2]
ctcurve = open(readdlm,"others/CT1.txt")[:,2]

# Otimização de Layout
# layout = Opt_PSO(43,rsf,true)
 #layout = Powell(43,2000,1E-5,rsf,true)
 @time layout = CMAES(43,2000,1000,100,rsf,true)

# Cálcula o layout no RSF já carregado
AEP = Calc_Park(layout,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve);

# Grava layout no disco
writedlm(string(runname,"_lay.txt"),[layout[:,1]+minx layout[:,2]+miny])

# Faz o display e grava o resultado no disco
Disp_Res(runname,rsf,AEP)
clf()
Plot_Grid(p_grid);
Plot_WT(layout,gridsize);
