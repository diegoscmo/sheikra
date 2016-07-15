################################################################################
#                         Sheikra's Crazy Joe Module                           #
################################################################################
function Crazy_Joe(numturb,numero_chutes,rsf,toplot)

  # Loads chosen RSF
  p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,minx,miny = Load_RSF(string("maps/",rsf,".rsf"));

  # Inclui os dados das regioes para limitar espaço de busca de cada turbina
  # regioes, centroides e reg_turbinas
  @everywhere  include("regioes_chico.jl")

  # Loads Linear Interpolation of the power curve (here for standard density)
  pcurve  = open(readdlm,"others/power_curve.txt")[:,2]
  ctcurve = open(readdlm,"others/CT1.txt")[:,2]

  # Armazena melhor chute e melhor valor - Inicializa tudo em zero.
  melhor_posicao = SharedArray(Float64,2*numturb,1);
  melhor_valor   = SharedArray(Float64,1,1)


  # Inicializa partículas, gbest e pbest
  @sync @parallel for p = 1:numero_chutes

      # Gera as coordenadas para cada particula, sem pegar locais sem vento
        xp = Inicializa_Particula(numturb,A_grid,gridsize,regioes,centroides,reg_turb)

        # Calcula a função objetivo e o valor das restrições
        obj,restr = Fun_Obj(xp',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

        # Se não tiver violação, então vamos ser felizes
        # if restr==0.0
           if obj<melhor_valor[1]
              println("\n Improving... ",obj)
              melhor_valor[1] = obj
              melhor_posicao[:,1] = xp
           end
        #end

  end #p

  #println(melhor_posicao')

  # Escreve a melhor solução
  #writedlm("melhor_crazy_joe.dat",sdata(melhor_posicao))


  #if toplot
  #    Atualiza_Display(melhor_posicao',numturb,melhor_valor[1],1,1,p_grid,gridsize)
  #end

  return vec(melhor_posicao')

end
