################################################################################
#                         Sheikra's OPT module                                 #
################################################################################

function Opt_PSO(numturb,rsf,toplot)

# Loads chosen RSF
  p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,minx,miny = Load_RSF(string("maps/",rsf,".rsf"));

# Inclui os dados das regioes para limitar espaço de busca de cada turbina
# regioes, centroides e reg_turbinas
@everywhere  include("regioes_chico.jl")

##### Entrada de Parâmetros #####

    # Parâmetros básicos
    nump = 50*numturb              # Número de partículas # 100/turbine
    nit  = 100                     # Número de iterações # 100
    iner = 0                       # 0. w=wi 1. Dim. c/T; 2. Dim. c/P
    wi   = 1.0                     # Inércia inicial
    wf   = 0.1                     # Inércia final
    C1   = 2.0                     # Coeficiente de aceleração 1
    C2   = 2.0                     # Coeficiente de aceleração 2

    # Loads Linear Interpolation of the power curve (here for standard density)
    pcurve  = open(readdlm,"others\\power_curve.txt")[:,2]
    ctcurve = open(readdlm,"others\\CT1.txt")[:,2]


##### PARTICLE SWARM OPTIMIZATION #####

    #Inicializações
    const nd    = 2*numturb                  # numero de dimensoes
    const xp    = zeros(nump,nd)             # Posições
    const ofxp  = zeros(nump)                # Valor da função não-penalizada
    const fxp   = zeros(nump)                # Valor da função penalizada
    const pbest = zeros(nump,nd)             # Posições Partic. Best
    const fpbest= zeros(nump)                # Valor dos Partic. Bbest
    const gbest = zeros(nd)                  # Posição Group Best
    const fgbest= 0.0                        # Valor dos Group Best
    const rxp   = zeros(nump)                # Valor das restrições da partícula
    const fpiorviavel = 0.0                 # Valor da Pior posição viável de todas


    # Inicializa as velocidades usando uma distribuição normal centrada em zero
    d = Normal()
    vp = rand(d,nump,nd)

    # Converte para shared array
    xp   = convert(SharedArray,xp)
    vp   = convert(SharedArray,vp)
    ofxp = convert(SharedArray,ofxp)
    rxp  = convert(SharedArray,rxp)


    # Inicializa partículas, gbest e pbest
    @sync @parallel for p = 1:nump

        # Gera as coordenadas para cada particula, sem pegar locais sem vento
        xp[p,:] = Inicializa_Particula(numturb,A_grid,gridsize,regioes,
                                                   centroides,reg_turb)

        # Calcula a função objetivo e o valor das restrições
        ofxp[p],rxp[p] = Fun_Obj(xp[p,:],numturb,f_grid,A_grid,k_grid,z_grid,
                                       p_grid,numsec,gridsize,pcurve,ctcurve)

    end #p

    # Agora que acabou o prmeiro loop de avaliação das partículas, podemos obter
    # o worsfeasible e dai penalizar os inviáveis
    @inbounds for p=1:nump

        # Verifica se é viável (rxp==0.0) e se o valor do objetivo é maior do que
        # o atual do fpiorviavel
        if rxp[p]==0.0 && ofxp[p]>fpiorviavel
            fpiorviavel = ofxp[p]
        end
    end #p

    # Se não temos nenhuma partícula viável, então
    # consideramos fpiorviavel como o maior
    # valor de objetivo entre as partículas
    if fpiorviavel == 0.0
        fpiorviavel = maximum(ofxp)
    end

    # Para cada partícula inviável, somamos a sua violação com o
    # valor do fpiorviavel
    @inbounds for p=1:nump
        if rxp[p]>0.0
            fxp[p] = fpiorviavel + rxp[p]
        else
            fxp[p] = ofxp[p]
        end
    end #p

    # Busca pelo globalbest
    p_min  = indmin(fxp)
    gbest  = xp[p_min,:]
    fgbest = fxp[p_min]

    # Primeira posição será a melhor inicial
    @inbounds for p=1:nump
       fpbest[p] = fxp[p]
       @inbounds for j=1:nd
            pbest[p,j] = xp[p,j]
       end
    end

    # Cálculo para fazer o índice de dispersão depois
    const disp_ini = 0.0
    @inbounds for p=1:nump
        disp_ini = sqrt(sum(sqrt(gbest.^2) - sqrt(xp[p,:].^2))^2) + disp_ini
    end

    if toplot
      Atualiza_Display(gbest,numturb,fgbest,0,nit,p_grid,gridsize)
    end

    #################################################################
    #                       LOOP PRINCIPAL
    #################################################################
    update = false
    for t = 1:nit

        #println(" Geração ",t) FIXME fisplay desligado

        # Opções de inercia
        if iner == 1                  # Diminui w em t de wi até wf
            w = wi + (t/nit)*(wf-wi)
        elseif iner == 2              # Diminui w em p de wi até wf
            w = wi + (p/nump)*(wf-wi)
        else
            w = wi                    # w fixo
        end #if iner

        @sync @parallel for p = 1:nump
          @inbounds for turb = 1:numturb

                # Localiza x e y da turbina no array
                d1 = 2*numturb-1
                d2 = 2*numturb

                # Atualiza posicao de uma turbina de modo que ela fique na regiao correta
                xp[p,d1:d2],vp[p,d1:d2] = Atualiza_Turbina(xp[p,d1:d2],vp[p,d1:d2],w,C1,C2,
                            pbest[p,d1:d2],gbest[d1:d2],turb,reg_turb,regioes,centroides)

            end #for turb

            # Calcula a função objetivo e o valor das restrições
            ofxp[p],rxp[p] = Fun_Obj(xp[p,:],numturb,f_grid,A_grid,k_grid,z_grid,
                                           p_grid,numsec,gridsize,pcurve,ctcurve)

        end #p

        # Agora que acabou o prmeiro loop de avaliação das partículas, podemos obter
        # o pior viavel e dai penalizar os inviáveis
        fpiorviavel = 0.0
        for p=1:nump

            # Verifica se é viável (rxp==0.0) e se o valor do objetivo é maior do que
            # o atual do fpiorviavel
            if rxp[p]==0.0 && ofxp[p]>fpiorviavel
                fpiorviavel = ofxp[p]
            end

        end #p

        # Se não temos nenhuma partícula viável, então
        # consideramos fpiorviavel como o maior
        # valor de objetivo entre as partículas
        if fpiorviavel == 0.0
            fpiorviavel = maximum(ofxp)
        end

        # Para cada partícula inviável, somamos a sua violação com o
        # valor do fpiorviavel
        @inbounds for p=1:nump
            if rxp[p]>0.0
                fxp[p] = fpiorviavel + rxp[p]
            else
                fxp[p] = ofxp[p]
            end
        end #p

        # Busca pelo globalbest
        p_min_atual  = indmin(fxp)
        bestinho = fxp[p_min_atual]
        if bestinho < fgbest
           gbest  = xp[p_min_atual,:]
           fgbest = bestinho
           update = true
        end

        # Atualização dos melhores valores de cada
        # partícula
        @inbounds for p=1:nump
            if ofxp[p]<fpbest[p]
                #println("\n Melhorou $p $(ofxp[p])") FIXME display desligado
                fpbest[p]  = ofxp[p]
                pbest[p,:] = xp[p,:]
            end
        end

        # índice de dispesão
        const disp_atual=0.0
        @inbounds for p=1:nump
            disp_atual = sqrt(sum(sqrt(gbest.^2) - sqrt(xp[p,:].^2))^2)+disp_atual
        end
        const disp_index = disp_atual/disp_ini

        # Se houverem mudanças no best, dá o display do valor
        if toplot #& update
           Atualiza_Display(gbest,numturb,fgbest,t,nit,p_grid,gridsize)
           update = false
        end #if

    end #for t


##### Display dos resultados finais #####

    if toplot
      Atualiza_Display(gbest,numturb,fgbest,nit,nit,p_grid,gridsize)
    end

    return  Array_Layout(gbest)

end
