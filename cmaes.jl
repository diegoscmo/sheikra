#
# Baseado em
# [1] Hansen, N. (2011). The CMA Evloution Strategy: A Tutorial
#
function CMAES(numturb,nc,maxiter,stopestag,rsf,toplot,xmedia=[])

  # numturb   número de turbinas
  # nc        número de chutes do Crazy_Joe
  # maxiter   máximo de iterações
  # stopestag máximo de iterações estagnado
  # rsf       nome do resource files
  # toplot    to plot, or not to plot...


  # Carrega o resource file
  p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,minx,miny = Load_RSF(string("maps/",rsf,".rsf"));

  # Inclui os dados das regioes para limitar espaço de busca de cada turbina
  # regioes, centroides e reg_turbinas
  @everywhere include("regioes_chico_esc.jl")

  # Loads Linear Interpolation of the power curve (here for standard density)
  pcurve  = open(readdlm,"others/power_curve.txt")[:,2]
  ctcurve = open(readdlm,"others/CT1.txt")[:,2]

  # Número de variáveis / dimensão do problema
  const ndim = 2*numturb

  # Define fora do If para ficar visível
  melhor_posicao = zeros(ndim,1) 
  melhor_valor   = zeros(1,1)    SharedArray(Float64,1,1)

  if length(xmedia)==0

    # Converte para shared array dentro do if...
    melhor_posicao = convert(SharedArray,melhor_posicao)
    melhor_valor   = convert(SharedArray,melhor_valor)
    
    # Inicializa partículas, gbest e pbest
    @sync @parallel for k = 1:nc

        # Gera as coordenadas para cada particula, sem pegar locais sem vento
          xp = Inicializa_Particula(numturb,A_grid,gridsize,regioes,centroides,reg_turb)

          # Calcula a função objetivo e o valor das restrições
          obj = Fun_Obj_Powell(xp,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)

          # Se melhorou, então guarda
          if obj<melhor_valor[1]
                println(" Crazy_Joe::Improving... ",obj)
                melhor_valor[1] = obj
                for j=1:ndim
                    melhor_posicao[j,1] = xp[j]
                end
          end

    end #p

    # Melhor chute é o nosso ponto incial no CMAES
    xmedia = vec(sdata(melhor_posicao)')

    # Libera a memória ...
    @everywhere gc()

  else

    # Recebemos um ponto inicial...
    # Vamos verificar a consistência
    if !size(xmedia,1)==mum2
        error("\n CMAES::Ponto inicial não tem a dimensão correta. ",size(xmedia))
     end
  end

  # Desvio padrao das coordenadas (step-size inicial) FIXME, estudar valor
  sigma = 50.0

  # Parámetros de seleção:
    # Tamanho da populacao >=2 [1-eq44]
  lambda = convert(Int64,4.0+floor(3.0*log(ndim)))  #FIXME isso pode ser um input

    # Numero de pts selec p/ recombinação [1-eq44]
  mu = convert(Int64,floor(lambda/2))               #FIXME isto pode ser outro input

    # mu x 1 pesos de recombinação normalizado [1-eq45]
  pesos = log(mu+1.0/2.0) - log(1.0:mu)
  pesos = pesos/sum(pesos)

    # Variancia efetiva do tamanho de mu [1-eq8]
  mueff = sum(pesos)^2/sum(pesos.^2)

  # Parametros de Adaptação:
    # Const de tempo p/acumulacao de C [1-eq47] !
  cc  = (4.0+mueff/ndim)/(ndim+4.0+2.0*mueff/ndim)

    # Const de tempo p/acumulacao do sigma [1-eq46] !
  cs  = (mueff+2.0)/(ndim+mueff+5.0)

    # Taxa de aprendizado para atualização rank-1 [1-eq48]
  c1  = 2.0 / ((ndim+1.3)^2 + mueff)

    # Atualização rank-mu [1-eq49 au=2]
  cmu = min(1.0-c1, 2.0*(mueff-2.0+1.0/mueff)/((ndim+2.0)^2.0+mueff))

    # Amortecimento para sigma [1-eq46]
  damps = 1.0 + 2.0*max(0.0,sqrt((mueff-1.0)/(ndim+1.0))-1.0) + cs

    # Expectativa de ||N(0,I)|| == norm(randn(N,1)) [1-pg25]
  chiN = ndim^0.5*(1.0-1.0/(4.0*ndim)+1.0/(21.0*ndim^2))

  # Aloca matrizes principais do CMA
  B = eye(ndim)                     # Matriz B define o sistema de coordenadas
  D = eye(ndim)                     # Matriz diagonal D, define a escala
  C = B*D*(B*D)'                    # Matriz de Covariância

  # Aloca demais variáveis
  pc = zeros(ndim)                  # Caminho de evolução para C
  ps = zeros(ndim)                  # Caminho de evolução para sigma
  arz = zeros(lambda,ndim)          # Cada linha é uma solução
  arx = zeros(lambda,ndim)
  arfitness = zeros(lambda)
  arindice  = zeros(lambda)

  iter       = 0                    # Contador de iterações
  eigeneval  = 0                    # do artigo original, removerei no futuro
  contaeval  = 0                    # Número de avaliações da Fun_Obj
  contaestag = 0                    # Número de iterações estagnado

  # Para processamento paralelo
  #arz = convert(SharedArray,arz)
  #arx = convert(SharedArray,arx)
  #arfitness = convert(SharedArray,arfitness)


  # Primeiro display
  if toplot
     Atualiza_Display(xmedia,numturb,melhor_valor[1],iter,maxiter,p_grid,gridsize)
  end

  ######### Loop do CMA-ES ############

  while true

    # Incrementa iterações
    iter += 1

    # Gera a avalia descendentes de lambda
    #@sync @parallel
    for k = 1:lambda
      #FIXME quando o código é paralelizado a minimização não acontece

      # Distribuição normal
      arz[k,:] = randn(ndim,1)

      # Adiciona mutação
      arx[k,:] = xmedia + sigma * (B*D*arz[k,:]')

      # Chamada da função objetivo
      arfitness[k] = Fun_Obj_Powell(arx[k,:]',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)

    end # for k

    # Incrementa avaliações de fobj
    contaeval += lambda

    # Organiza pontos pelo fitness
    arindice  = sortperm(arfitness)
    arfitness = arfitness[arindice]

    # Salva media antiga computa a nova poderada em xmedia
    xantig = copy(xmedia)
    xmedia = arx[arindice[1:mu],:]'*pesos           # Recombinacao dos selecionados
    zmedia = arz[arindice[1:mu],:]'*pesos           # == D^-1*B'*(xmedia-xantig)/sigma

    # Acumulação: Atualiza os caminhos de evolução [1-eq40, 1-pg25 e 1-eq42]
    ps = (1.0-cs)*ps + (sqrt(cs*(2.0-cs)*mueff))*(B*zmedia)

    hs = norm(ps)/sqrt(1.0-(1.0-cs)^(2*contaeval/lambda))/chiN < 1.4+2.0/(ndim+1.0)

    pc = (1.0-cc)*pc + hs*(sqrt(cc*(2.0-cc)*mueff)/sigma) * (xmedia-xantig)

    # Adapta a matriz de Covariância C
    C = ( (1.0-c1-cmu)*C + (1.0-hs)*cc*(2.0-cc) * C   # Mantido da iteração anterior
    + c1 * pc*pc'                                     # Atualização Rank-Um
    + cmu * (B*D*arz[arindice[1:mu],:]')              # Atualização Rank-mu
    * diagm(pesos)*(B*D*arz[arindice[1:mu],:]')')     # FIXME dá pra deixar mais rápido isso..

    # Adapta o tamanho de passo sigma [1-eq41]
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN-1.0))

    # Atualiza B e D utilizando C (Para chegar em O(N^2)) FIXME esse eigenval poderia sair
    if contaeval-eigeneval > lambda/(c1+cmu)/ndim/10.0

      eigeneval = contaeval
      C = triu(C) + triu(C,1)'                   # Força a simetria
      (tmp,B) = eig(C)                           # Decomposição de autovalores e avetores
      diagD   = sqrt(tmp)
      D = diagm(diagD)                           # D contém os desvios padrões

    end #if

    # Salva o melhor valor encontrado até então e dá o display
    if arfitness[1] < melhor_valor[1]

      melhor_valor = copy(arfitness[1])
      melhor_posicao = copy(arx[1,:])
      contaestag = 0
      if toplot
         Atualiza_Display(arx[1,:],numturb,melhor_valor[1],iter,maxiter,p_grid,gridsize)
      end

      # Como o metodo pode demorar, grava layout no disco
      layout = Array_Layout(melhor_posicao')
      writedlm(string("cmaes_lay.txt"),[layout[:,1]+minx layout[:,2]+miny])

      # Caso não achar, incrementa o contador de estagnação e sai se atingir stopestag
    else

      contaestag += 1
      if contaestag >= stopestag
        println("\nFim: Atingiu o número de iterações estagnado sem melhorias.")
        break
      end #contaestag

    end # if arfitness

    # Display a cada 10 iterações
    if iter%10 == 0
      @printf("iter: %d \t fitness: %.2f \t\n",iter, melhor_valor[1])
    end

    # Finaliza se passar o número máximo de iterações
    if iter > maxiter
      println("\nFim: Atingiu o máximo de iterações.")
      break
    end

  end #while

  # Adquire o melhor layout
  layout = Array_Layout(melhor_posicao')

  # Imprime se há restrições
  Verifica_Penalizacao(melhor_posicao',numturb,f_grid,A_grid,k_grid,
                               z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

  # Imprime resultado
  @printf("Fim. \t fitness: %.2f \t\n", melhor_valor[1])
  # e retorna com o layout encontrado
  return layout

end

######### Funções para testes #############

function Fun_ObjCMA(x)

  # Parabólica
  #valor = x[1]^2 + x[2]^2

  # Vale de Rosenbrock
  #valor =  100.0*(x[2]-x[1]^2)^2 + (1.0-x[1])^2

  # Rastringin
  valor = 20.0 + x[1]^2 - 10.0*cos(2.0*pi*x[1]) + x[2]-10.0*cos(2.0*pi*x[2])

  # Ackley
  valor = ( -20.0*exp(-0.2*sqrt(0.5*(x[1]^2+x[2]^2)))
  - exp(0.5*(cos(2.0*pi*x[1])+cos(2.0*pi*x[2])))
  + exp(1) + 20.0 )

  return valor
end
