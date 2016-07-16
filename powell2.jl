#
# Baseado em
# http://www.aip.de/groups/soe/local/numres/bookcpdf/c10-5.pdf, página 419
#
function Powell(numturb,nc,tol,rsf,toplot)

  # Tolerância para a parada (variação do objetivo)
  const ftol = tol

  # Loads chosen RSF
  p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,minx,miny = Load_RSF(string("maps/",rsf,".rsf"));

  # Inclui os dados das regioes para limitar espaço de busca de cada turbina
  # regioes, centroides e reg_turbinas
  @everywhere include("regioes_chico.jl")

  # Loads Linear Interpolation of the power curve (here for standard density)
  pcurve  = open(readdlm,"others/power_curve.txt")[:,2]
  ctcurve = open(readdlm,"others/CT1.txt")[:,2]

  # Usamos bastante este cara aqui
  const num2 = 2*numturb

  # Faz um início aleatório com nc tentativas. SE nc==1, então equivale ao
  # uso do Inicializa_Particulas.
  # Armazena melhor chute e melhor valor - Inicializa tudo em zero.
  melhor_posicao = SharedArray(Float64,num2,1);
  melhor_valor   = SharedArray(Float64,1,1)

  # Inicializa partículas, gbest e pbest
  @sync @parallel for k = 1:nc

      # Gera as coordenadas para cada particula, sem pegar locais sem vento
        xp = Inicializa_Particula(numturb,A_grid,gridsize,regioes,centroides,reg_turb)

        # Calcula a função objetivo e o valor das restrições
        obj,restr = Fun_Obj(xp',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

        # Se melhorou, então guarda
        if obj<melhor_valor[1]
              println(" Crazy_Joe::Improving... ",obj)
              melhor_valor[1] = obj
              melhor_posicao[:,1] = xp
        end


  end #p

  # Melhor chute é o nosso ponto incial no Powell
  const p = vec(sdata(melhor_posicao)')

  # Libera a memória ...
  @everywhere gc()

  # Já que o resto não está paralelizado, vamos desconectar os processos
  np = nprocs()
  if np>1
     for i=np:-1:2
         rmprocs(i)
    end
  end

  # Gera a base inicial de busca
  const L = eye(num2,num2)

  # Uma direação de busca (coluna de L)
  const xit = Array(Float64,num2)

  # Indicam a direção e o valor da maior mudança no valor da função objetivo
  const ibig = 0
  const del = 0.0

  # Vamos abrir um arquivo para escrita, pois o método é lento
  # e vale a pena acompanhar o andamento
  arquivo = open("convergencia_powell.txt","w")

  # Faz um primeiro calculo do objetivo na entrada
  fret,violret = Fun_Obj(p',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

  if toplot
     Atualiza_Display(p',numturb,fret,0,num2^2,p_grid,gridsize)
  end

  println(" Partindo da função objetivo ",fret)
  println(arquivo,0," ",fret," ",0," ",0.0," Ponto atual ",p')
  flush(arquivo)

   # Copia o ponto atual para a variável pt
  const pt = copy(p)

  # Vetores auxiliares
  const ptt = copy(p)
  const ppp = copy(p)
  
  # Loop principal
  for iter = 1:(10*num2)

    # Copia fret para fp
    const fp = fret

    # Incializa a variável que aponta para a posição de maior descrescimo do objetivo
    ibig = 0

    # Valor da maior variação
    del = 0.0

    # Loop pelas direções
    @showprogress 1 "Iteracao $iter || $fret " for i=1:num2

      # Extrai uma coluna da matriz de bases
      @inbounds for j=1:num2
          xit[j]=L[j,i]
      end

      # Copia do valor atual da função, já que vamos calcular um novo
      const fptt = fret

      # Line Search na direção de xit, partindo de p
      # A rotina limin também devolve o novo p
      fret, p = LS(p, xit, i, numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

      # Testa se esta direção é a de maior variação
      if (fptt-fret) > del
         del = fptt-fret
         ibig = i
      end #if

   end #i

   #  Critério de parada por variação do valor da função
   if abs(fp-fret) <= ftol*abs(fp)
      println(" Powell:: Tolerância atingida ",fp-fret)
      close(arquivo)
      return Array_Layout(p')
   end


   # Calcula os pontos para a frente (extrapolados) e translada o ponto atual
   @inbounds for j=1:num2
      ptt[j] = 2.0*p[j]-pt[j]
      xit[j] = p[j]-pt[j]
       pt[j] = p[j]
   end #j

   # Calcula o valor da função objetivo no ponto ptt
   fptt,violptt = Fun_Obj(ptt',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)


   # Critério do Powell
   if (fptt < fp)
      t = 2.0*(fp-2.0*fret+fptt)*sqrt(fp-fret-del)-del*sqrt(fp-fptt)
      if (t < 0.0)

          # Faz um line search na direção xit, partindo de p
          # Faz um line search na direção xit, partindo de p
          @inbounds for j=1:num2
              ppp[j] = p[j]
          end
          fret, p = LS(p, xit, ibig, numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

          # Copia o deslocamento de p para xit
          @simd for j=1:num2
             @inbounds xit[j] = p[j] - ppp[j]
          end

          # Move par o mínimo da nova direção e salva a nova direção
          @inbounds for j=1:num2
              L[j,ibig]=L[j,num2]
              L[j,num2]=xit[j]
          end #j
      end #if t
  end #if

  #println("\n Obj ",fret," na iteração externa ",iter)
  println(arquivo,iter," ",fret," ",ibig," ",del," Ponto atual ",p')
  flush(arquivo)

  if toplot
     Atualiza_Display(p',numturb,fret,iter,num2^2,p_grid,gridsize)
  end


end #iter

# Fecha o arquivo e retorna a solução
close(arquivo)
return Array_Layout(p')

end # Rotina..que só deve sair por tolerância ou por excesso de iterações


#
# Interface para fazer LS nos dois sentidos
#
function  LS(x, P, direcao, numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

        # Line search para frente (sentido = 1)
        flag, x_next,fN = Line_Search_Backtracing(x, P, 1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

        # Se o flag for -1, sentido contrário
        if flag==-1
           flag, x_next,fN = Line_Search_Backtracing(x, P, -1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)
        end

        if flag==-1
            println("\n LS::Powell::direção $direcao não mehorou em nenhum sentido.")
        end

        return fN, x_next
end



function Line_Search_Backtracing(x, d, sentido,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)


  # Sentido = 1 (frente) ou (-1) para tras...

  # Relaxação do alfa
    const tau = 0.5

  # Relaxação da inclinação inicial [0,0.25] (0.1)
    const cc = 0.1

  # Define um valor minimo de passo
    const minimo = 1E-12

  # Dimensao do vetor
    const nx = size(x,1)

  # d é um vetor
    d = vec(d)

  # Calcula o valor do custo no ponto atual
    f0,viol0 =  Fun_Obj(x',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

  # Normaliza a direção de busca, só para garantir...
    const nd = norm(d)
    @simd for j=1:nx
        @inbounds d[j] = d[j]/nd
    end

  # Verifica se não temos um alfa limite. Do contrário, utilizamos o máximo
  # FIXME -> adaptar ao tamanho do grid
    const alfa = 10.0

  # Define a posição futura
    const xf = Array(Float64,nx)

  # Flag de convergência
    const flag = 1

  # Produto interno D*d;..tem um -1 na frente porque a direção de minimização é D=-d...
    const direita1 = -1.0

  # Vou precisar deste cara fora do loop
    const fu = f0

  #println("\n ******************************************************************************************************* ")
  # Loop do Método
  @inbounds  for i=1:100

        # Posição futura
        @simd for j=1:nx
           @inbounds xf[j] =  x[j] + sentido*alfa*d[j]
        end


        # Valor na posição futura
        fu,violu = Fun_Obj(xf',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)


        #println(f0," ",fu," ",f0-fu," ", cc*alfa*direita1," ", cc*direita2," ",alfa)

        # Condição de Armijo
        if f0-fu < cc*(alfa*direita1)

            # Do contrário corta o passo
            alfa = alfa * tau

            # Se o passo for menor do que o mínimo, sai do algoritmo
            if alfa<minimo
                flag=-1
                break
            end
        else
            break
        end

    end #i
    #println("\n ======================================================================================================= ")
    return  flag,xf,fu

end # Armijo
