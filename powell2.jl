#
# Baseado em
# http://www.aip.de/groups/soe/local/numres/bookcpdf/c10-5.pdf, página 419
#
function Powell(numturb,nc,tol,rsf,toplot,p=[])

  # Tolerância para a parada (variação do objetivo)
  const ftol = tol

  # Loads chosen RSF
  p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,minx,miny = Load_RSF(string("maps/",rsf,".rsf"));

  # Inclui os dados das regioes para limitar espaço de busca de cada turbina
  # regioes, centroides e reg_turbinas
  @everywhere include("regioes_chico_esc.jl")

  # Loads Linear Interpolation of the power curve (here for standard density)
  pcurve  = open(readdlm,"others/power_curve.txt")[:,2]
  ctcurve = open(readdlm,"others/CT1.txt")[:,2]

  # Usamos bastante este cara aqui
  const num2 = 2*numturb

  # Se recebeu um ponto inicial, então não precisamos buscar um ponto
  # aleatório
  if length(p)==0


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
          obj = Fun_Obj_Powell(xp,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)

          # Se melhorou, então guarda
          if obj<melhor_valor[1]
                println(" Crazy_Joe::Improving... ",obj)
                melhor_valor[1] = obj
                melhor_posicao[:,1] = xp
          end


    end #p

    # Melhor chute é o nosso ponto incial no Powell
    p = vec(sdata(melhor_posicao)')

    # Libera a memória ...
    @everywhere gc()

    # Já que o resto não está paralelizado, vamos desconectar os processos
    np = nprocs()
    if np>1
       for i=np:-1:2
           rmprocs(i)
      end
    end

  else
      # Recebemos um ponto inicial...
      # Vamos verificar a consistência
      if size(p,1)!=num2
         error("\n Powell2::Ponto inicial não tem a dimensão correta. ",size(p))
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
  fret = Fun_Obj_Powell(p,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)

  if toplot
     Atualiza_Display(p,numturb,fret,0,num2^2,p_grid,gridsize)
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
    const ibig = 0

    # Valor da maior variação
    const del = 0.0

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
      fret, p = LS(p, xit, i, numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)

      # Testa se esta direção é a de maior variação
      if (fptt-fret) > del
         del = fptt-fret
         ibig = i
      end #if

   end #i

   #  Critério de parada por variação do valor da função
   #if 2.0*(fp-fret) <= ftol*(abs(fp)+abs(fret))+1E-10
   #if abs(fp-fret) <= ftol*abs(fret)
   #   println(" Powell:: Tolerância atingida ",(fp-fret))
   #   close(arquivo)
   #   return Array_Layout(p')
   #end


   # Calcula os pontos para a frente (extrapolados) e translada o ponto atual
   @inbounds for j=1:num2
      ptt[j] = 2.0*p[j]-pt[j]
      xit[j] = p[j]-pt[j]
       pt[j] = p[j]
   end #j

   # Calcula o valor da função objetivo no ponto ptt
   fptt = Fun_Obj_Powell(ptt,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)


   # Critério do Powell
   if (fptt < fp)
      t = 2.0*(fp-2.0*fret+fptt)*sqrt(fp-fret-del)-del*sqrt(fp-fptt)
      if (t < 0.0)

          # Faz um line search na direção xit, partindo de p
          @inbounds for j=1:num2
              ppp[j] = p[j]
          end
          fret, p = LS(p, xit, ibig, numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)

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

  # Como o metodo pode demorar, grava layout no disco
  layout = Array_Layout(p)
  writedlm(string("powellrun_lay.txt"),[layout[:,1]+minx layout[:,2]+miny])

  if toplot
     Atualiza_Display(p,numturb,fret,iter,num2^2,p_grid,gridsize)
  end


end #iter

# Fecha o arquivo
close(arquivo)

# Dá o display sobre a penalização do resultado
Verifica_Penalizacao(p,numturb,f_grid,A_grid,k_grid,
                             z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

# e retorna a solução
return Array_Layout(p)

end # Rotina..que só deve sair por tolerância ou por excesso de iterações


#
# Interface para fazer LS nos dois sentidos
#
function  LS(x, P, direcao, numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)

        # Line search para frente (sentido = 1)
        flag, x_next,fN,f0 = Line_Search_Backtracing(x, P, 1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)


        # Se o flag for -1, sentido contrário
        if flag==-1
           flag, x_next,fN,f0 = Line_Search_Backtracing(x, P, -1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)
        end

        if flag==-1
            println("\n LS::Powell::direção $direcao não mehorou em nenhum sentido. ")
            if f0<fN
               println("..original...")
               return f0,x
            else
               return fN,x_next
            end
        end

        return fN, x_next
end



function Line_Search_Backtracing(x, d, sentido,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)


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
    f0 =  Fun_Obj_Powell(x,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)

  # Normaliza a direção de busca, só para garantir...
    const nd = norm(d)
    @simd for j=1:nx
        @inbounds d[j] = d[j]/nd
    end

  # Verifica se não temos um alfa limite. Do contrário, utilizamos o máximo
  # FIXME -> adaptar ao tamanho do grid
    const alfa = 100.0

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
    @inbounds for i=1:100

        # Posição futura
        @simd for j=1:nx
           @inbounds xf[j] =  x[j] + sentido*alfa*d[j]
        end

        # Valor na posição futura
        fu = Fun_Obj_Powell(xf,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,regioes,centroides,reg_turb)


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
    return  flag,xf,fu,f0

end # Armijo

#
# A partir de uma posição retorna o valor da função e da soma das restrições
# fique dentro da sua região
#
@everywhere function Fun_Obj_Powell(x,numturb,f_grid,A_grid,k_grid,
                                    z_grid,p_grid,numsec,gridsize,pcurve,ctcurve,
                                    regioes,centroides,reg_turb)

    # Converte o array posicao para um layout
    const layout = Array_Layout(x)

    # Calcular o valor da funcao objetivo (-AEP total)
    # (Recupera somente a ultima posição da matriz)
    valor = -1.0 * Calc_Park(layout,f_grid,A_grid,k_grid,z_grid,
                             p_grid,numsec,gridsize,pcurve,ctcurve)[end,4]

    # Calcula o valor das restricoes

    # Verifica uma area livre ao redor da turbina definida por uma elipse
    # com 4D de largura e 16D de altura
    const nr1 = 0.0
    const altura  = 110.0*15.70
    const largura = 110.0*4.0
    for i=1:numturb
        for j=1:numturb
            if i!=j
                xi = layout[i,1]
                yi = layout[i,2]
                xj = layout[j,1]
                yj = layout[j,2]

                norma = Elipse(xj,yj,xi,yi,altura,largura,angulo)

                if norma < 1.0
                    nr1 = nr1 + (1.0 - norma)
                end #if norma
            end # if i!
        end #for j
    end # for i


    # Também devemos calcular, quando estamos utilizando o Powell,
    # Se cada turbina está no seu devido setor. Do contrário, temos que
    # Aplicar uma penalidade...para isto vamos precisar de novas informações,
    # tais como reg_turb, regioes e centróides.
    const nr2 = 0.0
    @inbounds for i=1:numturb

           # Busca regiao onde deveria estar a turbina
           const r = reg_turb[i,2]

           # Se
           if !To_Dentro_Regiao(regioes[r],centroides[r],layout[i,:])
               nr2 +=1
           end
    end



  # Aplica a violação das restrições no cálculo da função objetivo
  # onde um fator de penalização será utilizado para ajustar a dimensão
  # dos termos. Como a máxima violação que pode ocorrer em cada turbina
  # é de 200 e temos numturb turbinas, o caso mais crítico seria viol0
  # igual a 200*(numturb-1)^2. No entanto, o valor do objetivo é mumericamente bem maior
  # tal que um fator fixo de 100 pode ser usado nos testes iniciais.
    valor = valor + PENAL*nr1 + PENAL2*nr2

    return valor
end #FunObjPowell
