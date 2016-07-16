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
  include("regioes_chico.jl")

  # Loads Linear Interpolation of the power curve (here for standard density)
  pcurve  = open(readdlm,"others/power_curve.txt")[:,2]
  ctcurve = open(readdlm,"others/CT1.txt")[:,2]

  # Faz um início aleatório com nc tentativas. SE nc==1, então equivale ao
  # uso do Inicializa_Particulas.
  const p = Crazy_Joe(numturb,nc,rsf,toplot)

  # Usamos bastante este cara aqui
  const num2 = 2*numturb

  # Uma direação de busca
  const xit = Array(Float64,num2)

  # Libera a memória ...
  @everywhere gc()


  # Gera a base inicial de busca
  const L = eye(num2,num2)


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
  const ptt = copy(p)

  # Loop principal
  for iter = 1:(num2)^2

    # Copia fret para fp
    const fp = fret

    # Incializa a variável que aponta para a posição de maior descrescimo do objetivo
    ibig = 0

    # Valor da maior variação
    del = 0.0

    # Loop pelas direções
    @showprogress 1 "Iteracao $iter || $fret " for i=1:num2

      # Extrai uma coluna da matriz de bases
      for j=1:num2
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
   for j=1:num2
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
          ppp = copy(p)
          fret, p = LS(p, xit, ibig, numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

          # Copia o deslocamento de p para xit
          xit = p - ppp

          # Move par o mínimo da nova direção e salva a nova direção
          for j=1:num2
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

  # Calcula o valor do custo no ponto atual
    f0,viol0 =  Fun_Obj(x',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

  # Normaliza a direção de busca, só para garantir...
    d = vec(d/norm(d))

  # Verifica se não temos um alfa limite. Do contrário, utilizamos o máximo
  # FIXME -> adaptar ao tamanho do grid
    const alfa = 10.0

  # Copia para a saida
    const xf = copy(x)

  # Flag de convergência
    const flag = 1

  # Produto interno D*d;..tem um -1 na frente porque a direção de minimização é D=-d...
    const direita1 = -1.0

  # Vou precisar deste cara fora do loop
    const fu = f0

  #println("\n ******************************************************************************************************* ")
  # Loop do Método
    for i=1:100

        # Posição futura
        xf =  x + sentido*alfa*d

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
