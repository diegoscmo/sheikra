#
# Método de Powell
# http://www.personal.psu.edu/cxg286/Math555.pdf, página 99/100
# Artigo original: http://comjnl.oxfordjournals.org/content/7/2/155.full.pdf+html
#
function Powell(numturb,nc,tol,rsf,toplot)

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
  const x_now = Crazy_Joe(numturb,nc,rsf,toplot)
    
  # Agora vamos lá
  const xdiff = 1.0
  const count = 0
  const best_index = 0

  # Tolerância e numero máximo de iterações
  const Nmaxiter= (2*numturb)^2

  # Gera a base inicial de busca
  const L = eye(2*numturb,2*numturb)

  # Vamos abrir um arquivo para escrita, pois o método é lento
  # e vale a pena acompanhar o andamento
  arquivo = open("convergencia_powell.txt","w")

  # Faz um primeiro calculo do objetivo na entrada
  fN,violN = Fun_Obj(x_now',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

  println(arquivo,count," ",fN," ",0," ",0.0," Ponto atual ",x_now')
  flush(arquivo)

  # Loop principal
  while xdiff>tol && count < Nmaxiter

    # Incrementa o contador
    count +=1

    # Variáveis internas do while
    const Df = 0.0
    const x0 = copy(x_now)
    const r = Array(Float64,2*numturb)
    const P = Array(Float64,2*numturb)

    # Preciso deste cara depois do loop de direções
    const x_next = zeros(x_now)

    #println("\n Avaliando as direções para esta iteração")

    # Loop das direções
    @showprogress 10 "Iteracao $count / $Nmaxiter..." for i=1:2*numturb

        # Vetor direção atual
        @inbounds for j=1:2*numturb
            P[j] = L[j,i]
        end

        # Guarda o valor de fN para podermos calcular a diferença
        const valor_anterior = fN

        # Line search para frente (sentido = 1)
        flag, x_next,fN = Line_Search_Backtracing(x_now, P, 1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

        # Se o flag for -1, sentido contrário
        if flag==-1
           flag, x_next,fN = Line_Search_Backtracing(x_now, P, -1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)
        end

        if flag==-1
            println("\n Powell::direção $i não mehorou em nenhum sentido do LS")
        end

        # Vetor de diferença e sua norma
        @simd for j=1:2*numturb
           @inbounds r[j] = x_next[j] - x_now[j]
        end
        const Nr = norm(r)

        # Seleciona a direção com a mehor descida
        if (valor_anterior - fN)>Df
          Df = (valor_anterior - fN)
          #println(Df," ",i," ",Nr)
          best_index = i
        end

        # Rola o ponto atual
        @inbounds for j=1:2*numturb
           x_now[j] = x_next[j]
        end

    end #i

    # Ao terminarmos de varrer todas as direções, temos que nos preocupar com a base
    # pois os vetores podem ficar L.I
    # Vamos usar a Heurística do Powell...
    # Calcula a função em diferentes pontos
    f0,viol0 = Fun_Obj(x0',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)
    fN,violN = Fun_Obj(x_next',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)
    fE,violE = Fun_Obj((2.0*x_next-x0)',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

    # Segundo Powell, esta estratégia é melhor do que a tradicional se
    # a função não for quadrática.
    if ! ( (fE>=f0) || (2*(f0-2*fN+fE)*(f0-fN-Df)^2 >= Df*(f0-fE)^2) )

         println("\n Powel::Resetando a direção $best_index")

         # Arruma a base para não perder a independência linear.
         @simd for j=1:2*numturb
             @inbounds P[j] = x_next[j] - x0[j]
         end

         # Grava a nova direção
         # http://www.aip.de/groups/soe/local/numres/bookcpdf/c10-5.pdf, página 419
         @inbounds for j=1:2*numturb
             L[j,best_index] = L[j,end]
             L[j,end]        = P[j]
          end

        # Calcula o passo na direção deste novo P, e assume como ponto inicial
        # para a próxima iteração.
        flag, x_now,fN = Line_Search_Backtracing(x0, P, 1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

        # Se o flag for -1, sentido contrário
        if flag==-1
           flag, x_now,fN = Line_Search_Backtracing(x0, P, -1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)
        end


    end #if !

    println("\n Obj ",fN," na iteração externa ",count)

    println(arquivo,count," ",fN," ",best_index," ",Df," Ponto atual ",x_now')
    flush(arquivo)

    # Criterio de convergência
    xdiff = norm(P)

    if toplot
      Atualiza_Display(x_now',numturb,fN,count,Nmaxiter,p_grid,gridsize)
    end


  end #While

  close(arquivo)
  return  Array_Layout(x_now)


end # Powell



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
