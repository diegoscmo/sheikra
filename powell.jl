#
# Método de Powell
# http://www.personal.psu.edu/cxg286/Math555.pdf, página 99/100
#
function Powell(numturb,tol,rsf,toplot)

  # Loads chosen RSF
  p_grid,z_grid,f_grid,A_grid,k_grid,numsec,gridsize,minx,miny = Load_RSF(string("maps/",rsf,".rsf"));

  # Inclui os dados das regioes para limitar espaço de busca de cada turbina
  # regioes, centroides e reg_turbinas
  include("regioes_chico.jl")

  # Loads Linear Interpolation of the power curve (here for standard density)
  pcurve  = open(readdlm,"others/power_curve.txt")[:,2]
  ctcurve = open(readdlm,"others/CT1.txt")[:,2]

  # Gera um ponto inicial para o método
  x_now = Inicializa_Particula(numturb,A_grid,gridsize,regioes,centroides,reg_turb)

  #println("\n Saindo do ponto ", x_now')

  # Agora vamos lá
  const xdiff = 1.0
  const count = 0
  const best_index = 0

  # Tolerância e numero máximo de iterações
  const Nmaxiter= (2*numturb)^2

  # Gera a base inicial de busca
  const L = speye(2*numturb,2*numturb)

  # Loop principal
  while xdiff>tol && count < Nmaxiter


    # Incrementa o contador
    count +=1

    # Variáveis internas do while
    const Df = 0.0
    const x0 = copy(x_now)

    # Preciso deste cara depois do loop de direções
    const x_next = zeros(x_now)
    const fN = 0.0

    println("\n Avaliando as direções para esta iteração")

    # Loop das direções
    for i=1:2*numturb

        # Vetor direção atual
        P = full(L[:,i])

        # Line search para frente (sentido = 1)
        flag, x_next = Line_Search_Backtracing(x_now, P, 1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

        # Se o flag for -1, sentido contrário
        if flag==-1
           flag, x_next = Line_Search_Backtracing(x_now, P, -1.0,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)
        end

        # Vetor de diferença e sua norma
        const r = x_next - x_now
        const Nr = norm(r)

        #println("\n residuo $Nr na dir $i")

        # Melhor Indice até agora
        if Df<Nr
          Df = Nr
          best_index = i
        end

        # Rola o ponto atual
        x_now = copy(x_next)

    end #i

    # Ao terminarmos de varrer todas as direções, temos que nos preocupar com a base
    # pois os vetores podem ficar L.I
    # Vamos usar a Heurística do Powell...

    P = x_next - x0
    f0,viol0 = Fun_Obj(x0',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)
    fN,violN = Fun_Obj(x_next',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)
    fE,violE = Fun_Obj((2*x_next-x0)',numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

    if !(fE<=f0) || (2*(2*fN+fE-f0)*((fN-f0)-Df)^2 >= Df*(fE-f0)^2)

         # Arruma a base para não perder a IL.
         L[:,best_index] = P
    end

    println("\n Obj ",fN," ",count)

    # Criterio de convergência
    xdiff = norm(P)

    if toplot
      Atualiza_Display(x_now',numturb,fN,count,Nmaxiter,p_grid,gridsize)
    end


  end #While

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
    return  flag,xf

end # Armijo
