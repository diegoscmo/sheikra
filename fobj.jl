#
# A partir de uma posição retorna o valor da função e da soma das restrições
# fique dentro da sua região
#
@everywhere function Fun_Obj(x,numturb,f_grid,A_grid,k_grid,
                             z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

    # Converte o array posicao para um layout
    layout = Array_Layout(x)

    # Calcular o valor da funcao objetivo (-AEP total)
    # (Recupera somente a ultima posição da matriz)
    valor = -1.0 * Calc_Park(layout,f_grid,A_grid,k_grid,z_grid,
                             p_grid,numsec,gridsize,pcurve,ctcurve)[end,4]

    # Calcula o valor das restricoes
    # Aqui é verificado o espacamento minimo de 2*D entre as turbinas
    const nr = 0.0
    @inbounds for i=1:numturb
        @inbounds for j=1:numturb
            if i!=j
                dis = norm(layout[i,:]-layout[j,:])
                if dis < 2.0*(110.0)
                    nr = nr + (220-dis) #(1.0/dis)
                end
            end
        end
    end

  # Aplica a violação das restrições no cálculo da função objetivo
  # onde um fator de penalização será utilizado para ajustar a dimensão
  # dos termos. Como a máxima violação que pode ocorrer em cada turbina
  # é de 200 e temos numturb turbinas, o caso mais crítico seria viol0
  # igual a 200*(numturb-1)^2. No entanto, o valor do objetivo é mumericamente bem maior
  # tal que um fator fixo de 100 pode ser usado nos testes iniciais.
    valor = valor + PENAL*nr

    return valor, nr
end #FunObj


#
# Gera a posicao de uma particula de modo que cada turbina
# fique dentro da sua região
#
@everywhere function Inicializa_Particula(numturb,anygrid,gridsize,regioes,centroides,reg_turb)

  # Inicializa array e valor maximo de posicao
  layout = Array(Float64,numturb,2)
  xmax   = size(anygrid,1)*gridsize
  ymax   = size(anygrid,2)*gridsize

  @inbounds for i=1:numturb
    # Busca regiao onde deveria estar a turbina
    const r = reg_turb[i,2]

    # Insere a turbina aleatoriamente ate que ela esteja na sua devida regiao
    while To_Dentro_Regiao(regioes[r],centroides[r],layout[i,:]) == false

      layout[i,1] = xmax*rand()
      layout[i,2] = ymax*rand()

    end #while

  end # for i

  # Converte a saida para array 1D
  array = Layout_Array(layout)

  return array
end # Inicializa particula



#
# Atualiza a visualizacao dos resultados.
#
function Atualiza_Display(gbest,numturb,fgbest,t,nit,p_grid,gridsize)

    const pt = Array(Float64,numturb,2)
    for i=1:numturb
       pt[i,1] = gbest[i*2-1]
       pt[i,2] = gbest[i*2]
    end #

    clf()
    Plot_Grid(p_grid)
    Plot_WT(pt,gridsize)
    @printf "Geração total: %.4f [MWh] - IT: " -fgbest
    @printf "%d/%d\n" t nit

end

#
# Transforma o vetor posicao em um layout
#
@everywhere function Array_Layout(array)

  # Determina o numero de turbinas
  const numturb = convert(Int64,size(array,2)/2)

  # Inicializa o vetor Layout
  const layout = Array(Float64,numturb,2)

  # Converte array para layout x,y
  @inbounds for i=1:numturb
    layout[i,1] = array[i*2-1]
    layout[i,2] = array[i*2]
  end

  return layout

end # Array_Layout

#
# Transforma o layout em um array 1D
#
@everywhere function Layout_Array(layout)

  # Determina o numero de turbinas
  const numturb = size(layout,1)

  # Inicializa o vetor Layout
  const array = Array(Float64,numturb*2)

  # Converte layout x,y para array 1D
  @inbounds for i=1:numturb
    array[i*2-1] = layout[i,1]
    array[i*2]   = layout[i,2]
  end

  return array

end # Layout_Array

