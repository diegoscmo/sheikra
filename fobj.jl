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

  # Verifica uma area livre ao redor da turbina definida por uma elipse
  # com 4D de largura (2D entre turbinas) e 16D de altura (8D entre turbinas)
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
              end  #if norma
          end #if i!
      end # for j
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
      nr2 += 1.0
    end
  end

  # Aplica a violação das restrições no cálculo da função objetivo
  # onde um fator de penalização será utilizado para ajustar a dimensão
  # dos termos. Como a máxima violação que pode ocorrer em cada turbina
  # é de 220 e temos numturb turbinas, o caso mais crítico seria viol0
  # igual a 220*(numturb-1)^2. No entanto, o valor do objetivo é mumericamente bem maior
  # tal que um fator fixo de 100 pode ser usado nos testes iniciais.
  nr = nr1 + nr2

  return valor, nr
end #FunObj

#
# Testa se um ponto está dentro da elipse e retorna a distância normalizada
#
@everywhere function Elipse(xp,yp,xc,yc,altura,largura,angulo)

  c = cos(deg2rad(180.0+angulo))
  s = sin(deg2rad(180.0+angulo))

  xce = xp - xc
  yce = yp - yc

  xct = xce * c - yce * s
  yct = xce * s + yce * c

  dist_n = (xct^2/(largura/2.)^2) + (yct^2/(altura/2.)^2)

  return dist_n
end

#
# Devolve angulo em graus do setor predominante do vento
# a partir do grid de frequências do RSF
#
@everywhere function Ang_Pred(f_grid,numsec)

  pred = 0
  f_soma = 0.0

  # Compara a soma de frequencias de cada setor
  for i=1:numsec
    nov_soma = sum(f_grid[:,:,i])
    if  nov_soma > f_soma
      f_soma = nov_soma
      pred = i
  end
end

  # Devolve o setor como angulo e retorna
  ang = rad2deg(((pred-1)*360.0/numsec))
  return ang
end


#
#
# Verifica se o resultado está penalizado
#
function Verifica_Penalizacao(x,numturb,f_grid,A_grid,k_grid,
  z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

  valor,nr = Fun_Obj(x,numturb,f_grid,A_grid,k_grid,z_grid,p_grid,numsec,gridsize,pcurve,ctcurve)

  if nr > 0.0
    @printf("\nCUIDADO! Valor da função objetivo penalizado em %.3f.\n",nr)
  else
    @printf("\nValor obtido sem penalização.\n")
  end

end

#
# Gera a posicao de uma particula de modo que cada turbina
# fique dentro da sua região
#
@everywhere function Inicializa_Particula(numturb,anygrid,gridsize,regioes,centroides,reg_turb)

  # Inicializa array e valor maximo de posicao
  layout = zeros(Float64,numturb,2)
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
  return Layout_Array(layout)

end # Inicializa particula



#
# Atualiza a visualizacao dos resultados.
#
function Atualiza_Display(gbest,numturb,fgbest,t,nit,p_grid,gridsize,mostra_texto=false)

  const pt = Array(Float64,numturb,2)
  @inbounds for i=1:numturb
    pt[i,1] = gbest[i*2-1]
    pt[i,2] = gbest[i*2]
  end #

  clf()
  Plot_Grid(p_grid)
  Plot_WT(pt,gridsize)
  if mostra_texto
    @printf "Geração total: %.4f [MWh] - IT: " -fgbest
    @printf "%d/%d\n" t nit
  end

end

#
# Transforma o vetor posicao em um layout
#
@everywhere function Array_Layout(array)

  # Determina o numero de turbinas
  const numturb = convert(Int64,size(array,1)/2)

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
