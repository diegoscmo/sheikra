
#
# Minha implementação do Jarvis March para determinar o Convex Hull
# Obviamnte, está bizonha e pode ser melhorada...mas  eu entendo :0)
#
# Testando se funciona ... diferença somente nesta linha...
@everywhere function CH_Lenz(points)


  # A ideia é a seguinte:
  # 1) Iniciamos com o ponto de menor coordenada y
  #    Este é o ponto atual e o ângulo inicial é zero.
  # 2) Varremos todos os pontos, calculando o ângulo
  #    entre o ponto atual e todos os demais. O menor dos ângulos positivos
  #    entre o ângulo de referência e os demais indica quem é o próximo ponto do CH.
  #    Este é o ponto atual e o novo ângulo.
  # 3) Fica fazendo isto até que o primeiro ponto seja o próximo ponto
  #
  # Complexidade: O(N) por iteração, totalizando h*O(N), onde h é o número de pontos
  # no convex Hull.


  # Numero de pontos na lista de coordenadas
  const N = length(points)

  # Acha o ponto com a menor coordenada y
  const menor_y       = Inf
  const ponto_inicial = 0
  @inbounds for i=1:N

      # Recupera a coordenada y
      const y = points[i][2]

      # Compara
      if y<menor_y
         menor_y = y
         ponto_inicial = i
      end

  end #i

  # Inicializa o convex Hull
  const ch = []
  push!(ch,points[ponto_inicial])

  #println("Ponto inicial ", ponto_inicial)

  # Inicia o algoritmo
  const ponto_atual  = ponto_inicial
  const angulo_atual = 0.0

  # Loop principal..pior caso
  @inbounds for ext=1:N

     # Coordenadas do ponto atual
     const xa = points[ponto_atual][1]
     const ya = points[ponto_atual][2]

     # Calcula o ângulo entre o ponto atual e os
     # demais pontos
     const menor_angulo_positivo = 2*pi
     const candidato = 0
     @inbounds for i = 1:N

         if i!=ponto_atual

            # Coordenadas do ponto em estudo
            const xc = points[i][1]
            const yc = points[i][2]

            # Monta o vetor posição
            const vpx =  xc-xa
            const vpy =  yc-ya

            # Normaliza
            norma = sqrt(vpx^2 + vpy^2)

            # Calcula o ângulo entre o vetor de referencia
            # dado por [cos(angulo_atual) sin(angulo_atual)]
            # e o vetor atual [vpx vpy]/norma
            # se ele for positivo, entra na comparação (sentido anti-horário)
            # do contrário, nem olhamos...
            const entrada  = cos(angulo_atual)*vpx/norma + sin(angulo_atual)*vpy/norma

            # Evita que a entrada do acos seja menor do que -1.0 e maior do que 1.0
            entrada = max(-1.0, min(entrada, 1.0))
            const angulo_teste = acos(entrada)

            if angulo_teste >= 0.0 && angulo_teste < menor_angulo_positivo
               candidato = i
               menor_angulo_positivo = angulo_teste
            end

         end # if i..

     end #i

     # Terminamos de varrer todos os pontos e temos um candidato.
     push!(ch,points[candidato])

     # Partimos deste ponto e deste ângulo atualizado
     ponto_atual   = candidato
     angulo_atual += menor_angulo_positivo

     # Agora, temos que verificar se ele não é o primeiro ponto e, se não for, seguir o baile...
     if (candidato==ponto_inicial)
         break
     end

  end # While...

  # Devolve o Convex Hull
  return ch

end


#
# Enconta o centróide de um CH::
#
# Estratégia simples. Partindo da origem,
# fazemos uma média ponderada pela distância
# de cada ponto. Esta posição será o centróide.
#
#
@everywhere function Centroide_CH(CH)

  # Coordenadas do centróide
  const cx = 0.0
  const cy = 0.0

  # Loop pelos pontos do Convex Hull
  const NP = 0
  for ponto in CH

    # Recupera as coordenadas
    const x = ponto[1]
    const y = ponto[2]

    # Acumula as coordenadas
    cx  += x
    cy  += y

    # Acumula o numero de pontos
    NP += 1

  end # ponto

  # Pondera pelo número de coordenadas no CH
  cx = cx / NP
  cy = cy / NP


  # Retorna o centróide
  return [cx cy]
end



#
# Detecta se um ponto está dentro de um triângulo:: MUITO IDIOTA....
#
# Ponto: par de coordenadas do ponto que está sendo testado [x y]
# v1: par de coordenadas de um vértice (ponto do CH) [x y]
# v2: par de coordenadas de um vértice (próximo ponto do CH) [x y]
# c3: par de coordenadas do centróide do CH [x y]
#
@everywhere function Dentro_Triangulo(ponto,v1,v2,ce)


  # Bem idiota....soma os ângulos (produtos internos) dos vetores
  # X1 = v1 - ponto -> vetor posição entre o ponto e o vértice 1
  # X2 = v2 - ponto -> vetor posição entre o ponto e o vértice 2
  # X3 = v3 - ponto -> vetor posição entre o ponto e o vértice 3
  #
  # Se o ponto estiver dentro do triângulo, então a soma
  # dos ângulos vai dar 2pi :0)
  const X1 = [(v1[1] - ponto[1]); (v1[2] - ponto[2])]
  const X2 = [(v2[1] - ponto[1]); (v2[2] - ponto[2])]
  const X3 = [(ce[1] - ponto[1]); (ce[2] - ponto[2])]

  # Se de uma cagada inominável e o ponto coincidir
  # com um dos vértices, então estamos dentro :0)
  const DX1 = norm(X1)
  const DX2 = norm(X2)
  const DX3 = norm(X3)

  DX1 == 0.0 && return true
  DX2 == 0.0 && return true
  DX3 == 0.0 && return true

  # Normaliza
  X1 = X1 / DX1
  X2 = X2 / DX2
  X3 = X3 / DX3

  # Soma os ângulos (supostamente) internos,
  # lembrando que dot(a,b) = ||a|| ||b|| cos (teta)
  # e, se ||a||=||b||=1.0, então a soma dos produtos
  # internos será cos(teta1) + cos(teta2) + cos(teta3)
  # tal que acos(cos(.......)) = 2pi
  #
  const argumento1 = max(-1.0,min(1.0,dot(X1,X2)))
  const argumento2 = max(-1.0,min(1.0,dot(X2,X3)))
  const argumento3 = max(-1.0,min(1.0,dot(X3,X1)))
  const a1 = acos(argumento1)
  const a2 = acos(argumento2)
  const a3 = acos(argumento3)

  const soma = a1 + a2 + a3

  # Se for aproximadamente igual a 2pi, devolve verdadeiro,
  # do contrário, false
  return isapprox(soma,2*pi,rtol=1E-4)

end

#
# Estratégia bobinha para detectar se um ponto
# está dentro de uma das regiões de projeto
#
# A lista de regiões contém matrizes em que cada linha é um
# ponto e cada coluna tem as coordenadas x e y do convex Hull (CH)
# da região, respectivamente.
#
@everywhere function To_Dentro(regioes,centroides,ponto)

   # Flag que indica se está dentro de uma região
   const dentro = false

   # Indicador de região
   const dentro_de_quem = 0

   # Loop por cada região
   @inbounds for r in 1:length(regioes)

       # Recupera a região
       const regiao = regioes[r]

       # Recupera o centróide desta região
       const cent = centroides[r]

       # Chama a rotina para uma região
       const dentro = To_Dentro_Regiao(regiao,cent,ponto)

       # Se dentro triangulo for verdadeiro, então podemos
       # sair do loop de regioes
       if dentro
          dentro_de_quem = r
          break
       end

   end #regiao

  # Devolve a região em que o ponto está dentro...
  return dentro_de_quem

end



#
# Estratégia bobinha para detectar se um ponto
# está dentro de uma das regiões de projeto
#
# A lista de regiões contém matrizes em que cada linha é um
# ponto e cada coluna tem as coordenadas x e y do convex Hull (CH)
# da região, respectivamente.
#
@everywhere function To_Dentro_Regiao(regiao,centroide,ponto)

   # Flag que indica se está dentro de uma região
   const dentro = false

  # A estratégia é a seguinte:
  # Se estiver dentro de algum triangulo entre o ponto e
  # os vértices do CH, então está dentro da região
  @inbounds for i=1:length(regiao)-1

    # ponto 1 e ponto 2
    const p1 = regiao[i]
    const p2 = regiao[i+1]

    # Testa este triangulo..se for true,
    # marca o flag e sai do primeiro loop
    if Dentro_Triangulo(ponto,p1,p2,centroide)
      dentro = true
      break
    end
  end

  # Devolve a região em que o ponto está dentro...
  return dentro

end






function TESTA_AQUI()

  # Define duas regiões ...aqui são os pontos do domínio....
  r1 = [(3.0,-1.0)  (1.0,1.0)  (2.0,3.0)  (3.0,2.0)  (5.0,3.0)  (5.0,1.0)  (6.0,-1.0)]
  r2 = [(6.0, 2.0)  (7.0,4.0)  (9.0,2.0)  (7.0, 0.0) ]

  # Calcula os Convex Hull dos conjuntos de pontos
  ch1 = CH_Lenz(r1)
  ch2 = CH_Lenz(r2)

  # Calcula o centroide da região
  c1 = Centroide_CH(ch1)
  c2 = Centroide_CH(ch2)

  # Monta as listas
  regioes    = []
  push!(regioes,ch1)
  push!(regioes,ch2)

  centroides = []
  push!(centroides,c1)
  push!(centroides,c2)

  println(regioes)
  println(centroides)

  # Agora podemos testar aonde está um dado ponto

  # Tem que dar 1
  rp1 = To_Dentro(regioes,centroides,[4.0 1.0])

  # Tem que dar 2
  rp2 = To_Dentro(regioes,centroides,[7.1 1.9])

  # Tem que dar 0
  rp3 = To_Dentro(regioes,centroides,[0.0 0.0])

  # Se algum der problema, mostramos um erro
  @assert rp1==1 "Erro ao testar o primeiro ponto"
  @assert rp2==2 "Erro ao testar o segundo ponto"
  @assert rp3==0 "Erro ao testar o terceiro ponto"

  # Todos os testes deram certo
  println("\n Tudo certo...")


  # Agora vamos testar o update na primeira região (só velocidade)
  #v, p = Atualiza_Posicao(ch1,c1,[4 2],[2 2],1.0,1.0,1.0,[4 2],[4 2])

end
