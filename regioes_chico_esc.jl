# Para regioes segmentadas do mapa, cada regiao tem um numero
# de turbinas na solucao de micrositting feita pelo Marco:
# Regiao r1 - 4t ; r2 - 5t ; r3 - 3t ; r4 - 4t ; r5 - 6t ;
#        r6 - 8t ; r7 - 8t ; r8 - 5t ;

r = []

r1 = [  (0.0,0.0),(140.0*50.0,0.0),(140.0*50.0,93.0*50.0),(0.0,93.0*50.0) ]
push!(r,r1)

centroides = []
regioes = []

for i = 1:size(r,1)

  push!(regioes,CH_Lenz(r[i]))
  push!(centroides,Centroide_CH(regioes[i]))

end

# Relaciona cada turbina a sua regiao
# nturb reg
reg_turb = [
1     1
2     1
3     1
4     1
5     1
6     1
7     1
8     1
9     1
10    1
11    1
12    1
13    1
14    1
15    1
16    1
17    1
18    1
19    1
20    1
21    1
22    1
23    1
24    1
25    1
26    1
27    1
28    1
29    1
30    1
31    1
32    1
33    1
34    1
35    1
36    1
37    1
38    1
39    1
40    1
41    1
42    1
43    1 ]
