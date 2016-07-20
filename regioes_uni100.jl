# Para regioes segmentadas do mapa, cada regiao tem um numero
# de turbinas na solucao de micrositting feita pelo Marco:
# Regiao r1 - 4t ; r2 - 5t ; r3 - 3t ; r4 - 4t ; r5 - 6t ;
#        r6 - 8t ; r7 - 8t ; r8 - 5t ;

r = []

r1 = [  (0.0,0.0),(101.0*50.0,0.0),(101.0*50.0,201.0*50.0),(0.0,201.0*50.0)]
push!(r,r1)

centroides = []
regioes = []

for i = 1:size(r,1)

  push!(regioes,CH_Lenz(r[i]))
  push!(centroides,Centroide_CH(regioes[i]))

end

# Relaciona cada turbina a sua regiao
# nturb reg
reg_turb = zeros(Int64,100,2)
reg_turb[:,1] = [1:100]
reg_turb[:,2] = 1
