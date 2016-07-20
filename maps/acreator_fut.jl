#cd(string(ENV["SHEIKRA_DIR"],"\\map_creator"))

function Modulator(x,y,xmin,ymin,xsize,ysize)

# grampeia a funcao escolhida no retangulo definido

# modula o valor de A

return A

end
      # "Map name.rsf"
      mapname = "map_test.rsf"

      # Geographic location, use at least 100000 and 6000000
      xmin = 500000.0
      ymin = 6000000.0

      xsize = 1000.0
      ysize = 1000.0
      gridsize = 50.0
      numsec = 12

      xpts = convert(Int64,floor(xsize/gridsize))
      ypts = convert(Int64,floor(ysize/gridsize))

      text=Array(ASCIIString,(xpts+1)*(ypts+1),1)

      xpos = 0.0
      ypos = 0.0
      count = 1
      for i=1:xpts+1
        for j=1:ypts+1

          for s=1:numsec

          end

          xstring = @sprintf("%.1f ",xmin+xpos)
          ystring = @sprintf("%.1f     ",ymin+ypos)

          Amedstring = @sprintf("%.1f ",Amed) #11.0
          kmedstring = @sprintf("%.3f        ",kmed) #3.000
          pstring    = @sprintf("%.3f ",p) #400.431
          numsecstring = @sprintf("%.0f ",numsec)

          posicao = ["            "
          xstring
          ystring
          "0.0 "   # altura do relevo
          "120.0 " #hub height
          Amedstring #alpha de Weibull medio
          kmedstring #k de Weibull medio
          pstring #potencia disponivel
          numsecstring #numero de setores
          ]

          valor = 3.20
          fstring = @sprintf("%.0f",f)
          Astring = @sprintf("%.0f",A)
          kstring = @sprintf("%.0f",k)

          direcoes = [
          fstring #frequencia da direcao i #990
          " "
          Astring #weibull alpha da direcao #110
          "  "
          kstring #weibull k da direcao #300
          " "]

          #direcoes
          linha=""
          for i=1:15
            linha = linha*posicao[i]
          end
          for i=1:12
            for j=1:6
            linha= linha*direcoes[j]
            end
          end

          text[count,1] = linha*"\r\n"

          xpos = xpos+gridsize
          count += 1
        end
        xpos = 0.0
        ypos = ypos+gridsize
      end

      path = string(mapname)
      fhand=open(path,"w")
      write(fhand,text)
      close(fhand)
      println("\nAll done, "mapname," created.")
