function modulator(x,y,xmin,ymin,xsize,ysize)

# grampeia a funcao escolhida no retangulo definido

# modula o valor de A

return A

end

mapname = "map_test"

      xmin = 550906.0
      ymin = 6676455.0

      xsize = 5000.0
      ysize = 10000.0
      gridsize = 50.0
      numsec = 12

      xpts = convert(Int64,floor(xsize/gridsize))
      ypts = convert(Int64,floor(ysize/gridsize))

      text=Array(ASCIIString,(xpts+1)*(ypts+1),1)

      xpos = 0.0
      ypos = 0.0
      count = 1
      for i=1:ypts+1
        for j=1:xpts+1

          for s=1:numsec

          end

          posicao = ["            "
          string(convert(Int64,xmin+xpos))
          ".0 "   # fim do x
          string(convert(Int64,ymin+ypos))
          ".0     "    # fim do y
          "0.0"   # altura do relevo
          " 100.0 " #hub height
          string(9.19) #alpha de Weibull medio
          " "
          string("3.000") #k de Weibull medio
          "        "
          string(400.431) #potencia disponivel
          " "
          string(numsec) #numero de setores
          " "]



          #direcoes
          linha=""
          for i=1:15
            linha = linha*posicao[i]
          end
          for i=1:12

            if i==1
              freq = 999
            else
              freq = "  0"
            end

            direcoes = [
            string(freq) #frequencia da direcao i
            " "
            string(110) #weibull alpha da direcao
            "  "
            string(300) #weibull k da direcao
            ""]


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

      path = string("\\maps",mapname,".rsf")
      fh=open(path,"w")
      write(fh,text)
      close(fh)
      println("\nAll done, "mapname," created.")
