#-----------------------------------------------------#
# Diretivas de compilação para o programa SIMS        #
#                                                     #
#             Escrito inicialmente por                #
#               Rafael Gabler Gontijo                 #
#                        em                           #
#                    01/12/2012                       #
#                                                     #
# Objetivo: Compilar uma determinada versão do código #
#                                                     #
# Versão atual: - Aceita compilação para GNU/Linux ou #
#               - Compilador: GNU gfortran.           #
#                                                     #
#                                                     #
# Última modificação por: Rafael Gabler Gontijo       #
# Última modificação em: 01/12/2012                   #
#-----------------------------------------------------#
#
# Arquivos fonte padrão (linux) e outros sistemas
#
SRC-LNX = funcoes.f90 variaveis.f90 entrada.f90  \
          principal.f90 saida.f90 sims.f90  \
          
 OBJ-LNX = $(SRC-LNX:.f90=.o)

#
# Definição do compilador específico e das flags

#ENGINE-LNX = gfortran
ENGINE-LNX = ifort
FLAGS-LNX = -m64 -O2
#FLAGS-LNX = -m64 -O2 -ldislin 
#FLAGS-LNX = -m64 -O2 -openmp -traceback -fpe0 -g

#
# Objetivo principal: geração de executável em linux


sims.ex : $(SRC-LNX) 
	$(ENGINE-LNX) $(FLAGS-LNX) -o sims.ex $(SRC-LNX)


#
# Limpeza
clean :
	rm sims.ex
