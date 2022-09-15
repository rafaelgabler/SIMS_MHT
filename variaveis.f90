module variaveis

!******************************** Definindo variaveis **************************************!
integer i,j,k ,q,s, a,b,c,d,g
integer N, auxper(5), auxper2(5), auxper3(5)
integer razao, iter
real razao2, diferenca_interpol1, diferenca_interpol2
real diferenca, ajuda, lambda, kreal, aux_real, bifmax
real derivada1, derivada2, derivada3, trapezio, intervalo
real, allocatable ::gpvetor(:)
real, allocatable ::campo(:)
real, allocatable :: y(:)
real, allocatable :: trap(:)
real, allocatable :: cof1(:,:)
real, allocatable :: cof2(:,:)
real, allocatable :: cof3(:,:)
real, allocatable :: cof4(:,:)
real, allocatable :: cof5(:,:)
real, allocatable :: cof6(:,:)
real, allocatable :: cof7(:,:)
real, allocatable :: cof8(:,:)
real, allocatable :: Di(:,:,:)
real, allocatable :: centro_massa(:,:)
real, allocatable :: aux1(:,:),aux2(:,:),aux3(:,:), aux4(:,:)
real, allocatable :: auxt(:,:),auxf(:,:)
real, allocatable :: torquereal(:,:),torquereciproco(:,:)
real, allocatable :: forcareal(:,:),forcareciproca(:,:)
real, allocatable :: contribuicao_self(:,:)
real, allocatable :: contribuicao_fisico(:,:)
real, allocatable :: contribuicao_reciproco(:,:)
real, allocatable :: contribuicoes(:,:)
real, allocatable :: FORCAS(:,:,:,:)
real, allocatable :: TORQUES(:,:,:,:)
real, allocatable :: FT(:,:,:)
real, allocatable :: Tt(:,:,:)
real, allocatable :: X(:,:,:)
real, allocatable :: W(:,:,:)
real, allocatable :: XI(:,:,:,:)
real, allocatable :: U(:,:,:)
real, allocatable :: V(:,:,:)
real, allocatable :: var(:,:)
real, allocatable :: errovar(:,:)
real, allocatable :: flutmag(:,:)
real, allocatable :: auxiliarcor(:,:)
real, allocatable :: nr(:)
real, allocatable :: tempototal(:)
real, allocatable :: aux_erro_vel(:,:,:)
real, allocatable :: aux_erro_var(:,:,:)
real, allocatable :: variancia(:,:,:)
real, allocatable :: flut(:,:,:,:)
real, allocatable :: auxcor(:,:,:)
real, allocatable :: errocor(:,:)
real, allocatable :: funcaor(:,:)
real, allocatable :: dif(:,:)
real, allocatable :: velmedia(:,:,:)
real, allocatable :: vmedia(:,:)
real, allocatable :: errovmedia(:,:)
real, allocatable :: fmedia(:,:)
real, allocatable :: errofmedia(:,:)
real, allocatable :: flut_tempo(:,:)
real, allocatable :: usistema(:,:)
real, allocatable :: ILF(:,:)
real, allocatable :: ILR(:,:)
real, allocatable :: DIAM(:,:)
real, allocatable :: auxflut(:,:,:)
real, allocatable :: dt(:,:)
real, allocatable :: diarand(:)
real, allocatable :: hidrodinamica_aux1(:,:)
real, allocatable :: hidrodinamica_aux2(:,:)
real, allocatable :: hidro1(:,:)
real, allocatable :: hidro2(:,:)
real, allocatable :: magtempo(:,:)
real, allocatable :: beta(:,:)
real, allocatable :: potencial(:,:)
real, allocatable :: energia(:,:)
real, allocatable :: energiaantes(:,:)
real, allocatable :: dedt(:,:)
real, allocatable :: heatcurrent(:,:,:)
real, allocatable :: heataux(:,:,:)
real, allocatable :: auxj(:,:)
real, allocatable :: intj(:,:)
real, allocatable :: jmedio(:,:)
real coeficiente3, kr2, coeficiente7, coeficiente8
real mobilidade(3,3)
real resistencia(3,3), teste_identidade(3,3), coeficiente4, coeficiente5, coeficiente6
real mobilidade_self(3,3), mobilidade1(3,3), mobilidade2(3,3), coeficiente1, coeficiente2
real resistencia_self(3,3), resistencia1(3,3), resistencia2(3,3), pi
integer aux_periodico, posicao_campo
real dt_inicial, percentual
real contri1,contri2,contri3,contri4,contri5,contri6
integer nb,nbr
real rij(3), rn(3)
real l, h
integer npast
integer ncor
real phi
real nr1, nr2, nr3, Pe, Per, ragreg
integer rea
integer semente(1)
integer auxint1,auxint2,auxint3
integer auxint4
integer nnr
real r, alpha2,alpha
real xcentro, ycentro, zcentro
real xmin, xmax, ymin, ymax, zmin, zmax
real modrij
real modrand
real termo1,termo2
real termo3,termo4,termo5,eps
real auxcont, freq, shearratei
real tempo
real ti, tf, tpros
real reale, reale2, Str, freqcampo, freqbeat
integer inteiro, inteiro2, loop, loop2
integer auxiliar1
logical beating
logical posicao, fator
logical velocidade
logical estatistica, dipolo_ordenado
logical estatica
logical ordenado
logical torque
logical externo
logical gravidade
logical leito
logical mistura
logical gravadipolo
logical grafmag
logical browniano
logical tmagper
logical fmagper
logical ligaih
logical periodicidade
logical shear
logical oscillatory
logical agregado_inicial
logical continua
logical polidispersidade
logical printphi
logical oscilacampo
logical inertia
logical rotating
logical greenkubo
logical duffing
logical bifurcation, bifshear
 character(3) rea_char
 character(23) linha1
 character(20) linha2
 character(25) linha3
integer n2, auxiliar_continua
real n3, modip, shearrate, St
real teste1, frequencia
integer teste2, nfreq, contfreqinteiro1, contfreqinteiro2, multiplofreq
real h2, freqmax
real UMEDIA(3)
real SIGMA(3)
real, allocatable :: difusaoaux(:,:)
real difusao
real umsobreene
integer DeAllocateStatus
real dist
real qsi, konda(3), knormal(3), modk, kr
real C1, C2, C3, C4
real k1,k2,k3,k4,g1,g2,g3,g4

!****************************** Legenda das variaveis utilizadas ***************************!

! Variaveis inteiras

! i,j,k,q,s = variaveis inteiras utilizadas em loops
! N = numero de particulas
! nr(:) = vetor que armazena numeros randomicos necessarios para expressar movimento browniano,
! posicao inical e momento de dipolo inicial
! npast = numero de iteracoes e numero de passos de tempo
! rea = numero de realizacoes
! semente(1) = semente randomica do programa para auxilio na subrotina geradora de numeros aleatorios
! auxint1, auxint2, auxint3 e auxint4 =  Inteiros auxiliares para garantir que cada numero aleatorio
! seja gerado a partir de uma semente randomica
! nnr = quantidade de numeros randomicos gerados
! fundo = conta quantas particulas ja chegaram ao fundo do box


! Variaveis reais

! X(:,:) = Matriz posicao de todas as particulas
! U(:,:) = Matriz velocidade de todas as particulas
! W(:,:) = Matriz velocidade angular de todas as particulas
! flut1(:,:,:) = Matriz que contabiliza a flutuacao de velocidade de cada particula em cada time step
! e em cada realizacao (matriz tri-dimensional)
! flut2(:,:,:) = Matriz que contabiliza a flutuacao de velocidade de cada particula em cada time step
! e em cada realizacao (matriz tri-dimensional)
! flut3(:,:,:) = Matriz que contabiliza a flutuacao de velocidade de cada particula em cada time step
! e em cada realizacao (matriz tri-dimensional)
! flutrea1(:,:) = Matriz que conta a flutuacao media da suspensao em cada realizacao (direcao 1)
! flutrea2(:,:) = Matriz que conta a flutuacao media da suspensao em cada realizacao (direcao 1)
! flutrea3(:,:) = Matriz que conta a flutuacao media da suspensao em cada realizacao (direcao 1)
! flutfinal1(:) = Matriz que conta a flutuacao final da suspensao (media de realizacoes) (direcao 1)
! flutfinal2(:) = Matriz que conta a flutuacao final da suspensao (media de realizacoes) (direcao 2)
! flutfinal3(:) = Matriz que conta a flutuacao final da suspensao (media de realizacoes) (direcao 3)
! flutmedia1(:) = velocidade media da suspensao naquele instante de tempo (direcao 1)
! flutmedia2(:) = velocidade media da suspensao naquele instante de tempo (direcao 2)
! flutmedia3(:) = velocidade media da suspensao naquele instante de tempo (direcao 3)
! autocor1(:) = funcao autocorrelacao das flutuacoes de velocidade na direcao 1
! autocor2(:) = funcao autocorrelacao das flutuacoes de velocidade na direcao 2
! autocor3(:) = funcao autocorrelacao das flutuacoes de velocidade na direcao 3
! veli1(:,:) = velocidade de cada particula em cada instante de tempo na direcao 1
! veli2(:,:) = velocidade de cada particula em cada instante de tempo na direcao 2
! veli3(:,:) = velocidade de cada particula em cada instante de tempo na direcao 3
! velmedia1(:) = velocidade media de cada particula naquela realizacao (direcao 1)
! velmedia2(:) = velocidade media de cada particula naquela realizacao (direcao 2)
! velmedia3(:) = velocidade media de cada particula naquela realizacao (direcao 3)
! dif(:,:) = coeficiente de difusao nas 3 direcoes
! magnetizacao(:) = magnetizacao ao final de cada realizacao da suspensao
! Fp(:,:) = Matriz forcas relacionadas a influencia das paredes fisicas do box nas particulas
! Fr(:,:) = Matriz forcas de repulsao entre particulas do box
! Fc(:,:) = Matriz forcas de contato entre particulas do box
! Fcp(:,:) = Matriz forcas de contato entre particulas do box e as paredes fisicas do recepiente
! Fb(:,:) = Matriz forcas brownianas de todas as particulas
! Ft(:,:) = Matriz forca total de todas as particulas
! Fm(:,:) = Matriz forcas mageticas
! Fce(:,:) = Matriz forcas mageticas devido a um campo externo aplicado
! Tb(:,:) = Matriz torques brownianos
! Tm(:,:) = Matriz torques magneticos
! Tme(:,:) = Matriz torques magneticos externos
! Tt(:,:) = Matriz torques totais
! Di(:,:) = Matriz momento de dipolo magnetico de todas as particulas
! aux1(:) = Vetor auxiliar 1 no processo de implementacao de forcas magneticas
! aux2(:) = Vetor auxiliar 2 no processo de implementacao de forcas magneticas
! aux3(:) = Vetor auxiliar 3 no processo de implementacao de forcas magneticas
! T(:) = Vetor temporal
! rij(3) = Vetor R_{ij} unitario que aponta da particula i na direcao da particula j
! agregtotal(:) = Matriz que informa o numero de agregados de x particulas das realizacoes acumuladas

! dt = passo de tempo
! l e h = dimensoes do box (adimensionais)


! St = Numero de Stokes
! Pe = Numero de Peclet
! Arm = Numero de Stokes magnetico
! Str = Numero de Stokes rotacional
! Per = Numero de Peclet rotacional
! Armr = Numero de Stokes magnetico rotacional
! alpha = Energia browniana sobre energia magnetica
! Pc = Parametro de contato, numero adimensional utilizado na implementacao de forcas de contato
! phi = Fracao volumetrica de particulas
! magnetizacao_media = magnetizacao media em cima de varias realizacoes

! nr1,nr2,nr3 = numeros randomicos gerados para implementar forca browniana
! r = Distancia entre duas particulas quaisquer
! modip = modulo do vetor momento de dipolo de uma particula
! modrij = modulo do vetor R_{ij} utilizado para torna-lo unitario
! termo1, termo 2, temor 3 e termo 4 = termos utilizados para implementacao de forcas magneticas
! Eij = parametro utilizado para implementacao de forcas de contato entre particulas
! auxcont = constante auxiliar na implementacao das forcas de contato entre particulas
! tempo = tempo de simulacao
! ti,tf e tpros = tempo inicial, tempo final e tempo de processamento
! modrand = modulo do vetor randomico usado para implementar movimento browniano

! Variaveis logicas

! torque  = variavel logicas que determina se sera resolvida a equacao do torque
! estatistica = variavel logica que determina se sera feita analise estatistica em cima de varias realizacoes
! inercia = variavel logicas que determina se a particula possui inercia
! externo = variavel logica que verifica se existe ou nao um campo externo aplicado
! gravidade = variavel logica que liga ou desliga a gravidade
! browniano = variavel logica que liga ou desliga o movimento browniano
! agregado = variavel logica que determina se a partir do momento em que duas particulas
! formam um agregado elas podem se separar ou nao

end module variaveis
