!********************************************************************************************************!
!		            SUBROTINA PRINCIPAL DO PROGRAMA MOBILITY				         !
!											  		 !
! Ultima atualizacao em 14/01/2016							    		 !
! 											    		 !
!											    		 !
! Escrito por: Rafael Gabler Gontijo, PhD 						    		 !
!********************************************************************************************************!


! Esta subrotina tem como objetivo resolver as equacoes do momento linear e angular para uma
! suspensao de particulas rigidas monodispersas, na qual as seguintes forcas estao presentes:

! - Forcas de repulsao entre particulas (lubrificacao);
! - Forcas de contato entre particulas (Hertz);
! - Forcas por interacao magnetica entre os momentos de dipolo das particulas;
! - Forcas por interacao hidrodinamicas - sistemas periodicos

! Para a equacao do momento angular, consideram-se os seguintes torques presentes:

! - Torques magneticos (por interacao entre particula e por campo externo);
! - Torque viscoso;

! As escalas caracteristicas para o processo de adimensionalizacao sao:

! - Velocidade de Stokes
! - Tempo (a/U_s) ou (a²/D)

! As forças e torques por interacoes magneticas podem ser computadas de forma periodica ou nao


subroutine principal

use variaveis
use funcoes

! Numero de passos de tempo e passo de tempo

pi = acos(-1.0)
dt_inicial=min(0.01, pi/freqcampo)


! Se a simulacao e estatica (Monte Carlo) o numero de passos de tempo e 2, caso contrario depende do
! tempo total de simulacao escolhido pelo usuario

if(estatica) then
npast=2
else
npast=tempo/dt_inicial
end if

! Quantidade de numeros randomicos necessarios em cada timestep

nnr=3*N*rea

!********************************************************************************************************!

! Alocando variaveis na memoria

allocate(trap(npast))
allocate(X(rea,N,3))
allocate(U(rea,N,3))
allocate(W(rea,N,3))
if(browniano)then
allocate(FORCAS(6,rea,N,3))
allocate(TORQUES(3,rea,N,3))
else
allocate(FORCAS(5,rea,N,3))
allocate(TORQUES(2,rea,N,3))
end if
allocate(FT(rea,N,3))
allocate(Tt(rea,N,3))
allocate(Di(rea,N,3))
allocate(aux1(rea,N))
allocate(aux2(rea,N))
allocate(aux3(rea,N))
allocate(aux4(rea,N))
allocate(nr(nnr))
allocate(dt(rea,N))
allocate(hidrodinamica_aux1(N,3))
allocate(hidrodinamica_aux2(N,3))
allocate(contribuicao_self(rea,N))
allocate(contribuicao_fisico(rea,N))
allocate(contribuicao_reciproco(rea,N))
allocate(hidro1(nb,3))
allocate(hidro2(nbr,3))
if(tmagper)then
allocate(auxt(N,3))
allocate(torquereal(nb,3))
allocate(torquereciproco(nbr,3))
allocate(cof4(2,10000))
allocate(cof5(2,10000))
allocate(cof7(2,10000))
end if
if(fmagper) then
allocate(cof6(2,10000))
allocate(cof8(2,10000))
allocate(auxf(N,3))
allocate(forcareal(nb,3))
allocate(forcareciproca(nbr,3))
end if
allocate(ILF(nb,3))
allocate(ILR(nbr,3))
allocate(XI(nb,rea,N,3))
allocate(cof1(2,10000))
allocate(cof2(2,10000))
allocate(cof3(2,10000))
if(leito)then
allocate(usistema(N,3))
end if
if(grafmag)then
allocate(magtempo(3,npast))
allocate(flutmag(N,rea))
end if
allocate(tempototal(npast))
if(agregado_inicial) then
allocate(centro_massa(rea,3))
end if
allocate(DIAM(rea,N))
allocate(beta(rea,N))
allocate(diarand(rea*N))
if(greenkubo)then
allocate(potencial(rea,N))
allocate(energia(rea,N))
allocate(energiaantes(rea,N))
allocate(heatcurrent(rea,N,3))
allocate(dedt(rea,N))
end if
allocate(campo(npast))
allocate(y(npast))


!********************************************************************************************************!
!********************************************** ZERANDO TUDO ********************************************!
!********************************************************************************************************!

509 FORMAT(F30.4,F30.4,F30.4,F30.4)
510 FORMAT(F30.4,F30.4,F30.4)
666 FORMAT(F30.4,F30.4,F30.4,F30.4,F30.4,F30.4,F30.4)
1012 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)
2024 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x)


X=0.0
U=0.0
aux1=0.0
aux2=0.0
aux3=0.0
Di=0.0
nr=0.0
FORCAS=0.0
FT=0.0
hidrodinamica_aux1=0.0
hidrodinamica_aux2=0.0
hidro1=0.0
hidro2=0.0
contribuicao_self=0.0
contribuicao_fisico=0.0
contribuicao_reciproco=0.0
torquereal=0.0
torquereciproco=0.0
auxt=0.0
shearratei=shearrate

! Definindo o numero de Peclet rotacional em termos do translacional utilizando a relacao entre os coe-
! ficientes de difusao de Stokes-Einstein translacionais e rotacionais

Per=(4.0/3.0)*Pe

! Definindo o passo de tempo de todas as particulas igual ao inicial

do j=1,rea
do i=1,N
dt(j,i)=dt_inicial
end do
end do

! Definindo o tamanho das particulas

if(polidispersidade)then

call randomica(1.5,2.5,diarand,(N*rea),8)

do j=1,rea
do i=1,N
DIAM(j,i)=diarand((j-1)*N + i)
end do
end do

else

do j=1,rea
do i=1,N
DIAM(j,i)=2.0
end do
end do

end if

do j=1,rea
do i=1,n
beta(j,i) = DIAM(j,i)/DIAM(j,1)
end do
end do


! Calculando o tamanho do box para atender a fracao volumetrica desejada

if(agregado_inicial) then

ragreg=(N/phi)**(1.0/3.0)
l=100.0*ragreg
h=l

else
 l=((N/(razao*phi))*(4.0*3.1416)/(3.0))**(1.0/3.0)
 razao2=razao
 h=razao2*l
end if


! Criacao dos arquivos necessarios para armazenamento de resultados

if(.not.continua) then
 call gera_arquivos(posicao,velocidade,rea)
end if

! Calculando o numero pi

 pi=acos(-1.0)

! Definindo o parametro de convergencia da soma periodica qsi

qsi=1.0*((pi**0.5)/((l*l*h)**(1.0/3.0)))

! Tabelando as funcoes utilizadas na computacao das interacoes periodicas e criando os indices das 
! Lattices periodicas (caso algum tipo de interacao considere sistemas periodicos)

if(periodicidade) then
 call tabelagreen(qsi,l,nb,nbr,h)
 call estrutura_periodica
end if

! Distribuindo agora os momentos de dipolo iniciais das particulas

if(.not.continua) then
 call distribui_dipolo(Di,rea,N)
end if

! Vamos agora realizar uma simulacao paralela para todas as realizacoes

! Distribuindo inicialmente as particulas

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                         GENERATING INITIAL CONDITIONS                      *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

 call condicao_inicial

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                   INITIAL CONDITIONS SUCCESSFULLY GENERATED                *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

! Caso estejamos trabalhando com uma excitação de campo oriunda da solução do oscilador não-linear de Duffing:

if(duffing) then

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*             DUFFING HARMONIC EXCITATION SUCCESSFULLY PRE-CALCULATED        *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

! Condições iniciais para o campo H e sua derivada H'=y

y(1)=0.0
campo(1)=1.0
tempototal(1)=0.0

write(300*rea,510) tempototal(1),campo(1),y(1)

! Solução por meio do método de Ruge-Kutta de 4 ordem da excitação H(t) e de sua derivada H'(t)

do k=2,npast

k1=C4*cos(freqcampo*tempototal(k-1)) - (C1*y(k-1) + C2*campo(k-1) + C3*(campo(k-1)**3.0))
k2=C4*cos(freqcampo*(tempototal(k-1)+ 0.5*dt_inicial))-(C1*(y(k-1)+ 0.5*dt_inicial*k1)+C2*campo(k-1)+C3*(campo(k-1)**3.0))
k3=C4*cos(freqcampo*(tempototal(k-1)+ 0.5*dt_inicial))-(C1*(y(k-1)+ 0.5*dt_inicial*k2)+C2*campo(k-1)+C3*(campo(k-1)**3.0))
k4=C4*cos(freqcampo*(tempototal(k-1)+dt_inicial))-(C1*(y(k-1)+ dt_inicial*k3)+C2*campo(k-1)+C3*(campo(k-1)**3.0))

g1=y(k-1)
g2=y(k-1) + dt_inicial*0.5*g1
g3=y(k-1) + dt_inicial*0.5*g2
g4=y(k-1) + dt_inicial*g3

y(k) = y(k-1) + (dt_inicial/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
campo(k) = campo(k-1) + (dt_inicial/6.0)*(g1 + 2.0*g2 + 2.0*g3 + g4)
tempototal(k) = tempototal(k-1) + dt_inicial 

write(300*rea,510) tempototal(k),campo(k),y(k)

end do
end if


! Caso estejamos trabalhando com uma excitação de campo do tipo batimento:

if(beating) then

! Condições iniciais para o campo H e sua derivada H'=y

y(1)=0.0
campo(1)=2.0
tempototal(1)=0.0

do k=2,npast
tempototal(k) = tempototal(k-1) + dt_inicial 
campo(k) = cos(freqcampo*tempototal(k)) + cos(freqbeat*tempototal(k))
y(k) = -freqcampo*sin(freqcampo*tempototal(k)) - freqbeat*sin(freqbeat*tempototal(k))
end do

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*           BEATING PATTERN FIELD EXCITATION SUCCESSFULLY PRE-CALCULATED     *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

end if

! No caso em que temos um campo oscilatorio e não estamos lidando com excitações do tipo
! duffing ou por padrões de batimento, estamos então no contexto de uma excitação harmônica do tipo seno

if(oscilacampo) then
if(.not.duffing) then
if(.not.beating)then

! Condições iniciais para o campo H e sua derivada H'=y

y(1)=freqcampo
campo(1)=0.0
tempototal(1)=0.0

! Calculando o intervalo de tempo de atuação de cada frequência (para construção dos diagramas de bifurcação)

intervalo=npast*dt_inicial/nfreq
multiplofreq=0

do k=2,npast

if(bifurcation) then
tempototal(k) = tempototal(k-1) + dt_inicial 
contfreqinteiro1=tempototal(k-1)/intervalo
contfreqinteiro2=tempototal(k)/intervalo
if(contfreqinteiro1.ne.contfreqinteiro2) then
multiplofreq=multiplofreq+1 
end if
campo(k)=sin((freqcampo+(freqmax-freqcampo)*multiplofreq)*tempototal(k)) 
y(k)=(freqcampo + (freqmax-freqcampo)*multiplofreq)*cos((freqcampo + (freqmax-freqcampo)*multiplofreq)*tempototal(k))
else
tempototal(k) = tempototal(k-1) + dt_inicial 
campo(k) = sin(freqcampo*tempototal(k)) 
y(k) = freqcampo*cos(freqcampo*tempototal(k))
end if

end do
multiplofreq=0

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                HARMONIC FIELD EXCITATION SUCCESSFULLY PRE-CALCULATED       *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

end if
end if
end if


!********************************************************************************************************!
!********************************* INICIO DO PROCESSO NUMERICO DE SIMULACAO *****************************!
!********************************************************************************************************!


if(continua)then
iter=iter
auxiliar_continua=npast
else
iter=1
auxiliar_continua=npast-1
end if


aux_real=auxiliar_continua

! Iniciando o processo global numerico

if(printphi)then
call campo_phi(rea,k)   
end if

! Aqui é o looping que faz o tempo evoluir

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                                SIMULATING                                  *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

do k=iter, auxiliar_continua

k_real=k

if(shear) then

if(oscillatory) then
shearrate=shearratei*sin(freq*dt_inicial*k)
end if

end if

! Calculando a energia das partículas antes de evoluir os potenciais e velocidades das particulas

if(greenkubo) then
do j=1,rea
do i=1,N
energiaantes(j,i) = potencial(j,i) !+ 0.5*(U(j,i,1)**2.0 + U(j,i,2)**2.0 + U(j,i,3)**2.0)**0.5
end do
end do
end if



! Primeira coisa de tudo: determinar o passo de tempo de cada particula antes de comecar a simulacao

do j=1,rea
do i=1,N
do q=1,N
if(i.ne.q) then
! Calcula-se a distancia entre a particula em questao e as outras particulas
r=(((X(j,i,1)-X(j,q,1))**2.0)+((X(j,i,2)-X(j,q,2))**2.0)+((X(j,i,3)-X(j,q,3))**2.0))**0.5
dt(j,i)=dt_inicial
dt(j,q)=dt_inicial
end if
end do
end do
end do

!********************************************************************************************************!
!********************************* Determinando as forcas Brownianas  ***********************************!
!********************************************************************************************************!

if(browniano)then
 call brownian
end if

!********************************************************************************************************!
!****************** Determinando as forcas de repulsao entre as particulas do box  **********************!
!********************************************************************************************************!

 call repulsao

!********************************************************************************************************!
!******************** Determinando as forcas de contato entre as particulas do box  *********************!
!********************************************************************************************************!

 call contato

!********************************************************************************************************!
!************** Determinando as forcas por interacao magnetica entre as particulas do box ***************!
!********************************************************************************************************!

if(.not.fmagper) then
 call forca_magnetica
end if

!********************************************************************************************************!
!************** Determinando as forcas por interacao magnetica entre as particulas e campo***************!
!********************************************************************************************************!

!if(externo) then
! call campo_externo
!else
do j=1,rea
do i=1,N
FORCAS(5,j,i,1)=0.0
FORCAS(5,j,i,2)=0.0
FORCAS(5,j,i,3)=0.0
end do
end do
!end if

! Impondo gravidade ou condicao de particulas neutrally buoyant

if(gravidade)then
do j=1,rea
do i=1,N
FORCAS(3,j,i,1)=0.0
FORCAS(3,j,i,2)=0.0
FORCAS(3,j,i,3)=-beta(j,i)**3.0
end do
end do
else
do j=1,rea
do i=1,N
FORCAS(3,j,i,1)=0.0
FORCAS(3,j,i,2)=0.0
FORCAS(3,j,i,3)=0.0
end do
end do
end if

! Calculando todas as forcas que atuam nas particulas

do j=1,rea
do i=1,N
FT(j,i,1)=FORCAS(1,j,i,1)+FORCAS(2,j,i,1)+FORCAS(3,j,i,1)+FORCAS(4,j,i,1)+FORCAS(5,j,i,1)+FORCAS(6,j,i,1)
FT(j,i,2)=FORCAS(1,j,i,2)+FORCAS(2,j,i,2)+FORCAS(3,j,i,2)+FORCAS(4,j,i,2)+FORCAS(5,j,i,2)+FORCAS(6,j,i,2)
FT(j,i,3)=FORCAS(1,j,i,3)+FORCAS(2,j,i,3)+FORCAS(3,j,i,3)+FORCAS(4,j,i,3)+FORCAS(5,j,i,3)+FORCAS(6,j,i,3)
end do
end do

!********************************************************************************************************!
!******************************* DETERMINANDO AS INTERACOES PERIODICAS  *********************************!
!********************************************************************************************************!

if(periodicidade) then
call intper
end if

if(.not.ligaih) then
if(inertia) then

do j=1,rea
do i=1,N
 call resvel(U(j,i,1),dt_inicial,St,FT(j,i,1))
 call resvel(U(j,i,2),dt_inicial,St,FT(j,i,2))
 call resvel(U(j,i,3),dt_inicial,St,FT(j,i,3))
end do
end do
else
U=FT
endif
end if

! Determinando a velocidade em caso de cisalhamento simples

if(shear)then
do j=1,rea
do i=1,N
U(j,i,2) = U(j,i,2) + shearrate*X(j,i,3)
end do
end do
end if


!***************************** MODELO DE LEITO FLUIDIZADO MAGNETICO ***********************************!

if(leito)then
do i=1,N
usistema(i,1)=sum(U(:,i,1))/rea
usistema(i,2)=sum(U(:,i,2))/rea
usistema(i,3)=sum(U(:,i,3))/rea
end do


do q=1,rea
do i=1,N
U(q,i,1)=U(q,i,1)-usistema(i,1)
U(q,i,2)=U(q,i,2)-usistema(i,2)
U(q,i,3)=U(q,i,3)-usistema(i,3)
end do
end do
end if

! Calculando a energia total (cinética + potencial) das partículas para uso em Green-Kubo 

if(greenkubo) then
do j=1,rea
do i=1,N
energia(j,i) = potencial(j,i) !+ 0.5*(U(j,i,1)**2.0 + U(j,i,2)**2.0 + U(j,i,3)**2.0)**0.5

dedt(j,i)= (energiaantes(j,i)-energia(j,i))/dt_inicial

heatcurrent(j,i,1) = energia(j,i)*U(j,i,1) + X(j,i,1)*dedt(j,i)
heatcurrent(j,i,2) = energia(j,i)*U(j,i,2) + X(j,i,2)*dedt(j,i)
heatcurrent(j,i,3) = energia(j,i)*U(j,i,3) + X(j,i,3)*dedt(j,i)
end do
end do
end if


!******************************************************************************************************!

! Calculando a posicao atual das particulas (Utilizando Euler)

do j=1,rea
 do i=1,N
 call respos(X(j,i,1),dt(j,i),U(j,i,1))
 call respos(X(j,i,2),dt(j,i),U(j,i,2))
 call respos(X(j,i,3),dt(j,i),U(j,i,3))
 end do
end do


!******************************************************************************************************!
! Implementando a condicao de contorno de periodicidade

!if(periodicidade) then
do j=1,rea
do i=1,N
if(X(j,i,1).gt.l) then
X(j,i,1)=X(j,i,1)-l
end if
if(X(j,i,1).lt.0.0)then
X(j,i,1)=l-X(j,i,1)
end if
if(X(j,i,2).gt.l) then
X(j,i,2)=X(j,i,2)-l
end if
if(X(j,i,2).lt.0.0)then
X(j,i,2)=l-X(j,i,2)
end if
if(X(j,i,3).gt.h) then
X(j,i,3)=X(j,i,3)-h
end if
if(X(j,i,3).lt.0.0)then
X(j,i,3)=h-X(j,i,3)
end if
end do
end do
!end if

!********************************************************************************************************!

! Escrevendo em um arquivo de saida a posicao e velocidade de cada particula em cada instante de tempo


if(continua)then

teste1=k/n3
teste2=k/n2

if(teste1.eq.teste2) then


if(posicao)then
do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt_inicial

write(*,*) 'Percentual:', (k_real/aux_real)*100, '%'

if(agregado_inicial) then
write(666,*) k*dt_inicial, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if


if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if


if(greenkubo)then
do j=1,rea
write(2012*j,*) 'zone t="',k,'"'
do i=1,N
write(2012*j,1012) potencial(j,i),energia(j,i),dedt(j,i),heatcurrent(j,i,1),heatcurrent(j,i,2),heatcurrent(j,i,3) 
end do
end do
end if


end if



else

if(k.eq.1) then

if(posicao)then
do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
if(gravadipolo)then
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
else
write(j,509)X(j,i,1),X(j,i,2),X(j,i,3),DIAM(j,i)
end if
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt_inicial

write(*,*) 'Percentual:', (k_real/aux_real)*100, '%'

if(agregado_inicial) then
write(666,*) k*dt_inicial, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if


if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if

end if

teste1=k/n3
teste2=k/n2

if(k.ne.1) then

if(teste1.eq.teste2) then
if(posicao)then
do j=1,rea
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
if(gravadipolo)then
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
else
write(j,509)X(j,i,1),X(j,i,2),X(j,i,3),DIAM(j,i)
end if
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt_inicial

write(*,*) 'Percentual:', (k_real/aux_real)*100, '%'  

if(agregado_inicial) then
write(666,*) k*dt_inicial, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if

if(greenkubo)then
do j=1,rea
write(2012*j,*) 'zone t="',k,'"'
do i=1,N
write(2012*j,1012) potencial(j,i),energia(j,i),dedt(j,i),heatcurrent(j,i,1),heatcurrent(j,i,2),heatcurrent(j,i,3)
end do
end do
end if



if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if
end if

end if
end if

 
 

!********************************************************************************************************!
!********************************************************************************************************!
!***************************** SOLUCAO DO MOVIMENTO ROTACIONAL DAS PARTICULAS ***************************!
!********************************************************************************************************!
!********************************************************************************************************!

if(torque) then

!*******************************************************************************************!
!**************** Determinando os torques magneticos sobre as particulas *******************!
!*******************************************************************************************!

 if(.not.tmagper)then
 call torque_magnetico
 end if

if(externo) then

if(rotating)then

call rotating_field(alpha, freqcampo*k*dt_inicial)

else

if(oscilacampo) then
 call torque_externo(alpha*campo(k*dt_inicial))
else
 call torque_externo(alpha)
end if
end if

end if

! Caso existam torques brownianos eles ja foram computados acima, na subrotina "brownian"

! Computando os torques totais que atuam sobre as particulas

do j=1,rea
do i=1,N
if(browniano) then
Tt(j,i,1)= TORQUES(1,j,i,1) + TORQUES(2,j,i,1) + TORQUES(3,j,i,1)
Tt(j,i,2)= TORQUES(1,j,i,2) + TORQUES(2,j,i,2) + TORQUES(3,j,i,2)
Tt(j,i,3)= TORQUES(1,j,i,3) + TORQUES(2,j,i,3) + TORQUES(3,j,i,3)
else
Tt(j,i,1)= TORQUES(1,j,i,1) + TORQUES(2,j,i,1) 
Tt(j,i,2)= TORQUES(1,j,i,2) + TORQUES(2,j,i,2) 
Tt(j,i,3)= TORQUES(1,j,i,3) + TORQUES(2,j,i,3) 
end if
end do
end do


! Resolvendo a velocidade angular

if(mistura)then
do j=1,rea
do i=(percentual*N)+1,N
 call resomega(W(j,i,1),dt(j,i),Str,Tt(j,i,1))
 call resomega(W(j,i,2),dt(j,i),Str,Tt(j,i,2))
 call resomega(W(j,i,3),dt(j,i),Str,Tt(j,i,3))
end do
end do
else
do j=1,rea
do i=1,N
 call resomega(W(j,i,1),dt(j,i),Str,Tt(j,i,1))
 call resomega(W(j,i,2),dt(j,i),Str,Tt(j,i,2))
 call resomega(W(j,i,3),dt(j,i),Str,Tt(j,i,3))
! call resomega_sem_inercia(W(j,i,1),Tt(j,i,1))
! call resomega_sem_inercia(W(j,i,2),Tt(j,i,2))
! call resomega_sem_inercia(W(j,i,3),Tt(j,i,3))
end do
end do
end if

if(shear)then
do j=1,rea
do i=1,N
W(j,i,1) = W(j,i,1) - shearrate*0.5
end do
end do
end if

! Evoluindo o vetor momento de dipolo magneticos das particulas

if(mistura)then
do j=1,rea
do i=(N*percentual)+1,N
 call evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,3),W(j,i,2),W(j,i,3),dt(j,i))
 call evoldip(Di(j,i,2),Di(j,i,3),Di(j,i,1),W(j,i,3),W(j,i,1),dt(j,i))
 call evoldip(Di(j,i,3),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt(j,i))
end do
end do
else
do j=1,rea
do i=1,N
 call evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,3),W(j,i,2),W(j,i,3),dt(j,i))
 call evoldip(Di(j,i,2),Di(j,i,3),Di(j,i,1),W(j,i,3),W(j,i,1),dt(j,i))
 call evoldip(Di(j,i,3),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt(j,i))
end do
end do
end if

! Normalizando os dipolos

do j=1,rea
if(mistura)then
do i=1,(percentual*N)
Di(j,i,1)=0.0
Di(j,i,2)=0.0
Di(j,i,3)=0.0
end do
do i=(percentual*N)+1,N
Di(j,i,1)=Di(j,i,1)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
Di(j,i,2)=Di(j,i,2)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
Di(j,i,3)=Di(j,i,3)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
end do
else
do i=1,N
Di(j,i,1)=Di(j,i,1)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
Di(j,i,2)=Di(j,i,2)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
Di(j,i,3)=Di(j,i,3)/((Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5)
end do
end if
end do

! Determinando a magnetizacao da suspensao, caso seja solicitado pelo usuario

if(grafmag) then
 call media_ativa(Di,N,rea,magtempo(1,k),1)
 call media_ativa(Di,N,rea,magtempo(2,k),2)
 call media_ativa(Di,N,rea,magtempo(3,k),3)
! Calculando a barra de erro da magnetizacao de equilibrio
do j=1,rea
do i=1,N
flutmag(i,j)=(Di(j,i,3)-magtempo(3,k))**2.0
end do
end do
erromag=((1.0/(N*rea))*sum(flutmag))**0.5


!write(5*rea,*)magtempo(3,k),k*dt_inicial, ((magtempo(k)-magtempo(k-1))/dt_inicial),  erromag

derivada1=(magtempo(1,k)-magtempo(1,k-1))/dt_inicial
derivada2=(magtempo(2,k)-magtempo(2,k-1))/dt_inicial
derivada3=(magtempo(3,k)-magtempo(3,k-1))/dt_inicial

write(5*rea,2024) campo(k),y(k), magtempo(1,k),magtempo(2,k),magtempo(3,k), derivada1, derivada2, derivada3, k*dt_inicial

contfreqinteiro1= ((k-1)*dt_inicial)/intervalo
contfreqinteiro2= (k*dt_inicial)/intervalo

if(contfreqinteiro1.ne.contfreqinteiro2) then
multiplofreq=multiplofreq+1
frequencia=(freqcampo+(freqmax-freqcampo)*multiplofreq)
end if

! Escrever aqui na unidade certa um arquivo mag_tempo para cada frequencia
 
write(400*rea+multiplofreq+1,2024) campo(k),y(k), magtempo(1,k),magtempo(2,k),magtempo(3,k), derivada1, derivada2, derivada3, k*dt_inicial

end if

tempototal(k)=k*dt_inicial
end if
end do

do k=2,npast
trap(k-1)= ((magtempo(3,k-1)*y(k-1)) + (magtempo(3,k)*y(k)))*dt_inicial*0.5
end do

trapezio=sum(trap)

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                                 MHT REPORT	                              *'
print *,'*                                                                            *'
print *,'******************************************************************************'
write(*,*)''
write(*,*)'INTEGRAL OF M.dH:',abs(trapezio)



if(printphi) then
call campo_phi(rea,k)
end if



if(fator)then
 call fator_estrutura(X,N,l,h,dt_inicial,rea)
end if



! Fechando os arquivos abertos
do j=1,2*rea
 close(j)
end do

 close(100*rea)
 close(300*rea)

if(greenkubo) then
do j=1,rea
 close(2012*j)
end do
end if

! Desalocando todas as matrizes e vetores

deallocate(X, STAT = DeAllocateStatus)
deallocate(U, STAT = DeAllocateStatus)
deallocate(FORCAS, STAT = DeAllocateStatus)
deallocate(FT, STAT = DeAllocateStatus)
deallocate(nr, STAT = DeAllocateStatus)
deallocate(dt, STAT = DeAllocateStatus)
deallocate(hidrodinamica_aux1, STAT = DeAllocateStatus)
deallocate(hidrodinamica_aux2, STAT = DeAllocateStatus)
deallocate(hidro1, STAT = DeAllocateStatus)
deallocate(hidro2, STAT = DeAllocateStatus)
deallocate(ILF, STAT = DeAllocateStatus)
deallocate(ILR, STAT = DeAllocateStatus)
deallocate(XI, STAT = DeAllocateStatus)
deallocate(Tt, STAT = DeAllocateStatus)
deallocate(Di, STAT = DeAllocateStatus)
deallocate(aux1, STAT = DeAllocateStatus)
deallocate(aux2, STAT = DeAllocateStatus)
deallocate(aux3, STAT = DeAllocateStatus)
deallocate(aux4, STAT = DeAllocateStatus)
deallocate(contribuicao_self, STAT = DeAllocateStatus)
deallocate(contribuicao_fisico, STAT = DeAllocateStatus)
deallocate(contribuicao_reciproco, STAT = DeAllocateStatus)
if(tmagper)then
deallocate(auxt, STAT = DeAllocateStatus)
deallocate(torquereal, STAT = DeAllocateStatus)
deallocate(torquereciproco, STAT = DeAllocateStatus)
deallocate(cof4, STAT = DeAllocateStatus)
deallocate(cof5, STAT = DeAllocateStatus)
deallocate(cof7, STAT = DeAllocateStatus)
end if
if(fmagper) then
deallocate(cof6, STAT = DeAllocateStatus)
deallocate(cof8, STAT = DeAllocateStatus)
deallocate(auxf, STAT = DeAllocateStatus)
deallocate(forcareal, STAT = DeAllocateStatus)
deallocate(forcareciproca, STAT = DeAllocateStatus)
end if
deallocate(ILF, STAT = DeAllocateStatus)
deallocate(ILR, STAT = DeAllocateStatus)
deallocate(XI, STAT = DeAllocateStatus)
deallocate(cof1, STAT = DeAllocateStatus)
deallocate(cof2, STAT = DeAllocateStatus)
deallocate(cof3, STAT = DeAllocateStatus)
if(leito)then
deallocate(usistema, STAT = DeAllocateStatus)
end if
if(grafmag)then
deallocate(magtempo, STAT = DeAllocateStatus)
deallocate(flutmag, STAT = DeAllocateStatus)
end if
deallocate(tempototal, STAT = DeAllocateStatus)
if(agregado_inicial) then
deallocate(centro_massa, STAT = DeAllocateStatus)
end if
deallocate(DIAM, STAT = DeAllocateStatus)
deallocate(beta, STAT = DeAllocateStatus)
deallocate(diarand, STAT = DeAllocateStatus)
if(greenkubo)then
deallocate(potencial, STAT = DeAllocateStatus)
deallocate(energia, STAT = DeAllocateStatus)
deallocate(energiaantes, STAT = DeAllocateStatus)
deallocate(heatcurrent, STAT = DeAllocateStatus)
deallocate(dedt, STAT = DeAllocateStatus)
end if

write(*,*) ''

write(*,*) 'Finalizando o modulo de processamento...'

write(*,*) ''

end subroutine principal
