subroutine saida
use variaveis
use funcoes

write(*,*) 'Iniciando o modulo de pos-processamento de dados...'
write(*,*) ''
! Alocando algumas variaveis

npast=npast/n2

ncor=((1+(N-1))*(N-1))/2

allocate(U(rea,N,3))
allocate(V(rea,N,6))
allocate(X(rea,N,3))
allocate(flut(N,npast,rea,3))
allocate(velmedia(npast,rea,3))
allocate(vmedia(npast,3))
allocate(errovmedia(npast,3))
allocate(auxcor(rea,npast,3))
allocate(auxiliarcor(rea,3))
allocate(errocor(npast,3))
allocate(funcaor(npast,3))
allocate(dif(npast,3))
allocate(aux_erro_vel(npast,rea,3))
allocate(aux_erro_var(npast,rea,6))
allocate(variancia(npast,rea,6))
allocate(var(npast,6))
allocate(errovar(npast,6))
if(greenkubo)then
allocate(potencial(rea,N))
allocate(energia(rea,N))
allocate(energiaantes(rea,N))
allocate(heatcurrent(rea,N,3))
allocate(heataux(rea,npast,3))
allocate(dedt(rea,N))
allocate(jmedio(npast,3))
allocate(auxj(npast,3))
allocate(intj(npast,3))
end if

509 FORMAT(F30.4,F30.4,F30.4)
510 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x)
513 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)
514 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x)
1012 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)

!********************************** VELOCIDADE MEDIA ***************************************!

! Para uma suspensao com varias particulas, temos que a velocidade media em cada time-step e dada
! como uma media da velocidade de todas particulas e em todas as realizacoes para cada passo de
! tempo...

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocidade'//rea_char//'.plt',STATUS='OLD')
end do

do j=1,rea
read (rea+j,'(A)') linha1
do k=1,npast-1
read(rea+j,'(A)') linha2
do i=1,N
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)
end do

! Determinando a velocidade media em cada time-step e em cada realizacao (para cada direcao)
! considerando somente as particulas da parte ativa da suspensao

!if(mistura)then

!do i=N*percentual+1,N
!U(j,i,1)=0.0
!U(j,i,2)=0.0
!U(j,i,3)=0.0
!end do


!velmedia(k,j,1)= sum(U(j,:,1))/(N*percentual)
!velmedia(k,j,2)= sum(U(j,:,2))/(N*percentual)
!velmedia(k,j,3)= sum(U(j,:,3))/(N*percentual)


!else
velmedia(k,j,1)= sum(U(j,:,1))/N
velmedia(k,j,2)= sum(U(j,:,2))/N
velmedia(k,j,3)= sum(U(j,:,3))/N
!end if

end do
end do



! Calculando agora a velocidade media em cada time-step, tirando uma media em cima das
! realizacoes do valor de velmedia(k,j,direcao)

do k=1,npast-1
vmedia(k,1)=sum(velmedia(k,:,1))/rea
vmedia(k,2)=sum(velmedia(k,:,2))/rea
vmedia(k,3)=sum(velmedia(k,:,3))/rea
end do


! Calculando o erro-bar da velocidade media

do k=1,npast-1
do j=1,rea
aux_erro_vel(k,j,1)=(velmedia(k,j,1)-vmedia(k,1))**2.0
aux_erro_vel(k,j,2)=(velmedia(k,j,2)-vmedia(k,2))**2.0
aux_erro_vel(k,j,3)=(velmedia(k,j,3)-vmedia(k,3))**2.0
end do
end do

do k=1,npast-1
errovmedia(k,1)=((1.0/(rea))*sum(aux_erro_vel(k,:,1)))**0.5
errovmedia(k,2)=((1.0/(rea))*sum(aux_erro_vel(k,:,2)))**0.5
errovmedia(k,3)=((1.0/(rea))*sum(aux_erro_vel(k,:,3)))**0.5
end do


! Fechando os arquivos
do j=1,rea
 close(rea+j)
end do



! Escrevendo em um arquivo de saida a velocidade media e seu errobar

i=2*rea+1
open (i,file='velocidade_media.plt')
write(i,*) 'Variables="U","V","W","UMEDIA","DU","DV","DW","T"'

do k=1,npast-1
write(i,510)vmedia(k,1),vmedia(k,2),vmedia(k,3),   &
(((vmedia(k,1)**2.0)+(vmedia(k,2)**2.0)+(vmedia(k,3)**2.0))**0.5),   &
errovmedia(k,1),errovmedia(k,2),errovmedia(k,3),k*dt_inicial*n2
end do


write(*,*) 'Analise estatistica em cima da velocidade media - OK'

!*******************************************************************************************!


!********************************************** VARIANCIA **********************************!
do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocidade'//rea_char//'.plt',STATUS='OLD')
end do


do j=1,rea
read (rea+j,'(A)') linha1
do k=1,npast-1
read(rea+j,'(A)') linha2
do i=1,N
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)

! Determinando a variancia em cada time-step e em cada realizacao (para cada direcao)
! considerando somente as particulas da parte ativa da suspensao


if(mistura)then

if(i.le.(N*percentual)) then

V(j,i,1) = (U(j,i,1)-velmedia(k,j,1))**2.0
V(j,i,2) = (U(j,i,1)-velmedia(k,j,1))*(U(j,i,2)-velmedia(k,j,2))
V(j,i,3) = (U(j,i,1)-velmedia(k,j,1))*(U(j,i,3)-velmedia(k,j,3))

V(j,i,4) = (U(j,i,2)-velmedia(k,j,2))**2.0
V(j,i,5) = (U(j,i,2)-velmedia(k,j,2))*(U(j,i,3)-velmedia(k,j,3))
V(j,i,6) = (U(j,i,3)-velmedia(k,j,3))**2.0

else

V(j,i,1) = 0.0
V(j,i,2) = 0.0
V(j,i,3) = 0.0

V(j,i,4) = 0.0
V(j,i,5) = 0.0
V(j,i,6) = 0.0

end if

else

V(j,i,1) = (U(j,i,1)-velmedia(k,j,1))**2.0
V(j,i,2) = (U(j,i,1)-velmedia(k,j,1))*(U(j,i,2)-velmedia(k,j,2))
V(j,i,3) = (U(j,i,1)-velmedia(k,j,1))*(U(j,i,3)-velmedia(k,j,3))

V(j,i,4) = (U(j,i,2)-velmedia(k,j,2))**2.0
V(j,i,5) = (U(j,i,2)-velmedia(k,j,2))*(U(j,i,3)-velmedia(k,j,3))
V(j,i,6) = (U(j,i,3)-velmedia(k,j,3))**2.0


end if


end do

if(mistura) then

variancia(k,j,1)=sum(V(j,:,1))/(N*percentual)
variancia(k,j,2)=sum(V(j,:,2))/(N*percentual)
variancia(k,j,3)=sum(V(j,:,3))/(N*percentual)
variancia(k,j,4)=sum(V(j,:,4))/(N*percentual)
variancia(k,j,5)=sum(V(j,:,5))/(N*percentual)
variancia(k,j,6)=sum(V(j,:,6))/(N*percentual)

else

variancia(k,j,1)=sum(V(j,:,1))/N
variancia(k,j,2)=sum(V(j,:,2))/N
variancia(k,j,3)=sum(V(j,:,3))/N
variancia(k,j,4)=sum(V(j,:,4))/N
variancia(k,j,5)=sum(V(j,:,5))/N
variancia(k,j,6)=sum(V(j,:,6))/N


end if


V=0.0

end do
end do

! Calculando agora a variancia em cada time-step, tirando uma media em cima das
! realizacoes do valor de variancia(k,j,direcao)

do k=1,npast-1
var(k,1)=sum(variancia(k,:,1))/rea
var(k,2)=sum(variancia(k,:,2))/rea
var(k,3)=sum(variancia(k,:,3))/rea
var(k,4)=sum(variancia(k,:,4))/rea
var(k,5)=sum(variancia(k,:,5))/rea
var(k,6)=sum(variancia(k,:,6))/rea
end do


! Calculando o erro-bar da variancia

do k=1,npast-1
do j=1,rea
aux_erro_var(k,j,1)=(variancia(k,j,1)-var(k,1))**2.0
aux_erro_var(k,j,2)=(variancia(k,j,2)-var(k,2))**2.0
aux_erro_var(k,j,3)=(variancia(k,j,3)-var(k,3))**2.0
aux_erro_var(k,j,4)=(variancia(k,j,4)-var(k,4))**2.0
aux_erro_var(k,j,5)=(variancia(k,j,5)-var(k,5))**2.0
aux_erro_var(k,j,6)=(variancia(k,j,6)-var(k,6))**2.0
end do
end do

do k=1,npast-1
errovar(k,1)=((1.0/(rea))*sum(aux_erro_var(k,:,1)))**0.5
errovar(k,2)=((1.0/(rea))*sum(aux_erro_var(k,:,2)))**0.5
errovar(k,3)=((1.0/(rea))*sum(aux_erro_var(k,:,3)))**0.5
errovar(k,4)=((1.0/(rea))*sum(aux_erro_var(k,:,4)))**0.5
errovar(k,5)=((1.0/(rea))*sum(aux_erro_var(k,:,5)))**0.5
errovar(k,6)=((1.0/(rea))*sum(aux_erro_var(k,:,6)))**0.5
end do

! Fechando os arquivos
do j=1,rea
 close(rea+j)
end do


! Escrevendo em um arquivo de saida a variancia e seu errobar

i=2*rea+2
open (i,file='variancia.plt')
write(i,*) 'Variables="V11","V12","V13","V22","V23","V33","DV11","DV12","DV13","DV22","DV23","DV33","T"'


do k=2,npast-1
write(i,510)var(k,1),var(k,2),var(k,3),   &
var(k,4),var(k,5),var(k,6),   &
errovar(k,1),errovar(k,2),errovar(k,3),errovar(k,4),errovar(k,5),errovar(k,6),k*dt_inicial*n2
end do


write(*,*) 'Analise estatistica em cima da variancia - OK'

! Abrindo novamente os arquivos

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocidade'//rea_char//'.plt',STATUS='OLD')
end do

! Calculando a flutuacao em cada time-step e realizacao para cada particula

do j=1,rea
read (rea+j,'(A)') linha1
do k=1,npast-1
read(rea+j,'(A)') linha2

do i=1,N
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)

if(mistura)then

if(i.le.(N*percentual))then

flut(i,k,j,1)= U(j,i,1)-velmedia(k,j,1)
flut(i,k,j,2)= U(j,i,2)-velmedia(k,j,2)
flut(i,k,j,3)= U(j,i,3)-velmedia(k,j,3)

else

flut(i,k,j,1)= 0.0
flut(i,k,j,2)= 0.0
flut(i,k,j,3)= 0.0

end if

else
flut(i,k,j,1)= U(j,i,1)-velmedia(k,j,1)
flut(i,k,j,2)= U(j,i,2)-velmedia(k,j,2)
flut(i,k,j,3)= U(j,i,3)-velmedia(k,j,3)

end if

end do

end do
end do

! Fechando os arquivos
do j=1,rea
 close(rea+j)
end do

do j=1,rea
do i=1,npast

do k=1,npast-i
auxcor(j,k,1)=sum((flut(:,k,j,1)*flut(:,k+i-1,j,1)))/sum((flut(:,k,j,1)**2.0))
auxcor(j,k,2)=sum((flut(:,k,j,2)*flut(:,k+i-1,j,2)))/sum((flut(:,k,j,2)**2.0))
auxcor(j,k,3)=sum((flut(:,k,j,3)*flut(:,k+i-1,j,3)))/sum((flut(:,k,j,3)**2.0))
end do

! Essa funcao "r" já é a função autocorrelação naquele tempo R(t)

funcaor(i,1)=sum(auxcor(:,:,1))/((npast-i+1))
funcaor(i,2)=sum(auxcor(:,:,2))/((npast-i+1))
funcaor(i,3)=sum(auxcor(:,:,3))/((npast-i+1))

! Em seguida a gente estima a barra de erro (desvio padrão) da função autocorrelação

auxiliarcor(j,1)=(funcaor(i,1)-(auxcor(j,i,1)/(npast-i+1)))**2.0
auxiliarcor(j,2)=(funcaor(i,2)-(auxcor(j,i,2)/(npast-i+1)))**2.0
auxiliarcor(j,3)=(funcaor(i,3)-(auxcor(j,i,3)/(npast-i+1)))**2.0

errocor(i,1)=((1.0/rea)*sum(auxiliarcor(:,1)))**0.5
errocor(i,2)=((1.0/rea)*sum(auxiliarcor(:,2)))**0.5
errocor(i,3)=((1.0/rea)*sum(auxiliarcor(:,3)))**0.5


auxcor=0.0
end do
end do


! Escrevendo em um arquivo a funcao autocorrelacao das flutuacoes

i=2*rea+3
open (i,file='autocorrelacao.plt')
write(i,*) 'Variables="R1","R2","R3","DR1","DR2","DR3","T"'

do k=1,npast-1
write(i,510)funcaor(k,1),funcaor(k,2),funcaor(k,3),errocor(k,1),errocor(k,2),errocor(k,3),k*dt_inicial*n2
end do

write(*,*) 'Analise estatistica em cima da funcao autocorrelacao - OK'


!*******************************************************************************************!

!**************************************** COEFICIENTE DE DIFUSAO ***************************!


! Determinando o coeficiente de difusao

dif(1,1)=0.0
dif(1,2)=0.0
dif(1,3)=0.0

do k=2,npast-3
dif(k,1)=dif(k-1,1)+((funcaor(k,1)+funcaor(k-1,1))*(n2*dt_inicial)/2.0)
dif(k,2)=dif(k-1,2)+((funcaor(k,2)+funcaor(k-1,2))*(n2*dt_inicial)/2.0)
dif(k,3)=dif(k-1,3)+((funcaor(k,3)+funcaor(k-1,3))*(n2*dt_inicial)/2.0)
end do


! Escrevendo em um arquivo o coeficiente de difusao

i=2*rea+4
open (i,file='difusao.plt')
write(i,*) 'Variables="D1","D2","D3","T"'

do k=1,npast-1
write(i,510)dif(k,1),dif(k,2),dif(k,3),k*dt_inicial*n2
end do


write(*,*) 'Analise estatistica em cima do coeficiente de difusao - OK'

!*******************************************************************************************!

!**************************************** CONDUTIVIDADE TERMICA DE PARTICULA ***************************!

if(greenkubo) then
do i=1,rea
write(rea_char, '(I3)') i
open (2012*i,file='energia'//rea_char//'.plt',STATUS='OLD')
end do
do j=1,rea
read (2012*j,'(A)') linha1
do k=1,npast-1
read(2012*j,'(A)') linha2
do i=1,N
read(2012*j,1012) potencial(j,i),energia(j,i),dedt(j,i),heatcurrent(j,i,1),heatcurrent(j,i,2),heatcurrent(j,i,3) 
end do
jmedio(k,1)= sum(heatcurrent(:,:,1))/(N*rea)
jmedio(k,2)= sum(heatcurrent(:,:,2))/(N*rea)
jmedio(k,3)= sum(heatcurrent(:,:,3))/(N*rea)

heataux(j,k,1)=sum(heatcurrent(j,:,1))/N
heataux(j,k,2)=sum(heatcurrent(j,:,2))/N
heataux(j,k,3)=sum(heatcurrent(j,:,3))/N

end do
end do

do i=1,npast
auxj(i,1)=sum(heataux(:,i,1)*heataux(:,1,1))/rea
auxj(i,2)=sum(heataux(:,i,2)*heataux(:,1,2))/rea
auxj(i,3)=sum(heataux(:,i,3)*heataux(:,1,3))/rea
end do

! Essa funcao "auxj" já é o kernel da integral que levará ao cálculo da condutividade térmica. Vamos escrever esse kernel em um arquivo de saída

i=2*rea+5001
open (i,file='heatcurrent_kernel.plt')
write(i,*) 'Variables="K1","K2","K3","t"'

do k=1,npast-1
write(i,510)auxj(k,1),auxj(k,2),auxj(k,3),k*dt_inicial*n2
end do

! Vamos agora calcular a integral desse sinal que no limite em que a nossa estatística é feita em cima de um long-time deve saturar para o valor da condutividade térmica de partícula da suspensão 

intj(1,1)=0.0
intj(1,2)=0.0
intj(1,3)=0.0

do k=2,npast-3
intj(k,1)=intj(k-1,1)+((auxj(k,1)+auxj(k-1,1))*(n2*dt_inicial)/2.0)
intj(k,2)=intj(k-1,2)+((auxj(k,2)+auxj(k-1,2))*(n2*dt_inicial)/2.0)
intj(k,3)=intj(k-1,3)+((auxj(k,3)+auxj(k-1,3))*(n2*dt_inicial)/2.0)
end do

i=2*rea+5002
open (i,file='integral_heat_current.plt')
write(i,*) 'Variables="I1","I2","I3","t"'

do k=1,npast-1
write(i,510)intj(k,1),intj(k,2),intj(k,3),k*dt_inicial*n2
end do


end if

! Escrevendo em um arquivo os valores medios das componentes do vetor J(t)

i=2*rea+5000
open (i,file='heatcurrent.plt')
write(i,*) 'Variables="J1","J2","J3","t"'

do k=1,npast-1
write(i,510)jmedio(k,1),jmedio(k,2),jmedio(k,3),k*dt_inicial*n2
end do


!*******************************************************************************************!



write(*,*) 'Geracao dos arquivos de saida - OK'



if(grafmag)then
 CALL SYSTEM('gnuplot -persist "script4.gnu"')
else
 CALL SYSTEM('gnuplot -persist "script5.gnu"')
end if


! Desalocando as variaveis para liberar espaco na memoria

write(*,*) ''
write(*,*) 'Desalocando as matrizes criadas no modulo de pos-processamento...'
write(*,*) ''


deallocate(U)
write(*,*) 'Desalocando matriz 1 - OK'
deallocate(V)
write(*,*) 'Desalocando matriz 2 - OK'
deallocate(X)
write(*,*) 'Desalocando matriz 3 - OK'
deallocate(flut)
write(*,*) 'Desalocando matriz 4 - OK'
deallocate(velmedia)
write(*,*) 'Desalocando matriz 5 - OK'
deallocate(vmedia)
write(*,*) 'Desalocando matriz 6 - OK'
deallocate(errovmedia)
write(*,*) 'Desalocando matriz 7 - OK'
deallocate(auxcor)
write(*,*) 'Desalocando matriz 8 - OK'
deallocate(funcaor)
write(*,*) 'Desalocando matriz 9 - OK'
deallocate(dif)
write(*,*) 'Desalocando matriz 10 - OK'
deallocate(aux_erro_vel)
write(*,*) 'Desalocando matriz 11- OK'
deallocate(aux_erro_var)
write(*,*) 'Desalocando matriz 12- OK'
deallocate(variancia)
write(*,*) 'Desalocando matriz 13- OK'
deallocate(var)
write(*,*) 'Desalocando matriz 14- OK'
deallocate(errovar)
write(*,*) 'Desalocando matriz 15- OK'



end subroutine saida
