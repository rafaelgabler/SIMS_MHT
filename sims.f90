program magpart

use variaveis

! Titulo do programa e apresentacao do mesmo

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                  SIMS - SIMULATION OF MAGNETIC SUSPENSIONS	              *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                    NON-LINEAR MAGNETIC HYPERTHERMIA MODULE                 *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                     PROF. RAFAEL GABLER GONTIJO, PhD                       *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                        IN DEVELOPMENT SINCE 2009                           *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                         LAST UPDATE: 24/06/2022                            *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''
print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*       Numerical simulation of magnetic suspensions of hard spheres         *'  
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                       Langevin and Stokesian Dynamics                      *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''
!**************************************************************************************************************************************!
 call entrada

! Começa a contar o tempo de simulacao (ativa o cronometro)

 call cpu_time(ti)

! Chama a subrotina principal que executa de fato a simulacao

 call principal

if(estatistica) then
 call saida
end if

! Para o cronometro

 call cpu_time(tf)

! Calcula o tempo de processamento

 tpros=tf-ti


print *, 'TOTAL SIMULATION TIME:',tpros,'SECONDS'

write(*,*) ''

write(*,*) ''


end program magpart
