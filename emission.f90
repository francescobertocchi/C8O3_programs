program aggregatex

implicit none

character jobz,uplo
integer i,lwork,info,n,lda,j,k
real*8 theta,z,pi,mu0,rr(3),wt,work(6000),w(2000),h(2000,2000),step,erre
real*8 mu(2000,2000),mut(2000,2000,3),rz(2000),ry(2000),rs(2000),g(2000),ma(2000,2000),mat(2000,2000,3),ma0,ospes(2000),rspes(2000)
real*8 refr_ind,f_mon,sum_g,sum_rs,mut_agg(2000,3),os(2000),rx(2000),x,d,cd,s,omega,gamma,sigma, wave, wave2, sig
real*8 ptot,kb,T

pi=dacos(-1.d0)

open (1,file='CPL.dat')
open (2,file='coordinates.dat')
open (3,file='fluorescence.dat')
open(4,file="positions.dat")
open(23,file="electric-dipoles.dat")
open(24,file="magnetic-dipoles.dat")
open(33,file="Eigenval.dat")
open(37,file="Strengths.dat")
open(26,file="Boltzmann.dat")
write(*,*) 'number of molecules and refractive index, wt e mu0 (D)'
read(*,*) n, refr_ind,wt,mu0


write(*,*) "magnetic dipole moment amplitude (D) and spectral bandwidth"
read (*,*) ma0,sigma

!!!boltzmann constant (eV/K)
kb=8.39d-5

write(*,*) "temperature (in K)"
read (*,*) T

!!!mu is the electric dipole moment
!!!r is the position vector
mu=0.d0
rx=0.d0
ry=0.d0
rz=0.d0

!reads geometries and dipoles from input file
do i=1,n

   read(4,*) rx(i),ry(i),rz(i)
   read(23,*) mu(i,1),mu(i,2),mu(i,3)
   read(24,*) ma(i,1),ma(i,2),ma(i,3)

   !!!mu0 is the electric dipole moment amplitude in D. mu are unit vectors.
   
   
   mu(i,1)=mu0*mu(i,1)
   mu(i,2)=mu0*mu(i,2)
   mu(i,3)=mu0*mu(i,3)

!!!ma0 is the magnetic dipole moment amplitude in D. ma are unit vectors.
   

ma(i,1)=ma(i,1)*ma0
ma(i,2)=ma(i,2)*ma0
ma(i,3)=ma(i,3)*ma0
!!!allineati
   enddo

 !!!build hamiltonian matrix 
h=0.d0

do i=1,n
   !!diagonal term: exciton energy
   h(i,i)=wt !!!wt in eV
  
   do j=i+1,n !!!upper triangle
      d=sqrt((rz(j)-rz(i))**2+(rx(j)-rx(i))**2+(ry(j)-ry(i))**2)!!!distance between the i-th and j-th molecules

      !!non diagonal terms: interaction energies (dipolar approximation)
 

       h(i,j)=(0.625/refr_ind)*(((mu(i,1)*mu(j,1)+mu(i,2)*mu(j,2)+mu(i,3)*mu(j,3))/d**3)&
            -3*((mu(i,1)*(rx(j)-rx(i))+&
            mu(i,2)*(ry(j)-ry(i))+mu(i,3)*(rz(j)-rz(i)))*(mu(j,1)*(rx(j)-rx(i))+mu(j,2)*(ry(j)-ry(i))+mu(j,3)*(rz(j)-rz(i))))/d**5)!
      write(*,*) i,j,h(i,j)
      
   enddo

enddo

jobz='v'
uplo='u'
lda=2000
lwork=6000

call dsyev (jobz,uplo,n,h,lda,w,work,lwork,info)
!!!hamiltonian matrix has been diagonalized

do i=1,n
   write(33,*) i,w(i) !!!energie del modello eccitonico tucur, non riscalate
enddo


!transition dipole moments
mut=0.d0
mat=0.d0
do i=1,n !!!excited states
   do k=1,n !!!molecules
     do j=1,3 !! x,y,z component
        mut(k,i,j)=mut(k,i,j)+mu(k,j)*h(k,i)
        mat(k,i,j)=mat(k,i,j)+ma(k,j)*h(k,i)
     enddo
  enddo
enddo

!compute rotational strengths
!!!
rs=0 !!!rotational strength
g=0  !!!geometrical factor
os=0 !!!oscillator strength
sum_rs=0
sum_g=0
mut_agg=0

do i=1,n !!!excited states
   do j=1,n !!!molecules
      do k=j,n !!!compute mu-mu coupling contributions
         rr(3)=(mut(j,i,1)*mut(k,i,2)-mut(k,i,1)*mut(j,i,2))

         rr(2)=(mut(j,i,3)*mut(k,i,1)-mut(k,i,3)*mut(j,i,1))

         rr(1)=(mut(j,i,2)*mut(k,i,3)-mut(k,i,2)*mut(j,i,3))

         
if(j.eq.k) then !!!intramoelcular mu-m coupling
rs(i)=rs(i)+mat(k,i,1)*mut(j,i,1)+mat(k,i,2)*mut(j,i,2)+mat(k,i,3)*mut(j,i,3)
endif
 if(j.ne.k) then !!!intermolecular mu-m coupling
rs(i)=rs(i)+mat(k,i,1)*mut(j,i,1)+mat(k,i,2)*mut(j,i,2)+mat(k,i,3)*mut(j,i,3)

            rs(i)=rs(i)+mut(k,i,1)*mat(j,i,1)+mut(k,i,2)*mat(j,i,2)+mut(k,i,3)*mat(j,i,3)
           endif

           !!!convert rotational strength (from mu-mu coupling) in D**2
         g(i)=g(i)-(rz(k)-rz(j))*rr(3)-(rx(k)-rx(j))*rr(1)-(ry(k)-ry(j))*rr(2)
          rs(i)=rs(i)-w(i)*((rz(k)-rz(j))*rr(3)+(rx(k)-rx(j))*rr(1)+(ry(k)-ry(j))*rr(2))*2.532*0.0001


      enddo
   enddo

   !!!compute oscillator strengths
   do k=1,n
      do j=1,3 !!!cartesian axes
         mut_agg(i,j)=mut_agg(i,j)+mut(k,i,j)
      enddo
   enddo
   os(i)=mut_agg(i,1)**2+mut_agg(i,2)**2+mut_agg(i,3)**2
   sum_rs=sum_rs+rs(i)
   sum_g=sum_g+g(i)
   write(37,22) i*1.d0,w(i),rs(i),g(i),os(i)

enddo

write(*,*) 'sum rotational strengths', sum_rs
write(*,*) 'sum rotational g factor', sum_g

!write(*,*) "kbT", (kb*T)
ptot=0.d0
do i=1,n

   !!!strengths are weighted by Boltzmann distribution
ptot=ptot+exp(-(w(i)-w(1))/(kb*T))

ospes(i)=os(i)*exp(-(w(i)-w(1))/(kb*T))
rspes(i)=rs(i)*exp(-(w(i)-w(1))/(kb*T))
enddo

ospes=ospes/ptot

do i=1,n
write(26,*) i,ospes(i),rspes(i)

enddo
!stop
!write (*,*) "ptot", ptot
omega=50
do i=1,1000
   omega=omega+1.0d0
   wave=1/omega
   cd=0.d0
   s=0.d0
   do j=1,n
      wave2=w(j)/1240 
      sig=sigma/1240
      cd=cd+5.742*1d-4*((w(j))**3.d0)*rspes(j)*(1.d0/(sig*sqrt(pi)))*exp(-1.d0*((wave-wave2)/sig)**2) 

      s=s+3.2082*1d6*((w(j))**3.d0)*ospes(j)*(1.d0/sig)*exp(-1.d0*((wave-wave2)/sig)**2) 
   enddo
   write(1,*) omega, cd
   write(3,*) omega, s

  !!!cd and s are CPL and emission, in arbitrary units
   
enddo

22 format (9(1x,e10.3))
endprogram aggregatex
