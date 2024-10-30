program kmc
implicit none
! KMC computes the population vs time of a number of species involved in several dynamical processes
! rate: rate constant of a given processes  
! p:    population of a given species (p0 is its initial value)
! re (pr): reactant (product) for a given process
character*80 :: title
integer, parameter :: dp = selected_real_kind(15,307)
integer :: m,nran,nr,nesp,inran,i,j,k,l,mu,pd,kk,ia,iz,ijk
real(dp) :: t,tmax,tprint,tint,a0,r2a0,suma,vol,avog,conv,kc,kp,R,temp,alpha,pressure,voll,alpha_eq
real(dp) :: rnd(2)
integer,dimension(:),allocatable :: re1,re2,pr1,pr2,n,p,p0
real (dp) ,dimension(:),allocatable :: a,rate,c0,c0read,cont
avog=6.022d23
R=0.08205
read(*,"(a80)") title
print "(t3,a80)",title
print "(t3,a,/)","+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
nran=1
read(*,*) m,nesp,temp
allocate(re1(m),re2(m),pr1(m),pr2(m),a(m),rate(m),p0(nesp),p(nesp),n(nesp),c0(nesp),c0read(nesp),cont(m))
n=(/ (l,l=1,nesp) /)
! The vol relates to population of species i p_i and concentration c_i as: vol=p_i/c_i/avog
! So, choose vol for a expected population of a given species i (p_i)
! vol  is volume in cm^3
! voll is volume in L
! conc is read in mol/L but in the code is in mol/cm^3
read(*,*) vol
voll=vol*1e-3
print "(t3,a,1p,e10.2,a,1p,e10.2,a,1p,e10.2,a)","Nanovessel of:",vol," cm^3",voll," L",vol*1e21," nm^3"
conv=vol*avog
print*, ""
print "(t3,a)"," Rxn  Rate(s-1) re1 re2 pr1 pr2"
do i=1,m
   read(*,*) rate(i),re1(i),re2(i),pr1(i),pr2(i)
   if(re1(i)/=0.and.re2(i)/=0) rate(i)=rate(i)/vol
   print "(t3,i4,1p,e10.2,4(i4))",i,rate(i),re1(i),re2(i),pr1(i),pr2(i)
   if(re1(i)==re2(i)) rate(i)=2*rate(i)
enddo
print*, ""
print "(t3,a)","Species C0(M)       molecules     moles"
do ijk=1,m
   cont(ijk)=0
enddo
do i=1,nesp
   read(*,*) c0read(i)
   c0(i)=c0read(i)*1e-3
   p0(i)=int(c0(i)*conv+0.5)
   print "(t3,i4,1p,e10.2,i15,1p,e10.2)",i,c0read(i),p0(i),p0(i)/avog
enddo
read(*,*) tmax,tint
ia=1
iz=nesp
open(unit=7,file="pop.csv")
open(unit=8,file="rate.csv")
print "(/,t3,a,1p,e10.2,a,/,t3,a,1p,e10.2,a,/)","Total time=",tmax," seconds","Step size =",tint," seconds"
big: do inran=1,nran
   print "(/t3,a/)","SIMULATINS START HERE"
   p=p0
   t=0.d0
   tprint=0.d0
   do while(tprint<tmax)
     do j=1,m
        if(re1(j)==0) then
! unimolecular reaction
          a(j)=rate(j)*p(re2(j))
        else if(re2(j)==0) then
          a(j)=rate(j)*p(re1(j))
        else
          if(re1(j)/=re2(j)) a(j)=rate(j)*p(re1(j))*p(re2(j))
          if(re1(j)==re2(j)) a(j)=rate(j)*p(re1(j))*(p(re2(j))-1)/2
        endif
     enddo
     a0=sum(a)
     call random_number(rnd)
     t=t-log(rnd(1))/a0      
     do while (t>=tprint) 
        print "(t3,a,1p,e10.2)","TIME=",tprint
        do ijk=1,nesp 
           if(p(ijk)>0) print "(t3,a,i4,a,i15)","Molecules of:",ijk,":",p(ijk)
        enddo
        write(7,*) tprint,",",p(1),",",p(2)
        print "(t3,a)","EQUILIBRIUM PROPERTIES"
        print "(t3,a)","    Kc        Kp         alpha   p(atm)  alpha(eq)"
        kc=real(p(2),dp)*real(p(2),dp)/real(p(1),dp)/vol
        kp=kc/avog*R*temp*1e3
        alpha=(real(p0(1),dp)-real(p(1),dp))/real(p0(1),dp)
        pressure=(real(p(1),dp)+real(p(2),dp))/avog*R*temp/voll
        alpha_eq=(kp/(kp+4*pressure))**0.5
        if(p(1)>0) print "(t3,(1p,e10.4),(1p,e10.4),e10.4,e10.4,e10.4)",kc,kp,alpha,pressure,alpha_eq
        print*, ""
        if(p(1)>0) write(8,*) tprint,",",a(1),",",a(2)
        tprint=tprint+tint
        if(tprint>tmax) cycle big
     enddo
     r2a0=rnd(2)*a0
     suma=0.d0
     s1: do mu=1,m
       suma=suma+a(mu)
       if(suma>=r2a0) exit s1
     enddo s1
     if(re1(mu)==0) then
! unimolecular reaction
        p(re2(mu))=p(re2(mu))-1
     else if(re2(mu)==0) then
        p(re1(mu))=p(re1(mu))-1
! unimolecular reaction
     else
        p(re1(mu))=p(re1(mu))-1
        p(re2(mu))=p(re2(mu))-1
     endif
     if(pr1(mu)==0) then
! only one product
        p(pr2(mu))=p(pr2(mu))+1
     else if(pr2(mu)==0) then
        p(pr1(mu))=p(pr1(mu))+1
     else 
        p(pr1(mu))=p(pr1(mu))+1
        p(pr2(mu))=p(pr2(mu))+1
     endif
     cont(mu)=cont(mu)+1  
   enddo
enddo big

close(unit=7)
close(unit=8)

end program kmc
