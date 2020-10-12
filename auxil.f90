!------------------------------------------------------------------------
! Some Auxiliary subroutines:
!   
!      print_orbs
!      print_NO
!      zbisect
!      zbrent_so
!      zbrent_gr
!      times_init
!      print_times
!      consistency
!      MathOut
!      normal_test
!      eigenvector_test
!      orthog_vec
!      symmet_orthog
!      Akima
! zbrent... is from the Numerical recipes. Akima is public domain
!           The diagon_NAGrest created by N. N. Lathiotakis and I.Theophilou
!------------------------------------------------------------------------
subroutine print_orbs(nb,nat,vecnat)

!..Global
   use params_general
   implicit none

!..Arguments
   integer, intent(in) :: nb,nat
   integer, intent(in) :: vecnat(lnbasis,lnnatorb)

!..Local variables
   integer :: i,ia

   do ia=1,nat
      write(50+ia,'(e20.10)')(vecnat(i,ia),i=1,nb)
      write(50+ia,*)
!     write(6,'(f20.10)')(vecnat(i,ia),i=1,nb)
!     write(6,*)
   enddo

end subroutine print_orbs

!------------------------------------------------------------------------
subroutine print_NO(nnatorb)

!..Global
   use matrices
   implicit none

!..Arguments
   integer, intent(in) :: nnatorb 

!..Local Variables
   integer :: ia,ib

   print*,'-------------------------------------------'
   Print*, 'HcoreNO:'
   do ia=1,nnatorb
      do ib=1,nnatorb
         print*,ia,ib,HcoreNO(ia,ib)
      enddo
   enddo
   print*,'-------------------------------------------'
   Print*, 'CoulNO:'
   do ia=1,nnatorb
      do ib=1,nnatorb
         print*,ia,ib,CoulNO(ia,ib)
      enddo
   enddo
   print*,'-------------------------------------------'
   Print*, 'ExchNO:'
   do ia=1,nnatorb
      do ib=1,nnatorb
         print*,ia,ib,ExchNO(ia,ib)
      enddo
   enddo
   print*,'-------------------------------------------'
end subroutine print_NO
!------------------------------------------------------------------------

function zbisect(func,a1,a2,ITMAX,tol,it_done,occnumtmp)
!------------------------------------------------------------------------
! BISECTION METHOD FOR ROOT FINDING
!    func: function the root of which will be found.
!  a1, a2: lower, upper bound.
!   ITMAX: Maximum number of iterations.
!     tol: tolerance.
! it_done: final number of iterations.
! The bisection or false position methods ca be chosen by 
! commenting in/out the eppropriate lines below.
!------------------------------------------------------------------------

!..Global
   use params_general
   implicit none

!..The function
   real(dp) :: zbisect

!..Arguments
   integer, intent(in) :: ITMAX
   integer, intent(out) :: it_done
   real(dp), intent(in) :: a1, a2, tol
   real(dp), external :: func
   real(dp) :: occnumtmp(lnbasis,3)

!..Local
   integer :: iter
   real(dp) :: x1, x2, g1, g2, xm, gm

   x1 = a1
   x2 = a2

   g1 = func(x1,occnumtmp)
   g2 = func(x2,occnumtmp)

   iter=1
   do iter=1, ITMAX

!..Bisection rule:
   xm = (x1 + x2)/2._dp

!..False position method:
!  xm = (x1*g2 - x2*g1)/(g2 - g1)

   gm=func(xm,occnumtmp)
      if(gm*g1.lt.0._dp) then
         x2 = xm
         g2 = gm
      else
         x1 = xm
         g1 = gm
      endif
      if(abs((x2-x1)/x2) < tol) exit
!     if(abs(gm) < tol) goto 148
   enddo !iter  
   it_done = iter
   zbisect = xm
end function zbisect

!------------------------------------------------------------------------
function zbisect_pot(func,a1,a2,ITMAX,tol,it_done,X_0)
!------------------------------------------------------------------------
! BISECTION METHOD FOR ROOT FINDING
!    func: function the root of which will be found.
!  a1, a2: lower, upper bound.
!   ITMAX: Maximum number of iterations.
!     tol: tolerance.
! it_done: final number of iterations.
! The bisection or false position methods ca be chosen by 
! commenting in/out the eppropriate lines below.
!------------------------------------------------------------------------

!..Global
   use params_general
   implicit none

!..The function
   real(dp) :: zbisect_pot

!..Arguments
   integer, intent(in) :: ITMAX
   integer, intent(out) :: it_done
   real(dp), intent(in) :: a1, a2, tol
   real(dp), external :: func
   real(dp) :: X_0(*)

!..Local
   integer :: iter
   real(dp) :: x1, x2, g1, g2, xm, gm

   x1 = a1
   x2 = a2

   g1 = func(x1, X_0)
   g2 = func(x2, X_0)

   iter=1
   do iter=1, ITMAX

!..Bisection rule:
   xm = (x1 + x2)/2._dp

!..False position method:
!  xm = (x1*g2 - x2*g1)/(g2 - g1)

   gm=func(xm, X_0)
      if(gm*g1.lt.0._dp) then
         x2 = xm
         g2 = gm
      else
         x1 = xm
         g1 = gm
      endif
      if(abs((x2-x1)/x2) < tol) exit
!     if(abs(gm) < tol) goto 148
   enddo !iter  
   it_done = iter
   zbisect_pot = xm
end function zbisect_pot
!------------------------------------------------------------------------


function zbrent_so(f,x1,x2,ITMAX,tol,iter,occnumtmp)
!------------------------------------------------------------------------
! The brent method (From Numerical recipes f77 code adapted)
! Using the Brent method, find the root of a function func bracketed
! between x1 and x2. The root returned as zbrent, will be refined until
! its accuracy is tol. Parameters: Maximum allowed number of iterations,
! and machine floating-point precision.
!------------------------------------------------------------------------

!..Global
   use params_general
   implicit none
      
!..The function      
   real(dp) :: zbrent_so

!..Arguments
   integer, intent(in) :: ITMAX
   real(dp), intent(in) :: tol,x1,x2
   real(dp), external :: f
   real(dp), parameter :: EPS=zero
   integer, intent(inout) :: iter
   real(dp), intent(inout) :: occnumtmp(lnbasis,3)

!..Local
   real(dp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

   a=x1
   b=x2
   fa=f(a,occnumtmp)
   fb=f(b,occnumtmp)
   if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
      print*,'a,b',a,b
      print*,'fa,fb',fa,fb
      stop 'zbrent_so: root must be bracketed for zbrent_so'
   endif
   c=b
   fc=fb
   do 11 iter=1,ITMAX
      if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
         c=a
         fc=fa
         d=b-a
         e=d
      endif
      if(abs(fc).lt.abs(fb)) then
         a=b
         b=c
         c=a
         fa=fb
         fb=fc
         fc=fa
      endif
      tol1=2._dp*EPS*abs(b)+0.5_dp*tol
      xm=.5_dp*(c-b)
      if(abs(xm).le.tol1 .or. fb.eq.0.)then
         zbrent_so=b
!        print*,'To vrika metaxy:',c,b
!        print*,'solveocc:',fb,fc
         return
      endif
      if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
         s=fb/fa
         if(a.eq.c) then
            p=2._dp*xm*s
            q=1._dp-s
         else
            q=fa/fc
            r=fb/fc
            p=s*(2._dp*xm*q*(q-r)-(b-a)*(r-1._dp))
            q=(q-1._dp)*(r-1._dp)*(s-1._dp)
         endif
         if(p.gt.0._dp) q=-q
         p=abs(p)
         if(2._dp*p .lt. min(3._dp*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
         else
            d=xm
            e=d
         endif
      else
         d=xm
         e=d
      endif
      a=b
      fa=fb
      if(abs(d) .gt. tol1) then
         b=b+d
      else
         b=b+sign(tol1,xm)
      endif
      fb=f(b,occnumtmp)
11 continue
   print *, 'zbrent_so:zbrent exceeding maximum iterations'
   zbrent_so=b
END function zbrent_so
!  (C) Copr. 1986-92 Numerical Recipes Software 5.)2ptN75L:52.
!------------------------------------------------------------------------

function zbrent_gr(f,x1,x2,ITMAX,tol,iter,occnumtmp)

!..Global
   use params_general
   implicit none
  
!..Arguments     
   real(dp) :: zbrent_gr
   integer,intent(in) :: ITMAX
   real(dp), intent(in) :: tol,x1,x2
   real(dp), external :: f
   real(dp), parameter :: EPS=zero
   integer, intent(out) :: iter
   real(dp), intent(inout):: occnumtmp(lnbasis,3)

!..Local
   real(dp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
   a=x1
   b=x2
   fa=f(a,occnumtmp)
   fb=f(b,occnumtmp)
   if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
      print*,'a,b',a,b
      print*,'fa,fb',fa,fb
      stop 'zbrent_gr: root must be bracketed for zbrent'
   endif
   c=b
   fc=fb
   do 11 iter=1,ITMAX
      if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
         c=a
         fc=fa
         d=b-a
         e=d
      endif
      if(abs(fc).lt.abs(fb)) then
         a=b
         b=c
         c=a
         fa=fb
         fb=fc
         fc=fa
      endif
      tol1=2._dp*EPS*abs(b)+0.5_dp*tol
      xm=.5_dp*(c-b)
      if(abs(xm).le.tol1 .or. fb.eq.0.)then
         zbrent_gr=b
         return
      endif
      if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
         s=fb/fa
         if(a.eq.c) then
            p=2._dp*xm*s
            q=1._dp-s
         else
            q=fa/fc
            r=fb/fc
            p=s*(2._dp*xm*q*(q-r)-(b-a)*(r-1._dp))
            q=(q-1._dp)*(r-1._dp)*(s-1._dp)
         endif
         if(p.gt.0._dp) q=-q
         p=abs(p)
         if(2._dp*p .lt. min(3._dp*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
         else
            d=xm
            e=d
         endif
      else
         d=xm
         e=d
      endif
      a=b
      fa=fb
      if(abs(d) .gt. tol1) then
         b=b+d
      else
         b=b+sign(tol1,xm)
      endif
      fb=f(b,occnumtmp)
11 continue
   zbrent_gr=b
END function zbrent_gr

function zbrent_pot(f,x1,x2,ITMAX,tol,iter, &
                           X_0)

!..Global
   use params_general
   implicit none
  
!..Arguments     
   real(dp) :: zbrent_pot
   integer,intent(in) :: ITMAX
   real(dp), intent(in) :: tol,x1,x2
   real(dp), external :: f
   real(dp), parameter :: EPS=zero
   integer, intent(out) :: iter
   real(dp), intent(in):: X_0(*)

!..Local
   real(dp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
   a=x1
   b=x2
   fa=f(a,X_0)
   fb=f(b,X_0)
   if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
      print*,'a,b',a,b
      print*,'fa,fb',fa,fb
      stop 'zbrent_pot: root must be bracketed for zbrent'
   endif
   c=b
   fc=fb
   do 11 iter=1,ITMAX
      if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
         c=a
         fc=fa
         d=b-a
         e=d
      endif
      if(abs(fc).lt.abs(fb)) then
         a=b
         b=c
         c=a
         fa=fb
         fb=fc
         fc=fa
      endif
      tol1=2._dp*EPS*abs(b)+0.5_dp*tol
      xm=.5_dp*(c-b)
      if(abs(xm).le.tol1 .or. fb.eq.0.)then
         zbrent_pot=b
         return
      endif
      if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
         s=fb/fa
         if(a.eq.c) then
            p=2._dp*xm*s
            q=1._dp-s
         else
            q=fa/fc
            r=fb/fc
            p=s*(2._dp*xm*q*(q-r)-(b-a)*(r-1._dp))
            q=(q-1._dp)*(r-1._dp)*(s-1._dp)
         endif
         if(p.gt.0._dp) q=-q
         p=abs(p)
         if(2._dp*p .lt. min(3._dp*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
         else
            d=xm
            e=d
         endif
      else
         d=xm
         e=d
      endif
      a=b
      fa=fb
      if(abs(d) .gt. tol1) then
         b=b+d
      else
         b=b+sign(tol1,xm)
      endif
      fb=f(b,X_0)
11 continue
   zbrent_pot=b
END function zbrent_pot


!------------------------------------------------------------------------
subroutine times_init()
! Just initialize the times spend in each subroutine to zero

!..Global
   use global
   implicit none

   time_tot=0._dp 
   time_fock=0._dp 
   time_occ=0._dp 
   time_vec=0._dp
   time_nbas=0._dp 
   time_tote=0._dp 
   time_rdgm=0._dp 
   time_sumint=0._dp 

end subroutine times_init
!------------------------------------------------------------------------

subroutine print_times()
! Prints the time spent in each (important) subroutine

!..Global
   use params_general
   use global
   implicit none

   write(6,*)'-----------------------------------------------------'
   write(6,*)'TIME DISTRIBUTION INFORMATION:                       '
   write(6,'(a30,f10.2,a4)') '         Total time=', time_tot, &
         'sec'
   write(6,*)'-----------------------------------------------------'
   write(6,'(a30,f10.2,a4,f6.2,a1)') 'Read Gamess out   time=',&
         time_rdgm, 'sec',100._dp*time_rdgm/time_tot,'%'
   write(6,'(a30,f10.2,a4,f6.2,a1)') 'Fock construction time=',&
         time_fock, 'sec',100._dp*time_fock/time_tot,'%'
   write(6,'(a30,f10.2,a4,f6.2,a1)') 'K_ab, J_ab constr time=',&
         time_nbas, 'sec',100._dp*time_nbas/time_tot,'%'
   write(6,'(a30,f10.2,a4,f6.2,a1)') 'Orbital minimiz.  time=',&
         time_vec, 'sec',100._dp*time_vec/time_tot,'%'
   write(6,'(a30,f10.2,a4,f6.2,a1)') 'Occ.Nu. minimiz.  time=',&
         time_occ, 'sec',100._dp*time_occ/time_tot,'%'
   write(6,'(a30,f10.2,a4,f6.2,a1)') 'Integral summ.    time=',&
         time_sumint, 'sec',100._dp*time_sumint/time_tot,'%'
   write(6,'(a30,f10.2,a4,f6.2,a1)') 'Lagr. mat.  calc. time=',&
         time_lagr, 'sec',100._dp*time_lagr/time_tot,'%'
   write(6,'(a30,f10.2,a4,f6.2,a1)') 'Tot. Energy calc. time=',&
         time_tote, 'sec',100._dp*time_tote/time_tot,'%'
   write(6,*)'-----------------------------------------------------'
end subroutine print_times
!------------------------------------------------------------------------

subroutine normal_test(test,criter,iwork)
!------------------------------------------------------------------------
! Checks if vecnat is an orthonormal set of vectors. The basis set 
! is not orthogonal and the overlap matrix is provided by the
! 'matrices.com', while the orbitals are procided through the 
! 'orbocc.com'
! OUTPUT: 
!       test (logical): successful/unsuccessful
!      iwork (integer): 
!                           0: do not print, 
!                       other: print
!------------------------------------------------------------------------

!..Global
   use global
   use matrices
   use orbocc
   implicit none

!..Arguments
   logical :: test
   integer :: iwork
   real(dp) :: criter

!..Local
   integer :: ia,ib,i,j
   real(dp), allocatable :: Fsmall(:,:)
   real(dp) :: trial_mat, comp_0, comp_1, divmax

!..Local array allocation
   allocate ( Fsmall(lnbasis,lnbasis) )

   divmax=0._dp
   test=.true.
   if(iwork /= 0) print*,'-------------------------------------------------'
   if(iwork /= 0) print*,'Orthonormality test:'
   do ia=1,nnatorb
      do ib=1,nbasis
         Fsmall(ia,ib)=0._dp
         do i=1,nbasis
            do j=1,nbasis
!              trial_mat=F(i,j,ia)-F(i,j,ib) ! hermiticity of Lagr mult
               trial_mat=ovlap(i,j)
               Fsmall(ia,ib)= Fsmall(ia,ib) + &
               vecnat(i,ia)*trial_mat*vecnat(j,ib)
            enddo
         enddo
         comp_0 =  abs(Fsmall(ia,ib))
         comp_1 =  abs(Fsmall(ia,ib)-1._dp)
         divmax= max(divmax,min(comp_0, comp_1))
         if(comp_0>criter.and.comp_1>criter) then
            if(iwork /= 0) print*,'for ',ia,ib,'Ovlap=',Fsmall(ia,ib)
            test=.false.
         endif
      enddo
   enddo

   if(test) then
      if(iwork /= 0) print*,'The orthonormality test was PASSED'
   else
      if(iwork /= 0) print*,'The orthonormality test NOT PASSED'
   endif
   if(iwork /= 0) print*,'Maximum difference from zero or one:',divmax
   if(iwork /= 0) print*,'-------------------------------------------------'

   deallocate ( Fsmall )
end subroutine normal_test
!-------------------------------------------------------------------------

subroutine symmet_orthog()
!-----------------------------------------------------------------------
!  Symmetric orthogonalization of vecnat. The basis set 
! is not orthogonal and the overlap matrix is provided by the
! 'matrices.com', while the orbitals are provided through the 
! 'orbocc' module.
! The Formula used is C_orth = ((3/2) I - (1/2) D) C_north
! where C_north is the non orthogonal set of vectors,
!        C_orth is the final orthogonal set,
!        D is the overlap matrix of the C_orth set, and
!        I is the unit matrix.
! It works well for initial vectors not far from being orthogonal 
! since this formula is a taylor expansion for small D
!-----------------------------------------------------------------------

!..Global
   use global
   use orbocc
   use matrices
   implicit none
      
!..local variables:
   integer i,ia,ib
   real(dp), allocatable :: ovl(:,:) 
   real(dp), allocatable ::  vnew(:,:)
   real(dp) :: q_ovl

!..Local array allocation
   allocate( ovl(lnbasis,lnbasis), vnew(lnbasis,lnbasis) )

   do ia=1,nbasis
      do ib=1, ia
         ovl(ia,ib)=dot_product(vecnat(:,ia),(matmul(ovlap,vecnat(:,ib))))
!        do i=1,nbasis
!           do j=1,nbasis
!              ovl(ia,ib)=ovl(ia,ib) &
!                        +vecnat(i,ia)*ovlap(i,j)*vecnat(j,ib)
!           enddo
!        enddo
         ovl(ib,ia)= ovl(ia,ib)
      enddo
   enddo

   do ia=1,nbasis
      do i=1,nbasis
         vnew(i,ia)=0._dp
         do ib=1,nbasis
            q_ovl=-0.5_dp*ovl(ia,ib)
            if(ib == ia) q_ovl=q_ovl+1.5_dp
            vnew(i,ia)=vnew(i,ia)+vecnat(i,ib)*q_ovl
         enddo
      enddo
   enddo

   vecnat=vnew
   deallocate( ovl, vnew )

end subroutine symmet_orthog
!-----------------------------------------------------------------------

function BBC1_strength(ec_over_ek,MC)
!-----------------------------------------------------------------------
! This function calculates the strength of the BBC1 correction i.e.
! the factor multiplying the xc terms for both orbitals being weakly
! occupied as a function of the ratio of the correlation over the 
! kinetic energy 'ec_over_ek'. That is done by spline interpolation
! to the results of the Homogeneous electron gas. For the Homogeneous
! electron gas there are two sets of data. Those fitted to reproduce
! the Ceperley&Adler Monte-Carlo results for the correlation energy 
! and those fitted to reproduce the Ortiz&Balone Monte-Carlo results.
!
! INPUT:
!    ec_over_ek : The ratio E_correlation / E_kinetic
!            MC : 1: Ceperley&Adler, 2: Ortiz&Balone
! OUTPUT (function value): the strength of the xc term
!-----------------------------------------------------------------------

!..Global
   use global
   implicit none

!..The function
   real(dp) :: BBC1_strength
       
!..Arguments
   integer, intent(in) :: MC !1: Ceperley&Adler, 2: Ortiz&Balone
   real(dp),intent(in) :: ec_over_ek

!..Local
   integer, parameter :: SIZE=16
   real(dp) :: str1(0:SIZE), ecek1(0:SIZE)
   real(dp) :: str2(0:SIZE), ecek2(0:SIZE)
   integer :: ntable, i, ierr
   real(dp) :: alpha,str_new


   ntable=SIZE ! The number of str/ecek pairs below

!..Data for strengt/ratio for Ceperley-Alder
   str1(1)  =4.80_dp
   str1(2)  =2.6366_dp
   str1(3)  =1.7585_dp
   str1(4)  =0.981_dp
   str1(5)  =0.6253_dp
   str1(6)  =0.34203_dp
   str1(7)  =-0.0064_dp
   str1(8)  =-0.1176_dp
   str1(9)  =-0.202_dp
   str1(10) =-0.24_dp

   ecek1(1) =-1.0939794e-3_dp
   ecek1(2) =-3.66072678e-3_dp
   ecek1(3) =-7.331623097e-3_dp
   ecek1(4) =-1.73353975e-2_dp
   ecek1(5) =-3.02585e-2_dp
   ecek1(6) =-5.40964154e-2_dp
   ecek1(7) =-0.1620329_dp
   ecek1(8) =-0.3008926129_dp
   ecek1(9) =-0.638405507_dp
   ecek1(10)=-1.680826122_dp

!..Ortiz Ballone
   str2(1)  =4.913_dp !rs=0.1
   str2(2)  =2.7510_dp !rs=0.2
   str2(3)  =1.86703_dp !rs=0.3
   str2(4)  =1.389858_dp !rs=0.4
   str2(5)  =1.08743_dp !rs=0.5
   str2(6)  =0.87688_dp !rs=0.6
   str2(7)  =0.727_dp !rs=0.7
   str2(8)  =0.602166_dp !rs=0.8
   str2(9)  =0.4354_dp !rs=1.0
   str2(10) =0.190003_dp !rs=1.5
   str2(11) =0.0594_dp !rs=2.0
   str2(12) =-0.074_dp !rs=3.0
   str2(13) =-0.1461026_dp !rs=4.0
   str2(14) =-0.189_dp !rs=5.0
   str2(15) =-0.2338115_dp !rs=7.0
   str2(16) =-0.26305_dp !rs=10.0

   ecek2(1) =-1.086055e-3_dp !rs=0.1
   ecek2(2) =-3.616377e-3_dp !rs=0.2
   ecek2(3) =-7.21488e-3_dp !rs=0.3
   ecek2(4) =-1.170173e-2_dp !rs=0.4
   ecek2(5) =-1.6961415e-2_dp !rs=0.5
   ecek2(6) =-2.291107355e-2_dp !rs=0.6
   ecek2(7) =-2.948745e-2_dp !rs=0.7
   ecek2(8) =-3.6640197e-2_dp !rs=0.8
   ecek2(9) =-5.251647e-2_dp !rs=1.0
   ecek2(10) =-0.100056734_dp !rs=1.5
   ecek2(11) = -0.1569155_dp !rs=2.0
   ecek2(12) =-0.292846_dp !rs=3.0
   ecek2(13) =-0.45262455_dp !rs=4.0
   ecek2(14) =-0.631738_dp !rs=5.0
   ecek2(15) =-1.036181313_dp !rs=7.0
   ecek2(16) =-1.728741821_dp !rs=10.0

   do i=1,ntable
     ecek1(i)=-ecek1(i)
     ecek2(i)=-ecek2(i)
   enddo
   alpha=-ec_over_ek

   if(MC==1) then
      alpha=max(min(alpha,ecek1(SIZE)-0.000001_dp),ecek1(1)+0.000001_dp)
      call Akima(ntable, ierr, alpha, str_new, ecek1, str1)
   elseif(MC==2) then
      alpha=max(min(alpha,ecek2(SIZE)-0.000001_dp),ecek2(1)+0.000001_dp)
      call Akima(ntable, ierr, alpha, str_new, ecek2, str2)
   else
      print*,'Function BBC1_strength: MC =', MC
      stop 'BBC1_strength: MC must be 1 or 2 '
   endif
      
   BBC1_strength = str_new

end function BBC1_strength
function BBC1_strength_b(ec_over_ex,MC)
!LIKE BBC1_strength above but the paramere is Ec/Ex
!-----------------------------------------------------------------------
! This function calculates the strength of the BBC1 correction i.e.
! the factor multiplying the xc terms for both orbitals being weakly
! occupied as a function of the ratio of the correlation over the 
! exchange energy 'ec_over_ex'. That is done by spline interpolation
! to the results of the Homogeneous electron gas. For the Homogeneous
! electron gas there are two sets of data. Those fitted to reproduce
! the Ceperley&Adler Monte-Carlo results for the correlation energy 
! and those fitted to reproduce the Ortiz&Balone Monte-Carlo results.
!
! INPUT:
!    ec_over_ek : The ratio E_correlation / E_kinetic
!            MC : 1: Ceperley&Adler, 2: Ortiz&Balone
! OUTPUT (function value): the strength of the xc term
!-----------------------------------------------------------------------

!..Global
   use global
   implicit none

!..The function
   real(dp) :: BBC1_strength_b
       
!..Arguments
   integer, intent(in) :: MC !1: Ceperley&Adler, 2: Ortiz&Balone
   real(dp),intent(in) :: ec_over_ex

!..Local
   integer, parameter :: SIZE=16
   real(dp) :: str1(0:SIZE), ecex1(0:SIZE)
   real(dp) :: str2(0:SIZE), ecex2(0:SIZE)
   integer :: ntable, ierr
   real(dp) :: alpha,str_new


   ntable=SIZE ! The number of str/x pairs below

   if( MC == 1) then
      stop 'auxil:  BBC1_strength_b works only with Ortiz Balone, MC=2'
   endif

!..Ortiz Ballone
   str2(1)  =4.913_dp !rs=0.1
   str2(2)  =2.7510_dp !rs=0.2
   str2(3)  =1.86703_dp !rs=0.3
   str2(4)  =1.389858_dp !rs=0.4
   str2(5)  =1.08743_dp !rs=0.5
   str2(6)  =0.87688_dp !rs=0.6
   str2(7)  =0.727_dp !rs=0.7
   str2(8)  =0.602166_dp !rs=0.8
   str2(9)  =0.4354_dp !rs=1.0
   str2(10) =0.190003_dp !rs=1.5
   str2(11) =0.0594_dp !rs=2.0
   str2(12) =-0.074_dp !rs=3.0
   str2(13) =-0.1461026_dp !rs=4.0
   str2(14) =-0.189_dp !rs=5.0
   str2(15) =-0.2338115_dp !rs=7.0
   str2(16) =-0.26305_dp !rs=10.0

   ecex2(1) =0.26192229e-1_dp !rs=0.1
   ecex2(2) =0.43607816e-1_dp !rs=0.2
   ecex2(3) =0.58000046e-1_dp !rs=0.3
   ecex2(4) =0.70552246e-1_dp !rs=0.4
   ecex2(5) =0.81811200e-1_dp !rs=0.5
   ecex2(6) =0.92090506e-1_dp !rs=0.6
   ecex2(7) =0.10159207_dp  !rs=0.7
   ecex2(8) =0.11045579_dp  !rs=0.8
   ecex2(9) = 0.12665320_dp !rs=1.0
   ecex2(10) = 0.16087025_dp!rs=1.5
   ecex2(11) = 0.18921539_dp!rs=2.0
   ecex2(12) = 0.23541737_dp!rs=3.0
   ecex2(13) = 0.27289701_dp!rs=4.0
   ecex2(14) = 0.30471095_dp!rs=5.0
   ecex2(15) = 0.35699192_dp!rs=7.0
   ecex2(16) = 0.41691815_dp!rs=10.0

!  do i=1,ntable
!    ecex1(i)=-ecex1(i)
!    ecex2(i)=-ecex2(i)
!  enddo
!  alpha=-ec_over_ex
  
   if( ec_over_ex < ecex2(1) .or. ec_over_ex > ecex2(16) ) then
      print*,'!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print*,'!!!!!!!!!a=',ec_over_ex,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print*,'!!!!!!!!!!!!!alpha out of bounds!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   endif


   if(MC==1) then
      alpha=max(min(alpha,ecex1(SIZE)-0.000001_dp),ecex1(1)+0.000001_dp)
      call Akima(ntable, ierr, alpha, str_new, ecex1, str1)
   elseif(MC==2) then
      alpha=max(min(alpha,ecex2(SIZE)-0.000001_dp),ecex2(1)+0.000001_dp)
      call Akima(ntable, ierr, alpha, str_new, ecex2, str2)
   else
      print*,'Function BBC1_strength: MC =', MC
      stop 'BBC1_strength: MC must be 1 or 2 '
   endif
      
   BBC1_strength_b = str_new

end function BBC1_strength_b
!-----------------------------------------------------------------------

!********************************************************
!*          Akima spline fitting subroutine             *
!* ---------------------------------------------------- *
!* The input table is X(i), Y(i), where Y(i) is the     *
!* dependant variable. The interpolation point is x_int,*
!* which is assumed to be in the interval of the table  *
!* with at least one table value to the left, and three *
!* to the right. The interpolated returned value is     *
!* y_int.                                               *
!* ierr is returned as an error check (ierr=0 implies   *
!* error).                                              *
!* It is also assumed that the X(i) are in ascending    *
!* order.                                               *
!********************************************************
subroutine Akima(npoints,ierr,x_int,y_int,X,Y)  
   use global; implicit none
   integer SIZE
   parameter(SIZE=16)
   integer i,npoints,ierr
   real*8 x_int,y_int
   real*8 X(0:SIZE), Y(0:SIZE)
   real*8 XM(0:SIZE+3)
   real*8 Z (0:SIZE)
   real*8 a,b

   ierr=1
!..special case x_int=0
   if (x_int.eq.0.0) then
      y_int=0._dp
      return
   end if
!..Check to see if interpolation point is correct
   if (x_int.lt.X(1).or.x_int.ge.X(npoints)) then
      ierr=0 
      return
   end if
   X(0)=2._dp*X(1)-X(2)
!..Calculate Akima coefficients, a and b
   do i = 1, npoints-1
!.....Shift i to i+2
      XM(i+2)=(Y(i+1)-Y(i))/(X(i+1)-X(i))
   end do
   XM(npoints+2)=2._dp*XM(npoints+1)-XM(npoints)
   XM(npoints+3)=2._dp*XM(npoints+2)-XM(npoints+1)
   XM(2)=2._dp*XM(3)-XM(4)
   XM(1)=2._dp*XM(2)-XM(3)
   do i = 1, npoints
      a=dabs(XM(i+3)-XM(i+2))
      b=dabs(XM(i+1)-XM(i))
      if (a+b.ne.0._dp) goto 10
      Z(i)=(XM(i+2)+XM(i+1))/2._dp
      goto 20
10    Z(i)=(a*XM(i+1)+b*XM(i+2))/(a+b)
20 end do
!..Find relevant table interval
   i=0
30 i=i+1
   if (x_int.gt.X(i)) goto 30
   i=i-1
!..Begin interpolation
   b=X(i+1)-X(i)
   a=x_int-X(i)
   y_int=Y(i)+Z(i)*a+(3._dp*XM(i+2)-2._dp*Z(i)-Z(i+1))*a*a/b
   y_int=y_int+(Z(i)+Z(i+1)-2._dp*XM(i+2))*a*a*a/(b*b)
end subroutine Akima
!-----------------------------------------------------------------------
      
function s_rho(rho)
!-----------------------------------------------------------------------
! INPUT:
!    rho : The local density
! OUTPUT (function value): the strength of the xc term
! The Ortiz Ballone data are used
!-----------------------------------------------------------------------

!..Global
   use global; implicit none

!..The function
   real(dp) :: s_rho
       
!..Arguments
   real(dp),intent(in) :: rho

!..Local
   integer, parameter :: SIZE=16
   real(dp) :: str2(0:SIZE), rstab2(0:SIZE), r_s
   integer :: ntable, ierr
   real(dp) :: str_new


   ntable=SIZE ! The number of str/ecek pairs below

!..Ortiz Ballone
   str2(1)  =4.913_dp;      rstab2(1)=0.1_dp
   str2(2)  =2.7510_dp;     rstab2(2)=0.2_dp
   str2(3)  =1.86703_dp;    rstab2(3)=0.3_dp
   str2(4)  =1.389858_dp;   rstab2(4)=0.4_dp
   str2(5)  =1.08743_dp;    rstab2(5)=0.5_dp
   str2(6)  =0.87688_dp;    rstab2(6)=0.6_dp
   str2(7)  =0.727_dp;      rstab2(7)=0.7_dp
   str2(8)  =0.602166_dp;   rstab2(8)=0.8_dp
   str2(9)  =0.4354_dp;     rstab2(9)=1.0_dp
   str2(10) =0.190003_dp;   rstab2(10)=1.5_dp
   str2(11) =0.0594_dp;     rstab2(11)=2.0_dp
   str2(12) =-0.074_dp;     rstab2(12)=3.0_dp
   str2(13) =-0.1461026_dp; rstab2(13)=4.0_dp
   str2(14) =-0.189_dp;     rstab2(14)=5.0_dp
   str2(15) =-0.2338115_dp; rstab2(15)=7.0_dp
   str2(16) =-0.26305_dp;   rstab2(16)=10.0_dp

!   ecek2(1) =-1.086055d-3 !rs=0.1
!   ecek2(2) =-3.616377d-3 !rs=0.2
!   ecek2(3) =-7.21488d-3 !rs=0.3
!   ecek2(4) =-1.170173d-2 !rs=0.4
!   ecek2(5) =-1.6961415d-2 !rs=0.5
!   ecek2(6) =-2.291107355d-2 !rs=0.6
!   ecek2(7) =-2.948745d-2 !rs=0.7
!   ecek2(8) =-3.6640197d-2 !rs=0.8
!   ecek2(9) =-5.251647d-2 !rs=1.0
!   ecek2(10) =-0.100056734d0 !rs=1.5
!   ecek2(11) = -0.1569155d0 !rs=2.0
!   ecek2(12) =-0.292846d0 !rs=3.0
!   ecek2(13) =-0.45262455d0 !rs=4.0
!   ecek2(14) =-0.631738d0 !rs=5.0
!   ecek2(15) =-1.036181313d0 !rs=7.0
!   ecek2(16) =-1.728741821d0 !rs=10.0

   r_s = (3._dp/(4._dp*PI*rho))**(1._dp/3._dp)

   r_s=max(min(r_s,rstab2(SIZE)-0.000001_dp),rstab2(1)+0.000001_dp)
   call Akima(ntable, ierr, r_s, str_new, rstab2, str2)
      
   s_rho = str_new

end function s_rho
!-----------------------------------------------------------------------

subroutine diagon_lapack(ln,Hmat,Bmat,Emat,Psimat)
!-----------------------------------------------------------------------
! Wrapper for the diagonalizing routine of LAPACK
! Diagonalizer of hermitian matrix Hmat(ln,ln)
!-----------------------------------------------------------------------

!..Global
   use global
   use matrices
   use orbocc
   implicit none

!..Arguments
   integer :: ln
   real(dp) :: Hmat(ln, ln), Bmat(ln, ln),  Emat(ln), Psimat(ln, ln)

!..Local
   real(dp) :: A(ln,ln)
   real(dp) :: B(ln,ln)
   real(dp) :: vectr(ln,ln), eigs(ln)

!..Definitions for LAPACK routine
   character(1) :: JOBZ, RANGE, UPLO
   INTEGER :: IL, INFO, ITYPE,  IU,  LDA,  LDB,  LDZ
   INTEGER LWORK, M, N
   real(dp) :: VL,VU
   real(dp) :: ABSTOL
   real(dp) :: WORK(10*ln)
   integer :: IWORK(5*ln)
   integer :: ifail(ln)

   real(dp) :: DLAMCH
   external DLAMCH
      
   ABSTOL = 2*DLAMCH('S')

   ITYPE=1; JOBZ='V'; RANGE='I'; UPLO='L'; N = ln; LDA = ln
   B = Bmat; LDB = ln; VL = 0._dp; VU = 0._dp; IL = 1; IU = ln
   LDZ = ln; LWORK=10*ln
      
   A = Hmat
   call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
                LDB,  VL,  VU, IL, IU, ABSTOL, M, eigs, vectr, &
                LDZ, WORK, LWORK, IWORK, IFAIL, INFO )
   if ( INFO .ne. 0 ) then 
      print*, 'ERROR!', INFO
      STOP
   endif
!..Now store eigs --> ennat, vectr --> vecnat
   Emat=eigs; Psimat=vectr

end subroutine diagon_lapack

!--------------------------------------------------------------
subroutine diagon_NAG(ln,Hmat,Bmat,Emat,Psimat)
! Wrapper for the diagonalizing routine of NAG
!-----------------------------------------------------------------------

!..Global
   use global
   use matrices
   use orbocc
   implicit none

!..Arguments
   integer :: ln
   real(dp) :: Hmat(ln,ln), Bmat(ln,ln), Emat(ln), Psimat(ln,ln)

!..For NAG routine
   integer :: ITYPE=1, LWORK, IFAIL 
   character(1) :: JOB='V', UPLO='U'
   real(dp), allocatable :: A(:,:), B(:,:), W(:), WORK(:)

   LWORK = 64*ln

   allocate( A(ln,ln), B(ln,ln), W(ln), WORK(LWORK) )

   IFAIL=0

   A=Hmat; B=Bmat; 

   call F02FDF( ITYPE, JOB, UPLO, ln, A, ln, B, ln, W, &
                WORK, LWORK, IFAIL)

   if(ifail/=0) then
      print*,'ifail',ifail
      stop 'auxil:diagon_NAG: ifail not zero'
   endif

!..Now store W --> ennat, A --> vecnat
   Emat=W; Psimat=A   

   deallocate ( A, B, W, WORK )

end subroutine diagon_NAG

!-----------------------------------------------------------------------

subroutine trace_bond_anti(occ)
! Traces the bonding and anti-bonding for BBC3
   use global
   use functional_m
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)

!..Local variables:
   integer :: ia
   real(dp) :: occbond1, occanti1, occbond2, occanti2

   ibond(1)=1; ianti(1)=1
   ibond(2)=1; ianti(2)=1
   occbond1=2; occanti1=2
   occbond2=2; occanti2=2
   do ia=1,nnatorb
      if(occ(ia,1) > 0.5_dp) then
         if(occ(ia,1)-0.5_dp <= occbond1) then
            ibond(1) = ia
            occbond1 = occ(ia,1)-0.5_dp 
         endif
      endif
      if(occ(ia,1) < 0.5_dp ) then
         if(0.5_dp-occ(ia,1) < occanti1) then
            ianti(1) = ia
            occanti1 = 0.5_dp-occ(ia,1)
         endif
      endif
      if(.not.cl_shell) then
         if(occ(ia,2) > 0.5_dp) then
            if(occ(ia,2)-0.5_dp <= occbond2) then
               ibond(2) = ia
               occbond2 = occ(ia,2)-0.5_dp
            endif
         endif
         if(occ(ia,2) < 0.5_dp) then
            if(0.5_dp-occ(ia,2) < occanti2) then
               ianti(2) = ia
               occanti2 = 0.5_dp-occ(ia,2)
            endif
         endif
      endif
   enddo
   if(cl_shell) then
      ibond(2)=ibond(1)
      ianti(2)=ianti(1)
   endif

!Not trace at all keep it fixed!
   ibond(1)=nele(1); ibond(2)=nele(2)
   ianti(1)=nele(1)+1; ianti(2)=nele(2)+1
   
!  print*,'TRACE BONDING ANTIBONDING FOR BBC3:'
!  print*,'ibond(1)=',ibond(1),'ianti(1)=',ianti(1)
!  print*,'ibond(2)=',ibond(2),'ianti(2)=',ianti(2)
!  print*,nele
!  print*,'ibond(2)=',ibond(2)

!  print*,'ianti(2)=',ianti(2)
!  stop 'BBC3 stop'
!   ibond(1)=2; ibond(2)=2; ianti(1)=3; ianti(2)=3
!   ibond(1)=7; ibond(2)=7; ianti(1)=8; ianti(2)=8
   
end subroutine trace_bond_anti
   
!---------------------------------------------------------------------
subroutine BB_orb_kinds()
   use global; use functional_m; use orbocc
   implicit none

!..Local variables
   integer :: ia, ispin, nsp, ideg
   real(dp) :: del
   real(dp), parameter :: deg_crit=1.e-5_dp

!  ideg=0 ! No degeneracies: kind_BBC3 = 1 1...1 1 2 3 4 4...4
   ideg=1 ! respect degeneracies of bonding-antibonding

   do ispin=1,2
      nsp = xnele(ispin)+small; if(xnele(ispin) - nsp >= 0.009_dp) nsp=nsp+1
      ibond(ispin) = nsp
      if (.not. (Functional == 'RHF'.or. Functional == 'RHA'.or.Functional == 'DFT')) then
         if(nnatorb-nsp < 1 ) then
            stop 'BB_orb_kinds: too few natural orbitals'
         endif
      endif
      do ia=1,nnatorb
         if(ia < nsp) then
            kind_BBC3(ia,ispin) = 1 !Strongly occupied
         else if (ia == nsp) then
            kind_BBC3(ia,ispin) = 2 !Bonding
         else if (ia == nsp+1 ) then
            kind_BBC3(ia,ispin) = 3 !Antibonding
         else
            kind_BBC3(ia,ispin) = 4 !Weakly occupied
         endif
      enddo

!.....Degenerate bonding
      if(ideg == 1) then
         do ia=nsp,1,-1
            del=abs(ennat(ia)-ennat(nsp))
            if(del < deg_crit )  kind_BBC3(ia,ispin) = 2
         enddo

!.....Degenerate anti-bonding
         do ia=nsp+1,nnatorb
            del=abs(ennat(ia)-ennat(nsp+1))
            if(del < deg_crit ) kind_BBC3(ia,ispin) = 3
         enddo
      endif !(ideg == 1)
   enddo
 
!   ibond(1) = 2
!   ibond(2) = 2

!  kind_BBC3 = 4
!  kind_BBC3(1,1) = 1
!  kind_BBC3(2,1) = 1
!  kind_BBC3(3,1) = 2
!  kind_BBC3(4,1) = 2
!  kind_BBC3(5,1) = 3
!  kind_BBC3(6,1) = 3

!  kind_BBC3(1,2) = 1
!  kind_BBC3(2,2) = 1
!  kind_BBC3(3,2) = 2
!  kind_BBC3(4,2) = 2
!  kind_BBC3(5,2) = 3
!  kind_BBC3(6,2) = 3

!..Original BBC3 functional: TO BE COMMENTED OUT!
!  UNCOMMENT ONLY FOR TESTING PURPOSES
!  do ispin = 1,2
!     nsp = xnele(ispin)+small; if(xnele(ispin) - nsp > 0.5d0) nsp=nsp+1
!     do ia=1,nsp-1
!        kind_BBC3(ia,ispin) = 1
!     enddo
!     kind_BBC3(nsp,ispin) = 2
!     kind_BBC3(nsp+1,ispin) =3
!     do ia=nsp+2,nnatorb
!        kind_BBC3(ia,ispin) = 4
!     enddo
!  enddo

end subroutine BB_orb_kinds


!=======================================================================
subroutine  invert_S()
!..Global
   use global; use matrices; use orbocc; implicit none

!..Local
   real(dp), allocatable :: A(:,:)
   real(dp), allocatable :: B(:,:)
   real(dp), allocatable :: vectr(:,:), eigs(:)
   real(dp), allocatable :: Zmat(:)
   integer :: lndim, ndim, mu, nu, ku, info
!..Definitions for LAPACK routine
   character(1) :: JOBZ, RANGE, UPLO
   INTEGER :: IL, ITYPE,  IU,  LDA,  LDB,  LDZ
   INTEGER LWORK, M, N
   real(dp) :: VL,VU
   real(dp) :: ABSTOL
   real(dp), allocatable :: WORK(:)
   integer, allocatable :: IWORK(:)
   integer, allocatable :: ifail(:)
   
   real(dp) :: DLAMCH
   external DLAMCH

!..Local array allocation
   ndim=nbasis
   lndim=ndim
   allocate( A(lndim,lndim), B(lndim,lndim), Zmat(ndim), &
               vectr(lndim,lndim), eigs(lndim), &
               WORK(10*lndim), IWORK(5*lndim), ifail(lndim) )

   ITYPE=1; JOBZ='V'; RANGE='I'; UPLO='L'; N = ndim; LDA = lndim
   LDB = lndim; VL = 0._dp; VU = 0._dp; IL = 1; IU = ndim
   LDZ = lndim; LWORK=10*ndim

   ABSTOL = 2._dp*DLAMCH('S')
!  ABSTOL = 1.d-30
      
   do mu=1, ndim   
      do nu=1, mu
         A(mu,nu) = ovlap(mu,nu)
         B(mu,nu) = 0._dp
         B(nu,mu) = 0._dp
         A(nu,mu) = A(mu,nu) 
      enddo
      B(mu,mu)=1._dp
   enddo

   WORK=0._dp
   IWORK=0
   call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
                LDB,  VL,  VU, IL, IU, ABSTOL, M, eigs, vectr, &
                LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

   if (info /= 0) then
      print*,'Info:',info
      if(info < 0) print*,info, '-th argument of DSYGVX has wrong value'
      if(info > 0 .and. info <= lndim ) print*,info, '-th eigenvalue not converged'
      if(info > lndim .and. info <= 2*lndim ) &
         print*,info-lndim, '-th leading minor of ovlap not positive definite'
      stop 'diagon_Ssq:Diagonalization failed!'
   endif

   do ku=1,nbasis
      eigs(ku)=1._dp/max(eigs(ku),zero)
   enddo

   do mu=1,nbasis
      do nu=1,mu
         ov_inv(mu,nu)=0._dp
         do ku=1,nbasis
            ov_inv(mu,nu) = ov_inv(mu,nu) + vectr(mu,ku)*eigs(ku)*vectr(nu,ku)
         enddo
         ov_inv(nu,mu)= ov_inv(mu,nu)
      enddo
   enddo
end subroutine  invert_S
!=====================================================================================
subroutine Calc_Scr_dens(Nst, x0, y0, z0, x1, y1, z1)

   use global; use grid_params; use DFT_par; implicit none
   real(dp) :: x0, y0, z0, x1, y1, z1

   real(dp) :: xstep, ystep, zstep, r_d, r_screening, rho_xyz, Vxc_Laplace
   real(dp) :: Vxc_xyz, Q_scr_pos, Q_scr_neg, x, y, z
   integer :: Nst, ig, ist

   xstep=(x1-x0)/real(Nst-1); ystep=(y1-y0)/real(Nst-1); zstep=(z1-z0)/real(Nst-1)
   x=x0; y=y0; z=z0
   open(unit=51, file='Scr_dens',status='unknown')
   do ist = 1, Nst
      r_d = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      call Calc_Vxc_Laplace(x, y, z, rho_xyz, Vxc_xyz, Vxc_Laplace)
      r_screening = rho_xyz - one_ov_4pi*Vxc_Laplace
      write(51,*) r_d, r_screening, rho_xyz

! Reconstruct Vxc and Vha from rho_xc and rho
!     Vreco=0.d0; Vha=0.d0
!     do ig=1,ngrid
!        rrr=sqrt((x-x_grid(ig))**2+(y-y_grid(ig))**2+(z-z_grid(ig))**2)
!        if(rrr>1.d-6 ) then
!           rrr=1.d0/rrr
!        else
!           rrr=0.d0
!        endif
!        call Calc_Vxc_Laplace(x_grid(ig), y_grid(ig), z_grid(ig), rho_xyz, &
!                              Vdumm, Vxc_Laplace)
!        rrr=rrr* w_grid(ig)
!        Vreco= Vreco + (- one_ov_4pi*Vxc_Laplace)*rrr
!        Vha = Vha + rho_xyz*rrr
!     enddo
!     write(71,'(5e20.10)') r_d, Vha, Vreco, Vxc_xyz, Vha + Vreco
!     print*,'ist',ist

      write(71,'(5e20.10)') r_d, Vxc_xyz

 
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo
   close(51)

!..Find matrix elements of reconstructed Vxc(r) for comparison
!  allocate (Vtest(lnbasis,lnbasis), Vtest_ig(ngrid))
!  do ig=1,ngrid
!     call Calc_Vxc_Laplace(x_grid(ig), y_grid(ig), z_grid(ig), rho_xyz, Vxc_xyz, Vxc_Laplace)
!     Vtest_ig(ig)=Vxc_xyz
!  enddo
!  do l=1,nbasis
!     do k=1,l
!        Vtest(l,k)=0.d0
!        do ig=1,ngrid
!           Vtest(l,k)=Vtest(l,k)+w_grid(ig)*bas_f_grid(ig,l)*Vtest_ig(ig)*bas_f_grid(ig,k)
!        enddo
!        Vtest(k,l)=Vtest(l,k)
!     enddo
!  enddo

!  open(unit=121,file='Vtest', status='unknown')
!     write(121,'(2i4,f20.10)') ((l,k,Vtest(k,l),k=1,l),l=1,nbasis)
!  close(121)

   Q_scr_pos=0._dp; Q_scr_neg=0._dp
   do ig=1,ngrid
      call Calc_Vxc_Laplace(x_grid(ig),y_grid(ig),z_grid(ig),rho_xyz,Vxc_xyz,Vxc_Laplace)
      r_screening=rho_xyz-one_ov_4pi*Vxc_Laplace
      if(r_screening > 0._dp) then
         Q_scr_pos=Q_scr_pos+w_grid(ig)*r_screening
      else
         Q_scr_neg=Q_scr_neg+w_grid(ig)*r_screening
      endif
   enddo
   write(6,*)'----------------------------------------'
   write(6,'("Total Possitive Screen. Charge:",F20.12,"Total Negative Screen. Charge:",F20.12)')Q_scr_pos, Q_scr_neg
   write(6,'("Total Screen. Charge:",F20.12)')Q_scr_pos+Q_scr_neg
   write(6,*)'----------------------------------------'

end subroutine Calc_Scr_dens
!=====================================================================================

subroutine Calc_Vxc_Laplace(x,y,z,rho_xyz,Vxc_xyz, Vxc_Laplace)
   use global; use grid_params; use DFT_par; implicit none

   real(dp) :: Vxc_Laplace, x,y,z, rho_xyz, Vxc_xyz

   real(dp) :: x_temp(7), y_temp(7), z_temp(7), rho_temp(7), Vxc_temp(7), ddd=1d-4, ddd1
   real(dp) :: Exc_DFT_temp(7), Vc_temp(7), Ec_DFT_temp(7), G_temp(7), Vxc_sig_temp(7), Vc_sig_temp(7)
   real(dp) :: dvdx2, dvdy2, dvdz2
   real(dp) :: x_temp_p(7), x_temp_m(7), rho_temp_m(7)
   real(dp) :: y_temp_p(7), y_temp_m(7)
   real(dp) :: z_temp_p(7), z_temp_m(7), G_temp_p(7), G_temp_m(7)
   real(dp) ::G_temp_c(7,3), G_temp_c_p(7,3), G_temp_c_m(7,3), rho_temp_p(7) 
   real(dp) :: Vxc_temp_p(7), Vxc_sig_temp_p(7), Vc_temp_p(7), Vc_sig_temp_p(7)
   real(dp) :: Vxc_temp_m(7), Vxc_sig_temp_m(7), Vc_temp_m(7), Vc_sig_temp_m(7)
   real(dp),external :: electron_density
   integer :: i


   x_temp(1)=x; y_temp(1)=y; z_temp(1)=z
   x_temp(2)=x+ddd; y_temp(2)=y; z_temp(2)=z
   x_temp(3)=x-ddd; y_temp(3)=y; z_temp(3)=z
   x_temp(4)=x; y_temp(4)=y+ddd; z_temp(4)=z
   x_temp(5)=x; y_temp(5)=y-ddd; z_temp(5)=z
   x_temp(6)=x; y_temp(6)=y; z_temp(6)=z+ddd
   x_temp(7)=x; y_temp(7)=y; z_temp(7)=z-ddd

   do i=1,7
      rho_temp(i) = electron_density(x_temp(i), y_temp(i), z_temp(i),3)
   enddo
   rho_xyz=rho_temp(1)

!..Exchange Correlation Functional
   select case (xc_f90_info_family(xc_info))
   case(XC_FAMILY_LDA)
      call xc_f90_lda_exc_vxc(xc_func, 7, rho_temp(1), Exc_DFT_temp(1), Vxc_temp(1))
   case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call density_gradient_select(7, x_temp, y_temp, z_temp, G_temp, G_temp_c)
      call xc_f90_gga_exc_vxc(xc_func, 7, rho_temp(1), G_temp(1), Exc_DFT_temp(1), Vxc_temp(1), Vxc_sig_temp(1) )

      ddd1=0.5_dp*ddd
      x_temp_p = x_temp + ddd1
      x_temp_m = x_temp - ddd1
      y_temp_p = y_temp + ddd1
      y_temp_m = y_temp - ddd1
      z_temp_p = z_temp + ddd1
      z_temp_m = z_temp - ddd1

      do i=1,7
         rho_temp_p(i) = electron_density(x_temp_p(i), y_temp(i), z_temp(i),3)
         rho_temp_m(i) = electron_density(x_temp_m(i), y_temp(i), z_temp(i),3)
      enddo
      call density_gradient_select(7, x_temp_p, y_temp, z_temp, G_temp_p, G_temp_c_p)
      call density_gradient_select(7, x_temp_m, y_temp, z_temp, G_temp_m, G_temp_c_m)
      call xc_f90_gga_exc_vxc(xc_func, 7, rho_temp_p(1), G_temp_p(1), Exc_DFT_temp(1), Vxc_temp_p(1), Vxc_sig_temp_p(1) )
      call xc_f90_gga_exc_vxc(xc_func, 7, rho_temp_m(1), G_temp_m(1), Exc_DFT_temp(1), Vxc_temp_m(1), Vxc_sig_temp_m(1) )
      do i=1,7
         Vxc_temp(i) = Vxc_temp(i)-2._dp*(Vxc_sig_temp_p(i)* G_temp_c_p(i,1) - Vxc_sig_temp_m(i)* G_temp_c_m(i,1))/ddd
      enddo

      Vxc_temp_p=0._dp; Vxc_sig_temp_p=0._dp; Vxc_temp_m=0._dp; Vxc_sig_temp_m=0._dp
      do i=1,7
         rho_temp_p(i) = electron_density(x_temp(i), y_temp_p(i), z_temp(i),3)
         rho_temp_m(i) = electron_density(x_temp(i), y_temp_m(i), z_temp(i),3)
      enddo
      call density_gradient_select(7, x_temp, y_temp_p, z_temp, G_temp_p, G_temp_c_p)
      call density_gradient_select(7, x_temp, y_temp_m, z_temp, G_temp_m, G_temp_c_m)
      call xc_f90_gga_exc_vxc(xc_func, 7, rho_temp_p(1), G_temp_p(1), Exc_DFT_temp(1), Vxc_temp_p(1), Vxc_sig_temp_p(1) )
      call xc_f90_gga_exc_vxc(xc_func, 7, rho_temp_m(1), G_temp_m(1), Exc_DFT_temp(1), Vxc_temp_m(1), Vxc_sig_temp_m(1) )
      do i=1,7
         Vxc_temp(i) = Vxc_temp(i)-2._dp*(Vxc_sig_temp_p(i)* G_temp_c_p(i,2) - Vxc_sig_temp_m(i)* G_temp_c_m(i,2))/ddd
      enddo

      Vxc_temp_p=0._dp; Vxc_sig_temp_p=0._dp; Vxc_temp_m=0._dp; Vxc_sig_temp_m=0._dp
      do i=1,7
         rho_temp_p(i) = electron_density(x_temp(i), y_temp(i), z_temp_p(i),3)
         rho_temp_m(i) = electron_density(x_temp(i), y_temp(i), z_temp_m(i),3)
      enddo
      call density_gradient_select(7, x_temp, y_temp, z_temp_p, G_temp_p, G_temp_c_p)
      call density_gradient_select(7, x_temp, y_temp, z_temp_m, G_temp_m, G_temp_c_m)
      call xc_f90_gga_exc_vxc(xc_func, 7, rho_temp_p(1), G_temp_p(1), Exc_DFT_temp(1), Vxc_temp_p(1), Vxc_sig_temp_p(1) )
      call xc_f90_gga_exc_vxc(xc_func, 7, rho_temp_m(1), G_temp_m(1), Exc_DFT_temp(1), Vxc_temp_m(1), Vxc_sig_temp_m(1) )
      do i=1,7
         Vxc_temp(i) = Vxc_temp(i)-2._dp*(Vxc_sig_temp_p(i)*G_temp_c_p(i,3) - Vxc_sig_temp_m(i)* G_temp_c_m(i,3))/ddd
      enddo

      gga_x=.true.
   case default
      Exc_DFT_temp=0._dp; Vxc_temp=0._dp
   end select

!..Correlation Functional for seperate additional correlation functional
   select case (xc_f90_info_family(c_info))
   case(XC_FAMILY_LDA)
      call xc_f90_lda_exc_vxc(c_func, 7, rho_temp(1), Ec_DFT_temp(1), Vc_temp(1))
   case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA) 
      if (.not. gga_x) call density_gradient_select(7, x_temp, y_temp, z_temp, G_temp, G_temp_c )
      call xc_f90_gga_exc_vxc(c_func, 7, rho_temp(1), G_temp(1), Ec_DFT_temp(1), Vc_temp(1), Vc_sig_temp(1) )

      x_temp_p = x_temp + ddd1
      x_temp_m = x_temp - ddd1
      y_temp_p = y_temp + ddd1
      y_temp_m = y_temp - ddd1
      z_temp_p = z_temp + ddd1
      z_temp_m = z_temp - ddd1
      do i=1,7
         rho_temp_p(i) = electron_density(x_temp_p(i), y_temp(i), z_temp(i),3)
         rho_temp_m(i) = electron_density(x_temp_m(i), y_temp(i), z_temp(i),3)
      enddo
      call density_gradient_select(7, x_temp_p, y_temp, z_temp, G_temp_p, G_temp_c_p)
      call density_gradient_select(7, x_temp_m, y_temp, z_temp, G_temp_m, G_temp_c_m)
      call xc_f90_gga_exc_vxc(c_func, 7, rho_temp_p(1), G_temp_p(1), Ec_DFT_temp(1), Vc_temp_p(1), Vc_sig_temp_p(1) )
      call xc_f90_gga_exc_vxc(c_func, 7, rho_temp_m(1), G_temp_m(1), Ec_DFT_temp(1), Vc_temp_m(1), Vc_sig_temp_m(1) )
      do i=1,7
         Vc_temp(i) = Vc_temp(i)-2._dp*(Vc_sig_temp_p(i)* G_temp_c_p(i,1) - Vc_sig_temp_m(i)* G_temp_c_m(i,1))/ddd
      enddo

      Vc_temp_p=0._dp; Vc_sig_temp_p=0._dp; Vc_temp_m=0._dp; Vc_sig_temp_m=0._dp
      do i=1,7
         rho_temp_p(i) = electron_density(x_temp(i), y_temp_p(i), z_temp(i),3)
         rho_temp_m(i) = electron_density(x_temp(i), y_temp_m(i), z_temp(i),3)
      enddo
      call density_gradient_select(7, x_temp, y_temp_p, z_temp, G_temp_p, G_temp_c_p)
      call density_gradient_select(7, x_temp, y_temp_m, z_temp, G_temp_m, G_temp_c_m)
      call xc_f90_gga_exc_vxc(c_func, 7, rho_temp_p(1), G_temp_p(1), Ec_DFT_temp(1), Vc_temp_p(1), Vc_sig_temp_p(1) )
      call xc_f90_gga_exc_vxc(c_func, 7, rho_temp_m(1), G_temp_m(1), Ec_DFT_temp(1), Vc_temp_m(1), Vc_sig_temp_m(1) )
      do i=1,7
         Vc_temp(i) = Vc_temp(i)-2._dp*(Vc_sig_temp_p(i)* G_temp_c_p(i,2) - Vc_sig_temp_m(i)* G_temp_c_m(i,2))/ddd
      enddo

      Vc_temp_p=0._dp; Vc_sig_temp_p=0._dp; Vc_temp_m=0._dp; Vc_sig_temp_m=0._dp
      do i=1,7
         rho_temp_p(i) = electron_density(x_temp(i), y_temp(i), z_temp_p(i),3)
         rho_temp_m(i) = electron_density(x_temp(i), y_temp(i), z_temp_m(i),3)
      enddo
      call density_gradient_select(7, x_temp, y_temp, z_temp_p, G_temp_p, G_temp_c_p)
      call density_gradient_select(7, x_temp, y_temp, z_temp_m, G_temp_m, G_temp_c_m)
      call xc_f90_gga_exc_vxc(c_func, 7, rho_temp_p(1), G_temp_p(1), Ec_DFT_temp(1), Vc_temp_p(1), Vc_sig_temp_p(1) )
      call xc_f90_gga_exc_vxc(c_func, 7, rho_temp_m(1), G_temp_m(1), Ec_DFT_temp(1), Vc_temp_m(1), Vc_sig_temp_m(1) )
      do i=1,7
         Vc_temp(i) = Vc_temp(i) -2._dp*(Vc_sig_temp_p(i)* G_temp_c_p(i,3) - Vc_sig_temp_m(i)* G_temp_c_m(i,3))/ddd
      enddo

      gga_c=.true.
   case default
      Ec_DFT_temp(1)=0._dp; Vc_temp=0._dp
   end select

   Vxc_temp=Vxc_temp +Vc_temp
   Vxc_xyz=Vxc_temp(1) 
   
   dvdx2=(Vxc_temp(3)+Vxc_temp(2)-2._dp*Vxc_temp(1))/(ddd*ddd)
   dvdy2=(Vxc_temp(5)+Vxc_temp(4)-2._dp*Vxc_temp(1))/(ddd*ddd)
   dvdz2=(Vxc_temp(7)+Vxc_temp(6)-2._dp*Vxc_temp(1))/(ddd*ddd)

   Vxc_Laplace=dvdx2+dvdy2+dvdz2

end subroutine Calc_Vxc_Laplace
!=====================================================================================
subroutine Plot_Vxc_DFT_r(Nst, x0, y0, z0, x1, y1, z1)
!..Global
   use global; use grid_params; use DFT_par; implicit none

!..Arguments
   real(dp) :: x0, y0, z0, x1, y1, z1
   integer :: Nst

!..Local
   real(dp) :: x(Nst), y(Nst), z(Nst)
   real(dp) :: xstep, ystep, zstep, r_d(Nst), rho_xyz(Nst), ddd=1e-4_dp, ddd1
   real(dp) :: Vxc_xyz(Nst), Exc_DFT_xyz(Nst), G_xyz(Nst), G_xyz_c(Nst,3)
   real(dp) :: Vxc_sig_xyz(Nst), Vc_xyz(Nst), Ec_DFT_xyz(Nst), Vc_sig_xyz(Nst)

   real(dp) :: rho_p(Nst), rho_m(Nst),  G_p(Nst), G_m(Nst)
   real(dp) :: G_c_p(Nst,3), G_c_m(Nst,3)
   real(dp) :: x_p(Nst), x_m(Nst)
   real(dp) :: y_p(Nst), y_m(Nst)
   real(dp) :: z_p(Nst), z_m(Nst) 
   real(dp) :: Vxc_p(Nst), Vxc_sig_p(Nst), Vc_p(Nst), Vc_sig_p(Nst)
   real(dp) :: Vxc_m(Nst), Vxc_sig_m(Nst), Vc_m(Nst), Vc_sig_m(Nst)
   real(dp),external :: electron_density
   integer :: i, ist

   if ( xc_f90_info_family(xc_info) == XC_FAMILY_HYB_GGA) then
      print*,'Plot_Vxc_DFT_r: Potential for Hybrid functionals are not local'
      print*,'Plot_Vxc_DFT_r: Exits...'
      return
   endif

   xstep=(x1-x0)/real(Nst-1); ystep=(y1-y0)/real(Nst-1); zstep=(z1-z0)/real(Nst-1)

   do ist = 1, Nst   
     x(ist)=x0+real(ist-1)*xstep
     y(ist)=y0+real(ist-1)*ystep
     z(ist)=z0+real(ist-1)*zstep
     r_d(ist) = sqrt((x(ist)-x0)**2+(y(ist)-y0)**2+(z(ist)-z0)**2)
     rho_xyz(ist) = electron_density(x(ist), y(ist), z(ist), 3)
   enddo
   if(xc_f90_info_family(xc_info) == XC_FAMILY_GGA .or. xc_f90_info_family(c_info) == XC_FAMILY_GGA) then
      call density_gradient_select(Nst, x, y, z, G_xyz, G_xyz_c)
   endif

!..Exchange Correlation Functional
   select case (xc_f90_info_family(xc_info))
   case(XC_FAMILY_LDA)
      call xc_f90_lda_exc_vxc(xc_func, Nst, rho_xyz(1), Exc_DFT_xyz(1), Vxc_xyz(1))
   case(XC_FAMILY_GGA)
      call xc_f90_gga_exc_vxc(xc_func, Nst, rho_xyz(1), G_xyz(1), Exc_DFT_xyz(1), Vxc_xyz(1), Vxc_sig_xyz(1) )
      
      ddd1=0.5_dp*ddd
      x_p = x + ddd1
      x_m = x - ddd1
      y_p = y + ddd1
      y_m = y - ddd1
      z_p = z + ddd1
      z_m = z - ddd1

      do i=1,nst
         rho_p(i) = electron_density(x_p(i), y(i), z(i),3)
         rho_m(i) = electron_density(x_m(i), y(i), z(i),3)
      enddo
      call density_gradient_select(Nst, x_p, y, z, G_p, G_c_p)
      call density_gradient_select(Nst, x_m, y, z, G_m, G_c_m)
      call xc_f90_gga_exc_vxc(xc_func, Nst, rho_p(1), G_p(1), Exc_DFT_xyz(1), Vxc_p(1), Vxc_sig_p(1) )
      call xc_f90_gga_exc_vxc(xc_func, Nst, rho_m(1), G_m(1), Exc_DFT_xyz(1), Vxc_m(1), Vxc_sig_m(1) )
      do i=1,Nst
         Vxc_xyz(i) = Vxc_xyz(i)-2._dp*(Vxc_sig_p(i)* G_c_p(i,1) - Vxc_sig_m(i)* G_c_m(i,1))/ddd
      enddo

      Vxc_p=0._dp; Vxc_sig_p=0._dp; Vxc_m=0._dp; Vxc_sig_m=0._dp
      do i=1,Nst
         rho_p(i) = electron_density(x(i), y_p(i), z(i),3)
         rho_m(i) = electron_density(x(i), y_m(i), z(i),3)
      enddo
      call density_gradient_select(Nst, x, y_p, z, G_p, G_c_p)
      call density_gradient_select(Nst, x, y_m, z, G_m, G_c_m)
      call xc_f90_gga_exc_vxc(xc_func, Nst, rho_p(1), G_p(1), Exc_DFT_xyz(1), Vxc_p(1), Vxc_sig_p(1) )
      call xc_f90_gga_exc_vxc(xc_func, Nst, rho_m(1), G_m(1), Exc_DFT_xyz(1), Vxc_m(1), Vxc_sig_m(1) )
      do i=1,Nst
         Vxc_xyz(i) = Vxc_xyz(i)-2._dp*(Vxc_sig_p(i)* G_c_p(i,2) - Vxc_sig_m(i)* G_c_m(i,2))/ddd
      enddo

      Vxc_p=0._dp; Vxc_sig_p=0._dp; Vxc_m=0._dp; Vxc_sig_m=0._dp
      do i=1,Nst
         rho_p(i) = electron_density(x(i), y(i), z_p(i),3)
         rho_m(i) = electron_density(x(i), y(i), z_m(i),3)
      enddo
      call density_gradient_select(Nst, x, y, z_p, G_p, G_c_p)
      call density_gradient_select(Nst, x, y, z_m, G_m, G_c_m)
      call xc_f90_gga_exc_vxc(xc_func, Nst, rho_p(1), G_p(1), Exc_DFT_xyz(1), Vxc_p(1), Vxc_sig_p(1) )
      call xc_f90_gga_exc_vxc(xc_func, Nst, rho_m(1), G_m(1), Exc_DFT_xyz(1), Vxc_m(1), Vxc_sig_m(1) )
      do i=1,Nst
         Vxc_xyz(i) = Vxc_xyz(i)-2._dp*(Vxc_sig_p(i)*G_c_p(i,3) - Vxc_sig_m(i)* G_c_m(i,3))/ddd
      enddo

      gga_x=.true.
   case default
   end select

!..Correlation Functional for seperate additional correlation functional
   select case (xc_f90_info_family(c_info))
   case(XC_FAMILY_LDA)
      call xc_f90_lda_exc_vxc(c_func, Nst, rho_xyz(1), Ec_DFT_xyz(1), Vc_xyz(1))
   case(XC_FAMILY_GGA) 
      call xc_f90_gga_exc_vxc(c_func, Nst, rho_xyz(1), G_xyz(1), Ec_DFT_xyz(1), Vc_xyz(1), Vc_sig_xyz(1) )

      x_p = x + ddd1
      x_m = x - ddd1
      y_p = y + ddd1
      y_m = y - ddd1
      z_p = z + ddd1
      z_m = z - ddd1
      do i=1,Nst
         rho_p(i) = electron_density(x_p(i), y(i), z(i),3)
         rho_m(i) = electron_density(x_m(i), y(i), z(i),3)
      enddo
      call density_gradient_select(Nst, x_p, y, z, G_p, G_c_p)
      call density_gradient_select(Nst, x_m, y, z, G_m, G_c_m)
      call xc_f90_gga_exc_vxc(c_func, Nst, rho_p(1), G_p(1), Ec_DFT_xyz(1), Vc_p(1), Vc_sig_p(1) )
      call xc_f90_gga_exc_vxc(c_func, Nst, rho_m(1), G_m(1), Ec_DFT_xyz(1), Vc_m(1), Vc_sig_m(1) )
      do i=1,Nst
         Vc_xyz(i) = Vc_xyz(i)-2._dp*(Vc_sig_p(i)* G_c_p(i,1) - Vc_sig_m(i)* G_c_m(i,1))/ddd
      enddo

      Vc_p=0._dp; Vc_sig_p=0._dp; Vc_m=0._dp; Vc_sig_m=0._dp
      do i=1,Nst
         rho_p(i) = electron_density(x(i), y_p(i), z(i),3)
         rho_m(i) = electron_density(x(i), y_m(i), z(i),3)
      enddo
      call density_gradient_select(Nst, x, y_p, z, G_p, G_c_p)
      call density_gradient_select(Nst, x, y_m, z, G_m, G_c_m)
      call xc_f90_gga_exc_vxc(c_func, Nst, rho_p(1), G_p(1), Ec_DFT_xyz(1), Vc_p(1), Vc_sig_p(1) )
      call xc_f90_gga_exc_vxc(c_func, Nst, rho_m(1), G_m(1), Ec_DFT_xyz(1), Vc_m(1), Vc_sig_m(1) )
      do i=1,Nst
         Vc_xyz(i) = Vc_xyz(i)-2._dp*(Vc_sig_p(i)* G_c_p(i,2) - Vc_sig_m(i)* G_c_m(i,2))/ddd
      enddo

      Vc_p=0._dp; Vc_sig_p=0._dp; Vc_m=0._dp; Vc_sig_m=0._dp
      do i=1,Nst
         rho_p(i) = electron_density(x(i), y(i), z_p(i),3)
         rho_m(i) = electron_density(x(i), y(i), z_m(i),3)
      enddo
      call density_gradient_select(Nst, x, y, z_p, G_p, G_c_p)
      call density_gradient_select(Nst, x, y, z_m, G_m, G_c_m)
      call xc_f90_gga_exc_vxc(c_func, Nst, rho_p(1), G_p(1), Ec_DFT_xyz(1), Vc_p(1), Vc_sig_p(1) )
      call xc_f90_gga_exc_vxc(c_func, Nst, rho_m(1), G_m(1), Ec_DFT_xyz(1), Vc_m(1), Vc_sig_m(1) )
      do i=1,Nst
         Vc_xyz(i) = Vc_xyz(i) -2._dp*(Vc_sig_p(i)* G_c_p(i,3) - Vc_sig_m(i)* G_c_m(i,3))/ddd
      enddo

      gga_c=.true.
   case default
   end select

   Vxc_xyz = Vxc_xyz + Vc_xyz

   open(unit=91, file='Vxc_DFT.dat',status='unknown')
   do i=1,Nst
      write(91,'(2e20.10)') r_d(i), Vxc_xyz(i)
   enddo
   close(91)  

end subroutine Plot_Vxc_DFT_r
!=====================================================================================
function electron_density(x,y,z,ispin)
!..ispin: 1: spin up 2: spin down 3: total
!..Global
   use global; use orbocc, ONLY: occnum, vecnat; use functional_m, only:maxorb ;implicit none

!..The function
   real(dp) :: electron_density

!..Arguments
   real(dp), intent(in) :: x, y, z
   integer, intent(in) :: ispin

!..Local
   integer :: mu, nu, ia   
   real(dp) :: P(lnbasis,lnbasis), chi(lnbasis)
   real(dp), external :: f_bas

!..Construct so called density matrix by quantum chemists 
   do mu=1,nbasis
      do nu=1,nbasis
         P(mu,nu) = 0._dp
         do ia=1,maxorb
            P(mu,nu) = P(mu,nu) + occnum(ia,ispin)*vecnat(mu,ia)*vecnat(nu,ia) 
         enddo
      enddo
   enddo

!..Construct the basis functions vectors
   do mu=1,nbasis
      chi(mu) = f_bas(mu, x, y, z)
   enddo
   electron_density = 0._dp 
   do mu=1,nbasis
      do nu=1,nbasis
         electron_density = electron_density + P(mu,nu)*chi(mu)*chi(nu)
      enddo
   enddo

   electron_density = max(electron_density,1e-256_dp)
   return
end function electron_density
!======================================================================
function sph_average_dens(x,y,z,ispin)
!..Global
   use global; use grid_params
   implicit none

!..The function
   real(dp) :: sph_average_dens

!..Arguments
   real(dp), intent(in) :: x, y, z
   integer, intent(in) :: ispin

!..Local 
   integer :: ig, imin, imax
   real(dp) :: r, x1, y1, z1, dd, ss, dmax, dmin
   real(dp) :: electron_density


   if ( Natom > 1) stop 'auxil.f90: sph_average_dens: Natom > 1! What do do?'

   r=sqrt(x*x+y*y+z*z)

   ss=0._dp; dmax=0._dp; dmin=1.e10_dp
   do ig=1,nang
      x1=r*xang(ig); y1=r*yang(ig); z1=r*zang(ig)
      dd=electron_density(x1,y1,z1,ispin)
      if ( dd > dmax ) then
         imax=ig
         dmax=dd
      endif
      if ( dd < dmin ) then
         imin=ig
         dmin=dd
      endif
      ss=ss+dd*wang(ig)
   enddo
   sph_average_dens = ss/4._dp/pi

end function sph_average_dens
!======================================================================

subroutine plot_orb_eff(x0,y0,z0,x1,y1,z1)
   use global; use matrices; use orbocc; use basis_set; implicit none

!..Arguments
   real(dp) :: x0,y0,z0,x1,y1,z1

!..Local
   integer :: ist, i, N
   real(dp) :: x,y,z,d, orb, xstep,ystep,zstep
   real(dp), external :: f_bas

   N=500

   xstep=(x1-x0)/real(N-1); ystep=(y1-y0)/real(N-1); zstep=(z1-z0)/real(N-1)
   x=x0; y=y0; z=z0

   open(unit=53,file='screen_orb_of_r',status='unknown')
   do ist=1,N
      d=sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      orb = 0._dp
      do i=1,nbasis
         orb=orb+ vec_eff(i)*f_bas(i, x, y, z)
      enddo
      write(53,*) d, orb
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo
   close(53)

end subroutine plot_orb_eff

!======================================================================

subroutine plot_orbital(istart, istop,x0,y0,z0,x1,y1,z1)
   use global; use matrices; use orbocc; use basis_set; implicit none

!..Arguments
   integer :: istart, istop
   real(dp) :: x0,y0,z0,x1,y1,z1

!..Local
   integer :: i, ist, N, iorb, irnk
   real(dp) :: x,y,z,d, xstep, ystep, zstep, orb
   real(dp), external :: f_bas
   
   N=3000

   xstep=(x1-x0)/real(N-1); ystep=(y1-y0)/real(N-1); zstep=(z1-z0)/real(N-1)

!  if(ia == 0) then 
!     istart=1; istop=nnatorb
!     occ_tmp(:) = occnum(:,3)
!     do iorb=1,nnatorb
!        occmax=0.d0
!        imax = 1
!        do jorb=1,nnatorb
!           if ( occmax < occ_tmp(jorb) ) then
!              imax = jorb; occmax=occ_tmp(jorb)
!           endif
!        enddo
!        irank(iorb) = imax
!        occ_tmp(imax) =0.d0 
!        print*,'irank',irank(iorb)
!     enddo
!  endif

   open(unit=53,file='orbitals',status='unknown')
   do iorb=istart,istop 
      x=x0; y=y0; z=z0; irnk=iorb !irank(iorb)
      do ist=1,N
         orb = 0.d0
         do i=1,nbasis
            orb=orb+ vecnat(i,irnk)*f_bas(i, x, y, z)
         enddo
         d=sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
         write(53,*) d, orb 
         x=x+xstep; y=y+ystep; z=z+zstep
      enddo
      write(53,*)' '
   enddo
   close(53)

end subroutine plot_orbital
!========================================================================

!subroutine plot_potential_DFT((x0,y0,z0,x1,y1,z1,pfile)
!   use global; use matrices; use orbocc;
!   use DFT_par;implicit none

!   if (.not.allocated(chdens_gr)) then
!      allocate ( chdens_gr(ngrid,3) )
!
!      do ig=1,ngrid
!         chdens_gr(ig,1)=0.d0; chdens_gr(ig,2)=0.d0
!         do ia=1,nnatorb
!            chdens_gr(ig,1) = chdens_gr(ig,1) + occnum(ia,1)*vec_nat_grid(ig,ia)**2
!            chdens_gr(ig,2) = chdens_gr(ig,2) + occnum(ia,2)*vec_nat_grid(ig,ia)**2
!         enddo
!         chdens_gr(ig,3) = chdens_gr(ig,1)+chdens_gr(ig,2)
!      enddo
!   endif

   

!end subroutine plot_orbital_DFT

!========================================================================
subroutine plot_potential(x0,y0,z0,x1,y1,z1,pfile)
   use global; use matrices; use orbocc; use basis_set
   use grid_params; use energies, ONLY: E_HOMO 
   use functional_m, ONLY: scr_ch_const, do_EFO; use DFT_par;implicit none

!..Arguments
   character(10) :: pfile
   real(dp) :: x0, y0, z0, x1, y1, z1

   if ( do_EFO ) then
      call plot_potential_EFO(x0,y0,z0,x1,y1,z1,pfile)
   else
      call plot_potential_EFD(x0,y0,z0,x1,y1,z1,pfile)
   endif
end subroutine plot_potential

!========================================================================
subroutine plot_potential_EFD(x0,y0,z0,x1,y1,z1,pfile)
   use global; use matrices; use orbocc; use basis_set
   use grid_params; use energies, ONLY: E_HOMO,E_HOMO_HF; 
   use functional_m, ONLY: Functional, int_pot_basis, scr_ch_const; use DFT_par;implicit none

!..Arguments
   character(10) :: pfile
   real(dp) :: x0, y0, z0, x1, y1, z1

!..Local
   integer :: ia, ig, ist, iat, N=500
   integer :: n_exp_x, n_exp_y, n_exp_z, iga
   real(dp) :: x,y,z, xstep, ystep, zstep, expob, Xinteg
   real(dp), external :: prod_funct, electron_density
   real(dp), allocatable :: f_bas_pot_plot(:,:), r_d(:), rho(:), dyn(:), dyn_Ha(:), dyn_XC(:), ig_cl(:),scr_charg_plot(:)
   real(dp), allocatable :: dyn_XC_DFT(:)
   real(dp) :: Scharge, Tcharge, x_dg, y_dg, z_dg, rr, Sch, r_min, E_shift
   real(dp) :: rho1, rho2, vx_1, vx_2, ex, vc_1, vc_2, ec
   real(dp),external :: f_bas_pot, X_integr

   allocate ( f_bas_pot_plot(lnbasis_pot, N), r_d(N), rho(N) )
   allocate ( dyn(N), dyn_Ha(N), dyn_XC(N), ig_cl(N), scr_charg_plot(N))
   if ( Functional == 'DFT' ) then 
      allocate ( dyn_XC_DFT(N) )
   endif
   if (.not.allocated(chdens_gr)) then
      allocate ( chdens_gr(ngrid,3) )

      do ig=1,ngrid
         chdens_gr(ig,1)=0._dp; chdens_gr(ig,2)=0._dp
         do ia=1,nnatorb
            chdens_gr(ig,1) = chdens_gr(ig,1) + occnum(ia,1)*vec_nat_grid(ig,ia)**2
            chdens_gr(ig,2) = chdens_gr(ig,2) + occnum(ia,2)*vec_nat_grid(ig,ia)**2
         enddo
         chdens_gr(ig,3) = chdens_gr(ig,1)+chdens_gr(ig,2)
      enddo
   endif

   Tcharge=0._dp; Scharge=0._dp
   do ig=1,ngrid
      Tcharge=Tcharge+ w_grid(ig)*chdens_gr(ig,3)
      Sch=0._dp
      if(int_pot_basis) then
         do ia=1,nbasis_pot
            Sch= Sch + V_bs(ia)*charg_pot(ig,ia)
         enddo
         Scharge= Scharge +  w_grid(ig)*Sch
      endif
   enddo

   print*,'Total     Charge: ', Tcharge
   if(int_pot_basis) then
      print*,'Screening Charge: ', Scharge
   else
      print*,'Screening Charge cannot be calculated for int_pot_basis=.false.'
   endif
   print*,'-------------------------------------------------'

   xstep=(x1-x0)/real(N-1); ystep=(y1-y0)/real(N-1); zstep=(z1-z0)/real(N-1)
   x=x0; y=y0; z=z0
   do ist=1,N
      ig_cl(ist) = 1; r_min=1e30_dp
      do ig=1,ngrid
         x_dg = abs(x-x_grid(ig)); y_dg=abs(y-y_grid(ig)); z_dg=abs(z-z_grid(ig)) 
         rr=sqrt( x_dg*x_dg + y_dg*y_dg + z_dg*z_dg)
         if(rr < r_min) then
            r_min = rr; ig_cl(ist)=ig
         endif
      enddo
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo

   x=x0; y=y0; z=z0
   do ist=1,N
      dyn_Ha(ist) = 0._dp
      f_bas_pot_plot(:,ist) = 0._dp
      if(.not.int_pot_basis) then
         do ia=1,nbasis_pot
            f_bas_pot_plot(ia,ist) = f_bas_pot( ia, x, y, z)
         enddo
      else
         if ( xi_analytic ) then
            do ia=1,nbasis_pot
               iat=iat_bas_pot(ia)
               x_dg = xcoo(iat) - x
               y_dg = ycoo(iat) - y
               z_dg = zcoo(iat) - z
               n_exp_x=Int(expo_x_pot(ia)+0.001_dp)
               n_exp_y=Int(expo_y_pot(ia)+0.001_dp)
               n_exp_z=Int(expo_z_pot(ia)+0.001_dp)
               f_bas_pot_plot(ia,ist)=0._dp
               do iga=1,nga_bas_pot(ia)
                  expob=expo_bas_pot(iga,ia)
                  Xinteg=X_integr(n_exp_x, n_exp_y, n_exp_z, expob, x_dg, y_dg, z_dg)
                  f_bas_pot_plot(ia,ist)=f_bas_pot_plot(ia,ist)+coef_bas_pot(iga,ia)*fact_norm_pot(iga,ia)*Xinteg
               enddo
            enddo
         endif
         scr_charg_plot(ist)=0._dp
         do ia=1,nbasis_pot
            scr_charg_plot(ist) = scr_charg_plot(ist) + V_bs(ia)*f_bas_pot( ia, x, y, z)
         enddo
         rho1=electron_density(x,y,z,1)
         rho2=electron_density(x,y,z,2)
         rho(ist) = rho1+rho2
         if ( Functional == 'DFT' ) then 
            call lsd_slater(rho1,rho2,vx_1,vx_2, ex)
            vc_1=0._dp; vc_2=0._dp; ec=0._dp
            call lsd_VWN_cor(rho1,rho2,vc_1,vc_2,ec)
            dyn_XC_DFT(ist) = vx_1+vc_1 ! For open shell Vx_DFT, Vc_DFT need a second index
         endif
      endif
      do ig=1,ngrid
         x_dg = abs(x-x_grid(ig)); y_dg=abs(y-y_grid(ig)); z_dg=abs(z-z_grid(ig)) 
!        if((2._dp*x_dg > dx) .or. (2._dp*y_dg > dy) .or. (2._dp*z_dg > dz)) then
         if(ig /= ig_cl(ist) ) then
            rr= sqrt( x_dg*x_dg + y_dg*y_dg + z_dg*z_dg)
            rr = w_grid(ig)/rr
            dyn_Ha(ist) = dyn_Ha(ist) + rr*chdens_gr(ig,3)
            if(int_pot_basis .and. (.not.xi_analytic)) then
               do ia=1,nbasis_pot
                  f_bas_pot_plot(ia,ist) = f_bas_pot_plot(ia,ist) + charg_pot(ig,ia)*rr
               enddo
            endif
         endif
      enddo
      r_d(ist) = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo

   E_shift=-(E_HOMO-E_HOMO_HF)
!  E_shift=0.d0
!  E_shift=6.875706d0
   if(int_pot_basis.and.scr_ch_const) E_shift=0._dp
   print*,'Potential Shifted by:', E_shift, ' (a.u.)'

   print*,'-------------------------------------------------'
   print*,'Eff. Potential expansion (in aux basis) coefficients:'
   write(6,'(5f20.10)')V_bs

   open(unit=53, file=pfile, status='unknown')
   open(unit=57, file='potf', status='unknown')
   if(Functional == 'DFT') open(unit=61, file='pot_DFT_XC', status='unknown')
   do ist=1,N
      dyn_XC(ist) = 0._dp; dyn(ist) = 0._dp
      do ia=1,nbasis_pot
         dyn(ist) = dyn(ist) + V_bs(ia)*f_bas_pot_plot(ia,ist)
      enddo
      dyn(ist) = dyn(ist) + E_shift
      dyn_XC(ist) = dyn(ist) - dyn_Ha(ist) 
      write(53,'(f10.3,3e20.10)') r_d(ist), dyn(ist), dyn_Ha(ist), dyn_XC(ist)
      write(57,'(f10.3,16e20.10)') r_d(ist), (f_bas_pot_plot(ia, ist),ia=1,15)
      if(Functional == 'DFT') write(61,'(f10.3,e20.10)') r_d(ist), dyn_XC_DFT(ist)
   enddo
   close(53)
   close(57)
   if(Functional == 'DFT') close(61)

   open(unit=191, file='density', status='unknown')
      do ist=1,N
         write(191,'(3f20.10)') r_d(ist), r_d(ist)*r_d(ist)*rho(ist), r_d(ist)*r_d(ist)*scr_charg_plot(ist)
!        write(191,'(3f20.10)') r_d(ist), rho(ist), scr_charg_plot(ist)
      enddo
   close(191)

   deallocate (f_bas_pot_plot, r_d, rho )

end subroutine plot_potential_EFD
!=================================================================

subroutine pot_nik_plot(x0,y0,z0,x1,y1,z1,pfile)
   use global; use orbocc; use grid_params; use basis_set
   implicit none

!..Arguments
   character(10) :: pfile
   real(dp) :: x0, y0, z0, x1, y1, z1

!..Local
   integer :: ist, ia, l, n_occ, iat, N=500
   real(dp) :: xstep, ystep, zstep, x, y, z
   real(dp) :: x_dg, y_dg, z_dg, Xinteg, expob
   real(dp), allocatable :: dyn_Nik_plot(:), xi_plot(:, :), r_d(:)
   real(dp), allocatable :: vecnat_plot(:, :)
   integer :: n_exp_x, n_exp_y, n_exp_z, iga
   real(dp), external :: X_integr, f_bas

   allocate ( dyn_Nik_plot(N), xi_plot(lnbasis, N), vecnat_plot(lnbasis, N), r_d(N) ) 

   xstep=(x1-x0)/real(N-1); ystep=(y1-y0)/real(N-1); zstep=(z1-z0)/real(N-1)

   n_occ=ibond(1)
   print*,'n_occ',n_occ

   open(unit=153, file=pfile, status='unknown') 

   x=x0; y=y0; z=z0
   do ist=1,N
      xi_plot(:,ist) = 0._dp
      do ia=1,nbasis
         iat=iat_bas(ia)
         x_dg = xcoo(iat) - x
         y_dg = ycoo(iat) - y
         z_dg = zcoo(iat) - z
         n_exp_x=Int(expo_x(ia)+0.001_dp)
         n_exp_y=Int(expo_y(ia)+0.001_dp)
         n_exp_z=Int(expo_z(ia)+0.001_dp)
         xi_plot(ia,ist)=0._dp
         do iga=1,nga_bas(ia)
            expob=expo_bas(iga,ia)
            Xinteg=X_integr(n_exp_x, n_exp_y, n_exp_z, expob, x_dg, y_dg, z_dg)
            xi_plot(ia,ist)=xi_plot(ia,ist)+coef_bas(iga,ia)*fact_norm(iga,ia)*Xinteg
         enddo
      enddo
      do ia=1,n_occ
         vecnat_plot(ia,ist) = 0._dp
         do l=1,nbasis
            vecnat_plot(ia,ist) = vecnat_plot(ia,ist) + vecnat(l,ia)*f_bas(l, x, y, z)
         enddo
      enddo
      dyn_Nik_plot(ist)=0._dp
      do ia=1,n_occ
         do l=1,nbasis
            dyn_Nik_plot(ist) = dyn_Nik_plot(ist) & 
                              + vecnat_plot(ia,ist)*xi_plot(l,ist)*vecnat(l,ia)
         enddo
      enddo
      dyn_Nik_plot(ist) = -0.5_dp*dyn_Nik_plot(ist)
      r_d(ist) = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      write(153,*) r_d(ist), dyn_Nik_plot(ist)
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo

   close(153)

end subroutine pot_nik_plot

!=================================================================
subroutine plot_potential_EFQ(x0,y0,z0,x1,y1,z1)
   use global; use matrices; use orbocc; use basis_set
   use grid_params; use functional_m, only:maxorb; implicit none

!..Arguments
   real(dp) :: x0, y0, z0, x1, y1, z1

!..Local
   integer :: N, ig, ist, m,l, ia, ib
   real(dp) :: x,y,z, xstep, ystep, zstep, rr
   real(dp) :: Scharge, Tcharge, x_dg, y_dg, z_dg, corr_0, rstep
   real(dp), external :: f_bas, electron_density
   real(dp), allocatable :: r_d(:), rho(:), dyn(:), dyn_Ha(:), dyn_XC(:), effscr(:)
   real(dp), allocatable :: effscr_gr(:), EderV(:), Aden(:,:) 

   N=500

   call grid_setup
   call set_grid_volume(corr_0)
   call orbs_on_grid(maxorb)
   allocate ( effscr_gr(ngrid), Aden(lnbasis, lnbasis) )
   allocate ( r_d(N), rho(N), dyn(N), dyn_Ha(N), dyn_XC(N), effscr(N), EderV(N) )

   Tcharge=0._dp; Scharge=0._dp
   do ig=1,ngrid
      effscr_gr(ig) = 0._dp; chdens_gr(ig,1) = 0._dp; chdens_gr(ig,2) = 0._dp
      do m=1,nbasis
         effscr_gr(ig) = effscr_gr(ig) + vec_eff(m)*bas_f_grid(ig,m)
      enddo
      effscr_gr(ig) = veffnorm*effscr_gr(ig)**2
      do ia=1,nnatorb
          chdens_gr(ig,1) = chdens_gr(ig,1) + occnum(ia,1)*vec_nat_grid(ig,ia)**2
          chdens_gr(ig,2) = chdens_gr(ig,2) + occnum(ia,2)*vec_nat_grid(ig,ia)**2
      enddo
      chdens_gr(ig,3) = chdens_gr(ig,1)+chdens_gr(ig,2)
      Tcharge=Tcharge+ w_grid(ig)*chdens_gr(ig,3)
      Scharge=Scharge+ w_grid(ig)*effscr_gr(ig)
   enddo
   print*,'Total     Charge: ', Tcharge
   print*,'Screening Charge: ', Scharge

   print*,'dx,dy,dz', dx,dy,dz
   xstep=(x1-x0)/real(N-1); ystep=(y1-y0)/real(N-1); zstep=(z1-z0)/real(N-1)
   rstep=sqrt(xstep*xstep+ystep*ystep+zstep*zstep)

   if(.not.allocated(dEdVab)) then
      allocate ( dEdVab(lnbasis,lnbasis) )
      dEdVab = 0._dp
   endif

   do m=1,nbasis
      do l=1,nbasis
         Aden(m,l)=0._dp
         do ia=1,nnatorb
            do ib=1,nnatorb
               Aden(m,l) = Aden(m,l) + vecnat(m,ia)*dEdVab(ia,ib)*vecnat(l,ib)
            enddo
         enddo
      enddo
   enddo

   x=x0; y=y0; z=z0

   do ist=1,N
      dyn(ist) = 0._dp;  dyn_Ha(ist) = 0._dp
      rho(ist) = electron_density(x,y,z,3)
      EderV(ist) = 0._dp
      do m=1,nbasis
         do l=1,nbasis
            EderV(ist) = EderV(ist) + Aden(m,l)*f_bas(m,x,y,z)*f_bas(l,x,y,z)
         enddo
      enddo

      effscr(ist)=0._dp
      do m=1,nbasis
         effscr(ist)= effscr(ist)+vec_eff(m)*f_bas(m,x,y,z)
      enddo
      effscr(ist) = veffnorm*effscr(ist)**2
      r_d(ist) = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      do ig=1,ngrid
         x_dg = abs(x-x_grid(ig)); y_dg=abs(y-y_grid(ig)); z_dg=abs(z-z_grid(ig))
         if((2._dp*x_dg > dx) .or. (2._dp*y_dg > dy) .or. (2._dp*z_dg > dz)) then
            rr= sqrt( x_dg*x_dg + y_dg*y_dg + z_dg*z_dg)
            dyn(ist) = dyn(ist) + w_grid(ig)*effscr_gr(ig)/(rr+1.e-8_dp)
            dyn_Ha(ist) = dyn_Ha(ist) + w_grid(ig)*chdens_gr(ig,3)/(rr+1.e-8_dp)
         endif
      enddo
      dyn_XC(ist) = dyn(ist) - dyn_Ha(ist)
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo

   open(unit=191, file='density', status='unknown')
   open(unit=192, file='potential', status='unknown')
      do ist=1,N
         write(191,'(3f20.10)') r_d(ist), r_d(ist)*r_d(ist)*rho(ist), r_d(ist)*r_d(ist)*effscr(ist)
         write(192,'(4f20.10)') r_d(ist), dyn(ist), dyn_Ha(ist), dyn_XC(ist)
         write(193,'(2f20.10)') r_d(ist), EderV(ist)
      enddo
   close(191)
   close(192)

   deallocate ( r_d, rho, dyn, dyn_Ha, dyn_XC, effscr )
   deallocate ( effscr_gr, chdens_gr, Aden, EderV )

end subroutine plot_potential_EFQ

subroutine spin_ensemble()
   use global; use functional_m; use matrices; use orbocc; use energies
   implicit none 
   
!..Local
   real(dp) :: occ_1(lnnatorb,2), occ_2(lnnatorb,2), s, spin, omega
   real(dp) :: xnele_1(2), xnele_2(2), spin_stp
   integer :: i, Nst

   occ_1(:,1)=occnum(:,1); occ_1(:,2)=occnum(:,2); xnele_1(1) = xnele(1); xnele_1(2) = xnele(2) 
   occ_2(:,2)=occnum(:,1); occ_2(:,1)=occnum(:,2); xnele_2(1) = xnele(2); xnele_2(2) = xnele(1) 

   spin=abs(xnele(1)-xnele(2))/2._dp
   if ( spin < 1.e-3_dp) then
      print*,'spin_ensemble: It is pointless to do this for closed shells'
      return
   endif

   open(unit=97,file='ensemble.dat',status='unknown')
   spin_stp = 2.e-2_dp
   Nst = 2._dp*spin/spin_stp + 1
   s=spin
   do i=1,Nst
      omega = (spin+s)/(2._dp*spin)
      occnum(:,1:2) = omega*occ_1(:,1:2) + (1._dp-omega)*occ_2(:,1:2)
      occnum(:,3) = occnum(:,1)+occnum(:,2)
      xnele(1) = omega*xnele_1(1) + (1._dp-omega)*xnele_2(1)
      xnele(2) = omega*xnele_1(2) + (1._dp-omega)*xnele_2(2)
      if(functional == 'BB3' .or. functional == 'AC3' ) then
         call BB_orb_kinds      
      endif
!     print*,'SPIN= ',s
!     print*,occnum(:,1)
!     print*,' '
!     print*,occnum(:,2)
!     print*,' '
!     print*,xnele
!     print*,' '
!     print*,kind_BBC3(:,1)
!     print*,' '
!     print*,kind_BBC3(:,2)
!     print*,'--------------------------------------------'
      call calc_occ_factors(occnum, fac_h, fac_c, fac_e)
      call construct_f(0)
      call total_energy_MO(1) 
      write(97,'(2f20.10)') s, Tot_energy
      s=s-spin_stp
   enddo
   close(97)
end subroutine spin_ensemble
subroutine Lagrange_Hermiticity(Lcond)
!..Global
   use global; use matrices; use orbocc; use functional_m, only:maxorb
   implicit none

!..Arguments
   logical :: Lcond

!..Local
   integer :: ia, ib
   real(dp) :: Fa_fb, sum_f, f_max

   call calc_Lagrange_mult_all()

!  print*,'--------------------------------------------------------------------------'
!  print*,'--------------------------------------------------------------------------'
!  print*,'Divergence from Hermiticity of the Lagrangian Matrix:'
!  print*,'          a           b       F_ab-Fba           (F_ab-Fba)/(F_ab+Fba)'
!  print*,'--------------------------------------------------------------------------'
   f_max=0._dp
   sum_f=0._dp
   do ia=2,maxorb
      do ib=1,ia-1
         Fa_fb=e_ab(ia,ib)-e_ab(ib,ia)
         if(abs(Fa_fb) > f_max) f_max=abs(Fa_fb)
         sum_f = sum_f + abs(Fa_fb)
!        Fa_fb=2.d0*Fa_fb/(e_ab(ia,ib)+e_ab(ib,ia))
!        if(abs(Fa_fb) > 1.d-2.and.abs(e_ab(ia,ib)) > 1d-5) then
!           write(6,*)ia,ib,Fa_fb, e_ab(ib,ia)
!        endif
      enddo
   enddo
   sum_f=2._dp*sum_f/maxorb/(maxorb-1)

   Lcond=(sum_f < 1.e-4_dp)
   Lcond= Lcond.and.( f_max < 5.e-4_dp)
   print*,'Hermiticity test:', 'Average:', sum_f, 'Max:',f_max

end subroutine Lagrange_Hermiticity

!================================================================================
subroutine three_overlap(fname)
!..Global
   use global; use matrices; use grid_params; implicit none

!..Argument
  character(60) :: fname 

!..Local
   integer :: i,j,k, ig, iint, nnb, npb
   integer(8) :: ipack
   real(dp) :: s
   logical :: ffile

   if(.not.grid_stat) stop 'three_overlap: grid_setup must be called first'
   if(.not.allocated(ovlap3)) allocate(ovlap3(lnbasis,lnbasis,lnbasis_pot))

   ovlap3=0._dp

   inquire(file=fname,exist=ffile)
   if(ffile) then
      open(unit=77,file=fname,status='old',form='unformatted')
      read(77) nnb, npb
      if(nnb /= nbasis .and. npb /= nbasis_pot) then
         print*,'auxil:three_overlap: Wrong dimensions in ',fname
         stop 'Please remove or rename and restart'
      endif
      print*,'Reading three function overlaps ...'
      do iint=1,nbasis*nbasis*nbasis_pot
         read(77,end=12) ipack, s
         call unpack3(i,j,k,ipack)
         ovlap3(i,j,k)=s
         ovlap3(j,i,k)=s
      enddo
 12   continue
      close(77)
      return
   else
   open(unit=77,file=fname,status='unknown',form='unformatted')
   write(77) nbasis, nbasis_pot
   print*,'Calculate three function overlaps ...'
!.....Calculate overlap elements
   do i=1,nbasis
      do j=1,i
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,s,ig,ipack)
!$OMP DO
         do k=1,nbasis_pot
            s=0._dp
            do ig=1,ngrid
                if (abs(bas_f_grid(ig,i)) > zero .and. abs(bas_f_grid(ig,j)) > zero &
                    .and. abs(bas_f_grid_pot(ig,k))>zero) then
                    s = s + w_grid(ig)*bas_f_grid(ig,i)*bas_f_grid(ig,j)*&
                        bas_f_grid_pot(ig,k)
                endif
            enddo
            if( abs(s) > zero) then
               ovlap3(i,j,k) = s; ovlap3(j,i,k) = s
               call pack3(i,j,k,ipack)
!$OMP CRITICAL
               write(77) 'i',ipack,s
!$OMP END CRITICAL
            endif      
         enddo 
!$OMP END DO
!$OMP END PARALLEL
      enddo          
   enddo
   close(77)
   endif
   print*,'... Done'
end subroutine three_overlap

!================================================================================
subroutine three_xi_ovlap(fname)
!..Global
   use global; use matrices; use grid_params; implicit none

!..Argument
  character(60) :: fname 

!..Local
   integer :: i,j,k, ig, iint, nnb, npb
   integer(8) :: ipack
   real(dp) :: s
   logical :: ffile

   if(.not.grid_stat) stop 'three_xi_ovlap: grid_setup must be called first'
   if(.not.allocated(ovlap3_xi)) allocate(ovlap3_xi(lnbasis_pot,lnbasis_pot,lnbasis_pot))

   ovlap3_xi=0._dp

   inquire(file=fname,exist=ffile)
   if(ffile) then
      ovlap3_xi = 0._dp
      open(unit=77,file=fname,status='old',form='unformatted')
      read(77) nnb, npb
      if(nnb /= nbasis_pot .and. npb /= nbasis_pot) then
         print*,'auxil:three_overlap: Wrong dimensions in ',fname
         stop 'Please remove or rename and restart'
      endif
      print*,'Reading three function overlaps ...'
      do iint=1,nbasis_pot*nbasis_pot*nbasis_pot
         read(77,end=12) ipack, s
         call unpack3(i,j,k,ipack)
         ovlap3_xi(i,j,k)=s
         ovlap3_xi(j,i,k)=s
      enddo
 12   continue
      close(77)
      return
   else
   open(unit=77,file=fname,status='unknown',form='unformatted')
   write(77) nbasis_pot, nbasis_pot
   print*,'Calculate three function overlaps ...'
!.....Calculate overlap elements
   do i=1,nbasis_pot
      do j=1,i
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,s,ig,ipack)
!$OMP DO
         do k=1,nbasis_pot
            s=0._dp
            do ig=1,ngrid
                if (abs(bas_f_grid_pot(ig,i)) > zero .and. abs(bas_f_grid_pot(ig,j)) > zero &
                    .and. abs(bas_f_grid_pot(ig,k))>zero) then
                    s = s + w_grid(ig)*bas_f_grid_pot(ig,i)*bas_f_grid_pot(ig,j)*&
                        bas_f_grid_pot(ig,k)
                endif
            enddo
            if( abs(s) > zero) then
               ovlap3_xi(i,j,k) = s; ovlap3_xi(j,i,k) = s
               call pack3(i,j,k,ipack)
!$OMP CRITICAL
               write(77) ipack,s
!$OMP END CRITICAL
            endif      
         enddo 
!$OMP END DO
!$OMP END PARALLEL
      enddo          
   enddo
   close(77)
   endif
   print*,'... Done'

end subroutine three_xi_ovlap

!================================================================================
subroutine correl_entropy()
   use global  
   use orbocc
   implicit none

   real(dp) :: c_entropy, ent, q_c, occ
   integer :: ia, ispin

   print*,zero,one,real(nele(3))

   c_entropy=0._dp; q_c=0._dp
   do ispin=1,2
      do ia=1,nnatorb
         occ=min(max(occnum(ia,ispin),zero),one)
         ent=-occ*log(occ)
         if(occ < 0.46_dp) q_c=q_c+occ
!        print*,ia,ent
         c_entropy=c_entropy+ent
      enddo
   enddo

   print*,'Correlation entropy:', c_entropy
   print*,'Charge on weakly occupied orbitals:',q_c

end subroutine correl_entropy
!=====================================================================================
subroutine Make_projector(AtA, Pro_null, zz)
!..Global
   use global; use matrices; use orbocc; implicit none
!..Arguments
   real(dp) :: AtA(lnbasis_pot, lnbasis_pot), Pro_null(lnbasis_pot, lnbasis_pot)
   real(dp) :: zz

!..Local
   real(dp), allocatable :: A(:,:)
   real(dp), allocatable :: B(:,:)
   real(dp), allocatable :: vectr(:,:), eigs(:)
   real(dp), allocatable :: Zmat(:)
   integer :: lndim, ndim, mu, nu, ku, lu, info
!..Definitions for LAPACK routine
   character(1) :: JOBZ, RANGE, UPLO
   INTEGER :: IL, ITYPE,  IU,  LDA,  LDB,  LDZ
   INTEGER LWORK, M, N
   real(dp) :: VL,VU
   real(dp) :: ABSTOL
   real(dp), allocatable :: WORK(:)
   integer, allocatable :: IWORK(:)
   integer, allocatable :: ifail(:)
   
   real(dp) :: DLAMCH
   external DLAMCH

!..Local array allocation
   ndim=nbasis_pot
   lndim=ndim
   allocate( A(lndim,lndim), B(lndim,lndim), Zmat(ndim), &
               vectr(lndim,lndim), eigs(lndim), &
               WORK(10*lndim), IWORK(5*lndim), ifail(lndim) )

   ITYPE=1; JOBZ='V'; RANGE='I'; UPLO='L'; N = ndim; LDA = lndim
   LDB = lndim; VL = 0._dp; VU = 0._dp; IL = 1; IU = ndim
   LDZ = lndim; LWORK=10*ndim

   ABSTOL = 2._dp*DLAMCH('S')
!  ABSTOL = 1.d-30
      
   do mu=1, ndim   
      do nu=1, mu
         A(mu,nu) = AtA(mu,nu)
         B(mu,nu) = ovlap_pot(mu,nu)
         B(nu,mu) = B(mu,nu)
         A(nu,mu) = A(mu,nu) 
      enddo
!     B(mu,mu)=1.d0
   enddo

   WORK=0._dp
   IWORK=0
   call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
                LDB,  VL,  VU, IL, IU, ABSTOL, M, eigs, vectr, &
                LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

   if (info /= 0) then
      print*,'Info:',info
      if(info < 0) print*,info, '-th argument of DSYGVX has wrong value'
      if(info > 0 .and. info <= lndim ) print*,info, '-th eigenvalue not converged'
      if(info > lndim .and. info <= 2*lndim ) &
         print*,info-lndim, '-th leading minor of matrix not positive definite'
      stop 'Make_projector: Diagonalization failed!'
   endif

   print*,'Eigs proj:'
   print*,eigs
   print*,'-----------------------------------'
 
   do ku=1,ndim
      if (abs(eigs(ku)) > zz ) then 
         eigs(ku)=0._dp
      else
         eigs(ku)= 1._dp
      endif
   enddo

!  do mu=1,ndim
!     do nu=1,mu
!        Pro_null(mu,nu)=0.d0
!        do ku=1,ndim
!           Pro_null(mu,nu) = Pro_null(mu,nu) + vectr(mu,ku)*eigs(ku)*vectr(nu,ku)
!        enddo
!        Pro_null(nu,mu)= Pro_null(mu,nu)
!     enddo
!  enddo

   do mu=1,ndim
      do nu=1,ndim
         Pro_null(mu,nu)=0._dp
         do ku=1,ndim
            do lu=1,ndim
            Pro_null(mu,nu) = Pro_null(mu,nu)+eigs(ku)*vectr(mu,ku)*ovlap_pot(lu,nu)*vectr(lu,ku)
            enddo
         enddo
      enddo
   enddo  

end subroutine Make_projector
!=====================================================================================

!---------------------------------------------------------------
function Boys(n,r)
use global; implicit none
!..The function
   real (dp) :: Boys

!..Arguments
   integer, intent(IN) :: n 
   real (dp), intent(IN) :: r

!..Local
   integer, parameter :: N_max=250
   integer :: i
   real(dp) :: fc, term, xna, qq
   real(dp),external :: S14AAF
   external :: S14BAF

!..Use two approximations:
!  1) Expansion (not valid for large r)
!  2) asymprotic form for large r
!  For rmax=15 Nmax=50 accuracy up to 10^-8 

   if ( r <= 25._dp ) then
      term=1.0d0
      xna=real(n)+0.5_dp
      Boys=term/xna
      do i=1,N_max
         term=-term*r/real(i)
         qq=term/(xna+real(i))
         if (abs(qq) .lt. 1.e-10_dp) exit
         Boys=Boys+qq
      enddo
   else 
      qq=real(2*n+1)
      fc=1._dp/2._dp**(n+1)
      do i=2*n-1,1,-2
         fc=fc*real(i)
      enddo
      Boys=2._dp*fc*sqrt(pi/r**qq) 
   endif
   

!  TOL=1.d-16
!  xna=n+0.5d0
!  
!  term= S14AAF(xna, ifail)
!  call S14BAF (xna, r, TOL, dum, gma, IFAIL)
!  
!  boys=term*(1.d0-gma)/r**xna

end function Boys

!---------------------------------------------------------------
function Gma_mine(n)
!..Global
   use global; implicit none

!..The function
   real(dp) :: Gma_mine

!..Argument
   integer, intent(IN) :: n

   select case(n)
      case(0) 
         Gma_mine=sqpi
      case(1)
         Gma_mine=1._dp
      case(2)
         Gma_mine=sqpi*0.5_dp
      case(3)
         Gma_mine=1.d0
      case(4)
         Gma_mine=sqpi*0.75_dp
      case(5)
         Gma_mine=2._dp
      case(6)
         Gma_mine=1.875_dp*sqpi
      case(7) 
         Gma_mine=6._dp
      case(8)
         Gma_mine=6.5625_dp*sqpi
     case(9)
         Gma_mine=24._dp
      case default
         stop 'Gma_mine in auxil.F90: Gamma((n+1)/2) function not defined please define'
   end select
      
end function Gma_mine

!------------------------------------------------------------------------------------------------
function F_N_ov_m(N, m)
!..N over m: N!/[(N-m)!m!]
   use global; implicit none

!..The function
   real(dp) :: F_N_ov_m

!..Arguments
   integer, intent(IN) :: N, m

!..Local
   integer :: i, Nmm

   if( N<m .or. N<0 .or. m<0 ) stop 'F_N_ov_m in oep.F90: N must be larger or equal to m'
   
   F_N_ov_m=1._dp
   if (N == 0) then 
!     return
   elseif ( m==0 ) then
!     return
   else
      Nmm=N-m
      do i=1,Nmm
         F_N_ov_m=F_N_ov_m*real(m+i)/real(i)   
      enddo
   endif

end function F_N_ov_m

!----------------------------------------------------------------------------------------------
function X_integr(ne_x, ne_y, ne_z, expob, xd, yd, zd)
!..Global
   use global; implicit none

!..The function
   real(dp) :: X_integr

!..Arguments
   integer :: ne_x, ne_y, ne_z
   real(dp) :: expob, xd, yd, zd

!..Local
   Integer :: N_tot, k1, k2, k3, K_tot, kz, kt2, ip
   real(dp) :: x_k, y_k, z_k, ff, rbz, sg_ip, X_i
   real(dp), external :: F_N_ov_m, Gma_mine, Boys

!  rbz=1000.d0/99999; X_i=0.d0
!  do ip=1,100000
!     write(9,'(2e20.10)')X_i,Boys(3,X_i)
!     X_i=X_i+rbz
!  enddo
!  stop 'aaa'
   
   N_tot =  ne_x + ne_y + ne_z

   rbz=expob*(xd*xd+yd*yd+zd*zd)
   X_i=0._dp
   do k1=0,ne_x,2
      x_k= F_N_ov_m(ne_x, k1)*Gma_mine(k1)*xd**(ne_x-k1)
      do k2=0,ne_y,2
         y_k= F_N_ov_m(ne_y, k2)*Gma_mine(k2)*yd**(ne_y-k2)
         do k3=0,ne_z,2
            z_k= F_N_ov_m(ne_z, k3)*Gma_mine(k3)*zd**(ne_z-k3)
            K_tot = k1+k2+k3;  kt2=K_tot/2; kz=1+kt2
            ff=0._dp
            sg_ip=1._dp
            do ip=0,kt2
               ff=ff + sg_ip*F_N_ov_m(kt2,ip)*Boys(N_tot-K_tot+ip, rbz)
               sg_ip=-sg_ip
            enddo
            X_i=X_i+x_k*y_k*z_k*ff/expob**kz
         enddo
      enddo
   enddo

   X_integr = X_i*(-1._dp)**N_tot/sqpi
   return
end function X_integr

!=====================================================================================
subroutine diagon_H(info)
   implicit none
   integer :: info
!  call diagon_H_NAG(info)
   call diagon_H_lapack(info)
end subroutine diagon_H
!=====================================================================================
subroutine diagon_H_lapack(info)
!-----------------------------------------------------------------------
! Wrapper for the diagonalizing routine of LAPACK
! Diagonalizer of the effective Hamiltonian
!-----------------------------------------------------------------------

!..Global
   use global
   use matrices
   use orbocc
   implicit none

!..Argument
   integer :: info

!..Local
   real(dp), allocatable :: A(:,:)
   real(dp), allocatable :: B(:,:)
   real(dp), allocatable :: vectr(:,:), eigs(:)

!..Definitions for LAPACK routine
   character(1) :: JOBZ, RANGE, UPLO
   INTEGER :: IL, ITYPE,  IU,  LDA,  LDB,  LDZ
   INTEGER LWORK, M, N
   real(dp) :: VL,VU
   real(dp) :: ABSTOL
   real(dp), allocatable :: WORK(:)
   integer, allocatable :: IWORK(:)
   integer, allocatable :: ifail(:)
   
   real(dp) :: DLAMCH
   external DLAMCH

!..Local array allocation
   allocate( A(lnbasis,lnbasis), B(lnbasis,lnbasis), vectr(lnbasis,lnbasis),eigs(lnbasis) ) 
   allocate( WORK(10*lnbasis), IWORK(5*lnbasis), ifail(lnbasis) )

   ITYPE=1; JOBZ='V'; RANGE='I'; UPLO='L'; N = nbasis; LDA = lnbasis
   LDB = lnbasis; VL = 0._dp; VU = 0._dp; IL = 1; IU = nbasis
   LDZ = lnbasis; LWORK=10*nbasis

   ABSTOL = 2._dp*DLAMCH('S')
!  ABSTOL = 1.e-30_dp
      
! The matrices that define the eigenvalue problem A x = B lambda x
   A = H_eff
   B = ovlap
   WORK=0._dp
   IWORK=0
   call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
                LDB,  VL,  VU, IL, IU, ABSTOL, M, eigs, vectr, &
                LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

   if (info /= 0) then
      print*,'Info:',info
      if(info < 0) print*,info, '-th argument of DSYGVX has wrong value'
      if(info > 0 .and. info <= lnbasis ) print*,info, '-th eigenvalue not converged'
      if(info > lnbasis .and. info <= 2*lnbasis ) &
         print*,info-lnbasis, '-th leading minor of matrix not positive definite'
      stop 'diagon_H:Diagonalization failed!'
   endif

!..Now store eigs --> ennat, vectr --> vecnat
   ennat=eigs; vecnat=vectr
  
   deallocate ( A, B, vectr, eigs, WORK, IWORK, ifail )

end subroutine diagon_H_lapack

subroutine diagon_H_NAG(info)
!-----------------------------------------------------------------------
! Wrapper for the diagonalizing routine of NAG
! Diagonalizer of the effective Hamiltonian
!-----------------------------------------------------------------------

!..Global
   use global
   use matrices
   use orbocc
   implicit none

!..Argument
   integer :: info

!..For NAG routine
   integer :: ITYPE=1, LWORK, IFAIL
   character(1) :: JOB='V', UPLO='U'
   real(dp), allocatable :: A(:,:), B(:,:), WORK(:), W(:)

   LWORK = 64*lnbasis

   allocate( A(lnbasis,lnbasis), B(lnbasis,lnbasis), W(lnbasis), WORK(LWORK) )

   IFAIL=0

   A=H_eff; B=ovlap; 

   call F02FDF( ITYPE, JOB, UPLO, lnbasis, A, lnbasis, B, lnbasis, W, &
                WORK, LWORK, IFAIL)

   if(ifail/=0) then
      print*,'ifail',ifail
      stop 'auxil:diagon_H_NAG: ifail not zero'
   endif
   info=ifail
   ennat=W; vecnat=A   

   deallocate ( A, B, W, WORK )

end subroutine diagon_H_NAG


subroutine MATMUL_N(A,B,C,i1,j1)
!..Global
   use global; implicit none

   integer :: i,j,i1,j1
   real(dp) :: A(i1,j1), B(j1), C(i1)

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i, j)
!$OMP DO
   do i=1,i1
      C(i)=0._dp
      do j=1,j1
         C(i)=C(i)+A(i,j)*B(j)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine MATMUL_N
subroutine pack3(i,j,k,ipack)
!..Global
   use global; implicit none
!..Arguments
   integer :: i,j,k
   integer(8) :: ipack

!..Local
   integer, parameter :: ibase=lbsa
   integer(8), parameter :: ibase2=lbsa*lbsa

   ipack=ibase2*i+ibase*j+k

   return
end subroutine pack3
subroutine unpack3(i,j,k,ipack)
!..Global
   use global; implicit none
!..Arguments
   integer :: i,j,k
   integer(8) :: ipack

!..Local
   integer, parameter :: ibase=lbsa
   integer(8), parameter :: ibase2=lbsa*lbsa

   i=ipack/ibase2
   ipack=ipack-i*ibase2
   j=ipack/ibase
   k=ipack-j*ibase

   return
end subroutine unpack3
subroutine pack4(i,j,k,l,ipack)
!..Global
   use global; implicit none
!..Arguments
   integer :: i,j,k,l
   integer(8) :: ipack

!..Local
   integer, parameter :: ibase=lbsa
   integer(8), parameter :: ibase2=lbsa*lbsa
   integer(8), parameter :: ibase3=ibase2*lbsa

   ipack=ibase3*i+ibase2*j+ibase*k+l

   return
end subroutine pack4
subroutine unpack4(i,j,k,l,ipack)
!..Global
   use global; implicit none
!..Arguments
   integer :: i,j,k,l
   integer(8) :: ipack

!..Local
   integer, parameter :: ibase=lbsa
   integer(8), parameter :: ibase2=lbsa*lbsa
   integer(8), parameter :: ibase3=ibase2*lbsa

   i=ipack/ibase3
   ipack=ipack-i*ibase3
   j=ipack/ibase2
   ipack=ipack-j*ibase2
   k=ipack/ibase
   l=ipack-k*ibase

   return
end subroutine unpack4

!=====================================================================================
!Check convergence routine
subroutine check_conv(iscf, conv_E1, conv_n, DE, dn, Lconv)
!..Global
   use global; use check_conv_m; use energies; use orbocc 
   implicit none

!..Arguments
   integer, INTENT(in) :: iscf
   logical, INTENT(out) :: Lconv
   real(dp), INTENT(in) :: conv_E1, conv_n

!..Local
   real(dp) :: DE, dn

   integer :: ia, l, m

   DE=(Tot_energy-Tot_old)/abs( Tot_old )

   if(iscf==1) then
      allocate(Dmat(lnbasis,lnbasis),Dmat_old(lnbasis,lnbasis))
      Dmat_old=0._dp
   endif

   do l=1,nbasis
      do m=1,l
         Dmat(l,m)=0._dp
         do ia=1,nnatorb
            Dmat(l,m)=Dmat(l,m)+occnum(ia,3)*vecnat(l,ia)*vecnat(m,ia)
         enddo
      enddo
   enddo

   dn=0._dp
   do l=1,nbasis
      do m=1,l
         dn=dn+abs(Dmat(l,m)-Dmat_old(l,m))
      enddo
   enddo
   dn=dn/float(nbasis*nbasis)

   if(iscf<3) then
      Lconv=.false.
   else
!     Lconv = ( abs(DE) < conv_E1 )
!     Lconv = Lconv .and. (dn < conv_n )
      Lconv = ( abs(DE) < conv_E1 .and. abs(DE_old) < conv_E1 )
      Lconv = Lconv .and. (dn < conv_n .and. dn_old < conv_n)
   endif

!  print*,'CONV:',DE, DE_old
!  print*,'CONV:',Tot_energy,dE, dn
!  print*,'CONV:'

   Tot_old = Tot_energy; DE_old=DE
   Dmat_old=Dmat; dn_old=dn

end subroutine check_conv



subroutine diagon_ovlap_pot()

!..Global
   use global
   use matrices
   use orbocc
   implicit none

!..Local
   real(dp), allocatable :: A(:,:)
   real(dp), allocatable :: B(:,:)
   real(dp), allocatable :: vectr(:,:), eigs(:)
   real(dp) :: eig_min, eig_max

   integer :: ia

!..Definitions for LAPACK routine
   character(1) :: JOBZ, RANGE, UPLO
   INTEGER :: IL, INFO, ITYPE,  IU,  LDA,  LDB,  LDZ
   INTEGER LWORK, M, N
   real(dp) :: VL,VU
   real(dp) :: ABSTOL
   real(dp), allocatable :: WORK(:)
   integer, allocatable :: IWORK(:)
   integer, allocatable :: ifail(:)
      
   real(dp) :: DLAMCH
   external DLAMCH
      
   ABSTOL = 2._dp*DLAMCH('S')

!..Local array allocation
   allocate( A(lnbasis_pot,lnbasis_pot), B(lnbasis_pot,lnbasis_pot), &
               vectr(lnbasis_pot,lnbasis_pot), eigs(lnbasis_pot), &
               WORK(10*lnbasis_pot), IWORK(5*lnbasis_pot), ifail(lnbasis_pot) )

   ITYPE=1; JOBZ='V'; RANGE='I'; UPLO='L'; N = nbasis_pot; LDA = lnbasis_pot
   B = 0._dp; LDB = lnbasis_pot; VL = 0._dp; VU = 0._dp; IL = 1; IU = nbasis_pot
   LDZ = lnbasis_pot; LWORK=10*nbasis_pot
      
   A=ovlap_pot

   do ia=1,nbasis_pot
      B(ia,ia)=1._dp
   enddo
   
   call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
                LDB,  VL,  VU, IL, IU, ABSTOL, M, eigs, vectr, &
                LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

!..Now store eigs --> ennat, vectr --> vecnat
!  ennat=eigs; vecnat=vectr

   print*,'Eigenvalues of auxiliary basis overlap:'
   eig_max=-1.e20_dp; eig_min=1.e20_dp
   do ia=1,nbasis_pot
      if(eigs(ia) > eig_max) eig_max=eigs(ia)
      if(eigs(ia) < eig_min) eig_min=eigs(ia)
      print*,ia,eigs(ia)
   enddo
   print*,'Maximum ovlap_pot eigenvalue:', eig_max
   print*,'Minimum ovlap_pot eigenvalue:', eig_min
   print*,'ovlap_pot matrix condition:',abs(eig_max/eig_min)
   
   stop 'In diagon_ovlap_pot'

   deallocate ( A, B, vectr, eigs, WORK, IWORK, ifail )
end subroutine diagon_ovlap_pot

!---------------------------------------------------------------------------
! Calculates the average non diagonal matrix element of gamma
subroutine average_nondiagonal()
   use global; use orbocc; use matrices, ONLY:ovlap
   implicit none

   real(dp) :: TT(lnbasis, lnnatorb), D1nm, gamma_12
   integer :: k, m, n, ia

!..Use this TT to calculate the matrix element
   do ia=1,nnatorb
      do m=1,nbasis
         TT(m,ia)=0.0_dp
         do k=1,nbasis
            TT(m,ia) = TT(m,ia) + ovlap(m,k)*vecnat(k,ia) 
         enddo
      enddo
   enddo

! Set TT= vecnat instead to use the spectral form
!  TT=vecnat

   gamma_12=0.0_dp
   do m = 2, nbasis
!     do n = 1, m-1
         n=m-1
         D1nm = 0.0_dp
         do ia= 1, nnatorb
            D1nm = D1nm + occnum(ia,3)*TT(m,ia)*TT(n,ia)
         enddo
         print*,m,D1nm
         gamma_12=gamma_12 + D1nm*D1nm
!     enddo
   enddo 

   print*,xnele(3)
   gamma_12 = sqrt(2.0_dp*gamma_12/(xnele(3)*(xnele(3)-1.0_dp)))

   print*,'Average non-diagonal element of Gamma:', gamma_12

!  stop 'In average_nondiagonal'

end subroutine average_nondiagonal

!-------------------------------------------------------------------------------
subroutine read_nat_orbs(filename)
!..Global
   use global; use orbocc; implicit none

!..Arguments
   character(60) :: filename

!..Local
   integer :: limit=10000000, linear_ind, nblock, i, j, ia, ib, iamin, iamax
   logical :: linear_depend
   character(20) :: line
   character(45) :: dumm
   character(30) :: dumm4

   open(unit=11, file=filename, status='old')
   !..Check for Linearly dependent orbitals (less orbitals printed)
   rewind(11)
   linear_depend=.false.
   do i=1, limit
      read(11,'(a20)',end=135) line
        if(line == ' NUMBER OF LINEARLY ') then
           linear_depend=.true.
           exit
        endif
   enddo
 135  continue
    rewind(11)
    if (linear_depend) then
       do i=1,limit
          read(11,'(a20)') line
          if(line == ' NUMBER OF LINEARLY ') exit
       enddo
       read(11,'(a45,i5)')dumm, linear_ind
    end if ! end if linar_depend=true
    rewind(11)

!..linear_ind is the number of linearly independent orbitals to be read
   if (.not.linear_depend)  linear_ind = nbasis
 
   do i=1,limit
      read(11,'(a26)',end=55) dumm4
      if(dumm4 == '          NATURAL ORBITALS' .or. &
         dumm4 == '          MCSCF NATURAL OR' ) then
         read(11,*)
         goto 56
      endif
   enddo
55 stop 'gamess:string: NATURAL ORBITALS not found in log file'

56 continue
   print*,'NATURAL ORBITALS FOUND!'

  nblock = linear_ind/5
     do ib=0,nblock
        read(11,*)
        read(11,*)
        read(11,*)
        read(11,'(a22)')dumm4
        do j=1,nbasis
           iamin= ib*5 +1
           iamax= min(iamin+4,linear_ind)
           read(11,'(15x,5e21.12)') (vecnat(j,i),i=iamin,iamax)
        enddo
     enddo

   print*,'NATURAL ORBITALS READ!'
   rewind(11)
   
   do i=1,limit
      read(11,'(a26)',end=55) dumm4
      if(dumm4 == '          NATURAL ORBITALS' .or. &
         dumm4 == '          MCSCF NATURAL OR' ) then
         read(11,*)
         goto 156
      endif
   enddo
   stop 'gamess:string: NATURAL ORBITALS not found in log file'

156 continue
   print*,'READ OCCUPATIONS'

   iamin=1
   iamax=min(5,linear_ind)
   do ib=0,linear_ind/5
      read(11,*)
      read(11,*)
      print*,'aaaa',ib, iamin, iamax
      read(11,*) (occnum(ia,3),ia=iamin,iamax)
      read(11,*)
      do j=1,nbasis
         read(11,*) 
      enddo
      iamin=iamin+5
      iamax=iamax+5
      iamax=min(iamax,linear_ind)
   enddo

   print*, 'OCCUPATIONS READ'
   occnum(:,1)=0.5d0*occnum(:,3)
   occnum(:,2)=occnum(:,1)

   close(11)
end subroutine read_nat_orbs

