subroutine vary_occ_th_n(xmu)
!------------------------------------------------------------------------
! This subroutine finds the optimal occ. numbers. The constraint of 
! \sum n_i = N is enforced during the variation. mu is also a parameter
! that converges to the correct value at the end
!
! INPUT:
!   xmu(2)              : on input contains the initial values of \mu
!                         for spin_up and spin_down
!   xmu_range           : a range for the search of \mu. 
!
! In common:
!   occnum(lnnatorb)    : on input contains the starting occ numbers
!   occnum(lnnatorb)    : on output the refined occ numbers
!
! OUTPUT:
!   xmu(2)              : on output the refined value of \mu
! 
!         created by N. N. Lathiotakis and I.Theophilou
! 
!------------------------------------------------------------------------

! Global
   use global; use functional_m; use orbocc; use vary_occ_par
   implicit none

!..Aguments
   real(dp), intent(inout) :: xmu(2)

!..Local variables
   external solve_occ_th_n
   real(dp) :: solve_occ_th_n
   real(dp), allocatable :: occnumtmp(:,:)

   real(dp) :: xmui, sum_a
   real(dp) :: sum_occ(3)
   integer :: ia, iterb
   integer :: nspin, isc, nsc
   character(4) :: spin_nam(2)

!..Local array allocation
   allocate (occnumtmp(lnnatorb,3))

   spin_nam(1) = ' UP '
   spin_nam(2) = 'DOWN'

   if(cl_shell) then
      nsc=1
      nspin = 1
   else
      nspin = 2
      nsc=10
   endif
   
   occnumtmp(:,:2) = occnum(:,:2)
   where(occnumtmp(:,:2) < smallocc) occnumtmp(:,:2)=smallocc
   where(occnumtmp(:,:2) > 1.d0-smallocc) occnumtmp(:,:2)=1.d0-smallocc
   occnumtmp(:,3) = occnumtmp(:,1) +  occnumtmp(:,2)
   
   do isc=1,nsc
      print*,'======================================'
      print*,'isc=',isc
      do ispin = 1, nspin
!     if(nele(ispin) /= 0) then !to run He triplet
!          Define the opposite spin:
         if(ispin == 1) then
            ispin_op = 2
         else
            ispin_op = 1
         endif
 
         sum_a = solve_occ_th_n(xmu(ispin),occnumtmp)
         print*,'sum occ - nel for optimal mu:',sum_a+nele(ispin)
         print*,'mu of spin ',spin_nam(ispin),' :',xmu(ispin)
         xmui=xmu(ispin)
 
         occnum(:,ispin) = occnumtmp(:,ispin)
         occnum(:,3) = occnumtmp(:,3)
         sum_occ(ispin) = sum(occnum(:,ispin))
 
         print*,'------------------------------------------------------'
         print*,'SPIN: ',spin_nam(ispin),' :'
         print*,'Optimal mu:',xmui,' Nu of mu iter.:',iterb
         print*,'New occ numbers:'
         write(6,'("NEWOCC",5f12.8)') (occnum(ia,ispin), ia=1,5)
         write(6,'("      ",5f12.8)') (occnum(ia,ispin), ia=6,nnatorb)
         print*,'Sum of occ numbers:', sum_occ(ispin)

!     endif !nele(ispin) /= 0
      enddo !ispin = 1,nspin
   enddo !isc

   if(cl_shell) then
      xmu(2) = xmu(1)
      sum_occ(2) = sum_occ(1)
      occnum(:,2)=occnum(:,1)
      print*,'For closed shell singlet spin DOWN quantities'
      print*,'are the same as spin UP.'
   endif
   sum_occ(3) = sum_occ(1) + sum_occ(2)
   print*,'TOTAL NUMBER OF ELECTRONS',sum_occ(3) 
   print*,'------------------------------------------------------'

!..After all occ. numbers are redefined call 
   call calc_occ_factors(occnum, fac_h, fac_c, fac_e)

   deallocate (occnumtmp)

end subroutine vary_occ_th_n
!-----------------------------------------------------------------------

function solve_occ_th_n(xmu,occnumtmp)
!----------------------------------------------------------------------
! This subroutine solves the minimization problem with respect to the 
! occ numbers. Method: steepest descent.
!
! On output solve_occ_th_n = sum of the occupation numbers * 2 - nele
!
! INPUT:
!         xmu:  \mu initial for the appropriate spin
! Also from the common:
!      occnum: Old occupation numbers
!    niter_on: Max. number of directions in steepest descent method
!    nbracket: max. number of bracketing steps in each direction
!    niteron1: max. number of divisions in each direction
!     crit_on:   conv criterion for occ numbers (all derivatives must be
!               smaller than this number).
!    crit_on1:  linear search conv. criterion
!        step: length of steepest descent segment
!       ispin: The spin index 1,2 for up down
!
! OUTPUT: (through common)
!   occnumtmp: new occupation numbers
!   solve_occ_th_n: sum_occ - float(nele(ispin))
!---------------------------------------------------------------------

!..Global
   use global; use functional_m; use energies; use orbocc; use vary_occ_par
   use func_der_args
   implicit none

!..The function
   real(dp) :: solve_occ_th_n

!..Arguments
   real(dp), intent(inout) :: xmu
   real(dp), intent(inout) :: occnumtmp(lnnatorb,3)

!..Local variables
   real(dp), allocatable :: DE_Dn(:), DE_Dtheta(:)
   real(dp), allocatable :: delta_occ(:)
   real(dp) :: stpa, stpb, stp1, stplim, grada, gradb
   real(dp) :: fac,stpdif, sum_occ, hnorm, sum_ab_var
   real(dp) :: arg, smu
   real(dp), external :: grad_th_n, zbisect, zbrent_gr
   integer :: icount, ia, iter, iterb, isp, nfrac
   logical :: crit
   real(dp) :: enum, denom, dndthsq, fr

!..Local array allocation
   allocate ( DE_Dn(lnnatorb), DE_Dtheta(lnnatorb), delta_occ(lnnatorb) )
   xmu1 = xmu
   ispin1 = ispin
   ispin1_op = ispin_op
   fac = 0.1d0
   stplim = 0.01d0
   
   occnumtmp(:,ispin)=occnum(:,ispin)
   where(occnumtmp(:,ispin) < smallocc) occnumtmp(:,ispin)=smallocc
   where(occnumtmp(:,ispin) > 1.d0-smallocc) &
         occnumtmp(:,ispin) = 1.d0-smallocc
   occnumtmp(:,3)=occnumtmp(:,1) + occnumtmp(:,2)

!..Initialize the direction 
   call Func_der_n_n(occnumtmp, DE_Dn, ispin)

!..Project out the perpendicular component:
   enum=0.d0
   denom=0.d0
   do ia=1,nnatorb
      arg = sqrt(occnumtmp(ia,ispin))
      if(abs(arg) > 1.d0) print*,'arg1=',arg
      theta1(ia) = tovpi*asin(arg) 

      dndthsq=(0.5d0*pi*sin(pi*theta1(ia)))**2
      enum=enum+DE_Dn(ia)*dndthsq
      denom=denom+dndthsq
      DE_Dtheta(ia) = DE_Dn(ia)*piovtwo*sin(pi*theta1(ia))
   enddo
   fr=enum/denom

   h=DE_Dtheta-fr*0.5d0*pi*sin(pi*theta1)
   hnorm=1.d0/sqrt(dot_product(h,h))

   h=hnorm*h

!..h is the direction in theta space that I minimize along

   do iter = 1, niter_on !loop on the directions (segments)
      stpa = -zero
      stpb =  step
      stpdif = stpb - stpa
      icount = 0

!.....Bracket the minmum on the direction h
      grada = grad_th_n(stpa,occnumtmp)
      gradb = grad_th_n(stpb,occnumtmp)

  10  continue

      if(grada*gradb > 0.d0) then

!........Geometric increase of the window
         if(abs(grada) >= abs(gradb)) then
            stpb = stpb + fac*stpdif
            gradb = grad_th_n(stpb,occnumtmp)
         else
            stpa = stpa - fac*stpdif
            grada = grad_th_n(stpa,occnumtmp)
         endif
         stpdif = stpb - stpa
         if(stpdif > stplim) then
            stp1=stplim
            goto 30
         endif

         icount = icount + 1
         if(icount < nbracket) then
            go to 10
         else
            stp1=stplim
            go to 30
!           stop 'vary_occ: Line minimization problem plot debug'
         endif
      endif !grada*gradb > 0.d0
  
!.....Bisection rule
      stp1 = zbisect(grad_th_n,stpa,stpb,niter_on1,crit_on1,iterb, &
                     occnumtmp)

!.....Brent method:
!     stp1 = zbrent_gr(grad_th_n,stpa,stpb,niter_on1,crit_on1,iterb,&
!                      occnumtmp)

 30   continue

      grada=grad_th_n(stp1,occnumtmp)

!.....The average absolute change in occupation numbers:
      nfrac=0
      smu=0.d0
      do ia=1,nnatorb
         if(occnum1(ia,ispin) < 0.99d0) then
            if(occnum1(ia,ispin) > smallocc) then
               smu=smu+DE_Dn1(ia)
               nfrac=nfrac+1
            endif
         endif
      enddo
      delta_occ=abs(occnum1(:,ispin)-occnumtmp(:,ispin))
      sum_ab_var = sum(delta_occ)/float(nnatorb)
      xmu=smu/float(nfrac)

!.....Redefine direction h:
      enum=0.d0
      denom=0.d0
      DE_Dn=DE_Dn1
      DE_Dtheta=DE_Dtheta1
      occnumtmp(:,ispin) = occnum1(:,ispin)
      if(cl_shell) occnumtmp(:,2) = occnum1(:,2)
      occnumtmp(:,3) = occnum1(:,3)
      do ia =1,nnatorb
         dndthsq=(0.5d0*pi*sin(pi*theta1(ia)))**2
         enum=enum+DE_Dn(ia)*dndthsq
         denom=denom+dndthsq
      enddo
      fr=enum/denom

      h = DE_Dtheta -fr*0.5d0*pi*sin(pi*theta1)
      hnorm=1.d0/sqrt(dot_product(h,h))
      h = hnorm*h

      crit=sum_ab_var < crit_on

      if(crit) then
         goto 14
      endif

   enddo !iter
14 continue
 
   if((.not.cl_shell).and.common_mu) then 
      isp=3
   else
      isp=ispin
   endif
   sum_occ= sum(occnumtmp(:,isp))
   solve_occ_th_n = sum_occ - xnele(isp)

   deallocate ( DE_Dn, DE_Dtheta, delta_occ )

end function solve_occ_th_n
!-----------------------------------------------------------------------

subroutine Func_der_n_n(occ, DE_Dn, ispin)
!-----------------------------------------------------------------------
! This subroutine calculates the derivatives with respect
! to the occ numbers. 
! Input: 
!     occ: occupation numbers
!   ispin: the spin index
! OUTPUT:
!       DE_Dn : Derivative of energy w.r.t. the occ. numbers.
!-----------------------------------------------------------------------

!..Global
   use global; use functional_m; use matrices
   implicit none

!..Arguments
   integer, intent(in) :: ispin
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: DE_Dn(lnnatorb)

!..Local variables
   integer :: ia, ib
   real(dp), allocatable :: fac_c_d(:), fac_e_d(:)
   real(dp), external :: TS_nlnn_d

!..Local array allocation
   allocate (fac_c_d(lnnatorb), fac_e_d(lnnatorb))

   do ia=1,nnatorb
      DE_Dn(ia) = HcoreNO(ia,ia) 
      call calc_der_occ_factors(ia, ispin, occ, fac_c_d, fac_e_d)

!.....Now calculate the derivative DE/Dn_a
      do ib = 1, nnatorb
         DE_Dn(ia) = DE_Dn(ia) + fac_c_d(ib)*CoulNO(ib,ia) &
                               + fac_e_d(ib)*ExchNO(ib,ia)
      enddo
   enddo

!..Add the entropy term:
   if (temp_dep) then
      do ia=1,nnatorb
         DE_Dn(ia) = DE_Dn(ia) + TS_nlnn_d(ia,ispin)     
      enddo
   endif

   deallocate (fac_c_d, fac_e_d)

end subroutine Func_der_n_n

!----------------------------------------------------------------------
function grad_th_n(stp,occnumtmp)
! This subroutine is just a wrapper so we have the function
! grad_th_n with only one argument to be used in root finding method
! occnum1, De_Dn1, sum_occ are  by-product outputs
! needed only when convergence is achieved

!..Global
   use global; use functional_m; use func_der_args
   implicit none
   
!..The function
   real(dp) :: grad_th_n

!..Arguments
   real(dp), intent(in) :: stp
   real(dp), intent(in) :: occnumtmp(lnnatorb,3)

!..Local variables
   real(dp) :: xnel

   theta1 = tovpi*asin(sqrt(occnumtmp(:,ispin1)))

!..Vary the occ numbers (thetas) on h direction by stp
   theta1 = theta1 - stp*h
   occnum1(:,ispin1) = sin(piovtwo*theta1)**2
   if(cl_shell) then
      occnum1(:,ispin1_op)=occnum1(:,ispin1)
   else
      occnum1(:,ispin1_op)=occnumtmp(:,ispin1_op)
   endif
   occnum1(:,3) = occnum1(:,ispin1) + occnum1(:,ispin1_op)

!..Impose the constraint again:
   xnel=xnele(ispin1)
   call impose_constr(occnum1,theta1,ispin1,ispin1_op,xnel)

!..Recalculate the derivative:
   call Func_der_n_n(occnum1, DE_Dn1, ispin1)

   DE_Dtheta1 = piovtwo*DE_Dn1*sin(pi*theta1)

   grad_th_n = dot_product(DE_Dtheta1, h)/sqrt(dot_product(h,h))

end function grad_th_n
!-----------------------------------------------------------------------

subroutine impose_constr(occnum1,theta1,ispin,ispin_op,xnel)
!-----------------------------------------------------------------------
! Imposes the number of particle conservation constraint by multiplying
! the fractional occupancies by an appropriate factor.
! occnum1: occupation to be corrected
! theta1: the corresponding theta variables 
! ispin: the spin index
! ispin_op: the oposite spin index
! xnel: the number of electrons for ispin
!----------------------------------------------------------------------


!..Global
   use global; use functional_m
   implicit none

!..Argumens
   integer, intent(in) :: ispin, ispin_op
   real(dp), intent(inout) :: theta1(lnnatorb), occnum1(lnnatorb,3)
   real(dp), intent(in) ::  xnel

!..Local
   real(dp) :: sum1, sum2, xnel2, arg, f, sumnew
   integer :: ia
   
!..Sum all fractional occupancies
   sum1=0.d0
   sum2=0.d0
   xnel2=xnel
   do ia=1,nnatorb
      sum1=sum1+occnum1(ia,ispin)
      if(occnum1(ia,ispin) < 1.d0-small) then
         sum2=sum2+occnum1(ia,ispin)
      else
         xnel2=xnel2-1.d0 
      endif
   enddo

   if(abs(xnel-sum1) > small) then
      if(xnel > sum1) then
         f=(xnel-sum1)/nnatorb
         do ia=1,nnatorb
            occnum1(ia,ispin)=occnum1(ia,ispin)+f
            if(occnum1(ia,ispin) > 1.d0) then
               if(ia /= nnatorb) then
                  f=f+(occnum1(ia,ispin)-1.d0)/dfloat(nnatorb-ia)
               endif
               occnum1(ia,ispin)=one
            endif
            arg = sqrt(occnum1(ia,ispin))
            theta1(ia) = tovpi*asin(arg)
         enddo
      else
         f=xnel2/sum2
         do ia=1,nnatorb
            if(occnum1(ia,ispin) < 1.d0-small) then
               occnum1(ia,ispin)=occnum1(ia,ispin)*f
               if(occnum1(ia,ispin) < zero) occnum1(ia,ispin)=zero
            endif
            arg = sqrt(occnum1(ia,ispin))
            theta1(ia) = tovpi*asin(arg)
         enddo
      endif
   endif

   sumnew=sum(occnum1(:,ispin))

   if(cl_shell) occnum1(:,ispin_op)=  occnum1(:,ispin)
   occnum1(:,3) = occnum1(:,ispin) + occnum1(:,ispin_op)

end subroutine impose_constr
!----------------------------------------------------------------------
