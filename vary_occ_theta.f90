subroutine vary_occ_th(xmu, xmu_range)
!---------------------------------------------------------------------
! This subroutine finds the correct value of mu that gives the correct
! number of electrons. The occupation numbers for each mu are optimized
! using the solve_occ routine below.
! The substitution n=sin^2 (theta*pi/2) is used to take care
! of the constraint 0<=n<=1
!
! INPUT:
!   xmu(2)              : on input contains the initial values of \mu
!                         for spin_up and spin_down
!   xmu_range           : a range for the search of \mu. 
!
! In Module orbocc:
!   occnum(lnnatorb)    : on input contains the starting occ numbers
!   occnum(lnnatorb)    : on output the refined occ numbers
!
! OUTPUT:
!   xmu(2)              : on output the refined value of \mu
! 
!         created by N. N. Lathiotakis and I.Theophilou
!
!----------------------------------------------------------------------

!..Global
   use global; use functional_m; use orbocc; use vary_occ_par
   implicit none

!..Arguments
   real(dp), intent(in) :: xmu_range
   real(dp), intent(inout) :: xmu(2)

!..Local variables
   real(dp), external :: solve_occ_th
   real(dp), external :: zbisect, zbrent_so
   real(dp) :: xmui, xmu_a, xmu_b, xmu_stp, sum_a, sum_b
   real(dp) :: dif, sum_occ(3)
   real(dp), allocatable :: occnumtmp(:,:)
   integer :: icount, ia, iterb, nspin, isc, nsc
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
      nsc=2
   endif
      
!..Initialize the occnumtmp. Smallocc should be no less than 1d-8
!..for numerical stability
   occnumtmp(:,:2) = occnum(:,:2)
   where(occnumtmp(:,:2) < smallocc) occnumtmp(:,:2)=smallocc
   where(occnumtmp(:,:2) > 1.d0-smallocc) occnumtmp(:,:2)=1.d0-smallocc
   occnumtmp(:,3) = occnumtmp(:,1) +  occnumtmp(:,2)

!..Up and down occ nums are varied independently/alternatively nsc times
   do isc=1,nsc 
      print*,'================================================================='
      print*,'ISC=',isc
      do ispin = 1, nspin
!  if(nele(ispin) /= 0) then ! If you want to run He triplet uncomment

!........Define the opposite spin
         print*,'-----------------------------------------------------------------'
         print*,'SPIN: ',spin_nam(ispin),' :'
         if(ispin == 1) then
            ispin_op = 2
         else
            ispin_op = 1
         endif
 
         if(common_mu.and.ispin == 2) then
            xmu(2)=xmu(1)
            goto 111
         endif

         xmui = xmu(ispin)

         xmu_stp=min(max(xmu_range * abs(xmui),0.01d0),0.2d0)
         xmu_a = xmui - xmu_stp
         xmu_b = xmui + xmu_stp

         icount=0

!........Calculate occ. nums at the extremal mu values
         sum_a = solve_occ_th(xmu_a,occnumtmp)
         sum_b = solve_occ_th(xmu_b,occnumtmp)

!........Check if variation is necessary.
         if(abs(sum_a) < crit_mu.and.abs(sum_b) < crit_mu) then
            print*,'crit_mu satisfied within the entire range'
            xmu(ispin) = 0.5d0*(xmu_a + xmu_b)
            goto 111
         endif

!........Find a mu range (xmu_a, xmu_b) that brackets the correct mu
  5      continue
         print*,'(',xmu_a,'< mu < ',xmu_b,')'

         if(sum_a*sum_b >= 0.d0) then
            print*,'For mu=',xmu_a,'sum_a=',sum_a
            print*,'For mu=',xmu_b,'sum_a=',sum_b
            print*,'Redefining xmu, xmu_a, xmu_b ...'
            if(sum_b <= 0.d0) then
               dif   = xmu_b - xmu_a
               xmu_a = xmu_b
               sum_a = sum_b 
               xmu_b = xmu_a+dif
               sum_b = solve_occ_th(xmu_b,occnumtmp)
            elseif(sum_a >= 0.d0) then
               dif   = xmu_b - xmu_a
               xmu_b = xmu_a
               sum_b = sum_a
               xmu_a = xmu_b-dif
               sum_a = solve_occ_th(xmu_a,occnumtmp)
            endif
            icount= icount+1
            if(icount < nbracket) then
               goto 5
            else
               print*, 'Mu can not be bracketed'
               go to 11
               stop 'vary_occ_theta:Try to increase the range of mu'
            endif
         else
!           print*,'For mu=',xmu_a,'sum_a=',sum_a
!           print*,'For mu=',xmu_b,'sum_a=',sum_b
            print*,'Mu Bracketed!!!'
         endif
!........END: Find a mu range (xmu_a, xmu_b) that brackets the correct mu

         sum_a=solve_occ_th(xmu_a,occnumtmp)
         sum_b=solve_occ_th(xmu_b,occnumtmp)

!........Find xmu with Bisection method
!        xmui = zbisect(solve_occ_th, xmu_a, xmu_b, nmu, crit_mu,iterb, &
!                       occnumtmp)

!........Find xmu with Brent method
         xmui = zbrent_so(solve_occ_th,xmu_a, xmu_b, nmu, crit_mu, iterb,&
                          occnumtmp)
         xmu(ispin) = xmui

 111     continue
         sum_a = solve_occ_th(xmu(ispin),occnumtmp)
!        print*,'sum occ - nel for optimal mu:',sum_a

         occnum(:,ispin) = occnumtmp(:,ispin)
         occnum(:,3) = occnumtmp(:,3)
         sum_occ(ispin) = sum(occnum(:,ispin))

  11     continue
         print*,'Optimal mu:',xmui,' Nu of mu iter.:',iterb
         print*,'New occ numbers:'
         write(6,'("NEWOCC",5f12.8)') (occnum(ia,ispin), ia=1,5)
         write(6,'("      ",5f12.8)') (occnum(ia,ispin), ia=6,nnatorb)
         print*,'Sum of occ numbers:', sum_occ(ispin)

!        endif !nele(ispin) /= 0
      enddo !ispin = 1,nspin
   enddo !isc

   if(cl_shell) then
      xmu(2) = xmu(1)
      sum_occ(2) = sum_occ(1)
      occnum(:,2) = occnum(:,1)
      print*,'For closed shell singlet spin DOWN quantities'
      print*,'are the same as spin UP.'
   endif
   sum_occ(3) = sum_occ(1) + sum_occ(2)
   print*,'------------------------------------------------------'

!..After all occ. numbers are redefined: 
   call calc_occ_factors(occnum, fac_h, fac_c, fac_e)

   deallocate (occnumtmp)


end subroutine vary_occ_th

!-----------------------------------------------------------------------

function solve_occ_th(xmu,occnumtmp)
!------------------------------------------------------------------------
! This subroutine solves the minimization problem with respect to the 
! occ numbers for fixed mu (input). 
!
! On output solve_occ_th = sum of the occupation numbers *  - nele(ispin)
! for the right mu and occnum it should be zero.
!
!
! INPUT:
!         xmu:  \mu
!
! Also from the vary_occ_par module:
!      occnum: Old occupation numbers
!    niter_on: Max. number of directions in conjugate gradient method
!    nbracket: max. number of bracketing steps in each direction
!   niteron_1: max. number of divisions in each direction
!     crit_on:   conv criterion for occ numbers (all derivatives must be
!               smaller than this number).
!    crit_on1:  linear search conv. criterion
!        step: length of steepest descent segment
!       ispin: The spin index 1,2 for up down
!
! OUTPUT: 
!   occnumtmp: new occupation numbers
!   solve_occ_th: sum_occ - float(nele(ispin))
!---------------------------------------------------------------------------

!..Global
   use global; use functional_m; use energies; use orbocc; use vary_occ_par
   use func_der_args
   implicit none

!..The function
   real(dp) :: solve_occ_th

!..Arguments
   real(dp), intent(in) :: xmu
   real(dp), intent(out) :: occnumtmp(lnnatorb,3)

!..Local variables
   real(dp), allocatable :: DE_Dn(:), delta_occ(:)
   real(dp), allocatable :: DE_Dtheta(:), DE_Dtheta_old(:)
   real(dp) :: stpa, stpb, stp1, stplim
   real(dp) :: grada, gradb, prod, prod_old, gamma_f, prod1
   real(dp) :: fac,stpdif, sum_occ, hnorm, sum_ab_var
   real(dp), external :: grad_th, zbisect, zbrent_gr
   integer :: icount, iter, iterb, isp
   logical :: crit

!..Local array allocation
   allocate ( DE_Dn(lnnatorb), delta_occ(lnnatorb),&
              DE_Dtheta(lnnatorb), DE_Dtheta_old(lnnatorb) )

   xmu1 = xmu
   ispin1 = ispin
   ispin1_op = ispin_op
   fac = 0.1d0
   stplim = 5.d2
   
   occnumtmp(:,ispin) = occnum(:,ispin)
   where(occnumtmp(:,ispin) < smallocc) occnumtmp(:,ispin)=smallocc
   where(occnumtmp(:,ispin) > 1.d0-smallocc) &
         occnumtmp(:,ispin) = 1.d0-smallocc
   occnumtmp(:,3) = occnumtmp(:,1) + occnumtmp(:,2)

!..Initialize the direction 
   call Func_der_n(xmu, occnumtmp, DE_Dn, ispin)

   theta1=tovpi*asin(sqrt(occnumtmp(:,ispin)))
   DE_Dtheta = piovtwo*DE_Dn*sin(pi*theta1)

   h=DE_Dtheta
   hnorm=1.d0/sqrt(dot_product(h,h))
   h=hnorm*h

   prod=dot_product(DE_Dtheta, DE_Dtheta)
      
   do iter = 1, niter_on !loop on the directions (segments)
      prod_old = prod !For the Flecher-Reeves gamma
      DE_Dtheta_old = DE_Dtheta !For the Pollak-Rivere gamma
      stpa = -zero
      stpb =  step
      icount = 0

!.....Bracket the minmum on the direction h
      grada = grad_th(stpa,occnumtmp)
      gradb = grad_th(stpb,occnumtmp)

  10  continue

      if(grada*gradb > 0.d0) then
         if(icount == 1.and.abs(grada-gradb) > small) then
            stp1=(grada*stpb-stpa*gradb)/(grada-gradb)
            if(stp1 > stpb) then
               stpb=stp1
               gradb=grad_th(stp1,occnumtmp)
            elseif(stp1 < stpa) then
               stpa=stp1
               grada=grad_th(stp1,occnumtmp)
            endif
         endif
!........Linear increase of the window
!        if(abs(grada) > abs(gradb)) then
!           stpb = stpb + step
!        else
!           stpa = stpa - step
!        endif

!........Geometric increase of the window
         stpdif = stpb - stpa
         if(abs(grada) >= abs(gradb)) then
            stpb = stpb + fac*stpdif
            gradb = grad_th(stpb,occnumtmp)
         else
            stpa = stpa - fac*stpdif
            grada = grad_th(stpa,occnumtmp)
         endif

         icount = icount + 1
         if(icount < nbracket) then
            go to 10
         else
            go to 40 !Line minimization problem: abort
         endif
      endif !grada*gradb > 0.d0

!.....Bisection rule
!     stp1 = zbisect(grad_th,stpa,stpb,niter_on1,crit_on1,iterb,occnumtmp)

!.....Brent method:
      stp1 = zbrent_gr(grad_th,stpa,stpb,niter_on1,crit_on1,iterb,occnumtmp)

      if(stp1 > stplim) then
         stp1=stplim
         grada=grad_th(stp1,occnumtmp)
      endif

!.....The average absolute change in occupation numbers:
      delta_occ = abs(occnum1(:,ispin)-occnumtmp(:,ispin))
!     where (delta_occ < 0.d0) delta_occ=-delta_occ
      sum_ab_var= sum(delta_occ)/float(nnatorb)

      DE_Dn=DE_Dn1
      DE_Dtheta=DE_Dtheta1
      occnumtmp(:,ispin) = occnum1(:,ispin)
      if(cl_shell) occnumtmp(:,2) = occnum1(:,2)
      occnumtmp(:,3) = occnum1(:,3)

 40   continue

!.....Redefine gamma and h (see numer. recipes sect 10.6)
      select case (n_vary_meth)
!........Steepest descent: 
         case (1)
            gamma_f=0.d0

!........Conjugate gradient with
!........Fletcher-Reeves formula for gamma:
         case (2)
            prod = dot_product(DE_Dtheta, DE_Dtheta)
            if(iter == 1) prod_old=1d3
            gamma_f = prod / max(prod_old, small)
            prod_old=prod

!........Conjugate gradient with
!........Polak-Ribiere formula for gamma:
         case (3)
            prod = dot_product(De_Dtheta-De_Dtheta_old, De_Dtheta)
            prod1= dot_product(De_Dtheta_old, De_Dtheta_old)
            gamma_f = prod/max(prod1,small)

!........Unkknown method!
         case default
            stop 'vary_occ_theta: unknown variation method!'
      end select
        
      h = DE_Dtheta + gamma_f*h
      hnorm = 1.d0/sqrt(dot_product(h,h))
      h= hnorm*h

!....CONVERGENCE CRITERIA:
!....(1)The derivative equal to zero:
      crit=.true.
!     do ia=1, nnatorb
!        crit=crit.and.(abs(DE_Dtheta(ia)) < crit_on)
!     enddo

!....(2)The change in occ numbers is small:
      crit=crit.and.(sum_ab_var < crit_on)

      if(crit) then
         exit
      endif
   enddo !iter
    
   if((.not.cl_shell).and.common_mu) then 
      isp=3
   else
      isp=ispin
   endif
   sum_occ = sum(occnumtmp(:,isp))
   solve_occ_th = sum_occ - xnele(isp)

   deallocate ( DE_Dn, delta_occ, DE_Dtheta, DE_Dtheta_old )

end function solve_occ_th
!-----------------------------------------------------------------------

subroutine Func_der_n(xmu, occ, DE_Dn, ispin)
!-----------------------------------------------------------------------
! This subroutine calculates the derivatives with respect
! to the occ numbers. The constraint -2\mu (\sum n_i - N/2) is added 
! also.
! INPUT: 
!    xmu : \mu
!    occ : occupation numbers
!  ispin : spin index
! 
! OUTPUT:
!       DE_Dn : Derivative of energy w.r.t. the occ. numbers.
!-----------------------------------------------------------------------

!..Global
   use global; use functional_m; use matrices 
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3), xmu
   real(dp), intent(out) :: DE_Dn(lnnatorb)
   integer, intent(in) :: ispin

!..Local variables
   integer :: ia, ib
   real(dp), allocatable :: fac_c_d(:), fac_e_d(:)

   real(dp) :: TS_nlnn_d

!..Local array allocation
   allocate( fac_c_d(lnnatorb), fac_e_d(lnnatorb) )

   do ia=1,nnatorb
      DE_Dn(ia) = HcoreNO(ia,ia) - xmu
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

   deallocate( fac_c_d, fac_e_d )

end subroutine Func_der_n

!-----------------------------------------------------------------------

function grad_th(stp,occnumtmp)
!-----------------------------------------------------------------------
! This subroutine is just a wrapper so we have the function
! grad_th with only one argument to be used in root finding method
! occnum1, De_Dn1, sum_occ are  by-product outputs
! needed only when convergence is achieved
!-----------------------------------------------------------------------

!..Global
   use global; use functional_m; use func_der_args      
   implicit none

!..The function
   real(dp) :: grad_th

!..Arguments
   real(dp), intent(in) :: stp
   real(dp), intent(in) :: occnumtmp(lnnatorb,3)

   theta1=tovpi*asin(sqrt(occnumtmp(:,ispin1)))
   
!..Vary along the direction
   theta1 = theta1- stp*h
   occnum1(:,ispin1) = sin(piovtwo*theta1)**2
   occnum1(:,ispin1) = max(min(occnum1(:,ispin1), one), small)
   if(cl_shell) then
      occnum1(:,ispin1_op) = occnum1(:,ispin1)
   else
      occnum1(:,ispin1_op) = occnumtmp(:,ispin1_op)
   endif
   occnum1(:,3) = occnum1(:,ispin1) + occnum1(:,ispin1_op)

   call Func_der_n(xmu1, occnum1, DE_Dn1, ispin1)
   DE_Dtheta1 = piovtwo*DE_Dn1*sin(pi*theta1)

   grad_th = dot_product(DE_Dtheta1, h)/sqrt(dot_product(h,h))

end function grad_th
