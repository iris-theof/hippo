subroutine vary_occ(xmu, xmu_range)
!----------------------------------------------------------------------
! This subroutine finds the correct value of mu that gives the correct
! number of electrons. For each mu the occ numbers are optimized using
! the solve_occ function.
!
! INPUT:
!   xmu(2)              : on input contains the initial values of \mu
!                         for spin_up and spin_down
!   xmu_range           : a range for the search of \mu. 
!
! In vary_occ_par module:
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
   real(dp), external :: solve_occ
   real(dp), external :: zbisect, zbrent_so
   real(dp), allocatable :: occnumtmp(:,:)
   real(dp) :: xmui, xmu_a, xmu_b, xmu_stp, sum_a, sum_b
   real(dp) :: dif, sum_occ(3)
   integer :: icount, ia, iterb
   integer :: nspin, isc, nsc
   character(4) spin_nam(2)

!..Local array allocation
   allocate( occnumtmp(1:lnnatorb,1:3))
   spin_nam(1) = ' UP '
   spin_nam(2) = 'DOWN'

   if(cl_shell) then
      nsc=1
      nspin = 1
   else
      nspin = 2
      nsc=10
   endif
   
!..Initialize the occnumtmp. Smallocc should be no less than 1d-8
!..for numerical stability
   occnumtmp(:,:2) = occnum(:,:2)
   where(occnumtmp(:,:2) < smallocc) occnumtmp(:,:2)=smallocc
   where(occnumtmp(:,:2) > 1.d0-smallocc) occnumtmp(:,:2)=1.d0-smallocc
   occnumtmp(:,3) = occnumtmp(:,1) +  occnumtmp(:,2)

!..Up and down occ nums are varied independently/alternatively nsc times
   do isc=1,nsc
      print*,'======================================'
      print*,'isc=',isc
      do ispin = 1, nspin
!        if(nele(ispin).ne.0) then !To run He triplet for instance
!........Define the opposite spin
         if(ispin.eq.1) then
            ispin_op = 2
         else
            ispin_op = 1
         endif

         if(common_mu.and.ispin.eq.2) then
            xmu(2)=xmu(1)
            goto 111
         endif

         xmui = xmu(ispin)

         xmu_stp=min(max(xmu_range * abs(xmui),0.01d0),0.2d0)
         xmu_a = xmui - xmu_stp
         xmu_b = xmui + xmu_stp

         icount=0

!.......Calculate occ. nums at the extremal mu values
         sum_a = solve_occ(xmu_a,occnumtmp)
         sum_b = solve_occ(xmu_b,occnumtmp)

!........Check if variation is necessary.
         if(abs(sum_a) < crit_mu.and.abs(sum_b) < crit_mu) then
            print*,'crit_mu satisfied within the entire range'
            xmu(ispin) = 0.5d0*(xmu_a + xmu_b)
            goto 111
         endif

!.......Find a mu range (xmu_a, xmu_b) that brackets the correct mu
  5      continue
         print*,'(',xmu_a,'< mu < ',xmu_b,')'

         if(sum_a*sum_b.ge.0.d0) then
            print*,'For mu=',xmu_a,'sum_a=',sum_a
            print*,'For mu=',xmu_b,'sum_a=',sum_b
            print*,'Redefining xmu, xmu_a, xmu_b ...'
            if(sum_b.le.0.d0) then
               dif   = xmu_b - xmu_a
               xmu_a = xmu_b
               sum_a = sum_b 
               xmu_b = xmu_a+dif
               sum_b = solve_occ(xmu_b,occnumtmp)
            elseif(sum_a.ge.0.d0) then
               dif   = xmu_b - xmu_a
               xmu_b = xmu_a
               sum_b = sum_a
               xmu_a = xmu_b-dif
               sum_a = solve_occ(xmu_a,occnumtmp)
            endif
            icount= icount+1
            if(icount.lt.200) then
               goto 5
            else
               print*, 'Mu can not be bracketed'
               go to 11
               stop 'vary_occ:Try to increase the range of mu'
            endif
         else
            print*,'For mu=',xmu_a,'sum_a=',sum_a
            print*,'For mu=',xmu_b,'sum_a=',sum_b
            print*,'Mu Bracketed!!!'
         endif
!.....END: Find a mu range (xmu_a, xmu_b) that brackets the correct mu

         sum_a=solve_occ(xmu_a,occnumtmp)
         sum_b=solve_occ(xmu_b,occnumtmp)

!.....Find xmu with Bisection method
!        xmui = zbisect(solve_occ, xmu_a, xmu_b, nmu, crit_mu,iterb,&
!                       occnumtmp)

!.....Find xmu with Brent method
         xmui = zbrent_so(solve_occ, xmu_a, xmu_b, nmu, crit_mu,iterb,&
                          occnumtmp)
         xmu(ispin) = xmui

 111     continue
         sum_a = solve_occ(xmu(ispin),occnumtmp)
         print*,'sum occ - nel for optimal mu:',sum_a

         occnum(:,ispin) = occnumtmp(:,ispin)
         occnum(:,3) = occnumtmp(:,3)
         sum_occ(ispin) = sum(occnum(:,ispin))

  11     continue
         print*,'------------------------------------------------------'
         print*,'SPIN: ',spin_nam(ispin),' :'
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

!..After all occ. numbers are redefined call 
   call calc_occ_factors(occnum, fac_h, fac_c, fac_e)

   deallocate( occnumtmp)

end subroutine vary_occ

!----------------------------------------------------------------------

function solve_occ(xmu, occnumtmp)
!----------------------------------------------------------------------
! This subroutine solves the minimization problem with respect to the 
! occ numbers for fixed mu (input). 
!
! On output solve_occ = sum of the occupation numbers * nele
! for the right mu and occnum it should be zero.
!
!
! INPUT:
!         xmu:  \mu
!
! Also from module vary_occ_par:
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
! OUTPUT: (through common)
!   occnumtmp: new occupation numbers
!   solve_occ: sum_occ - float(nele(ispin))
!------------------------------------------------------------------------

!..Global
   use global; use functional_m; use energies; use orbocc; use vary_occ_par
   use func_der_args
   implicit none

!..The function
   real(dp) solve_occ

!..Arguments
   real(dp), intent(in) :: xmu
   real(dp), intent(inout) :: occnumtmp(lnnatorb,3)

!..Local variables
   real(dp), allocatable :: DE_Dn(:), DE_Dn_old(:)
   real(dp) :: Delta_occ(lnnatorb)
   real(dp) :: stpa, stpb, stp1, stplim
   real(dp) :: grada, gradb
   real(dp) :: prod, prod_old, gamma_f, prod1
   real(dp) :: fac,stpdif, sum_occ, hnorm, sum_ab_var
   real(dp) :: En_con_old
   real(dp), external :: grad, zbisect, zbrent_gr
   integer :: icount, ia, iter, iterb, isp
   logical :: crit

!..Local array allocation
   allocate ( DE_Dn(lnnatorb), DE_Dn_old(lnnatorb) )

   nstart = 1 !The lower index of the occ. numbers to vary
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

   h = DE_Dn
   hnorm=1.d0/sqrt(dot_product(h,h))
   h=hnorm*h

!..Define lower index of occ. numbers to vary
   do ia=1,nnatorb
      occ_var(ia)=.true.
      if(occnumtmp(ia,ispin).ge.one-small.and.h(ia).lt.der_min)  then
         occ_var(ia)=.false.
      endif
   enddo

   nstart=1
   do ia=1,nnatorb
      if(occ_var(ia)) then
         nstart=nstart+1
         goto 99
      endif
   enddo
99 continue
   if(nstart.gt.nnatorb) stop 'solve_occ: nstart.ge.nnatorb'

   prod=dot_product(DE_Dn, DE_Dn)
      
   En_con_old=0.d0

   do iter = 1, niter_on !loop on the directions (segments)
      prod_old = prod !For the Flecher-Reeves gamma
      DE_Dn_old=DE_Dn
      stpa = -zero
      stpb =  step
      icount = 0

!..Bracket the minmum on the direction h
      grada = grad(stpa,occnumtmp)
      gradb = grad(stpb,occnumtmp)

  10  continue

      if(grada*gradb.gt.0.d0) then
         if(icount.eq.1.and.abs(grada-gradb).gt.1.d-12) then
            stp1=(grada*stpb-stpa*gradb)/(grada-gradb)
            if(stp1.gt.stpb) then
               stpb=stp1
               gradb=grad(stp1,occnumtmp)
            elseif(stp1.lt.stpa) then
               stpa=stp1
               grada=grad(stp1,occnumtmp)
            endif
         endif
!..Linear increase of the window
!        if(abs(grada).gt.abs(gradb)) then
!           stpb = stpb + step
!        else
!           stpa = stpa - step
!        endif

!..Geometric increase of the window
         stpdif = stpb - stpa
         if(abs(grada).ge.abs(gradb)) then
            stpb = stpb + fac*stpdif
            gradb = grad(stpb,occnumtmp)
         else
            stpa = stpa - fac*stpdif
            grada = grad(stpa,occnumtmp)
         endif

         icount = icount + 1
         if(icount.lt.nbracket) then
            go to 10
         else
            prod_old=1.d3*prod_old
            go to 40 !Line minimization problem: abort
         endif
      endif !grada*gradb.gt.0.d0

!..Bisection rule
!     stp1 = zbisect(grad,xa,xb,niter_on1,crit_on1,iterb,occnumtmp)

!..Brent method:
      stp1 = zbrent_gr(grad,stpa,stpb,niter_on1,crit_on1,iterb,occnumtmp)

  40  continue

      if(stp1.gt.stplim) then
         stp1=stplim
         grada=grad(stp1,occnumtmp)
      endif

! The average absolute change in occupation numbers:
      delta_occ=abs(occnum1(:,ispin)-occnumtmp(:,ispin))
!     where (delta_occ < 0.d0) delta_occ=-delta_occ
      sum_ab_var = sum(delta_occ)/float(nnatorb)

      DE_Dn = DE_Dn1
      occnumtmp(:,ispin) = occnum1(:,ispin)
      if(cl_shell) occnumtmp(:,2) = occnum1(:,2)
      occnumtmp(:,3) = occnum1(:,3)

! 40  continue

!..Redefine gamma and h (see numer. recipes sect 10.6)
      select case(n_vary_meth)
!........Steepest descent: 
         case (1) 
            gamma_f=0.d0

!........Conjugate gradient with Fletcher-Reeves formula for gamma
         case (2) 
            prod = dot_product(DE_Dn,DE_Dn)
            if(iter == 1) prod_old=1.d3
            gamma_f = prod / max(prod_old, small)
            prod_old=prod

!........Conjugate gradient with Polak-Ribiere formula for gamma
         case (3) 
            prod = dot_product(De_Dn-De_Dn_old, De_Dn)
            prod1 = dot_product(De_Dn_old, De_Dn_old)
            gamma_f = prod/max(prod1,small)

!........Unkknown method!
         case default 
            stop 'vary_occ: unknown variation method!'
      end select
        
      h=DE_Dn + gamma_f * h
      hnorm=1.d0/sqrt(dot_product(h,h))
      h=hnorm*h

      do ia=1,nnatorb
         occ_var(ia)=.true.
         if(occnumtmp(ia,ispin).ge.one-small.and.h(ia).lt.der_min) then
            occ_var(ia)=.false.
         endif
      enddo
      nstart=1
      do ia=1,nnatorb
         if(occ_var(ia)) then
            nstart=nstart+1
            goto 101
         endif
      enddo
 101  continue
      if(nstart.gt.nnatorb) stop 'solve_occ: nstart.ge.nnatorb'

! CONVERGENCE CRITERIA:
! (1)The derivative equal to zero:
!     crit=.true.
!     do ia=1, nnatorb
!        if(occ_var(ia)) then
!           crit=crit.and.(abs(DE_Dn(ia)).lt.crit_on)
!        endif
!     enddo

! (2)The change in occ numbers is small:
      crit=sum_ab_var.lt.crit_on

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
   solve_occ = sum_occ - xnele(isp)
   deallocate ( DE_Dn, DE_Dn_old )

end function solve_occ

!----------------------------------------------------------------------
! The next function projects vector vec on the direction of vec1
! vec1 need not be normalized as that is taken care of:
!       vec_proj = vec \dot vec1 / (vec1 \dot vec1)
! ACHTUNG: only if occ_var(ia) there is contribution
function vec_proj(vec, vec1, n, occ_var)

!..Global
   use params_general
   implicit none

   real(dp) :: vec_proj
!..Arguments
   real(dp), intent(in) :: vec(lnnatorb), vec1(lnnatorb)
   logical, intent(in) :: occ_var(lnnatorb)
   integer, intent(in) :: n

!..Local variables
   integer :: ia
   real(dp) :: sum_v

   vec_proj = 0.d0
   sum_v = 0.d0
   do ia=1, n
      if(occ_var(ia)) then
         vec_proj = vec_proj + vec(ia)*vec1(ia)
         sum_v = sum_v + vec1(ia)*vec1(ia)
      endif
   enddo
   vec_proj = vec_proj / sqrt(sum_v)
end function vec_proj

!-------------------------------------------------------------------------

function grad(stp,occnumtmp)
!-------------------------------------------------------------------------
! This subroutine is just a wrapper so we have the function
! grad with only one argument to be used in root finding method
! occnum1, De_Dn1, sum_occ are  by-product outputs
! needed only when convergence is achieved
!-------------------------------------------------------------------------

!..Global
   use global; use functional_m; use func_der_args
   implicit none
      
!..The function
   real(dp) :: grad
   
!..Arguments
   real(dp), intent(in) :: stp
   real(dp), intent(in) :: occnumtmp(lnnatorb,3)

!..Local variables
   integer :: ia
   real(dp) :: tmpa
   real(dp) :: vec_proj

   do ia = 1, nnatorb
      if(.not.occ_var(ia)) then
         occnum1(ia,ispin1)=occnumtmp(ia,ispin1)
      else
         tmpa = occnumtmp(ia,ispin1) - stp*h(ia)
         occnum1(ia,ispin1) = max(min(tmpa, one),small)
      endif

      if(cl_shell) then
         occnum1(ia,ispin1_op)=occnum1(ia,ispin1)
      else
         occnum1(ia,ispin1_op)=occnumtmp(ia,ispin1_op)
      endif
      occnum1(ia,3) = occnum1(ia,ispin1) + occnum1(ia,ispin1_op)
   enddo

   call Func_der_n(xmu1, occnum1, DE_Dn1, ispin1)

   grad = vec_proj(DE_Dn1, h ,nnatorb, occ_var)
end function grad
