subroutine Optimize_FD(iscf)
!.. Optimizes the total energy with respect to beta and mu
!.. under the constraint that the sum of occupations per spin is N/2
!..Global
   use global; use orbocc, ONLY: occnum, ennat, beta_occ, xmu_occ
   implicit none

!..Argument
   integer :: iscf

!..Local
   real(dp) :: sum_occ, sum_occ1, sum_occ2, Sgra, exp_mine
   integer :: ia

!..FOR NAG MINIMIZATION
   integer :: Nvar=2, nclin=0, ncnlin=1, lda=1, ldcj=1, ldr=2, iter, istate(3)
   integer :: liwork=8, lwork=100, iuser(1), ifail
   real(dp) :: A(1,1), BL(3), BU(3), C(1), cjac(1,2), clamda(3), objf, objgrd(2)
   real(dp) :: R(2,2), X(2), iwork(8), work(100), user(1) 

   external :: E_OBJ_FERMI, confun_FD, e04uef, e04ucf !, e04djf, e04dgf

   print*,'---------------------------------------------------------'
   print*,'Minimizing occupations assuming Fermi Dirac distribution:'
   
!..Optional parameters for the NAG routine

!..Initialize mu and beta for first scf step
   if(iscf == 1) then
!     xmu_occ(1) = 0.5d0*(ennat(ibond(1)) + ennat(1+ibond(1)))
      xmu_occ(1) = ennat(ibond(1))
      Temperat= 1.e6_dp
      beta_occ(1) = 1._dp/(Boltz_konst*Temperat)
      Print*, 'Initial beta, mu:', beta_occ(1), xmu_occ(1)
   endif
   X(1) = beta_occ(1); X(2) = xmu_occ(1)

   call E_OBJ_FERMI(2, 2, X, OBJF, OBJGRD, 0, IUSER, USER)
   Sgra=0._dp
   do ia=1,Nvar
      Sgra = Sgra+ OBJGRD(ia)*OBJGRD(ia)
   enddo

   print*,'Sgra: ',Sgra
!  if (Sgra > 1.d-8) then
      print*,'Optimizing beta and mu...'
  
      BL(1) = 1.d-10
      BU(1) = 1d20
      
!     BL(2) = ennat(ibond(1)-1)
!     BU(2) = ennat(ibond(1)+2)
      BL(2) = -1d20
      BU(2) = 1d20

      BL(3) = xnele(1)
      BU(3) = xnele(1)

!  e04ucf settings
!  call e04uef("Derivative level = 3")
!  call e04uef("Verify Level = 3")
! NOTE: difference interval *must* be smaller than 'rdmepsocc'
      call e04uef("Difference interval = 1d-12")
      call e04uef("Major Iteration Limit = 1000")
      call e04uef("Minor Iteration Limit = 5000")
      call e04uef("Central difference interval = 1d-5")

!..accuracy of solution 
      call e04uef("Function Precision = 1d-9")
      call e04uef("Optimality Tolerance = 1d-8")

!..How accurately a variable should obey a constraint:
      call e04uef("Feasibility Tolerance = 1d-12")

!..print level
      call e04uef("Major print level = 0")
      call e04uef("Minor print level = 0")


      call e04ucf(Nvar, nclin, ncnlin, lda, ldcj, ldr, A, &
               bl, bu, confun_FD, E_OBJ_FERMI, iter, istate, c, &
               cjac, clamda, objf, objgrd, r, X, iwork, liwork, &
               work, lwork, iuser, user, ifail)

      beta_occ(1)=X(1)
      xmu_occ(1)=X(2)
      beta_occ(2)=X(1)
      xmu_occ(2)=X(2)

      do ia=1,nnatorb
         occnum(ia,1) = 1._dp/(1._dp+exp_mine(X(1)*(ennat(ia)-X(2))))
         occnum(ia,2) = occnum(ia,1)
         occnum(ia,3) = occnum(ia,1)+occnum(ia,2)
      enddo
!  else
!     print*,'Beta and mu already optimal no need to optimize further'
!  endif

   open(unit=121,file='n_vs_e',status='unknown')
   do ia=1,nnatorb
      write(121,'(2f25.10)')ennat(ia),occnum(ia,1)
   enddo
   close(121)

   Temperat=1._dp/Boltz_konst/beta_occ(1)
   print*,'Optimals:'
   print*,'Effective    T:', Temperat, 'Kelvin'
   print*,'Effective beta:', beta_occ(1), '(a.u.)'
   print*,'Effective   mu:', xmu_occ(1), '(a.u.)'
   print*,'-----------------------------------------------------'

   sum_occ1=sum(occnum(:,1))
   sum_occ2=sum(occnum(:,2))
   sum_occ=sum_occ1+sum_occ2

   print*,'Occupation numbers for spin up:'
   write(6,'(5f12.8)')(occnum(ia,1),ia=1,nnatorb)
   print*,'Mu for spin up (From Fermi Dirac):',xmu_occ(1)
   print*,'Total spin up electrons: ', sum_occ1
   print*,'Occupation numbers for spin down:'
   write(6,'(5f12.8)')(occnum(ia,2),ia=1,nnatorb)
   print*,'Mu for spin dn (From Fermi Dirac):',xmu_occ(2)
   print*,'Total spin down electrons: ', sum_occ2
   print*,'Total number of electrons: ', sum_occ
   print*,'******************************************************'

end subroutine Optimize_FD

!-----------------------------------------------------------------------
subroutine E_OBJ_FERMI(MODE, N, X, OBJF, OBJGRD, NSTATE, IUSER, RUSER)
!..Global
   use global; use functional_m; use energies; use orbocc 
   implicit none

!..Arguments
   integer :: MODE, N, NSTATE, IUSER(*)
   real(dp) :: X(N), OBJF, OBJGRD(N), RUSER(*) 

!..Local
   real(dp), allocatable :: occ(:,:), DE_Dn(:)
   real(dp) :: dm, drvia, fcdrv, exp_mine, rdum
   integer :: ia, i, idum

   allocate (occ(lnnatorb,3), DE_Dn(lnnatorb))

   idum=nstate; idum=iuser(1); rdum=ruser(1)
   
   do ia=1,nnatorb
      dm = X(1)*(ennat(ia)-X(2))
      occ(ia,1) = 1._dp/(1._dp+exp_mine(dm))
      occ(ia,1) = max(occ(ia,1),zero)
      occ(ia,2) = occ(ia,1)
      occ(ia,3) = occ(ia,1)+occ(ia,2)
   enddo
   call calc_occ_factors(occ, fac_h, fac_c, fac_e)

   if (MODE == 0 .or. mode==2) then
      call total_energy_MO(1)
      OBJF = Tot_energy
   endif
   if (mode ==1 .or. mode==2) then
      call Func_der_n_NAG(occ, DE_Dn, 1)

      OBJGRD = 0._dp
      do ia = 1,nnatorb
         dm = exp_mine(X(1)*(ennat(ia)-X(2)))
         drvia = dm/(1._dp+dm)**2
         do i = 1,2 
            if(i == 1) then
               fcdrv = (X(2) - ennat(ia))*drvia
            else
               fcdrv = X(1)*drvia
            endif
            OBJGRD(i) = OBJGRD(i) + 2._dp*fcdrv*DE_Dn(ia) !2 is for spin
         enddo
      enddo
   endif

end subroutine E_OBJ_FERMI

!-----------------------------------------------------------------------
subroutine confun_FD(MODE,ncnln,N,LDCJ,needc,X,C,CJAC,nstate,iuser,ruser)
!..Global
   use global; use functional_m; use energies; use orbocc
   implicit none

!..Arguments
   integer :: MODE, ncnln, n, ldcj, needc(ncnln), nstate, iuser(*)
   real(dp) :: X(n), C(ncnln), CJAC(ldcj,n), ruser(*)

!..Local
   real(dp) :: occ(lnnatorb)
   real(dp) :: dm, drvia, exp_mine
   integer :: ia

   do ia=1,nnatorb
      dm = exp_mine(X(1)*(ennat(ia)-X(2)))
      occ(ia) = 1._dp/(1._dp+dm)
   enddo

   if(mode == 0 .or. mode == 2 ) then !Calculate c
      if(needc(1) > 0) then
         c(1) = 0._dp
         do ia=1,nnatorb
            c(1)=c(1) + occ(ia)
         enddo
      endif
   endif
   if(mode == 1 .or. mode == 2) then !calculate cjac
      if(needc(1) > 0) then
         cjac(1,1) = 0._dp; cjac(1,2) = 0._dp
         do ia=1,nnatorb
            dm = exp_mine(X(1)*(ennat(ia)-X(2)))
            drvia = dm/(1._dp+dm)**2
            cjac(1,1) = cjac(1,1) + (X(2) - ennat(ia))*drvia
            cjac(1,2) = cjac(1,2) + X(1)*drvia
         enddo
      endif
   endif

end subroutine confun_FD
!-----------------------------------------------------------------------
subroutine Fit_Fermi(iscf,xmu)
!..Given a set of occupation numbers occnum and a set of energies ennat
!..It fits a 2 Fermi functions (1 per spin) yielding the effecive beta and
!..mu per spin.
!..Global
   use global; use orbocc, ONLY: occnum, ennat, beta_occ, xmu_occ
   implicit none

!..Arguments
   integer :: iscf
   real(dp) :: xmu(2)

!..Local
   real(dp) :: xin(2), Sgra, exp_mine
   integer :: ia

!..FOR NAG MINIMIZATION
   integer :: IWORK(3), IUSER(1), IFAIL, INFORM
   real(dp) :: OBJF, OBJGRD(2), WORK(26), RUSER(1)

   external :: OBJ_FERMI, e04djf, e04dgf

   print*,'-----------------------------------------------------'
   print*,'Effective Fermi Dirac Fitting:'
   
!..Optional parameters for the NAG routine
   open(unit=8, file='optinals_nag', status='unknown')
      write(8,'(a5)')'Begin'
      write(8,*)'  Print Level = 0'
      write(8,*)'  Optimality Tolerance = 1D-9'
      write(8,*)'  Verify = YES'
      write(8,*)'  Verify Level = 1'
      write(8,'(a3)')'End'
      rewind(8)
      call e04djf(8,INFORM)
   close(8)

!..Initialize mu and beta for first scf step
   if(iscf == 1) then
      xmu_occ=xmu
      beta_occ(1) = log((1.d0/occnum(1+ibond(1),1))-1.d0) / (ennat(1+ibond(1))-xmu_occ(1))
   endif
   xin(1) = beta_occ(1); xin(2) = xmu_occ(1)

   call OBJ_FERMI(2, 2, Xin, OBJF, OBJGRD, 0, IUSER, RUSER)
   Sgra=0._dp
   do ia=1,nnatorb
      Sgra = Sgra+ OBJGRD(ia)*OBJGRD(ia)
   enddo

   print*,'Sgra: ',Sgra
   if (Sgra > 1.e-8_dp) then
      print*,'Optimizing beta and mu...'
!     call e04dgf(2, OBJ_FERMI, ITER, OBJF, OBJGRD, Xin, IWORK, WORK, IUSER, RUSER, IFAIL)
      beta_occ(1)=xin(1)
      xmu_occ(1)=xin(2)
   else
       Print*,'Beta and mu already optimal no need to optimize further!'
   endif

   open(unit=121,file='n_vs_e',status='unknown')
   do ia=1,nnatorb
      write(121,'(2f25.10)')ennat(ia),occnum(ia,1)
   enddo
   close(121)

   Temperat=1._dp/(Boltz_konst*  beta_occ(1))

   Temperat=1.e-3_dp
   beta_occ(1) = 1.0_dp/(Boltz_konst*Temperat)
   xmu_occ=xmu

   do ia=1,nnatorb
      occnum(ia,1) = 1.0_dp/(1.0_dp + exp_mine(beta_occ(1)*(ennat(ia)-xmu_occ(1))))
      occnum(ia,2) = occnum(ia,1); occnum(ia,3)=2.0_dp*occnum(ia,1)
   enddo

   print*,'Effective    T:', Temperat, 'Kelvin'
   print*,'Effective beta:', beta_occ(1), '(a.u.)'
   print*,'Effective   mu:', xmu_occ(1), '(a.u.)'
   print*,'-----------------------------------------------------'

end subroutine Fit_Fermi

!-----------------------------------------------------------------------
subroutine OBJ_FERMI(MODE, N, X, OBJF, OBJGRD, NSTATE, IUSER, RUSER)
!..Global
   use global
   implicit none

!..Arguments
   integer :: MODE, N, NSTATE, IUSER(*)
   real(dp) :: X(N), OBJF, OBJGRD(N), RUSER(*) 
   real(dp), external :: Fermi_LSQ

   if (MODE == 0) then
      OBJF = Fermi_LSQ(X)   
   elseif (MODE == 2) then
      OBJF = Fermi_LSQ(X) 
      call DFermi_LSQ(X, OBJGRD)
   else
      stop 'OBJ_FERMI:wrong MODE' 
   endif
end subroutine OBJ_FERMI

!-----------------------------------------------------------------------
function Fermi_LSQ(xin)
!..Global
   use global; use orbocc, ONLY: occnum, ennat
   implicit none

!..Arguments
   real(dp) :: xin(2)

!..The function
   real(dp) :: Fermi_LSQ

!..Local
   real(dp) :: sum_d, dm
   integer :: ia

   sum_d=0.0_dp
   do ia=1,nnatorb
      dm = xin(1)*(ennat(ia) - xin(2))
      sum_d = sum_d + (occnum(ia,1) - 1.0_dp/(1.0_dp + exp(dm)))**2
   enddo
   Fermi_LSQ = sum_d
   return 
end function Fermi_LSQ
!-----------------------------------------------------------------------

subroutine DFermi_LSQ(xin, Der_x)
!..Global
   use global; use orbocc, ONLY: occnum, ennat
   implicit none

!..Arguments
   real(dp) :: xin(2), Der_x(2)
!..Local
   integer :: ia
   real(dp) :: dm, ff
   
   Der_x=0.0_dp
   do ia=1,nnatorb
      dm = xin(1)*(ennat(ia) - xin(2))
      ff = (1.0_dp/(1.0_dp+exp(dm))-occnum(ia,1))/(cosh(dm)+1.0_dp)
      Der_x(1) = Der_x(1) + ff*(xin(2)-ennat(ia))
      Der_x(2) = Der_x(2) + ff*xin(1)
   enddo

   return
end subroutine DFermi_LSQ

!----------------------------------------------------------------------------
function Dn_De_FD(epsi)
! The derivative of Fermi dirac with respect to epsilon
!..Global
   use global; use orbocc, ONLY: beta_occ, xmu_occ
   implicit none
   
   real(dp) :: Dn_De_FD, epsi
    
   Dn_De_FD = -0.5_dp*beta_occ(1)/(cosh(beta_occ(1)*(epsi-xmu_occ(1))) + 1.0_dp)

end function Dn_De_FD
!----------------------------------------------------------------------------

function exp_mine(x)
   use global
   real(dp), intent(IN) :: x
   real(dp) :: exp_mine, y
   
   y=min(max(x,-200.0_dp),200.0_dp)
   exp_mine=exp(y)
end function exp_mine
