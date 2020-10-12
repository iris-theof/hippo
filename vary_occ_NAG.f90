subroutine vary_occ_NAG(xmu)
!-----------------------------------------------------------------------
! This subroutine varies the occ numbers using the NAG optimization
! routines. mu(output) is calculated at the end as the derivative 
! DE/Dn (n fractional).
! 
!         created by N. N. Lathiotakis and I.Theophilou
! 
!-----------------------------------------------------------------------
   use global; use functional_m; use orbocc 
   implicit none

!..Argument
   real(dp), intent(inout) :: xmu(2)

!..Parameters for the NAG routine:
   integer,parameter :: ldcj=1
   integer :: lda, liwork, lwork, ldr, iuser(1)
   integer, allocatable :: istate(:), iwork(:)
   integer :: ia, ib, nvar, nvar1, nvar2, nclin, ifail, ncnlin, iter

   real(dp), allocatable :: A(:, :), bl(:), bu(:)
   real(dp), allocatable :: X(:), work(:)
   real(dp), allocatable :: clamda(:) 
   real(dp), allocatable :: objgrd(:), r(:, :)
   real(dp), allocatable :: c(:), cjac(:,:)
   real(dp) :: user(1), objf

!..Local variables
   integer :: icnt, icnt2, mode, ihomo, imax
   real(dp) :: sumgrd, sumgrd2
   real(dp) :: sum_occ1, sum_occ2, sum_occ, s2sum, s2cor, dviol(lnnatorb,2)
   real(dp), allocatable :: occ_sorted(:)
   real(dp) :: occmax, occtmp, occ_str, occ_weak

   external Energy_FUN, e04udm, e04uef, e04ucf, confun

!..Local array allocation
   ihomo=xnele(1)

   lda=2; if( Functional == 'PN5' ) lda=ihomo
   ldr=2*lnnatorb+ihomo
   liwork=6*lnnatorb+2*lda
   lwork=55*lnnatorb+8*lnnatorb*lnnatorb+11*lda
   allocate ( istate(2*lnnatorb+lda+1), iwork(liwork), &
              A(lda,2*lnnatorb), bl(2*lnnatorb+ihomo), bu(2*lnnatorb+ihomo), &
              X(2*lnnatorb), work(lwork), &
              clamda(2*lnnatorb+ihomo+1), &
              c(ldcj), cjac(ldcj,2*lnnatorb), &
              objgrd(2*lnnatorb), r(ldr, ldr) )

   allocate( occ_sorted(lnnatorb) )

   A=0.d0

   where(occnum(:,:2) < smallocc) occnum(:,:2)=smallocc
   where(occnum(:,:2) > 1._dp-smallocc) occnum(:,:2)=1._dp-smallocc
   occnum(:,3) = occnum(:,1) +  occnum(:,2)

!  e04ucf settings
!  call e04uef("Derivative level = 3")
!  call e04uef("Verify Level = 3")
! NOTE: difference interval *must* be smaller than 'rdmepsocc'
   call e04uef("Difference interval = 1d-12")
   call e04uef("Major Iteration Limit = 10000")
   call e04uef("Minor Iteration Limit = 5000")
   call e04uef("Central difference interval = 1d-5")

!..accuracy of solution 
   call e04uef("Function Precision = 1d-16")
   call e04uef("Optimality Tolerance = 1d-12")

!..How accurately a variable should obey a constraint:
   call e04uef("Feasibility Tolerance = 1d-14")
!  call e04uef("Nonlinear feasibility = 1d-6")

!..print level
   call e04uef("Major print level = 0")
   call e04uef("Minor print level = 0")

!  Npin1=int(xnele(1))-2 !Number of states kept pinned for majority spin
!  Npin1=0
!  Npin2=0 !Number of states kept pinned to 1 for minority spin

!  Nzer1=16 !Number of states kept pinned to 0 for minority spin
!  Nzer2=16


   if ( Functional == 'PN5' ) then
      Nzer1=nnatorb-2*ihomo + Npin1
      if ( .not. cl_shell ) stop 'vary_occ_NAG: Do not know how to do open shells for PNOF5'
   endif 

   if(cl_shell) then
      Npin2=Npin1
      Nzer2=Nzer1
   endif

!..parameters of the NAG routine:
   if(cl_shell) then
      nvar1=nnatorb-Npin1-Nzer1
      nvar=nvar1
      nclin=1
   else
      nvar1= nnatorb - Npin1 - Nzer1
      nvar2= nnatorb - Npin2 - Nzer2
      nvar= nvar1 + Nvar2 !For spin up and spin down
      if(common_mu) then 
         nclin=1 !One linear constraint
      else
         nclin=2 !Two linear constraints
      endif
   endif
   if(xnele(2) <= 1.e-5_dp) then !all electrons in one spin channel
      nvar=nvar1
      nclin=1
   endif

   ncnlin=0  ! Number of non linear constraints
!  if(.not.cl_shell) ncnlin=1


   bl(1:nvar)=small !Lower bound of occ numbers
   bu(1:nvar)=1._dp-zero !Upper bound of occ numbers

   bl(nvar+1)=xnele(1)-Npin1 !The first linear constraint
   bu(nvar+1)=xnele(1)-Npin1 !The first linear constraint
   if (cl_shell .and. ncnlin>0) then
      bl(nvar+2)=0._dp      
      bu(nvar+2)=0._dp      
   endif

   if((.not.cl_shell).and.(.not.common_mu).and.(xnele(2)>=1.e-5_dp)) then
      bl(nvar+2)=xnele(2)-Npin2 !The second linear constraint
      bu(nvar+2)=xnele(2)-Npin2 !The second linear constraint
   endif
   if((.not.cl_shell).and.common_mu) then
      if(xnele(2) < 1.e-5_dp) &
         stop 'No electrons in minority chanel. common_mu MUST BE false'
      bl(nvar+1) = xnele(1)+xnele(2)-Npin1-Npin2
      bu(nvar+1) = xnele(1)+xnele(2)-Npin1-Npin2
   endif 
   mode=0

   A(1,1:nvar1) = 1._dp; A(1,nvar1+1:nvar)=0._dp
   if ((.not.cl_shell).and.common_mu) then
      A(1,nvar1+1:nvar)=1._dp
   endif
   if (.not.common_mu) then
      A(2,1:nvar1) = 0._dp
      A(2,nvar1+1:nvar) = 1._dp
   endif

   if ( Functional == 'PN5' ) then
      A=0._dp
      nclin=ihomo
      do ia=1,nclin
         bl(ia + nvar) = 1._dp
         bu(ia + nvar) = 1._dp
         A(ia,ia) = 1._dp
         A(ia,nvar-ia+1) = 1._dp
      enddo
   endif

   X(1:nvar1)=occnum(npin1+1:nnatorb-Nzer1,1)
   if((.not.cl_shell).and.(xnele(2)>1.e-5_dp)) X(nvar1+1:nvar)=occnum(npin2+1:nnatorb-Nzer2,2)

!  call TotS2violation(mode, occnum, s2sum, dviol)
   if ((.not.cl_shell) .and. ncnlin>0) then
      s2sum=max(abs(s2sum)/10,small)
      bl(nvar+3)=-1.e-3_dp
      bu(nvar+3)=1.e-3_dp
   endif

   ifail=-1

   if (Functional == 'HDR' ) call orbs_on_grid(maxorb)

 
!..CALL OF NAG OPTIMIZATION ROUTINE:
   if ( ncnlin>0 ) then
      call e04ucf(nvar, nclin, ncnlin, lda, ldcj, ldr, A, &
                  bl, bu, confun, Energy_FUN, iter, istate, c, &
                  cjac, clamda, objf, objgrd, r, X, iwork, liwork, &
                  work, lwork, iuser, user, ifail)
   else
      call e04ucf(nvar, nclin, ncnlin, lda, ldcj, ldr, A, &
                  bl, bu, e04udm, Energy_FUN, iter, istate, c, &
                  cjac, clamda, objf, objgrd, r, X, iwork, liwork, &
                  work, lwork, iuser, user, ifail)
   endif

   print*,'obj:'
   print*,(objgrd(ia),ia=1,nvar)

   occnum(Npin1+1:nnatorb-Nzer1,1)=X(1:nvar1)
   if(cl_shell) then
      occnum(Npin2+1:nnatorb-Nzer2,2)=X(1:nvar1)
   else
      if (xnele(2)>1.e-5_dp) then
         occnum(Npin2+1:nnatorb-Nzer2,2)=X(nvar1+1:nvar)
      else
         occnum(Npin2+1:nnatorb-Nzer2,2)=zero
      endif
   endif

! A small noisy negative value can crash the program. To avoid:
   occnum(Npin1+1:nnatorb-Nzer1,1) = max(occnum(Npin1+1:nnatorb-Nzer1,1),zero)
   occnum(Npin1+1:nnatorb-Nzer1,1) = min(occnum(Npin1+1:nnatorb-Nzer1,1),1._dp-zero)

   occnum(Npin1+1:nnatorb-Nzer1,3)=occnum(Npin1+1:nnatorb-Nzer1,2)+occnum(Npin1+1:nnatorb-Nzer1,1)

   icnt=0; sumgrd=0._dp; icnt2=0; sumgrd2=0._dp
   do ia=1,nnatorb
      if(1._dp-occnum(ia,1) > 1.e-7_dp.and.occnum(ia,1) > 1.e-7_dp) then
         icnt=icnt+1
         sumgrd=sumgrd+objgrd(ia-Npin1)
      endif
      if(xnele(2) >1.e-5_dp) then
         if(1._dp-occnum(ia,2) > 1e-7_dp.and.occnum(ia,2) > 1.e-7_dp) then 
            icnt2=icnt2+1
            sumgrd2=sumgrd2+objgrd(nvar1+ia-Npin2)
         endif
      endif
   enddo
   call TotS2violation(0, occnum, s2sum, dviol)

!  print*,'aaaa',xlambda_sp,s2sum
!  enddo
   xmu(1)=sumgrd/dfloat(icnt)
   if(common_mu) then
      xmu(2)=xmu(1)
   else
      xmu(2)=sumgrd2/dfloat(icnt2)
   endif

   sum_occ1=sum(occnum(:,1))
   sum_occ2=sum(occnum(:,2))

   sum_occ=sum_occ1+sum_occ2

!  call orbital_sort()

!..Enforces the aufbau principle 
!  call aufbau()

   mode=0
   if(.not.cl_shell) then
      call TotS2violation(mode, occnum, s2sum, dviol)
      s2cor=abs(xnele(1)-xnele(2))/2._dp
      s2cor=s2cor*(s2cor+1._dp)
      print*,'Error <S^2>-<S^2>_cor:',s2sum
      print*,'<S^2>_cor:', s2cor
   endif
   call calc_occ_factors(occnum, fac_h, fac_c, fac_e)

   print*,'******************************************************'
   print*,'NAG minimization routine for occ numbers:'
   print*,'Occupation numbers for spin up:'
   write(6,'(5f12.8)')(occnum(ia,1),ia=1,nnatorb)
   print*,'Mu for spin up:',xmu(1)
   print*,'Total spin up electrons: ', sum_occ1
   print*,'Occupation numbers for spin down:'
   write(6,'(5f12.8)')(occnum(ia,2),ia=1,nnatorb)
   print*,'Mu for spin dn:',xmu(2)
   print*,'Total spin down electrons: ', sum_occ2
   print*,'Total number of electrons: ', sum_occ
   print*,'******************************************************'

   occ_sorted(:)=occnum(:,1)

   do ia=1,nnatorb
      occmax=occ_sorted(ia)
      imax = ia
      do ib = ia+1, nnatorb
         if( occ_sorted(ib) > occmax ) then
            imax=ib; occmax= occ_sorted(ib)
         endif
      enddo
      occtmp=occ_sorted(ia)
      occ_sorted(ia)=occ_sorted(imax)
      occ_sorted(imax)=occtmp
   enddo

   occ_str=0._dp
   do ia=1,ihomo
       occ_str= occ_str + occ_sorted(ia)
   enddo

   occ_weak=0._dp
   do ia=ihomo+1,nnatorb
      occ_weak = occ_weak + occ_sorted(ia)
   enddo

   write(6,'("N_H_L ",5f12.8)') occ_sorted(ihomo), occ_sorted(ihomo+1), &
                                occ_sorted(ihomo)-occ_sorted(ihomo+1), &
                                occ_str, occ_weak

   deallocate (occ_sorted)

   deallocate ( istate, iwork, A, bl, bu, X, work, clamda, objgrd, r )

end subroutine vary_occ_NAG

!------------------------------------------------------------------------

subroutine Energy_FUN(mode,n,x,objf,objgrd,nstate,iuser,user)
!-----------------------------------------------------------------------
! Wrapper function for calculating the objective function
! and its gradient. It is required by NAG routine e04ucf to 
! have this particular form.
!
! it calls total_energy_MO and Func_der_n_NAG for the objective 
! function and the gradient respectively
!-----------------------------------------------------------------------


!..Global
   use global; use functional_m; use energies; use orbocc; use grid_params
   implicit none

!..Arguments
   integer :: mode, n, nstate, iuser(*)
   real(dp) :: x(n), objf, objgrd(n), user(*)

!..Local
   real(dp), allocatable :: occ(:,:), DE_Dn(:), dviol(:,:)
   real(dp) :: viol
   integer :: nspin, ispin, nvar1, nvar2, nvar

!..Local array allocation
   allocate (occ(lnnatorb,3), DE_Dn(lnnatorb), dviol(lnnatorb,2))

!!..Calculate Coul and Exch integrals over Nat orbitals J_ab and K_ab
!   if(nstate == 1) then
!      call nat_orb_basis()
!   endif
   
!  xlamspin = user(1)
   if(cl_shell.or.xnele(2)<=1.e-5_dp) then
      nspin=1
   else
      nspin=2
   endif

   nvar1=nnatorb-Npin1-Nzer1 
   nvar2=nnatorb-Npin2-Nzer2
   nvar=nvar1+nvar2

   occ(1:Npin1,1)=occnum(1:Npin1,1)
   occ(Npin1+1:nnatorb-Nzer1,1)=X(1:nvar1)
   occ(nnatorb-Nzer1+1:nnatorb,1)=occnum(nnatorb-Nzer1+1:nnatorb,1)
   if(cl_shell) then
     occ(:,2)=occ(:,1)
   else
!    if(xnele(2)>1.d-5) then
     if(xnele(2)>1.e-5_dp.and.(nvar2 >0 )) then
        occ(1:Npin2,2)=occnum(1:Npin2,2)
        occ(Npin2+1:nnatorb-Nzer2,2)=X(nvar1+1:nvar)
        occ(nnatorb-Nzer2+1:nnatorb,2)=occnum(nnatorb-Nzer2+1:nnatorb,2)
     else
        occ(:,2)=zero
     endif
   endif
   occ(:,3)=occ(:,1)+occ(:,2)
      
!..The mode cases are described in the manual of NAG e04ucf
   if(mode == 1.or.mode == 2) then !Calculate the gradient
      call TotS2violation(mode, occ, viol, dviol)
      do ispin=1,nspin
         call Func_der_n_NAG(occ, DE_Dn, ispin)
         if (ispin == 1) then
            objgrd(1:nvar1)=DE_Dn(Npin1+1:nnatorb-Nzer1)  !-xlamspin*dviol(:,ispin)
         else
            objgrd(nvar1+1:nvar)=DE_Dn(Npin2+1:nnatorb-Nzer2)  !-xlamspin*dviol(:,ispin)
         endif
      enddo
   endif
   if(mode == 0.or.mode == 2) then !Calculate the obj function
      call calc_occ_factors(occ, fac_h, fac_c, fac_e)
      if(Functional == 'HDR'.and. mode == 0) then
         call DFT_functional_libxc(occ,orbgridsq)
      endif
      call total_energy_MO(0)
      call TotS2violation(mode, occ, viol, dviol)
      objf=Tot_energy !-xlamspin*viol
   endif

   deallocate (occ, DE_Dn, dviol)
end subroutine Energy_FUN

!-----------------------------------------------------------------------
subroutine Func_der_n_NAG(occ, DE_Dn, ispin)
! This subroutine calculates the derivatives with respect
! to the occ numbers. 
! INPUT:
!     occ: occupation numbers
!   ispin: spin index
!
! OUTPUT:
!       DE_Dn : Derivative of energy w.r.t. the occ. numbers.
!       
!-----------------------------------------------------------------------

!..global
   use global; use functional_m; use matrices; use grid_params; use orbocc; use DFT_par
   implicit none

!..Arguments
   integer, intent(in) :: ispin
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: DE_Dn(lnnatorb)

!..Local variables
   integer :: ia, ib, ig
   real(dp), allocatable :: fac_c_d(:), fac_e_d(:), Vxc_tot(:), occ_pow(:,:)
   real(dp), external :: TS_nlnn_d

!..Local array allocation
   allocate ( fac_c_d(lnnatorb), fac_e_d(lnnatorb), Vxc_tot(ngrid) )
   allocate ( occ_pow(lnnatorb,3) )

   do ia=1,nnatorb
      DE_Dn(ia) = HcoreNO(ia,ia) 
      call calc_der_occ_factors(ia, ispin, occ, fac_c_d, fac_e_d)

!.....Now calculate the derivative DE/Dn_a
      do ib=1,nnatorb
         DE_Dn(ia) = DE_Dn(ia) + fac_c_d(ib)*CoulNO(ib,ia) &
                               + fac_e_d(ib)*ExchNO(ib,ia)
      enddo
   enddo

   if (( Functional == 'GPC').and.ispin ==2) then
      DE_Dn(3) = DE_dn(3)+occ(1,2)*CoulNO(1,1)
      DE_Dn(1) = DE_dn(1)+occ(3,2)*CoulNO(3,3)
      
      if ( l_non_JK ) then
         DE_Dn(1) = DE_dn(1)+ fd_non_JK1*q_non_JK 
         DE_Dn(2) = DE_dn(2)+ fd_non_JK2*q_non_JK 
         DE_Dn(3) = DE_dn(3)+ fd_non_JK3*q_non_JK 
      endif
   end if

!..Add DFT term for hybrid DFT/RDMFT
   if ( Functional == 'HDR' ) then

      call DFT_functional_libxc(occ,orbgridsq)

!$omp parallel default(shared) &
!$omp private(ig)
!$omp do
      do ig=1,ngrid
         Vxc_tot(ig) = w_grid(ig)*xmix_DF*(Vlxc(1,ig)+Vlc(1,ig))
      enddo
!$omp end do
!$omp end parallel

!$omp parallel default(shared) &
!$omp private(ig,ia)
!$omp do
     do ia=1, nnatorb
        do ig=1,ngrid
           DE_Dn(ia) = DE_Dn(ia) + Vxc_tot(ig)*orbgridsq(ig,ia)
        enddo
     enddo
!$omp end do
!$omp end parallel

   endif !( Functional == 'HDR' )

!..SLR Slater with power
   if ( Functional == 'SLR' ) then

      do ia=1,nnatorb
         occ_pow(ia,1)=occnum(ia,1)**alpha
         occ_pow(ia,2)=occnum(ia,2)**alpha
      enddo

!$omp parallel default(shared) &
!$omp private(ig,ia)
!$omp do
      do ig=1,ngrid
         chdens_gr(ig,1) = 0._dp; chdens_gr(ig,2) = 0._dp
         do ia=1,maxorb
            chdens_gr(ig,1) = chdens_gr(ig,1) + occ_pow(ia,1)*vec_nat_grid(ig,ia)**2
            chdens_gr(ig,2) = chdens_gr(ig,2) + occ_pow(ia,2)*vec_nat_grid(ig,ia)**2
         enddo
         chdens_gr(ig,3)=chdens_gr(ig,1)+chdens_gr(ig,2)
      enddo
!$omp end do
!$omp end parallel

      call DFT_functional_libxc

   print*,'--------------------------------------------------------'
   print*,'--------------------------------------------------------'
   print*,'DE_Dn',DE_Dn
   print*,'--------------------------------------------------------'
!$omp parallel default(shared) &
!$omp private(ig,ia)
!$omp do
      do ia=1,nnatorb
         do ig=1,ngrid
            DE_Dn(ia) = DE_Dn(ia) + w_grid(ig)*0.5_dp*fac_h(ia)*Vlxc(1,ig)*vec_nat_grid(ig,ia)**2
         enddo
      enddo
!$omp end do
!$omp end parallel

   endif


!..Add the entropy term:
   if(temp_dep) then
      do ia=1,nnatorb
         if(temp_dep) DE_Dn(ia) = DE_Dn(ia) + TS_nlnn_d(ia,ispin)     
      enddo
   endif

   deallocate ( fac_c_d, fac_e_d, Vxc_tot )

end subroutine Func_der_n_NAG

!-----------------------------------------------------------------------
subroutine confun(MODE,ncnln,N,LDCJ,needc,X,C,CJAC,nstate,iuser,ruser)
!..Global
   use global; use functional_m; use orbocc
   implicit none

!..Arguments
   integer :: MODE, ncnln, n, ldcj, needc(ncnln), nstate, iuser(*)
   real(dp) :: X(n), C(ncnln), CJAC(ldcj,n), ruser(*)

!..Local
   integer :: ia, nvar1, nvar2, nvar
   real(dp) :: occ(lnnatorb,3), viol, dviol(lnnatorb,2)

   nvar1=nnatorb-Npin1-Nzer1 
   nvar2=nnatorb-Npin2-Nzer2
   nvar=nvar1+nvar2

   occ(1:Npin1,1)=occnum(1:Npin1,1)
   occ(Npin1+1:nnatorb-Nzer1,1)=X(1:nvar1)
   occ(nnatorb-Nzer1+1:nnatorb,1)=occnum(nnatorb-Nzer1+1:nnatorb,1)
   if(cl_shell) then
     occ(:,2)=occ(:,1)
   else
!    if(xnele(2)>1.d-5) then
     if(xnele(2)>1.e-5_dp.and.(nvar2 >0 )) then
        occ(1:Npin2,2)=occnum(1:Npin2,2)
        occ(Npin2+1:nnatorb-Nzer2,2)=X(nvar1+1:nvar)
        occ(nnatorb-Nzer2+1:nnatorb,2)=occnum(nnatorb-Nzer2+1:nnatorb,2)
     else
        occ(:,2)=zero
     endif
   endif
   occ(:,3)=occ(:,1)+occ(:,2)

   if(mode == 0 .or. mode == 2 ) then !Calculate c
      if(needc(1) > 0) then
         call TotS2violation(mode, occ, viol, dviol)
         c(1) = viol
      endif
   endif
   if(mode == 1 .or. mode == 2) then !calculate cjac
      if(needc(1) > 0) then
         call TotS2violation(mode, occ, viol, dviol)
         do ia=1,nnatorb
            CJAC(1,ia) = dviol(ia,1)
            CJAC(1,nnatorb+ia) = dviol(ia,2)
         enddo
      endif
   endif
   CJAC=1e-5_dp*CJAC
   c(1)=1e-5_dp*c(1)
end subroutine confun

!---------------------------------------------------------------------
subroutine orbital_sort
!..Global
   use global; use orbocc
   implicit none

   integer :: ia, ib, imax
   real(dp) :: otemp, omax
   real(dp) :: vtemp(lnnatorb)

   print*,'*********BEFORE SORTING*******************************'
   print*,'NAG minimization routine for occ numbers:'
   print*,'Occupation numbers for spin up:'
   write(6,'(5f12.8)')(occnum(ia,1),ia=1,nnatorb)
   print*,'******************************************************'


   do ia=1,nnatorb-1

      imax=ia; omax=occnum(ia,1)
      do ib=ia+1,nnatorb
         if (occnum(ib,1) > omax) then
            omax=occnum(ib,1); imax=ib
         endif
      enddo

      if(imax /= ia) then
         vtemp(:)=vecnat(:,ia); otemp=occnum(ia,1)
         vecnat(:,ia) = vecnat(:,imax); occnum(ia,1) = occnum(imax,1)
         vecnat(:,imax) = vtemp(:); occnum(imax,1) = otemp
      endif

   enddo

   occnum(:,2) = occnum(:,1)
   occnum(:,3) = 2._dp*occnum(:,1)
   
end subroutine orbital_sort
!---------------------------------------------------------------------
subroutine aufbau()
!..Global
   use global; use orbocc
   implicit none

   integer :: ia, ib, imax, imin
   real(dp) :: otemp, omax, etemp, emin
   real(dp) :: vtemp(lnnatorb)
   logical :: f_aufbau


!..First order the orbitals and occ numbers according to orbital energies
   do ia=1,nnatorb-1
      imin=ia; emin=ennat(ia)
      do ib=ia+1,nnatorb
         if (ennat(ib) < emin) then
            emin=ennat(ib); imin=ib
         endif
      enddo

      if(imin /= ia) then
         vtemp(:)=vecnat(:,ia); otemp=occnum(ia,1); etemp=ennat(ia)
         vecnat(:,ia) = vecnat(:,imin); occnum(ia,1) = occnum(imin,1); ennat(ia)=emin
         vecnat(:,imin) = vtemp(:); occnum(imin,1) = otemp; ennat(imin)=etemp
      endif
   enddo

   occnum(:,2)=occnum(:,1); occnum(:,3)=occnum(:,1)+occnum(:,2)

   f_aufbau=.true.
!..Now just sort occupations 
   do ia=1,nnatorb-1
      imax=ia; omax=occnum(ia,1)
      do ib=ia+1,nnatorb
         if(occnum(ib,1) > omax) then
            omax=occnum(ib,1); imax=ib
         endif
      enddo

      if(imax /= ia) then 
         f_aufbau=.false.
         otemp=occnum(ia,1); occnum(ia,1)=omax;  occnum(imax,1) = otemp
      endif
   enddo
   occnum(:,2)=occnum(:,1); occnum(:,3)=occnum(:,1)+occnum(:,2)

   if(.not.f_aufbau) print*,'AUFBAU was corrected'

!  print*,'----------------------------------------------------------------'
!  do ia=1, nnatorb
!     write(6,'(i3,"  Energy:",f20.15,"  Occupation:",f20.15)')ia, ennat(ia), occnum(ia,3)   
!  enddo
!  print*,'----------------------------------------------------------------'

!  do ia=1, nnatorb-1
!     do ib=ia+1,nnatorb
!        print*,ia,ib,(occnum(ia,3)-occnum(ib,3))/(ennat(ib)-ennat(ia))
!     enddo
!  enddo

end subroutine aufbau
!---------------------------------------------------------------------
subroutine vary_occ_NAG_36(xmu)
!-----------------------------------------------------------------------
! This subroutine varies the occ numbers using the NAG optimization
! routines. mu(output) is calculated at the end as the derivative 
! DE/Dn (n fractional).
! 
!         created by N. N. Lathiotakis
! 
!-----------------------------------------------------------------------
   use global; use functional_m; use orbocc 
   implicit none

!..Argument
   real(dp), intent(inout) :: xmu(2)

!..Parameters for the NAG routine:
   integer,parameter :: ldcj=1
   integer :: lda, liwork, lwork, ldr, iuser(1)
   integer, allocatable :: istate(:), iwork(:)
   integer :: ia, ib, nvar, nvar1, nvar2, nclin, ifail, ncnlin, iter

   real(dp), allocatable :: A(:, :), bl(:), bu(:)
   real(dp), allocatable :: X(:), work(:)
   real(dp), allocatable :: clamda(:) 
   real(dp), allocatable :: objgrd(:), r(:, :)
   real(dp), allocatable :: c(:), cjac(:,:)
   real(dp) :: user(1), objf

!..Local variables
   integer :: icnt, icnt2, mode, ihomo, imax
   real(dp) :: sumgrd, sumgrd2
   real(dp) :: sum_occ1, sum_occ2, sum_occ, s2sum, s2cor, dviol(lnnatorb,2)
   real(dp), allocatable :: occ_sorted(:)
   real(dp) :: occmax, occtmp, occ_str, occ_weak

   external Energy_FUN, e04udm, e04uef, e04ucf, confun

!..Local array allocation
   ihomo=xnele(1)
   nclin=5
   lda=nclin; if( Functional == 'PN5' ) lda=ihomo
   ldr=2*lnnatorb+nclin
   liwork=6*lnnatorb+2*lda
   lwork=55*lnnatorb+8*lnnatorb*lnnatorb+11*lda
   allocate ( istate(2*lnnatorb+lda+1), iwork(liwork), &
              A(lda,2*lnnatorb), bl(2*lnnatorb+ihomo), bu(2*lnnatorb+ihomo), &
              X(2*lnnatorb), work(lwork), &
              clamda(2*lnnatorb+ihomo+1), &
              c(ldcj), cjac(ldcj,2*lnnatorb), &
              objgrd(2*lnnatorb), r(ldr, ldr) )

   allocate( occ_sorted(lnnatorb) )

   A=0.d0



   where(occnum(:,:2) < smallocc) occnum(:,:2)=smallocc
   where(occnum(:,:2) > 1._dp-smallocc) occnum(:,:2)=1._dp-smallocc
   occnum(:,3) = occnum(:,1) +  occnum(:,2)

!  e04ucf settings
!  call e04uef("Derivative level = 3")
!  call e04uef("Verify Level = 3")
! NOTE: difference interval *must* be smaller than 'rdmepsocc'
   call e04uef("Difference interval = 1d-12")
   call e04uef("Major Iteration Limit = 10000")
   call e04uef("Minor Iteration Limit = 5000")
   call e04uef("Central difference interval = 1d-5")

!..accuracy of solution 
   call e04uef("Function Precision = 1d-12")
   call e04uef("Optimality Tolerance = 1d-10")

!..How accurately a variable should obey a constraint:
   call e04uef("Feasibility Tolerance = 1d-12")
!  call e04uef("Nonlinear feasibility = 1d-6")

!..print level
   call e04uef("Major print level = 0")
   call e04uef("Minor print level = 0")
!  Npin1=int(xnele(1))-2 !Number of states kept pinned for majority spin
!  Npin1=0
!  Npin2=0 !Number of states kept pinned to 1 for minority spin

!  Nzer1=16 !Number of states kept pinned to 0 for minority spin
!  Nzer2=16
!..parameters of the NAG routine:
     nvar1=3
     nvar2=3
     nvar=nvar1+nvar2

   if ( Functional == 'PN5' ) then
      Nzer1=nnatorb-2*ihomo + Npin1
      if ( .not. cl_shell ) stop 'vary_occ_NAG: Do not know how to do open shells for PNOF5'
   endif 

   if(cl_shell) then
      Npin2=Npin1
      Nzer2=Nzer1
   endif


   ncnlin=0  ! Number of non linear constraints
!  if(.not.cl_shell) ncnlin=1


   bl(1:nvar)=small !Lower bound of occ numbers
   bu(1:nvar)=1._dp-zero !Upper bound of occ numbers

   bl(nvar+1)=xnele(1)-Npin1 !The first linear constraint
   bu(nvar+1)=xnele(1)-Npin1 !The first linear constraint
   if (cl_shell .and. ncnlin>0) then
      bl(nvar+2)=0._dp      
      bu(nvar+2)=0._dp      
   endif

   if((.not.cl_shell).and.(.not.common_mu).and.(xnele(2)>=1.e-5_dp)) then
      bl(nvar+2)=xnele(2)-Npin2 !The second linear constraint
      bu(nvar+2)=xnele(2)-Npin2 !The second linear constraint
   endif
   if((.not.cl_shell).and.common_mu) then
      if(xnele(2) < 1.e-5_dp) &
         stop 'No electrons in minority chanel. common_mu MUST BE false'
      bl(nvar+1) = xnele(1)+xnele(2)-Npin1-Npin2
      bu(nvar+1) = xnele(1)+xnele(2)-Npin1-Npin2
   endif 
   mode=0

   A(1,1:nvar1) = 1._dp; A(1,nvar1+1:nvar)=0._dp
   if ((.not.cl_shell).and.common_mu) then
      A(1,nvar1+1:nvar)=1._dp
   endif
   if (.not.common_mu) then
      A(2,1:nvar1) = 0._dp
      A(2,nvar1+1:nvar) = 1._dp
   endif

   if ( Functional == 'PN5' ) then
      A=0._dp
      nclin=ihomo
      do ia=1,nclin
         bl(ia + nvar) = 1._dp
         bu(ia + nvar) = 1._dp
         A(ia,ia) = 1._dp
         A(ia,nvar-ia+1) = 1._dp
      enddo
   endif
   !lambda1up+lambda3down=1
     A(nclin-2,1)=1._dp
     A(nclin-2,nvar1+3)=1._dp

   !lambda2up+lambda2down=1
     A(nclin-1,2)=1._dp
     A(nclin-1,nvar1+2)=1._dp

   !lambda3up+lambda1down=1
     A(nclin,3)=1._dp
     A(nclin,nvar1+1)=1._dp

   bl(nvar+nclin-2)=1._dp
   bu(nvar+nclin-2)=1._dp
   bl(nvar+nclin-1)=1._dp
   bu(nvar+nclin-1)=1._dp
   bl(nvar+nclin)=1._dp
   bu(nvar+nclin)=1._dp

   X(1:nvar1)=occnum(npin1+1:nnatorb-Nzer1,1)
   if((.not.cl_shell).and.(xnele(2)>1.e-5_dp)) X(nvar1+1:nvar)=occnum(npin2+1:nnatorb-Nzer2,2)

!  call TotS2violation(mode, occnum, s2sum, dviol)
   if ((.not.cl_shell) .and. ncnlin>0) then
      s2sum=max(abs(s2sum)/10,small)
      bl(nvar+3)=-1.e-3_dp
      bu(nvar+3)=1.e-3_dp
   endif

   ifail=-1

   if (Functional == 'HDR' ) call orbs_on_grid


!..CALL OF NAG OPTIMIZATION ROUTINE:
   if ( ncnlin>0 ) then
      call e04ucf(nvar, nclin, ncnlin, lda, ldcj, ldr, A, &
                  bl, bu, confun, Energy_FUN, iter, istate, c, &
                  cjac, clamda, objf, objgrd, r, X, iwork, liwork, &
                  work, lwork, iuser, user, ifail)
   else
      call e04ucf(nvar, nclin, ncnlin, lda, ldcj, ldr, A, &
                  bl, bu, e04udm, Energy_FUN, iter, istate, c, &
                  cjac, clamda, objf, objgrd, r, X, iwork, liwork, &
                  work, lwork, iuser, user, ifail)
   endif

   print*,'obj:'
   print*,(objgrd(ia),ia=1,nvar)

   occnum(Npin1+1:nnatorb-Nzer1,1)=X(1:nvar1)
   if(cl_shell) then
      occnum(Npin2+1:nnatorb-Nzer2,2)=X(1:nvar1)
   else
      if (xnele(2)>1.e-5_dp) then
         occnum(Npin2+1:nnatorb-Nzer2,2)=X(nvar1+1:nvar)
      else
         occnum(Npin2+1:nnatorb-Nzer2,2)=zero
      endif
   endif

! A small noisy negative value can crash the program. To avoid:
   occnum(Npin1+1:nnatorb-Nzer1,1) = max(occnum(Npin1+1:nnatorb-Nzer1,1),zero)
   occnum(Npin1+1:nnatorb-Nzer1,1) = min(occnum(Npin1+1:nnatorb-Nzer1,1),1._dp-zero)

   occnum(Npin1+1:nnatorb-Nzer1,3)=occnum(Npin1+1:nnatorb-Nzer1,2)+occnum(Npin1+1:nnatorb-Nzer1,1)

   icnt=0; sumgrd=0._dp; icnt2=0; sumgrd2=0._dp
   do ia=1,nnatorb
      if(1._dp-occnum(ia,1) > 1.e-7_dp.and.occnum(ia,1) > 1.e-7_dp) then
         icnt=icnt+1
         sumgrd=sumgrd+objgrd(ia-Npin1)
      endif
      if(xnele(2) >1.e-5_dp) then
         if(1._dp-occnum(ia,2) > 1e-7_dp.and.occnum(ia,2) > 1.e-7_dp) then 
            icnt2=icnt2+1
            sumgrd2=sumgrd2+objgrd(nvar1+ia-Npin2)
         endif
      endif
   enddo
   call TotS2violation(0, occnum, s2sum, dviol)

!  print*,'aaaa',xlambda_sp,s2sum
!  enddo
   xmu(1)=sumgrd/dfloat(icnt)
   if(common_mu) then
      xmu(2)=xmu(1)
   else
      xmu(2)=sumgrd2/dfloat(icnt2)
   endif

   sum_occ1=sum(occnum(:,1))
   sum_occ2=sum(occnum(:,2))

   sum_occ=sum_occ1+sum_occ2

!  call orbital_sort()

!..Enforces the aufbau principle 
!  call aufbau()

   mode=0
   if(.not.cl_shell) then
      call TotS2violation(mode, occnum, s2sum, dviol)
      s2cor=abs(xnele(1)-xnele(2))/2._dp
      s2cor=s2cor*(s2cor+1._dp)
      print*,'Error <S^2>-<S^2>_cor:',s2sum
      print*,'<S^2>_cor:', s2cor
   endif
   call calc_occ_factors(occnum, fac_h, fac_c, fac_e)

   print*,'******************************************************'
   print*,'NAG minimization routine for occ numbers:'
   print*,'Occupation numbers for spin up:'
   write(6,'(5f12.8)')(occnum(ia,1),ia=1,nnatorb)
   print*,'Mu for spin up:',xmu(1)
   print*,'Total spin up electrons: ', sum_occ1
   print*,'Occupation numbers for spin down:'
   write(6,'(5f12.8)')(occnum(ia,2),ia=1,nnatorb)
   print*,'Mu for spin dn:',xmu(2)
   print*,'Total spin down electrons: ', sum_occ2
   print*,'Total number of electrons: ', sum_occ
   print*,'******************************************************'

   occ_sorted(:)=occnum(:,1)

   do ia=1,nnatorb
      occmax=occ_sorted(ia)
      imax = ia
      do ib = ia+1, nnatorb
         if( occ_sorted(ib) > occmax ) then
            imax=ib; occmax= occ_sorted(ib)
         endif
      enddo
      occtmp=occ_sorted(ia)
      occ_sorted(ia)=occ_sorted(imax)
      occ_sorted(imax)=occtmp
   enddo

   occ_str=0._dp
   do ia=1,ihomo
       occ_str= occ_str + occ_sorted(ia)
   enddo

   occ_weak=0._dp
   do ia=ihomo+1,nnatorb
      occ_weak = occ_weak + occ_sorted(ia)
   enddo

   write(6,'("N_H_L ",5f12.8)') occ_sorted(ihomo), occ_sorted(ihomo+1), &
                                occ_sorted(ihomo)-occ_sorted(ihomo+1), &
                                occ_str, occ_weak

   deallocate (occ_sorted)

   deallocate ( istate, iwork, A, bl, bu, X, work, clamda, objgrd, r )

end subroutine vary_occ_NAG_36


