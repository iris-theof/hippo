!--------------------------------------------------------------------------
! This subroutine constructs the Fock-matrices F^(ia), In case of HF F is 
! independent of ia. In case of DFT it constructs the Konh Sham Hamiltonian.
!
! INPUT (through common variables):
!        occnum: the occupation numbers
!        vecnat: the natural orbitals
!        nbasis: nu of basis functions
!        nnatorb: nu of natural orbitals
!        nintgr: number of non zero two integrals
!        hcore: The 1-e integrals (kinetic+external potential)
!        twoin: the 2-e integrals
!        mu,nu,lambda,sigma: integer matrices for indirect indexing of
!                            2-e integrals
! OUTPUT:
!        F: The fock matrix
!        Coul: The coulomb matrix
!        Exch: The exchange matrix
!  created by N. N. Lathiotakis and I.Theophilou
!--------------------------------------------------------------------------
subroutine construct_f(iaa)

!..Global      
   use global; use functional_m; use matrices; use orbocc; use DFT_par; use grid_params
   implicit none

!..Arguments
   integer, intent(in) :: iaa

!..Local variables
   integer :: i, ia, ib, m, n, ig
   integer :: k,l,ic
   real(dp) :: DXX(lnbasis, lnbasis,3), DM(lnbasis, lnbasis), DMX(3), ss, dd
   real(dp) :: occ_pow(lnnatorb,3), X_int, vv(lnbasis)
   real(dp) :: fac_cM(lnnatorb,lnnatorb), fac_eM(lnnatorb,lnnatorb)

   if(iaa<0.or.iaa>1) stop 'In construct_f: Invalid iaa'
!..iaa: =0 construct the HcoreNO, CoulNO, ExchNO
!       =1 Do not construct HcoreNO, CoulNO, ExchNO
   
   !..Sum the 2-integrals
   call sum_intg()

   if ( functional /= 'DBF' ) then
      
      if ( l_hybc ) then
         if ( .not. allocated( VHF ) ) allocate( VHF(nbasis, nbasis))
         do m=1,nbasis
            do n=1,m
               VHF(m,n)=0._dp
               do ib=1,maxorb
                  VHF(m,n)= VHF(m,n) + 2._dp*Coul(m,n,ib) - Exch(m,n,ib)
               enddo
               VHF(n,m)= VHF(m,n)
            enddo
         enddo
      endif

   do ia=1,nnatorb
      do ib=1,ia
         fac_cM(ia,ib) = (occnum(ia,1)+occnum(ia,2))*(occnum(ib,1)+occnum(ib,2))
         fac_eM(ia,ib) = -occnum(ia,1)*occnum(ib,1)-occnum(ia,2)*occnum(ib,2)
!........The double counting 1/2 factor
         fac_cM(ia,ib) = 0.5_dp*fac_cM(ia,ib)
         fac_eM(ia,ib) = 0.5_dp*fac_eM(ia,ib)
         fac_cM(ib,ia) = fac_cM(ia,ib)
         fac_eM(ib,ia) = fac_eM(ia,ib)
      enddo
   enddo

!..The generalized Fock Matrix
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ib,m,n)
!$OMP DO
      do ia=1, maxorb
         do m=1,nbasis
            do n=1,m
               F_xc(m,n,ia)=0._dp; F_ha(m,n,ia)=0._dp; F_haM(m,n,ia)=0._dp;  F_x(m,n,ia)=0._dp
               F(m,n,ia)=Hcore(m,n)*0.5_dp*fac_h(ia) 
               do ib=1,maxorb
                  F(m,n,ia)= F(m,n,ia) + fac_c(ia,ib)*Coul(m,n,ib) + fac_e(ia,ib)*Exch(m,n,ib)
                  F_xc(m,n,ia)= F_xc(m,n,ia) + fac_e(ia,ib)*Exch(m,n,ib)
                  F_x(m,n,ia)= F_x(m,n,ia) + fac_eM(ia,ib)*Exch(m,n,ib)
                  F_ha(m,n,ia)= F_ha(m,n,ia) + fac_c(ia,ib)*Coul(m,n,ib)
                  F_haM(m,n,ia)= F_haM(m,n,ia) + fac_cM(ia,ib)*Coul(m,n,ib)
               enddo
               F(n,m,ia)=F(m,n,ia)
               F_xc(n,m,ia)=F_xc(m,n,ia)
               F_x(n,m,ia)=F_x(m,n,ia)
               F_ha(n,m,ia)=F_ha(m,n,ia)
               F_haM(n,m,ia)=F_haM(m,n,ia)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
   else
!..The generalized Fock Matrix
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ib,m,n)
!$OMP DO
      do ia=1, maxorb
         do m=1,nbasis
            do n=1,m
               F(m,n,ia)=Hcore(m,n)*0.5_dp*fac_h(ia) 
               do ib=1,nnatorb
                  F(m,n,ia)= F(m,n,ia) + fac_c(ia,ib)*Coul(m,n,ib) ! Exc from DFT
               enddo
               F(n,m,ia)=F(m,n,ia)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
   endif !( functional /= 'DBF' )

   if(iaa==0) then
!.....Calculate CoulNO, ExchNO, HcoreNO
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ib,m,n,dd)
!$OMP DO
      do ia=1,maxorb
         HcoreNO(ia,ia) = 0.0_dp
         do m =1, nbasis
            do n = 1, nbasis
               HcoreNO(ia,ia) = HcoreNO(ia,ia) + vecnat(m,ia)*hcore(m,n)*vecnat(n,ia)
            enddo
         enddo
         do ib = 1, ia
            CoulNO(ia,ib)=0.0_dp; ExchNO(ia,ib)=0.0_dp
            do m =1, nbasis
               do n = 1, nbasis
                  dd=vecnat(m,ia)*vecnat(n,ia)
                  CoulNO(ia,ib) = CoulNO(ia,ib) + Coul(m,n,ib) * dd
                  ExchNO(ia,ib) = ExchNO(ia,ib) + Exch(m,n,ib) * dd 
               enddo ! m
            enddo ! n
            CoulNO(ib,ia) = CoulNO(ia,ib); ExchNO(ib,ia) = ExchNO(ia,ib)
         enddo ! ib
      enddo !ia
!$OMP END DO
!$OMP END PARALLEL
   endif !(iaa==0)

   if (l_non_JK) then
      call psi4_int(vecnat(:,1),vecnat(:,2),vecnat(:,2),vecnat(:,3),X_int)
      q_non_JK=X_int ! \int\int dr dr\prime \frac{phi_1(r) phi_2(r) phi_2(r\prime) phi_3(r\prime)}{|r-r\prime|}
      if ( .not. allocated(XI1) ) allocate ( XI1(lnnatorb), XI2(lnnatorb), XI3(lnnatorb) )
      do i=1,nbasis
         vv=0._dp; vv(i)=1._dp
         call psi4_int(vv,vecnat(:,2),vecnat(:,2),vecnat(:,3),X_int)
         XI1(i) = f_non_JK*X_int
         call psi4_int(vv,vecnat(:,1),vecnat(:,2),vecnat(:,3),X_int)
         XI2(i) = f_non_JK*X_int
         call psi4_int(vv,vecnat(:,3),vecnat(:,1),vecnat(:,2),X_int)
         XI2(i) = f_non_JK*X_int
         call psi4_int(vv,vecnat(:,2),vecnat(:,2),vecnat(:,1),X_int)
         XI3(i) = f_non_JK*X_int
      enddo
   endif !(l_non_JK)
   
   if(Functional == 'SLR') then
      call orbs_on_grid(nnatorb)

      do ia=1,nnatorb
         occ_pow(ia,1) = occnum(ia,1)**alpha
         occ_pow(ia,2) = occnum(ia,2)**alpha
      enddo

      if ( .not.allocated(chdens_gr) ) allocate(chdens_gr(ngrid,3))

!$omp parallel default(shared) &
!$omp private(ig,ia)
!$omp do
       do ig=1,ngrid
         chdens_gr(ig,1) = 0.0_dp; chdens_gr(ig,2) = 0.0_dp
         do ia=1,maxorb
            chdens_gr(ig,1) = chdens_gr(ig,1) + occ_pow(ia,1)*vec_nat_grid(ig,ia)**2
            chdens_gr(ig,2) = chdens_gr(ig,2) + occ_pow(ia,2)*vec_nat_grid(ig,ia)**2
         enddo
         chdens_gr(ig,3)=chdens_gr(ig,1)+chdens_gr(ig,2)
      enddo
!$omp end do
!$omp end parallel

      call DFT_functional_libxc(occnum,orbgridsq)
!     call DFT_functional_mine

      if(.not.allocated(Vxc_mn_DFT)) allocate (Vxc_mn_DFT(nbasis,nbasis))

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,n,ig)
!$OMP DO
      do m=1,nbasis
         do n=1,m
            Vxc_mn_DFT(m,n)=0.0_dp
            do ig=1,ngrid
               Vxc_mn_DFT(m,n) = Vxc_mn_DFT(m,n) + w_grid(ig)*bas_f_grid(ig,m)*Vlxc(1,ig)*bas_f_grid(ig,n)
            enddo
            Vxc_mn_DFT(n,m)=Vxc_mn_DFT(m,n)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      do ia=1,maxorb
         do m=1,nbasis
            do n=1,m
               F(m,n,ia) = F(m,n,ia)+ 0.5_dp*fac_h(ia)*Vxc_mn_DFT(m,n)
            enddo
         enddo
      enddo

   endif !(Functional == 'SLR')


   if((Functional == 'DFT') .or.  (Functional == 'HDR') .or. (Functional == 'DBF') ) then
      call orbs_on_grid(maxorb)


      call DFT_functional_libxc(occnum,orbgridsq)

      if(.not.allocated(Vxc_mn_DFT)) allocate (Vxc_mn_DFT(nbasis,nbasis))
      if(.not.allocated(Vc_mn_DFT)) allocate (Vc_mn_DFT(nbasis,nbasis))

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,n,ig)
!$OMP DO
      do m=1,nbasis
         do n=1,m
            Vxc_mn_DFT(m,n)=0.0_dp; Vc_mn_DFT(m,n)=0.0_dp
            do ig=1,ngrid
               Vxc_mn_DFT(m,n) = Vxc_mn_DFT(m,n) + w_grid(ig)*bas_f_grid(ig,m)*Vlxc(1,ig)*bas_f_grid(ig,n)
               Vc_mn_DFT(m,n) = Vc_mn_DFT(m,n) + w_grid(ig)*bas_f_grid(ig,m)*Vlc(1,ig)*bas_f_grid(ig,n)
            enddo
            Vxc_mn_DFT(n,m)=Vxc_mn_DFT(m,n); Vc_mn_DFT(n,m)=Vc_mn_DFT(m,n)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      if(gga_x .or. gga_c) then
         do m=1,nbasis
            do n=1,nbasis
               DM(m,n)=0.0_dp
               do ia=1,nnatorb
                  DM(m,n)=DM(m,n)+occnum(ia,3)*vecnat(m,ia)*vecnat(n,ia) 
               enddo
            enddo
         enddo

         do ig=1,ngrid
            do m=1,nbasis
               do n=1,nbasis
                  do ic=1,3
                     DXX(m,n,ic)=grad_bas_f_grid(ig,m,ic)*bas_f_grid(ig,n)
                  enddo
               enddo
            enddo

            do ic=1,3
               DMX(ic) = 0.d0
               do m=1,nbasis
                  do n=1,nbasis
                     DMX(ic)=DMX(ic) + DM(m,n)*DXX(m,n,ic)
                  enddo
               enddo
               DMX(ic) = 4.0_dp*DMX(ic) 
            enddo

            do k=1,nbasis
               do l=1,k
                  ss=0.0_dp
                  do ic=1,3
                     ss=ss+DMX(ic)*(DXX(k,l,ic)+DXX(l,k,ic))
                  enddo
                  if(gga_x) Vxc_mn_DFT(k,l) = Vxc_mn_DFT(k,l)+w_grid(ig)*Vlxc_sig(1,ig)*ss
                  if(gga_c) Vc_mn_DFT(k,l) = Vc_mn_DFT(k,l)+w_grid(ig)*Vlc_sig(1,ig)*ss
               enddo
            enddo
         enddo
         do k=1,nbasis
            do l=1,k
               Vxc_mn_DFT(l,k) = Vxc_mn_DFT(k,l); Vc_mn_DFT(l,k) = Vc_mn_DFT(k,l)
            enddo
         enddo
      endif !(gga_x .or. gga_c)

      open(unit=121,file='Vxc', status='unknown')
         write(121,'(2i4,f20.10)') ((l,k,Vxc_mn_DFT(k,l)+Vc_mn_DFT(k,l),k=1,l),l=1,nbasis)
      close(121)

      if (Functional/='HDR') then
         do ia=1,maxorb
            do m=1,nbasis
               do n=1,m
                 F(m,n,ia) = F(m,n,ia)+ 0.5_dp*fac_h(ia)*(Vxc_mn_DFT(m,n)+Vc_mn_DFT(m,n))
                 F(n,m,ia) = F(m,n,ia)
               enddo
            enddo
         enddo

!        if (.not. do_oep) then
!           if(Fexist) then
!              F=xmix_dft*F+(1.d0-xmix_dft)*F_old !Mixing. For OEP its already done
!           else
!              if ( .not. allocated(F_old) ) allocate ( F_old(lnbasis,lnbasis,lnnatorb) )
!              Fexist=.true.
!           endif
!           F_old=F
!        endif
      else
         do ia=1,maxorb
            do m=1,nbasis
               do n=1,m
                 F(m,n,ia) = F(m,n,ia)+ 0.5_dp*fac_h(ia)*xmix_DF*(Vxc_mn_DFT(m,n) + Vc_mn_DFT(m,n))
               enddo
            enddo
         enddo
      endif

   endif !((Functional == 'DFT' .or.  (Functional == 'HDR'))

end subroutine construct_f
!--------------------------------------------------------------------------

      
subroutine sum_intg()
!--------------------------------------------------------------------------
! This subroutine does the dirty and time consuming job of suming over 
! the 2-e integrals for the calculation of Matrix elements of coulomb 
! and exchange kind. Most of the computer time is 'wasted' in this routine.
! Although the best is done to make it optimal there might still  be 
! space of improvement. BUT DO WITH CARE AND ONLY IF YOU KNOW WHAT YOU
! ARE DOING.
!--------------------------------------------------------------------------

!..Global
   use global; use integrals; use matrices 
   use orbocc, ONLY: vecnat
   use functional_m, ONLY: maxorb
   implicit none

!..Local variables
   integer :: m, n, l, s
   integer(8) :: ind_p
   integer :: ia, ib, nrec, nchunk, nrec_last, ichunk, irec
   logical :: inv_pairs, read_int
   real(dp) :: two_int, wmn, wls
   real(dp), allocatable :: DM(:,:,:)

   if(.not.allocated(DM) ) allocate( DM(lnbasis,lnbasis,lnnatorb) )
   
   
   do ia=1,maxorb
      do m = 1,nbasis
         do n =1, m
            DM(m,n,ia)=vecnat(m,ia)*vecnat(n,ia)
            DM(n,m,ia)=DM(m,n,ia)
         enddo
      enddo
   enddo

   nrec=50000000 ! x20 = 1GB 
   nchunk=nintgr/nrec
   nchunk=nchunk+1
   nrec_last=mod(nintgr,nrec)

   read_int=.true.
   if(nchunk==1) then
      nrec=nrec_last
      read_int=.false.
   endif
   
   if(.not.allocated(twoin)) then
      allocate( indpacked(nrec), mu(nrec), nu(nrec), lambda(nrec), sigma(nrec), twoin(nrec) )
      read_int=.true.
   endif

   Coul=0.0_dp; Exch=0.0_dp

   rewind(42)
!..Do loop over the nonzero and unique orbitals
   do ichunk=1,nchunk
      if(ichunk==nchunk) nrec=nrec_last

      if ( read_int ) then
         read(42,err=140) (indpacked(irec), twoin(irec), irec=1,nrec)
      endif
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(irec,ind_p,m,n,l,s)
!$OMP DO
      do irec = 1, nrec
         ind_p=indpacked(irec)
         call unpack4(m,n,l,s,ind_p)
         mu(irec)=m
         nu(irec)=n
         lambda(irec)=l 
         sigma(irec)=s 
      enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ib,irec,ind_p,m,n,l,s,two_int,wmn,wls,inv_pairs)
!$OMP DO
   do ib=1,maxorb
      do irec = 1, nrec
         m=mu(irec)
         n=nu(irec)
         l=lambda(irec)
         s=sigma(irec)

         two_int = twoin(irec)

!.....Create the weights due to symmetry of (mn|ls)
         if(m==n) then
            wmn = 1.0_dp
         else
            wmn = 2.0_dp
         endif
         if(l==s) then
            wls = 1.0_dp
         else
            wls = 2.0_dp
         endif

!.....The coulomb and exchange should be added to different F elements:
!.....Coulomb
         inv_pairs = (l/=m.or.n/=s)
         Coul(m,n,ib) = Coul(m,n,ib) + DM(l,s,ib) * wls * two_int

         if(inv_pairs) Coul(l,s,ib) = Coul(l,s,ib) + DM(m,n,ib) * wmn * two_int

!.....Exchange
         if(m>=l) Exch(m,l,ib) = Exch(m,l,ib) + two_int * DM(n,s,ib)
         if(inv_pairs.and.l>=m) Exch(l,m,ib) = Exch(l,m,ib) + two_int * DM(n,s,ib)
         if(m/=n) then
            if(n>=l) Exch(n,l,ib) = Exch(n,l,ib) + two_int * DM(m,s,ib)
            if(inv_pairs.and.l>=n) Exch(l,n,ib) = Exch(l,n,ib) + two_int * DM(m,s,ib)

            if(l/=s) then
               if(n>=s) Exch(n,s,ib) = Exch(n,s,ib) + two_int * DM(m,l,ib)
               if(inv_pairs.and.s>=n) Exch(s,n,ib)=Exch(s,n,ib) + two_int * DM(m,l,ib)
            endif !(l/=s)
         endif !(m/=n)

         if(l/=s) then
            if(m>=s) Exch(m,s,ib) = Exch(m,s,ib) + two_int * DM(n,l,ib)
            if(inv_pairs.and.s>=m) Exch(s,m,ib) = Exch(s,m,ib) + two_int * DM(n,l,ib)
         endif !l/=s
            
      enddo ! irec
   enddo !ib
!$OMP END DO
!$OMP END PARALLEL

   enddo ! ichunk

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,m,n)
!$OMP DO
   do ia=1,maxorb
      do m=1,nbasis
         do n=1,m
            Coul(n,m,ia)=Coul(m,n,ia)
            Exch(n,m,ia)=Exch(m,n,ia)
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   deallocate( DM )

   return

 140 stop 'construct_f:sum_intg: error reading 2-e intergral file'

end subroutine sum_intg

!=====================================================================================
subroutine DFT_functional_libxc(occ,orbgridsq)
   use global; use grid_params, ONLY: ngrid, chdens_gr
   use DFT_par; use functional_m, ONLY:Functional, cl_shell
!  use xc_f90_types_m; use xc_f90_lib_m
   implicit none

!..Arguments
   real(dp) :: occ(lnnatorb,3), orbgridsq(ngrid,lnnatorb)
   
!..Local
   real(dp), allocatable :: rho(:,:), grad_n(:,:)
   integer :: ig, ia, i1

   if(.not.allocated(chdens_gr)) allocate (chdens_gr(ngrid,3))
   
   if ( cl_shell ) then
      i1=1
   else
      i1=2
   endif
   allocate ( rho(i1,ngrid), grad_n(i1,ngrid) )
   rho=0.0_dp

!$omp parallel default(shared) &
!$omp private(ig,ia)
!$omp do
    do ig=1,ngrid
      do ia=1,nnatorb
         if (cl_shell) then
            rho(1,ig) = rho(1,ig) + occ(ia,3)*orbgridsq(ig,ia)
         else
            rho(1,ig) = rho(1,ig) + occ(ia,1)*orbgridsq(ig,ia)
            rho(2,ig) = rho(2,ig) + occ(ia,2)*orbgridsq(ig,ia)
         endif
      enddo
   enddo
!$omp end do
!$omp end parallel

   if ( cl_shell ) then
      chdens_gr(:,3)=rho(1,:)
      chdens_gr(:,1)=0.5_dp*chdens_gr(:,3)
      chdens_gr(:,2)=chdens_gr(:,1)
   else
      chdens_gr(:,1)=rho(1,:)
      chdens_gr(:,2)=rho(2,:)
      chdens_gr(:,3)=chdens_gr(:,1)+chdens_gr(:,2)
   endif

   if (.not. allocated(Exc_DFT)) allocate(Exc_DFT(ngrid))
   if (.not. allocated(Ec_DFT)) allocate(Ec_DFT(ngrid))

   if (.not. allocated(Vlxc)) allocate ( Vlxc(i1,ngrid) )
   if (.not. allocated(Vlc)) allocate ( Vlc(i1,ngrid) )
   if (.not. allocated(Vlxc_sig)) allocate ( Vlxc_sig(i1,ngrid) )
   if (.not. allocated(Vlc_sig)) allocate ( Vlc_sig(i1,ngrid) )

!  allocate (Fxc_DFT(ngrid), Fc_DFT(ngrid))

!..Exchange Correlation Functional
   select case (xc_f90_info_family(xc_info))
   case(XC_FAMILY_LDA)
      call xc_f90_lda_exc_vxc(xc_func, ngrid, rho(1,1), Exc_DFT(1), Vlxc(1,1))
   case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call density_gradient( i1, grad_n )
      call xc_f90_gga_exc_vxc(xc_func, ngrid, rho(1,1), grad_n(1,1), Exc_DFT(1), Vlxc(1,1), Vlxc_sig(1,1) )
      gga_x=.true.
   case default
      Exc_DFT=0.0_dp; Vlxc=0.0_dp; Vlxc_sig=0.0_dp
   end select

!..Correlation Functional for seperate additional correlation functional
   if (id_func_c/=0 .and. Functional /= 'SLR') then
      select case (xc_f90_info_family(c_info))
      case(XC_FAMILY_LDA)
         call xc_f90_lda_exc_vxc(c_func, ngrid, rho(1,1), Ec_DFT(1), Vlc(1,1))
      case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA) 
         if (.not. gga_x) call density_gradient ( i1, grad_n )
         call xc_f90_gga_exc_vxc(c_func, ngrid, rho(1,1), grad_n(1,1), Ec_DFT(1), Vlc(1,1), Vlc_sig(1,1) )
         gga_c=.true.
      case default
         Ec_DFT=0.0_dp; Vlc=0.0_dp; Vlc_sig=0.0_dp
      end select
   else
      Ec_DFT=0.0_dp; Vlc=0.0_dp; Vlc_sig=0.0_dp 
   endif

   

end subroutine DFT_functional_libxc
!==========================================================================================
!-----------------------------------------------------------------------------------
  subroutine density_gradient (i1, top_orb, grad_n)
   use global; use orbocc; use grid_params; implicit none 

!..Arguments
   integer, intent(in) :: i1,top_orb
   real(dp),intent(out) :: grad_n(2*i1-1,ngrid)

!..Local 
   integer :: ibs, ig, k, m, ia, ispin, ic
   real(dp) :: gx, gy, gz, DM(lnbasis, lnbasis, 2), g_ig(3,3)

!..Gradient of Basis Functions on grid
   if (.not.allocated(grad_bas_f_grid)) then
      allocate(grad_bas_f_grid(ngrid,nbasis,3))
   endif
   if (.not. lgradset) then
      do ibs=1,nbasis
         do ig=1,ngrid
            call gradient_f_bas(ibs, x_grid(ig), y_grid(ig), z_grid(ig), gx, gy, gz)
            grad_bas_f_grid(ig,ibs,1)=gx
            grad_bas_f_grid(ig,ibs,2)=gy
            grad_bas_f_grid(ig,ibs,3)=gz
         enddo
      enddo
      lgradset=.true.
   endif

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,m,ispin,ia)
!$OMP DO
   do k=1,nbasis
      do m=1,nbasis
         do ispin=1,2
            DM(k,m,ispin) = 0.0_dp
            do ia=1,top_orb
               if( occnum(ia,ispin) > small) then
                  DM(k,m,ispin) = DM(k,m,ispin) + occnum(ia,ispin)*vecnat(k,ia)*vecnat(m,ia)
               endif
            enddo
         enddo  
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,ic,ispin,g_ig,k,m)
!$OMP DO
   do ig=1,ngrid
      do ic=1,3
         do ispin=1,2
            g_ig(ic,ispin) = 0.0_dp  
            do k=1,nbasis
               do m=1,nbasis
                  g_ig(ic,ispin) = g_ig(ic,ispin) + DM(k,m,ispin)*grad_bas_f_grid(ig,k,ic)*bas_f_grid(ig,m)
               enddo
            enddo
            g_ig(ic,ispin) = 2.0_dp*g_ig(ic,ispin) 
         enddo
      enddo

      !Construct the 'contracted gradients' required for libxc input
      if (i1==1) then
         g_ig(:,3) = g_ig(:,1) + g_ig(:,2)
         grad_n(1,ig) = dot_product(g_ig(:,3),g_ig(:,3))
      else
        grad_n(1,ig) = dot_product(g_ig(:,1),g_ig(:,1))
        grad_n(2,ig) = dot_product(g_ig(:,1),g_ig(:,2))
        grad_n(3,ig) = dot_product(g_ig(:,2),g_ig(:,2))
      endif
   enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine density_gradient

!--------------------------------------------------------------------------------------------------------------
subroutine density_gradient_OLD (i1, grad_n )
   use global; use orbocc; use grid_params; use functional_m, ONLY:cl_shell; implicit none 

!..Arguments
   integer, intent(in) :: i1
   real(dp),intent(out) :: grad_n(i1,ngrid)

!..Local 
   integer :: ibs, ig, k, m, ia, ispin, ic
   real(dp) :: gx, gy, gz, DM(lnbasis, lnbasis, 2), g_ig(3,2)

!..Gradient of Basis Functions on grid
   if (.not.allocated(grad_bas_f_grid)) then
      allocate(grad_bas_f_grid(ngrid,nbasis,3))
   endif
   if (.not. lgradset) then
      do ibs=1,nbasis
         do ig=1,ngrid
            call gradient_f_bas(ibs, x_grid(ig), y_grid(ig), z_grid(ig), gx, gy, gz)
            grad_bas_f_grid(ig,ibs,1)=gx
            grad_bas_f_grid(ig,ibs,2)=gy
            grad_bas_f_grid(ig,ibs,3)=gz
         enddo
      enddo
      lgradset=.true.
   endif

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,m,ispin,ia)
!$OMP DO
   do k=1,nbasis
      do m=1,nbasis
         do ispin=1,2
            DM(k,m,ispin) = 0.0_dp
            do ia=1,nnatorb
               if( occnum(ia,ispin) > small) then
                  DM(k,m,ispin) = DM(k,m,ispin) + occnum(ia,ispin)*vecnat(k,ia)*vecnat(m,ia)
               endif
            enddo
         enddo  
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,ic,ispin,g_ig,k,m)
!$OMP DO
   do ig=1,ngrid
      do ic=1,3
         do ispin=1,2
            g_ig(ic,ispin) = 0.0_dp  
            do k=1,nbasis
               do m=1,nbasis
                  g_ig(ic,ispin) = g_ig(ic,ispin) + DM(k,m,ispin)*grad_bas_f_grid(ig,k,ic)*bas_f_grid(ig,m)
               enddo
            enddo
            g_ig(ic,ispin) = 4.0_dp*g_ig(ic,ispin) 
         enddo
      enddo
      if (cl_shell) then
         grad_n(1,ig) = g_ig(1,1)*g_ig(1,1)+g_ig(2,1)*g_ig(2,1)+g_ig(3,1)*g_ig(3,1)
      else
!   To be checked
         stop 'construct f: check first'
!        grad_n(1,ig) = g_ig(1,1)*g_ig(1,2)+g_ig(2,1)*g_ig(2,2)+g_ig(3,1)*g_ig(3,2)
!        grad_n(2,ig) = g_ig(1,2)*g_ig(1,2)+g_ig(2,2)*g_ig(2,2)+g_ig(3,2)*g_ig(3,2)
      endif
   enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine density_gradient_OLD

!==========================================================================================
subroutine density_gradient_select(Npoints,Xpoint,Ypoint,Zpoint,Gpoint,G_c)
   use global; use orbocc; use grid_params; implicit none 

!..Arguments
   integer :: Npoints
   real(dp) :: Xpoint(Npoints), Ypoint(Npoints), Zpoint(Npoints), Gpoint(Npoints), G_c(Npoints,3)

!..Local 
   integer :: ibs, ig, k, m, ia, ispin, ic
   real(dp) :: gx, gy, gz, DM(lnbasis, lnbasis, 2), g_ig(3,2)
   real(dp) :: gbf(nbasis,3), bf(nbasis)
   real(dp), external :: f_bas


!..Gradient of Basis Functions 

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,m,ispin,ia)
!$OMP DO
   do k=1,nbasis
      do m=1,nbasis
         do ispin=1,2
            DM(k,m,ispin) = 0.0_dp
            do ia=1,nnatorb
               if( occnum(ia,ispin) > small) then
                  DM(k,m,ispin) = DM(k,m,ispin) + occnum(ia,ispin)*vecnat(k,ia)*vecnat(m,ia)
               endif
            enddo
         enddo  
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,ic,ibs,bf,gbf,gx,gy,gz,ispin,g_ig,k,m)
!$OMP DO
   do ig=1,Npoints
      do ibs=1,nbasis
         call gradient_f_bas(ibs, Xpoint(ig), Ypoint(ig), Zpoint(ig), gx, gy, gz)
         gbf(ibs,1)=gx
         gbf(ibs,2)=gy
         gbf(ibs,3)=gz
         bf(ibs) = f_bas(ibs, Xpoint(ig), Ypoint(ig), Zpoint(ig))
      enddo
      grad_bas_f_grid(ig,:,:) = gbf

      do ic=1,3
         do ispin=1,2
            g_ig(ic,ispin) = 0.0_dp  
            do k=1,nbasis
               do m=1,nbasis
                  g_ig(ic,ispin) = g_ig(ic,ispin) + DM(k,m,ispin)*gbf(k,ic)*bf(m)
               enddo
            enddo
            g_ig(ic,ispin) = 2.0_dp*g_ig(ic,ispin) 
         enddo
      enddo
!....The factors of 2.d0 and 4.0 in the following is because we want total gradient and there is spin degeneracy
      G_c(ig,1)=2.0_dp*g_ig(1,1); G_c(ig,2)=2.0_dp*g_ig(2,1); G_c(ig,3)=2.0_dp*g_ig(3,1)
      Gpoint(ig) = 4.0_dp*(g_ig(1,1)*g_ig(1,1)+g_ig(2,1)*g_ig(2,1)+g_ig(3,1)*g_ig(3,1))
!     grad_n(2,ig) = g_ig(1,1)*g_ig(1,2)+g_ig(2,1)*g_ig(2,2)+g_ig(3,1)*g_ig(3,2)
!     grad_n(3,ig) = g_ig(1,2)*g_ig(1,2)+g_ig(2,2)*g_ig(2,2)+g_ig(3,2)*g_ig(3,2)
   enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine density_gradient_select
!=====================================================================================
subroutine DFT_functional_mine(occ,orbgridsq)
!..Made for testing only 
   use global; use grid_params, ONLY:ngrid; use DFT_par
   implicit none

!..Arguments
   real(dp) :: occ(lnnatorb), orbgridsq(ngrid,lnnatorb)

   integer :: ig, ia
   real(dp) :: rho, rho1, rho2, vx_1, vx_2, vc_1, vc_2, ex, ec, zero2

   zero2=zero*zero

   if (.not. allocated(Exc_DFT)) allocate(Exc_DFT(ngrid))
   if (.not. allocated(Ec_DFT)) allocate(Ec_DFT(ngrid))
   if (.not. allocated(Vlxc)) allocate(Vlxc(1,ngrid))
   if (.not. allocated(Vlc)) allocate(Vlc(1,ngrid))

   do ig=1,ngrid
      rho = 0.0_dp
      do ia=1,nnatorb
         rho = rho + occ(ia)*orbgridsq(ig,ia)
      enddo
      rho1=rho/2.0_dp
      rho2=rho1
      rho1=max(zero2,rho1)
      rho2=max(zero2,rho2)

      call lsd_slater(rho1,rho2,vx_1,vx_2, ex)
      vc_1=0.0_dp; vc_2=0.0_dp; ec=0.0_dp
      call lsd_VWN_cor(rho1,rho2,vc_1,vc_2,ec)

      Exc_DFT(ig) = ex/rho; Ec_DFT(ig)=ec/rho
      Vlxc(1,ig) = vx_1; Vlc(1,ig)=vc_1 ! For open shell Vx_DFT, Vc_DFT need a second index
   enddo

end subroutine DFT_functional_mine

subroutine lsd_slater(rho1, rho2, vx_1, vx_2, ex)
!..Slater Exchange for testing
   use global; implicit none
   real(dp) :: rho1, rho2, vx_1, vx_2, ex
   
   real(dp) :: rhoth1, rhoth2, fpot, fexc, third, alp

   alp=2.0_dp/3.0_dp
   third=1.0_dp/3.0_dp
   rhoth1=rho1**third; rhoth2=rho2**third
   fpot=(3.0_dp/4.0_dp/pi)**third
   fexc=-(9.0_dp/4.0_dp)*alp*fpot
   fpot=-3.0_dp*alp*fpot
   ex=fexc*(rhoth1**4.0_dp + rhoth2**4.0_dp)
   vx_1 = fpot*rhoth1
   vx_2 = fpot*rhoth2
end subroutine lsd_slater

subroutine lsd_VWN_cor(rho1,rho2,vc_1,vc_2,ec)
!..LDA VWN correlation for testing 
!..Global
   use global; implicit none
!..Arguments
   real(dp), intent(IN) :: rho1, rho2
   real(dp), intent(OUT) :: vc_1, vc_2, ec
 
!..Parameters
   real(dp), parameter :: facr = 0.78762331789974325053_dp, &
                                 fach = 1.70992093416136561756_dp,&
                                 facg = 1.125d0, facdg=1.5d0,&
                                 VWNAP= 0.0621814_dp/2.0_dp, VWNBP= 3.7274400_dp, &
                                 VWNCP=12.9352000D+00, VWNX0P=-0.1049800D+00, &
                                 VWNAF= VWNAP/2.0_dp, VWNBF= 7.0604200_dp, &
                                 VWNCF=18.0578000_dp, VWNX0F=-0.3250000_dp, &
                                 VWNAA=-0.03377372_dp/2.0_dp, VWNBA= 1.1310700_dp, &
                                 VWNCA=13.0045000_dp, VWNX0A=-0.0047584_dp

!..Local
   real(dp) :: rhoa, rhob, rhot, rho16, qp, qf, qa, xx0p, xx0f, xx0a
   real(dp) :: coef1p, coef1f, coef1a, coef2p, coef2f, coef2a
   real(dp) :: coef3p, coef3f, coef3a, x, xxp, xxf, xxa, t1p, t1f, t1a 
   real(dp) :: t2p, t2f, t2a, t3p, t3f, t3a, epscp, epscf, epsca
   real(dp) :: xbp, xbf, xba
   real(dp) :: xbq2p,xbq2f,xbq2a,brackp,brackf,bracka,curlyp,curlyf,curlya,depscp,depscf,depsca
   real(dp) :: hx,dhx, zeta, zeta2, zeta3, zeta4, gzeta, dgzeta, dpot1, dpot2, duma, dumb, pot

   rhoa  =rho1
   rhob  =rho2
   rhot  =rhoa+rhob
   rho16 =rhot**(1._dp/6._dp)

   qp    = sqrt(4._dp*vwncp-vwnbp*vwnbp)
   qf    = sqrt(4._dp*vwncf-vwnbf*vwnbf)
   qa    = sqrt(4._dp*vwnca-vwnba*vwnba)
   xx0p  = VWNX0P*VWNX0P+VWNX0P*VWNBP+VWNCP
   xx0f  = VWNX0F*VWNX0F+VWNX0F*VWNBF+VWNCF
   xx0a  = VWNX0A*VWNX0A+VWNX0A*VWNBA+VWNCA
 
   COEF1P= VWNAP
   COEF1F= VWNAF
   COEF1A= VWNAA
   COEF2P=-VWNAP*VWNBP*VWNX0P/XX0P
   COEF2F=-VWNAF*VWNBF*VWNX0F/XX0F
   COEF2A=-VWNAA*VWNBA*VWNX0A/XX0A
   COEF3P= VWNAP*(2._dp*VWNBP/QP)*((VWNCP-VWNX0P*VWNX0P)/XX0P)
   COEF3F= VWNAF*(2._dp*VWNBF/QF)*((VWNCF-VWNX0F*VWNX0F)/XX0F)
   COEF3A= VWNAA*(2._dp*VWNBA/QA)*((VWNCA-VWNX0A*VWNX0A)/XX0A)
 
   X     =FACR/RHO16
   XXP   =X*X+X*VWNBP+VWNCP
   XXF   =X*X+X*VWNBF+VWNCF
   XXA   =X*X+X*VWNBA+VWNCA
   T1P   =LOG(X*X/XXP)
   T1F   =LOG(X*X/XXF)
   T1A   =LOG(X*X/XXA)
   T2P   =LOG((X-VWNX0P)*(X-VWNX0P)/XXP)
   T2F   =LOG((X-VWNX0F)*(X-VWNX0F)/XXF)
   T2A   =LOG((X-VWNX0A)*(X-VWNX0A)/XXA)
   T3P   =ATAN(QP/(2._dp*X+VWNBP))
   T3F   =ATAN(QF/(2._dp*X+VWNBF))
   T3A   =ATAN(QA/(2._dp*X+VWNBA))
 
   EPSCP =COEF1P*T1P+COEF2P*T2P+COEF3P*T3P
   EPSCF =COEF1F*T1F+COEF2F*T2F+COEF3F*T3F
   EPSCA =COEF1A*T1A+COEF2A*T2A+COEF3A*T3A
 
   XBP    = 2._dp*X+VWNBP
   XBF    = 2._dp*X+VWNBF
   XBA    = 2._dp*X+VWNBA
   XBQ2P  = XBP*XBP + QP*QP
   XBQ2F  = XBF*XBF + QF*QF
   XBQ2A  = XBA*XBA + QA*QA
   BRACKP = 2._dp/(X-VWNX0P) -XBP/XXP -4._dp*(2._dp*VWNX0P+VWNBP)/XBQ2P
   BRACKF = 2._dp/(X-VWNX0F) -XBF/XXF -4._dp*(2._dp*VWNX0F+VWNBF)/XBQ2F
   BRACKA = 2._dp/(X-VWNX0A) -XBA/XXA -4._dp*(2._dp*VWNX0A+VWNBA)/XBQ2A
   CURLYP = 2._dp/X - XBP/XXP - 4._dp*VWNBP/XBQ2P
   CURLYF = 2._dp/X - XBF/XXF - 4._dp*VWNBF/XBQ2F
   CURLYA = 2._dp/X - XBA/XXA - 4._dp*VWNBA/XBQ2A
   DEPSCP = VWNAP*(CURLYP - VWNBP*VWNX0P*BRACKP/XX0P)
   DEPSCF = VWNAF*(CURLYF - VWNBF*VWNX0F*BRACKF/XX0F)
   DEPSCA = VWNAA*(CURLYA - VWNBA*VWNX0A*BRACKA/XX0A)
 
   HX  = FACH*(EPSCF-EPSCP)/EPSCA - 1._dp
   DHX = FACH*(DEPSCF-DEPSCP-(EPSCF-EPSCP)*DEPSCA/EPSCA)/EPSCA

   ZETA  = (RHOA-RHOB)/RHOT
   ZETA2 = ZETA *ZETA
   ZETA3 = ZETA2*ZETA
   ZETA4 = ZETA2*ZETA2
   GZETA = ((1._dp+ZETA)**(4._dp/3._dp)+(1._dp-ZETA)**(4._dp/3._dp)-2._dp)*FACG
   DGZETA= ((1._dp+ZETA)**(1._dp/3._dp) -(1._dp-ZETA)**(1._dp/3._dp) )*FACDG
 
!          TOTAL POTENTIAL, AND CONTRIBUTION TO ENERGY
 
   POT   = EPSCP + EPSCA*GZETA*(1._dp+HX*ZETA4)
   ec  = POT*RHOT
 
!          CONTRIBUTION TO FUNCTIONAL DERIVATIVE
 
   DPOT1 = -(X/6._dp)*(DEPSCP + DEPSCA*GZETA*(1._dp+HX*ZETA4)&
                            +  EPSCA*GZETA*DHX*ZETA4)
   DPOT2 = EPSCA*(DGZETA*(1._dp+HX*ZETA4) + 4._dp*GZETA*HX*ZETA3)
   DUMA  = POT + (DPOT1 + DPOT2*(1._dp-ZETA))
   DUMB  = POT + (DPOT1 - DPOT2*(1._dp+ZETA))
 
   vc_1  = vc_1 + DUMA
   vc_2  = vc_2 + DUMB

end subroutine lsd_VWN_cor

!-------------------------------------------------------------------------------------
! Calculates the integral \int d^3x d^3y \frac{\phi_a(x)\phi_b(x)\phi_c(y)\phi_d(y)}
! {|x-y|} where x,y are 3d space vectors and \phi_a, \phi_b ... are molecular orbitals
! 
subroutine psi4_int(phi_a,phi_b,phi_c,phi_d,X_int)
!..Global
   use global; use integrals
   implicit none

!..Arguments
   real(dp), intent(in) :: phi_a(lnbasis), phi_b(lnbasis), phi_c(lnbasis), phi_d(lnbasis)
   real(dp), intent(out) :: X_int

!..Local 
   integer :: k,l,m,n, nrec, nchunk, nrec_last, ichunk, irec
   integer(8) :: ind_p
   logical :: inv_pairs, read_int
   real(dp) :: two_int

   X_int=0._dp

   nrec=50000000 ! x20 = 1GB 
   nchunk=nintgr/nrec
   nchunk=nchunk+1
   nrec_last=mod(nintgr,nrec)

   read_int=.true.
   if(nchunk==1) then
      nrec=nrec_last
      read_int=.false.
   endif
   
   if(.not.allocated(twoin)) then
      allocate( indpacked(nrec), mu(nrec), nu(nrec), lambda(nrec), sigma(nrec), twoin(nrec) )
      read_int=.true.
   endif

   rewind(42)
!..Do loop over the nonzero and unique orbitals
   do ichunk=1,nchunk
      if(ichunk==nchunk) nrec=nrec_last

      if ( read_int ) then
         read(42,err=149) (indpacked(irec), twoin(irec), irec=1,nrec)
      endif
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(irec,ind_p,k,l,m,n)
!$OMP DO
      do irec = 1, nrec
         ind_p=indpacked(irec)
         call unpack4(k,l,m,n,ind_p)
         mu(irec)=k
         nu(irec)=l
         lambda(irec)=m 
         sigma(irec)=n 
      enddo
!$OMP END DO
!$OMP END PARALLEL

      do irec = 1, nrec
         k=mu(irec)
         l=nu(irec)
         m=lambda(irec)
         n=sigma(irec)

         two_int = twoin(irec)

!        X_int=X_int+two_int*vecnat(k,ia)*vecnat(l,ib)*vecnat(m,ic)*vecnat(n,id)

!        if(k /= l) X_int=X_int+two_int*vecnat(l,ia)*vecnat(k,ib)*vecnat(m,ic)*vecnat(n,id)
!        if(m /= n) X_int=X_int+two_int*vecnat(k,ia)*vecnat(l,ib)*vecnat(n,ic)*vecnat(m,id)
!        if(k /= l .and. m /= n) X_int=X_int+two_int*vecnat(l,ia)*vecnat(k,ib)*vecnat(n,ic)*vecnat(m,id)

!        inv_pairs= (k/=m.or.l/=n)
!        if (inv_pairs) then
!           X_int=X_int+two_int*vecnat(m,ia)*vecnat(n,ib)*vecnat(k,ic)*vecnat(l,id)
!           if(m /= n) X_int=X_int+two_int*vecnat(n,ia)*vecnat(m,ib)*vecnat(k,ic)*vecnat(l,id)
!           if(k /= l) X_int=X_int+two_int*vecnat(m,ia)*vecnat(n,ib)*vecnat(l,ic)*vecnat(k,id)
!           if(k /= l .and. m /= n) X_int=X_int+two_int*vecnat(n,ia)*vecnat(m,ib)*vecnat(l,ic)*vecnat(k,id)
!        endif

         X_int=X_int+two_int*phi_a(k)*phi_b(l)*phi_c(m)*phi_d(n)
     
         if(k /= l) X_int=X_int+two_int*phi_a(l)*phi_b(k)*phi_c(m)*phi_d(n)
         if(m /= n) X_int=X_int+two_int*phi_a(k)*phi_b(l)*phi_c(n)*phi_d(m)
         if(k /= l .and. m /= n) X_int=X_int+two_int*phi_a(l)*phi_b(k)*phi_c(n)*phi_d(m)

         inv_pairs= (k/=m.or.l/=n)
         if (inv_pairs) then
            X_int=X_int+two_int*phi_a(m)*phi_b(n)*phi_c(k)*phi_d(l)
            if(m /= n) X_int=X_int+two_int*phi_a(n)*phi_b(m)*phi_c(k)*phi_d(l)
            if(k /= l) X_int=X_int+two_int*phi_a(m)*phi_b(n)*phi_c(l)*phi_d(k)
            if(k /= l .and. m /= n) X_int=X_int+two_int*phi_a(n)*phi_b(m)*phi_c(l)*phi_d(k)
         endif
      enddo ! irec
   enddo ! ichunk

   return

 149 stop 'psi4_int: error reading 2-e intergral file'

end subroutine psi4_int

