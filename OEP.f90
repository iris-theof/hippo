subroutine analytic_Veff(iscf)
!..Global
   use global; use grid_params; use basis_set; use orbocc; use matrices
   use energies, ONLY:E_HOMO; use functional_m; use DFT_par
   implicit none

!..Argument
   integer :: iscf

!..Local
   logical :: automatic_limit
   logical,allocatable  :: neg_scdens(:)
   integer :: k, l, m, ig, igp, ia, ii, n_occ, info
   real(dp),allocatable :: Pkn(:,:), A_kn(:,:), Qkn(:,:), DMMM(:,:) 
   real(dp) :: B_kn(lnbasis_pot), R_kn(lnbasis_pot), S_kn(lnbasis_pot)
   real(dp),allocatable :: Akn_inv(:,:), FF(:,:), Fii(:,:), FFp(:,:)
   real(dp) :: SAinvS, SAinvB, AinvS(lnbasis_pot), xmu_oep, fc_mix, rr
   real(dp),allocatable :: ovlap3nat(:,:,:), A_tild(:,:)
   real(dp) :: DDD(lnbasis,lnbasis)
   real(dp) :: B_tild(lnbasis_pot), B_kn_in(lnbasis_pot)
   real(dp) :: dm_gr, pos, pos_old, fac_dif
   real(dp), allocatable :: Schdens(:)
   real(dp) :: Uhx(ngrid), Uxx(ngrid)
   real(dp) :: Delta_e, Delta_n, AinvB(lnbasis_pot)
   real(dp) :: X_bar(lnbasis_pot), uxxtmp, uhxtmp, ss, volume_neg, fc_pen_v
   real(dp), external :: Dn_De_FD
   real(dp) :: DV_bs(lnbasis_pot), V_bs_in(lnbasis_pot)
   real(dp), save :: pos_mix1
   real(dp) :: fc_neg
   integer :: ipos, Nipos=6000

   automatic_limit = .false. !Automatic lambda -> 0
   do_ceda = .false. !Common energy denomicator approximation

   n_occ=max(ibond(1),ibond(2)) ! HOMO index

   if (automatic_limit .and. do_ceda) then
      svd_cut=1d10
   endif

   if (.not.allocated(vec_nat_grid)) allocate ( vec_nat_grid(ngrid,lnbasis) )
   if(.not.allocated(chdens_gr)) allocate(chdens_gr(ngrid,3))
   if ( .not.allocated(Vcoul) ) allocate( Vcoul(lnbasis,lnbasis) )
   if ( .not.allocated(Vexch) ) allocate( Vexch(lnbasis,lnbasis) )
   allocate (ovlap3nat(lnbasis, lnbasis, lnbasis_pot ), DMMM(nbasis,nbasis) )
   allocate (A_kn(lnbasis_pot,lnbasis_pot), A_tild(lnbasis_pot,lnbasis_pot))
   allocate (Akn_inv(lnbasis_pot,lnbasis_pot), FF(lnbasis, lnbasis), FFp(lnbasis, lnbasis) )
   allocate (Fii(lnbasis, lnbasis), Schdens(ngrid), neg_scdens(ngrid) ) 

   call orbs_on_grid(maxorb)

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,ii,ia,l,m,DDD)
!$OMP DO
   do k=1,nbasis_pot
      do ia=1,nbasis
         do l=1,nbasis
            DDD(l,ia)=0.d0
            do m=1,nbasis
               DDD(l,ia) = DDD(l,ia) + vecnat(m,ia)* ovlap3(l,m,k)
            enddo
         enddo
      enddo
      do ii=1,nbasis
         do ia=1,ii
            ovlap3nat(ii,ia,k)=0.d0
            do l=1,nbasis
               ovlap3nat(ii,ia,k)=ovlap3nat(ii,ia,k) + vecnat(l,ii)*DDD(l,ia)
            enddo
            ovlap3nat(ia,ii,k)=ovlap3nat(ii,ia,k)
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   if (Functional == 'RHF' .or. Functional == 'DFT' .or. Functional == 'RHA' ) then

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,l,ii,ia)
!$OMP DO
      do k=1,nbasis_pot
         do l=1,nbasis_pot
            A_kn(k,l)=0.d0
            do ii=1, n_occ
               do ia=n_occ+1,nbasis
                  A_kn(k,l) = A_kn(k,l) + occnum(ii,3)*ovlap3nat(ii,ia,k)*ovlap3nat(ia,ii,l) &
                            / max(ennat(ia)-ennat(ii), zero)
               enddo
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
   else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,l,ii,ia,Delta_n,Delta_e)
!$OMP DO
      do k=1,nbasis_pot
         do l=1,nbasis_pot
            A_kn(k,l)=0.d0
            do ii=1, nbasis
               do ia=1, nbasis; if (ia /= ii ) then
                  Delta_n=occnum(ii,3) - occnum(ia,3)
                  if ( abs(Delta_n) > 2.d0*Dn_lim ) then
                     Delta_e=ennat(ia) - ennat(ii) 
                     if(abs(Delta_e) < small_e) then
                     if ( Delta_e >= 0 ) then 
                         Delta_e=small_e
                     else
                         Delta_e=-small_e
                     endif
                     endif
                     A_kn(k,l) = A_kn(k,l) + ovlap3nat(ii,ia,k)*ovlap3nat(ia,ii,l) &
                               * Delta_n/Delta_e
                  endif
               endif; enddo
            enddo
            A_kn(k,l)=0.5d0*A_kn(k,l)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
   endif

   if (xld > zero .or. automatic_limit) call invert_Akn(A_kn, Akn_inv, svd_cut)

   if (Functional == 'RHF' .or. Functional == 'DFT' .or. Functional == 'RHA' ) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ii,ia,k,l)
!$OMP DO
     do ii=1, n_occ
        do ia=n_occ+1,nbasis
!     do ii=1, nbasis
!        do ia=1, nbasis
            FF(ia,ii) = 0.d0
            do k=1,nbasis
               do l=1,nbasis
                  if (l < k) then
                     FF(ia,ii) = FF(ia,ii) + vecnat(k,ia)*(F(k,l,1)-Hcore(k,l))*vecnat(l,ii)
                  else
                     FF(ia,ii) = FF(ia,ii) + vecnat(k,ia)*(F(l,k,1)-Hcore(l,k))*vecnat(l,ii)
                  endif
               enddo
            enddo
         enddo
      enddo   
!$OMP END DO
!$OMP END PARALLEL
   else !RDMFT
      do ii=1, nbasis
         fac_dif=0.5d0*fac_h(ii)
         do k=1,nbasis
            do l=1,k
               Fii(k,l)=F(k,l,ii)-fac_dif*hcore(k,l)-F_ha(k,l,ii)-F_x(k,l,ii)
            enddo
         enddo
         do ia=1, nbasis
            FFp(ia,ii)=0.d0
            do k=1,nbasis
               do l=1,nbasis
                  if(l<k) then
                     FFp(ia,ii) = FFp(ia,ii) &
                                + vecnat(k,ia)*Fii(k,l)*vecnat(l,ii)
                  else
                     FFp(ia,ii) = FFp(ia,ii) &
                                + vecnat(k,ia)*Fii(l,k)*vecnat(l,ii)
                  endif
               enddo
            enddo
         enddo
      enddo
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ii,ia)
!$OMP DO
      do ii=1,nbasis
         do ia=1,nbasis
            FF(ia,ii) = FFp(ia,ii)-FFp(ii,ia)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
   endif !(Functional == 'RHF' .or. Functional == 'DFT' )

print*, 'hrsreghgksrgkusregksrueyghfksehrgflserhglsriughsrileu'


   if ( Functional == 'DFT' .or. Functional == 'RHF' .or. Functional == 'RHA' ) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,ii,ia)
!$OMP DO
      do k=1,nbasis_pot
         B_kn(k)=0.d0
         do ii=1, n_occ
            do ia=n_occ+1,nbasis
               B_kn(k) = B_kn(k) + occnum(ii,3)*ovlap3nat(ii,ia,k)*FF(ia,ii)&
                                 / max(ennat(ia)-ennat(ii), zero)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
   else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,ii,ia,Delta_n,Delta_e)
!$OMP DO
      do k=1,nbasis_pot
         B_kn(k)=0.d0
         do ii=1, nbasis
            do ia=1, nbasis
            if (ia /= ii ) then 
               Delta_n=occnum(ii,3) - occnum(ia,3)
               if ( abs(Delta_n) > 2.d0*Dn_lim ) then
                  Delta_e=ennat(ia)-ennat(ii)
                  if(abs(Delta_e) < small_e) then
                     if ( Delta_e >= 0 ) then 
                        Delta_e=small_e
                     else
                        Delta_e=-small_e
                     endif
                  endif
                  B_kn(k) = B_kn(k) + ovlap3nat(ii,ia,k)*FF(ia,ii)/Delta_e
               endif
            endif
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
   endif

!.....Density on Grid
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,ia,dm_gr)
!$OMP DO
   do ig=1,ngrid
      chdens_gr(ig,1) = 0.d0
      do ia=1,nbasis
         dm_gr = vec_nat_grid(ig,ia)*vec_nat_grid(ig,ia)
         chdens_gr(ig,1) = chdens_gr(ig,1) + occnum(ia,1)*dm_gr
      enddo
      chdens_gr(ig,3) = 2.d0*chdens_gr(ig,1); chdens_gr(ig,2)=chdens_gr(ig,1)
   enddo
!$OMP END DO
!$OMP END PARALLEL

   grad_S=.false.

!..Calculate the CEDA part of A and B (CPU transistor exterminator!)
   if ((.not. grad_S) &!.and.(Functional == 'DFT' .or. Functional == 'RHF') &
        .and. (xld > zero .or. do_ceda)) then

      allocate (Pkn(lnbasis_pot,lnbasis_pot), Qkn(lnbasis_pot,lnbasis_pot))

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,igp,rr,ia,dm_gr,uxxtmp)
!$OMP DO
      do ig=1,ngrid
         uxxtmp=0.d0
         if ( Functional /= 'DFT' .or. hyb_mix > small) then
            do igp=1,ngrid
               dm_gr=0.d0
               do ia=1,n_occ
                  dm_gr= dm_gr + occnum(ia,3)*vec_nat_grid(ig,ia)*vec_nat_grid(igp,ia)
               enddo
               if(ig /= igp) then
                  rr=sqrt(&
                          (x_grid(ig)-x_grid(igp))**2+&
                          (y_grid(ig)-y_grid(igp))**2+&
                          (z_grid(ig)-z_grid(igp))**2)
                  uxxtmp=uxxtmp + w_grid(igp)*dm_gr**2 / rr
               endif
            enddo
         endif
         Uxx(ig)=uxxtmp
      enddo
!$OMP END DO
!$OMP END PARALLEL

   if (Functional == 'DFT' .or. Functional == 'RHF' .or. Functional == 'RHA' ) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,ii,ia,k,l)
!$OMP DO
      do k=1,nbasis_pot
         do l=1,k
            Pkn(k,l) = 0.d0
            Qkn(k,l) = 0.d0
            do ig=1,ngrid
               Pkn(k,l) = Pkn(k,l) + w_grid(ig)*bas_f_grid_pot(ig,k)*chdens_gr(ig,3)*bas_f_grid_pot(ig,l)
            enddo
            do ii=1,n_occ
               do ia=1,n_occ
                  Qkn(k,l) = Qkn(k,l) + occnum(ii,3)*ovlap3nat(ii,ia,k)*ovlap3nat(ia,ii,l)
               enddo
            enddo
            Pkn(l,k) = Pkn(k,l)
            Qkn(l,k) = Qkn(k,l)
         enddo   
      enddo
!$OMP END DO
!$OMP END PARALLEL
   else !RDMFT
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,ii,ia,k,l)
!$OMP DO
      do k=1,nbasis_pot
         do l=1,k
            Pkn(k,l) = 0.d0
            Qkn(k,l) = 0.d0
            do ig=1,ngrid
               Pkn(k,l) = Pkn(k,l) + w_grid(ig)*bas_f_grid_pot(ig,k)*chdens_gr(ig,3)*bas_f_grid_pot(ig,l)
            enddo
            do ii=1, nbasis
               do ia=1, nbasis
                  Qkn(k,l) = Qkn(k,l) + occnum(ii,3)*ovlap3nat(ii,ia,k)*ovlap3nat(ia,ii,l)
               enddo
            enddo
            Pkn(l,k) = Pkn(k,l)
            Qkn(l,k) = Qkn(k,l)
         enddo   
      enddo
!$OMP END DO
!$OMP END PARALLEL
   endif

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig, igp, rr, uhxtmp )
!$OMP DO
      do ig=1,ngrid
         uhxtmp =0.d0
         do igp=1,ngrid
            if (ig /= igp) then
               rr=sqrt(&
                       (x_grid(ig)-x_grid(igp))**2+&
                       (y_grid(ig)-y_grid(igp))**2+&
                       (z_grid(ig)-z_grid(igp))**2)
               uhxtmp = uhxtmp + w_grid(igp)*chdens_gr(igp,3) / rr
            endif
         enddo
         Uhx(ig) = uhxtmp
         if ( Functional == 'DFT') then
            Uhx(ig) = Uhx(ig) + Vlc(1,ig) + Vlxc(1,ig)
         endif
         Uhx(ig) = chdens_gr(ig,3)*Uhx(ig)
      enddo
!$OMP END DO
!$OMP END PARALLEL

     fc_mix=0.5d0
     if ( Functional == 'DFT' .and. hyb_mix > small) fc_mix=0.5d0*hyb_mix

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k, ig)
!$OMP DO
      do k=1,nbasis_pot
         R_kn(k) = 0.d0
         do ig=1,ngrid
            R_kn(k)=R_kn(k)+ w_grid(ig)*bas_f_grid_pot(ig,k)*(Uhx(ig)-fc_mix*Uxx(ig))
         enddo
      enddo 
!$OMP END DO
!$OMP END PARALLEL
 
   if (Functional == 'DFT' .or. Functional == 'RHF' .or. Functional == 'RHA') then
      do k=1,nbasis_pot
         S_kn(k)=0.d0
         do ii=1,n_occ
            do ia=1,n_occ
               S_kn(k) = S_kn(k) + occnum(ii,3)*ovlap3nat(ii,ia,k)*FF(ii,ia)
            enddo
         enddo
      enddo
   else !RDMFT
      do k=1,nbasis_pot
         S_kn(k)=0.d0
         do ii=1,nbasis
            do ia=1,nbasis
               S_kn(k) = S_kn(k) + ovlap3nat(ii,ia,k)*FF(ii,ia)
            enddo
         enddo
      enddo
   endif

      A_tild=Pkn-Qkn; B_tild=R_kn-S_kn
      deallocate (Pkn, Qkn)
   else !(... .and. (xld > zero ))
      A_tild=0.d0; B_tild=0.d0
   endif

   if(grad_S) then
      A_tild=Ovchixi 
      B_tild=0.d0
   endif

   if(.not.allocated(V_bs)) then 
      allocate ( V_bs(lnbasis_pot),  V_bs_o(lnbasis_pot) )
      V_bs_o=0.d0
   endif

   if ( automatic_limit ) then
      call do_automatic_limit(A_kn, svd_cut) 
   else !( .not. automatic_limit )
      if (do_ceda) then
         A_kn = A_tild
         B_kn = B_tild
      else
         if(xld > small) then
            A_kn = A_kn + xld*A_tild
            B_kn = B_kn + xld*B_tild
         endif
      endif
 
      call invert_Akn(A_kn, Akn_inv, svd_cut)

!     AinvB = MATMUL(B_kn, Akn_inv) 

      do k=1,nbasis_pot
         ss=0.d0
         do l=1,nbasis_pot
            ss=ss+Akn_inv(k,l)*B_kn(l)
         enddo
         AinvB(k)=ss
      enddo

      if( (.not.int_pot_basis) .or. (.not.scr_ch_const) ) then
         xmu_oep=0.d0 ! Constraint applied only to integral basis xi
         V_bs = AinvB
      else
!..Enable screening charge and positivity constraint
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,l,ss)
!$OMP DO
         do k=1,nbasis_pot
            ss=0.d0
            do l=1,nbasis_pot
               ss=ss+Akn_inv(k,l)*X_bas_pot(l)
            enddo
            AinvS(k)=ss
         enddo
!$OMP END DO
!$OMP END PARALLEL

         SAinvS = DOT_PRODUCT(X_bas_pot, AinvS)
         SAinvB = DOT_PRODUCT(X_bas_pot, AinvB)

         xmu_oep=(veffnorm-SAinvB)/SAinvS

         V_bs = AinvB + xmu_oep*AinvS

         pos_mix1=pos_mix
         pos_old=1.d20

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,ig)
!$OMP DO
         do ig=1,ngrid
            Schdens(ig)=0.d0
            do k=1,nbasis_pot
               Schdens(ig)=Schdens(ig)+V_bs(k)*charg_pot(ig,k)
            enddo
            if( Schdens(ig) > -zero) then
               neg_scdens(ig)=.false.
            else
               neg_scdens(ig)=.true.
            endif
         enddo
!$OMP END DO
!$OMP END PARALLEL

         volume_neg=0.d0
         do ig=1,ngrid
            if ( neg_scdens(ig) ) volume_neg=volume_neg+w_grid(ig)
         enddo
         fc_pen_v=pos_penalty !/volume_neg

!........Positivity criterion
         pos=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig)
!$OMP DO REDUCTION (+:pos)
         do ig=1,ngrid
            if ( neg_scdens(ig) )&
               pos=pos-w_grid(ig)*Schdens(ig)
         enddo
!$OMP END DO
!$OMP END PARALLEL

         write(6,'("Initial Positivity test:", f20.10, " electrons")')pos

         B_kn_in = B_kn
         V_bs_in=V_bs
!........Lool to enforce positivity of screeing charge
         if(posit_const) then
Pos_l:   do ipos=1,Nipos

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,ig,ss,fc_neg)
!$OMP DO
            do k=1,nbasis_pot
               ss=X_bas_pot(k)
               do ig = 1,ngrid
                  if(neg_scdens(ig))&
                     ss=ss+w_grid(ig)*charg_pot(ig,k)

!                    fc_neg=1._dp/(1._dp+exp(Schdens(ig)/0.000000001_dp))
!                    ss=ss+w_grid(ig)*fc_neg*charg_pot(ig,k)

!                    ss=ss-2.d0*w_grid(ig)*charg_pot(ig,k) ! Does not work

               enddo
               X_bar(k)=ss
            enddo
!$OMP END DO
!$OMP END PARALLEL
         
!...........Add penalty to B
            B_kn = B_kn_in + fc_pen_v * X_bar 

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,l,ss)
!$OMP DO
            do k=1,nbasis_pot
               ss=0.d0
               do l=1,nbasis_pot
                  ss=ss+Akn_inv(k,l)*B_kn(l)
               enddo
               AinvB(k)=ss
            enddo
!$OMP END DO
!$OMP END PARALLEL

            SAinvB = DOT_PRODUCT(X_bas_pot, AinvB) 

            xmu_oep=(veffnorm-SAinvB)/SAinvS

!...........New Potential
            DV_bs = AinvB + xmu_oep*AinvS

            V_bs = (1.d0-pos_mix1)*V_bs +pos_mix1*DV_bs

!..........positivity
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,k,ss)
!$OMP DO
            do ig=1,ngrid
               ss=0.d0
               do k=1,nbasis_pot
                  ss=ss+V_bs(k)*charg_pot(ig,k)
               enddo
               Schdens(ig)=ss
               if( Schdens(ig) > -zero) then
                  neg_scdens(ig)=.false.
               else
                  neg_scdens(ig)=.true.
               endif
            enddo
!$OMP END DO
!$OMP END PARALLEL

             pos=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ig,k)
!$OMP DO REDUCTION (+:pos)
             do ig=1,ngrid
                if ( neg_scdens(ig) )&
                   pos=pos-w_grid(ig)*Schdens(ig)
             enddo
!$OMP END DO
!$OMP END PARALLEL

            if(mod(ipos,50)==0 .or. ipos == 1) &
               write(6,'(i5,"  Positivity test:", f20.10, "  Pos mix:",e15.5)')ipos,pos,pos_mix1

!...........Convergence check
            if(pos < pos_small) then
               write(6,'(i5,"  Positivity test:", f20.10, "  Pos mix:",e15.5)')ipos,pos,pos_mix1
               exit Pos_l
            endif

            if( pos < pos_old ) then
                pos_mix1=pos_mix1*1.01d0
            else
                pos_mix1=pos_mix1*0.98d0
            endif
            pos_old=pos
            if(pos_mix1 < 1.d-12) pos_mix1=pos_mix
 
            if(ipos == Nipos) print*,"OEP:find positive scr charge: All iterations done", ipos
         enddo Pos_l
         print*,"Positivity test:",pos, 'After:', ipos, 'cycles'
         endif !posit_const
      endif !(.not.int_pot_basis) .or. (.not.scr_ch_const) )
   endif !( automatic_limit )

!..Mix for more scf iterations
   if (iscf >1 ) V_bs = xmix_OEP * V_bs + (1.d0-xmix_OEP) * V_bs_o
   V_bs_o = V_bs

!..Calculate Veff, i.e. V_bs matrix elements
   if(.not.allocated(Veff)) allocate ( Veff(lnbasis, lnbasis)  )
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(l,m,k,ss)
!$OMP DO
   do l=1,nbasis
      do m=1,l
         ss=0.d0
         do k=1, nbasis_pot
            ss = ss + ovlap3(m,l,k)*V_bs(k)
         enddo
         Veff(m,l)= ss
         Veff(l,m) = ss
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!..The effective Hamiltonian
!   H_eff = Hcore + Veff

! print*,' Veff:'
! write(6,'(15F20.10)')((Veff(k,l),k=1,15),l=1,15)



!   do kk=1,1 !20
   do k=1,nbasis
      do l=1,nbasis
         DMMM(k,l)=0._dp
         do ia=1,nbasis
!           do ia=1,n_occ
        DMMM(k,l)=DMMM(k,l) + occnum(ia,3)*vecnat(k,ia)*vecnat(l,ia)     
!         DMMM(k,l)=DMMM(k,l) + 2*vecnat(k,ia)*vecnat(l,ia)
         enddo
      enddo
   enddo
print*, 'lalalallalalallllllllllllllllllllllllllllllllllllllllllllllllllllllllllll'
!   call sum_intg_22(DMMM,Vcoul)
!   call sum_intg_11(DMMM,Vexch)

 H_eff = Hcore + Vcoul - 0.5*Vexch + Veff

!..Calculation of the objective function
!  do ia=1,nbasis
!     do ii=1,nbasis
!        Veff_aa(ia,ii)=0.d0
!        do k=1, nbasis_pot
!           Veff_aa(ia,ii) = Veff_aa(ia,ii) + ovlap3nat(ia,ii,k)*V_bs(k)
!        enddo
!     enddo
!  enddo

!  if (xld > 1.d-24 ) then
!     Eob = 0.d0
!     do ii=1, n_occ
!        do ia=n_occ+1,nbasis
!           Eob = Eob + (FF(ia,ii)-Veff_aa(ia,ii))**2/max(ennat(ia)-ennat(ii),zero)
!        enddo 
!     enddo
!     SEob=0.d0
!     Eob=0.d0
!     do k=1, n_occ
!        do l=1, n_occ
!           SEob = SEob + A_tild(k,l)*V_bs(k)*V_bs(l)
!           Eob = Eob + A_kn(k,l)*V_bs(k)*V_bs(l)
!        enddo
!        SEob = SEob-2.d0*B_tild(l)*V_bs(l)
!        Eob = Eob-2.d0*B_kn(l)*V_bs(l)
!     enddo
!     print*,'===============For T vs S (lambda curve):'
!     write(6,'("lambda, T0, S",3f20.10)')xld, Eob, SEob
!     print*,'========================================='
!  endif

   call diagon_H(info)

   E_HOMO=ennat(max(ibond(1),ibond(2)))

   deallocate ( ovlap3nat, A_kn, Akn_inv, FF, A_tild, Fii, FFp )
   deallocate ( Schdens, neg_scdens )
end subroutine analytic_Veff

