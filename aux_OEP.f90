!==========================================================================================
subroutine set_up_pot_basis()
   use global, ONLY:xi_analytic; implicit none
   if ( xi_analytic ) then
      call set_up_pot_basis_analytic()
   else
      call set_up_pot_basis_grid()
   endif
end subroutine set_up_pot_basis
!==========================================================================================

subroutine set_up_pot_basis_analytic()
!..Global
   use global; use grid_params; use matrices; use basis_set
   use functional_m, only:int_pot_basis,grad_S; implicit none

!..Local
   integer :: ibs, jbs, igr, iga, iat, n_exp_x, n_exp_y, n_exp_z
   real(dp) :: xx, yy, zz 
   real(dp) :: xd, yd, zd, expob, Xinteg
   real(dp), external :: f_bas_pot, X_integr

   if (.not. allocated(bas_f_grid_pot) ) allocate( bas_f_grid_pot(ngrid, lnbasis_pot) )
   if (.not. allocated(charg_pot) )      allocate( charg_pot(ngrid, lnbasis_pot) )
   if (.not. allocated(X_bas_pot) )      allocate( X_bas_pot(lnbasis_pot) )

   print*,'setup pot basis...'

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ibs, igr, xx, yy, zz)
!$OMP DO
   do ibs=1,nbasis_pot
      X_bas_pot(ibs)=0.0_dp
      do igr=1,ngrid
         xx=x_grid(igr); yy=y_grid(igr); zz=z_grid(igr)
         charg_pot(igr,ibs) = f_bas_pot( ibs, xx, yy, zz)
         if(.not.int_pot_basis) bas_f_grid_pot(igr,ibs) = charg_pot(igr,ibs)
         X_bas_pot(ibs) = X_bas_pot(ibs) + w_grid(igr)*charg_pot(igr,ibs)
      enddo
   enddo 
!$OMP END DO
!$OMP END PARALLEL

   if (int_pot_basis) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ibs, igr, iat, xd, yd, zd, n_exp_x, n_exp_y, n_exp_z, iga, expob, Xinteg )
!$OMP DO
   do igr=1,ngrid
      do ibs=1,nbasis_pot
         iat=iat_bas_pot(ibs)
         xd = xcoo(iat) - x_grid(igr)
         yd = ycoo(iat) - y_grid(igr)
         zd = zcoo(iat) - z_grid(igr)
         n_exp_x=Int(expo_x_pot(ibs)+0.001_dp)
         n_exp_y=Int(expo_y_pot(ibs)+0.001_dp)
         n_exp_z=Int(expo_z_pot(ibs)+0.001_dp)
         bas_f_grid_pot(igr,ibs)=0.0_dp
         do iga=1,nga_bas_pot(ibs)
            expob=expo_bas_pot(iga,ibs)
            Xinteg=X_integr(n_exp_x, n_exp_y, n_exp_z, expob, xd, yd, zd)
            bas_f_grid_pot(igr,ibs)=bas_f_grid_pot(igr,ibs)+coef_bas_pot(iga,ibs)*fact_norm_pot(iga,ibs)*Xinteg
         enddo
!        print*,igr,ibs,bas_f_grid_pot(igr,ibs)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL


   if(.not.allocated(ovlap_pot_chi)) then
      allocate(ovlap_pot_chi(lnbasis_pot, lnbasis_pot))
      ovlap_pot_chi=ovlap_pot
   endif

!.. Overlap of xi basis: 
!.. xi basis is not normalizable and overlap is either infinite or 0 
!.. depending on symmetry.
   
   if ( grad_S) then
      if(.not.allocated(Ovchixi)) allocate(Ovchixi(lnbasis_pot, lnbasis_pot))
      do ibs=1,nbasis_pot
         do jbs=1, ibs
            Ovchixi(ibs,jbs)=0.0_dp
            do igr=1,ngrid
               Ovchixi(ibs,jbs)=Ovchixi(ibs,jbs)+0.5_dp*w_grid(igr)*(bas_f_grid_pot(igr,ibs)*charg_pot(igr,jbs)+ &
                                                                    bas_f_grid_pot(igr,jbs)*charg_pot(igr,ibs))
            enddo
            Ovchixi(jbs,ibs)=Ovchixi(ibs,jbs)
         enddo
      enddo
   endif

   endif

   print*,'... Done'
end subroutine set_up_pot_basis_analytic

!==========================================================================================
subroutine set_up_pot_basis_grid()
!..Global
   use global; use grid_params; use matrices; use basis_set
   use functional_m, ONLY:int_pot_basis; implicit none

!..Local
   integer :: ibs, igr, jgr
   real(dp) :: xx, yy, zz, rr, wrr
   real(dp), external :: f_bas_pot
   
   if (.not. allocated(bas_f_grid_pot) ) allocate( bas_f_grid_pot(ngrid, lnbasis_pot) )
   if (.not. allocated(charg_pot) )      allocate( charg_pot(ngrid, lnbasis_pot) )
   if (.not. allocated(X_bas_pot) )      allocate( X_bas_pot(lnbasis_pot) )

   print*,'setup pot basis...'

!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(ibs, igr, xx, yy, zz)
!!$OMP DO
   do ibs=1,nbasis_pot
      X_bas_pot(ibs)=0.0_dp
      do igr=1,ngrid
         xx=x_grid(igr); yy=y_grid(igr); zz=z_grid(igr)
         charg_pot(igr,ibs) = f_bas_pot( ibs, xx, yy, zz)
         if(.not.int_pot_basis)   bas_f_grid_pot(igr,ibs) = charg_pot(igr,ibs)
         X_bas_pot(ibs) = X_bas_pot(ibs) + w_grid(igr)*charg_pot(igr,ibs)
      enddo
   enddo 
!!$OMP END DO
!!$OMP END PARALLEL

   if (int_pot_basis) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ibs, igr, jgr, rr)
!$OMP DO
   do igr=1,ngrid
      bas_f_grid_pot(igr,:)=0.0_dp
      do jgr=1,ngrid
         if(igr /= jgr) then
            rr=sqrt(&
                     (x_grid(igr)-x_grid(jgr))**2+&
                     (y_grid(igr)-y_grid(jgr))**2+&
                     (z_grid(igr)-z_grid(jgr))**2)
            wrr=w_grid(jgr)/rr
            do ibs=1,nbasis_pot
               bas_f_grid_pot(igr,ibs) = bas_f_grid_pot(igr,ibs) + wrr*charg_pot(jgr,ibs)
            enddo
         endif
      enddo
   enddo 
!$OMP END DO
!$OMP END PARALLEL
   endif

   if(.not.allocated(ovlap_pot_chi)) then
      allocate(ovlap_pot_chi(lnbasis_pot, lnbasis_pot))
      ovlap_pot_chi=ovlap_pot
   endif

   print*,'...done'

end subroutine set_up_pot_basis_grid
!-----------------------------------------------------------------------------------------------------------------

subroutine invert_Akn(A_kn, Akn_inv, zz)
!..Global
   use global; use matrices; use orbocc; implicit none
!..Arguments
   real(dp) :: A_kn(lnbasis_pot, lnbasis_pot), Akn_inv(lnbasis_pot, lnbasis_pot)
   real(dp) :: zz

!..Local
   real(dp), allocatable :: A(:,:)
   real(dp), allocatable :: B(:,:)
   real(dp), allocatable :: vectr(:,:), eigs(:)
   real(dp), allocatable :: Zmat(:)
   integer :: lndim, ndim, mu, nu, ku, info, inull, inotn, k
!..Definitions for LAPACK routine
   character(1) :: JOBZ, RANGE, UPLO
   INTEGER :: IL, ITYPE,  IU,  LDA,  LDB,  LDZ
   INTEGER LWORK, M, N
   real(dp) :: VL,VU
   real(dp) :: ABSTOL
   real(dp), allocatable :: WORK(:)
   integer, allocatable :: IWORK(:)
   integer, allocatable :: ifail(:)
   logical :: l_S_diag
   
   real(dp) :: DLAMCH
   external DLAMCH

!..Local array allocation
   ndim=nbasis_pot
   lndim=ndim
   allocate( A(lndim,lndim), B(lndim,lndim), Zmat(ndim), &
               vectr(lndim,lndim), eigs(lndim), &
               WORK(10*lndim), IWORK(5*lndim), ifail(lndim) )

   ITYPE=1; JOBZ='V'; RANGE='I'; UPLO='L'; N = ndim; LDA = lndim
   LDB = lndim; VL = 0.0_dp; VU = 0.0_dp; IL = 1; IU = ndim
   LDZ = lndim; LWORK=10*ndim

!  ABSTOL = 2*DLAMCH('S')
   ABSTOL = 0.0_dp
!  ABSTOL = 1.e-30_dp
      
   do mu=1, ndim   
      do nu=1, mu
         A(mu,nu) = A_kn(mu,nu)
         A(nu,mu) = A(mu,nu) 
      enddo
   enddo

   l_S_diag=.false.   
   if(l_S_diag) then
      do mu=1, ndim   
         do nu=1, mu
            B(mu,nu)=ovlap_pot(mu,nu)
            B(nu,mu)=B(mu,nu)
         enddo
      enddo
   else
      B=0.0_dp
      do mu=1, ndim   
         B(mu,mu)=1.0_dp
      enddo
   endif

   WORK=0.0_dp
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
      stop 'invert_Akn: Diagonalization failed!'
   endif

   print*,'Eigenvalues of A:'
   write(6,'(5e12.5)') (eigs(k),k=1,nbasis_pot)
   print*,'-----------------------------------'

   n_null_eiA=0; n_notn_eiA=0
   do ku=1,ndim
      if (abs(eigs(ku)) < zz ) then 
         n_null_eiA = n_null_eiA+1
      else
         n_notn_eiA = n_notn_eiA+1
      endif
   enddo
   ln_null_eiA = max(n_null_eiA,1); ln_notn_eiA = max(n_notn_eiA,1)
   zero_null = (n_null_eiA == 0 ); zero_notn = (n_notn_eiA == 0 )

   if (allocated(vec_A_null)) deallocate(vec_A_null)
      allocate(vec_A_null(lnbasis_pot,ln_null_eiA))
   if (allocated(vec_A_notn)) deallocate(vec_A_notn)
      allocate(vec_A_notn(lnbasis_pot,ln_notn_eiA))
   if (allocated(enn_A_notn)) deallocate(enn_A_notn)
      allocate(enn_A_notn(ln_notn_eiA))
   inull=1; inotn=1
   do ku=1,ndim
      if (abs(eigs(ku)) < zz ) then 
         vec_A_null(:,inull)=vectr(:,ku)
         eigs(ku)= 0.d0
!        eigs(ku)=1.0_dp/zz
         inull=inull+1
      else
         vec_A_notn(:,inotn)=vectr(:,ku)
         enn_A_notn(inotn)=eigs(ku)
         eigs(ku)=1.0_dp/eigs(ku)
         inotn=inotn+1
      endif 
   enddo
   
!..Reconstruct Akn_inv from eigenvalues and eigenvectors
   do mu=1,ndim
      do nu=1,mu
         Akn_inv(mu,nu)=0.0_dp
         do ku=1,ndim
            Akn_inv(mu,nu) = Akn_inv(mu,nu) + vectr(mu,ku)*eigs(ku)*vectr(nu,ku)
         enddo
         Akn_inv(nu,mu)= Akn_inv(mu,nu)
      enddo
   enddo

end subroutine  invert_Akn
!=====================================================================================

subroutine invert_Akn_tilde(A_kn, Akn_inv, ndim)
!..Global
   use global; use matrices; use orbocc; implicit none
!..Arguments
   integer :: ndim
   real(dp) :: A_kn(ndim, ndim), Akn_inv(ndim, ndim)

!..Local
   real(dp), allocatable :: A(:,:)
   real(dp), allocatable :: B(:,:)
   real(dp), allocatable :: vectr(:,:), eigs(:)
   real(dp), allocatable :: Zmat(:)
   integer :: mu, nu, ku, info
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
   allocate( A(ndim,ndim), B(ndim,ndim), Zmat(ndim), &
               vectr(ndim,ndim), eigs(ndim), &
               WORK(10*ndim), IWORK(5*ndim), ifail(ndim) )

   ITYPE=1; JOBZ='V'; RANGE='I'; UPLO='L'; N = ndim; LDA = ndim
   LDB = ndim; VL = 0.0_dp; VU = 0.0_dp; IL = 1; IU = ndim
   LDZ = ndim; LWORK=10*ndim

   ABSTOL = 2*DLAMCH('S')
   ABSTOL = 0.0_dp
!  ABSTOL = 1.e-30_dp
      
   do mu=1, ndim   
      do nu=1, mu
         A(mu,nu) = A_kn(mu,nu)
         B(mu,nu) = ovlap_pot(mu,nu)
         A(nu,mu) = A(mu,nu) 
         B(nu,mu) = B(mu,nu)
      enddo
   enddo

   WORK=0.0_dp
   IWORK=0
   call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
                LDB,  VL,  VU, IL, IU, ABSTOL, M, eigs, vectr, &
                LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

   if (info /= 0) then
      print*,'Info:',info
      if(info < 0) print*,info, '-th argument of DSYGVX has wrong value'
      if(info > 0 .and. info <= ndim ) print*,info, '-th eigenvalue not converged'
      if(info > ndim .and. info <= 2*ndim ) &
         print*,info-ndim, '-th leading minor of matrix not positive definite'
      stop 'invert_Akn_tilde: Diagonalization failed!'
   endif

   do ku=1,ndim
      if (abs(eigs(ku)) > zero ) then 
         eigs(ku)=1.0_dp/eigs(ku)
      else
         eigs(ku)= 0.0_dp
      endif
   enddo

   do mu=1,ndim
      do nu=1,mu
         Akn_inv(mu,nu)=0.0_dp
         do ku=1,ndim
            Akn_inv(mu,nu) = Akn_inv(mu,nu) + vectr(mu,ku)*eigs(ku)*vectr(nu,ku)
         enddo
         Akn_inv(nu,mu)= Akn_inv(mu,nu)
      enddo
   enddo

end subroutine  invert_Akn_tilde

!=======================================================================================
subroutine do_automatic_limit(A_kn, B_kn, A_tild, B_tild, svd_cut, do_ceda)
!..Global
   use global; use matrices; use orbocc
   implicit none

!..Arguments
   real(dp) :: A_kn(lnbasis_pot, lnbasis_pot), B_kn(lnbasis_pot)
   real(dp) :: A_tild(lnbasis_pot, lnbasis_pot), B_tild(lnbasis_pot), svd_cut
   logical :: do_ceda
    
!..Local
   real(dp) :: Akn_inv(lnbasis_pot, lnbasis_pot), sss, V_bs0(lnbasis_pot)
   real(dp), allocatable :: Atpp(:,:), Atp(:,:), btp(:), Atpp_inv(:,:), bi(:), v0i(:), vbp(:)
   integer :: k, l, m, n

   call invert_Akn(A_kn, Akn_inv, svd_cut)
   allocate (Atpp(ln_null_eiA, ln_null_eiA), Atp(ln_null_eiA, ln_notn_eiA),&
             bi(ln_notn_eiA), btp(ln_null_eiA), Atpp_inv(ln_null_eiA, ln_null_eiA),&
             vbp(ln_null_eiA),v0i(ln_notn_eiA) )

      print*,'Dimension of Zero eigenvalue space:',n_null_eiA
      print*,'Dimension of non Zero eigenvalue space:',n_notn_eiA

   v0i = 0.0_dp
   do k=1,n_notn_eiA
      bi(k) = 0.0_dp
      do l=1,nbasis_pot
         bi(k)=bi(k)+B_kn(l)*vec_A_notn(l,k)
      enddo
      v0i(k) = bi(k)/max(enn_A_notn(k), zero)
   enddo

   Atpp=0.0_dp; Atp=0.0_dp; btp=0.0_dp; vbp=0.0_dp
   do k=1,n_null_eiA
      do l=1,n_null_eiA
         do m=1, nbasis_pot
            do n=1, nbasis_pot
               Atpp(k,l) = Atpp(k,l) + vec_A_null(m,k)*A_tild(m,n)*vec_A_null(n,l)
            enddo
            enddo
         enddo
      enddo

         do k=1,n_null_eiA
            do l=1,n_notn_eiA
               do m=1, nbasis_pot
                  do n=1, nbasis_pot
                     Atp(k,l) = Atp(k,l) + vec_A_null(m,k)*A_tild(m,n)*vec_A_notn(n,l)
                  enddo
               enddo
            enddo
         enddo

      do k=1,n_null_eiA
         do l=1,nbasis_pot
            btp(k)= btp(k)+ B_tild(l)*vec_A_null(l,k)
         enddo
      enddo

      call invert_Akn_tilde(Atpp, Atpp_inv, n_null_eiA)  

      do k=1,n_null_eiA
         do l=1, n_null_eiA
            sss=0.0_dp
            if(.not. do_ceda) then
            do m=1,n_notn_eiA         
               sss=sss+ Atp(l,m)*v0i(m)
            enddo
            endif
            vbp(k)=vbp(k)+Atpp_inv(k,l)*(-sss+btp(l))
         enddo
      enddo

   do l=1,nbasis_pot
      V_bs0(l)=0.0_dp
      do k=1, n_notn_eiA
         V_bs0(l)=V_bs0(l)+ v0i(k)*vec_A_notn(l,k)
      enddo
      V_bs(l)=V_bs0(l)
      do k=1, n_null_eiA
         V_bs(l)=V_bs(l)+ vbp(k)*vec_A_null(l,k)
      enddo
   enddo

   deallocate (Atpp, Atp, btp, Atpp_inv, bi, v0i, vbp)

end subroutine do_automatic_limit

!-----------------------------------------------------------------

subroutine grad_E_vec_eff_R(x0,y0,z0,x1,y1,z1,pfile)
!..Calculates the derivatives delta E / delta f(r) and
!  delta E / delta v(r), f: effective oorbital, v: effective potential
!  in real space line segment: (x0,y0,z0) ---> (x1,y1,z1) and outputs in pfile

   use global; use orbocc; use matrices; use energies
   implicit none

!..Arguments
   real(dp) :: x0,y0,z0,x1,y1,z1
   character(10) :: pfile

   integer :: m,ia,k,l,ii, n_occ, Nstep, ist
   real(dp) :: ss, ss1, dE_df, dE_dv, ff, xstep, ystep, zstep, xint, xx, yy, zz,dd
   real, allocatable :: St(:,:), PP(:,:), QQ(:,:), fb(:)
   real(dp),external :: f_bas, f_bas_pot

   Nstep=200

   allocate (St(nbasis,nbasis), PP(nbasis,nbasis), QQ(nbasis,nbasis), fb(nbasis) ) 

   n_occ=max(ibond(1),ibond(2))

   call make_S_klm()

   do k=1,nbasis
      do l=1,k
         ss=0._dp
         do m=1,nbasis_pot
            ss=ss+S3(k,l,m)*vec_eff(m)
         enddo
         St(k,l)=ss
         St(l,k)=ss
      enddo
   enddo

   do ii=1,n_occ
      do ia=n_occ+1,nbasis
         ss=0._dp; ss1=0._dp
         do m=1,nbasis
            do l=1,nbasis
               ss=ss+St(m,l)*vecnat(m,ia)*vecnat(l,ii)
               ss1=ss1+(F(m,l,1)-Hcore(m,l))*vecnat(m,ia)*vecnat(l,ii)
            enddo
         enddo
         PP(ii,ia)=-ss-ss1/max(ennat(ia)-ennat(ii), small) 
      enddo
   enddo

   do k=1,nbasis
      do l=1,nbasis
         ss=0._dp
         do ii=1,n_occ
            do ia=n_occ+1,nbasis
               ss=ss+vecnat(k,ii)*vecnat(l,ia)*PP(ii,ia)
            enddo
         enddo
         QQ(k,l) = ss       
      enddo
   enddo
 
   open(unit=937,file=pfile,status='unknown')
   xstep=(x1-x0)/float(Nstep-1); ystep=(y1-y0)/float(Nstep-1); zstep=(z1-z0)/float(Nstep-1)
   xx=x0; yy=y0; zz=z0
   do ist=1,Nstep
      do k=1,nbasis
         fb(k)=f_bas(k,xx,yy,zz)
      enddo
      dE_df=0._dp; dE_dv=0._dp; ff=0._dp
      do k=1,nbasis
         do l=1,nbasis
            call Xintegral2(k,l,xx,yy,zz,xint)
            ss=fb(k)*fb(l)
            dE_df=dE_df + QQ(k,l)*xint
            dE_dv=dE_dv + QQ(k,l)*ss
         enddo
         ff=ff+vec_eff(k)*f_bas_pot(k,xx,yy,zz)
      enddo
      dd=dE_df
      dE_df=dE_df*ff
      xx=xx+xstep; yy=yy+ystep; zz=zz+zstep
      write(937,'(3f10.6,3e20.10)')xx,yy,zz,dE_df,dd,dE_dv
   enddo

   close(937)
   return
end subroutine grad_E_vec_eff_R
