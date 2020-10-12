!Subroutine for performing the energy minimisation of the EFO 
!Objective functional by Nikitas and Toms method
subroutine eff_orbital_NIK_TOM(iscf,x0,y0,z0,x1,y1,z1,pfile)
!..Global
   use global; use orbocc; use matrices; use functional_m, only: xmix_OEP, ahyb, VHF, l_hybc, Functional
   use energies; use DFT_par, only: Vc_mn_DFT ; implicit none

!..Arguments
   real(dp) :: x0,y0,z0,x1,y1,z1
   character(10) :: pfile
   integer :: iscf

!..Local
   integer :: i, a, m, l, info, nvar, meth_case, ist
   real(dp) :: ss, xnorm, objf, objf_old=1.D5
   real(dp) :: obj1, obj2, obj3, obj4, mid, eps
   logical :: l_new_obj
   real(dp), allocatable :: X(:), K(:,:), B(:), C(:), U(:,:), XF(:), BU(:), DF(:) 
!..Vectors for Conjugate gradient 
   real(dp), allocatable :: GG(:), HH(:) 
   integer :: icall=0, Ndiis=6, Nd
   real(dp),allocatable :: f_out(:), f_old(:,:), df_old(:,:), ov(:,:), XF_o(:)
   real(dp) :: emax

   external En_FUN_vec, make_mats_TOM

   nvar=nbasis

   Allocate( X(nvar), K(nvar,nvar), B(nvar), C(nvar), U(nvar,nvar), XF(nvar), BU(nvar), DF(nvar) ) 
   Allocate( GG(nvar), HH(nvar) )

   allocate ( f_out(nvar), f_old(nvar,0:Ndiis), df_old(nvar,Ndiis), XF_o(Nvar), ov(nvar, nvar) )
!..Initialize effective orbital (first call)   
   if (.not.l_init_efo) then   
      print*,'Effective orbital initialization'
      if( .not. allocated(vec_eff) ) allocate (vec_eff(lnbasis))
      call expand_sqdensity(vec_eff)
      call norm_vec_eff(vec_eff, xnorm)
      l_init_efo=.true.
   endif

!..Initial value
   X=vec_eff

   ov=0._dp
   do l=1,nvar
      ov(l,l)=1._dp
   enddo

!..Loop goes here 
   do i=1,2000

      !..Initialise K, B

      if ( Functional == 'DFT' .or. Functional == 'RHF' ) then
         call make_mats_TOM(K,B)
      else
         call make_mats_LRDMFT(K,B)
      endif

      call diagon_lapack(nvar,K,ovlap,C,U)

      !..Loop to transform X to new basis
      XF=0._dp
      do a=1,nvar
         do m=1,nvar
            do l=1,nvar
               XF(a) = XF(a) + X(m)*ovlap(m,l)*U(l,a)
            enddo
         enddo
      enddo

      !..Transform B to new basis 
      BU=0._dp
      do a=1,nvar
         do m=1,nvar
            BU(a) = BU(a) + B(m) * U(m,a)
         enddo
      enddo


      !..Initial objective function calculation
      objf=0._dp   
      do a=1,nvar
         objf = objf + XF(a) * BU(a) - 0.5_dp * C(a) * XF(a)**2
      enddo

      meth_case=1 ! direct solving for f
!     meth_case=2 ! iterative solving with steepest decent (Nikitas)
!     meth_case=3 ! iterative solving with DIIS

      select case (meth_case)
      case(1)
         !..Interval bisection initialisation 
         m=0
         obj1=abs(C(nvar))
         obj2=-1._dp
         obj3=1._dp 
         obj4=0._dp

         ss=0._dp
         do a=1,nvar
            ss = ss + (BU(a)/min( C(a) + obj2,small ) )**2
         enddo
         obj4=ss- veffnorm

         if ( (obj3.lt.0._dp .and. obj4.lt.0._dp ) .or. ( obj3.gt.0._dp .and. obj4.gt.0._dp )) then 
            print*, 'Lambda incorrectly bracketed' 
            stop 
         endif

         if ( obj3 .lt. 0._dp ) then 
            mid=obj1
            ss=obj3
            obj1=obj2
            obj3=obj4
            obj2=mid
            obj4=ss
         endif
         !..Interval bisection loop 
         do 
            m=m+1

            mid = 0.5 * ( obj1 + obj2 )

            ss=0._dp
            do a=1,nvar
               ss = ss + (BU(a)/( C(a) + mid ))**2
            enddo
            ss=ss - veffnorm
     
            if ( ss .gt. 0._dp ) then 
               obj1=mid
               obj3 = ss
            else
               obj2=mid
               obj4=ss
            endif

            if(abs(obj1-obj2).lt.1.D-16) then 
               ss=0.5*(obj1 + obj2) 
               exit 
            endif

         enddo
         !..End of Interval Bisection 
         obj4=ss

         do a=1,nvar
            DF(a) = BU(a) / (ss + C(a))
         enddo

         !..Mixing
!         DF = 0.05_dp *  DF + 0.95_dp * XF
!         call norm_vec_eff(DF, ss)
         XF=DF

!         ss=DOT_PRODUCT(DF,DF)
         obj4=ss
!XF=veffnorm/ss * DF
         objf=0._dp
         do a=1,nvar
            objf = objf + XF(a) * BU(a) - 0.5_dp * C(a) * XF(a)**2
         enddo

         !..End of method
         !...............................

         !..Loop calculating the minimum with B&K held constant
         !..Via Conjugate gradient  method 
      case(2)
         !..initialise the CG method with the gradient DF
         ss=0._dp
         do a=1,nvar
            ss = ss + XF(a) * BU(a) - C(a) * XF(a)**2
         enddo

         !..Calculates df                                                                   
         DF = -( BU - C * XF ) + XF * ss/veffnorm
         !..Vectors for conjugate gradient 
         GG=DF 
         HH=GG
         m=0
         !..Conjugate gradient loop 
         do
            m=m+1
            !..Algerbraic minimum finding in direciton HH 
            obj2=0._dp
            obj3=0._dp 
            ss=0._dp
            do a=1,nvar
               obj2 = obj2 + HH(a) * BU(a) - C(a) * HH(a) * XF(a)
               obj3 = obj3 + C(a) * HH(a)**2
            enddo

            ss =  obj2/obj3
            DF = XF + ss * HH      

            !..Calculate the objective functional 
            XF=DF
            obj1=objf
            objf=0._dp 
            do a=1,nvar
               objf = objf + XF(a) * BU(a) - 0.5_dp * C(a) * XF(a)**2
            enddo
            !..Escape if convergence is achieved 

            if (  abs(objf-obj1).lt.1.D-16 ) exit 

            !..Calculation of the gradient at the new minimum 
            !..determine lambda
            ss=0._dp
            do a=1,nvar
               ss = ss + XF(a) * BU(a) - C(a) * XF(a)**2
            enddo
            !..Calculates gradient 
            DF = -( BU - C * XF ) + XF * ss/veffnorm

            !..Conjugate gradient calculation 
            obj2=0._dp 
            obj3=0._dp 
            do a=1,nvar
               obj2=obj2 + DF(a)**2!( DF(a) - GG(a) ) * DF(a) 
               obj3=obj3 + GG(a)**2
            enddo

            ss = obj2/obj3
            GG=DF
            HH=GG + ss * HH

         enddo
         !..End of CG loop 
         !..Transform f back into original basis 
      case(3)
         do ist=1,123456789
            ss=0._dp
            do a=1,nvar
               ss = ss + XF(a) * BU(a) - C(a) * XF(a)**2
            enddo
            DF = -( BU - C * XF ) + XF * ss/veffnorm
            eps=0._dp; ss=0._dp
            do a=1,nvar
               eps = eps + DF(a)*BU(a) - C(a) * XF(a) *DF(a)
               ss = ss + C(a) * DF(a)*DF(a)
            enddo
            eps=eps/ss
            ss=0._dp
            do a=1,nvar
               ss = ss+DF(a)*DF(a)
            enddo
            ss=sqrt(1._dp+eps*eps*ss/veffnorm)
            XF=(XF+eps*DF)/ss

!           print*,'a'
!           print*,XF
            if (ist == 100 ) f_old(:,0)=XF_o
            if ( ist > 100) then
               call DIIS_f(icall, Ndiis, nvar, XF, f_out, f_old, df_old, ov, emax)
               XF=f_out
            endif
!           print*,'b'
!           print*,XF

            objf=0._dp
            do a=1,nvar
               objf = objf + XF(a) * BU(a) - 0.5_dp * C(a) * XF(a)**2
            enddo
   xnorm=0._dp
   emax=0._dp
   do a=1,nvar
      emax=emax+(XF(a)-XF_o(a))**2
      xnorm=xnorm+XF(a)**2
   enddo
   print*, ist, 'emax',sqrt(emax), sqrt(xnorm), objf
 
!           print*,ist, objf
            XF_o = XF
         enddo
   
      case default
         stop 'Effective orbital; unknown case'
      end select

      !stop 
      !..End of steepest descent
      
      DF=X
      X=0._dp 
      do l=1,nvar
         do a=1,nvar
            X(l) = X(l) + XF(a) * U(l,a)
         enddo
      enddo
      !..mixing 
      X=0.5*X+0.5*DF
      !..Now normalize optimal X and vec_eff 
      call norm_vec_eff(X, ss)
      vec_eff = X

      !..Print info for each step 
      print*, i, m, objf , objf - objf_old, C(nvar)
!      print*, i, m, C(nvar) , obj4,  objf - objf_old

      if ( abs(objf - objf_old) < 1.D-14 ) then 
         print*, 'convergence achieved, Objective function: ', objf 
         exit
      else if ( i == 5000 ) then 
         print*, 'convergence not achieved after 5000 cycles!' 
         stop
      else
         objf_old=objf
      endif

   enddo


   print*,' Objective function: ',objf
   print*,'---------------------------------------------------------'
   !..prints eigenvalues of K 
   if (.true.) then 
      print*, '----------K EIGENVALUES-----------'
      print*, C
      print*, '----------------------------------'
   endif

   if( .not. allocated(vec_eff0) ) allocate(vec_eff0(nbasis))
   if (iscf >1 ) vec_eff = xmix_OEP * vec_eff + (1.d0-xmix_OEP) * vec_eff0
!   if (iscf >1 ) vec_eff = 0.5_dp * (vec_eff +  vec_eff0)
!   if (iscf >1 ) vec_eff =0.6_dp * vec_eff0 +  0.4_dp * vec_eff
   call norm_vec_eff(vec_eff, ss)
   vec_eff0 = vec_eff

   call build_Veff()
   if ( Functional=='DFT' ) then
      if ( .not. l_hybc ) then 
         H_eff = Hcore + Veff + Vc_mn_DFT
      else 
         !      H_eff = Hcore + ahyb*Veff + (1._dp - ahyb)*VHF !..original hybrid
         H_eff = Hcore + ahyb*Veff + (1._dp - ahyb)*VHF  +  Vc_mn_DFT !..Hybrid + correlation
      endif

   else
      H_eff = Hcore + Veff
   endif

   call diagon_H(info)    

   call construct_f(0)

!   call plot_potential_EFO(x0,y0,z0,x1,y1,z1,'potentialf_EFO')
!  call grad_E_vec_eff_R(x0,y0,z0,x1,y1,z1,'func_ders')
   return
   Deallocate( X, K, B, XF, U, C, BU )
end subroutine eff_orbital_NIK_TOM 

!..Calculates the Matrix A, and the vector B 
subroutine make_mats_TOM(A,B)
   use global; use orbocc; use matrices
   use DFT_par, only: Vc_mn_DFT
   use functional_m, only: Functional
   use energies; implicit none

!..Arguments
   real(dp) :: A(lnbasis, lnbasis), B(lnbasis) 

   integer :: m,ia,k,l,ii, n_occ
   real(dp) :: ss
   real(dp), allocatable :: QQ(:,:,:), ST(:,:,:), FF(:,:), VV(:,:) 

   allocate ( QQ(nbasis,nbasis,nbasis), ST(nbasis,nbasis,nbasis), &
              FF(nbasis,nbasis), VV(nbasis,nbasis) )

   n_occ=max(ibond(1),ibond(2))

   call make_S_klm()

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,k,ss,l) 
!$OMP DO
   do m=1,nbasis
      do ia=n_occ+1,nbasis
         do k=1,nbasis
            ss=0._dp
            do l=1,nbasis
               ss = ss + S3(k,l,m)*vecnat(l,ia)
            enddo
            QQ(k,ia,m) = ss
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,ii,k,ss) 
!$OMP DO
   do m=1,nbasis
      do ia=n_occ+1,nbasis
         do ii=1,n_occ
            ss = 0._dp
            do k=1,nbasis
               ss = ss + QQ(k,ia,m)*vecnat(k,ii)
            enddo
            ST(ii,ia,m) = ss
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ii,k,l) 
!$OMP DO
   do k=1,nbasis
      do l=1,nbasis
         A(k,l) = 0.d0
         do ia=n_occ+1,nbasis
            do ii=1, n_occ
               A(k,l) = A(k,l) + (ST(ii,ia,k) * ST(ii,ia,l))/max(ennat(ia)-ennat(ii), small)
            enddo
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL
   A=-2._dp * A
!..V_HX
if ( Functional=='RHF' ) then 
   VV(:,:)=F(:,:,1)-Hcore(:,:)
else
   VV(:,:)=F(:,:,1)-Hcore(:,:)-Vc_mn_DFT(:,:)
endif 
!..V_HXC
!   VV(:,:)=F(:,:,1)-Hcore(:,:)

   FF=0._dp

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ii,k,l) 
!$OMP DO
   do ia=n_occ+1,nbasis
      do ii=1, n_occ
         FF(ii,ia) = 0.d0
         do k=1,nbasis
            do l=1,nbasis
                  FF(ii,ia) = FF(ii,ia) + vecnat(k,ia)*VV(k,l)*vecnat(l,ii)
            enddo
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL
!  write(6,'(15f10.6)')((FF(ia,ii),ia=1,nbasis),ii=1,nbasis)

   do m=1,nbasis
      ss=0._dp
      do ii=1,n_occ
         do ia=n_occ+1,nbasis
            ss = ss + FF(ii,ia)*ST(ii,ia,m)/max(ennat(ia)-ennat(ii), small)
         enddo
      enddo
      B(m)=-2._dp*ss
   enddo

   deallocate ( QQ, ST, FF )
   return
end subroutine make_mats_TOM

!========================================
!........>>END OF TOMS CHANGES<<.........
!========================================


!..Calculates the Matrix A, and the vector B 
subroutine make_mats_LRDMFT(A,B)
   use global; use orbocc; use matrices
   use DFT_par, only: Vc_mn_DFT
   use functional_m, only: Functional, Dn_lim, small_e
   use energies; implicit none

!..Arguments
   real(dp) :: A(lnbasis, lnbasis), B(lnbasis) 

   integer :: m,ia,k,l,ii, n_occ
   real(dp) :: ss, Delta_e, Delta_n, fac_dif
   real(dp), allocatable :: QQ(:,:,:), ST(:,:,:), FF(:,:), VV(:,:), FFp(:,:) 

   allocate ( QQ(nbasis,nbasis,nbasis), ST(nbasis,nbasis,nbasis), &
              FF(nbasis,nbasis), VV(nbasis,nbasis), FFp(nbasis,nbasis))

   n_occ=max(ibond(1),ibond(2))

   call make_S_klm()

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,k,ss,l) 
!$OMP DO
   do m=1,nbasis
      do ia=1,nbasis
         do k=1,nbasis
            ss=0._dp
            do l=1,nbasis
               ss = ss + S3(k,l,m)*vecnat(l,ia)
            enddo
            QQ(k,ia,m) = ss
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,ii,k,ss) 
!$OMP DO
   do m=1,nbasis
      do ia=1,nbasis
         do ii=1,nbasis; if(ia /= ii ) then
            ss = 0._dp
            do k=1,nbasis
               ss = ss + QQ(k,ia,m)*vecnat(k,ii)
            enddo
            ST(ii,ia,m) = ss
         endif; enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ii,k,l) 
!$OMP DO
   do k=1,nbasis
      do l=1,nbasis
         A(k,l) = 0.d0
         do ia=1,nbasis
            do ii=1,nbasis; if(ia /= ii ) then
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
                  A(k,l) = A(k,l) - (ST(ii,ia,k)*ST(ii,ia,l))*Delta_n/Delta_e
               endif
            endif; enddo
         enddo
         A(k,l) = 0.5_dp * A(k,l)
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ii,fac_dif,k,l,VV,ia) 
!$OMP DO
   do ii=1, nbasis
      fac_dif=0.5d0*fac_h(ii)
      do k=1,nbasis
         do l=1,k
            VV(k,l)=F(k,l,ii)-fac_dif*hcore(k,l)
            VV(l,k)=VV(k,l)
         enddo
      enddo

      do ia=1, nbasis
         FFp(ii,ia) = 0.d0
         do k=1,nbasis
            do l=1,nbasis
                  FFp(ii,ia) = FFp(ii,ia) &
                             + vecnat(k,ia)*VV(k,l)*vecnat(l,ii)
            enddo
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

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

   do m=1,nbasis
      ss=0._dp
      do ii=1,nbasis
         do ia=1,nbasis; if ( ia /= ii ) then
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
               ss = ss + FF(ii,ia)*ST(ii,ia,m)/Delta_e
            endif
         endif; enddo
      enddo
      B(m)=-ss
   enddo

   deallocate ( QQ, ST, FF, VV, FFp )
   return
end subroutine make_mats_LRDMFT

! For RDMFT matrices A and B shoule be defined differently:
subroutine make_mats_RDMFT(A,B)
   use global; use orbocc; use matrices
   use DFT_par, only: Vc_mn_DFT
   use functional_m, only: Functional, Dn_lim, small_e
   use energies; implicit none

!..Arguments
   real(dp) :: A(lnbasis, lnbasis), B(lnbasis) 

   integer :: m,ia,k,l,ii, n_occ
   real(dp) :: ss, fac_dif, Delta_n, Delta_e
   real(dp), allocatable :: QQ(:,:,:), ST(:,:,:), FF(:,:), VV(:,:), FFp(:,:), Fii(:,:) 

   allocate ( QQ(nbasis,nbasis,nbasis), ST(nbasis,nbasis,nbasis), &
              FF(nbasis,nbasis), VV(nbasis,nbasis), FFp(nbasis,nbasis),Fii(nbasis,nbasis) )

   n_occ=max(ibond(1),ibond(2))

   call make_S_klm()

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,k,ss,l) 
!$OMP DO
   do m=1,nbasis
      do ia=1,nbasis
         do k=1,nbasis
            ss=0._dp
            do l=1,nbasis
               ss = ss + S3(k,m,l)*vecnat(l,ia)
            enddo
            QQ(k,m,ia) = ss
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,ii,k,ss) 
!$OMP DO
   do m=1,nbasis
      do ia=1,nbasis
         do ii=1,nbasis
            ss = 0._dp
            do k=1,nbasis
               ss = ss + QQ(m,k,ia)*vecnat(k,ii)
            enddo
            ST(ii,ia,m) = ss
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ii,k,l) 
!$OMP DO
   do k=1,nbasis
      do l=1,nbasis
         A(k,l) = 0.d0
         do ia=1, nbasis
            do ii=1, nbasis; if ( ia /= ii ) then
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
                  A(k,l) = A(k,l) - (ST(ii,ia,l)*ST(ii,ia,k))*Delta_n/Delta_e
               endif
            endif; enddo
         enddo
         A(k,l) = 0.5_dp * A(k,l)
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ii,fac_dif,k,l,Fii,ia) 
!$OMP DO
   do ii=1, nbasis
      fac_dif=0.5d0*fac_h(ii)
      do k=1,nbasis
         do l=1,k
            Fii(k,l)=F(k,l,ii)-fac_dif*hcore(k,l)
            Fii(l,k)=Fii(k,l)
         enddo
      enddo
      do ia=1, nbasis
         FFp(ia,ii)=0.d0
         do k=1,nbasis
            do l=1,nbasis
               FFp(ia,ii) = FFp(ia,ii) &
                          + vecnat(k,ia)*Fii(k,l)*vecnat(l,ii)
            enddo
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ii,ia)
!$OMP DO
   do ii=1,nbasis
      do ia=1,nbasis
         FF(ia,ii) = FFp(ii,ia)-FFp(ia,ii)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   do m=1,nbasis
      ss=0._dp
      do ii=1,nbasis
         do ia=1,nbasis; if (ii /= ia) then
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
               ss = ss + FF(ii,ia)*ST(ii,ia,m)/Delta_e
            endif
         endif; enddo
      enddo
      B(m)=-ss
   enddo

   deallocate ( QQ, ST, FF, Fii, FFp )
   return
end subroutine make_mats_RDMFT


!========================================
!........>>END OF TOMS CHANGES<<.........
!========================================

! Find the effective orbital using NAG minimization routines. 
! Screening charge constraint already incorporated so minimization
! is unconstrained
subroutine eff_orbital_NAG(iscf,x0,y0,z0,x1,y1,z1,pfile)
!..Global
   use global; use orbocc; use matrices; use functional_m, only: xmix_OEP 
   use energies; implicit none

!..Arguments
   real(dp) :: x0,y0,z0,x1,y1,z1
   character(10) :: pfile
   integer :: iscf

!..Parameters for the NAG routine:
   integer :: nclin, ncnln
   integer :: lda, liwork, lwork, ldr, iuser(1), ldcj 
   integer :: ifail=0, nvar, nstate, iter=0
   real, parameter :: BIGBND=1.e20_dp
   integer, allocatable :: istate(:), iwork(:)
   real(dp), allocatable :: A(:, :), bl(:), bu(:)
   real(dp), allocatable :: X(:), work(:), grd(:)
   real(dp), allocatable :: clamda(:) 
   real(dp), allocatable :: objgrd(:), r(:, :)
   real(dp), allocatable :: c(:), cjac(:,:), Xs(:)
   real(dp), allocatable :: Bmat(:,:),Emat(:),Psimat(:,:),Dveceff(:)
   real(dp) :: user(1), objf

!..Local
   integer :: k, info
   real(dp) :: ss, xnorm
   logical :: l_new_obj
!..TOMS THING
   character(len=30) :: Func_pres , format_string
   integer :: pres

   external En_FUN_vec, e04udm, e04uef, e04ucf, OBJ_LRDMFT, CONFUN1, OBJ_NIK

   nvar=nbasis

   nclin=0
   ncnln=0
!----TOMS CHANGE 
!   ncnln=1
!   ifail=1

   lda=max(1,nclin)
   ldr=nvar
   ldcj=max(1,ncnln)
   nstate=1

   liwork=3*nvar+nclin+2*ncnln
   lwork=2*nvar*nvar+nvar*nclin+2*nvar*ncnln+20*nvar+11*nclin+21*ncnln+1

   allocate ( A(lda,1), bl(nvar+nclin+ncnln), bu(nvar+nclin+ncnln), &
              X(nvar), istate(nvar+nclin+ncnln), c(ldcj), cjac(ldcj,nvar),&
              clamda(nvar+nclin+ncnln), objgrd(nvar), r(ldr,nvar), iwork(liwork),&
              work(lwork), grd(nvar), Xs(Nvar) )
   allocate ( Bmat(nvar,nvar), Emat(nvar), Psimat(nvar,nvar), Dveceff(nvar) )

!..Unconstrained minimization:
   do k=1,nvar
      bl(k)=-BIGBND; bu(k)=BIGBND
   enddo
!....TOMS CHANGE
!   bl(nvar+1)=veffnorm ; bu(nvar+1)=veffnorm
pres=6
if ( pres .gt. 9 ) then 
   format_string= "(A24,I2)"
else 
   format_string = "(A24,I1)"
endif

write (Func_pres,format_string) 'Function Precision = 1d-',pres
print*, 'LOOK HERE'
print*, Func_pres

!  call e04uef("Difference interval = 1d-6")
   call e04uef("Major Iteration Limit = 2000")
   call e04uef("Minor Iteration Limit = 1000")
!  call e04uef("Central difference interval = 1d-4")

!..accuracy of solution 
   call e04uef("Function Precision = 1d-12")
!   call e04uef(Func_pres)
   call e04uef("Optimality Tolerance = 1d-10")

!..How accurately a variable should obey a constraint:
!  call e04uef("Feasibility Tolerance = 1d-10")
!  call e04uef("Nonlinear feasibility = 1d-12")

!..print level
   call e04uef("Major print level = 5")
   call e04uef("Minor print level = 0")
   
!..Initialize effective orbital (first call)   
   if (.not.l_init_efo) then   
      print*,'Effective orbital initialization'
      if( .not. allocated(vec_eff) ) allocate (vec_eff(lnbasis))
      call expand_sqdensity(vec_eff)
      call norm_vec_eff(vec_eff, xnorm)
      l_init_efo=.true.
      call build_Veff()
      H_eff = Hcore + Veff
      call diagon_H(info)    
      call construct_f(0)
      call total_energy_MO(0)
      Print*,'------------- ENERGIES with INITIAL Vec_eff  ---------'
      write(6,113) Bare_energy, Coul_energy, Exch_energy, &
                Totel_energy, Rep_energy, TS_ener, Tot_energy

      Print*,'-------------END INITIAL ENERGIES --------------------'
   endif

!..Initial value
   X=vec_eff

   l_new_obj=.true.

   print*,X
   print*,'-----------------------------'
   if ( .not.l_new_obj ) then !Energy is the objective function

      call e04ucf(nvar, nclin, ncnln, lda, ldcj, ldr, A, &
                  bl, bu, e04udm, En_FUN_vec, iter, istate, c, &
!                  bl, bu, CONFUN1, En_FUN_vec, iter, istate, c, &
                  cjac, clamda, objf, objgrd, r, X, iwork, liwork, &
                  work, lwork, iuser, user, ifail)

      print*,'objgrd', (objgrd(k),k=1,Nvar)
      print*,'----------------------------------------'

   else !Nikitas objective function frozen orbitals

      call e04ucf(nvar, nclin, ncnln, lda, ldcj, ldr, A, &
                  bl, bu, e04udm, OBJ_LRDMFT, iter, istate, c, &
!                 bl, bu, e04udm, OBJ_NIK, iter, istate, c, &
!                  bl, bu, CONFUN1, OBJ_LRDMFT, iter, istate, c, &
                  cjac, clamda, objf, objgrd, r, X, iwork, liwork, &
                  work, lwork, iuser, user, ifail)
   endif

!..Now normalize optimal X
   call norm_vec_eff(X, ss)

!..Final optimal vec_eff
   vec_eff=X
   print*,X
   print*,'-----------------------------'

   if( .not. allocated(vec_eff0) ) allocate(vec_eff0(nbasis))
   if (iscf >1 ) vec_eff = xmix_OEP * vec_eff + (1.d0-xmix_OEP) * vec_eff0
   call norm_vec_eff(vec_eff, ss)
   vec_eff0 = vec_eff

   print*,' Objective function: ',objf
   print*,'---------------------------------------------------------'

  
   call build_Veff()
   H_eff = Hcore + Veff
   call diagon_H(info)    
   call construct_f(0)
   call total_energy_MO(0)
   Print*,'------------- ENERGIES  after veff step ---------'
   write(6,113) Bare_energy, Coul_energy, Exch_energy, &
                Totel_energy, Rep_energy, TS_ener, Tot_energy


113 format('   One-electron     Energy: ',f18.10,/, &
           '   Hartree          Energy: ',f18.10,/, &
           '   Exch./Corr.      Energy: ',f18.10,/, &
           '   Total Electronic Energy: ',f18.10,/, &
           '   Repulsive nuclei Energy: ',f18.10,/, &
           '   Entropy       term (TS): ',f18.10,/, &
           '   *****TOTAL  ENERGY*****: ',f18.10,/) 

   call plot_potential_EFO(x0,y0,z0,x1,y1,z1,'potentialf_EFO')
!  call grad_E_vec_eff_R(x0,y0,z0,x1,y1,z1,'func_ders')
   return
end subroutine eff_orbital_NAG

!------TOMS NONLINEAR CONSTRAINT FOR THE SCREENING DENSITY 
!------------------------------------------------------------------------------
subroutine CONFUN1(MODE,NCNLN,N,LDCJ,NEEDC,X,C,CJAC,NSTATE,IUSER,USER)

use global; use matrices; 
integer :: LDCJ, MODE, N, NCNLN, NSTATE, NEEDC(*), IUSER(*)
real(dp) :: X(N), C(NCNLN), CJAC(LDCJ,N), USER(*)

integer i,j

if (NSTATE.eq.1) then 
   do i=1,NCNLN
      do j=1,N 
         CJAC(i,j) = 0.d0
      enddo
   enddo
endif

if (NEEDC(1) .gt. 0 ) then 
   if (MODE.eq.0 .or. MODE.eq.2) then 
      C(1)=0.d0
      do i=1,N
         do j=1,N
            C(1)=C(1) + X(i) * X(j) * ovlap(i,j)
         enddo
      enddo
   endif
   if (MODE.eq.1 .or. MODE.eq.2 ) then 
      CJAC=0.d0
      do i=1,N
         do j=1,N
            CJAC(1,i)=CJAC(1,i) + X(j) * ovlap(i,j)
         enddo
      enddo
   endif
endif

return
end subroutine CONFUN1


!------------------------------------------------------------------------------
subroutine expand_sqdensity(vec_eff)
!..Global
   use global; use matrices, only:ovlap; use orbocc; use grid_params; use functional_m, only:maxorb; implicit none

!..Arguments
   real(dp),intent(out) :: vec_eff(lnbasis)

!..Local
   real(dp) :: delta(lnbasis,lnbasis), oinv(lnbasis,lnbasis), eigs(lnbasis),vectr(lnbasis,lnbasis)
   real(dp), allocatable :: sqrho(:), b(:)
   real(dp) :: fc, sqrho_e, Obj, x_m, y_m, z_m, sqdif, sqdif_m
   real(dp),external :: electron_density, f_bas
   integer :: k,l,m,n, ig, ia 

   call orbs_on_grid(maxorb)

   allocate ( sqrho(Ngrid), b(lnbasis) )

   fc=(real(nele(3))-1_dp)/real(nele(3))

   print*,'--------------------------------------------------------------------'
   print*,'Initializing Effective Orbital:'
   b=0._dp
   do ig=1, ngrid
      sqrho(ig) = 0._dp
      do ia=1, maxorb
         sqrho(ig) = sqrho(ig) + occnum(ia,3)*orbgridsq(ig,ia)
      enddo
      sqrho(ig) = sqrt(fc*sqrho(ig))
      do k=1,nbasis
         b(k)=b(k) + w_grid(ig)*sqrho(ig)*charg_pot(ig,k)
      enddo
   enddo
   delta=0_dp
   do k=1,nbasis
      delta(k,k)=1_dp
   enddo

   call diagon_lapack(lnbasis,ovlap,delta,eigs,vectr) 

   do k=1,nbasis
      if(abs(eigs(k)) > zero) then
         eigs(k)=1_dp/eigs(k)
      else
         eigs(k)=zero
      endif
   enddo

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,n,k)
!$OMP DO
   do m=1,nbasis
      do n=1,m
         oinv(m,n)=0._dp
         do k=1,nbasis
            oinv(m,n) = oinv(m,n) + vectr(m,k)*eigs(k)*vectr(n,k)
         enddo
         oinv(n,m)= oinv(m,n)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   do k=1,nbasis
      vec_eff(k)=0_dp
      do l=1,nbasis
         vec_eff(k) = vec_eff(k) + oinv(k,l)*b(l)
      enddo
   enddo
   sqdif_m=0_dp
   Obj=0_dp
   do ig=1, ngrid
      sqrho_e=0_dp
      do k=1,nbasis
         sqrho_e=sqrho_e + vec_eff(k)*charg_pot(ig,k)
      enddo
      sqdif=sqrho(ig)-sqrho_e
      if(abs(sqdif_m) < abs(sqdif)) then
         sqdif_m = sqdif; x_m=x_grid(ig); y_m=y_grid(ig); z_m=z_grid(ig)
      endif
!     if( sq < -zero ) print*, x_grid(ig), y_grid(ig), z_grid(ig), sqrho_e
      Obj=Obj+ w_grid(ig)*(sqrho(ig)-sqrho_e)**2
   enddo  
   print*,'Maximum diff from rho:',sqdif_m,'at:',x_m,y_m,z_m

   print*,'Objective Function: ', Obj
   print*,'--------------------------------------------------------------------'

end subroutine expand_sqdensity

!-------------------------------------------------------------------------------
!..Normalize Effective Orbital
subroutine norm_vec_eff(vec_eff, xnorm)
!..Global
   use global; use matrices, only:ovlap; implicit none

!..Arguments
   real(dp) :: vec_eff(lnbasis)

!..Local
   real(dp) :: xnorm,fno
   integer :: i,j

   xnorm=0._dp
   do i=1,nbasis
      do j=1,nbasis
         xnorm=xnorm+vec_eff(i)*ovlap(i,j)*vec_eff(j)
      enddo
   enddo
   fno=sqrt(veffnorm/max(xnorm,zero))
   vec_eff=fno*vec_eff

   return
end subroutine norm_vec_eff

!------------------------------------------------------------------
subroutine build_Veff()
   use global; use matrices, only:vec_eff,Veff; use integrals
   implicit none

!..Local variables
   integer :: m, n, l, s
   integer(8) :: ind_p
   integer :: nrec, nchunk, nrec_last, ichunk, irec
   logical :: inv_pairs, read_int
   real(dp) :: two_int, wmn, wls

   real(dp), allocatable :: DM(:,:) 

   if ( .not.allocated(Veff) ) allocate( Veff(nbasis,nbasis) )
   allocate(DM(lnbasis,lnbasis))

   do m = 1,nbasis
      do n =1, m
         DM(m,n)=vec_eff(m)*vec_eff(n)
         DM(n,m)=DM(m,n)
      enddo
   enddo

   Veff = 0._dp

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

!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(irec,m,n,l,s,two_int,wmn,wls,inv_pairs)
!!$OMP DO
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

!.....Coulomb
      inv_pairs = (l/=m.or.n/=s)
      Veff(m,n) = Veff(m,n) + DM(l,s) * wls * two_int

      if(inv_pairs) Veff(l,s) = Veff(l,s) + DM(m,n) * wmn * two_int

      enddo ! irec
!!$OMP END DO
!!$OMP END PARALLEL

   enddo ! ichunk

   do m=1,nbasis
      do n=1,m
         Veff(n,m)=Veff(m,n)
      enddo
   enddo

   deallocate(DM)
   return

 140 stop 'Effective_orbital: build_Veff: error reading 2-e intergral file'

end subroutine build_Veff
!-------------------------------------------------------------------------------------------
subroutine build_Veff0()
   use global; use matrices, only:vec_eff,Veff; use integrals_mix
   implicit none

   integer :: iint, m, n, l, s
   real(dp) :: twint

   if ( .not.allocated(Veff) ) allocate( Veff(nbasis,nbasis) )

   Veff = 0._dp

   do iint=1, nint_mix
      m=m_mix(iint); n=n_mix(iint); l=l_mix(iint); s=s_mix(iint)
      twint=twin_mix(iint)
      if ( m /= n ) twint=2._dp*twint
      Veff(l,s) = Veff(l,s) + vec_eff(m)*vec_eff(n)*twint
   enddo

   do l=2,nbasis
      do s=1,l-1
         Veff(s,l)=Veff(l,s)
      enddo
   enddo

   return
end subroutine build_Veff0

!----------------------------------------------------------------------------
subroutine En_FUN_vec(mode,n,X,objf,objgrd,nstate,iuser,user)
   use global; use matrices; use energies, ONLY:Tot_energy
   implicit none

!..Arguments
   integer :: mode, n, nstate, iuser(*)
   real(dp) :: X(n), objf, objgrd(n), user(*)

!..Local
   integer :: info, k, l
   real(dp) :: grd(n), fac_o=1._dp
   real(dp) :: qq(n), ss, x1, x3, xnorm

   vec_eff=X
   call norm_vec_eff(vec_eff, xnorm)
   x1=1._dp/sqrt(xnorm)
   x3=x1*x1*x1

   call build_Veff()
   H_eff = Hcore + Veff
   call diagon_H(info) 
   call construct_f(0)

!..The mode cases are described in the manual of NAG e04ucf
   if(mode == 1.or.mode == 2) then !Calculate the gradient
      do k=1,n
         ss=0._dp
         do l=1,n
            ss=ss+ovlap(k,l)*X(l)
         enddo
         qq(k)=ss
      enddo 
      call grd_vec_eff(grd)
      grd=8._dp*grd*fac_o
      ss=0._dp
      do k=1,n
         ss=ss+grd(k)*X(k)
      enddo
      ss=ss*x3
      do k=1,n
         objgrd(k)=veffnorm*(x1*grd(k) - ss*qq(k))
      enddo 
   endif
   if(mode == 0.or.mode == 2) then !Calculate the obj function
      call total_energy_MO(0)
      objf=Tot_energy
!     print*,'E',Tot_energy
   endif
   return
end subroutine En_FUN_vec

!----------------------------------------------------------------------------
subroutine OBJ_NIK(mode,n,X,objf,objgrd,nstate,iuser,user)
   use global; use matrices; implicit none

!..Arguments
   integer :: mode, n, nstate, iuser(*)
   real(dp) :: X(n), objf, objgrd(n), user(*)

!..Local
   integer :: k, l
   real(dp) :: grd(n), qq(n), ss, x1, x3, xnorm, tt2

   vec_eff=X
   call norm_vec_eff(vec_eff, xnorm)
   x1=1._dp/sqrt(xnorm)
   x3=x1*x1*x1

   call build_Veff()

!..The mode cases are described in the manual of NAG e04ucf
   if(mode == 1.or.mode == 2) then !Calculate the gradient
      do k=1,n
         ss=0._dp
         do l=1,n
            ss=ss+ovlap(k,l)*X(l)
         enddo
         qq(k)=ss
      enddo 
      call grd_vec_eff_nik(grd)
      grd=grd
      ss=0._dp
      do k=1,n
         ss=ss+grd(k)*X(k)
      enddo
      ss=ss*x3
      do k=1,n
         objgrd(k)=veffnorm*(x1*grd(k) - ss*qq(k))
      enddo 
   endif
   if(mode == 0.or.mode == 2) then !Calculate the obj function
      call T_2(tt2)
      objf=tt2 
   endif
   return
end subroutine OBJ_NIK

!--------------------------------------------------------------------------------
subroutine T_2(tt2)
   use global; use orbocc; use matrices
   use energies
   implicit none

!..Arguments
   real(dp) :: tt2

!..Local
   integer :: ia,k,l,ii, n_occ
   real(dp) :: ss1, ss2
   real(dp), allocatable :: VV1(:,:),  VV2(:,:)

   allocate ( VV1(nbasis,nbasis), VV2(nbasis,nbasis) )

   n_occ=max(ibond(1),ibond(2))

   VV1=F(:,:,1)-Hcore(:,:)
   VV2=Veff(:,:) 

   tt2=0._dp
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ii,ss1, ss2,k,l) REDUCTION(+:tt2)
!$OMP DO
   do ia=n_occ+1,nbasis
      do ii=1, n_occ
         ss1=0._dp; ss2=0._dp
         do k=1,nbasis
            do l=1,nbasis
                  ss1=ss1+vecnat(k,ia)*VV1(k,l)*vecnat(l,ii)
                  ss2=ss2+vecnat(k,ia)*VV2(k,l)*vecnat(l,ii)
            enddo
         enddo
         tt2=tt2+(ss2**2 - 2*ss1*ss2)/max(ennat(ia)-ennat(ii), small)
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

    print*,'tt2',tt2
    stop

   deallocate ( VV1, VV2 )
   return
end subroutine T_2
!----------------------------------------------------------------------------
subroutine OBJ_LRDMFT(mode,n,X,objf,objgrd,nstate,iuser,user)
   use global; use matrices; implicit none

!..Arguments
   integer :: mode, n, nstate, iuser(*)
   real(dp) :: X(n), objf, objgrd(n), user(*)

!..Local
   integer :: k, l
   real(dp) :: grad(n), qq(n), ss, x1, x3, xnorm, tt2
   real(dp),allocatable :: A(:,:), B(:)

   allocate(A(nbasis,nbasis), B(nbasis))
   vec_eff=X
   call norm_vec_eff(vec_eff, xnorm)
   x1=1._dp/sqrt(xnorm)
   x3=x1*x1*x1

   call make_mats_LRDMFT(A,B)
!..The mode cases are described in the manual of NAG e04ucf
   if(mode == 1.or.mode == 2) then !Calculate the gradient
      do k=1,n
         ss=0._dp
         do l=1,n
            ss=ss+ovlap(k,l)*X(l)
         enddo
         qq(k)=ss
      enddo

      do k=1,nbasis
         grad(k)=0._dp
         do l=1,nbasis
            grad(k)=grad(k)+A(k,l)*vec_eff(l)
         enddo
         grad(k)=-grad(k)+B(k)
      enddo
!     grad=grad
      ss=0._dp
      do k=1,n
         ss=ss+grad(k)*X(k)
      enddo
      ss=ss*x3
      do k=1,n
         objgrd(k)=veffnorm*(x1*grad(k) - ss*qq(k))
      enddo
   endif
   if(mode == 0.or.mode == 2) then !Calculate the obj function
      objf=0._dp
      do k=1,nbasis
         do l=1,nbasis
            objf=objf+vec_eff(k)*A(k,l)*vec_eff(l)
         enddo
      enddo
      objf=-0.5_dp*objf
      do k=1,nbasis
         objf=objf+vec_eff(k)*B(k)
      enddo
   endif
   deallocate(A,B)
   return
end subroutine OBJ_LRDMFT


!----------------------------------------------------------------------------
subroutine OBJ_LRDMFT_old(mode,n,X,objf,objgrd,nstate,iuser,user)
   use global; use matrices; implicit none

!..Arguments
   integer :: mode, n, nstate, iuser(*)
   real(dp) :: X(n), objf, objgrd(n), user(*)

!..Local
   integer :: k, l
   real(dp) :: grad(n), qq(n), ss, x1, x3, xnorm, tt2

   vec_eff=X
   call norm_vec_eff(vec_eff, xnorm)
   x1=1._dp/sqrt(xnorm)
   x3=x1*x1*x1

   call build_Veff()

!..The mode cases are described in the manual of NAG e04ucf
   if(mode == 1.or.mode == 2) then !Calculate the gradient
      do k=1,n
         ss=0._dp
         do l=1,n
            ss=ss+ovlap(k,l)*X(l)
         enddo
         qq(k)=ss
      enddo
      call grd_vec_eff_LRDMFT(grad)
      print*,'grd',grad
!     grad=grad
      ss=0._dp
      do k=1,n
         ss=ss+grad(k)*X(k)
      enddo
      ss=ss*x3
      do k=1,n
         objgrd(k)=veffnorm*(x1*grad(k) - ss*qq(k))
      enddo
   endif
   if(mode == 0.or.mode == 2) then !Calculate the obj function
      call T_LRDMFT(tt2)
      objf=tt2
      print*,'TT2',TT2
   endif
   stop 'a'
   return
end subroutine OBJ_LRDMFT_old

!--------------------------------------------------------------------------------
subroutine T_LRDMFT(tt2)
   use global; use orbocc; use matrices
   use functional_m, only: Dn_lim, small_e
   use energies
   implicit none

!..Arguments
   real(dp) :: tt2

!..Local
   integer :: ia,k,l,ii, n_occ
   real(dp) :: ss, ssc, Delta_n, Delta_e, fac_dif
   real(dp), allocatable :: VV(:,:), FF(:,:) 

   allocate ( VV(nbasis,nbasis), FF(nbasis,nbasis) )

   n_occ=max(ibond(1),ibond(2))

   call make_FF(FF)
 
   tt2=0._dp

!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(ia,ii,ss,k,l,Delta_n,Delta_e) REDUCTION(+:tt2)
!!$OMP DO
   do ia=1,nbasis
      do ii=1, nbasis; if (ii /= ia) then
         ss=0._dp
         do k=1,nbasis
            do l=1,nbasis
                  ss=ss+vecnat(k,ia)*(Veff(k,l))*vecnat(l,ii)
            enddo
         enddo

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
            tt2=tt2+(Delta_n*ss**2 - 2*ss*FF(ia,ii))/Delta_e
         endif
      endif; enddo
   enddo
!!$OMP END DO
!!$OMP END PARALLEL

!    print*,'A_tt2',tt2

   deallocate ( VV, FF )
   return
end subroutine T_LRDMFT

!----------------------------------------------------------------------------------------

subroutine make_FF(FF)
   use global; use orbocc; use matrices
   implicit none

!..Arguments
   real(dp), intent(OUT) :: FF(nbasis,nbasis)
 
   integer :: m,ia,k,l,ii, n_occ
   real(dp) :: ss, Delta_e, Delta_n, fac_dif
   real(dp), allocatable :: QQ(:,:,:), ST(:,:,:), VV(:,:), FFp(:,:) 

   allocate ( VV(nbasis,nbasis), FFp(nbasis,nbasis))

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ii,fac_dif,k,l,VV,ia) 
!$OMP DO
   do ii=1, nbasis
      fac_dif=0.5d0*fac_h(ii)
      do k=1,nbasis
         do l=1,k
            VV(k,l)=F(k,l,ii)-fac_dif*hcore(k,l)
            VV(l,k)=VV(k,l)
         enddo
      enddo

      do ia=1, nbasis
         FFp(ii,ia) = 0.d0
         do k=1,nbasis
            do l=1,nbasis
                  FFp(ii,ia) = FFp(ii,ia) &
                             + vecnat(k,ia)*VV(k,l)*vecnat(l,ii)
            enddo
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

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

   deallocate ( VV, FFp )
   return
end subroutine make_FF

!-------------------------------------------------------------------------------
subroutine grd_vec_eff_LRDMFT(grad)
   use global; use orbocc; use matrices
   use functional_m, only: Functional, Dn_lim, small_e 
   use energies; implicit none

!..Arguments
   real(dp) :: grad(lnbasis), gradcheck(nbasis)

   integer :: m,ia,k,l,ii,n
   real(dp) :: ss, gr, Delta_n, Delta_e, fac_dif, t_minus, t_plus, xlim, tt2
   real(dp), allocatable :: QQ(:,:,:), ST(:,:,:), FF(:,:) 

   allocate ( QQ(nbasis,nbasis,nbasis), ST(nbasis,nbasis,nbasis), &
              FF(nbasis,nbasis))

   xlim=1.e-6_dp

   call make_S_klm()

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,k,ss,l) 
!$OMP DO
   do m=1,nbasis
      do ia=1,nbasis
         do k=1,nbasis
            ss=0._dp
            do l=1,nbasis
               ss = ss + S3(k,l,m)*vecnat(l,ia)
            enddo
            QQ(k,ia,m) = ss
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,ii,k,ss) 
!$OMP DO
   do m=1,nbasis
      do ia=1,nbasis
         do ii=1,nbasis
            ss = 0._dp
            do k=1,nbasis
               ss = ss + QQ(k,ia,m)*vecnat(k,ii)
            enddo
            ST(ii,ia,m) = ss
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

call make_FF(FF)

!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(ia,ii,ss,n,k,l,gr,Delta_n,Delta_e) 
!!$OMP DO
   do n=1,nbasis
      gr=0._dp
      do ia=1,nbasis
         do ii=1, nbasis; if (ii /= ia) then
            ss=0._dp
            do k=1,nbasis
               do l=1,nbasis
                  ss=ss+vecnat(k,ia)*Veff(k,l)*vecnat(l,ii)
               enddo
            enddo
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
              gr=gr+(Delta_n*ss*ST(ii,ia,n)-FF(ia,ii)*ST(ii,ia,n))/Delta_e
            endif
         endif; enddo
      enddo
      grad(n)=4.0_dp*gr
   enddo

!!$OMP END DO
!!$OMP END PARALLEL

! print*,'grad',grad

!   do k=1, nbasis
!      vec_eff(k)=vec_eff(k)-xlim
!      call T_LRDMFT(tt2)
!      t_minus=tt2
!      vec_eff(k)=vec_eff(k)+2*xlim
!      call T_LRDMFT(tt2)
!      t_plus=tt2
!      vec_eff(k)=vec_eff(k)-xlim
!     gradcheck(k)=(t_plus-t_minus)/(2*xlim)
!   enddo

!  print*,'gradcheck',gradcheck


   deallocate ( QQ, ST, FF )

   return
end subroutine grd_vec_eff_LRDMFT

!-------------------------------------------------------------------------------------

subroutine grd_vec_eff_nik(dvec)
   use global; use orbocc; use matrices
   use energies; implicit none

!..Arguments
   real(dp) :: dvec(lnbasis)

   integer :: m,ia,k,l,ii, n_occ
   real(dp) :: ss
   real(dp), allocatable :: QQ(:,:,:), ST(:,:,:), FF(:,:), VV(:,:)

   allocate ( QQ(nbasis,nbasis,nbasis), ST(nbasis,nbasis,nbasis), &
              FF(nbasis,nbasis), VV(nbasis,nbasis) )

   n_occ=max(ibond(1),ibond(2))

   call make_S_klm()

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,k,ss,l) 
!$OMP DO
   do m=1,nbasis
      do ia=n_occ+1,nbasis
         do k=1,nbasis
            ss=0._dp
            do l=1,nbasis
               ss = ss + S3(k,l,m)*vecnat(l,ia)
            enddo
            QQ(k,ia,m) = ss
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,ii,k,ss) 
!$OMP DO
   do m=1,nbasis
      do ia=n_occ+1,nbasis
         do ii=1,n_occ
            ss = 0._dp
            do k=1,nbasis
               ss = ss + QQ(k,ia,m)*vecnat(k,ii)
            enddo
            ST(ii,ia,m) = ss
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   VV(:,:)=F(:,:,1)-Hcore(:,:)-Veff(:,:)

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ii,k,l) 
!$OMP DO
   do ia=n_occ+1,nbasis
      do ii=1, n_occ
         FF(ii,ia) = 0.d0
         do k=1,nbasis
            do l=1,nbasis
                  FF(ii,ia) = FF(ii,ia) + vecnat(k,ia)*VV(k,l)*vecnat(l,ii)
            enddo
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   do m=1,nbasis
      ss=0._dp
      do ii=1,n_occ
         do ia=n_occ+1,nbasis
            ss = ss + FF(ii,ia)*ST(ii,ia,m)/max(ennat(ia)-ennat(ii), small)
         enddo
      enddo
      dvec(m)=-4._dp*ss
   enddo

   deallocate ( QQ, ST, FF )
   return
end subroutine grd_vec_eff_nik

!========================================================================
subroutine plot_potential_EFO(x0,y0,z0,x1,y1,z1,pfile)
   use global; use matrices; use orbocc; use basis_set; use grid_params
   use functional_m, ONLY: maxorb
   use DFT_par; implicit none

!..Arguments
   character(14) :: pfile
   real(dp) :: x0, y0, z0, x1, y1, z1

!..Local
   integer :: ia, ib, ig, ist, N=500
   real(dp) :: x,y,z, xstep, ystep, zstep
   real(dp), external :: electron_density, f_bas
   real(dp), allocatable :: r_d(:), rho(:), dyn(:), dyn_Ha(:), dyn_XC(:), ig_cl(:),scr_charg_plot(:)
   real(dp), allocatable ::  vec_eff_plot(:)
   real(dp), allocatable ::  Sch(:)
   real(dp) :: Scharge, Tcharge, x_dg, y_dg, z_dg, rr, r_min

   allocate ( r_d(N), rho(N) )
   allocate ( dyn(N), dyn_Ha(N), dyn_XC(N), ig_cl(N), scr_charg_plot(N), vec_eff_plot(N) )
   if (.not.allocated(chdens_gr)) then
      allocate ( chdens_gr(ngrid,3) )

      if (.not.allocated(vec_nat_grid)) call orbs_on_grid(maxorb)
      do ig=1,ngrid
         chdens_gr(ig,1)=0._dp; chdens_gr(ig,2)=0._dp
         do ia=1,maxorb
            chdens_gr(ig,1) = chdens_gr(ig,1) + occnum(ia,1)*vec_nat_grid(ig,ia)**2
            chdens_gr(ig,2) = chdens_gr(ig,2) + occnum(ia,2)*vec_nat_grid(ig,ia)**2
         enddo
         chdens_gr(ig,3) = chdens_gr(ig,1)+chdens_gr(ig,2)
      enddo
   endif

   allocate ( Sch(ngrid) )
   Tcharge=0._dp
   do ig=1,ngrid
      Tcharge=Tcharge+ w_grid(ig)*chdens_gr(ig,3)
      Sch(ig)=0._dp
      do ia=1,nbasis
         do ib=1,nbasis
            Sch(ig)= Sch(ig) + vec_eff(ia)*vec_eff(ib)*charg_pot(ig,ib)*charg_pot(ig,ia)
         enddo
      enddo
   enddo

!..Analytic Screening charge
   Scharge=0._dp
   do ia=1,nbasis
      do ib=1,nbasis
         Scharge=Scharge  + vec_eff(ia)*vec_eff(ib)*ovlap_pot(ia,ib)
      enddo
   enddo

   print*,'Screening Charge: ', Scharge
   print*,'    Total Charge: ',Tcharge
   print*,'-------------------------------------------------'

   xstep=(x1-x0)/real(N-1); ystep=(y1-y0)/real(N-1); zstep=(z1-z0)/real(N-1)
   x=x0; y=y0; z=z0
   do ist=1,N
      ig_cl(ist) = 1; r_min=1.e30_dp
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
      dyn(ist) = 0._dp
      dyn_Ha(ist) = 0._dp
      scr_charg_plot(ist)=0._dp
      vec_eff_plot(ist)=0._dp
      do ia=1,nbasis
         vec_eff_plot(ist)=vec_eff_plot(ist) + vec_eff(ia)*f_bas( ia, x, y, z)
      enddo
      scr_charg_plot(ist) = vec_eff_plot(ist)**2
      do ig=1,ngrid
         x_dg = abs(x-x_grid(ig)); y_dg=abs(y-y_grid(ig)); z_dg=abs(z-z_grid(ig)) 
!        if((2.d0*x_dg > dx) .or. (2.d0*y_dg > dy) .or. (2.d0*z_dg > dz)) then
         if(ig /= ig_cl(ist) ) then
            rr= sqrt( x_dg*x_dg + y_dg*y_dg + z_dg*z_dg)
            rr = w_grid(ig)/rr
            dyn_Ha(ist) = dyn_Ha(ist) + rr*chdens_gr(ig,3)
            dyn(ist) = dyn(ist) + rr*Sch(ig)
         endif
      enddo
      r_d(ist) = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo

   open(unit=53, file=pfile, status='unknown')
   do ist=1,N
      dyn_XC(ist) = dyn(ist) - dyn_Ha(ist) 
      write(53,'(f10.3,3e20.10)') r_d(ist), dyn(ist), dyn_Ha(ist), dyn_XC(ist)
   enddo
   close(53)

   open(unit=191, file='density_EFO', status='unknown')
      do ist=1,N
         write(191,'(3f20.10)') r_d(ist), r_d(ist)*r_d(ist)*rho(ist), r_d(ist)*r_d(ist)*scr_charg_plot(ist)
      enddo
   close(191)

   open(unit=191, file='Eff_orbital', status='unknown')
      do ist=1,N
         write(191,'(3f20.10)') r_d(ist), vec_eff_plot(ist)
      enddo
   close(191)


   deallocate (r_d, rho, Sch )

end subroutine plot_potential_EFO

!-------------------------------------------------------------------------------------------
subroutine make_S_klm()
   use global; use matrices, only:vec_eff,S3; use integrals
   implicit none

   integer :: m, n, l, s
   integer(8) :: ind_p
   integer :: nrec, nchunk, nrec_last, ichunk, irec
   logical :: inv_pairs, read_int
   real(dp) :: two_int

   if ( .not.allocated(S3) ) allocate( S3(nbasis,nbasis,nbasis) )
   S3 = 0._dp

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

!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(irec,m,n,l,s,two_int,wmn,wls,inv_pairs)
!!$OMP DO
   do irec = 1, nrec
      m=mu(irec)
      n=nu(irec)
      l=lambda(irec)
      s=sigma(irec)

      two_int = twoin(irec)

      S3(l,s,m) = S3(l,s,m) + two_int*vec_eff(n)
      if ( m /= n ) S3(l,s,n) = S3(l,s,n) + two_int*vec_eff(m)

      inv_pairs = (l/=m.or.n/=s)
      if ( inv_pairs ) then
         S3(m,n,l) = S3(m,n,l) + two_int*vec_eff(s)
         if ( l /= s) S3(m,n,s) = S3(m,n,s) + two_int*vec_eff(l)
      endif
      
   enddo ! irec
!!$OMP END DO
!!$OMP END PARALLEL

   enddo ! ichunk

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(l,m,n)
!$OMP DO
   do l=1,nbasis
      do m=1,nbasis
         do n=1,m
            S3(n,m,l)=S3(m,n,l)
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   return

 140 stop 'Effective_orbital: build_Veff: error reading 2-e intergral file'

end subroutine make_S_klm

!-----------------------------------------------------------------------------------
subroutine grd_vec_eff(dvec)
   use global; use orbocc; use matrices
   use energies
   implicit none

!..Arguments
   real(dp) :: dvec(lnbasis)

   integer :: m,ia,k,l,ii, n_occ
   real(dp) :: ss
   real(dp), allocatable :: QQ(:,:,:), ST(:,:,:), Am(:,:), bm(:), FF(:,:) 

   allocate ( QQ(nbasis,nbasis,nbasis), ST(nbasis,nbasis,nbasis), &
              Am(nbasis, nbasis), bm(nbasis), FF(nbasis,nbasis) )

   n_occ=max(ibond(1),ibond(2))

   call make_S_klm()

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,k,l,ss) 
!$OMP DO
   do m=1,nbasis
      do ia=n_occ+1,nbasis
         do k=1,nbasis
            ss=0._dp
            do l=1,nbasis
               ss = ss + S3(k,l,m)*vecnat(l,ia)
            enddo
            QQ(k,ia,m) = ss
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,ii,k,ss) 
!$OMP DO
   do m=1,nbasis
      do ia=n_occ+1,nbasis
         do ii=1,n_occ
            ss = 0._dp
            do k=1,nbasis
               ss = ss + QQ(k,ia,m)*vecnat(k,ii)
            enddo
            ST(ii,ia,m) = ss
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,l,ss,ia,ii) 
!$OMP DO
   do k=1,nbasis
      do l=1,nbasis
         ss = 0._dp
         do ia=n_occ+1,nbasis
            do ii=1, n_occ
               ss = ss + ST(ii,ia,k)*ST(ii,ia,l) / max(ennat(ia)-ennat(ii), small)
            enddo
         enddo
         Am(k,l) = ss
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ii,k,l) 
!$OMP DO
   do ia=n_occ+1,nbasis
      do ii=1, n_occ
         FF(ii,ia) = 0.d0
         do k=1,nbasis
            do l=1,nbasis
!                 FF(ii,ia) = FF(ii,ia) + vecnat(k,ia)*(F(k,l,1)-Hcore(k,l)-Veff(k,l))*vecnat(l,ii)
                  FF(ii,ia) = FF(ii,ia) + vecnat(k,ia)*(F(k,l,1)-Hcore(k,l))*vecnat(l,ii)
!                 FF(ii,ia) = FF(ii,ia) + vecnat(k,ia)*(F(k,l,1)-H_eff(k,l))*vecnat(l,ii)
            enddo
         enddo
      enddo
   enddo   
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(m,ia,ii,ss) 
!$OMP DO
   do m=1,nbasis
      ss=0._dp
      do ii=1,n_occ
         do ia=n_occ+1,nbasis
            ss = ss + FF(ii,ia)*ST(ii,ia,m)/max(ennat(ia)-ennat(ii), small)
         enddo
      enddo
      Bm(m)=ss
   enddo
!$OMP END DO
!$OMP END PARALLEL

   do m=1,nbasis
      dvec(m)=0._dp
      do l=1,nbasis
         dvec(m)=dvec(m) + Am(m,l)*vec_eff(l)
      enddo
      dvec(m)= dvec(m) - Bm(m)
   enddo

   deallocate ( QQ, ST, Am, bm, FF )
   return
end subroutine grd_vec_eff

!--------------------------------------------------------------------------------------------------
! The improved DIIS extrapolation method: P. Pulay, J. Comput. Chem. 3, 556 (1982).
! icall (in/out) : how many times DIIS_extrapolation has been called before/ number of present call
! Ndiis (in) : number of previous runs to be used / dimension of B matrix
! H (in): the hamiltonian to be extrapolated
! H_extrap (out): the extrapolated Hamiltonian
! emax : the maximum error of the e vector
subroutine DIIS_f(icall, Ndiis, Nd, f_in, f_out, f_old, df_old, ov, emax)
!..Global
   use global; implicit none

!..Arguments
   integer, intent(in) :: Ndiis, Nd
   integer, intent(inout) :: icall
   real(dp), intent(in) :: f_in(Nd)
   real(dp), intent(in) :: ov(Nd,Nd)
   real(dp), intent(out) :: f_out(Nd)
   real(dp), intent(inout):: f_old(Nd,0:Ndiis), df_old(Nd,Ndiis)
   real(dp), intent(out) :: emax

!..Local
   integer Np, ip, jp, k, m, l
   real(dp), allocatable :: B(:,:), C(:) !, U(:,:), eigB(:), Binv(:,:), delta_mat(:,:)
   real(dp) :: xnorm, Bs

   icall=icall+1
   Np=min(icall,Ndiis) 
   
   allocate( B(Np+1, Np+1), C(Np+1) )
!  allocate( B(Np, Np), U(Np, Np), C(Np), eigB(Np), Binv(Np, Np), delta_mat(Np,Np) )

!  delta_mat=0._dp
!  do k=1,Np
!     delta_mat(k,k)=1._dp
!  enddo

!..Save present f_in in f_old, create df_old, new B
   if( icall > Ndiis) then  !Shift f_old by -1 to empty f_old(:,ip) 
      do ip=1,Np-1
         f_old(:,ip)=f_old(:,ip+1)
         df_old(:,ip)=df_old(:,ip+1)
      enddo
   endif
 
!..Save f_in in table 
   f_old(:,Np) = f_in(:) 
   df_old(:,Np) = f_old(:,Np) - f_old(:,Np-1)

   xnorm=0._dp
   do l=1,Np
      xnorm=xnorm+df_old(l,Np)*df_old(l,Np)
   enddo
   xnorm = sqrt(veffnorm/xnorm)
!  df_old(:,Np) = xnorm*df_old(:,Np)

   if (mod(icall,ndiis+1) == 0) then 
   do ip=1,Np
      do jp=1,ip
         B(jp,ip) = 0._dp
         do m=1,Nd
            do l=1,Nd
               B(jp,ip) = B(jp,ip)+df_old(m,jp)*ov(m,l)*df_old(l,ip)
            enddo
         enddo
         B(ip,jp) = B(jp,ip)
      enddo
   enddo

   do k=1,Np
      B(k,Np+1) = -1._dp
      B(Np+1,k) = -1._dp
      C(k) =0._dp
   enddo
   B(Np+1, Np+1) = 0._dp
   C(Np+1)=-1._dp
 
!  print*,'Matrix B:'
!  write(6,'(7f20.10)')((B(m,l),m=1,Np+1), l=1,Np+1)

!  call diagon_NAG(Np,B,delta_mat,eigB,U)
!  print*,'Matrix B eigenvalues:'
!  write(6,'(6e20.10)')(eigB(m),m=1,Np)

!  do k=1,Np
!     if(abs(eigB(k)) > zero ) then
!        eigB(k)=1._dp/eigB(k)
!     else
!        if(eigB(k) > 0 )  then 
!           eigB(k)=1._dp/zero
!        else
!           eigB(k)=-1._dp/zero
!        endif
!     endif
!  enddo

!  do m=1,Np
!     do k=1,Np
!        Binv(m,k)=0.0_dp
!        do l=1, Np
!           Binv(m,k) = Binv(m,k) + U(m,l)*eigB(l)*U(k,l)
!        enddo
!     enddo
!  enddo

!   do m=1,Np
!     do k=1,Np
!        Bs=0._dp
!        do l=1,Np
!           Bs=Bs+B(m,l)*Binv(l,k)
!        enddo
!        print*,m,k,Bs
!     enddo
!  enddo  

!  Bs=0._dp
!  do m=1,Np
!     do k=1,Np
!        Bs=Bs+Binv(m,k)
!     enddo
!  enddo

!  do k=1,Np
!     C(k)=0._dp
!     do l=1,Np
!        C(k) = C(k) + Binv(k,l)
!     enddo
!     C(k) = C(k) /  Bs
!  enddo 

!..Solve linear system B X = C 
   call lin_sqsym_solve(Np+1, B, C) 
!  print*,'C:'
!  write(6,'(7e20.10)')(C(m),m=1,Np+1)

!  xnorm=0._dp
!  do k=1,Np
!     xnorm=xnorm+C(k) 
!  enddo
!  print*,'Sum of C(k):',xnorm

   do m=1,Nd
      f_out(m)= 0._dp
      do k=1,Np
         f_out(m) = f_out(m) + f_old(m,k)*C(k) 
      enddo
   enddo   

   print*,' '
!  print*,'f_old:'
!  do k=1,ndiis
!  write(6,'(7e20.10)')(f_old(m,k),m=1,Nd)
!  enddo

!  print*,'fout:'
!  write(6,'(7e20.10)')(f_out(m),m=1,Nd)


!  xnorm=0._dp
!  do m=1,Nd
!     xnorm=xnorm+f_out(m)**2
!  enddo
!  print*,'intermediate norm', sqrt(xnorm)
!  xnorm=sqrt(veffnorm/xnorm)
!  f_out(:) = xnorm*f_out(:)

   else
      f_out = f_in
   endif
!  emax=0._dp
!  do k=1,Nd
!     emax=emax+f_out(k)-f_in(k)
!  enddo

end subroutine DIIS_f
