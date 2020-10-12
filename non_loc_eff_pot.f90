! Find the effective orbital using NAG minimization routines. 
subroutine eff_pot_NAG(iscf)
!..Global
   use global; use orbocc; use matrices; use functional_m, only:maxorb, svd_cut, xmix_OEP 
   use energies
   implicit none

!..Arguments
   integer :: iscf

!..Parameters for the NAG routine:
   integer :: nclin, ncnln
   integer :: lda, liwork, lwork, ldr, iuser(1), ldcj 
   integer :: ifail=0, nvar, nstate, ia, iter=0
   real, parameter :: BIGBND=1.e20_dp
   integer, allocatable :: istate(:), iwork(:)
   real(dp), allocatable :: A(:, :), bl(:), bu(:)
   real(dp), allocatable :: X(:), work(:), grd(:)
   real(dp), allocatable :: clamda(:) 
   real(dp), allocatable :: objgrd(:), r(:, :), Q1(:,:)
   real(dp), allocatable :: c(:), cjac(:,:), H_eff_t(:,:)
   real(dp) :: user(1), objf, dl, sav

!..Local
   integer :: mode, k, m, n, l, info, iii, ivec, mm,mm2
   real(dp) :: xlam, ss, xnorm, DE, E_0
   real(dp), allocatable :: Veff_t(:,:)
   logical :: l_new_obj, ffile

   external En_FUN_vec, e04udm, e04uef, e04ucf, OBJ_NL

   nvec_eff=1 
   nvar=nbasis*nvec_eff

   nclin=0
   ncnln=0

   lda=max(1,nclin)
   ldr=nvar
   ldcj=max(1,ncnln)
   nstate=1

   liwork=3*nvar+nclin+2*ncnln
   lwork=2*nvar*nvar+nvar*nclin+2*nvar*ncnln+20*nvar+11*nclin+21*ncnln+1

   allocate ( A(lda,1), bl(nvar+nclin+ncnln), bu(nvar+nclin+ncnln), &
              X(nvar), istate(nvar+nclin+ncnln), c(ldcj), cjac(ldcj,nvar),&
              clamda(nvar+nclin+ncnln), objgrd(nvar), r(ldr,nvar), iwork(liwork),&
              work(lwork), grd(nvar), Q1(nbasis,nbasis), Veff_t(nbasis,nbasis) )

!..Unconstrained minimization:
   do k=1,nvar
      bl(k)=-BIGBND; bu(k)=BIGBND
   enddo

!  call e04uef("Difference interval = 1d-6")
   call e04uef("Major Iteration Limit = 5000")
   call e04uef("Minor Iteration Limit = 2000")
!  call e04uef("Central difference interval = 1d-4")

!..accuracy of solution 
  call e04uef("Function Precision = 1d-2")
   call e04uef("Optimality Tolerance = 1d-3")

!..How accurately a variable should obey a constraint:
   call e04uef("Feasibility Tolerance = 5d-8")
!  call e04uef("Nonlinear feasibility = 1d-12")

!..print level
   call e04uef("Major print level = 0")
   call e04uef("Minor print level = 0")
   
!..Initialize effective orbital (first call)   
   if (.not.l_init_efo) then   
      if ( .not. allocated(vec_eff_n) ) allocate(vec_eff_n(lnbasis,nvec_eff) )
!     vec_eff_n(:,1)=  vecnat(:,2) !+vecnat(:,3)+vecnat(:,4)+vecnat(:,5)
!     vec_eff_n(:,2)=  vecnat(:,1) 
!     vec_eff_n(:,3)=  vecnat(:,3) 
!     vec_eff_n(:,4)=  vecnat(:,4) 
!     vec_eff_n(:,5)=  vecnat(:,6) 
!     vec_eff_n(:,6)=  vecnat(:,7) 
        vec_eff_n(:,1:nvec_eff)=vecnat(:,1:nvec_eff)
      do ivec=1,nvec_eff
         call  vec_eff_norm(ivec, vec_eff_n, xnorm)
         print*,'Normalization of vec_eff_n(:,ivec):',ivec,xnorm
      enddo

!     inquire ( file='vec_eff.dat', exist=ffile )
!     if ( ffile ) then
!        open(unit=971, file='vec_eff.dat', status='unknown' )
!        read(971,*)vec_eff_n
!        close(971)
!        print*,vec_eff_n
!     endif
      if ( .not.allocated(DMM) ) allocate( DMM(lnbasis,lnbasis) )
      if ( .not.allocated(Vcoul) ) allocate( Vcoul(lnbasis,lnbasis) )
      if ( .not.allocated(Vexch) ) allocate( Vexch(lnbasis,lnbasis) )
      if ( .not. allocated(Veff) ) allocate( Veff(lnbasis,lnbasis) ) 
      if ( .not. allocated(H_eff_t) ) allocate(H_eff_t (lnbasis,lnbasis) ) 

      Veff=0._dp
      do ivec=1,nvec_eff
         do k=1,nbasis
            do l=1,nbasis
               Q1(k,l)=vec_eff_n(k,ivec)*vec_eff_n(l,ivec)
            enddo
         enddo
         call sum_intg_1(Q1,Veff_t)
         Veff= Veff+ Veff_t
      enddo

      do iii=1,nconv
         do k=1,nbasis
            do l=1,nbasis
               DMM(k,l)=0._dp
               do ia=1,nbasis
                  DMM(k,l)=DMM(k,l) + occnum(ia,3)*vecnat(k,ia)*vecnat(l,ia)
               enddo
            enddo
         enddo
         call sum_intg_2(DMM,Vcoul)
         call sum_intg_1(DMM,Vexch)
         H_eff = Hcore + Vcoul - 0.5_dp*Vexch - Veff 
         if (iii .gt. 1) H_eff = x_mix_1*H_eff + (1._dp-x_mix_1)*H_eff_t
         H_eff_t = H_eff
         call diagon_H(info)
         call construct_f(0)
         call total_energy_MO(0)
         if (iii .gt. 1) then
            DE=abs(Tot_energy-E_0)
            if ( DE < de_conv ) goto 5
         endif
         E_0=Tot_energy
      enddo
5     print*,'1',iii,'DE=',DE,'E_0=',E_0
      Print*,'------------- ENERGIES with INITIAL Vec_eff  ---------'
      write(6,113) Bare_energy, Coul_energy, Exch_energy, &
                Totel_energy, Rep_energy, TS_ener, Tot_energy

      Print*,'-------------END INITIAL ENERGIES --------------------'
      l_init_efo=.true.
   endif


!..Initial value
   do ivec=1,nvec_eff
      mm=(ivec-1)*nbasis+1
      mm2=ivec*nbasis 
      X(mm:mm2)=vec_eff_n(:,ivec)
   enddo

   call e04ucf(nvar, nclin, ncnln, lda, ldcj, ldr, A, &
               bl, bu, e04udm, OBJ_NL, iter, istate, c, &
               cjac, clamda, objf, objgrd, r, X, iwork, liwork, &
               work, lwork, iuser, user, ifail)

!..Final optimal vec_eff
   do ivec=1,nvec_eff
      mm=(ivec-1)*nbasis+1
      mm2=ivec*nbasis 
      vec_eff_n(:,ivec)=X(mm:mm2)
   enddo

   do ivec=1,nvec_eff
      call  vec_eff_norm(ivec,vec_eff_n, xnorm)
      print*,'Normalization of vec_eff_n:',ivec,xnorm
   enddo

   print*,' Objective function: ',objf
   print*,'---------------------------------------------------------'

!  if( .not. allocated(vec_eff0_n) ) allocate(vec_eff0_n(nbasis_pot,nvec_eff))
!  if (iscf >1 ) vec_eff = xmix_OEP * vec_eff + (1.d0-xmix_OEP) * vec_eff0
!  call norm_vec_eff(vec_eff, ss)
!  vec_eff0 = vec_eff

   Veff=0._dp
   do ivec=1,nvec_eff
      do k=1,nbasis
         do l=1,nbasis
            Q1(k,l)=vec_eff_n(k,ivec)*vec_eff_n(l,ivec)
         enddo
      enddo
      call sum_intg_1(Q1,Veff_t)
      Veff= Veff+ Veff_t
   enddo

   do iii=1,nconv
   do k=1,nbasis
      do l=1,nbasis
         DMM(k,l)=0._dp
         do ia=1,nbasis
            DMM(k,l)=DMM(k,l) + occnum(ia,3)*vecnat(k,ia)*vecnat(l,ia)
         enddo
      enddo
   enddo

   call sum_intg_2(DMM,Vcoul)
   call sum_intg_1(DMM,Vexch)

   H_eff = Hcore + Vcoul - 0.5_dp*Vexch - Veff
   if (iii .gt. 1) H_eff = x_mix_1*H_eff + (1._dp-x_mix_1)*H_eff_t
   H_eff_t = H_eff

   call diagon_H(info)    

      call construct_f(0)
      call total_energy_MO(0)
      if (iii .gt. 1) then
         DE=abs(Tot_energy-E_0)
         if ( DE < de_conv ) goto 15
      endif
      E_0=Tot_energy
!        print*,'b',Tot_energy

   enddo !iii
15 print*,'End Iter',iii,'DE=',DE,'E_0=',E_0
   Print*,'------------- ENERGIES with FINAL Vec_eff  ---------'
   write(6,113) Bare_energy, Coul_energy, Exch_energy, &
             Totel_energy, Rep_energy, TS_ener, Tot_energy

   Print*,'-------------END FINAL ENERGIES --------------------'

   open(unit=971, file='vec_eff.dat', status='unknown' )
      write(971,*) vec_eff_n
      close(971)


113 format('   One-electron     Energy: ',f18.10,/, &
           '   Hartree          Energy: ',f18.10,/, &
           '   Exch./Corr.      Energy: ',f18.10,/, &
           '   Total Electronic Energy: ',f18.10,/, &
           '   Repulsive nuclei Energy: ',f18.10,/, &
           '   Entropy       term (TS): ',f18.10,/, &
           '   *****TOTAL  ENERGY*****: ',f18.10,/) 
   return
end subroutine eff_pot_NAG

!----------------------------------------------------------------------------
subroutine OBJ_NL(mode,n,X,objf,objgrd,nstate,iuser,user)
   use global; use matrices; use orbocc; use energies, ONLY:Tot_energy
   implicit none

!..Arguments
   integer :: mode, n, nstate, iuser(*)
   real(dp) :: X(n), objf, objgrd(n), user(*)

!..Local
   integer :: info, k, l,m, ia, iii, ivec, mm, mm2
   real(dp) :: grd(n), dl=1e-8_dp, fac_o=1._dp, fac_g=1._dp
   real(dp) :: ss, x1, x3, xnorm, Q1(nbasis,nbasis), Veff_t(nbasis,nbasis,nvec_eff)
   real(dp) :: Veff_n(nbasis,nbasis), grd_t(nbasis)
   real(dp) :: ddd=1.e-5_dp, objf0

   real(dp) :: H_eff_t(nbasis,nbasis), DE, E_0

   do ivec=1,nvec_eff
      mm=(ivec-1)*nbasis+1
      mm2=ivec*nbasis 
      vec_eff_n(:,ivec)=X(mm:mm2)
   enddo

   Veff=0._dp
   do ivec=1,nvec_eff
      do k=1,nbasis
         do l=1,nbasis
            Q1(k,l)=vec_eff_n(k,ivec)*vec_eff_n(l,ivec)
         enddo
      enddo
      call sum_intg_1(Q1,Veff_n)
      Veff= Veff+ Veff_n
   enddo

   do iii=1,nconv
   do k=1,nbasis
      do l=1,nbasis
         DMM(k,l)=0._dp
         do ia=1,nbasis
            DMM(k,l)=DMM(k,l) + occnum(ia,3)*vecnat(k,ia)*vecnat(l,ia)
         enddo
      enddo
   enddo

   call sum_intg_2(DMM,Vcoul)
   call sum_intg_1(DMM,Vexch)

   H_eff = Hcore + Vcoul - 0.5_dp*Vexch - Veff
   if (iii .gt. 1) H_eff = x_mix_1*H_eff + (1._dp-x_mix_1)*H_eff_t
   H_eff_t = H_eff

   call diagon_H(info) 
   call construct_f(0)
   call total_energy_MO(0)
   if (iii .gt. 1) then
      DE=abs(Tot_energy-E_0)
      if ( DE < de_conv ) goto 5
   endif
   E_0=Tot_energy

   enddo !iii
 5 print*,'3',iii,'DE=',DE,'E_0=',E_0

!..The mode cases are described in the manual of NAG e04ucf
   if(mode == 1.or.mode == 2) then !Calculate the gradient
      call total_energy_MO(0)
      objf0=Tot_energy
!     print*,'objf0',objf0

!      goto 131
!.....Numerical gradient
!      do ivec=1,nvec_eff
!      do m=1,nbasis
!         mm=(ivec-1)*nbasis+m
!         X(mm)=X(mm)+ddd
!         vec_eff_n(m,ivec)=X(mm)
!         do k=1,nbasis
!            do l=1,nbasis
!               Q1(k,l)=vec_eff_n(k,ivec)*vec_eff_n(l,ivec)
!            enddo
!         enddo
!         call sum_intg_1(Q1,Veff_n)
!         Veff= Veff-Veff_t(:,:,ivec)+Veff_n

!         do iii=1,nconv
!         do k=1,nbasis
!            do l=1,nbasis
!               DMM(k,l)=0._dp
!               do ia=1,nbasis
!                  DMM(k,l)=DMM(k,l) + occnum(ia,3)*vecnat(k,ia)*vecnat(l,ia)
!               enddo
!            enddo
!         enddo

!         call sum_intg_2(DMM,Vcoul)
!         call sum_intg_1(DMM,Vexch)
!         H_eff = Hcore + Vcoul - 0.5_dp*Vexch - Veff
!         if (iii .gt. 1) H_eff = x_mix_1*H_eff + (1._dp-x_mix_1)*H_eff_t
!         H_eff_t = H_eff
       
!         call diagon_H(info) 
!         call construct_f(0)
!         call total_energy_MO(0)
!         if (iii .gt. 1) then
!            DE=abs(Tot_energy-E_0)
!            if ( DE < de_conv ) goto 15
!         endif
!         E_0=Tot_energy
!         enddo !iii
!15       continue !print*,'i',iii,'DE=',DE,'E_0=',E_0
!         objf=Tot_energy
!         Veff= Veff-Veff_n+Veff_t(:,:,ivec)
         
!         objgrd(m)=(objf-objf0)/ddd
!         X(mm)=X(mm)- ddd
!         vec_eff_n(m,ivec)=X(mm)
!      enddo !m
!      enddo !ivec

!      Veff=0._dp
!      do ivec=1,nvec_eff
!         do k=1,nbasis
!            do l=1,nbasis
!               Q1(k,l)=vec_eff_n(k,ivec)*vec_eff_n(l,ivec)
!            enddo
!         enddo
!         call sum_intg_1(Q1,Veff_n)
!         Veff= Veff+ Veff_n
!      enddo

!      do iii=1,nconv
!      do k=1,nbasis
!         do l=1,nbasis
!            DMM(k,l)=0._dp
!            do ia=1,nbasis
!               DMM(k,l)=DMM(k,l) + occnum(ia,3)*vecnat(k,ia)*vecnat(l,ia)
!            enddo
!         enddo
!      enddo

!      call sum_intg_2(DMM,Vcoul)
!      call sum_intg_1(DMM,Vexch)

!      H_eff = Hcore + Vcoul - 0.5_dp*Vexch - Veff
!      if (iii .gt. 1) H_eff = x_mix_1*H_eff + (1._dp-x_mix_1)*H_eff_t
!      H_eff_t = H_eff
      
!      call diagon_H(info) 
!      call construct_f(0)
!      call total_energy_MO(0)
!      if (iii .gt. 1) then
!         DE=abs(Tot_energy-E_0)
!         if ( DE < de_conv ) goto 25
!      endif
!      E_0=Tot_energy
!     print*,'5',Tot_energy
!      enddo !iii
!25 continue !print*,'5',iii,'DE=',DE,'E_0=',E_0
!....End numerical gradient

!131   continue

      call calc_Lagrange_mult_all

      do ivec=1,nvec_eff
         call GRAD_OBJ(ivec,grd_t)
         mm=(ivec-1)*nbasis+1
         mm2=ivec*nbasis
         objgrd(mm:mm2)=grd_t
      enddo 
!     do k=1,nbasis
!        print*,k,objgrd(k),grd(k), objgrd(k)/grd(k)
!     enddo
!      stop
   endif
   if(mode == 0.or.mode == 2) then !Calculate the obj function
      call total_energy_MO(0)
      objf=Tot_energy
!     print*,'E',Tot_energy
   endif

   return
end subroutine OBJ_NL

!-------------------------------------------------------------------------------

subroutine GRAD_OBJ(ivec,dk)

   use global; use matrices; use orbocc; use functional_m, only:small_e; implicit none

   integer :: ia, ib, l, k, ivec
   real(dp), allocatable :: Q(:,:), T(:,:)
   real(dp) :: dk(nbasis), Delta_e

   allocate( Q(lnbasis, lnbasis), T(lnbasis, lnbasis) )
!Q calculation

   do l=1,nbasis
      do k=1,l
         Q(l,k)=0._dp
         do ia=1,nbasis
            do ib=1,nbasis
               if (ia/=ib) then
                  Delta_e=ennat(ia) - ennat(ib) 
                  if(abs(Delta_e) < small_e) then
                     if ( Delta_e >= 0 ) then 
                         Delta_e=small_e
                     else
                         Delta_e=-small_e
                     endif
                  endif
                  Q(l,k)=Q(l,k)+vecnat(l,ia)*vecnat(k,ib)*(e_ab(ia,ib)-e_ab(ib,ia))/Delta_e
               endif
            end do 
         end do
         Q(l,k)=4._dp*Q(l,k)
         Q(k,l)=Q(l,k)
      end do
   end do

   call sum_intg_1(Q,T) 

!dk calculation
   do k=1,nbasis
       dk(k)=0._dp
       do l=1,nbasis
          dk(k)= dk(k)+T(k,l)*vec_eff_n(l,ivec)
       end do
!      print*,'a',k,dk(k)
   end do

   deallocate ( Q, T )

end subroutine GRAD_OBJ
subroutine sum_intg_1(Q,T)
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
   real(dp) :: Q(lnbasis, lnbasis),  T(lnbasis, lnbasis)

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

   T=0._dp

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

      do irec = 1, nrec
         m=mu(irec)
         n=nu(irec)
         l=lambda(irec)
         s=sigma(irec)

         two_int = twoin(irec)

!.....The exchange should be added to different T elements:
!.....Exchange
    
         inv_pairs = (l/=m.or.n/=s)
         if(m>=l) T(m,l) = T(m,l) + two_int * Q(n,s)
         if(inv_pairs.and.l>=m) T(l,m) = T(l,m) + two_int * Q(n,s)
         if(m/=n) then
            if(n>=l) T(n,l) = T(n,l) + two_int * Q(m,s)
            if(inv_pairs.and.l>=n) T(l,n) = T(l,n) + two_int * Q(m,s)

            if(l/=s) then
               if(n>=s) T(n,s) = T(n,s) + two_int * Q(m,l)
               if(inv_pairs.and.s>=n) T(s,n)=T(s,n) + two_int * Q(m,l)
            endif !(l/=s)
         endif !(m/=n)

         if(l/=s) then
            if(m>=s) T(m,s) = T(m,s) + two_int * Q(n,l)
            if(inv_pairs.and.s>=m) T(s,m) = T(s,m) + two_int * Q(n,l)
         endif !l/=s
            
      enddo ! irec

   enddo ! ichunk

   do m=1,nbasis
      do n=1,m
         T(n,m)=T(m,n)
      enddo
   enddo

   return

 140 stop 'construct_f:sum_intg: error reading 2-e intergral file'

end subroutine sum_intg_1

subroutine sum_intg_2(Q,T)
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
   real(dp) :: Q(lnbasis, lnbasis),  T(lnbasis, lnbasis)

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

   T=0._dp

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
         T(m,n) = T(m,n) + Q(l,s) * wls * two_int

         if(inv_pairs) T(l,s) = T(l,s) + Q(m,n) * wmn * two_int

      enddo ! irec

   enddo ! ichunk

   do m=1,nbasis
      do n=1,m
         T(n,m)=T(m,n)
      enddo
   enddo

   return

 140 stop 'construct_f:sum_intg: error reading 2-e intergral file'

end subroutine sum_intg_2

subroutine non_loc_eff_pot_HF(iscf)
!..Global
   use global; use matrices; use orbocc; use functional_m; implicit none

!..Arguments
   integer :: iscf
   real(dp) :: xmix=0.1_dp

!..Local
   integer :: i, ia, ib
   real(dp) :: Heff(lnnatorb,lnnatorb)
   real(dp), allocatable, save :: Heff_old(:,:)
   real(dp) :: vectr(lnnatorb,lnnatorb), diag(lnnatorb)

   if(.not.allocated(Heff_old) ) allocate( Heff_old(lnnatorb,lnnatorb) )

   call calc_Lagrange_mult_all()   

   call RHF(occnum, fac_h, fac_c, fac_e)

   call construct_f(1)
   Heff = F(:,:,1)
   call diagon_lapack(lnbasis,Heff,ovlap,ennat,vecnat)

   do ia=1,nbasis
     do ib=1,nbasis
 if( ia == ib ) then
           Heff(ia,ib)=ennat(ia)
        else
           if(ia>ib) then 
              Heff(ia,ib)= xa*(e_ab(ia, ib)-e_ab(ib, ia))
           else
              Heff(ia,ib)= xa*(e_ab(ib, ia)-e_ab(ia, ib))
           endif
        endif 
     end do
   end do

!   Heff=Heff/real(xnele(3))
!   if( iscf > 1) Heff=xmix*Heff+(1.d0-xmix)*Heff_old
!   Heff_old=Heff
   print*,'Non local potential HF'

!  call diagon_Vab(Heff, vectr, diag)
   call diagon_lapack(lnbasis,Heff,ovlap,diag,vectr) 
   call assign_eigfunctions(vectr, diag)   

   call calc_occ_factors(occnum, fac_h, fac_c, fac_e)
   call construct_f(0)

end subroutine non_loc_eff_pot_HF

! Baldsiefen Gross effective Hamiltonian
subroutine non_loc_eff_pot_BG(Temper,shift,iscf)
!..Global
   use global; use matrices; use orbocc; use functional_m; use energies; implicit none

!..Arguments
   integer :: iscf
   real(dp) :: Temper, shift

!..Local
   integer :: ia, ib, ispin, inn
   real(dp) :: Heff(lnnatorb,lnnatorb), arg, diff_inv, diff
   real(dp) :: occ_ia, de, e_old, DE_Dn(lnnatorb)
   real(dp) :: vec_old(lnbasis,lnnatorb), en_old(lnnatorb), DE_Dn_1(lnnatorb)
   real(dp) :: vectr(lnnatorb,lnnatorb), diag(lnnatorb), occ_t(lnnatorb,3)
   real(dp), allocatable, save :: Heff_old(:,:)
   real(dp) :: dde, aexp, occ_ib, tau, tt, dnn, fr
   real(dp) :: DE_Dn_plus, DE_Dn_minus, xpos, xneg, scale_f
   logical :: orb_temp
   save e_old, scale_f

   tau=2.0_dp
   orb_temp=.false.

   if(iscf == 1) then
      e_old=1d20
      scale_f=1.e2_dp
   endif
   print*,'iscf',iscf,'scale_f',scale_f

   print*,'Baldsiefen effective potential:'
    
   if(.not.allocated(Heff_old)) allocate(Heff_old(lnnatorb, lnnatorb))

!..Aplies to analytic calculation of second derivative
   if(orb_temp) then
      if(Functional == 'BB0' ) then
         aexp=0.5_dp
      elseif(Functional == 'POW' ) then
         aexp=alpha
      else
!
      endif
   endif
!  xmix=1.0d0
   ispin=1

   do inn=1,15
      call calc_Lagrange_mult_all()

      call Func_der_n_NAG(occnum, DE_Dn, ispin)

      if (orb_temp) then
         dnn=small*0.9_dp
         do ia=1,nnatorb
            occ_t(ia,1)=max(min(occnum(ia,1),1._dp-small),small)
            occ_t(ia,2)=max(min(occnum(ia,2),1._dp-small),small)
            occ_t(ia,3)= occ_t(ia,1) + occ_t(ia,2)
         enddo
      endif

!!..The generalized Fock Matrix
!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(ia,ib,occ_ia,arg,dde,occ_ib,occ_t,DE_Dn_1,ispin,DE_Dn_plus,DE_Dn_minus,tt,diff,diff_inv)
!!$OMP DO
      do ia=1, nnatorb
         occ_ia=max(min(occnum(ia,1),1.d0-small),small)

         if(orb_temp) then
            if(Functional == 'POW' .or. Functional == 'BB0') then
!...........Calculation of dde (d^2E / dn_a^2 ) analytically (Muller, Power)
               dde=(1._dp-aexp*aexp*occ_ia**(2._dp*(aexp-1._dp)))*CoulNO(ia,ia) 
               do ib=1,nnatorb
                  occ_ib=max(min(occnum(ib,1),1._dp-small),small)
                  dde=dde-aexp*(aexp-1._dp)*occ_ia**(aexp-2._dp)*occ_ib**aexp*ExchNO(ia,ib)
               enddo
            else
!...........Calculation of dde (d^2E / dn_a^2 ) numerically (valid for any functional)
               occ_t(ia,1) = occ_t(ia,1) + dnn        
               occ_t(ia,3) = occ_t(ia,1) + occ_t(ia,2)
               call Func_der_n_NAG(occ_t, DE_Dn_1, ispin)
               DE_Dn_plus=DE_Dn_1(ia)
               occ_t(ia,1) = occ_t(ia,1) - 2._dp*dnn
               occ_t(ia,3) = occ_t(ia,1) + occ_t(ia,2)
               call Func_der_n_NAG(occ_t, DE_Dn_1, ispin)
               DE_Dn_minus=DE_Dn_1(ia)
               dde= (DE_Dn_plus-DE_Dn_minus)/(dnn*2._dp)
               occ_t(ia,1) = occ_t(ia,1) + dnn
               occ_t(ia,3) = occ_t(ia,1) + occ_t(ia,2)
            endif
             tt=occ_ia*(1._dp-occ_ia)*dde/tau
         else
            tt=Temper
         endif

         arg=occ_ia/(1._dp-occ_ia)
         Heff(ia,ia) = shift-tt*log(arg) + DE_Dn(ia) 
         do ib=1, ia-1
            diff=occnum(ia,1)-occnum(ib,1)
            if (abs(diff) >= 1.e-8_dp) then
               diff_inv=1._dp/diff
            else
               diff_inv=0._dp
            endif
            Heff(ia,ib) = (e_ab(ib,ia) - e_ab(ia,ib))*diff_inv
!...........Now scale down off-diagonal elements
            if(abs(Heff(ia,ib)) > scale_f) then
               print*,'big:',ia,ib,Heff(ia,ib)
               Heff(ia,ib)=Heff(ia,ib)*scale_f/abs(Heff(ia,ib))
            endif
            Heff(ib,ia) =  Heff(ia,ib)
         enddo !ib
      enddo !ia=1,nnatorb
!!$OMP END DO
!!$OMP END PARALLEL

!     if ( iscf > 1) then 
!        Heff = xmix*Heff + (1.d0-xmix)*Heff_old
!     endif

      Heff_old=Heff; en_old = ennat; vec_old = vecnat
      call diagon_Vab(Heff, vectr, diag)
      call assign_eigfunctions(vectr, diag)
      call construct_f(0)
      call total_energy_MO(1)
      de=Tot_energy-e_old
      print*,inn,' Energy:',Tot_energy, 'DE:',de
      e_old=Tot_energy

      if ( de < 0 ) then
         xneg = xneg + 1._dp
      else
         xpos = xpos + 1._dp
      endif

   enddo ! inn=1,15

   fr = xneg/max(xpos,small)
   if (fr > 1.5_dp ) then
      scale_f = 1.01_dp*scale_f
!     Temper=0.99*Temper
   elseif ( fr > 1.1_dp ) then
!     Temper=Temper
      scale_f = scale_f
   else
!     Temper=1.01*Temper
      scale_f = 0.95_dp* scale_f 
   endif

end subroutine non_loc_eff_pot_BG


!=====================================================================================
!  PIRIS non local potential method
subroutine non_loc_eff_pot(iscf)
!..Global
   use global; use matrices; use orbocc; use functional_m; use energies; implicit none

!..Arguments
   integer :: iscf

!..Local
   real(dp) :: F_ab(lnnatorb,lnnatorb), v(lnbasis,lnbasis), v_extrap(lnbasis,lnbasis)
   real(dp) :: vectr(lnnatorb,lnnatorb), diag(lnnatorb)
   real(dp) :: vec_old(lnbasis,lnnatorb), en_old(lnnatorb)
   real(dp) :: DD(lnbasis)
   real(dp),save :: scale_f, e_old
   real(dp) :: de, xneg, xpos
   real(dp), save :: fr
   real(dp), external :: Theta
   integer :: ia, ib, k, l, inn, icall

   en_old = ennat; vec_old = vecnat

   if(iscf == 1) then
      e_old=1e20_dp
      scale_f=1e-2_dp
!  elseif(mod(iscf,1) == 0) then
!     scale_f = 0.9*scale_f
   endif

   print*,'ZETA:',scale_f, 'fr=', fr

   if(iscf == 1 .and. iguess /= 3) then
      F_ab=0._dp
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ib,k,l,DD)
!$OMP DO
      do ia=1,nnatorb
         do k=1, nbasis
            DD(k)=0._dp
            do l=1, nbasis
               DD(k)=DD(k)+F(k,l,ia)*vecnat(l,ia)
            enddo
         enddo
         do ib=1,nnatorb
            F_ab(ia,ib) = 0._dp
            do k=1, nbasis
               F_ab(ia,ib) = F_ab(ia,ib) + vecnat(k,ib)*DD(k)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      do ia=1,nnatorb
         do ib=1, ia
            v(ia,ib)=0.5_dp*(F_ab(ia,ib) + F_ab(ib,ia))
            v(ib,ia)=v(ia,ib)
         enddo
      enddo
      call diagon_Vab(v, vectr, diag)
      call assign_eigfunctions(vectr,diag)   
      call construct_f(0)
   else !if(iscf == 1 .and. iguess /= 3)
      xpos=0._dp; xneg=0._dp
      icall=0
      do inn=1,15
         F_ab = 0._dp
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ib,k,l,DD)
!$OMP DO
         do ia=1,nnatorb
            do k=1, nbasis
               DD(k)=0._dp
               do l=1, nbasis
                  DD(k)=DD(k)+F(k,l,ia)*vecnat(l,ia)
               enddo
            enddo
            do ib=1,nnatorb
               do k=1, nbasis
                   F_ab(ia,ib) = F_ab(ia,ib) + vecnat(k,ib)*DD(k)
               enddo
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL

         do ia=1, nnatorb
            v(ia,ia) = ennat(ia)
            do ib=1, nnatorb; if(ia /= ib) then
               v(ia,ib) = Theta(ia-ib)*(F_ab(ia,ib) - F_ab(ib,ia))&
                        + Theta(ib-ia)*(F_ab(ib,ia) - F_ab(ia,ib))
               if( abs(v(ia,ib)) > scale_f) then
!                    print*,'For ',ia,ib,'V=',v(ia,ib)
                     v(ia,ib)=v(ia,ib)*scale_f/abs(v(ia,ib))
               endif
            endif; enddo
         enddo
         v_extrap=v
         call diagon_Vab(v_extrap, vectr, diag)
         call assign_eigfunctions(vectr,diag)   
         call construct_f(0)
         call total_energy_MO(1)
         de=Tot_energy-e_old
         print*,inn,' Energy:',Tot_energy, 'DE:',de
         e_old=Tot_energy
         if (de < 0._dp) then
            xneg = xneg + 1._dp
         else
            xpos = xpos + 1._dp
         endif
      enddo
      fr = xneg/max(xpos,small)
      if (fr > 1.5d0 ) then
         scale_f = 1.01*scale_f
      elseif ( fr > 1.1_dp ) then
         scale_f = scale_f
      else
         scale_f = 0.95_dp* scale_f 
      endif
   endif


!  print*,'Eigenvalues:'
!  write(6,'(10f10.5)')(diag(ia),ia=1,nnatorb)
!  print*,'Eigenvectors:'
!  do ia=1,nnatorb
!     print*,'ia=',ia
!     write(6,'(10f10.5)')(vectr(ia,ib),ib=1,nnatorb)
!  enddo
   

end subroutine non_loc_eff_pot
!=====================================================================================
function Theta(i)
use global; implicit none
integer :: i
real(dp) :: Theta
   if(i>=0) then
      Theta = 0._dp
   else
      Theta = 1._dp
   endif 
end function Theta
!======================================================================================

subroutine assign_eigfunctions(vectr, diag)
   use global; use orbocc; implicit none

!..Arguments 
   real(dp) :: vectr(lnnatorb, lnnatorb)
   real(dp) :: diag(lnnatorb)

!..Local 
   real(dp) :: vecnat_new(lnbasis,lnnatorb)
   integer :: i, ia, ib

!..Expand eigenvectors in original fixed basis
   do ia=1,nnatorb
      do i=1,nbasis
         vecnat_new(i,ia) = 0._dp
         do ib=1,nnatorb
            vecnat_new(i,ia)=vecnat_new(i,ia)+ vecnat(i,ib)* vectr(ib,ia)
         enddo
      enddo
   enddo

!..Assign eigenvectors according to overlap with old ones

   vecnat = vecnat_new
   ennat= diag

end subroutine assign_eigfunctions

!======================================================================================
subroutine diagon_Vab(V_ab, vectr, diag)
   use global; implicit none

!..Arguments 
   real(dp), intent(in) :: V_ab(lnnatorb,lnnatorb)
   real(dp), intent(out) :: vectr(lnnatorb,lnnatorb), diag(lnnatorb)

!..Local
   integer :: ia, ib

!..For NAG routine
   integer :: ITYPE=1, LWORK, IFAIL 
   character(1) :: JOB='V', UPLO='U'
   real(dp), allocatable :: A(:,:), B(:,:), W(:), WORK(:)

   LWORK = 64*lnnatorb+1

   allocate( A(lnnatorb,lnnatorb), B(lnnatorb,lnnatorb), W(lnnatorb), WORK(LWORK) )

   do ia=1,nnatorb
      do ib=1,nnatorb
         B(ia,ib) = 0._dp
      enddo
      B(ia,ia) = 1._dp
   enddo

   A = V_ab

   IFAIL = 0
!  print*,nnatorb

   call F02FDF( ITYPE, JOB, UPLO, lnnatorb, A, lnnatorb, B, lnnatorb, W, &
                WORK, LWORK, IFAIL)

   if(ifail/=0) then
      print*,'ifail',ifail
      stop 'non_loc_eff_pot_piris:diagon_Vab: ifail not zero'
   endif

   vectr = A
   diag = W

   deallocate ( A, B, W, WORK )

end subroutine diagon_Vab

!.................................................
!..Norm of Effective Orbital
subroutine vec_eff_norm(ivec, vec_eff_n, xnorm)
!..Global
   use global; use matrices, only:ovlap,nvec_eff; implicit none

!..Arguments
   real(dp) :: vec_eff_n(lnbasis,nvec_eff)
   integer :: ivec

!..Local
   real(dp) :: xnorm,fno
   integer :: i,j

   xnorm=0._dp
   do i=1,nbasis
      do j=1,nbasis
         xnorm=xnorm+vec_eff_n(i,ivec)*ovlap(i,j)*vec_eff_n(j,ivec)
      enddo
   enddo
   fno=sqrt(veffnorm/max(xnorm,zero))
!  vec_eff_n=fno*vec_eff_n

   return
end subroutine vec_eff_norm

!------------------------------------------------------------------

