!-------------------------------------------------------------------------------
! Implements in RDMFT the model potential:
! Ryabinkin, Kananenka, Staroverov,  Phys. Rev. Lett 111, 013001 (2013) (for HF)
! Ryabinkin, Kohut, Staroverov,  Phys. Rev. Lett 115, 083001 (2015) (for any WF) 
!-------------------------------------------------------------------------------
subroutine model_local_pot(x0,y0,z0,x1,y1,z1)
   use global; use orbocc; use grid_params; use matrices; use functional_m; use energies; implicit none

!..Arguments
   real(dp), intent(in) :: x0,y0,z0,x1,y1,z1

!..Local
   integer :: ig, ia, ib,k, l, iter, Niter, Nist=1500, ist, icall=0, Ndiis, Nmix, mode, maxorb_0, IP_meth, id, jd
   real(dp) :: x, y, z, delt=1.e10_dp, delta_rho, emax=1.e10_dp
   real(dp), allocatable :: rho(:), pot_grid_0(:), pot_grid(:), Vloc_0(:,:), Vloc(:,:), rho_WF(:)
   real(dp), allocatable :: DD(:,:,:), H_0(:,:), H_eff_0(:,:), I_aver_WF(:), tau_WF(:), I_aver(:), tau(:)
   real(dp), allocatable :: rho_p(:), pot_p(:), pot_p_0(:), rho_p_WF(:), H_extrap(:,:)
   real(dp), allocatable :: I_aver_p_WF(:), tau_p_WF(:), r_d(:)
   real(dp), allocatable :: I_aver_p(:), tau_p(:), x_p(:), y_p(:), z_p(:), vec_can(:,:), ennat_can(:)
   real(dp) :: xstep, ystep, zstep, xmix_sta, conv, conv_diis, Total_Energy, e_IP
   real(dp) :: x_infty(1),y_infty(1),z_infty(1),rho_infty(1),I_aver_infty(1),tau_infty(1), E_kin
   real(dp), external :: electron_density
   logical :: print_in_all_steps

!..Parameters
   Niter=500 ! Number of iterations for Staroverov potential
   Nmix=50! Number of steps with only mixing (before DIIS)
   xmix_sta=0.2_dp !Mixing for convergence of Staroverov potential
   print_in_all_steps=.false. !Print Staroverov potential in all steps
   Ndiis=6 ! Number of previous Hamiltonians used in the extrapolation
   conv=1e-9_dp ! convergence in potential change
   conv_diis=1e-12_dp ! convergence in DIIS error

   IP_meth=1 ! 1: usual EKT, 2: asymptotic behavior of epsilon 3: shift
   
   call read_basis()
   call orbs_on_grid(maxorb) 
  
   allocate ( rho(ngrid), rho_WF(ngrid), pot_grid_0(ngrid), pot_grid(ngrid) )
   allocate ( I_aver_WF(ngrid), tau_WF(ngrid), I_aver(ngrid), tau(ngrid) )
   allocate ( Vloc_0(lnbasis,lnbasis), Vloc(lnbasis,lnbasis), DD(lnbasis,lnbasis, lnnatorb) )
   allocate ( H_0(lnbasis, lnbasis), H_eff_0(lnbasis,lnbasis), H_extrap(lnbasis,lnbasis) )
   allocate ( rho_p(Nist), pot_p_0(Nist), pot_p(Nist), rho_p_WF(Nist) )
   allocate ( I_aver_p(Nist), tau_p(Nist), x_p(Nist), y_p(Nist), z_p(Nist) )
   allocate ( I_aver_p_WF(Nist), tau_p_WF(Nist), r_d(Nist), vec_can(lnbasis, lnnatorb), ennat_can(lnnatorb) )

!..The density matrix
   call construct_DD(nbasis,maxorb,vecnat,DD)

!..Generate initial effective Hamiltonian (except xc part)
   call Calc_H0(DD,H_0) 

!..The density 
   do ig=1,ngrid
!     rho_WF(ig) = 0._dp
!     do ia=1,maxorb
!        rho_WF(ig) = rho_WF(ig) + occnum(ia,3)*orbgridsq(ig,ia)
!     enddo
!     rho_WF(ig) = max(rho_WF(ig),small)
      rho_WF(ig) = electron_density(x_grid(ig),y_grid(ig),z_grid(ig),3)
   enddo
 
!..Find the Ionization Potential
   x_infty=0._dp; y_infty=0._dp; z_infty=10._dp
   if (functional == 'RHF' .or. functional == 'RHA' ) then
      mode=0 ! Use the KS/HF definition of I_aver
      rho_infty(1)=electron_density(x_infty(1),y_infty(1),z_infty(1),3)  
      call average_local_IE_tau(mode, 1,x_infty,y_infty,z_infty,rho_infty,DD,I_aver_infty,tau_infty) 
      print*,'*********************************************************'
      print*,'       E_HOMO HF = ', ennat(nele(1))
      print*,'Iaver at infinity= ', I_aver_infty, ' a.u.;  Tau at infinity= ', tau_infty(1)
      print*,'*********************************************************'
      if (IP_meth == 1) then
         e_IP=ennat(nele(1))
         print*,'IP target selected from Koopmans:', e_IP, ' a.u.'
         print*,'*********************************************************'
      else
         e_IP=I_aver_infty(1)
         print*,'IP target selected from Iaver at Infinity:', e_IP, ' a.u.'
         print*,'*********************************************************'
      endif
   else !(functional /= 'RHF')
      mode=1 ! Use the WF for the calculation of I_aver 
      call canonical_orbitals(IP_meth, e_IP) ! EKT 
      rho_infty(1)=electron_density(x_infty(1),y_infty(1),z_infty(1),3)  
      call average_local_IE_tau(mode, 1,x_infty,y_infty,z_infty,rho_infty,DD,I_aver_infty,tau_infty) 
      print*,'*********************************************************'
      print*,'  E_HOMO from EKT=', e_IP, ' a.u.'
      print*,'Iaver at infinity= ', I_aver_infty(1), ' a.u.;  Tau at infinity= ', tau_infty(1)
      print*,'*********************************************************'
      if (IP_meth == 1) then
         print*,'IP target selected from EKT=', e_IP, ' a.u.'
         print*,'*********************************************************'
      else
         e_IP=I_aver_infty(1)
         print*,'IP target selected from Iaver at Infinity:', e_IP, ' a.u.'
         print*,'*********************************************************'
      endif
   endif !(functional == 'RHF')

!..Slater potential, I_aver_WF, tau_WF on the grid
   pot_grid=0._dp
   call Slater_potential(Ngrid,x_grid,y_grid,z_grid,pot_grid,rho_WF) 
   call average_local_IE_tau(mode, Ngrid,x_grid,y_grid,z_grid,rho_WF,DD,I_aver_WF,tau_WF)

   do ig=1,ngrid
      pot_grid_0(ig) = pot_grid(ig) - I_aver_WF(ig) + tau_WF(ig)
   enddo
!..Now: pot_grid is the Slater potential; pot_grid_0 is a constant that we add quntities below

!..FOR PLOTTING 
!..Slater potential, I_aver_WF, tau_WF on the plotting line
   xstep=(x1-x0)/real(Nist-1); ystep=(y1-y0)/real(Nist-1); zstep=(z1-z0)/real(Nist-1)
   x=x0; y=y0; z=z0
   do ist=1,Nist
      x_p(ist)=x; y_p(ist)=y; z_p(ist)=z
      rho_p_WF(ist) = electron_density(x,y,z,3)
!     rho_p_WF(ist) = max(small,rho_p_WF(ist))
      r_d(ist) = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo
   pot_p=0._dp
   call Slater_potential(Nist,x_p,y_p,z_p,pot_p,rho_p_WF)
   call average_local_IE_tau(mode, Nist,x_p,y_p,z_p,rho_p_WF,DD,I_aver_p_WF,tau_p_WF)
   open(unit=861,file='Staroverov.dat',status='unknown')
   open(unit=860,file='Slater.dat',status='unknown')
   do ist=1,Nist
      pot_p_0(ist)=pot_p(ist) - I_aver_p_WF(ist) + tau_p_WF(ist)
      write(860,'(4e20.10)')r_d(ist), pot_p(ist), I_aver_p_WF(ist), tau_p_WF(ist)
!     rho_p_WF(ist)*4.d0*Pi*r_d(ist)**2, tau_p_WF(ist)*rho_p_WF(ist)*4.d0*Pi*r_d(ist)**2
   enddo
   close(860)
!..END PLOTTING

!..Matrix elements of the model potential (grid integration) / add to effective Hamiltonian
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,l,ig)
!$OMP DO 
   do k=1,nbasis
      do l=1,k
         Vloc(k,l)=0._dp
         do ig=1,ngrid
            Vloc(k,l)=Vloc(k,l)+w_grid(ig)*pot_grid(ig)*bas_f_grid(ig,k)*bas_f_grid(ig,l)
         enddo
         Vloc(l,k)=Vloc(k,l)
         H_eff_0(k,l) = H_0(k,l)+Vloc(k,l)
         H_eff_0(l,k) = H_eff_0(k,l)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!..Now read the initial guess (starting orbitals for iteration 1)
!..Orbitals and energies will replace WF values
   call read_vec_e()

!..Idempotent KS system:
   maxorb_0=maxorb
   mode=0; maxorb=nele(1); occnum(1:maxorb,1)=1._dp; occnum(1+maxorb:nnatorb,1)=zero
   occnum(:,2)=occnum(:,1); occnum(:,3)=occnum(:,1)+occnum(:,2)
   call calc_occ_factors(occnum, fac_h, fac_c, fac_e)

!..Iterate the potential 
   do iter=1, Niter
      print*,'----------------------------------------------------------------------------'
      print*,'Iteration:', iter

!.....Redefine DD
      call construct_DD(nbasis,maxorb,vecnat,DD)

!.....Generate effective Hamiltonian (except xc part)
      call Calc_H0(DD,H_0) 
   
!.....Shift orbital energies
      if ( iter > 1) then
         do ia=1,maxorb
            ennat(ia)=ennat(ia)-ennat(maxorb) + e_IP! Shift energies to match EKT IP
         enddo
      endif

!.....Recalculate the density, I_aver, tau, for the new orbitals, ennat
      do ig=1, ngrid
         rho(ig) = electron_density(x_grid(ig),y_grid(ig),z_grid(ig),3)
      enddo
      
      call average_local_IE_tau(mode,Ngrid,x_grid,y_grid,z_grid,rho,DD,I_aver,tau)
         
      delta_rho=0._dp
      do ig=1, ngrid
         pot_grid(ig) = pot_grid_0(ig) + I_aver(ig) - tau(ig) 
         delta_rho=delta_rho+w_grid(ig)*abs(rho(ig)-rho_WF(ig))
      enddo
      delta_rho=delta_rho/real(nele(3))

!.....Print Staroverov potential in every step (1 set for xmgrace)
      if ( print_in_all_steps ) then
         do ist=1,Nist
            rho_p(ist) = electron_density(x_p(ist),y_p(ist),z_p(ist),3)
         enddo
         call average_local_IE_tau(mode,Nist,x_p,y_p,z_p,rho_p,DD,I_aver_p,tau_p)
         do ist=1,Nist
            pot_p(ist)=pot_p_0(ist) + I_aver_p(ist) - tau_p(ist) 
            write(861,'(5e20.10)')r_d(ist), pot_p(ist), I_aver_p(ist), tau_p(ist)
         enddo
         write(861,*) ' '
      endif

!..Matrix elements of the model potential (grid integration) / add to effective Hamiltonian
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,l,ig)
!$OMP DO 
      do k=1,nbasis
         do l=1,k
            Vloc(k,l)=0._dp
            do ig=1,ngrid
               Vloc(k,l)=Vloc(k,l)+w_grid(ig)*pot_grid(ig)*bas_f_grid(ig,k)*bas_f_grid(ig,l)
            enddo
            Vloc(l,k)=Vloc(k,l)
            H_eff(k,l) = H_0(k,l)+Vloc(k,l)
            H_eff(l,k) = H_eff(k,l)
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      if (iter < Nmix) then
!.....Just Mixing
         if (iter > 1) then
            H_eff=xmix_sta*H_eff+(1._dp - xmix_sta)*H_eff_0
         endif
         H_eff_0=H_eff
      else
!.....Apply DIIS convergence extrapolation
         call DIIS_extrapolation(icall, Ndiis, H_eff, H_extrap, emax) 
         H_eff_0=H_eff
         H_eff=H_extrap
      endif

!.....The effective KS Hamiltonian/diagonalization   
      call diagon_lapack(lnbasis,H_eff,ovlap,ennat,vecnat)
      
      print*,'e_IP',e_IP
      print*,'Eigenvalues:'
      write(6,'(5F20.12)')(ennat(ia),ia=1,maxorb)
      print*, ' '
      write(6,'(" E_HOMO: ",F20.12," a.u.,",F20.12," eV")')ennat(maxorb), ennat(maxorb)*ha_ev
      print*, ' '

!.....Check convergence
      if(iter > 1) then
         delt=0._dp
         do k=1,nbasis
            do l=1,k
               delt=delt+(Vloc(k,l)-Vloc_0(k,l))**2
            enddo
         enddo
         delt=2._dp*sqrt(delt)/real(nbasis*(nbasis-1))
         print*,'Delt:',delt,'   Delta_rho:',delta_rho
         if(iter> Nmix) then
            print*,'Maximum DIIS error:',emax
            print*,' '
            if ( iter > Nmix) print*,'Emax,conv_diis:',emax , conv_diis
         endif
         if ( delt < conv .or. emax < conv_diis ) exit
      endif
      Vloc_0 = Vloc
      print*,'............................................................................'
   enddo ! iter=1,Niter

!..Print final Staroverov potential in unit 861 (Staroverov.dat):
   do ist=1,Nist
      rho_p(ist) = electron_density(x_p(ist),y_p(ist),z_p(ist),3)
   enddo
   call average_local_IE_tau(mode,Nist,x_p,y_p,z_p,rho_p,DD,I_aver_p,tau_p)
   do ist=1,Nist
      pot_p(ist)=pot_p_0(ist) + I_aver_p(ist) - tau_p(ist)
      write(861,'(5e20.10)')r_d(ist), pot_p(ist), I_aver_p(ist), tau_p(ist)
   enddo
   call Kinetic_energy(E_kin)
   Print*,'T_s: ',E_kin

   close(861)
   maxorb=maxorb_0
   deallocate ( rho, rho_WF, pot_grid_0, pot_grid, I_aver_WF, tau_WF, I_aver, tau, Vloc_0)
   deallocate ( Vloc, DD, H_0, H_eff_0, H_extrap, rho_p, pot_p_0, pot_p, rho_p_WF)
   deallocate ( I_aver_p, tau_p, I_aver_p_WF, tau_p_WF, r_d, x_p, y_p, z_p )

   print*,'Target IP=', e_IP

end subroutine model_local_pot
!-----------------------------------------------------------------------------------------------
! Calculates the slater potential (pot) at the point (x,y,z) in real space
subroutine Slater_potential(Np,x,y,z, pot, rho)
   use global; use grid_params; use orbocc; use functional_m, ONLY:maxorb;implicit none

!..Arguments
   integer, intent(in) :: Np
   real(dp), intent(in) :: x(Np),y(Np),z(Np)
   real(dp), intent(in) :: rho(Np)
   real(dp), intent(out) :: pot(Np)

!..Local
   integer :: ia, ib, ig , k, l, ip
   real(dp) :: xint, rr, ff, ss, xx, yy, zz
   real(dp) :: vc(lnnatorb), vc2(lnnatorb), r(ngrid), x_int_kl(lnbasis, lnbasis), x_int_ka(lnbasis, lnnatorb)
   real(dp) :: xint_b(lnnatorb), xint_c, fac_h_HF(lnnatorb), fac_c_HF(lnnatorb, lnnatorb)
   real(dp) ::  fac_e_HF(lnnatorb,lnnatorb), fac_n(lnnatorb,lnnatorb)
   real(dp), external :: f_bas, electron_density
   logical :: analytic_integral

!..To do analytically the integrals I_kl(r)= \int dr_1 \frac{\chi_k(r_1) \chi_l(r_1)}{|r-r_1|}
!..where chi_k and chi_l are basis functions.
   analytic_integral=.true.

   call RHF(occnum, fac_h_HF, fac_c_HF, fac_e_HF)
   fac_n=fac_c-fac_c_HF
   if ( analytic_integral ) then

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip, l, k, ia, ib, ff, vc, vc2, xint, rr, xx, yy, zz, ss, r, x_int_kl, x_int_ka, xint_b, xint_c)
!$OMP DO 
      do ip=1,Np
         pot(ip)=0._dp
         xx=x(ip); yy=y(ip); zz=z(ip)

         vc=0._dp

            do l=1,nbasis
               ff=f_bas(l, xx, yy, zz)
               do ia=1,maxorb
                  vc(ia)=vc(ia)+vecnat(l,ia)*ff
               enddo
            enddo

            do ia=1,maxorb
               vc2(ia)=vc(ia)*vc(ia)
            enddo

            do k=1,nbasis
               do l=1,k
                  xint=0._dp
                  call Xintegral2(k,l,xx,yy,zz,xint)
                  x_int_kl(k,l)=xint
                  x_int_kl(l,k)=xint      
               enddo
            enddo
            do ia=1,maxorb
               do k=1,nbasis
                  ss=0._dp
                  do l=1,nbasis
                     ss=ss+vecnat(l,ia)*x_int_kl(k,l)
                  enddo
                  x_int_ka(k,ia)=ss
               enddo
            enddo
            do ia=1,maxorb
               xint_b(ia)=0._dp
               do k=1,nbasis
                  xint_b(ia)=xint_b(ia)+x_int_ka(k,ia)*vecnat(k,ia)
               enddo
            enddo
            do ia=1,maxorb
               do ib=1,maxorb
                  xint=0._dp
                  do k=1,nbasis
                     xint=xint+vecnat(k,ia)*x_int_ka(k,ib)
                  enddo
                  xint = fac_e(ia,ib)*vc(ia)*vc(ib)*xint
                  xint_c= fac_n(ia,ib)*vc2(ia)*xint_b(ib)
                  pot(ip)=pot(ip) + xint + xint_c
               enddo
            enddo
         pot(ip)=2._dp*pot(ip)/rho(ip)
      enddo
!$OMP END DO
!$OMP END PARALLEL

   else ! not analytic_integral

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip, l, k, ia, ib, ig, ff, vc, xint, rr, xx, yy, zz, ss, r, vc2, xint_b, xint_c )
!$OMP DO 
      do ip=1,Np
         pot(ip)=0._dp
         xx=x(ip); yy=y(ip); zz=z(ip)

         vc=0._dp
   
            do l=1,nbasis
               if (Np == Ngrid) then
                  ff = bas_f_grid(ip,l) 
               else
                  ff=f_bas(l, xx, yy, zz)
               endif
               do ia=1,maxorb
                  vc(ia)=vc(ia)+vecnat(l,ia)*ff
               enddo
            enddo

            do ia=1,maxorb
               vc2(ia)=vc(ia)*vc(ia)
            enddo

            do ig=1,ngrid
               r(ig)=0._dp
               rr=sqrt( (x_grid(ig)-xx)**2&
                      + (y_grid(ig)-yy)**2&
                      + (z_grid(ig)-zz)**2)
               if (rr > small) r(ig)=w_grid(ig)/rr
            enddo
            do ia=1,maxorb
               xint_b(ia)=0._dp
               do ig=1,ngrid
                  if(r(ig)>zero) xint_b(ia)=xint_b(ia)+r(ig)*vec_nat_grid(ig,ia)*vec_nat_grid(ig,ia)
               enddo
            enddo
            do ia=1,maxorb
               do ib=1,maxorb
                  xint=0._dp
                  do ig=1,ngrid
                     if(r(ig) > zero) xint=xint+r(ig)*vec_nat_grid(ig,ia)*vec_nat_grid(ig,ib)
                  enddo
!                 xint= occnum(ia,3)*occnum(ib,3)*vc(ia)*vc(ib)*xint
                  xint= fac_e(ia,ib)*vc(ia)*vc(ib)*xint
                  xint_c= fac_n(ia,ib)*vc2(ia)*xint_b(ib)
                  pot(ip) = pot(ip) + xint + xint_c
               enddo
            enddo
!        pot(ip)=-pot(ip)/(2.d0*rho(ip))
         pot(ip)=2._dp*pot(ip)/rho(ip)
      enddo
!$OMP END DO
!$OMP END PARALLEL
   endif

   return
end subroutine Slater_potential

!----------------------------------------------------------------------------------------------
subroutine average_local_IE_tau(mode,Np,x,y,z,rho,DD,I_aver,tau)
   use global; use files; use grid_params; use matrices, only:e_ab
   use orbocc, only:occnum, vecnat, ennat, fac_h;use functional_m, only:maxorb; implicit none
!..Arguments
   integer, intent(in) :: Np, mode
   real(dp), INTENT(in) :: x(Np),y(Np),z(Np), rho(Np)
   real(dp), INTENT(out) :: I_aver(Np)
   real(dp), INTENT(out) :: tau(Np)
   real(dp), INTENT(in) :: DD(lnbasis,lnbasis, lnnatorb)
!..Local
   integer :: ip, ia, ib, k, l
   real(dp) ::  BB(lnbasis,lnbasis), Bs(Np,lnbasis)
   real(dp) ::  gx(lnbasis), gy(lnbasis), gz(lnbasis)
   real(dp) :: G_lap, psum, rr
   real(dp) , external :: f_bas_Laplacian, f_bas
   real(dp) :: pr, Q, ss, rrn(Np)
   real(dp) :: phi_ia, phi_ib, g_ia_x, g_ia_y, g_ia_z, g_ib_x, g_ib_y, g_ib_z,rho_orb

   integer :: I_kinet

   integer :: i, j, nblock, imin, imax
   real(dp) :: vecekt(lnbasis,lnnatorb), e_a(lnbasis), vecmcscf(lnbasis,lnnatorb), den
   real(dp), allocatable :: H_lang(:,:), delta_mat(:,:), en_lang(:), vec_lang(:, :)
   character(22) :: dumm
   logical :: ffile

   I_aver=0._dp
   tau =0._dp

!..constructs the average local ionization potential on the grid
   do ip=1,Np
      do k=1,nbasis
         Bs(ip,k)= f_bas(k,x(ip),y(ip),z(ip))
      enddo
   enddo
     
   if (mode == 0) then !I_KS (using ennat)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip, ia, l, k, psum, BB)
!$OMP DO 
      do ip=1,Np
         do l=1,nbasis
            do k=1,nbasis
               BB(k,l)= Bs(ip,k)*Bs(ip,l)
            enddo
         enddo
         do ia = 1, maxorb  !only for closed shell implemented
            psum=0._dp
            do l=1, nbasis
               do k=1, nbasis
                  psum = psum + DD(k,l,ia)*BB(k,l)
               enddo
            enddo
            psum = psum*ennat(ia)
            I_aver(ip) = I_aver(ip) + psum
         enddo
         I_aver(ip)=2._dp*I_aver(ip)/rho(ip)
      enddo
!$OMP END DO
!$OMP END PARALLEL

   elseif (mode == 1) then !I_WF (using non diagonal lagrange multipliers)
!.... Diagonalize lagrangian matrix
      allocate ( H_lang(maxorb,maxorb), delta_mat(maxorb,maxorb), en_lang(maxorb), vec_lang(maxorb,maxorb) )
      do ia=1,maxorb
         do ib=1,ia
            delta_mat(ia,ib)=0._dp
!           den=max(sqrt(0.5d0*(occnum(ia,1)*occnum(ib,1)+occnum(ia,2)*occnum(ib,2))),zero)
            H_lang(ia,ib)=0.5_dp*(e_ab(ia,ib)+e_ab(ib,ia))     !/den
            H_lang(ib,ia)=H_lang(ia,ib)
         enddo
         delta_mat(ia,ia)=1._dp
      enddo
      call diagon_lapack(maxorb,H_lang,delta_mat,en_lang,vec_lang)
!     print*,'--------------------------------------------------------------'
!     print*,'Eigenvalues of Lagrangian:'
!     print*,en_lang
!     print*,'--------------------------------------------------------------'

!.....Convert vec_lang to vecekt
      do ia=1,maxorb
         do k=1,nbasis
            vecekt(k,ia)=0._dp
            do ib=1,maxorb
               vecekt(k,ia)=vecekt(k,ia) + vec_lang(ib,ia)*vecnat(k,ib)
            enddo
         enddo
      enddo
   
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip, ia, ib, l, k, psum, BB)
!$OMP DO 
      do ip=1,Np
         do l=1,nbasis
            do k=1,nbasis
               BB(k,l)= Bs(ip,k)*Bs(ip,l)
            enddo
         enddo
         do ia=1, maxorb
            psum=0._dp
            do l=1, nbasis
               do k=1, nbasis
                  psum = psum + vecekt(k,ia)*vecekt(l,ia)*BB(k,l)
               enddo
            enddo
            psum = psum*en_lang(ia)
            I_aver(ip) = I_aver(ip) + psum
         enddo
         I_aver(ip)=2._dp*I_aver(ip)/rho(ip)
      enddo
!$OMP END DO
!$OMP END PARALLEL

   else !(mode ne 1,2)
      stop 'model_OEP: Wrong mode'
   endif !mode = 0

   I_kinet = 2 !recommended
   if ( I_kinet == 1 ) then !I_kinet == 1: Never worked properly

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip,ia, l, k, G_lap, BB, psum)
!$OMP DO 
      do ip=1,Np
         do k=1,nbasis
            G_lap = f_bas_Laplacian(k,x(ip),y(ip),z(ip))
            do l=1,nbasis
               BB(k,l)= Bs(ip,l)*G_lap
            enddo
         enddo
         do ia = 1 , maxorb  !only for closed shell implemented
            psum=0._dp
            do k=1,nbasis
               do l=1, nbasis
                  psum = psum+BB(k,l)*DD(k,l,ia)
               enddo
            enddo
            psum=psum*occnum(ia,3)
            tau(ip) = tau(ip) + psum
         enddo
         tau(ip)=-0.5_dp*tau(ip)/rho(ip)
      enddo
!$OMP END DO
!$OMP END PARALLEL

   elseif ( I_kinet == 2 ) then !use grad**2 : tested and working

!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(ip,ia, l, k, psum, gx, gy, gz, BB, rr)
!!$OMP DO REDUCTION(+:tau)
      do ip=1,Np
         do k=1,nbasis
            gx(k)=0._dp; gy(k)=0._dp; gz(k)=0._dp
            call gradient_f_bas(k, x(ip), y(ip), z(ip), gx(k), gy(k), gz(k))
         enddo
         do k=1,nbasis
            do l=1,k
               BB(k,l)=gx(k)*gx(l)+gy(k)*gy(l)+gz(k)*gz(l)
               BB(l,k)=BB(k,l)
            enddo
         enddo

         do ia = 1, maxorb  !only for closed shell implemented
            psum = 0._dp
            do k=1,nbasis
               do l=1, nbasis
                  psum = psum + DD(k,l,ia)*BB(k,l)
               enddo
            enddo
            psum = 0.5_dp*fac_h(ia)*psum
            tau(ip) = tau(ip) + psum
         enddo
         tau(ip)=tau(ip)/rho(ip)
      enddo
!!$OMP END DO
!!$OMP END PARALLEL

   else ! Use t_P for modified RKS method : To be corrected! Never worked properly
!!$OMP PARALLEL DEFAULT(SHARED) &
!!$OMP PRIVATE(ip,ia, ib, k, psum, gx, gy, gz, Q, phi_ia, phi_ib, g_ia_x, g_ia_y, g_ia_z, g_ib_x, g_ib_y, g_ib_z)
      do ip=1,Np
         do k=1,nbasis
            gx(k)=0._dp; gy(k)=0._dp; gz(k)=0._dp
            call gradient_f_bas(k, x(ip), y(ip), z(ip), gx(k), gy(k), gz(k))
         enddo

         psum=0._dp
         do ia=2,maxorb
            phi_ia=0._dp; g_ia_x=0._dp; g_ia_y=0._dp; g_ia_z=0._dp
            do k=1,nbasis
               phi_ia=phi_ia+vecnat(k,ia)*Bs(ip,k)
               g_ia_x=g_ia_x+vecnat(k,ia)*gx(k)
               g_ia_y=g_ia_y+vecnat(k,ia)*gy(k)
               g_ia_z=g_ia_z+vecnat(k,ia)*gz(k)
            enddo
            
            do ib=1,ia-1
               phi_ib=0._dp; g_ib_x=0._dp; g_ib_y=0._dp; g_ib_z=0._dp
               do k=1,nbasis
                  phi_ib=phi_ib+vecnat(k,ib)*Bs(ip,k)
                  g_ib_x=g_ib_x+vecnat(k,ib)*gx(k)
                  g_ib_y=g_ib_y+vecnat(k,ib)*gy(k)
                  g_ib_z=g_ib_z+vecnat(k,ib)*gz(k)
               enddo

               Q=(phi_ia*g_ib_x-phi_ib*g_ia_x)**2 &
                +(phi_ia*g_ib_y-phi_ib*g_ia_y)**2 &
                +(phi_ia*g_ib_z-phi_ib*g_ia_z)**2 
               
               psum=psum+0.25_dp*fac_h(ia)*fac_h(ib)*Q
            enddo
         enddo
         tau(ip)=0.5_dp*psum/rho(ip)
      enddo
!!$OMP END DO
!!$OMP END PARALLEL

   endif ! use_laplacian

   if ( allocated(H_lang) ) deallocate ( H_lang, delta_mat, en_lang, vec_lang )
end subroutine average_local_IE_tau

subroutine Calc_H0(DD,H_0) 
! Calculates the local Hamiltonian Kinetic/External/Hartree 

!..Global
   use global; use integrals; use matrices 
   use orbocc, ONLY: fac_c, fac_h
   use functional_m, ONLY: maxorb
   implicit none

!..Arguments
   real(dp), intent(in) :: DD(lnbasis,lnbasis,lnnatorb)
   real(dp), intent(out) :: H_0(lnbasis,lnbasis)
      
!..Local variables
   integer :: m, n, l, s
   integer(8) :: ind_p
   integer :: ia, ib, nrec, nchunk, nrec_last, ichunk, irec
   logical :: inv_pairs, read_int
   real(dp) :: two_int, wmn, wls

 
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

   Coul=0._dp

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
            wmn = 1._dp
         else
            wmn = 2._dp
         endif
         if(l==s) then
            wls = 1._dp
         else
            wls = 2._dp
         endif

!.....The coulomb and exchange should be added to different F elements:
!.....Coulomb
         inv_pairs = (l/=m.or.n/=s)
         Coul(m,n,ib) = Coul(m,n,ib) + DD(l,s,ib) * wls * two_int

         if(inv_pairs) Coul(l,s,ib) = Coul(l,s,ib) + DD(m,n,ib) * wmn * two_int

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
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   do m=1,nbasis
      do n=1,m
         H_0(m,n)=Hcore(m,n)
         do ib=1,maxorb
            H_0(m,n)= H_0(m,n) + 2._dp*Coul(m,n,ib) 
         enddo
         H_0(n,m)=H_0(m,n)
      enddo
   enddo

   return

 140 stop 'Staroverov:sum_intg_sta: error reading 2-e intergral file'

end subroutine Calc_H0
!..................................................................................

subroutine read_vec_e()
!..Global
   use global; use files; use orbocc; implicit none

!..Local
   integer :: limit=10000000, linear_ind, nblock, i, j, ia, ib, iamin, iamax
   logical :: linear_depend
   character(20) :: line
   character(45) :: dumm
   character(30) :: dumm4

   print*,'Reading orbilats from ',gam_out_pot
   open(unit=11, file=gam_out_pot, status='old')
   !..Check for Linearly dependent orbitals (less orbitals printed)
   rewind(11)
   linear_depend=.false.
   do i=1, limit
      read(11,'(a20)',end=135) line
        if(line == ' NUMBER OF LINEARLY ') then
           linear_depend=.true.
           exit
        endif
    enddo
 135  continue
    rewind(11)
    if (linear_depend) then
       do i=1,limit
          read(11,'(a20)') line
          if(line == ' NUMBER OF LINEARLY ') exit
       enddo
       read(11,'(a45,i5)')dumm, linear_ind
    end if ! end if linar_depend=true
    rewind(11)

!..linear_ind is the number of linearly independent orbitals to be read
   if (.not.linear_depend)  linear_ind = nbasis
 
   do i=1,limit
      read(11,'(a22)',end=55) dumm4
      if(dumm4 == '          EIGENVECTORS') then
         read(11,*)
         goto 56
      endif
   enddo
55 stop  'Model_OEP:read_vec_e: string:          EIGENVECTORS not found in log file'

56 continue

  nblock = linear_ind/5
     do ib=0,nblock
        read(11,*)
        read(11,*)
        read(11,*)
        read(11,'(a22)')dumm4
        do j=1,nbasis
           iamin= ib*5 +1
           iamax= min(iamin+4,linear_ind)
           read(11,'(15x,5e21.12)') (vecnat(j,i),i=iamin,iamax)
        enddo
     enddo

   rewind(11)
   
   do i=1,limit
      read(11,'(a22)') dumm4
      if(dumm4 == '          EIGENVECTORS') goto 156
   enddo
   stop 'gamess:string:          EIGENVECTORS not found in log file'

156 continue
   read(11,*)

   iamin=1
   iamax=min(5,linear_ind)
   do ib=0,linear_ind/5
      read(11,*)
      read(11,*)
      read(11,*) (ennat(ia),ia=iamin,iamax)
      read(11,*)
      do j=1,nbasis
         read(11,*) 
      enddo
      iamin=iamin+5
      iamax=iamax+5
      iamax=min(iamax,linear_ind)
   enddo

   close(11)
end subroutine read_vec_e

!--------------------------------------------------------------------------------------------------
! The improved DIIS extrapolation method: P. Pulay, J. Comput. Chem. 3, 556 (1982).
! icall (in/out) : how many times DIIS_extrapolation has been called before/ number of present call
! Ndiis (in) : number of previous runs to be used / dimension of B matrix
! H (in): the hamiltonian to be extrapolated
! H_extrap (out): the extrapolated Hamiltonian
! emax : the maximum error of the e vector
subroutine DIIS_extrapolation(icall, Ndiis, H, H_extrap, emax)
!..Global
   use global; use orbocc; use matrices, only: ovlap; use functional_m, only: maxorb; implicit none

!..Arguments
   integer, intent(in) :: Ndiis
   integer :: icall
   real(dp), intent(in) :: H(lnbasis, lnbasis)
   real(dp), intent(out) :: H_extrap(lnbasis, lnbasis), emax

!..Local
   real(dp) :: Dm(lnbasis,lnbasis), abe
   real(dp) :: Ds(lnbasis,lnbasis), FDS(lnbasis,lnbasis), e_t(lnbasis, lnbasis)
   real(dp), allocatable :: B(:,:), C(:)
   real(dp), allocatable, save :: e(:,:,:), Fs(:,:,:), sqov(:,:), B_o(:,:)

   integer :: i, k, l, m, n, Np, ip

   if (.not.allocated(e) ) allocate( e(lnbasis, lnbasis, Ndiis) ) 
   if (.not.allocated(Fs) ) allocate( Fs(lnbasis, lnbasis, Ndiis) ) 
   if (.not.allocated(sqov) ) allocate( sqov(lnbasis, lnbasis) )
   if (.not.allocated(B_o) ) allocate( B_o(Ndiis, Ndiis) )

   icall=icall+1
   if (icall==1) then
      call sqrt_mat(nbasis,ovlap,sqov)
   endif

!..The density matrix
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(l, k, i)
!$OMP DO
   do k=1,nbasis
      do l=1,k
         Dm(k,l)=0._dp
         do i=1,maxorb
            Dm(k,l) = Dm(k,l) + vecnat(k,i)*vecnat(l,i)
         enddo
         Dm(l,k)=Dm(k,l)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!..DS
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(l, k, m)
!$OMP DO
   do k=1,nbasis
      do l=1,nbasis
         Ds(k,l)=0._dp
         do m=1,nbasis
            Ds(k,l)=Ds(k,l)+Dm(k,m)*ovlap(m,l)
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!..FDS
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(l, k, m)
!$OMP DO
   do k=1,nbasis
      do l=1,nbasis
         FDS(k,l)=0._dp
         do m=1,nbasis
            FDS(k,l)=FDS(k,l)+H(k,m)*Ds(m,l)
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL


!..For testing
!  do i=1,nbasis
!     do j=1, nbasis
!        Ds(i,j) = 0.d0
!        do ia=1,maxorb
!           do k=1,nbasis
!              do l=1, nbasis
!                 Ds(i,j)=Ds(i,j) + ovlap(i,k)*vecnat(k,ia)*vecnat(l,ia)*ovlap(l,j)*ennat(ia)
!              enddo
!           enddo
!        enddo
!        print*,i,j,Ds(i,j),FDS(i,j)
!     enddo
!  enddo
!  stop 'a'

   Np=Ndiis
   if ( icall < Ndiis) Np=icall

   allocate( B(Np+1, Np+1), C(Np+1) )

!..Save present hamiltonian
   if( icall > Ndiis) then
      do ip=1,Np-1
         Fs(:,:,ip)=Fs(:,:,ip+1)
         e(:,:,ip)=e(:,:,ip+1)
      enddo
!     Fs = CSHIFT(Fs, 1, 3)
!     e  = CSHIFT(e, 1, 3)
   endif
   Fs(:,:,Np) = H

!..FDS-SDF
   do k=1,nbasis
      do l=1,nbasis
         e_t(k,l) = FDS(k,l)-FDS(l,k)
      enddo
   enddo      

!.. A+ e_t A (A=S^-/2)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(l, k, m)
!$OMP DO
    do k=1,nbasis
       do l=1,nbasis
          Ds(k,l)=0._dp
          do m=1,nbasis
             Ds(k,l)= Ds(k,l) + e_t(k,m)*sqov(m,l)
          enddo
       enddo
    enddo 
!$OMP END DO
!$OMP END PARALLEL

    emax=0._dp
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(l, k, m, abe)
!$OMP DO
    do k=1,nbasis
       do l=1,nbasis
          e(k,l,Np)=0._dp
          do m=1,nbasis
             e(k,l,Np)= e(k,l,Np) + sqov(m,k)*Ds(m,l)
          enddo
          abe=abs(e(k,l,Np))
!$OMP CRITICAL
          if ( abe > emax ) then
              emax=abe
          endif
!$OMP END CRITICAL 
       enddo
    enddo 
!$OMP END DO
!$OMP END PARALLEL

    print*,'DIIS emax:',emax

    do k=1,np-1
       do l=1,np-1
          if(icall > Ndiis) B_o(k,l) = B_o(k+1,l+1)
          B(k,l)=B_o(k,l)
       enddo
    enddo
   
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k, m, n)
!$OMP DO
    do k=1,np
       B_o(k,Np) = 0._dp 
       do m=1,nbasis
          do n=1,nbasis 
             B_o(k,Np) = B_o(k,Np) + e(m,n,k)*e(m,n,Np)
          enddo
       enddo
       B(k,np)=B_o(k,Np)
       B_o(Np,k)=B_o(k,Np)
       B(Np,k)=B(k,Np)
    enddo    
!$OMP END DO
!$OMP END PARALLEL
 
   do k=1,Np
      B(k,Np+1) = -1._dp
      B(Np+1,k) = -1._dp
      C(k) =0._dp
   enddo
   B(Np+1, Np+1) = 0._dp
   C(Np+1)=-1._dp

!..Singular value decomposition/ invertion
!  call invert_B(Np+1, B)

!  s=0.d0
!  do k=1,Np
!     C(k)=-B(Np+1,k)
!     print*,'C',k,c(k)
!     s=s+c(k)
!  enddo
!  print*,'s',s

!..Solve linear system B X = C 
   call lin_sqsym_solve(Np+1, B, C) 

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k, l, ip)
!$OMP DO
   do k=1,nbasis
      do l=1,k
         H_extrap(k,l)=0._dp
         do ip=1,Np
            H_extrap(k,l)= H_extrap(k,l) + C(ip)*Fs(k,l,ip)
         enddo
         H_extrap(l,k)=H_extrap(k,l)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine DIIS_extrapolation

!----------------------------------------------------
! Just invert B matrix
subroutine invert_B(N, B)
!..Global
   use params_general; implicit none

!..Arguments
   integer, intent(in) :: N 
   real(dp), intent(in) :: B(N, N)

!..Local 
   real(dp) :: B_i(N, N)
   real(dp) :: WORK(N*N)
   integer :: INFO
   integer :: IPIV(N)

!  LAPACK ROUTINES
   external DGETRF 
   external DGETRI

   B_i=B
!..DGETRF computes an LU factorization of a general M-by-N matrix A
   CALL DGETRF( N, N, B, N, IPIV, INFO )

   IF(INFO.EQ.0)THEN
!     PRINT '(" LU decomposition successful ")'
   ENDIF
   IF(INFO.LT.0)THEN
      PRINT '(" LU decomposition:  illegal value ")'
      stop 'DIIS:invert_B:DGETRF'
   ENDIF
   IF(INFO.GT.0)THEN
      WRITE(*,35)INFO,INFO
35    FORMAT( 'LU decomposition: U(',I4,',',I4,') = 0 ')
   ENDIF

!..DGETRI computes the inverse of a matrix using the LU factorization
   CALL DGETRI(N, B, N, IPIV, WORK, N*N, INFO)
   IF (info.NE.0) THEN
      stop 'DIIS:invert_B:DGETRI Matrix inversion failed!'
   ELSE
!     PRINT '(" Inverse of B Successful ")'
   ENDIF

end subroutine invert_B

!=======================================================================
subroutine  sqrt_mat(n,X,Xsq)
!..Global
   use params_general; implicit none

   integer :: n
   real(dp), intent(in) :: X(n,n)
   real(dp), intent(out) :: Xsq(n,n)

!..Local
   real(dp) :: A(n,n)
   real(dp) :: B(n,n)
   real(dp) :: vectr(n,n), eigs(n)
   integer :: mu, nu, ku, info
!..Definitions for LAPACK routine
   character(1) :: JOBZ, RANGE, UPLO
   INTEGER :: IL, ITYPE,  IU,  LDA,  LDB,  LDZ
   INTEGER LWORK, M
   real(dp) :: VL,VU
   real(dp) :: ABSTOL
   real(dp) :: WORK(10*n)
   integer IWORK(5*N)
   integer :: ifail(N)
   
   real(dp) :: DLAMCH
   external DLAMCH

   ITYPE=1; JOBZ='V'; RANGE='I'; UPLO='L'; LDA = N
   LDB = N; VL = 0._dp; VU = 0._dp; IL = 1; IU = N
   LDZ = N; LWORK=10*N

   ABSTOL = 2*DLAMCH('S')
!  ABSTOL = 1.d-30
      
   do mu=1, n   
      do nu=1, mu
         A(mu,nu) = X(mu,nu)
         B(mu,nu) = 0._dp
         B(nu,mu) = 0._dp
         A(nu,mu) = A(mu,nu) 
      enddo
      B(mu,mu)=1._dp
   enddo

   WORK=0._dp
   IWORK=0
   call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, &
                LDB,  VL,  VU, IL, IU, ABSTOL, M, eigs, vectr, &
                LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

   if (info /= 0) then
      print*,'Info:',info
      if(info < 0) print*,info, '-th argument of DSYGVX has wrong value'
      if(info > 0 .and. info <= n ) print*,info, '-th eigenvalue not converged'
      if(info > n .and. info <= 2*n ) &
         print*,info-n, '-th leading minor of ovlap not positive definite'
      stop 'sqrt_mat:Diagonalization failed!'
   endif

   do ku=1,N
      eigs(ku)=1.d0/sqrt(max(eigs(ku),zero))
   enddo

   do mu=1,N
      do nu=1,mu
         Xsq(mu,nu)=0._dp
         do ku=1,N
            Xsq(mu,nu) = Xsq(mu,nu) + vectr(mu,ku)*eigs(ku)*vectr(nu,ku)
         enddo
         Xsq(nu,mu)= Xsq(mu,nu)
      enddo
   enddo
end subroutine sqrt_mat 

subroutine lin_sqsym_solve(N, B, C)
!..Global
   use params_general; implicit none

!..Arguments
   integer, intent(in) :: N 
   real(dp), intent(in) :: B(N, N)
   real(dp), intent(inout) :: C(N,1)

!..Local 
   integer, parameter :: LWMAX=100, NRHS=1
   real(dp) :: A(N,N)
   real(dp), allocatable :: WORK(:)
   integer :: INFO
   integer :: IPIV(N), LWORK

!  LAPACK ROUTINES
   external DSYSV

   allocate( WORK(LWMAX) )

   A=B

   LWORK = -1
   CALL DSYSV( 'Lower', N, NRHS, A, N, IPIV, C, N, WORK, LWORK, INFO )
   LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

   deallocate ( work); allocate( work(lwork) )

   CALL DSYSV( 'Lower', N, NRHS, A, N, IPIV, C, N, WORK, LWORK, INFO )

   IF( INFO.GT.0 ) THEN
      WRITE(*,*)'The element of the diagonal factor '
      WRITE(*,*)'D(',INFO,',',INFO,') is zero, so that'
      WRITE(*,*)'D is singular; the solution could not be computed.'
      stop 'In DIIS_extrapolation: lin_sqsym_solve, Non invertible matrix'
   END IF
   
end subroutine lin_sqsym_solve

!-------------------------------------------------------------------------------------
subroutine construct_DD(nbasis,maxorb,vecnat,DD)

!..Global
   use params_general; implicit none

!..Arguments
   integer, intent(in) :: nbasis, maxorb
   real(dp), intent(in) :: vecnat(lnbasis, lnbasis)
   real(dp), intent(out) :: DD(lnbasis, lnbasis, lnnatorb)

!..Local
   integer :: ia, l, k

!..The density matrix
!..Redefine the density matrix
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,l,k)
!$OMP DO
   do ia=1,maxorb
      do l=1,nbasis
         do k=1,l
            DD(k,l,ia)=vecnat(k, ia)*vecnat(l, ia)
            DD(l,k,ia)=DD(k,l,ia)
         enddo
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

   return
end subroutine construct_DD

subroutine canonical_orbitals(IP_meth, e_IP)
!..Global
   use global; use matrices; use orbocc; use functional_m, only:maxorb; implicit none

!..Arguments
   integer, intent(in) :: IP_meth
   real(dp), intent(inout) :: e_IP

!..Local
   real(dp) :: vec_can(lnbasis, lnbasis), ennat_can(lnbasis)
   real(dp), allocatable :: xlam(:,:), smat(:,:), eigs1(:), Psimat(:,:), xpro(:)
   real(dp) :: den, e_shift, e_IP_0
   integer :: ia, ib, ic, ie, ih, indx(lnbasis)
   logical :: use_projection

   e_IP_0=e_IP

   call calc_Lagrange_mult_all()
!..Local array allocation
   allocate ( xlam(maxorb,maxorb), smat(maxorb,maxorb), eigs1(maxorb), Psimat(maxorb,maxorb), xpro(maxorb) ) 

   do ia=1,maxorb
      do ib=1,maxorb
         den=max(sqrt(0.5_dp*(occnum(ia,1)*occnum(ib,1)+occnum(ia,2)*occnum(ib,2))),zero)
         xlam(ia,ib)=0.5_dp*(e_ab(ia,ib)+e_ab(ib,ia)) /den
!        if(abs(xlam(ia,ib)) > 1.e2_dp)  xlam(ia,ib)=0._dp
         smat(ia,ib)=0._dp
      enddo
      smat(ia,ia)=1._dp
   enddo
!  write(22,*) ((ia,ib, xlam(ib,ia), xlam(ia,ib), ib=1,ia), ia=1,lnbasis)

   call diagon_lapack(maxorb, xlam, smat, eigs1, Psimat)
!  write(23,'(15f10.4)') ((Psimat(ia,ib), ia=1,lnbasis), ib=1,lnbasis)

   use_projection=.true.
   if ( use_projection ) then
      ih=nele(1)
      do ie=1,maxorb
         indx(ie)=1; xpro(ie)=Psimat(1,ie)**2
         do ia=1,maxorb
            if (xpro(ie) <= Psimat(ia,ie)**2) then
               xpro(ie)=Psimat(ia,ie)**2; indx(ie)=ia
            endif
         enddo
         if(indx(ie)==ih) e_IP=eigs1(ie)
      enddo
   else
      e_IP=-1.e20_dp
      do ie=1,maxorb
         if(eigs1(ie).lt.0._dp .and. eigs1(ie).gt.e_IP ) e_IP=eigs1(ie)
      enddo
   endif

! write(*,'(5f20.8)')((xlam(ia,ib),ia=1,nnatorb),ib=1,nnatorb)

   if (IP_meth == 3) then
      do ia=1,maxorb
         do ib=1,maxorb
            xlam(ia,ib)=0.5_dp*(e_ab(ia,ib)+e_ab(ib,ia))
            smat(ia,ib)=0._dp
         enddo
         smat(ia,ia)=1._dp
      enddo
      call diagon_lapack(maxorb, xlam, smat, eigs1, Psimat)
      e_shift=e_IP-e_IP_0
      eigs1=eigs1+e_shift
      do ia=1,maxorb
         do ib=1,maxorb
            e_ab(ib,ia)=0.d0
            do ic=1,maxorb
               e_ab(ib,ia) = e_ab(ib,ia) + eigs1(ic)*Psimat(ic,ib)*Psimat(ic,ia)
            enddo
         enddo
      enddo 
   endif

!  print*,'Extended Koopman`s theorem eigenvalues: (a.u.) Starov'
!  print*,'  Nu    Eigenv.    Indx       Max Pro'
!  write(6,"(i4,':',f20.8,i5,f10.4)"),(ia, ha_ev*eigs1(ia) ,indx(ia), xpro(ia),ia=1,maxorb)
!  print*,'------------------------------------------------------'
   deallocate ( xlam, eigs1, Psimat ) 

end subroutine canonical_orbitals

!-----------------------------------------------------------------------------------------
subroutine test_Xintegral2
!..Global
   use global; use basis_set; use matrices; implicit none

!..Local
   real(dp) :: test_V(lnbasis, lnbasis), xint, x, y, z
   integer :: k,l,iat

   test_V=0._dp

   do iat=1,natom
      x=xcoo(iat); y=ycoo(iat); z=zcoo(iat)
      do k=1,nbasis
         do l=1,k
            call Xintegral2(k,l,x,y,z,xint)
            test_V(k,l)=test_V(k,l)-xint*charge(iat)
         enddo
      enddo
   enddo

      do k=1,nbasis
         do l=1,k
!            if ( abs(test_V(k,l)-(hcore(k,l)-kin(k,l))) > 1e-8_dp ) then
               print*,k,l,abs(test_V(k,l)-(hcore(k,l)-kin(k,l)))
               print*,k,l,test_V(k,l),hcore(k,l)-kin(k,l)
               print*,' '
!            endif
         enddo
      enddo

end subroutine test_Xintegral2

!-----------------------------------------------------------------------------------------
subroutine Xintegral2(k,l,x,y,z,xint)
!  This routine calculates the integral I(r) = \int d^3r_1 \frac{chi_k(r_1) chi_l(r_1)}{|r-r_1|)
!  where r and r_1 are vectors, r=(x,y,z) and chi_k and chi_l are basis functions. 
!  The matrix element of the Coulomb ionic potential on k and l basis functions is then given 
!  by V_{kl}= \sum_a I(R_a)*Z_a where index a counts the atoms, R_a=(X_a, Y_a, Z_a) is the
!  atomic position and Z_a is the atomic number.
 
!..Global
   use global; use basis_set; implicit none

!..Arguments
   real(dp), intent(in) :: x,y,z
   integer, intent(in) :: k,l
   real(dp), intent(out) :: xint

!..Local
   integer :: k_at, l_at, k_ga, l_ga, k_g, l_g, n_xk, n_yk, n_zk, n_xl, n_yl, n_zl
   integer :: i_xk, i_yk, i_zk, i_xl, i_yl, i_zl, ne_x, ne_y, ne_z
   real(dp) :: P_x, P_y, P_z, expo_c, xd, yd, zd, cc_x(0:10,0:10), cc_y(0:10,0:10), cc_z(0:10,0:10)
   real(dp) ::  r_ab2, expo_k, expo_l, exp_1, coef_1, coef_c, xx
   real(dp), external :: X_integr, F_N_ov_m
     
   k_at=iat_bas(k)
   l_at=iat_bas(l)

   k_ga=nga_bas(k)
   l_ga=nga_bas(l)

   r_ab2=R_AB(k_at,l_at)**2

   xint=0._dp

   do k_g=1,k_ga
      do l_g=1,l_ga
         expo_k = expo_bas(k_g,k)
         expo_l=expo_bas(l_g,l)
         expo_c=expo_k+expo_l
         P_x=(expo_k*xcoo(k_at)+expo_l*xcoo(l_at))/expo_c
         P_y=(expo_k*ycoo(k_at)+expo_l*ycoo(l_at))/expo_c
         P_z=(expo_k*zcoo(k_at)+expo_l*zcoo(l_at))/expo_c
         xd = P_x - x 
         yd = P_y - y 
         zd = P_z - z 
         n_xk=Int(expo_x(k)+0.001_dp)
         n_yk=Int(expo_y(k)+0.001_dp)
         n_zk=Int(expo_z(k)+0.001_dp)
         n_xl=Int(expo_x(l)+0.001_dp)
         n_yl=Int(expo_y(l)+0.001_dp)
         n_zl=Int(expo_z(l)+0.001_dp)

         exp_1=expo_k*expo_l*r_ab2/expo_c
         coef_1 = exp(-exp_1)*fact_norm(k_g,k)*fact_norm(l_g,l)*coef_bas(k_g,k)*coef_bas(l_g,l)

         do i_xk=0,n_xk
         do i_xl=0,n_xl
            cc_x(i_xk, i_xl) = xnovm(n_xk,i_xk)*xnovm(n_xl,i_xl)* &
!           cc_x(i_xk, i_xl) = F_N_ov_m(n_xk,i_xk)*F_N_ov_m(n_xl,i_xl)* &
                               (P_x-xcoo(k_at))**(n_xk-i_xk)* &
                               (P_x-xcoo(l_at))**(n_xl-i_xl)
         enddo
         enddo
         
         do i_yk=0,n_yk
         do i_yl=0,n_yl
            cc_y(i_yk, i_yl) = xnovm(n_yk,i_yk)*xnovm(n_yl,i_yl)* &
!           cc_y(i_yk, i_yl) = F_N_ov_m(n_yk,i_yk)*F_N_ov_m(n_yl,i_yl)* &
                               (P_y-ycoo(k_at))**(n_yk-i_yk)* &
                               (P_y-ycoo(l_at))**(n_yl-i_yl)
         enddo
         enddo
         
         do i_zk=0,n_zk
         do i_zl=0,n_zl
            cc_z(i_zk, i_zl) = xnovm(n_zk,i_zk)*xnovm(n_zl,i_zl)* &
!           cc_z(i_zk, i_zl) = F_N_ov_m(n_zk,i_zk)*F_N_ov_m(n_zl,i_zl)* &
                               (P_z-zcoo(k_at))**(n_zk-i_zk)* &
                               (P_z-zcoo(l_at))**(n_zl-i_zl)
         enddo
         enddo
         
         do i_xk=0,n_xk
         do i_xl=0,n_xl
         do i_yk=0,n_yk
         do i_yl=0,n_yl
         do i_zk=0,n_zk
         do i_zl=0,n_zl
            ne_x=i_xk+i_xl
            ne_y=i_yk+i_yl
            ne_z=i_zk+i_zl
            coef_c = coef_1*cc_x(i_xk, i_xl)*cc_y(i_yk, i_yl)*cc_z(i_zk, i_zl)
            xx=X_integr(ne_x, ne_y, ne_z, expo_c, xd, yd, zd)
            xint=xint + coef_c*xx
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      enddo
   enddo

end subroutine Xintegral2

subroutine Kinetic_energy(E_kin)
   use global; use matrices; use orbocc; use functional_m, only:maxorb; implicit none

!..Arguments
   real(dp), intent(out) :: E_kin

!..Local variables
   real(dp) :: ss
   integer :: ia, l, m

   do l=2,nbasis
      do m=1,l-1
         kin(m,l)=kin(l,m)
      enddo
   enddo

   print*,'maxorb',maxorb
   E_kin=0._dp
   do ia=1,maxorb
      ss=0._dp
      do l=1,nbasis
         do m=1,nbasis
            ss=ss+vecnat(l,ia)*kin(l,m)*vecnat(m,ia)
         enddo
      enddo
      E_kin=E_kin+occnum(ia,3)*ss
   enddo
   
end subroutine Kinetic_energy

!-----------------------------------------------------------------------------------------
!subroutine total_energy_Veff(rho, pot_grid, Total_Energy)
!   use global; use matrices; use orbocc; use grid_params
!   use functional_m, only:maxorb; implicit none
!
!!..Arguments
!   real(dp), intent(in) :: rho(ngrid), pot_grid(ngrid)
!   real(dp), intent(out) :: Total_Energy
!
!!..Local
!   real(dp) :: Bare_energy, Coul_energy, Exch_energy
!   integer :: ia, ib, ig
!
!   call construct_f(0)
!
!   Bare_energy = 0.d0
!   do ia=1, nnatorb
!      Bare_energy = Bare_energy + HcoreNO(ia,ia)*fac_h(ia)
!   enddo
!
!   Coul_energy = 0.d0
!   do ia=1, maxorb
!      do ib=1, maxorb
!         Coul_energy = Coul_energy + fac_c(ia,ib)*CoulNO(ia,ib)
!      enddo
!   enddo
!
!   Exch_energy = 0.d0
!   do ig=1,ngrid
!      Exch_energy=Exch_energy + w_grid(ig)*rho(ig)*pot_grid(ig)
!   enddo
!
!   Total_energy= Bare_energy + Coul_energy + Exch_energy
!   print*,Bare_energy,Coul_energy, Exch_energy
!end subroutine total_energy_Veff
