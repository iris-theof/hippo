!--------------------------------------------------------------------------
! Contains routines for the calculation of the total energy 
! or one electron energies 
!
!   Created by N. N. Lathiotakis and I.Theophilou
!
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine total_energy_MO(ii)
!--------------------------------------------------------------------------
!
! This subroutine calculates the Electronic energy components and sums
! them together with repulsive energy to fine the total energy.
! 
! the formula (3.181) of Attila Szabo, N. S. Ostlund: Modern Quantum
! Chemistry is used
!
! Input: 
!       hcore: 1-e hamiltonian matrix
!       Coul : the Coulomb Interaction matrix(ces)
!       Exch : the Exchange Interaction matrix(ces)
!      vecnat: the 1-e orbitals
!          ii: for BDF functional ii=0 : RDMFT functional for the total energy
!                                 ii=1 : DFT functional for the total energy
!
! Output:
!     Bare_energy: Contribution of 1-e terms to the total energy 
!     Coul_energy: Contribution of e-e Coulomb interaction to the total
!                  energy
!     Exch_energy: Contribution of electron Exchange to the total energy
!
!--------------------------------------------------------------------------

!..Global
   use global; use functional_m; use matrices; use orbocc; use energies
   use grid_params; use DFT_par
   implicit none

!..Arguments
   integer :: ii
!..For functional == 'DBF' ii=0 use the RDM functional ii=1 use KS

!..Local variables
   integer :: ia, ib, ig
   real(dp) :: E_non_JK
   real(dp), external :: TS_nlnn

   Bare_energy = 0.d0
   do ia=1, nnatorb
      Bare_energy = Bare_energy + HcoreNO(ia,ia)*fac_h(ia)
   enddo

   Coul_energy = 0.d0
   Exch_energy = 0.d0
   do ia=1, maxorb
      do ib=1, maxorb
         Coul_energy = Coul_energy + fac_c(ia,ib)*CoulNO(ia,ib)
         if ( functional /= 'DBF' ) then
!.... TOMS energy change to add EX energy to HYB 
            if ( l_hybc ) then 
               Exch_energy = Exch_energy - ExchNO(ia,ib)
            else
               Exch_energy = Exch_energy + fac_e(ia,ib)*ExchNO(ia,ib)
            endif
         else
            if ( ii == 0 ) Exch_energy = Exch_energy + fac_e(ia,ib)*ExchNO(ia,ib)
         endif
      enddo
   enddo

   if (Functional == 'DFT' .or. (Functional == 'DBF' .and. ii /= 0) ) then
      ExcDFT_ener=0.d0; EcDFT_ener=0.d0
      do ig=1,ngrid
         ExcDFT_ener=ExcDFT_ener+ w_grid(ig)*chdens_gr(ig,3)*Exc_DFT(ig)
         EcDFT_ener=EcDFT_ener+ w_grid(ig)*chdens_gr(ig,3)*Ec_DFT(ig)
      enddo
      ExcDFT_ener=ExcDFT_ener+EcDFT_ener
!     print*,'******** DFT XC Energies: ********'
!     print*,'  HF Exchange: ', Exch_energy
!....TOMS energy change to add EX energy to HYB
      if ( l_hybc ) then 
         Exch_energy = (1-ahyb) * (Exch_energy+EcDFT_ener) + ahyb* ExcDFT_ener !new hyb
!         Exch_energy = (1-ahyb) * (Exch_energy) + ahyb* ExcDFT_ener ! original Hyb
      else 
         Exch_energy = Exch_energy + ExcDFT_ener
      endif
!     print*,' DFT Exc/Corr: ', ExcDFT_ener
!     if (abs(EcDFT_ener) > small) print*,'  Energy from DFT Corr. Functional: ', EcDFT_ener
!     print*,'**********************************'
   elseif ( Functional == 'HDR' ) then
      ExcDFT_ener=0.d0; EcDFT_ener=0.d0
      do ig=1,ngrid
         ExcDFT_ener=ExcDFT_ener+ w_grid(ig)*chdens_gr(ig,3)*Exc_DFT(ig)*xmix_DF
         EcDFT_ener=EcDFT_ener+ w_grid(ig)*chdens_gr(ig,3)*Ec_DFT(ig)*xmix_DF
      enddo
      ExcDFT_ener=ExcDFT_ener+EcDFT_ener
      Exch_energy = Exch_energy + ExcDFT_ener
   endif

   Totel_energy = Bare_energy + Coul_energy + Exch_energy

   if ( l_non_JK ) then
      E_non_JK=f_non_JK*q_non_JK
!     print*,'Energy contribution from non_JK terms:',f_non_JK, q_non_JK
      Totel_energy = Totel_energy + E_non_JK
   endif
   

   TS_ener = 0.d0
   if(temp_dep) TS_ener = TS_nlnn()

   Totel_energy= Totel_energy + TS_ener
      
   Tot_energy = Rep_energy + Totel_energy

end subroutine total_energy_MO

!--------------------------------------------------------------------------
! The nuclear repulsive Energy:
   subroutine repulsive_energy()

!.....Global
   use global
   use energies
   implicit none

!.....Local variables
   integer :: iatom, jatom
   real(dp) :: xdist, ydist, zdist, rdist

   Rep_energy = 0.d0
   do iatom = 2, natom
      do jatom = 1, iatom - 1
         xdist = xcoo(iatom) - xcoo(jatom)
         ydist = ycoo(iatom) - ycoo(jatom)
         zdist = zcoo(iatom) - zcoo(jatom)
         rdist = sqrt(xdist*xdist + ydist*ydist + zdist*zdist)
         Rep_energy = Rep_energy + charge(iatom)*charge(jatom)/rdist
      enddo
   enddo
end subroutine repulsive_energy

!--------------------------------------------------------------------------
subroutine orb_energies(xmu_energies)
! Find Orbital energies as expectation values of nat orbs on
! Fock operators:

!..Global
   use global; use functional_m; use matrices; use orbocc
   implicit none

!..Argument
   real(dp), intent(out) :: xmu_energies(lnnatorb)

!..Local variables
   integer :: ia,ib

   do ia= 1, nnatorb
      xmu_energies(ia) = HcoreNO(ia,ia)*fac_h(ia)*0.5d0
      do ib = 1, nnatorb
         xmu_energies(ia) = xmu_energies(ia) + fac_e(ia,ib)*ExchNO(ia,ib)
         xmu_energies(ia) = xmu_energies(ia) + fac_c(ia,ib)*CoulNO(ia,ib)
      enddo
   enddo

!   print*,'E_aa (Lagrange multiplier diagonal elements:'
!   print*,(xmu_energies(ia),ia=1,nnatorb)
      
end subroutine orb_energies

!-------------------------------------------------------------------------
subroutine energy_spectrum(meth_sp, xmu)
! This routine calculates the energy spectrum using 3 proposed methods:
! meth_sp=1 : Diagonal of Lagrange Multipliers
! meth_sp=2 : Pernal-Cioslowski Extended-Koopmans theorem (CPL 412, 71)
! meth_sp=3 : Sharma et al Method (arXiv cond-mat: 0912.1118v1)

!..Global
   use global; use functional_m; use matrices; use orbocc; use energies; use grid_params; use DFT_par
   implicit none

!..Arguments
   real(dp),intent(in) :: xmu(2)
   integer :: meth_sp

!..Local variables
   integer :: ia, ib, ie, indx(lnnatorb), leig, iho, ig
   real(dp), allocatable :: occnum_temp(:,:)
   real(dp), allocatable :: DE_Dn(:)
   real(dp), allocatable :: xlam(:,:), eigs1(:)
   real(dp), allocatable :: Psimat(:,:), smat(:,:)
   real(dp) :: Evalues(lnnatorb,2), den
   real(dp) :: xpro(lnnatorb)
   real(dp) :: e_2, e_3, e_4, e_5, rho(ngrid), dens_homo_gr, eho
   character(4) :: orb_t(lnnatorb)

   orb_t='    '
   orb_t(nele(1))='HOMO'
   orb_t(nele(1)+1)='LUMO'
   print*,'------------------------------------------------------'
   print*,'------------------------------------------------------'
   print*,'||                ENERGY SPECTRUM                   ||'
   print*,'------------------------------------------------------'
   print*,'------------------------------------------------------'

   if((Functional == 'RHF') .or. (Functional == 'DFT') .or. (Functional == 'DBF') .or. (Functional == 'RHA') .or. do_oep ) then
      print*,'------------------------------------------------------'
      if (.not. do_oep) then
         select case (functional)
            case('RHF') 
               print*,'RHF eigenvalues:'
            case('RHA') 
               print*,'RHA eigenvalues:'
            case('DFT')
               print*,'DFT eigenvalues:'
            case('DBF')
               print*,'DFT eigenvalues:'
            case default
         end select
      else
         print*,'OEP effective Hamiltonian eigenvalues:'
      endif
      print*,'------------------------------------------------------'
      write(6,'(a4,"  ",e20.8," a.u.   ", e20.8," eV")') &
           (orb_t(ia),ennat(ia),ha_ev*ennat(ia),ia=1,nnatorb)
      print*,'------------------------------------------------------'
      if(Functional=='DFT' .and. (.not.do_oep)) then
      if( &
         xc_f90_info_family(xc_info) == XC_FAMILY_LDA .and. &
         xc_f90_info_family(c_info) == XC_FAMILY_LDA ) & 
      then
         iho=nele(1)
         
         e_2=0.5d0*CoulNO(iho,iho)
         
         e_3=0.d0; e_4=0.d0
         do ig=1,ngrid
            dens_homo_gr = occnum(ia,1)*vec_nat_grid(ig,iho)**2 
            rho(ig)=chdens_gr(ig,3)-dens_homo_gr
            e_3=e_3+w_grid(ig)*dens_homo_gr*(Vlxc(1,ig)+Vlc(1,ig))
            e_4=e_4+w_grid(ig)*chdens_gr(ig,3)*(Exc_DFT(ig)+Ec_DFT(ig))
         enddo

         call xc_f90_lda_exc_vxc(xc_func, ngrid, rho(1), Exc_DFT(1), Vlxc(1,1))
         call xc_f90_lda_exc_vxc(c_func, ngrid, rho(1), Ec_DFT(1), Vlc(1,1))
         e_5=0.d0
         do ig=1,ngrid
            dens_homo_gr = occnum(ia,1)*vec_nat_grid(ig,iho)**2
            e_5=e_5+w_grid(ig)*rho(ig)*(Exc_DFT(ig)+Ec_DFT(ig))
         enddo
         print*,'e_2',e_2
         print*,'e_3',e_3
         print*,'e_4',e_4
         print*,'e_5',e_5
         print*,'- e_3 + e_4 - e_5:',- e_3 + e_4 - e_5
         eho=ennat(iho) - e_2 - e_3 + e_4 - e_5

         write(6,'("Nikitas corrected (-IP): ",e20.8," a.u.   ", e20.8," eV")') eho, ha_ev*eho
!        print*,'Nikitas corrected LUMO: '
!        write(6,'(e20.8," a.u.   ", e20.8," eV")') &
!        (ennat(iho+1)+0.5d0*CoulNO(iho+1,iho+1)),ha_ev*(ennat(iho+1)+0.5d0*CoulNO(iho+1,iho+1))
         print*,'------------------------------------------------------'
      endif; endif
      return
   endif

!..1. Find the Orbital energies from the Lagrange Multiplier Matrix

   if (meth_sp == 1) then
      call calc_Lagrange_mult_all()
      do ia=1,nnatorb
         ennat(ia)=e_ab(ia,ia)
      enddo
      print*,'From the diagonals of Lagrange multipliers:'
      write(6,"(i4,':',f15.8)")(ia,ha_ev*ennat(ia),ia=1,nnatorb)
      print*,'------------------------------------------------------'

   elseif(meth_sp == 2) then
!..2. Pernal-Cioslowski Extended-Koopmans theorem (CPL 412, 71):
      leig=nnatorb !max(ibond(1),ibond(2))+4

      call calc_Lagrange_mult_all()
!..Local array allocation
      allocate ( xlam(leig,leig), smat(leig,leig), eigs1(leig), Psimat(leig,leig) ) 

      do ia=1,leig
         do ib=1,leig
            den=max(sqrt(0.5d0*(occnum(ia,1)*occnum(ib,1)+occnum(ia,2)*occnum(ib,2))),zero)
            xlam(ia,ib)=0.5d0*(e_ab(ia,ib)+e_ab(ib,ia)) /den
!           if(abs(xlam(ia,ib)) > 1.d2)  xlam(ia,ib)=0.d0
            smat(ia,ib)=0.d0
         enddo
         smat(ia,ia)=1.d0
      enddo
      write(22,*) ((ia,ib, xlam(ib,ia), xlam(ia,ib), ib=1,ia), ia=1,leig)

      call diagon_lapack(leig, xlam, smat, eigs1, Psimat)
      write(23,'(15f10.4)') ((Psimat(ia,ib), ia=1,leig), ib=1,leig)

      do ie=1,leig
         indx(ie)=1; xpro(ie)=Psimat(1,ie)**2
         do ia=1,leig
            if (xpro(ie) <= Psimat(ia,ie)**2) then
               xpro(ie)=Psimat(ia,ie)**2; indx(ie)=ia
            endif
         enddo
      enddo

!  write(*,'(5f20.8)')((xlam(ia,ib),ia=1,nnatorb),ib=1,nnatorb)

      print*,'Extended Koopman`s theorem eigenvalues: (eV)'
      print*,'  Nu    Eigenv.    Indx       Max Pro'
      write(6,"(i4,':',f20.8,i5,f10.4)") (ia, ha_ev*eigs1(ia) ,indx(ia), xpro(ia),ia=1,leig)
      print*,'------------------------------------------------------'
      deallocate ( xlam, eigs1, Psimat ) 

   elseif(meth_sp==3) then
!..3. Sharma et al Method (arXiv cond-mat: 0912.1118v1):
      allocate ( occnum_temp(lnnatorb, 3),DE_Dn(lnnatorb) )
      occnum_temp=occnum

      do ia=1,nnatorb
         occnum(ia,1) = 0.5d0
         occnum(ia,3) = occnum(ia,1) + occnum(ia,2) 
         call Func_der_n(xmu(1), occnum, DE_Dn, 1)
         Evalues(ia,1)=DE_Dn(ia)+xmu(1) !Func_der_n finds DE_Dn-xmu
  
         occnum(ia,1) = 1.d0
         occnum(ia,3) = occnum(ia,1) + occnum(ia,2)
         call calc_occ_factors(occnum, fac_h, fac_c, fac_e)
         call construct_f(0)
         call total_energy_MO(1)
         Evalues(ia,2) = Tot_energy

         occnum(ia,1) = 0.d0
         occnum(ia,3) = occnum(ia,1) + occnum(ia,2)
         call calc_occ_factors(occnum, fac_h, fac_c, fac_e)
         call construct_f(0)
         call total_energy_MO(1)
         Evalues(ia,2) = Evalues(ia,2)-Tot_energy

         occnum(ia,1) = occnum_temp(ia,1)
         occnum(ia,3) = occnum_temp(ia,3)
      enddo
      call calc_occ_factors(occnum, fac_h, fac_c, fac_e)
      

!     do ia=1,nnatorb-1
!        emin=Evalues(ia,1); imin=ia
!        do ib=ia+1, nnatorb
!           if( emin > Evalues(ib,1) ) then
!              emin=Evalues(ib,1)
!              imin=ib
!           endif
!        enddo
!        emin = Evalues(ia,1)
!        Evalues(ia,1) = Evalues(imin,1)
!        Evalues(imin,1) = emin

!        emin = Evalues(ia,2)
!        Evalues(ia,2) = Evalues(imin,2)
!        Evalues(imin,2) = emin
!     enddo

      print*,'------------------------------------------------------'
      print*,'Sharma et al Method (arXiv cond-mat: 0912.1118v1)'
      print*,'ENERGIES (eV):'
      print*,'        dE/dn |(1/2)       DE (1->0)'
      write(6,"(i4,':',2f16.7)") (ia,ha_ev*Evalues(ia,1),ha_ev*Evalues(ia,2), ia=1,nnatorb)
      print*,'------------------------------------------------------'
      deallocate ( occnum_temp, DE_Dn)

   else
      stop 'energy_spectrum: meth_sp out of range, no such spectrum method.'
   endif

end subroutine energy_spectrum
!--------------------------------------------------------------------------

!subroutine total_energy_AO() !This is obsolete! Replaced by total_energy_MO
!--------------------------------------------------------------------------
!
! This subroutine calculates the Electronic energy components and sums
! them together with repulsive energy to fine the total energy.
! THE EXPANSION IN BASIS FUNCTIONS is used, i.e.
! the formula (3.184) of Attila Szabo, N. S. Ostlund Modern Quantum
! Chemistry is used
!
! Input: 
!       hcore: 1-e hamiltonian matrix
!       Coul : the Coulomb Interaction matrix(ces)
!       Exch : the Exchange Interaction matrix(ces)
!      vecnat: the 1-e orbitals
!
! Output:
!     Bare_energy: Contribution of 1-e terms to the total energy 
!     Coul_energy: Contribution of e-e Coulomb interaction to the total
!                  energy
!     Exch_energy: Contribution of electron Exchange to the total energy
!
!---------------------------------------------------------------------------

!   use global; use functional_m; use matrices; use orbocc; use energies
!   implicit none

!!.....Local variables
!   integer :: ia, ib, mu, nu
!   real(dp) :: hbare, hcoul, hexch, orb_fac
!   real(dp) :: TS_nlnn

!   Bare_energy = 0.d0
!   Coul_energy = 0.d0
!   Exch_energy = 0.d0

!   do ia=1, nnatorb
!!  if (occnum(ia,3).gt.smallocc) then
!      do mu = 1,nbasis
!         do nu = 1,nbasis
!            if(mu >= nu) then !Take care of the upper part of Hcore
!               hbare = hcore(mu,nu) * fac_h(ia)
!            else
!               hbare = hcore(nu,mu) * fac_h(ia)
!            endif ! mu >= nu
!            orb_fac = vecnat(mu,ia) * vecnat(nu,ia) !P_nm
!            Bare_energy = Bare_energy + hbare * orb_fac
!            do ib=1,nnatorb
!               if(mu >= nu) then !Take care of the upper part of Hcore
!                  hcoul = Coul(mu,nu,ib)
!                  hexch = Exch(mu,nu,ib)
!               else
!                  hcoul = Coul(nu,mu,ib)
!                  hexch = Exch(nu,mu,ib)
!               endif ! mu >= nu
!               Coul_energy = Coul_energy + hcoul*orb_fac*fac_c(ia,ib)
!               Exch_energy = Exch_energy + hexch*orb_fac*fac_e(ia,ib)
!            enddo
!         enddo
!      enddo
!!  endif !(occnum(ia,3).gt.smallocc) 
!   enddo
!
!   Totel_energy = Bare_energy + Coul_energy + Exch_energy
!
!   TS_ener = 0.d0
!   if(temp_dep) TS_ener = TS_nlnn()
!   
!   Totel_energy= Totel_energy + TS_ener
!   
!   Tot_energy = Rep_energy + Totel_energy
!
!end subroutine total_energy_AO

!------------------------------------------------------------------
! 1RDM on the basis functions:
! \gamma_\sigma(r,r\prime)=\sum \gamma_{kl} \chi_k^*(r\prime)*\chi_k(r)
subroutine gamma_on_basis
!..Global
   use global; use functional_m; use matrices; use orbocc
   implicit none

!..Local
   integer :: ia, k, l
   real(dp) :: vecpro

   if ( .not. allocated(gamma_b) ) allocate( gamma_b(nbasis,nbasis,3) )

   do k=1,nbasis
      do l=1,k
         gamma_b(k,l,1)=0.d0; gamma_b(k,l,2)=0.d0
         do ia=1, nnatorb
            vecpro=vecnat(k,ia)*vecnat(l,ia)
            gamma_b(k,l,1)=gamma_b(k,l,1)+vecpro*occnum(ia,1)
            gamma_b(k,l,2)=gamma_b(k,l,1)+vecpro*occnum(ia,2)
         enddo
         gamma_b(k,l,3)=gamma_b(k,l,1)+gamma_b(k,l,2)
         gamma_b(l,k,1)=gamma_b(k,l,1)
         gamma_b(l,k,2)=gamma_b(k,l,2)
         gamma_b(l,k,3)=gamma_b(k,l,3)
      enddo
   enddo
   return
end subroutine gamma_on_basis

subroutine off_diag_average()
! Same as in Mazzioti s paper JCP 133, 014104 (2010)
!..Global
   use global; use files; use functional_m; use matrices; use orbocc; use basis_set
   implicit none

!..Local
   integer :: k,l, nu
   real (dp) :: d
 
   call read_basis()

   call gamma_on_basis

   d = 0.d0
   do k = 1,nbasis
      do l = 1, k-1
         d=d+gamma_b(k,l,3)**2 
      enddo
   enddo

   d=sqrt(d/real(nbasis*(nbasis-1)))
   print*,'Harmonic average of Gamma off-diagonal:',d

   nu=0
   d = 0.d0
   do k = 1,nbasis
      do l = 1, k-1
         if ( iat_bas(k) /= iat_bas(l) ) then
            nu=nu+1
            d=d+gamma_b(k,l,3)**2 
         endif
      enddo
   enddo
   if ( nu > 0 ) then
      d=sqrt(d/real(nu))
      print*,'Harmonic average of Gamma off-diagonal (diff atoms):',d
   else
      print*,'For 1 atom systems Gamma off-diagonal (diff atoms) does not exist'
   endif

return
end subroutine off_diag_average
