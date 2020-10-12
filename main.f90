program hippo
!-----------------------------------------------------------------
!
!            Main of the RDMFT code
!
! The code minimizes the Energy w.r.t the natural orbitals
! and occupation numbers for functionals of 1-RDM
! It handles atomic and molecular systems.
! It relies on GAMESS for the calculation of 1 and 2 electron 
! integrals. In the routine readgamess.f90 there are details
! on how to modify GAMESS source code and how to use it in order
! to get all the information in a format readable by this program
!
!         created by N. N. Lathiotakis and I.Theophilou
!
!-----------------------------------------------------------------
! Conventions:
!  1) 3 space identation 
!    
!  2) Although everywhere 'implicit none' exists the following 
!     is strictly followed for readability:
!
!     a) variables starting with i,j,k,l,m,n are integer
!     b) variables starting with any other letter are usually
!        real. However they can be also logical or character. 
!-----------------------------------------------------------------

!..Global variables
   use global; use files; use functional_m; use matrices; use orbocc; use DFT_par
   use energies; use vary_occ_par; use basis_set; use grid_params
   use check_conv_m; use ovlap_four
   implicit none

!..Local variables
   integer :: nscf, iscf=0, nupocc, ncy_nvar, icy_nvar, icall_diis=0, Ndiis=20
   real(dp) :: xmu_in(2), xmu_range
   real(dp), allocatable :: xmu_energies(:)
   real(dp) :: rel_ene, conv_E1, conv_DD, ddmat
   integer ::  vmajor, vminor, vmicro

   logical :: ffile, conv, restart_save, test_orth, MLP_defs, Lcond
   real(dp) :: rat_ec_ek, rat_ec_ex, BBC1_strength
   real(dp) :: EHF

   integer :: ia, iorth, n

   real(dp) :: x0, y0, z0, x1, y1, z1
   integer :: iele1, iele2, n_occ

   real(dp) :: Pot_funct_bas
   external :: Pot_funct_bas

   real(dp), allocatable :: occnum_old(:,:), H_d(:,:), F1(:,:)
   real(dp) :: xmix_occ, Temper, shift, emax

   l_non_JK=.false.
   invert_dens=.true.

   ! initialize basis sets
   call basis_set_init_constants()

!..Read input file
   call read_inp(restart_save, nscf, &
                 nupocc, conv_E1, xmu_range)

!..Check the consistency of some input parameters
   if(iguess == 3) then
      inquire(file=rest_file, exist=ffile)
      if(.not.ffile) stop 'main: restart file does not exist!'
   endif

   if(cl_shell.and.(.not.common_mu)) then
      print*,'Main: You choose closed shell but with'
      print*,'not common_mu. That can be done using'
      print*,'cl_shell=.false. spin_state="multlet"'
      print*,'and all the occupancies as in closed shell'
      print*,'I presume thats not what you want and set '
      print*,'          common_mu=.true.'
      common_mu=.true.
   endif

   if (cl_shell) then
      if ( Npin2 /= Npin1 ) then
           Npin2 = Npin1
           print*,'Main: Warning: closed shell but Npin2 /= Npin1.'
           print*,'I set Npin2=Npin1 and assume Npin1 is correct!'
      endif
      if ( Nzer2 /= Nzer1 ) then
           Nzer2 = Nzer1
           print*,'Main: Warning: closed shell but Nzer2 /= Nzer1.'
           print*,'I set Nzer2=Nzer1 and assume Npin1 is correct!'
      endif
   endif
   if((Nzer1.lt.0) .or. (Nzer2.lt.0)) stop 'main: Nzer1 or Nzer2 < 0'

!..Read the gamess output
   call readgamess()

 
   if ( Npin1 >= nele(1) .and. .not.(Functional == 'RHF' .or. Functional == 'DFT' .or. Functional == 'RHA' ) ) then
        print*,' You can not pin ',Npin1, ' spin-orbitals if you have only ', nele(1), ' spin-up electrons!'
        stop 'main: insufficient variational freedom: Npin1 >= nele(1) !!'
   endif 
   if ( Npin2 > nele(2) .and. .not.(Functional == 'RHF' .or. Functional == 'DFT' .or. Functional == 'RHA' ) ) then
        print*,' You can not pin ',Npin2, ' spin-orbitals if you have only ', nele(2), ' spin-down electrons!'
        stop 'main: insufficient variational freedom: Npin2 >= nele(2) !!' 
   endif 
   if ( Nzer1 >= nnatorb-nele(1) .and. .not.(Functional == 'RHF' .or. Functional == 'DFT' .or. Functional == 'RHA') ) then
        print*,' You can not empty ',Nzer1, ' spin-orbitals if you have ', nele(1), ' spin-up electrons!'
        stop 'main:  Insufficiant variational freedom: Nzer1 >= nnatorb-nele(1) !!'
   endif 
   if ( Nzer2 >= nnatorb-nele(2) .and. .not.(Functional == 'RHF' .or. Functional == 'DFT' .or. Functional == 'RHA') ) then
        print*,' You can not empty ',Nzer2, ' spin-orbitals if you have ', nele(2), ' spin-down electrons!'
        stop 'main:  Insufficiant variational freedom: Nzer2 >= nnatorb-nele(2) !!'
   endif 
   if((Npin1.lt.0) .or. (Npin2.lt.0)) stop 'main: Npin1 or Npin2 < 0'
   if((Nzer1.lt.0) .or. (Nzer2.lt.0)) stop 'main: Nzer1 or Nzer2 < 0'

   SA1_fac = 1.d0

   if(rdnel.and.cl_shell.and.abs(xnele(1)-xnele(2))>1d-8) then
      stop 'main: for closed shell xnele(1) has to be equal to xnele(2)'
   endif

   if(.not.rdnel) then
      xnele(1)=nele(1)
      xnele(2)=nele(2)
   endif
   xnele(3) = xnele(1)+xnele(2)


!  print*,'Read the natural orbitals from ',gam_out
!  call read_nat_orbs()
!  print*,'Now the matrix element'
!  call average_nondiagonal
!  stop 'main:after average_nondiagonal'

   if (iguess /= 4) then
      maxorb=nnatorb-min(Nzer1,Nzer2)
      if ( Functional == 'RHF' .or. Functional == 'DFT' .or. Functional == 'RHA' ) maxorb = max(nele(1),nele(2))
   endif

   if ( Functional=='DBF' ) DBF_fun_RDM='BB3'

   call BB_orb_kinds()
   print*,'======================================================================='
   print*,'For spin_up:',xnele(1),'electrons and',ibond(1),'Strongly occ. orbitals'
   print*,'For spin_dn:',xnele(2),'electrons and',ibond(2),'Strongly occ. orbitals'
   print*,'======================================================================='

   if(functional == 'BB3' .or. (functional == 'DBF'.and. DBF_fun_RDM=='BB3')) then
      print*,'==============================================================='
      print*,'Orbital kinds for BBC3:'
      print*,'1: strongly occ, 2: bonding 3: antibonding 4: weakly occ.'
      print*,'Spin up:'
      write(*,'(20i3)')(kind_BBC3(ia,1),ia=1,nnatorb)
      print*,'Spin down:'
      write(*,'(20i3)')(kind_BBC3(ia,2),ia=1,nnatorb)
      print*,'==============================================================='
   endif

   if(functional == 'GPC') then
      l_non_JK=.true.
      f_GPC_13 = 1._dp
      fc12=1._dp
      fc23=1._dp
   endif

   if (functional=='LSH' .or. functional=='LSM'.or.functional=='LST') then
      if ( .not.allocated(f_LSH)) allocate (f_LSH(lnnatorb,2))
      do ia=1,nnatorb
         if(ia <= nele(1) ) then
            f_LSH(ia,1) = 1.d0
         else
            f_LSH(ia,1) = -1.d0
         endif 
         if(ia <= nele(2) ) then
            f_LSH(ia,2) = 1.d0
         else
            f_LSH(ia,2) = -1.d0
         endif 
         if (nele(2) < nele(1)) then
           f_LSH(ia,2) =-1.d0
         end if
      enddo
   endif

   if (functional=='LSH') then
      if ( nele(3) /= 2) stop 'Main: Lowdin Schull functional is for 2 electron systems only'
   endif
   if (functional=='LSM') then
      if ( nele(1) /= nele(2)) stop 'Main: Lowdin Schull Modified  functional is for closed shell systems only'
   endif
   if (functional=='LST') then
      if ( nele(1) /= (nele(2)+2)) stop 'Main: Lowdin Schull Modified functional is for triplets only'
   endif

!..Calculate the nuclear repulsive energy
   call repulsive_energy()

!..Parameters for the non local potential
   inquire(file='nl_params', exist=ffile)
   if(ffile) then
      open(unit=451,file='nl_params',status='unknown')
         read(451,*) do_nonloc; print*,do_nonloc
         read(451,*) non_loc_pot; print*,non_loc_pot
         read(451,*) Temper; print*,Temper
         read(451,*) shift; print*, shift
         read(451,*) xa
      close(451)
   else
      do_nonloc=.true. ! Do the non local potential
      non_loc_pot='P'  ! B: Baldsiefen P: Piris
!      non_loc_pot='N'  
      Temper=5.d0      ! Temperature for Baldsiefen
      shift=0.d0       ! Shift for Baldsiefen
      xa=0.d0
   endif

!..Read DFT and OEP parameters
   open(unit=731, file='DFT_OEP_pars', status='old')
   read(731,*) id_func_xc, id_func_c
   read(731,*) do_oep, int_pot_basis, scr_ch_const, posit_const, do_staroverov, do_EFO, l_hybc
   read(731,*) xld, svd_cut, xmix_OEP, xmix_DFT 
   read(731,*) pos_small, pos_penalty, pos_mix, small_e, Dn_lim, ahyb
   read(731,*) nrad, nang
   read(731,*) x0,y0,z0
   read(731,*) x1,y1,z1

   if(Functional=='DFT' .or. ((Functional=='RHF'.or. Functional == 'RHA').and.cl_shell) .or. Functional=='DBF') then
      print*, 'Main: Non local effective potential selected for Functional=',Functional,'... Swicthing to do_nonloc=.false.'
      do_nonloc=.false.
   endif

   
!..................DENSITY INVERSION.............................!
!    call read_basis()
!    call test_Xintegral2 
!      stop


   if (invert_dens) then

      print*, 'DENSITY INVERSION'

      call read_basis()
      call grid_setup
      call Invert_density(x0, y0, z0, x1, y1, z1)
      stop
   endif
     
!..................................................................!



   if(do_EFO) l_init_efo=.false.

!..Read basis functions
!..For usual coulomb and HF exchange terms they are not required. Used in DFT and OEP.
   if(do_oep .or. (Functional == 'DFT') .or. (Functional=='HDR') &
             .or. (Functional == 'DBF') .or. (Functional=='SLR') ) use_grid=.true.
   if(use_grid) then 
      call read_basis()
      call grid_setup
      if( do_oep ) then
         grad_S=.false.
         call read_basis_pot()
         call set_up_pot_basis
         open(unit=251,file="ovname",status='old')
         if ( .not.do_EFO ) then
            read(251,'(a60)')ov3_file
            call three_overlap(ov3_file)
         endif
      endif
   endif

   if(Functional == 'DFT' .or. Functional=='HDR' .or. Functional=='DBF' ) then
 !     call xc_f90_version(vmajor, vminor, vmicro)
!      write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro
!     id_func_xc = XC_LDA_X
!     id_func_c = XC_LDA_C_VWN

!     id_func_xc = XC_GGA_X_PBE
!     id_func_c = XC_GGA_C_PBE

!     id_func_xc = XC_GGA_X_B88
!     id_func_c = XC_GGA_C_LYP

!     id_func_xc = XC_GGA_X_PW91
!     id_func_c = XC_GGA_C_PW91

!     id_func_xc = XC_GGA_X_PW91
!     id_func_c = 0 !XC_GGA_C_P86
!.....Exchange Correlation Functional
      if ( cl_shell ) then
         call xc_f90_func_init(xc_func, xc_info, id_func_xc, XC_UNPOLARIZED)
!........Seperate Correlation Functional
         call xc_f90_func_init(c_func, c_info, id_func_c, XC_UNPOLARIZED)
      else
         call xc_f90_func_init(xc_func, xc_info, id_func_xc, XC_POLARIZED)
!........Seperate Correlation Functional
         call xc_f90_func_init(c_func, c_info, id_func_c, XC_POLARIZED)
      endif

      hyb_mix=0.d0
      if(xc_f90_info_family(xc_info) == XC_FAMILY_HYB_GGA) then 
!........For old libxc
!        call XC_F90_hyb_gga_exx_coef(xc_func, hyb_mix)
!........New Libxc
         call xc_f90_hyb_exx_coef(xc_func, hyb_mix)
         print*,'HYBRID Exx MIXING COEF: ',hyb_mix
      endif
   endif

   if (Functional == 'SLR') then
      id_func_xc=XC_LDA_X
      call xc_f90_func_init(xc_func, xc_info, id_func_xc, XC_UNPOLARIZED)
   endif

   if ( do_oep ) then
      if ( do_nonloc ) then
         print*,'main: Warning: You can not do both OEP and non local pot...'
         print*,'main: OEP with local potential will be done instead'
      endif

!.....N-1 electrons for the screening charge:
      veffnorm=xnele(1)+xnele(2)-1.0d0 

      n_occ=max(ibond(1), ibond(2))
      E_HOMO_HF=ennat_HF(n_occ)
   endif !do_oep

    if(functional == 'BB1' .or. functional == 'BBp' .or. functional == 'BBL' ) then
       strength=1.d0 !Fixed for BB1, BBL initial value for BBp
    endif

   if(functional == "POW" .or. functional == "SA1" .or. functional == "HDR" .or. &
      functional == "NEK" .or. functional == "SLR" .or. &
      (functional == "DBF" .and. DBF_fun_RDM=="POW") ) then
      if (functional == "POW") alpha=0.578d0 !Default value
      if (functional == "SLR") alpha=0.578d0 !Default value
      if (functional == "DBF") alpha=0.578d0 !Default value
      if (functional == "NEK") alpha=0.125d0 !Default value
      inquire(file='alpha', exist=ffile)
      if(.not.ffile) then 
         print*,'main: file "alpha" not found. a=0.578 will be used for POW functional'
      else
         open(unit=193, file='alpha', status='old')
            read(193,*) alpha
         close(193)
         print*,'main: file "alpha" found. New value for alpha:', alpha
      endif
   endif

   if(functional == "MPF") then
      inquire(file='P_mull', exist=ffile)
      if(.not.ffile) stop 'main: file "P_mull" required for MPF functional does not exist!'
      open(unit=193, file='P_mull', status='old')
         read(193,*) P_mull
         if ((P_mull <= 0.d0) .or. (P_mull >= 1.d0) ) then
            stop 'P_mull should be strictly between 0 and 1'
         endif
      close(193)
   endif

   if(functional == "MLP") then
      inquire(file='MLP_params', exist=ffile)
      if(.not.ffile) then
         print*,'main: file "MLP_params" required for MLP functional does not exist!'
         MLP_defs=.true.
      else
         open(unit=193, file='MLP_params', status='old')
            read(193,*,err=143) fit_a1, fit_a2, fit_b1
         close(193)
      endif
      goto 144
 143     print*,'main: Wrong data in file MLP_params'
         MLP_defs=.true.
 144  continue
      if ( MLP_defs ) then
         print*, 'The default MPL params will be used: 126.3101 2213.33 2338.64'
         fit_a1=126.3101d0;  fit_a2=2213.33d0; fit_b1=2338.64d0
      endif
   endif

   if(nnatorb > nbasis) then
      print*, 'nnatorb:',nnatorb,'nbasis:',nbasis
      stop 'main:nnatorb cannot be larger than nbasis!'
   endif

   if ( Functional == 'HDR' ) then
      open(unit=53, file='HDR.dat', status='old')
         read(53,*) xmix_DF
      close(53)
!     xmix_DF= 0.5d0
      xmix_RD= 1.d0 !- xmix_DF
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!  x_t=0.d0
!  do igrid=1, ngrid
!     x=x_grid(igrid)
!     y=y_grid(igrid)
!     z=z_grid(igrid)
!     call C_density(x,y,z,rhoA,rhoB,rho)
!     x_t = x_t + rhoB*w_grid(igrid)
!  enddo
!  print*,'Total number of electrons:',x_t
!  stop 'a'

!..Construct the initial Fock matrices
   call calc_occ_factors(occnum, fac_h, fac_c, fac_e)

   call construct_f(0)

!  call test_Xintegral2
!  stop 'ab'
   
!..Calculate the initial Total Energy
   call total_energy_MO(1)
   write(6,113) Bare_energy, Coul_energy, Exch_energy, &
                Totel_energy, Rep_energy, TS_ener,iscf,  Tot_energy
   if(rdnel) then
      print*,'This total energy is obtained by using the occupancies'
      print*,'n_el(i) read from input file. The xnele(1), xnele(2) will'
      print*,'taken into account in the minimization w.r.t. the occ. numbers'
   endif
   HF_exc_ener= Exch_energy

!..Check orbital orthonormality:
   call normal_test(test_orth,1d-8,1)
   if(test_orth) then
      print*,'Initial orbitals are orthonormal!'
   else
      print*,'Initial orbitals are NOT orthonormal!'
      do iorth=1,10
         call symmet_orthog()
         call normal_test(test_orth,1.d-10,0)
         if(test_orth) exit
      enddo
      if(iorth == 10) stop 'main: Symmetric orthogonalization failed!'
!.....Construct the initial Fock matrices
      call calc_occ_factors(occnum, fac_h, fac_c, fac_e)
      call construct_f(0)

!.....Calculate the initial Total Energy
      call total_energy_MO(1)
      write(6,113) Bare_energy, Coul_energy, Exch_energy, &
                Totel_energy, Rep_energy, TS_ener,iscf,  Tot_energy

   endif

   tot_old = Tot_energy !Previous total energy value for check convergence
   EHF = tot_old


!.. START OF SCF CYCLE
   print*,'   SELF CONSISTENT ITERATIONS FOR ORBITALS:'
   print*,'----------------------------------------------------------------'
   conv=.false.

   if ( do_oep ) then
      allocate ( Inc_par(lnnatorb, lnnatorb) )
!     open(unit=23,file='Inc_par',status='old')
      Inc_par=.false.
!     do i=1,200
!        read(23,*,end=127) iii, jjj
!        Inc_par(iii,jjj)=.true.
!     enddo
!127   continue
!     close(23)
    endif
   
    call correl_entropy()

!   call average_nondiagonal()

SC: do iscf = 1, nscf
      print*, 'ITERATION NU:',iscf
!     if ( do_oep .and. iscf > 5 ) nupocc=5

      if((iscf == 1.or.mod(iscf,nupocc) == 0)) then
!........Initialize mu
         if(iscf == 1) then
            allocate ( xmu_energies(lnnatorb))
            call orb_energies(xmu_energies)
            if(iguess == 3) then
               xmu_in(1)= xmu_read(1)
               xmu_in(2)= xmu_read(2)
            else !(iguess < 3)
               xmu_in(1) = 2.d0*xmu_energies(nele(1))/fac_h(nele(1))
               if(common_mu) then
                  xmu_in(2) = xmu_in(1)
               else
                  if(nele(2) > 0) xmu_in(2) = 2.d0*xmu_energies(nele(2))/fac_h(nele(2))
               endif
            endif !(iguess == 3)
            deallocate (xmu_energies)
         endif !(iscf == 1) 

!........Variation of Occupation numbers
         if(functional /= 'RHF' .and. functional /= 'DFT' .and. functional /= 'RHA' ) then
            ncy_nvar=1
            if(functional == 'SA1') then
               ncy_nvar=5
            endif
            do icy_nvar=1, ncy_nvar
               select case(meth_occ)
                  case(1) !Very Reliable
                     call vary_occ_th(xmu_in, xmu_range)
                  case(2) !Reliable
                     call vary_occ(xmu_in, xmu_range)
                  case(3) !Very Reliable
                      if ((.not.do_oep).or.(iscf>0)) then
                         if(.not.allocated(occnum_old)) allocate(occnum_old(lnnatorb,3))
                         xmix_occ=1.0d0 !Mixing parameter for occupation numbers
                         occnum_old=occnum
                         if (constr == 0 ) then
                            print*,'calculation of occupation numbers with only ensemble constraints'
!                           call BB_orb_kinds()
                            call vary_occ_NAG(xmu_in)
                         else if (constr == 1.and. functional /='PN5') then
                            print*,'Occupation number minimization with Borland-Dennis conditions' 
                            call vary_occ_NAG_36(xmu_in)
                         endif
!                        call Optimize_FD(iscf)
                         occnum=(1.d0-xmix_occ)*occnum_old+xmix_occ*occnum
                         call calc_occ_factors(occnum, fac_h, fac_c, fac_e)
!                         call construct_f(0)
!                        if ( do_oep ) call Fit_Fermi(iscf,xmu_in)
                      endif
                      call correl_entropy()
                  case(4) !Reliable
                     call vary_occ_th_n(xmu_in)
                  case default
                     print*,'Wrong meth_occ: ', meth_occ
                     print*,'use: 1 for vary_occ_th,' 
                     print*,'     2 for vary_occ,'
                     print*,'     3 for vary_occ_NAG'
                     print*,'     4 for vary_occ_th with strict N conservation'
               end select
               SA1_fac=sum(occnum(1:nnatorb,3)**(2.d0*alpha))
               print*,'sum n_i^2a = ',SA1_fac
               SA1_fac= xnele(3) / SA1_fac
               print*,'new beta: ', SA1_fac
            enddo

!......Calculate total energy after variation of occ numbers
            call construct_f(0)
!           if(iscf == 1) then
               print*,'Energy after occ num variation:'
               call total_energy_MO(1) !Called after construct_f!!!!
               write(6,113) Bare_energy, Coul_energy, Exch_energy, &
                         Totel_energy, Rep_energy, TS_ener,iscf,  Tot_energy
!           endif
         else
            print*,'For RHF (ROHF) or RHA or DFT occ. numbers are not varied'
            if ( rdnel .and. iscf==1 ) then
               iele1 = xnele(1)
               iele2 = xnele(2)
               occnum=0
               if (iele1 >= 1) occnum(1:iele1,1) = 1.d0; occnum(iele1+1,1)=xnele(1)-iele1 
               if (iele2 >= 1) occnum(1:iele2,1) = 1.d0; occnum(iele2+1,2)=xnele(2)-iele2 
               occnum(:,3)=occnum(:,1)+occnum(:,2)
               call calc_occ_factors(occnum, fac_h, fac_c, fac_e)
            endif
         endif
      endif

!     goto 238
      if(do_oep) then
         if ( do_EFO ) then
!           call eff_orbital_STPD(x0,y0,z0,x1,y1,z1,'potentialf')
!           call eff_orbital_L_BFGS_B(iscf,x0,y0,z0,x1,y1,z1,'potentialf')
!           call eff_orbital_NIK_TOM(iscf,x0,y0,z0,x1,y1,z1,'potentialf')
             call eff_orbital_NAG(iscf,x0,y0,z0,x1,y1,z1,'potentialf')
         else   
            call analytic_Veff(iscf)
         endif
         call construct_F(0)
      elseif(do_nonloc) then
!        call Vary_NLP_cg(iscf)
         if (non_loc_pot=='P') then
            print*,'Piris Minimization for the orbitals'
            call non_loc_eff_pot(iscf)
         elseif (non_loc_pot=='B') then
            call non_loc_eff_pot_BG(Temper,shift,iscf)
         elseif (non_loc_pot=='H') then
            call non_loc_eff_pot_HF(iscf)
         elseif (non_loc_pot=='N') then
!           call  NONLOCAL_RDMFT_SIMPLE(iscf)
!           call  NONLOCAL_RDMFT33(iscf)
!           call  NONLOCAL_RDMFT(iscf)
         else
            stop 'main: unknown orbital variation method'
         endif
         
!        call Vary_NLP_NAG(iscf)
      elseif(((Functional == 'RHF' .or. Functional == 'RHA') .and.cl_shell) .or. (Functional == 'DFT' .and. &
           cl_shell) .or. Functional == 'DBF' ) then

!........Find orbitals through diagonalization of Hamiltonian
         if (.not.allocated(H_d)) allocate( H_d(lnbasis, lnbasis) , F1(lnbasis, lnbasis) )
         if ( iscf > 1) then  ! F mixing or Enable DIIS extrapolation
            if ( iscf <= 5 ) then
               if ( functional == 'DFT' .or. functional == 'RHA' ) F(:,:,1)=xmix_DFT*F(:,:,1)+(1.d0-xmix_DFT)*F_old(:,:,1)
            else
               F1(:,:)=F(:,:,1)
               call DIIS_extrapolation(icall_diis, Ndiis, F1, H_d, emax)
               F(:,:,1)=H_d
            endif
            F_old(:,:,1) = F(:,:,1)
         else
            if (.not. allocated(F_old) ) allocate( F_old(lnbasis,lnbasis,lnnatorb) )
            F_old(:,:,1) = F(:,:,1)
         endif
         call diagon_lapack(lnbasis,F(:,:,1),ovlap,ennat,vecnat)
         call construct_F(0) !Rebuild Hamiltonian using new orbitals
!     elseif (do_staroverov) then
!       call model_local_pot(x0, y0, z0, x1, y1, z1)  
!       call construct_F(0)
      else ! Do the notorious procedure of orbital minimization in RDMFT
!........Vary the orbitals (Direct minimization)
         print*,'Direct Minimization for the orbitals'
         call vary_orbs_direct
!        call vary_orbs_all(iscf)
      endif

!238  continue

!....Calculate the total energy again:
      call total_energy_MO(1) !Called after construct_f!!!!
      write(6,113) Bare_energy, Coul_energy, Exch_energy, &
                   Totel_energy, Rep_energy, TS_ener, iscf, Tot_energy

!....Vary the strength of the BBp functional:
      if(functional == 'BBp') then
         rat_ec_ek=(Tot_energy-HF_tot_ener)/HF_kin_ener
         rat_ec_ex=(Tot_energy-HF_tot_ener)/HF_exc_ener
         print*,'*****Strength Calculation:******'
         print*,'Old strength: ',strength
         strength = BBC1_strength(rat_ec_ek,2)
         print*,'New strength: ',strength,' (ratio: ',rat_ec_ek,')'
      endif
          
!....Convergence check
      conv_DD=1.d-5
      if (Functional == 'RHF' .or. Functional == 'RHA' .or. Functional == 'DFT' .or. Functional == 'DBF' &
          .or.  do_OEP) then
          Lcond = .true. !Disable Largange Matrix Hermiticity condition
      else
          call Lagrange_Hermiticity(Lcond)
      endif
      call check_conv(iscf, conv_E1, conv_DD, rel_ene, ddmat,  conv)
      conv=conv.and.Lcond
      write(6,110) Tot_energy,rel_ene,conv_E1,ddmat,conv_DD

      print*,'------------------------------------------------------------'
      print*,'Iteration',iscf

!....Save restart information (in every scf cycle)
      if(restart_save) then
         call save_restart_info(xmu_in)
      endif

      if ( conv ) exit SC
      call energy_spectrum(2, xmu_in)
 
   enddo SC !iscf
!..END OF SCF CYCLE

  print*,'------------------------------------------------------------'
   if(conv) then
      print*, 'Convergence achieved in',iscf,' cycles'
   else
      print*, 'Convergence not achieved after ',iscf-1, ' cycles'
   endif

  print*,'------------------------------------------------------------'
  print*,'------------------------------------------------------------'
  print*,'||||||||             END OF SCF CYCLE             ||||||||||'
  print*,'------------------------------------------------------------'
  print*,'------------------------------------------------------------'

   print*,'ENERGY FINAL RESULTS: '
   write(6,113) Bare_energy, Coul_energy, Exch_energy, &
           Totel_energy, Rep_energy, TS_ener, iscf, Tot_energy 
   if(iguess == 2) then
      Print*,'Correlation Energy: ', Tot_energy-EHF
   else
      Print*,'DEenergy (from initial): ', Tot_energy-EHF
   endif
   call off_diag_average()

!..Find Orbital energies
   call energy_spectrum(2, xmu_in)
   call energy_spectrum(3, xmu_in)

110 format('RESULT: Total Energy: ',f18.10,/,13x,'RESULT:  Relative Energy Change:',e13.6,' (',e13.6,')'/,&
                                     13x,'RESULT:  Density Matrix  Change:',e13.6,' (',e13.6,')')
113 format('   One-electron     Energy: ',f18.10,/, &
           '   Hartree          Energy: ',f18.10,/, &
           '   Exch./Corr.      Energy: ',f18.10,/, &
           '   Total Electronic Energy: ',f18.10,/, &
           '   Repulsive nuclei Energy: ',f18.10,/, &
           '   Entropy       term (TS): ',f18.10,/, &
    i3,')'    '   *****TOTAL  ENERGY*****: ',f18.10,/) 

   if((Functional /= 'RHF') .and. (Functional /= 'DFT') .and. (Functional /= 'DBF') .and. (Functional /= 'RHA') ) then 
!.....Orbitals orthonormality test
      call normal_test(test_orth,1d-8,1)

!.....Check for the hermiticity of Lagrange multipliers
   endif

   if(nscf /= 0 .and. do_oep) then
      call plot_potential(x0,y0,z0,x1,y1,z1,'potentialf')
   endif
   call save_restart_info(xmu_in)
      if ( Functional /= 'DFT' .and. (.not.do_OEP) .and. do_staroverov ) then
          call energy_spectrum(2, xmu_in)
          call model_local_pot(x0, y0, z0, x1, y1, z1)  
          call construct_F(0)
      endif
!..Check Positivity of the screening density in LDA
!  if (Functional == 'DFT') then
!     call  Calc_Scr_dens(5000, x0, y0, z0, x1, y1, z1)
!     call Plot_Vxc_DFT_r(5000, x0, y0, z0, x1, y1, z1)
!  endif

!  call plot_orbital(1,7,x0,y0,z0,x1,y1,z1)

!  if(functional /= 'DFT' .and. functional /= 'RHF' ) then
!     call correl_entropy()
!  endif

!  call spin_ensemble()!occupations messed up after this call


!   close(111); close(222); close(333); 
    close(42)

   print*,'*********************************'
   print*,'*     TELOS KALO, OLA KALA      *'
   print*,'*********************************'

end program hippo

!---------------------------------------------------------------------------
! Reads the input file 'inp'
!--------------------------------------------------------------------------
subroutine read_inp(restart_save, nscf,&
                    nupocc, conv_E1, xmu_range)
!..Global
   use global; use files; use functional_m; use vary_occ_par; use orbocc
   implicit none

!..Arguments
   integer, intent(out) :: nscf, nupocc
   logical, intent(out) :: restart_save
   real(dp), intent(out) :: conv_E1, xmu_range

!..Local
   character(60) :: inpfile
   integer :: i

   read*,inpfile
!..Read the input
   open(unit=2,file=inpfile,status='old')
   read(2,*); read(2,'(a60/a60/a60/a60)',err=1)gam_out, gam_out_pot, int_file, rest_file; read(2,*)
   read(2,*); read(2,*,err=2)Functional, meth_occ, vary_str, strength; read(2,*)
   read(2,*); read(2,*,err=3)nnatorb, Npin1, Npin2, Nzer1, Nzer2; read(2,*)
   read(2,*); read(2,*,err=4)nscf, conv_E1, nupocc; read(2,*)
   read(2,*); read(2,*,err=5)iguess, restart_save; read(2,*)
   read(2,*); read(2,*,err=6); read(2,*)niter_on, nbracket, niter_on1, n_vary_meth; read(2,*)
   read(2,*); read(2,*,err=7)crit_on, crit_on1, step; read(2,*)
   read(2,*); read(2,*,err=8) xmu_range,nmu,crit_mu,rdnel,xnele(1),xnele(2),common_mu; read(2,*)
   read(2,*); read(2,*,err=9)step_dir_init, scale_up, scale_down,ntry_dir,nst_orb_dir; read(2,*)
   read(2,*); read(2,*,err=10)cl_shell; read(2,*,err=11)spin_st

   allocate(n_el(nnatorb))
   read(2,*,err=12)(n_el(i),i=1,nnatorb); read(2,*)
   read(2,*); read(2,*,err=13)temp_dep,temperat; read(2,*)
   read(2,*); read(2,*,err=14)constr
   close(2)

   print*,'Successful reading of file ', inpfile
   return
1  stop 'Error reading gam_out, gam_out_pot, int_file, rest_file.'
2  stop 'Error reading Functional, meth_occ, vary_str, strength.'
3  stop 'Error reading nnatorb, Npin1, Npin2, Nzer1, Nzer2.'
4  stop 'Error reading nscf, conv_E1, nupocc.'
5  stop 'Error reading iguess, restart_save.'
6  stop 'Error reading niter_on, nbracket, niter_on1, n_vary_meth.'
7  stop 'Error reading crit_on, crit_on1, step.'
8  stop 'Error reading xmu_range,nmu,crit_mu,rdnel,xnele(1),xnele(2),common_mu.'
9  stop 'Error reading step_dir_init, scale_up, scale_down,ntry_dir,nst_orb_dir.'
10 stop 'Error reading cl_shell'
11 stop 'Error reading spin_st'
12 stop 'Error reading n_el(i)'
13 stop 'Error reading temp_dep,temperat'
14 stop 'Error reading constr'

end subroutine read_inp

subroutine save_restart_info(xmu_in)
use global; use files; use orbocc; use energies; implicit none; 
   real(dp), intent(IN) :: xmu_in(2)

!..Local
   integer :: ia, i

   open(unit=10,file=rest_file,action='write')
      write(10,*,err=59)'NBASIS, NNATORB:'
      write(10,'(2i5)',err=59) nbasis, nnatorb
      write(10,*,err=59)'MU FOR SPIN up, down:'
      write(10,*,err=59)xmu_in(1), xmu_in(2)
      write(10,*,err=59)'OCC NUMBERS FOR SPIN UP:'
      write(10,'(4f20.16)',err=59) (occnum(ia,1),ia=1,nnatorb)
      write(10,*,err=59)'OCC NUMBERS FOR SPIN DOWN:'
      write(10,'(4f20.16)',err=59) (occnum(ia,2),ia=1,nnatorb)
      write(10,*,err=59)'NATURAL ORBITALS:'
      write(10,'(4e25.16)',err=59) ((vecnat(i,ia),i=1,nbasis),ia=1,nbasis)
      write(10,*,err=59)'ORBITAL MINIMIZATION QUANTITIES'
      write(10,'(4e20.10)',err=59) (step_dir(ia),ia=1,nnatorb)
      write(10,*,err=59)'ENERGIES:'
      write(10,113,err=59) Bare_energy, Coul_energy, Exch_energy, &
                           Totel_energy, Rep_energy, TS_ener, Tot_energy 
      write(10,'(4e20.12)',err=59) (ennat(ia),ia=1,nnatorb)
   close(10)
   goto 60
59 print*,'main: WARNING: The restart file can not be writen!'
60 continue
113 format('One-electron     Energy: ',f18.10,/, &
           'Hartree          Energy: ',f18.10,/, &
           'Exchange         Energy: ',f18.10,/, &
           'Total Electronic Energy: ',f18.10,/, &
           'Repulsive nuclei Energy: ',f18.10,/, &
           'Entropy       term (TS): ',f18.10,/, &
           '*****TOTAL  ENERGY*****: ',f18.10,/) 

end subroutine save_restart_info

