!----------------------------------------------------------------------
!   Module for very general parameters 
!----------------------------------------------------------------------
module params_general

!  lnatom : maximum number of atoms. 
   integer :: lnatom=200

   integer :: lnnatorb, lnbasis, lnbasis_pot, lnintgr, ln4ov=10000
!     lnbasis : Leading dimesion of number of basis functions
!     lnnatorb :   "        "     of number of natural orbitals 
!     lnintgr : Leading dimension of number of integrals stored 
!               All the above are assigned a value later and the
!               matrices are allocated

!..Basis set dimensions max number of gaussians per basis, max number of basis
   integer, parameter :: lnga=50, lbsa=10000 !lsbsa should be a power of 10 If modified
! modification of gamess should be done accordingly in pack routine in int2a.src

! Global definition of real type
   integer, parameter :: dp=selected_real_kind(p=14, R=30)

   real(dp), parameter :: small = 1.e-12_dp, &
                              small_de = 1.e5_dp, &
                                  zero = 1.e-24_dp, &
                                   one = 1.0_dp-zero, &
                                   two = 2.0_dp-zero-zero, &
                              smallocc = 1.e-12_dp, &
                               der_min = -1.e-4_dp, &
                                    pi = 3.1415926535897932385_dp, &
                                 tovpi = 2.0_dp/pi, & 
                               piovtwo = pi/2.0_dp, &
                                 tovsqpi = 1.1283791670955126_dp, & !2/sqrt(pi)
                                    sqpi = 1.7724538509055160_dp, & !sqrt(pi)
                                 ha_ev = 27.211396_dp, & 
                             one_ov_4pi = 0.079577471545947667884_dp !1/(4pi)
   real(dp), parameter :: xnovm(0:5,0:5) = reshape( (/ 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &
                                              1._dp, 1._dp, 2._dp, 3._dp, 4._dp, 5._dp, &
                                              1._dp, 2._dp, 1._dp, 3._dp, 6._dp,10._dp, &
                                              1._dp, 3._dp, 3._dp, 1._dp, 4._dp,10._dp, &
                                              1._dp, 4._dp, 6._dp, 4._dp, 1._dp, 5._dp, &
                                              1._dp, 5._dp,10._dp,10._dp, 5._dp, 1._dp /) ,  (/ 6,6/))
!    small : something positive and small
!     zero : something positive close to zero (SMALLER than small)
! smallocc : some small value for occ. numbers (initial)
!  der_min : lowest acceptable negative derivative for varying occ
!            numbers. Purpose: to keep the iner shells at occ. nu. = 1.
!    xnovm : N over M N!/(M-N)!/M!: xnofm(n,m) is max(m,n) over min(m,n) i.e. order does not matter. 

!..The Boltzman constant in Hartree/K
   real(dp), parameter :: Boltz_konst=3.166815342e-6_dp

end module params_general

!----------------------------------------------------------------------
!  Global module
!----------------------------------------------------------------------
module global
   use params_general
   save
! Global variables
   integer :: nbasis, nbasis_pot, nnatorb, natom, nintgr, nele(3), iguess 

   integer :: meth_occ

   integer :: iocc, iunocc, il_pot, iu_pot, jl_pot, ju_pot, nn_pr_occ

! Atomic data
   character(12), allocatable   :: name_ch(:)
   real(dp), allocatable :: xcoo(:), ycoo(:), zcoo(:), &
                                   charge(:)
   real(dp), allocatable :: x_AB(:,:), y_AB(:,:), z_AB(:,:), R_AB(:,:)

! Initial mu data
   real(dp) :: xmu_read(2)

! Parameters related to fractional total number of particles
   logical rdnel ! if true then read number of electrons: xnele(nspin)
   real(dp) :: xnele(3)

! Orbital variation parameters
   real(dp), allocatable :: step_dir(:) 
   real(dp) :: step_dir_init, scale_up, scale_down
   integer :: ntry_dir, nst_orb_dir

! Check if effective orbital is initialized
   logical :: l_init_efo

! Spin variables
   character(7) :: spin_st
   integer, allocatable :: n_el(:)

! Timing variables
   real(dp) :: time_tot, time_fock, time_occ, time_vec, &
                      time_nbas, time_tote, time_rdgm, &
                      time_sumint, time_lagr

! The temperature
   logical  :: temp_dep !t: add entropy term f: plain T=0 calculation
   real(dp) :: Temperat, entropy

! Bonding and antibonding for two spins (used by BBC3)
   integer :: ibond(2), ianti(2)
   integer, allocatable :: kind_BBC3(:,:)

! Powers of nbasis
   integer :: nbasis2, nbasis3

   real(dp) :: veffnorm

   integer :: n4ov
! The restart file

   real(dp) :: ggama=0.0_dp

   real(dp) :: exp_th=1.0_dp

   logical :: xi_analytic=.true.

   logical, allocatable :: Inc_par(:,:)

! For the non local effective potential
   real(dp) :: x_mix_1=0.1_dp, de_conv=1.e-12_dp
   integer :: nconv=500


   interface
     function oct_secnds(time)
       use params_general
       real(dp), intent(in) :: time
       real(dp) :: oct_secnds
     end function oct_secnds
   end interface
end module global

module files
   use params_general
   save
    
   character(60) :: gam_out, gam_out_pot, int_file, ov3_file, xi3_file, rest_file, gam_out_mix, int_file_mix
   
end module files

!----------------------------------------------------------------------
!  Energies module: Energies and components
!----------------------------------------------------------------------
module energies
   use params_general
   save
   real(dp) :: Rep_energy, Tot_energy
   real(dp) :: Bare_energy, Coul_energy, Exch_energy, Totel_energy
   real(dp) :: HF_kin_ener, HF_exc_ener, HF_tot_ener, TS_ener
   real(dp) :: ExcDFT_ener, EcDFT_ener
   real(dp) :: E_HOMO_HF, E_HOMO
end module energies

!----------------------------------------------------------------------
!  Functional_m module: global parameters related to the functional
!----------------------------------------------------------------------
module functional_m
   use params_general
   save
   character(3) :: Functional, DBF_fun_RDM
   real(dp) :: strength
   real(dp) :: alpha, xa
   real(dp) :: SA1_fac
   real(dp) :: P_mull
   real(dp) :: fit_a0, fit_a1, fit_a2, fit_a3, fit_a4, fit_b1, fit_b2 !For pade function
   real(dp) :: fit_a5, fit_a6, fit_a7 !For pade function
   real(dp) :: alpf, bta, gma !For pade function
   real(dp) :: xi_fit, xli_fit, q_fit, xg_fit !For pade function
   logical, allocatable :: si_term(:)
   logical :: cl_shell, common_mu
   logical :: vary_str, field_calc
   logical :: lex_projection
   logical :: do_OEP, do_EFO, do_ceda, int_pot_basis, scr_ch_const, posit_const, do_staroverov,invert_dens
   logical :: do_nonloc
   logical :: grad_S, l_non_JK
   integer :: maxorb
   integer :: constr
   character(1) :: non_loc_pot='P'
   real(dp) :: xld, xmix_OEP, svd_cut, pos_small, pos_mix, pos_penalty, small_e, Dn_lim
   real(dp) :: xmix_DF, xmix_RD
   real(dp), allocatable :: f_LSH(:,:)

!..Alpha parameter for hybrid constrained potential
   real(dp) :: ahyb
   logical :: l_hybc
   real(dp), allocatable :: VHF(:,:)
end module functional_m

!----------------------------------------------------------------------
!  Matrices module: Fock and other big matrices 
!----------------------------------------------------------------------
module matrices
   use params_general
   save
   
! The fock matrix:
   real(dp), allocatable :: F(:,:,:) !The Fock Matrix
   real(dp), allocatable :: F_xc(:,:,:), F_ha(:,:,:), F_haM(:,:,:), F_x(:,:,:)!The Fock Matrix
   real(dp), allocatable :: Coul(:,:,:), Exch(:,:,:) 

! The overlap of basis functions:
   real(dp), allocatable :: ovlap(:,:), ov_inv(:,:)
   real(dp), allocatable :: ovlap_pot(:,:), ovlap_pot_chi(:,:)

! The 1-e Hamiltonian
   real(dp), allocatable :: kin(:,:) !The kinetic term
   real(dp), allocatable :: hcore(:,:) ! kinetic + v_external

! The x,y,z integrals
   real(dp), allocatable :: x_integ(:,:), &
                                   y_integ(:,:), &
                                   z_integ(:,:)
   real(dp) :: EField_x, EField_y, EField_z

! Hcore, Coulomb, Exchange on the natural orbitals basis
   real(dp), allocatable :: CoulNO(:,:) 
   real(dp), allocatable :: ExchNO(:,:)
   real(dp), allocatable :: HcoreNO(:,:)
   
! Four atom basis overlap:\int phi_i(r) \phi_j(r) \phi_k(r) \phi_l(r)
   real(dp), allocatable :: FourOvlap(:)
   integer, allocatable :: i_ov(:), j_ov(:), k_ov(:), l_ov(:) 

! Three atom basis overlap:
  real(dp), allocatable :: ovlap3(:,:,:), ovlap3_xi(:,:,:), Tr_prod(:,:) 

! The <a|F^(b)|b> matrix:
   real(dp), allocatable :: e_ab(:,:), ee_ab(:,:),  eee_ab(:,:)

! The effective Hamiltonian
   real(dp), allocatable :: H_eff(:,:), H_eff_c(:,:), H_eff_xc(:,:), H_eff_HF(:,:), V_ij_HF(:,:)
   real(dp), allocatable :: dV(:), V_bs(:), V_bs_o(:), Veff(:,:), Ovchixi(:,:), Veff2(:,:)
   real(dp), allocatable :: dVNL(:,:)
   real(dp), allocatable :: d_x(:), dX_i(:), dX_o(:)
   real(dp), allocatable :: dEdVab(:,:)
   real(dp), allocatable :: vec_eff(:), vec_eff0(:)
   real(dp), allocatable :: vec_eff_n(:,:), vec_eff0_n(:,:), Vef(:,:), Vef2(:), vec_hf(:,:)
   integer :: nvec_eff

   real(dp), allocatable :: H_scr(:,:)

! The hessian 
   real(dp), allocatable :: Hess(:,:,:)

! For the non JK terms
   real(dp), allocatable :: XI1(:), XI2(:), XI3(:)
   real(dp) :: q_non_JK

   real(dp), allocatable :: DMM(:,:)
   real(dp), allocatable :: Vcoul(:,:)
   real(dp), allocatable :: Vexch(:,:)

   real(dp), allocatable :: S3(:,:,:)

end module matrices

!-------------------------------------------------------------------
! orbocc module: occuparion numbers, natural orbitals and related
!-------------------------------------------------------------------
module orbocc
   use params_general
   save
! Occupation numbers and natural orbitals
   real(dp), allocatable :: occnum(:,:) !The occupation numbers
   real(dp), allocatable :: vecnat(:,:) !The natural orbitals
   real(dp), allocatable :: vecnat_targ(:,:), vecnat0(:,:)  
   real(dp), allocatable :: ennat(:), ennat_HF(:)

   real(dp), allocatable :: vec_A_null(:,:), vec_A_notn(:,:), enn_A_notn(:)
   integer :: Npin1, Npin2, Nzer1, Nzer2
   integer :: n_null_eiA, n_notn_eiA, ln_null_eiA, ln_notn_eiA
   logical :: zero_null, zero_notn

! Occupation number factors depending on the functional.
! See functionals.f
   real(dp), allocatable :: fac_c(:,:), fac_e(:,:)
   real(dp), allocatable :: fac_h(:)
   real(dp) :: f_GPC_13, f_non_JK, fd_non_JK1, fd_non_JK2, fd_non_JK3, fc12, fc23

! Spectral decomposition of Gamma on the original basis:
   real(dp), allocatable :: gamma_b(:,:,:)

   real(dp):: beta_occ(2), xmu_occ(2)
end module orbocc

!--------------------------------------------------------------------
! integrals module: 2-electron integral parameters
!-------------------------------------------------------------------
module integrals
   use params_general
   save
!..The two-e integrals (indirect indexing)
   real(dp), allocatable :: twoin(:)

!..The indirect index matrices
   integer, allocatable :: mu(:), nu(:), lambda(:), sigma(:)
   integer(8), allocatable :: indpacked(:)

end module integrals

module integrals_mix
   use params_general
   save
!..The two-e integrals between orbital and auxiliary basis
   real(dp), allocatable :: twin_mix(:)

!..Number of non zero mixed orbitals
   integer :: nint_mix

!..The indirect index matrices
   integer, allocatable :: m_mix(:), n_mix(:), l_mix(:), s_mix(:)

end module integrals_mix

!------------------------------------------------------------------------
! four overlap module
!------------------------------------------------------------------------
module ovlap_four
   use params_general
   save
   character(60) :: ov4_file
   real(dp), allocatable :: ovlap4(:)
   integer, allocatable :: mu_o4(:), nu_o4(:), lambda_o4(:), sigma_o4(:)
   integer(8), allocatable :: indpacked_o4(:)
end module ovlap_four

!----------------------------------------------------------------------
! basis set parameters module for the analytic reconstruction of basis functions
!----------------------------------------------------------------------
module basis_set
   use global
   integer, allocatable :: iat_bas(:) ! iat_bas(i): which atom nu the i basis corresponds
   integer, allocatable :: iat_bas_pot(:) ! iat_bas(i): which atom nu the i basis corresponds
   integer, allocatable :: ityp_bas(:) ! ityp_bas(i): type of i basis element, 0=s, 1=p, 2=d ...
   integer, allocatable :: ityp_bas_pot(:) ! ityp_bas(i): type of i basis element, 0=s, 1=p, 2=d ...
   integer, allocatable :: nga_bas(:) ! Number of gaussians in i basis element
   integer, allocatable :: nga_bas_pot(:) ! Number of gaussians in i basis element
   real(dp), allocatable :: expo_bas(:,:) ! expo_bas(iga,ibas) exponent of iga gaussian in ibas element
   real(dp), allocatable :: expo_bas_pot(:,:) ! expo_bas(iga,ibas) exponent of iga gaussian in ibas element
   real(dp), allocatable :: coef_bas(:,:) ! coef_bas(iga,ibas) coefficient of .... 
   real(dp), allocatable :: coef_bas_pot(:,:) ! coef_bas(iga,ibas) coefficient of .... 
   real(dp), allocatable :: expo_x(:), expo_y(:), expo_z(:) !exponent in the angular momentum prefactor
   real(dp), allocatable :: expo_x_pot(:), expo_y_pot(:), expo_z_pot(:) !exponent in the angular momentum prefactor
   real(dp), allocatable :: fact_norm(:,:) !Normalization factor for primitive Gaussians 
   real(dp), allocatable :: fact_norm_pot(:,:) !Normalization factor for primitive Gaussians 
   real(dp) :: s_fun_norm, p_fun_norm, d_fun_norm, f_fun_norm, g_fun_norm, h_fun_norm, expo_bas_min
 contains

!..Some constants dor the normalization of basis functions
   subroutine basis_set_init_constants()
     s_fun_norm=(8.0_dp/pi/pi/pi)**0.25_dp
     p_fun_norm=(128.0_dp/pi/pi/pi)**0.25_dp
     d_fun_norm=(2048.0_dp/pi/pi/pi/9.0_dp)**0.25_dp
     f_fun_norm=(32768.0_dp/pi/pi/pi/225.0_dp)**0.25_dp
     g_fun_norm=((2._dp/pi)**0.75_dp)*((256._dp/105._dp)**0.5_dp)
     h_fun_norm=((2._dp/pi)**0.75_dp)*((32765._dp/30240._dp)**0.5_dp)
   end subroutine basis_set_init_constants

end module basis_set

!----------------------------------------------------------------------------------------------
!  Parameters related to the grid implementation
!----------------------------------------------------------------------------------------------
module grid_params
   use global

   logical :: grid_stat=.false.
   logical :: use_grid=.false.
   logical :: lgradset=.false.
   integer :: ngrid ! number of grid points
   real(dp), allocatable :: x_grid(:),& ! x coord of grid point
                                   y_grid(:),& ! y coord of grid point
                                   z_grid(:),& ! z coord of grid point
                                   w_grid(:),& ! weight of grid point 
                                   v_grid(:)   ! volume corresponding to grid point

!..Parameters for the new spherical grid
   integer :: nrad=150, nang=350 ! Radial and angular mesh sizes
!..Choices for nang: 6, 14, 26, 38,  50, 74,  86, 110, 146, 170, 194, 230, 266,  302,&
!                    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074,&
!                    3470, 3890, 4334, 4802, 5294, 5810
   real(dp) :: tol_wg=1.e-12_dp ! weight tolerance for accepting grid points
   real(dp) :: Del_fuz=0.13_dp ! fuzzy parameter in Voronoi faces
   real(dp) :: rsb_cut_off=20.0_dp ! Cutoff in radial integration (in Bragg-Slater radii)
   real(dp), allocatable :: xang(:), yang(:), zang(:), wang(:) ! Angular grid xang, yang, zang on unitary sphere

!..Parameters for the old rectangular grid
   real(dp) :: coo_extend=4.0_dp ! space added at ends
   real(dp) :: dx=2.e-1_dp, dy=2.e-1_dp,  dz=2.0e-1_dp ! grid size 
   real(dp) :: grid_volume
   real(dp) :: small_f=1.e-12_dp
   integer :: n_x, n_y, n_z

!..Quantities on the grid
   real(dp), allocatable :: bas_f_grid(:,:)
   real(dp), allocatable :: charg_pot(:,:), bas_f_grid_pot(:,:), X_bas_pot(:), C_bas_pot(:)
   real(dp), allocatable :: vec_nat_grid(:,:), orbgridsq(:,:), vec_nat_targ_grid(:,:), orbgridsq_targ(:,:)
   real(dp), allocatable :: xnorm_pot(:)
   real(dp), allocatable :: chdens_gr(:,:)
   real(dp), allocatable :: grad_bas_f_grid(:,:,:)
contains
  
   subroutine calc_grid_volume
      grid_volume = SUM(w_grid)
   end subroutine calc_grid_volume

   subroutine set_grid_volume(corr_0)
       real(dp) :: corr_0
       grid_volume = dx*dy*dz
       corr_0 = 4.0_dp*pi* (grid_volume*3.0_dp/(4.0_dp*pi))**(2.0_dp/3.0_dp)
   end subroutine set_grid_volume
end module grid_params

!----------------------------------------------------------------------
! vary_occ_par : parameters vor vary_occ.... subroutines
!----------------------------------------------------------------------
module vary_occ_par
   use global
   save
   integer :: niter_on, niter_on1, nbracket, nmu, ispin, ispin_op
   integer :: n_vary_meth
   real(dp) :: crit_on, crit_on1, crit_mu, step  
end module vary_occ_par

module DFT_par
   use global; use xc_f90_types_m; use xc_f90_lib_m
   real(dp), allocatable :: Exc_DFT(:), Ec_DFT(:), Vlxc(:,:), Vlc(:,:)
   real(dp), allocatable :: Vlxc_sig(:,:), Vlc_sig(:,:)
   real(dp), allocatable :: Vxc_mn_DFT(:,:), Vc_mn_DFT(:,:), F_old(:,:,:)
   logical :: Fexist=.false.
   real(dp) :: xmix_dft=0.4_dp
   logical :: gga_x=.false., gga_c=.false.

!..For Libcx
   integer :: id_func_xc=XC_LDA_X, id_func_c=XC_LDA_C_VWN
   TYPE(xc_f90_pointer_t) :: xc_func, c_func
   TYPE(xc_f90_pointer_t) :: xc_info, c_info
   real(dp) :: hyb_mix=0.0_dp

end module DFT_par

module HFD_par
   use global;
   real(dp), allocatable :: V_hfd(:)
end module 

!----------------------------------------------------------------------
!.....For the wrapper function grad, grad_th, grad_th_n:
!----------------------------------------------------------------------
module func_der_args
   use global
   save
   integer :: nstart, ispin1, ispin1_op
   real(dp) :: xmu1
   real(dp), allocatable :: occnum1(:,:) 
   real(dp), allocatable :: DE_Dn1(:), DE_Dtheta1(:)
   real(dp), allocatable :: h(:), theta1(:)
   logical, allocatable :: occ_var(:)
end module func_der_args

module check_conv_m
   use global
   save
   real(dp) :: Tot_old, DE_old, dn_old
   real(dp), allocatable ::  Dmat(:,:), Dmat_old(:,:)
end module check_conv_m

!----------------------------------------------------------------------
! To allocate the memory for most important global objects
!----------------------------------------------------------------------
subroutine array_allocate()

!..Global
   use global
   use matrices
   use functional_m
   use orbocc
   use func_der_args
   use basis_set
   implicit none

!..FROM GLOBAL:
   allocate (step_dir(lnnatorb)) 
   allocate (kind_BBC3(lnnatorb,2))


!..FROM FUNCTIONAL_M
   allocate (si_term(lnnatorb))

!..FROM MATRICES:
!..Fock Matrices and components
   allocate (  F(lnbasis, lnbasis, lnnatorb), &
               F_xc(lnbasis, lnbasis, lnnatorb), &
               F_x(lnbasis, lnbasis, lnnatorb), &
               F_ha(lnbasis, lnbasis, lnnatorb), &
               F_haM(lnbasis, lnbasis, lnnatorb), &
               Coul(lnbasis, lnbasis, lnnatorb), &
               Exch(lnbasis, lnbasis, lnnatorb) )

!..One electron matrices
   allocate (  ovlap(lnbasis, lnbasis), &
                 kin(lnbasis, lnbasis), &
               hcore(lnbasis, lnbasis)  )

!..Polarization matrices
   allocate (  x_integ(lnbasis, lnbasis), &
               y_integ(lnbasis, lnbasis), &
               z_integ(lnbasis, lnbasis)  )

!..The J_ab, K_ab, H_ab
   allocate (  CoulNO(lnnatorb, lnnatorb), & 
               ExchNO(lnnatorb, lnnatorb), &
              HcoreNO(lnnatorb, lnnatorb) )

!..Orbital variation Lagrange multipliers
   allocate (e_ab(lnbasis, lnbasis), & 
             ee_ab(lnbasis, lnbasis), &
             eee_ab(lnbasis, lnbasis) )

!..Effective Hamiltonian
   allocate (   H_eff(lnbasis,lnbasis), &
              H_eff_c(lnbasis,lnbasis), &
             H_eff_xc(lnbasis,lnbasis)  )

!..FROM ORBOCC
!..Occupation numbers, natural orbitals
   allocate(occnum(lnnatorb,3),vecnat(lnbasis,lnbasis),vecnat0(lnbasis,lnbasis),&
             ennat(lnbasis),ennat_HF(lnbasis), vec_hf(lnbasis,lnbasis))

!..Occupation related factors
   allocate(fac_c(lnnatorb,lnnatorb), &
            fac_e(lnnatorb,lnnatorb), &
            fac_h(lnnatorb) )

!..FROM basis_set
   allocate(  iat_bas(lnbasis), ityp_bas(lnbasis), nga_bas(lnbasis), &
              expo_bas(lnga,lnbasis), coef_bas(lnga,lnbasis),&
              fact_norm(lnga,lnbasis),&
              expo_x(lnbasis), expo_y(lnbasis), expo_z(lnbasis) ) 

! FROM func_der_args
   allocate(occnum1(lnnatorb,3), DE_Dn1(lnnatorb), DE_Dtheta1(lnnatorb), &
            h(lnnatorb), theta1(lnnatorb))
   allocate(occ_var(lnnatorb))

end subroutine array_allocate

