!----------------------------------------------------------
! THIS IS A FILE THAT THE FUNCTIONALS ARE DEFINED
! Entropy functionals are also included at the end
!        created by N. N. Lathiotakis and I.Theophilou
!----------------------------------------------------------

subroutine calc_occ_factors(occ, fac_h, fac_c, fac_e)
! Just a wrapper calling the appropriate functional routine

!..Global
   use global; use functional_m
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

   select case (functional)
      case("RHF") 
         call RHF(occ, fac_h, fac_c, fac_e) ! Restricted Hartree-Fock
      case("RHA") 
         call RHA(occ, fac_h, fac_c, fac_e) ! Restricted Hartree-Fock
      case("GUM") 
         call GUM(occ, fac_h, fac_c, fac_e) ! Goedecker Umrigar
      case("BB0") 
         call BB0(occ, fac_h, fac_c, fac_e) ! Mueller, Buijse Baerends
      case("BB1")
         call BB1(occ, fac_h, fac_c, fac_e) ! BBC1
      case("BB2")
         call BB2(occ, fac_h, fac_c, fac_e) ! BBC2 
      case("BB3")
         call BB3(occ, fac_h, fac_c, fac_e) ! BBC3
      case("BBp")
         call BB1(occ, fac_h, fac_c, fac_e) ! BBC++
      case("CHF")
         call CHF(occ, fac_h, fac_c, fac_e) ! Corected Hartree Fock (Csanyi-Arias)
      case("CGA")
         call CGA(occ, fac_h, fac_c, fac_e) ! Csanyi-Goedecker-Arias
      case("POW")
         call POW(occ, fac_h, fac_c, fac_e) ! Power functional
      case("PNO")
         call PNO(occ, fac_h, fac_c, fac_e) ! Piris NO functional
      case("PNP")
         call PNP(occ, fac_h, fac_c, fac_e) ! Piris NO functional with pinning
      case("MPF")
         call MPF(occ, fac_h, fac_c, fac_e) ! Mueller general-p functional
      case("MLP")
         call MLP(occ, fac_h, fac_c, fac_e) ! Marques Lathiotakis Pade
      case("SA1")
         call SA1(occ, fac_h, fac_c, fac_e) ! New Sangeeta idea
      case("AC3")
         call AC3(occ, fac_h, fac_c, fac_e) ! JCP129, 164105 (2009)
      case("TST") 
          call TST(occ, fac_h, fac_c, fac_e)
      case("OSM")
         call OSM(occ, fac_h, fac_c, fac_e) 
      case("DFT") ! For DFT calculation
         call DFT(occ, fac_h, fac_c, fac_e) 
      case("SLR") ! Slater version of Power functional
         call DFT(occ, fac_h, fac_c, fac_e)
      case("HDR") 
         call HDR(occ, fac_h, fac_c, fac_e) 
      case("DBF")
         call DBF(occ, fac_h, fac_c, fac_e)
      case("NEK")
         call NEK(occ, fac_h, fac_c, fac_e)
      case("LSH")
         call LSH(occ, fac_h, fac_c, fac_e)
      case("PN5") ! PNOF-5
         call PN5(occ, fac_h, fac_c, fac_e)
      case("BBN") 
         call BBN(occ, fac_h, fac_c, fac_e) ! Mueller modified
      case("LSM")
         call LSM(occ, fac_h, fac_c, fac_e) !Modified LSH
      case("LST")
         call LST(occ, fac_h, fac_c, fac_e) !Modified LSH for triplet
      case("GPC") 
          call GPC(occ, fac_h, fac_c, fac_e) !Functional for three electrons based on the expansion of a 3 determinental wavefunction
      case default
         print*,'Functional:',functional
         stop 'calc_occ_factors: not implemented functional'
   end select 

end subroutine calc_occ_factors

!----------------------------------------------------------------
subroutine calc_der_occ_factors(ia, ispin, occ, fac_c_d, fac_e_d)
! Just a wrapper calling the appropriate functional routine

!..Global
   use global; use functional_m
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb)
   real(dp), intent(out) :: fac_e_d(lnnatorb)
      
   select case (functional)
      case ("RHF") !Restricted Hartree Fock
         call RHF_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("RHA") !Restricted Hartree Fock
         call RHA_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("GUM") !Goedecker Umrigar
         call GUM_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("BB0") !Mueller functional
         call BB0_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("BB1") !BBC1
         call BB1_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("BB2") !BBC2
         call BB2_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("BB3") !BBC3
         call BB3_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("BBp")
         call BB1_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("CHF") 
         call CHF_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("CGA")
         call CGA_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("POW")
         call POW_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("PNO")
         call PNO_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("PNP")
         call PNP_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("MPF")
         call MPF_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("MLP")
         call MLP_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("SA1")
         call SA1_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("AC3")
         call AC3_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("TST")
         call BB0_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("OSM")
         call OSM_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("DFT")
         call DFT_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("SLR")
         call DFT_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("HDR") 
         call HDR_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("DBF") 
         call DBF_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("NEK") 
         call NEK_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("LSH") 
         call LSH_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("PN5") 
         call PN5_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("BBN") !Mueller modified
         call BBN_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case ("LSM") 
         call LSM_d(ia, ispin, occ, fac_c_d, fac_e_d)!LSH modified
      case ("LST") 
         call LST_d(ia, ispin, occ, fac_c_d, fac_e_d)!LSH modified for triplets
      case ("GPC") 
         call GPC_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case default
         stop 'calc_der_occ_factors: not implemented functional'
   end select

end subroutine calc_der_occ_factors


!-------------------------------------------------------------------
! RESTRICTED (OPEN-SHELL) HARTREE FOCK (RHF)
subroutine RHF(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!.....Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!.....Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -occ(ia,1)*occ(ib,1)-occ(ia,2)*occ(ib,2)
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine RHF
!---------------AND THE DERIVATIVE----------------------------
subroutine RHF_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)
      
!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = -ocib
   enddo

end subroutine RHF_d

!-------------------------------------------------------------------
! RESTRICTED (OPEN-SHELL) HARTREE (NO FOCK!!!!!)
subroutine RHA(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!.....Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!.....Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = 0._dp
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0._dp

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine RHA
!---------------AND THE DERIVATIVE----------------------------
subroutine RHA_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)
      
!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = 0._dp
   enddo

end subroutine RHA_d


!--------------------------------------------------------------------
! THE MUELLER FUNCTIONAL (or Buisze-Baerents) Phys Lett. 105A, 446 (1984), 
! Mol Phys 100, 401 (2002)
subroutine BB0(occ, fac_h, fac_c, fac_e)!
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -sqrt(occ(ia,1)*occ(ib,1)) &
                        -sqrt(occ(ia,2)*occ(ib,2))

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine BB0
!---------------AND THE DERIVATIVE----------------------------
subroutine BB0_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)
   enddo

end subroutine BB0_d

!--------------------------------------------------------------
! The Goedecker Umrigar functional (PRL81 866, 1998)
! see also PRA 72, 030501(R) (2005) for open shells
subroutine GUM(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)
      
!..Local
   integer :: ia,ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = occ(ia,1)*occ(ib,2) &
                      + occ(ia,2)*occ(ib,1)
         fac_e(ia,ib) = 0.0_dp
         if(ia/=ib) then
            fac_c(ia,ib) = fac_c(ia,ib) + occ(ia,1)*occ(ib,1) &
                                        + occ(ia,2)*occ(ib,2)
            fac_e(ia,ib) =  - sqrt(occ(ia,1)*occ(ib,1)) &
                            - sqrt(occ(ia,2)*occ(ib,2))
            fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)
         endif
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
!         if(abs(fac_c(ia,ib))<1d-8) fac_c(ia,ib)=0.d0
!         if(abs(fac_e(ia,ib))<1d-8) fac_e(ia,ib)=0.d0

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine GUM
!---------------AND THE DERIVATIVE----------------------------
subroutine GUM_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op 
      fac_e_d(ib) = 0._dp  
      if(ib/=ia) then
         fac_c_d(ib) = fac_c_d(ib) + ocib
         fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)
      endif
   enddo

end subroutine GUM_d

!-------------------------------------------------------------
! BBC1,  Journal of Chem Phys, 122, 204102-1 (2005).
subroutine BB1(occ, fac_h, fac_c, fac_e)
   use global; use functional_m
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib
   real(dp) :: exch1, exch2

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         exch1 = -sqrt(max(occ(ia,1)*occ(ib,1),small))
         exch2 = -sqrt(max(occ(ia,2)*occ(ib,2),small))

!........BBC1 correction (strength: PRB 75, 195120, (2007))
         if(ia/=ib) then
            if(ia>ibond(1).and.ib>ibond(1)) exch1 = -strength*exch1 
            if(ia>ibond(2).and.ib>ibond(2)) exch2 = -strength*exch2 
         endif

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*(exch1+exch2)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine BB1
!---------------AND THE DERIVATIVE----------------------------
subroutine BB1_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use functional_m
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = max(occ(ia, ispin),small)
   do ib = 1, nnatorb
      ocib= max(occ(ib,ispin),small)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)

!.....BBC1 correction
      if(ia/=ib) then
         if(ia>ibond(ispin).and.ib>ibond(ispin)) fac_e_d(ib)=-strength*fac_e_d(ib) 
      endif
   enddo

end subroutine BB1_d

!--------------------------------------------------------------------
! BBC2, Journal of Chem Phys, 122, 204102-1 (2005).
subroutine BB2(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib
   real(dp) :: exch1, exch2

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         exch1 = -sqrt(occ(ia,1)*occ(ib,1))
         exch2 = -sqrt(occ(ia,2)*occ(ib,2))

         if(ia/=ib) then !Here go the BBC1, BBC2 corrections
            if(ia>ibond(1).and.ib>ibond(1)) exch1 = -exch1 
            if(ia>ibond(2).and.ib>ibond(2)) exch2 = -exch2 
            if(ia<=ibond(1).and.ib<=ibond(1)) exch1 = -occ(ia,1)*occ(ib,1)
            if(ia<=ibond(2).and.ib<=ibond(2)) exch2 = -occ(ia,2)*occ(ib,2)
         endif

!........And the factor of a half for the double counting
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*(exch1+exch2)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine BB2
!---------------AND THE DERIVATIVE----------------------------
subroutine BB2_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)
      
!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)
      if(ia/=ib) then !BBC1, BBC2 corrections
         if(ia>ibond(ispin).and.ib>ibond(ispin)) fac_e_d(ib)=-fac_e_d(ib) 
         if(ia<=ibond(ispin).and.ib<=ibond(ispin)) fac_e_d(ib)=-ocib
      endif
   enddo

end subroutine BB2_d

!-------------------------------------------------------------
! BBC3_old (Journal of Chem Phys, 122, 204102-1 (2005).
! BBC3_old uses one bonding and one anti-bonding orbital
! An improved version respecting degeneracies of bonding
! and antibonding is given below
subroutine BB3_old(occ, fac_h, fac_c, fac_e) !spin dependent
   use global
   use functional_m
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib
   integer :: ilim1, ilim2
   real(dp) :: exch1, exch2

   ilim1=xnele(1); if(xnele(1)-ilim1 > 5*small) ilim1=ilim1+1
   ilim2=xnele(2); if(xnele(2)-ilim2 > 5*small) ilim2=ilim2+1

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia !NOTE ia >= ib 
!........The non SI part added for all:
         fac_c(ia,ib) = occ(ia,1)*occ(ib,2) &
                      + occ(ia,2)*occ(ib,1)
         exch1=0.0_dp
         exch2=0.0_dp

!........Add SI terms for bonding antibonding (spin resolved):
         if(((ia==ianti(1)).or.(ia==ibond(1))).and.(ib==ia)) then
            fac_c(ia,ib)= fac_c(ia,ib) + occ(ia,1)*occ(ib,1)
            exch1 = -sqrt(occ(ia,1)*occ(ib,1))
         endif
         if(((ia==ianti(2)).or.(ia==ibond(2))).and.(ib==ia)) then
            fac_c(ia,ib)= fac_c(ia,ib) + occ(ia,2)*occ(ib,2)
            exch2 = -sqrt(occ(ia,2)*occ(ib,2))
         endif
          
         if(ia/=ib) then
            fac_c(ia,ib)= fac_c(ia,ib) + occ(ia,1)*occ(ib,1) &
                                       + occ(ia,2)*occ(ib,2)
            exch1 = -sqrt(occ(ia,1)*occ(ib,1))
            exch2 = -sqrt(occ(ia,2)*occ(ib,2))

!...........BBC1 correction:
            if(ia>ilim1.and.ib>ilim1) exch1 = -exch1 
            if(ia>ilim2.and.ib>ilim2) exch2 = -exch2 

!...........BBC2 correction:
            if(ia<=ilim1.and.ib<=ilim1) exch1 = -occ(ia,1)*occ(ib,1)
            if(ia<=ilim2.and.ib<=ilim2) exch2 = -occ(ia,2)*occ(ib,2)

!...........Include antibonding in BBC2 correction
            if(ia == ianti(1) .and. ib <= ilim1) then
               if(ib /= ibond(1)) exch1 = -occ(ia,1)*occ(ib,1)
            endif
            if(ia == ianti(2) .and. ib <= ilim2) then
               if(ib /= ibond(2)) exch2 = -occ(ia,2)*occ(ib,2)
            endif
         endif

!........And the factor of a half for the double counting
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*(exch1+exch2)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine BB3_old
!---------------AND THE DERIVATIVE----------------------------
subroutine BB3_d_old(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op
   integer :: ilim

   ilim=xnele(ispin); if(xnele(ispin)-ilim > 5*small) ilim=ilim+1

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)

      fac_c_d(ib) = ocib_op 
      fac_e_d(ib) = 0._dp
!.....Add SI terms for bonding-antibonding
      if(((ia==ianti(ispin)).or.(ia==ibond(ispin))).and.(ib==ia)) then
         fac_c_d(ib) = fac_c_d(ib) + ocib
         fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)
      endif
      if(ia/=ib) then 
         fac_c_d(ib) = fac_c_d(ib) + ocib
         fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)

!........BBC1 correction
         if(ia > ilim.and.ib > ilim) fac_e_d(ib)=-fac_e_d(ib)

!........BBC2 correction
         if(ia <= ilim .and. ib <= ilim) fac_e_d(ib)=-ocib

!........Add BBC2 for the antibonding
         if(ia == ianti(ispin) .and. ib <= ilim) then
            if(ib /= ibond(ispin)) fac_e_d(ib)=-ocib
         endif
         if(ib == ianti(ispin) .and. ia <= ilim) then
            if(ia /= ibond(ispin)) fac_e_d(ib)=-ocib
         endif
      endif
   enddo
end subroutine BB3_d_old
!-------------------------------------------------------------------

!-------------------------------------------------------------
! BBC3 (Journal of Chem Phys, 122, 204102-1 (2005).
! Degeneracies of bonding/antibonding are taken into account
! as in JCP 128, 184103 (2008)
subroutine BB3(occ, fac_h, fac_c, fac_e) !spin dependent
   use global
   use functional_m
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib
   real(dp) :: exch1, exch2

! See auxil.f90 subroutine BB_orb_kinds. 
! kind_BBC3(ia,ispin) = 1   : orbital 'ia' is strongly occupied for ispin
! kind_BBC3(ia,ispin) = 2   : orbital 'ia' is bonding for ispin
! kind_BBC3(ia,ispin) = 3   : orbital 'ia' is anti-bonding for ispin
! kind_BBC3(ia,ispin) = 4   : orbital 'ia' is weakly occupied for ispin

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia !NOTE ia >= ib 
!........The non SI part added for all:
         fac_c(ia,ib) = occ(ia,1)*occ(ib,2) &
                      + occ(ia,2)*occ(ib,1)
         exch1=0.0_dp
         exch2=0.0_dp

!........Add SI terms for bonding antibonding (spin resolved):
         if((kind_BBC3(ia,1)==2.or.kind_BBC3(ia,1)==3).and.(ib==ia)) then
            fac_c(ia,ib)= fac_c(ia,ib) + occ(ia,1)*occ(ib,1)
            exch1 = -sqrt(occ(ia,1)*occ(ib,1))
         endif
         if((kind_BBC3(ia,2)==2.or.kind_BBC3(ia,2)==3).and.(ib==ia)) then
            fac_c(ia,ib)= fac_c(ia,ib) + occ(ia,2)*occ(ib,2)
            exch2 = -sqrt(occ(ia,2)*occ(ib,2))
         endif
          
         if(ia/=ib) then
            fac_c(ia,ib)= fac_c(ia,ib) + occ(ia,1)*occ(ib,1) &
                                       + occ(ia,2)*occ(ib,2)
            exch1 = -sqrt(occ(ia,1)*occ(ib,1))
            exch2 = -sqrt(occ(ia,2)*occ(ib,2))

!...........BBC1 correction:
            if(kind_BBC3(ia,1) >= 3 .and. kind_BBC3(ib,1) >= 3) exch1 = -exch1 
            if(kind_BBC3(ia,2) >= 3 .and. kind_BBC3(ib,2) >= 3) exch2 = -exch2 

!...........BBC2 correction:
            if(kind_BBC3(ia,1) <= 2 .and. kind_BBC3(ib,1) <= 2) exch1 = -occ(ia,1)*occ(ib,1)
            if(kind_BBC3(ia,2) <= 2 .and. kind_BBC3(ib,2) <= 2) exch2 = -occ(ia,2)*occ(ib,2)

!...........Include antibonding in BBC2 correction
            if(kind_BBC3(ia,1)==3 .and. kind_BBC3(ib,1)==1) exch1 = -occ(ia,1)*occ(ib,1)
            if(kind_BBC3(ia,2)==3 .and. kind_BBC3(ib,2)==1) exch2 = -occ(ia,2)*occ(ib,2)
         endif

!........And the factor of a half for the double counting
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*(exch1+exch2)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine BB3
!---------------AND THE DERIVATIVE----------------------------
subroutine BB3_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)

      fac_c_d(ib) = ocib_op 
      fac_e_d(ib) = 0._dp
!.....Add SI terms for bonding-antibonding
!     if(((ia==ianti(ispin)).or.(ia==ibond(ispin))).and.(ib==ia)) then
      if((kind_BBC3(ia,ispin)==2 .or. kind_BBC3(ia,ispin)==3) &
                                                .and. ia==ib) then
         fac_c_d(ib) = fac_c_d(ib) + ocib
         fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)
      endif
      if(ia/=ib) then 
         fac_c_d(ib) = fac_c_d(ib) + ocib
         fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)

!........BBC1 correction
         if(kind_BBC3(ia,ispin)>=3 .and. kind_BBC3(ib,ispin)>=3) &
                    fac_e_d(ib)=-fac_e_d(ib)

!........BBC2 correction
         if(kind_BBC3(ia,ispin)<=2 .and. kind_BBC3(ib,ispin)<=2) &
                    fac_e_d(ib)=-ocib

!........Add BBC2 for the antibonding
         if(kind_BBC3(ia,ispin)==3 .and. kind_BBC3(ib,ispin)==1) &
                    fac_e_d(ib)=-ocib

         if(kind_BBC3(ib,ispin)==3 .and. kind_BBC3(ia,ispin)==1) &
                    fac_e_d(ib)=-ocib
      endif
   enddo
end subroutine BB3_d
!-------------------------------------------------------------------

!--------------------------------------------------------------------
! THE CHF FUNCTIONAL (Corrected-Hartree-Fock Csanyi-Arias, PRB 61, 7348 (2000))
subroutine CHF(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      if(occ(ia,1) < 0._dp) stop 'CHF:1 negative occ'
      if(occ(ia,2) < 0._dp) stop 'CHF:2 negative occ'
      if(occ(ia,1) > 1._dp) stop 'CHF:3 occ >1 '
      if(occ(ia,2) > 1._dp) stop 'CHF:4 occ >1 '
      do ib=1,ia
         if(occ(ib,1) < 0._dp) stop 'CHF:5 negative occ'
         if(occ(ib,2) < 0._dp) stop 'CHF:6 negative occ'
         if(occ(ib,1) > 1._dp) stop 'CHF:7 occ >1 '
         if(occ(ib,2) > 1._dp) stop 'CHF:8 occ >1 '
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -occ(ia,1)*occ(ib,1) &
                        -occ(ia,2)*occ(ib,2) &
             -sqrt(occ(ia,1)*(1._dp-occ(ia,1))*occ(ib,1)*(1._dp-occ(ib,1)))&
             -sqrt(occ(ia,2)*(1._dp-occ(ia,2))*occ(ib,2)*(1._dp-occ(ib,2)))

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine CHF
!---------------AND THE DERIVATIVE----------------------------
subroutine CHF_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = min(occ(ia, ispin),1.d0-small)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = -ocib - 0.5_dp*(1._dp-2._dp*ocia)*  &
                     sqrt((ocib*(1._dp-ocib))/(ocia*(1._dp-ocia)))
   enddo

end subroutine CHF_d

!--------------------------------------------------------------------
! THE CGA FUNCTIONAL (Csanyi-Goedecker-Arias, PRA 65, 032510 (2002))
subroutine CGA(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -0.5_dp*(occ(ia,1)*occ(ib,1) &
                        +occ(ia,2)*occ(ib,2) &
             +sqrt(occ(ia,1)*(2._dp-occ(ia,1))*occ(ib,1)*(2._dp-occ(ib,1)))&
             +sqrt(occ(ia,2)*(2._dp-occ(ia,2))*occ(ib,2)*(2._dp-occ(ib,2))))

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine CGA
!---------------AND THE DERIVATIVE----------------------------
subroutine CGA_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = -0.5_dp*(ocib + (1._dp-1._dp*ocia)* &
                     sqrt((ocib*(2._dp-ocib))/(ocia*(2._dp-ocia))))
   enddo

end subroutine CGA_d

!--------------------------------------------------------------------
! POWER FUNCTIONAL, Phys. Rev. B, 78, 201103(R) (2008), Lathiotakis et al PRA, 79 (2009)
! (modification of the Mueller functional: alpha exponent instead of 1/2 )
subroutine POW(occ, fac_h, fac_c, fac_e)
   use global; use functional_m
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -(occ(ia,1)*occ(ib,1))**alpha &
                        -(occ(ia,2)*occ(ib,2))**alpha

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine POW
!---------------AND THE DERIVATIVE----------------------------
subroutine POW_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use functional_m
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = - alpha*(ocib**alpha)*(ocia**(alpha-1._dp))
   enddo

end subroutine POW_d

!--------------------------------------------------------------
! The Piris Natural Orbital Functional (PNOF) 
! see IJQC 106, 1093 (2006); JCP 123, 214102 (2005);
subroutine PNO(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)
      
!..Local
   integer :: ia,ib
   real(dp) :: occia1, occia2, occib1, occib2

!..External
   real(dp), external :: xLambda
   
   do ia=1,nnatorb
      occia1=min(occ(ia,1),1._dp-zero)
      occia2=min(occ(ia,2),1._dp-zero)
      fac_h(ia) = occia1 + occia2
      do ib=1,ia
         occib1=min(occ(ib,1),1._dp-zero)
         occib2=min(occ(ib,2),1._dp-zero)
         fac_c(ia,ib) = occia1*occib2 &
                      + occia2*occib1
         fac_e(ia,ib) = 0.0_dp
         if(ia/=ib) then
            fac_c(ia,ib) = fac_c(ia,ib) + occia1*occib1 &
                                        + occia2*occib2
            fac_e(ia,ib) =  - xLambda(occia1,occib1) &
                            - xLambda(occia2,occib2)
            fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)
         endif
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine PNO
!---------------AND THE DERIVATIVE----------------------------
subroutine PNO_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

!..External
   real(dp), external :: dxLambda

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = min(occ(ia, ispin), 1._dp-zero)

   do ib = 1, nnatorb
      ocib= min(occ(ib,ispin), 1._dp-zero)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op 
      fac_e_d(ib) = 0._dp 
      if(ib/=ia) then
         fac_c_d(ib) = fac_c_d(ib) + ocib
         fac_e_d(ib) = -dxLambda(ocia,ocib)
      endif
   enddo

end subroutine PNO_d
!---------------Auxiliary functions-------------------
function xLambda(ocia,ocib)
!..Global variables
   use global; use orbocc
   implicit none

!..The function
   real(dp) :: xLambda
!..Arguments
   real(dp) :: ocia,ocib

!..External
   real(dp), external :: H_theta

!..Local
   real(dp) :: xnocc=0.45_dp

   xLambda = (1._dp - 2._dp*H_theta(xnocc-ocib)*H_theta(xnocc-ocia))* &
             sqrt(ocib*ocia) &
           + H_theta(ocib-xnocc)* H_theta(ocia-xnocc)* &
             sqrt((1._dp-ocib)*(1._dp-ocia))
   
end function xLambda
!---------------------------------------------------------
function dxLambda(ocia,ocib)
!..Global variables
   use global; use orbocc
   implicit none

!..The function
   real(dp) :: dxLambda
!..Arguments
   real(dp) :: ocia,ocib

!..Local
   real(dp) :: xnocc=0.45_dp

!..External
   real(dp), external :: H_theta

   dxLambda = (1._dp - 2._dp*H_theta(xnocc-ocib)*H_theta(xnocc-ocia))* &
             sqrt(ocib/ocia)/2._dp &
           - H_theta(ocib-xnocc)* H_theta(ocia-xnocc)* &
             sqrt((1._dp-ocib)/(1._dp-ocia))/2._dp
   
end function dxLambda
!---------------------------------------------------------
function H_theta(x) ! Heavyside theta
!..Global variables
   use global;
   implicit none

   real(dp) :: H_theta

!..Argument
   real(dp) :: x

   if(x >= 0.d0) then
      H_theta = 1._dp 
   else
      H_theta = 0._dp
   endif
end function H_theta
!-------------------------------------------------------------------
!--------------------------------------------------------------
! The Piris Natural Orbital Functional (PNOF) 
! see IJQC 106, 1093 (2006); JCP 123, 214102 (2005);
! VERSION WITHOUT THE PINNED STATES CORRECTION
subroutine PNP(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)
      
!..Local
   integer :: ia,ib

!..External
   real(dp), external :: xLambdaP
   
   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = occ(ia,1)*occ(ib,2) &
                      + occ(ia,2)*occ(ib,1)
         fac_e(ia,ib) = 0._dp
         if(ia/=ib) then
            fac_c(ia,ib) = fac_c(ia,ib) + occ(ia,1)*occ(ib,1) &
                                        + occ(ia,2)*occ(ib,2)
            fac_e(ia,ib) =  - xLambdaP(occ(ia,1),occ(ib,1)) &
                            - xLambdaP(occ(ia,2),occ(ib,2))
            fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)
         endif
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine PNP
!---------------AND THE DERIVATIVE----------------------------
subroutine PNP_d(ia, ispin, occ, fac_c_d, fac_e_d)
! VERSION WITHOUT THE PINNED STATES CORRECTION
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

!..External
   real(dp), external :: dxLambdaP

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = min(occ(ia, ispin),1._dp-small)

   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op 
      fac_e_d(ib) = 0._dp  
      if(ib/=ia) then
         fac_c_d(ib) = fac_c_d(ib) + ocib
         fac_e_d(ib) = -dxLambdaP(ocia,ocib)
      endif
   enddo

end subroutine PNP_d
!---------------Auxiliary functions-------------------
function xLambdaP(ocia,ocib)
!..Global variables
   use global; use orbocc
   implicit none

!..The function
   real(dp) :: xLambdaP
!..Arguments
   real(dp) :: ocia,ocib

!..External
   real(dp), external :: H_theta

   xLambdaP = (1._dp - 2._dp*H_theta(0.5_dp-ocib)*H_theta(0.5_dp-ocia))* &
             sqrt(ocib*ocia) 
   
end function xLambdaP
!---------------------------------------------------------
function dxLambdaP(ocia,ocib)
!..Global variables
   use global; use orbocc
   implicit none

!..The function
   real(dp) :: dxLambdaP
!..Arguments
   real(dp) :: ocia,ocib

!..External
   real(dp), external :: H_theta

   dxLambdaP = (1._dp - 2._dp*H_theta(0.5_dp-ocib)*H_theta(0.5_dp-ocia))* &
             sqrt(ocib/ocia)/2._dp 
   
end function dxLambdaP
!------------------------------------------------------------------------

!--------------------------------------------------------------------
! THE MUELLER GENERAL-p FUNCTIONAL: for p=1/2 equivalent to BB0
! See AMK Mueller, Phys Lett 105A, 446 (1984).
! Modifying the FunP and FunP_d functions it can be used as a 
! general test functional. The function of n_i and n_j DOES not need
! to be symmetric in i and j! 
subroutine MPF(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

!..External
   real(dp), external :: funP

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -0.5_dp*(funP(occ(ia,1),occ(ib,1)) &
                               +funP(occ(ib,1),occ(ia,1))) &
                        -0.5_dp*(funP(occ(ia,2),occ(ib,2)) &
                               +funP(occ(ib,2),occ(ia,2)))

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine MPF
!---------------AND THE DERIVATIVE----------------------------
subroutine MPF_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

!..External
   real(dp), external :: funP_d

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = -0.5_dp*(funP_d(ocia,ocib,1)+funP_d(ocib,ocia,2))
   enddo

end subroutine MPF_d
! Auxiliary functions: 
function funP(ocia,ocib) !The function replacing the sq root
   use params_general, ONLY:dp
   use functional_m, ONLY:P_mull
   implicit none
!..The function
   real(dp) :: funP
!..Arguments
   real(dp) :: ocia,ocib

   funP = ocia**P_mull * ocib**(1._dp-P_mull)

end function funP

function funP_d(ocia,ocib,ind_d) !The derivatives of the function
   use params_general, ONLY:dp
   use functional_m, ONLY:P_mull
   implicit none
!..The function
   real(dp) :: funP_d
!..Arguments
   real(dp) :: ocia,ocib
   integer :: ind_d !see bellow:

   if (ind_d == 1) then ! derivative with respect to ocia
      funP_d = P_mull*ocia**(P_mull-1._dp) * ocib**(1._dp-P_mull)
   else ! derivative with respect to ocib
      funP_d = ocia**P_mull * (1._dp-P_mull)*ocib**(-P_mull)
   endif

end function funP_d
!------------------------------------------------------------------------

!--------------------------------------------------------------------
! Replace the square root with a PADE function 
! Marques, Lathiotakis PRA 2008
! f(x) = \frac{A_1 x + A_2 x^2}{1+ B_1 x} where x=n_a^\sigma n_b^\sigma
subroutine MLP(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

!..External
   real(dp), external :: FunFit

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = occ(ia,1)*occ(ib,2) &
                      + occ(ia,2)*occ(ib,1)
         fac_e(ia,ib) = 0._dp

! The folowing if block ads the extra terms which are excluded 
! if SIC is assumed. Thus for including SI terms comment out 
! the next line (and the endif).
         !if(ia/=ib) then
            fac_c(ia,ib) = fac_c(ia,ib) + occ(ia,1)*occ(ib,1) &
                                        + occ(ia,2)*occ(ib,2)
            fac_e(ia,ib) = -FunFit(occ(ia,1),occ(ib,1)) &
                           -FunFit(occ(ia,2),occ(ib,2))
            fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)
         !endif

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine MLP
!---------------AND THE DERIVATIVE----------------------------
subroutine MLP_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..External
   real(dp), external :: FunFit_d

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
!.....If SI terms are included remove the two lines bellow
!.....as well the 'if' and 'endif' lines bellow them.
      fac_c_d(ib) = ocib_op
      fac_e_d(ib) = 0._dp
      !if(ib/=ia) then
         fac_c_d(ib) = fac_c_d(ib) + ocib
         fac_e_d(ib) = -FunFit_d(ocia,ocib)*ocib
      !endif
   enddo

end subroutine MLP_d

function FunFit(ocia,ocib)
   use params_general, ONLY: dp
   use functional_m, ONLY: fit_a1, fit_a2, fit_b1
   implicit none
!..The function
   real(dp) :: FunFit
!..Arguments
   real(dp) :: ocia,ocib
!..Local
   real(dp) :: ff

   ff = ocia * ocib

   FunFit = (fit_a1*ff + fit_a2*ff*ff)/(1.0_dp + fit_b1*ff)

end function FunFit

function FunFit_d(ocia,ocib)
   use params_general, ONLY: dp
   use functional_m, ONLY: fit_a1, fit_a2, fit_b1
   implicit none
!..The function
   real(dp) :: FunFit_d
!..Arguments
   real(dp) :: ocia,ocib
!..Local
   real(dp) :: ff, ff1

   ff = ocia * ocib
   ff1 = (1.0_dp + fit_b1*ff)
   FunFit_d = ((fit_a1 + 2.0_dp*fit_a2*ff)*ff1 - fit_b1*(fit_a1*ff + fit_a2*ff*ff)) / &
           (ff1*ff1)

end function FunFit_d
!------------------------------------------------------------------------

!--------------------------------------------------------------------
! SA1 FUNCTIONAL (Sangeeta s new idea of modified Mueler functional
! with a beta to fix the sum rule).
subroutine SA1(occ, fac_h, fac_c, fac_e)
   use global; use functional_m
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

!   alpha=0.55d0

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -(occ(ia,1)*occ(ib,1))**alpha &
                        -(occ(ia,2)*occ(ib,2))**alpha

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*SA1_fac*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine SA1
!---------------AND THE DERIVATIVE----------------------------
subroutine SA1_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use functional_m
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = - alpha*SA1_fac*(ocib**alpha)*(ocia**(alpha-1._dp))
   enddo
! ACHTUNG!!! for the derivative SA1 is assumed constant. Thats for the
! sake of uniformity with other functionals SA1 is readjusted at the
! end of n_i minimization. 

end subroutine SA1_d
!--------------------------------------------------------------

!-------------------------------------------------------------
! AC3 (Journal of Chem Phys, 129, 164105 (2008).
subroutine AC3(occ, fac_h, fac_c, fac_e)
   use global; use functional_m
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib
   real(dp), external :: D_rohr
   real(dp) :: exch1, exch2, Drohr1, Drohr2

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         exch1 = sqrt(max(occ(ia,1)*occ(ib,1),zero))
         exch2 = sqrt(max(occ(ia,2)*occ(ib,2),zero))

!........BBC1 correction (PRB 75, 195120, (2007))
         if(ia/=ib) then
            if(ia>ibond(1).and.ib>ibond(1)) exch1 = -exch1 
            if(ia>ibond(2).and.ib>ibond(2)) exch2 = -exch2 
         endif

!........Multiply with AC3 factor
         Drohr1=D_rohr(ia,ib,1,occ(ia,1), occ(ib,1))
         Drohr2=D_rohr(ia,ib,2,occ(ia,2), occ(ib,2))
         exch1 = exch1*(1._dp-Drohr1)+occ(ia,1)*occ(ib,1)*Drohr1
         exch2 = exch2*(1._dp-Drohr2)+occ(ia,2)*occ(ib,2)*Drohr2

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = -0.5_dp*(exch1+exch2)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine AC3
!---------------AND THE DERIVATIVE----------------------------
subroutine AC3_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use functional_m
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op
   real(dp) :: f_ij
   real(dp) :: Drohr, Drohrder
   real(dp), external :: D_rohr, D_rohr_der

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = max(occ(ia, ispin),small)
   do ib = 1, nnatorb
      ocib= max(occ(ib,ispin),small)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib

!.....BBC1 factor
      f_ij=1._dp
      if(ia/=ib) then
         if(ia>ibond(ispin).and.ib>ibond(ispin)) f_ij=-1._dp
      endif

      Drohr=D_rohr(ia,ib,ispin,ocia, ocib)
      Drohrder=D_rohr_der(ia,ib,ispin,ocia, ocib)
      fac_e_d(ib) = -0.5_dp*f_ij*sqrt(max(ocib,zero))*(1._dp-Drohr)/sqrt(max(ocia,zero)) &
                  + f_ij*sqrt(max(ocia*ocib,zero))*Drohrder &
                  - ocib*Drohr - ocia*ocib*Drohrder
   enddo

end subroutine AC3_d
!----Auxiliary functions for AC3--------------------------------------
!--------D_ij---------------------------------------------------------
function D_rohr(i,j, ispin, occi, occj)
!..Global
   use global; use orbocc; use functional_m
   implicit none
   real(dp) :: D_rohr

!..Arguments
   integer, intent(in) :: i,j,ispin
   real(dp), intent(in)  :: occi, occj
   
!..Local
   real(dp) :: thi, thj, dij
   real(dp), external :: D_d, D_o

   if(ibond(ispin) >= i) then; thi=1._dp; else; thi=0._dp; endif !theta(ibond(ispin)-i)
   if(ibond(ispin) >= j) then; thj=1._dp; else; thj=0._dp; endif !theta(ibond(ispin)-j)
   if(i /= j) then; dij=0._dp; else; dij=1._dp; endif !delta_ij

   D_rohr = (1._dp - dij)*thi*thj &
          + dij*D_d(2._dp*occi-1._dp) &
          + (thi*(1._dp-thj)+thj*(1._dp-thi))*D_o(2._dp*(occi+occj-1._dp))

end function D_rohr

!--------D_ij derivative--------------------------------------------------------
function D_rohr_der(i,j, ispin, occi, occj)
!..Global
   use global; use orbocc; use functional_m
   implicit none
   real(dp) :: D_rohr_der

!..Arguments
   integer, intent(in) :: i,j,ispin
   real(dp), intent(in)  :: occi, occj
   
!..Local
   real(dp) :: thi, thj, dij
   real(dp), external :: D_d_der, D_o_der

   if(ibond(ispin) >= i) then; thi=1._dp; else; thi=0._dp; endif !theta(ibond(ispin)-i)
   if(ibond(ispin) >= j) then; thj=1._dp; else; thj=0._dp; endif !theta(ibond(ispin)-j)
   if(i /= j) then; dij=0._dp; else; dij=1._dp; endif !delta_ij

   D_rohr_der =  dij*2._dp*D_d_der(2._dp*occi-1._dp)  &
              + (thi*(1._dp-thj)+thj*(1._dp-thi))*2._dp*D_o_der(2._dp*(occi+occj-1._dp))
!..Factors of 2 in front of D_d_der, D_o_der are due to spin resolved version
end function D_rohr_der

!--------D_d--------------------------------------------------------
function D_d(x)
   use global; implicit none
   real(dp) :: D_d
   real(dp), intent(in) :: x
   real(dp) :: a_d, P_d, P_d2
   a_d = 1.4423_dp
   P_d = a_d * x*x * (x*x -2._dp); P_d2=P_d*P_d
   D_d = P_d2 / (1._dp + P_d2)
end function D_d
!--------D_o--------------------------------------------------------
function D_o(x)
   use global; implicit none
   real(dp) :: D_o
   real(dp), intent(in) :: x
   real(dp) :: a_o, P_o, P_o2
   a_o = 1.5552_dp
   P_o = a_o * x*x * (x*x -2._dp); P_o2=P_o*P_o
   D_o = P_o2 / (1._dp + P_o2)
end function D_o
!--------D_d-derivative-----------------------------------------
function D_d_der(x)
   use global; implicit none
   real(dp) :: D_d_der
   real(dp), intent(in) :: x
   real(dp) :: a_d, P_d, P_d2
   a_d = 1.4423_dp
   P_d = a_d * x*x * (x*x -2._dp); P_d2=P_d*P_d
   D_d_der = 8._dp*a_d*P_d*(x**3-x)/(1._dp+P_d2)**2
end function D_d_der
!--------D_o-derivative-----------------------------------------
function D_o_der(x)
   use global; implicit none
   real(dp) :: D_o_der
   real(dp), intent(in) :: x
   real(dp) :: a_o, P_o, P_o2
   a_o = 1.5552_dp
   P_o = a_o * x*x * (x*x -2._dp); P_o2=P_o*P_o
   D_o_der = 8._dp*a_o*P_o*(x**3-x)/(1._dp+P_o2)**2
end function D_o_der
!--------------------------------------------------------------

!--------------------------------------------------------------------
! Mueller for open shell with Necs new idea 
subroutine TST(occ, fac_h, fac_c, fac_e)!
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib
   real(dp) :: ocia_av, ocib_av, ocia_df, ocib_df

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      ocia_av=0.5_dp*(occ(ia,1) + occ(ia,2))
      ocia_df=0.5_dp*(occ(ia,1) - occ(ia,2))
      do ib=1,ia
         ocib_av=0.5_dp*(occ(ib,1) + occ(ib,2))
         ocib_df=0.5_dp*(occ(ib,1) - occ(ib,2))
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -(ocia_df*ocib_df + sqrt(ocia_av*ocib_av))

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo
end subroutine TST
!---------------AND THE DERIVATIVE----------------------------
subroutine TST_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op
   real(dp) :: ocia_av, ocib_av, ocib_df

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   ocia_av = 0.5_dp*(occ(ia, 1) + occ(ia, 2))
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_av = 0.5_dp*(occ(ib, 1) + occ(ib, 2))
      ocib_op=occ(ib, ispin_op)
      ocib_df = 0.5_dp*(occ(ib, 1) - occ(ib, 2))
      fac_c_d(ib) = ocib_op + ocib
      if(ispin == 1) then
         fac_e_d(ib) = -ocib_df - 0.5_dp*sqrt(ocib_av)/sqrt(ocia_av)
      else
         fac_e_d(ib) =  ocib_df - 0.5_dp*sqrt(ocib_av)/sqrt(ocia_av)
      endif
   enddo

end subroutine TST_d

!--------------------------------------------------------------------
! Entropy term functionals
function TS_nlnn()

!..Global variables
   use global; use orbocc
   implicit none

   real(dp) ::  TS_nlnn

!..Local variables
   integer :: ispin,ia
   real(dp) ::  xn, x1_n

   entropy=0._dp
   do ispin=1,2
      do ia=1,nnatorb
         xn=max(occnum(ia,ispin), small)
         x1_n=1._dp-xn
         entropy=entropy+xn*log(xn)+x1_n*log(x1_n)
      enddo
   enddo
   entropy=Boltz_konst*entropy
      
   TS_nlnn= -Temperat*entropy

end function TS_nlnn
!---------------AND THE DERIVATIVE----------------------------
function TS_nlnn_d(ia,ispin)
      
!..Global variables
   use global; use orbocc
   implicit none

   real(dp) :: TS_nlnn_d

!..Arguments
   integer :: ia, ispin

!..Local variables
   real(dp) :: xn, x1_n

   xn=max(occnum(ia,ispin), small)
   x1_n=1._dp-xn

   TS_nlnn_d=-Boltz_konst*Temperat*log(xn/x1_n)
      
end

!-------------------------------------------------------------------
! COULOMB and HF EXCHANGE FOR HYBRID FUNCTIONALS
subroutine DFT(occ, fac_h, fac_c, fac_e)
   use global; use DFT_par
   implicit none

!.....Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!.....Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -occ(ia,1)*occ(ib,1)-occ(ia,2)*occ(ib,2)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)*hyb_mix
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine DFT
!---------------AND THE DERIVATIVE----------------------------
subroutine DFT_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use DFT_par
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)
      
!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = -ocib*hyb_mix
   enddo

end subroutine DFT_d

!--------------------------------------------------------------------
! THE MUELLER FUNCTIONAL for open shell attempt 1
subroutine OSM(occ, fac_h, fac_c, fac_e)!
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = sqrt(occ(ia,1)*occ(ib,1))+occ(ia,1)*occ(ib,1)&
                        +sqrt(occ(ia,2)*occ(ib,2))+occ(ia,2)*occ(ib,2)&
                        +sqrt(occ(ia,1)*occ(ib,2))-occ(ia,1)*occ(ib,2)&
                        +sqrt(occ(ia,2)*occ(ib,1))-occ(ia,2)*occ(ib,1)

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = -0.25_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine OSM
!---------------AND THE DERIVATIVE----------------------------
subroutine OSM_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = - 0.25_dp*((sqrt(ocib)/sqrt(ocia) + 2._dp*ocib)&
                             +(sqrt(ocib_op)/sqrt(ocia) - 2._dp*ocib_op))
   enddo

end subroutine OSM_d

!--------------------------------------------------------------------
! HYBRID Power/DFT-LDA
subroutine HDR(occ, fac_h, fac_c, fac_e)!
   use global; use functional_m, ONLY: alpha, xmix_RD
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fac_e(ia,ib) = -(occ(ia,1)*occ(ib,1))**alpha &
                        -(occ(ia,2)*occ(ib,2))**alpha

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)*xmix_RD

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine HDR
!---------------AND THE DERIVATIVE----------------------------
subroutine HDR_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use functional_m, ONLY: alpha, xmix_RD
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_c_d(ib) = ocib_op + ocib
      fac_e_d(ib) = - xmix_RD*alpha*(ocib**alpha)*(ocia**(alpha-1._dp))
   enddo

end subroutine HDR_d

!--------------------------------------------------------------------
subroutine DBF(occ, fac_h, fac_c, fac_e)!
   use global; use functional_m, only: DBF_fun_RDM
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

   select case (DBF_fun_RDM)
      case("BB0") 
         call BB0(occ, fac_h, fac_c, fac_e)
      case("BB3") 
         call BB3(occ, fac_h, fac_c, fac_e)
      case default
         print*,'Functional DBF:', DBF_fun_RDM, ' UKNOWN'
         stop 
   end select 

end subroutine DBF

!---------------AND THE DERIVATIVE----------------------------
subroutine DBF_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use functional_m, only: DBF_fun_RDM
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

   select case (DBF_fun_RDM)
      case("BB0") 
         call BB0_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case("BB3") 
         call BB3_d(ia, ispin, occ, fac_c_d, fac_e_d)
      case default
         print*,'Functional DBF: ',DBF_fun_RDM, ' UKNOWN'
         stop 
   end select 

end subroutine DBF_d

!--------------------------------------------------------------------
! Hybrid Mueller - HF: f = a f_HF + (1-a) f_BB0
subroutine NEK(occ, fac_h, fac_c, fac_e)
   use global; use functional_m, only: alpha
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb) 
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local
   real(dp) :: fac_h_HF(lnnatorb), fac_c_HF(lnnatorb,lnnatorb), fac_e_HF(lnnatorb,lnnatorb)

   call RHF(occ, fac_h_HF, fac_c_HF, fac_e_HF)

   call BB0(occ, fac_h, fac_c, fac_e)

   fac_h = alpha*fac_h_HF + (1._dp -alpha)*fac_h
   fac_c = alpha*fac_c_HF + (1._dp -alpha)*fac_c
   fac_e = alpha*fac_e_HF + (1._dp -alpha)*fac_e

end subroutine NEK

!---------------AND THE DERIVATIVE----------------------------
subroutine NEK_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use functional_m, only: alpha
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   real(dp) :: fac_c_d_HF(lnnatorb), fac_e_d_HF(lnnatorb)
 
   call RHF_d(ia, ispin, occ, fac_c_d_HF, fac_e_d_HF)

   call BB0_d(ia, ispin, occ, fac_c_d, fac_e_d)

   fac_c_d= alpha*fac_c_d_HF + (1._dp -alpha)*fac_c_d
   fac_e_d= alpha*fac_e_d_HF + (1._dp -alpha)*fac_e_d

end subroutine NEK_d

!--------------------------------------------------------------------
! THE FUNCTIONAL of H. Shull and P. Lo¨wdin, J. Chem. Phys. 30, 617 (1959).
! See also W. Kutzelnigg, Theor. Chim. Acta 1, 327 (1963).
! Exact for 2 electron systems
subroutine LSH(occ, fac_h, fac_c, fac_e)!
   use global; use functional_m, ONLY:f_LSH
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib

   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
         fac_c(ia,ib) = 0._dp
         fac_e(ia,ib) =  f_LSH(ia,1)*f_LSH(ib,1)*sqrt(occ(ia,1)*occ(ib,1)) &
                        +f_LSH(ia,2)*f_LSH(ib,2)*sqrt(occ(ia,2)*occ(ib,2))

!........The double counting 1/2 factor
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine LSH
!---------------AND THE DERIVATIVE----------------------------
subroutine LSH_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use functional_m, ONLY:f_LSH
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      fac_c_d(ib) = 0._dp
      fac_e_d(ib) = 0.5_dp*f_LSH(ia,ispin)*f_LSH(ib,ispin)*sqrt(ocib)/sqrt(ocia)
   enddo

end subroutine LSH_d

!--------------------------------------------------------------------
! PNOF 5
subroutine PN5(occ, fac_h, fac_c, fac_e)!
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib, iap, Nmax

   Nmax=2*nele(1)

   fac_h = 0._dp; fac_c=0._dp; fac_e=0._dp
   do ia=1,Nmax
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      iap=Nmax-ia+1
      do ib=1,ia
         if ( ia /= ib .and. ia+ib /= Nmax+1 ) then
            fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
            fac_e(ia,ib) = -occ(ia,1)*occ(ib,1)-occ(ia,2)*occ(ib,2)

!........The double counting 1/2 factor
            fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
            fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

            fac_c(ib,ia) = fac_c(ia,ib)
            fac_e(ib,ia) = fac_e(ia,ib)
         endif
      enddo
      fac_c(ia,ia) = 0.5_dp*(occ(ia,1)+occ(ia,2))
      fac_e(ia,ia) = 0._dp
      fac_c(ia,iap) = 0._dp; fac_c(iap,ia)=0._dp
      fac_e(ia,iap) = -0.5_dp*(sqrt(occ(ia,1)*occ(iap,1))+sqrt(occ(ia,2)*occ(iap,2)))
      fac_e(iap,ia) = fac_e(ia,iap)
   enddo

end subroutine PN5
!---------------AND THE DERIVATIVE----------------------------
subroutine PN5_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op, Nmax, iap
   real(dp) :: ocia, ocib, ocib_op, ociap

   Nmax=2*nele(1)
   iap=Nmax-ia+1

   ispin_op=2
   if(ispin==2) ispin_op=1

   fac_e_d=0._dp; fac_c_d=0._dp

   ocia = max(occ(ia, ispin),zero)
   do ib = 1, Nmax
      if ( ia /= ib .and. ia+ib /= Nmax+1 ) then
         ocib= occ(ib,ispin)
         ocib_op=occ(ib, ispin_op)
         fac_c_d(ib) = ocib_op + ocib
         fac_e_d(ib) = -ocib
      endif
   enddo
   fac_c_d(ia) = 0.5_dp
   ociap=occ(iap, ispin)
   fac_e_d(iap) = -0.5_dp*sqrt(ociap/ocia)

end subroutine PN5_d

!--------------------------------------------------------------------
! Modification of the Mueller functional where in the Direct Coulomb term 
! the dependence on n_i is dropped.
subroutine BBN(occ, fac_h, fac_c, fac_e)!
   use global
   implicit none

!..Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!..Local variables
   integer :: ia, ib
   real(dp) :: fca1,fca2,fcb1,fcb2


   do ia=1,nnatorb
      fac_h(ia) = occ(ia,1) + occ(ia,2)
      do ib=1,ia
!        fac_c(ia,ib) = (occ(ia,1)+occ(ia,2))*(occ(ib,1)+occ(ib,2))
         fca1=0._dp; if ( occ(ia,1) > 0.5_dp ) fca1=1._dp 
         fca2=0._dp; if ( occ(ia,2) > 0.5_dp ) fca2=1._dp 
         fcb1=0._dp; if ( occ(ib,1) > 0.5_dp ) fcb1=1._dp 
         fcb2=0._dp; if ( occ(ib,2) > 0.5_dp ) fcb2=1._dp 
         fac_c(ia,ib) = (fca1+fca2)*(fcb1+fcb2)
         fac_e(ia,ib) = -sqrt(occ(ia,1)*occ(ib,1)) &
                        -sqrt(occ(ia,2)*occ(ib,2))

!........The double counting 1/2 factor
         fac_c(ia,ib) = 0.5_dp*fac_c(ia,ib)
         fac_e(ia,ib) = 0.5_dp*fac_e(ia,ib)

         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

end subroutine BBN
!---------------AND THE DERIVATIVE----------------------------
subroutine BBN_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)

!..Local
   integer :: ib, ispin_op
   real(dp) :: ocia, ocib, ocib_op

   ispin_op=2
   if(ispin==2) ispin_op=1

   ocia = occ(ia, ispin)
   do ib = 1, nnatorb
      ocib= occ(ib,ispin)
      ocib_op=occ(ib, ispin_op)
      fac_e_d(ib) = - 0.5_dp*sqrt(ocib)/sqrt(ocia)
      if(ocib>0.5_dp) then
         ocib=1._dp
      else
         ocib=0._dp
      endif 
      if(ocib_op>0.5_dp) then
         ocib_op=1._dp
      else
         ocib_op=0._dp
      endif 
      if(ocia>0.5_dp) then
         ocia=1._dp
      else
         ocia=0._dp
      endif 
     
      fac_c_d(ib) = ocib_op + ocib
   enddo

end subroutine BBN_d


!-------------------------------------------------------------------
! Modified LSH it treats the inner electrons with HF and the two outer ones with LSH
! The mixed terms interaction is treated with HF
subroutine LSM(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!.....Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!.....Local variables
   integer :: ia, ib
   real(dp) :: fac_e_LSH(lnnatorb,lnnatorb), fac_e_RHF(lnnatorb,lnnatorb), fac_c_RHF(lnnatorb,lnnatorb)
    
   call RHF(occ, fac_h, fac_c_RHF, fac_e_RHF)
   call LSH(occ, fac_h, fac_c, fac_e_LSH)
   fac_c =0._dp
   fac_e =0._dp

   do ia=1, nele(1)-1 !for closed shell only
      do ib=1, nnatorb
!       do ib=1, ia
         fac_c(ia,ib) = fac_c_RHF(ia, ib)
         fac_e(ia,ib) = fac_e_RHF(ia, ib)
         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

   do ia = nele(1), nnatorb
     do ib = ia, nnatorb
       fac_e(ia,ib) = fac_e_LSH(ia, ib)
       fac_e(ib,ia) = fac_e(ia,ib)
     end do
   end do

end subroutine LSM
!---------------AND THE DERIVATIVE----------------------------
subroutine LSM_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)


!.. Local
   integer :: ib
   real(dp) ::  fac_c_d_RHF(lnnatorb), fac_e_d_RHF(lnnatorb), fac_e_d_LSH(lnnatorb)

   call RHF_d(ia, ispin, occ, fac_c_d_RHF, fac_e_d_RHF)
   call LSH_d(ia, ispin, occ, fac_c_d, fac_e_d_LSH)
   
   fac_c_d = 0._dp
   fac_e_d = 0._dp

   if (ia.lt.nele(1) ) then 
     do ib = 1, nnatorb  
       fac_c_d(ib) =  fac_c_d_RHF(ib)
       fac_e_d(ib) =  fac_e_d_RHF(ib)
     end do
   else
     do ib = 1, nele(1)-1 
       fac_c_d(ib) =  fac_c_d_RHF(ib)
       fac_e_d(ib) =  fac_e_d_RHF(ib)
     end do
     do ib = nele(1), nnatorb 
       fac_e_d(ib) =  fac_e_d_LSH(ib)
     end do
   end if

end subroutine LSM_d
!-------------------------------------------------------------------
! Modified LSH for triplets it treats the inner electrons with HF and the two outer ones with LSH
! The mixed terms interaction is treated with HF
subroutine LST(occ, fac_h, fac_c, fac_e)
   use global
   implicit none

!.....Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!.....Local variables
   integer :: ia, ib
   real(dp) :: fac_e_LSH(lnnatorb,lnnatorb), fac_e_RHF(lnnatorb,lnnatorb), fac_c_RHF(lnnatorb,lnnatorb)
    
   call RHF(occ, fac_h, fac_c_RHF, fac_e_RHF)
   call LSH(occ, fac_h, fac_c, fac_e_LSH)
   fac_c =0._dp
   fac_e =0._dp

   do ia=1, nele(2) !for open shell triplets only
      do ib=1, nnatorb
         fac_c(ia,ib) = fac_c_RHF(ia, ib)
         fac_e(ia,ib) = fac_e_RHF(ia, ib)
         fac_c(ib,ia) = fac_c(ia,ib)
         fac_e(ib,ia) = fac_e(ia,ib)
      enddo
   enddo

   do ia = nele(2)+1, nnatorb
     do ib = ia, nnatorb
       fac_e(ia,ib) = fac_e_LSH(ia, ib)
       fac_e(ib,ia) = fac_e(ia,ib)
     end do
   end do

end subroutine LST
!---------------AND THE DERIVATIVE----------------------------
subroutine LST_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)


!.. Local
   integer :: ib
   real(dp) ::  fac_c_d_RHF(lnnatorb), fac_e_d_RHF(lnnatorb), fac_e_d_LSH(lnnatorb)

   call RHF_d(ia, ispin, occ, fac_c_d_RHF, fac_e_d_RHF)
   call LSH_d(ia, ispin, occ, fac_c_d, fac_e_d_LSH)
   
   fac_c_d = 0._dp
   fac_e_d = 0._dp

   if (ia.lt.(nele(2) + 1)) then 
     do ib = 1, nnatorb  
       fac_c_d(ib) =  fac_c_d_RHF(ib)
       fac_e_d(ib) =  fac_e_d_RHF(ib)
     end do
   else
     do ib = 1, nele(2) 
       fac_c_d(ib) =  fac_c_d_RHF(ib)
       fac_e_d(ib) =  fac_e_d_RHF(ib)
     end do
     do ib = nele(2)+1, nnatorb 
       fac_e_d(ib) =  fac_e_d_LSH(ib)
     end do
   end if

end subroutine LST_d
!Functional for 3 electrons with spin 1/2 based on the wavefunction expansion in 3 determinants
subroutine GPC(occ, fac_h, fac_c, fac_e)
   use global
   use functional_m, only:l_non_JK; use orbocc, only:f_non_JK, fc12, fc23, f_GPC_13
   implicit none

!.....Arguments
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_h(lnnatorb)
   real(dp), intent(out) :: fac_c(lnnatorb,lnnatorb)
   real(dp), intent(out) :: fac_e(lnnatorb,lnnatorb)

!.....Local variables
   real(dp) ::  fac_e_RHF(lnnatorb,lnnatorb), fac_c_RHF(lnnatorb,lnnatorb)
   

   call RHF(occ, fac_h, fac_c_RHF, fac_e_RHF)
   fac_c = fac_c_RHF
   fac_e = fac_e_RHF

         fac_c(1,1) = fac_c(1,1) + occ(1,2)*occ(3,2)
         fac_c(2,2) = fac_c(2,2) - occ(2,1)*occ(2,2)
         fac_c(1,3) = fac_c(1,3) - 0.5*(occ(1,2)*occ(3,1)-occ(3,2)*occ(1,1)-2*occ(1,2)*occ(3,2))
         fac_c(3,1) = fac_c(1,3) 
         fac_c(3,3) = fac_c(3,3) + occ(1,2)*occ(3,2)
         fac_e(1,2) = fac_e(1,2) + 0.5*occ(2,1)*occ(2,2)
         fac_e(2,1) = fac_e(1,2) 
         fac_e(1,3) = fac_e(1,3) + occ(1,2)*occ(3,2)+f_GPC_13*sqrt(occ(1,2)*occ(3,2))
         fac_e(3,1) = fac_e(1,3) 
         fac_e(2,3) = fac_e(2,3) + 0.5*occ(2,1)*occ(2,2)
         fac_e(3,2) = fac_e(3,2)

   if ( l_non_JK ) then
      f_non_JK = 2._dp*(fc12*sqrt(occ(1,2)*occ(2,2)) + fc23*sqrt(occ(2,2)*occ(3,2)))
   endif
end subroutine GPC
!---------------AND THE DERIVATIVE----------------------------
subroutine GPC_d(ia, ispin, occ, fac_c_d, fac_e_d)
   use global; use orbocc, only:fd_non_JK1, fd_non_JK2, fd_non_JK3, fc12, fc23, f_GPC_13
   use functional_m, only:l_non_JK
   implicit none

!..Arguments
   integer, intent(in) :: ia, ispin 
   real(dp), intent(in) :: occ(lnnatorb,3)
   real(dp), intent(out) :: fac_c_d(lnnatorb), fac_e_d(lnnatorb)


!.. Local
   integer :: ib
   real(dp) ::  fac_c_d_RHF(lnnatorb), fac_e_d_RHF(lnnatorb), fac_e_d_LSH(lnnatorb)

   call RHF_d(ia, ispin, occ, fac_c_d_RHF, fac_e_d_RHF)
   fac_c_d = fac_c_d_RHF
   fac_e_d = fac_e_d_RHF
   
  if ((ia == 1).and.( ispin==2))  then
    fac_c_d(1) = fac_c_d(1) + occ(3,2)
    fac_c_d(3) = fac_c_d(3) -0.5*occ(3,1) - occ(3,2)
    fac_e_d(3) = fac_e_d(3) + occ(3,2) + 0.5*f_GPC_13*sqrt(occ(3,2)/occ(1,2))
  end if

  if ((ia == 3).and.( ispin==2)) then
    fac_c_d(1) = fac_c_d(1) - 0.5*occ(1,1) -occ(1,2)
    fac_c_d(3) = fac_c_d(3) - occ(1,2)
    fac_e_d(1) = fac_e_d(1) + occ(1,2) +0.5*f_GPC_13*sqrt(occ(1,2)/occ(3,2))
  end if 

  if ((ia == 1).and.(ispin == 1)) then
   fac_c_d(3)=fac_c_d(3)-0.5*occ(3,2)
  end if

  if ((ia == 2).and.(ispin == 2)) then
   fac_c_d(2) =fac_c_d(2) - occ(2,1)
   fac_e_d(1) =fac_e_d(1) + 0.5*occ(2,1)
   fac_e_d(3) =fac_e_d(3) + 0.5*occ(2,1)
  end if

  if ((ia == 2).and.(ispin == 1)) then
   fac_c_d(2) =fac_c_d(2) - occ(2,2)
   fac_e_d(1) =fac_e_d(1) + 0.5*occ(2,2)
   fac_e_d(3) =fac_e_d(3) + 0.5*occ(2,2)
  end if


  if ((ia ==3).and.(ispin == 1)) then
    fac_c_d(1) = fac_c_d(1)-0.5*occ(1,2)
  end if

  if ( l_non_JK .and. ispin == 2 ) then
     fd_non_JK1 = fc12*sqrt(occ(2,2)/max(small,occ(1,2))) 
     fd_non_JK2 = fc12*sqrt(occ(1,2)/max(small,occ(2,2))) + fc23*sqrt(occ(3,2)/max(small,occ(2,2)))
     fd_non_JK3 = fc23*sqrt(occ(2,2)/max(small,occ(3,2)))
  endif

end subroutine GPC_d
