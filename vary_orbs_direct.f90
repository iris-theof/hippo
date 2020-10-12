subroutine vary_orbs_direct()
!-----------------------------------------------------------------------
! Optimizes the orbitals by minimizing the original energy functional
! under the orthonormality constraint.
! A brute force algorithm is implemented. The gradient direction
! is chosen for every orbital and at each step an attempt to lower the
! energy is made. The attempt determines the step size. The step
! does not lead to the minimum in the particular direction. It just 
! has to lower the energy more than 'thresh'. The step size is 
! adjusted in every step with a scale_up and scale_down factor
! This procedure is very slow, with slow convergence but is 
! safe. 
!
!         Created by N. N. Lathiotakis and I.Theophilou
!
!-----------------------------------------------------------------------

   use global; use matrices; use orbocc; use functional_m; use energies
   implicit none

!..Local variables
   integer :: ia, itry, ntry, istep, nstep
   real(dp) :: vec_old(lnbasis,lnbasis), E_deriv(lnbasis,lnnatorb)
   real(dp) :: h(lnbasis,lnnatorb)
   real(dp), allocatable :: F_old(:,:,:)
   real(dp) :: old_energy, old_energy1, Delta_E
   real(dp) :: small_step, thresh
   integer :: iconstr, icnt, natlim
   logical :: test_orth

   print*,'In vary_orbs_direct'

!..Local array allocation
   allocate ( F_old(lnbasis,lnbasis,lnnatorb) )

   iconstr=0
   icnt=0
!..Parameters
!  natlim = min(max(nele(1),nele(2))+5,nbasis-1)
   natlim = min(nnatorb,nbasis-1)
   ntry=ntry_dir
   nstep=nst_orb_dir
   small_step=1.d-8
   thresh = 1.d-16 !threshold for succesful step

!..Initialize steps 
!  if(iscf == 1.or.mod(iscf,50) == 0) then 
!     step_dir= step_dir_init
!  endif

!..Give failed variation a new try
   do ia=1,nnatorb
      if(abs(step_dir(ia)) < small_step) &
             step_dir(ia) = sqrt(small_step*step_dir_init)
   enddo

!..save old values
   old_energy = Tot_energy
   vec_old=vecnat
   F_old=F

   do istep=1,nstep
!........Save old orbitals/ Fock matrices:
      old_energy1=old_energy
      do ia=1,natlim 
         if(abs(step_dir(ia)) >= small_step) then

!...........Calculate the ia-th gradient
            call calc_Lagrange_mult(ia)
            call E_deriv_phi(ia, E_deriv)
    
            h=E_deriv 

!...........Normalize the ia-th gradient
            call orthog_h_ia(1,ia,h)

            do itry=1,ntry

!..............Change orbitals along the conjugate gradient by step:
               vecnat(:,ia) = vec_old(:,ia)-h(:,ia)*step_dir(ia)

!..............Orthogonalize new orbitals
               call orthog_vec_ia(ia)

!.................Calculate Fock/Exch/Coul/Energies
               iconstr=iconstr+1
               icnt = icnt + (nnatorb - ia + 1)
               call construct_f(0)
               call total_energy_MO(1)

!.................Check if Total Energy is lower
               Delta_E=Tot_energy-old_energy
               if(Delta_E < -thresh) then 
!.................Successful step
                  step_dir(ia)=step_dir(ia)*scale_up !Increase step
                  old_energy=Tot_energy
                  vec_old(:,ia:nbasis) = vecnat(:,ia:nbasis)
                  F_old=F
                  exit
               else 
!.................Not Successful step
                  step_dir(ia)=-step_dir(ia)*scale_down !Decrease step size
                  if(abs(step_dir(ia)) <= small_step.or.itry == ntry) then
!....................Successful step not possible 
                     exit
                  endif !(abs(step_dir(ia)) <= small_step.or.itry == ntry)
               endif !(Delta_E < -thresh)
            enddo !itry=1,ntry
            Tot_energy=old_energy
            vecnat(:,ia:nbasis) = vec_old(:,ia:nbasis)
            F = F_old
  
         endif !(step_dir(ia) >= small_step
      enddo !ia
      write(6,111)Istep, old_energy, old_energy-old_energy1  
      print*,'Nu of Fock calculations:',iconstr, icnt
   enddo !istep

111 format(' Istep:',i3,' Energy:',f20.12,' Diff: ',e20.12) 
  
   call normal_test(test_orth, small, 1)

   deallocate ( F_old )


end subroutine vary_orbs_direct

!-----------------------------------------------------------------------
   subroutine orthog_vec_ia(ia)
! Othogonalization of occupied / unoccupied orbitals.
! for ib >= ia: orthonormalize ib to all with smaller index

!..Global
   use global; use matrices; use orbocc
   implicit none

!..Arguments
   integer, intent(in) :: ia

!..Local variables
   integer :: ib, ic
   real(dp) :: prod, xnorm

!..Remove from ia orbital the projection of ib < ia orbitals 
   do ib = 1, ia-1
      prod = dot_product(vecnat(:,ia), matmul(ovlap, vecnat(:,ib)))
      vecnat(:,ia) = vecnat(:,ia) - prod *vecnat(:,ib)
   enddo

!..Normalize vecnat(ia)
   xnorm = 1.d0/sqrt(dot_product(vecnat(:,ia),matmul(ovlap, vecnat(:,ia))))
   vecnat(:,ia) = xnorm*vecnat(:,ia)

!..For orbs with index ib > ia remove the projection of those 
!..ic: ia <= ic < ib
   do ib = ia+1, nbasis
      do ic = ia,ib-1
         prod = dot_product(vecnat(:,ic), matmul(ovlap, vecnat(:,ib)))
         vecnat(:,ib) = vecnat(:,ib) - prod * vecnat(:,ic)
      enddo !ic

!.....Normalize ib orbital
      xnorm = 1.d0/sqrt(dot_product(vecnat(:,ib),matmul(ovlap,vecnat(:,ib))))
      vecnat(:,ib) = xnorm*vecnat(:,ib) 
   enddo !ib = ia+1, nbasis

   return
   end
!-----------------------------------------------------------------------

subroutine orthog_h_ia(mode,ia,h)
!-----------------------------------------------------------------------
! Normalize h and orthogonalize to the orbital with index ia
! mode=0 only normalize h
! mode>0 also orthogonalize
!-----------------------------------------------------------------------

!..Global
   use global; use matrices; use orbocc
   implicit none

!..Arguments
   integer, intent(in) :: mode,ia
   real(dp), intent(out) :: h(lnbasis,lnnatorb)

!..Local variables
   integer :: ib
   real(dp) :: prod, xnorm

!..Remove from ia-th gradient the projection of ib=ia orbital 
   if(mode > 0) then
      do ib = ia, ia !alternatively ib=ia,nnatorb
         prod = dot_product(vecnat(:,ib),matmul(ovlap,h(:,ia)))
         h(:,ia) = h(:,ia) - prod * vecnat(:,ib)
      enddo
   endif

!..Now normalize h(ia)
   xnorm = dot_product(h(:,ia),matmul(ovlap,h(:,ia)))
   if (xnorm > zero) then 
      xnorm = 1.d0/sqrt(xnorm)
      h(:,ia) = xnorm*h(:,ia)
   endif
end subroutine orthog_h_ia
!--------------------------------------------------------------------

subroutine calc_Lagrange_mult_all() !ALL LAGRANGE MULTIPLIERS
!--------------------------------------------------------------------
! Calculates ALL the quantities e_{ab} = <a|F^a|b> = <b|F^a|a> 
! Note! e_{ab} = <a|F^b|b> = <b|F^b|a> i.e. the second index
! of e corresponds to the index of the F matrix involved.
!--------------------------------------------------------------------

!..Global
   use global; use matrices; use orbocc; use functional_m, only:l_non_JK, maxorb
   implicit none

!..Local variables
   integer :: ia,ib
   real(dp) :: X_int

  F_xc=F_xc+(F_ha-F_haM) -F_x

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ia,ib)
!$OMP DO
   do ib = 1,nbasis
      do ia = 1,nbasis
         e_ab(ia,ib) = 0.d0
         ee_ab(ia,ib) = 0.d0
         if(ib <= nnatorb) then
            e_ab(ia,ib) = dot_product(vecnat(:,ia),matmul(F(:,:,ib),vecnat(:,ib)))
            ee_ab(ia,ib) = dot_product(vecnat(:,ia),matmul(F_xc(:,:,ib),vecnat(:,ib)))
!           do k=1,nbasis
!              do l=1,nbasis
!                 e_ab(ia,ib) = e_ab(ia,ib) + vecnat(k,ia)*F(k,l,ib)*vecnat(l,ib)
!              enddo
!           enddo
         endif !(ib <= nnatorb) 
      enddo !ib
   enddo !ia
!$OMP END DO
!$OMP END PARALLEL

   if ( l_non_JK ) then
      do ib = 1,maxorb
         call psi4_int(vecnat(:,ib),vecnat(:,2),vecnat(:,2),vecnat(:,3),X_int)
         e_ab(ib,1) = e_ab(ib,1)+f_non_JK*X_int
         call psi4_int(vecnat(:,ib),vecnat(:,1),vecnat(:,2),vecnat(:,3),X_int)
         e_ab(ib,2) = e_ab(ib,2)+f_non_JK*X_int
         call psi4_int(vecnat(:,ib),vecnat(:,3),vecnat(:,1),vecnat(:,2),X_int)
         e_ab(ib,2) = e_ab(ib,2)+f_non_JK*X_int
         call psi4_int(vecnat(:,ib),vecnat(:,2),vecnat(:,2),vecnat(:,1),X_int)
         e_ab(ib,3) = e_ab(ib,3)+f_non_JK*X_int
      enddo
   endif

   return
   end
!--------------------------------------------------------------------

subroutine calc_Lagrange_mult(iaa) !ONLY E_AB INVOLVING iaa
!--------------------------------------------------------------------
! Calculates only the quantities e_{ba} = <a|F^a|b> = <b|F^a|a>
! where the either of a,b is equal to the argument index iaa.
! Note! e_{ab} = <a|F^b|b> = <b|F^b|a> i.e. the second index
! of e corresponds to the index of the F matrix involved.
!--------------------------------------------------------------------

!..Global
   use global; use matrices; use orbocc
   implicit none

!..Arguments
   integer, intent(in) :: iaa

!..Local variables
   integer :: ib
   real(dp) :: X_int

   do ib = 1,nbasis
      e_ab(ib,iaa) = 0.d0
      if(iaa <= nnatorb) then
         e_ab(ib,iaa) = dot_product(vecnat(:,ib),matmul(F(:,:,iaa),vecnat(:,iaa)))
      endif !(iaa <= nnatorb) 
   enddo !ib

   do ib = 1,nbasis
      e_ab(iaa,ib) = 0.d0
      if(ib <= nnatorb) then
         e_ab(iaa,ib) = dot_product(vecnat(:,iaa),matmul(F(:,:,ib),vecnat(:,ib)))
      endif !(ib <= nnatorb) 
   enddo !ib

!  if ( l_non_JK ) then
!     do ib = 1,maxorb
!        if ( iaa == 1) then
!        call psi4_int(vecnat(:,ib),vecnat(:,2),vecnat(:,2),vecnat(:,3),X_int)
!        e_ab(ib,1) = e_ab(ib,1)+f_non_JK*X_int
!        elseif( iaa == 2) then
!        call psi4_int(vecnat(:,ib),vecnat(:,1),vecnat(:,2),vecnat(:,3),X_int)
!        e_ab(ib,2) = e_ab(ib,2)+f_non_JK*X_int
!        call psi4_int(vecnat(:,ib),vecnat(:,3),vecnat(:,1),vecnat(:,2),X_int)
!        e_ab(ib,2) = e_ab(ib,2)+f_non_JK*X_int
!        elseif ( iaa==3) then
!        call psi4_int(vecnat(:,ib),vecnat(:,2),vecnat(:,2),vecnat(:,1),X_int)
!        e_ab(ib,3) = e_ab(ib,3)+f_non_JK*X_int
!        endif
!     enddo
!  endif


!  if ( l_non_JK ) then
!     do ib = 1,nbasis
!        if ( iaa == 1) then
!           call psi4_int(vecnat(:,ib),vecnat(:,2),vecnat(:,2),vecnat(:,3),X_int)
!           e_ab(ib,iaa) = e_ab(ib,iaa)+f_non_JK*X_int
!        elseif ( iaa == 2 ) then
!           call psi4_int(vecnat(:,ib),vecnat(:,1),vecnat(:,2),vecnat(:,3),X_int)
!           e_ab(ib,iaa) = e_ab(ib,iaa)+f_non_JK*X_int
!           call psi4_int(vecnat(:,ib),vecnat(:,3),vecnat(:,1),vecnat(:,2),X_int)
!           e_ab(ib,iaa) = e_ab(ib,iaa)+f_non_JK*X_int
!        elseif ( iaa == 3 ) then
!           call psi4_int(vecnat(:,ib),vecnat(:,2),vecnat(:,2),vecnat(:,1),X_int)
!           e_ab(ib,iaa) = e_ab(ib,iaa)+f_non_JK*X_int
!        endif
!     enddo
!     if (iaa == 1) then
!        call psi4_int(vecnat(:,1),vecnat(:,1),vecnat(:,2),vecnat(:,3),X_int)
!        e_ab(1,2) = e_ab(1,2) + f_non_JK*X_int
!        call psi4_int(vecnat(:,1),vecnat(:,3),vecnat(:,1),vecnat(:,2),X_int)
!        e_ab(1,2) = e_ab(1,2) + f_non_JK*X_int
!        call psi4_int(vecnat(:,1),vecnat(:,2),vecnat(:,2),vecnat(:,1),X_int)
!        e_ab(1,3) = e_ab(1,3) + f_non_JK*X_int
!     elseif ( iaa == 2) then
!        call psi4_int(vecnat(:,2),vecnat(:,2),vecnat(:,2),vecnat(:,3),X_int)
!        e_ab(2,1) = e_ab(2,1) + f_non_JK*X_int
!        call psi4_int(vecnat(:,2),vecnat(:,2),vecnat(:,2),vecnat(:,1),X_int)
!        e_ab(2,3) = e_ab(2,3) + f_non_JK*X_int
!     elseif ( iaa == 3) then
!        call psi4_int(vecnat(:,3),vecnat(:,2),vecnat(:,2),vecnat(:,3),X_int)
!        e_ab(3,1) = e_ab(3,1) + f_non_JK*X_int
!        call psi4_int(vecnat(:,3),vecnat(:,1),vecnat(:,2),vecnat(:,3),X_int)
!        e_ab(3,2) = e_ab(3,2) + f_non_JK*X_int
!        call psi4_int(vecnat(:,3),vecnat(:,3),vecnat(:,1),vecnat(:,2),X_int)
!        e_ab(3,2) = e_ab(3,2) + f_non_JK*X_int
!     endif
!  endif

end subroutine calc_Lagrange_mult
!---------------------------------------------------------------------------
  
subroutine  E_deriv_phi(ia, E_deriv)
!---------------------------------------------------------------------------
! Calculates the gradient with respect to the ia-orbital:
!---------------------------------------------------------------------------

!..Global
   use global; use matrices; use orbocc; use functional_m, only:l_non_JK
   implicit none

!..Arguments
   integer,  intent(in) :: ia
   real(dp),  intent(out) :: E_deriv(lnbasis,lnnatorb)

!..Local variables
   real(dp) :: weight
   integer :: i,j

   weight=0.5d0

   do i=1,nbasis
      E_deriv(i,ia) = dot_product(F(i,:,ia), vecnat(:,ia)) 
!.....Now add the gradient correction:
!     (Eq. (12) Chem. Phys. Lett. 364, 409 (2002); PRA 55, 1765 (1997))
      do j=1,nbasis
         E_deriv(i,ia)=E_deriv(i,ia)-weight*ovlap(i,j)*&
             dot_product(e_ab(ia,:)+e_ab(:,ia),vecnat(j,:))
      enddo
   enddo

   if (l_non_JK) then
      do i=1,nbasis
         E_deriv(i,1) =E_deriv(i,1)+XI1(i)
         E_deriv(i,2) =E_deriv(i,2)+XI2(i)
         E_deriv(i,3) =E_deriv(i,3)+XI3(i)
      enddo
   endif
!  print*,ia,(E_deriv(i,ia),i=1,nbasis)

end subroutine E_deriv_phi
!---------------------------------------------------------------------------
