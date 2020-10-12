subroutine readgamess()
!----------------------------------------------------------------------
! This subroutine reads the output of gamess 
! Its a bit dirty! but works ...
! To work, the GAMESS source needs a little tweek.
! These modifications refer to the GAMESS version of November 2004
!---------------------------------------------------------------------
! 1) In the file int1.src:
!   After line 1343 i.e. after the line:
!     DBUG = MASWRK  .AND.
!    *       (NPRINT.EQ.3 .OR. EXETYP.EQ.DEBUG .OR. EXETYP.EQ.DBUGME)
!     add a line:
!         dbug = dbug.or.(nprint.eq.4)
! --> The purpose is to print both 1-e and 2-e integrals when NPRINT=4
!----------------------------------------------------------------------
! 2) In the file int2a.src at line 2947 it reads:
!     IF(NPRINT.NE.-5  .AND.  .NOT.CMBDIR .AND. MASWRK) THEN
!        IF(ICOUNT.LE.NINTIC) THEN
!           WRITE(IW,9015) II,JST,KST,LST,ICOUNT
!        ELSE
!           WRITE(IW,9010) II,JST,KST,LST,NREC,ICOUNT-NINTIC
!        ENDIF
!     ENDIF
!     comment out the two write statements:
!           WRITE(IW,9015) II,JST,KST,LST,ICOUNT
!           WRITE(IW,9010) II,JST,KST,LST,NREC,ICOUNT-NINTIC
! --> The purpose is not to polute the output of the integrals
!     with lines: II,JST,KST ...
!--------------------------------------------------------------------
! 3) In the file mthlib.src:
!     in the lines 3713, and 1991 increase the accuracy from f11.6 to 
!     f16.11:
!             9028 FORMAT(15X,10(6X,I4,6X))
!             9048 FORMAT(I5,2X,A8,10E21.11)
! --> The purpose is to increase the accuracy in the 1-e integrals
!--------------------------------------------------------------------
! 4) In the file mthlib.src the lines 1835-1838:
!             9028 FORMAT(15X,10(7X,I4,4X))
!             9048 FORMAT(I5,2X,A8,10F17.12)
!             9068 FORMAT(15X,10F17.10)
!             9078 FORMAT(16X,10(8X,A4,5X))
! --> The purpose is to increase the accuracy in the output orbitals
!    
!    Also at line 4295 (approx) it should read:
! 9048 FORMAT(I5,2X,A8,10F17.12)
!--------------------------------------------------------------------
! 5) In the gamess.src: 
!
! a) The variable INTFNAM should be declared e.g. around line 50:
!           character*50 INTFNAM
!
! b) After line 431 add:
!     call GETENV('INTGRLS',INTFNAM) 
!     open(unit=42, file=INTFNAM, status='unknown')
! 
! c) Then in the rungms script add at around 150 line:
!    setenv INTGRLS $SCR/$JOB.int
!
! --> The purpose is to write 2-e integrals in an unformatted file
!     which will be opened as $SCR/$JOB.int
!-------------------------------------------------------------------
! 6) In the int2a.src replace the lines nu 495(approx)
!     IF(OUT .AND. JC.GT.0) 
!    *WRITE(IW,9000) (J1(M),J2(M),J3(M),J4(M),Q(M),N1(M),V(M),M=1,JC)
! with
!     IF(OUT .AND. JC.GT.0) then
!     do m=1,jc
!       if(abs(V(M)).gt.1.d-16) then
!         WRITE(42,'(4i4,f26.16)') J1(M),J2(M),J3(M),J4(M),V(M)
!       endif
!     enddo
!     endif
! 
! And the 1492:
!      IF (MASWRK) 
!     *WRITE (IW,9008)(J1(M),J2(M),J3(M),J4(M),Q(M),N1(M),V(M),M = 1,3)
! with:
!      IF (MASWRK) then 
!       do m = 1,3
!         if(abs(V(M)).gt.1.d-16) then
!           WRITE (42,'(4i4,f26.16)')J1(M),J2(M),J3(M),J4(M),V(M)
!         endif
!       enddo
!     endif
!
! --> The purpose is to write 2-e integrals in an unformatted file
!     (unit=42)
!-----------------------------------------------------------------------
! GAMESS COMPILATION:
!    a) amd64 everything straightforward with PGI compiler (excellent!!!)
!    b) pc-linux with Intel 8 compiler for AMD athlon: 
!         i) Remove all loop unrolling leave only -O3 optimization.
!        ii) Compile DDI with -D_UNDERSCORES=1
!-----------------------------------------------------------------------
! RUN GAMESS:
! Then in the INPUT FILE:
! 1) Gamess should be run with NPRINT=4 in $CONTRL 
! 2) also symmetry C1 0 (in the $DATA block)
! 3) and remove the next empty line so atom coordinates start right after
! 4) Add $GUESS  GUESS=HUCKEL PRTMO=.true. $END
!    HUCKEL could be anything, but the important is PRTMO=.true.
! 5) To run a finite field calculation, use: RUNTYP=FFIELD in $CONTRL,
!    $FFCALC ONEFLD=.true. EONE(1)=0.000,0.000,0.005 $END,
!    $ELMOM IEMINT=1 $END. The EONE(1) contains the input field
!    components.
!
!              Created by N. N. Lathiotakis and I.Theophilou
!
!--------------------------------------------------------------------------

!..Global
   use global; use files; use matrices; use orbocc; use energies; use functional_m
   use integrals
   implicit none

!..Local
   character(20) :: line
   character(50) :: dumm
   character(30) :: dumm1
   character(22) :: dumm4
   character(40) :: dumm5
   character(23) :: dumm8!iris
   character(3)  :: dumm9!iris
   character(10) :: dumm10!iris
   integer :: limit,i,ia,iat,jat,it,j,ilin2
   integer :: n_filled, n_active !iris
   integer :: nblock, ib, imin,imax, idum
   integer :: ibs, iorb, ichrg_tot, linear_ind
   integer :: ialpha, ibeta, ispst
   integer :: iamin, iamax
   logical :: ffile, exist_file, exist_int, onefld, linear_depend
   real(dp) :: repenergy
   real(dp),allocatable :: occ_r(:,:), step_dir_t(:)

   print*,'Reading from file: ',gam_out
   inquire(file=gam_out,exist=ffile)
   if(.not.ffile) stop 'Readgamess: Gamess output doesnt exist'
   open(unit=1,file=gam_out,status='old')

   print*,'Reading the integrals from file: ',int_file
   inquire(file=int_file,exist=ffile)
   if(.not.ffile) stop 'Readgamess: Integrals file doesnt exist'
!  open(unit=42,file=int_file,status='old',form='binary')
   open(unit=42,file=int_file,status='old',form='unformatted', access='stream')
   limit=10000000

!.....START READING:
!.....read the atomic coordinates
   do i=1,limit
      read(1,'(a20)') line
      if(line == ' TOTAL NUMBER OF ATO') exit
   enddo
   backspace(1)
   read(1,'(48x,i4)')natom
   if( natom < 1 ) stop 'readgamess:unaccepted number of atoms'
   lnatom=natom

   allocate (name_ch(lnatom), charge(lnatom), xcoo(lnatom), ycoo(lnatom), zcoo(lnatom))

   rewind(1)
   do i=1,limit
      read(1,'(a20)') line
      if(line == '           CHARGE   ') exit
   enddo

   read(1,121,err=21) (name_ch(iat), charge(iat), xcoo(iat), &
                       ycoo(iat), zcoo(iat), iat=1,lnatom)
      
21 continue
   
   allocate (x_AB(lnatom, lnatom), y_AB(lnatom, lnatom), z_AB(lnatom, lnatom), R_AB(lnatom, lnatom) )

!..Calculate Atomic vector matrices:
!..x_AB, y_AB, z_AB: unit vectors along A-B, R_AB: distance of A,B
   do iat=1,natom
      do jat=1,iat-1
         x_AB(iat,jat) = xcoo(jat)-xcoo(iat)
         y_AB(iat,jat) = ycoo(jat)-ycoo(iat)
         z_AB(iat,jat) = zcoo(jat)-zcoo(iat)

         R_AB(iat,jat) = sqrt(x_AB(iat,jat)*x_AB(iat,jat)& 
                             +y_AB(iat,jat)*y_AB(iat,jat)&
                             +z_AB(iat,jat)*z_AB(iat,jat))

         x_AB(iat,jat)=x_AB(iat,jat)/R_AB(iat,jat)
         y_AB(iat,jat)=y_AB(iat,jat)/R_AB(iat,jat)
         z_AB(iat,jat)=z_AB(iat,jat)/R_AB(iat,jat)

         x_AB(jat,iat)=-x_AB(iat,jat); y_AB(jat,iat)=-y_AB(iat,jat)
         z_AB(jat,iat)=-z_AB(iat,jat); R_AB(jat,iat)=R_AB(iat,jat)
      enddo
   enddo

   print*,'NUMBER OF ATOMS: ',natom
   print*,'ATOM       NUC-CHARGE      X         Y         Z'
   do it=1,natom
      write(6,'(a12,4f10.5)')  name_ch(it), charge(it), xcoo(it), &
                                ycoo(it), zcoo(it)
   enddo
121 format(a12,f4.2, f17.10,3x,f17.10,3x,f17.10)
!..read parameters
   do i=1,limit
      read(1,'(a20)') line
      if(line == ' TOTAL NUMBER OF BAS') exit
   enddo
   read(1,'(a47,i5)')dumm,nbasis
   read(1,'(a47,i5)')dumm,nele(3)
   print*,'NUMBER OF ELECTRONS (GAMESS)  : ',nele(3)
   print*,'NUMBER OF  BASIS-SET  ELEMENTS: ',nbasis
   read(1,'(a47,i5)')dumm, ichrg_tot !The total charge
   read(1,'(a47,i5)')dumm, ispst !The spin multiplicity
   read(1,'(a47,i5)')dumm, ialpha !Nu of spin up orbitals
   read(1,'(a47,i5)')dumm, ibeta !Nu of spin down orbitals
   read(1,*)   ! Total number of atoms
   read(1,'(a35,f17.10)')dumm,repenergy

   print*,'INITIAL NUMBER  OF  NATURAL   ORBITALS: ',nnatorb

   if ( Functional == 'PN5' ) then 
      nnatorb=nele(3)
      print*,'PNOF5: NUMBER OF NATURAL ORBITALS RESET TO: ',nnatorb
   endif

!..Now all dimensions are known lets allocate
   lnnatorb=max(nnatorb,1)
   lnbasis=max(nbasis,1) 
   nbasis2= nbasis*nbasis
   nbasis3=nbasis2*nbasis
   call array_allocate()
   step_dir= step_dir_init

   if (nbasis > lbsa) stop 'readgamess: System too big for hippo! Crash!'

!..Initialize the occupation numbers:
   call initialize_occ(ialpha, ibeta)

!..read 1-e integrals
   do i=1,limit
      read(1,'(a30)') dumm1
      if(dumm1 == '          1 ELECTRON INTEGRALS') exit
   enddo
   read(1,*)
   print*,'read the overlap...'
   call read_bl_m(ovlap,nbasis)
   print*,'read the hcore...'
   call read_bl_m(Hcore,nbasis)
   print*,'read the kinetic...'
   call read_bl_m(kin,nbasis)

!..Check for Linearly dependent orbitals (less orbitals printed)
   rewind(1)
   linear_depend=.false.
   do i=1, limit
     read(1,'(a20)',end=135) line
      if(line == ' NUMBER OF LINEARLY ') then
         linear_depend=.true.
       exit
       endif
   enddo
 135  continue
   rewind(1)
   if (linear_depend) then
     do i=1,limit
       read(1,'(a20)') line
       if(line == ' NUMBER OF LINEARLY ') exit
     enddo
     read(1,'(a45,i5)')dumm, linear_ind
   end if ! end if linar_depend=true
   rewind(1)

!..linear_ind is the number of linearly independent orbitals to be read
   if (.not.linear_depend)  linear_ind = nbasis

!..read dipole moment integrals / Electric fields:
   rewind(1)
   field_calc=.false.
   do i=1, limit
      read(1,'(a30)',end=137) dumm1
      if(dumm1 == '      *  FINITE FIELD CALCULAT') then
         field_calc=.true.
         exit
      endif
   enddo
137 continue
   if(field_calc) then
      onefld=.false.
      do i=1,50
         read(1,'(a22)') dumm4
         if(dumm4 == '    SYM= F  ONEFLD= T') then
            onefld=.true.
            exit
         endif
      enddo
      if (.not.onefld) then
         print*,'The finite field mode works only for user defined'
         print*,'Electric field... see ONEFLD=.true. option'
         stop 'readgamess: ONEFLD in gamess input should be .true.'
      else
         read(1,'(34x,3f9.5)') EField_x, EField_y, EField_z
         print*,'External Electric Field Present with components:'
         write(*,"(' Ex= ',f8.4,' Ey= ',f8.4,' Ez= ',f8.4)") EField_x, &
               EField_y, EField_z
      endif
        
      exist_int=.false.
      do i=1, limit
         read(1,'(a22)',end=139) dumm4
         if(dumm4 == '             X INTEGRA') then
            exist_int=.true.
            exit
         endif
      enddo
 139  continue
      if(.not.exist_int) then
         print*,'Is this a finite Field Calculation? If yes then'
         print*,'set  IEMINT=1 in $ELMOM group in Gamess input.'
         print*,'If not, then set RUNTYP=ENERGY and remove the'
         print*,'$FFCALC and $ELMOM fields from Gamess input.'
         stop 'readgamess: Can not find the X,Y,Z integrals'
      else
         print*,'read the X integrals...'
         call read_bl_m(x_integ, nbasis)
         print*,'read the Y integrals...'
         read(1,*)
         read(1,*)
         call read_bl_m(y_integ, nbasis)
         print*,'read the Z integrals...'
         read(1,*)
         read(1,*)
         call read_bl_m(z_integ, nbasis)
      endif !(.not.exist_int) 

!....Add to the 1-e integrals:
      do i=1,nbasis
         do j=1,i
           hcore(i,j) = hcore(i,j) + EField_x * x_integ(i,j) &
                                   + EField_y * y_integ(i,j) &
                                   + EField_z * z_integ(i,j)
         enddo
      enddo
   endif !field_calc

!..Hermiticity of hcore:
   do i=1,nbasis
      do j=i,nbasis
         hcore(i,j)=hcore(j,i)
         ovlap(i,j)=ovlap(j,i)
      enddo
   enddo

!..read 2-e integrals
   rewind(1)
   do i=1,limit
      read(1,'(a50)')dumm 
      if(dumm == ' TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS = ') &
            exit
   enddo
   ilin2=i
   rewind(1)
   do i=1,ilin2-1
      read(1,*)
   enddo
   read(1,'(a50, i19)') dumm,nintgr
   print*,'NUMBER OF UNIQUE 2-e INTEGRALS: ',nintgr
   print*,'Reading Integrals ...'
! Reading the integrals from the log file...
!  rewind(1)
!     do i=1,limit
!        read(1,'(a30)') dumm1
!        if(dumm1 == '          2 ELECTRON INTEGRALS') exit
!     enddo
!     do i=1,5
!        read(1,'(a30)')dumm1
!     enddo
!     read(1,55)(mu(intg),nu(intg),lambda(intg), &
!          sigma(intg),q(intg),nn(intg),twoin(intg),intg=1,nintgr)
! Reading the integrals from a separate file opened as 'binary'
! for the intel fortran compiler

!..Allocate integrals space
!  lnintgr=nintgr+1
!  allocate(twoin(lnintgr))
!  allocate(mu(lnintgr), nu(lnintgr), lambda(lnintgr), sigma(lnintgr))

!!..Reading the integrals from int file:
!!  small_in=small !reject all integrals < small_in
!   small_in=1.d-9 !reject all integrals < small_in
!   twoin_min=1.d+2
!   intg1=1
!   print*,'Number of integrals before:',nintgr
!   do intg=1,nintgr
!      read(42,'(4i4,f26.16)',err=140) &
!         mu(intg1),nu(intg1),lambda(intg1), sigma(intg1), twoin(intg1)
!         absint=abs(twoin(intg1))
!      if(absint > small_in) then 
!         if(absint < twoin_min) twoin_min=absint
!         intg1=intg1+1
!      endif
!   enddo
!   go to 141
!140 print*,'Cant read integral nu:',intg
!   print*,'Out of', nintgr, ' in total...'
!   stop 'readgamess:cant read integral file properly'
!141 nintgr=intg1-1
!   print*,'Number of integrals exceeding smallin:',nintgr
!   print*,'Smallest absolute value for 2-e integr:',twoin_min

!..read the initial guess for orbitals and occ. numbers

!..Find the appropriate point in the input file to read the orbitals:
   rewind(1)
   if(iguess == 1) then ! Initial Gamess orbitals
      do i=1,limit
         read(1,'(a40)') dumm5
         if(dumm5 == '                              INITIAL GU') then
            print*,'Using the initial orbitals of gammess as initial'
            read(1,*)
            read(1,*)
            read(1,*)
            read(1,*)
            idum=nbasis/5 
            if(idum*5.ne.nbasis) idum=idum+1
            do j=1,idum
               read(1,*)
            enddo
            exit
         endif
      enddo
   elseif (iguess == 2) then !final Gamess orbitals
      do i=1,limit
         read(1,'(a22)',end=51) dumm4
         if(dumm4 == '          EIGENVECTORS') then
            print*,'Using the converged vectors of gamess as initial'
            read(1,*)
            goto 56
         endif
      enddo
51    stop 'readgamess: string:          EIGENVECTORS not found in log file'
   elseif (iguess == 4) then ! MCSCF natural orbitals   !iris
      do i=1,limit
         read(1,'(a20)')  line
         if (line == ' FORMING THE "STANDA') then
            read(1,'(a3,i3,a10,i3)')dumm9, n_filled ,dumm10, n_active
            maxorb = n_filled + n_active
            read(1,*)
            read(1,*)
            read(1,'(a23)') dumm8
            if (dumm8 == '          MCSCF NATURAL') then
               print*,'Using the MCSCF natural orbitals of gamess'
               read(1,*)
               goto 56
            endif
         endif
      enddo
   endif

56 continue

!..Now the orbitals will be read. For iguess=3 they will be
!..overwritten partially up to nnatorb later. 

   if(iguess == 0) then ! Delta matrix for orbitals
      do ia=1,nnatorb
         do i=1,nbasis
            vecnat(i,ia)=0
            if(i == ia) vecnat(i,ia)=1.d0
         enddo
      enddo
   elseif(iguess == 1 .or. iguess == 2) then
      nblock = linear_ind/5
      do ib=0,nblock
         read(1,*)
         read(1,*)
         read(1,*)
         read(1,'(a22)')dumm4
         do j=1,nbasis
            imin= ib*5 +1
            imax= min(imin+4,min(linear_ind,149))
            read(1,'(15x,5e21.12)') (vecnat(j,i),i=imin,imax)
         enddo
      enddo
!   else 
 !     nblock = linear_ind/5
 !     if(iguess == 4) nblock = maxorb/5
 !     do ib=0,nblock
 !        read(1,*)
 !        read(1,*)
 !        read(1,*)
 !        read(1,'(a22)')dumm4
 !        do j=1,nbasis
 !           imin= ib*5 +1
 !           imax= min(imin+4,min(linear_ind,149))
 !           imax= min(imin+4,min(maxorb,149))
 !           read(1,'(15x,5e21.12)') (vecnat(j,i),i=imin,imax)
 !        enddo
 !     enddo
!      do ia=1,nbasis
!         do j=1,nbasis
!           vec_hf(j,ia)=vecnat(j,ia) 
!         enddo
!      enddo     
 !  endif !if(iguess == 0)
 
   elseif (iguess == 3) then ! Restart read from unit 10
      print*,'Initial guess for natural orbitals and' 
      print*,'Occupation numbers will be read from',rest_file
      inquire(file=rest_file, exist=exist_file)
      if(.not.exist_file) stop 'readgamess: restart file do not exist'
      open(unit=10,file=rest_file,status='old')
      read(10,*)
      read(10,*)ibs,iorb
      read(10,*)
      read(10,*,err=211) xmu_read(1), xmu_read(2)
      goto 212
 211  xmu_read(1)=0.d0; xmu_read(2)=0.d0
 212  continue
      if(ibs /= nbasis) then 
         print*,'Restart file is for:'
         print*,ibs,' basis functions'
         print*,'while the calculation concerns'
         print*,nbasis,' basis functions'
         stop 'readgammess: Inconsistency in restart file'
      endif
      allocate (occ_r(iorb,2))
      allocate (step_dir_t(iorb))
      read(10,*)
      read(10,'(4f20.16)') (occ_r(ia,1),ia=1,iorb)
      read(10,*)
      read(10,'(4f20.16)') (occ_r(ia,2),ia=1,iorb)
      read(10,*)
      read(10,'(4e25.16)') ((vecnat(i,ia),i=1,nbasis), ia=1,nbasis)
      read(10,*)
      read(10,'(4e20.10)') (step_dir_t(ia),ia=1,iorb)
      if(iorb < nnatorb) then
         print*,'nnatorb > iorb, so the missing occ. numbers'
         print*,'will be initialized to zero'
         occnum(1:iorb,1:2) = occ_r
         step_dir(1:iorb) = step_dir_t
      else !iorb >= nnatorb
         if(iorb > nnatorb) then
            print*,'nnatorb < iorb, so only a subset of the'
            print*,'read-in occupation numbers will be used'
         else
            print*,'iorb=nnatorb. Restart file exactly matches'
         endif
         occnum(1:nnatorb,1:2) = occ_r
         if(iorb < nbasis)  occnum(iorb:nbasis,1:2)=0._dp
         step_dir(1:nnatorb) = step_dir_t
      endif
      occnum(:,3)=occnum(:,1) + occnum(:,2)
      deallocate (occ_r)
  elseif (iguess == 4) then !iris
     print*,'reading MCSCF natural orbitals'
     nblock = maxorb/5
     do ib=0,nblock!this should be changed
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,'(a22)')dumm4
        do j=1,nbasis
           imin= ib*5 +1
           imax= min(imin+4,min(maxorb,149))
           read(1,'(15x,5e21.12)') (vecnat(j,i),i=imin,imax)
        enddo
     enddo
  else
      print*,'readgamess: orbitals not read by readgamess'
   endif!iguess 

!...Now read the Hartree  Fock kinetic energy and total energy. 
!...For the strength fit. Its better to read since my run might 
!...be a restart not doing HF at all
   rewind(1) !For some unknown reason I have to rewind unit 1
   do i=1, limit
      read(1,'(a40)') dumm5
      if(dumm5 == '           NUCLEAR REPULSION ENERGY =   ') then
         read(1,*)
         read(1,'(a40,f16.10)') dumm5, HF_tot_ener
      endif
        
      if(dumm5 == '             TOTAL POTENTIAL ENERGY =   ') then
          read(1,'(a40,f16.10)') dumm5, HF_kin_ener
          exit
      endif
   enddo
   print*,'Hartree-Fock Kin Energy:',HF_kin_ener
   print*,'Hartree-Fock Exc Energy:',HF_exc_ener
   print*,'Hartree-Fock Tot Energy:',HF_tot_ener

   rewind(1)
!  !Read the EKT IP
!  if (iguess.eq.4) then !iris
!      do i=1, limit 
!         read(1,'(a22)',end=170) dumm4
!         if(dumm4 == '          EKT ORBITALS') goto 155
!      enddo

!155   continue
!     read(1,*)
!     nblock = maxorb/5
!     iamin=1
!     iamax=min(5,maxorb)
!     do ib=0,nblock
!        read(1,*)
!        read(1,*)
!        read(1,*) (ennat(ia),ia=iamin,iamax) ! reading EKT eigenvalues
!        read(1,*)
!        do j=1,nbasis
!           read(1,*) 
!        enddo
!        occnum=0.0_dp
!        iamin= iamin+5
!        iamax= min(iamin+4,min(maxorb,149))
!     end do
!     E_HOMO=ennat(maxorb)
!     print*,'The EKT IP is=', E_HOMO
!  end if
!170   print*, 'readgamess: EKT orbitals not found'
     
!..Read the Gamess orbital energies or occupation numbers

   rewind(1) 
   if (iguess == 1) then 
         print*, 'readgamess: for iguess == 1 no initial orbital energies can be trusted' 
   elseif (iguess == 2 .or. iguess == 3) then 
      do i=1,limit
         read(1,'(a22)',end=171) dumm4
         if(dumm4 == '          EIGENVECTORS') goto 156
      enddo
171   stop 'readgamess: EIGENVECTORS not found in gamess output'

156    continue
       read(1,*)

      iamin=1
      iamax=min(5,nbasis)
      nblock=linear_ind/5
      do ib=0,nblock
         read(1,*)
         read(1,*)
         read(1,*) (ennat(ia),ia=iamin,iamax)
         read(1,*)
         do j=1,nbasis
            read(1,*) 
         enddo
         iamin=iamin+5
         iamax=iamax+5
         iamax=min(iamax,min(linear_ind,149))
      enddo
      ennat_HF=ennat
      if (iguess == 3) then !overwrite ennat
         do ia=1,9
            read(10,*)
         enddo
         read(10,*) (ennat(ia),ia=1,iorb)
      endif
   else !if iguess=4 iris
      do i=1,limit
         read(1,'(a23)',end=172) dumm8
         if (dumm8 == '          MCSCF NATURAL') goto 157
      enddo
172   stop 'readgamess: NATURAL ORBITALS not found in gamess output'

157    continue
      
      print*,'reading occupation numbers from gamess MCSCF'
      read(1,*)

      nblock = maxorb/5
      iamin=1
      iamax=min(5,maxorb)
      do ib=0,nblock
         read(1,*)
         read(1,*)
         read(1,*) (occnum(ia,3),ia=iamin,iamax) ! reading occupation numbers
         read(1,*)
         do j=1,nbasis
            read(1,*) 
         enddo
         iamin= iamin+5
         iamax= min(iamin+4,min(maxorb,149))
      end do
      occnum(1:n_filled,3)=2_dp
      occnum(:,1)=0.5_dp*occnum(:,3)
      occnum(:,2)=occnum(:,1)
   endif

!  else
!     stop 'readgamess:read orbital energies for iguess=1 not implemented'
!  endif
 
!ta ebgala gia na diabastoun swsta sto invert density 5/3!   occnum=0._dp
!   occnum(1:nele(1),1:2) = 1._dp
!   occnum(1:nele(1),3) = 2._dp
!   if(iguess == 4) maxorb=nbasis

   close(1); close(10)
end subroutine readgamess

!--------------------------------------------------------------------
subroutine read_bl_m(onemat,n_n)
! This routine reads the structure gamess uses to print out the 
! 1-e integrals

!..Global
   use params_general
   implicit none

!..Arguments
   integer :: n_n
   real(dp) :: onemat(n_n,n_n)
    
!..Local
   integer :: nblocks, ii, imin, i, j, jmax 
   character(15) :: dumm2
   character(43) :: dumm3

   read(1,*)dumm3

   nblocks = n_n/5 
   if(nblocks*5==n_n) nblocks=nblocks-1
   do ii=0,nblocks
      read(1,*)
      read(1,*)
      read(1,*)
      imin = ii * 5 + 1
      do i=imin,n_n
         jmax= min(imin+4,i) 
         read(1,'(a15,5e21.12)') dumm2,(onemat(i,j),j=imin,jmax)
      enddo
   enddo

end subroutine read_bl_m

!---------------------------------------------------------------------
subroutine initialize_occ(ialpha, ibeta )
! Initializes the occupation numbers according to the spin 
! state of the system

!..Global
   use global; use functional_m; use orbocc
   implicit none

!..Arguments
   integer :: ialpha, ibeta

!..Local variables:
   integer norb, i, nsum

   if(cl_shell.and.(spin_st /= 'singlet')) then
      print*,'Warning closed shell calculation should have'
      print*,'            spin_st=singlet'
      print*,'   I assign spin_st=singlet and continue'
      spin_st='singlet'
   endif
   if(cl_shell.and.mod(nele(3),2) /= 0) then
      print*,'Closed shell calculation with odd number of electrons?'
   endif

!..Initialize to zero
   do i=1,nnatorb
      occnum(i,1) = zero; occnum(i,2) = zero; occnum(i,3) = zero+zero
   enddo

   if(cl_shell) then
!.....Singlet closed shell
      nele(1) = nele(3)/2
      nele(2) = nele(1)
      print*,'nele(1)', nele(1)
      print*,'nele(2)', nele(2)

      do i = 1, nele(1)
         occnum(i,1)= 1._dp-zero*zero; occnum(i,2)= 1._dp-zero*zero; occnum(i,3)= 2._dp-2._dp*zero*zero
      enddo
      do i = nele(1)+1, nnatorb
         occnum(i,1)= zero*zero; occnum(i,2)= zero*zero; occnum(i,3)= zero*zero+zero*zero
      enddo
   else  !not closed shell 
      if(spin_st == 'doublet') then
!........Doublet open shell
         if(mod(nele(3),2) == 0) then
            print*,'doublet is compatible with odd number of electrons'
            stop 'initialize_occ:mod(nele,2) == 0 and doublet'
         endif
         nele(1) = (nele(3)-1)/2 + 1
         nele(2) = (nele(3)-1)/2 
         do i=1,nele(2) ! Doubly occupied states
            occnum(i,1) = 1._dp-zero*zero; occnum(i,2) = 1._dp-zero*zero; occnum(i,3) = 2._dp-2._dp*zero*zero
         enddo
         occnum(nele(1),1)=1._dp-zero*zero; occnum(nele(1),2)=zero*zero; occnum(nele(1),3)=1._dp-zero*zero
      elseif(spin_st == 'triplet') then
!........Triplet open shell
         if(mod(nele(3),2) == 1) then
            print*,'triplet is compatible with even number of electrons'
            stop 'initialize_occ:mod(nele,2) == 1 and triplet'
         endif
         nele(1) = (nele(3)-2)/2 + 2
         nele(2) = (nele(3)-2)/2 
         do i=1,nele(2) ! Doubly occupied states
            occnum(i,1) = 1._dp-zero*zero; occnum(i,2) = 1._dp-zero*zero; occnum(i,3) = 2._dp-2._dp*zero*zero
         enddo
         do i=nele(2)+1,nele(1) ! Singly occupied states
            occnum(i,1) = 1._dp-zero*zero; occnum(i,2) = zero*zero; occnum(i,3) = 1._dp-zero*zero
         enddo
      elseif(spin_st == 'singlet') then
!........You can not have cl_shell==.false. and singlet
         print*, 'You choose singlet for an open shell?'
         stop 'initialize_occ: open shell and singlet'
      elseif(spin_st == 'multlet') then ! Use read in occupations
!........Spin state given as input
         nsum=0
         do i=1,nnatorb
            if(n_el(i) == 0) then ! Empty state
               occnum(i,1) = zero; occnum(i,2) = zero; occnum(i,3) = zero+zero
            elseif(n_el(i) == 2) then ! Doubly occ. state
               occnum(i,1) = one; occnum(i,2) = one; occnum(i,3) = two
               nele(1)=nele(1)+1
               nele(2)=nele(2)+1
               norb = i
            elseif(n_el(i) == 1) then ! Singly (spin up) occ. state
               occnum(i,1) = one; occnum(i,2) = zero; occnum(i,3) = one
               nele(1)=nele(1)+1
               norb = i
            elseif(n_el(i) == -1) then ! Singly (spin down) occ. state
               occnum(i,1) = zero; occnum(i,2) = one; occnum(i,3) = one
               nele(2)=nele(2)+1
               norb = i
            else                       ! Wrong: not defined
               print*,'n_el(',i,') is ',n_el(i), ' and'
               print*,' not one of (-1,0,1,2) check input!'
               stop 'initialize_occ: Not defined initial occupancy'
            endif
         enddo
         nsum=nele(1)+nele(2)
         print*,'As given in the input file we have:'
         print*,nele(1), 'Spin_UP electrons'
         print*,nele(2), 'Spin_DOWN electr.'
         if(nsum /= nele(3)) then
            print*,'Summing the occupancies we have:', nsum, 'electr.'
            print*,'While in readgamess we have:', nele(3), 'electrons'
            print*,'I assume the value:', nsum
            nele(3) = nsum
         endif
      else ! Unknown spin state
         print*,'spin state: ', spin_st, ' not defined'
         stop 'initialize_occ: Not defined spin state!'
      endif !spin_st
   endif !cl_shell
   if(nele(1) /= ialpha) then
      print*,'Warning: iniatialize_occ in readgamess.f:'
      print*,'Different number of spin up electrons '
      print*,'in input file:', nele(1), ' and gamess output:', ialpha
   endif
   if(nele(2) /= ibeta) then
      print*,'Warning: iniatialize_occ in readgamess.f:'
      print*,'Different number of spin down electrons '
      print*,'in input file:', nele(2), ' and gamess output:', ibeta 
   endif
end subroutine initialize_occ
!===============================================================================

subroutine read_basis()
!..Global      
   use global; use files; use basis_set
   implicit none

!..Local
   logical :: ffile
   character(80) :: line   
   character(70) :: str_1
   character(1) :: angm, chardum
   character(11) :: at_name, dumm_1
   integer :: iat, nga, i, iga, idum, idum1, id, ibs, ibsa
   integer :: ityp
   real(dp) :: xdum, ydum, coef_L(lnga), expo_bas_max

   print*,'Reading basis set from file: ',gam_out
   inquire(file=gam_out,exist=ffile)
   if(.not.ffile) stop 'Readgamess: Gamess output doesnt exist'
   open(unit=1,file=gam_out,status='old')

   do i=1,1000000
      read(1,'(a70)') str_1
      if ( str_1 == '  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIE' ) exit
      if ( str_1 == ' SHELL TYPE PRIMITIVE    EXPONENT          CONTRACTION COEFFICIENTS   ' ) exit
   enddo 
   read(1,*) ! Skip a line
   ibs=0
   read(1,'(a80)') line
   do iat=1,natom
      read(line,*) at_name !; print*, 'ATOM:', iat, at_name
      read(1,*)
      do ibsa=1, lbsa
         read(1,'(a80)') line
         read(line,*,err=17) idum
          dumm_1 = line
         if( dumm_1 == ' TOTAL NUMB') goto 17
         ibs=ibs+1
         iat_bas(ibs) = iat
         do iga = 1, lnga 
            nga_bas(ibs) = iga
            read(line,*,err=13) idum, angm, nga, expo_bas(iga,ibs), coef_bas(iga,ibs)
            if(angm == 'L' .or. angm == 'l') then
               read(line,*) idum, chardum, idum1, xdum, ydum, coef_L(iga) 
            endif
            read(1,'(a80)') line
            read(line,'(a11)') dumm_1
            if ( dumm_1 == '           ') goto 13
         enddo!iga=1,lnga
 13      continue
!        print*,angm
         if ( angm == 'S' .or. angm == 's' ) then
            ityp_bas(ibs) = 0
            do iga = 1, nga_bas(ibs)
               fact_norm(iga,ibs)= s_fun_norm*expo_bas(iga,ibs)**0.75_dp
            enddo
         elseif ( angm == 'P' .or. angm == 'p' ) then
            ityp_bas(ibs) = 1 !p_x
            do iga = 1, nga_bas(ibs)
               fact_norm(iga,ibs)= p_fun_norm*expo_bas(iga,ibs)**1.25_dp
            enddo
            do id = 2, 3 !p_y, p_z
               ibs = ibs + 1
               ityp_bas(ibs) = id
               iat_bas(ibs) = iat
               nga_bas(ibs) = nga_bas(ibs-1)
               do iga = 1, nga_bas(ibs)
                  expo_bas(iga,ibs) = expo_bas(iga, ibs-1)
                  coef_bas(iga,ibs) = coef_bas(iga, ibs-1)
                  fact_norm(iga,ibs) = fact_norm(iga, ibs-1)
               enddo!iga 
            enddo
         elseif ( angm == 'L' .or. angm == 'l' ) then
            ityp_bas(ibs) = 0
            do iga = 1, nga_bas(ibs)
               fact_norm(iga,ibs)= s_fun_norm*expo_bas(iga,ibs)**0.75_dp
            enddo
            do id = 1, 3 !p_x, p_y, p_z
               ibs = ibs + 1
               ityp_bas(ibs) = id
               iat_bas(ibs) = iat
               nga_bas(ibs) = nga_bas(ibs-1)
               do iga = 1, nga_bas(ibs)
                  expo_bas(iga,ibs) = expo_bas(iga, ibs-1)
                  coef_bas(iga,ibs) = coef_L(iga)
                  fact_norm(iga,ibs) = p_fun_norm*expo_bas(iga,ibs)**1.25_dp
               enddo !iga
            enddo !id
         elseif ( angm == 'D' .or. angm == 'd' ) then
            ityp_bas(ibs) = 4 !d_xx
            do iga = 1, nga_bas(ibs)
               fact_norm(iga,ibs)= d_fun_norm*expo_bas(iga,ibs)**1.75_dp
            enddo
            do id = 5, 9 !dyy, dzz, dxy, dxz, dyz 
               ibs = ibs + 1
               ityp_bas(ibs) = id
               iat_bas(ibs) = iat
               nga_bas(ibs) = nga_bas(ibs-1)
               do iga = 1, nga_bas(ibs)
                  expo_bas(iga,ibs) = expo_bas(iga, ibs-1)
                  coef_bas(iga,ibs) = coef_bas(iga, ibs-1)
                  if ( id < 7 ) then
                     fact_norm(iga,ibs) = d_fun_norm*expo_bas(iga,ibs)**1.75_dp
                  else
                     fact_norm(iga,ibs) = sqrt(3.d0)*d_fun_norm*expo_bas(iga,ibs)**1.75_dp
                  endif
               enddo!iga 
            enddo
         elseif ( angm == 'F' .or. angm == 'f' ) then!!!!CAUTION HAS NOT BEEN CHECKED!!!!!!!!
            ityp_bas(ibs) = 10 !f_xxx
            do iga = 1, nga_bas(ibs) 
               fact_norm(iga,ibs)= f_fun_norm*expo_bas(iga,ibs)**2.25_dp
            enddo
            do id = 11, 19 !fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fzzx, fzzzy, fxyz
               ibs = ibs + 1
               ityp_bas(ibs) = id
               iat_bas(ibs) = iat
               nga_bas(ibs) = nga_bas(ibs-1)
               do iga = 1, nga_bas(ibs)
                  expo_bas(iga,ibs) = expo_bas(iga, ibs-1)
                  coef_bas(iga,ibs) = coef_bas(iga, ibs-1)
                  if( id < 13) then
                     fact_norm(iga,ibs)= f_fun_norm*expo_bas(iga,ibs)**2.25_dp
                  elseif( id < 19) then
                     fact_norm(iga,ibs)= sqrt(5._dp)*f_fun_norm*expo_bas(iga,ibs)**2.25_dp
                  else
                     fact_norm(iga,ibs)= sqrt(15._dp)*f_fun_norm*expo_bas(iga,ibs)**2.25
                  endif
               enddo!iga 
            enddo !id
         elseif ( angm == 'G' .or. angm == 'g' ) then!!!!CAUTION HAS NOT BEEN CHECKED!!!!!!!!
           ityp_bas(ibs) = 20 !g_xxxx
           do iga = 1, nga_bas(ibs)
              fact_norm(iga,ibs)= g_fun_norm*expo_bas(iga,ibs)**2.75_dp
           enddo
           do id = 21, 34 ! gxxxx, gyyyy, gzzzz, gxxxy, gxxxz, gyyyx, gyyyx, etc
              ibs = ibs + 1
              ityp_bas(ibs) = id
              iat_bas(ibs) = iat
              nga_bas(ibs) = nga_bas(ibs-1)
              do iga = 1, nga_bas(ibs)
                 expo_bas(iga,ibs) = expo_bas(iga, ibs-1)
                 coef_bas(iga,ibs) = coef_bas(iga, ibs-1)
                 if( id < 23) then
                    fact_norm(iga,ibs)= g_fun_norm*expo_bas(iga,ibs)**2.75_dp
                 elseif( id < 29) then
                    fact_norm(iga,ibs)= sqrt(7._dp)*g_fun_norm*expo_bas(iga,ibs)**2.75_dp
                 elseif (id < 32) then
                    fact_norm(iga,ibs)= sqrt(35._dp/3._dp)*g_fun_norm*expo_bas(iga,ibs)**2.75_dp
                 else
                    fact_norm(iga,ibs)= sqrt(35._dp)*g_fun_norm*expo_bas(iga,ibs)**2.75_dp
                 endif
              enddo!iga 
           enddo !id
         elseif (angm =='H' .or. angm == 'h' ) then
           ityp_bas(ibs) = 35 !h_x^5
           do iga=1,nga_bas(ibs)
              fact_norm(iga,ibs) = h_fun_norm*expo_bas(iga,ibs)**3.25_dp
           enddo
           do id=36,55
              ibs=ibs+1
              ityp_bas(ibs) = id
              iat_bas(ibs) = iat
              nga_bas(ibs) = nga_bas(ibs-1)
              do iga = 1, nga_bas(ibs)
                 expo_bas(iga,ibs) = expo_bas(iga, ibs-1)
                 coef_bas(iga,ibs) = coef_bas(iga, ibs-1)
                 if (id < 38) then
                    fact_norm(iga,ibs)= h_fun_norm*expo_bas(iga,ibs)**3.25_dp
                 else if (id<44) then
                    fact_norm(iga,ibs) = 3._dp*h_fun_norm*expo_bas(iga,ibs)**3.25_dp
                 else if (id<50) then
                    fact_norm(iga,ibs) = sqrt(21._dp)*h_fun_norm*expo_bas(iga,ibs)**3.25_dp
                 else if (id<53) then
                    fact_norm(iga,ibs) = sqrt(63._dp)*h_fun_norm*expo_bas(iga,ibs)**3.25_dp
                 else
                    fact_norm(iga,ibs) = sqrt(105._dp)*h_fun_norm*expo_bas(iga,ibs)**3.25_dp
                 endif
              enddo
           enddo
         else
            print*,angm,': not implemented angular momentum'
            stop 'read_basis'
         endif
      enddo 
 17   continue
   enddo !iat=1,natom 
   close(1)

!Set the exponents for the Basis reconstruction
   do ibs=1,nbasis
      ityp = ityp_bas(ibs)
      if(ityp == 0) then ! S
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 1) then ! P_x
         expo_x(ibs) = 1._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 2) then ! P_y
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 1._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 3) then ! P_z
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 1._dp
      elseif(ityp == 4) then ! D_xx
         expo_x(ibs) = 2._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 5) then ! D_yy
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 2._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 6) then ! D_zz
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 2._dp
      elseif(ityp == 7) then ! D_xy
         expo_x(ibs) = 1._dp
         expo_y(ibs) = 1._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 8) then ! D_xz
         expo_x(ibs) = 1._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 1._dp
      elseif(ityp == 9) then ! D_yz
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 1._dp
         expo_z(ibs) = 1._dp
      elseif(ityp == 10) then ! F_xxx
         expo_x(ibs) = 3._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 11) then ! F_yyy
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 3._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 12) then ! F_zzz
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 3._dp
      elseif(ityp == 13) then ! F_xxy
         expo_x(ibs) = 2._dp
         expo_y(ibs) = 1._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 14) then ! F_xxz
         expo_x(ibs) = 2._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 1._dp
      elseif(ityp == 15) then ! F_xyy
         expo_x(ibs) = 1._dp
         expo_y(ibs) = 2._dp
         expo_z(ibs) = 0._dp
      elseif(ityp == 16) then ! F_yyz
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 2._dp
         expo_z(ibs) = 1._dp
      elseif(ityp == 17) then ! F_xzz
         expo_x(ibs) = 1._dp
         expo_y(ibs) = 0._dp
         expo_z(ibs) = 2._dp
      elseif(ityp == 18) then ! F_yzz
         expo_x(ibs) = 0._dp
         expo_y(ibs) = 1._dp
         expo_z(ibs) = 2._dp
      elseif(ityp == 19) then ! F_xyz
         expo_x(ibs) = 1._dp
         expo_y(ibs) = 1._dp
         expo_z(ibs) = 1._dp
     elseif(ityp == 20) then ! g_xxxx
        expo_x(ibs) = 4._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 0._dp
     elseif(ityp == 21) then ! g_yyyy
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 4._dp
        expo_z(ibs) = 0._dp
     elseif(ityp == 22) then ! g_zzzz
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 4._dp
     elseif(ityp == 23) then ! g_xxxy
        expo_x(ibs) = 3._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 0._dp
     elseif(ityp == 24) then !gxxxz
        expo_x(ibs) = 3._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 1._dp
     elseif(ityp == 25) then !gyyyx
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 3._dp
        expo_z(ibs) = 0._dp
     elseif(ityp == 26) then !gyyyz
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 3._dp
        expo_z(ibs) = 1._dp
     elseif(ityp == 27) then !gzzzx
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 3._dp
     elseif(ityp == 28) then !gzzzy
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 3._dp
     elseif(ityp == 29) then !gxxyy
        expo_x(ibs) = 2._dp
        expo_y(ibs) = 2._dp
        expo_z(ibs) = 0._dp
     elseif(ityp == 30) then !gxxzz
        expo_x(ibs) = 2._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 2._dp
     elseif(ityp == 31) then !gyyzz
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 2._dp
        expo_z(ibs) = 2._dp
     elseif(ityp == 32) then !g_xxyz
        expo_x(ibs) = 2._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 1._dp
     elseif(ityp == 33) then ! gyyxz
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 2._dp
        expo_z(ibs) = 1._dp
     elseif(ityp == 34) then ! gzzxy
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 2._dp
     elseif(ityp == 35) then !h x^5
        expo_x(ibs) = 5._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 0._dp
     else if (ityp == 36) then
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 5._dp
        expo_z(ibs) = 0._dp
     else if (ityp == 37) then
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 5._dp
     else if (ityp == 38) then
        expo_x(ibs) = 4._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 0._dp
     else if (ityp == 39) then
        expo_x(ibs) = 4._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 1._dp
     else if (ityp ==40) then
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 4._dp
        expo_z(ibs) = 0._dp
     else if (ityp == 41) then
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 4._dp
        expo_z(ibs) = 1._dp
     else if (ityp == 42) then
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 4._dp
     else if (ityp == 43) then
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 4._dp
     else if (ityp == 44) then
        expo_x(ibs) = 3._dp
        expo_y(ibs) = 2._dp
        expo_z(ibs) = 0._dp
     else if (ityp == 45) then
        expo_x(ibs) = 3._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 2._dp
     else if (ityp == 46) then
        expo_x(ibs) = 2._dp
        expo_y(ibs) = 3._dp
        expo_z(ibs) = 0._dp
     else if (ityp == 47) then
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 3._dp
        expo_z(ibs) = 2._dp
     else if (ityp == 48) then
        expo_x(ibs) = 2._dp
        expo_y(ibs) = 0._dp
        expo_z(ibs) = 3._dp
     else if (ityp == 49) then
        expo_x(ibs) = 0._dp
        expo_y(ibs) = 2._dp
        expo_z(ibs) = 3._dp
     else if (ityp == 50) then
        expo_x(ibs) = 3._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 1._dp
     else if (ityp ==51) then
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 3._dp
        expo_z(ibs) = 1._dp
     else if (ityp == 52) then
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 3._dp
     else if (ityp == 53) then
        expo_x(ibs) = 2._dp
        expo_y(ibs) = 2._dp
        expo_z(ibs) = 1._dp
     else if (ityp == 54) then
        expo_x(ibs) = 2._dp
        expo_y(ibs) = 1._dp
        expo_z(ibs) = 2._dp
     else if (ityp == 55) then
        expo_x(ibs) = 1._dp
        expo_y(ibs) = 2._dp
        expo_z(ibs) = 2._dp
      else
         print*,'ityp= ', ityp
         stop 'read_basis: wrong ityp'
      endif
   enddo

!  print*,'-----------------------------------------------------'
!  print*,'BASIS SET: '
!  do ibs=1,nbasis
!     print*,'ATOM: ', iat_bas(ibs)
!     print*,'TYPE: ', ityp_bas(ibs)
!     do iga=1, nga_bas(ibs)
!        write(6,'(i3,3f20.5)')iga, expo_bas(iga, ibs),  coef_bas(iga, ibs)
!     enddo
!  enddo
   expo_bas_min=1e10_dp
   do ibs=1,nbasis
      expo_bas_max=-1e10_dp
      do iga=1,nga_bas(ibs)
         expo_bas_max=max(expo_bas(iga, ibs),expo_bas_max)
      enddo
  expo_bas_min=min(expo_bas_min,expo_bas_max)
   enddo
end subroutine read_basis
!===================================================================================
subroutine read_basis_pot()
!..Global      
   use global; use files; use matrices, ONLY:ovlap_pot; use basis_set;
   use functional_m, ONLY:do_EFO;
   implicit none

!..Local
   logical :: ffile
   character(80) :: line
   character(20) :: line1
   character(70) :: str_1
   character(1) :: angm, chardum
   character(11) :: at_name, dumm_1
   character(30) :: dumm1
   character(47) :: dumm
   integer :: iat, nga, i,j, iga, idum, idum1, id, ibs, ibsa
   integer :: ityp, limit
   real(dp) :: xdum, ydum, coef_L_pot(lnga)

    if(do_EFO) gam_out_pot=gam_out

   print*,'Reading the potential basis set from file: ',gam_out_pot
   inquire(file=gam_out_pot,exist=ffile)
   if(.not.ffile) stop 'Readgamess: Gamess output for the potential basis doesnt exist'
   open(unit=1,file=gam_out_pot,status='old')

   limit=10000000
   do i=1,limit
      read(1,'(a20)',end=11) line1
      if(line1 == ' TOTAL NUMBER OF BAS') exit
   enddo
   goto 12
11 stop 'read_basis_pot:error in input file'
12 read(1,'(a47,i5)')dumm,nbasis_pot
   print*, 'nbasis_pot', nbasis_pot
   rewind(1)
   lnbasis_pot=max(nbasis_pot,1)
   allocate(  iat_bas_pot(lnbasis_pot), ityp_bas_pot(lnbasis_pot), nga_bas_pot(lnbasis_pot),&
              expo_bas_pot(lnga,lnbasis_pot), coef_bas_pot(lnga,lnbasis_pot),&
              fact_norm_pot(lnga,lnbasis_pot),&
              expo_x_pot(lnbasis_pot), expo_y_pot(lnbasis_pot), expo_z_pot(lnbasis_pot) ) 

   do i=1,1000000
      read(1,'(a70)') str_1
      if ( str_1 == '  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIE' ) exit
      if ( str_1 == ' SHELL TYPE PRIMITIVE    EXPONENT          CONTRACTION COEFFICIENTS   ' ) exit
   enddo 
   read(1,*) ! Skip a line
   ibs=0
   read(1,'(a80)') line
   do iat=1,natom
      read(line,*) at_name !; print*, 'ATOM:', iat, at_name
      read(1,*)
      do ibsa=1, lbsa
         read(1,'(a80)') line
         read(line,*,err=17) idum
          dumm_1 = line
         if( dumm_1 == ' TOTAL NUMB') goto 17
         ibs=ibs+1
         iat_bas_pot(ibs) = iat
         do iga = 1, lnga
            nga_bas_pot(ibs) = iga
            read(line,*,err=13) idum, angm, nga, expo_bas_pot(iga,ibs), coef_bas_pot(iga,ibs)
            if(angm == 'L' .or. angm == 'l') then
               read(line,*) idum, chardum, idum1, xdum, ydum, coef_L_pot(iga) 
            endif
            read(1,'(a80)') line
            read(line,'(a11)') dumm_1
            if ( dumm_1 == '           ') goto 13
         enddo!iga=1,lnga
 13      continue
!        print*,angm
         if ( angm == 'S' .or. angm == 's' ) then
            ityp_bas_pot(ibs) = 0
            do iga = 1, nga_bas_pot(ibs)
               fact_norm_pot(iga,ibs)= s_fun_norm*expo_bas_pot(iga,ibs)**0.75
            enddo
         elseif ( angm == 'P' .or. angm == 'p' ) then
            ityp_bas_pot(ibs) = 1 !p_x
            do iga = 1, nga_bas_pot(ibs)
               fact_norm_pot(iga,ibs)= p_fun_norm*expo_bas_pot(iga,ibs)**1.25
            enddo
            do id = 2, 3 !p_y, p_z
               ibs = ibs + 1
               ityp_bas_pot(ibs) = id
               iat_bas_pot(ibs) = iat
               nga_bas_pot(ibs) = nga_bas_pot(ibs-1)
               do iga = 1, nga_bas_pot(ibs)
                  expo_bas_pot(iga,ibs) = expo_bas_pot(iga, ibs-1)
                  coef_bas_pot(iga,ibs) = coef_bas_pot(iga, ibs-1)
                  fact_norm_pot(iga,ibs) = fact_norm_pot(iga, ibs-1)
               enddo!iga 
            enddo
         elseif ( angm == 'L' .or. angm == 'l' ) then
            ityp_bas_pot(ibs) = 0
            do iga = 1, nga_bas_pot(ibs)
               fact_norm_pot(iga,ibs)= s_fun_norm*expo_bas_pot(iga,ibs)**0.75_dp
            enddo
            do id = 1, 3 !p_x, p_y, p_z
               ibs = ibs + 1
               ityp_bas_pot(ibs) = id
               iat_bas_pot(ibs) = iat
               nga_bas_pot(ibs) = nga_bas_pot(ibs-1)
               do iga = 1, nga_bas_pot(ibs)
                  expo_bas_pot(iga,ibs) = expo_bas_pot(iga, ibs-1)
                  coef_bas_pot(iga,ibs) = coef_L_pot(iga)
                  fact_norm_pot(iga,ibs) = p_fun_norm*expo_bas_pot(iga,ibs)**1.25_dp
               enddo !iga
            enddo !id
         elseif ( angm == 'D' .or. angm == 'd' ) then
            ityp_bas_pot(ibs) = 4 !d_xx
            do iga = 1, nga_bas_pot(ibs)
               fact_norm_pot(iga,ibs)= d_fun_norm*expo_bas_pot(iga,ibs)**1.75_dp
            enddo
            do id = 5, 9 !dyy, dzz, dxy, dxz, dyz 
               ibs = ibs + 1
               ityp_bas_pot(ibs) = id
               iat_bas_pot(ibs) = iat
               nga_bas_pot(ibs) = nga_bas_pot(ibs-1)
               do iga = 1, nga_bas_pot(ibs)
                  expo_bas_pot(iga,ibs) = expo_bas_pot(iga, ibs-1)
                  coef_bas_pot(iga,ibs) = coef_bas_pot(iga, ibs-1)
                  if ( id < 7 ) then
                     fact_norm_pot(iga,ibs) = d_fun_norm*expo_bas_pot(iga,ibs)**1.75_dp
                  else
                     fact_norm_pot(iga,ibs) = sqrt(3.d0)*d_fun_norm*expo_bas_pot(iga,ibs)**1.75_dp
                  endif
               enddo!iga 
            enddo
         elseif ( angm == 'F' .or. angm == 'f' ) then!!!!CAUTION HAS NOT BEEN CHECKED!!!!!!!!
            ityp_bas_pot(ibs) = 10 !f_xxx
            do iga = 1, nga_bas_pot(ibs) 
               fact_norm_pot(iga,ibs)= f_fun_norm*expo_bas_pot(iga,ibs)**2.25_dp
            enddo
            do id = 11, 19 !fyyy, fzzz, fxxy, fxxz, fyyx, fyyz, fzzx, fzzzy, fxyz
               ibs = ibs + 1
               ityp_bas_pot(ibs) = id
               iat_bas_pot(ibs) = iat
               nga_bas_pot(ibs) = nga_bas_pot(ibs-1)
               do iga = 1, nga_bas_pot(ibs)
                  expo_bas_pot(iga,ibs) = expo_bas_pot(iga, ibs-1)
                  coef_bas_pot(iga,ibs) = coef_bas_pot(iga, ibs-1)
                  if( id < 13) then
                     fact_norm_pot(iga,ibs)= f_fun_norm*expo_bas_pot(iga,ibs)**2.25_dp
                  elseif( id < 19) then
                     fact_norm_pot(iga,ibs)= sqrt(5.d0)*f_fun_norm*expo_bas_pot(iga,ibs)**2.25_dp
                  else
                     fact_norm_pot(iga,ibs)= sqrt(15.d0)*f_fun_norm*expo_bas_pot(iga,ibs)**2.25_dp
                  endif
               enddo!iga 
            enddo !id 
        elseif ( angm == 'G' .or. angm == 'g' ) then!!!!CAUTION HAS NOT BEEN CHECKED!!!!!!!!
           ityp_bas_pot(ibs) = 20 !g_xxxx
           do iga = 1, nga_bas_pot(ibs) 
              fact_norm_pot(iga,ibs)= g_fun_norm*expo_bas_pot(iga,ibs)**2.75_dp
           enddo
           do id = 21, 34 ! gxxxx, gyyyy, gzzzz, gxxxy, gxxxz, gyyyx, gyyyx, etc
              ibs = ibs + 1
              ityp_bas_pot(ibs) = id
              iat_bas_pot(ibs) = iat
              nga_bas_pot(ibs) = nga_bas_pot(ibs-1)
              do iga = 1, nga_bas_pot(ibs)
                 expo_bas_pot(iga,ibs) = expo_bas_pot(iga, ibs-1)
                 coef_bas_pot(iga,ibs) = coef_bas_pot(iga, ibs-1)
                 if( id < 23) then
                    fact_norm_pot(iga,ibs)= g_fun_norm*expo_bas_pot(iga,ibs)**2.75_dp
                 elseif( id < 29) then
                    fact_norm_pot(iga,ibs)= sqrt(7._dp)*g_fun_norm*expo_bas_pot(iga,ibs)**2.75_dp
                 elseif (id < 32) then
                    fact_norm_pot(iga,ibs)= sqrt(35._dp/3._dp)*g_fun_norm*expo_bas_pot(iga,ibs)**2.75_dp
                 else
                    fact_norm_pot(iga,ibs)= sqrt(35._dp)*g_fun_norm*expo_bas_pot(iga,ibs)**2.75_dp
                 endif
              enddo!iga 
           enddo !id
        elseif (angm =='H' .or. angm == 'h' ) then
           ityp_bas_pot(ibs) = 35 !h_x^5
           do iga=1,nga_bas_pot(ibs)
              fact_norm_pot(iga,ibs) = h_fun_norm*expo_bas_pot(iga,ibs)**3.25_dp
           enddo
           do id=36,55
              ibs=ibs+1
              ityp_bas_pot(ibs) = id
              iat_bas_pot(ibs) = iat
              nga_bas_pot(ibs) = nga_bas_pot(ibs-1)
              do iga = 1, nga_bas_pot(ibs)
                 expo_bas_pot(iga,ibs) = expo_bas_pot(iga, ibs-1)
                 coef_bas_pot(iga,ibs) = coef_bas_pot(iga, ibs-1)
                 if (id < 38) then
                    fact_norm_pot(iga,ibs)= h_fun_norm*expo_bas_pot(iga,ibs)**3.25_dp
                 else if (id<44) then
                    fact_norm_pot(iga,ibs) = 3._dp*h_fun_norm*expo_bas_pot(iga,ibs)**3.25_dp
                 else if (id<50) then
                    fact_norm_pot(iga,ibs) = sqrt(21._dp)*h_fun_norm*expo_bas_pot(iga,ibs)**3.25_dp
                 else if (id<53) then
                    fact_norm_pot(iga,ibs) = sqrt(63._dp)*h_fun_norm*expo_bas_pot(iga,ibs)**3.25_dp
                 else
                    fact_norm_pot(iga,ibs) = sqrt(105._dp)*h_fun_norm*expo_bas_pot(iga,ibs)**3.25_dp
                 endif
              enddo
           enddo

         else
            print*,angm,': not implemented angular momentum'
            stop 'read_basis'
         endif
      enddo 
 17   continue
   enddo !iat=1,natom 

   rewind(1) ! Now read the overlap of the potential basis elements
     
  do i=1,limit
     read(1,'(a30)') dumm1
     if(dumm1 == '          1 ELECTRON INTEGRALS') exit
  enddo
  read(1,*)
  print*,'read the overlap of aux. basis functions...'
  if(.not.allocated(ovlap_pot)) allocate(ovlap_pot(lnbasis_pot,lnbasis_pot))
  call read_bl_m(ovlap_pot,nbasis_pot)
  do i=1,nbasis_pot
     do j=1,i
        ovlap_pot(j,i)=ovlap_pot(i,j)
     enddo
  enddo

  close(1)

!Set the exponents for the Basis reconstruction
   do ibs=1,nbasis_pot
      ityp = ityp_bas_pot(ibs)
      if(ityp == 0) then ! S
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 1) then ! P_x
         expo_x_pot(ibs) = 1._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 2) then ! P_y
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 1._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 3) then ! P_z
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 1._dp
      elseif(ityp == 4) then ! D_xx
         expo_x_pot(ibs) = 2._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 5) then ! D_yy
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 2._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 6) then ! D_zz
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 2._dp
      elseif(ityp == 7) then ! D_xy
         expo_x_pot(ibs) = 1._dp
         expo_y_pot(ibs) = 1._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 8) then ! D_xz
         expo_x_pot(ibs) = 1._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 1._dp
      elseif(ityp == 9) then ! D_yz
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 1._dp
         expo_z_pot(ibs) = 1._dp
      elseif(ityp == 10) then ! F_xxx
         expo_x_pot(ibs) = 3._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 11) then ! F_yyy
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 3._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 12) then ! F_zzz
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 3._dp
      elseif(ityp == 13) then ! F_xxy
         expo_x_pot(ibs) = 2._dp
         expo_y_pot(ibs) = 1._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 14) then ! F_xxz
         expo_x_pot(ibs) = 2._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 1._dp
      elseif(ityp == 15) then ! F_xyy
         expo_x_pot(ibs) = 1._dp
         expo_y_pot(ibs) = 2._dp
         expo_z_pot(ibs) = 0._dp
      elseif(ityp == 16) then ! F_yyz
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 2._dp
         expo_z_pot(ibs) = 1._dp
      elseif(ityp == 17) then ! F_xzz
         expo_x_pot(ibs) = 1._dp
         expo_y_pot(ibs) = 0._dp
         expo_z_pot(ibs) = 2._dp
      elseif(ityp == 18) then ! F_yzz
         expo_x_pot(ibs) = 0._dp
         expo_y_pot(ibs) = 1._dp
         expo_z_pot(ibs) = 2._dp
      elseif(ityp == 19) then ! F_xyz
         expo_x_pot(ibs) = 1._dp
         expo_y_pot(ibs) = 1._dp
         expo_z_pot(ibs) = 1._dp
    elseif(ityp == 20) then ! g_xxxx
        expo_x_pot(ibs) = 4._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 0._dp
     elseif(ityp == 21) then ! g_yyyy
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 4._dp
        expo_z_pot(ibs) = 0._dp
     elseif(ityp == 22) then ! g_zzzz
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 4._dp
     elseif(ityp == 23) then ! g_xxxy
        expo_x_pot(ibs) = 3._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 0._dp
     elseif(ityp == 24) then !gxxxz
        expo_x_pot(ibs) = 3._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 1._dp
     elseif(ityp == 25) then !gyyyx
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 3._dp
        expo_z_pot(ibs) = 0._dp
     elseif(ityp == 26) then !gyyyz
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 3._dp
        expo_z_pot(ibs) = 1._dp
     elseif(ityp == 27) then !gzzzx
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 3._dp
     elseif(ityp == 28) then !gzzzy
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 3._dp
     elseif(ityp == 29) then !gxxyy
        expo_x_pot(ibs) = 2._dp
        expo_y_pot(ibs) = 2._dp
        expo_z_pot(ibs) = 0._dp
     elseif(ityp == 30) then !gxxzz
        expo_x_pot(ibs) = 2._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 2._dp
     elseif(ityp == 31) then !gyyzz
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 2._dp
        expo_z_pot(ibs) = 2._dp
     elseif(ityp == 32) then !g_xxyz
        expo_x_pot(ibs) = 2._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 1._dp
     elseif(ityp == 33) then ! gyyxz
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 2._dp
        expo_z_pot(ibs) = 1._dp
     elseif(ityp == 34) then ! gzzxy
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 2._dp
     elseif(ityp == 35) then !h x^5
        expo_x_pot(ibs) = 5._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 0._dp
     else if (ityp == 36) then
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 5._dp
        expo_z_pot(ibs) = 0._dp
     else if (ityp == 37) then
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 5._dp
     else if (ityp == 38) then
        expo_x_pot(ibs) = 4._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 0._dp
     else if (ityp == 39) then
        expo_x_pot(ibs) = 4._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 1._dp
     else if (ityp ==40) then
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 4._dp
        expo_z_pot(ibs) = 0._dp
     else if (ityp == 41) then
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 4._dp
        expo_z_pot(ibs) = 1._dp
     else if (ityp == 42) then
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 4._dp
     else if (ityp == 43) then
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 4._dp
     else if (ityp == 44) then
        expo_x_pot(ibs) = 3._dp
        expo_y_pot(ibs) = 2._dp
        expo_z_pot(ibs) = 0._dp
     else if (ityp == 45) then
        expo_x_pot(ibs) = 3._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 2._dp
     else if (ityp == 46) then
        expo_x_pot(ibs) = 2._dp
        expo_y_pot(ibs) = 3._dp
        expo_z_pot(ibs) = 0._dp
     else if (ityp == 47) then
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 3._dp
        expo_z_pot(ibs) = 2._dp
     else if (ityp == 48) then
        expo_x_pot(ibs) = 2._dp
        expo_y_pot(ibs) = 0._dp
        expo_z_pot(ibs) = 3._dp
     else if (ityp == 49) then
        expo_x_pot(ibs) = 0._dp
        expo_y_pot(ibs) = 2._dp
        expo_z_pot(ibs) = 3._dp
     else if (ityp == 50) then
        expo_x_pot(ibs) = 3._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 1._dp
     else if (ityp ==51) then
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 3._dp
        expo_z_pot(ibs) = 1._dp
     else if (ityp == 52) then
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 3._dp
     else if (ityp == 53) then
        expo_x_pot(ibs) = 2._dp
        expo_y_pot(ibs) = 2._dp
        expo_z_pot(ibs) = 1._dp
     else if (ityp == 54) then
        expo_x_pot(ibs) = 2._dp
        expo_y_pot(ibs) = 1._dp
        expo_z_pot(ibs) = 2._dp
     else if (ityp == 55) then
        expo_x_pot(ibs) = 1._dp
        expo_y_pot(ibs) = 2._dp
        expo_z_pot(ibs) = 2._dp
      else
         print*,'ityp= ', ityp
         stop 'f_bas: wrong ityp'
      endif
   enddo

!  print*,'-----------------------------------------------------'
!  print*,'BASIS SET: '
!  do ibs=1,nbasis_pot
!     print*,'ATOM: ', iat_bas_pot(ibs)
!     print*,'TYPE: ', ityp_bas_pot(ibs)
!     do iga=1, nga_bas_pot(ibs)
!        write(6,'(i3,3f20.5)')iga, expo_bas_pot(iga, ibs),  coef_bas_pot(iga, ibs)
!     enddo
!  enddo
end subroutine read_basis_pot
!===================================================================================
subroutine construct_orbs_ongrid
!..Global
   use global; use basis_set;use  orbocc
    implicit none
!Local
    real(dp),allocatable :: x(:), y(:), z(:) 
    integer :: ip, npoints,ibs,ia
    real(dp), external :: f_bas
    real(dp) :: ff


    open(unit=800, file='coords', status='unknown')
    open(unit=802, file='numpoints', status='unknown')
    open(unit=803, file='wf', status='unknown')
    read(802,*) npoints
    allocate(x(1:npoints),y(1:npoints),z(1:npoints))

    do ip = 1, npoints
      read(800,'(E14.7,x,E14.7,x,E14.7)') x(ip), y(ip), z(ip)
    end do
    
  
    do ia = 1, nnatorb
      do ip = 1, npoints
        ff=0.d0
        do ibs = 1, nbasis
          ff=ff+vecnat(ibs,ia)*f_bas(ibs,x(ip),y(ip),z(ip))
        end do
       write(803,'(E14.7)')ff
      end do
    end do
    
    deallocate(x,y,z)
    close(801)
    close(802)
    close(803)
end subroutine construct_orbs_ongrid
!===================================================================================
function f_bas(ibs, x, y, z)
!..Global
   use global
   use basis_set
   implicit none

!..The function
   real(dp) :: f_bas

!..Arguments
   real(dp), intent(in) :: x, y, z
   integer, intent(in) :: ibs

!..Local
   integer :: iat,  iga
   real(dp) :: xd, yd, zd, rdis2
   real(dp) :: factx, facty, factz, fact

   iat=iat_bas(ibs) 
   xd = x - xcoo(iat) 
   yd = y - ycoo(iat) 
   zd = z - zcoo(iat)

   factx=xd**expo_x(ibs)
   facty=yd**expo_y(ibs)
   factz=zd**expo_z(ibs)

   fact = factx*facty*factz

   rdis2 = xd*xd+yd*yd+zd*zd

   f_bas = 0._dp
   do iga=1,nga_bas(ibs)
      f_bas = f_bas +  coef_bas(iga,ibs)* fact_norm(iga,ibs)*exp(-expo_bas(iga,ibs)*rdis2)
   enddo
   f_bas = fact * f_bas
end function f_bas
!===============================================================================================

function f_bas_pot(ibs, x, y, z)
!..Global
   use global
   use basis_set
   implicit none

!..The function
   real(dp) :: f_bas_pot

!..Arguments
   real(dp), intent(in) :: x, y, z
   integer, intent(in) :: ibs

!..Local
   integer :: iat, iga
   real(dp) :: xd, yd, zd, rdis2
   real(dp) :: factx, facty, factz, fact

   iat=iat_bas_pot(ibs) 
   xd = x - xcoo(iat) 
   yd = y - ycoo(iat) 
   zd = z - zcoo(iat)

   factx=xd**expo_x_pot(ibs)
   facty=yd**expo_y_pot(ibs)
   factz=zd**expo_z_pot(ibs)

   fact = factx*facty*factz

   rdis2 = xd*xd+yd*yd+zd*zd

   f_bas_pot = 0._dp
   do iga=1,nga_bas_pot(ibs)
      f_bas_pot = f_bas_pot +  coef_bas_pot(iga,ibs)* fact_norm_pot(iga,ibs)*exp(-expo_bas_pot(iga,ibs)*rdis2)
   enddo
   f_bas_pot = fact * f_bas_pot
end function f_bas_pot
!===============================================================================================

subroutine gradient_f_bas(ibs, x, y, z, gx, gy, gz)
!  The gradient of a basis function with index ibs at the point (x,y,z). 
!  The x-, y-, z-components of the gradient are gx, gy, gz
!..Global
   use global
   use basis_set
   implicit none

!..Arguments
   real(dp), intent(in) :: x, y, z
   real(dp), intent(out) :: gx, gy, gz
   integer, intent(in) :: ibs

!..Local
   integer :: iat, iga
   real(dp) :: xd, yd, zd, rdis2
   real(dp) :: factx, facty, factz, fact, ff, x2, y2, z2, fx, fy, fz

   iat=iat_bas(ibs) 
   xd = x - xcoo(iat) 
   yd = y - ycoo(iat) 
   zd = z - zcoo(iat)

   factx=xd**expo_x(ibs)
   facty=yd**expo_y(ibs)
   factz=zd**expo_z(ibs)

   fact = factx*facty*factz

   rdis2 = xd*xd+yd*yd+zd*zd

   x2 = 2._dp*xd
   y2 = 2._dp*yd
   z2 = 2._dp*zd

   fx=expo_x(ibs)*xd**max(expo_x(ibs)-1._dp,0._dp)*facty*factz
   fy=expo_y(ibs)*yd**max(expo_y(ibs)-1._dp,0._dp)*factx*factz
   fz=expo_z(ibs)*zd**max(expo_z(ibs)-1._dp,0._dp)*factx*facty

   gx = 0._dp; gy=0._dp; gz=0._dp
   do iga = 1,nga_bas(ibs)
      ff = coef_bas(iga,ibs)* fact_norm(iga,ibs)*exp(-expo_bas(iga,ibs)*rdis2)
      gx  = gx + ff*(fx - fact * expo_bas(iga,ibs)*x2)
      gy  = gy + ff*(fy - fact * expo_bas(iga,ibs)*y2)
      gz  = gz + ff*(fz - fact * expo_bas(iga,ibs)*z2)
   enddo

end subroutine gradient_f_bas
!===============================================================================================

function f_bas_gradsq(ibs, x, y, z)
!  The square of the gradient at point x,y,z of the basis function with index ibs at (x,y,z)
!..Global
   use global
   use basis_set
   implicit none

!..The function
   real(dp) :: f_bas_gradsq

!..Arguments
   real(dp), intent(in) :: x, y, z
   integer, intent(in) :: ibs

!..Local 
   real(dp) :: gx, gy, gz

   call gradient_f_bas(ibs, x, y, z, gx, gy, gz)

   f_bas_gradsq=gx*gx+gy*gy+gz*gz

   return
end function f_bas_gradsq
!===============================================================================================

function f_bas_Laplacian(ibs, x, y, z)
!  The Laplacian of the basis function with index ibs at the real space point (x,y,z).
!..Global
   use global
   use basis_set
   implicit none

!..The function
   real(dp) :: f_bas_Laplacian

!..Arguments
   real(dp), intent(in) :: x, y, z
   integer, intent(in) :: ibs

!..Local
   integer :: iat,  iga
   real(dp) :: xd, yd, zd, rdis2
   real(dp) :: factx, facty, factz, fact, ff, f_L
   real(dp) :: x2, y2, z2, fx, fy, fz, gx, gy, gz

   iat=iat_bas(ibs) 
   xd = x - xcoo(iat) 
   yd = y - ycoo(iat) 
   zd = z - zcoo(iat)

   factx=xd**expo_x(ibs)
   facty=yd**expo_y(ibs)
   factz=zd**expo_z(ibs)

   fact = factx*facty*factz

   rdis2 = xd*xd+yd*yd+zd*zd

   fx = 0._dp; fy=0._dp; fz=0._dp
   gx = 0._dp; gy=0._dp; gz=0._dp
   if( xd /= 0.d0) then 
!  if ( abs(xd) > zero ) then
      fx = expo_x(ibs)/xd; gx=fx/xd
   endif
   if( yd /= 0._dp) then
!  if ( abs(yd) > zero ) then
      fy = expo_y(ibs)/yd; gy=fy/yd
   endif
   if( zd /= 0._dp) then
!  if ( abs(zd) > zero ) then
      fz = expo_z(ibs)/zd; gz=fz/zd
   endif
   x2 = 2._dp*xd; y2 = 2._dp*yd; z2 = 2._dp*zd
   
   f_L=0._dp
   do iga = 1,nga_bas(ibs)
      ff =  coef_bas(iga,ibs)* fact_norm(iga,ibs)*fact*exp(-expo_bas(iga,ibs)*rdis2)
      f_L = f_L + ((fx - expo_bas(iga,ibs)*x2)**2 -gx &
                +  (fy - expo_bas(iga,ibs)*y2)**2 -gy &
                +  (fz - expo_bas(iga,ibs)*z2)**2 -gz &
                -  6._dp*expo_bas(iga,ibs) ) * ff
   enddo

   f_bas_Laplacian=f_L
   return
end function f_bas_Laplacian

!===============================================================================================
subroutine check_Laplacian
! Check the laplacian by calculating the matrix elements of the kinetic energy <phi_i | T | phi_j>
! where |phi_i> are basis functions. Also calculates the total kinetic energy. 
!..Global
   use global; use orbocc; use matrices; use grid_params; use functional_m
   implicit none

!..Local
   integer :: ig, i, j, k, l
   real(dp), allocatable :: test_kin(:,:), G_lap(:,:), G_x(:,:), G_y(:,:), G_z(:,:)
   real(dp), external :: f_bas, f_bas_Laplacian, f_bas_gradsq
   real(dp) :: x,y,z, del, gx, gy, gz, TT
   logical :: use_laplacian

   allocate ( test_kin(lnbasis, lnbasis), G_lap(ngrid,lnbasis) ) 
   allocate ( G_x(ngrid,lnbasis), G_y(ngrid,lnbasis), G_z(ngrid,lnbasis) )

   use_laplacian=.false.
!..If true it uses the laplacian for the kinetic energy if false
!  it uses the gradients

!..The Laplacian or the gradient of a basis function on the grid
   do ig=1,ngrid
      x=x_grid(ig); y=y_grid(ig); z=z_grid(ig)
      do i=1,nbasis
         if (use_laplacian) then
            G_lap(ig,i) = f_bas_Laplacian(i,x,y,z)
!...........G_lap(ig, i) the Laplacian of the basis function with index i at the point 
!           of the grid with index ig
         else !use gradient
            call gradient_f_bas(i, x, y, z, gx, gy, gz)
            G_x(ig,i) =  gx
            G_y(ig,i) =  gy
            G_z(ig,i) =  gz
!...........G_x(ig, i), G_y(ig,i),  G_z(ig,i) the gradient components of the basis function with index i at the point 
!           of the grid with index ig
         endif
      enddo
   enddo

!..Calculation of the kinetic energy matrix element on the grid
!..test_kin(i,j) the matrix element of the kinetic energy < phi_i | T | phi_j >
   test_kin=0._dp
   do ig=1,ngrid

      x=x_grid(ig); y=y_grid(ig); z=z_grid(ig)
      do i=1,nbasis
         do j=1,nbasis
            if (use_laplacian) then
!              test_kin(i,j)=test_kin(i,j) - 0.25d0*w_grid(ig)*&
!                 (bas_f_grid(ig,j)*G_lap(ig,i)+G_lap(ig,j)*bas_f_grid(ig,i))
               test_kin(i,j)=test_kin(i,j) - 0.5_dp*w_grid(ig)*&
                  bas_f_grid(ig,i)*G_lap(ig,j)
            else !use gradients
               test_kin(i,j)=test_kin(i,j) + 0.5_dp*w_grid(ig)* &
                  (G_x(ig,i)*G_x(ig,j)+G_y(ig,i)*G_y(ig,j)+G_z(ig,i)*G_z(ig,j) )
            endif
         enddo
      enddo
   enddo

!..Comparison with analytic kinetic energy matrix elements kin(i,j)
   do i=1,nbasis
      do j=1,i
         del=test_kin(i,j)-kin(i,j)
         
!        if(abs(del) > 1.d-5 ) 
         print*,i,j,test_kin(i,j),kin(i,j)
      enddo
   enddo

!..The kinetic energy numerical: TT
   TT=0._dp
   do i=1,nbasis
      do k=1,nbasis
         do l=1,nbasis
            TT=TT+occnum(i,3)*vecnat(k,i)*vecnat(l,i)*test_kin(k,l)
         enddo
      enddo
   enddo
   print*,'Num Kinetic energy:', TT

   stop 'at Check_Laplacian'

   deallocate ( test_kin, G_lap, G_x, G_y, G_z )

end subroutine check_Laplacian


