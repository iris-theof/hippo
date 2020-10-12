
 subroutine Invert_density(x0, y0, z0, x1, y1, z1)

!..Global
   use global; use grid_params; use basis_set; use orbocc; use matrices
   use energies; use functional_m; use DFT_par
   implicit none

!..Argument
   real(dp),intent (IN) :: x0, y0, z0, x1, y1, z1
   integer :: iscf, ig,ia, iag, info, k, l, m, n, ii, iter,nocc
   real(dp) :: alphaa, dm, dm_t, paramE, rr, x, y, z, xint, factor, utest,rho_int,obj_fun,rho_int2,Tr_sc,rho_xc 
   real(dp) :: Qneg, Tr_xc,deltaQneg, Qt, xnn
   real(dp), allocatable :: targ_dens(:), dens(:), U(:,:), Delta_rho(:),DM1(:,:),DM2(:,:)
   real(dp), allocatable :: rho_sc(:),U_tar(:,:),U_sc(:,:),D_sc(:,:), Dzer(:,:)

   if ( .not.allocated(rho_sc) ) allocate( rho_sc(ngrid) )
   if ( .not.allocated(U) ) allocate( U(nbasis, nbasis) )
   if ( .not.allocated(U_sc) ) allocate( U_sc(nbasis, nbasis) )
   if ( .not.allocated(U_tar) ) allocate( U_tar(nbasis, nbasis) )
   if ( .not.allocated(dens) ) allocate( dens(ngrid) )
   if ( .not.allocated(Delta_rho) ) allocate( Delta_rho(ngrid) )
   if ( .not.allocated(targ_dens) ) allocate( targ_dens(ngrid) )
   if ( .not.allocated(DM1) ) allocate( DM1(nbasis, nbasis) )
   if ( .not.allocated(DM2) ) allocate( DM2(nbasis, nbasis) )
   if ( .not.allocated(D_sc) ) allocate( D_sc(nbasis, nbasis) )
   if ( .not.allocated(Dzer) ) allocate( Dzer(nbasis, nbasis) )

  paramE=4.d-1
  alphaa= 1.d0
  iter = 1000
  nocc=nele(1)
  Qt=0.d0
 factor= (xnele(3)-alphaa)/(xnele(3))


 !....................Upologismos target density and target potential .......................!

    call orbs_on_grid(nnatorb)

       do ig= 1,ngrid
           targ_dens(ig)=0.d0
           do ia=1,nbasis
             targ_dens(ig)= targ_dens(ig)+occnum(ia,3)*vec_nat_grid(ig,ia)*vec_nat_grid(ig,ia)
           enddo
        enddo
  
        do k=1,nbasis
            do l=1,nbasis
               DM1(k,l)=0._dp
               do ia=1,nbasis
                  DM1(k,l)=DM1(k,l) + occnum(ia,3)*vecnat(k,ia)*vecnat(l,ia)
               enddo
            enddo
         enddo

!....NEKTARIOS
     call e2_zeros( Dzer )
!    xnn=0._dp
!    do k=1,nbasis
!       do l=1,nbasis
!          xnn=xnn+DM1(k,l)*Dzer(k,l)*ovlap(k,l)
!       enddo
!    enddo

!Disable
!   Dzer=1._dp

        do k=1,nbasis
            do l=1,nbasis
               DM1(k,l)=DM1(k,l)*Dzer(k,l)
            enddo
        enddo
        
     call  sum_intg_2(DM1,U_tar)


 !..................Initialize density and potential..................................!

     rho_sc=factor*targ_dens
     D_sc=factor*DM1 
     U=factor*U_tar

     U_sc=U 
 !......................Enter loop.....................................................! 
 
 print*, 'ENTER LOOP'
 

  open(unit=665, file='Obj_ii', status='unknown')

 do ii= 1, iter 

   print*, '--------------ITERATION No', ii ,'------------------------------'

   !........Effective Hamiltonian..........................................!

     H_eff = Hcore + U

     call diagon_H(info)
  !...........................Update......................................!


!   print*, 'CALC NEW DENSITY'

     call orbs_on_grid(nocc)


     do ig= 1,ngrid
       dens(ig)=0.d0
         do ia=1,nocc
            dens(ig)= dens(ig)+ vec_nat_grid(ig,ia)*vec_nat_grid(ig,ia)
        enddo
        dens(ig)=2*dens(ig)
     enddo

   do k=1,nbasis
     do l=1,nbasis
        DM2(k,l)=0._dp
          do ia=1,nocc
             DM2(k,l)=DM2(k,l) + vecnat(k,ia)*vecnat(l,ia)
          enddo
           DM2(k,l)=2*DM2(k,l)
     enddo
   enddo
 
        do k=1,nbasis
            do l=1,nbasis
               DM2(k,l)=DM2(k,l)*Dzer(k,l)
            enddo
        enddo

   call  sum_intg_2(DM2,U)


  !...................Check convergence...................................!

   Delta_rho=targ_dens-dens

   rho_int=0.d0
    do ig=1,ngrid
      rho_int= rho_int+ w_grid(ig)*abs(Delta_rho(ig))
     enddo


 print*, 'RHO_INT', rho_int

 !....objective function.......

  obj_fun=0.d0
  do k=1,nbasis
     do l=1,nbasis
       obj_fun=obj_fun + DM1(k,l)*U_tar(k,l) + DM2(k,l)*U(k,l) - 2*DM1(k,l)*U(k,l)
     enddo
  enddo

  print*, 'OBJECTIVE FUNCTION', obj_fun

!  if ((obj_fun .le. 1.d-8) .or. (rho_int .le. 1.d-5)) then 
!    print*, '---------------CONVERGENCE ACHIEVED-------------------------'
!     go to 1
!  endif 


!  call determine_paramE(U_tar,DM1,U_sc,U,ii,paramE)

!  print*,'BGHKE TO PARAME',paramE

  D_sc=D_sc-paramE*(DM1-DM2)
  rho_sc=rho_sc-paramE*(targ_dens-dens)
  U= U_sc - paramE*(U_tar-U)
  
  U_sc=U

!  print*,'RHO_SC',rho_sc


!----------------------CHECK INTEGRATION--------------------------------------!
!    rho_int2=0.d0
!    rho_xc=0.d0
!    do ig=1,ngrid
!      rho_int2= rho_int2+ w_grid(ig)*rho_sc(ig)
!      rho_xc=rho_xc+w_grid(ig)*(rho_sc(ig)-dens(ig))
!    enddo

  Tr_sc=0.d0
  Tr_xc=0.d0
  do k=1,nbasis
    do l=1,nbasis
  Tr_sc= Tr_sc + D_sc(k,l)*ovlap(k,l)
  Tr_xc= Tr_xc + (D_sc(k,l)-DM2(k,l))*ovlap(k,l)
    enddo
  enddo

  print*,'SCREENING DENSITY INTEGRATES TO AND RHO_xc', Tr_sc,Tr_xc
!-----------------------------------------------------------------------------!


 !----------------------Q CRITERION--------------------------------------------!

     Qneg=0.d0
    do ig=1,ngrid
      Qneg= Qneg + w_grid(ig)*0.5*(abs(rho_sc(ig))-rho_sc(ig))
     enddo

    deltaQneg=Qneg-Qt

    Qt=Qneg

  print*,'------Qneg, DeltaQneg----------',Qneg, deltaQneg
!-----------------------------------------------------------------------------!

!if (((Qneg .gt. 1.d-2) .and. (deltaQneg .gt. 5.d-3) ) .or. (Qneg .gt. 5.d-2)) then
!     print*, '---------------CONVERGENCE ACHIEVED-------------------------'
!     go to 1
!  endif


      write(665,*) ii, obj_fun, rho_int, Qneg

    print*,'I.P=', -27.2114*ennat(nele(1)), 'e.V'

 enddo!......ii...........
 
  close(665)      

 print*,'I.P=', -27.2114*ennat(nele(1)), 'e.V'

 call plot_U(dens,rho_sc,DM2,D_sc,x0,y0,z0,x1,y1,z1)

end subroutine Invert_density



  subroutine plot_U(dens,rho_sc,DM2,D_sc,x0,y0,z0,x1,y1,z1)

  use global; use matrices; use orbocc; use basis_set
   use grid_params; implicit none

!..Arguments
   real(dp),intent (IN) :: x0, y0, z0, x1, y1, z1
   real(dp),intent (IN) :: rho_sc(ngrid), dens(ngrid),DM2(nbasis,nbasis),D_sc(nbasis,nbasis)
   real(dp),external :: f_bas

!..Local
   integer :: N, ig, ist, k, l
   real(dp) :: x,y,z, xstep, ystep, zstep, rr, r_d
   real(dp) :: x_dg, y_dg, z_dg, rstep, xint
   real(dp) :: dyn, dyn_Ha, dyn_XC, rho_sc_plot, rho_xc_plot, rho_u, r_xc_u
   logical :: plotanalyt=.true.
   N=500
  

   print*,'PLOT POTENTIAL' 
!------------------------Plot Potential-------------------------------!

 
   xstep=(x1-x0)/real(N-1); ystep=(y1-y0)/real(N-1); zstep=(z1-z0)/real(N-1)
   rstep=sqrt(xstep*xstep+ystep*ystep+zstep*zstep)

   x=x0; y=y0; z=z0

   open(unit=192, file='Potential_plot', status='unknown')
  if (plotanalyt) then
   do ist=1,N
      dyn = 0._dp;  dyn_Ha= 0._dp
      do k=1,nbasis
        do l=1,nbasis
         call Xintegral2(k,l,x,y,z,xint)
          dyn=dyn + D_sc(k,l)*xint
          dyn_Ha = dyn_Ha + DM2(k,l)*xint
        enddo
      enddo
      r_d= sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      dyn_XC = dyn - dyn_Ha
      write(192,'(4f20.10)') r_d, dyn_XC, dyn_Ha, dyn
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo      
 else        
   do ist=1,N
      dyn = 0._dp;  dyn_Ha= 0._dp
     do ig=1,ngrid
         x_dg = abs(x-x_grid(ig)); y_dg=abs(y-y_grid(ig)); z_dg=abs(z-z_grid(ig))
            rr= sqrt( x_dg*x_dg + y_dg*y_dg + z_dg*z_dg)
            dyn = dyn + w_grid(ig)*rho_sc(ig)/(rr+1.e-8_dp)
            dyn_Ha = dyn_Ha + w_grid(ig)*dens(ig)/(rr+1.e-8_dp)
     enddo
      r_d= sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      dyn_XC = dyn - dyn_Ha
      write(192,'(2f20.10)') r_d, dyn_XC, dyn_Ha, dyn
      x=x+xstep; y=y+ystep; z=z+zstep
  enddo
  endif
  close(192)

 x=x0; y=y0; z=z0

  open(unit=666, file='Dens_plot', status='unknown')
   do ist=1,N
      rho_sc_plot = 0._dp;  rho_xc_plot= 0._dp; rho_u= 0._dp
      do k=1,nbasis
        do l=1,nbasis
          rho_sc_plot= rho_sc_plot + D_sc(k,l)*f_bas(k,x,y,z)*f_bas(l,x,y,z)
          rho_xc_plot=rho_xc_plot + (D_sc(k,l)-DM2(k,l))*f_bas(k,x,y,z)*f_bas(l,x,y,z)
          rho_u=rho_u + DM2(k,l)*f_bas(k,x,y,z)*f_bas(l,x,y,z)
        enddo
      enddo
      r_xc_u= rho_u- rho_sc_plot
      r_d= sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      write(666,'(5f20.10)') r_d, (r_d)**2*rho_sc_plot, rho_xc_plot, rho_u, r_xc_u
      x=x+xstep; y=y+ystep; z=z+zstep
   enddo
  close(666)

  
 end subroutine plot_U




 subroutine determine_paramE(U_tar,DM1,U_sc,U,ii,paramE)

!..Global
   use global; use grid_params; use basis_set; use orbocc; use matrices
   use energies; use functional_m; use DFT_par
   implicit none

!..Argument
   real(dp),intent (IN) :: DM1(nbasis,nbasis), U_tar(nbasis,nbasis), U_sc(nbasis,nbasis), U(nbasis,nbasis)
   integer,intent (IN) :: ii
   real(dp),intent (OUT) :: paramE
   real(dp) :: DMM1(nbasis,nbasis), DMM2(nbasis,nbasis), DMM3(nbasis,nbasis)
   real(dp) :: U1(nbasis,nbasis), U2(nbasis,nbasis), U3(nbasis,nbasis)
   integer :: k, l, ia,nocc, info
   real(dp) :: testE1,  testE2,  testE3, Rho_cb, Rho_ba, S_cb, S_ba
   real(dp):: objfun1,objfun2,objfun3,aaa,aaaa,aa

   nocc=nele(1)

!  print*,'MPHKE PARAME', paramE

  if (ii == 1) then
     testE1=1.d-4
     testE2=1.d-3
     testE3=1.d-2
  else
     testE1=paramE+0.1*paramE
     testE2=paramE
     testE3=paramE-0.1*paramE
!  print*,'testE1,E2,E3',testE1,testE2,testE3
  endif

!........No1..............

U1=U_sc-testE1*(U_tar-U)


  H_eff = Hcore + U1

  call diagon_H(info)


   do k=1,nbasis
     do l=1,nbasis
        DMM1(k,l)=0._dp
          do ia=1,nocc
             DMM1(k,l)=DMM1(k,l) + vecnat(k,ia)*vecnat(l,ia)
          enddo
           DMM1(k,l)=2*DMM1(k,l)
     enddo
   enddo

 
call  sum_intg_2(DMM1,U1)


  objfun1=0.d0
  do k=1,nbasis
     do l=1,nbasis
       objfun1 =objfun1 + DM1(k,l)*U_tar(k,l) + DMM1(k,l)*U1(k,l) - 2*DM1(k,l)*U1(k,l)
     enddo
  enddo



!......No2.................


U2=U_sc-testE2*(U_tar-U)


  H_eff = Hcore + U2

  call diagon_H(info)


  do k=1,nbasis
     do l=1,nbasis
        DMM2(k,l)=0._dp
          do ia=1,nocc
             DMM2(k,l)=DMM2(k,l) + vecnat(k,ia)*vecnat(l,ia)
          enddo
           DMM2(k,l)=2*DMM2(k,l)
     enddo
   enddo

call  sum_intg_2(DMM2,U2)

  objfun2=0.d0
  do k=1,nbasis
     do l=1,nbasis
       objfun2 =objfun2 + DM1(k,l)*U_tar(k,l) + DMM2(k,l)*U2(k,l) - 2*DM1(k,l)*U2(k,l)
     enddo
  enddo


!......No3....................


U3=U_sc-testE3*(U_tar-U)

  H_eff = Hcore + U3

  call diagon_H(info)


   do k=1,nbasis
     do l=1,nbasis
        DMM3(k,l)=0._dp
          do ia=1,nocc
             DMM3(k,l)=DMM3(k,l) + vecnat(k,ia)*vecnat(l,ia)
          enddo
           DMM3(k,l)=2*DMM3(k,l)
     enddo
   enddo

call  sum_intg_2(DMM3,U3)

  objfun3=0.d0
  do k=1,nbasis
     do l=1,nbasis
       objfun3 =objfun3 + DM1(k,l)*U_tar(k,l) + DMM3(k,l)*U3(k,l) - 2*DM1(k,l)*U3(k,l)
     enddo
  enddo


!......optimalE.............


  Rho_cb=(testE3**2-testE2**2)/(objfun3-objfun2)

  Rho_ba=(testE2**2-testE1**2)/(objfun2-objfun1)

  S_cb=(testE3-testE2)/(objfun3-objfun2)

  S_ba=(testE2-testE1)/(objfun2-objfun1)


  aaa= objfun1*(testE2-testE1)+ objfun2*(testE3-testE2) +objfun3*(testE1-testE2) 
  aaaa=testE1**2*(testE2-testE3)+testE3**2*(testE1-testE2)+testE2**2*(testE3-testE1)

  aa=aaa/aaaa
!  print*,'ALPHA',aa

  paramE=(Rho_cb-Rho_ba)/(2*(S_cb-S_ba))


  print*,'THE OPTIMAL PARAMETER E IS', paramE
  
 end subroutine determine_paramE        

!-----------------------------------------------------------------------------------------
subroutine e2_zeros( Dzer )
!..Global
   use global; use integrals; use matrices
   implicit none

!..Arguments
   real(dp), intent(OUT) :: Dzer(nbasis,nbasis)

!..Local 
   integer :: j, k
   
!..Local variables
   integer :: m, n, l, s
   integer(8) :: ind_p
   integer :: ia, ib, nrec, nchunk, nrec_last, ichunk, irec
   logical :: inv_pairs, read_int
   real(dp) :: two_int

   Dzer=0._dp

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

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(irec,m,n,l,s,two_int)
!$OMP DO
      do irec = 1, nrec
         m=mu(irec)
         n=nu(irec)
         l=lambda(irec)
         s=sigma(irec)

         two_int = twoin(irec)
         if ( abs(two_int) .gt. 1e-1_dp ) then
            Dzer(m,n)=1._dp
            Dzer(l,s)=1._dp
            Dzer(n,m)=1._dp
            Dzer(s,l)=1._dp
         endif
!        if ( ( m .eq. l) .and. ( n.eq.s ) ) then
!            if ( .not. l_Dz(m,n) ) then
!               l_Dz(m,n)=.true.
!               print*,m,n, '[mn|mn]= ', two_int, '<m|n>= ',ovlap(m,n)
!            else
!               print*,m,n, ' l_Dz(m,n) already true!'
!            endif
!        endif
!        if ( ( m .eq. s) .and. ( n.eq.l ) ) then
!            if ( .not. l_Dz(m,n) ) then
!               l_Dz(m,n)=.true.
!               print*,m,n, '[mn|mn]= ', two_int, '<m|n>= ',ovlap(m,n)
!            else
!               print*,m,n, ' l_Dz(m,n) already true!'
!            endif
!        endif
      enddo ! irec
!$OMP END DO
!$OMP END PARALLEL

   enddo ! ichunk

   return

 140 stop 'construct_f:sum_intg: error reading 2-e intergral file'

end subroutine e2_zeros
