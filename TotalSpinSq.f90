subroutine TotS2violation(mode, occ, viol, dviol)
!..Global
   use global; use functional_m
   implicit none

   integer :: mode
   real(dp) :: occ(lnnatorb,3), viol, dviol(lnnatorb,2)

         call BB0_viol(mode, occ, viol, dviol) ! Mueller
!  select case (functional)
!     case("RHF") 
!     case("GUM") 
!        call BB0_viol(mode, occ, viol, dviol) ! Mueller
!        call GUM_viol(mode, occ, viol, dviol) !Goedecker Umrigar
!     case("BB0") 
!        call BB0_viol(mode, occ, viol, dviol) ! Mueller
!     case("BB1")
!        call BB0_viol(mode, occ, viol, dviol) ! Mueller
!        call BB1_viol(mode, occ, viol, dviol) ! BBC1
!     case("BB2")
!        call BB0_viol(mode, occ, viol, dviol) ! Mueller
!        call BB2_viol(mode, occ, viol, dviol) ! BBC2
!     case("BB3")
!        call BB0_viol(mode, occ, viol, dviol) ! Mueller
!        call BB3_viol(mode, occ, viol, dviol) ! BBC3
!     case("CHF")
!     case("CGA")
!     case("POW")
!     case("PNO")
!     case("PNP")
!     case("MPF")
!     case("MLP")
!     case default
!        stop 'TotS2violation: not implemented functional'
!  end select 
end subroutine TotS2violation

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine BB0_viol(mode, occ, viol, dviol) 
!..Global
   use global; use functional_m
   implicit none

!..Arguments
   real(dp) occ(lnnatorb,3), viol, dviol(lnnatorb,2)
   integer :: mode

!..Local
   integer :: ia
   real(dp) :: qaa, qbb, qab, qab1, qab2, fac

   fac=1.0

   if (mode==0 .or. mode==2) then
      viol=0.d0
      do ia=1,nnatorb
         qaa=occ(ia,1)
         qbb=occ(ia,2)
         qab=2.d0*sqrt(occ(ia,1)*occ(ia,2))
         viol=viol+qaa+qbb+qab
      enddo
      if(xnele(1) > xnele(2)) then
         viol=viol-(xnele(1)+3.d0*xnele(2))
      else 
         viol=viol-(xnele(2)+3.d0*xnele(1))
      endif
      viol=0.5d0*viol*fac
   endif
  
   if(mode==1 .or. mode==2) then
      do ia=1,nnatorb
         qaa=1
         qbb=1
         qab1=sqrt(occ(ia,2)/max(occ(ia,1),zero*zero))
         qab2=sqrt(occ(ia,1)/max(occ(ia,2),zero*zero))
         dviol(ia,1)=0.5d0*fac*(qaa+qab1)
         dviol(ia,2)=0.5d0*fac*(qbb+qab2)
      enddo
   endif
end subroutine BB0_viol

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GUM_viol(mode, occ, viol, dviol) 
!..Global
   use global; use functional_m
   implicit none

!..Arguments
   real(dp) occ(lnnatorb,3), viol, dviol(lnnatorb,2)
   integer :: mode

!..Local
   integer :: ia
   real(dp) :: qaa, qbb, qab, qab1, qab2, fac

   fac=1.0

   if (mode==0 .or. mode==2) then
      viol=0.d0
      do ia=1,nnatorb
         qaa=occ(ia,1)*occ(ia,1)
         qbb=occ(ia,2)*occ(ia,2)
         qab=2.d0*sqrt(max(occ(ia,1)*occ(ia,2),small))
         viol=viol+qaa+qbb+qab
      enddo
      if(xnele(1) > xnele(2)) then
         viol=viol-(xnele(1)+3.d0*xnele(2))
      else 
         viol=viol-(xnele(2)+3.d0*xnele(1))
      endif
      viol=-0.5d0*viol*fac
   endif

   if(mode==1 .or. mode==2) then
      do ia=1,nnatorb
         qaa=2.d0*occ(ia,1)
         qbb=2.d0*occ(ia,2)
         qab1=sqrt(max(occ(ia,2)/max(occ(ia,1),small_de),small))
         qab2=sqrt(max(occ(ia,1)/max(occ(ia,2),small_de),small))
         dviol(ia,1)=-0.5d0*fac*(qaa+qab1)
         dviol(ia,2)=-0.5d0*fac*(qbb+qab2)
      enddo
   endif
end subroutine GUM_viol

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine BB1_viol(mode, occ, viol, dviol) 
!..Global
   use global; use functional_m
   implicit none

!..Arguments
   real(dp) occ(lnnatorb,3), viol, dviol(lnnatorb,2)
   integer :: mode

!..Local
   integer :: ilim1, ilim2, ia
   real(dp) :: qaa, qbb, qab, qab1, qab2, fac

   fac=1.0

   ilim1=xnele(1)+small
   ilim2=xnele(2)+small

   if (mode==0 .or. mode==2) then
      viol=0.d0
      do ia=1,nnatorb
         qaa=occ(ia,1)
         qbb=occ(ia,2)
         qab=2.d0*sqrt(max(small,occ(ia,1)*occ(ia,2)))
         if(ia > ilim1 .and. ia > ilim2) qab=-qab 
         viol=viol+qaa+qbb+qab
      enddo
      if(xnele(1) > xnele(2)) then
         viol=viol-(xnele(1)+3.d0*xnele(2))
      else 
         viol=viol-(xnele(2)+3.d0*xnele(1))
      endif
      viol=-0.5d0*viol*fac
   endif

   if(mode==1 .or. mode==2) then
      do ia=1,nnatorb
         qaa=1
         qbb=1
         if (ia > ilim1 .and. ia > ilim2) then
            qab1=-sqrt(max(occ(ia,2)/max(occ(ia,1),small_de),small))
            qab2=-sqrt(max(occ(ia,1)/max(occ(ia,2),small_de),small))
         else
            qab1=sqrt(max(occ(ia,2)/max(occ(ia,1),small_de),small))
            qab2=sqrt(max(occ(ia,1)/max(occ(ia,2),small_de),small))
         endif
         dviol(ia,1)=-0.5d0*fac*(qaa+qab1)
         dviol(ia,2)=-0.5d0*fac*(qbb+qab2)
      enddo
   endif
end subroutine BB1_viol

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine BB2_viol(mode, occ, viol, dviol) 
!..Global
   use global; use functional_m
   implicit none

!..Arguments
   real(dp) occ(lnnatorb,3), viol, dviol(lnnatorb,2)
   integer :: mode

!..Local
   integer :: ilim1, ilim2, ia
   real(dp) :: qaa, qbb, qab, qab1, qab2, fac

   fac=1.0

   ilim1=xnele(1)+small
   ilim2=xnele(2)+small

   if (mode==0 .or. mode==2) then
      viol=0.d0
      do ia=1,nnatorb
         qaa=occ(ia,1)
         qbb=occ(ia,2)
         qab=2.d0*sqrt(max(small,occ(ia,1)*occ(ia,2)))
         if(ia > ilim1 .and. ia > ilim2) qab=-qab 
         if(ia <= ilim1 .and. ia <= ilim2) then
            qab=2.d0*occ(ia,1)*occ(ia,2)
         endif
         viol=viol+qaa+qbb+qab
      enddo
      if(xnele(1) > xnele(2)) then
         viol=viol-(xnele(1)+3.d0*xnele(2))
      else 
         viol=viol-(xnele(2)+3.d0*xnele(1))
      endif
      viol=-viol*fac
   endif

   if(mode==1 .or. mode==2) then
      do ia=1,nnatorb
         qaa=1
         qbb=1
         if (ia > ilim1 .and. ia > ilim2) then
            qab1=-sqrt(max(occ(ia,2)/max(occ(ia,1),small_de),small_de))
            qab2=-sqrt(max(occ(ia,1)/max(occ(ia,2),small_de),small_de))
         elseif ( ia <= ilim1 .and. ia <= ilim2 ) then 
            qab1=2.d0*occ(ia,2)
            qab2=2.d0*occ(ia,1)
         else
            qab1=sqrt(max(occ(ia,2)/max(occ(ia,1),small_de),small_de))
            qab2=sqrt(max(occ(ia,1)/max(occ(ia,2),small_de),small_de))
         endif
         dviol(ia,1)=-fac*(qaa+qab1)
         dviol(ia,2)=-fac*(qbb+qab2)
      enddo
   endif
end subroutine BB2_viol

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine BB3_viol(mode, occ, viol, dviol) 
!..Global
   use global; use functional_m
   implicit none

!..Arguments
   real(dp) occ(lnnatorb,3), viol, dviol(lnnatorb,2)
   integer :: mode

!..Local
   integer :: ia
   real(dp) :: qaa, qbb, qab, qab1, qab2, fac

   fac=1.0

   if (mode==0 .or. mode==2) then
      if(xnele(1) > xnele(2)) then
         viol=viol+(xnele(1)+3.d0*xnele(2))
      else 
         viol=viol+(xnele(2)+3.d0*xnele(1))
      endif
      do ia=1,nnatorb
         if(kind_BBC3(ia,1)==1 .or. kind_BBC3(ia,1)==4) then
            qaa=occ(ia,1)*occ(ia,1)
         else
            qaa=occ(ia,1)
         endif
         if(kind_BBC3(ia,2)==1 .or. kind_BBC3(ia,2)==4) then
            qbb=occ(ia,2)*occ(ia,2)
         else
            qbb=occ(ia,2)
         endif
         if(kind_BBC3(ia,1)==4 .and. kind_BBC3(ia,2)==4) then 
            qab=-2.d0*sqrt(occ(ia,1)*occ(ia,2))
         elseif(kind_BBC3(ia,1)<=2 .and. kind_BBC3(ia,2)<=2) then
            qab=2.d0*occ(ia,1)*occ(ia,2)
         elseif(kind_BBC3(ia,1)<=2 .and. kind_BBC3(ia,2)==3) then
            qab=2.d0*occ(ia,1)*occ(ia,2)
         elseif(kind_BBC3(ia,1)==3 .and. kind_BBC3(ia,2)<=2 ) then
            qab=2.d0*occ(ia,1)*occ(ia,2)
         else
            qab= 2.d0*sqrt(occ(ia,1)*occ(ia,2))
         endif
         viol=viol-qaa-qbb-qab
      enddo
      viol=0.5d0*viol*fac
   endif

   if(mode==1 .or. mode==2) then
      do ia=1,nnatorb
         if(kind_BBC3(ia,1)==1 .or. kind_BBC3(ia,1)==4) then
            qaa=2.d0*occ(ia,1)
         else
            qaa=1.d0
         endif
         if(kind_BBC3(ia,2)==1 .or. kind_BBC3(ia,2)==4) then
            qbb=2.d0*occ(ia,2)
         else
            qbb=1.d0
         endif
         if(kind_BBC3(ia,1)==4 .and. kind_BBC3(ia,2)==4) then 
            qab1=-sqrt(max(occ(ia,2)/max(occ(ia,1),small_de),small))
            qab2=-sqrt(max(occ(ia,1)/max(occ(ia,2),small_de),small))
         elseif(kind_BBC3(ia,1)<=2 .and. kind_BBC3(ia,2)<=2) then
            qab1=2.d0*occ(ia,2)
            qab2=2.d0*occ(ia,1)
         elseif(kind_BBC3(ia,1)<=2 .and. kind_BBC3(ia,2)==3) then
            qab1=2.d0*occ(ia,2)
            qab2=2.d0*occ(ia,1)
         elseif(kind_BBC3(ia,1)==3 .and. kind_BBC3(ia,2)<=2 ) then
            qab1=2.d0*occ(ia,2)
            qab2=2.d0*occ(ia,1)
         else
            qab1=sqrt(max(occ(ia,2)/max(occ(ia,1),small_de),small))
            qab2=sqrt(max(occ(ia,1)/max(occ(ia,2),small_de),small))
         endif
         dviol(ia,1)=-0.5d0*fac*(qaa+qab1)
         dviol(ia,2)=-0.5d0*fac*(qbb+qab2)
      enddo
   endif
end subroutine BB3_viol


