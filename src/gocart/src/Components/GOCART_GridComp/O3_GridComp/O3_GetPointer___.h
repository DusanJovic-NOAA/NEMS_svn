!                          -------------------
!                          W  A  R  N  I  N  G
!                          -------------------
!
!   This code fragment is automatically generated by mapl_acg.
!   Please DO NOT edit it. Any modification made in here will be overwritten
!   next time this file is auto-generated. Instead, enter your additions
!   or deletions in the .rc file. 
!


!     Bin sizes
!     ---------
      integer, parameter              :: NBIN_O3CM = 2 ! Total ozone 

!     Bin-indexed Chem Arrays
!     -----------------------
      type(Chem_Array), target        ::    O3CM(NBIN_O3CM) ! Export: Total ozone 
      type(Chem_Array), pointer       :: ptrO3CM(:)  ! Export: Total ozone 

!     Local array referencing the Import/Export states
!     ------------------------------------------------
      type(Chem_Array), target        ::    O3 ! Export: Ozone from GOCART
      type(Chem_Array), pointer       :: ptrO3 ! Export: Ozone from GOCART
      type(Chem_Array), target        ::    O3PPMV ! Export: Ozone from GOCART
      type(Chem_Array), pointer       :: ptrO3PPMV ! Export: Ozone from GOCART

!     Get pointers to data in state
!     -----------------------------

      ptrO3 => O3   ! Ozone from GOCART
      call MAPL_GetPointer ( EXPORT, O3%data3d,  'O3', RC=STATUS )
      VERIFY_(STATUS)

      ptrO3PPMV => O3PPMV   ! Ozone from GOCART
      call MAPL_GetPointer ( EXPORT, O3PPMV%data3d,  'O3PPMV', RC=STATUS )
      VERIFY_(STATUS)

      ptrO3CM => O3CM   ! Total ozone Bin 001
      call MAPL_GetPointer ( EXPORT, O3CM(1)%data2d,  'O3CM001', RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GetPointer ( EXPORT, O3CM(2)%data2d,  'O3CM002', RC=STATUS )
      VERIFY_(STATUS)
