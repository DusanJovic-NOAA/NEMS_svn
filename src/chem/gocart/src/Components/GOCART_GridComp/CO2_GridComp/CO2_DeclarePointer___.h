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
      integer, parameter              :: NBIN_CO2SC = 4 ! CO2 Surface Concentration 
      integer, parameter              :: NBIN_CO2EM = 4 ! CO2 Emission 
      integer, parameter              :: NBIN_CO2CL = 4 ! CO2 Column Load 

!     Bin-indexed Chem Arrays
!     -----------------------
      type(Chem_Array), target        ::    CO2SC(NBIN_CO2SC) ! EXPORT: CO2 Surface Concentration 
      type(Chem_Array), pointer       :: ptrCO2SC(:)  ! EXPORT: CO2 Surface Concentration 
      type(Chem_Array), target        ::    CO2EM(NBIN_CO2EM) ! EXPORT: CO2 Emission 
      type(Chem_Array), pointer       :: ptrCO2EM(:)  ! EXPORT: CO2 Emission 
      type(Chem_Array), target        ::    CO2CL(NBIN_CO2CL) ! EXPORT: CO2 Column Load 
      type(Chem_Array), pointer       :: ptrCO2CL(:)  ! EXPORT: CO2 Column Load 

!     Local array referencing the Import/Export states
!     ------------------------------------------------
      type(Chem_Array), target        ::    CO2 ! EXPORT: Carbon Dioxide
      type(Chem_Array), pointer       :: ptrCO2 ! EXPORT: Carbon Dioxide
      type(Chem_Array), target        ::    CO2NAMER ! EXPORT: North American Carbon Dioxide
      type(Chem_Array), pointer       :: ptrCO2NAMER ! EXPORT: North American Carbon Dioxide
      type(Chem_Array), target        ::    CO2SAMER ! EXPORT: South American Carbon Dioxide
      type(Chem_Array), pointer       :: ptrCO2SAMER ! EXPORT: South American Carbon Dioxide
      type(Chem_Array), target        ::    CO2AFRIC ! EXPORT: African Carbon Dioxide
      type(Chem_Array), pointer       :: ptrCO2AFRIC ! EXPORT: African Carbon Dioxide
