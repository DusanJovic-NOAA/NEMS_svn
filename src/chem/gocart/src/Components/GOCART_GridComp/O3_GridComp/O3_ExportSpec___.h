!                          -------------------
!                          W  A  R  N  I  N  G
!                          -------------------
!
!   This code fragment is automatically generated by mapl_acg.
!   Please DO NOT edit it. Any modification made in here will be overwritten
!   next time this file is auto-generated. Instead, enter your additions
!   or deletions in the .rc file. 
!


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'O3',  &
        LONG_NAME          = 'Ozone from GOCART',  &
        UNITS              = 'kg/kg', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
        PRECISION          = KIND(0.0),&
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'O3PPMV',  &
        LONG_NAME          = 'Ozone from GOCART',  &
        UNITS              = 'ppmv', &
        DIMS               = MAPL_DimsHorzVert,    &
        VLOCATION          = MAPL_VLocationCenter,    &
        PRECISION          = KIND(0.0),&
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'O3CM001',  &
        LONG_NAME          = 'Total ozone Bin 001',  &
        UNITS              = 'Dobsons', &
        DIMS               = MAPL_DimsHorzOnly,    &
        VLOCATION          = MAPL_VLocationNone,    &
        PRECISION          = KIND(0.0),&
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,  &
        SHORT_NAME         = 'O3CM002',  &
        LONG_NAME          = 'Total ozone Bin 002',  &
        UNITS              = 'Dobsons', &
        DIMS               = MAPL_DimsHorzOnly,    &
        VLOCATION          = MAPL_VLocationNone,    &
        PRECISION          = KIND(0.0),&
                                                       RC=STATUS  )
     VERIFY_(STATUS)

