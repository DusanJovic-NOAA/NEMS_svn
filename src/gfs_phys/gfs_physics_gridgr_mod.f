!
! !MODULE: gfs_physics_gridgr_mod  ---      Definition of the atmospheric 
!                                           fields in the ESMF internal state.
!
! !DESCRIPTION: gfs_physics_gridgr_mod ---    Define the atmospheric states
!                                             in the ESMF internal state.
!---------------------------------------------------------------------------
! !REVISION HISTORY:
!
!  2009/03/10      Sarah Lu,  Initial code.
!  2009/05/20      Sarah Lu,  Updated to the latest trunk
!  2009/08/09      Sarah Lu,  Add tracer field
!  2009/10/12      Sarah Lu,  Port to the latest trunk
!
! !INTERFACE:
!
 MODULE gfs_physics_gridgr_mod

 use machine , only : kind_grid

 IMPLICIT none

 INTEGER, PARAMETER            :: MAX_NTRAC=100

 TYPE PHY_R3D
    real(kind=kind_grid), dimension(:,:,:), pointer :: flds
 END TYPE PHY_R3D

 TYPE Grid_Var_Data
! 2D fields
    real(kind=kind_grid),  pointer:: z(:,:)
    real(kind=kind_grid),  pointer:: ps(:,:)

! 3D fields
    real(kind=kind_grid),  pointer:: u(:,:,:)
    real(kind=kind_grid),  pointer:: v(:,:,:)
    real(kind=kind_grid),  pointer:: t(:,:,:)
    real(kind=kind_grid),  pointer:: p(:,:,:)
    real(kind=kind_grid),  pointer:: dp(:,:,:)
    real(kind=kind_grid),  pointer:: dpdt(:,:,:)

! tracer fields
    real(kind=kind_grid),  pointer:: q(:,:,:)
    real(kind=kind_grid),  pointer:: oz(:,:,:)
    real(kind=kind_grid),  pointer:: cld(:,:,:)

! chemical tracer fields
    TYPE (PHY_R3D), DIMENSION (MAX_NTRAC)  :: ChemTracers
   
 end type Grid_Var_Data

 contains

! !IROUTINE: gridvar_aldata ---

!---------------------------------------------------------------------------
    subroutine gridvar_aldata(dim1, dim2, dim3, dim4, grid_fld, iret)

    TYPE(Grid_Var_Data), INTENT(out)   :: grid_fld
    integer, intent(in)                :: dim1, dim2, dim3, dim4
    implicit none

    integer, intent(out)             :: iret
!
allocate(                                    &
           grid_fld%z      (dim1,dim2),      &
           grid_fld%ps     (dim1,dim2),      &
           grid_fld%u      (dim1,dim2,dim3), &
           grid_fld%v      (dim1,dim2,dim3), &
           grid_fld%t      (dim1,dim2,dim3), &
           grid_fld%p      (dim1,dim2,dim3), &
           grid_fld%dp     (dim1,dim2,dim3), &
           grid_fld%dpdt   (dim1,dim2,dim3), &
           grid_fld%q      (dim1,dim2,dim3), &
           grid_fld%oz     (dim1,dim2,dim3), &
           grid_fld%cld    (dim1,dim2,dim3), &
           stat=iret)
    if(iret.ne.0) iret=-3

    return
  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 END MODULE gfs_physics_gridgr_mod
