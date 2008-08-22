!
! !description: gfs physics gridded component error messages
!
! !revision history:
!
!  July 2007 		Shrinivas Moorti
!  November 2007 	Hann-Ming Henry Juang
!
!
! !interface:
!
      module gfs_physics_err_msg_mod

!
!!uses:
!
      use esmf_mod, ONLY: esmf_logmsgfounderror, esmf_failure, &
                          esmf_success

      implicit none

      logical, parameter :: lprint=.false.

      contains

      subroutine gfs_physics_err_msg_var(rc1,msg,var,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      character (len=*), intent(in) :: var
      if(esmf_logmsgfounderror(rc1, msg)) then
          rcfinal = esmf_failure
          print*, 'error happened in physics for ',msg,' ',var,' rc = ', rc1
          rc1     = esmf_success
      else
          if(lprint) print*, 'pass in physics for ',msg,' ',var
      end if
      end subroutine gfs_physics_err_msg_var

      subroutine gfs_physics_err_msg(rc1,msg,rc)
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rc
      character (len=*), intent(in) :: msg
      if(esmf_logmsgfounderror(rc1, msg)) then
          rc  = esmf_failure
          print*, 'error happened in physics for ',msg, ' rc = ', rc1
          rc1 = esmf_success
      else
          if(lprint) print*, 'pass in physics for ',msg
      end if
      end subroutine gfs_physics_err_msg

      subroutine gfs_physics_err_msg_final(rcfinal,msg,rc)
      integer, intent(inout)        :: rcfinal
      integer, intent(inout)        :: rc
      character (len=*), intent(in) :: msg
      if(rcfinal == esmf_success) then
         if(lprint) print*, "final pass in physics for ",msg
      else
         print*, "final fail in physics for",msg
      end if
      rc = rcfinal
      end subroutine gfs_physics_err_msg_final

      end module gfs_physics_err_msg_mod
