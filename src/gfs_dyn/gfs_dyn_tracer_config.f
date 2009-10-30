!
!! ! Module: gfs_dyn_tracer_config
!
! ! Description: gfs dynamics tracer configuration module
!
! ! Revision history:
!   Aug 08 2009   Sarah Lu, initial code
!   Aug 09 2009   Sarah Lu, add ntrac_chem, ntrac_met
!   Aug 10 2009   Sarah Lu, gfs_dyn_tracer is determined from ChemRegistry
!   Oct 16 2009   Sarah Lu, remove ChemRegistry; hardwire tracer specification 
!                           for testing; port to the latest trunk
! -------------------------------------------------------------------------
!
      module gfs_dyn_tracer_config
      use gfs_dyn_machine , only : kind_grid
      implicit none
      SAVE
!
! tracer specification
!
      type    gfs_dyn_tracer_type
        character*20,    pointer     :: vname(:)    ! variable name
        integer                  :: ntrac
        integer                  :: ntrac_met
        integer                  :: ntrac_chem
      endtype gfs_dyn_tracer_type

      type (gfs_dyn_tracer_type), save     ::  gfs_dyn_tracer
!
! misc tracer options
!
      logical, save                  :: glbsum  = .true.
!

! --- public interface
      public     tracer_config_init

      contains

! -------------------------------------------------------------------   
      subroutine tracer_config_init (gfs_dyn_tracer,ntrac,
     &                               ntoz,ntcw,ncld,me)

c  
c  This subprogram sets up gfs_dyn_tracer
c 
      implicit none
! input
      integer, intent(in)    ::  me, ntoz,ntcw,ncld
! output
      type (gfs_dyn_tracer_type), intent(out)    ::  gfs_dyn_tracer
!
      integer, intent(inout)  :: ntrac
! local
      integer                 :: i, status, ierr

! ntrac_chem = number of chem tracers

      gfs_dyn_tracer%ntrac_chem = 0                        

! ntrac_met = number of met tracers
      if ( ntoz < ntcw ) then                       
        gfs_dyn_tracer%ntrac_met = ntcw + ncld - 1   
      else                                                           
        gfs_dyn_tracer%ntrac_met = ntoz                              
      endif                                          
      if ( gfs_dyn_tracer%ntrac_met /= ntrac ) then
        print *,'LU_TRC: ERROR ! inconsistency in ntrac:',
     &           ntrac, gfs_dyn_tracer%ntrac_met
        stop     
      endif

! update ntrac = total number of tracers
      gfs_dyn_tracer%ntrac = gfs_dyn_tracer%ntrac_met +     
     &                       gfs_dyn_tracer%ntrac_chem
      ntrac = gfs_dyn_tracer%ntrac

! Set up tracer name
      if ( gfs_dyn_tracer%ntrac > 0 ) then      
       allocate(gfs_dyn_tracer%vname(ntrac), stat=status)
           if( status .ne. 0 ) go to 999         
      gfs_dyn_tracer%vname(1) = 'spfh'   
      gfs_dyn_tracer%vname(ntoz) = 'o3mr'  
      gfs_dyn_tracer%vname(ntcw) = 'clwmr' 
      endif

      print *,'LU_TRC: exit tracer_config_init'
      return

999   print *,'LU_TRC: error in allocate gfs_dyn_tracer :',status,me

      end subroutine tracer_config_init

! ========================================================================= 

      end module gfs_dyn_tracer_config
