!------------------------------------------------------------------------------
! Subset of the music module defaults.
!------------------------------------------------------------------------------
Module defaults
  Implicit None
  Save

  Public 

  !--------------------------------------------
  ! Type Definitions
  !--------------------------------------------
  Integer, Parameter   :: RDbl = Selected_real_kind(10, 50)
  !--------------------------------------------
  ! Mathematical Constants
  !--------------------------------------------
  Real(kind=RDbl), Parameter  :: pi    = 3.1415927_RDbl
  Real(kind=RDbl), Parameter  :: twopi = 2.0_RDbl*pi
  Real(kind=RDbl), Parameter  :: degTorad = pi/180.0_RDbl

End Module defaults




