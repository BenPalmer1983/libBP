! --------------------------------------------------------------!
! Materials module
! materialsTypes, materials
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 26th June 2016
! ----------------------------------------


Module materialsTypes
! Setup Modules
  Use kinds
  Type :: materialObj
    Integer(kind=StandardInteger) :: isotopeCount
    Integer(kind=StandardInteger), Dimension(1:100) :: isotopeP
    Integer(kind=StandardInteger), Dimension(1:100) :: isotopeN
    Real(kind=DoubleReal), Dimension(1:100) :: ratioByMass
    Real(kind=DoubleReal), Dimension(1:100) :: percentageByMass
    Real(kind=DoubleReal), Dimension(1:100) :: percentageByNumber
  End Type materialObj

End Module materialsTypes


Module materials
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------!
! Vector Functions
! ------------------------------------------------------------------------!



End Module materials

!-----------------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------------------
