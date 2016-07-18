! --------------------------------------------------------------!
! Static Caculations module
! staticCalcsTypes, staticCalcs
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calculates the forces between atoms
! Calculates the energies of atoms/total energy
! Calculates stresses
!
! ----------------------------------------
! Updated: 21st May 2016
! ----------------------------------------

Module reactionRateTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_
! Make private
  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived data types
!  Public :: o

!  Type :: o

End Module reactionRateTypes


Module reactionRate
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
  Use reactionRateTypes
! Force declaration of all variables
  Implicit None
! Make private
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines


!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------



End Module reactionRate
