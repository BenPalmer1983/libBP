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

Module staticCalcsTypes
! Setup Modules
  Use kinds




End Module staticCalcsTypes


Module staticCalcs
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
  Use strings
  Use constants
  Use matrix
  Use basicMaths
  Use rng
  Use linearAlgebra
  Use coordFunctions
  Use geomTypes
  Use potentialsTypes
  Use staticCalcsTypes
! Force declaration of all variables
  Implicit None
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: atomLabelIDs


!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------

! ----------------
! Prep
! ----------------

  Subroutine atomLabelIDs(potential, coords)
! Look through the atom labels from the potential and coordinates and assign unique
! atom IDs - case insensitive
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialType) :: potential
    Type(coordsType) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Character(Len=16), Dimension(1:32) :: atomLabelsPotential_A
    Character(Len=16), Dimension(1:32) :: atomLabelsPotential_B
    Character(Len=16), Dimension(1:1024) :: atomLabelsCoords
    Character(Len=16), Dimension(1:32) :: tempIDs
! Temp store
    Do i=1,32
      atomLabelsPotential_A(i) = StrToUpper(Trim(Adjustl(potential%atomLabel_A(i))))
      atomLabelsPotential_B(i) = StrToUpper(Trim(Adjustl(potential%atomLabel_B(i))))
    End Do
    Do i=1,1024
      atomLabelsCoords(i) = StrToUpper(Trim(Adjustl(coords%label(i))))
    End Do
    tempIDs = BlankStringArray(tempIDs)
 ! Check potential labels A
    Do i=1,32  ! Loop through atomLabelsPotential_A
      If(IfStringEmpty(atomLabelsPotential_A(i)))Then
        Exit
      End If
      Do j=1,32  ! Loop through tempIDs
        If(IfStringEmpty(tempIDs(j)))Then
          tempIDs(j) = atomLabelsPotential_A(i)  ! Store temp label
          potential%atomID_A(i) = j              ! Store ID
          Exit
        Else
          If(atomLabelsPotential_A(i).eq.tempIDs(j))Then
            potential%atomID_A(i) = j  ! Store ID
            Exit
          End If
        End If
      End Do
    End Do
 ! Check potential labels B
    Do i=1,32  ! Loop through atomLabelsPotential_B
      If(IfStringEmpty(atomLabelsPotential_B(i)))Then
        Exit
      End If
      Do j=1,32  ! Loop through tempIDs
        If(IfStringEmpty(tempIDs(j)))Then
          tempIDs(j) = atomLabelsPotential_B(i)  ! Store temp label
          potential%atomID_B(i) = j              ! Store ID
          Exit
        Else
          If(atomLabelsPotential_B(i).eq.tempIDs(j))Then
            potential%atomID_B(i) = j  ! Store ID
            Exit
          End If
        End If
      End Do
    End Do
! Check coords labels
    Do i=1,1024  ! Loop through coords labels
      If(IfStringEmpty(atomLabelsCoords(i)))Then
        Exit
      End If
      Do j=1,32  ! Loop through tempIDs
        If(IfStringEmpty(tempIDs(j)))Then
          tempIDs(j) = atomLabelsCoords(i)  ! Store temp label
          coords%labelID(i) = j             ! Store ID
          Exit
        Else
          If(atomLabelsCoords(i).eq.tempIDs(j))Then
            coords%labelID(i) = j  ! Store ID
            Exit
          End If
        End If
      End Do
    End Do
    Do i=1,32
      If(IfStringEmpty(tempIDs(i)))Then
        j = i - 1
        Exit
      End If
    End Do
    potential%IDcount = j
    coords%IDcount = j


  End Subroutine atomLabelIDs



  Subroutine calcEFS(coords, nl)
! Make neighbour list for atoms
! The input coords and alat must be large enough so the same atom does not interact
! with a copy in a periodic cell surrounding the original
! e.g. if the rVerlet cutoff is 5, the alat must be greater than 5
!
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(nlType) :: coords
    Type(nlType) :: nl
! Vars:  Private



  End Subroutine calcEFS


! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------













End Module staticCalcs
