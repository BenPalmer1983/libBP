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
  Use printModTypes
  Use printMod
  Use matrix
  Use mpiSubsTypes
  Use mpiSubs
  Use basicMaths
  Use keysMod
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
  Public :: printAtomLabelIDs
  Public :: calcEFS


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
    Type(coordsType), Dimension(:) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, cKey
    Character(Len=16), Dimension(1:32) :: atomLabelsPotential_A
    Character(Len=16), Dimension(1:32) :: atomLabelsPotential_B
    Character(Len=16), Dimension(1:Size(coords),1:1024) :: atomLabelsCoords
    Character(Len=16), Dimension(1:32) :: tempIDs
    !Do cKey=1,Size(coords,1)
! Temp store
    Do i=1,32
      atomLabelsPotential_A(i) = StrToUpper(Trim(Adjustl(potential%atomLabel_A(i))))
      atomLabelsPotential_B(i) = StrToUpper(Trim(Adjustl(potential%atomLabel_B(i))))
    End Do
    Do cKey=1,Size(coords,1)
      Do i=1,1024
        atomLabelsCoords(cKey,i) = StrToUpper(Trim(Adjustl(coords(cKey)%label(i))))
      End Do
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
    Do cKey=1,Size(coords,1)
      Do i=1,1024  ! Loop through coords labels
        If(IfStringEmpty(atomLabelsCoords(cKey,i)))Then
          Exit
        End If
        Do j=1,32  ! Loop through tempIDs
          If(IfStringEmpty(tempIDs(j)))Then
            tempIDs(j) = atomLabelsCoords(cKey, i)  ! Store temp label
            coords(cKey)%labelID(i) = j             ! Store ID
            Exit
          Else
            If(atomLabelsCoords(cKey,i).eq.tempIDs(j))Then
              coords(cKey)%labelID(i) = j  ! Store ID
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
    End Do
    potential%IDcount = j
    Do cKey=1,Size(coords,1)
      coords%IDcount = j
    End Do
! Store full key in coords obj and potential obj
    Do i=1,32
      potential%atomIDs(i) = tempIDs(i)
      Do cKey=1,Size(coords,1)
        coords(cKey)%atomIDs(i) = tempIDs(i)
      End Do
    End Do
  End Subroutine atomLabelIDs


  Subroutine printAtomLabelIDs(coords)
! Print Label
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(coordsType), Dimension(:) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Character(Len=64) :: tempLine
! Print out atom ids
    Do i=1,32
      If(IfStringEmpty(coords(1)%atomIDs(i)))Then
        Exit
      Else
        Write(tempLine,*) i,coords(1)%atomIDs(i)
        Call addLinePage(tempLine)
        print *,tempLine
      End If
    End Do
  End Subroutine printAtomLabelIDs



  Subroutine calcEFS(coords, nl, potential, cKeyIn, resetForcesIn)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger), Optional :: cKeyIn
    Logical, Optional :: resetForcesIn
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey
    Logical :: resetForces
! Optional Arguments
    resetForces = .false.
    cKey = 0
    If(Present(resetForcesIn))Then
      resetForces = resetForcesIn
    End If
    If(Present(cKeyIn))Then
      cKey = cKeyIn
    End If
! Calc EFS
    Call calcEFS_MPI(coords, nl, potential, cKey, resetForces)
  End Subroutine calcEFS

  Subroutine calcEFS_MPI(coords, nl, potential, cKey, resetForces)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
    Logical :: resetForces
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, loopID
    Type(mpiObj) :: efsMPI
! One config
    If(cKey.gt.0)Then
! Reset forces if required (across all processes)
      If(resetForces)Then
        Do i=1,p_cMax
          Do j=1,3
            coords(cKey)%forces(i,j) = 0.0D0
          End Do
        End Do
      End If
    Else
! Init mpi
      Call m_initMPI(efsMPI)
      loopID = 0
      Do cKey=1,p_confs
! Reset forces if required (across all processes)
        If(resetForces)Then
          Do i=1,p_cMax
            Do j=1,3
              coords(cKey)%forces(i,j) = 0.0D0
            End Do
          End Do
        End If
        If((coords(cKey)%points.gt.0).and.(nl(cKey)%length.gt.0))Then
          loopID = loopID + 1
          If(M_Loop(efsMPI, loopID))Then
            Call calcEFS_Action(coords, nl, potential, cKey)
          End If
        End If
      End Do
    End If
  End Subroutine calcEFS_MPI

  Subroutine calcEFS_Action(coords, nl, potential, cKey)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, pKey, dKey
    Integer(kind=StandardInteger), Dimension(1:32,0:32) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Dimension(1:32,0:32) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    Logical :: pairPot, densPot


!-------------------------
! Pair Potential Keys
!-------------------------
    pairKeyArray = 0
! Loop through potentials
    Do i=1,potential%fCount
      pairPot = .false.
      If(StrMatch(potential%fPotential(i),"PAIR"))Then
        pairPot = .true.
      End If
      If(StrMatch(potential%fPotential(i),"MORSE"))Then
        pairPot = .true.
      End If
      If(pairPot)Then
        pKey = DoubleKey(potential%atomID_A(i), potential%atomID_B(i))
        print *,potential%atomID_A(i),potential%atomID_B(i),i,j,pKey
        Do j=1,32
          If(pairKeyArray(pKey,j).eq.0)Then
            pairKeyArray(pKey,j) = i
            pairKeyArray(pKey,0) = j
            Exit
          End If
        End Do
      End If
    End Do
!-------------------------
! Density Keys  (EAM, 2BMEAM)
!-------------------------
    densityKeyArray = 0
! Loop through potentials
    Do i=1,potential%fCount
      pairPot = .false.
      If(StrMatch(potential%fPotential(i),"DENS"))Then
        densPot = .true.
      End If
      If(StrMatch(potential%fPotential(i),"DDEN"))Then
        densPot = .true.
      End If
      If(pairPot)Then
        dKey = DoubleKey(potential%atomID_A(i), potential%atomID_B(i))
        Do j=1,32
          If(densityKeyArray(dKey,j).eq.0)Then
            densityKeyArray(dKey,j) = i
            densityKeyArray(dKey,0) = j
            Exit
          End If
        End Do
      End If
    End Do




  !coords(cKey)%forces(nl(cKey)%atomA_ID(i),1) = 0.0D0
!------------------------------
! First neighbour list loop
!------------------------------

    Do i=1,nl(cKey)%length
! Get pair key
      pKey = DoubleKey(nl(cKey)%atomA_Type(i), nl(cKey)%atomB_Type(i))


      If(i.lt.5)Then
        print *,nl(cKey)%atomA_Type(i), nl(cKey)%atomB_Type(i)
        print *,i,nl(cKey)%atomPairKey(i)
        Do j=1,pairKeyArray(pKey,0)
          print *,"   pair ",j
        End Do
      End If
    End Do
    print *,"NL complete:"
    print *,nl(cKey)%length


! PAIR Potentials



! EAMPAIR Potentials




  End Subroutine calcEFS_Action


! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------













End Module staticCalcs
