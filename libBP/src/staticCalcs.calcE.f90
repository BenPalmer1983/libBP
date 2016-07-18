! -------------------------------------------------
!  Include File:   Calc energy only
!
! -------------------------------------------------
  Subroutine calcE(nl, potential, cKeyIn)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger), Optional :: cKeyIn
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey
! Optional Arguments
    cKey = 0
    If(Present(cKeyIn))Then
      cKey = cKeyIn
    End If
! Calc EFS
    Call calcE_MPI(nl, potential, cKey)
  End Subroutine calcE

  Subroutine calcE_MPI(nl, potential, cKey)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: loopID
    Type(mpiObj) :: eMPI
! One config
    If(cKey.gt.0)Then
! Run E calculation
      If(nl(cKey)%length.gt.0)Then  ! Only run if there are coord points and neighbour list > 0
        Call calcE_Action(nl, potential, cKey)
      End If
    Else
! Init mpi
      Call m_initMPI(eMPI)
      loopID = 0
      Do cKey=1,p_confs
! Run EFS calculation
        If(nl(cKey)%length.gt.0)Then  ! Only run if there are coord points and neighbour list > 0
          loopID = loopID + 1
          If(M_Loop(eMPI, loopID))Then
            Call calcE_Action(nl, potential, cKey)
          End If
        End If
      End Do
    End If
  End Subroutine calcE_MPI

  Subroutine calcE_Action(nl, potential, cKey)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, fN, pKey, dKey, eKey, nlKey, cKey_Loop
    Logical :: pairPot
! Vars:  Private - Search
    Type(potentialSearchType) :: searchObj
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
!-------------------------
! Potential Keys
!-------------------------
    If(.not.nl(cKey)%keysSet)Then
      Call nlPotentialKeys(nl, potential, cKey)
    End If
!------------------------------
! Init variables
!------------------------------
    nl(cKey)%pairEnergy = 0.0D0
    nl(cKey)%embeddingEnergy = 0.0D0
    nl(cKey)%electronDensity = 0.0D0
    nl(cKey)%atomEnergy = 0.0D0
!------------------------------
! First neighbour list loop
! Pair functions and density functions
!------------------------------
    Do nlKey=1,nl(cKey)%length
!-------------------------------------------------------------------------------
! Get pair key
      pKey = DoubleKey(nl(cKey)%atomA_Type(nlKey), nl(cKey)%atomB_Type(nlKey))
! Loop through pair potentials between atoms
      Do j=1,nl(cKey)%pairKeyArray(pKey,0)
! Set search object
        searchObj%Fn = nl(cKey)%pairKeyArray(pKey,j)
        searchObj%x = nl(cKey)%rD(nlKey)
        yArray = SearchPotential (searchObj, potential)
! Store Pair Energy (Individual)
        nl(cKey)%atomEnergy(nl(cKey)%atomA_ID(nlKey),1) = &
          nl(cKey)%atomEnergy(nl(cKey)%atomA_ID(nlKey),1) + yArray(1)  ! Atom A pair energy
        nl(cKey)%atomEnergy(nl(cKey)%atomB_ID(nlKey),1) = &
          nl(cKey)%atomEnergy(nl(cKey)%atomB_ID(nlKey),1) + yArray(1)  ! Atom B pair energy
! Store energy (total pair)
        nl(cKey)%pairEnergy = nl(cKey)%pairEnergy + yArray(1)
      End Do
!-------------------------------------------------------------------------------
! Electron Density - DENS  (EAM Density) - A electron density at B
      dKey = nl(cKey)%atomA_Type(nlKey)
! Loop through pair potentials between atoms
      Do j=1,nl(cKey)%densityKeyArray(dKey,0,1)
! Set search object
        searchObj%Fn = nl(cKey)%densityKeyArray(dKey,j,1)
        searchObj%x = nl(cKey)%rD(nlKey)
        yArray = SearchPotential (searchObj, potential)
! Store density
        nl(cKey)%electronDensity(nl(cKey)%atomB_ID(nlKey),1) = &
          nl(cKey)%electronDensity(nl(cKey)%atomB_ID(nlKey),1) + yArray(1)
      End Do
!-------------------------------------------------------------------------------
! Electron Density - DENS  (EAM Density) - B electron density at A
      dKey = nl(cKey)%atomB_Type(nlKey)
! Loop through pair potentials between atoms
      Do j=1,nl(cKey)%densityKeyArray(dKey,0,1)
! Set search object
        searchObj%Fn = nl(cKey)%densityKeyArray(dKey,j,1)
        searchObj%x = nl(cKey)%rD(nlKey)
        yArray = SearchPotential (searchObj, potential)
! Store density
        nl(cKey)%electronDensity(nl(cKey)%atomA_ID(nlKey),1) = &
          nl(cKey)%electronDensity(nl(cKey)%atomA_ID(nlKey),1) + yArray(1)
      End Do
    End Do
!------------------------------
! First coords Loop
! Embedding energy
!------------------------------
    Do cKey_Loop=1,nl(cKey)%coordsLength
      eKey = nl(cKey)%labelID(cKey_Loop)
! Loop through embedding energies for atom
      Do j=1,nl(cKey)%embeddingKeyArray(eKey,0,1)
! Set search object
        searchObj%Fn = nl(cKey)%embeddingKeyArray(eKey,j,1)
        searchObj%x = nl(cKey)%electronDensity(cKey_Loop,1)  ! DENS electron density, stored in (cKey_Loop,1)
        yArray = SearchPotential (searchObj, potential)
! Store energy (Individual)
        nl(cKey)%atomEnergy(cKey_Loop,2) = &
          nl(cKey)%atomEnergy(cKey_Loop,2) + yArray(1)  ! Atom A embedding energy
        nl(cKey)%atomEnergy(cKey_Loop,3) = &
          nl(cKey)%atomEnergy(cKey_Loop,1) + nl(cKey)%atomEnergy(cKey_Loop,2)  ! Atom A total energy
! Store energy (total embedding)
        nl(cKey)%embeddingEnergy = nl(cKey)%embeddingEnergy + yArray(1)
      End Do
    End Do
!-------------------------------------------------------------------------------
!
!
!------------------------------
! Total Energy
!------------------------------
    nl(cKey)%totalEnergy = nl(cKey)%pairEnergy + nl(cKey)%embeddingEnergy
!
  End Subroutine calcE_Action
