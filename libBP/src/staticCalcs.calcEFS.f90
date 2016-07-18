! -------------------------------------------------
!  Include File:   Calc energy, force, stress
!
! -------------------------------------------------
  Subroutine calcEFS(nl, potential, cKeyIn, resetForcesIn)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
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
    Call calcEFS_MPI(nl, potential, cKey, resetForces)
  End Subroutine calcEFS

  Subroutine calcEFS_MPI(nl, potential, cKey, resetForces)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
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
            nl(cKey)%forces(i,j) = 0.0D0
          End Do
        End Do
      End If
! Run EFS calculation
      If(nl(cKey)%length.gt.0)Then  ! Only run if there are coord length and neighbour list > 0
        Call calcEFS_Action(nl, potential, cKey)
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
              nl(cKey)%forces(i,j) = 0.0D0
            End Do
          End Do
        End If
! Run EFS calculation
        If(nl(cKey)%length.gt.0)Then  ! Only run if there are coord length and neighbour list > 0
          loopID = loopID + 1
          If(M_Loop(efsMPI, loopID))Then
            Call calcEFS_Action(nl, potential, cKey)
          End If
        End If
      End Do
    End If
  End Subroutine calcEFS_MPI

  Subroutine calcEFS_Action(nl, potential, cKey)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, fN, pKey, dKey, eKey, nlKey, cKey_Loop
! Vars:  Private - Search
    Type(potentialSearchType) :: searchObj
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), Dimension(1:3) :: forceArray
    Real(kind=DoubleReal) :: forceMagnitude
    Real(kind=DoubleReal) :: embeDerivA, embeDerivB, densDerivBA, densDerivAB
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
          nl(cKey)%atomEnergy(nl(cKey)%atomA_ID(nlKey),1) + yArray(1)          ! Atom A pair energy (1)
        nl(cKey)%atomEnergy(nl(cKey)%atomB_ID(nlKey),1) = &
          nl(cKey)%atomEnergy(nl(cKey)%atomB_ID(nlKey),1) + yArray(1)          ! Atom B pair energy (1)
! Store energy (total pair)
        nl(cKey)%pairEnergy = nl(cKey)%pairEnergy + yArray(1)
! Force from pair potential
        forceArray(1) = yArray(2) * nl(cKey)%vecAB(nlKey,1)
        forceArray(2) = yArray(2) * nl(cKey)%vecAB(nlKey,2)
        forceArray(3) = yArray(2) * nl(cKey)%vecAB(nlKey,3)
! Force on atom A
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),1) = &
          nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),1) - forceArray(1)
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),2) = &
          nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),2) - forceArray(2)
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),3) = &
          nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),3) - forceArray(3)
! Force on atom B
        nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),1) = &
          nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),1) + forceArray(1)
        nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),2) = &
          nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),2) + forceArray(2)
        nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),3) = &
          nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),3) + forceArray(3)
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
          nl(cKey)%atomEnergy(cKey_Loop,2) + yArray(1)                          ! Atom A embedding energy (2)
        nl(cKey)%atomEnergy(cKey_Loop,3) = &
          nl(cKey)%atomEnergy(cKey_Loop,1) + nl(cKey)%atomEnergy(cKey_Loop,2)   ! Atom A total energy (3)
! Store energy (total embedding)
        nl(cKey)%embeddingEnergy = nl(cKey)%embeddingEnergy + yArray(1)
      End Do
    End Do
!------------------------------
! Second neighbour list loop
! Embedding function force
!------------------------------
    Do nlKey=1,nl(cKey)%length
!-------------------------------------------------------------------------------
      embeDerivA = 0.0D0
      embeDerivB = 0.0D0
      densDerivBA = 0.0D0
      densDerivAB = 0.0D0
! @Fi(p)/@p
      eKey = nl(cKey)%atomA_Type(nlKey)
      Do j=1,nl(cKey)%embeddingKeyArray(eKey,0,1)
        searchObj%Fn = nl(cKey)%embeddingKeyArray(eKey,j,1)     ! Set search object
        searchObj%x = nl(cKey)%electronDensity(nl(cKey)%atomA_ID(nlKey),1)  ! DENS electron density, stored in (nlKey,1)
        yArray = SearchPotential (searchObj, potential)
        embeDerivA = embeDerivA + yArray(2)
      End Do
! @Fj(p)/@p
      eKey = nl(cKey)%atomB_Type(nlKey)
      Do j=1,nl(cKey)%embeddingKeyArray(eKey,0,1)
        searchObj%Fn = nl(cKey)%embeddingKeyArray(eKey,j,1)     ! Set search object
        searchObj%x = nl(cKey)%electronDensity(nl(cKey)%atomB_ID(nlKey),1)  ! DENS electron density, stored in (nlKey,1)
        yArray = SearchPotential (searchObj, potential)
        embeDerivB = embeDerivB + yArray(2)
      End Do
! @Pij(r)/@r
      dKey = nl(cKey)%atomA_Type(nlKey)
      Do j=1,nl(cKey)%densityKeyArray(dKey,0,1)
        searchObj%Fn = nl(cKey)%densityKeyArray(dKey,j,1)     ! Set search object
        searchObj%x = nl(cKey)%rD(nlKey)  ! DENS electron density, stored in (nlKey,1)
        yArray = SearchPotential (searchObj, potential)
        densDerivAB = densDerivAB + yArray(2)
      End Do
! @Pji(r)/@r
      dKey = nl(cKey)%atomB_Type(nlKey)
      Do j=1,nl(cKey)%densityKeyArray(dKey,0,1)
        searchObj%Fn = nl(cKey)%densityKeyArray(dKey,j,1)     ! Set search object
        searchObj%x = nl(cKey)%rD(nlKey)  ! DENS electron density, stored in (nlKey,1)
        yArray = SearchPotential (searchObj, potential)
        densDerivBA = densDerivBA + yArray(2)
      End Do
! Force from embedding functional
      forceMagnitude = (embeDerivA*densDerivBA+embeDerivB*densDerivAB)
      forceArray(1) = forceMagnitude * nl(cKey)%vecAB(nlKey,1)
      forceArray(2) = forceMagnitude * nl(cKey)%vecAB(nlKey,2)
      forceArray(3) = forceMagnitude * nl(cKey)%vecAB(nlKey,3)
! Force on atom A
      nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),1) = &
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),1) + forceArray(1)
      nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),2) = &
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),2) + forceArray(2)
      nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),3) = &
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),3) + forceArray(3)
! Force on atom B
      nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),1) = &
        nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),1) - forceArray(1)
      nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),2) = &
        nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),2) - forceArray(2)
      nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),3) = &
        nl(cKey)%forcesMD(nl(cKey)%atomB_ID(nlKey),3) - forceArray(3)
    End Do
!-------------------------------------------------------------------------------
!
!
!------------------------------
! Total Energy
!------------------------------
    nl(cKey)%totalEnergy = nl(cKey)%pairEnergy + nl(cKey)%embeddingEnergy
!
  End Subroutine calcEFS_Action
