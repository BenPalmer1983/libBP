! -------------------------------------------------
!  Include File:   Calc energy, force - Individual Atoms
!                  Used by geom optimise subroutines
!
! -------------------------------------------------
  Subroutine calcEF(nl, potential, cKeyIn, atomKeyIn, energyOnlyIn, resetForcesIn)      ! Make a data type to clean this up
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger), Optional :: cKeyIn
    Integer(kind=StandardInteger), Optional :: atomKeyIn
    Logical, Optional :: energyOnlyIn
    Logical, Optional :: resetForcesIn
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey, atomKey
    Logical :: energyOnly
    Logical :: resetForces
! Optional Arguments
    resetForces = .false.
    cKey = 0
    atomKey = 0
    energyOnly = .true.
    If(Present(resetForcesIn))Then
      resetForces = resetForcesIn
    End If
    If(Present(energyOnlyIn))Then
      energyOnly = energyOnlyIn
    End If
    If(Present(cKeyIn))Then
      cKey = cKeyIn
    End If
    If(Present(atomKeyIn))Then
      atomKey = atomKeyIn
    End If
! Calc EF
    Call calcEF_LKeys(nl, potential, cKey, atomKey, energyOnly, resetForces) ! Loop through cKey config keys
  End Subroutine calcEF


  Subroutine calcEF_LKeys(nl, potential, cKey, atomKey, energyOnly, resetForces)
! Loop through cKey config keys
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: atomKey
    Logical :: resetForces
    Logical :: energyOnly
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey_Loop, i
! Loop through nl keys (or specific key)
    If(cKey.eq.0)Then
      Do cKey_Loop=1,size(nl)
        Call calcEF_LAtoms(nl, potential, cKey_Loop, atomKey, energyOnly, resetForces) ! Loop through atomKey atom keys
      End Do
    Else
      Call calcEF_LAtoms(nl, potential, cKey, atomKey, energyOnly, resetForces) ! Loop through atomKey atom keys
    End If
  End Subroutine calcEF_LKeys

  Subroutine calcEF_LAtoms(nl, potential, cKey, atomKey, energyOnly, resetForces)
! Loop through atom IDs
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger), Optional :: cKey
    Integer(kind=StandardInteger), Optional :: atomKey
    Logical :: resetForces
    Logical :: energyOnly
! Vars:  Private
    Integer(kind=StandardInteger) :: atomKey_Loop, i
! Reset forces
    If(atomKey.eq.0)Then
      Do atomKey_Loop=1,nl(cKey)%coordsLength
        Call calcEF_Action(nl, potential, cKey, atomKey_Loop, energyOnly, resetForces)
      End Do
    Else
      Call calcEF_Action(nl, potential, cKey, atomKey, energyOnly, resetForces)
    End If
  End Subroutine calcEF_LAtoms


! -----------------------------------------
!



  Subroutine calcEF_Action(nl, potential, cKey, atomKey, energyOnly, resetForces)
! Select energy only or energy + force
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger), Optional :: cKey
    Integer(kind=StandardInteger), Optional :: atomKey
    Logical :: resetForces
    Logical :: energyOnly
! Choose subroutine
    If(energyOnly)Then
      Call calcEF_Action_E(nl, potential, cKey, atomKey)
    Else

    End If
  End Subroutine calcEF_Action


! -----------------------------------------
!    Energy Only
! -----------------------------------------

  Subroutine calcEF_Action_E(nl, potential, cKey, atomKey)
! Calculates energy of a single atom
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: atomKey, atomKeyBs
    Logical :: resetForces
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, k, fN, pKey, dKey, eKey, nlKey, cKey_Loop
    Logical :: pairPot
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
      Call nlPotentialKeys_opt(nl, potential, cKey)
    End If
!------------------------------
! Init variables
!------------------------------
    nl(cKey)%pairEnergy = 0.0D0
    nl(cKey)%embeddingEnergy = 0.0D0
    Do j=1,3
      nl(cKey)%electronDensity(atomKey,j) = 0.0D0
      nl(cKey)%atomEnergy(atomKey,j) = 0.0D0
    End Do
! reset electron density at surrounding atoms
    Do i=1,nl(cKey)%nlOverKey(atomKey,0)
      nlKey = nl(cKey)%nlOverKey(atomKey,i)
      atomKeyBs = nl(cKey)%atomB_Type(nlKey)
      Do j=1,3
        nl(cKey)%electronDensity(atomKeyBs,j) = 0.0D0
      End Do
    End Do
!------------------------------
! First neighbour list loop
! Pair functions and density functions
!------------------------------
    Do i=1,nl(cKey)%nlOverKey(atomKey,0)
      nlKey = nl(cKey)%nlOverKey(atomKey,i)
!-------------------------------------------------------------------------------
! Get pair key
      pKey = DoubleKey(nl(cKey)%atomA_Type(nlKey), nl(cKey)%atomB_Type(nlKey))
! Loop through pair potentials between atoms
      Do j=1,nl(cKey)%pairKeyArray(pKey,0)
! Set search object
        searchObj%Fn = nl(cKey)%pairKeyArray(pKey,j)
        searchObj%x = nl(cKey)%rD(nlKey)
        yArray = SearchPotential (searchObj, potential)
! Store Pair Energy - Atom A only
        nl(cKey)%atomEnergy(nl(cKey)%atomA_ID(nlKey),1) = &
          nl(cKey)%atomEnergy(nl(cKey)%atomA_ID(nlKey),1) + yArray(1)           ! Atom A pair energy
! Store energy (total pair)
        nl(cKey)%pairEnergy = nl(cKey)%pairEnergy + yArray(1)
      End Do
!-------------------------------------------------------------------------------
! Electron Density - DENS  (EAM Density) - Sum of Bs electron density at A
      dKey = nl(cKey)%atomB_Type(nlKey)
! Loop through pair potentials between atoms
      Do j=1,nl(cKey)%densityKeyArray(dKey,0,1)
! Set search object
        searchObj%Fn = nl(cKey)%densityKeyArray(dKey,j,1)
        searchObj%x = nl(cKey)%rD(nlKey)
        yArray = SearchPotential (searchObj, potential)
! Store density
        nl(cKey)%electronDensity(atomKey,1) = &
          nl(cKey)%electronDensity(atomKey,1) + yArray(1)
      End Do
    End Do
!------------------------------
! First coords Loop
! Embedding energy
!------------------------------
    eKey = nl(cKey)%labelID(atomKey)
! Loop through embedding energy functions for atom
    Do j=1,nl(cKey)%embeddingKeyArray(eKey,0,1)
! Set search object
      searchObj%Fn = nl(cKey)%embeddingKeyArray(eKey,j,1)
      searchObj%x = nl(cKey)%electronDensity(atomKey,1)  ! DENS electron density, stored in (i,1)
      yArray = SearchPotential (searchObj, potential)
! Store energy (Individual)
      nl(cKey)%atomEnergy(atomKey,2) = &
        nl(cKey)%atomEnergy(atomKey,2) + yArray(1)                      ! Atom A embedding energy
      nl(cKey)%atomEnergy(atomKey,3) = &
        nl(cKey)%atomEnergy(atomKey,1) + nl(cKey)%atomEnergy(atomKey,2) ! Atom A total energy
! Store energy (total embedding)
      nl(cKey)%embeddingEnergy = nl(cKey)%embeddingEnergy + yArray(1)
    End Do
!-------------------------------------------------------------------------------
!
!
!------------------------------
! Total Energy
!------------------------------
    nl(cKey)%totalEnergy = nl(cKey)%pairEnergy + nl(cKey)%embeddingEnergy

!
  End Subroutine calcEF_Action_E





  Subroutine calcEF_Action_EF(nl, potential, cKey, atomKey, resetForces)
! Calculates energy/force of a collection of atoms
! **** not working 18/07/2016 ****
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: atomKey, atomKeyBs
    Logical :: resetForces
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, k, fN, pKey, dKey, eKey, nlKey, cKey_Loop
    Logical :: pairPot
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
      Call nlPotentialKeys_opt(nl, potential, cKey)
    End If
!------------------------------
! Init variables
!------------------------------
    nl(cKey)%pairEnergy = 0.0D0
    nl(cKey)%embeddingEnergy = 0.0D0
    Do j=1,3
      nl(cKey)%electronDensity(atomKey,j) = 0.0D0
      nl(cKey)%atomEnergy(atomKey,j) = 0.0D0
    End Do
! reset electron density at surrounding atoms
    Do i=1,nl(cKey)%nlOverKey(atomKey,0)
      nlKey = nl(cKey)%nlOverKey(atomKey,i)
      atomKeyBs = nl(cKey)%atomB_Type(nlKey)
      Do j=1,3
        nl(cKey)%atomEnergy(atomKeyBs,j) = 0.0D0
      End Do
    End Do
! Calculate electron density to use for small perturbations of atoms
    !If(.not.fixedDensityCalculated)Then
      !electronDensityFixed
    !End If

!------------------------------
! First neighbour list loop
! Pair functions and density functions
!------------------------------
    Do i=1,nl(cKey)%nlOverKey(atomKey,0)
      nlKey = nl(cKey)%nlOverKey(atomKey,i)
!-------------------------------------------------------------------------------
! Get pair key
      pKey = DoubleKey(nl(cKey)%atomA_Type(nlKey), nl(cKey)%atomB_Type(nlKey))
! Loop through pair potentials between atoms
      Do j=1,nl(cKey)%pairKeyArray(pKey,0)
! Set search object
        searchObj%Fn = nl(cKey)%pairKeyArray(pKey,j)
        searchObj%x = nl(cKey)%rD(nlKey)
        yArray = SearchPotential (searchObj, potential)
! Store Pair Energy - Atom A only
        nl(cKey)%atomEnergy(nl(cKey)%atomA_ID(nlKey),1) = &
          nl(cKey)%atomEnergy(nl(cKey)%atomA_ID(nlKey),1) + yArray(1)           ! Atom A pair energy
! Store energy (total pair)
        nl(cKey)%pairEnergy = nl(cKey)%pairEnergy + yArray(1)
! Force from pair potential
        forceArray(1) = yArray(2) * nl(cKey)%vecAB(nlKey,1)
        forceArray(2) = yArray(2) * nl(cKey)%vecAB(nlKey,2)
        forceArray(3) = yArray(2) * nl(cKey)%vecAB(nlKey,3)
! Force on atom A Only
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),1) = &
          nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),1) - forceArray(1)
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),2) = &
          nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),2) - forceArray(2)
        nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),3) = &
          nl(cKey)%forcesMD(nl(cKey)%atomA_ID(nlKey),3) - forceArray(3)
      End Do
!-------------------------------------------------------------------------------
! Electron Density - DENS  (EAM Density) - Sum of Bs electron density at A
      dKey = nl(cKey)%atomB_Type(nlKey)
! Loop through pair potentials between atoms
      Do j=1,nl(cKey)%densityKeyArray(dKey,0,1)
! Set search object
        searchObj%Fn = nl(cKey)%densityKeyArray(dKey,j,1)
        searchObj%x = nl(cKey)%rD(nlKey)
        yArray = SearchPotential (searchObj, potential)
! Store density
        nl(cKey)%electronDensity(atomKey,1) = &
          nl(cKey)%electronDensity(atomKey,1) + yArray(1)
      End Do
    End Do
!------------------------------
! Second neighbour list loop
! Pair functions and density functions
!------------------------------
    Do i=1,nl(cKey)%nlOverKey(atomKey,0)  ! Loop through B atoms that surround the A atom "atomKey"
      nlKey = nl(cKey)%nlOverKey(atomKey,i)
      atomKeyBs = nl(cKey)%atomB_Type(nlKey)
      Do k=1,nl(cKey)%nlOverKey(atomKeyBs,0)  ! Loop through C atoms that surround each B atom "atomKeyBs"
        nlKey = nl(cKey)%nlOverKey(atomKey,i)
!-------------------------------------------------------------------------------
! Electron Density - DENS  (EAM Density) - Sum of Cs electron density at Bs (Bs surround the original A atom)
        dKey = nl(cKey)%atomB_Type(nlKey)
! Loop through pair potentials between atoms
        Do j=1,nl(cKey)%densityKeyArray(dKey,0,1)
! Set search object
          searchObj%Fn = nl(cKey)%densityKeyArray(dKey,j,1)
          searchObj%x = nl(cKey)%rD(nlKey)
          yArray = SearchPotential (searchObj, potential)
! Store density
          nl(cKey)%electronDensity(atomKeyBs,1) = &
            nl(cKey)%electronDensity(atomKeyBs,1) + yArray(1)
        End Do
      End Do
    End Do
!------------------------------
! First coords Loop
! Embedding energy
!------------------------------
    eKey = nl(cKey)%labelID(atomKey)
! Loop through embedding energy functions for atom
    Do j=1,nl(cKey)%embeddingKeyArray(eKey,0,1)
! Set search object
      searchObj%Fn = nl(cKey)%embeddingKeyArray(eKey,j,1)
      searchObj%x = nl(cKey)%electronDensity(atomKey,1)  ! DENS electron density, stored in (i,1)
      yArray = SearchPotential (searchObj, potential)
! Store energy (Individual)
      nl(cKey)%atomEnergy(atomKey,2) = &
        nl(cKey)%atomEnergy(atomKey,2) + yArray(1)                      ! Atom A embedding energy
      nl(cKey)%atomEnergy(atomKey,3) = &
        nl(cKey)%atomEnergy(atomKey,1) + nl(cKey)%atomEnergy(atomKey,2) ! Atom A total energy
! Store energy (total embedding)
      nl(cKey)%embeddingEnergy = nl(cKey)%embeddingEnergy + yArray(1)
    End Do
!------------------------------
! Third neighbour list loop
! Embedding function force
!------------------------------
    Do i=1,nl(cKey)%nlOverKey(atomKey,0)
      nlKey = nl(cKey)%nlOverKey(atomKey,i)
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
    End Do
!-------------------------------------------------------------------------------
!
!
!------------------------------
! Total Energy
!------------------------------
    nl(cKey)%totalEnergy = nl(cKey)%pairEnergy + nl(cKey)%embeddingEnergy
  End Subroutine calcEF_Action_EF





! -----------------------------------------
!    Calculate Densities Only
! -----------------------------------------

  Subroutine calcDensities(nl, potential, cKey_In)
! Calculates Electron density at each atom
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger), Optional :: cKey_In
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey, cKey_Loop
! Optional arguments
    cKey = 0
    If(Present(cKey_In))Then
      cKey = cKey_In
    End If
    If(cKey.eq.0)Then
      Do cKey_Loop=1,size(nl)
        Call calcDensities_Action(nl, potential, cKey_Loop) !
      End Do
    Else
      Call calcDensities_Action(nl, potential, cKey)        !
    End If
  End Subroutine calcDensities
!-----------------------------------------------------------
  Subroutine calcDensities_Action(nl, potential, cKey)
! Calculates Electron density at each atom
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_opt), Dimension(:) :: nl
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
      Call nlPotentialKeys_opt(nl, potential, cKey)  ! Make keys if not set
    End If
!------------------------------
! Init variables
!------------------------------
    nl(cKey)%fixedDensityCalculated = .true.
    nl(cKey)%electronDensityFixed = 0.0D0
!------------------------------
! First neighbour list loop
! Pair functions and density functions
!------------------------------
    Do nlKey=1,nl(cKey)%length
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
        nl(cKey)%electronDensityFixed(nl(cKey)%atomB_ID(nlKey),1) = &
          nl(cKey)%electronDensityFixed(nl(cKey)%atomB_ID(nlKey),1) + yArray(1)
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
        nl(cKey)%electronDensityFixed(nl(cKey)%atomA_ID(nlKey),1) = &
          nl(cKey)%electronDensityFixed(nl(cKey)%atomA_ID(nlKey),1) + yArray(1)
      End Do
    End Do
  End Subroutine calcDensities_Action










!-------------------------------------------------------------------------------
