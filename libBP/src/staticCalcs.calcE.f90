  Subroutine calcE(coords, nl, potential, cKeyIn)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
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
    Call calcE_MPI(coords, nl, potential, cKey)
  End Subroutine calcE

  Subroutine calcE_MPI(coords, nl, potential, cKey)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: loopID
    Type(mpiObj) :: eMPI
! One config
    If(cKey.gt.0)Then
! Run E calculation
      If((coords(cKey)%points.gt.0).and.(nl(cKey)%length.gt.0))Then  ! Only run if there are coord points and neighbour list > 0
        Call calcE_Action(coords, nl, potential, cKey)
      End If
    Else
! Init mpi
      Call m_initMPI(eMPI)
      loopID = 0
      Do cKey=1,p_confs
! Run EFS calculation
        If((coords(cKey)%points.gt.0).and.(nl(cKey)%length.gt.0))Then  ! Only run if there are coord points and neighbour list > 0
          loopID = loopID + 1
          If(M_Loop(eMPI, loopID))Then
            Call calcE_Action(coords, nl, potential, cKey)
          End If
        End If
      End Do
    End If
  End Subroutine calcE_MPI

  Subroutine calcE_Action(coords, nl, potential, cKey)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, fN, pKey, dKey, eKey
    Integer(kind=StandardInteger), Dimension(1:128,0:32) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Dimension(1:128,0:32,1:3) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Dimension(1:128,0:32,1:3) :: embeddingKeyArray  ! Total potentials for each combination pair stored in "0"
    Logical :: pairPot
! Vars:  Private - Search
    Type(potentialSearchType) :: searchObj
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
!-------------------------
! Pair Potential Keys
!-------------------------
! Fill an array that links the double key of atom A and atom B to the (pair) potential key for quicker searching
    pairKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        pairPot = .true.
! If density function or embedding functional
        If(StrMatch(potential%fPotential(fN),"DENS").or.StrMatch(potential%fPotential(fN),"DDEN").or.&
        StrMatch(potential%fPotential(fN),"SDEN"))Then
          pairPot = .false.
        End If
        If(StrMatch(potential%fPotential(fN),"EMBE").or.StrMatch(potential%fPotential(fN),"SEMB").or.&
        StrMatch(potential%fPotential(fN),"DEMB"))Then
          pairPot = .false.
        End If
! If pair pot
        If(pairPot)Then
          pKey = DoubleKey(potential%atomID_A(fN), potential%atomID_B(fN))
          Do j=1,32
            If(pairKeyArray(pKey,j).eq.0)Then
              pairKeyArray(pKey,j) = fN
              pairKeyArray(pKey,0) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!-------------------------
! Density Keys  (EAM, 2BMEAM)
!-------------------------
! Fill an array that links the key of each atom to the (density) potential key for quicker searching
    densityKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        dKey = potential%atomID_A(fN)
! DENS function
        If(StrMatch(potential%fPotential(fN),"DENS"))Then
          Do j=1,32
            If(densityKeyArray(dKey,j,1).eq.0)Then
              densityKeyArray(dKey,j,1) = fN
              densityKeyArray(dKey,0,1) = j
              Exit
            End If
          End Do
        End If
! DDEN function
        If(StrMatch(potential%fPotential(fN),"DDEN"))Then
          Do j=1,32
            If(densityKeyArray(dKey,j,2).eq.0)Then
              densityKeyArray(dKey,j,2) = fN
              densityKeyArray(dKey,0,2) = j
              Exit
            End If
          End Do
        End If
! SDEN function
        If(StrMatch(potential%fPotential(fN),"SDEN"))Then
          Do j=1,32
            If(densityKeyArray(dKey,j,3).eq.0)Then
              densityKeyArray(dKey,j,3) = fN
              densityKeyArray(dKey,0,3) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!-------------------------
! Embedding Keys  (EAM, 2BMEAM)
!-------------------------
! Fill an array that links the key of each atom to the (density) potential key for quicker searching
    embeddingKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        eKey = potential%atomID_A(fN)
! DENS function
        If(StrMatch(potential%fPotential(fN),"EMBE"))Then
          Do j=1,32
            If(embeddingKeyArray(eKey,j,1).eq.0)Then
              embeddingKeyArray(eKey,j,1) = fN
              embeddingKeyArray(eKey,0,1) = j
              Exit
            End If
          End Do
        End If
! DDEN function
        If(StrMatch(potential%fPotential(fN),"DEMB"))Then
          Do j=1,32
            If(embeddingKeyArray(eKey,j,2).eq.0)Then
              embeddingKeyArray(eKey,j,2) = fN
              embeddingKeyArray(eKey,0,2) = j
              Exit
            End If
          End Do
        End If
! SDEN function
        If(StrMatch(potential%fPotential(fN),"SEMB"))Then
          Do j=1,32
            If(embeddingKeyArray(eKey,j,3).eq.0)Then
              embeddingKeyArray(eKey,j,3) = fN
              embeddingKeyArray(eKey,0,3) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!------------------------------
! Init variables
!------------------------------
    coords(cKey)%pairEnergy = 0.0D0
    coords(cKey)%electronDensity = 0.0D0
    coords(cKey)%atomEnergy = 0.0D0
!------------------------------
! First neighbour list loop
! Pair functions and density functions
!------------------------------
    Do i=1,nl(cKey)%length
!-------------------------------------------------------------------------------
! Get pair key
      pKey = DoubleKey(nl(cKey)%atomA_Type(i), nl(cKey)%atomB_Type(i))
! Loop through pair potentials between atoms
      Do j=1,pairKeyArray(pKey,0)
! Set search object
        searchObj%Fn = pairKeyArray(pKey,j)
        searchObj%x = nl(cKey)%rD(i)
        yArray = SearchPotential (searchObj, potential)
! Store Pair Energy (Individual)
        coords(cKey)%atomEnergy(nl(cKey)%atomA_ID(i),1) = &
          coords(cKey)%atomEnergy(nl(cKey)%atomA_ID(i),1) + yArray(1)
        coords(cKey)%atomEnergy(nl(cKey)%atomB_ID(i),1) = &
          coords(cKey)%atomEnergy(nl(cKey)%atomB_ID(i),1) + yArray(1)
! Store energy (total pair)
        coords(cKey)%pairEnergy = coords(cKey)%pairEnergy + yArray(1)
      End Do
!-------------------------------------------------------------------------------
! Electron Density - DENS  (EAM Density) - A electron density at B
      dKey = nl(cKey)%atomA_Type(i)
! Loop through pair potentials between atoms
      Do j=1,densityKeyArray(dKey,0,1)
! Set search object
        searchObj%Fn = densityKeyArray(dKey,j,1)
        searchObj%x = nl(cKey)%rD(i)
        yArray = SearchPotential (searchObj, potential)
! Store density
        coords(cKey)%electronDensity(nl(cKey)%atomB_ID(i),1) = &
          coords(cKey)%electronDensity(nl(cKey)%atomB_ID(i),1) + yArray(1)
      End Do
!-------------------------------------------------------------------------------
! Electron Density - DENS  (EAM Density) - B electron density at A
      dKey = nl(cKey)%atomB_Type(i)
! Loop through pair potentials between atoms
      Do j=1,densityKeyArray(dKey,0,1)
! Set search object
        searchObj%Fn = densityKeyArray(dKey,j,1)
        searchObj%x = nl(cKey)%rD(i)
        yArray = SearchPotential (searchObj, potential)
! Store density
        coords(cKey)%electronDensity(nl(cKey)%atomA_ID(i),1) = &
          coords(cKey)%electronDensity(nl(cKey)%atomA_ID(i),1) + yArray(1)
      End Do
    End Do
!------------------------------
! First coords Loop
! Embedding energy
!------------------------------
    Do i=1,coords(cKey)%points
      eKey = coords(cKey)%labelID(i)
! Loop through embedding energies for atom
      Do j=1,embeddingKeyArray(eKey,0,1)
! Set search object
        searchObj%Fn = embeddingKeyArray(eKey,j,1)
        searchObj%x = coords(cKey)%electronDensity(i,1)  ! DENS electron density, stored in (i,1)
        yArray = SearchPotential (searchObj, potential)
! Store energy (Individual)
        coords(cKey)%atomEnergy(i,2) = &
          coords(cKey)%atomEnergy(i,2) + yArray(1)
        coords(cKey)%atomEnergy(i,3) = &
          coords(cKey)%atomEnergy(i,1) + coords(cKey)%atomEnergy(i,1)
! Store energy (total embedding)
        coords(cKey)%embeddingEnergy = coords(cKey)%embeddingEnergy + yArray(1)
      End Do
    End Do
!-------------------------------------------------------------------------------
!
!
!------------------------------
! Total Energy
!------------------------------
    coords(cKey)%totalEnergy = coords(cKey)%pairEnergy + coords(cKey)%embeddingEnergy
  End Subroutine calcE_Action
