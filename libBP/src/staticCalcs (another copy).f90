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
  Use potentials
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
    Integer(kind=StandardInteger) :: i, j, cKey, k
    Character(Len=16) :: tempLabel
    Character(Len=16), Dimension(1:64) :: tempLabels
! Blank tempLabels array
    tempLabels = BlankStringArray(tempLabels)
! Loop through potential functions
    Do i=1,potential%fCount
      tempLabel = StrToUpper(potential%atomLabel_A(i))
      If(IsBlank(tempLabel))Then
        print *,"Blank"
        potential%atomID_A(i) = 0
      Else
        Do j=1,Size(tempLabels,1)
          If(StrMatch(tempLabel,tempLabels(j)))Then
            potential%atomID_A(i) = j
            Exit
          End If
          If(IfStringEmpty(tempLabels(j)))Then
            tempLabels(j) = tempLabel
            potential%atomID_A(i) = j
            Exit
          End If
        End Do
      End If
      print *,i,"|",tempLabel,"|",potential%atomID_A(i)
      tempLabel = StrToUpper(potential%atomLabel_B(i))
      If(IsBlank(tempLabel))Then
        potential%atomID_B(i) = 0
      Else
        Do j=1,Size(tempLabels,1)
          If(StrMatch(tempLabel,tempLabels(j)))Then
            potential%atomID_B(i) = j
            Exit
          End If
          If(IfStringEmpty(tempLabels(j)))Then
            tempLabels(j) = tempLabel
            potential%atomID_B(i) = j
            Exit
          End If
        End Do
      End If
      print *,i,tempLabel,potential%atomID_B(i)
    End Do
! Loop through coords
    Do cKey=1,Size(coords,1)
      Do i=1,coords(cKey)%points
        tempLabel = StrToUpper(coords(cKey)%label(i))
        Do j=1,Size(tempLabels,1)
          If(StrMatch(tempLabel,tempLabels(j)))Then
            coords(cKey)%labelID(i) = j
            Exit
          End If
          If(IfStringEmpty(tempLabels(j)))Then
            tempLabels(j) = tempLabel
            coords(cKey)%labelID(i) = j
            Exit
          End If
        End Do
      End Do
    End Do
! Store IDs
    k = 0
    Do j=1,64
      If(tempLabels(j)(1:1).eq." ")Then
        Exit
      End If
      k = k + 1
      potential%atomIDs(j) = tempLabels(j)
    End Do
    potential%atomID_Count = k
    Do cKey=1,Size(coords,1)
      Do j=1,64
        coords(cKey)%atomIDs(j) = tempLabels(j)
      End Do
      coords(cKey)%atomID_Count = k
    End Do
! Fill in uKey for potentials
    Do i=1,potential%fCount
      If(potential%atomID_B(i).eq.0)Then
        potential%uKey(i) = potential%atomID_A(i)
      Else
        potential%uKey(i) = DoubleKey(potential%atomID_A(i),potential%atomID_B(i))
      End If
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
    Call addLinePage("Atom Labels and IDs","T")
    Do i=1,32
      If(IfStringEmpty(coords(1)%atomIDs(i)))Then
        Exit
      Else
        Write(tempLine,*) i,coords(1)%atomIDs(i)
        Call addLinePage(tempLine)
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
! Run EFS calculation
      If((coords(cKey)%points.gt.0).and.(nl(cKey)%length.gt.0))Then  ! Only run if there are coord points and neighbour list > 0
        Call calcEFS_Action(coords, nl, potential, cKey)
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
! Run EFS calculation
        If((coords(cKey)%points.gt.0).and.(nl(cKey)%length.gt.0))Then  ! Only run if there are coord points and neighbour list > 0
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
    Integer(kind=StandardInteger) :: i, j, fN, pKey, dKey, eKey
    Integer(kind=StandardInteger), Dimension(1:128,0:32) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Dimension(1:128,0:32,1:3) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Dimension(1:128,0:32,1:3) :: embeddingKeyArray  ! Total potentials for each combination pair stored in "0"
    Logical :: pairPot
! Vars:  Private - Search
    Type(potentialSearchType) :: searchObj
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), Dimension(1:3) :: forceArray
    Real(kind=DoubleReal) :: forceMagnitude
    Real(kind=DoubleReal) :: embeDerivA, embeDerivB, densDerivBA, densDerivAB

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
! Force from pair potential
        forceArray(1) = yArray(2) * nl(cKey)%vecAB(i,1)
        forceArray(2) = yArray(2) * nl(cKey)%vecAB(i,2)
        forceArray(3) = yArray(2) * nl(cKey)%vecAB(i,3)
! Force on atom A
        coords(cKey)%forces(nl(cKey)%atomA_ID(i),1) = &
          coords(cKey)%forces(nl(cKey)%atomA_ID(i),1) - forceArray(1)
        coords(cKey)%forces(nl(cKey)%atomA_ID(i),2) = &
          coords(cKey)%forces(nl(cKey)%atomA_ID(i),2) - forceArray(2)
        coords(cKey)%forces(nl(cKey)%atomA_ID(i),3) = &
          coords(cKey)%forces(nl(cKey)%atomA_ID(i),3) - forceArray(3)
! Force on atom B
        coords(cKey)%forces(nl(cKey)%atomB_ID(i),1) = &
          coords(cKey)%forces(nl(cKey)%atomB_ID(i),1) + forceArray(1)
        coords(cKey)%forces(nl(cKey)%atomB_ID(i),2) = &
          coords(cKey)%forces(nl(cKey)%atomB_ID(i),2) + forceArray(2)
        coords(cKey)%forces(nl(cKey)%atomB_ID(i),3) = &
          coords(cKey)%forces(nl(cKey)%atomB_ID(i),3) + forceArray(3)
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

        !If(i.le.16)Then
          !print *,i, nl(cKey)%rD(i), yArray(1)
        !End If
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
!------------------------------
! Second neighbour list loop
! Embedding function force
!------------------------------
    Do i=1,nl(cKey)%length
!-------------------------------------------------------------------------------
      embeDerivA = 0.0D0
      embeDerivB = 0.0D0
      densDerivBA = 0.0D0
      densDerivAB = 0.0D0
! @Fi(p)/@p
      eKey = nl(cKey)%atomA_Type(i)
      Do j=1,embeddingKeyArray(eKey,0,1)
        searchObj%Fn = embeddingKeyArray(eKey,j,1)     ! Set search object
        searchObj%x = coords(cKey)%electronDensity(nl(cKey)%atomA_ID(i),1)  ! DENS electron density, stored in (i,1)
        yArray = SearchPotential (searchObj, potential)
        embeDerivA = embeDerivA + yArray(2)
      End Do
! @Fj(p)/@p
      eKey = nl(cKey)%atomB_Type(i)
      Do j=1,embeddingKeyArray(eKey,0,1)
        searchObj%Fn = embeddingKeyArray(eKey,j,1)     ! Set search object
        searchObj%x = coords(cKey)%electronDensity(nl(cKey)%atomB_ID(i),1)  ! DENS electron density, stored in (i,1)
        yArray = SearchPotential (searchObj, potential)
        embeDerivB = embeDerivB + yArray(2)
      End Do
! @Pij(r)/@r
      dKey = nl(cKey)%atomA_Type(i)
      Do j=1,densityKeyArray(dKey,0,1)
        searchObj%Fn = densityKeyArray(dKey,j,1)     ! Set search object
        searchObj%x = nl(cKey)%rD(i)  ! DENS electron density, stored in (i,1)
        yArray = SearchPotential (searchObj, potential)
        densDerivAB = densDerivAB + yArray(2)
      End Do
! @Pji(r)/@r
      dKey = nl(cKey)%atomB_Type(i)
      Do j=1,densityKeyArray(dKey,0,1)
        searchObj%Fn = densityKeyArray(dKey,j,1)     ! Set search object
        searchObj%x = nl(cKey)%rD(i)  ! DENS electron density, stored in (i,1)
        yArray = SearchPotential (searchObj, potential)
        densDerivBA = densDerivBA + yArray(2)
      End Do
! Force from embedding functional
      forceMagnitude = (embeDerivA*densDerivBA+embeDerivB*densDerivAB)
      forceArray(1) = forceMagnitude * nl(cKey)%vecAB(i,1)
      forceArray(2) = forceMagnitude * nl(cKey)%vecAB(i,2)
      forceArray(3) = forceMagnitude * nl(cKey)%vecAB(i,3)
! Force on atom A
      coords(cKey)%forces(nl(cKey)%atomA_ID(i),1) = &
        coords(cKey)%forces(nl(cKey)%atomA_ID(i),1) + forceArray(1)
      coords(cKey)%forces(nl(cKey)%atomA_ID(i),2) = &
        coords(cKey)%forces(nl(cKey)%atomA_ID(i),2) + forceArray(2)
      coords(cKey)%forces(nl(cKey)%atomA_ID(i),3) = &
        coords(cKey)%forces(nl(cKey)%atomA_ID(i),3) + forceArray(3)
! Force on atom B
      coords(cKey)%forces(nl(cKey)%atomB_ID(i),1) = &
        coords(cKey)%forces(nl(cKey)%atomB_ID(i),1) - forceArray(1)
      coords(cKey)%forces(nl(cKey)%atomB_ID(i),2) = &
        coords(cKey)%forces(nl(cKey)%atomB_ID(i),2) - forceArray(2)
      coords(cKey)%forces(nl(cKey)%atomB_ID(i),3) = &
        coords(cKey)%forces(nl(cKey)%atomB_ID(i),3) - forceArray(3)
    End Do
!-------------------------------------------------------------------------------
!
!
!------------------------------
! Total Energy
!------------------------------
    coords(cKey)%totalEnergy = coords(cKey)%pairEnergy + coords(cKey)%embeddingEnergy



! Print

    print *,"NL complete: ",nl(cKey)%length
    print *,coords(cKey)%pairEnergy,(coords(cKey)%pairEnergy/coords(cKey)%points)
    print *,coords(cKey)%embeddingEnergy,(coords(cKey)%embeddingEnergy/coords(cKey)%points)
    print *,coords(cKey)%totalEnergy,(coords(cKey)%totalEnergy/coords(cKey)%points)


    Do i=1,coords(cKey)%points
      !print *,coords(cKey)%forces(i,1),coords(cKey)%forces(i,2),coords(cKey)%forces(i,3)
    End Do


  End Subroutine calcEFS_Action


! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------













End Module staticCalcs
