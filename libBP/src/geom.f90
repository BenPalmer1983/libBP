! --------------------------------------------------------------!
! Geometry module
! geomTypes, geom
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Makes configuration of atoms
! Calculates neighbour list for a collection of atoms
!
! ----------------------------------------
! Updated: 4th November 2015
! ----------------------------------------

Module geomTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
  Integer(kind=StandardInteger), Parameter :: p_confs = 10
  Integer(kind=StandardInteger), Parameter :: p_cMax = 1024        ! Coords max length
  Integer(kind=StandardInteger), Parameter :: p_cuMax = 128        ! Coords unit type max length
  Integer(kind=StandardInteger), Parameter :: p_nlMax = 50000      ! Neighbour list length - no longer used
! Make private
  Private
! Public Variables and Parameters
  Public :: p_confs, p_cuMax, p_cMax, p_nlMax
! Public derived types
  Public :: coordsUnitType, coordsType, nlType, nlType_Opt

  Type :: coordsUnitType
    Integer(kind=StandardInteger) :: points = 0
    Integer(kind=StandardInteger) :: xCopy = 1
    Integer(kind=StandardInteger) :: yCopy = 1
    Integer(kind=StandardInteger) :: zCopy = 1
    Real(kind=DoubleReal) :: aLat = 1.0D0
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Character(Len=16), Dimension(1:p_cuMax) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:p_cuMax) :: labelID
    Real(kind=DoubleReal), Dimension(1:p_cuMax,1:3) :: unitCoords   ! Unit coords - fractional
    Real(kind=DoubleReal), Dimension(1:p_cuMax,1:3) :: unitForces   ! Unit coords - fractional
  End Type coordsUnitType

  Type :: coordsType
    Integer(kind=StandardInteger) :: length = 0
    Real(kind=DoubleReal) :: aLat = 0.0D0
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Character(Len=16), Dimension(1:p_cMax) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:p_cMax) :: labelID
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: coords       ! coords (fractional)
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: forces       ! forces
    Integer(kind=StandardInteger) :: IDcount
    Character(Len=16), Dimension(1:p_cuMax) :: atomIDs
    Integer(kind=StandardInteger) :: atomID_Count
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: electronDensity   ! Moved to nl object
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: atomEnergy        ! Moved to nl object
    Real(kind=DoubleReal) :: pairEnergy
    Real(kind=DoubleReal) :: embeddingEnergy
    Real(kind=DoubleReal) :: totalEnergy
  End Type coordsType

  Type :: nlType
! 6.2MB per neighbour list
    Integer(kind=StandardInteger) :: length = 0
    Integer(kind=StandardInteger) :: arrayLength = 0
! Coord details
    Integer(kind=StandardInteger) :: coordsLength
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Character(Len=16), Dimension(1:p_cMax) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:p_cMax) :: labelID
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: coords   ! Initial Coords
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: forces   ! Initial Forces
    Character(Len=16), Dimension(1:p_cuMax) :: atomIDs
    Integer(kind=StandardInteger) :: atomID_Count
! Potential Keys
    Logical :: keysSet = .false.
    !Integer(kind=StandardInteger), Dimension(1:128,0:32) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    !Integer(kind=StandardInteger), Dimension(1:128,0:32,1:3) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    !Integer(kind=StandardInteger), Dimension(1:128,0:32,1:3) :: embeddingKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: embeddingKeyArray  ! Total potentials for each combination pair stored in "0"
! Atom details - EFS Calculation
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: electronDensity   ! Moved to nl object
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: atomEnergy        ! Moved to nl object
    Real(kind=DoubleReal) :: pairEnergy
    Real(kind=DoubleReal) :: embeddingEnergy
    Real(kind=DoubleReal) :: totalEnergy
! MD forces/coord - updated as they change
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: coordsMD   ! Initial Coords
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: forcesMD   ! Initial Forces
! Neighbour List Details
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomB_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_Type
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomB_Type
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomPairKey
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: inCell
! Position Details
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: rD
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: vecAB
! Subcell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: subCell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: cell
! Subcell
    Integer(kind=StandardInteger) :: scCount = 0
    Real(kind=DoubleReal) :: scALat
    Integer(kind=StandardInteger) :: maxAtomsPerSC = 0
! Misc
    Real(kind=DoubleReal) :: rMin
    Real(kind=DoubleReal) :: rMax
    Real(kind=DoubleReal) :: rVerlet
    Real(kind=DoubleReal) :: totalRD
    Real(kind=DoubleReal) :: totalRDSq
! NL Update Parameters
    Real(kind=DoubleReal) :: aLat_Original    ! Lattice parameter when NL built
  End Type nlType

  Type :: nlType_Opt
! 6.2MB per neighbour list
    Integer(kind=StandardInteger) :: length = 0
    Integer(kind=StandardInteger) :: arrayLength = 0
! Coord details
    Integer(kind=StandardInteger) :: coordsLength
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Character(Len=16), Dimension(1:p_cMax) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:p_cMax) :: labelID
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: coords   ! Initial Coords
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: forces   ! Initial Forces
    Character(Len=16), Dimension(1:p_cuMax) :: atomIDs
    Integer(kind=StandardInteger) :: atomID_Count
! Potential Keys
    Logical :: keysSet = .false.
    !Integer(kind=StandardInteger), Dimension(1:128,0:32) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    !Integer(kind=StandardInteger), Dimension(1:128,0:32,1:3) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    !Integer(kind=StandardInteger), Dimension(1:128,0:32,1:3) :: embeddingKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: pairKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: densityKeyArray  ! Total potentials for each combination pair stored in "0"
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: embeddingKeyArray  ! Total potentials for each combination pair stored in "0"
! Atom details - EFS Calculation
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: electronDensity
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: atomEnergy             ! pair, embedding, total
    Logical :: fixedDensityCalculated = .false.
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: electronDensityFixed   ! used to approximate electron density without recalculating for small perturbations of atoms
    Real(kind=DoubleReal) :: pairEnergy
    Real(kind=DoubleReal) :: embeddingEnergy
    Real(kind=DoubleReal) :: totalEnergy
! MD forces/coord - updated as they change
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: coordsMD   ! Initial Coords
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: forcesMD   ! Initial Forces
! Neighbour List Details
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: nlOverKey   ! Used to isolate individual atoms and their local neighbour list
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomB_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_Type
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomB_Type
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomPairKey
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: inCell
! Position Details
    Real(kind=DoubleReal), Allocatable, Dimension(:) :: rD
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: vecAB
! Subcell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: subCell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: cell
! Subcell
    Integer(kind=StandardInteger) :: scCount = 0
    Real(kind=DoubleReal) :: scALat
    Integer(kind=StandardInteger) :: maxAtomsPerSC = 0
! Misc
    Real(kind=DoubleReal) :: rMin
    Real(kind=DoubleReal) :: rMax
    Real(kind=DoubleReal) :: rVerlet
    Real(kind=DoubleReal) :: totalRD
    Real(kind=DoubleReal) :: totalRDSq
! NL Update Parameters
    Real(kind=DoubleReal) :: aLat_Original    ! Lattice parameter when NL built
  End Type nlType_Opt

End Module geomTypes


Module geom
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
  Use basicMaths
  Use keysMod
  Use rng
  Use rngDist
  Use linearAlgebra
  Use coordFunctions
  Use geomTypes
! Force declaration of all variables
  Implicit None
! Public variables
!  Type(nlType) :: nl
!  Real(kind=DoubleReal),Dimension(1:3) :: coords
! Make private
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: initUnitCoords
  Public :: initCoords
  Public :: standardCoords
  Public :: expandUnitCoords
  Public :: heatCoords
  Public :: zeroForces
  Public :: makeNL
  Public :: updateNL
  Public :: makeNL_opt
  Public :: printCoords
  Public :: printNLSummary
  Public :: Tensor_Homogeneous
  Public :: Tensor_Orthorhombic
  Public :: Tensor_Tetragonal


! Interfaces
!  Interface makeNL
!    Module Procedure makeNL_A, makeNL_B
!  End Interface makeNL
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------


! ------------------------------------------------------------
!               Coordinates
! ------------------------------------------------------------

  Subroutine initUnitCoords(coords)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(coordsUnitType), Dimension(:) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, cKey
! Initialise data type
    Do cKey =1,size(coords)
      coords(cKey)%xCopy = 1
      coords(cKey)%yCopy = 1
      coords(cKey)%zCopy = 1
      coords(cKey)%aLat = 1.0D0
      Do i=1,size(coords(cKey)%label,1)
        coords(cKey)%label(i) = BlankString(coords(cKey)%label(i))
        coords(cKey)%labelID(i) = 0
        Do j=1,3
          coords(cKey)%unitCoords(i,j) = 0.0D0
          coords(cKey)%unitForces(i,j) = 0.0D0
        End Do
      End Do
      Do i=1,3
        Do j=1,3
          If(i.eq.j)Then
            coords(cKey)%unitCell(i,j) = 1.0D0
          Else
            coords(cKey)%unitCell(i,j) = 0.0D0
          End If
        End Do
      End Do
      coords(cKey)%points = 0
    End Do
  End Subroutine initUnitCoords

  Subroutine initCoords(coords)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(coordsType), Dimension(:) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, cKey
! Initialise data type
    Do cKey=1,size(coords)
      Do i=1,size(coords(cKey)%label,1)
        coords(cKey)%label(i) = BlankString(coords(cKey)%label(i))
        coords(cKey)%labelID(i) = 0
        Do j=1,3
          coords(cKey)%coords(i,j) = 0.0D0
          coords(cKey)%forces(i,j) = 0.0D0
        End Do
      End Do
      coords(cKey)%length = 0
    End Do
  End Subroutine initCoords

  Subroutine standardCoords(typeCell, coords, cKey, atomLabels, atomIDs_in)
! Loads standard unit coords
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Character(*) :: typeCell
    Type(coordsUnitType), Dimension(:) :: coords
    Integer(kind=StandardInteger) :: cKey
    Character(*), Dimension(:) :: atomLabels
    Integer(kind=StandardInteger), Dimension(:), Optional :: atomIDs_in
! Vars:  Private
    Character(Len(typeCell)) :: typeCellUpper
    Integer(kind=StandardInteger), Dimension(1:size(atomLabels,1)) :: atomIDs
! Optional argument
    atomIDs = 0
    If(Present(atomIDs_in))Then
      atomIDs = atomIDs_in
    End If
! Read in standard type
    typeCellUpper = StrToUpper(typeCell)
    If(typeCellUpper(1:3).eq."FCC")Then
! Unit cell
      coords(cKey)%unitCell(1,1) = 1.0D0
      coords(cKey)%unitCell(1,2) = 0.0D0
      coords(cKey)%unitCell(1,3) = 0.0D0
      coords(cKey)%unitCell(2,1) = 0.0D0
      coords(cKey)%unitCell(2,2) = 1.0D0
      coords(cKey)%unitCell(2,3) = 0.0D0
      coords(cKey)%unitCell(3,1) = 0.0D0
      coords(cKey)%unitCell(3,2) = 0.0D0
      coords(cKey)%unitCell(3,3) = 1.0D0
! Points
      coords(cKey)%points = 4
      coords(cKey)%label(1) = atomLabels(1)
      coords(cKey)%labelID(1) = atomIDs(1)
      coords(cKey)%unitCoords(1,1) = 0.0D0
      coords(cKey)%unitCoords(1,2) = 0.0D0
      coords(cKey)%unitCoords(1,3) = 0.0D0
      coords(cKey)%label(2) = atomLabels(2)
      coords(cKey)%labelID(2) = atomIDs(2)
      coords(cKey)%unitCoords(2,1) = 0.5D0
      coords(cKey)%unitCoords(2,2) = 0.5D0
      coords(cKey)%unitCoords(2,3) = 0.0D0
      coords(cKey)%label(3) = atomLabels(3)
      coords(cKey)%labelID(3) = atomIDs(3)
      coords(cKey)%unitCoords(3,1) = 0.5D0
      coords(cKey)%unitCoords(3,2) = 0.0D0
      coords(cKey)%unitCoords(3,3) = 0.5D0
      coords(cKey)%label(4) = atomLabels(4)
      coords(cKey)%labelID(4) = atomIDs(4)
      coords(cKey)%unitCoords(4,1) = 0.0D0
      coords(cKey)%unitCoords(4,2) = 0.5D0
      coords(cKey)%unitCoords(4,3) = 0.5D0
    End If
    If(typeCellUpper(1:3).eq."BCC")Then
! Unit cell
      coords(cKey)%unitCell(1,1) = 1.0D0
      coords(cKey)%unitCell(1,2) = 0.0D0
      coords(cKey)%unitCell(1,3) = 0.0D0
      coords(cKey)%unitCell(2,1) = 0.0D0
      coords(cKey)%unitCell(2,2) = 1.0D0
      coords(cKey)%unitCell(2,3) = 0.0D0
      coords(cKey)%unitCell(3,1) = 0.0D0
      coords(cKey)%unitCell(3,2) = 0.0D0
      coords(cKey)%unitCell(3,3) = 1.0D0
! Points
      coords(cKey)%points = 2
      coords(cKey)%label(1) = atomLabels(1)
      coords(cKey)%labelID(1) = atomIDs(1)
      coords(cKey)%unitCoords(1,1) = 0.0D0
      coords(cKey)%unitCoords(1,2) = 0.0D0
      coords(cKey)%unitCoords(1,3) = 0.0D0
      coords(cKey)%label(2) = atomLabels(2)
      coords(cKey)%labelID(2) = atomIDs(2)
      coords(cKey)%unitCoords(2,1) = 0.5D0
      coords(cKey)%unitCoords(2,2) = 0.5D0
      coords(cKey)%unitCoords(2,3) = 0.5D0
    End If
  End Subroutine standardCoords

  Subroutine expandUnitCoords(coordsUnit, coords)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsUnitType), Dimension(:) :: coordsUnit
    Type(coordsType), Dimension(:) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, k, m, n, cKey, cKeyE
    Real(kind=DoubleReal), Dimension(1:3) :: xCoords
    Real(kind=DoubleReal) :: xMult, yMult, zMult
! Loop through configs
    Do cKey=1,size(coordsUnit,1)
      If((cKey.le.size(coords,1)).and.(coordsUnit(cKey)%points.gt.0))Then
        cKeyE = cKey
! Adjust unit vector
        xMult = coordsUnit(cKey)%xCopy/coordsUnit(cKey)%xCopy
        yMult = coordsUnit(cKey)%yCopy/coordsUnit(cKey)%xCopy
        zMult = coordsUnit(cKey)%zCopy/coordsUnit(cKey)%xCopy
        Do i=1,3
          coords%unitCell(1,i) = xMult*coords%unitCell(1,i)
          coords%unitCell(2,i) = xMult*coords%unitCell(2,i)
          coords%unitCell(3,i) = xMult*coords%unitCell(3,i)
        End Do
        m = 0
        Do i=1,coordsUnit(cKey)%xCopy
          Do j=1,coordsUnit(cKey)%yCopy
            Do k=1,coordsUnit(cKey)%zCopy
              Do n=1,coordsUnit(cKey)%points
                m = m + 1
                coords(cKeyE)%label(m) = coordsUnit(cKey)%label(n)
                coords(cKeyE)%labelID(m) = coordsUnit(cKey)%labelID(n)
! Fractional coords
                coords(cKeyE)%coords(m,1) = (coordsUnit(cKey)%unitCoords(n,1)+1.0D0*(i-1))/(1.0D0*coordsUnit(cKey)%xCopy)
                coords(cKeyE)%coords(m,2) = (coordsUnit(cKey)%unitCoords(n,2)+1.0D0*(j-1))/(1.0D0*coordsUnit(cKey)%yCopy)
                coords(cKeyE)%coords(m,3) = (coordsUnit(cKey)%unitCoords(n,3)+1.0D0*(k-1))/(1.0D0*coordsUnit(cKey)%zCopy)
! Forces
                coords(cKeyE)%forces(m,1) = coordsUnit(cKey)%unitForces(n,1)
                coords(cKeyE)%forces(m,2) = coordsUnit(cKey)%unitForces(n,2)
                coords(cKeyE)%forces(m,3) = coordsUnit(cKey)%unitForces(n,3)
              End Do
            End Do
          End Do
        End Do
        coords(cKeyE)%length = m
        coords(cKeyE)%aLat = coordsUnit(cKey)%xCopy*coordsUnit(cKey)%aLat
        coords(cKeyE)%unitCell = coordsUnit(cKey)%unitCell
      End If
    End Do
! Initialise data type
  End Subroutine expandUnitCoords



  Subroutine heatCoords(coords, maxVar, cKey_In)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Real(kind=DoubleReal) :: maxVar
    Integer(kind=StandardInteger), Optional :: cKey_In
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: i, j, k
! Optional Arguments
    cKey = 0
    If(Present(cKey_In))Then
      cKey = cKey_In
    End If
! Heat co-ordinates
    Do k=1,size(coords,1)
      If((cKey.eq.0).or.(k.eq.cKey))Then
        Do i=1,coords(cKey)%length
          Do j=1,3
            coords(cKey)%coords(i,j) = coords(cKey)%coords(i,j) + &
              2.0D0 * (RandomDist("G") - 0.5D0) *  maxVar
          End Do
        End Do
      End If
    End Do
  End Subroutine heatCoords


  Subroutine printCoords(coords, cKey)
! Print coords
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Integer(kind=StandardInteger) :: cKey
! Print summary
    If(cKey.eq.0)Then
      Do cKey=1,size(coords,1)
        Call printCoordsSR(coords,cKey)
      End Do
    Else
      Call printCoordsSR(coords,cKey)
    End If
  End Subroutine printCoords

  Subroutine printCoordsSR(coords, cKey)
! Print coords
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i
! Print coords
    Do i=1,coords(cKey)%length
      print "(I6,A4,A4,A1,I2,A2,E14.5,E14.5,E14.5)",i,"    ",coords(cKey)%label(i),"(",coords(cKey)%labelID(i),") ",&
        coords(cKey)%coords(i,1),coords(cKey)%coords(i,2),coords(cKey)%coords(i,3)
    End Do

  End Subroutine printCoordsSR



  Subroutine zeroForces(coords, cKey)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Dimension(:) :: coords
    Integer(kind=StandardInteger) :: cKey
! Zero
    If(cKey.eq.0)Then
      Do cKey=1,size(coords,1)
        coords(cKey)%forces = 0.0D0
      End Do
    Else
      coords(cKey)%forces = 0.0D0
    End If
  End Subroutine zeroForces



! ------------------------------------------------------------
!               Neighbour List
! ------------------------------------------------------------

  Include "geom.nl.f90"

! ------------------------------------------------------------
!               Neighbour List (for structure optimisation)
! ------------------------------------------------------------

  Include "geom.nl_opt.f90"

! ------------------------------------------------------------
!               Distortion Tensors
! ------------------------------------------------------------

  Include "geom.tensors.f90"



! -----------------------------------------------
!        Print Summaries
! -----------------------------------------------

  Subroutine printNLSummary(nl, cKeyIn)
! Print neighbour list summary
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Integer(kind=StandardInteger), Optional :: cKeyIn
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey
! Optional Arguments
    cKey = 0
    If(Present(cKeyIn))Then
      cKey = cKeyIn
    End If
! Start page
    Call addLinePage("Neighbour List Summary","T")
! Loop through configs
    If(cKey.eq.0)Then
      Do cKey=1,size(nl,1)
        If(nl(cKey)%length.gt.0)Then
          Call printNLSummary_Individual(nl, cKey)
        End If
      End Do
    Else
      If(nl(cKey)%length.gt.0)Then
        Call printNLSummary_Individual(nl, cKey)
      End If
    End If
    !Call printPage(newPage)
  End Subroutine printNLSummary

  Subroutine printNLSummary_Individual(nl, cKey)
! Print neighbour list summary
    Implicit None ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Character(Len=64) :: tempLine
! Output
    Write(tempLine,*) "Config ",cKey
    Call addLinePage(tempLine)
    Write(tempLine,*) "R Verlet: ",nl(cKey)%rVerlet
    Call addLinePage(tempLine)
    Write(tempLine,*) "R Min: ",nl(cKey)%rMin
    Call addLinePage(tempLine)
    Write(tempLine,*) "R Max: ",nl(cKey)%rMax
    Call addLinePage(tempLine)
    Write(tempLine,*) "NL Length: ",nl(cKey)%length
    Call addLinePage(tempLine)


  End Subroutine printNLSummary_Individual





! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------



! Basic energy-force functions

  Function ljEnergy(sigma, r) Result (vr)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal) :: sigma
    Real(kind=DoubleReal) :: r
! Out
    Real(kind=DoubleReal) :: vr
! Calc
    vr = 1.0D-5*4.0D0*((sigma/r)**12.0D0-(sigma/r)**6.0D0)
  End Function ljEnergy

  Function ljForce(sigma, r) Result (fr)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal) :: sigma
    Real(kind=DoubleReal) :: r
! Out
    Real(kind=DoubleReal) :: fr
! Calc
    fr = 1.0D-5*(48.0D0/r)*((sigma/r)**12.0D0-(sigma/r)**6.0D0)
  End Function ljForce













End Module geom
