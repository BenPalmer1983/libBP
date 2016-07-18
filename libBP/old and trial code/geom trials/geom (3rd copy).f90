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
  Integer(kind=StandardInteger), Parameter :: p_nlMax = 12300      ! Neighbour list length
! Make private
  Private
! Public Variables and Parameters
  Public :: p_confs, p_cuMax, p_cMax, p_nlMax
! Public derived types
  Public :: coordsUnitType, coordsType, nlType

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
    Integer(kind=StandardInteger) :: points = 0
    Real(kind=DoubleReal) :: aLat = 0.0D0
    Character(Len=16), Dimension(1:p_cMax) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:p_cMax) :: labelID
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: fracCoords   ! coords - fractional
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: coords   ! coords - fractional
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: forces   ! coords - fractional
    Integer(kind=StandardInteger) :: IDcount
    Character(Len=16), Dimension(1:64) :: atomIDs
    Integer(kind=StandardInteger) :: atomID_Count
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: electronDensity
    Real(kind=DoubleReal), Dimension(1:p_cMax,1:3) :: atomEnergy
    Real(kind=DoubleReal) :: pairEnergy
    Real(kind=DoubleReal) :: embeddingEnergy
    Real(kind=DoubleReal) :: totalEnergy
  End Type coordsType

  Type :: nlType
! 6.2MB per neighbour list
! Atom details
    !Integer(kind=StandardInteger), Dimension(1:p_nlMax) :: atomA_ID
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: atomA_ID
    Integer(kind=StandardInteger), Dimension(1:p_nlMax) :: atomB_ID
    Integer(kind=StandardInteger), Dimension(1:p_nlMax) :: atomA_Type
    Integer(kind=StandardInteger), Dimension(1:p_nlMax) :: atomB_Type
    Integer(kind=StandardInteger), Dimension(1:p_nlMax) :: atomPairKey
    Integer(kind=StandardInteger), Dimension(1:p_nlMax) :: inCell

! Position Details
    Real(kind=DoubleReal), Dimension(1:p_nlMax) :: rD
    Real(kind=DoubleReal), Dimension(1:p_nlMax,1:3) :: vecAB
    Real(kind=DoubleReal), Dimension(1:p_nlMax,1:3) :: coordsA
    Real(kind=DoubleReal), Dimension(1:p_nlMax,1:3) :: coordsB
! Subcell
    Integer(kind=StandardInteger), Dimension(1:p_nlMax,1:3) :: subCell
    Integer(kind=StandardInteger), Dimension(1:p_nlMax,1:3) :: cell
! Subcell
    Integer(kind=StandardInteger) :: scCount = 0
    Real(kind=DoubleReal) :: scALat
    Integer(kind=StandardInteger) :: maxAtomsPerSC = 0
! Misc
    Integer(kind=StandardInteger) :: length = 0
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: rMin
    Real(kind=DoubleReal) :: rMax
    Real(kind=DoubleReal) :: rVerlet
    Real(kind=DoubleReal) :: totalRD
    Real(kind=DoubleReal) :: totalRDSq
  End Type



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
  Public :: zeroForces
  Public :: initNL
  Public :: makeNL
  Public :: printCoords
  Public :: printNLSummary



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
          coords(cKey)%fracCoords(i,j) = 0.0D0
          coords(cKey)%coords(i,j) = 0.0D0
          coords(cKey)%forces(i,j) = 0.0D0
        End Do
      End Do
      coords(cKey)%points = 0
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
! Vars:  In
    Type(coordsUnitType), Dimension(:) :: coordsUnit
    Type(coordsType), Dimension(:) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, k, m, n, cKey
    Real(kind=DoubleReal), Dimension(1:3) :: xCoords
! Loop through
    Do cKey=1,size(coordsUnit,1)
      If((cKey.le.size(coords,1)).and.(coordsUnit(cKey)%points.gt.0))Then
        m = 0
        Do i=1,coordsUnit(cKey)%xCopy
          Do j=1,coordsUnit(cKey)%yCopy
            Do k=1,coordsUnit(cKey)%zCopy
              Do n=1,coordsUnit(cKey)%points
                m = m + 1
                coords(cKey)%label(m) = coordsUnit(cKey)%label(n)
                coords(cKey)%labelID(m) = coordsUnit(cKey)%labelID(n)
! Fractional coords
                coords(cKey)%fracCoords(m,1) = (coordsUnit(cKey)%unitCoords(n,1)+1.0D0*(i-1))/(1.0D0*coordsUnit(cKey)%xCopy)
                coords(cKey)%fracCoords(m,2) = (coordsUnit(cKey)%unitCoords(n,2)+1.0D0*(j-1))/(1.0D0*coordsUnit(cKey)%yCopy)
                coords(cKey)%fracCoords(m,3) = (coordsUnit(cKey)%unitCoords(n,3)+1.0D0*(k-1))/(1.0D0*coordsUnit(cKey)%zCopy)
! Real coords
                xCoords(1) = coordsUnit(cKey)%aLat*(coordsUnit(cKey)%unitCoords(n,1)+1.0D0*(i-1))
                xCoords(2) = coordsUnit(cKey)%aLat*(coordsUnit(cKey)%unitCoords(n,2)+1.0D0*(j-1))
                xCoords(3) = coordsUnit(cKey)%aLat*(coordsUnit(cKey)%unitCoords(n,3)+1.0D0*(k-1))
                xCoords = TransformCoords(xCoords, coordsUnit(cKey)%unitCell)
                coords(cKey)%coords(m,1) = xCoords(1)
                coords(cKey)%coords(m,2) = xCoords(2)
                coords(cKey)%coords(m,3) = xCoords(3)
! Forces
                coords(cKey)%fracCoords(m,1) = coordsUnit(cKey)%unitForces(n,1)
                coords(cKey)%fracCoords(m,2) = coordsUnit(cKey)%unitForces(n,2)
                coords(cKey)%fracCoords(m,3) = coordsUnit(cKey)%unitForces(n,3)
              End Do
            End Do
          End Do
        End Do
        coords(cKey)%points = m
        coords(cKey)%aLat = coordsUnit(cKey)%xCopy*coordsUnit(cKey)%aLat
      End If
    End Do
! Initialise data type
  End Subroutine expandUnitCoords



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
    Do i=1,coords(cKey)%points
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

  Subroutine initNL(nl)
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey

    Do cKey = 1,size(nl,1)
      Allocate(nl(cKey)%atomA_ID(1:12000))
    End Do

  End Subroutine initNL


  Subroutine makeNL(nl, coords, rVerlet, cKey_NL_in)
! Make neighbour list for atoms
! The input coords and alat must be large enough so the same atom does not interact
! with a copy in a periodic cell surrounding the original
! e.g. if the rVerlet cutoff is 5, the alat must be greater than 5
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(coordsType), Dimension(:)  :: coords
    Real(kind=DoubleReal) :: rVerlet
    Integer(kind=StandardInteger), Optional :: cKey_NL_in
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey_NL
    Integer(kind=StandardInteger) :: cKey
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: l, m ,n
    Integer(kind=StandardInteger) :: scX, scY, scZ
    Integer(kind=StandardInteger) :: atomA, atomB, atomA_ID, atomB_ID
    Integer(kind=StandardInteger) :: coordLength, scW, scCount, scKey, maxAtomsPerSC
    Integer(kind=StandardInteger) :: scKeyA, scKeyB
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: rVerletSQ, scAlat
    Integer(kind=StandardInteger) :: nlKey
    Integer(kind=StandardInteger) :: xKey, yKey, zKey
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB
    Real(kind=DoubleReal) :: xD, yD, zD, rD
    Real(kind=DoubleReal) :: xDsq, yDsq, zDsq, rDsq
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Integer(kind=StandardInteger), Dimension(1:3) :: shiftArr
    Integer(kind=StandardInteger) :: inCell
! Vars:  Allocatable arrays
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: scAtomCount     ! Number of atoms per sub cell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: scKeyArr        ! array of scKey and x,y,z position
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: scCoordsI
    Real(kind=DoubleReal), Allocatable, Dimension(:,:,:) :: scCoordsR
! Allocate if not allocated
    Do cKey=1,size(nl,1)
      Allocate(nl(cKey)%atomA_ID(1:16000))
    End Do
! Optional arguments
    cKey_NL = 0
    If(Present(cKey_NL_in))Then
      cKey_NL = cKey_NL_in
    End If
! Loop through configs
    Do cKey=1,size(coords,1)
      If((cKey.le.size(nl,1)).and.(coords(cKey)%points.gt.0))Then  ! Check that there is space in the nl object and there are coordinates
        If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
    print *,cKey,coords(cKey)%label(1),coords(cKey)%labelID(1)
    print *,cKey,coords(cKey)%label(2),coords(cKey)%labelID(2)
    print *,cKey,coords(cKey)%label(3),coords(cKey)%labelID(3)
!----------------------------------------------
! Init config specific variables
    coordLength = coords(cKey)%points
    rVerletSQ = rVerlet**2
    aLat = coords(cKey)%aLat
    nl(cKey)%aLat = coords(cKey)%aLat
    nl(cKey)%totalRD = 0.0D0
    nl(cKey)%totalRDSq = 0.0D0
    nl(cKey)%rVerlet = rVerlet
    nl(cKey)%rMin = rVerlet
    nl(cKey)%rMax = 0.0D0
    nlKey = 0
! calculate sub cell parameters
    scW = floor(aLat/(1.0D0*rVerlet))
    If(scW.lt.1)Then
      scW = 1
    End If
    scCount = scW**3
    scAlat = aLat/(1.0D0*scW)
    maxAtomsPerSC = 5*ceiling(coordLength/(1.0D0*scCount))
    nl(cKey)%scCount = scCount
    nl(cKey)%scAlat = scAlat
    nl(cKey)%maxAtomsPerSC = maxAtomsPerSC
    If(Allocated(scAtomCount))Then
      Deallocate(scAtomCount)
    End If
    If(Allocated(scKeyArr))Then
      Deallocate(scKeyArr)
    End If
    If(Allocated(scCoordsI))Then
      Deallocate(scCoordsI)
    End If
    If(Allocated(scCoordsR))Then
      Deallocate(scCoordsR)
    End If
! Allocate arrays
    Allocate(scAtomCount(1:scCount))
    Allocate(scKeyArr(1:scCount,1:3))
    Allocate(scCoordsI(1:scCount,1:maxAtomsPerSC,1:2))
    Allocate(scCoordsR(1:scCount,1:maxAtomsPerSC,1:3))
! Init arrays
    scAtomCount = 0
    Do k=1,scW
      Do j=1,scW
        Do i=1,scW
          scKey = i+scW*(j-1)+scW**2*(k-1)
          scKeyArr(scKey,1) = i
          scKeyArr(scKey,2) = j
          scKeyArr(scKey,3) = k
        End Do
      End Do
    End Do
! Loop through atoms
    Do i=1,coordLength
      !scKey = SubCellKey(scAlat,scW,atomCoords(i,1),atomCoords(i,2),atomCoords(i,3))
      scKey = SubCellKey(scAlat,scW,coords(cKey)%coords(i,1),coords(cKey)%coords(i,2),coords(cKey)%coords(i,3))
      scAtomCount(scKey) = scAtomCount(scKey) + 1
! Store in sub cell arrays
      scCoordsI(scKey,scAtomCount(scKey),1) = i                 ! unique atom id
      scCoordsI(scKey,scAtomCount(scKey),2) = coords(cKey)%labelID(i)      ! atom type
      Do j=1,3
        !scCoordsR(scKey,scAtomCount(scKey),j) = atomCoords(i,j)
        scCoordsR(scKey,scAtomCount(scKey),j) = coords(cKey)%coords(i,j)
      End Do
    End Do
! Make neighbour list
! Loop through subcells
    Do scKeyA=1,scCount
! loop through surrounding and image sub cells
      Do l=-1,1
        Do m=-1,1
          Do n=-1,1
! Atom B subcell key
            xKey = scKeyArr(scKeyA,1)+l
            yKey = scKeyArr(scKeyA,2)+m
            zKey = scKeyArr(scKeyA,3)+n
! PBC subcell key
            scX = Modulus(xKey-1,scW)+1
            scY = Modulus(yKey-1,scW)+1
            scZ = Modulus(zKey-1,scW)+1
            scKeyB = scX+scW*(scY-1)+scW**2*(scZ-1)
! set default coord shift
            xShift = 0.0D0
            yShift = 0.0D0
            zShift = 0.0D0
            shiftArr = 0
! adjust shift
            inCell = 1
            If(xKey.lt.1)Then
              xShift = -1.0D0 * aLat
              shiftArr(1) = -1
              inCell = 0
            End If
            If(yKey.lt.1)Then
              yShift = -1.0D0 * aLat
              shiftArr(2) = -1
              inCell = 0
            End If
            If(zKey.lt.1)Then
              zShift = -1.0D0 * aLat
              shiftArr(3) = -1
              inCell = 0
            End If
            If(xKey.gt.scW)Then
              xShift = 1.0D0 * aLat
              shiftArr(1) = 1
              inCell = 0
            End If
            If(yKey.gt.scW)Then
              yShift = 1.0D0 * aLat
              shiftArr(2) = 1
              inCell = 0
            End If
            If(zKey.gt.scW)Then
              zShift = 1.0D0 * aLat
              shiftArr(3) = 1
              inCell = 0
            End If
! loop through A atoms
            Do atomA =1,scAtomCount(scKeyA)
! loop through B atoms
              Do atomB =atomA+1,scAtomCount(scKeyB)
                atomA_ID = scCoordsI(scKeyA,atomA,1)
                atomB_ID = scCoordsI(scKeyB,atomB,1)
                xA = scCoordsR(scKeyA,atomA,1)
                xB = scCoordsR(scKeyB,atomB,1)+xShift
                xD = xA-xB
                xDsq = xd**2
                If(xDsq.le.rVerletSQ)Then
                  yA = scCoordsR(scKeyA,atomA,2)
                  yB = scCoordsR(scKeyB,atomB,2)+yShift
                  yD = yA-yB
                  yDsq = yd**2
                  If(yDsq.le.rVerletSQ)Then
                    zA = scCoordsR(scKeyA,atomA,3)
                    zB = scCoordsR(scKeyB,atomB,3)+zShift
                    zD = zA-zB
                    zDsq = zd**2
                    If(zDsq.le.rVerletSQ)Then
                      rdSq = xDsq + yDsq + zDsq
                      If(rdSq.le.rVerletSq)Then
                        rd = sqrt(rdSq)
                        ! Key
                        nlKey = nlKey + 1
                        ! Key/Type
                        nl(cKey)%atomA_ID(nlKey) = atomA_ID                   ! A ID
                        nl(cKey)%atomB_ID(nlKey) = atomB_ID                   ! B ID
                        nl(cKey)%atomA_Type(nlKey) = scCoordsI(scKeyA,atomA,2)  ! A type
                        nl(cKey)%atomB_Type(nlKey) = scCoordsI(scKeyB,atomB,2)  ! B type
                        nl(cKey)%inCell(nlKey) = inCell                     ! 1 if in cell, 0 if out of cell
                        nl(cKey)%atomPairKey = DoubleKey(scCoordsI(scKeyA,atomA,2),scCoordsI(scKeyB,atomB,2))
                        ! Displacement/Direction
                        nl(cKey)%rD(nlKey) = rD
                        nl(cKey)%vecAB(nlKey,1) = xD/rD                ! Vector from B to A (x)
                        nl(cKey)%vecAB(nlKey,2) = yD/rD                ! Vector from B to A (y)
                        nl(cKey)%vecAB(nlKey,3) = zD/rD                ! Vector from B to A (z)
                        nl(cKey)%coordsA(nlKey,1) = xA
                        nl(cKey)%coordsA(nlKey,2) = yA
                        nl(cKey)%coordsA(nlKey,3) = zA
                        nl(cKey)%coordsB(nlKey,1) = xB
                        nl(cKey)%coordsB(nlKey,2) = yB
                        nl(cKey)%coordsB(nlKey,3) = zB
                        ! Cell
                        nl(cKey)%subCell(nlKey,1) = shiftArr(1)
                        nl(cKey)%subCell(nlKey,2) = shiftArr(2)
                        nl(cKey)%subCell(nlKey,3) = shiftArr(3)
                        If(rD.gt.nl(cKey)%rMax)Then
                          nl(cKey)%rMax = rD
                        End If
                        If(rD.lt.nl(cKey)%rMin)Then
                          nl(cKey)%rMin = rD
                        End If
                      End If
                    End If
                  End If
                End If
              End Do
            End Do
          End Do
        End Do
      End Do
    End Do
    nl(cKey)%length = nlKey
!----------------------------------------------
        End If
      End If
    End Do

  End Subroutine makeNL



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


  Function SubCellKey (scAlat, scW, x, y, z) Result (key)
! Return sub cell key
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal) :: scAlat, x, y, z
    Integer(kind=StandardInteger) :: scW
! Out
    Integer(kind=StandardInteger) :: key
! Private
    Integer(kind=StandardInteger) :: i, j, k
    i = floor(x/(1.0D0*scAlat))+1
    j = floor(y/(1.0D0*scAlat))+1
    k = floor(z/(1.0D0*scAlat))+1
    key = i+scW*(j-1)+scW**2*(k-1)
  End Function SubCellKey


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
