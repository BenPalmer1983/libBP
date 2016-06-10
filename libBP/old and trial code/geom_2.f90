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

  Type :: coordsUnitType
    Integer(kind=StandardInteger) :: xCopy = 1
    Integer(kind=StandardInteger) :: yCopy = 1
    Integer(kind=StandardInteger) :: zCopy = 1
    Real(kind=DoubleReal) :: aLat = 1.0D0
    Character(Len=16), Dimension(1:128) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:128) :: labelID
    Real(kind=DoubleReal), Dimension(1:128,1:3) :: unitCoords   ! Unit coords - fractional
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: unitCell
    Integer(kind=StandardInteger) :: points = 0
  End Type coordsUnitType

  Type :: coordsType
    Real(kind=DoubleReal) :: aLat = 0.0D0
    Character(Len=16), Dimension(1:1024) :: label           ! Unit labels
    Integer(kind=StandardInteger), Dimension(1:1024) :: labelID
    Real(kind=DoubleReal), Dimension(1:1024,1:3) :: fracCoords   ! coords - fractional
    Real(kind=DoubleReal), Dimension(1:1024,1:3) :: coords   ! coords - fractional
    Integer(kind=StandardInteger) :: points = 0
  End Type coordsType

  Type :: nlType
! 6.2MB per neighbour list
    Integer(kind=StandardInteger), Dimension(1:50000,1:5) :: i        ! A ID, B ID, A Type, B Type
    Real(kind=DoubleReal), Dimension(1:50000,1:10) :: r                ! rD, xD/rD, yD/rD, zD/rD
    Integer(kind=StandardInteger), Dimension(1:50000,1:3) :: subCell
    Integer(kind=StandardInteger), Dimension(1:50000,1:3) :: cell
    Integer(kind=StandardInteger) :: length = 0
    Integer(kind=StandardInteger) :: scCount = 0
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
  Use matrix
  Use basicMaths
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
  Public :: makeNL



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
    Type(coordsUnitType) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
! Initialise data type
    coords%xCopy = 1
    coords%yCopy = 1
    coords%zCopy = 1
    coords%aLat = 1.0D0
    Do i=1,size(coords%label,1)
      coords%label(i) = BlankString(coords%label(i))
      coords%labelID(i) = 0
      Do j=1,3
        coords%unitCoords(i,j) = 0.0D0
      End Do
    End Do
    Do i=1,3
      Do j=1,3
        If(i.eq.j)Then
          coords%unitCell(i,j) = 1.0D0
        Else
          coords%unitCell(i,j) = 0.0D0
        End If
      End Do
    End Do
    coords%points = 0
  End Subroutine initUnitCoords

  Subroutine initCoords(coords)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(coordsType) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
! Initialise data type
    Do i=1,size(coords%label,1)
      coords%label(i) = BlankString(coords%label(i))
      coords%labelID(i) = 0
      Do j=1,3
        coords%fracCoords(i,j) = 0.0D0
        coords%coords(i,j) = 0.0D0
      End Do
    End Do
    coords%points = 0
  End Subroutine initCoords

  Subroutine standardCoords(typeCell, coords)
! Loads standard unit coords
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Character(*) :: typeCell
    Type(coordsUnitType) :: coords
! Vars:  Private
    Character(Len(typeCell)) :: typeCellUpper
! Read in standard type
    typeCellUpper = StrToUpper(typeCell)
    If(typeCellUpper(1:3).eq."FCC")Then
      coords%points = 4
      coords%label(1) = "A"
      coords%labelID(1) = 1
      coords%unitCoords(1,1) = 0.0D0
      coords%unitCoords(1,2) = 0.0D0
      coords%unitCoords(1,3) = 0.0D0
      coords%label(2) = "B"
      coords%labelID(2) = 2
      coords%unitCoords(2,1) = 0.5D0
      coords%unitCoords(2,2) = 0.5D0
      coords%unitCoords(2,3) = 0.0D0
      coords%label(3) = "C"
      coords%labelID(3) = 3
      coords%unitCoords(3,1) = 0.5D0
      coords%unitCoords(3,2) = 0.0D0
      coords%unitCoords(3,3) = 0.5D0
      coords%label(4) = "D"
      coords%unitCoords(4,1) = 0.0D0
      coords%unitCoords(4,2) = 0.5D0
      coords%unitCoords(4,3) = 0.5D0
    End If
  End Subroutine standardCoords

  Subroutine expandUnitCoords(coordsUnit, coords)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(coordsUnitType) :: coordsUnit
    Type(coordsType) :: coords
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, k, m, n
    Real(kind=DoubleReal), Dimension(1:3) :: xCoords
! Loop through
    m = 0
    Do i=1,coordsUnit%xCopy
      Do j=1,coordsUnit%yCopy
        Do k=1,coordsUnit%zCopy
          Do n=1,coordsUnit%points
            m = m + 1
            coords%label(m) = coordsUnit%label(n)
            coords%labelID(m) = coordsUnit%labelID(n)
! Fractional coords
            coords%fracCoords(m,1) = (coordsUnit%unitCoords(n,1)+1.0D0*(i-1))/(1.0D0*coordsUnit%xCopy)
            coords%fracCoords(m,2) = (coordsUnit%unitCoords(n,2)+1.0D0*(j-1))/(1.0D0*coordsUnit%yCopy)
            coords%fracCoords(m,3) = (coordsUnit%unitCoords(n,3)+1.0D0*(k-1))/(1.0D0*coordsUnit%zCopy)
! Real coords
            xCoords(1) = coordsUnit%aLat*(coordsUnit%unitCoords(n,1)+1.0D0*(i-1))
            xCoords(2) = coordsUnit%aLat*(coordsUnit%unitCoords(n,2)+1.0D0*(j-1))
            xCoords(3) = coordsUnit%aLat*(coordsUnit%unitCoords(n,3)+1.0D0*(k-1))
            xCoords = TransformCoords(xCoords, coordsUnit%unitCell)
            coords%coords(m,1) = xCoords(1)
            coords%coords(m,2) = xCoords(2)
            coords%coords(m,3) = xCoords(3)
          End Do
        End Do
      End Do
    End Do
    coords%points = m
    coords%aLat = coordsUnit%xCopy*coordsUnit%aLat
! Initialise data type
  End Subroutine expandUnitCoords

! ------------------------------------------------------------
!               Neighbour List
! ------------------------------------------------------------

  Subroutine makeNL(nl, coords, rVerlet)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    !Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    !Real(kind=DoubleReal), Dimension(:,:) :: atomCoords
    Type(coordsType) :: coords
    Real(kind=DoubleReal) :: rVerlet, aLat
! Private Variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: l, m ,n
    Integer(kind=StandardInteger) :: scX, scY, scZ
    Integer(kind=StandardInteger) :: atomA, atomB, atomA_ID, atomB_ID
    Integer(kind=StandardInteger) :: coordLength, scW, scCount, scKey, maxAtomsPerSC
    Integer(kind=StandardInteger) :: scKeyA, scKeyB
    Real(kind=DoubleReal) :: rVerletSQ, scAlat
    Integer(kind=StandardInteger), &
    Dimension(1:100000) :: nlUniqueKeyArr
    Integer(kind=StandardInteger) :: nlKey, uKey
    Integer(kind=StandardInteger) :: xKey, yKey, zKey
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB
    Real(kind=DoubleReal) :: xD, yD, zD, rD
    Real(kind=DoubleReal) :: xDsq, yDsq, zDsq, rDsq
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Integer(kind=StandardInteger), Dimension(1:3) :: shiftArr
    Integer(kind=StandardInteger) :: inCell
! Allocatable arrays
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: scAtomCount     ! Number of atoms per sub cell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: scKeyArr        ! array of scKey and x,y,z position
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: scCoordsI
    Real(kind=DoubleReal), Allocatable, Dimension(:,:,:) :: scCoordsR
! Init config specific variables
    coordLength = coords%points
    rVerletSQ = rVerlet**2
    aLat = coords%aLat
    nl%aLat = coords%aLat
    nl%totalRD = 0.0D0
    nl%totalRDSq = 0.0D0
    nl%rVerlet = rVerlet
    nlUniqueKeyArr = 0
    nlKey = 0
! calculate sub cell parameters
    scW = floor(aLat/(1.0D0*rVerlet))
    If(scW.lt.1)Then
      scW = 1
    End If
    scCount = scW**3
    scAlat = aLat/(1.0D0*scW)
    maxAtomsPerSC = 5*ceiling(coordLength/(1.0D0*scCount))
    nl%scCount = scCount
    print *,aLat,rVerlet,scW,scCount,scAlat,maxAtomsPerSC
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
      scKey = SubCellKey(scAlat,scW,coords%coords(i,1),coords%coords(i,2),coords%coords(i,3))
      scAtomCount(scKey) = scAtomCount(scKey) + 1
! Store in sub cell arrays
      scCoordsI(scKey,scAtomCount(scKey),1) = i                 ! unique atom id
      scCoordsI(scKey,scAtomCount(scKey),2) = coords%labelID(i)      ! atom type
      Do j=1,3
        !scCoordsR(scKey,scAtomCount(scKey),j) = atomCoords(i,j)
        scCoordsR(scKey,scAtomCount(scKey),j) = coords%coords(i,j)
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
              Do atomB =1,scAtomCount(scKeyB)
                If(l.eq.0.and.m.eq.0.and.n.eq.0.and.atomA.eq.atomB)Then
                  ! skip
                Else
                  atomA_ID = scCoordsI(scKeyA,atomA,1)
                  atomB_ID = scCoordsI(scKeyB,atomB,1)
                  If(atomA_ID.gt.atomB_ID)Then
                    uKey = (atomA_ID-1)*(atomA_ID)/2+atomB_ID
                  Else
                    uKey = (atomB_ID-1)*(atomB_ID)/2+atomA_ID
                  End If
                  If(nlUniqueKeyArr(uKey).eq.0)Then
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
                            nlUniqueKeyArr(uKey) = 1
                            rd = sqrt(rdSq)
                            ! Key
                            nlKey = nlKey + 1
                            ! Key/Type
                            nl%i(nlKey,1) = atomA_ID                   ! A ID
                            nl%i(nlKey,2) = atomB_ID                   ! B ID
                            nl%i(nlKey,3) = scCoordsI(scKeyA,atomA,2)  ! A type
                            nl%i(nlKey,4) = scCoordsI(scKeyB,atomB,2)  ! B type
                            nl%i(nlKey,5) = inCell                     ! 1 if in cell, 0 if out of cell
                            ! Displacement/Direction
                            nl%r(nlKey,1) = rD
                            nl%r(nlKey,2) = xD/rD                ! Vector from B to A (x)
                            nl%r(nlKey,3) = yD/rD                ! Vector from B to A (y)
                            nl%r(nlKey,4) = zD/rD                ! Vector from B to A (z)
                            nl%r(nlKey,5) = xA
                            nl%r(nlKey,6) = yA
                            nl%r(nlKey,7) = zA
                            nl%r(nlKey,8) = xB
                            nl%r(nlKey,9) = yB
                            nl%r(nlKey,10) = zB
                            ! Cell
                            nl%subCell(nlKey,1) = shiftArr(1)
                            nl%subCell(nlKey,2) = shiftArr(2)
                            nl%subCell(nlKey,3) = shiftArr(3)
                            If(nlKey.eq.1)Then
                              nl%rMin = rD
                              nl%rMax = rD
                            Else
                              If(rD.gt.nl%rMax)Then
                                nl%rMax = rD
                              End If
                              If(rD.lt.nl%rMin)Then
                                nl%rMin = rD
                              End If
                            End If
                          End If
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
    nl%length = nlKey
    print *,"NL Length: ",nl%length
  End Subroutine makeNL










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
