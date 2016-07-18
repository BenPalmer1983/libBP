! ------------------------------------------------------------
!               GEOM: Neighbour List
! ------------------------------------------------------------

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
    Real(kind=DoubleReal), Dimension(1:3) :: position, positionT
    Integer(kind=StandardInteger), Dimension(1:3) :: shiftArr
    Integer(kind=StandardInteger) :: inCell
! Vars: Neighbour list estimate
    Integer(kind=StandardInteger) :: nlArrayLength
! Vars:  Allocatable arrays
    Integer(kind=StandardInteger), Allocatable, Dimension(:) :: scAtomCount     ! Number of atoms per sub cell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: scKeyArr        ! array of scKey and x,y,z position
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:,:) :: scCoordsI
    Real(kind=DoubleReal), Allocatable, Dimension(:,:,:) :: scCoordsR
!----------------------------------------------
! Optional arguments
!----------------------------------------------
    cKey_NL = 0
    If(Present(cKey_NL_in))Then
      cKey_NL = cKey_NL_in
    End If
!
!----------------------------------------------
! Array Management
!----------------------------------------------
!
! Deallocate if allocated
    Do cKey=1,size(nl,1)
      If((cKey.le.size(nl,1)).and.(coords(cKey)%length.gt.0))Then  ! Check that there is space in the nl object and there are coordinates
        If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!-------------------
    If(Allocated(nl(cKey)%atomA_ID))Then
      Deallocate(nl(cKey)%atomA_ID)
    End If
    If(Allocated(nl(cKey)%atomB_ID))Then
      Deallocate(nl(cKey)%atomB_ID)
    End If
    If(Allocated(nl(cKey)%atomA_Type))Then
      Deallocate(nl(cKey)%atomA_Type)
    End If
    If(Allocated(nl(cKey)%atomB_Type))Then
      Deallocate(nl(cKey)%atomB_Type)
    End If
    If(Allocated(nl(cKey)%atomPairKey))Then
      Deallocate(nl(cKey)%atomPairKey)
    End If
    If(Allocated(nl(cKey)%inCell))Then
      Deallocate(nl(cKey)%inCell)
    End If
    If(Allocated(nl(cKey)%rD))Then
      Deallocate(nl(cKey)%rD)
    End If
    If(Allocated(nl(cKey)%vecAB))Then
      Deallocate(nl(cKey)%vecAB)
    End If
    If(Allocated(nl(cKey)%subCell))Then
      Deallocate(nl(cKey)%subCell)
    End If
    If(Allocated(nl(cKey)%cell))Then
      Deallocate(nl(cKey)%cell)
    End If
!-------------------
        End If
      End If
    End Do
! Allocate if not allocated
    Do cKey=1,size(nl,1)
      If((cKey.le.size(nl,1)).and.(coords(cKey)%length.gt.0))Then  ! Check that there is space in the nl object and there are coordinates
        If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!-------------------
! Estimate NL size
      nlArrayLength = 2000+&
        ceiling(0.35D0*(coords(cKey)%length)**2*(8.0D0*rVerlet**3)/((coords(cKey)%aLat)**3))
      nl(cKey)%arrayLength = nlArrayLength
! Atom details
      Allocate(nl(cKey)%atomA_ID(1:nlArrayLength))
      Allocate(nl(cKey)%atomB_ID(1:nlArrayLength))
      Allocate(nl(cKey)%atomA_Type(1:nlArrayLength))
      Allocate(nl(cKey)%atomB_Type(1:nlArrayLength))
      Allocate(nl(cKey)%atomPairKey(1:nlArrayLength))
      Allocate(nl(cKey)%inCell(1:nlArrayLength))
! Position Details
      Allocate(nl(cKey)%rD(1:nlArrayLength))
      Allocate(nl(cKey)%vecAB(1:nlArrayLength,1:3))
! Subcell
      Allocate(nl(cKey)%subCell(1:nlArrayLength,1:3))
      Allocate(nl(cKey)%cell(1:nlArrayLength,1:3))
!-------------------
        End If
      End If
    End Do
!
!----------------------------------------------
! Build NL
!----------------------------------------------
!
! Loop through configs
!
    Do cKey=1,size(coords,1)
      If((cKey.le.size(nl,1)).and.(coords(cKey)%length.gt.0))Then  ! Check that there is space in the nl object and there are coordinates
        If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!----------------------------------------------
! Store from coords
    nl(cKey)%coordsLength = coords(cKey)%length
    nl(cKey)%aLat = coords(cKey)%aLat
    nl(cKey)%unitCell = coords(cKey)%unitCell
    nl(cKey)%label = coords(cKey)%label
    nl(cKey)%labelID = coords(cKey)%labelID
    nl(cKey)%coords = coords(cKey)%coords  ! Store initial coords
    nl(cKey)%forces = coords(cKey)%forces  ! Store initial forces
    nl(cKey)%coordsMD = coords(cKey)%coords  ! Store initial coords as first set of MD coords
    nl(cKey)%forcesMD = coords(cKey)%forces  ! Store initial forces as first set of MD forces
    nl(cKey)%atomIDs = coords(cKey)%atomIDs
    nl(cKey)%atomID_Count = coords(cKey)%atomID_Count
! Original NL parameters
    nl(cKey)%aLat_Original = nl(cKey)%aLat
!----------------------------------------------
! Init config specific variables
    coordLength = coords(cKey)%length
    rVerletSQ = rVerlet**2
    aLat = coords(cKey)%aLat
    nl(cKey)%totalRD = 0.0D0
    nl(cKey)%totalRDSq = 0.0D0
    nl(cKey)%rVerlet = rVerlet
    nl(cKey)%rMin = rVerlet
    nl(cKey)%rMax = 0.0D0
    nlKey = 0
!----------------------------------------------
! Keep coord within 1x1x1 box
    Do i = 1,nl(cKey)%coordsLength
      Do j=1,3
        nl(cKey)%coords(i,j) = Modulus(nl(cKey)%coords(i,j),1.0D0)
      End Do
    End Do
!----------------------------------------------
! calculate sub cell parameters
    scW = floor(aLat/(1.0D0*rVerlet))   ! Number of subcells (1D)
    If(scW.lt.1)Then
      scW = 1
    End If
    scCount = scW**3                                        ! Total number of subcells in 3D
    scAlat = aLat/(1.0D0*scW)                               ! Subcell width
    maxAtomsPerSC = 5*ceiling(coordLength/(1.0D0*scCount))  ! Estimate max atoms per subcell
    nl(cKey)%scCount = scCount
    nl(cKey)%scAlat = scAlat
    nl(cKey)%maxAtomsPerSC = maxAtomsPerSC
! Deallocate if allocated
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
! Init subcell arrays
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
      position(1) = aLat*nl(cKey)%coords(i,1)   ! x coord
      position(2) = aLat*nl(cKey)%coords(i,2)   ! y coord
      position(3) = aLat*nl(cKey)%coords(i,3)   ! z coord
      positionT = MatMul(nl(cKey)%unitCell,position)  ! To fix/implement
      scKey = SubCellKey(scAlat,scW,position(1),position(2),position(3))
      scAtomCount(scKey) = scAtomCount(scKey) + 1
! Store in sub cell arrays
      scCoordsI(scKey,scAtomCount(scKey),1) = i                          ! unique atom id
      scCoordsI(scKey,scAtomCount(scKey),2) = nl(cKey)%labelID(i)        ! atom type
     Do j=1,3
        scCoordsR(scKey,scAtomCount(scKey),j) = position(j) ! atom coords
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
                        ! super cell coords of atom B
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
!
!
  End Subroutine makeNL


  Subroutine updateNL(nl, cKey_NL_in, distortion_in)
! Update neighbour list without rebuilding
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Integer(kind=StandardInteger), Optional :: cKey_NL_in
    Real(kind=DoubleReal), Dimension(1:3,1:3), Optional :: distortion_in
! Vars:  Private
    Integer(kind=StandardInteger) :: n, cKey, cKey_NL
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: xD, yD, zD
    Real(kind=DoubleReal), Dimension(1:3) :: aVec, bVec
    Logical :: cellDistortion
!----------------------------------------------
! Optional arguments
!----------------------------------------------
    cKey_NL = 0
    If(Present(cKey_NL_in))Then
      cKey_NL = cKey_NL_in
    End If
    distortion = IdentityMatrix(distortion)
    cellDistortion = .false.
    If(Present(distortion_in))Then
      distortion = distortion_in
      cellDistortion = .true.
    End If
!----------------------------------------------
! Loop through neighbour lists
!----------------------------------------------
    Do cKey=1,size(nl,1)
      If((cKey_NL.eq.0).or.(cKey_NL.eq.cKey))Then
!------------------------
    aLat = nl(cKey)%aLat
! Loop through pairs
    Do n=1,nl(cKey)%length
      aVec(1) = nl(cKey)%coords(nl(cKey)%atomA_ID(n),1)
      aVec(2) = nl(cKey)%coords(nl(cKey)%atomA_ID(n),2)
      aVec(3) = nl(cKey)%coords(nl(cKey)%atomA_ID(n),3)
      bVec(1) = (nl(cKey)%coords(nl(cKey)%atomB_ID(n),1)+nl(cKey)%subCell(n,1))
      bVec(2) = (nl(cKey)%coords(nl(cKey)%atomB_ID(n),2)+nl(cKey)%subCell(n,2))
      bVec(3) = (nl(cKey)%coords(nl(cKey)%atomB_ID(n),3)+nl(cKey)%subCell(n,3))
! Distort position vectors
      If(cellDistortion)Then
        aVec = MatMul(distortion, aVec)
        bVec = MatMul(distortion, bVec)
      End If
! Update atom seperation
      xD = aLat*(aVec(1)-bVec(1)) ! xA - yA
      yD = aLat*(aVec(2)-bVec(2)) ! xA - yA
      zD = aLat*(aVec(3)-bVec(3)) ! xA - yA
      nl(cKey)%rD(n) = sqrt(xD**2+yD**2+zD**2)
! Update atom A to atom B vector
      nl(cKey)%vecAB(n,1) = xD/nl(cKey)%rD(n)
      nl(cKey)%vecAB(n,2) = yD/nl(cKey)%rD(n)
      nl(cKey)%vecAB(n,3) = zD/nl(cKey)%rD(n)
    End Do
!------------------------
      End If
    End Do
  End Subroutine updateNL


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
