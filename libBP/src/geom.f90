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
  
  Type :: nlType    
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: i        ! A ID, B ID, A Type, B Type 
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: r                ! rD, xD/rD, yD/rD, zD/rD
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: subCell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: cell    
    Integer(kind=StandardInteger) :: length = 0
    Integer(kind=StandardInteger) :: scCount = 0
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: rMin
    Real(kind=DoubleReal) :: rMax
    Real(kind=DoubleReal) :: rVerlet 
    Real(kind=DoubleReal) :: totalRD
    Real(kind=DoubleReal) :: totalRDSq
  End Type  
  
  Type :: mdType  
    Real(kind=DoubleReal) :: atomEnergy
    Real(kind=DoubleReal),Dimension(1:3) :: coords
    Real(kind=DoubleReal),Dimension(1:3) :: force
    Real(kind=DoubleReal),Dimension(1:3) :: velocity
  End Type    

End Module geomTypes


Module geom
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use basicMaths
  Use rng
  Use linearAlgebra
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
  Public :: makeCoords
  Public :: initNL
  Public :: makeNL
  Public :: updateNL
  Public :: makeNLOriginal
  !Public :: mdRun
  !Public :: configOpt
  !Public :: bfgsRelax
  !Public :: simpleRelax
  
! Interfaces  
  Interface makeNL
    Module Procedure makeNL_A, makeNL_B
  End Interface makeNL
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  
! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------  
  
  Subroutine makeCoords(coords, typeCell, copiesX, copiesY, copiesZ, aLat)
    Implicit None   ! Force declaration of all variables 
    
    Real(kind=DoubleReal), Dimension(:,:) :: coords
    Integer(kind=StandardInteger) :: copiesX, copiesY, copiesZ, typeCell
    Real(kind=DoubleReal) :: aLat
    
    Real(kind=DoubleReal), Dimension(1:2,1:3) :: bccUnit 
    Real(kind=DoubleReal), Dimension(1:4,1:3) :: fccUnit         
! ----------
! Init
! ----------
    coords = 0    
! -------------------
! Structure atom arrangements       
! -------------------       
! FCC atom 1
    fccUnit(1,1) = 0.0D0
    fccUnit(1,2) = 0.0D0
    fccUnit(1,3) = 0.0D0    
! FCC atom 1
    fccUnit(2,1) = 0.5D0
    fccUnit(2,2) = 0.5D0
    fccUnit(2,3) = 0.0D0    
! FCC atom 1
    fccUnit(3,1) = 0.5D0
    fccUnit(3,2) = 0.0D0
    fccUnit(3,3) = 0.5D0    
! FCC atom 1
    fccUnit(4,1) = 0.0D0
    fccUnit(4,2) = 0.5D0
    fccUnit(4,3) = 0.5D0   
! BCC atom 1
    bccUnit(1,1) = 0.0D0
    bccUnit(1,2) = 0.0D0
    bccUnit(1,3) = 0.0D0    
! BCC atom 1
    bccUnit(2,1) = 0.5D0
    bccUnit(2,2) = 0.5D0
    bccUnit(2,3) = 0.5D0     

    If(typeCell.eq.1)Then
      Call makeCoordsProcess(coords, fccUnit, copiesX, copiesY, copiesZ, aLat)    
    End If
    If(typeCell.eq.2)Then
      Call makeCoordsProcess(coords, bccUnit, copiesX, copiesY, copiesZ, aLat)    
    End If
  End Subroutine makeCoords
  
  Subroutine makeCoordsProcess(coords, unitCellIn, copiesX, copiesY, copiesZ, aLat)
    Implicit None   ! Force declaration of all variables
! In/Out
    Real(kind=DoubleReal), Dimension(:,:) :: coords
    Real(kind=DoubleReal), Dimension(:,:) :: unitCellIn
    Integer(kind=StandardInteger) :: copiesX, copiesY, copiesZ
    Real(kind=DoubleReal) :: aLat
! Private variables  
    Real(kind=DoubleReal), Dimension(1:size(unitCellIn,1),1:size(unitCellIn,2)) :: unitCell
    Integer(kind=StandardInteger) :: coordID, i
    Integer(kind=StandardInteger) :: xLoop, yLoop, zLoop
! Init
    coordID = 1
    unitCell = unitCellIn
    Do xLoop=1,copiesX
      Do yLoop=1,copiesY
        Do zLoop=1,copiesZ 
          Do i=1,size(unitCell,1)   
! Fractional coords
            coords(coordID,1) = (xLoop + unitCell(i,1) - 1.0D0)*aLat
            coords(coordID,2) = (yLoop + unitCell(i,2) - 1.0D0)*aLat
            coords(coordID,3) = (zLoop + unitCell(i,3) - 1.0D0)*aLat
! Increment
            coordID = coordID + 1 
          End Do  
        End Do
      End Do   
    End Do
  End Subroutine makeCoordsProcess
  
  
  
  Subroutine initNL(nl, allocatedLengthIn)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    Integer(kind=StandardInteger), Optional :: allocatedLengthIn
! Private    
    Integer(kind=StandardInteger) :: allocatedLength
! Optional argument
    allocatedLength = 25000    
    If(Present(allocatedLengthIn))Then
      allocatedLength = allocatedLengthIn
    End If    
! Deallocate arrays
    If(Allocated(nl%i))Then
      Deallocate(nl%i)
    End If
    If(Allocated(nl%r))Then
      Deallocate(nl%r)
    End If
    If(Allocated(nl%cell))Then
      Deallocate(nl%cell)
    End If
! Allocate arrays
    Allocate(nl%i(1:allocatedLength,1:5))
    Allocate(nl%r(1:allocatedLength,1:10))
    Allocate(nl%subCell(1:allocatedLength,1:3))
    Allocate(nl%cell(1:allocatedLength,1:3))    
  End Subroutine initNL
  
  
  Subroutine makeNL_A(nl, atomTypesIn, atomCoordsIn, coordStart, coordEnd, rVerlet, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    Integer(kind=StandardInteger), Dimension(:) :: atomTypesIn
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoordsIn  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rVerlet, aLat   
! Private Variables
    Integer(kind=StandardInteger) :: i, j, coordLength
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomCoords 
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1)) :: atomTypes      
! Init config specific variables
    coordLength = coordEnd-coordStart+1    
! Sort out coords so they fall within alatxalatxalat box
    Do i=1,coordLength
      atomTypes(i) = atomTypesIn(coordStart+i-1)
      Do j=1,3
        atomCoords(i,j) = Modulus(atomCoordsIn(coordStart+i-1,j),aLat)
      End Do
    End Do    
! Call process subroutine
    Call makeNLProcess(nl, atomTypes, atomCoords, coordStart, coordEnd, rVerlet, aLat)    
  End Subroutine makeNL_A 
  
  Subroutine makeNL_B(nl, atomTypesIn, atomCoordsIn, coordStart, coordEnd, rVerlet, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    Integer(kind=StandardInteger), Dimension(:,:) :: atomTypesIn
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoordsIn  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rVerlet, aLat   
! Private Variables
    Integer(kind=StandardInteger) :: i, j, coordLength
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomCoords 
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1)) :: atomTypes      
! Init config specific variables
    coordLength = coordEnd-coordStart+1    
! Sort out coords so they fall within alatxalatxalat box
    Do i=1,coordLength
      atomTypes(i) = atomTypesIn(coordStart+i-1,1)
      Do j=1,3
        atomCoords(i,j) = Modulus(atomCoordsIn(coordStart+i-1,j),aLat)
      End Do
    End Do    
! Call process subroutine
    Call makeNLProcess(nl, atomTypes, atomCoords, coordStart, coordEnd, rVerlet, aLat)    
  End Subroutine makeNL_B 
  
  Subroutine makeNLProcess(nl, atomTypes, atomCoords, coordStart, coordEnd, rVerlet, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords
    Integer(kind=StandardInteger) :: coordStart, coordEnd
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
    Dimension(1:((coordEnd-coordStart+2)*(coordEnd-coordStart+1))/2) :: nlUniqueKeyArr
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
    coordLength = coordEnd-coordStart+1
    rVerletSQ = rVerlet**2
    nl%aLat = aLat
    nl%totalRD = 0.0D0
    nl%totalRDSq = 0.0D0
    nl%rVerlet = rVerlet
    nlUniqueKeyArr = 0
    nlKey = 0
! Sort out coords so they fall within alatxalatxalat box
    !Do i=1,coordLength
    !  atomTypes(i) = atomTypesIn(coordStart+i-1)
    !  Do j=1,3
    !    atomCoords(i,j) = Modulus(atomCoordsIn(coordStart+i-1,j),aLat)
    !  End Do
    !End Do    
! calculate sub cell parameters
    scW = floor(aLat/(1.0D0*rVerlet))
    scCount = scW**3
    scAlat = aLat/(1.0D0*scW)
    maxAtomsPerSC = 5*ceiling(coordLength/(1.0D0*scCount))    
    nl%scCount = scCount
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
      scKey = SubCellKey(scAlat,scW,atomCoords(i,1),atomCoords(i,2),atomCoords(i,3))
      scAtomCount(scKey) = scAtomCount(scKey) + 1
! Store in sub cell arrays
      scCoordsI(scKey,scAtomCount(scKey),1) = i                 ! unique atom id 
      scCoordsI(scKey,scAtomCount(scKey),2) = atomTypes(i)      ! atom type                          
      Do j=1,3
        scCoordsR(scKey,scAtomCount(scKey),j) = atomCoords(i,j)
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
  End Subroutine makeNLProcess
  
  
  
  
  Subroutine updateNL(nl, atomCoords, coordStart)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart
! Private
    Integer(kind=StandardInteger) :: i, atomA_ID, atomB_ID
    Real(kind=DoubleReal) :: xD, yD, zD, rD, rDSq
! Update
    Do i=1,nl%length
! ID    
      atomA_ID = nl%i(i,1)+coordStart-1
      atomB_ID = nl%i(i,2)+coordStart-1   
! Make change in coord  
      nl%r(i,5) = atomCoords(atomA_ID,1)
      nl%r(i,6) = atomCoords(atomA_ID,2)
      nl%r(i,7) = atomCoords(atomA_ID,3)
      nl%r(i,8) = atomCoords(atomB_ID,1)+1.0D0*nl%subCell(i,1)*nl%aLat
      nl%r(i,9) = atomCoords(atomB_ID,2)+1.0D0*nl%subCell(i,2)*nl%aLat
      nl%r(i,10) = atomCoords(atomB_ID,3)+1.0D0*nl%subCell(i,3)*nl%aLat
! Seperation
      xD = nl%r(i,5)-nl%r(i,8)
      yD = nl%r(i,6)-nl%r(i,9)
      zD = nl%r(i,7)-nl%r(i,10)
      rDSq = xD**2+yD**2+zD**2
      rD = rDSq**0.5
! Store
      nl%r(i,1) = rD  
      nl%r(i,2) = xD/rD  
      nl%r(i,3) = yD/rD  
      nl%r(i,4) = zD/rD   
    End Do
  End Subroutine updateNL
  
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
! -----------------------------------------------
!        Module Test Subroutines
!
! -----------------------------------------------  
  
  
  
    Subroutine makeNLOriginal(nl, atomTypesIn, atomCoordsIn, coordStart, coordEnd, rcut, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    Integer(kind=StandardInteger), Dimension(:) :: atomTypesIn
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoordsIn  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
! Private variables  
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomCoords 
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1)) :: atomTypes
    Real(kind=DoubleReal) :: rcutsq
    Integer(kind=StandardInteger), &
    Dimension(1:(coordEnd-coordStart)*(coordEnd-coordStart)) :: &
    nlUniqueKeyArr
    Integer(kind=StandardInteger) :: i, j, l, m, n
    Integer(kind=StandardInteger) :: atomA, atomB, uKey, nlKey
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xD, yD, zD, rD
    Real(kind=DoubleReal) :: x1sq, x2sq, x3sq, rsq
    Real(kind=DoubleReal) :: rMinSq, rMaxSq
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Integer(kind=StandardInteger) :: coordLength    
! Allocate nl arrays    
    If(.not.Allocated(nl%i))Then
      Allocate(nl%i(1:100000,1:4))
    End If
    If(.not.Allocated(nl%r))Then
      Allocate(nl%r(1:100000,1:10))
    End If
    If(.not.Allocated(nl%cell))Then
      Allocate(nl%cell(1:100000,1:3))
    End If
! Init config specific variables
    rMinSq = 2.0D21
    rMaxSq = -2.0D21
    coordLength = coordEnd-coordStart+1
    rcutsq = rcut**2
    nl%aLat = aLat
    nl%totalRD = 0.0D0
    nl%totalRDSq = 0.0D0
! Sort out coords so they fall within alatxalatxalat box
    Do i=1,coordLength
      atomTypes(i) = atomTypesIn(coordStart+i-1)
      Do j=1,3
        atomCoords(i,j) = Modulus(atomCoordsIn(coordStart+i-1,j),aLat)
      End Do
    End Do    
! Loop through Atom B 3x3x3
    nlKey = 0
    Do l=-1,1
      Do m=-1,1
        Do n=-1,1
! Set co-ordinate shift
          xShift = aLat * l
          yShift = aLat * m
          zShift = aLat * n
! Reset unique key list
          nlUniqueKeyArr = 0
! Loop through atom pairs
          Do atomA=1,coordLength
            Do atomB=1,coordLength
              If(l.eq.0.and.m.eq.0.and.n.eq.0.and.atomA.eq.atomB)Then  ! Don't self count atom
              Else
! calculate the key of the atom A atom B combination
! half length list A=i,B=j == A=j,B=i
                If(atomA.gt.atomB)Then
                  uKey = (atomA-1)*(atomA)/2+atomB
                Else
                  uKey = (atomB-1)*(atomB)/2+atomA
                End If!
                If(nlUniqueKeyArr(uKey).eq.0)Then
                  nlUniqueKeyArr(uKey) = 1
                  xA = 1.0D0*atomCoords(atomA,1)
                  xB = 1.0D0*(xshift + atomCoords(atomB,1))
                  yA = 1.0D0*atomCoords(atomA,2)
                  yB = 1.0D0*(yshift + atomCoords(atomB,2))
                  zA = 1.0D0*atomCoords(atomA,3)
                  zB = 1.0D0*(zshift + atomCoords(atomB,3))
                  xD = xA-xB
                  x1sq = xD**2
                  If(x1sq.le.rcutsq)Then
                    yD = yA-yB
                    x2sq = yD**2
                    If(x2sq.le.rcutsq)Then
                      zD = zA-zB
                      x3sq = zD**2
                      If(x3sq.le.rcutsq)Then
                        rsq = x1sq + x2sq + x3sq
                        If(rsq.le.rcutsq)Then
                          rD = rsq**0.5
                          nl%totalRD = nl%totalRD + rD
                          nl%totalRDSq = nl%totalRDSq + rsq
                          nlKey = nlKey + 1
                          ! Key/Type
                          nl%i(nlKey,1) = atomA                ! A ID
                          nl%i(nlKey,2) = atomB                ! B ID
                          nl%i(nlKey,3) = atomTypes(atomA)  ! A type
                          nl%i(nlKey,4) = atomTypes(atomB)  ! B type
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
                          nl%cell(nlKey,1) = l
                          nl%cell(nlKey,2) = m
                          nl%cell(nlKey,3) = n
                          If(rsq.lt.rMinSq)Then
                            rMinSq = rsq
                          End If
                          If(rsq.gt.rMaxSq)Then
                            rMaxSq = rsq
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
    nl%rMin = rMinSq**0.5D0
    nl%rMax = rMaxSq**0.5D0
    nl%Length = nlKey
  End Subroutine makeNLOriginal
  
End Module geom