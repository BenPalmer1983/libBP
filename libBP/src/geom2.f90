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
    !Integer(kind=StandardInteger), Dimension(1:100000,1:4) :: i
    !Real(kind=DoubleReal), Dimension(1:100000,1:10) :: r
    !Integer(kind=StandardInteger), Dimension(1:100000,1:3) :: cell
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: i
    Real(kind=DoubleReal), Allocatable, Dimension(:,:) :: r
    Integer(kind=StandardInteger), Allocatable, Dimension(:,:) :: cell
    
    Integer(kind=StandardInteger) :: length = 0
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal) :: rMin
    Real(kind=DoubleReal) :: rMax
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
  Public :: makeNL
  Public :: updateNL
  Public :: mdRun
  !Public :: configOpt
  !Public :: bfgsRelax
  !Public :: simpleRelax
  
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  
  
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
! 
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
  
  
  
  
  
  
  
  
  
  !makeNL(atomID, coords, 1, 32, 6.5D0,2*4.04D0, nl)
  
  Subroutine makeNL(nl, atomTypesIn, atomCoordsIn, coordStart, coordEnd, rcut, aLat)
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
    Integer(kind=StandardInteger) :: atomA, atomB, uKey, nlKey, atomA_ID, atomB_ID
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
    !print *,nl%Length
  End Subroutine makeNL
  
  
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
    nl%totalRD = 0.0D0
    nl%totalRDSq = 0.0D0
    Do i=1,nl%length
! ID    
      atomA_ID = nl%i(i,1)+coordStart-1
      atomB_ID = nl%i(i,2)+coordStart-1       
! Make change in coord  
      nl%r(i,5) = atomCoords(atomA_ID,1)
      nl%r(i,6) = atomCoords(atomA_ID,2)
      nl%r(i,7) = atomCoords(atomA_ID,3)
      nl%r(i,8) = atomCoords(atomB_ID,1)+1.0D0*nl%cell(i,1)*nl%aLat
      nl%r(i,9) = atomCoords(atomB_ID,2)+1.0D0*nl%cell(i,2)*nl%aLat
      nl%r(i,10) = atomCoords(atomB_ID,3)+1.0D0*nl%cell(i,3)*nl%aLat
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
      nl%totalRD = nl%totalRD + rD
      nl%totalRDSq = nl%totalRDSq + rDSq
    End Do
  End Subroutine updateNL
  
  
  Subroutine efsCalc(nl,forceArr)
! Calc energy/stress/force using lj potential - no units
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    Real(kind=DoubleReal), Dimension(:,:) :: forceArr
! Private    
    Real(kind=DoubleReal) :: ljSigma, energyVal, forceVal, rd
    Real(kind=DoubleReal) :: fX, fY, fZ
    Integer(kind=StandardInteger) :: i, nlKey, atomA, atomB
    ljSigma = 4.0D0
    forceArr = 0.0D0
    Do nlKey=1,nl%length
      atomA = nl%i(nlKey,1)
      atomB = nl%i(nlKey,2)      
      rd = nl%r(nlKey,1)
      !energyVal = ljEnergy(ljSigma,rd)
      forceVal = ljForce(ljSigma, rd)
      fX = forceVal*nl%r(nlKey,2)
      fY = forceVal*nl%r(nlKey,3)
      fZ = forceVal*nl%r(nlKey,4)
! Store forces
      forceArr(atomA,1) = forceArr(atomA,1)+fX
      forceArr(atomA,2) = forceArr(atomA,2)+fY
      forceArr(atomA,3) = forceArr(atomA,3)+fZ
      forceArr(atomB,1) = forceArr(atomB,1)-fX
      forceArr(atomB,2) = forceArr(atomB,2)-fY
      forceArr(atomB,3) = forceArr(atomB,3)-fZ      
    End Do
  
  End Subroutine efsCalc
  
  Subroutine mdMove(nl, atomCoords, forceArr, velocityArr, timeInc)
    Implicit None   ! Force declaration of all variables
! In/Out
    Type(nlType) :: nl
    Real(kind=DoubleReal), Dimension(:,:) :: forceArr
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Real(kind=DoubleReal), Dimension(:,:) :: velocityArr
    Real(kind=DoubleReal) :: timeInc
    Real(kind=DoubleReal) :: dV, v0, vAvg, s0, dS
    Real(kind=DoubleReal) :: vFactor
! Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal), Dimension(1:size(forceArr,1),1:size(forceArr,2)) :: forceArrLast
! Init    
    forceArrLast = forceArr
! Calculate force on atoms
    Call efsCalc(nl,forceArr)
! Mass = 1, acceleration = force
    Do i=1,size(atomCoords,1)
      Do j=1,3
        vFactor = 1.0D0
        If(CompareSign(forceArr(i,j),forceArrLast(i,j)))Then
          vFactor = 0.5D0
        End If
! Velocity
        v0 = velocityArr(i,j)
        dV = timeInc * forceArr(i,j)
        vAvg = v0 + 0.5D0*dV        
        velocityArr(i,j) =  vFactor*(v0 + dV)
! Displacement        
        s0 = atomCoords(i,j)
        dS = 0.5D0 * forceArr(i,j) * timeInc**2 + vAvg * timeInc
        atomCoords(i,j) = Modulus(s0 + dS,nl%alat)
      End Do
    End Do
! Update neighbour list
    !Call updateNL(nl, atomCoords, 1)
  End Subroutine mdMove
  
  
  Subroutine mdRun(atomTypesIn, atomCoordsIn, coordStart, coordEnd, rcut, aLat, timeSteps, timeInc)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypesIn
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoordsIn  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
    Integer(kind=StandardInteger) :: timeSteps
    Real(kind=DoubleReal) :: timeInc
! Private
    Type(nlType) :: nl
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1)) :: atomTypes
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomCoords  
    Integer(kind=StandardInteger) :: coordCount
    Integer(kind=StandardInteger) :: i, j, k
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: forceArr
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: velocityArr
 
    print *,size(atomTypesIn,1)
    print *,size(atomCoordsIn,1)
    print *,nl%length
    
    print *, atomCoordsIn(1,1), atomCoordsIn(1,2), atomCoordsIn(1,3)
     
! Init      
    velocityArr = 0.0D0
    coordCount = (coordEnd-coordStart+1)
! Load atom coords from input array
    k = coordStart   
    Do i=1,coordCount
      atomTypes(i) = atomTypesIn(k)
      Do j=1,3
        atomCoords(i,j) = atomCoordsIn(k,j)
      End Do
      k = k + 1
    End Do
    
    Call makeNL(nl, atomTypes, atomCoords, 1, coordCount, rcut, aLat)   
    Call efsCalc(nl,forceArr)
    
    Do i=1,coordCount
      print *,i,atomCoords(i,1),atomCoords(i,2),atomCoords(i,3),&
        forceArr(i,1),forceArr(i,2),forceArr(i,3),nl%length
    End Do

    
    Do j=1,timeSteps
      Call makeNL(nl, atomTypes, atomCoords, 1, coordCount, rcut, aLat)    
      Call mdMove(nl, atomCoords, forceArr, velocityArr, timeInc)   
      !If(j.ge.8.and.j.le.16)Then
      !print *,""
      !Do i=1,nl%length
      !  If(nl%i(i,1).eq.1.or.nl%i(i,2).eq.1)Then
      !    print *,nl%i(i,1),nl%i(i,2),nl%r(i,1)
      !  End If  
      !End Do
      !print *,""
      !End If
    End Do
    
    print *,""
    print *,""
    print *,""
    
    Do i=1,coordCount
      print *,i,atomCoords(i,1),atomCoords(i,2),atomCoords(i,3),&
      forceArr(i,1),forceArr(i,2),forceArr(i,3),nl%length
    End Do
    
    
    
    
    
  
  End Subroutine mdRun
  
  
  
  
  
  
  
  
  
  
! Functions
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  Subroutine makeNLTrial(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat)
! Make neighbour list for atoms
! Assumes no interferrence i.e. rcut<0.5*alat
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
    !Type(nlType) :: nl
! Private variables  
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: coordsA
    Real(kind=DoubleReal), Dimension(1:27*(coordEnd-coordStart+1),1:3) :: coordsB
    Real(kind=DoubleReal), Dimension(1:27*(coordEnd-coordStart+1),1:3) :: cell
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1),1:2) :: coordsA_I 
    Integer(kind=StandardInteger), Dimension(1:27*(coordEnd-coordStart+1),1:2) :: coordsB_I 
    Integer(kind=StandardInteger) :: aCount, bCount
    Integer(kind=StandardInteger) :: i, j, k
    
    Real(kind=DoubleReal) :: xMin, xMax    
    Integer(kind=StandardInteger) :: atomA, atomB
    Real(kind=DoubleReal) :: rcutsq
    Real(kind=DoubleReal) :: xB, yB, zB
    Real(kind=DoubleReal) :: xD, yD, zD
    Real(kind=DoubleReal) :: xDsq, yDsq, zDsq, rDsq
    Integer(kind=StandardInteger) :: uKey
    Logical :: savePair
    

    !Integer(kind=StandardInteger), &
    !Dimension(1:(((coordEnd-coordStart)*(coordEnd-coordStart-1))/2+&
    !   (coordEnd-coordStart+1))) :: nlUniqueKeyArr
    Integer(kind=StandardInteger), Dimension(1:100000) :: nlU_In
    Integer(kind=StandardInteger), Dimension(1:100000) :: nlU_Out
    Integer(kind=StandardInteger) :: l, m, n
    Integer(kind=StandardInteger) :: atomB_ID
    !Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xD, yD, zD, rD
    !Real(kind=DoubleReal) :: x1sq, x2sq, x3sq, rsq
    !Real(kind=DoubleReal) :: rMinSq, rMaxSq
    !Real(kind=DoubleReal) :: xShift, yShift, zShift
    Integer(kind=StandardInteger), Dimension(1:100) :: neighbourTally
    
! Init    
    xMin = (-1.0D0)*rcut
    xMax = (1.0D0)*rcut+aLat
    rcutsq = rcut**2
    neighbourTally = 0
    nlU_In = 0
    nlU_Out = 0
! Cube A
    i = 0
    Do k=coordStart,coordEnd
      i = i+1
      Do j=1,3
        coordsA_I(i,1) = atomTypes(k)
        coordsA(i,j) = atomCoords(k,j)
      End Do
    End Do
    aCount = i    
    print *,aCount
    
! Cube B
    k = 0
    Do l=-1,1
      Do m=-1,1
        Do n=-1,1
          Do i=1,aCount 
            xB = coordsA(i,1) + l*aLat
            If(xB.ge.xMin.and.xB.le.xMax)Then
              yB = coordsA(i,2) + m*aLat
              If(yB.ge.xMin.and.yB.le.xMax)Then
                zB = coordsA(i,3) + n*aLat
                If(zB.ge.xMin.and.zB.le.xMax)Then
                  k = k + 1
                  coordsB(k,1) = xB
                  coordsB(k,2) = yB
                  coordsB(k,3) = zB
                  cell(k,1) = l
                  cell(k,2) = m
                  cell(k,3) = n
                  coordsB_I(k,1) = coordsA_I(i,1)
                  coordsB_I(k,2) = i
                End If
              End If    
            End If                
          End Do  
        End Do  
      End Do  
    End Do  
    bCount = k    
    print *,bCount
    
    k = 0
    Do atomA=1,aCount    
      Do atomB=1,bCount      
        If(cell(atomB,1).eq.0.and.cell(atomB,2).eq.0.and.&
        cell(atomB,3).eq.0.and.atomA.eq.coordsB_I(atomB,2))Then
          !Skip
        Else        
          savePair = .false.
          If(atomA.gt.atomB)Then
            uKey = (atomA-1)*(atomA)/2+coordsB_I(atomB,2)
          Else
            uKey = (coordsB_I(atomB,2)-1)*(coordsB_I(atomB,2))/2+atomA
          End If
          If(cell(atomB,1).eq.0.and.cell(atomB,2).eq.0.and.cell(atomB,3).eq.0)Then
            If(nlU_In(uKey).eq.0)Then
              savePair = .true.
              nlU_In(uKey) = 1
            End If
          Else  
            If(nlU_Out(uKey).eq.0)Then
              savePair = .true.
              nlU_Out(uKey) = 1
            End If
          End If
          
          
          !If(savePair)Then
            xD = coordsA(atomA,1)-coordsB(atomB,1)
            xDsq = xD**2
            If(xDsq.le.rcutsq)Then
              yD = coordsA(atomA,2)-coordsB(atomB,2)
              yDsq = yD**2
              If(yDsq.le.rcutsq)Then
                zD = coordsA(atomA,3)-coordsB(atomB,3)
                zDsq = zD**2
                If(zDsq.le.rcutsq)Then
                  rDsq = xDsq+yDsq+zDsq
                  If(rDsq.le.rcutsq)Then
                    k = k + 1
                    
                    atomB_ID = coordsB_I(atomB,2)
                    neighbourTally(atomA) = neighbourTally(atomA) + 1
                    neighbourTally(atomB_ID) = neighbourTally(atomB_ID) + 1
                    
                  End If
                End If
              End If
            End If
          !End If         
          
        End If
      End Do
    End Do
    print *,k
    
    Do i=1,aCount
      print *,i,neighbourTally(i)
    End Do
  
  End Subroutine makeNLTrial
  
  
  
  
  
  
End Module geom