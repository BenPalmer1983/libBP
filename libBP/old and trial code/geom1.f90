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
    Integer(kind=StandardInteger), Dimension(1:100000,1:4) :: i
    Real(kind=DoubleReal), Dimension(1:100000,1:10) :: r
    Integer(kind=StandardInteger), Dimension(1:100000,1:3) :: cell
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
  Use rng
  Use linearAlgebra
  Use geomTypes
! Force declaration of all variables
  Implicit None
! Public variables  
  Type(nlType) :: nl
  Real(kind=DoubleReal),Dimension(1:3) :: coords
! Make private
  Private
! ---- Variables
  Public :: nl
! ---- Subroutines
  Public :: makeCoords
  Public :: makeNL
  Public :: updateNL
  Public :: configOpt
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
  
  Subroutine makeNL(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
    !Type(nlType) :: nl
! Private variables  
    Real(kind=DoubleReal) :: rcutsq
    Integer(kind=StandardInteger), &
    Dimension(1:(((coordEnd-coordStart)*(coordEnd-coordStart-1))/2+&
       (coordEnd-coordStart+1))) :: nlUniqueKeyArr
    Integer(kind=StandardInteger) :: l, m, n
    Integer(kind=StandardInteger) :: atomA, atomB, uKey, nlKey, atomA_ID, atomB_ID
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xD, yD, zD, rD
    Real(kind=DoubleReal) :: x1sq, x2sq, x3sq, rsq
    Real(kind=DoubleReal) :: rMinSq, rMaxSq
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Integer(kind=StandardInteger) :: coordLength    
    Integer(kind=StandardInteger), Dimension(1:100) :: neighbourTally
! Init config specific variables
    rMinSq = 2.0D21
    rMaxSq = -2.0D21
    coordLength = coordEnd-coordStart+1
    rcutsq = rcut**2
    nl%aLat = aLat
    nl%totalRD = 0.0D0
    nl%totalRDSq = 0.0D0
    neighbourTally = 0
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
                If(atomA.lt.atomB)Then
                  uKey = (atomB-1)*(atomB-2)/2+atomA
                Else
                  uKey = (atomA-1)*(atomA-2)/2+atomB
                End If
!
                If(nlUniqueKeyArr(uKey).eq.0)Then
                  nlUniqueKeyArr(uKey) = 1
                  atomA_ID = coordStart+atomA-1
                  atomB_ID = coordStart+atomB-1
                  xA = 1.0D0*atomCoords(atomA_ID,1)
                  xB = 1.0D0*(xshift + atomCoords(atomB_ID,1))
                  yA = 1.0D0*atomCoords(atomA_ID,2)
                  yB = 1.0D0*(yshift + atomCoords(atomB_ID,2))
                  zA = 1.0D0*atomCoords(atomA_ID,3)
                  zB = 1.0D0*(zshift + atomCoords(atomB_ID,3))
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
                          neighbourTally(atomA_ID) = neighbourTally(atomA_ID) + 1
                          neighbourTally(atomB_ID) = neighbourTally(atomB_ID) + 1
                          rD = rsq**0.5
                          nl%totalRD = nl%totalRD + rD
                          nl%totalRDSq = nl%totalRDSq + rsq
                          nlKey = nlKey + 1
                          ! Key/Type
                          nl%i(nlKey,1) = atomA                ! A ID
                          nl%i(nlKey,2) = atomB                ! B ID
                          nl%i(nlKey,3) = atomTypes(atomA_ID)  ! A type
                          nl%i(nlKey,4) = atomTypes(atomB_ID)  ! B type
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
    Do n=1,coordEnd
      print *,neighbourTally(n)
    End Do  
  End Subroutine makeNL
    
  
  Subroutine updateNL(atomCoords, coordStart)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
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
    
    
  Subroutine efCalc()
! Energy Force calculation using L-J potential
    Implicit None   ! Force declaration of all variables  
! Private
    Integer(kind=StandardInteger) :: i
    Do i=1,nl%length
      !nl%r(i,1)
      !ljEnergy
      
      
    End Do  
    
    
    
  End Subroutine efCalc
    
    
    
    
    
    
  
  Subroutine configOpt(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
! Private
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomCoordsIn 
    Integer(kind=StandardInteger) :: i, j, k
    Real(kind=DoubleReal) :: y, yRef, h, dY
    Real(kind=DoubleReal), Dimension(1:3*(coordEnd-coordStart+1)) :: parameters, parametersIn  
    Integer(kind=StandardInteger), Dimension(1:3*(coordEnd-coordStart+1)) :: pMap, pMapR
    Integer(kind=StandardInteger) :: pSize    
    Real(kind=DoubleReal), Dimension(1:3*(coordEnd-coordStart+1)) :: gradArr, diffArr
! Init  
    i = 0  
    Do k=coordStart,coordEnd
      i = i + 1
      Do j=1,3
        atomCoordsIn(i,j) = atomCoords(k,j)
      End Do
    End Do
    h = 0.1D0
! Parameters array
    Call coordsToParameters(atomCoordsIn, parameters)
    parametersIn = parameters
! Make NL
    Call makeNL(atomTypes, atomCoordsIn, coordStart, coordEnd, rcut, aLat)
    yRef = nl%totalRDSq     
    print *,nl%length, nl%totalRD, nl%totalRDSq     
    
    Do k=1,30
    Do i=1,size(parameters,1)
      parameters = parametersIn    
      parameters(i) = parameters(i)+h
      Call parametersToCoords(parameters, atomCoordsIn)
      Call updateNL(atomCoordsIn, 1)
      !print *,i,nl%length, nl%totalRD, nl%totalRDSq, (nl%totalRDSq-yRef)
      If(nl%totalRDSq.lt.yRef)Then
        parametersIn = parameters
        yRef = nl%totalRDSq
      End If
    End Do
    End Do
    
    parameters = parametersIn   
    Call parametersToCoords(parameters, atomCoordsIn)
    Call updateNL(atomCoordsIn, 1)
    print *,i,nl%length, nl%totalRD, nl%totalRDSq, (nl%totalRDSq-yRef)

    
    
    


  End Subroutine configOpt
  
  Subroutine configOptProcess(pSize, parameters, pmr, gradIn, diffIn) 
! 
    Implicit None   ! Force declaration of all variables
! In/Out    
    Integer(kind=StandardInteger) :: pSize
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Integer(kind=StandardInteger), Dimension(:) :: pmr
    Real(kind=DoubleReal), Dimension(:) :: gradIn
    Real(kind=DoubleReal), Dimension(:) :: diffIn
! Private
    Real(kind=DoubleReal), Dimension(1:pSize) :: P
    Real(kind=DoubleReal), Dimension(1:pSize) :: grad
    Real(kind=DoubleReal), Dimension(1:pSize) :: diff
    Real(kind=DoubleReal), Dimension(1:1,1:pSize) :: J
    Real(kind=DoubleReal), Dimension(1:pSize) :: R
    Real(kind=DoubleReal), Dimension(1:pSize,1:1) :: JT
    Real(kind=DoubleReal), Dimension(1:pSize,1:pSize) :: JTJ
    Real(kind=DoubleReal), Dimension(1:pSize) :: JTR
    Integer(kind=StandardInteger) :: i
    
    Do i=1,pSize
      P(i) =  parameters(pmr(i))
      grad(i) =  gradIn(pmr(i)) 
      diff(i) =  diffIn(pmr(i)) 
      print *,i,P(i),grad(i),diff(i)      
      JT(i,1) = grad(i)
      J(1,i) = grad(i)
      R(i) = diff(i)
    End Do  
    
    !JTJ = matmul(JT,J)
    !JTR = matmul(JT,R)
    !JTR = -1.0D0*JTR
  
  End Subroutine configOptProcess
  
  Subroutine configOptCalc(parameters, y)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Real(kind=DoubleReal), Dimension(:) :: parameters  
    Real(kind=DoubleReal) :: y
! Private
    Real(kind=DoubleReal), Dimension(1:size(parameters)/3,1:3) :: atomCoords
    Integer(kind=StandardInteger) :: i, j, k
    k=0
    Do i=1,size(parameters)/3
      Do j=1,3
        k = k + 1
        atomCoords(i,j) = parameters(k)
      End Do
    End Do      
! Update neighbour list
    Call updateNL(atomCoords, 1)
    y = nl%totalRDSq
  End Subroutine configOptCalc
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  Subroutine coordsToParameters(coords, parameters)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Real(kind=DoubleReal), Dimension(:,:) :: coords  
    Real(kind=DoubleReal), Dimension(1:3*size(coords)) :: parameters
! Private
    Integer(kind=StandardInteger) :: i, j, k
! Convert
    k=0
    Do i=1,size(coords,1)
      Do j=1,3
        k = k + 1
        parameters(k) = coords(i,j)
      End Do
    End Do      
  End Subroutine coordsToParameters

  Subroutine parametersTocoords(parameters, coords)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal), Dimension(1:size(parameters)/3,1:3) :: coords  
! Private
    Integer(kind=StandardInteger) :: i, j, k
! Convert
    k=0
    Do i=1,size(parameters,1)/3
      Do j=1,3
        k = k + 1
        coords(i,j) = parameters(k)
      End Do
    End Do      
  End Subroutine parametersTocoords
  
  
  
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
    vr = 4*((sigma/r)**12-(sigma/r)**6)
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
    fr = (48/r)*((sigma/r)**12-(sigma/r)**6)
  End Function ljForce
  
End Module geom