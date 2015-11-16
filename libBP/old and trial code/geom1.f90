Module geom
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: makeCoords
  Public :: makeNL
  Public :: makeNLTrial
  Public :: makeNLOld
  
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
  
  Subroutine makeNL(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat, nl_I, nl_R, nlLength)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
    Integer(kind=StandardInteger), Dimension(:,:) :: nl_I
    Real(kind=DoubleReal), Dimension(:,:) :: nl_R
    Integer(kind=StandardInteger) :: nlLength
! Private variables  
    Real(kind=DoubleReal) :: rcutsq
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomsA_I
    Integer(kind=StandardInteger), Dimension(1:9*(coordEnd-coordStart+1),1:3) :: atomsB_I  
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomsA_R
    Real(kind=DoubleReal), Dimension(1:9*(coordEnd-coordStart+1),1:3) :: atomsB_R
    Integer(kind=StandardInteger), &
    Dimension(1:(((coordEnd-coordStart)*(coordEnd-coordStart-1))/2+&
       (coordEnd-coordStart+1))) :: nlUniqueKeyArr
    Integer(kind=StandardInteger) :: i, j, coordID, atoms, atomACount, atomBCount
    Integer(kind=StandardInteger) :: k, l, m, n
    Integer(kind=StandardInteger) :: atomA, atomB, uKey, nlKey, atomA_ID, atomB_ID
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xD, yD, zD, rD
    Real(kind=DoubleReal) :: xMin, xMax
    Real(kind=DoubleReal) :: x1sq, x2sq, x3sq, rsq
    Real(kind=DoubleReal) :: rMin, rMax, rMinSq, rMaxSq
    
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Integer(kind=StandardInteger) :: coordLength
    
! Init config specific variables
    rMin = 2.0D21
    rMax = -2.0D21
    rMinSq = 2.0D21
    rMaxSq = -2.0D21
    coordLength = coordEnd-coordStart+1
    rcutsq = rcut**2
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
                          rD = rsq**0.5
                          nlKey = nlKey + 1
                          nl_I(nlKey,1) = atomTypes(atomA_ID)  ! A type
                          nl_I(nlKey,2) = atomTypes(atomB_ID)  ! B type
                          nl_I(nlKey,3) = atomA                ! A ID
                          nl_I(nlKey,4) = atomB                ! B ID
                          nl_I(nlKey,5) = atomA_ID             ! A ID in coords array
                          nl_I(nlKey,6) = atomB_ID             ! B ID in coords array
                          nl_R(nlKey,1) = rD
                          nl_R(nlKey,2) = xD/rD
                          nl_R(nlKey,3) = yD/rD
                          nl_R(nlKey,4) = zD/rD
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
    nlLength = nlKey
  End Subroutine makeNL
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
! Test functions/subroutines  
  
  
  Subroutine makeNLTrial(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
! Private variables  
    Real(kind=DoubleReal) :: rcutsq
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomsA_I
    Integer(kind=StandardInteger), Dimension(1:8*(coordEnd-coordStart+1),1:3) :: atomsB_I  
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomsA_R
    Real(kind=DoubleReal), Dimension(1:8*(coordEnd-coordStart+1),1:3) :: atomsB_R
    Integer(kind=StandardInteger), &
    Dimension(1:(((coordEnd-coordStart)*(coordEnd-coordStart-1))/2+&
       (coordEnd-coordStart+1))) :: uniqueIn
    Integer(kind=StandardInteger), &
    Dimension(1:(50*(coordEnd-coordStart)*(coordEnd-coordStart))) :: uniqueOut
    Integer(kind=StandardInteger) :: i, j, coordID, atoms, atomACount, atomBCount
    Integer(kind=StandardInteger) :: aID, bID
    Integer(kind=StandardInteger) :: k, l, m, n
    Integer(kind=StandardInteger) :: atomA, atomB, nlKey
    Real(kind=DoubleReal) :: x1, x2, x3
    Real(kind=DoubleReal) :: xMin, xMax
    Real(kind=DoubleReal) :: x1sq, x2sq, x3sq, rsq
    
    print *,size(uniqueIn,1)
    print *,size(uniqueOut,1)
    
! Init
    xMin = -1.0D0*rcut 
    xMax = aLat+1.0D0*rcut
    uniqueIn = 0
    uniqueOut = 0
    rcutsq = rcut**2
! Atom block A
    atoms = (coordEnd-coordStart+1)
    aID = 0
    Do i=coordStart,coordEnd
      aID = aID + 1
      atomsA_I(aID,1) = aID        ! Atom ID
      atomsA_I(aID,2) = atomTypes(i)   ! Atom type
      atomsA_I(aID,3) = i              ! Atom key in coord array
      atomsA_R(aID,1) = atomCoords(i,1)
      atomsA_R(aID,2) = atomCoords(i,2)
      atomsA_R(aID,3) = atomCoords(i,3)
    End Do
    atomACount = aID
! Halo Atom block B    
    bID = 0
    Do l=-1,1
      Do m=-1,1
        Do n=-1,1
          Do i=1,atomACount
            If(l.eq.0.and.m.eq.0.and.n.eq.0)Then
            Else
              x1 = atomsA_R(i,1)+1.0D0*l*aLat
              If(x1.ge.xMin.and.x1.le.xMax)Then
                x2 = atomsA_R(i,1)+1.0D0*m*aLat
                If(x2.ge.xMin.and.x2.le.xMax)Then
                  x3 = atomsA_R(i,1)+1.0D0*n*aLat
                  If(x3.ge.xMin.and.x3.le.xMax)Then
                    bID = bID + 1
                    atomsB_I(bID,1) = atomsA_I(i,1)   ! Atom ID
                    atomsB_I(bID,2) = atomsA_I(i,2)   ! Atom type
                    atomsB_I(bID,3) = 0               ! Outside A box
                    atomsB_R(bID,1) = x1
                    atomsB_R(bID,2) = x2
                    atomsB_R(bID,3) = x3
                  End If
                End If
              End If
            End If
          End Do 
        End Do
      End Do
    End Do
    atomBCount = bID
    
! Pair search in box
    k = 0
    Do i=1,atomACount
      Do j=1,atomACount
        If(atomsA_I(i,1).ne.atomsA_I(j,1))Then
          If(atomsA_I(i,1).lt.atomsA_I(j,1))Then
            nlKey = (atomsA_I(j,1)-1)*(atomsA_I(j,1)-2)/2+atomsA_I(i,1)
          Else
            nlKey = (atomsA_I(i,1)-1)*(atomsA_I(i,1)-2)/2+atomsA_I(j,1)
          End If    
          If(uniqueIn(nlKey).eq.0)Then
            uniqueIn(nlKey) = 1          
            x1sq = (atomsA_R(i,1)-atomsA_R(j,1))**2
            If(x1sq.le.rcutsq)Then       
              x2sq = (atomsA_R(i,2)-atomsA_R(j,2))**2
              If(x2sq.le.rcutsq)Then  
                x3sq = (atomsA_R(i,3)-atomsA_R(j,3))**2
                If(x3sq.le.rcutsq)Then
                  rsq = (x1sq+x2sq+x3sq)
                  If(rsq.le.rcutsq)Then
                    k = k + 1  
                  End If
                End If
              End If
            End If
          End If
        End If
      End Do
    End Do
    print *,k
! Pair search out of box
    Do i=1,atomACount
      Do j=1,atomBCount   
          If(i.lt.j)Then
            nlKey = (j-1)*(j-2)/2+i
          Else
            nlKey = (i-1)*(i-2)/2+j
          End If 
          If(uniqueOut(nlKey).eq.0)Then
            uniqueOut(nlKey) = 1     
            x1sq = (atomsA_R(i,1)-atomsB_R(j,1))**2
            If(x1sq.le.rcutsq)Then       
              x2sq = (atomsA_R(i,2)-atomsB_R(j,2))**2
              If(x2sq.le.rcutsq)Then  
                x3sq = (atomsA_R(i,3)-atomsB_R(j,3))**2
                If(x3sq.le.rcutsq)Then
                  rsq = (x1sq+x2sq+x3sq)
                  If(rsq.le.rcutsq)Then
                    k = k + 1   
                  End If
                End If
              End If
            End If
          End If
      End Do
    End Do
    
    
    
    !atomA
    print *,k,atomACount,atomBCount
  
  End Subroutine makeNLTrial
  
  
  
  
  Subroutine makeNLOld(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
! Private variables  
    Real(kind=DoubleReal) :: rcutsq
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomsA_I
    Integer(kind=StandardInteger), Dimension(1:9*(coordEnd-coordStart+1),1:3) :: atomsB_I  
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomsA_R
    Real(kind=DoubleReal), Dimension(1:9*(coordEnd-coordStart+1),1:3) :: atomsB_R
    Integer(kind=StandardInteger), &
    Dimension(1:(((coordEnd-coordStart)*(coordEnd-coordStart-1))/2+&
       (coordEnd-coordStart+1))) :: nlUniqueKeyArr
    Integer(kind=StandardInteger) :: i, j, coordID, atoms, atomACount, atomBCount
    Integer(kind=StandardInteger) :: k, l, m, n
    Integer(kind=StandardInteger) :: atomA, atomB, uKey, nlKey
    Real(kind=DoubleReal) :: x1, x2, x3
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB
    Real(kind=DoubleReal) :: xMin, xMax
    Real(kind=DoubleReal) :: x1sq, x2sq, x3sq, rsq
    Real(kind=DoubleReal) :: rMin, rMax, rMinSq, rMaxSq
    
    Real(kind=DoubleReal) :: xShift, yShift, zShift
    Integer(kind=StandardInteger) :: coordLength

    
! Init config specific variables
    rMin = 2.0D21
    rMax = -2.0D21
    rMinSq = 2.0D21
    rMaxSq = -2.0D21
    coordLength = coordEnd-coordStart+1
    rcutsq = rcut**2
    
    ! loop through Atom B 3x3x3
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
                  xA = 1.0D0*atomCoords(coordStart+atomA-1,1)
                  xB = 1.0D0*(xshift + atomCoords(coordStart+atomB-1,1))
                  yA = 1.0D0*atomCoords(coordStart+atomA-1,2)
                  yB = 1.0D0*(yshift + atomCoords(coordStart+atomB-1,2))
                  zA = 1.0D0*atomCoords(coordStart+atomA-1,3)
                  zB = 1.0D0*(zshift + atomCoords(coordStart+atomB-1,3))
                  x1sq = (xA-xB)**2
                  If(x1sq.le.rcutsq)Then
                    x2sq = (yA-yB)**2
                    If(x2sq.le.rcutsq)Then
                      x3sq = (zA-zB)**2
                      If(x3sq.le.rcutsq)Then
                        rsq = x1sq + x2sq + x3sq
                        If(rsq.le.rcutsq)Then
                          nlKey = nlKey + 1
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
    print *,nlKey,rMinSq,rMaxSq,(rMinSq**0.5),(rMaxSq**0.5)
  End Subroutine makeNLOld
  
  
  
  
  Subroutine makeNLA(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
! Private variables  
    Real(kind=DoubleReal) :: rcutsq
    Integer(kind=StandardInteger), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomsA_I
    Integer(kind=StandardInteger), Dimension(1:9*(coordEnd-coordStart+1),1:3) :: atomsB_I  
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: atomsA_R
    Real(kind=DoubleReal), Dimension(1:9*(coordEnd-coordStart+1),1:3) :: atomsB_R
    Integer(kind=StandardInteger), &
    Dimension(1:(((coordEnd-coordStart)*(coordEnd-coordStart-1))/2+&
       (coordEnd-coordStart+1))) :: uniqueIn, uniqueOut
    Integer(kind=StandardInteger), &
    Dimension(1:(((coordEnd-coordStart)*(coordEnd-coordStart-1))/2+&
       (coordEnd-coordStart+1))) :: nlUniqueKeyArr
    Integer(kind=StandardInteger) :: i, j, coordID, atoms, atomACount, atomBCount
    Integer(kind=StandardInteger) :: k, l, m, n
    Integer(kind=StandardInteger) :: atomA, atomB, nlKey
    Real(kind=DoubleReal) :: x1, x2, x3
    Real(kind=DoubleReal) :: xMin, xMax
    Real(kind=DoubleReal) :: x1sq, x2sq, x3sq, rsq
    print *,size(nlUniqueKeyArr,1)
! Init
    xMin = -1.0D0*rcut 
    xMax = aLat+1.0D0*rcut
    nlUniqueKeyArr = 0
    rcutsq = rcut**2
! Atom block A
    atoms = (coordEnd-coordStart+1)
    coordID = 0
    Do i=coordStart,coordEnd
      coordID = coordID + 1
      atomsA_I(coordID,1) = coordID        ! Atom ID
      atomsA_I(coordID,2) = atomTypes(i)   ! Atom type
      atomsA_I(coordID,3) = i              ! Atom key in coord array
      atomsA_R(coordID,1) = atomCoords(i,1)
      atomsA_R(coordID,2) = atomCoords(i,2)
      atomsA_R(coordID,3) = atomCoords(i,3)
    End Do
    atomACount = coordID
! Overlay Atom block B    
    coordID = 0
    Do l=-1,1
      Do m=-1,1
        Do n=-1,1
          Do i=1,atomACount
            If(l.eq.0.and.m.eq.0.and.n.eq.0)Then
              coordID = coordID + 1
              atomsB_I(coordID,1) = atomsA_I(i,1)     ! Atom ID
              atomsB_I(coordID,2) = atomsA_I(i,2)     ! Atom type
              atomsB_I(coordID,3) = 1                 ! Inside A box
              atomsB_R(coordID,1) = atomsA_R(i,1)
              atomsB_R(coordID,2) = atomsA_R(i,1)
              atomsB_R(coordID,3) = atomsA_R(i,1)
            Else
              x1 = atomsA_R(i,1)+1.0D0*l*aLat
              If(x1.ge.xMin.and.x1.le.xMax)Then
                x2 = atomsA_R(i,1)+1.0D0*m*aLat
                If(x2.ge.xMin.and.x2.le.xMax)Then
                  x3 = atomsA_R(i,1)+1.0D0*n*aLat
                  If(x3.ge.xMin.and.x3.le.xMax)Then
                    coordID = coordID + 1
                    atomsB_I(coordID,1) = atomsA_I(i,1)   ! Atom ID
                    atomsB_I(coordID,2) = atomsA_I(i,2)   ! Atom type
                    atomsB_I(coordID,3) = 0               ! Outside A box
                    atomsB_R(coordID,1) = x1
                    atomsB_R(coordID,2) = x2
                    atomsB_R(coordID,3) = x3
                  End If
                End If
              End If
            End If
          End Do 
        End Do
      End Do
    End Do
    atomBCount = coordID
    !Do i=1,atomACount
      !print *,i,atomsA_I(i,1),atomsA_I(i,2),atomsA_R(i,1),atomsA_R(i,2),atomsA_R(i,3)
    !End Do
    !Do i=1,atomBCount
      !print *,i,atomsB_I(i,1),atomsB_I(i,2),atomsB_R(i,1),atomsB_R(i,2),atomsB_R(i,3)
    !End Do
    
    
! Pair search
    k = 0
    Do i=1,atomACount
      Do j=1,atomBCount
        If(atomsB_I(j,3).eq.1.and.atomsA_I(i,1).eq.atomsB_I(j,1))Then
          ! Skip - same atom
        Else
        
          If(atomsA_I(i,1).lt.atomsB_I(j,1))Then
            nlKey = (atomsB_I(j,1)-1)*(atomsB_I(j,1)-2)/2+atomsA_I(i,1)
          Else
            nlKey = (atomsA_I(i,1)-1)*(atomsA_I(i,1)-2)/2+atomsB_I(j,1)
          End If            
          
          If(i.ge.250.and.j.gt.250)Then
            print *,"...",i,j,nlKey,nlUniqueKeyArr(nlKey),atomsB_I(j,1)
          End If
          
          If(nlUniqueKeyArr(nlKey).eq.0)Then
            nlUniqueKeyArr(nlKey) = 1          
            x1sq = (atomsA_R(i,1)-atomsB_R(j,1))**2
            If(x1sq.le.rcutsq)Then       
              x2sq = (atomsA_R(i,2)-atomsB_R(j,2))**2
              If(x2sq.le.rcutsq)Then  
                x3sq = (atomsA_R(i,3)-atomsB_R(j,3))**2
                If(x3sq.le.rcutsq)Then
                  rsq = (x1sq+x2sq+x3sq)
                  If(rsq.le.rcutsq)Then
                    k = k + 1     
                    If(i.ge.200.and.j.gt.200)Then                    
                      print *,i,j,atomsA_I(i,1),atomsB_I(j,1),rsq**0.5    
                    End If  
                  End If
                End If
              End If
            End If
          End If            
          
          !x1sq, x2sq, x3sq, rsq
          
          !atomA = atomsA_I(i)
          !atomB = atomsB_I(j)
          
          !nlUniqueKeyArr
        End If
      End Do
    End Do      
    !atomA
    print *,k,atomACount,atomBCount
  
  End Subroutine makeNLA
  
  
  
  
  
  
  
  
  
  

End Module geom