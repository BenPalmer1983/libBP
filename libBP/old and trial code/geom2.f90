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
  
  End Type  
  

End Module geomTypes


Module geom
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use rng
  Use geomTypes
! Force declaration of all variables
  Implicit None
! Public variables  
  !Integer(kind=LongInteger) :: randomLCG_n=0
! Make private
  Private
! ---- Variables
  Public :: randomLCG_n
! ---- Subroutines
  Public :: makeCoords
  Public :: makeNL
  Public :: updateNL
  Public :: bfgsRelax
  Public :: simpleRelax
  
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
  
  Subroutine makeNL(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat, nl_I, nl_R, nlLength, totalRD)
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
    Real(kind=DoubleReal) :: totalRD
! Private variables  
    Real(kind=DoubleReal) :: rcutsq
    Integer(kind=StandardInteger), &
    Dimension(1:(((coordEnd-coordStart)*(coordEnd-coordStart-1))/2+&
       (coordEnd-coordStart+1))) :: nlUniqueKeyArr
    Integer(kind=StandardInteger) :: l, m, n
    Integer(kind=StandardInteger) :: atomA, atomB, uKey, nlKey, atomA_ID, atomB_ID
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xD, yD, zD, rD
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
    totalRD = 0.0D0
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
                          totalRD = totalRD + rD
                          nlKey = nlKey + 1
                          nl_I(nlKey,1) = atomTypes(atomA_ID)  ! A type
                          nl_I(nlKey,2) = atomTypes(atomB_ID)  ! B type
                          nl_I(nlKey,3) = atomA                ! A ID
                          nl_I(nlKey,4) = atomB                ! B ID
                          nl_I(nlKey,5) = atomA_ID             ! A ID in coords array
                          nl_I(nlKey,6) = atomB_ID             ! B ID in coords array
                          nl_R(nlKey,1) = rD                          
                          nl_R(nlKey,2) = xD/rD                ! Vector from B to A (x)
                          nl_R(nlKey,3) = yD/rD                ! Vector from B to A (y)
                          nl_R(nlKey,4) = zD/rD                ! Vector from B to A (z)
                          nl_R(nlKey,5) = xA
                          nl_R(nlKey,6) = yA
                          nl_R(nlKey,7) = zA
                          nl_R(nlKey,8) = xB
                          nl_R(nlKey,9) = yB
                          nl_R(nlKey,10) = zB
                          
                          If(nlKey.lt.5)Then
                            !print *,nlKey,nl_I(nlKey,1),nl_I(nlKey,2),nl_R(nlKey,1),nl_R(nlKey,2),&
                            !nl_R(nlKey,3),nl_R(nlKey,4)
                            !print *,nlKey,nl_I(nlKey,1),nl_I(nlKey,2),nl_R(nlKey,5),nl_R(nlKey,6),&
                            !nl_R(nlKey,7),nl_R(nlKey,8),nl_R(nlKey,9),nl_R(nlKey,10)
                            !print *,""
                          End If  
                          
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
  
  
  Subroutine updateNL(atomCoordsChange, nl_I, nl_R, nlLength, totalRD)
! Make neighbour list for atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoordsChange  
    Integer(kind=StandardInteger), Dimension(:,:) :: nl_I
    Real(kind=DoubleReal), Dimension(:,:) :: nl_R
    Integer(kind=StandardInteger) :: nlLength
    Real(kind=DoubleReal) :: totalRD
! Private
    Integer(kind=StandardInteger) :: i, atomA_ID, atomB_ID
    Real(kind=DoubleReal) :: xD, yD, zD, rD
! Update
    totalRD = 0.0D0
    Do i=1,nlLength
! ID    
      atomA_ID = nl_I(i,5)
      atomB_ID = nl_I(i,6)       
! Make change in coord      
      nl_R(i,5) = nl_R(i,5)+atomCoordsChange(atomA_ID,1)
      nl_R(i,6) = nl_R(i,6)+atomCoordsChange(atomA_ID,2)
      nl_R(i,7) = nl_R(i,7)+atomCoordsChange(atomA_ID,3)
      nl_R(i,8) = nl_R(i,8)+atomCoordsChange(atomB_ID,1)
      nl_R(i,9) = nl_R(i,9)+atomCoordsChange(atomB_ID,2)
      nl_R(i,10) = nl_R(i,10)+atomCoordsChange(atomB_ID,3)
      xD = nl_R(i,5)-nl_R(i,8)
      yD = nl_R(i,6)-nl_R(i,9)
      zD = nl_R(i,7)-nl_R(i,10)      
      rD = (xD**2+yD**2+zD**2)**0.5
      !print *,rD,nl_R(i,1) 
      nl_R(i,1) = rD  
      nl_R(i,2) = xD/rD  
      nl_R(i,3) = yD/rD  
      nl_R(i,4) = zD/rD   
      totalRD = totalRD + rD
    End Do
  End Subroutine updateNL
  
  
  Subroutine bfgsRelax(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat)
! simple relaxation of atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
! Private
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger), Dimension(1:20000,1:6) :: nl_I
    Real(kind=DoubleReal), Dimension(1:20000,1:10) :: nl_R
    Integer(kind=StandardInteger) :: nlLength
    !Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: forces
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: velocity
    Real(kind=DoubleReal) :: timeInc, totalRD, randVar
    Real(kind=DoubleReal), &
    Dimension(1:size(atomCoords,1),1:size(atomCoords,2)) :: atomCoordsChange  
    
    
    timeInc = 0.1D0
    velocity = 0.0D0
    atomCoordsChange = 0.0D0
    Call makeNL(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat, nl_I, nl_R, nlLength, totalRD) 
    
    Do i=coordStart,coordEnd
      print *,i,atomCoords(coordStart+i-1,1),atomCoords(coordStart+i-1,2),atomCoords(coordStart+i-1,3)
    End Do
    
    
    
    
    print *,"0",nlLength, totalRD
    Do i=1,4
      Do j=coordStart,coordEnd
        Do k=1,3
          randVar = RandomLCG()
          atomCoordsChange(j,k) = 1.50D0*randVar
          !print *,j,atomCoordsChange(j,1),atomCoordsChange(j,2),atomCoordsChange(j,3)
        End Do 
      End Do
      
      Call updateNL(atomCoordsChange, nl_I, nl_R, nlLength, totalRD)
      print *,i ,nlLength, totalRD
      
      
      
      !Call simpleRelaxForce(forces, nl_I, nl_R, nlLength)
    End Do    
    
    !Function LMA(points, calcFunction, parametersIn, weightingIn, limitsLowerIn, limitsUpperIn) &
    !RESULT (parametersOut)
    
    !fx = calcFunction(x,parameters,size(parameters,1))
    
    
  End Subroutine bfgsRelax
    
  
  
  
  
  
  
  
  
  Subroutine simpleRelax(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat)
! simple relaxation of atoms
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal) :: rcut, aLat
! Private
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger), Dimension(1:20000,1:6) :: nl_I
    Real(kind=DoubleReal), Dimension(1:20000,1:10) :: nl_R
    Integer(kind=StandardInteger) :: nlLength
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: forces
    Real(kind=DoubleReal), Dimension(1:(coordEnd-coordStart+1),1:3) :: velocity
    Real(kind=DoubleReal) :: timeInc, totalRD
    
    timeInc = 0.1D0
    velocity = 0.0D0
    Do i=1,100
      Call makeNL(atomTypes, atomCoords, coordStart, coordEnd, rcut, aLat, nl_I, nl_R, nlLength, totalRD)    
      Call simpleRelaxForce(forces, nl_I, nl_R, nlLength)
      Call simpleRelaxMove(atomCoords, coordStart, coordEnd, forces, velocity, timeInc)
      !Call simpleRelaxOutputCoords(atomTypes, atomCoords, coordStart, coordEnd)
    End Do    
    
  End Subroutine simpleRelax
  
  
  Subroutine simpleRelaxForce(forces, nl_I, nl_R, nlLength)
! simple relaxation of atoms - force calculation subroutine
    Implicit None   ! Force declaration of all variables
! In/Out
    Real(kind=DoubleReal), Dimension(:,:) :: forces
    Integer(kind=StandardInteger), Dimension(:,:) :: nl_I
    Real(kind=DoubleReal), Dimension(:,:) :: nl_R
    Integer(kind=StandardInteger) :: nlLength
! Private
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: force, forceTotal
    Real(kind=DoubleReal) :: r 
! Init
    
    forces = 0.0D0
    
    Do i=1,nlLength
      r = nl_R(i,1)
      force = (2.0D0/(r**2)+1.0D0/(r**3))
      forces(nl_I(i,3),1) = forces(nl_I(i,3),1) + nl_R(i,2) * (-1.0D0) * force
      forces(nl_I(i,3),2) = forces(nl_I(i,3),2) + nl_R(i,3) * (-1.0D0) * force
      forces(nl_I(i,3),3) = forces(nl_I(i,3),3) + nl_R(i,4) * (-1.0D0) * force
      forces(nl_I(i,4),1) = forces(nl_I(i,4),1) + nl_R(i,2) * force
      forces(nl_I(i,4),2) = forces(nl_I(i,4),2) + nl_R(i,3) * force
      forces(nl_I(i,4),3) = forces(nl_I(i,4),3) + nl_R(i,4) * force
      If(i.lt.20)Then
        !-1.0D0*
        !print *,i,r,force
        !print *,i,nl_I(i,1),nl_I(i,2),nl_R(i,1),nl_R(i,2),&
        !nl_R(i,3),nl_R(i,4)
        !print *,i,nl_I(i,1),nl_I(i,2),nl_R(i,5),nl_R(i,6),&
        !nl_R(i,7),nl_R(i,8),nl_R(i,9),nl_R(i,10)
        !print *,""
      End If
    End Do    
    forceTotal = 0.0D0
    Do i=1,255
      !print *,i,forces(i,1),forces(i,2),forces(i,3)
      forceTotal = forceTotal + abs(forces(i,1)) + abs(forces(i,2)) + abs(forces(i,3))
    End Do   
    print *, forceTotal
  End Subroutine simpleRelaxForce
  
  
  Subroutine simpleRelaxMove(atomCoords, coordStart, coordEnd, forces, velocity, timeInc)
! simple relaxation of atoms - force calculation subroutine
    Implicit None   ! Force declaration of all variables
! In/Out
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
    Real(kind=DoubleReal), Dimension(:,:) :: forces    
    Real(kind=DoubleReal), Dimension(:,:) :: velocity    
    Real(kind=DoubleReal) :: timeInc
! Private
    Integer(kind=StandardInteger) :: i, j, k
    Real(kind=DoubleReal) :: maxForce
    
    Do i=1,size(forces,1)
      Do j=1,3 
        velocity(i,j) = velocity(i,j) + timeInc - forces(i,j)
      End Do
    End Do
    
    j = 0
    Do i=coordStart,coordEnd
      j = j + 1
      Do k=1,3
        atomCoords(i,k) = atomCoords(i,k) - forces(j,k) * timeInc**2 + velocity(j,k) * timeInc
      End Do   
    End Do
    
    
    
    
    maxForce = 0.0D0
    Do i=1,size(forces,1)
      Do j=1,3 
        If(i.eq.1.and.j.eq.1)Then
          maxForce = abs(forces(i,j))
        End If
        If(abs(forces(i,j)).gt.maxForce)Then
          maxForce = abs(forces(i,j))
        End If  
      End Do
    End Do      
    !print *,maxForce
    
    j = 0
    Do i=coordStart,coordEnd
      j = j + 1
      Do k=1,3
        !atomCoords(i,k) = atomCoords(i,k) - 0.2D0 * forces(j,k)
      End Do   
    End Do
  
  End Subroutine simpleRelaxMove
  
  Subroutine simpleRelaxOutputCoords(atomTypes, atomCoords, coordStart, coordEnd)
! simple relaxation of atoms - force calculation subroutine
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger), Dimension(:) :: atomTypes
    Real(kind=DoubleReal), Dimension(:,:) :: atomCoords  
    Integer(kind=StandardInteger) :: coordStart, coordEnd
! Private
    Integer(kind=StandardInteger) :: i
    
    print *,(coordEnd-coordStart+1)
    print *,"Comment"
    Do i=coordStart,coordEnd
      print *,atomTypes(i), atomCoords(i,1), atomCoords(i,2), atomCoords(i,3)
    End Do
  
  End Subroutine simpleRelaxOutputCoords
  
  

End Module geom