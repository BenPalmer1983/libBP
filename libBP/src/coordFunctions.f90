Module coordFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use rngDist
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: TransformCoords
  Public :: RdCoords
  Public :: HeatCoords
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------


  Function TransformCoords (xVect, tVect) RESULT (xPVect)
! Transform coords with transformation matrix
    Implicit None  ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal), Dimension(1:3) :: xVect, xPVect
    Real(kind=DoubleReal), Dimension(1:3, 1:3) :: tVect
    xPVect(1) = xVect(1)*tVect(1,1)+&
    xVect(2)*tVect(1,2)+&
    xVect(3)*tVect(1,3)
    xPVect(2) = xVect(1)*tVect(2,1)+&
    xVect(2)*tVect(2,2)+&
    xVect(3)*tVect(2,3)
    xPVect(3) = xVect(1)*tVect(3,1)+&
    xVect(2)*tVect(3,2)+&
    xVect(3)*tVect(3,3)
  End Function TransformCoords

  Function RdCoords (xVect, yVect) RESULT (rd)
! Transform coords with transformation matrix
    Implicit None  ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal), Dimension(1:3) :: xVect, yVect
    Real(kind=DoubleReal) :: rd
    rd = sqrt((xVect(1)-yVect(1))**2+(xVect(2)-yVect(2))**2+(xVect(3)-yVect(3))**2)
  End Function RdCoords

  Function HeatCoords (inCoords, maxVar) RESULT (outCoords)
! Transform coords with transformation matrix
    Implicit None  ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal), Dimension(:,:) :: inCoords
    Real(kind=DoubleReal), Dimension(1:size(inCoords,1),1:size(inCoords,2)) :: outCoords
    Real(kind=DoubleReal) :: maxVar
    Integer(kind=StandardInteger) :: i, j
! Vary Coords
    Do i=1,size(inCoords,1)
      Do j=1,size(inCoords,2)
        outCoords(i,j) = inCoords(i,j) + 2.0D0 * (RandomDist("G") - 0.5D0) *  maxVar
      End Do
    End Do
  End Function HeatCoords

  
End Module coordFunctions