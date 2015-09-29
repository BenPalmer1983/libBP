PROGRAM testProg
! University of Birmingham
! Ben Palmer
  Use maths
  Real(kind=DoubleReal), Dimension(1:2,1:2) :: aMatrix
  Real(kind=DoubleReal) :: aTr
  Integer(kind=StandardInteger), Dimension(1:2,1:2) :: bMatrix
  Integer(kind=StandardInteger) :: bTr
! Test str upper lower
  print *,StrToLower("Apples")
  print *,StrToUpper("Bananas")
! Test trace real 
  aMatrix(1,1) = 1.2D0
  aMatrix(1,2) = 2.1D0
  aMatrix(2,1) = 3.4D0
  aMatrix(2,2) = 2.6D0  
  aTr = Trace(aMatrix)
  print *,aTr  
! Test trace int 
  bMatrix(1,1) = 1
  bMatrix(1,2) = 2
  bMatrix(2,1) = 3
  bMatrix(2,2) = 4  
  bTr = Trace(bMatrix)
  print *,bTr
  
  print *, UnitConvert(5.2D0, "ANG", "M") 
  print *, UnitConvert(1.2D0, "ANG", "NM") 
  print *, UnitConvert(2.2D0, "ANG", "M") 
  print *, UnitConvert(2.2D0, "ANGS", "M") 
  print *, UnitConvert(2.2D0, "ANGS", "ANG") 
End Program testProg