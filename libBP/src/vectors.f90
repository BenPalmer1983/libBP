Module vectors
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
  Public :: CrossProduct
  Public :: DotProduct
  Public :: TripleProduct
  Public :: TripleProductSq
  Public :: ColToSquare
  Public :: SquareToCol
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------!
! Vector Functions
! ------------------------------------------------------------------------!

  Function CrossProduct(VectorA, VectorB) RESULT (VectorC)
! Calculates cross product of two vectors
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB, VectorC
! Calculate cross product
    VectorC(1) = VectorA(2)*VectorB(3)-VectorA(3)*VectorB(2)
    VectorC(2) = VectorA(3)*VectorB(1)-VectorA(1)*VectorB(3)
    VectorC(3) = VectorA(1)*VectorB(2)-VectorA(2)*VectorB(1)
  End Function CrossProduct

  Function DotProduct(VectorA, VectorB) RESULT (DotProductResult)
! Calculates dot product of two vectors
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB
    Real(kind=DoubleReal) :: DotProductResult
! Calculate dot product
    DotProductResult = VectorA(1)*VectorB(1)+VectorA(2)*VectorB(2)+VectorA(3)*VectorB(3)
  End Function DotProduct

  Function TripleProduct(VectorA, VectorB, VectorC) RESULT (TripleProductResult)
! Calculates (scalar) triple product of three vectors (resulting in volume of 3 vectors)
    Implicit None  ! Force declaration of all variables
! declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB, VectorC
    Real(kind=DoubleReal) :: TripleProductResult
! Calculate cross product
    TripleProductResult = DotProduct(VectorA,CrossProduct(VectorB,VectorC))
  End Function TripleProduct

  Function TripleProductSq(VectorIn) RESULT (TripleProductResult)
! Calculates (scalar) triple product of three vectors input in 3x3 matrix (resulting in volume of 3 vectors)
    Implicit None  ! Force declaration of all variables
! declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: VectorIn
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB, VectorC
    Real(kind=DoubleReal) :: TripleProductResult
! Calculate cross product
    Do i=1,3
      VectorA(i) = VectorIn(1,i)
      VectorB(i) = VectorIn(2,i)
      VectorC(i) = VectorIn(3,i)
    End Do
    TripleProductResult = DotProduct(VectorA,CrossProduct(VectorB,VectorC))
  End Function TripleProductSq  
   
! ------------------------------------------------------------------------!
! Unit Vector Functions
! ------------------------------------------------------------------------!

  Function ColToSquare(columnMatrix) RESULT (squareMatrix)
! Converts from PWscf 6 row unit vector to standard 3x3 unit vector
    Implicit None  ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(1:6) :: columnMatrix
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: squareMatrix
    Real(kind=DoubleReal) :: a, b, c, BC, AC, AB, cosBC, &
    cosAC, cosAB, sinBC, sinAC, sinAB
! Initialise variables
    squareMatrix = 0.0D0
! Input variables
    a = columnMatrix(1)
    b = a * columnMatrix(2)
    c = a * columnMatrix(3)
    cosBC = columnMatrix(4)
    cosAC = columnMatrix(5)
    cosAB = columnMatrix(6)
    BC = acos(cosBC)
    AC = acos(cosAC)
    AB = acos(cosAB)
    sinBC = sin(BC)
    sinAC = sin(AC)
    sinAB = sin(AB)
! Fill square matrix
    squareMatrix(1,1) = 1.0D0 * a
    squareMatrix(1,2) = 0.0D0
    squareMatrix(1,3) = 0.0D0
    squareMatrix(2,1) = 1.0D0*b*cosAB
    squareMatrix(2,2) = 1.0D0*b*sinAB
    squareMatrix(2,3) = 0.0D0
    squareMatrix(3,1) = 1.0D0*c*cosAC
    squareMatrix(3,2) = 1.0D0*c*(cosBC-cosAC*cosAB)/sinAB
    squareMatrix(3,3) = 1.0D0*((1.0D0+2.0D0*cosBC*cosAC*cosAB-&
    cosBC**2-cosAC**2-cosAB**2)**0.5)/sinAB
  End Function ColToSquare

  Function SquareToCol(squareMatrix) RESULT (columnMatrix)
! Convert from square matrix to PWscf unit vector
    Implicit None  ! Force declaration of all variables
! declare variables
    Real(kind=DoubleReal), Dimension(1:6) :: columnMatrix
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: squareMatrix
    Real(kind=DoubleReal) :: a, b, c, cosBC, cosAC, cosAB
! Convert
    a = sqrt(squareMatrix(1,1)**2+squareMatrix(1,2)**2+squareMatrix(1,3)**2)
    b = sqrt(squareMatrix(2,1)**2+squareMatrix(2,2)**2+squareMatrix(2,3)**2)/a
    c = sqrt(squareMatrix(3,1)**2+squareMatrix(3,2)**2+squareMatrix(3,3)**2)/a
    cosBC = squareMatrix(2,1)*squareMatrix(3,1)+squareMatrix(2,2)*squareMatrix(3,2)+&
    squareMatrix(2,3)*squareMatrix(3,3)
    cosAC = squareMatrix(1,1)*squareMatrix(3,1)+squareMatrix(1,2)*squareMatrix(3,2)+&
    squareMatrix(1,3)*squareMatrix(3,3)
    cosAB = squareMatrix(1,1)*squareMatrix(2,1)+squareMatrix(1,2)*squareMatrix(2,2)+&
    squareMatrix(1,3)*squareMatrix(2,3)
! Store to PWscf type column vector
    columnMatrix(1) = a
    columnMatrix(2) = b
    columnMatrix(3) = c
    columnMatrix(4) = cosBC
    columnMatrix(5) = cosAC
    columnMatrix(6) = cosAB
  End Function SquareToCol   

End Module vectors