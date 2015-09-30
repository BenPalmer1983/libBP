Module functionPoints
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use calcFunctions
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: MakeFunctionPoints
  Public :: PolyPoints
  Public :: BirchMurnPoints
  
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

  Function MakeFunctionPoints(calcFunction,parameters,xStart,xEnd,pointCount) RESULT (dataPoints)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: xStart, xEnd
    Integer(kind=StandardInteger) :: pointCount
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:pointCount,1:2) :: dataPoints
! Private: Declare variables
    Real(kind=DoubleReal), External :: calcFunction
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: xInc, xVal, yVal
! Set values
    xInc = (xEnd-xStart)/(1.0D0*(pointCount-1.0D0))
! Loop through points
    Do i=1,pointCount
      xVal = xStart+1.0D0*(i-1)*xInc
      yVal = calcFunction(parameters, xVal)
      dataPoints(i,1) = xVal
      dataPoints(i,2) = yVal
    End Do
  End Function MakeFunctionPoints
  
  Function PolyPoints(parameters,xStart,xEnd,pointCount) RESULT (dataPoints)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: xStart, xEnd
    Integer(kind=StandardInteger) :: pointCount
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:pointCount,1:2) :: dataPoints
! Private: Declare variables
    Real(kind=DoubleReal), External :: calcFunction
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: xInc, xVal, yVal
! Set values
    xInc = (xEnd-xStart)/(1.0D0*(pointCount-1.0D0))
! Loop through points
    Do i=1,pointCount
      xVal = xStart+1.0D0*(i-1)*xInc
      yVal = CalcPolynomial(parameters, xVal)
      dataPoints(i,1) = xVal
      dataPoints(i,2) = yVal
    End Do
  End Function PolyPoints
  
  Function BirchMurnPoints(parameters,xStart,xEnd,pointCount) RESULT (dataPoints)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: xStart, xEnd
    Integer(kind=StandardInteger) :: pointCount
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:pointCount,1:2) :: dataPoints
! Private: Declare variables
    Real(kind=DoubleReal), External :: calcFunction
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: xInc, xVal, yVal
! Set values
    xInc = (xEnd-xStart)/(1.0D0*(pointCount-1.0D0))
! Loop through points
    Do i=1,pointCount
      xVal = xStart+1.0D0*(i-1)*xInc
      yVal = BirchMurnCalc(xVal, parameters)
      dataPoints(i,1) = xVal
      dataPoints(i,2) = yVal
    End Do
  End Function BirchMurnPoints


! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------  
  

End Module functionPoints