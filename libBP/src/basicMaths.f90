Module basicMaths
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
  Public :: RoundDP
  Public :: Factorial, FactorialDP, FactorialQ
  Public :: BinomialCoefficient, BinomialCoefficientDP, BinomialCoefficientQ
  Public :: Odd
  Public :: Even
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

  Function RoundDP(dpIn) RESULT (intOut)
! Round DP to nearest int
    Implicit None ! Force declaration of all variables
    Real(kind=DoubleReal) :: dpIn
    Integer(kind=StandardInteger) :: intOut
    intOut = Floor(dpIn+0.5D0)
  End Function RoundDP 
 
  Function Factorial(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=StandardInteger) :: output
! calculate factorial
    output = 1
    Do i=1,input
      output = i * output
    End Do
  End Function Factorial
  
  Function FactorialDP(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=VeryLongInteger) :: tempInt
    Real(kind=DoubleReal) :: output
! calculate factorial
    tempInt = 1
    Do i=1,input
      tempInt = i * tempInt
    End Do
    output = 1.0D0*tempInt
  End Function FactorialDP
  
  Function FactorialQ(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=VeryLongInteger) :: tempInt
    Real(kind=QuadrupoleReal) :: tempQ
    Real(kind=QuadrupoleReal) :: output
! calculate factorial
    tempInt = 1
    tempQ = 1
    Do i=1,input
      If(i.le.33)Then
        tempInt = i * tempInt
        tempQ = 1.0D0*tempInt
      End If
      If(i.eq.34)Then
        tempQ = 1.0D0*i*tempInt
      End If
      If(i.ge.35)Then
        tempQ = 1.0D0*i*tempQ
      End If
    End Do
    output = tempQ
  End Function FactorialQ

  Function BinomialCoefficient(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: c,n,k
! calculate factorial
    c = Factorial(n)/(Factorial(n-k)*Factorial(k))
  End Function BinomialCoefficient
  
  Function BinomialCoefficientDP(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: n,k
    Real(kind=DoubleReal) :: c, nDP, nkDP, kDP
! calculate factorial
    nDP = FactorialDP(n)
    nkDP = Factorial(n-k)
    kDP = FactorialDP(k)
    c = 1.0D0*nDP/(nkDP*kDP)
  End Function BinomialCoefficientDP

  Function BinomialCoefficientQ(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: n,k
    Real(kind=QuadrupoleReal) :: c, nDP, nkDP, kDP
! calculate factorial
    nDP = FactorialDP(n)
    nkDP = Factorial(n-k)
    kDP = FactorialDP(k)
    c = 1.0D0*nDP/(nkDP*kDP)
  End Function BinomialCoefficientQ
  
  Function Odd(input) RESULT (output)
! Returns true if odd, false if even
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: input
    Real(kind=DoubleReal) :: dpA, dpB
    Logical :: output
    output = .false.
    dpA = 1.0D0*input
    dpB = 2.0D0*ceiling(dpA/2.0D0)
    If(dpB.gt.dpA)Then
      output = .true.
    End If
  End Function Odd

  Function Even(input) RESULT (output)
! Returns true if even, false if odd
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: input
    Real(kind=DoubleReal) :: dpA, dpB
    Logical :: output
    output = .true.
    dpA = 1.0D0*input
    dpB = 2.0D0*ceiling(dpA/2.0D0)
    If(dpB.gt.dpA)Then
      output = .false.
    End If
  End Function Even
  

End Module basicMaths