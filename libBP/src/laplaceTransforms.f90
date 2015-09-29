Module laplaceTransforms
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use basicMaths
! Force declaration of all variables
  Implicit None
! Public variables  
! Make private
  Private
! Public
! --functions--!
  Public :: GaverStehfestWeighting, GaverStehfestWeightingQ
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

  Function GaverStehfestWeighting(N, weightingIn) RESULT (weighting)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: N
    Integer(kind=StandardInteger) :: j, k, jStart, jEnd
    Real(kind=DoubleReal) :: factor, wSum
! Real(kind=DoubleReal), Dimension(1:2*N) :: weighting
    Real(kind=DoubleReal), Dimension(:) :: weightingIn
    Real(kind=DoubleReal), Dimension(1:size(weightingIn)) :: weighting
! Init array
    weighting = 0.0D0
! k loop
    Do k=1,2*N
      factor = (-1)**(k+N)/(1.0D0*FactorialDP(N))
      jStart = Floor((k+1)/2.0D0)
      jEnd = min(k,N)
      wSum = 0.0D0
! j loop
      Do j=jStart,jEnd
        wSum = wSum + 1.0D0*(j**(N+1))*BinomialCoefficientDP(N,j)*&
        BinomialCoefficientDP(2*j,j)*BinomialCoefficientDP(j,k-j)
      End Do
      weighting(k) = factor*wSum
    End Do
  End Function GaverStehfestWeighting
  
  Function GaverStehfestWeightingQ(N, weightingIn) RESULT (weighting)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: N
    Integer(kind=StandardInteger) :: j, k, jStart, jEnd
    Real(kind=QuadrupoleReal) :: factor, wSum
    Real(kind=QuadrupoleReal), Dimension(:) :: weightingIn
    Real(kind=QuadrupoleReal), Dimension(1:size(weightingIn)) :: weighting
! Init array
    weighting = 0.0D0
! k loop
    Do k=1,2*N
      factor = (-1)**(k+N)/(1.0D0*FactorialQ(N))
      jStart = Floor((k+1)/2.0D0)
      jEnd = min(k,N)
      wSum = 0.0D0
! j loop
      Do j=jStart,jEnd
        wSum = wSum + 1.0D0*(j**(N+1))*BinomialCoefficientQ(N,j)*&
        BinomialCoefficientQ(2*j,j)*BinomialCoefficientQ(j,k-j)
      End Do
      weighting(k) = factor*wSum
    End Do
  End Function GaverStehfestWeightingQ


End Module laplaceTransforms