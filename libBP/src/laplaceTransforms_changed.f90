Module laplaceTransforms
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use basicMaths
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --functions--!
  Public :: GaverStehfestWeighting
  Public :: GaverStehfestWeightingQ
  Public :: GaverStehfest

  Public :: ltDecay
  Public :: ltExp

! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -------------------------------------------------
! FUNCTIONS
! -------------------------------------------------

  Function GaverStehfestWeighting(N, weightingIn) RESULT (weighting)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: N
    Integer(kind=StandardInteger) :: j, k, jStart, jEnd
    Real(kind=DoubleReal) :: factor, wSum
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
    weighting = 0.0_QuadrupoleReal
    factor = 0.0_QuadrupoleReal
    wSum = 0.0_QuadrupoleReal
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


  Function GaverStehfest(funcInS, t, p, mIn) RESULT (ft)
! Gaver Stehfest
    Implicit None ! Force declaration of all variables
! Vars In
    Real(kind=DoubleReal), External :: funcInS
    Real(kind=DoubleReal) :: t
    Real(kind=DoubleReal), Dimension(1:10) :: p
    Integer(kind=StandardInteger), Optional :: mIn
! Vars Out
    Real(kind=DoubleReal) :: ft
! Vars Private
    Integer(kind=StandardInteger) :: k, m
    Real(kind=DoubleReal) :: s, fs
    Real(kind=DoubleReal), Dimension(1:20) :: w
! Make weighting array
    m = 7
    If(Present(mIn))Then
      m = mIn
    End If
    w = GaverStehfestWeighting(m, w)
! Loop through and approximate ft
    ft = 0.0D0
    Do k=1,2*m
      s = (1.0D0*k*lnTwo)/t
      fs = funcInS(s,p)
      ft = ft + w(k)*fs
    End Do
    ft = (lnTwo/t)*ft
  End Function GaverStehfest

! -------------------------------------------------
! SUBROUTINES
! -------------------------------------------------










! -------------------------------------------------
! FUNCTIONS for GS
! -------------------------------------------------

  Function ltDecay(s, p) RESULT (fs)
! First decay calculation function
    Implicit None ! Force declaration of all variables
! Vars In
    Real(kind=DoubleReal) :: s
    Real(kind=DoubleReal), Dimension(1:100) :: p
! Vars Out
    Real(kind=DoubleReal) :: fs
! Vars Private
    Integer(kind=StandardInteger) :: i, n, k
! p(1)             isotope to calculate
! p(2)             total isotopes in chain
! p(2+4*(n-1)+1)   production rate
! p(2+4*(n-1)+2)   decay constant
! p(2+4*(n-1)+3)   Branch Factor
! p(2+4*(n-1)+4)   Start Amount
    print *,"ltDecay"
! Calculate F(s) for parent
    fs = (1.0D0/(s+p(4)))*(p(3)/s+p(6))
    print *,p(4)
    print *,p(3)
    print *,p(6)
    print *,""
! Calculate F(s) until desired isotope (parent/nth daughter)
    n = nint(p(1))
    Do i=2,n
      k = 2+4*(i-1)
      !fs = (1.0D0/(s+p(k+2)))*(p(k+1)/s+p(k+2)*p(k-2)*fs+p(k+3))
      fs = (1.0D0/(s+p(k+2)))*(p(k+3)*p(k-3)*fs+p(k+4))

      print *,i,(k+1),p(k+1)
      print *,i,(k+2),p(k+2)
      print *,i,(k+3),p(k+3)
      print *,i,(k+4),p(k+4)
      print *,""
    End Do

    !fs = p(1)*(1.0D0/(s*(s+p(2))))+p(3)*(1.0D0/(s+p(2)))
  End Function ltDecay





  Function ltExp(s, p) RESULT (fs)
! Exp
    Implicit None ! Force declaration of all variables
! Vars In
    Real(kind=DoubleReal) :: s
    Real(kind=DoubleReal), Dimension(1:10) :: p
    !Integer(kind=StandardInteger) :: nParameters
! Vars Out
    Real(kind=DoubleReal) :: fs
! Calculate
    fs = (1.0D0)/(s-p(1))
  End Function ltExp














End Module laplaceTransforms
