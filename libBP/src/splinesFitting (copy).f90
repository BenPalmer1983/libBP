Module splinesFitting
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Fitting module: BM, exp etc
! Splines: Poly splines, exp(poly) splines
! --------------------------------------------------------------!
  Use kinds
  Use strings
  Use constants
  Use matrix
  Use linearAlgebra
  Use rng
  Use calcFunctions
  Use regression
  Use interpolation
  Use lmaM
! Force declaration of all variables
  Implicit None
! Public variables
  Character(Len=128), Dimension(1:20) :: fittingReport
  Real(kind=DoubleReal), Dimension(1:1000,1:10) :: splineCoeffs
! Make private
  Private
! Public
! --variables--!
  Public :: fittingReport
  Public :: splineCoeffs
! --- Functions - Fitting
  Public :: BirchMurnFit
  Public :: ExpFit
  Public :: SingleDecayFit
  Public :: DoubleDecayFit
  Public :: TripleDecayFit
  Public :: FittingPoints
  Public :: FitEmbeddingA
  Public :: FitEmbeddingB
  Public :: FitEmbeddingC
  Public :: FitDensity
! --- Functions - Spline
  Public :: SplineAB
  Public :: SplineExpThird
  Public :: SplineExpFifth
  Public :: SplineNodes
  Public :: SplineComplete
  Public :: VaryNode
  Public :: SplinePoints
! --- Subroutines
  Public :: CompleteNodeData
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------
! Birch-Murnaghan Fitting with LMA
! ----------------------------------------------------------------------------------

  Function BirchMurnFit(points, bp0Lower_In, bp0Upper_In) RESULT (coefficients)
! Fit Murnaghan EoS to data
! Fitting method adapted from http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/appendices/appendix-eos.html
! Birch-Murnaghan equation described by Hebbache 2004
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Optional :: bp0Lower_In
    Real(kind=DoubleReal), Optional :: bp0Upper_In
! Out
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients   ! E0, V0, B0, B'0
! Private variables
    Real(kind=DoubleReal) :: randDouble
    Real(kind=DoubleReal), Dimension(1:3) :: coefficientsQ
    Real(kind=DoubleReal) :: bp0Lower
    Real(kind=DoubleReal) :: bp0Upper
! LMA Vars
    Real(kind=DoubleReal), Dimension(1:4) :: parameters
    Real(kind=DoubleReal), Dimension(1:4) :: upper, lower
! Optional arguments
    bp0Lower = 2.0D0
    If(Present(bp0Lower_In))Then
      bp0Lower = bp0Lower_In
    End If
    bp0Upper = 9.0D0
    If(Present(bp0Upper_In))Then
      bp0Upper = bp0Upper_In
    End If
! Quadratic fit  (as a starting point)
    coefficientsQ = PolyFit(points,2)
! Random number
    randDouble = RandomLCG()
! Starting values for fit
    coefficients(2) = (-1.0D0*coefficientsQ(2))/(2.0D0*coefficientsQ(3))    !V0
    coefficients(1) = coefficientsQ(3)*coefficients(2)**2+&                 !E0
    coefficientsQ(2)*coefficients(2)+&
    coefficientsQ(1)
    coefficients(3) = 2.0D0 * coefficientsQ(3) * coefficients(2)            !B0
    coefficients(4) = 4.0D0 + 2.0D0 * randDouble                            !B'0
! --------------------------------------------------
! LMA
! --------------------------------------------------
    coefficientsQ(1) = coefficients(1)
    coefficientsQ(2) = coefficients(2)
    coefficientsQ(3) = coefficients(3)
! Set limits on bulk property pressure derivative
    lower = 1.0D0
    upper = -1.0D0
    lower(4) = bp0Lower
    upper(4) = bp0Upper
! Fit parameters
    parameters = LMA(points, LMA_BirchMurn, coefficients,.true.,lower,upper)
! store LMA if better
    coefficients = parameters
  End Function BirchMurnFit

! ----------------------------------------------------------------------------------
! EXP fitting - single term, double terms, triple terms
! ----------------------------------------------------------------------------------

  Function ExpFit(dataPoints,terms) RESULT (parameters)
    Implicit None   ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Integer(kind=StandardInteger) :: terms
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:(2*terms)) :: parameters
    parameters = 0.0D0
    If(terms.eq.1)Then
      parameters = SingleDecayFit(dataPoints)
    End If
    If(terms.eq.2)Then
      parameters = DoubleDecayFit(dataPoints)
    End If
    If(terms.eq.3)Then
      parameters = TripleDecayFit(dataPoints)
    End If
  End Function ExpFit

  Function SingleDecayFit(dataPoints) RESULT (parameters)
! Fit aexp(mx) to data
    Implicit None   ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:2) :: parameters
! Private: Declare variables
    Integer(kind=StandardInteger) :: n
    Real(kind=DoubleReal) :: lnY
    Real(kind=DoubleReal) :: sumY, sumX_Y, sumX_X_Y
    Real(kind=DoubleReal) :: sumY_LnY, sumX_Y_LnY
! Linear regression
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:2) :: yMatrix, cMatrix
! LMA
    Integer(kind=StandardInteger) :: i
! --------------------
! Linear regression
! --------------------
    sumY = 0.0D0
    sumX_Y = 0.0D0
    sumX_X_Y = 0.0D0
    sumY_LnY = 0.0D0
    sumX_Y_LnY = 0.0D0
    Do n=1,size(dataPoints,1)
      lnY = log(dataPoints(n,2))
      sumY = sumY + dataPoints(n,2)
      sumX_Y = sumX_Y + dataPoints(n,1)*dataPoints(n,2)
      sumX_X_Y = sumX_X_Y + dataPoints(n,1)*dataPoints(n,1)*dataPoints(n,2)
      sumY_LnY = sumY_LnY + dataPoints(n,2)*lnY
      sumX_Y_LnY = sumX_Y_LnY + dataPoints(n,1)*dataPoints(n,2)*lnY
    End Do
! x matrix
    xMatrix(1,1) = sumY
    xMatrix(1,2) = sumX_Y
    xMatrix(2,1) = sumX_Y
    xMatrix(2,2) = sumX_X_Y
! y matrix
    yMatrix(1) = sumY_LnY
    yMatrix(2) = sumX_Y_LnY
! solve
    xMatrix = InvertMatrix(xMatrix)
    cMatrix = matmul(xMatrix,yMatrix)
! --------------------
! LMA
! --------------------
! set parameters
    parameters(1) = exp(cMatrix(1))
    parameters(2) = cMatrix(2)
    Do i=1,2
      If(isnan(parameters(i)))Then  ! If fitting fails
        parameters(i) = 1.0D0
      End If
    End Do
! Fit parameters
    parameters = LMA(dataPoints, LMA_Exp, parameters)
  End Function SingleDecayFit

  Function DoubleDecayFit(dataPoints) RESULT (parameters)
! Fit Murnaghan EoS to data
! Fitting method adapted from http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/appendices/appendix-eos.html
! Birch-Murnaghan equation described by Hebbache 2004
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
! Out
    Real(kind=DoubleReal), Dimension(1:4) :: parameters
! Private
    Integer(kind=StandardInteger) :: n
! Lin Reg approx vars
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: S
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: SS
    Real(kind=DoubleReal), Dimension(1:4) :: cMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: xMatrix
    Real(kind=DoubleReal) :: sumX, sumY, sumX_X, sumX_Y, sumS_X, sumS_Y, sumSS_X, sumSS_Y
    Real(kind=DoubleReal) :: sumS, sumSS, sumS_S, sumS_SS, sumSS_SS
    Real(kind=DoubleReal) :: p1, q1
    Real(kind=DoubleReal), Dimension(1:2) :: cbMatrix, ybMatrix
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: xbMatrix
    Real(kind=DoubleReal) :: sumB_B, sumB_N, sumN_N, sumB_Y, sumN_Y
    Real(kind=DoubleReal) :: beta, eta
! LMA
    Integer(kind=StandardInteger) :: i
! -----------------------------------------
! Approximate linear regression
! -----------------------------------------
! Advice from Claude Leibovici
! Book by Jean Jacquelin:  https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
!
! Init vars
    sumX = 0.0D0
    sumY = 0.0D0
    sumX_X = 0.0D0
    sumX_Y = 0.0D0
    sumS_X = 0.0D0
    sumS_Y = 0.0D0
    sumSS_X = 0.0D0
    sumSS_Y = 0.0D0
    sumS = 0.0D0
    sumSS = 0.0D0
    sumS_S = 0.0D0
    sumS_SS = 0.0D0
    sumSS_SS = 0.0D0
! Numeric integration to calc S and SS array
    n = size(dataPoints,1)      ! number of data points
    Do i=1,n
      If(i.eq.1)Then
        S(i) = 0.0D0
        SS(i) = 0.0D0
      Else
        S(i) = S(i-1) + 0.5D0*(dataPoints(i,2)+dataPoints(i-1,2))*&
        (dataPoints(i,1)-dataPoints(i-1,1)) ! Numeric integration
        SS(i) = SS(i-1) + 0.5D0*(S(i)+S(i-1))*&
        (dataPoints(i,1)-dataPoints(i-1,1)) ! Numeric integration
      End If
    End Do
! Sum
    Do i=1,n
      sumX = sumX + dataPoints(i,1)
      sumY = sumY + dataPoints(i,2)
      sumX_X = sumX_X + (dataPoints(i,1)*dataPoints(i,1))
      sumX_Y = sumX_Y + (dataPoints(i,1)*dataPoints(i,2))
      sumS_X = sumS_X + S(i)*dataPoints(i,1)
      sumS_Y = sumS_Y + S(i)*dataPoints(i,2)
      sumSS_X = sumSS_X + SS(i)*dataPoints(i,1)
      sumSS_Y = sumSS_Y + SS(i)*dataPoints(i,2)
      sumS = sumS + S(i)
      sumSS = sumSS + SS(i)
      sumS_S = sumS_S + S(i)*S(i)
      sumS_SS = sumS_SS + S(i)*SS(i)
      sumSS_SS = sumSS_SS + SS(i)*SS(i)
    End Do
! Make y matrix
    yMatrix(1) = sumSS_Y
    yMatrix(2) = sumS_Y
    yMatrix(3) = sumX_Y
    yMatrix(4) = sumY
! Make xMatrix
    xMatrix(1,1) = sumSS_SS
    xMatrix(1,2) = sumS_SS
    xMatrix(1,3) = sumSS_X
    xMatrix(1,4) = sumSS
    xMatrix(2,1) = sumS_SS
    xMatrix(2,2) = sumS_S
    xMatrix(2,3) = sumS_X
    xMatrix(2,4) = sumS
    xMatrix(3,1) = sumSS_X
    xMatrix(3,2) = sumS_X
    xMatrix(3,3) = sumX_X
    xMatrix(3,4) = sumX
    xMatrix(4,1) = sumSS
    xMatrix(4,2) = sumS
    xMatrix(4,3) = sumX
    xMatrix(4,4) = n
! Solve set of equations
    xMatrix = InvertMatrix(xMatrix)
    cMatrix = matmul(xMatrix,yMatrix)
! calculate P and Q for next regression
    p1 = 0.5D0*(cMatrix(2)+sqrt(cMatrix(2)*cMatrix(2)+4*cMatrix(1)))
    q1 = 0.5D0*(cMatrix(2)-sqrt(cMatrix(2)*cMatrix(2)+4*cMatrix(1)))
! Sum
    sumB_B = 0.0D0
    sumB_N = 0.0D0
    sumN_N = 0.0D0
    sumB_Y = 0.0D0
    sumN_Y = 0.0D0
    Do i=1,n
      beta = exp(p1*dataPoints(i,1))
      eta = exp(q1*dataPoints(i,1))
      sumB_B = sumB_B + beta*beta
      sumB_N = sumB_N + beta*eta
      sumN_N = sumN_N + eta*eta
      sumB_Y = sumB_Y + beta*dataPoints(i,2)
      sumN_Y = sumN_Y + eta*dataPoints(i,2)
    End Do
! Make next x matrix
    xbMatrix(1,1) = sumB_B
    xbMatrix(1,2) = sumB_N
    xbMatrix(2,1) = sumB_N
    xbMatrix(2,2) = sumN_N
! Make next y matrix
    ybMatrix(1) = sumB_Y
    ybMatrix(2) = sumN_Y
! Calc cb
    xbMatrix = InvertMatrix(xbMatrix)
    cbMatrix = matmul(xbMatrix,ybMatrix)
! --------------------
! LMA
! --------------------
! set parameters
    parameters(1) = cbMatrix(1)
    parameters(2) = p1
    parameters(3) = cbMatrix(2)
    parameters(4) = q1
    Do i=1,4
      If(isnan(parameters(i)))Then  ! If fitting fails
        parameters(i) = 1.0D0
      End If
    End Do
! Fit parameters
    parameters = LMA(dataPoints, LMA_Exp, parameters)
  End Function DoubleDecayFit

  Function TripleDecayFit(dataPoints) RESULT (parameters)
! Fits double exponential to data
! f(x) = a exp(lA) + b exp(lB)
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
! Out
    Real(kind=DoubleReal), Dimension(1:6) :: parameters
! Private
    Integer(kind=StandardInteger) :: i, n
! Approx regression
    Integer(kind=StandardInteger) :: iA, iB, iC, gridSize, searchBetterFailCount
    Real(kind=DoubleReal) :: maxLambda, minLambda
    Real(kind=DoubleReal) :: lambdaInc
    Real(kind=DoubleReal) :: lA, lB, lC
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: yReg
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:3) :: xReg
    Real(kind=DoubleReal), Dimension(1:3) :: linCoeffs
    Real(kind=DoubleReal) :: linRSS, linBestRSS, bestRSSLastSearch
! Init vars
    parameters = 0.0D0
    linBestRSS = 0.0D0
    bestRSSLastSearch = 0.0D0
! Set grid size
    gridSize = 15
! Make y array
    Do i=1,size(dataPoints,1)
      yReg(i) = dataPoints(i,2)
    End Do
! Loop through combinations
    maxLambda = 0.5D0
    searchBetterFailCount = 0
    Do n=1,20   ! expand search "area"
      maxLambda = 2.0D0*maxLambda   ! grid from -maxLambda to maxLambda
      minLambda = -1.0D0 * maxLambda
      lambdaInc = (2*maxLambda)/(gridSize-1)
      Do iA=1,gridSize-2
        lA = minLambda + (iA-1)*lambdaInc
        Do i=1,size(dataPoints,1)
          xReg(i,1) = exp(lA*dataPoints(i,1))  ! Make x1 array (for the A exp(lA x) function)
        End Do
        Do iB=iA+1,gridSize-1
          lB = minLambda + (iB-1)*lambdaInc
          Do i=1,size(dataPoints,1)
            xReg(i,2) = exp(lB*dataPoints(i,1))  ! Make x1 array (for the A exp(lA x) function)
          End Do
          Do iC=iB+1,gridSize
            lC = minLambda + (iC-1)*lambdaInc
            Do i=1,size(dataPoints,1)
              xReg(i,3) = exp(lC*dataPoints(i,1))
            End Do
            linCoeffs = LinearRegression(xReg, yReg)
            linRSS = TripleDecayFitRSS(dataPoints, linCoeffs(1), lA, linCoeffs(2), lB, linCoeffs(3), lC)
            If(iA.eq.1.and.iB.eq.2.and.iC.eq.3)Then
              linBestRSS = linRSS
              parameters(1) = linCoeffs(1)
              parameters(2) = lA
              parameters(3) = linCoeffs(2)
              parameters(4) = lB
              parameters(5) = linCoeffs(3)
              parameters(6) = lC
            Else
              If(linRSS.lt.linBestRSS)Then
                linBestRSS = linRSS
                parameters(1) = linCoeffs(1)
                parameters(2) = lA
                parameters(3) = linCoeffs(2)
                parameters(4) = lB
                parameters(5) = linCoeffs(3)
                parameters(6) = lC
              End If
            End If
          End Do
        End Do
      End Do
      If(linBestRSS.gt.bestRSSLastSearch)Then
        searchBetterFailCount = searchBetterFailCount + 1
        If(searchBetterFailCount.eq.2)Then ! If two successive fails, break out
          Exit
        End If
      Else
        searchBetterFailCount = 0 ! reset fail count
        bestRSSLastSearch = linBestRSS
      End If
    End Do
! --------------------------------------------------
! LMA
! --------------------------------------------------
    Do i=1,6
      If(isnan(parameters(i)))Then  ! If fitting fails
        parameters(i) = 1.0D0
      End If
    End Do
! Fit parameters
    parameters = LMA(dataPoints, LMA_Exp, parameters)
  End Function TripleDecayFit














  Function FitEmbeddingA(dataPoints, startPoint, endPoint) Result (coefficients)
! Fit embedding functional to data points
! E(p) = a + bp^0.5 + cp^2
!
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Integer(kind=StandardInteger) :: startPoint, endPoint
! Out
    Real(kind=DoubleReal), Dimension(1:2) :: coefficients
! Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1),1:2) :: X
    Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1)) :: Y
! Adjust to fit in data array size
    If(startPoint.lt.1)Then
      startPoint = 1
    End If
    If(endPoint.gt.size(dataPoints,1).or.endPoint.lt.1)Then
      endPoint = size(dataPoints,1)
    End If
! Prepare matrices
    j = 0
    Do i=startPoint,endPoint
      j = j + 1
      X(j,1) = 1.0D0
      X(j,2) = 1.0D0*dataPoints(i,1)**0.5D0
      Y(j) = 1.0D0*dataPoints(i,2)
    End Do
! Calculate
    coefficients = LinearRegression(X,Y)
  End Function FitEmbeddingA
! ----------------------------------------------------------------------------------
  Function FitEmbeddingB(dataPoints, startPoint, endPoint) Result (coefficients)
! Fit embedding functional to data points
! E(p) = a + bp^0.5 + cp^2
!
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Integer(kind=StandardInteger) :: startPoint, endPoint
! Out
    Real(kind=DoubleReal), Dimension(1:3) :: coefficients
! Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1),1:3) :: X
    Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1)) :: Y
! Adjust to fit in data array size
    If(startPoint.lt.1)Then
      startPoint = 1
    End If
    If(endPoint.gt.size(dataPoints,1).or.endPoint.lt.1)Then
      endPoint = size(dataPoints,1)
    End If
! Prepare matrices
    j = 0
    Do i=startPoint,endPoint
      j = j + 1
      X(j,1) = 1.0D0
      X(j,2) = 1.0D0*dataPoints(i,1)**0.5D0
      X(j,3) = 1.0D0*dataPoints(i,1)**2.0D0
      Y(j) = 1.0D0*dataPoints(i,2)
    End Do
! Calculate
    coefficients = LinearRegression(X,Y)
  End Function FitEmbeddingB
! ----------------------------------------------------------------------------------
  Function FitEmbeddingC(dataPoints, startPoint, endPoint) Result (coefficients)
! Fit embedding functional to data points
! E(p) = a + bp^0.5 + cp^2 + dp^4
!
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Integer(kind=StandardInteger) :: startPoint, endPoint
! Out
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
! Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1),1:4) :: X
    Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1)) :: Y
! Adjust to fit in data array size
    If(startPoint.lt.1)Then
      startPoint = 1
    End If
    If(endPoint.gt.size(dataPoints,1).or.endPoint.lt.1)Then
      endPoint = size(dataPoints,1)
    End If
! Prepare matrices
    j = 0
    Do i=startPoint,endPoint
      j = j + 1
      X(j,1) = 1.0D0
      X(j,2) = 1.0D0*dataPoints(i,1)**0.5D0
      X(j,3) = 1.0D0*dataPoints(i,1)**2.0D0
      X(j,4) = 1.0D0*dataPoints(i,1)**4.0D0
      Y(j) = 1.0D0*dataPoints(i,2)
    End Do
! Calculate
    coefficients = LinearRegression(X,Y)
  End Function FitEmbeddingC


  Function FitDensity(dataPoints, startPointIn, endPointIn) Result (coefficients)
! Fit embedding functional to data points
!
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Integer(kind=StandardInteger), Optional :: startPointIn, endPointIn
    Integer(kind=StandardInteger) :: startPoint, endPoint
! Out
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
! Private
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: pointCount
! optional arguments
    startPoint = 1
    endPoint = size(dataPoints,1)
    If(Present(startPointIn))Then
      startPoint = startPointIn
    End If
    If(Present(endPointIn))Then
      endPoint = endPointIn
    End If
! Adjust to fit in data array size
    If(startPoint.lt.1)Then
      startPoint = 1
    End If
    If(endPoint.gt.size(dataPoints,1))Then
      endPoint = size(dataPoints,1)
    End If
! pointCount
    pointCount = 0
    Do i=startPoint,endPoint
      If(dataPoints(i,1).ne.0.0D0)Then
        pointCount = pointCount + 1
      End If
    End Do
! run subroutine
    Call FitDensity_Process(dataPoints, pointCount, startPoint, endPoint, coefficients)
  End Function FitDensity
! ----------
  Subroutine FitDensity_Process(dataPoints, pointCount, startPoint, endPoint, coefficients)
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Integer(kind=StandardInteger) :: pointCount, startPoint, endPoint
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
! Private
    Real(kind=DoubleReal), Dimension(1:pointCount,1:2) :: fitPoints, fitPointsLinReg
    Integer(kind=StandardInteger) :: i, n, k
    Real(kind=DoubleReal), Dimension(1:6) :: parameters, parameters_Last
! Lin reg
    Integer(kind=StandardInteger) :: iA, iB, gridSize, searchBetterFailCount
    Real(kind=DoubleReal) :: maxLambda, minLambda
    Real(kind=DoubleReal) :: lambdaInc
    Real(kind=DoubleReal) :: lA, lB
    Real(kind=DoubleReal), Dimension(1:pointCount) :: yReg
    Real(kind=DoubleReal), Dimension(1:pointCount,1:2) :: xReg
    Real(kind=DoubleReal), Dimension(1:2) :: linCoeffs
    Real(kind=DoubleReal) :: linRSS, linBestRSS, bestRSSLastSearch
! p(r) = a r^2 exp(-br^2) + c r^2 exp(-dr^2)
    i = 0
    Do k=startPoint,endPoint
      If(dataPoints(k,1).ne.0.0D0)Then
        i = i + 1
        fitPoints(i,1) = dataPoints(k,1)
        fitPoints(i,2) = dataPoints(k,2)
        fitPointsLinReg(i,1) = (dataPoints(k,1)**2)
        fitPointsLinReg(i,2) = dataPoints(k,2)/(dataPoints(k,1)**2)
      End If
    End Do
! Linear reg to find starting point
! fit p(r)/r^2 = a exp(-br^2) + c exp(-dr^2)
! Init vars
    parameters = 0.0D0
    parameters_Last = 0.0D0
    linBestRSS = 0.0D0
    bestRSSLastSearch = 0.0D0
! Set grid size
    gridSize = 15
! Make y array
    Do i=1,pointCount
      yReg(i) = fitPointsLinReg(i,2)
    End Do
! Loop through combinations
    maxLambda = 0.5D0
    searchBetterFailCount = 0
    Do n=1,20   ! expand search "area"
      maxLambda = 2.0D0*maxLambda   ! grid from -maxLambda to maxLambda
      minLambda = -1.0D0 * maxLambda
      lambdaInc = (2*maxLambda)/(gridSize-1)
      Do iA=1,gridSize-1
        lA = minLambda + (iA-1)*lambdaInc
        Do i=1,pointCount
          xReg(i,1) = exp(lA*fitPointsLinReg(i,1))    ! Make x1 array (for the A exp(lA x) function)
        End Do
        Do iB=iA+1,gridSize
          lB = minLambda + (iB-1)*lambdaInc
          Do i=1,pointCount
            xReg(i,2) = exp(lB*fitPointsLinReg(i,1))  ! Make x1 array (for the A exp(lA x) function)
          End Do
          linCoeffs = LinearRegression(xReg,yReg)
          linRSS = DoubleDecayFitRSS(fitPointsLinReg, linCoeffs(1), lA, linCoeffs(2), lB)
          If(iA.eq.1.and.iB.eq.2)Then
            linBestRSS = linRSS
            parameters(1) = linCoeffs(1)
            parameters(2) = lA
            parameters(3) = linCoeffs(2)
            parameters(4) = lB
          Else
            If(linRSS.lt.linBestRSS)Then
              linBestRSS = linRSS
              parameters(1) = linCoeffs(1)
              parameters(2) = lA
              parameters(3) = linCoeffs(2)
              parameters(4) = lB
            End If
          End If
        End Do
      End Do
      If(linBestRSS.gt.bestRSSLastSearch)Then
        searchBetterFailCount = searchBetterFailCount + 1
        If(searchBetterFailCount.eq.2)Then ! If two successive fails, break out
          Exit
        End If
      Else
        searchBetterFailCount = 0 ! reset fail count
        bestRSSLastSearch = linBestRSS
      End If
    End Do
! --------------------------------------------------
! LMA
! --------------------------------------------------
    Do i=1,4
      If(isnan(parameters(i)))Then  ! If fitting fails
        parameters(i) = 1.0D0
      End If
    End Do
! Fit parameters
    parameters = LMA(fitPoints, LMA_ExpDens, parameters)
! Store results
    Do i=1,4
      coefficients(i) = parameters(i)
    End Do
  End Subroutine FitDensity_Process



! ----------------------------------------------------------------------------------
! Fitting Points
! Use fitting functions, regression functions
!
!
!
!
! ----------------------------------------------------------------------------------

  Function FittingPoints(dataPointsIn, calcFunction, pointsOutCount, optArgA, optArgB) Result (dataPointsOut)
! Fit to input data, and give set of points for the fit function
    Implicit None   ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: dataPointsIn
    Character(*) :: calcFunction
    Integer(kind=StandardInteger) :: pointsOutCount
    Real(kind=DoubleReal), Optional :: optArgA, optArgB
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:pointsOutCount,1:2) :: dataPointsOut
! Private: Declare variables
    Integer(kind=StandardInteger) :: i, pointsInCount
    Character(Len=12) :: calcFunctionT
    Real(kind=DoubleReal), Dimension(1:2) :: parameters2
    Real(kind=DoubleReal), Dimension(1:3) :: parameters3
    Real(kind=DoubleReal), Dimension(1:4) :: parameters4
    Real(kind=DoubleReal), Dimension(1:5) :: parameters5
    Real(kind=DoubleReal), Dimension(1:6) :: parameters6
    Real(kind=DoubleReal) :: xStart, xEnd, xInc
    Real(kind=DoubleReal) :: argA, argB
! Optional arguments
    argA = 0.0D0
    argB = 0.0D0
    If(Present(optArgA))Then
      argA = optArgA
    End If
    If(Present(optArgB))Then
      argB = optArgB
    End If
! Init
    pointsInCount = size(dataPointsIn,1)
    calcFunctionT = "            "
    calcFunctionT = Trim(AdjustL(StrToUpper(calcFunction)))
    parameters2 = 0.0D0
    parameters3 = 0.0D0
    parameters4 = 0.0D0
    parameters5 = 0.0D0
    parameters6 = 0.0D0
    fittingReport = BlankStringArray(fittingReport)
! Start/End x
    xStart = dataPointsIn(1,1)
    xEnd = dataPointsIn(pointsInCount,1)
    xInc = (xEnd-xStart)/(pointsOutCount-1.0D0)
! Poly fit
    If(calcFunctionT(1:5).eq."POLY2")Then   ! Second order (3 params)
      parameters3 = PolyFit(dataPointsIn,2)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = CalcPolynomial(parameters3,dataPointsOut(i,1))
      End Do
      fittingReport(1) = "2nd Order Polynomial"
      fittingReport(2) = "P(x)="&
                         //trim(DpToStr(parameters3(1)))//"+"&
                         //trim(DpToStr(parameters3(2)))//"x+"&
                         //trim(DpToStr(parameters3(3)))//"x^2"
      fittingReport(3) = "--------------------------------------------------------------"
    End If
    If(calcFunctionT(1:5).eq."POLY3")Then   ! Third order (4 params)
      parameters4 = PolyFit(dataPointsIn,3)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = CalcPolynomial(parameters4,dataPointsOut(i,1))
      End Do
      fittingReport(1) = "3rd Order Polynomial"
      fittingReport(2) = "P(x)="&
                         //trim(DpToStr(parameters4(1)))//"+"&
                         //trim(DpToStr(parameters4(2)))//"x+"&
                         //trim(DpToStr(parameters4(3)))//"x^2"&
                         //trim(DpToStr(parameters4(4)))//"x^3"
      fittingReport(3) = "--------------------------------------------------------------"
    End If
    If(calcFunctionT(1:5).eq."POLY4")Then   ! Fourth order (5 params)
      parameters5 = PolyFit(dataPointsIn,4)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = CalcPolynomial(parameters5,dataPointsOut(i,1))
      End Do
      fittingReport(1) = "4th Order Polynomial"
      fittingReport(2) = "P(x)=("&
                         //trim(DpToStr(parameters5(1)))//")+("&
                         //trim(DpToStr(parameters5(2)))//")x+("&
                         //trim(DpToStr(parameters5(3)))//")x^2+("&
                         //trim(DpToStr(parameters5(4)))//")x^3+("&
                         //trim(DpToStr(parameters5(5)))//")x^4"
      fittingReport(3) = "--------------------------------------------------------------"
    End If
    If(calcFunctionT(1:5).eq."POLY5")Then   ! Fifth order (6 params)
      parameters6 = PolyFit(dataPointsIn,5)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = CalcPolynomial(parameters6,dataPointsOut(i,1))
      End Do
      fittingReport(1) = "5th Order Polynomial"
      fittingReport(2) = "P(x)=("&
                         //trim(DpToStr(parameters6(1)))//")+("&
                         //trim(DpToStr(parameters6(2)))//")x+("&
                         //trim(DpToStr(parameters6(3)))//")x^2+("&
                         //trim(DpToStr(parameters6(4)))//")x^3+("&
                         //trim(DpToStr(parameters6(5)))//")x^4+("&
                         //trim(DpToStr(parameters6(6)))//")x^5"
      fittingReport(3) = "--------------------------------------------------------------"
    End If
! EXP fit
    If(calcFunctionT(1:7).eq."EXPFIT1")Then ! 1 term
      parameters2 = ExpFit(dataPointsIn,1)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = ExpCalc(dataPointsOut(i,1),parameters2)
      End Do
      fittingReport(1) = "Exponential Fit, 1 term"
      fittingReport(2) = "P(x)=("&
                         //trim(DpToStr(parameters2(1)))//")exp[("&
                         //trim(DpToStr(parameters2(2)))//")x]"
      fittingReport(3) = "--------------------------------------------------------------"
    End If
    If(calcFunctionT(1:7).eq."EXPFIT2")Then ! 2 terms
      parameters4 = ExpFit(dataPointsIn,2)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = ExpCalc(dataPointsOut(i,1),parameters4)
      End Do
      fittingReport(1) = "Exponential Fit, 2 terms"
      fittingReport(2) = "P(x)=("&
                         //trim(DpToStr(parameters4(1)))//")exp[("&
                         //trim(DpToStr(parameters4(2)))//")x]+("&
                         //trim(DpToStr(parameters4(3)))//")exp[("&
                         //trim(DpToStr(parameters4(4)))//")x]"
      fittingReport(3) = "--------------------------------------------------------------"
    End If
    If(calcFunctionT(1:7).eq."EXPFIT3")Then ! 3 terms
      parameters6 = ExpFit(dataPointsIn,3)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = ExpCalc(dataPointsOut(i,1),parameters6)
      End Do
      fittingReport(1) = "Exponential Fit, 3 terms"
      fittingReport(2) = "P(x)=("&
                         //trim(DpToStr(parameters6(1)))//")exp[("&
                         //trim(DpToStr(parameters6(2)))//")x]+("&
                         //trim(DpToStr(parameters6(3)))//")exp[("&
                         //trim(DpToStr(parameters6(4)))//")x]+("&
                         //trim(DpToStr(parameters6(5)))//")exp[("&
                         //trim(DpToStr(parameters6(6)))//")x]"
      fittingReport(3) = "--------------------------------------------------------------"
    End If
! Bulk Modulus Fit
    If(calcFunctionT(1:3).eq."BM1")Then   ! BM Restrict BP0
      parameters4 = BirchMurnFit(dataPointsIn,1.0D0,9.0D0)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = BirchMurnCalc(dataPointsOut(i,1),parameters4)
      End Do
    End If
    If(calcFunctionT(1:3).eq."BM2")Then   ! BM
      parameters4 = BirchMurnFit(dataPointsIn)
      Do i=1,pointsOutCount
        dataPointsOut(i,1) = xStart+(i-1)*xInc
        dataPointsOut(i,2) = BirchMurnCalc(dataPointsOut(i,1),parameters4)
      End Do
    End If
! Spline Fit
    If(calcFunctionT(1:7).eq."SPLINE1")Then
      dataPointsOut = SplinePoints(dataPointsIn, pointsOutCount,1)
    End If
    If(calcFunctionT(1:7).eq."SPLINE3")Then
      dataPointsOut = SplinePoints(dataPointsIn, pointsOutCount,3)
    End If
    If(calcFunctionT(1:7).eq."SPLINE ".or.calcFunctionT(1:7).eq."SPLINE5")Then
      dataPointsOut = SplinePoints(dataPointsIn, pointsOutCount,5)
    End If
! Interp Fit
    If(calcFunctionT(1:7).eq."INTERP3")Then
      dataPointsOut = InterpPoints(dataPointsIn, pointsOutCount,3)
    End If
    If(calcFunctionT(1:7).eq."INTERP ".or.calcFunctionT(1:7).eq."INTERP4")Then
      dataPointsOut = InterpPoints(dataPointsIn, pointsOutCount,4)
    End If
    If(calcFunctionT(1:7).eq."INTERP5")Then
      dataPointsOut = InterpPoints(dataPointsIn, pointsOutCount,5)
    End If
    If(calcFunctionT(1:5).eq."EXACT")Then
      dataPointsOut = FullInterpPoints(dataPointsIn, pointsOutCount)
    End If
  End Function FittingPoints

! ----------------------------------------------------------------------------------
! RSS Functions
!
!
!
!
! ----------------------------------------------------------------------------------

  Function DoubleDecayFitRSS(dataPoints, a, b, lA, lB) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, lA, lB, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)+b*exp(lB*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
  End Function DoubleDecayFitRSS

  Function TripleDecayFitRSS(dataPoints, a, lA, b, lB, c, lC) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, c, lA, lB, lC, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)+b*exp(lB*x)+c*exp(lC*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
  End Function TripleDecayFitRSS

! ----------------------------------------------------------------------------------
! Spline Functions
!
!
!
!
! ----------------------------------------------------------------------------------

Function SplineAB(pointA, pointB, splineTypeIn) RESULT (coefficients)
! Polynomial to spline between points A and B - x, f(x), f'(x) and f''(x) supplied
  Implicit None  !Force declaration of all variables
! In:      Declare variables
  Real(kind=DoubleReal), Dimension(:) :: pointA  !1 x, 2 f(x), 3 f'(x), 4 f''(x)....
  Real(kind=DoubleReal), Dimension(:) :: pointB  !1 x, 2 f(x), 3 f'(x), 4 f''(x)....
  Integer(kind=StandardInteger), Optional :: splineTypeIn  ! 1 = P(x), 2 = exp
! Out:     Declare variables
  Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1))) :: coefficients
! Private: Declare variables
  Integer(kind=StandardInteger) :: splineType
  Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1))) :: yMatrix
  Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1)),1:(2*(size(pointA,1)-1))) :: xMatrix
  Integer(kind=StandardInteger) :: i, j, n, row, col, matrixSize, matrixHalfSize, expt
  Real(kind=DoubleReal) :: x, coeff
! Optional arguments
  splineType = 1
  If(Present(splineTypeIn))Then
    splineType = splineTypeIn
  End If
! Init variables
  coefficients = 0.0D0
  yMatrix = 0.0D0
  xMatrix = 0.0D0
  coeff = 0.0D0
  matrixHalfSize = (size(pointA,1)-1)
  matrixSize = 2*matrixHalfSize
! Make x-matrix
  row = 0
  Do i=1,matrixHalfSize
    row = row + 1
    n = 0
    Do col=1,matrixSize
      expt = col - i
      If(expt.lt.0)Then
        x = 0.0D0
        coeff = 0.0D0
      Else
        n = n + 1
        If(row.eq.1)Then
          coeff = 1.0D0
        Else
          coeff = 1.0D0
          Do j=1,(row-1)
            coeff = coeff * (n+j-1)
          End Do
        End If
      End If
      If(col.ge.row)Then
        xMatrix(2*row-1,col) = coeff*pointA(1)**(expt)
        xMatrix(2*row,col) = coeff*pointB(1)**(expt)
      Else
        xMatrix(2*row-1,col) = 0.0D0
        xMatrix(2*row,col) = 0.0D0
      End If
    End Do
  End Do
! make y-matrix
  row = 0
  Do i=1,matrixHalfSize
    row = row + 1
    yMatrix(row) = pointA(i+1)
    row = row + 1
    yMatrix(row) = pointB(i+1)
  End Do
! Spline Type 1
  If(splineType.eq.1)Then
! solve equation
    coefficients = SolveLinearSet(xMatrix,yMatrix)
  End If
! Spline Type 1
  If(splineType.eq.2)Then
! solve equation
    Call PositiveYMatrix(xMatrix,yMatrix)
    Call LnYMatrix(yMatrix)
    coefficients = SolveLinearSet(xMatrix,yMatrix)
  End If
End Function SplineAB


Function SplineExpThird(xA,fxA,fpxA,xB,fxB,fpxB) RESULT (coefficients)
! Third order polynomial
! Find parameters for spline between points A and B
! form of spline function:
!   f(x) = exp(a+bx+cx^2+dx^3)
!   f'(x) = (b+2cx+3dx^2)*exp(a+bx+cx^2+dx^3) = (b+2cx+3dx^2)*f(x)
! xA = Point A x value, fxA = point A f(x), fpxA = Point A f'(x)
! xB = Point B x value, fxB = point B f(x), fpxB = Point B f'(x)
  Implicit None   ! Force declaration of all variables
! In
  Real(kind=DoubleReal) :: xA,fxA,fpxA,xB,fxB,fpxB
! Out
  Real(kind=DoubleReal), Dimension(1:4) :: coefficients
! Private variables
  Real(kind=DoubleReal), Dimension(1:4,1:4) :: xMatrix
  Real(kind=DoubleReal), Dimension(1:4) :: yMatrix
! init matrices
  xMatrix = 0.0D0
  yMatrix = 0.0D0
! make x matrix
  xMatrix(1,1) = 1.0D0
  xMatrix(1,2) = xA
  xMatrix(1,3) = xA**2
  xMatrix(1,4) = xA**3
  xMatrix(2,1) = 1.0D0
  xMatrix(2,2) = xB
  xMatrix(2,3) = xB**2
  xMatrix(2,4) = xB**3
  xMatrix(3,1) = 0.0D0
  xMatrix(3,2) = 1.0D0
  xMatrix(3,3) = 2.0D0*xA
  xMatrix(3,4) = 3.0D0*xA**2
  xMatrix(4,1) = 0.0D0
  xMatrix(4,2) = 1.0D0
  xMatrix(4,3) = 2.0D0*xB
  xMatrix(4,4) = 3.0D0*xB**2
! make y matrix
  yMatrix(1) = log(fxA)
  yMatrix(2) = log(fxB)
  yMatrix(3) = fpxA/fxA
  yMatrix(4) = fpxB/fxB
! solve
  coefficients = SolveLinearSet(xMatrix,yMatrix) ! By LU decomposition
End Function SplineExpThird

Function SplineExpFifth(xA,fxA,fpxA,fppxA,xB,fxB,fpxB,fppxB) RESULT (coefficients)
! Fifth order polynomial
! Find parameters for spline between points A and B
! form of spline function:
!   f(x) = exp(a+bx+cx^2+dx^3+ex^4+fx^5)
!   f'(x) = 0a+(b+2cx+3dx^2+4ex^3+5fx^4)*f(x)
!   f''(x) = a*0+b*y'+c*2(y+xy')+d*3x(2y+xy')+e*4x^2*(3y+xy')+f*5x^3(4y+xy')
! xA = Point A x value, fxA = point A f(x), fpxA = Point A f'(x)
! xB = Point B x value, fxB = point B f(x), fpxB = Point B f'(x)
  Implicit None   ! Force declaration of all variables
! In
  Real(kind=DoubleReal) :: xA,fxA,fpxA,fppxA,xB,fxB,fpxB,fppxB
! Out
  Real(kind=DoubleReal), Dimension(1:6) :: coefficients
! Private variables
  Real(kind=DoubleReal), Dimension(1:6,1:6) :: xMatrix
  Real(kind=DoubleReal), Dimension(1:6) :: yMatrix
! init matrices
  xMatrix = 0.0D0
  yMatrix = 0.0D0
! make x matrix
  xMatrix(1,1) = 1.0D0
  xMatrix(1,2) = xA
  xMatrix(1,3) = xA**2
  xMatrix(1,4) = xA**3
  xMatrix(1,5) = xA**4
  xMatrix(1,6) = xA**5
  xMatrix(2,1) = 1.0D0
  xMatrix(2,2) = xB
  xMatrix(2,3) = xB**2
  xMatrix(2,4) = xB**3
  xMatrix(2,5) = xB**4
  xMatrix(2,6) = xB**5
  xMatrix(3,1) = 0.0D0
  xMatrix(3,2) = 1.0D0
  xMatrix(3,3) = 2.0D0*xA
  xMatrix(3,4) = 3.0D0*xA**2
  xMatrix(3,5) = 4.0D0*xA**3
  xMatrix(3,6) = 5.0D0*xA**4
  xMatrix(4,1) = 0.0D0
  xMatrix(4,2) = 1.0D0
  xMatrix(4,3) = 2.0D0*xB
  xMatrix(4,4) = 3.0D0*xB**2
  xMatrix(4,5) = 4.0D0*xB**3
  xMatrix(4,6) = 5.0D0*xB**4
  xMatrix(5,1) = 0.0D0
  xMatrix(5,2) = 1.0D0*fpxA
  xMatrix(5,3) = 2.0D0*(fxA+xA*fpxA)
  xMatrix(5,4) = 3.0D0*xA*(2.0D0*fxA+xA*fpxA)
  xMatrix(5,5) = 4.0D0*(xA**2)*(3.0D0*fxA+xA*fpxA)
  xMatrix(5,6) = 5.0D0*(xA**3)*(4.0D0*fxA+xA*fpxA)
  xMatrix(6,1) = 0.0D0
  xMatrix(6,2) = 1.0D0*fpxB
  xMatrix(6,3) = 2.0D0*(fxB+xB*fpxB)
  xMatrix(6,4) = 3.0D0*xB*(2.0D0*fxB+xB*fpxB)
  xMatrix(6,5) = 4.0D0*(xB**2)*(3.0D0*fxB+xB*fpxB)
  xMatrix(6,6) = 5.0D0*(xB**3)*(4.0D0*fxB+xB*fpxB)
! make y matrix
  yMatrix(1) = log(fxA)
  yMatrix(2) = log(fxB)
  yMatrix(3) = fpxA/fxA
  yMatrix(4) = fpxB/fxB
  yMatrix(5) = fppxA
  yMatrix(6) = fppxB
! solve
  coefficients = SolveLinearSet(xMatrix,yMatrix) ! By LU decomposition
End Function SplineExpFifth

Function SplineNodes(inputNodes,numDataPoints,startPoint,endPoint,&
dataSize,splineTypeIn,forceCalcDervIn,interpNodeIn) RESULT (dataPoints)
! Input nodes x,f(x) for each node, calculate f'(x) and f''(x) from the set of nodes, then spline
! inputNodes         array of nodes
! numDataPoints      total data points to output
! startPoint         starting node
! endPoint           ending node
! dataSize           size of data points array (ge numDataPoints)
! splineTypeIn       array of type of spline to use for each segment
! Variable length output
  Implicit None  !Force declaration of all variables
! Declare variables - arg
  Real(kind=DoubleReal), Dimension(:,:) :: inputNodes
  Integer(kind=StandardInteger) :: numDataPoints, startPoint, endPoint, dataSize
  Integer(kind=StandardInteger), Dimension(1:1000), Optional :: splineTypeIn
  Logical, Dimension(1:1000), Optional :: forceCalcDervIn, interpNodeIn
! Declare variables - priv
  Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1),1:4) :: splineNodeArr
  Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1),1:2) :: interpeNodesArr
  Integer(kind=StandardInteger) :: i, j, nodeKey, nodeCount, nodeCountI, interpKey
  Real(kind=DoubleReal), Dimension(1:dataSize,1:4) :: dataPoints ! output data points
  Real(kind=DoubleReal) :: x, xStart, xEnd, xIncrement
  Real(kind=DoubleReal), Dimension(1:6) :: coefficients
  Real(kind=DoubleReal), Dimension(1:4) :: expThird
  Real(kind=DoubleReal), Dimension(1:6) :: expFifth, polyFitCoeffs
  Real(kind=DoubleReal), Dimension(1:4) :: embFuncC, densFunc
  Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
  Real(kind=DoubleReal), Dimension(1:3) :: pointAA, pointBB
  Real(kind=DoubleReal), Dimension(1:3) :: yArray, embFuncB
  Real(kind=DoubleReal), Dimension(1:2) :: embFuncA
  Integer(kind=StandardInteger), Dimension(1:1000) :: splineType, nodeMap
  !Integer(kind=StandardInteger), Dimension(1:1000,1:2) :: nodeMap
  Logical, Dimension(1:1000) :: forceCalcDerv, interpNode
  Real(kind=DoubleReal), Dimension(1:4,1:2) :: interpPoints
! optional arguments
  splineType = 1   ! set to standard 5th order polynomial  1 a+bx+...  2 exp(a+bx+...)
  If(present(splineTypeIn))Then
    splineType = splineTypeIn
  End If
  forceCalcDerv = .false.
  If(present(forceCalcDervIn))Then  ! Interp between nodes to calc f'(x) and f''(x)
    forceCalcDerv = forceCalcDervIn
  End If
  interpNode = .true.
  If(present(interpNodeIn))Then  ! Interp between nodes to calc f'(x) and f''(x)
    interpNode = interpNodeIn
  End If
! Init Variables
  splineCoeffs = 0.0D0
  nodeMap = 0
  interpPoints = 0.0D0
  dataPoints = 0.0D0
  If(startPoint.eq.0)Then
    startPoint = 1
  End If
  If(endPoint.eq.0)Then
    endPoint = size(inputNodes,1)
  End If
! Make array of points used for interpolation
  j = 0
  nodeKey = 0
  nodeCountI = 0
  Do i=startPoint,endPoint
    nodeKey = nodeKey + 1
    If(interpNode(nodeKey))Then
      j = j + 1
      interpeNodesArr(j,1) = inputNodes(i,1)
      interpeNodesArr(j,2) = inputNodes(i,2)
      nodeMap(nodeKey) = j
    End If
  End Do
  nodeCountI = j
! Transfer nodes from inputNodes to splineNodeArr and calculate 1st/2nd derivatives
  nodeKey = 0
  nodeCount = 0
  Do i=startPoint,endPoint
    nodeKey = nodeKey + 1
! don't change x/y values
    splineNodeArr(nodeKey,1) = inputNodes(i,1)
    splineNodeArr(nodeKey,2) = inputNodes(i,2)
! If calculate 1st/2nd deriv
    If(forceCalcDerv(nodeKey))Then
      x = inputNodes(i,1)
      interpKey = nodeMap(nodeKey) - 1
      If(interpKey.lt.1)Then
        interpKey = 1
      End If
      If(interpKey.gt.(nodeCountI-3))Then
        interpKey = (nodeCountI-3)
      End If
      Do j=1,4
        interpPoints(j,1) = interpeNodesArr(interpKey+j-1,1)
        interpPoints(j,2) = interpeNodesArr(interpKey+j-1,2)
      End Do
      splineNodeArr(nodeKey,3) = InterpLagrange(x,interpPoints,1)
      splineNodeArr(nodeKey,4) = InterpLagrange(x,interpPoints,2)
    Else
! Else store input value
      splineNodeArr(nodeKey,3) = inputNodes(i,3)
      splineNodeArr(nodeKey,4) = inputNodes(i,4)
    End If
  End Do
  nodeCount = nodeKey
! Calculate spline data points
  xStart = splineNodeArr(1,1)
  xEnd = splineNodeArr(nodeCount,1)
  xIncrement = (xEnd-xStart)/(1.0D0*numDataPoints-1.0D0)
! Loop through data points
  nodeKey = 0
  x = xStart
  Do i=1,numDataPoints  ! loop through data points being output
    If((i.eq.1).or.(x.ge.splineNodeArr(nodeKey+1,1).and.(nodeKey+1).lt.nodeCount))Then
      nodeKey = nodeKey + 1
      pointA(1) = splineNodeArr(nodeKey,1)
      pointA(2) = splineNodeArr(nodeKey,2)
      pointA(3) = splineNodeArr(nodeKey,3)
      pointA(4) = splineNodeArr(nodeKey,4)
      pointB(1) = splineNodeArr(nodeKey+1,1)
      pointB(2) = splineNodeArr(nodeKey+1,2)
      pointB(3) = splineNodeArr(nodeKey+1,3)
      pointB(4) = splineNodeArr(nodeKey+1,4)
      If(splineType(nodeKey).eq.1)Then           ! Normal 5th order spline
        coefficients = SplineAB(pointA, pointB)
        Do j=1,6
          splineCoeffs(nodeKey,j) = coefficients(j)
        End Do
      End If
      If(splineType(nodeKey).eq.2)Then           ! exp(3rd order)
        pointAA(1) = pointA(1)
        pointAA(2) = pointA(2)
        pointAA(3) = pointA(3)
        pointBB(1) = pointB(1)
        pointBB(2) = pointB(2)
        pointBB(3) = pointB(3)
        expThird = SplineAB(pointAA, pointBB, 2)
        Do j=1,4
          splineCoeffs(nodeKey,j) = expThird(j)
        End Do
      End If
      If(splineType(nodeKey).eq.3)Then           ! exp(5th order)
        expFifth = SplineExpFifth(pointA(1),pointA(2),pointA(3),pointA(4),&
        pointB(1),pointB(2),pointB(3),pointB(4))
        Do j=1,6
          splineCoeffs(nodeKey,j) = expFifth(j)
        End Do
      End If
      If(splineType(nodeKey).eq.4.and.nodeKey.eq.1)Then           ! embedding function F(p) = a+bp^0.5+bp^2+dp^4
        embFuncC = FitEmbeddingC(splineNodeArr,1,nodeCount)
      End If
      If(splineType(nodeKey).eq.5.and.nodeKey.eq.1)Then           ! density function p(r) = a r^2 exp(b r^2) + c r^2 exp(d r^2)
        densFunc = FitDensity(splineNodeArr,1,nodeCount)
      End If
      If(splineType(nodeKey).eq.6.and.nodeKey.eq.1)Then           ! density function p(r) = 5th order poly
        polyFitCoeffs = PolyFit(splineNodeArr,5)
      End If
      If(splineType(nodeKey).eq.7.and.nodeKey.eq.1)Then           ! embedding function F(p) = a+bp^0.5+bp^2
        embFuncB = FitEmbeddingB(splineNodeArr,1,nodeCount)
      End If
      If(splineType(nodeKey).eq.8.and.nodeKey.eq.1)Then           ! embedding function F(p) = a+bp^0.5
        embFuncA = FitEmbeddingA(splineNodeArr,1,nodeCount)
      End If
    End If
    If(splineType(nodeKey).eq.1)Then
      dataPoints(i,1) = x
      dataPoints(i,2) = CalcPolynomial(coefficients, x, 0)
      dataPoints(i,3) = CalcPolynomial(coefficients, x, 1)
      dataPoints(i,4) = CalcPolynomial(coefficients, x, 2)
    End If
    If(splineType(nodeKey).eq.2)Then
      dataPoints(i,1) = x
      dataPoints(i,2) = CalcPolynomialExp(expThird, x, 0)
      dataPoints(i,3) = CalcPolynomialExp(expThird, x, 1)
      dataPoints(i,4) = CalcPolynomialExp(expThird, x, 2)
    End If
    If(splineType(nodeKey).eq.3)Then
      dataPoints(i,1) = x
      dataPoints(i,2) = CalcPolynomialExp(expFifth, x, 0)
      dataPoints(i,3) = CalcPolynomialExp(expFifth, x, 1)
      dataPoints(i,4) = CalcPolynomialExp(expFifth, x, 2)
    End If
    If(splineType(nodeKey).eq.4)Then
      dataPoints(i,1) = x
      yArray = EmbeddingC(x,embFuncC)
      dataPoints(i,2) = yArray(1)
      dataPoints(i,3) = yArray(2)
      dataPoints(i,4) = yArray(3)
    End If
    If(splineType(nodeKey).eq.5)Then
      dataPoints(i,1) = x
      dataPoints(i,2) = densFunc(1)*x**2*exp(densFunc(2)*x**2)+densFunc(3)*x**2*exp(densFunc(4)*x**2)
      dataPoints(i,3) = densFunc(1)*x**2*exp(densFunc(2)*x**2)+densFunc(3)*x**2*exp(densFunc(4)*x**2)
      dataPoints(i,4) = densFunc(1)*x**2*exp(densFunc(2)*x**2)+densFunc(3)*x**2*exp(densFunc(4)*x**2)
    End If
    If(splineType(nodeKey).eq.6)Then
      dataPoints(i,1) = x
      dataPoints(i,2) = CalcPolynomialExp(polyFitCoeffs, x, 0)
      dataPoints(i,3) = CalcPolynomialExp(polyFitCoeffs, x, 1)
      dataPoints(i,4) = CalcPolynomialExp(polyFitCoeffs, x, 2)
    End If
    If(splineType(nodeKey).eq.7)Then
      dataPoints(i,1) = x
      yArray = EmbeddingB(x,embFuncB)
      dataPoints(i,2) = yArray(1)
      dataPoints(i,3) = yArray(2)
      dataPoints(i,4) = yArray(3)
    End If
    If(splineType(nodeKey).eq.8)Then
      dataPoints(i,1) = x
      yArray = EmbeddingA(x,embFuncA)
      dataPoints(i,2) = yArray(1)
      dataPoints(i,3) = yArray(2)
      dataPoints(i,4) = yArray(3)
    End If
! Increment x
    x = x + xIncrement
  End Do
End Function SplineNodes

Function SplineComplete(inputPoints,interpSizeIn) RESULT (splinePoints)
! Complete missing points f'(x) f''(x) by interpolation
  Implicit None  !Force declaration of all variables
! Declare variables
  Real(kind=DoubleReal), Dimension(:,:) :: inputPoints
  Integer(kind=StandardInteger) :: interpSizeIn
  Real(kind=DoubleReal), Dimension(1:size(inputPoints,1),1:4) :: splinePoints
  Real(kind=DoubleReal), Dimension(1:3) :: yArray
  Real(kind=DoubleReal) :: x
  Integer(kind=StandardInteger) :: i
! Init variables
  splinePoints = 0.0D0
  Do i=1,size(inputPoints,1)
    x = inputPoints(i,1)
    yArray = PointInterp(inputPoints,x,interpSizeIn,2)
    splinePoints(i,1) = x
    splinePoints(i,2) = yArray(1)  ! f(x)
    splinePoints(i,3) = yArray(2)  ! f'(x)
    splinePoints(i,4) = yArray(3)  ! f''(x)
  End Do
End Function SplineComplete



! ----- Might not be useful

Function VaryNode(nodeValue, varyAmount) RESULT (outputValue)
! Used by VaryNode
  Implicit None   ! Force declaration of all variables
! Private variables
  Real(kind=DoubleReal) :: randDouble
  Real(kind=DoubleReal) :: nodeValue, varyAmount, outputValue
! Get rand number
  Call RANDOM_NUMBER(randDouble)
  outputValue = nodeValue + varyAmount*(randDouble-0.5D0)
End Function VaryNode

Function FillSplineResponse(dataPointsIn, startIn, endIn) RESULT (dataPointsOut)
! Fill in gaps where there was no response to adjusting the parameter
  Implicit None   ! Force declaration of all variables
! Private variables
  Real(kind=DoubleReal), Dimension(:,:) :: dataPointsIn
  Real(kind=DoubleReal), Dimension(1:size(dataPointsIn,1),1:size(dataPointsIn,2)) :: dataPointsOut
  Integer(kind=StandardInteger) :: i, j, k, startI, endI
  Real(kind=DoubleReal) :: x, y, xA, xB, yA, yB, grad
  Integer(kind=StandardInteger), optional :: startIn, endIn
! Init
  startI = 1
  endI = size(dataPointsIn,1)
  dataPointsOut = dataPointsIn
! optional
  If(Present(startIn))Then
    startI = startIn
  End If
  If(Present(endIn))Then
    endI = endIn
  End If
! Loop
  Do i=startI+1,endI-1
    If(dataPointsOut(i,2).eq.0.0D0)Then
      xA = dataPointsOut(i-1,1)
      yA = dataPointsOut(i-1,2)
! find next point that isn't 0
      k = 0
      Do j = i+1,endI-1
        k = k + 1
        If(dataPointsOut(j,2).ne.0.0D0)Then
          xB = dataPointsOut(j,1)
          yB = dataPointsOut(j,2)
          Exit
        End If
      End Do
      grad = (yB-yA)/(xB-xA)
      Do j=i,i+k
        x = dataPointsOut(j,1)
        y = yA+(x-xA)*grad
        dataPointsOut(j,2) = y
      End Do
    End If
  End Do
End Function FillSplineResponse

! ---------------------------------------------
! Spline Fitting
! ---------------------------------------------

Function SplinePoints(dataPointsIn, pointsOutCount,splineOrderIn) RESULT (dataPointsOut)
! Force declaration of all variables
  Implicit None
! In:      Declare variables
  Real(kind=DoubleReal), Dimension(:,:) :: dataPointsIn
  Integer(kind=StandardInteger) :: pointsOutCount
  Integer(kind=StandardInteger), Optional :: splineOrderIn
! Out:     Declare variables
  Real(kind=DoubleReal), Dimension(1:pointsOutCount,1:2) :: dataPointsOut
! Private: Declare variables
  Integer(kind=StandardInteger) :: i, j, n, k, pointsInCount, splineOrder
  Real(kind=DoubleReal) :: x
  Real(kind=DoubleReal) :: xStart,xEnd,xInc
  Real(kind=DoubleReal) :: xUpper,xLower
  Real(kind=DoubleReal) :: xA, xB
  Real(kind=DoubleReal), Dimension(1:4,1:2) :: interpArray
  Real(kind=DoubleReal), Dimension(1:2) :: A_A, B_A
  Real(kind=DoubleReal), Dimension(1:3) :: A_B, B_B
  Real(kind=DoubleReal), Dimension(1:4) :: A_C, B_C
  Real(kind=DoubleReal), Dimension(1:2) :: coefficients_A
  Real(kind=DoubleReal), Dimension(1:4) :: coefficients_B
  Real(kind=DoubleReal), Dimension(1:6) :: coefficients_C
! Optional
  splineOrder = 5
  If(Present(splineOrderIn))Then
    splineOrder = splineOrderIn
  End If
! Init
  pointsInCount = size(dataPointsIn,1)
  xStart = dataPointsIn(1,1)
  xEnd = dataPointsIn(pointsInCount,1)
  xInc = (xEnd-xStart)/(pointsOutCount-1.0D0)
! Loop through points to make
  n = 1
  Do i=1,pointsOutCount
! Output node x val
    x = xStart+(i-1)*xInc
! Input node x lower/upper
    xLower = dataPointsIn(n,1)
    xUpper = dataPointsIn(n+1,1)
! Get spline coefficients
    If(i.eq.1.or.x.gt.xUpper)Then
! Find start and end node
      Do k=n,(pointsInCount-1)
        xLower = dataPointsIn(k,1)
        xUpper = dataPointsIn(k+1,1)
        If(x.ge.xLower.and.x.le.xUpper)Then
          n = k
        End If
      End Do
      xA = dataPointsIn(n,1)
      xB = dataPointsIn(n+1,1)
! Interp node start (four point interp 1,(2),3,4)
      k = n - 1
      If(k.lt.1)Then
        k = 1
      End If
      If((k+3).gt.pointsInCount)Then
        k = pointsInCount-3
      End If
      interpArray = 0.0D0
      Do j=1,4
        interpArray(j,1) = dataPointsIn(k+(j-1),1)
        interpArray(j,2) = dataPointsIn(k+(j-1),2)
      End Do
! 1st order spline
      If(splineOrder.eq.1)Then
        A_A(1) = xA
        A_A(2) = InterpLagrange(xA,interpArray,0)
        B_A(1) = xB
        B_A(2) = InterpLagrange(xB,interpArray,0)
        coefficients_A = SplineAB(A_A,B_A)
      End If
! 3rd order spline
      If(splineOrder.eq.3)Then
        A_B(1) = xA
        A_B(2) = InterpLagrange(xA,interpArray,0)
        A_B(3) = InterpLagrange(xA,interpArray,1)
        B_B(1) = xB
        B_B(2) = InterpLagrange(xB,interpArray,0)
        B_B(3) = InterpLagrange(xB,interpArray,1)
        coefficients_B = SplineAB(A_B,B_B)
      End If
! 5th order spline
      If(splineOrder.eq.5)Then
        A_C(1) = xA
        A_C(2) = InterpLagrange(xA,interpArray,0)
        A_C(3) = InterpLagrange(xA,interpArray,1)
        A_C(4) = InterpLagrange(xA,interpArray,2)
        B_C(1) = xB
        B_C(2) = InterpLagrange(xB,interpArray,0)
        B_C(3) = InterpLagrange(xB,interpArray,1)
        B_C(4) = InterpLagrange(xB,interpArray,2)
        coefficients_C = SplineAB(A_C,B_C)
      End If
    End If
! store output points
    dataPointsOut(i,1) = x
    If(splineOrder.eq.1)Then
      dataPointsOut(i,2) = CalcPolynomial(coefficients_A,x)
    End If
    If(splineOrder.eq.3)Then
      dataPointsOut(i,2) = CalcPolynomial(coefficients_B,x)
    End If
    If(splineOrder.eq.5)Then
      dataPointsOut(i,2) = CalcPolynomial(coefficients_C,x)
    End If
  End Do
End Function SplinePoints



! ----------------------------------------------------------------------------------
! Module Subroutines
!
!
!
!
! ----------------------------------------------------------------------------------

Subroutine CompleteNodeData(splineNodes, startIn, endIn)
! Complete y'(x) and y''(x) values for a y(x) spline
! At least 4 data points required, otherwise exits
! At least 4 "columns" in splineNodes x, y(x), y'(x). y''(x)
  Implicit None  ! Force declaration of all variables
! Declare private variables
  Real(kind=DoubleReal), Dimension(:,:) :: splineNodes
  Integer(kind=StandardInteger) :: node, nodes, i, j
  Integer(kind=StandardInteger), optional :: startIn, endIn
  Integer(kind=StandardInteger) :: startNode, endNode, tempNode
  Integer(kind=StandardInteger) :: interpStart, interpEnd
  Real(kind=DoubleReal), Dimension(1:4,1:2) :: interpNodes
! optional arguments
  startNode = 1
  endNode = size(splineNodes,1)
  If(present(startIn))Then
    If(startIn.ge.1)Then
      startNode = startIn
    End If
  End If
  If(present(endIn))Then
    If(endIn.le.endNode)Then
      endNode = endIn
    End If
  End If
! Swap start/end if wrong way around
  If(startNode.gt.endNode)Then
    tempNode = startNode
    startNode = endNode
    endNode = tempNode
  End If
! Exit subroutine if too few points
  nodes = endNode-startNode+1
  If(nodes.lt.4)Then
    return
  End If
! Interp at each node
  Do node=startNode,endNode
! Set start-end nodes used for interpolation
    interpStart = node-2
    interpEnd = interpStart+3
! Adjust these nodes so they are within the bounds of the data
    If(interpStart.lt.startNode)Then
      interpStart = startNode
      interpEnd = interpStart + 3
    End If
    If(interpEnd.gt.endNode)Then
      interpEnd = endNode
      interpStart = interpEnd - 3
    End If
! make interpolation array
    j = 0
    Do i=interpStart,interpEnd
      j = j + 1
      interpNodes(j,1) = splineNodes(i,1)
      interpNodes(j,2) = splineNodes(i,2)
    End Do
! print *,node,splineNodes(node,1),splineNodes(node,2),splineNodes(node,3),splineNodes(node,4)
! lagrange interpolation
    splineNodes(node,2) = InterpLagrange(splineNodes(node,1), interpNodes, 0)
    splineNodes(node,3) = InterpLagrange(splineNodes(node,1), interpNodes, 1)
    splineNodes(node,4) = InterpLagrange(splineNodes(node,1), interpNodes, 2)
! print *,node,splineNodes(node,1),splineNodes(node,2),splineNodes(node,3),splineNodes(node,4)
! print *,""
  End Do
End Subroutine CompleteNodeData




End Module splinesFitting
