Module fitting
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use rng
  Use regression
  Use interpolation
  Use lmaM
! Force declaration of all variables
  Implicit None
! Public variables  
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: BirchMurnFit
  Public :: SingleDecayFit
  Public :: DoubleDecayFit
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  

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
    
  Function SingleDecayFit(dataPoints) RESULT (parameters)
! Fit Murnaghan EoS to data
! Fitting method adapted from http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/appendices/appendix-eos.html
! Birch-Murnaghan equation described by Hebbache 2004
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
! Out
    Real(kind=DoubleReal), Dimension(1:2) :: parameters
! Private
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
    Do i=1,2
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
            linCoeffs = LinearRegression(yReg, xReg)
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
  
  
End Module fitting