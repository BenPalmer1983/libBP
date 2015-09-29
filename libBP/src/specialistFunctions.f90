Module specialistFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use matrix
  Use regression
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: FitEmbedding
  Public :: FitDensity
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  
  Function FitEmbedding(dataPoints, startPointIn, endPointIn) Result (coefficients)
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
    pointCount = endPoint-startPoint+1
! run subroutine
    Call FitEmbedding_Process(dataPoints, pointCount, startPoint, endPoint, coefficients)
  End Function FitEmbedding
! ----------
  Subroutine FitEmbedding_Process(dataPoints, pointCount, startPoint, endPoint, coefficients)
    Implicit None  !Force declaration of all variables
! Declare variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Integer(kind=StandardInteger) :: pointCount, startPoint, endPoint
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
! Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal), Dimension(1:pointCount,1:4) :: X
    Real(kind=DoubleReal), Dimension(1:4,1:pointCount) :: XT
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: XTX
    Real(kind=DoubleReal), Dimension(1:pointCount) :: Y
    Real(kind=DoubleReal), Dimension(1:4) :: XTY
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
    XT = TransposeMatrix(X)
    XTX = matmul(XT,X)
    XTX = InvertMatrix(XTX)
    XTY = matmul(XT,Y)
    coefficients = matmul(XTX,XTY)
  End Subroutine FitEmbedding_Process

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
! LMA Vars
    Real(kind=DoubleReal) :: convergenceThreshold
    Real(kind=DoubleReal) :: rss, testRSS, bestRSS, convergence, lambda
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:4) :: J
    Real(kind=DoubleReal), Dimension(1:4,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:4) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:4) :: P      ! Change
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
    convergenceThreshold = 1.0D-8
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
          linCoeffs = LinearRegression(yReg, xReg)
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
! LMA Opt
    convergence = 1.0D0
    lambda = 1.0D0
! ----------------
! Start LMA loop
    Do n=1,100  ! maximum 100 loops
      rss = 0.0D0
! Make Jacobian and Residuals matrix
      Do i=1,pointCount
        R(i) = (parameters(1)*(fitPoints(i,1)**2)*&
        exp(parameters(2)*(fitPoints(i,1)**2))+&
        parameters(3)*(fitPoints(i,1)**2)*&
        exp(parameters(4)*(fitPoints(i,1)**2))&
        )-fitPoints(i,2)   ! f(x)-y
        J(i,1) = (fitPoints(i,1)**2)*exp(parameters(2)*(fitPoints(i,1)**2))                ! d/dx1
        J(i,2) = parameters(1)*(fitPoints(i,1)**4)*exp(parameters(2)*(fitPoints(i,1)**2))  ! d/dx2
        J(i,3) = (fitPoints(i,1)**2)*exp(parameters(4)*(fitPoints(i,1)**2))                ! d/dx3
        J(i,4) = parameters(3)*(fitPoints(i,1)**4)*exp(parameters(4)*(fitPoints(i,1)**2))  ! d/dx4
        rss = rss + R(i)**2
      End Do
      Do k=1,50 ! max 50
! calculate change matrix
! ***********
! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)
! ***********
! Transpose Jacobian
        JT = TransposeMatrix(J)
        JTJ = matmul(JT,J)
        JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
        JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
        JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)
        JTR = matmul(JT,R)
        JTR = -1.0D0*JTR ! Recycle JTR var
        P = matmul(JTJ,JTR)
! Store last loop values
        parameters_Last = parameters
! Update parameters
        Do i=1,size(P)
          parameters(i) = parameters(i) + P(i)
        End Do
! Calc RSS
        testRSS = 0.0D0
        Do i=1,pointCount
          testRSS = testRSS + &
          ((parameters(1)*(fitPoints(i,1)**2)*exp(parameters(2)*(fitPoints(i,1)**2))+&
          parameters(3)*(fitPoints(i,1)**2)*exp(parameters(4)*(fitPoints(i,1)**2))&
          )-fitPoints(i,2))**2
        End Do
! Delayed gratification scheme - 1.5*lambda or 0.2*lambda
        If(testRSS.gt.rss)Then  ! If worse
          lambda = lambda * 1.5D0
          parameters = parameters_Last
          bestRSS = rss
        Else  ! If better
          lambda = lambda * 0.2D0
          bestRSS = testRSS
          Exit
        End If
      End Do
      convergence = abs(testRSS-rss)
! Breakout if convergence threshold met
      If(convergence.lt.convergenceThreshold)Then
        Exit
      End If
! End LMA loop
! ----------------
    End Do
! Store results
    Do i=1,4
      coefficients(i) = parameters(i)
    End Do
  End Subroutine FitDensity_Process


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

  
End Module specialistFunctions