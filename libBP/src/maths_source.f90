
! ------------------------------------------------------------------------!
! General Maths Functions
! ------------------------------------------------------------------------!

 






! ------------------------------------------------------------------------!
! Polynomial Related Functions
! ------------------------------------------------------------------------!




  Function DerivativePolynomial (coefficientsIn) RESULT (coefficientsOut)
! Calculates p(x) by default, p'(x) for derivativeIn = 1 etc
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension( : ) :: coefficientsIn
    Real(kind=DoubleReal), Dimension(1:(size(coefficientsIn)-1)) :: coefficientsOut
    Integer(kind=StandardInteger) :: j
    Do j=1,size(coefficientsIn,1)-1
      coefficientsOut(j) = j * coefficientsIn(j+1)
    End Do
  End Function DerivativePolynomial


  Function MinimaPolynomial (coefficients, lower, upper) RESULT (x)
! Finds minima in section of polynomial
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal), Dimension(1:size(coefficients,1)-1) :: coefficientsD
    Real(kind=DoubleReal) :: upper, lower
    Real(kind=DoubleReal) :: x, y, xInc, minX, minY
! Find section with minima
    xInc = (upper-lower)/100.0D0
    Do i=0,100
      x = lower + i * xInc
      y = CalcPolynomial(coefficients,x)
      If(i.eq.0)Then
        minX = x
        minY = y
      Else
        If(y.lt.minY)Then
          minY = y
          minX = x
        End If
      End If
    End Do
    coefficientsD = DerivativePolynomial(coefficients)
! Find minimum
    x = SolvePolynomial(coefficientsD,minX-xInc,minX+xInc)
  End Function MinimaPolynomial


! -------------------------------------------------------------------------------------------------!
  Function PolyFitQ(points,order,extendedFitIn) RESULT (coefficientsOut)
! Fits a polynomial of order to the points input
! Input points and output coeffs DP, internal matrix QP
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: k,col,row,exponentValue
    Integer(kind=StandardInteger), Optional :: extendedFitIn
    Integer(kind=StandardInteger) :: extendedFit
    Integer(kind=StandardInteger) :: order
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficientsOut
! Real(kind=QuadrupleReal) :: vRSS, spRSS
    Real(kind=QuadrupleReal), Dimension(1:(order+1),1:1) :: coefficients
    Real(kind=QuadrupleReal), Dimension(1:(order+1),1:(order+1)) :: xMatrix
    Real(kind=QuadrupleReal), Dimension(1:(order+1),1:1) :: yMatrix
! Optional argument
    extendedFit = 0
    If(Present(extendedFitIn))Then
      extendedFit = extendedFitIn
    End If
! Step 1 - Standard fit with Vandermonde matrix
! Init matrices
    xMatrix = 0.0D0
    yMatrix = 0.0D0
    coefficientsOut = 0.0D0
! Build Least Squares Fitting Vandermonde matrix
    Do row=1,(order+1)
      Do col=1,(order+1)
        exponentValue = row+col-2
        Do k=1,size(points,1)
          xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*points(k,1)&
          **exponentValue
        End Do
      End Do
    End Do
    Do row=1,(order+1)
      exponentValue = row-1
      Do k=1,size(points,1)
        yMatrix(row,1) = 1.0D0*yMatrix(row,1)+1.0D0*points(k,2)*&
        points(k,1)**exponentValue
      End Do
    End Do
! invert xMatrix
    xMatrix = InvertMatrixQ(xMatrix)
! multiply inverse by y to get coefficients
    coefficients = MatMultQ(xMatrix,yMatrix)
! move values to output matrix
    Do row=1,(order+1)
      coefficientsOut(row) = Dble(coefficients(row,1))
    End Do
! Step 2 - extended fit
! If(extendedFit.ne.0)Then
! If required, use superposition fit and check if better than vandermonde
! vRSS = CalcResidualSquareSum(points,coefficients)
! coefficientsSP = SP_PolyFit(points,order,50)
! spRSS = CalcResidualSquareSum(points,coefficientsSP)
! If(spRSS.lt.vRSS)Then
!  coefficients = coefficientsSP
! End If
! End If
  End Function PolyFitQ

  Function PolyFitExact(points) RESULT (coefficients)
! Exactly fit polynomial to points, matrix method
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: yMatrix
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:size(points,1)) :: xMatrix
    Integer(kind=StandardInteger) :: i, j
! Initialise variables
    coefficients = 0.0D0
    xMatrix = 0.0D0
    yMatrix = 0.0D0
! Make xMatrix and yMatrix
    Do i=1,size(points,1)
      Do j=1,size(points,1)
        xMatrix(i,j) = points(i,1)**(j-1)
      End Do
      yMatrix(i) = points(i,2)
    End Do
    xMatrix = InvertMatrix(xMatrix)
    coefficients = MatMul(xMatrix,yMatrix)
  End Function PolyFitExact

  Function CalcResidualSquareSum(points,coefficients) RESULT (rss)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension( : ) :: coefficients
    Real(kind=DoubleReal) :: rss, x, y
    rss = 0.0D0
    Do i=1,size(points,1)
      x = 1.0D0*points(i,1)
      y = CalcPolynomial(coefficients,x)
      rss = rss + (y-points(i,2))**2
    End Do
  End Function CalcResidualSquareSum

  Function MinPolyFit(points,order) RESULT (x)
! Fit poly to points, then calculate minimum of the curve, assuming it's in the region of the points (+- 25%)
! and there isn't too much wobble
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Integer(kind=StandardInteger) :: order           ! Largest poly term e.g. x3+x2+x+1 3rd order
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficients
    Real(kind=DoubleReal) :: xMin, xMax, x
    Integer(kind=StandardInteger) :: i
! Init values
    x = 0.0D0
    Do i=1,size(points,1)
      If(i.eq.1)Then
        xMin = points(i,1)
        xMax = points(i,1)
      Else
        If(points(i,1).lt.xMin)Then
          xMin = points(i,1)
        End If
        If(points(i,1).gt.xMax)Then
          xMax = points(i,1)
        End If
      End If
    End Do
    coefficients = PolyFit(points,order)
    x = MinimaPolynomial (coefficients, xMin, xMax)
  End Function MinPolyFit
 
  Function InterpMatrix(points) RESULT (coefficients)
! Calculates coefficients of interpolation polynomial
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, j, matrixSize
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:size(points,1)) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: yMatrix
! Init
    matrixSize = size(points,1)
    coefficients = 0.0D0
    xMatrix = 0.0D0
    yMatrix = 0.0D0
! Build Y and X matrix
    Do i=1,matrixSize
      Do j=1,matrixSize
        xMatrix(i,j) = 1.0D0*points(i,1)**(j-1)
      End Do
      yMatrix(i) = 1.0D0*points(i,2)
    End Do
! Invert xMatrix (reuse xMatrix)
    xMatrix = InvertMatrix(xMatrix)
! calculate coefficients c = x^(-1)y
    coefficients = matmul(xMatrix,yMatrix)
  End Function InterpMatrix

  Function InterpMatrixPoint(x, points, derivativeIn) RESULT (y)
! Calculates values f(x), f'(x), f''(x) of a point using polynomial interpolation
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:size(points,2)) :: pointsTemp
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Integer(kind=StandardInteger) :: i, derivative
    Integer(kind=StandardInteger), Optional :: derivativeIn
    Real(kind=DoubleReal) :: x, y
! Initialise variables
    y = 0.0D0
! Handle optional argument
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! Get coefficients
    coefficients = InterpMatrix(points)
! x,f(x)
    If(derivative.lt.2)Then
! Calculate value of point
      y = CalcPolynomial(coefficients, x, derivative)
    Else
! Original polynomial not accurate enough - update points from x,f(x) to x,f'(x)
      Do i=1,size(points,1)
        pointsTemp(i,1) = points(i,1)
        pointsTemp(i,2) = CalcPolynomial(coefficients, points(i,1), 1)
      End Do
! Get coefficients
      coefficients = InterpMatrix(pointsTemp)
! Calculate value of point
      y = CalcPolynomial(coefficients, x, 1)
    End If
  End Function InterpMatrixPoint

  Function InterpMatrixPointInaccurate(x, points, derivativeIn) RESULT (y)
! Calculates values f(x), f'(x), f''(x) of a point using polynomial interpolation
! Purposeley inaccurate function
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Integer(kind=StandardInteger) :: derivative
    Integer(kind=StandardInteger), Optional :: derivativeIn
    Real(kind=DoubleReal) :: x, y
! Initialise variables
    y = 0.0D0
! Handle optional argument
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! Get coefficients
    coefficients = InterpMatrix(points)
! Calculate value of point
    y = CalcPolynomial(coefficients, x, derivative)
  End Function InterpMatrixPointInaccurate






  Function TripleDecayFit(dataPoints, convergenceThresholdIn) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA) + b exp(lB)
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Optional :: convergenceThresholdIn
    Real(kind=DoubleReal) :: convergenceThreshold
    Integer(kind=StandardInteger) :: i, n, k
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:7) :: output
! Approx regression
    Integer(kind=StandardInteger) :: iA, iB, iC, gridSize, searchBetterFailCount
    Real(kind=DoubleReal) :: maxLambda, minLambda
    Real(kind=DoubleReal) :: lambdaInc
    Real(kind=DoubleReal) :: lA, lB, lC
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: yReg
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:3) :: xReg
    Real(kind=DoubleReal), Dimension(1:3) :: linCoeffs
    Real(kind=DoubleReal) :: linRSS, linBestRSS, bestRSSLastSearch
! LMA Vars
    Real(kind=DoubleReal) :: rss, testRSS, bestRSS, convergence, lambda
    Real(kind=DoubleReal), Dimension(1:6) :: parameters, parameters_Last
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:6) :: J
    Real(kind=DoubleReal), Dimension(1:6,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:6,1:6) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:6,1:6) :: JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:6) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:6) :: P      ! Change
! Optional argument
    convergenceThreshold = 1.0D-8
    If(Present(convergenceThresholdIn))Then
      convergenceThreshold = convergenceThresholdIn
    End If
! Init vars
    parameters = 0.0D0
    parameters_Last = 0.0D0
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
    Do n=1,1000  ! maximum 1000 loops
      rss = 0.0D0
! Make Jacobian and Residuals matrix
      Do i=1,size(dataPoints,1)
        R(i) = (parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
        parameters(3)*exp(parameters(4)*dataPoints(i,1))+&
        parameters(5)*exp(parameters(6)*dataPoints(i,1))&
        )-dataPoints(i,2)   ! f(x)-y
        J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
        J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
        J(i,3) = exp(parameters(4)*dataPoints(i,1))  ! d/dx3
        J(i,4) = dataPoints(i,1)*parameters(3)*exp(parameters(4)*dataPoints(i,1))  ! d/dx4
        J(i,5) = exp(parameters(6)*dataPoints(i,1))  ! d/dx5
        J(i,6) = dataPoints(i,1)*parameters(5)*exp(parameters(6)*dataPoints(i,1))  ! d/dx6
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
        Do i=1,size(dataPoints,1)
          testRSS = testRSS + &
          ((parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
          parameters(3)*exp(parameters(4)*dataPoints(i,1))+&
          parameters(5)*exp(parameters(6)*dataPoints(i,1))&
          )-dataPoints(i,2))**2
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
! Output   f(x) = a exp(b x) + c exp(d x)
    output(1) = parameters(1)  ! a
    output(2) = parameters(2)  ! b
    output(3) = parameters(3)  ! c
    output(4) = parameters(4)  ! d
    output(5) = parameters(5)  ! c
    output(6) = parameters(6)  ! d
    output(7) = bestRSS
  End Function TripleDecayFit

  Function QuadDecayFit(dataPoints, convergenceThresholdIn) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA) + b exp(lB)
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Optional :: convergenceThresholdIn
    Real(kind=DoubleReal) :: convergenceThreshold
    Integer(kind=StandardInteger) :: i, n, k
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:9) :: output
! Approx regression
    Integer(kind=StandardInteger) :: iA, iB, iC, i_D, gridSize, searchBetterFailCount
    Real(kind=DoubleReal) :: maxLambda, minLambda
    Real(kind=DoubleReal) :: lambdaInc
    Real(kind=DoubleReal) :: lA, lB, lC, lD
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: yReg
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:4) :: xReg
    Real(kind=DoubleReal), Dimension(1:4) :: linCoeffs
    Real(kind=DoubleReal) :: linRSS, linBestRSS, bestRSSLastSearch
! LMA Vars
    Real(kind=DoubleReal) :: rss, testRSS, bestRSS, convergence, lambda
    Real(kind=DoubleReal), Dimension(1:8) :: parameters, parameters_Last
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:8) :: J
    Real(kind=DoubleReal), Dimension(1:8,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:8,1:8) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:8,1:8) :: JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:8) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:8) :: P      ! Change
! Optional argument
    convergenceThreshold = 1.0D-8
    If(Present(convergenceThresholdIn))Then
      convergenceThreshold = convergenceThresholdIn
    End If
! Init vars
    parameters = 0.0D0
    parameters_Last = 0.0D0
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
      Do iA=1,gridSize-3
        lA = minLambda + (iA-1)*lambdaInc
        Do i=1,size(dataPoints,1)
          xReg(i,1) = exp(lA*dataPoints(i,1))  ! Make x1 array (for the A exp(lA x) function)
        End Do
        Do iB=iA+1,gridSize-2
          lB = minLambda + (iB-1)*lambdaInc
          Do i=1,size(dataPoints,1)
            xReg(i,2) = exp(lB*dataPoints(i,1))  ! Make x1 array (for the A exp(lA x) function)
          End Do
          Do iC=iB+1,gridSize-1
            lC = minLambda + (iC-1)*lambdaInc
            Do i=1,size(dataPoints,1)
              xReg(i,3) = exp(lC*dataPoints(i,1))
            End Do
            Do i_D=iC+1,gridSize
              lD = minLambda + (i_D-1)*lambdaInc
              Do i=1,size(dataPoints,1)
                xReg(i,4) = exp(lD*dataPoints(i,1))
              End Do
              linCoeffs = LinearRegression(yReg, xReg)
              linRSS = QuadDecayFitRSS(&
              dataPoints, linCoeffs(1), lA, linCoeffs(2), lB,&
              linCoeffs(3), lC, linCoeffs(4), lD)
              If(iA.eq.1.and.iB.eq.2.and.iC.eq.3)Then
                linBestRSS = linRSS
                parameters(1) = linCoeffs(1)
                parameters(2) = lA
                parameters(3) = linCoeffs(2)
                parameters(4) = lB
                parameters(5) = linCoeffs(3)
                parameters(6) = lC
                parameters(7) = linCoeffs(4)
                parameters(8) = lD
              Else
                If(linRSS.lt.linBestRSS)Then
                  linBestRSS = linRSS
                  parameters(1) = linCoeffs(1)
                  parameters(2) = lA
                  parameters(3) = linCoeffs(2)
                  parameters(4) = lB
                  parameters(5) = linCoeffs(3)
                  parameters(6) = lC
                  parameters(7) = linCoeffs(4)
                  parameters(8) = lD
                End If
              End If
            End Do
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
    Do n=1,1000  ! maximum 1000 loops
      rss = 0.0D0
! Make Jacobian and Residuals matrix
      Do i=1,size(dataPoints,1)
        R(i) = (parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
        parameters(3)*exp(parameters(4)*dataPoints(i,1))+&
        parameters(5)*exp(parameters(6)*dataPoints(i,1))+&
        parameters(7)*exp(parameters(8)*dataPoints(i,1))&
        )-dataPoints(i,2)   ! f(x)-y
        J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
        J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
        J(i,3) = exp(parameters(4)*dataPoints(i,1))  ! d/dx3
        J(i,4) = dataPoints(i,1)*parameters(3)*exp(parameters(4)*dataPoints(i,1))  ! d/dx4
        J(i,5) = exp(parameters(6)*dataPoints(i,1))  ! d/dx5
        J(i,6) = dataPoints(i,1)*parameters(5)*exp(parameters(6)*dataPoints(i,1))  ! d/dx6
        J(i,7) = exp(parameters(8)*dataPoints(i,1))  ! d/dx7
        J(i,8) = dataPoints(i,1)*parameters(7)*exp(parameters(8)*dataPoints(i,1))  ! d/dx8
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
        Do i=1,size(dataPoints,1)
          testRSS = testRSS + &
          ((parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
          parameters(3)*exp(parameters(4)*dataPoints(i,1))+&
          parameters(5)*exp(parameters(6)*dataPoints(i,1))+&
          parameters(7)*exp(parameters(8)*dataPoints(i,1))&
          )-dataPoints(i,2))**2
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
! Output   f(x) = a exp(b x) + c exp(d x)
    output(1) = parameters(1)  ! a
    output(2) = parameters(2)  ! b
    output(3) = parameters(3)  ! c
    output(4) = parameters(4)  ! d
    output(5) = parameters(5)  ! c
    output(6) = parameters(6)  ! d
    output(7) = parameters(7)  ! c
    output(8) = parameters(8)  ! d
    output(9) = bestRSS
  End Function QuadDecayFit

  Function SingleDecayFitRSS(dataPoints, a, lA) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, lA, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
  End Function SingleDecayFitRSS

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
    ddfRssCount = ddfRssCount + 1
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

  Function QuadDecayFitRSS(dataPoints, a, lA, b, lB, c, lC, d, lD) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, c, d, lA, lB, lC, lD, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)+b*exp(lB*x)+c*exp(lC*x)+d*exp(lD*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
  End Function QuadDecayFitRSS

  Function FitExpPoly(dataPoints) Result (coefficients)
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
    dataPoints= 0.0D0
    coefficients= 0.0D0
  End Function FitExpPoly




! ------------------------------------------------------------------------!
! Matrix Functions
! ------------------------------------------------------------------------!

! Standard/Double Precision


!  Function LUDecomp(xMatrix) RESULT (lMatrix, uMatrix)

!  End Function LUDecomp





! ------------------------------------------------------------------------!
! Numbers
! ------------------------------------------------------------------------!


! ------------------------------------------------------------------------!
! Random Number Related Functions
! ------------------------------------------------------------------------!


! ------------------------------------------------------------------------!
! Co-ordinates
! ------------------------------------------------------------------------!


! ------------------------------------------------------------------------!
! Physical/Scientific functions
! ------------------------------------------------------------------------!


! ------------------------------------------------------------------------!
! Decay Functions
! ------------------------------------------------------------------------!


! ------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

! ------------------------------------------------------------------------!
! Swap Matrix Rows
! ------------------------------------------------------------------------!



! Subroutine SwapRows_Real(matrix,rowA,rowB)
!  Real(kind=SingleReal), Intent(in) :: matrix(:),rowA,rowB

! End Subroutine SwapRows_Real

! Subroutine SwapRows_Double(matrix,rowA,rowB)
!  Real(kind=DoubleReal), Intent(in) :: matrix(:),rowA,rowB

! End Subroutine SwapRows_Double

! ------------------------------------------------------------------------!
! Random Number Related Subroutines
! ------------------------------------------------------------------------!

  Subroutine SetRandomSeedArray()
    Integer(kind=StandardInteger) :: i,j,k
    Integer(kind=StandardInteger) :: clockReturn, gap, seedTemp, multiple
    Integer(kind=StandardInteger), Dimension(1:1000) :: randomSeed
    Integer(kind=StandardInteger), Dimension(0:9) :: randomIntegers
! random number seed from cpu
    randomIntegers(0) = 25441
    randomIntegers(1) = 37261
    randomIntegers(2) = 1622261
    randomIntegers(3) = 162982
    randomIntegers(4) = 72635
    randomIntegers(5) = 9927151
    randomIntegers(6) = 91
    randomIntegers(7) = 6452
    randomIntegers(8) = 448327
    randomIntegers(9) = 9253411
    j = 0
    Call SYSTEM_CLOCK(clockReturn)
    gap = clockReturn - 10*floor(1.0D0*clockReturn/10)
    Do i=1,1000
      k = j + gap
      If(k.gt.9)Then
        k = k - 10
      End If
      seedTemp = i * (randomIntegers(j)-randomIntegers(k))
      seedTemp = abs(seedTemp)
      multiple = floor(1.0D0 * seedTemp / 1.0D0 * clockReturn)
      seedTemp = seedTemp - multiple * clockReturn
      randomSeed(i) = abs(clockReturn - seedTemp)
      If(j.eq.9)Then
        j = 0
      End If
      j = j + 1
    End Do
    Call RANDOM_SEED(put=randomSeed)
  End Subroutine SetRandomSeedArray



! ------------------------------------------------------------------------!
! Miscellaneous Functions (Activity calc)
! ------------------------------------------------------------------------!


! ------------------------------------------------------------------------!
!                                                                         !
! MODULE [TEST] FUNCTIONS                                                 !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!


! ------------------------------------------------------------------------!
!                                                                         !
! MODULE [OLD] FUNCTIONS                                                 !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

  Function SplineExpOld(xA,fxA,fpxA,xB,fxB,fpxB) RESULT (coefficients)
! Find parameters for spline between points A and B
! form of spline function:  f(x) = exp(a+bx+cx^2+dx^3)
    Implicit None   ! Force declaration of all variables
! Private variables
! In
    Real(kind=DoubleReal) :: xA,fxA,fpxA,xB,fxB,fpxB
! Out
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
! other vars
    Real(kind=DoubleReal) :: a,b,c,d
    Real(kind=DoubleReal) :: a_last,b_last,c_last,d_last
    Real(kind=DoubleReal) :: ln_fxA, ln_fxB
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: xMat
    Real(kind=DoubleReal), Dimension(1:2) :: yMat, aMat
    Real(kind=DoubleReal) :: yA_A, yB_A, yA_B, yB_B
! LMA
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: J, JT, JTJ, JTJ_inv, JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:4) :: R, JTR, JTR_neg, P
    Real(kind=DoubleReal) :: rss, testRSS, lambda
    Integer(kind=StandardInteger) :: n, k
    Logical :: lambdaLoop, keepLooping
    coefficients = 0.0D0
! Starting point solve for f(x) = exp(a+bx)
    ln_fxA = log(fxA)
    ln_fxB = log(fxB)
! Set up matrices
    xMat(1,1) = 1.0D0
    xMat(1,2) = 1.0D0*xA
    xMat(2,1) = 1.0D0
    xMat(2,2) = 1.0D0*xB
    yMat(1) = ln_fxA
    yMat(2) = ln_fxB
! solve
    xMat = InvertMatrix(xMat)
    aMat = matmul(xMat,yMat)
! starting
    a = aMat(1)
    b = aMat(2)
    c = 0.0D0
    d = 0.0D0
    a = 6.42D0
    b = 1.79D0
    c = -4.51D0
    d = 1.08D0
! Starting rss
    yA_A = exp(a+b*xA+c*xA**2+d*xA**3)
    yB_A = (b+2.0D0*c*xA+3.0D0*d*xA**2)*yA_A
    yA_B = exp(a+b*xB+c*xB**2+d*xB**3)
    yB_B = (b+2*c*xB+3*d*xB**2)*yA_B
! starting lambda
    lambda = 1.000D0
    keepLooping = .true.
    n = 0
    Do while(keepLooping.and.n.lt.20)
      n = n + 1
! Make Jacobian and residue matrix
! y = exp(a+b*xA+c*xA**2+d*xA**3)   dy/dx = (b+2*c*xA+3*d*xA**2)*exp(a+b*xA+c*xA**2+d*xA**3)
      yA_A = exp(a+b*xA+c*xA**2+d*xA**3)
      yB_A = (b+2.0D0*c*xA+3.0D0*d*xA**2)*yA_A
      yA_B = exp(a+b*xB+c*xB**2+d*xB**3)
      yB_B = (b+2*c*xB+3*d*xB**2)*yA_B
      rss = 0.0D0
      rss = rss + (yA_A-fxA)**2+(yB_A-fpxA)**2+(yA_B-fxB)**2+(yB_B-fpxB)**2
! point A f(x)
      J(1,1) = yA_A
      J(1,2) = xA*yA_A
      J(1,3) = (xA**2)*yA_A
      J(1,4) = (xA**3)*yA_A
      R(1) = yA_A-fxA
! point B f(x)
      J(2,1) = yA_B
      J(2,2) = xB*yA_B
      J(2,3) = (xB**2)*yA_B
      J(2,4) = (xB**3)*yA_B
      R(2) = yA_B-fxB
! point A f'(x)
      J(3,1) = yB_A
      J(3,2) = yA_A+xA*yB_A
      J(3,3) = 2*xA*yA_A+xA*(xA**2)*yB_A
      J(3,4) = 3*(xA**2)*yA_A+(xA**3)*yB_A
      R(3) = yB_A-fpxA
! point B f'(x)
      J(4,1) = yB_B
      J(4,2) = yA_B+xB*yB_B
      J(4,3) = 2*xB*yA_B+xB*(xB**2)*yB_B
      J(4,4) = 3*(xB**2)*yA_B+(xB**3)*yB_B
      R(4) = yB_B-fpxB
      lambdaLoop = .true.
      k = 0
      Do While(lambdaLoop.and.k.le.10)
        k = k + 1
        JT = TransposeMatrix(J)
        JTJ = matmul(JT,J)
        JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
        JTJ = MatAdd(JTJ,JTJ_Diag) ! update JTJ
        JTJ_inv = InvertMatrix(JTJ)
        JTR = matmul(JT,R)
        JTR_neg = -1.0D0*JTR
        P = matmul(JTJ_inv,JTR_neg)
! store last values
        a_last = a
        b_last = b
        c_last = c
        d_last = d
! update
        a = a + P(1)
        b = b + P(2)
        c = c + P(3)
        d = d + P(4)
! Calc RSS
        yA_A = exp(a+b*xA+c*xA**2+d*xA**3)
        yB_A = (b+2.0D0*c*xA+3.0D0*d*xA**2)*yA_A
        yA_B = exp(a+b*xB+c*xB**2+d*xB**3)
        yB_B = (b+2*c*xB+3*d*xB**2)*yA_B
        testRSS = 0.0D0
        testRSS = testRSS + (yA_A-fxA)**2+(yB_A-fpxA)**2+(yA_B-fxB)**2+(yB_B-fpxB)**2
        print *,n,k,a,b,c,d, lambda, rss, testRSS
        If(testRSS.lt.rss)Then
          lambda = 0.2D0 * lambda
          lambdaLoop = .false.
        Else
          lambda = 1.5D0 * lambda
          a = a_last
          b = b_last
          c = c_last
          d = d_last
        End If
      End Do
      If(testRSS.lt.1.0D-5)Then
        keepLooping = .false.
      End If
    End Do
  End Function SplineExpOld

  Function Spline(inputNodes,numDataPoints,startPointIn,endPointIn) RESULT (dataPoints)
! NOT WORKING YET - 12092014
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: inputNodes
    Integer(kind=StandardInteger) :: i, numDataPoints, nodeKey
    Real(kind=DoubleReal), Dimension(1:numDataPoints,1:4) :: dataPoints
    Real(kind=DoubleReal) :: x, xStart, xEnd, xIncrement
    Integer(kind=StandardInteger), Optional :: startPointIn, endPointIn
    Integer(kind=StandardInteger) :: startPoint, endPoint
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Fill in
!
! Init Variables
    dataPoints = 0.0D0
    startPoint = 1
    endPoint = size(inputNodes,1)
    If(Present(startPointIn))Then
      startPoint = startPointIn
    End If
    If(Present(endPointIn))Then
      endPoint = endPointIn
    End If
    xStart = inputNodes(startPoint,1)
    xEnd = inputNodes(endPoint,1)
    xIncrement = (xEnd-xStart)/(1.0D0*numDataPoints-1)
! Loop through data points
    nodeKey = startPoint-1
    x = xStart
    Do i=1,numDataPoints
      If((i.eq.1).or.(x.ge.inputNodes(nodeKey+1,1).and.(nodeKey+1).lt.endPoint))Then
        nodeKey = nodeKey + 1
        pointA(1) = inputNodes(nodeKey,1)
        pointA(2) = inputNodes(nodeKey,2)
        pointA(3) = inputNodes(nodeKey,3)
        pointA(4) = inputNodes(nodeKey,4)
        pointB(1) = inputNodes(nodeKey+1,1)
        pointB(2) = inputNodes(nodeKey+1,2)
        pointB(3) = inputNodes(nodeKey+1,3)
        pointB(4) = inputNodes(nodeKey+1,4)
        coefficients = SplineAB(pointA, pointB)
      End If
      yArray = CalcPolynomial (coefficients, x, 2)
      dataPoints(i,1) = x
      dataPoints(i,2) = yArray(1)
      dataPoints(i,3) = yArray(2)
      dataPoints(i,4) = yArray(3)
! Increment x
      x = x + xIncrement
    End Do
  End Function Spline

  Function SplineABOld(pointA, pointB) RESULT (coefficients)
! Polynomial to spline between points A and B - x, f(x), f'(x) and f''(x) supplied
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:) :: pointA  !1 x, 2 f(x), 3 f'(x), 4 f''(x)....
    Real(kind=DoubleReal), Dimension(:) :: pointB  !1 x, 2 f(x), 3 f'(x), 4 f''(x)....
    Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1))) :: coefficients, yMatrix
    Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1)),1:(2*(size(pointA,1)-1))) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:(size(pointA,1)-1),1:(2*(size(pointA,1)-1))) :: &
    xMatrixCoeffs, xMatrixExp
    Integer(kind=StandardInteger) :: row, col, matrixSize, matrixHalfSize
! Init variables
    coefficients = 0.0D0
    yMatrix = 0.0D0
    xMatrix = 0.0D0
    matrixSize = 2*(size(pointA,1)-1)
    matrixHalfSize = (size(pointA,1)-1)
! Make xMatrixExp
    Do row=1,matrixHalfSize
      Do col=1,matrixSize
        xMatrixExp(row,col) = col-row
        If(xMatrixExp(row,col).lt.0)Then
          xMatrixExp(row,col) = 0
        End If
      End Do
    End Do
! Make xMatrixCoeffs
    Do col=1,matrixSize
      xMatrixCoeffs(1,col) = 1.0D0
    End Do
    Do row=2,matrixHalfSize
      Do col=1,matrixSize
        xMatrixCoeffs(row,col) = 1.0D0*xMatrixExp(row-1,col)*xMatrixCoeffs(row-1,col)
      End Do
    End Do
! Make xMatrix
    Do row=1,matrixHalfSize
      Do col=1,matrixSize
        xMatrix(row,col) = &
        1.0D0*xMatrixCoeffs(row,col)*pointA(1)**xMatrixExp(row,col)
      End Do
    End Do
    Do row=1,matrixHalfSize
      Do col=1,matrixSize
        xMatrix(row+matrixHalfSize,col) = &
        1.0D0*xMatrixCoeffs(row,col)*pointB(1)**xMatrixExp(row,col)
      End Do
    End Do
! Make yMatrix
    Do row=1,matrixHalfSize
      yMatrix(row) = pointA(row+1)
    End Do
    Do row=1,matrixHalfSize
      yMatrix(row+matrixHalfSize) = pointB(row+1)
    End Do
! Invert xMatrix
!    xMatrix = InvertMatrix(xMatrix)
! Solve equation
!    coefficients = matmul(xMatrix,yMatrix)
    coefficients = SolveLinearSet(xMatrix,yMatrix)
  End Function SplineABOld

End Module maths
