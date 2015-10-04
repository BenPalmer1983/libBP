Module fitting
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use strings
  Use constants
  Use matrix
  Use rng
  Use calcFunctions
  Use regression
  Use interpolation
  Use splines
  Use lmaM
! Force declaration of all variables
  Implicit None
! Public variables  
  Character(Len=128), Dimension(1:20) :: fittingReport
! Make private
  Private
! Public
! --variables--!
  Public :: fittingReport
! --functions--!
  Public :: BirchMurnFit
  Public :: ExpFit
  Public :: SingleDecayFit
  Public :: DoubleDecayFit
  Public :: TripleDecayFit
  Public :: FittingPoints    
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
  
    
  
! ----------------------------------------------------------------------------------
! Fitting Points
! Use fitting functions, regression functions
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
    If(calcFunctionT(1:6).eq."SPLINE")Then 
      dataPointsOut = SplinePoints(dataPointsIn, pointsOutCount)
    End If  
! Interp Fit
    If(calcFunctionT(1:6).eq."INTERP")Then 
      dataPointsOut = SplinePoints(dataPointsIn, pointsOutCount)
    End If  
    

  End Function FittingPoints
  
End Module fitting