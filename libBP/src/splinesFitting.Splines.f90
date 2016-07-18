! ----------------------------------------------------------------------------------
! Spline Functions!
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
