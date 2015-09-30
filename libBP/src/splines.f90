Module splines
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use linearAlgebra
  Use calcFunctions
  Use regression
  Use interpolation
  Use specialistFunctions
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
! --- Functions
  Public :: SplineAB
  Public :: SplineExpThird
  Public :: SplineExpFifth
  Public :: SplineNodes
  Public :: SplineNodesV
  Public :: SplineComplete
  Public :: VaryNode
! --- Subroutines  
  Public :: CompleteNodeData
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

  Function SplineAB(pointA, pointB) RESULT (coefficients)
! Polynomial to spline between points A and B - x, f(x), f'(x) and f''(x) supplied
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:) :: pointA  !1 x, 2 f(x), 3 f'(x), 4 f''(x)....
    Real(kind=DoubleReal), Dimension(:) :: pointB  !1 x, 2 f(x), 3 f'(x), 4 f''(x)....
    Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1))) :: coefficients, yMatrix
    Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1)),1:(2*(size(pointA,1)-1))) :: xMatrix
    Integer(kind=StandardInteger) :: i, j, n, row, col, matrixSize, matrixHalfSize, expt
    Real(kind=DoubleReal) :: x, coeff
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
! solve equation
! coefficients = SolveLinearSet(xMatrix,yMatrix)
    xMatrix = InvertMatrix(xMatrix)
    coefficients = matmul(xMatrix,yMatrix)
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
! xMatrix = InvertMatrix(xMatrix)               ! By matrix inversion
! coefficients = matmul(xMatrix,yMatrix)
  End Function SplineExpFifth

  Function SplineNodes(inputNodes,numDataPoints,startPointIn,endPointIn) RESULT (dataPoints)
! Input nodes x,f(x) for each node, calculate f'(x) and f''(x) from the set of nodes, then spline
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: inputNodes
    Real(kind=DoubleReal), Dimension(1:size(inputNodes,1),1:4) :: splinePoints
    Integer(kind=StandardInteger) :: i, numDataPoints, nodeKey
    Real(kind=DoubleReal), Dimension(1:numDataPoints,1:4) :: dataPoints
    Real(kind=DoubleReal) :: x, xStart, xEnd, xIncrement
    Integer(kind=StandardInteger), Optional :: startPointIn, endPointIn
    Integer(kind=StandardInteger) :: startPoint, endPoint
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
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
! set f'(x) and f''(x)
    Do i=startPoint,endPoint
      x = inputNodes(i,1)
      yArray = PointInterp(inputNodes,x,3,2,startPoint,endPoint)
      splinePoints(i,1) = x
      splinePoints(i,2) = yArray(1)
      splinePoints(i,3) = yArray(2)
      splinePoints(i,4) = yArray(3)
    End Do
! Calculate spline data points
    xStart = splinePoints(startPoint,1)
    xEnd = splinePoints(endPoint,1)
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
      dataPoints(i,1) = x
      dataPoints(i,2) = CalcPolynomial (coefficients, x, 0)
      dataPoints(i,3) = CalcPolynomial (coefficients, x, 1)
      dataPoints(i,4) = CalcPolynomial (coefficients, x, 2)
! Increment x
      x = x + xIncrement
    End Do
  End Function SplineNodes

  Function SplineNodesV(inputNodes,numDataPoints,startPoint,endPoint,dataSize,splineTypeIn,forceCalcDervIn) RESULT (dataPoints)
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
    Integer(kind=StandardInteger) :: numDataPoints, startPoint, endPoint, nodeCount, dataSize
! Declare variables - priv
    Real(kind=DoubleReal), Dimension(1:(endPoint-startPoint+1),1:4) :: splineNodeArr
    Integer(kind=StandardInteger) :: i, nodeKey
    Real(kind=DoubleReal), Dimension(1:dataSize,1:4) :: dataPoints
    Real(kind=DoubleReal) :: x, xStart, xEnd, xIncrement
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:4) :: expThird
    Real(kind=DoubleReal), Dimension(1:6) :: expFifth, polyFitCoeffs
    Real(kind=DoubleReal), Dimension(1:4) :: embFunc, densFunc
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Integer(kind=StandardInteger), Dimension(1:1000), Optional :: splineTypeIn
    Integer(kind=StandardInteger), Dimension(1:1000) :: splineType
    Logical, Dimension(1:1000), Optional :: forceCalcDervIn
    Logical, Dimension(1:1000) :: forceCalcDerv
! optional arguments
    splineType = 1   ! set to standard 5th order polynomial  1 a+bx+...  2 exp(a+bx+...)
    If(present(splineTypeIn))Then
      splineType = splineTypeIn
    End If
    forceCalcDerv = .false.
    If(present(forceCalcDervIn))Then  ! Interp between nodes to calc f'(x) and f''(x)
      forceCalcDerv = forceCalcDervIn
    End If
    nodeCount = endPoint - startPoint + 1
! Init Variables
    dataPoints = 0.0D0
    If(startPoint.eq.0)Then
      startPoint = 1
    End If
    If(endPoint.eq.0)Then
      endPoint = size(inputNodes,1)
    End If
! Transfer nodes from inputNodes to splineNodeArr
    nodeKey = 0
    Do i=startPoint,endPoint
      nodeKey = nodeKey + 1
      splineNodeArr(nodeKey,1) = inputNodes(i,1)
      splineNodeArr(nodeKey,2) = inputNodes(i,2)
      splineNodeArr(nodeKey,3) = inputNodes(i,3)
      splineNodeArr(nodeKey,4) = inputNodes(i,4)
    End Do
! If required, set f'(x) and f''(x)
    Do nodeKey=1,nodeCount
      If(forceCalcDerv(nodeKey))Then
        x = splineNodeArr(nodeKey,1)
        yArray = PointInterp(splineNodeArr,x,4,2,1,nodeCount)
        splineNodeArr(nodeKey,1) = x
        splineNodeArr(nodeKey,2) = yArray(1)
        splineNodeArr(nodeKey,3) = yArray(2)
        splineNodeArr(nodeKey,4) = yArray(3)
      End If
    End Do
! Calculate spline data points
    xStart = splineNodeArr(1,1)
    xEnd = splineNodeArr(nodeCount,1)
    xIncrement = (xEnd-xStart)/(1.0D0*numDataPoints-1.0D0)
! Loop through data points
    nodeKey = 0
    x = xStart
    Do i=1,numDataPoints
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
        End If
        If(splineType(nodeKey).eq.2)Then           ! exp(3rd order)
          expThird = SplineExpThird(pointA(1),pointA(2),pointA(3),&
          pointB(1),pointB(2),pointB(3))
        End If
        If(splineType(nodeKey).eq.3)Then           ! exp(5th order)
          expFifth = SplineExpFifth(pointA(1),pointA(2),pointA(3),pointA(4),&
          pointB(1),pointB(2),pointB(3),pointB(4))
        End If
        If(splineType(nodeKey).eq.4.and.nodeKey.eq.1)Then           ! embedding function F(p) = a+bp^0.5+bp^2+dp^4
          embFunc = FitEmbedding(splineNodeArr,1,nodeCount)
        End If
        If(splineType(nodeKey).eq.5.and.nodeKey.eq.1)Then           ! density function p(r) = a r^2 exp(b r^2) + c r^2 exp(d r^2)
          densFunc = FitDensity(splineNodeArr,1,nodeCount)
        End If
        If(splineType(nodeKey).eq.6.and.nodeKey.eq.1)Then           ! density function p(r) = 5th order poly
          polyFitCoeffs = PolyFit(splineNodeArr,5)
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
        dataPoints(i,2) = embFunc(1)+embFunc(2)*x**0.5D0+embFunc(3)*x**2.0D0+embFunc(4)*x**4.0D0
        dataPoints(i,3) = 0.5D0*embFunc(2)*x**(-0.5D0)+2.0D0*embFunc(3)*x+4.0D0*embFunc(4)*x**3.0D0
        dataPoints(i,4) = -0.25D0*embFunc(2)*x**(-1.5D0)+2.0D0*embFunc(3)+12.0D0*embFunc(4)*x**2.0D0
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
! Increment x
      x = x + xIncrement
    End Do
  End Function SplineNodesV

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
            exit
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

! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------  
  
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

End Module splines