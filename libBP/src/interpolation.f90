Module interpolation
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use calcFunctions
  Use matrix
  Use linearAlgebra
! Force declaration of all variables
  Implicit None
! Public variables  
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: InterpLagrange
  Public :: PointInterp
  Public :: InterpPoints
  Public :: FullInterp
  Public :: FullInterpPoints
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  
  Function InterpLagrange(x, points, derivativeIn) RESULT (output)
! Calculates y(x), y'(x) or y''(x) using Lagrange interpolation
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:2) :: pointsTemp
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Integer(kind=StandardInteger) :: j, i, n, k, derivative
    Integer(kind=StandardInteger), Optional :: derivativeIn
    Real(kind=DoubleReal) :: numerator, denominator, numeratorPart
    Real(kind=DoubleReal) :: x, y, dy, output, xTemp, dyy
! Initialise variables
    output = 0.0D0
! Handle optional argument
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! y(x)
    If(derivative.eq.0)Then
! Make coefficients
      Do n=1,size(points,1)
        numerator = 1.0D0
        denominator = 1.0D0
        Do k=1,size(points,1)
          If(k.ne.n)Then
            numerator=numerator*(x-points(k,1))
            denominator=denominator*(points(n,1)-points(k,1))
          End If
        End Do
        coefficients(n)=1.0D0*(numerator/denominator)
      End Do
! Calculate y
      y = 0.0D0
      Do n=1,size(points,1)
        y=y+points(n,2)*coefficients(n)
      End Do
      output = y
    End If
! y'(x)
    If(derivative.eq.1)Then
      Do n=1,size(points,1)
        numerator = 0.0D0
        denominator = 1.0D0
        Do k=1,size(points,1)
          If(k.ne.n)Then
            denominator=denominator*(points(n,1)-points(k,1))
            numeratorPart = 1.0D0
            Do i=1,size(points,1)
              If(i.ne.n.and.i.ne.k)Then
                numeratorPart=numeratorPart*(x-points(i,1))
              End If
            End Do
            numerator=numerator+numeratorPart
          End If
        End Do
        coefficients(n)=1.0D0*(numerator/denominator)
      End Do
! Calculate dy
      dy = 0.0D0
      Do n=1,size(points,1)
        dy = dy + points(n,2)*coefficients(n)
      End Do
      output = dy
    End If
! y''(x)
    If(derivative.eq.2)Then
! Could use recursive functions for higher orders, but y''(x) high enough for now
! Calculate y'(x) for each x input
      Do j = 1,size(points,1)
        xTemp = points(j,1)
        Do n=1,size(points,1)
          numerator = 0.0D0
          denominator = 1.0D0
          Do k=1,size(points,1)
            If(k.ne.n)Then
              denominator=denominator*(points(n,1)-points(k,1))
              numeratorPart = 1.0D0
              Do i=1,size(points,1)
                If(i.ne.n.and.i.ne.k)Then
                  numeratorPart=numeratorPart*(xTemp-points(i,1))
                End If
              End Do
              numerator=numerator+numeratorPart
            End If
          End Do
          coefficients(n)=1.0D0*(numerator/denominator)
        End Do
! Calculate dy
        dy = 0.0D0
        Do n=1,size(points,1)
          dy = dy + points(n,2)*coefficients(n)
        End Do
        pointsTemp(j,1) = xTemp
        pointsTemp(j,2) = dy
      End Do
! Use the x, y'(x) points to interpolate y''(x)
      Do n=1,size(points,1)
        numerator = 0.0D0
        denominator = 1.0D0
        Do k=1,size(points,1)
          If(k.ne.n)Then
            denominator=denominator*(pointsTemp(n,1)-pointsTemp(k,1))
            numeratorPart = 1.0D0
            Do i=1,size(pointsTemp,1)
              If(i.ne.n.and.i.ne.k)Then
                numeratorPart=numeratorPart*(x-pointsTemp(i,1))
              End If
            End Do
            numerator=numerator+numeratorPart
          End If
        End Do
        coefficients(n)=1.0D0*(numerator/denominator)
      End Do
! Calculate dy
      dyy = 0.0D0
      Do n=1,size(points,1)
        dyy = dyy + pointsTemp(n,2)*coefficients(n)
      End Do
      output = dyy
    End If
  End Function InterpLagrange

  Function PointInterp(points,x,subsetSize,derivativeIn,inputSetStartIn,inputSetLengthIn,verboseIn) RESULT (yArray)
! Takes large set of data points, finds region of points around the input "x", and interps with lagrange
! points - input data points
! x - value to interpolate y at
! subsetSize - number of points to use in interpolation
! derivativeIn - calc y, y and y' or y. y' and y''
! inputSetStartIn - data points input starts at this number
    Implicit None  !Force declaration of all variables
! Declare variables - In/Out
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger), optional :: derivativeIn, inputSetStartIn, inputSetLengthIn
    Integer(kind=StandardInteger) :: subsetSize, inputStart, inputLength, derivative
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Logical, optional :: verboseIn
    Logical :: verbose
! Declare variables
    Real(kind=DoubleReal), Dimension(1:subsetSize,1:2) :: pointsInterp
    Real(kind=DoubleReal) :: xLower, xUpper
    Integer(kind=StandardInteger) :: i, j, dataSetSize, xPos
    Integer(kind=StandardInteger) :: xPosOffset, xPosOffsetR, xPosStart
    Integer(kind=StandardInteger) :: inputEnd, xPosUpper, xPosLower
! Initialise Variables
    xPos = 1
    dataSetSize = size(points,1)
    yArray = 0.0D0
! Handle optional arguments
    inputStart = 1
    inputLength = dataSetSize
    derivative = 0
    If(Present(inputSetStartIn))Then
      inputStart = inputSetStartIn
    End If
    If(Present(inputSetLengthIn))Then
      inputLength = inputSetLengthIn
    End If
    inputEnd = inputStart+inputLength-1
    If(Present(derivativeIn))Then
      derivative = derivativeIn
      If(derivative.lt.0)Then
        derivative = 0
      End If
      If(derivative.gt.2)Then
        derivative = 2
      End If
    End If
    verbose = .false.
    If(Present(verboseIn))Then
      verbose = verboseIn
    End If
! Check data set size
    If(subsetSize.lt.2)Then
      subsetSize = 2
    End If
    If(subsetSize.gt.inputLength)Then
      subsetSize = inputLength
    End If
! Reduce set of data points
    If(subsetSize.eq.inputLength)Then
      j = 0
      Do i=inputStart,inputEnd
        j = j + 1
        pointsInterp(j,1) = points(i,1)
        pointsInterp(j,2) = points(i,2)
      End Do
    ElseIf(subsetSize.lt.inputLength)Then
! Reduce set of data points
      xLower = points(inputStart,1)
      xUpper = points(inputEnd,1)
      If(verbose)Then
        print *,inputStart,inputEnd
        print *,xLower,xUpper
      End If
! Find xPos
      If(x.lt.xLower)Then  !If x lower than data set, use lowest possible points
        xPos = inputStart
      ElseIf(x.gt.xUpper)Then  !If x higher than data set, use highest possible points
        xPos = inputEnd
      Else
! Estimate position
        xPos = INT(Floor(((x - xLower) / (xUpper - xLower)) * 1.0D0 * inputLength) + inputStart)
        If(xPos.lt.inputStart)Then
          xPos = inputStart
        End If
        If((xPos+1).gt.inputEnd)Then
          xPos = inputEnd-1
        End If
        xLower = points(xPos,1)
        xUpper = points(xPos+1,1)
! If estimate is incorrect, search for better value
        If(x.lt.xLower)Then
          xPosStart = xPos
          Do xPos=xPosStart,inputStart,-1    !Search down
            xLower = points(xPos,1)
            xUpper = points(xPos+1,1)
            If(x.le.xUpper.and.x.ge.xLower)Then
              Exit  !xPos found
            End If
          End Do
        End If
        If(x.gt.xUpper)Then
          xPosStart = xPos
          Do xPos=xPosStart,inputEnd,+1    !Search down
            xLower = points(xPos,1)
            xUpper = points(xPos+1,1)
            If(x.le.xUpper.and.x.ge.xLower)Then
              Exit  !xPos found
            End If
          End Do
        End If
      End If
! Adjust xPos to center of subset
      xPosOffset = INT(Floor(1.0D0*subsetSize/2))
      xPosOffsetR = subsetSize - xPosOffset
      xPosLower = xPos - xPosOffset
      xPosUpper = xPos + xPosOffsetR - 1
! Adjust xPos start, so it fits in the range of subset of selected data points
      If(xPosLower.lt.inputStart)Then
        xPosLower = inputStart
        xPosUpper = inputStart + subsetSize - 1
      End If
      If(xPosUpper.gt.inputEnd)Then
        xPosLower = inputEnd - subsetSize + 1
        xPosUpper = inputEnd
      End If
! Transfer data points to pointsInterp
      j = 0
      Do i=xPosLower,xPosUpper
        j = j + 1
        pointsInterp(j,1) = points(i,1)
        pointsInterp(j,2) = points(i,2)
      End Do
    End If
! If verbose
    If(verbose)Then
      Do j=1,subsetSize
        print *,pointsInterp(j,1),pointsInterp(j,2)
      End Do
    End If
! Store interpolation results
    If(derivative.ge.0)Then
      yArray(1) = InterpLagrange(x, pointsInterp)
    End If
    If(derivative.ge.1)Then
      yArray(2) = InterpLagrange(x, pointsInterp, 1)
    End If
    If(derivative.ge.2)Then
      yArray(3) = InterpLagrange(x, pointsInterp, 2)
    End If
  End Function PointInterp
  
  
  
! ---------------------------------------------
! Interpolation Fitting
! ---------------------------------------------
  
  Function InterpPoints(dataPointsIn, pointsOutCount, interpNodes) RESULT (dataPointsOut)
! Force declaration of all variables
    Implicit None
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: dataPointsIn
    Integer(kind=StandardInteger) :: pointsOutCount
    Integer(kind=StandardInteger) :: interpNodes
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:pointsOutCount,1:2) :: dataPointsOut
! Private: Declare variables
    Integer(kind=StandardInteger) :: i, j, n, k, pointsInCount
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal) :: xStart,xEnd,xInc
    Real(kind=DoubleReal) :: xUpper,xLower
    Real(kind=DoubleReal) :: xA, xB
    Real(kind=DoubleReal), Dimension(1:interpNodes,1:2) :: interpArray
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
! Interp node start (e.g. four point interp 1,(2),3,4)
        k = n - 1
        If(k.lt.1)Then
          k = 1
        End If
        If((k+(interpNodes-1)).gt.pointsInCount)Then
          k = pointsInCount-(interpNodes-1)
        End If
        interpArray = 0.0D0
        Do j=1,interpNodes
          interpArray(j,1) = dataPointsIn(k+(j-1),1)
          interpArray(j,2) = dataPointsIn(k+(j-1),2)      
        End Do
      End If
! store output points
      dataPointsOut(i,1) = x
      dataPointsOut(i,2) = InterpLagrange(x,interpArray,0)
    End Do
  End Function InterpPoints
  
  
  Function FullInterp(dataPoints) RESULT (coefficients)
! Ax = y
! Exact fit of (N-1) polynomial to N data points
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: coefficients
! Private: Declare variables
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:size(dataPoints,1)) :: A
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: Y
    Integer(kind=StandardInteger) :: row, col, matSize
! Init
    matSize = size(dataPoints,1)
! Make Y matrix and A matrix    
    Do row = 1,matSize
      Do col = 1,matSize
        A(row,col) = 1.0D0*dataPoints(row,1)**(col-1.0D0)         
      End Do  
      Y(row) = 1.0D0*dataPoints(row,2)
    End Do
! Solve
    coefficients = SolveLinearSet(A,Y) 
  End Function FullInterp
  
  Function FullInterpPoints(dataPoints, pointsOutCount) RESULT (dataPointsOut)
! Force declaration of all variables
    Implicit None
! In:      Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Integer(kind=StandardInteger) :: pointsOutCount
! Out:     Declare variables
    Real(kind=DoubleReal), Dimension(1:pointsOutCount,1:2) :: dataPointsOut
! Private: Declare variables
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: coefficients
    Integer(kind=StandardInteger) :: i, pointsInCount
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal) :: xStart,xEnd,xInc
! Init
    pointsInCount = size(dataPoints,1)
    xStart = dataPoints(1,1)    
    xEnd = dataPoints(pointsInCount,1)    
    xInc = (xEnd-xStart)/(pointsOutCount-1.0D0)
! Fit Points    
    coefficients = FullInterp(dataPoints)
! Loop through points to make
    Do i=1,pointsOutCount
! Output node x val    
      x = xStart+(i-1)*xInc      
! store output points
      dataPointsOut(i,1) = x
      dataPointsOut(i,2) = CalcPolynomial(coefficients,x)
    End Do  
  End Function FullInterpPoints

End Module interpolation