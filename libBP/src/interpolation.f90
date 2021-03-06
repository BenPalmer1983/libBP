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
  Public :: Lagrange_FX
  Public :: Lagrange_dFX
  Public :: Lagrange_nFX
  Public :: PointInterp
  Public :: InterpPoints
  Public :: FullInterp
  Public :: FullInterpPoints
  Public :: PointInterp3DArr
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

  Function InterpLagrange(x, points, derivativeIn) RESULT (output)
! Calculates y(x), y'(x) or y''(x) using Lagrange interpolation
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Integer(kind=StandardInteger), Optional :: derivativeIn
! Vars:  Out
    Real(kind=DoubleReal) :: output
! Vars:  Private
    Integer(kind=StandardInteger) :: derivative
! Initialise variables
    output = 0.0D0
! Handle optional argument
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! y(x)
    output = Lagrange_nFX(x, points, derivative)
  End Function InterpLagrange


  Function Lagrange_FX(x, points) RESULT (output)
! Calculates F(x) using Lagrange interpolation
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension( : , : ) :: points
! Vars:  Out
    Real(kind=DoubleReal) :: output
! Vars:  Private
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Integer(kind=StandardInteger) :: n, k
    Real(kind=DoubleReal) :: numerator, denominator
    Real(kind=DoubleReal) :: y
! Initialise variables
    output = 0.0D0
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
  End Function Lagrange_FX

  Function Lagrange_dFX(x, points) RESULT (output)
! Calculates F'(x) using Lagrange interpolation
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension( : , : ) :: points
! Vars:  Out
    Real(kind=DoubleReal) :: output
! Vars:  Private
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Integer(kind=StandardInteger) :: i, n, k
    Real(kind=DoubleReal) :: numerator, denominator, numeratorPart
    Real(kind=DoubleReal) :: dy
! Initialise variables
    output = 0.0D0
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
  End Function Lagrange_dFX

  Function Lagrange_nFX(x, points, n) RESULT (output)
! Calculates F^n(x) using Lagrange interpolation (nth derivaive)
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Integer(kind=StandardInteger) :: n
! Vars:  Out
    Real(kind=DoubleReal) :: output
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:2) :: pointsTempA
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:2) :: pointsTempB
! F(x)
    If(n.le.0)Then
      output = Lagrange_FX(x, points)
    End If
! F'(x)
    If(n.eq.1)Then
      output = Lagrange_dFX(x, points)
    End If
! Fn(x)
    If(n.ge.2)Then
      pointsTempA = points
      Do i=2,n
        Do j=1,size(pointsTempA,1)
          pointsTempB(j,1) = pointsTempA(j,1)
          pointsTempB(j,2) = Lagrange_dFX(points(j,1), pointsTempA)
        End Do
        pointsTempA = pointsTempB
      End Do
      output = Lagrange_dFX(x, pointsTempA)
    End If
  End Function Lagrange_nFX

  Function PointInterp(points,x,subsetSize,derivativeIn,inputSetStartIn,inputSetLengthIn,verboseIn) RESULT (yArray)
! Takes large set of data points, finds region of points around the input "x", and interps with lagrange
! points - input data points
! x - value to interpolate y at
! subsetSize - number of points to use in interpolation
! derivativeIn - calc y, y and y' or y. y' and y''
! inputSetStartIn - data points input starts at this number
! inputSetLengthIn - length/number of input data points
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger), optional :: derivativeIn, inputSetStartIn, inputSetLengthIn
    Integer(kind=StandardInteger) :: subsetSize, inputStart, inputLength, derivative
    Logical, optional :: verboseIn
    Logical :: verbose
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
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
        print *,"Input Start: ",inputStart,"Input End: ",inputEnd
        print *,"Start/Lower Val: ",xLower,"End/Upper Val: ",xUpper
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
      yArray(1) = Lagrange_FX(x, pointsInterp)
    End If
    If(derivative.ge.1)Then
      yArray(2) = Lagrange_dFX(x, pointsInterp)
    End If
    If(derivative.ge.2)Then
      yArray(3) = Lagrange_nFX(x, pointsInterp, 2)
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

! ---------------------------------------------
! 3D filled array
! ---------------------------------------------

  Function PointInterp3DArr(points,x,fN,subsetSize,derivativeIn) RESULT (yArray)
! Data held in a 3D array, first column being the key for each set of data, each data set is full and the same length
! points - array containing data points
! x - point being interpolated around
! fN - key for data array
! subsetSize - number of points used in interpolation
    Implicit None  !Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:,:,:) :: points
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger) :: fN
    Integer(kind=StandardInteger) :: subsetSize
    Integer(kind=StandardInteger), optional :: derivativeIn
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Real(kind=DoubleReal), Dimension(1:subsetSize,1:2) :: pointsInterp
    Real(kind=DoubleReal) :: xLower, xUpper
    Integer(kind=StandardInteger) :: i, dataSetSize, derivative
    Integer(kind=StandardInteger) :: xPos, xPos_New_U, xPos_New_L
    Integer(kind=StandardInteger) :: xPosUpper, xPosLower
    Logical :: loopXPos
! Optional arguments
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! Init vars
    dataSetSize = Size(points,2)
    xPos = 1
    xLower = points(fN,1,1)
    xUpper = points(fN,dataSetSize,1)
! Estimate xPos
    xPos = INT(Floor(((x - xLower) / (xUpper - xLower)) * 1.0D0 * dataSetSize))
    xPos = CheckXpos(1,dataSetSize,xPos,3)
    If(x.ge.(points(fN,xPos-1,1)).and.(x.le.points(fN,xPos+1,1)))Then
! xPos is fine, use this one
    Else
! Search for better
      loopXPos = .true.
      xPos_New_L = xPos
      xPos_New_U = xPos
      i=0
      Do While (loopXPos)
        i = i + 1
        If(i.gt.dataSetSize)Then
          Exit
        End If
        xPos_New_L = xPos_New_L - 1
        xPos_New_U = xPos_New_U + 1
        xPos_New_L = CheckXpos(1,dataSetSize,xPos_New_L,3)
        xPos_New_U = CheckXpos(1,dataSetSize,xPos_New_U,3)
        If(x.ge.(points(fN,xPos_New_U-1,1)).and.(x.le.points(fN,xPos_New_U+1,1)))Then
! xPos OK
          xPos = xPos_New_U
          Exit
        End If
        If(x.ge.(points(fN,xPos_New_U-1,1)).and.(x.le.points(fN,xPos_New_U+1,1)))Then
! xPos OK
          xPos = xPos_New_L
          Exit
        End If
      End Do
    End If
! get upper/lower
    Call xPosUpperLower(xPos, subsetSize, xPosUpper, xPosLower)
! Failsafe
    If(xPosLower.lt.1)Then
      xPosLower = 1
    End If
    If((xPosLower+subsetSize).gt.dataSetSize)Then
      xPosLower = dataSetSize-subsetSize+1
    End If
! Make set of points to use for interpolation
    Do i=1,subsetSize
      pointsInterp(i,1) = points(fN,xPosLower-1+i,1)
      pointsInterp(i,2) = points(fN,xPosLower-1+i,2)
    End Do
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
  End Function PointInterp3DArr

  Function CheckXpos(xStart, xEnd, xPos, subsetSize) Result (xPosNew)
! Check the xPos and it's surrounding set is in range of the full data set
    Implicit None  !Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger) :: xStart, xEnd, xPos, subsetSize
! Vars:  Out
    Integer(kind=StandardInteger) :: xPosNew
! Vars:  Private
    !Real(kind=DoubleReal) :: halfSubset
    Integer(kind=StandardInteger) :: xPosLower, xPosUpper
    !halfSubset = 0.5D0*(subsetSize-1)
    !xLower = xPos - floor(halfSubset)
    !xUpper = xPos + ceiling(halfSubset)
    Call xPosUpperLower(xPos, subsetSize, xPosUpper, xPosLower)
    xPosNew = xPos
    If(xPosLower.lt.xStart)Then
      xPosNew = xPosNew + xStart - xPosLower
    ElseIf(xPosUpper.gt.xEnd)Then
      xPosNew = xPosNew + xEnd - xPosUpper
    End If
  End Function CheckXpos

!-------------------------------------------------------------------------------
! Subroutines
!-------------------------------------------------------------------------------

  Subroutine xPosUpperLower(xPos, subsetSize, xPosUpper, xPosLower)
! Get upper and lower values of subset about xpos
    Implicit None  !Force declaration of all variables
! Vars:  In/Out
    Integer(kind=StandardInteger) :: xPos, subsetSize, xPosUpper, xPosLower
! Vars:  Private
    Real(kind=DoubleReal) :: halfSubset
! Calculate
    halfSubset = 0.5D0*(subsetSize-1)
    xPosUpper = xPos + ceiling(halfSubset)
    xPosLower = xPos - floor(halfSubset)
  End Subroutine xPosUpperLower







End Module interpolation
