Module rngDist
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use rng
  Use calcFunctions
  Use regression
  Use interpolation
! Force declaration of all variables
  Implicit None
! Public variables  
  Real(kind=DoubleReal) :: randomDist_randNumberG = -1.0D0
  Real(kind=DoubleReal) :: randomDist_randNumberH = -1.0D0
  Real(kind=DoubleReal), Dimension(0:100,1:2) :: randomDist_inverseInt = 0.0D0
! Make private
  Private
! Public
! --functions--!
  Public :: RandomDist
  Public :: RandomDist_GP
  Public :: RandomVaryPoint
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
 

  Function RandomDist(distTypeIn,setupDistIn,lowerIn,upperIn,sigmaIn) RESULT (output)
! Random float
! F flat distribution
! S square root (diagonal)
! G Gaussian - Box Muller - fits mu=2.5 sigma=0.5 multiplied by 0.05
! H Half Gaussian - Box Muller - fits mu=0.0 sigma=1.0 multiplied by 0.1
! P Test distribution
! M Maxwell-Boltzman distribution
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: lower, upper, sigma
    Character(len=1) :: distType, setupDist  !F flat, G Inverse gaussian distribution
    Real(kind=DoubleReal) :: output
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: randNumber, randNumberB, r, theta, x, y
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), optional :: lowerIn, upperIn, sigmaIn
    Character(len=1), optional :: distTypeIn, setupDistIn
    Real(kind=DoubleReal), Dimension(1:6,1:2) :: fitPoints
    Real(kind=DoubleReal), Dimension(1:21,1:2) :: mbPoints
! Optional arguments
    distType = "F"
    setupDist = "N"
    sigma = 0.4D0
    lower = 0.0D0
    upper = 1.0D0
    If(Present(distTypeIn))Then
      distType = distTypeIn
    End If
    If(Present(setupDistIn))Then
      setupDist = setupDistIn
    End If
    If(Present(sigmaIn))Then
      sigma = sigmaIn
    End If
    If(Present(lowerIn))Then
      lower = lowerIn
    End If
    If(Present(upperIn))Then
      upper = upperIn
    End If
! Get random number 0<= x <=1
    randNumber = RandomLCG()  ! maths.f90
! If Square Root Type
    If(distTypeIn.eq."S")Then
      randNumber = Sqrt(randNumber)
    End If
! If Gaussian Type - full curve, centered on 0.5
! Box Muller method - fits mu=2.5 sigma=0.5 multiplied by 0.05
    If(distTypeIn.eq."G")Then
      If(randomDist_randNumberG.gt.-1.0D0)Then
        randNumber = randomDist_randNumberG
        randomDist_randNumberG = -1.0D0
      Else
! Get second random number
        randNumberB = RandomLCG()
! Calculate x and y
        r = sqrt(-2.0D0*log(randNumber))
        theta = 2.0D0*pi*randNumberB
        x = r*cos(theta)
        y = r*sin(theta)
! Adjust values
        x = (x+5.0D0)/10.0D0
        y = (y+5.0D0)/10.0D0
! store second random number
        randNumber = x
        randomDist_randNumberG = y
      End If
      If(randNumber.lt.lower)Then
        randNumber = lower
      End If
      If(randNumber.gt.upper)Then
        randNumber = upper
      End If
    End If
! If Gaussian Type - half curve, centered on 0.0
! Box Muller method - fits mu=0.0 sigma=1.0 multiplied by 0.1
    If(distTypeIn.eq."H")Then
      If(randomDist_randNumberH.gt.-1.0D0)Then
        randNumber = randomDist_randNumberH
        randomDist_randNumberH = -1.0D0
      Else
! Get second random number
        randNumberB = RandomLCG()  ! maths.f90
! Calculate x and y
        r = sqrt(-2.0D0*log(randNumber))
        theta = 2.0D0*pi*randNumberB
        x = r*cos(theta)
        y = r*sin(theta)
! Adjust values
        x = abs(x)/5.0D0
        y = abs(y)/5.0D0
! store second random number
        randNumber = x
        randomDist_randNumberH = y
      End If
      If(randNumber.lt.lower)Then
        randNumber = lower
      End If
      If(randNumber.gt.upper)Then
        randNumber = upper
      End If
    End If
! Testing dist type - points
    If(distTypeIn.eq."P")Then
      If(setupDist.eq."Y")Then
! Example set of points
        fitPoints(1,1) = 0.0D0
        fitPoints(2,1) = 0.2D0
        fitPoints(3,1) = 0.4D0
        fitPoints(4,1) = 0.6D0
        fitPoints(5,1) = 0.8D0
        fitPoints(6,1) = 1.0D0
        fitPoints(1,2) = 0.0D0
        fitPoints(2,2) = 1.2D0
        fitPoints(3,2) = 0.8D0
        fitPoints(4,2) = 1.8D0
        fitPoints(5,2) = 1.9D0
        fitPoints(6,2) = 0.5D0
        randomDist_inverseInt = RandomDist_GP(fitPoints,"T")  ! maths.f90
      End If
      yArray = PointInterp(randomDist_inverseInt,randNumber,4)  ! maths.f90
      randNumber = yArray(1)
    End If
! Maxwell-Boltzman Distribution
! P(x) = srqt(2/pi)*(x^2*exp(...))
    If(distTypeIn.eq."M")Then
      If(setupDist.eq."Y")Then
! a = 0.25
        Do i=1,21
          mbPoints(i,1) = (i-0)/20.0D0
          mbPoints(i,2) = MaxwellBoltzman(mbPoints(i,1),0.25D0)
        End Do
        randomDist_inverseInt = RandomDist_GP(mbPoints,"T")
      End If
      yArray = PointInterp(randomDist_inverseInt,randNumber,4)  ! maths.f90
      randNumber = yArray(1)
    End If    
! Output (adjust to fall in range)
    output = lower + randNumber*(upper-lower)
  End Function RandomDist

  Function RandomDist_GP(inputPoints, integratorIn) RESULT (outputPoints)
! General purpose distribution function
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: inputPoints ! 6+ pairs
    Integer(kind=LongInteger) :: i
    Real(kind=DoubleReal), Dimension(0:100,1:2) :: outputPoints
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:7) :: coefficientsI
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal) :: maxY, dy, y, x, xLast, iVal, iValMax
    Character(len=1), optional :: integratorIn
    Character(len=1) :: integrator
! Init
    integrator = "P"  ! P Polynomial Fit  T Trapezoidal
    xLast = 0.0D0
    iVal = 0.0D0
! Optional arguments
    If(Present(integratorIn))Then
      integrator = integratorIn
    End If
! Polynomial fit + integration
    If(integrator.eq."P")Then
! Fit polynomial to the data points
      coefficients = PolyFit(inputPoints,5)  ! maths.f90
! Integral coefficients
      coefficientsI(1) = 0.0D0
      Do i=1,6
        coefficientsI(i+1) = (1.0D0/(1.0D0*i))*coefficients(i)
      End Do
! Output points
      maxY = CalcPolynomial(coefficientsI, 1.0D0)  ! maths.f90 ! f(x) should be always positive, Int(f(x)) always increasing
      Do i=0,100
        outputPoints(i,1) = CalcPolynomial(coefficientsI, (i/100.0D0)) / maxY
        outputPoints(i,2) = (i/100.0D0)
      End Do
    End If
! Trapezoidal integrate
    If(integrator.eq."T")Then
      dy = 1.0D0/100.0D0
      Do i=0,100
        y = (i/100.0D0)
        outputPoints(i,2) = y
        If(i.eq.0)Then
          iVal = 0.0D0
          outputPoints(i,1) = iVal
          x = inputPoints(1,1)
        Else
          yArray = PointInterp(inputPoints,y,3)  ! maths.f90
          x = yArray(1)
          If(x.lt.0.0D0)Then
            x = 0.0D0
          End If
          iVal = iVal + (0.5D0*dy*(x+xLast))
          outputPoints(i,1) = iVal
        End If
        xLast = x
      End Do
      iValMax = iVal
! Normalize so total integral = 1
      Do i=0,100
        outputPoints(i,1) = outputPoints(i,1)/iValMax
      End Do
    End If
  End Function RandomDist_GP

  Function RandomVaryPoint(pointValue, maxVariation, sigma) RESULT (output)
! Vary a point +/- a max amount using inverse Gaussian distribution
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: pointValue, maxVariation, sigma
    Real(kind=DoubleReal) :: randNumber, variation, output
! Initialise variables
    output = 0.0D0
! Make varied point
    variation = RandomDist("G","N",0.0D0,maxVariation,sigma)
    Call RANDOM_NUMBER(randNumber)
    If(randNumber.gt.0.5D0)Then
      variation = -1.0D0 * variation
    End If
    output = pointValue + variation
  End Function RandomVaryPoint

End Module rngDist