Module activityFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use calcFunctions
  Use solveFunctions
  Use laplaceTransforms
! Force declaration of all variables
  Implicit None
! Public variables  
! Make private
  Private
! Public
! --functions--!
  Public :: CalcIsotopeAmount
  Public :: CalcIsotopeAmountGS
  Public :: MaxTrajDepth
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------  
 
  Function CalcIsotopeAmount(w,decayDataArray,t,calcOptionIn) RESULT (isotopeChange)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j,decaySteps,decayStepCounter, noChanges
    Integer(kind=StandardInteger), optional :: calcOptionIn
    Integer(kind=StandardInteger) :: calcOption
    Real(kind=DoubleReal) :: halfLifeChange, randNumber, w, t
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange
    Real(kind=DoubleReal) :: stableLimit
! Quadrupole Reals
    Real(kind=QuadrupoleReal) :: resultQ, resultGS, tQ, tempQ
    Real(kind=QuadrupoleReal), Dimension(1:20) :: L ! Lambda
    Real(kind=QuadrupoleReal), Dimension(1:20) :: N ! Starting number of atoms
    Real(kind=QuadrupoleReal), Dimension(1:20) :: E ! Exp
    Real(kind=QuadrupoleReal), Dimension(1:19) :: B ! Exp
! -------------------------------------------------
! decaySteps really means decay isotopes in chain (steps = decaySteps-1)
! -------------------------------------------------
! Input decay chain array:
! decayDataArray(i,1) !Tally key
! decayDataArray(i,2) !No. Atoms
! decayDataArray(i,3) !Half life
! decayDataArray(i,4) !branching factor
! decayDataArray(i,5) !isotope Z
! decayDataArray(i,6) !isotope A
! -------------------------------------------------
! Output decay chain array:
! isotopeChange(i,1)    !Tally key
! isotopeChange(i,2)    !Change in isotope amount
! isotopeChange(i,3)    !Start amount
! isotopeChange(i,4)    !End amount
! isotopeChange(i,5)    !Isotope Z
! isotopeChange(i,6)    !Isotope A
! isotopeChange(i,7)    !T1/2
! isotopeChange(i,8)    !Decay constant
! isotopeChange(i,9)    !Branching factor
! isotopeChange(i,10)   !Parent production rate
! isotopeChange(i,11)   !Time
! isotopeChange(i,12)   !GS End
! -------------------------------------------------
! Optional arguments
    calcOption = 1  !(1) 1-4 analytic 5+ SG, (2)1+  SG,  (3) 1-4 analytic+GS 5+ SG
    If(Present(calcOptionIn))Then
      calcOption = calcOptionIn
    End If
! Init variables
    tQ = t
    resultQ = 0.0D0
    resultGS = 0.0D0
! -------------------------------------------------
! Alter decay chain
! -------------------------------------------------
! - If dTime * decay constant lt 1.0D-14 then assume stable for purposes of simulation
    decayStepCounter = 0
    Do i=1,size(decayDataArray,1)
      stableLimit = (log(2.0D0)/decayDataArray(i,3))*t
      decayStepCounter = decayStepCounter + 1
      If(stableLimit.lt.1.0D-14)Then
        decayDataArray(i,3) = -1    !set as stable
        Exit
      End If
    End Do
! Resize array
    decayDataArray = ArraySize2DDouble(decayDataArray,decayStepCounter)
! -------------------------------------------------
! Set stable isotope decay constant very small to avoid infinity error
! -------------------------------------------------
    Do i=1,size(decayDataArray,1)
      If(decayDataArray(i,3).eq.(-1))Then
        decayDataArray(i,3) = 1.0D100
      End If
    End Do
! -------------------------------------------------
! Break same decay constants by ~1E-3% to avoid singularities
! -------------------------------------------------
    noChanges = 0
    Do While(noChanges.eq.0)
      noChanges = 1
      Do i=1,size(decayDataArray,1)
        Do j=1,size(decayDataArray,1)
          If(i.ne.j)Then
            If(decayDataArray(i,3).eq.decayDataArray(j,3))Then
              Call RANDOM_NUMBER(randNumber)
              halfLifeChange = 0.1D0+randNumber*0.9D0
              halfLifeChange = decayDataArray(i,3)*1D-5*halfLifeChange
              decayDataArray(i,3) = decayDataArray(i,3)+halfLifeChange
              decayDataArray(j,3) = decayDataArray(j,3)-halfLifeChange
              noChanges = 0
            End If
          End If
        End Do
      End Do
    End Do
! set decay steps/isotopes
    decaySteps = size(decayDataArray,1)
! allocate isotopeChange array
    Allocate(isotopeChange(1:decaySteps,1:12))
! Fill with starting data
    Do i=1,decaySteps
      isotopeChange(i,1) = decayDataArray(i,1)
      isotopeChange(i,2) = 0.0D0          !default no change
      isotopeChange(i,3) = decayDataArray(i,2)
      isotopeChange(i,4) = decayDataArray(i,2)    !default no change
      isotopeChange(i,5) = decayDataArray(i,5)
      isotopeChange(i,6) = decayDataArray(i,6)
      isotopeChange(i,7) = decayDataArray(i,3)
      isotopeChange(i,8) = log(2.0D0)/decayDataArray(i,3)
      isotopeChange(i,9) = decayDataArray(i,4)
      isotopeChange(i,10) = w
      isotopeChange(i,11) = t
      isotopeChange(i,12) = 0.0D0          !default no change
    End Do
! Store lambda starting atom number data
    Do i=1,decaySteps
      If(decayDataArray(i,3).gt.9.9D99)Then
        L(i) = 0.0D0
      Else
        L(i) = lnTwoQ/isotopeChange(i,7)
      End If
      N(i) = isotopeChange(i,3)
      tempQ = -1.0D0*L(i)*tQ
      E(i) = exp(tempQ)
      B(i) = decayDataArray(i,4)
    End Do
!
! nP -> nA -> nB -> nC -> nD ...
!
! Set starting variables
    If(decaySteps.ge.1)Then
! calc nP
      If(calcOption.eq.1)Then
        resultQ = (w/L(1))*(1-E(1))+N(1)*E(1)
        isotopeChange(1,4) = dble(resultQ)
      End If
      If(calcOption.eq.2.or.ISNAN(resultQ))Then ! solve numerically
        resultGS = CalcIsotopeAmountGS(tQ,1,isotopeChange)
        isotopeChange(1,4) = dble(resultGS)
      End If
    End If
    If(decaySteps.ge.2)Then
! calc nA
      If(calcOption.eq.1)Then ! solve numerically
        resultQ = B(2)*L(1)*w*(1.0D0/(L(1)*L(2))+E(1)/(L(1)*(L(1)-L(2)))-&
        E(2)/(L(2)*(L(1)-L(2))))+&
        B(2)*L(1)*N(1)*(E(1)/(L(2)-L(1))+E(2)/(L(1)-L(2)))+&
        N(2)*E(2)
        isotopeChange(2,4) = dble(resultQ)
      End If
      If(calcOption.eq.2.or.ISNAN(resultQ))Then ! solve numerically
        resultGS = CalcIsotopeAmountGS(tQ,2,isotopeChange)
        isotopeChange(2,4) = dble(resultGS)
      End If
    End If
    If(decaySteps.ge.3)Then
! child B terms
      If(calcOption.eq.1)Then
        resultQ = &
        w*B(2)*B(3)*L(1)*L(2)*&                   ! Term 1
        (1.0D0/(L(1)*L(2)*L(3))-&
        E(1)/(L(1)*(L(1)-L(2))*(L(1)-L(3)))+&
        E(2)/(L(2)*(L(1)-L(2))*(L(2)-L(3)))+&
        E(3)/(L(3)*(L(1)-L(3))*(L(3)-L(2))))+&
        B(2)*B(3)*L(1)*L(2)*N(1)*&                ! Term 2
        (E(1)/((L(1)-L(2))*(L(1)-L(3)))-&
        E(2)/((L(1)-L(2))*(L(2)-L(3)))-&
        E(3)/((L(1)-L(3))*(L(3)-L(2))))+&
        B(3)*L(2)*N(2)*&                          ! Term 3
        (E(1)/(L(2)-L(1))+E(2)/(L(1)-L(2)))+&
        N(3)*E(3)
        isotopeChange(3,4) = dble(resultQ)
      End If
      If(calcOption.eq.2.or.ISNAN(resultQ))Then ! solve numerically
        resultGS = CalcIsotopeAmountGS(tQ,3,isotopeChange)
        isotopeChange(3,4) = dble(resultGS)
      End If
    End If
    If(decaySteps.ge.4)Then
! child C terms
      If(calcOption.eq.1)Then
        resultQ = &
        w*B(2)*B(3)*B(4)*L(1)*L(2)*L(3)*&          ! Term 1
        (&
        1.0D0/(L(1)*L(2)*L(3)*L(4))&
        +E(1)/(L(1)*(L(1)-L(2))*(L(1)-L(3))*(L(1)-L(4)))&
        -E(2)/(L(2)*(L(1)-L(2))*(L(1)-L(3))*(L(2)-L(4)))&
        -E(3)/(L(3)*(L(1)-L(3))*(L(3)-L(2))*(L(3)-L(4)))&
        -E(4)/(L(4)*(L(1)-L(4))*(L(4)-L(2))*(L(4)-L(3)))&
        )+&
        B(2)*B(3)*L(1)*L(2)*N(1)*&                  ! Term 2
        (&
        E(2)/((L(1)-L(2))*(L(2)-L(3))*(L(2)-L(4)))&
        -E(1)/((L(1)-L(2))*(L(1)-L(3))*(L(1)-L(4)))&
        +E(3)/((L(1)-L(3))*(L(3)-L(2))*(L(3)-L(4)))&
        +E(4)/((L(1)-L(4))*(L(4)-L(2))*(L(4)-L(3)))&
        )+&
        B(3)*B(4)*L(2)*L(3)*N(2)*&                   ! Term 3
        (&
        E(2)/((L(2)-L(3))*(L(2)-L(4)))&
        -E(3)/((L(2)-L(3))*(L(3)-L(4)))&
        -E(4)/((L(2)-L(4))*(L(4)-L(3)))&
        )+&
        B(4)*L(3)*N(3)*&                   ! Term 4
        (&
        E(3)/(L(4)-L(3))&
        +E(4)/(L(3)-L(4))&
        )+&
        E(4)*N(4)
        isotopeChange(4,4) = dble(resultQ)
      End If
      If(calcOption.eq.2.or.ISNAN(resultQ))Then ! solve numerically
        resultGS = CalcIsotopeAmountGS(tQ,4,isotopeChange)
        isotopeChange(4,4) = dble(resultGS)
      End If
    End If
! Numeric inverse laplace for remainder
    If(decaySteps.ge.5)Then
      Do i=4,decaySteps
        resultGS = CalcIsotopeAmountGS(tQ,i,isotopeChange)
        isotopeChange(i,4) = dble(resultGS)
        isotopeChange(i,12) = dble(resultGS)
      End Do
    End If
! Adjust the isotope values
    Do i=1,decaySteps
      If(isotopeChange(i,4).lt.0.0D0)Then
        isotopeChange(i,4) = 0.0D0
      End If
      If(isotopeChange(i,12).lt.0.0D0)Then
        isotopeChange(i,12) = 0.0D0
      End If
    End Do
! Store changes in isotope amounts
    Do i=1,size(isotopeChange,1)
      isotopeChange(i,2) = isotopeChange(i,4) - isotopeChange(i,3)
    End Do
  End Function CalcIsotopeAmount

  Function CalcIsotopeAmountGS(t,isotopeStep,isotopeChangeIn) RESULT (output)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i, isotopeStep, M, k
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChangeIn
    Real(kind=QuadrupoleReal), Dimension(1:50) :: weightingQ
    Real(kind=QuadrupoleReal), Dimension(1:20) :: L ! Lambda
    Real(kind=QuadrupoleReal), Dimension(1:20) :: N ! Starting number of atoms
    Real(kind=QuadrupoleReal), Dimension(1:20) :: B ! Starting number of atoms
    Real(kind=QuadrupoleReal) :: kQ, w, t, ft, s, FS, output
! -------------------------------------------------
! Output decay chain array:
! isotopeChange(i,1)    !Tally key
! isotopeChange(i,2)    !Change in isotope amount
! isotopeChange(i,3)    !Start amount
! isotopeChange(i,4)    !End amount
! isotopeChange(i,5)    !Isotope Z
! isotopeChange(i,6)    !Isotope A
! isotopeChange(i,7)    !T1/2
! isotopeChange(i,8)    !Decay constant
! isotopeChange(i,9)    !Branching factor
! isotopeChange(i,10)   !Parent production rate
! -------------------------------------------------
! Init variables
    M = 8
    weightingQ = GaverStehfestWeightingQ(M,weightingQ)
    w = isotopeChangeIn(1,10)
    output = 0.0D0
! Adjust the isotope values
    Do i=1,isotopeStep
      If(isotopeChangeIn(i,4).lt.0.0D0)Then
        isotopeChangeIn(i,4) = 0.0D0
      End If
    End Do
! Store lambda starting atom number data
    Do i=1,isotopeStep
      L(i) = lnTwoQ/isotopeChangeIn(i,7)
      N(i) = isotopeChangeIn(i,3)
      If(i.eq.1)Then
        B(i) = 1.0D0
      Else
        B(i) = isotopeChangeIn(i,9)
      End If
    End Do
! Perform calculation
    ft = 0.0D0
    Do k=1,2*M
      kQ = 1.0D0 * k
      s = (kQ*lnTwoQ)/t
! -----------------------
      FS = (1.0D0/(s+L(1)))*(w/s+N(1))
      Do i=2,isotopeStep
        FS = (1.0D0/(s+L(i)))*(B(i)*L(i-1)*FS+N(2))
      End Do
! FS = (1.0D0/(s+L(1)))*(w/s+N(1))
! -----------------------
      ft = ft + weightingQ(k)*FS
    End Do
    ft = (lnTwoQ/t)*ft
    output = Dble(ft)
! isotopeChangeOut(isotopeStep,4) = Dble(ft)
  End Function CalcIsotopeAmountGS

  Function MaxTrajDepth(coefficients, maxDepthIn) RESULT (maxDepth)
! Calc max depth (E = 0) for ion trajectory, described by polynomial
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, j, nMax
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal) :: c, x, y, maxDepth, maxDepthL
    Real(kind=DoubleReal), Optional :: maxDepthIn
! Set optional argument
    maxDepth = 1.0D10     ! 1m in ang
    If(Present(maxDepthIn))Then
      maxDepth = maxDepthIn
    End If
! Do three refinement loops
    Do i=1,4
      nMax = 20+(i*10)
      c = 10**(log10(maxDepth)/(1.0D0*nMax))
      Do j=1,50
        If(i.lt.3)Then
          x = 1.0D0*c**j
        Else
          x = (maxDepth+1.0D0)-c**(nMax-j)
        End If
        y = CalcPolynomial(coefficients, x)
        If(y.lt.0.0D0)Then
          maxDepth = x
          Exit
        End If
        If(i.eq.4)Then
          maxDepthL = x
        End If
      End Do
    End Do
    maxDepth = SolvePolynomial (coefficients, maxDepthL, maxDepth, 1.0D-6)
  End Function MaxTrajDepth  

! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------  

  Function ArraySize1DDouble (inputArray,arraySize) RESULT (outputArray)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: arraySize
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:arraySize))
! transfer data
    Do i=1,arraySize
      If(i.le.size(inputArray))Then
        outputArray(i) = inputArray(i)
      Else
        outputArray(i) = 0.0D0
      End If
    End Do
  End Function ArraySize1DDouble

  Function ArraySize2DDouble (inputArray,arraySizeHeight,arraySizeWidthIn) &
    RESULT (outputArray)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i, j
    Integer(kind=StandardInteger) :: arraySizeHeight
    Integer(kind=StandardInteger), optional :: arraySizeWidthIn
    Integer(kind=StandardInteger) :: arraySizeWidth
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: outputArray
! catch optional width
    If(present(arraySizeWidthIn))Then
      arraySizeWidth = arraySizeWidthIn
    Else
      arraySizeWidth = size(inputArray,2)
    End If
! Allocate output array
    Allocate(outputArray(1:arraySizeHeight,1:arraySizeWidth))
! transfer data
    Do i=1,arraySizeHeight
      Do j=1,arraySizeWidth
        If(i.le.size(inputArray,1).and.j.le.size(inputArray,2))Then
          outputArray(i,j) = inputArray(i,j)
        Else
          outputArray(i,j) = 0.0D0
        End If
      End Do
    End Do
  End Function ArraySize2DDouble



End Module activityFunctions