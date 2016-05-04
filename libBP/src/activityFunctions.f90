! --------------------------------------------------------------!
! Plot
! plotKinds, plotTypes
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Uses Matplotlib python library to build charts
! Requires matplotlib to be installed
!
! Two modules must be used - plotTypes to load the types, and plot to use the plot functions
! and subroutines.  It must be compiled as a part of this library of functions as there are
! quite a lot of modules it depends on.
!
! There is an option to perform one or multiple fits on the input data.
!
! ----------------------------------------
! Updated: 21st September 2015
! ----------------------------------------

Module activityFunctionsTypes
! Setup Modules
  Use kinds

  Type :: decayChainObj
    Real(kind=DoubleReal) :: time = 0.0D0
    !Real(kind=DoubleReal) :: productionRate = 0.0D0
    Integer(kind=StandardInteger) :: isotopes
    Character(Len=16), Dimension(1:100) :: label
    Real(kind=DoubleReal), Dimension(1:100) :: productionRate = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: branchFactor = 1.0D0    ! from isotope parent
    Real(kind=DoubleReal), Dimension(1:100) :: decayConstant = -1.0D0  ! negative for stable
    Real(kind=DoubleReal), Dimension(1:100) :: halfLife = -1.0D0       ! negative for stable
    Real(kind=DoubleReal), Dimension(1:100) :: amountStart = 0.0D0     ! N(0)
    Real(kind=DoubleReal), Dimension(1:100) :: amountEnd = 0.0D0       ! N(t)
    Real(kind=DoubleReal), Dimension(1:100) :: activity = 0.0D0        ! lambda N(t)
  End Type

  Type :: activityTimeObj
    Real(kind=DoubleReal) :: time = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: amount = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: activity = 0.0D0
  End Type

End Module activityFunctionsTypes


Module activityFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use calcFunctions
  Use solveFunctions
  Use laplaceTransforms
  Use activityFunctionsTypes
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
  Public :: ActivityCompareGS
  Public :: decayChainComplete
  Public :: decayChainPrint
  Public :: CalcIsotopeChainGS
  Public :: CalcIsotopeChain
  Public :: CalcActivities
  Public :: CalcActivitiesPrint
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
      If(calcOption.eq.1)Then ! solve analytically
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



! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------

  Subroutine ActivityCompareGS(lambda, dataPointsCount, endTime, nAtoms, productionRate)
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Real(kind=DoubleReal) :: lambda
    Integer(kind=StandardInteger) :: dataPointsCount
    Real(kind=DoubleReal) :: endTime, nAtoms, productionRate
! Private: Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: analyticResult, numericResult, numericResultB, t
    Real(kind=DoubleReal), Dimension(1:10) :: p
    p(1) = productionRate
    p(2) = lambda
    p(3) = nAtoms
!
    Do i=1,dataPointsCount
      t = (i-1)*(endTime/(dataPointsCount-1))
      analyticResult = productionRate*((1/lambda)-(exp(-1.0D0*lambda*t)/lambda))+nAtoms*exp(-1.0D0*lambda*t)
      If(t.eq.0.0D0)Then
        numericResult = nAtoms
      Else
        !numericResult = GaverStehfest(ltDecay,t,p,7)
        !numericResultB = GaverStehfest(ltDecay,t,p,8)
      End If
      print *,t,analyticResult,numericResult,numericResultB
    End Do
  End Subroutine ActivityCompareGS

  Subroutine decayChainComplete(decayChain)
! Completes the decay chain object
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Type(decayChainObj) :: decayChain
! Vars Private
    Integer(kind=StandardInteger) :: i, n
    n = 0
    Do i=1,100
      n = n + 1
      If(decayChain%decayConstant(i).eq.-1.0D0.and.decayChain%halfLife(i).gt.0.0D0)Then
! complete decay constant from half life
        decayChain%decayConstant(i) = lnTwo/decayChain%halfLife(i)
      End If
      If(decayChain%halfLife(i).le.0.0D0.and.decayChain%decayConstant(i).gt.0.0D0)Then
! complete decay constant from half life
        decayChain%halfLife(i) = lnTwo/decayChain%decayConstant(i)
      End If
! Adjust for stable isotope
      If(decayChain%decayConstant(i).lt.0.0D0)Then
        decayChain%halfLife(i) = -1.0D0
        decayChain%decayConstant(i) = 0.0D0
      End If
      If(decayChain%halfLife(i).lt.0.0D0)Then
        decayChain%halfLife(i) = -1.0D0
        decayChain%decayConstant(i) = 0.0D0
      End If
! Break out if stable
      If(decayChain%decayConstant(i).eq.0.0D0)Then
        Exit
      End If
    End Do
    decayChain%isotopes = n
  End Subroutine decayChainComplete


  Subroutine decayChainPrint(decayChain)
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Type(decayChainObj) :: decayChain
    Integer(kind=StandardInteger) :: i
! Output
    Print "(A96)","================================================================================================"
    Print "(A32,F16.8)","Time:                      ",decayChain%time
    Print "(A32,I8)","Isotopes in chain: ",decayChain%isotopes
    Print "(A8,A20,A8,A20,A20,A20,A20)","n","w_n","B_n-1,n","L_n","T1/2_n","nStart","nEnd"
    Do i=1,decayChain%isotopes
      Print "(I8,E20.8,F8.4,E20.8,E20.8,E20.8,E20.8)",&
      i,decayChain%productionRate(i),decayChain%branchFactor(i),decayChain%decayConstant(i),decayChain%halfLife(i),&
      decayChain%amountStart(i),decayChain%amountEnd(i)
    End Do
    Print "(A96)","================================================================================================"
    Print *,""
  End Subroutine decayChainPrint


  Subroutine CalcIsotopeChainGS(decayChain)
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Type(decayChainObj) :: decayChain
! Vars Private
    Integer(kind=StandardInteger) :: n, k
    Real(kind=DoubleReal) :: nEnd
    Real(kind=DoubleReal) :: t
    Real(kind=DoubleReal), Dimension(1:100) :: p
! Complete decay chain data
    Call decayChainComplete(decayChain)
! Store time
    t = decayChain%time
! Prepare parameters array
    p(1) = 0   ! Chain to calculate
    p(2) = decayChain%isotopes
    Do n=1,decayChain%isotopes
      k = 2+4*(n-1)+1
      p(k) = decayChain%productionRate(n)
      k = 2+4*(n-1)+2
      If(decayChain%decayConstant(n).le.0.0D0)Then  ! Set very small decay constant if stable
        p(k) = 1.0D-30
      Else
        p(k) = decayChain%decayConstant(n)
      End If
      k = 2+4*(n-1)+3
      p(k) = decayChain%branchFactor(n)
      k = 2+4*(n-1)+4
      p(k) = decayChain%amountStart(n)
    End Do
! Calculate isotope amounts
    Do n=1,decayChain%isotopes
      p(1) = n
      nEnd = GaverStehfest(ltDecay, t, p, 8)
      decayChain%amountEnd(n) = nEnd
    End Do
  End Subroutine CalcIsotopeChainGS







  Subroutine CalcIsotopeChain(decayChain)
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Type(decayChainObj) :: decayChain
! Vars Private
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: t, nEnd
    Real(kind=DoubleReal), Dimension(1:100) :: W ! Production Rate
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: N ! Starting number of atoms
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Exp
! Complete decay chain data
    Call decayChainComplete(decayChain)
! store input in shortned name arrays to make equations clearer
    t = decayChain%time
    Do i=1,decayChain%isotopes
      W(i) = decayChain%productionRate(i)
      L(i) = decayChain%decayConstant(i)
      B(i) = decayChain%branchFactor(i)
      N(i) = decayChain%amountStart(i)
    End Do
! Break infinities
    Call DecayBreakInfinities(L,decayChain%isotopes)
! Run analytic calculations
! Loop through isotopes
! Using L-1(1/(q+ps) = (1/p)*exp(-1*(q/p)*t) and partial fractions
    Do i=1,decayChain%isotopes
      nEnd = CalcIsotopeChainCalc(t,i,W,L,N,B)
      decayChain%amountEnd(i) = nEnd
      decayChain%activity(i) = nEnd*decayChain%decayConstant(i)
    End Do
  End Subroutine CalcIsotopeChain

  Subroutine DecayBreakInfinities(L,n)
!
    Implicit None ! Force declaration of all variables
! Vars In/Out
    Real(kind=DoubleReal), Dimension(:) :: L ! Lambda
    Integer(kind=StandardInteger) :: n
! Vars Private
    Integer(kind=StandardInteger) :: i,j
    Logical :: breaking
! Loop and alter matching decay constants slightly
    breaking = .true.
    Do While(breaking)
      breaking = .false.
      Do i=1,n-1
        Do j=i+1,n
          If(L(i).eq.L(j))Then
            breaking = .true.
            L(i) = L(i)*1.0000001D0  ! Vary by 0.00001%
          End If
        End Do
      End Do
    End Do
  End Subroutine DecayBreakInfinities

  Function CalcIsotopeChainCalc(t,m,W,L,N,B) Result (nEnd)
    Implicit None ! Force declaration of all variables
! Vars In
    Real(kind=DoubleReal) :: t
    Integer(kind=StandardInteger) :: m
    Real(kind=DoubleReal), Dimension(1:100) :: W ! Production Rate
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: N ! Starting number of atoms
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Branch factor
! Vars Out
    Real(kind=DoubleReal) :: nEnd
! Vars Private
    Integer(kind=StandardInteger) :: k
    Real(kind=DoubleReal) :: mult, multR
    Real(kind=DoubleReal) :: nChange
! Init
    nEnd = 0.0D0
! ---------------------------------------------
! UNSTABLE Isotopes
! ---------------------------------------------
    If(L(m).gt.0.0D0)Then
! Loop through terms
      Do k=1,m
        multR = CalcIsotopeChainMultR(k,m,L,B)
        mult = multR * N(k)
! decay of starting matter
        nChange = CalcIsotopeChainF_Unstable(k,m,t,mult,L)
        nEnd = nEnd + nChange
! production term
        mult = multR * W(k)
        nChange = CalcIsotopeChainG_Unstable(k,m,t,mult,L)
        nEnd = nEnd + nChange
      End Do
    Else
! ---------------------------------------------
! STABLE Isotopes
! ---------------------------------------------
! Loop through terms
      nEnd = N(m)+t*W(m)
      Do k=1,m-1
        multR = CalcIsotopeChainMultR(k,m,L,B)
        mult = multR * N(k)
! decay of starting matter
        nChange = CalcIsotopeChainF_Stable(k,m,t,mult,L)
        nEnd = nEnd + nChange
! production term
        mult = multR * W(k)
        nChange = CalcIsotopeChainG_Stable(k,m,t,mult,L)
        nEnd = nEnd + nChange
      End Do
    End If
  End Function CalcIsotopeChainCalc

  Function CalcIsotopeChainMultR(k,m,L,B) Result (multR)
! Vars In
    Integer(kind=StandardInteger) :: k, m
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Branching Factor
! Vars Out
    Real(kind=DoubleReal) :: multR
! Private
    Integer(kind=StandardInteger) :: i
! Result
    multR = 1.0D0
    If(k.lt.m)Then
      Do i=k,m-1
        multR = multR * B(i+1) * L(i)
      End Do
    End If
  End Function CalcIsotopeChainMultR

! --------------------------
! Unstable
! --------------------------

  Function CalcIsotopeChainF_Unstable(k,m,t,mult,L) Result (nChange)
! Vars In
    Integer(kind=StandardInteger) :: k, m
    Real(kind=DoubleReal) :: t, mult
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
! Vars Out
    Real(kind=DoubleReal) :: nChange
! Private
    Real(kind=DoubleReal) :: multP
    Real(kind=DoubleReal) :: r
    Integer(kind=StandardInteger) :: i, j
! Calculate isotope amount change
    nChange = 0.0D0
    multP = (-1.0D0)**(m-k)
! Loop through pfrac
    Do i=k,m
      r = 1.0D0
      Do j=k,m
        If(j.ne.i)Then
          r = r * (L(i)-L(j))
        End If
      End Do
      nChange = nChange + (1.0D0/r)*exp(-1.0D0*L(i)*t)*multP*mult
    End Do
  End Function CalcIsotopeChainF_Unstable

  Function CalcIsotopeChainG_Unstable(k,m,t,mult,L) Result (nChange)
! Vars In
    Integer(kind=StandardInteger) :: k, m
    Real(kind=DoubleReal) :: t, mult
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
! Vars Out
    Real(kind=DoubleReal) :: nChange
! Private
    Real(kind=DoubleReal) :: multP
    Real(kind=DoubleReal) :: r
    Integer(kind=StandardInteger) :: i, j
! Calculate isotope amount change
    nChange = 0.0D0
    multP = (-1.0D0)**(m-k+1)
  ! term A
    r = 1.0D0
    Do i=k,m
      r = r * L(i)
    End Do
    nChange = nChange + (1.0D0/r)*mult
! Loop through pfrac
    Do i=k,m
      r = 1.0D0*L(i)
      Do j=k,m
        If(j.ne.i)Then
          r = r * (L(i)-L(j))
        End If
      End Do
      nChange = nChange + (1.0D0/r)*exp(-1.0D0*L(i)*t)*multP*mult
    End Do
  End Function CalcIsotopeChainG_Unstable

! --------------------------
! Stable
! --------------------------

  Function CalcIsotopeChainF_Stable(k,mIn,t,mult,L) Result (nChange)
! Vars In
    Integer(kind=StandardInteger) :: k, mIn
    Real(kind=DoubleReal) :: t, mult
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
! Vars Out
    Real(kind=DoubleReal) :: nChange
! Private
    Integer(kind=StandardInteger) :: m
    Real(kind=DoubleReal) :: multP
    Real(kind=DoubleReal) :: r
    Integer(kind=StandardInteger) :: i, j
! In
    m = mIn-1
! Calculate isotope amount change
    nChange = 0.0D0
    multP = (-1.0D0)**(m-k+1)
! term A
    r = 1.0D0
    Do i=k,m
      r = r * L(i)
    End Do
    nChange = nChange + (1.0D0/r)*mult
! Loop through pfrac
    Do i=k,m
      r = 1.0D0*L(i)
      Do j=k,m
        If(j.ne.i)Then
          r = r * (L(i)-L(j))
        End If
      End Do
      nChange = nChange + (1.0D0/r)*exp(-1.0D0*L(i)*t)*multP*mult
    End Do
  End Function CalcIsotopeChainF_Stable


  Function CalcIsotopeChainG_Stable(k,mIn,t,mult,L) Result (nChange)
  ! Vars In
      Integer(kind=StandardInteger) :: k, mIn
      Real(kind=DoubleReal) :: t, mult
      Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
  ! Vars Out
      Real(kind=DoubleReal) :: nChange
  ! Private
      Integer(kind=StandardInteger) :: m
      Real(kind=DoubleReal) :: multP
      Real(kind=DoubleReal) :: p, q, r
      Integer(kind=StandardInteger) :: i, j
  ! In
      m = mIn-1
  ! Calculate isotope amount change
      nChange = 0.0D0
      multP = (-1.0D0)**(m-k)
  ! term A
      r = 1.0D0
      Do i=k,m
        r = r * L(i)
      End Do
      nChange = nChange + (1.0D0/r)*t*mult
  ! term B
      p = CalcIsotopeChainC(L,k,m)
      q = 1.0D0
      Do i=k,m
        q = q*L(i)*L(i)
      End Do
      r = (-1.0D0)*(p/q)
      nChange = nChange + r*mult
  ! Loop through pfrac
      Do i=k,m
        r = 1.0D0*L(i)*L(i)
        Do j=k,m
          If(j.ne.i)Then
            r = r * (L(i)-L(j))
          End If
        End Do
        nChange = nChange + (1.0D0/r)*exp(-1.0D0*L(i)*t)*multP*mult
      End Do
    End Function CalcIsotopeChainG_Stable

  Function CalcIsotopeChainC(L,kIn,mIn) Result (numerator)
! Calculates numerator in isotope activity function
    Implicit None ! Force declaration of all variables
! Vars In
    Real(kind=DoubleReal), Dimension(:) :: L
    Integer(kind=StandardInteger) :: kIn, mIn
    Real(kind=DoubleReal) :: numerator
! Vars Private
    Integer(kind=StandardInteger) :: i, j, k, m, n
    Integer(kind=StandardInteger), Dimension(1:(mIn-kIn)) :: combinationSet
    Logical :: loopCombinations
    Real(kind=DoubleReal) :: tempVal
! init
    n = mIn-kIn+1  ! Set size
    m = n-1
    If(n.eq.1)Then
      numerator = 1.0D0
    Else
      numerator = 0.0D0
! Set up starting combination
      Do i=1,m
        combinationSet(i) = i
      End Do
! Loop through all combinations (order important)
      loopCombinations = .true.
      k = 0
      Do while(loopCombinations)
        k = k + 1
        If(k.gt.1)Then
          j = m
          Do i=1,m
            loopCombinations = .false.
            If(combinationSet(j).lt.(n+1-i))Then
              combinationSet(j) = combinationSet(j) + 1
              loopCombinations = .true.
              Exit
            End If
            j = j - 1
          End Do
        End If
        If(loopCombinations)Then
          tempVal = 1.0D0
          Do i=1,m
            tempVal = tempVal * L(combinationSet(i)+kIn-1)
          End Do
          numerator = numerator + tempVal
        End If
      End Do
    End If
  End Function CalcIsotopeChainC

! --------------------------
! Calc Activities
! --------------------------

  Subroutine CalcActivities(decayChain,activityTime,endTime)
! calculate activities over time
! Vars In
    Type(decayChainObj) :: decayChain
    Type(activityTimeObj), Dimension(:) :: activityTime
    Real(kind=DoubleReal) :: endTime
! Vars Private
    Integer(kind=StandardInteger) :: i, j, timeSteps
    Real(kind=DoubleReal) :: t
! Time steps
    timeSteps = Size(activityTime,1)
! Loop through steps and calculate activity
    Do i=1,timeSteps
      t = (1.0D0*(i/(1.0D0*timeSteps)))*endTime
      activityTime(i)%time = t
      decayChain%time = t
      Call CalcIsotopeChain(decayChain)
      Do j=1,decayChain%isotopes
        activityTime(i)%amount(j) = decayChain%amountEnd(j)
        activityTime(i)%activity(j) = decayChain%activity(j)
      End Do
    End Do
  End Subroutine CalcActivities



  Subroutine CalcActivitiesPrint(decayChain,activityTime)
! calculate activities over time
! Vars In
    Type(decayChainObj) :: decayChain
    Type(activityTimeObj), Dimension(:) :: activityTime
! Vars Private
    Integer(kind=StandardInteger) :: i, j, k, m, n, timeSteps, stringLen
    Character(Len=14) :: tempLine
    Character(Len=2048) :: printLine, headerLine, breakLine
! Time steps
    timeSteps = Size(activityTime,1)
    stringLen = 16+16*decayChain%isotopes
! Loop through steps and calculate activity
    Do i=1,timeSteps
      Write(tempLine,"(E14.7)") activityTime(i)%time
      tempLine = Trim(Adjustl(tempLine))
      Call CalcIsotopeChain(decayChain)
      printLine(1:1) = " "
      printLine(2:15) = tempLine(1:14)
      printLine(16:16) = " "
      If(i.eq.1)Then
        headerLine = "    Time/s    "
        breakLine =  "=============="
      End If
      Do j=1,decayChain%isotopes
        Write(tempLine,"(E14.7)") activityTime(i)%amount(j)
        tempLine = Trim(Adjustl(tempLine))
        m = 16*j
        printLine(m+1:m+1) = " "
        If(i.eq.1)Then
          headerLine(m+1:m+1) = " "
          breakLine(m+1:m+1) = "="
        End If
        Do k=1,14
          n = k + m + 1
          printLine(n:n) = tempLine(k:k)
          If(i.eq.1)Then
            headerLine(n:n) = decayChain%label(j)(k:k)
            breakLine(n:n) = "="
          End If
        End Do
        printLine(m+16:m+16) = " "
        If(i.eq.1)Then
          headerLine(m+16:m+16) = " "
          breakLine(m+16:m+16) = "="
        End If
      End Do
      If(i.eq.1)Then
        print *,breakLine(1:stringLen)
        print *,headerLine(1:stringLen)
        print *,breakLine(1:stringLen)
      End If
      print *,printLine(1:stringLen)
    End Do
  End Subroutine CalcActivitiesPrint







End Module activityFunctions
