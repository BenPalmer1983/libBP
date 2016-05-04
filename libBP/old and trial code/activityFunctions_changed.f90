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
    Real(kind=DoubleReal), Dimension(1:100) :: amountStart = 0.0D0
    Real(kind=DoubleReal), Dimension(1:100) :: amountEnd = 0.0D0
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
    Real(kind=DoubleReal) :: dataPoints
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
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
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
    Real(kind=DoubleReal), Dimension(1:100) :: E ! Exp
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Exp
! Complete decay chain data
    Call decayChainComplete(decayChain)
! store input in shortned name arrays to make equations clearer
    t = decayChain%time
    !w = decayChain%productionRate
    print *,""
    print *,""
    Do i=1,decayChain%isotopes
      W(i) = decayChain%productionRate(i)
      L(i) = decayChain%decayConstant(i)
      E(i) = exp(-1.0D0*L(i)*t)
      B(i) = decayChain%branchFactor(i)
      N(i) = decayChain%amountStart(i)
    End Do
! Break infinities
    Call DecayBreakInfinities(L,decayChain%isotopes)
! Run analytic calculations
    print *,""
    print *,""
! Loop through isotopes
! Using L-1(1/(q+ps) = (1/p)*exp(-1*(q/p)*t) and partial fractions
    Do i=1,decayChain%isotopes
      nEnd = CalcIsotopeChainCalc(t,i,W,L,N,E,B)
      decayChain%amountEnd(i) = nEnd
    End Do
    print *,""
    print *,""
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

  Function CalcIsotopeChainCalc(t,z,W,L,N,E,B) Result (nEnd)
    Implicit None ! Force declaration of all variables
  ! Vars Private
    Integer(kind=StandardInteger) :: i, ii, termKey, j, jj, z
    Integer(kind=StandardInteger) :: termI, pfracI, bracI, nKey, lKey, lMin, startL
    Real(kind=DoubleReal) :: q, p, r, t, nEnd, pfMultiplier, termMultiplier, partVal
    Real(kind=DoubleReal), Dimension(1:100) :: W ! Production Rate
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: N ! Starting number of atoms
    Real(kind=DoubleReal), Dimension(1:100) :: E ! Exp
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Exp
! Calculate amounts from decaying isotopes
    If(L(z).gt.0.0D0)Then
!-------------------------------------------
! Unstable isotopes
!-------------------------------------------
      nEnd = 0.0D0
! Due to decaying isotopes
! Loop through terms
      Do termI=1,z
        nKey = z+1-termI
        pfMultiplier = (-1.0D0)**(termI+1)
        termMultiplier = 1.0D0
        Do i=2,termI
          termMultiplier = termMultiplier * L(z+1-i) * B(z+2-i)
        End Do
        lMin = z-termI+1
! Loop through pfrac term
        Do pfracI=1,termI
          lKey = z-termI+pfracI
          partVal = 1.0D0
          If(termI.ge.2)Then
            Do i=lMin,z
              If(i.ne.lKey)Then
                partVal = partVal * (L(lKey)-L(i))
              End If
            End Do
          End If
          p = 1.0D0 * pfMultiplier * partVal
          q = L(lKey)*p
          r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(nKey)
          nEnd = nEnd + r
        End Do
! Production term
        If(termI.eq.z)Then
          pfMultiplier = (-1.0D0)**(termI)
! 0 term
          p = 1.0D0
          Do i=1,z
            p = p * L(i)
          End Do
          q = 0.0D0
          r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*W(1)
          nEnd = nEnd + r
! Loop through pfrac terms
          lMin = 1
! Loop through pfrac term
          Do pfracI=1,z
             lKey = pfracI
             partVal = 1.0D0*L(lKey)
             If(termI.ge.2)Then
               Do i=lMin,z
                 If(i.ne.lKey)Then
                   partVal = partVal * (L(lKey)-L(i))
                 End If
               End Do
             End If
             p = 1.0D0 * pfMultiplier * partVal
             q = L(lKey)*p
             r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*W(1)
             nEnd = nEnd + r
           End Do
        End If
      End Do
    Else
!-------------------------------------------
! Stable isotopes
!-------------------------------------------
! First term (N_z(0)) stays the same (stable)
      nEnd = N(z)
! Loop through terms
      Do termI=2,z
        nKey = z+1-termI
        pfMultiplier = (-1.0D0)**(termI+1)
        termMultiplier = 1.0D0
        Do i=2,termI
          termMultiplier = termMultiplier * L(z+1-i) * B(z+2-i)
        End Do
        lMin = z-termI+1
! First term
        p = 1.0D0
        Do i=lMin,z-1
          p = p * L(i)
        End Do
        p = p * pfMultiplier
        q = 0
        r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(nKey)
        nEnd = nEnd + r
! Loop through pfrac term
        Do pfracI=2,termI
          lKey = z-termI+pfracI-1
          partVal = 1.0D0*L(lKey)
          If(termI.ge.2)Then
            Do i=lMin,z-1
              If(i.ne.lKey)Then
                partVal = partVal * (L(lKey)-L(i))
              End If
            End Do
          End If
          p = 1.0D0 * pfMultiplier * partVal
          q = L(lKey)*p
          r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(nKey)
          nEnd = nEnd + r
        End Do
      End Do
! Term due to parent production
      pfMultiplier = (-1.0D0)**(z)
      termMultiplier = 1.0D0
      Do i=2,z
        termMultiplier = termMultiplier * L(z+1-i) * B(z+2-i)
      End Do
! First pfrac 1/((abc...)s^2)
      partVal = 1.0D0
      Do i=1,z-1
        partVal = partVal * L(i)
      End Do
      partVal = 1.0D0/partVal
      r = partVal*t*termMultiplier*W(1)
      nEnd = nEnd + r
! Second pfrac (unrepeated combination)/((a^2b^2...z^2)s)
      q = 1.0D0  ! Numerator
      p = 1.0D0  ! Denominator
      q = CalcIsotopeChainNumerator(L,z-1)
      Do i=1,z-1
        p = p*L(i)*L(i)
      End Do
      r = (-1.0D0)*(q/p)*termMultiplier*W(1)
      nEnd = nEnd + r
! Remaining terms
      Do pfracI=1,z-1
        lKey = pfracI
        partVal = 1.0D0*L(lKey)*L(lKey)
        If(termI.ge.2)Then
          Do i=lMin,z-1
             If(i.ne.lKey)Then
               partVal = partVal * (L(lKey)-L(i))
             End If
           End Do
         End If
         p = 1.0D0 * pfMultiplier * partVal
         q = L(lKey)*p
         r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*W(1)
         nEnd = nEnd + r
      End Do
    End If
  End Function CalcIsotopeChainCalc


  Function CalcIsotopeChainNumerator(L,n) Result (numerator)
! Calculates numerator in isotope activity function
    Implicit None ! Force declaration of all variables
! Vars In
    Real(kind=DoubleReal), Dimension(:) :: L
    Integer(kind=StandardInteger) :: n
    Real(kind=DoubleReal) :: numerator
! Vars Private
    Integer(kind=StandardInteger) :: i, j, k, m
    Integer(kind=StandardInteger), Dimension(1:n-1) :: combinationSet
    Logical :: loopCombinations
    Real(kind=DoubleReal) :: tempVal
! init
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
            tempVal = tempVal * L(combinationSet(i))
          End Do
          numerator = numerator + tempVal
        End If
      End Do
    End If
  End Function CalcIsotopeChainNumerator







  Function CalcIsotopeChainCalcC(t,w,z,L,N,E,B) Result (nEnd)

    Implicit None ! Force declaration of all variables
! Vars Private
    Integer(kind=StandardInteger) :: i, ii, termKey, j, jj, z
    Integer(kind=StandardInteger) :: termI, pfracI, bracI, nKey, lKey, lMin, startL
    Real(kind=DoubleReal) :: q, p, r, w, t, nEnd, pfMultiplier, termMultiplier, partVal
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: N ! Starting number of atoms
    Real(kind=DoubleReal), Dimension(1:100) :: E ! Exp
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Exp
! Calculate amounts from decaying isotopes
    nEnd = 0.0D0

    If(z.eq.1)Then
      nEnd = 0.0D0
! Term 1 of 1
! pfrac 1 of 1
      termMultiplier = 1.0D0
      pfMultiplier = 1.0D0
      p = 1.0D0 * pfMultiplier
      q = L(1)
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(1)
      nEnd = nEnd + r
      print *,z,nEnd
    End If

    If(z.eq.2)Then
      nEnd = 0.0D0
    ! Term 1 of 2
      termMultiplier = 1.0D0
      pfMultiplier = 1.0D0
    ! pfrac 1 of 1
      p = 1.0D0 * pfMultiplier
      q = L(2)
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(2)
      nEnd = nEnd + r


    ! Term 2 of 2
      termMultiplier = 1.0D0*L(1)
      pfMultiplier = -1.0D0
    ! pfrac 1 of 2
      p = 1.0D0 * pfMultiplier * (L(1)-L(2))
      q = L(1)*p
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(1)
      nEnd = nEnd + r
    ! pfrac 2 of 2
      p = 1.0D0 * pfMultiplier * (L(2)-L(1))
      q = L(2)*p
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(1)
      nEnd = nEnd + r
      print *,z,nEnd
    End If

    If(z.eq.3)Then
      nEnd = 0.0D0
    ! Term 1 of 3
      termMultiplier = 1.0D0
      pfMultiplier = 1.0D0
    ! pfrac 1 of 1
      p = 1.0D0 * pfMultiplier
      q = L(3)
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(3)
      nEnd = nEnd + r


    ! Term 2 of 3
      termMultiplier = 1.0D0*L(2)
      pfMultiplier = -1.0D0
    ! pfrac 1 of 2
      p = 1.0D0 * pfMultiplier * (L(2)-L(3))
      q = L(2)*p
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(2)
      nEnd = nEnd + r
    ! pfrac 2 of 2
      p = 1.0D0 * pfMultiplier * (L(3)-L(2))
      q = L(3)*p
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(2)
      nEnd = nEnd + r


    ! Term 3 of 3
      termMultiplier = 1.0D0*L(1)*L(2)
      pfMultiplier = 1.0D0
    ! pfrac 1 of 3
      p = 1.0D0 * pfMultiplier * (L(1)-L(2)) * (L(1)-L(3))
      q = L(1)*p
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(1)
      nEnd = nEnd + r
    ! pfrac 2 of 3
      p = 1.0D0 * pfMultiplier * (L(2)-L(1)) * (L(2)-L(3))
      q = L(2)*p
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(1)
      nEnd = nEnd + r
    ! pfrac 3 of 3
      p = 1.0D0 * pfMultiplier * (L(3)-L(1)) * (L(3)-L(2))
      q = L(3)*p
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(1)
      nEnd = nEnd + r

      print *,z,nEnd
    End If

    print *,""
    print *,"------------------------"
    print *,z,""
    print *,"------------------------"
    nEnd = 0.0D0
    If(z.le.3)Then
    Do termI=1,z   ! Loop through terms
      print *,"term ",termI
      print "(A5,I2,A1)","   N(",nKey,")"
      print *,pfMultiplier,"  ",termMultiplier
      nKey = z+1-termI
      pfMultiplier = (-1.0D0)**(termI+1)
      termMultiplier = 1.0D0
      Do i=2,termI
        termMultiplier = termMultiplier * L(z+1-i)
      End Do
      lMin = z-termI+1
      Do pfracI=1,termI ! Loop through pfrac term
        lKey = z-termI+pfracI
        partVal = 1.0D0
        print *,"    pfrac ",pfracI,"  "
        print *,"    lKey  ",lKey
        print *,"    lMin  ",lMin
        If(termI.ge.2)Then
          Do i=lMin,z
            If(i.ne.lKey)Then
              print "(A10,I2,A4,I2,A1)","        (L",lKey," - L",i,")"
              partVal = partVal * (L(lKey)-L(i))
            End If
          End Do
        End If
        p = 1.0D0 * pfMultiplier * partVal
        q = L(lKey)*p
        r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(nKey)
        nEnd = nEnd + r
      End Do
    End Do
    print *,z,"___",nEnd
    End If
    print *,""
    print *,""

!termI, fracI








    If(z.eq.99)Then
    Do i=1,z  ! Loop through terms for decaying isotope terms
      termMultiplier = 1.0D0
      Do j=1,i-1
        termMultiplier = termMultiplier*B(j+1)*L(j)
      End Do
      termKey = z + 1 - i
      pfMultiplier = (-1)**(i+1)
      Do j=1,i  ! Loop through partial fraction terms
        p = 1.0D0
        partVal = 1.0D0
        Do jj=1,j  ! looping through bracket (a-b) etc in denominator
          If(jj.ne.j)Then
            partVal = L(j) - L(jj)
            p = p * partVal
          End If
        End Do
        p = p * pfMultiplier
        q = p * L(j)
        r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(termKey)
        nEnd = nEnd + r
        !print *,i,j,q,p,r,nEnd,termMultiplier
      End Do
    End Do
    print *,""
  End If


  End Function CalcIsotopeChainCalcC


  Function CalcIsotopeChainCalcB(t,w,z,L,N,E,B) Result (nEnd)
! t time
! w parent production rate
! z
    Implicit None ! Force declaration of all variables
  ! Vars Private
    Integer(kind=StandardInteger) :: i, j, ii, jj, z
    Real(kind=DoubleReal) :: q, p, r, w, t, nEnd, pfMultiplier, termMultiplier, partVal
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: N ! Starting number of atoms
    Real(kind=DoubleReal), Dimension(1:100) :: E ! Exp
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Exp
! Calculate amounts from decaying isotopes
    nEnd = 0.0D0
! --------------------
! Unstable Isotopes
! --------------------
!
    print *,z,L(z)
    If(L(z).gt.0.0D0)Then
! Unstable isotopes
      Do i=1,z
        pfMultiplier = (-1)**(i-1)    ! from partial fraction expansion
        termMultiplier = 1.0D0      ! from branching and lambda constants
        Do j=1,i-1
          termMultiplier = termMultiplier*B(j+1)*L(j)
        End Do
        If(i.eq.1)Then
          p = 1.0D0
          q = L(z-i+1)
          nEnd = nEnd + (1.0D0/p)*exp(-1.0D0*(q/p)*t)*N(z-i+1)
        Else  ! 2nd term and above
          Do j=1,i  ! loop through pFrac terms
            p = 1.0D0
            nEnd = nEnd + termMultiplier*N(z-i+1)
            Do jj=1,i  ! loop through lambda constants
              If(j.ne.jj)Then
                partVal = (L(j)-L(jj))
                p = p * partVal
              End If
            End Do
            p = p * pfMultiplier
            q = p * L(j)
            r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(z-i+1)
            nEnd = nEnd + r
          End Do
        End If
      End Do

    print *,"Isotope Z"
! Calculate contribution from production of parent (and extra atoms to decay)
    termMultiplier = 1.0D0      ! from branching and lambda constants
    Do i=1,z-1
      termMultiplier = termMultiplier*B(i+1)*L(i)
    End Do
! First pfrac term
    p = 1.0D0
    Do i=1,z
      p = p * L(i)
      print *,i
    End Do
    r = (1.0D0/p)*termMultiplier*w
    !nEnd = nEnd + r
! Remaining pfrac terms
    pfMultiplier = (-1)**(z)
    Do i=1,z  ! loop through pFrac terms
      p = 1.0D0
      Do jj=1,i  ! loop through lambda constants
        If(i.ne.jj)Then
          partVal = (L(i)-L(jj))
          p = p * partVal
          print *,i,jj,partVal
        End If
      End Do
      p = L(i) * p * pfMultiplier
      q = p * L(i)
      r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*w
      !nEnd = nEnd + r
      End Do



    Else
! Stable isotopes

    End If


! Calculates production terms

  End Function CalcIsotopeChainCalcB


  Function CalcIsotopeChainCalcA(t,w,z,L,N,E,B) Result (nEnd)
    Implicit None ! Force declaration of all variables
  ! Vars Private
    Integer(kind=StandardInteger) :: i, j, jj, z
    Real(kind=DoubleReal) :: q, p, r, w, t, nEnd, multiplier, termMultiplier, partVal
    Real(kind=DoubleReal), Dimension(1:100) :: L ! Lambda
    Real(kind=DoubleReal), Dimension(1:100) :: N ! Starting number of atoms
    Real(kind=DoubleReal), Dimension(1:100) :: E ! Exp
    Real(kind=DoubleReal), Dimension(1:100) :: B ! Exp
! Calculate amounts from decaying isotopes
    nEnd = 0.0D0
    Do i=1,z
      If(i.eq.1)Then
        p = 1.0D0
        q = L(z-i+1)
        nEnd = nEnd + (1.0D0/p)*exp(-1.0D0*(q/p)*t)*N(z-i+1)
      Else  ! 2nd term and above
        multiplier = (-1)**(i-1)    ! from partial fraction expansion
        termMultiplier = 1.0D0      ! from branching and lambda constants
        Do j=1,i-1
          termMultiplier = termMultiplier*B(j+1)*L(j)
        End Do
        Do j=1,i  ! loop through pFrac terms
          p = 1.0D0
          nEnd = nEnd + termMultiplier*N(z-i+1)
          Do jj=1,i  ! loop through lambda constants
            If(j.ne.jj)Then
              partVal = (L(j)-L(jj))
              p = p * partVal
            End If
          End Do
          p = p * multiplier
          q = p * L(j)
          r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*N(z-i+1)
          nEnd = nEnd + r
        End Do
      End If
    End Do
! Calculate contribution from production of parent (and extra atoms to decay)
    termMultiplier = 1.0D0      ! from branching and lambda constants
    Do i=1,z-1
      termMultiplier = termMultiplier*B(i+1)*L(i)
    End Do
! First term
    p = 1.0D0
    Do i=1,z
      p = p * L(i)
    End Do
    r = (1.0D0/p)*termMultiplier*w
    nEnd = nEnd + r
! Remaining pfrac terms
    multiplier = (-1)**(z)
    Do i=1,z  ! loop through pFrac terms
      p = 0.0D0




                If(L(j).le.0.0D0)Then
                  ! Stable, this fraction would cause infinity
                Else
                  Do jj=1,i  ! loop through lambda constants
                    If(jj.eq.1)Then
                      p = L(j)
                    End If
                    If(j.ne.jj)Then
                      partVal = (L(j)-L(jj))
                      p = p * partVal
                    End If
                  End Do
                  p = p * multiplier
                  q = p * L(j)
                  r = (1.0D0/p)*exp(-1.0D0*(q/p)*t)*termMultiplier*w
                  nEnd = nEnd + r
                End If
              End Do




! Calculates production terms

  End Function CalcIsotopeChainCalcA






End Module activityFunctions
