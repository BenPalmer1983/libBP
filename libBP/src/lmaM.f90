Module lmaM
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
  Use linearAlgebra
  Use calcFunctions
! Force declaration of all variables
  Implicit None
! Public variables  
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: LMA
  Public :: LMA_Polynomial
  Public :: LMA_BirchMurn
  Public :: LMA_Exp
  Public :: LMA_ExpDens
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
  
  Function LMA(points, calcFunction, parametersIn, weightingIn, limitsLowerIn, limitsUpperIn) &
  RESULT (parametersOut)
! --------------------------------------------------------------------------------------------------  
! Functions by Ben Palmer University of Birmingham 2015
!
! LMA algorithm
! Levenberg, K.: 1944, Quarterly of Applied Mathematics 2, 164:
!    Dampened Newton Gauss Algorithm
! Marquardt, D.: 1963, Journal of the Society for Industrial and Applied Mathematics 11(2), 431:
!    Introduction of Lambda (Marquardt Parameter)
! R. Fletcher 1971, A modified Marquardt Subroutine for Non-linear Least Squares
!    Starting Lambda and lower lambda cutoff (where it is replaced by 0)
! M. Transtrum, J. Sethna 2012, Improvements to the Levenberg-Marquardt algorithm for nonlinear
! least-squares minimization:
!    Delayed gratification scheme for varying lambda
! Jens Jessen-Hansen 2011, Levenberg-Marquardts algorithm Project in Numerical Methods:
!    Independent/Diagonal covariance matrix for weighting 
! Parameter limits implemented 
! Finite difference method used to calculate J
! ----------------------------
! Requires kinds to be specified, if not already:
! Integer, Parameter :: StandardInteger = Selected_Int_Kind(8) 
! Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307) 
! -------------------------------------------------------------------------------------------------- 
! points - array of input x,y data points the function is being fit to
! calcFunction - the name of the function being called to calculate f(x)
! parametersIn - the starting parameters that will be used and improved for calcFunction
! limitsLowerIn - array of lower limits (if required)
! limitsUpperIn - array of upper limits (if required)
! set all limitsUpperIn lower than limitsLowerIn to ignore limits on all parameters
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), External :: calcFunction
    Real(kind=DoubleReal), Dimension(:) :: parametersIn
    Logical, Optional :: weightingIn
    Real(kind=DoubleReal), Dimension(:), Optional :: limitsLowerIn, limitsUpperIn
! Out
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: parametersOut
! Private variables
    Logical :: weighting
    Integer(kind=StandardInteger) :: i, n, k
    Integer(kind=StandardInteger) :: nPoints, nParams
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: parameters, parametersT, parametersL
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: limitsLower, limitsUpper
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:size(parametersIn,1)) :: J
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: R
    Real(kind=DoubleReal) :: x, y, fx, fxP, dVal, dValFactor
    Real(kind=DoubleReal) :: lambda, lambdaCutoff
    Real(kind=DoubleReal) :: rss, rssTrial
    Real(kind=DoubleReal) :: convT, convV
    Logical :: converged
! Optional arguments
    weighting = .true.
    limitsLower = 1.0D0   
    limitsUpper = -1.0D0   
    If(Present(limitsLowerIn))Then
      limitsLower = limitsLowerIn
    End If
    If(Present(limitsUpperIn))Then
      limitsUpper = limitsUpperIn
    End If   
    If(Present(weightingIn))Then
      weighting = weightingIn
    End If     
! Init vars
    parameters = parametersIn
    nPoints = size(points,1)
    nParams = size(parametersIn,1)
    dValFactor = 1.0D-4
    converged = .false.
    convT = 1.0D-4
! ----------------
! Start LMA loop
! ----------------
    Do n=1,20  ! maximum 20 loops
      rss = 0.0D0
      Do i=1,nPoints
        x = points(i,1)
        y = points(i,2)
        fx = LMA_FWrapper(calcFunction,x,parameters,size(parameters),limitsLower,limitsUpper)
        R(i) = fx - y
        rss = rss + R(i)**2
        Do k=1,nParams
          dVal = abs(parametersT(k) * dValFactor)
          If(dVal.eq.0.0D0)Then
            dVal = 1.0D-5
          End If
          parametersT = parameters
          parametersT(k) = parameters(k) + dVal
          fxP = LMA_FWrapper(calcFunction,x,parametersT,size(parametersT),limitsLower,limitsUpper)
          J(i,k) = (fxP-fx)/dVal
        End Do
      End Do
! Starting lambda
      If(n.eq.1)Then
        lambdaCutoff = LMA_Lambda(J)
        lambda = lambdaCutoff
      End If
! Start loop - keep increasing lambda if not an improvement
      Do k=1,100 ! max 100
! Store last loop values
        parametersL = parameters
! calculate new parameters
        parameters = LMA_Calc(J,R,lambda,parameters,weighting)
! Check for parameter limits        
        Do i=1,size(parameters,1)
          If(limitsLower(i).lt.limitsUpper(i))Then
            If(parameters(i).lt.limitsLower(i))Then
              parameters(i) = limitsLower(i)
            End If
            If(parameters(i).gt.limitsUpper(i))Then
              parameters(i) = limitsUpper(i)
            End If
          End If
        End Do
! Calc RSS
        rssTrial = LMA_FunctionRSS(points, calcFunction, parameters)
! Breakout if NaN
        If(Isnan(rssTrial))Then
          parameters = parametersL 
          Exit          
        End If
! Check convergence
        convV = abs((rssTrial-rss)/rss)
        If(convV.lt.convT)Then
          converged = .true.
        End If
! Delayed gratification scheme - 1.5*lambda or 0.2*lambda
        If(converged)Then
          Exit
        Else
          If(rssTrial.gt.rss)Then  ! If worse
            If(lambda.lt.lambdaCutoff)Then
              lambda = 2.0D0*lambdaCutoff
            End If
            lambda = lambda * 1.5D0
            parameters = parametersL
          Else  ! If better
            lambda = lambda * 0.2D0
            Exit
          End If
        End If
      End Do
      If(converged)Then
        Exit
      End If
    End Do
    parametersOut = parameters
  End Function LMA

  Function LMA_Calc(J,R,lambda,parametersIn,weighting) RESULT (parametersOut)
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: J                                     ! Jacobian
    Real(kind=DoubleReal), Dimension(:) :: R                                       ! Residuals
    Real(kind=DoubleReal), Dimension(:) :: parametersIn
    Real(kind=DoubleReal) :: lambda
    Logical :: weighting
! Out
    Real(kind=DoubleReal), Dimension(1:size(parametersIn,1)) :: parametersOut
! Private variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:size(R),1:size(R)) :: W                ! 
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,1)) :: JT               ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,1)) :: JTW              ! 
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,2)) :: JTWJ
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,2)) :: JTWJ_Diag
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: JTWR 
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,2)) :: A  ! Left side
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: B              ! Right side
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: P                             ! Change
    Real(kind=DoubleReal) :: lambdaCutoff
! ***********
! P = (JTJ+Lambda*diag(JTJ))^(-1)(-1*JTR)
! ***********
! Check Marquardt parameter cutoff
    lambdaCutoff = LMA_Lambda(J)
! Check Marquardt parameter cutoff
    If(lambda.lt.lambdaCutoff)Then
      lambda = 0.0D0
    End If
! Prep covariance matrix    
    W = 0.0D0
    Do i=1,size(R)
      If(weighting)Then
        W(i,i) = 1.0D0/R(i)**2
      Else
        W(i,i) = 1.0D0
      End If
    End Do
! Matrix calculations
    JT = TransposeMatrix(J)
    JTW = matmul(JT,W)
    JTWJ = matmul(JTW,J)
    JTWR = matmul(JTW,R)
    JTWJ_Diag = lambda*DiagMatrix(JTWJ) ! Dampening Matrix
    A = MatAdd(JTWJ,JTWJ_Diag)
    B = -1.0D0*JTWR
    P = SolveLinearSet(A, B)
! Update parameters
    Do i=1,size(P)
      parametersOut(i) = parametersIn(i) + P(i)
    End Do
  End Function LMA_Calc

  Function LMA_Lambda(J) RESULT (lambda)
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: J                                     ! Jacobian
! Out
    Real(kind=DoubleReal) :: lambda
! Private variables
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,1)) :: JT                ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,2)) :: JTJ, JTJ_Inv
    Real(kind=DoubleReal) :: traceJTJ
! Calculate starting lambda/lambda cutoff
    JT = TransposeMatrix(J)
    JTJ = matmul(JT,J)
    JTJ_Inv = InvertMatrix(JTJ)
! Check Marquardt parameter cutoff
    traceJTJ = Trace(JTJ_Inv)
    lambda = traceJTJ**(-1)
  End Function LMA_Lambda

  Function LMA_FWrapper(calcFunction,x,parameters,pSize,limitsLower,limitsUpper) RESULT (fx)  
! Wrapper function around the function being called
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), External :: calcFunction
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters
    Real(kind=DoubleReal), Dimension(1:pSize) :: limitsLower, limitsUpper
    Real(kind=DoubleReal), Intent(IN) :: x
! Out    
    Real(kind=DoubleReal) :: fx
! Private    
    Integer(kind=StandardInteger) :: i
! Calculation    
    fx = calcFunction(x,parameters,pSize)  
    Do i=1,pSize
      If(limitsLower(i).lt.limitsUpper(i))Then  ! Limits have been set
        If(parameters(i).lt.limitsLower(i).or.parameters(i).gt.limitsUpper(i))Then
          If(fx.eq.0.0D0)Then
            fx = 1.0D0
          End If
          fx = fx * 1.0D5
        End If
      End If
    End Do
  End Function LMA_FWrapper    

  Function LMA_FunctionRSS(points, calcFunction, parameters) RESULT (rss)
! RSS between function and calculated points
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), External :: calcFunction
    Real(kind=DoubleReal), Dimension(:) :: parameters
! Out
    Real(kind=DoubleReal) :: rss
! Private variables
    Real(kind=DoubleReal) :: x, y, fx
    Integer(kind=StandardInteger) :: i
! RSS
    rss = 0.0D0
    Do i=1,size(points,1)
      x = points(i,1)
      y = points(i,2)
      fx = calcFunction(x,parameters,size(parameters,1))
      rss = rss + (y - fx)**2
    End Do
  End Function LMA_FunctionRSS
  
!---------------------------------------------------------------------------------------------------------------------------------------
! Functions for LMA 
!---------------------------------------------------------------------------------------------------------------------------------------

  Function LMA_Polynomial(x,parameters,pSize) RESULT (y)
    Implicit None  !Force declaration of all variables
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters
    Real(kind=DoubleReal), Intent(IN) :: x
    Real(kind=DoubleReal) :: y
    y = CalcPolynomial(parameters, x)
  End Function LMA_Polynomial

  Function LMA_BirchMurn(x,parameters,pSize) RESULT (y)
    Implicit None  !Force declaration of all variables
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters
    Real(kind=DoubleReal), Intent(IN) :: x
    Real(kind=DoubleReal) :: y
    y = BirchMurnCalc(x, parameters)
  End Function LMA_BirchMurn

  Function LMA_Exp(x,parameters,pSize) RESULT (y)
    Implicit None  !Force declaration of all variables
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters
    Real(kind=DoubleReal), Intent(IN) :: x
    Real(kind=DoubleReal) :: y
    y = ExpCalc(x,parameters)    
  End Function LMA_Exp

  Function LMA_ExpDens(x,parameters,pSize) RESULT (y)
! Density function  p(r) = ar^2 exp(br^2) + cr^2 exp(dr^2)  
    Implicit None  !Force declaration of all variables
    Integer(kind=StandardInteger), Intent(IN) ::  pSize
    Real(kind=DoubleReal), Dimension(1:pSize), Intent(IN) :: parameters
    Real(kind=DoubleReal), Intent(IN) :: x
    Real(kind=DoubleReal) :: y
    y = parameters(1)*x**2*exp(parameters(2)*x**2)+&
    parameters(3)*x**2*exp(parameters(4)*x**2)    
  End Function LMA_ExpDens


   

End Module lmaM