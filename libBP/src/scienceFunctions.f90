! --------------------------------------------------------------!
!Science Functions module
! scienceFunctions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 6th July 2016
! ----------------------------------------

!Module _Types
! Setup Modules
!  Use kinds
! Force declaration of all variables
!  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_ = 32
! Make private
!  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived types
!  Public ::
!End Module _Types




Module scienceFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Setup Modules
  Use mpi
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module scope parameters
! Vars:  Module scope variables
! Make private
  Private
! Public Variables and Parameters
! Public Subroutines and Functions
  Public :: F_ZBL
  Public :: F_ZblFull


!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------

  Function F_ZBL (parametersIn) RESULT (y)
! parametersIn(1) = x
! parametersIn(2) = qA
! parametersIn(3) = qB
! ZBL potential, separation x, charges qA and qB
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: parametersIn
! Vars:  Out
    Real(kind=DoubleReal) :: y
! Vars:  Private
    Real(kind=DoubleReal) :: xVal, xa, xs, exa
! Force none infinite result for 0
    If(parametersIn(1).eq.0.0D0)Then
      xVal = 0.00001D0
    Else
      xVal = parametersIn(1)
    End If
! Calculate y
    xs = 0.4683766 * (parametersIn(2)**(2.0D0/3.0D0)+parametersIn(3)**(2.0D0/3.0D0))**0.5
    xa = 1.0D0*xVal/xs
    exa = 0.1818D0*exp(-3.2D0*xa)+0.5099D0*exp(-0.9423D0*xa)+&
    0.2802D0*exp(-0.4029*xa)+0.02817*exp(-0.2016D0*xa)
    y = ((1.0D0*parametersIn(2)*parametersIn(3))/xVal)*exa
  End Function F_ZBL

  Function F_ZblFull (parametersIn) RESULT (yArray)
! parametersIn(1) = x
! parametersIn(2) = qA
! parametersIn(3) = qB
! ZBL potential, separation x, charges qA and qB
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: parametersIn
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Real(kind=DoubleReal) :: qA, qB, xVal, x, y, dy, ddy, xs
    Real(kind=DoubleReal) :: termFa, termFb, termFc, termGa, termGb, termGc
! Set vars
    x = parametersIn(1)
    qA = parametersIn(2)
    qB = parametersIn(3)
! Force none infinite result for 0
    If(x.eq.0.0D0)Then
      xVal = 0.00001D0
    Else
      xVal = x
    End If
    xs = 0.4683766 * (qA**(2.0D0/3.0D0)+qB**(2.0D0/3.0D0))**0.5
! Calculate y
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    y = termFa * termGa
    yArray(1) = y
! Calculate dy
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termFb = (1.0D0*qA*qB)*(xVal)**(-2.0D0)*(-1.0D0)                 !f'(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    termGb = (-3.2D0/xs)*0.1818D0*exp((-3.2D0/xs)*xVal)+&            !g'(x)
    (-0.9423D0/xs)*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)*0.02817*exp((-0.2016D0/xs)*xVal)
    dy = termFa*termGb+termFb*termGa
    yArray(2) = dy
! Calculate ddy
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termFb = (1.0D0*qA*qB)*(xVal)**(-2.0D0)*(-1.0D0)                        !f'(x)
    termFc = (1.0D0*qA*qB)*(xVal)**(-3.0D0)*(-1.0D0)*(-2.0D0)               !f''(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                             !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    termGb = (-3.2D0/xs)*0.1818D0*exp((-3.2D0/xs)*xVal)+&                 !g'(x)
    (-0.9423D0/xs)*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)*0.02817*exp((-0.2016D0/xs)*xVal)
    termGc = (-3.2D0/xs)**2*0.1818D0*exp((-3.2D0/xs)*xVal)+&                 !g''(x)
    (-0.9423D0/xs)**2*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)**2*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)**2*0.02817*exp((-0.2016D0/xs)*xVal)
    ddy = termFa*termGc+2*termFb*termGb+termFc*termGa
    yArray(3) = ddy
  End Function F_ZblFull


End Module scienceFunctions

!
