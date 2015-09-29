Module newtonGauss
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use constants
  Use matrix
! Force declaration of all variables
  Implicit None
! Public variables  
! Make private
  Private
! Public
! --variables--!
! --functions--!
  Public :: NewtonGaussOpt
! Interfaces  
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------

  Function NewtonGaussOpt(J,R) RESULT (P)!
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: J   ! Jacobian
    Real(kind=DoubleReal), Dimension(:) :: R      ! Residuals
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,2)) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: P      ! Change
! ***********
! P = (JTJ)^(-1)(-1*JTR)
! ***********
! Transpose Jacobian
    JT = TransposeMatrix(J)
    JTJ = matmul(JT,J)
    JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)
    JTR = matmul(JT,R)
    JTR = -1.0D0*JTR ! Recycle JTR var
    P = matmul(JTJ,JTR)
  End Function NewtonGaussOpt
End Module newtonGauss