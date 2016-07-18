! ------------------------------------------------------------
!               GEOM: Distortion Tensors
! ------------------------------------------------------------



! Elastic Constants
!
!      |   e1    e6/2    e5/2   |
!      |  e6/2    e2     e4/2   |
!      |  e5/2   e4/2     e3    |
!
!  E(ei) = E_0 - p(V)DV + V Sum i=1,6 Sum j=1,6 (c_ij e_i e_j)/2 + O(e_i^3)
!



  Function Tensor_Homogeneous (sigma) Result (distortion)
! Homogeneous strain
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: sigma
! Vars:  Out
    Real(kind=DoubleReal), Dimension (1:3,1:3) :: distortion
! Set up Matrix
    distortion = IdentityMatrix(distortion)
    distortion(1,1) = distortion(1,1) + sigma
    distortion(2,2) = distortion(2,2) + sigma
    distortion(3,3) = distortion(3,3) + sigma
  End Function Tensor_Homogeneous


  Function Tensor_Orthorhombic (sigma) Result (distortion)
! Makes volume conserving orthorhombic strain tensor
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: sigma
! Vars:  Out
    Real(kind=DoubleReal), Dimension (1:3,1:3) :: distortion
! Set up Matrix
    distortion = IdentityMatrix(distortion)
    distortion(1,1) = distortion(1,1) + sigma
    distortion(2,2) = distortion(2,2) - sigma
    distortion(3,3) = distortion(3,3) + (sigma**2/(1.0D0-sigma**2))
  End Function Tensor_Orthorhombic


  Function Tensor_Tetragonal (sigma) Result (distortion)
! Makes volume conserving orthorhombic strain tensor
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal) :: sigma
! Vars:  Out
    Real(kind=DoubleReal), Dimension (1:3,1:3) :: distortion
! Set up Matrix
    distortion = IdentityMatrix(distortion)
    distortion(1,2) = distortion(1,2) + sigma/2.0D0
    distortion(2,1) = distortion(2,1) + sigma/2.0D0
    distortion(3,3) = distortion(3,3) + (sigma**2/(4.0D0-sigma**2))
  End Function Tensor_Tetragonal



!
