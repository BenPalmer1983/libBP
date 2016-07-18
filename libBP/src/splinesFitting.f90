Module splinesFitting
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Fitting module: BM, exp etc
! Splines: Poly splines, exp(poly) splines
! --------------------------------------------------------------!
  Use kinds
  Use strings
  Use constants
  Use matrix
  Use linearAlgebra
  Use rng
  Use calcFunctions
  Use fitting
  Use regression
  Use interpolation
  Use lmaM
! Force declaration of all variables
  Implicit None
! Public variables
  Character(Len=128), Dimension(1:20) :: fittingReport
  Real(kind=DoubleReal), Dimension(1:1000,1:10) :: splineCoeffs
! Make private
  Private
! Public
! --variables--!
  Public :: fittingReport
  Public :: splineCoeffs
! --- Functions - Fitting
  Public :: BirchMurnFit
  Public :: ExpFit
  Public :: SingleDecayFit
  Public :: DoubleDecayFit
  Public :: TripleDecayFit
  Public :: FittingPoints
  Public :: FitEmbeddingA
  Public :: FitEmbeddingB
  Public :: FitEmbeddingC
  Public :: FitDensity
  Public :: BestFitLine
  Public :: BestFitExp
  Public :: BestFitExpSq
  Public :: MorseFit
  Public :: MorseFit_Extended
  Public :: LJFit
! --- Functions - Spline
  Public :: SplineAB
  Public :: SplineExpThird
  Public :: SplineExpFifth
  Public :: SplineNodes
  Public :: SplineComplete
  Public :: VaryNode
  Public :: SplinePoints
! --- Subroutines
  Public :: CompleteNodeData
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------
! Fitting
! ----------------------------------------------------------------------------------

  Include "splinesFitting.Fitting.f90"

! ----------------------------------------------------------------------------------
! Splines
! ----------------------------------------------------------------------------------

  Include "splinesFitting.Splines.f90"

! ----------------------------------------------------------------------------------
! Module Subroutines
!
!
!
!
! ----------------------------------------------------------------------------------

Subroutine CompleteNodeData(splineNodes, startIn, endIn)
! Complete y'(x) and y''(x) values for a y(x) spline
! At least 4 data points required, otherwise exits
! At least 4 "columns" in splineNodes x, y(x), y'(x). y''(x)
  Implicit None  ! Force declaration of all variables
! Declare private variables
  Real(kind=DoubleReal), Dimension(:,:) :: splineNodes
  Integer(kind=StandardInteger) :: node, nodes, i, j
  Integer(kind=StandardInteger), optional :: startIn, endIn
  Integer(kind=StandardInteger) :: startNode, endNode, tempNode
  Integer(kind=StandardInteger) :: interpStart, interpEnd
  Real(kind=DoubleReal), Dimension(1:4,1:2) :: interpNodes
! optional arguments
  startNode = 1
  endNode = size(splineNodes,1)
  If(present(startIn))Then
    If(startIn.ge.1)Then
      startNode = startIn
    End If
  End If
  If(present(endIn))Then
    If(endIn.le.endNode)Then
      endNode = endIn
    End If
  End If
! Swap start/end if wrong way around
  If(startNode.gt.endNode)Then
    tempNode = startNode
    startNode = endNode
    endNode = tempNode
  End If
! Exit subroutine if too few points
  nodes = endNode-startNode+1
  If(nodes.lt.4)Then
    return
  End If
! Interp at each node
  Do node=startNode,endNode
! Set start-end nodes used for interpolation
    interpStart = node-2
    interpEnd = interpStart+3
! Adjust these nodes so they are within the bounds of the data
    If(interpStart.lt.startNode)Then
      interpStart = startNode
      interpEnd = interpStart + 3
    End If
    If(interpEnd.gt.endNode)Then
      interpEnd = endNode
      interpStart = interpEnd - 3
    End If
! make interpolation array
    j = 0
    Do i=interpStart,interpEnd
      j = j + 1
      interpNodes(j,1) = splineNodes(i,1)
      interpNodes(j,2) = splineNodes(i,2)
    End Do
! print *,node,splineNodes(node,1),splineNodes(node,2),splineNodes(node,3),splineNodes(node,4)
! lagrange interpolation
    splineNodes(node,2) = InterpLagrange(splineNodes(node,1), interpNodes, 0)
    splineNodes(node,3) = InterpLagrange(splineNodes(node,1), interpNodes, 1)
    splineNodes(node,4) = InterpLagrange(splineNodes(node,1), interpNodes, 2)
! print *,node,splineNodes(node,1),splineNodes(node,2),splineNodes(node,3),splineNodes(node,4)
! print *,""
  End Do
End Subroutine CompleteNodeData




End Module splinesFitting
