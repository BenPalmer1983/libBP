Module rngFunctions
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use rng
  Use arrayFunctions
!  Use constants
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! ---- Variables
  Public :: IntegerList
! ---- Subroutines
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

  Function IntegerList(listStart,listEnd,shuffles) RESULT (list)
! Array filled with integers, possibly shuffled
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: listStart, listEnd, listSize, shuffles, rowA, rowB, shuffleCount
    Integer(kind=StandardInteger), Dimension(1:(listEnd-listStart+1)) :: list
! Initialise variables
    shuffleCount = 0
    listSize = listEnd-listStart+1
! Make list
    Do i=1,listSize
      list(i) = listStart+i-1
    End Do
! Shuffle list
    If(shuffles.gt.0)Then
      Do While(shuffleCount.lt.shuffles)
        rowA = RandomInteger(1,listSize)
        rowB = RandomInteger(1,listSize)
        If(rowA.ne.rowB)Then
          Call swapRows(list,rowA,rowB)
          shuffleCount = shuffleCount + 1
        End If
      End Do
    End If
  End Function IntegerList

End Module rngFunctions
