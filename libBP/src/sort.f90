! --------------------------------------------------------------!
! Array sorting module
! sort
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 2nd July 2016
! ----------------------------------------

Module sort
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
  Use arrayFunctions
! Force declaration of all variables
  Implicit None
! Make declared variables private
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: sortArray
! Interfaces
  Interface sortArray
    Module Procedure sortArray_DP_1D, sortArray_DP_2D
  End Interface sortArray

!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------


! ------------------
! DP 1D
! ------------------
  Subroutine sortArray_DP_1D(arrayIn,orderIn,startRowIn,endRowIn)
! Sort array
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Real(kind=DoubleReal), Dimension (:) :: arrayIn
    Character(*), Optional :: orderIn
    Integer(kind=StandardInteger), Optional :: startRowIn
    Integer(kind=StandardInteger), Optional :: endRowIn
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Character(Len=3) :: order
    Integer(kind=StandardInteger) :: startRow
    Integer(kind=StandardInteger) :: endRow
    Logical :: keepLooping
! Optional arguments
    order = "ASC"
    startRow = 1
    endRow = Size(arrayIn,1)
    If(Present(orderIn))Then
      If(orderIn(1:1).eq."d")Then
        order = "DSC"
      End If
      If(orderIn(1:1).eq."D")Then
        order = "DSC"
      End If
    End If
    If(Present(startRowIn))Then
      startRow = startRowIn
    End If
    If(Present(endRowIn))Then
      endRow = endRowIn
    End If
! Sort array: Ascending
    If(order.eq."ASC")Then
      keepLooping = .true.
      Do While(keepLooping)
        keepLooping = .false.
        Do i=startRow,endRow-1
          If(arrayIn(i).gt.arrayIn(i+1))Then
            keepLooping = .true.
            Call swapRows(arrayIn,i,i+1)
          End If
        End Do
      End Do
    Else
      keepLooping = .true.
      Do While(keepLooping)
        keepLooping = .false.
        Do i=startRow,endRow-1
          If(arrayIn(i).lt.arrayIn(i+1))Then
            keepLooping = .true.
            Call swapRows(arrayIn,i,i+1)
          End If
        End Do
      End Do
    End If
  End Subroutine sortArray_DP_1D



! ------------------
! DP 2D
! ------------------
  Subroutine sortArray_DP_2D(arrayIn,orderIn,sortColIn,startRowIn,endRowIn)
! Sort array
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Real(kind=DoubleReal), Dimension (:,:) :: arrayIn
    Character(*), Optional :: orderIn
    Integer(kind=StandardInteger), Optional :: sortColIn
    Integer(kind=StandardInteger), Optional :: startRowIn
    Integer(kind=StandardInteger), Optional :: endRowIn
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Character(Len=3) :: order
    Integer(kind=StandardInteger) :: sortCol
    Integer(kind=StandardInteger) :: startRow
    Integer(kind=StandardInteger) :: endRow
    Logical :: keepLooping
! Optional arguments
    order = "ASC"
    sortCol = 1
    startRow = 1
    endRow = Size(arrayIn,1)
    If(Present(orderIn))Then
      If(orderIn(1:1).eq."d")Then
        order = "DSC"
      End If
      If(orderIn(1:1).eq."D")Then
        order = "DSC"
      End If
    End If
    If(Present(sortColIn))Then
      sortCol = sortColIn
      If(sortCol.gt.size(arrayIn,2))Then
        sortCol = 1
      End If
    End If
    If(Present(startRowIn))Then
      startRow = startRowIn
    End If
    If(Present(endRowIn))Then
      endRow = endRowIn
    End If
! Sort array: Ascending
    If(order.eq."ASC")Then
      keepLooping = .true.
      Do While(keepLooping)
        keepLooping = .false.
        Do i=startRow,endRow-1
          If(arrayIn(i,sortCol).gt.arrayIn(i+1,sortCol))Then
            keepLooping = .true.
            Call swapRows(arrayIn,i,i+1)
          End If
        End Do
      End Do
    Else
      keepLooping = .true.
      Do While(keepLooping)
        keepLooping = .false.
        Do i=startRow,endRow-1
          If(arrayIn(i,sortCol).lt.arrayIn(i+1,sortCol))Then
            keepLooping = .true.
            Call swapRows(arrayIn,i,i+1)
          End If
        End Do
      End Do
    End If
  End Subroutine sortArray_DP_2D


! ------------------
! INT 1D
! ------------------
  Subroutine sortArray_Int_1D(arrayIn,orderIn,startRowIn,endRowIn)
! Sort array
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Integer(kind=StandardInteger), Dimension (:) :: arrayIn
    Character(*), Optional :: orderIn
    Integer(kind=StandardInteger), Optional :: startRowIn
    Integer(kind=StandardInteger), Optional :: endRowIn
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Character(Len=3) :: order
    Integer(kind=StandardInteger) :: startRow
    Integer(kind=StandardInteger) :: endRow
    Logical :: keepLooping
! Optional arguments
    order = "ASC"
    startRow = 1
    endRow = Size(arrayIn,1)
    If(Present(orderIn))Then
      If(orderIn(1:1).eq."d")Then
        order = "DSC"
      End If
      If(orderIn(1:1).eq."D")Then
        order = "DSC"
      End If
    End If
    If(Present(startRowIn))Then
      startRow = startRowIn
    End If
    If(Present(endRowIn))Then
      endRow = endRowIn
    End If
! Sort array: Ascending
    If(order.eq."ASC")Then
      keepLooping = .true.
      Do While(keepLooping)
        keepLooping = .false.
        Do i=startRow,endRow-1
          If(arrayIn(i).gt.arrayIn(i+1))Then
            keepLooping = .true.
            Call swapRows(arrayIn,i,i+1)
          End If
        End Do
      End Do
    Else
      keepLooping = .true.
      Do While(keepLooping)
        keepLooping = .false.
        Do i=startRow,endRow-1
          If(arrayIn(i).lt.arrayIn(i+1))Then
            keepLooping = .true.
            Call swapRows(arrayIn,i,i+1)
          End If
        End Do
      End Do
    End If
  End Subroutine sortArray_Int_1D







End Module sort

























!-----------------------------------
