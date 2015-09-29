Module general

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! ----------------------------------------
! Updated: 1st May 2014
! ----------------------------------------

! Setup Modules
  Use kinds

! force declaration of all variables
  Implicit None
! Include MPI header
! Include 'mpif.h'
! Privacy of functions/subroutines/variables
  Private
! Public subroutines
  Public :: swapArrayRows1D, swapArrayRows2D
  Public :: extractArrayColumnDP, extractArrayColumnInt
  Public :: makeDir, rmFile, rmDir, randFileName, tempFileName
  Public :: strToIntArr, strToDPArr, strToStrArr
  Public :: timeAcc
  Public :: readFile
! Public functions
  Public :: dpToString, intToString
  Public :: GetClockTime
  Public :: StrToUpper
  Public :: NumericOnly
  Public :: RemoveSpaces
  Public :: CorrectFilePath
  Public :: TrimSpaces
  Public :: BlankString
  Public :: BlankStringArray
  Public :: BlankString2DArray
  Public :: SpacesRight
  Public :: RemoveComments
  Public :: RemoveQuotes

  character( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  Contains

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

! List of Subroutines
! -------------------------------------------------------------------------
!

 


! ---------------------------------------------------------------------------------------------------


! ------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

! List of Functions
! -------------------------------------------------------------------------
!














  End Module general
