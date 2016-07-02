! --------------------------------------------------------------!
! File Output Module
! fileOutTypes, fileOut
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calculates the forces between atoms
! Calculates the energies of atoms/total energy
! Calculates stresses
!
! ----------------------------------------
! Updated: 21st May 2016
! ----------------------------------------


Module fileOutTypes
! Setup Modules
  Use kinds
  Type :: fileObj
    Character(Len=64) :: filePath = ""
    Character(Len=32) :: fileName = ""
    Character(Len=4) :: fileExtension = ""
    Integer(kind=StandardInteger) :: lineCount
    Character(Len=96), Dimension(1:10000) :: data
  End Type fileObj

End Module fileOutTypes


Module fileOut

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Read user input file

! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------

! Setup Modules
  Use kinds
  Use strings
  Use general
  Use mpiSubsTypes
  Use mpiSubs
  Use units
  Use fileOutTypes
! Force declaration of all variables
  Implicit None
! Define variables
  Type(fileObj) :: mainFile
! Privacy of variables/functions/subroutines
  Private
! Set Public Variables
  Public :: mainFile
! Public Subroutines - Print
  Public :: outputFile

! Interfaces
!  Interface subName
!    Module Procedure subA, subB
!  End Interface subName

  Contains

! --------------------------------------------------------------------------
!     Print Basic
! --------------------------------------------------------------------------

  Subroutine outputFile(fileOut_In)
! Add header row (to head each column)
    Implicit None   ! Force declaration of all variables
! Vars:  In/out
    Type(fileObj), Optional :: fileOut_In
! Vars:  Private
    Type(fileObj) :: fileOut
! Optional argument
    If(Present(fileOut_In))Then
      fileOut = fileOut_In
    Else
      fileOut = mainFile
    End If
  End Subroutine outputFile




End Module fileOut
