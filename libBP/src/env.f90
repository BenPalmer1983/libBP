Module env
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use strings
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! ---- Variables
! ---- Functions
!  Public :: GetClockTime
! ---- Subroutines
  Public :: systemR
  Public :: progInstalled

! ---- Subroutines
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------



! Run command and save  output
  Subroutine systemR(commandIn, commandOut)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
! In
    Character(*) :: commandIn
    Character(*) :: commandOut
! Private
    Character(Len=256) :: cwdVar
    Character(Len=256) :: tempFile
    Character(Len=1024) :: commandVar
    Character(len=255) :: fileRow
    Integer(kind=StandardInteger) :: i, ios
! get working directory
    Call getcwd(cwdVar)
! Temp file name
    tempFile = trim(cwdVar)//"/tempcmdresult"
! Run command and save output plus error stream
    commandVar = trim(commandIn)//" 2> "//trim(tempFile)
    Call system(commandVar)
! Read ouput file
    commandOut = BlankString(commandOut)
    Open(UNIT=9999,FILE=trim(tempFile),status='old',action='read')
    Do i=1,1000
      Read(9999,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
      commandOut = trim(commandOut)//trim(fileRow)
    End Do
! Close file
    close(9999)
! Remove temp file
    !Call system("rm -fR "//trim(tempFile))
  End Subroutine systemR

! Check program installed
  Subroutine progInstalled(progName)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
! In
    Character(*) :: progName

    progName = ""
    !Call system("rm -fR "//trim(path))


    !Real(kind=DoubleReal) :: time,timeStart,timeEnd
    !time = time + timeEnd - timeStart
  End Subroutine progInstalled







!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------
End Module env
