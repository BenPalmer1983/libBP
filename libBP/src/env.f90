! --------------------------------------------------------------!
! Array sorting module
! sort
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 4th July 2016
! ----------------------------------------

Module envTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter :: p_ = 32
! Make private
  Private
! Public Variables and Parameters
!  Public :: p_
! Public derived types
  Public :: envType
! Define types
  Type :: envType
    Character(Len=32) :: user
    Character(Len=64) :: cwd
    Character(Len=64) :: home
    !Character(Len=16), Dimension(1:p_potentials) :: atomLabel_B
    !Integer(kind=StandardInteger), Dimension(1:p_potentials) :: atomID_A
  End Type envType
End Module envTypes

Module env
  Use kinds
  Use strings
  Use envTypes
! Force declaration of all variables
  Implicit None
! Public variables
  Type(envType) :: envVars
! Make private
  Private
! Public
! ---- Variables
  Public :: envVars
! ---- Functions
!  Public ::
! ---- Subroutines
  Public :: loadVars
  Public :: systemR
  Public :: progInstalled

! ---- Subroutines
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------


  Subroutine loadVars()
! Run command and save  output
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
! Vars:  Private
    Character(Len=64) :: cmdLine
    Character(Len=128) :: tempLine
!-------------------------
! Store vars
!-------------------------
! current working directory
    Call getcwd(tempLine)
    envVars%cwd = trim(tempLine)
! user
    cmdLine = "whoami"
    Call systemR(cmdLine,tempLine)
    envVars%user = trim(tempLine)
! home
    envVars%home = "/home/"//trim(envVars%user)
  End Subroutine loadVars



  Subroutine systemR(commandIn, commandOut)
! Run command and save  output
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
    commandVar = trim(commandIn)//" 1> "//trim(tempFile)
    Call Execute_Command_Line(commandVar)
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
    Call Execute_Command_Line("rm "//trim(tempFile))
  End Subroutine systemR


  Subroutine progInstalled(progName)
! Check program installed
    Implicit None   ! Force declaration of all variables
! Declare private variables
! In
    Character(*) :: progName
    progName = ""
    !Call system("rm -fR "//trim(path))
    !Real(kind=DoubleReal) :: time,timeStart,timeEnd
    !time = time + timeEnd - timeStart
  End Subroutine progInstalled





  Function GetHomeDirText (input) RESULT (output)
! Get home directory from result text between : and :
    Implicit None  !Force declaration of all variables
! Vars:  In
    Character(*) :: input
! Vars:  Out
    Character(Len(input)) :: output
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Logical :: flag
! Read input text
    flag = .false.
    j = 0
    Do i=1,Len(input)
      If(input(i:i).eq.":")Then
        If(flag)Then
          Exit
        Else
          flag = .true.
        End If
      Else
        If(flag)Then
          j = j + 1
          output(j:j) = input(i:i)
        End If
      End If
    End Do
    output = TrimStr(output)
  End Function GetHomeDirText

!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------
End Module env
