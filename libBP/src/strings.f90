Module strings
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
! Force declaration of all variables
  Implicit None
! Public variables  
  Integer(kind=LongInteger) :: randomLCG_n_strings=0
  Integer(kind=LongInteger) :: randomLCG_xn_strings
  Integer(kind=LongInteger) :: randomLCG_R_n_strings=0
  Integer(kind=LongInteger) :: randomLCG_R_xn_strings
! Make private
  Private
! Public
! ---- Variables
  Public :: randomLCG_n_strings
  Public :: randomLCG_xn_strings
  Public :: randomLCG_R_n_strings
  Public :: randomLCG_R_xn_strings
! ---- Functions
  Public :: StrToUpper
  Public :: StrToLower
  Public :: NumericOnly
  Public :: RemoveSpaces
  Public :: TrimSpaces
  Public :: BlankString
  Public :: BlankStringArray
  Public :: BlankString2DArray
  Public :: SpacesRight
  Public :: RemoveComments
  Public :: RemoveQuotes
  Public :: IntToStr
  Public :: DpToStr
  Public :: StrToInt
  Public :: StrToDp
  Public :: RandName
  Public :: TempFileName
  Public :: CleanString
! ---- Subroutines  
  Public :: randCharacter
  
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
 
! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

  Function StrToUpper (input) RESULT (output)
! Returns factorial of input
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(IN) :: input
    Character(LEN(input)) :: output
! Private
    Integer(kind=StandardInteger) :: i, n
! Loop through characters
    output = "   "
    Do i=1,Len(input)
      n = Iachar(input(i:i))
      If(n.ge.97.and.n.le.122)Then
        n = n - 32
      End If
      output(i:i) = char(n)
    End Do
  End Function StrToUpper
!---------------------------------------------------------------------------------------------------------------------------------------
  Function StrToLower (input) RESULT (output)
! Returns factorial of input
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(IN) :: input
    Character(Len(input)) :: output
! Private
    Integer(kind=StandardInteger) :: i, n
! Loop through characters
    Do i=1,Len(input)
      n = Iachar(input(i:i))
      If(n.ge.65.and.n.le.90)Then
        n = n + 32
      End If
      output(i:i) = char(n)
    End Do
  End Function StrToLower
!---------------------------------------------------------------------------------------------------------------------------------------
  Function NumericOnly (input) RESULT (output)
! Returns factorial of input
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(IN) :: input
    Character(Len(input)) :: outputTemp
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i, n
! Remove characters
    outputTemp = input
    Do i = 1, Len( outputTemp )
      output( i:i ) = " "
    End Do
    n = 0
    Do i=1,Len(outputTemp)
      If(outputTemp( i:i ).eq.".".or.(iachar(outputTemp( i:i )).ge.48.and.iachar(outputTemp( i:i )).le.57))Then
        n = n + 1
        output( n:n ) = outputTemp( i:i )
      End If
    End Do
  End Function NumericOnly
!---------------------------------------------------------------------------------------------------------------------------------------
  Function RemoveSpaces (input) RESULT (output)
! Returns factorial of input
    Implicit None  !Force declaration of all variables
! Declare variables    
    Character(*), Intent(IN) :: input
    Character(Len(input)) :: outputTemp
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i, j
! Remove spaces
    outputTemp = input
! Blank output
    Do i=1,Len(outputTemp)
      output(i:i) = " "
    End Do
! transfer outputtemp to output without spaces
    j = 0
    Do i=1,Len(outputTemp)
      If(outputTemp( i:i ).ne." ")Then
        j = j + 1
        output(j:j) = outputTemp(i:i)
      End If
    End Do
  End Function RemoveSpaces
!---------------------------------------------------------------------------------------------------------------------------------------
    Function TrimSpaces (input) RESULT (output)
      Character(*), INTENT(IN) :: input
      Character(LEN(trim(adjustl(input)))) :: output
      output = trim(adjustl(input))
    End Function TrimSpaces

    Function BlankString (input) RESULT (output)
      Character(*), INTENT(IN) :: input
      Character(Len(input)) :: output
      Integer(kind=StandardInteger) :: i
      Do i=1,Len(input)
        output(i:i) = " "
      End Do
    End Function BlankString

    Function BlankStringArray (input) RESULT (output)
      Character(*), Dimension(:), INTENT(IN) :: input
      Character(Len(input)) :: line
      Character(Len(input)), Dimension(1:size(input,1)) :: output
      Integer(kind=StandardInteger) :: i
      Do i=1,Len(input)
        line(i:i) = " "
      End Do
      Do i=1,size(input,1)
        output(i) = line
      End Do
    End Function BlankStringArray

  Function BlankString2DArray (input) RESULT (output)
    Character(*), Dimension(:,:), INTENT(IN) :: input
    Character(Len(input)) :: line
    Character(Len(input)), Dimension(1:size(input,1),1:size(input,2)) :: output
    Integer(kind=StandardInteger) :: i, j
    Do i=1,Len(input)
      line(i:i) = " "
    End Do
    Do i=1,size(input,1)
      Do j=1,size(input,2)
        output(i,j) = line
      End Do
    End Do
  End Function BlankString2DArray
    
  Function SpacesRight (input) RESULT (output)
! Adds spaces to right of string
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: tempStr, output
    Integer(kind=StandardInteger) :: i
    tempStr = trim(adjustl(input))
    Do i=1,Len(tempStr)
      output(i:i) = " "
    End Do
    Do i=1,Len(trim(tempStr))
      output(i:i) = tempStr(i:i)
    End Do
  End Function SpacesRight

  Function RemoveComments (input) RESULT (output)
! Removes comments from
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i
    output = BlankString(output)
! Comment character !
    Do i=1,Len(input)
      If(input(i:i).eq."!")Then
        Exit
      Else
        output(i:i) = input(i:i)
      End If
    End Do
  End Function RemoveComments

  Function RemoveQuotes (input) RESULT (output)
! Removes comments from
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i, j
    output = BlankString(output)
! Comment character !
    j = 0
    Do i=1,Len(input)
      If(ichar(input(i:i)).eq.34.or.ichar(input(i:i)).eq.39)Then
! Do nothing
      Else
        j = j + 1
        output(j:j) = input(i:i)
      End If
    End Do
  End Function RemoveQuotes
  
  Function IntToStr (input) RESULT (output)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables  
  ! In:      Declare variables  
    Integer(kind=StandardInteger) :: input
  ! Out:     Declare variables    
    Character(Len=16) :: output
    Write(output,"(I16)") input
    output = trim(adjustl(output))
  End Function IntToStr
  
  Function DpToStr (input) RESULT (output)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables  
  ! In:      Declare variables  
    Real(kind=DoubleReal) :: input
  ! Out:     Declare variables    
    Character(Len=16) :: output
    Write(output,"(E12.5)") input
    output = trim(adjustl(output))
  End Function DpToStr
  
  Function StrToInt (input) RESULT (output)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables  
  ! In:      Declare variables  
    Character(*) :: input
  ! Out:     Declare variables    
    Integer(kind=StandardInteger) :: output
    output = 0
    Read(input,*) output
  End Function StrToInt
  
  Function StrToDp (input) RESULT (output)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables  
  ! In:      Declare variables  
    Character(*) :: input
  ! Out:     Declare variables    
    Real(kind=DoubleReal) :: output
    output = 0.0D0
    Read(input,*) output
  End Function StrToDp
  
  Function RandName(randSwitchIn, prefixIn) Result (randNameOut)
! Make a random 8 character "name"
    Implicit None  ! Force declaration of all variables  
! In:      Declare variables
    Logical, Optional :: randSwitchIn
    Character(*), Optional :: prefixIn
! Out:     Declare variables    
    Character(len=8) :: randNameOut
! Private: Declare variables
    Logical :: randSwitch
    Character(len=8) :: prefix
    Integer(kind=StandardInteger) :: i
    Character(len=1) :: randChar
    Logical :: writePrefix = .true.
! Optional
    randSwitch = .true.
    prefix = "        "
    If(Present(randSwitchIn))Then
      randSwitch = randSwitchIn
    End If
    If(Present(prefixIn))Then
      prefix = prefixIn
    End If
! Init output   
    randNameOut = "        "
! Loop through letters    
    Do i=1,8
      If(prefix(i:i).eq." ")Then
        writePrefix = .false.
      End If
      If(writePrefix)Then
        randNameOut(i:i) = prefix(i:i)
      Else
        Call randCharacter(randChar, randSwitchIn, 2)
        randNameOut(i:i) = randChar
      End If
    End Do
  End Function RandName  
    
  Function TempFileName(randSwitchIn) Result (fileNameOut)
! Make random name for temp file
    Implicit None  ! Force declaration of all variables  
! In:      Declare variables
    Logical, Optional :: randSwitchIn
! Out:     Declare variables    
    Character(len=8) :: fileNameOut  
! Private: Declare variables
    Logical :: randSwitch
! Optional
    randSwitch = .true.
    If(Present(randSwitchIn))Then
      randSwitch = randSwitchIn
    End If    
! Init
    fileNameOut = "        " 
! Make name    
    fileNameOut = RandName(randSwitch, "tmp")  
  End Function TempFileName
      
  Function CleanString(stringIn) Result (stringOut)
! Clean a string - a-z A-Z 0-9 space
    Implicit None  ! Force declaration of all variables 
! In:      Declare variables
    Character(*) :: stringIn  
! Out:     Declare variables  
    Character(Len(stringIn)) :: stringOut     
! Private: Declare variables
    Integer(kind=StandardInteger) :: i, n, charVal
! Blank output string
    stringOut = BlankString(stringOut)
    n = 0
    Do i=1,len(stringIn)
      charVal = ichar(stringIn(i:i))
      If(charVal.eq.32.or.&
      (charVal.ge.48.and.charVal.le.57).or.&
      (charVal.ge.65.and.charVal.le.90).or.&
      (charVal.ge.97.and.charVal.le.122))Then
        n = n + 1
        stringOut(n:n) = stringIn(i:i)
      End If
    End Do
  End Function CleanString
  
  
  
  
  
  
! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------  
  
  Subroutine randCharacter(letter, randSwitchIn, setIn)
! In/Out:      Declare variables
    Character(len=1) :: letter
    Logical, Optional :: randSwitchIn
    Integer(kind=StandardInteger), Optional :: setIn
! Private: Declare variables
    Logical :: randSwitch
    Integer(kind=StandardInteger) :: set
    Integer(kind=StandardInteger) :: characterNum
    Real(kind=DoubleReal) :: randNumber
    Character(len=52), Parameter :: alpha = &
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(len=62), Parameter :: alphaNum = &
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890'
! Optional arguments
    randSwitch = .false.   ! false use repeatable rng, true use random based off cpu clock
    set = 1                ! 1 = alpha, 2 = alpha numeric
    If(Present(randSwitchIn))Then
      randSwitch = randSwitchIn
    End If
    If(Present(setIn))Then
      set = setIn
    End If
    letter = " "
! Make character
    If(randSwitch)Then
      randNumber = RandomLCG_R_strings()
    Else
      randNumber = RandomLCG_strings()
    End If
    If(set.eq.1)Then
      characterNum = Ceiling(62.0D0*randNumber+1.0D0)
      If(characterNum.lt.1)Then
        characterNum = 1
      End If
      If(characterNum.gt.52)Then
        characterNum = 52
      End If
      letter = alpha(characterNum:characterNum)
    End If
    If(set.eq.2)Then
      characterNum = Ceiling(52.0D0*randNumber+1.0D0)
      If(characterNum.lt.1)Then
        characterNum = 1
      End If
      If(characterNum.gt.62)Then
        characterNum = 62
      End If
      letter = alphaNum(characterNum:characterNum)
    End If
  End Subroutine randCharacter
  
  
!---------------------------------------------------------------------------------------------------------------------------------------
! String random number functions (as these are loaded AFTER strings MOD with the rng MOD
!---------------------------------------------------------------------------------------------------------------------------------------

  Function RandomLCG_strings(seedIn) RESULT (output)
! Random number - linear congruential generator
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=LongInteger) :: m, a, c, clockTime
    Integer(kind=LongInteger) :: seed
    Integer(kind=StandardInteger), Optional :: seedIn
    Real(kind=DoubleReal) :: output
! init
    seed = 0
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! Read input, reset counter
    If(Present(seedIn))Then
      seed = seedIn
      If(seed.lt.0)Then ! If seed -1 (for example) get random seed
        Call SYSTEM_CLOCK(clockTime) ! "nano seconds" - well, an estimate
        seed = mod(clockTime,m) ! fit in m
      End If
      randomLCG_n_strings = 0
    End If
! If first iteration
    If(randomLCG_n_strings.eq.0)Then
      If(seed.eq.0)Then
        seed = 12791244+45778951 ! Use default seed
      End If
      randomLCG_n_strings = 0
      randomLCG_xn_strings = seed
    End If
! Increment counter
    randomLCG_n_strings = randomLCG_n_strings + 1
! calculate
    randomLCG_xn_strings = mod((a*randomLCG_xn_strings+c),m)
    output = (1.0D0*randomLCG_xn_strings)/(1.0D0*m)
  End Function RandomLCG_strings 
  
  Function RandomLCG_R_strings() RESULT (output)
! Random number - linear congruential generator
! This function starts with a random seed
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=LongInteger) :: m, a, c, clockTime
    Integer(kind=LongInteger) :: seed
    Real(kind=DoubleReal) :: output
! init
    seed = 0
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! If first iteration
    If(randomLCG_R_n_strings.eq.0)Then
! Make "random" seed
      Call SYSTEM_CLOCK(clockTime) ! "nano seconds" - well, an estimate
      seed = mod(clockTime,m)      
      Do i=1,10
        seed = mod((a*seed+c),m)
      End Do
      randomLCG_R_n_strings = 0
      randomLCG_R_xn_strings = seed
    End If
! Increment counter
    randomLCG_R_n_strings = randomLCG_R_n_strings + 1
! calculate
    randomLCG_R_xn_strings = mod((a*randomLCG_R_xn_strings+c),m)
    output = (1.0D0*randomLCG_R_xn_strings)/(1.0D0*m)
  End Function RandomLCG_R_strings

End Module strings