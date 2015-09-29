Module strings
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
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
  
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains 
!---------------------------------------------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------
End Module strings