! --------------------------------------------------------------!
! largeInt
! largeIntType
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 4th March 2016
! ----------------------------------------

Module largeIntTypes
! Setup Modules
  Use kinds

  Type :: liType
! CONTROL
    Integer(kind=StandardInteger) :: sign = 1
    Integer(kind=StandardInteger) :: iSize = 1
    Integer(kind=VeryLongInteger), Dimension(1:32768) ::   iVal=0_VeryLongInteger    ! 9999999999999999999999999999999999999   37 digits
  End Type




End Module largeIntTypes

Module largeInt
! Setup Modules
  Use kinds
  Use logicalMod
  Use strings
  Use largeIntTypes

! Force declaration of all variables
  Implicit None
! Variables
  Integer(kind=VeryLongInteger), Parameter :: vliMaxVal = 9999999999999999999999999999999999999_VeryLongInteger
  Integer(kind=VeryLongInteger), Parameter :: vliAdjust = 10000000000000000000000000000000000000_VeryLongInteger
! Privacy of variables/functions/subroutines
  Private
! Public Vars
  Public :: vliMaxVal
  Public :: vliAdjust
! Public Functions
  Public :: AddLargeInt
  Public :: SubtractLargeInt
! Public Subroutines
  Public :: storeLargeInt
  Public :: printLargeInt



  Contains
! ---------------------------------------------------------------------------------------------------



! -----------------------------------
! FUNCTIONS
! -----------------------------------

  Recursive Function AddLargeInt(liAIn, liBIn) RESULT (liC)
! Add large integers together
! liA + liB
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Type(liType) :: liAIn
    Type(liType) :: liBIn
! Out:     Declare variables
    Type(liType) :: liC
! Private
    Type(liType) :: liA
    Type(liType) :: liB
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: maxSections
    Integer(kind=VeryLongInteger) :: carry = 0_VeryLongInteger
    Integer(kind=VeryLongInteger) :: tempVLI = 0_VeryLongInteger
    Logical :: addFunction, reverseSign
! Store input
    liA = liAIn
    liB = liBIn
! get max sections
    maxSections = liA%iSize
    If(liB%iSize.gt.liA%iSize)Then
      maxSections = liB%iSize
    End If
! Carry on with adding, or use subtraction
    If(liA%sign.gt.0.and.liB%sign.gt.0)Then   ! |A| + |B| = C
      addFunction = .true.
      reverseSign = .false.
    End If
    If(liA%sign.lt.0.and.liB%sign.gt.0)Then   ! (-)|A| + |B| = C   |B|-|A| = C
      addFunction = .false.
      liA%sign = 1
      liB%sign = 1
      liC = SubtractLargeInt(liB,liA)
    End If
    If(liA%sign.gt.0.and.liB%sign.lt.0)Then   ! |A|+(-)|B| = C   |A| - |B| = C
      addFunction = .false.
      liA%sign = 1
      liB%sign = 1
      liC = SubtractLargeInt(liA,liB)
    End If
    If(liA%sign.lt.0.and.liB%sign.lt.0)Then   ! (-)|A| + (-)|B| = C   |A| + |B| = -C
      addFunction = .true.
      reverseSign = .true.
    End If
    If(addFunction)Then
      carry = 0_VeryLongInteger
      Do i=1,maxSections
        tempVLI = liA%iVal(i) + liB%iVal(i) + carry
        carry = 0_VeryLongInteger
        If(tempVLI.gt.vliMaxVal)Then
          tempVLI = tempVLI - (1_VeryLongInteger  + vliMaxVal)
          carry = 1_VeryLongInteger
        End If
        liC%iVal(i) = tempVLI
      End Do
      liC%iSize = maxSections
      If(carry.gt.0_VeryLongInteger)Then
        liC%iVal(maxSections+1) = carry
        liC%iSize = maxSections+1
      End If
      liC%sign = 1
      If(reverseSign)Then
        liC%sign = -1*liC%sign
      End If
    End If
  End Function AddLargeInt

  Recursive Function SubtractLargeInt(liAIn, liBIn) RESULT (liC)
! Subtract two large integers
! liA - liB
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Type(liType) :: liAIn
    Type(liType) :: liBIn
! Out:     Declare variables
    Type(liType) :: liC
! Private
    Type(liType) :: liA
    Type(liType) :: liB
    Type(liType) :: liTemp
    Integer(kind=StandardInteger) :: i, j
    Integer(kind=StandardInteger) :: maxSections
    Integer(kind=VeryLongInteger) :: carry = 0_VeryLongInteger
    Integer(kind=VeryLongInteger) :: tempVLI = 0_VeryLongInteger
    Logical :: subtractFunction, reverseSign, reverseValues
! Store input
    liA = liAIn
    liB = liBIn
! Max number of vli chunks
    maxSections = liA%iSize
    print *,"ok",liA%iSize,liB%iSize,liA%sign,liB%sign
    If(liB%iSize.gt.liA%iSize)Then
      maxSections = liB%iSize
    End If
! Check type of calculation
    If(liA%sign.gt.0.and.liB%sign.gt.0)Then   ! |A| - |B| = C
      subtractFunction = .true.
      reverseSign = .false.
      print *,"sub 1"
    End If
    If(liA%sign.lt.0.and.liB%sign.gt.0)Then   ! (-)|A| - |B| = C    |A| + |B| = -C
      subtractFunction = .false.
      liA%sign = 1
      liB%sign = 1
      liC = AddLargeInt(liA,liB)
      liC%sign = -1*liC%sign
      print *,"sub 2"
    End If
    If(liA%sign.gt.0.and.liB%sign.lt.0)Then   ! |A| - (-)|B| = C    |A| + |B| = C
      subtractFunction = .false.
      liA%sign = 1
      liB%sign = 1
      liC = AddLargeInt(liA,liB)
      liC%sign = -1*liC%sign
      print *,"sub 3"
    End If
    If(liA%sign.lt.0.and.liB%sign.lt.0)Then   ! (-)|A| - (-)|B| = C   |A| - |B| = -C     ( or |B| - |A| = C)
      subtractFunction = .true.
      reverseSign = .true.
      print *,"sub 4"
    End If
! Subtract
    If(subtractFunction)Then
! Force to subtract smallest value from largest, then reverse sign of result
      reverseValues = ALargerThanB(liB,liA)
      If(reverseValues)Then
        liTemp = liA
        liA = liB
        liB = liTemp
        reverseSign = FlipLogical(reverseSign)
      End If
! Set carry value
      carry = 0_VeryLongInteger
      j = 1
      Do i=1,maxSections
        tempVLI = liA%iVal(i) - liB%iVal(i) + carry
        !print *,tempVLI
        !print *,liA%iVal(i)
        !print *,liB%iVal(i)
        !print *,""
        carry = 0_VeryLongInteger
        If(tempVLI.lt.0_VeryLongInteger)Then
          tempVLI = tempVLI + vliAdjust
          carry = -1_VeryLongInteger
        End If
        liC%iVal(i) = tempVLI
        If(liC%iVal(i).ne.0_VeryLongInteger)Then
          j = i
        End If
      End Do
      liC%iSize = j
      liC%sign = 1
      If(reverseSign)Then
        liC%sign = -1*liC%sign
      End If
    End If
  End Function SubtractLargeInt

  Function ALargerThanB(liA, liB) RESULT (result)
! Returns TRUE if A larger than B
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Type(liType) :: liA
    Type(liType) :: liB
! Out:     Declare variables
    Logical :: result
! Private
! Default A>B
    result = .true.
! Max number of vli chunks
    If(liB%iSize.gt.liA%iSize)Then
      result = .false.
    End If
    If(liA%iSize.eq.liB%iSize)Then
      If(liA%iVal(liA%iSize).lt.liB%iVal(liB%iSize))Then
        result = .false.
      End If
    End If
  End Function ALargerThanB


! -----------------------------------
! SUBROUTINE
! -----------------------------------

  Subroutine storeLargeInt(objIn, valIn)
! Reset data objects for a plot
    Implicit None  ! Force declaration of all variables
! Input
    Type(liType) :: objIn
    Character(*) :: valIn
! Private
    Integer(kind=StandardInteger) :: i, j, k, l, valLen
    Character(Len=37) :: tempInt
! Clear out
    Do i=1,20
      objIn%iVal(i) = 0_VeryLongInteger
    End Do
!
    objIn%sign = 1
    valLen = len(valIn)
    k = 38
    l = 1
    tempInt = BlankString(tempInt)
    Do i=1,valLen
      j = valLen+1-i
      If(valIn(j:j).eq." ")Then
        Exit
      End If
      If(valIn(j:j).eq."-")Then
        objIn%sign = -1
        Exit
      End If
      k = k - 1
      tempInt(k:k) = valIn(j:j)
      If(k.eq.1)Then
        Read(tempInt,"(I37)"), objIn%iVal(l)
        objIn%iSize = l
        l = l + 1
        k = 38
        tempInt = BlankString(tempInt)
      End If
    End Do
    If(k.gt.1)Then
      Read(tempInt,"(I37)"), objIn%iVal(l)
      objIn%iSize = l
    End If
  End Subroutine storeLargeInt

  Subroutine printLargeInt(objIn)
! Reset data objects for a plot
    Implicit None  ! Force declaration of all variables
! Input
    Type(liType) :: objIn
! Private
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: leadingZeros
    Character(Len=1024) :: printStr
    Character(Len=37) :: chunk, leadingZerosStr
! Blank output string
    printStr = BlankString(printStr)
    j = objIn%iSize
    If(objIn%sign.lt.0)Then
      printStr(1:1) = "-"
    End If
    Do i=1,objIn%iSize
      Write(chunk,"(I37)") objIn%iVal(j)
      !print *,len(trim(adjustl(chunk))),chunk
      leadingZeros = 37-len(trim(adjustl(chunk)))
      leadingZerosStr = BlankString(leadingZerosStr)
      If(leadingZeros.gt.0.and.i.gt.1)Then
        Do k=1,leadingZeros
          leadingZerosStr(k:k) = "0"
        End Do
      End If
      !print *,objIn%iVal(j)
      printStr = trim(printStr)//trim(adjustl(leadingZerosStr))//trim(adjustl(chunk))
      j = j - 1
    End Do
! Output to terminal
    print *,trim(printStr)
  End Subroutine printLargeInt





  Subroutine multLargeInt(vli)
! Reset data objects for a plot
    Implicit None  ! Force declaration of all variables
! Input
    Integer(kind=StandardInteger) :: i, x
    Integer(kind=LongInteger) :: vli
! Private
    x = 2147483647_StandardInteger
    print *,x
    vli = 9999999999999999_LongInteger
    print *,vli

    vli = -2_LongInteger
    Do i=1,48
      print *,i,vli
      vli = vli * 2_LongInteger
    End Do



  End Subroutine multLargeInt


End Module largeInt
