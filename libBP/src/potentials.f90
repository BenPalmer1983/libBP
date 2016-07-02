! --------------------------------------------------------------!
! Static Caculations module
! staticCalcsTypes, staticCalcs
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calculates the forces between atoms
! Calculates the energies of atoms/total energy
! Calculates stresses
!
! ----------------------------------------
! Updated: 21st May 2016
! ----------------------------------------

Module potentialsTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
  Integer(kind=StandardInteger), Parameter :: p_potentials = 32
! Make private
  Private
! Public Variables and Parameters
  Public :: p_potentials
! Public derived types
  Public :: potentialType, potentialSearchType

  Type :: potentialType   ! approx 5MB
    Character(Len=64) :: label
    Integer(kind=StandardInteger) :: fCount   ! Number of potential functions
    Logical, Dimension(1:p_potentials) :: switchOn = .false.
    Character(Len=16), Dimension(1:p_potentials) :: atomLabel_A
    Character(Len=16), Dimension(1:p_potentials) :: atomLabel_B
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: atomID_A
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: atomID_B
    Integer(kind=StandardInteger) :: IDcount
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: uKey
    Character(Len=1), Dimension(1:p_potentials) :: fType         ! e.g. A(nalytic) N(umeric)
    Character(Len=16), Dimension(1:p_potentials) :: fPotential    ! MORSE, LJ, EAMPAIR, EAMDENS
!   Analytic vars
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:12) :: keyParameter = 0.0D0
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: parameterCount = 0
!   Numeric vars
    Character(Len=128), Dimension(1:p_potentials) :: tabFilePath
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: fPointCount = 0
    Real(kind=DoubleReal), Dimension(1:p_potentials) :: rMin = 0.0D0   ! Angstrom
    Real(kind=DoubleReal), Dimension(1:p_potentials) :: rMax = 6.5D0   ! Angstrom
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:1001,1:4) :: dataPoints = 0.0D0
!   Atom IDs
    Character(Len=16), Dimension(1:64) :: atomIDs
    Integer(kind=StandardInteger) :: atomID_Count
! Need a ZBL and spline setting

! analyticType:
! MORSE
  End Type potentialType

  Type :: potentialSearchType
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger):: atomID_A
    Integer(kind=StandardInteger):: atomID_B
    Integer(kind=StandardInteger):: fN   ! only used by plain search
    Character(Len=16) :: fPotential      ! only used by search by name
  End Type potentialSearchType

  Type :: potentialSearchResultType
    Integer(kind=StandardInteger):: resultCount
    Real(kind=DoubleReal), Dimension(1:100,1:3) :: yValues
  End Type potentialSearchResultType



End Module potentialsTypes


Module potentials
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
  Use strings
  Use general
  Use constants
  Use printModTypes
  Use printMod
  Use matrix
  Use basicMaths
  Use rng
  Use linearAlgebra
  Use interpolation
  Use coordFunctions
  Use geomTypes
  Use potentialsTypes
! Force declaration of all variables
  Implicit None
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: initPotential
  Public :: loadPotential
  Public :: printPotentialSummary
  Public :: SearchPotential
  Public :: SearchPotential_Name



!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------


  Subroutine initPotential(potential)
! Make neighbour list for atoms
! The input coords and alat must be large enough so the same atom does not interact
! with a copy in a periodic cell surrounding the original
! e.g. if the rVerlet cutoff is 5, the alat must be greater than 5
!
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialType) :: potential
! Vars:  Private
    potential%switchOn = .false.
    potential%label = BlankString(potential%label)
    potential%fCount = 0
    potential%atomLabel_A = WipeStringArray(potential%atomLabel_A)
    potential%atomLabel_B = WipeStringArray(potential%atomLabel_B)
    potential%atomID_A = 0
    potential%atomID_B = 0
    potential%IDcount = 0
    potential%uKey = 0
    potential%fType = WipeStringArray(potential%fType)
    potential%fPotential = WipeStringArray(potential%fPotential)
!   Analytic vars
    potential%keyParameter = 0.0D0
    potential%parameterCount = 0
!   Numeric vars
    potential%tabFilePath = WipeStringArray(potential%tabFilePath)
    potential%fPointCount = 0
    potential%rMin = 0.0D0
    potential%rMax = 6.5D0
    potential%dataPoints = 0.0D0
  End Subroutine initPotential





  Subroutine loadPotential(filePath, potential)
! Add a morse potential
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Character(*) :: filePath
    Type(potentialType) :: potential
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, n, fieldCount
    Character(Len=128), Dimension(1:1500) :: fileArray   ! 1MB
    Character(Len=128), Dimension(1:10) :: fieldArray
    Character(Len=128) :: tempField, filePathTemp
    Integer(kind=StandardInteger) :: fN ! function number
    Real(kind=DoubleReal), Dimension(1:5000,1:2) :: dataPointsIn
    Integer(kind=StandardInteger) :: pointCount
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal) :: xMin, xMax, xInc, x
! Init vars
    fN = 0
! Read in file
    Call readFile(filePath, fileArray, n)
! Read through potential file
    i = 0
    Do While (i.le.n)
      i = i + 1
! Read line
      Call readFieldsCharacter(fileArray(i),fieldArray,fieldCount)
      tempField = StrToUpper(fieldArray(1))
! Increment function counter
      If(tempField(1:8).eq."#F_START")Then
        fN = fN + 1
! Set any defaults here
        potential%switchOn(fN) = .true.
      End If
! Function On
      If(tempField(1:5).eq."#F_ON")Then
        tempField = StrToUpper(fieldArray(2))
        potential%switchOn(fN) = .false.
        If(TestBoolStr(tempField))Then
          potential%switchOn(fN) = .true.
        End If
      End If
! Function Type - Analytic/Numeric
      If(tempField(1:7).eq."#F_TYPE")Then
        tempField = StrToUpper(fieldArray(2))
        potential%fType(fN) = "N"
        If(tempField(1:1).eq."A")Then
          potential%fType(fN) = "A"
        End If
      End If
! Function Potential - MORSE/LJ/EAMPAIR etc
      If(tempField(1:12).eq."#F_POTENTIAL")Then
        tempField = StrToUpper(fieldArray(2))
! Store
        potential%fPotential(fN) = trim(tempField)
! Store parameter numbers
        If(potential%fPotential(fN)(1:2).eq."LJ")Then     ! Lennard-Jones
          potential%parameterCount(fN) = 2
        End If
        If(potential%fPotential(fN)(1:5).eq."MORSE")Then  ! Morse
          potential%parameterCount(fN) = 3
        End If
      End If
! Atom Label
      If(tempField(1:6).eq."#LABEL")Then
        tempField = StrToUpper(fieldArray(2))
        potential%atomLabel_A(fN) = Trim(tempField)
        tempField = StrToUpper(fieldArray(3))
        potential%atomLabel_B(fN) = Trim(tempField)
      End If
! Parameters
      If(tempField(1:11).eq."#PARAMETERS")Then
        Do j=1,potential%parameterCount(fN)
          i = i + 1
          Call readFieldsCharacter(fileArray(i),fieldArray,fieldCount)
          potential%keyParameter(fN,j) = StrToDp(fieldArray(1))
        End Do
      End If
! Tabulated Data File
      If(tempField(1:9).eq."#TAB_FILE")Then
        filePathTemp = fieldArray(2)
        Call completePath(filePathTemp)
        potential%tabFilePath(fN) = filePathTemp
      End If
    End Do
! Store function count
    potential%fCount = fN
! Load tabulated data
    Do i=1,fN
      If(potential%fType(i).eq."N")Then
        If(FileExists(potential%tabFilePath(i)))Then
          Call readFile(potential%tabFilePath(i), fileArray, n)
! Read each line of the file
          Do j=1,n
            Call readFieldsCharacter(fileArray(j),fieldArray,fieldCount)
! Read in data points
            dataPointsIn(j,1) = StrToDP(fieldArray(1))
            dataPointsIn(j,2) = StrToDP(fieldArray(2))
          End Do
          pointCount = n
          xMin = dataPointsIn(1,1)
          xMax = dataPointsIn(pointCount,1)
          xInc = (xMax - xMin)/1000.0D0
! Interpolate to 1001 points
          x = xMin
          Do j=1,1001
! x, f(x)
            yArray = PointInterp(dataPointsIn,x,4,2,1,pointCount)
            potential%dataPoints(i,j,1) = x
            potential%dataPoints(i,j,2) = yArray(1)
            potential%dataPoints(i,j,3) = yArray(2)
            potential%dataPoints(i,j,4) = yArray(3)
            !print *,x,yArray(1),yArray(2)
! Increment x
            x = x + xInc
          End Do
        End If
      End If
    End Do
  End Subroutine loadPotential



  Subroutine printPotentialSummary(potential)
  ! Add a morse potential
    Implicit None   ! Force declaration of all variables
  ! Vars:  In
    Type(potentialType) :: potential
  ! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Character(Len=128) :: tempLine
  !--------
    Write(tempLine,*) "Potential Functions Summary (",potential%fCount,")"
    Call addLinePage(tempLine,"T")
  ! Loop through functions
      Do i = 1,potential%fCount
        Write(tempLine,"(I2,A1,A1,A1,A10,A1,A4,I2,A4,I2,A2,I2,A1)") &
        i," ",potential%fType(i)," ",potential%fPotential(i)," ",&
        potential%atomLabel_A(i),potential%atomID_A(i),&
        potential%atomLabel_B(i),potential%atomID_B(i),&
        " (",potential%uKey(i),")"
        Call addLinePage(tempLine)
      End Do
    End Subroutine printPotentialSummary


  Subroutine printPotentialSummary_Full(potential)
! Add a morse potential
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialType) :: potential
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
    Character(Len=64) :: tempLine
!--------
    Write(tempLine,*) "Functions: ",potential%fCount
    Call addLinePage(tempLine,"T")
! Loop through functions
    Do i = 1,potential%fCount
      Write(tempLine,*) "Function:       ",i
      Call addLinePage(tempLine)
      Write(tempLine,*) "Type:           ",potential%fType(i)
      Call addLinePage(tempLine)
      Write(tempLine,*) "Potential:      ",potential%fPotential(i)
      Call addLinePage(tempLine)
      If(potential%fType(i).eq."A")Then
        Write(tempLine,*) "Parameters:     ",potential%parameterCount(i)
        Call addLinePage(tempLine)
      End If
      Write(tempLine,*) "Atom A:         ",potential%atomLabel_A(i),potential%atomID_A(i)
      Call addLinePage(tempLine)
      Write(tempLine,*) "Atom B:         ",potential%atomLabel_B(i),potential%atomID_B(i)
      Call addLinePage(tempLine)
      If(potential%fType(i).eq."A")Then
        Write(tempLine,*) "Parameters:     "
        Call addLinePage(tempLine)
        Do j = 1,potential%parameterCount(i)
          Write(tempLine,*) j,"  ",potential%keyParameter(i,j)
          Call addLinePage(tempLine)
        End Do
      Else
        Call addLinePage("Data File:     ")
        Call addLinePage(potential%tabFilePath(i))
      End If
    End Do
  End Subroutine printPotentialSummary_Full


! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------



  Function SearchPotential (searchObj, potential) Result (yArray)
! Add a morse potential
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialSearchType) :: searchObj
    Type(potentialType) :: potential
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:10) :: parameters
! Check in range
    If((searchObj%fN.ge.1).and.(searchObj%fN.le.p_potentials))Then
! Search and interpolate if numeric
      If(potential%fType(searchObj%fN).eq."N")Then
! Numeric - search through data points and interpolate
        yArray = PointInterp3DArr(potential%dataPoints,searchObj%x,searchObj%fN,4,1)
      End If
      If(potential%fType(searchObj%fN).eq."A")Then
! Load parameters
        Do i=1,potential%parameterCount(searchObj%fN)
          parameters(i) = potential%keyParameter(searchObj%fN,i)
        End Do
! Morse potential
        If(StrMatch(potential%fPotential(searchObj%fN),"MORSE"))Then
          yArray = Pot_Morse(parameters, searchObj%x, 2)
        End If
      End If
    End If
  End Function SearchPotential


  Function SearchPotential_Name (searchObj, potential) Result (yArray)
! Search potential by name
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialSearchType) :: searchObj
    Type(potentialType) :: potential
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Integer(kind=StandardInteger) :: i, fN
    Real(kind=DoubleReal), External :: potFunction
    Real(kind=DoubleReal), Dimension(1:10) :: parameters
! Loop through potential functions to find the right one
    fN = 1
    Do i=1,32
      If(CheckIDMatch(searchObj%atomID_A, searchObj%atomID_B, &
      potential%atomID_A(i), potential%atomID_B(i)))Then
        If(StrMatch(potential%fPotential(i),searchObj%fPotential,.false.))Then
          fN = i
          Exit
        End If
      End If
    End Do
! Search and interpolate if numeric
    If(potential%fType(fN).eq."N")Then
! Numeric - search through data points and interpolate
      yArray = PointInterp3DArr(potential%dataPoints,searchObj%x,fN,4,1)
    End If
    If(potential%fType(fN).eq."A")Then
! Load parameters
      Do i=1,10
        parameters(i) = potential%keyParameter(fN,i)
      End Do
! Morse potential
      If(StrMatch(potential%fPotential(fN),"MORSE"))Then
        yArray = Pot_Morse(parameters, searchObj%x, 2)
      End If
    End If
  End Function SearchPotential_Name


  Function CheckIDMatch (idSearch_A, idSearch_B, idPot_A, idPot_B) Result (result)
! Check if there is a match with the search and potential IDs
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger) :: idSearch_A, idSearch_B, idPot_A, idPot_B
! Vars:  Out
    Logical :: result
! Check
    result = .false.
    If((idSearch_A.eq.idPot_A).and.(idSearch_B.eq.idPot_B))Then
      result = .true.
    End If
    If((idSearch_A.eq.idPot_B).and.(idSearch_B.eq.idPot_A))Then
      result = .true.
    End If
  End Function CheckIDMatch





! -----------------------------------------------
!        Analytic functions
!
! -----------------------------------------------



  Function Pot_Morse (parameters, x, derivativeIn) Result (yArray)
! Calculate potential value - Morse potential
! param 1 = alpa
! param 2 = rc
! param 3 = D
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger), Optional :: derivativeIn
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Integer(kind=StandardInteger) :: derivative
    Real(kind=DoubleReal) :: A,R,D
! Optional Arguments
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! Init vars
    A = parameters(1)
    R = parameters(2)
    D = parameters(3)
! Morse Potential f(x)
    If(derivative.ge.1)Then
      yArray(1) = D*((1.0D0-exp(A*(R-x)))**2-1.0D0)
    End If
    If(derivative.ge.2)Then
      yArray(2) = 2.0D0*D*A*(exp(1.0D0*A*(R-x))-exp(2.0D0*A*(R-x)))
    End If
  End Function Pot_Morse



End Module potentials

























!-----------------------------------
