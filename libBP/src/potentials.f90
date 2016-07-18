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
  Integer(kind=StandardInteger), Parameter :: p_potPoints = 1001
! Make private
  Private
! Public Variables and Parameters
  Public :: p_potentials, p_potPoints
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
    Logical, Dimension(1:p_potentials) :: pairFunction = .false.
    Character(Len=1), Dimension(1:p_potentials) :: fType         ! e.g. A(nalytic) N(umeric)
    Character(Len=16), Dimension(1:p_potentials) :: fPotential    ! MORSE, LJ, EAMPAIR, EAMDENS
!   Analytic vars
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:12) :: keyParameter = 0.0D0
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: parameterCount = 0
!   Numeric vars
    Character(Len=128), Dimension(1:p_potentials) :: tabFilePath
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: fPointCount = 0
    Real(kind=DoubleReal), Dimension(1:p_potentials) :: xMin = 0.0D0   ! Angstrom
    Real(kind=DoubleReal), Dimension(1:p_potentials) :: xMax = 10.0D0   ! Angstrom
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:p_potPoints,1:4) :: dataPoints = 0.0D0
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:p_potPoints,1:4) :: dataPoints_C = 0.0D0  ! Calculated data points, from numeric and analytic, with ZBL, if switched ON
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:p_potPoints,1:2) :: splineNodes = 0.0D0
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:p_potPoints,1:4) :: splinePoints = 0.0D0
!   Atom IDs
    Character(Len=16), Dimension(1:64) :: atomIDs
    Integer(kind=StandardInteger) :: atomID_Count
! ZBL Hard Core Settings
    Logical, Dimension(1:p_potentials) :: zblOn = .false.
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: zbl_ZA
    Integer(kind=StandardInteger), Dimension(1:p_potentials) :: zbl_ZB
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:4) :: zbl_rA      !  A x,y,y',y''
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:4) :: zbl_rB      !  B x,y,y',y''
! ZBL to function spline
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:4) :: zblExpSpline
! Potential Fitting
    Real(kind=DoubleReal), Dimension(1:p_potentials) :: fitMin = 0.0D0
    Real(kind=DoubleReal), Dimension(1:p_potentials) :: fitMax = 10.0D0
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:3) :: morseFit = 0.0D0   ! De, a, rc
    Real(kind=DoubleReal), Dimension(1:p_potentials,1:4) :: zbl_rB_Morse      !  B x,y,y',y''  exp(p(x)) spline zbl to morse potential point B

! analyticType:
! MORSE
  End Type potentialType

  Type :: potentialSearchType
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger):: atomID_A
    Integer(kind=StandardInteger):: atomID_B
    Integer(kind=StandardInteger):: fN   ! only used by plain search
    Character(Len=16) :: fPotential      ! only used by search by name
    Logical :: zblOverride = .false.     ! Use to override any ZBL setting to get point if ZBL was ignored
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
  Use sort
  Use constants
  Use printModTypes
  Use printMod
  Use matrix
  Use basicMaths
  Use rng
  Use linearAlgebra
  Use calcFunctions
  Use scienceFunctions
  Use interpolation
  Use coordFunctions
  Use geomTypes
  Use potentialsTypes
  Use isotopesTypes
  Use isotopes
  Use splinesFitting
! Force declaration of all variables
  Implicit None
  Private
! ---- Variables
!  Public :: nl
! ---- Subroutines
  Public :: initPotential
  Public :: loadPotential
  Public :: fitStandardPotentials
  Public :: printPotentialSummary
  Public :: SearchPotential
  Public :: SearchPotential_Name
  Public :: outputPotential



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
    potential%xMin = 0.0D0
    potential%xMax = 10.0D0
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
    Call initPotential(potential)
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
! ZBL On
      If(tempField(1:7).eq."#ZBL_ON")Then
        tempField = StrToUpper(fieldArray(2))
        potential%zblOn(fN) = .false.
        If(TestBoolStr(tempField))Then
          potential%zblOn(fN) = .true.
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
! Set pair function
        potential%pairFunction(fN) = .true.
        If(StrMatch(potential%fPotential(fN),"DENS"))Then
          potential%pairFunction(fN) = .false.
        End If
        If(StrMatch(potential%fPotential(fN),"DDEN"))Then
          potential%pairFunction(fN) = .false.
        End If
        If(StrMatch(potential%fPotential(fN),"SDEN"))Then
          potential%pairFunction(fN) = .false.
        End If
        If(StrMatch(potential%fPotential(fN),"EMBE"))Then
          potential%pairFunction(fN) = .false.
        End If
        If(StrMatch(potential%fPotential(fN),"DEMB"))Then
          potential%pairFunction(fN) = .false.
        End If
        If(StrMatch(potential%fPotential(fN),"SEMB"))Then
          potential%pairFunction(fN) = .false.
        End If
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
! ZBL Parameters
      If(tempField(1:4).eq."#ZBL")Then
        potential%zbl_rA(fN,1) = StrToDp(fieldArray(2))
        potential%zbl_rB(fN,1) = StrToDp(fieldArray(3))
      End If
    End Do
! Store function count
    potential%fCount = fN
! Load tabulated data
    Do fN=1,potential%fCount
      If(potential%fType(fN).eq."N")Then
        If(FileExists(potential%tabFilePath(fN)))Then
          Call readFile(potential%tabFilePath(fN), fileArray, n)
! Read each line of the file
          Do j=1,n
            Call readFieldsCharacter(fileArray(j),fieldArray,fieldCount)
! Read in data points
            dataPointsIn(j,1) = StrToDP(fieldArray(1))
            dataPointsIn(j,2) = StrToDP(fieldArray(2))
          End Do
! Sort data points
          !Call sortArray(dataPointsIn,"ASC",1,n)
! min/max
          pointCount = n
          xMin = dataPointsIn(1,1)
          xMax = dataPointsIn(pointCount,1)
          xInc = (xMax - xMin)/1000.0D0
! Store xMin, xMax
          potential%xMin(fN) = xMin
          potential%xMax(fN) = xMax
! Interpolate to 1001 points for f(x) and f'(x)
          x = xMin
          Do j=1,1001
            yArray = PointInterp(dataPointsIn,x,4,2,1,pointCount)
            potential%dataPoints(fN,j,1) = x
            potential%dataPoints(fN,j,2) = yArray(1)
            potential%dataPoints(fN,j,3) = yArray(2)
            potential%dataPoints(fN,j,4) = yArray(3)
! Increment x
            x = x + xInc
          End Do
        End If
      End If
    End Do
! Update potential
    Call updatePotential(potential)
  End Subroutine loadPotential


  Subroutine updatePotential(potential)
! Update Potential
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialType) :: potential
! Vars:  Private
    Integer(kind=StandardInteger) :: fN
    Type(elementsObj) :: elementsList
    Type(elementObj) :: elementDetails
    Real(kind=DoubleReal), Dimension(1:3) :: zblParameters
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), Dimension(1:4) :: expCoefficients
    Type(potentialSearchType) :: searchObj
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: xInc, x
! Load an elements list object
    Call loadElements(elementsList)    !isotopes.f90
! Loop through potentials
    Do fN = 1,potential%fCount
! Search for matching elements from atom label inputs
      elementDetails = SearchElements(potential%atomLabel_A(fN), elementsList)
      potential%zbl_ZA(fN) = elementDetails%atomicNumber
      elementDetails = SearchElements(potential%atomLabel_B(fN), elementsList)
      potential%zbl_ZB(fN) = elementDetails%atomicNumber
! Default fit Min/Max
      potential%fitMin(fN) = potential%xMin(fN)
      potential%fitMax(fN) = potential%xMax(fN)
! Calculate ZBL points
      If((potential%zbl_ZA(fN).gt.0).and.(potential%zbl_ZB(fN).gt.0))Then
        If((potential%zbl_rA(fN,1).gt.0.0).and.(potential%zbl_rB(fN,1).gt.potential%zbl_rA(fN,1)))Then
! Only pair functions
          If(potential%pairFunction(fN))Then
! Store point A
            zblParameters(1) = potential%zbl_rA(fN,1)
            zblParameters(2) = potential%zbl_ZA(fN)
            zblParameters(3) = potential%zbl_ZB(fN)
            yArray = F_ZblFull(zblParameters)
            potential%zbl_rA(fN,2) = yArray(1)
            potential%zbl_rA(fN,3) = yArray(2)
            potential%zbl_rA(fN,4) = yArray(3)
! Store point B
            searchObj%fN = fN
            searchObj%x = potential%zbl_rB(fN,1)
            yArray = SearchPotential(searchObj, potential)
            potential%zbl_rB(fN,2) = yArray(1)
            potential%zbl_rB(fN,3) = yArray(2)
            potential%zbl_rB(fN,4) = yArray(3)
            print *,potential%zbl_rB(fN,1),potential%zbl_rB(fN,2),potential%zbl_rB(fN,3),potential%zbl_rB(fN,4)
! Exp(P(x)) spline from ZBL to potential function
            expCoefficients = SplineExpThird(potential%zbl_rA(fN,1),potential%zbl_rA(fN,2),&
            potential%zbl_rA(fN,3),potential%zbl_rB(fN,1),potential%zbl_rB(fN,2),potential%zbl_rB(fN,3))
            potential%zblExpSpline(fN,1) = expCoefficients(1)
            potential%zblExpSpline(fN,2) = expCoefficients(2)
            potential%zblExpSpline(fN,3) = expCoefficients(3)
            potential%zblExpSpline(fN,4) = expCoefficients(4)
! Change fit min
            potential%fitMin(fN) = potential%zbl_rB(fN,1)  ! Point b, x value
          End If
        End If
      End If
! Generate calculated data points
! from both ANALYTIC and NUMERIC potentials
      xInc = (potential%xMax(fN)-potential%xMin(fN))/1000
      x = potential%xMin(fN)
      searchObj%fN = fN
      Do i=1,1001
        searchObj%x = x
        !print *,i,x,potential%xMax(fN)
        yArray = SearchPotential(searchObj, potential)
! Store
        potential%dataPoints_C(fN,i,1) = x
        potential%dataPoints_C(fN,i,2) = yArray(1)
        potential%dataPoints_C(fN,i,3) = yArray(2)
        potential%dataPoints_C(fN,i,4) = yArray(3)
! Increment
        x = x + xInc
      End Do
    End Do
  End Subroutine updatePotential

! ----------------------------------------------------
! Potential fitting
! ----------------------------------------------------


  Subroutine fitStandardPotentials(potential)
! Fit potential to standard analytic potentials
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(potentialType) :: potential
! Vars:  Private
    Integer(kind=StandardInteger) :: fN, i, j, k, n
    Real(kind=DoubleReal), Dimension(1:p_potPoints,1:2) :: tempDataPoints
    Real(kind=DoubleReal), Dimension(1:20,1:2) :: tempDataPointsReduced
    Real(kind=DoubleReal), Dimension(1:3) :: morseParameters
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Loop through functions
    Do fN = 1,potential%fCount
      If(potential%pairFunction(fN))Then
        j = 0
        Do i=1,p_potPoints
          If(potential%dataPoints_C(fN,i,1).ge.potential%fitMin(fN))Then
            If(potential%dataPoints_C(fN,i,1).le.potential%fitMax(fN))Then
              j = j + 1
              tempDataPoints(j,1) = potential%dataPoints_C(fN,i,1)
              tempDataPoints(j,2) = potential%dataPoints_C(fN,i,2)
            End If
          End If
        End Do
        n = j
! Reduce points
        k = 1
        If(n.gt.10)Then
          k = ceiling(n/10.0D0)
        End If
        i = 1
        j = 0
        Do While(i.lt.n)
          j = j + 1
          tempDataPointsReduced(j,1) = tempDataPoints(i,1)
          tempDataPointsReduced(j,2) = tempDataPoints(i,2)
          i = i + k
        End Do
!
! Fit Morse Potential
!----------------------------
        Call fitPot_Morse(tempDataPointsReduced, j, morseParameters)
        potential%morseFit(fN,1) = morseParameters(1)
        potential%morseFit(fN,2) = morseParameters(2)
        potential%morseFit(fN,3) = morseParameters(3)
        potential%zbl_rB_Morse(fN,1) = potential%zbl_rB_Morse(fN,1)
        yArray = F_MorseFull(morseParameters,potential%zbl_rB_Morse(fN,1))
        potential%zbl_rB_Morse(fN,2) = yArray(1)
        potential%zbl_rB_Morse(fN,3) = yArray(2)
        potential%zbl_rB_Morse(fN,4) = yArray(3)
!
! Fit LJ Potential
!----------------------------

      End If
    End Do
  End Subroutine fitStandardPotentials
! -------------------------
  Subroutine fitPot_Morse(tempDataPoints, pointCount, parameters)
! Fit potential to standard analytic potentials
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Real(kind=DoubleReal), Dimension(:,:) :: tempDataPoints
    Integer(kind=StandardInteger) :: pointCount
    Real(kind=DoubleReal), Dimension(1:3) :: parameters
! Vars:  Private
    Real(kind=DoubleReal), Dimension(1:pointCount,1:2) :: dataPoints
    Integer(kind=StandardInteger) :: i
! Transfer data
    Do i=1,pointCount
      dataPoints(i,1) = tempDataPoints(i,1)
      dataPoints(i,2) = tempDataPoints(i,2)
    End Do
! Fit data points
    parameters = MorseFit(dataPoints)
  End Subroutine fitPot_Morse


! ----------------------------------------------------
! Printing
! ----------------------------------------------------


  Subroutine printPotentialSummary(potential)
! Add a morse potential
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
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



  Function SearchPotential (searchObj, potential, zblOverrideIn) Result (yArray)
! Search for y(x), y'(x), y''(x) from a potential
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialSearchType) :: searchObj
    Type(potentialType) :: potential
    Logical, Optional :: zblOverrideIn
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Logical :: zblOverride
! Optional arguments
    zblOverride = .false.
    If(Present(zblOverrideIn))Then
      zblOverride = zblOverrideIn
    End If
! Check in range
    If((searchObj%fN.ge.1).and.(searchObj%fN.le.p_potentials))Then
      If(zblOverride)Then
        Call SearchPotential_A(searchObj, potential, yArray)
      Else ! No ZBL override
        If(potential%zblOn(searchObj%fN))Then
! ZBL Core
          If(searchObj%x.le.potential%zbl_rA(searchObj%fN,1))Then
            yArray = GetPotentialZBL(searchObj, potential)
          End If
! EXP(P(x)) Spline
          If((searchObj%x.gt.potential%zbl_rA(searchObj%fN,1)).and.&
            (searchObj%x.lt.potential%zbl_rB(searchObj%fN,1)))Then
            yArray = GetPotentialExpP(searchObj, potential)
          End If
! Underlying potential function
          If(searchObj%x.ge.potential%zbl_rB(searchObj%fN,1))Then
            Call SearchPotential_A(searchObj, potential, yArray)
          End If
        Else
          Call SearchPotential_A(searchObj, potential, yArray)
        End If
      End If
    End If
  End Function SearchPotential
!---------------------
  Subroutine SearchPotential_A(searchObj, potential, yArray)
!search subroutine
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(potentialSearchType) :: searchObj
    Type(potentialType) :: potential
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:10) :: parameters
! Search and interpolate if numeric
!-----------------------
! NUMERIC POTENTIALS
!-----------------------
    If(potential%fType(searchObj%fN).eq."N")Then
! Numeric - search through data points and interpolate
      yArray = PointInterp3DArr(potential%dataPoints,searchObj%x,searchObj%fN,4,2)
    End If
!-----------------------
! ANALYTIC POTENTIALS
!-----------------------
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
  End Subroutine SearchPotential_A

  Function GetPotentialZBL (searchObj, potential) Result (yArray)
! Function to get ZBL point with minimal input
! Uses F_ZblFull function from scienceFunctions.f90
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialSearchType) :: searchObj
    Type(potentialType) :: potential
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Real(kind=DoubleReal), Dimension(1:3) :: zblParameters
! Set parameters
    zblParameters(1) = searchObj%x
    zblParameters(2) = potential%zbl_ZA(searchObj%fN)
    zblParameters(3) = potential%zbl_ZB(searchObj%fN)
    yArray = F_ZblFull (zblParameters)
  End Function GetPotentialZBL

  Function GetPotentialExpP (searchObj, potential) Result (yArray)
! Function to get ZBL point with minimal input
! Uses F_ZblFull function from scienceFunctions.f90
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialSearchType) :: searchObj
    Type(potentialType) :: potential
! Vars:  Out
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Vars:  Private
    Real(kind=DoubleReal), Dimension(1:4) :: expPCoeffs
! Set parameters
    expPCoeffs(1) = potential%zblExpSpline(searchObj%fN,1)
    expPCoeffs(2) = potential%zblExpSpline(searchObj%fN,2)
    expPCoeffs(3) = potential%zblExpSpline(searchObj%fN,3)
    expPCoeffs(4) = potential%zblExpSpline(searchObj%fN,4)
    yArray(1) = CalcPolynomialExp(expPCoeffs,searchObj%x,0)
    yArray(2) = CalcPolynomialExp(expPCoeffs,searchObj%x,1)
    yArray(3) = CalcPolynomialExp(expPCoeffs,searchObj%x,2)
  End Function GetPotentialExpP


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




! -----------------------------------------------
!        Output Potential Functions
!
! -----------------------------------------------


  Subroutine outputPotential(potential,fileDirectory)
! Output Potential to file
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialType) :: potential
    Character(*) :: fileDirectory
! Vars:  Private
    Integer(kind=StandardInteger) :: fN, i
    Character(Len=32) :: fileName
    Character(Len=128) :: filePath

    print *,"Output to ",trim(fileDirectory)

! Loop through potential functions
    Do fN = 1,potential%fCount
      fileName = "pot_out."//Trim(IntToStr(fN))//".pot"
      filePath = trim(fileDirectory)//"/"//fileName
      Open(UNIT=118,FILE=Trim(filePath))
! Loop through potential
      Do i=1,Size(potential%dataPoints,2)
        write(118,"(E17.10,A1,E17.10,A1,E17.10,A1,E17.10)") &
        potential%dataPoints(fN,i,1),",",&
        potential%dataPoints(fN,i,2),",",&
        potential%dataPoints(fN,i,3),",",&
        potential%dataPoints(fN,i,4)
      End Do
      Close(118)
    End Do
! Loop through potential functions
    Do fN = 1,potential%fCount
      fileName = "pot_out_c."//Trim(IntToStr(fN))//".pot"
      filePath = trim(fileDirectory)//"/"//fileName
      Open(UNIT=118,FILE=Trim(filePath))
! Loop through potential
      Do i=1,Size(potential%dataPoints_C,2)
        write(118,"(E17.10,A1,E17.10,A1,E17.10,A1,E17.10)") &
        potential%dataPoints_C(fN,i,1),",",&
        potential%dataPoints_C(fN,i,2),",",&
        potential%dataPoints_C(fN,i,3),",",&
        potential%dataPoints_C(fN,i,4)
      End Do
      Close(118)
    End Do

  End Subroutine outputPotential


End Module potentials

























!-----------------------------------
