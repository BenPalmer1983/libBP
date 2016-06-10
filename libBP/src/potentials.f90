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

  Type :: potentialType   ! approx 5MB
    Character(Len=64) :: label
    Integer(kind=StandardInteger) :: fCount   ! Number of potential functions
    Character(Len=16), Dimension(1:32) :: atomLabel_A
    Character(Len=16), Dimension(1:32) :: atomLabel_B
    Integer(kind=StandardInteger), Dimension(1:32) :: atomID_A
    Integer(kind=StandardInteger), Dimension(1:32) :: atomID_B
    Integer(kind=StandardInteger) :: IDcount
    Integer(kind=StandardInteger), Dimension(1:32) :: uKey
    Character(Len=1), Dimension(1:32) :: fType         ! e.g. A(nalytic) N(umeric)
    Character(Len=16), Dimension(1:32) :: fPotential    ! MORSE, LJ, EAMPAIR, EAMDENS
!   Analytic vars
    Real(kind=DoubleReal), Dimension(1:32,1:12) :: keyParameter = 0.0D0
    Integer(kind=StandardInteger), Dimension(1:32) :: parameterCount = 0
!   Numeric vars
    Character(Len=128), Dimension(1:32) :: tabFilePath
    Integer(kind=StandardInteger), Dimension(1:32) :: fPointCount = 0
    Real(kind=DoubleReal), Dimension(1:32) :: rMin = 0.0D0   ! Angstrom
    Real(kind=DoubleReal), Dimension(1:32) :: rMax = 6.5D0   ! Angstrom
    Real(kind=DoubleReal), Dimension(1:32,1:1001,1:4) :: dataPoints = 0.0D0

! Need a ZBL and spline setting

! analyticType:
! MORSE
  End Type potentialType



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
  Use matrix
  Use basicMaths
  Use rng
  Use linearAlgebra
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
    Character(Len=128), Dimension(1:1000) :: fileArray   ! 1MB
    Character(Len=128) :: tempLine
    Character(Len=128), Dimension(1:10) :: fieldArray
    Character(Len=128) :: tempField
    Integer(kind=StandardInteger) :: fN ! function number
    Integer(kind=StandardInteger) :: parameterCount
    Real(kind=DoubleReal) :: tempDP
! Init vars
    fN = 0
! Read in file
    Call readFile(filePath, fileArray, n)
    print *,trim(filePath)
    print *,n
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
          print *,"lj..."
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




    End Do
! Store function count
    potential%fCount = fN
  End Subroutine loadPotential



  Subroutine printPotentialSummary(potential)
! Add a morse potential
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(potentialType) :: potential
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j
!--------
    Print *,"--------------------------------------------------------------"
    Print *,"Functions: ",potential%fCount
    Print *,"--------------------------------------------------------------"
! Loop through functions
    Do i = 1,potential%fCount
      Print "(A16,I2)",       "Function:       ",i
      Print "(A16,A1)",       "Type:           ",potential%fType(i)
      Print "(A16,A16)",      "Potential:      ",potential%fPotential(i)
      If(potential%fType(i).eq."A")Then
        Print "(A16,I2)",     "Parameters:     ",potential%parameterCount(i)
      End If
      Print "(A16,A16,I4)",   "Atom A:         ",potential%atomLabel_A(i),potential%atomID_A(i)
      Print "(A16,A16,I4)",   "Atom B:         ",potential%atomLabel_B(i),potential%atomID_B(i)
      If(potential%fType(i).eq."A")Then
        Print "(A16)",        "Parameters:     "
        Do j = 1,potential%parameterCount(i)
          Print "(I2,A2,E16.8)",j,"  ",potential%keyParameter(i,j)
        End Do
      Else
        Print "(A16)",        "Data File:     "
        Print "(A64)",potential%tabFilePath(i)
      End If

      Print *,"--------------------------------------------------------------"

    End Do


  End Subroutine printPotentialSummary


! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------













End Module potentials
