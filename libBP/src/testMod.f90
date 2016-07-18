Module testMod
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
!
! Test functions for Library
! --------------------------------------------------------------!
  Use kinds
  Use logicalMod
  Use strings
  Use general
  Use envTypes
  Use env
  Use printModTypes
  Use printMod
  Use activityFunctionsTypes
  Use activityFunctions
  Use basicMaths
  Use keysMod
  Use laplaceTransforms
  Use interpolation
  Use regression
  Use splinesFitting
  Use geomTypes
  Use geom
  Use potentialsTypes
  Use potentials
  Use staticCalcsTypes
  Use staticCalcs
  Use isotopesTypes
  Use isotopes
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! --Subroutines--!
  Public :: testActivity
  Public :: testCalcIsotopeAmount
  Public :: testGaverStehfest
  Public :: testDecayChain
  Public :: testNeighbourList
  Public :: testStaticCalc
  Public :: testCombinations
  Public :: testLinearRegression
  Public :: testMorseFit

!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------


! ------------------------------------------------------------------------------------------
!  Test subroutines for activityFunctions.f90
! ------------------------------------------------------------------------------------------

  Subroutine testActivity()
! Tests activity function
    Implicit None ! Force declaration of all variables
! Vars: In/Out
!   None
! Vars: Private
    Print *,"Testing activity function"
! Polonium-218
    Call ActivityCompareGS(0.002236D0, 2000, 100000.0D0, 1.0D6, 1000.0D0)
  End Subroutine testActivity


  Subroutine testCalcIsotopeAmount()
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: n
    Real(kind=DoubleReal) :: w, t
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange

    t = 60.0D0
    w = 0.0D0
    Allocate(decayDataArray(1:5,1:6))
    decayDataArray(1,1) = 1.0D0
    decayDataArray(1,2) = 1.0D6
    decayDataArray(1,3) = 30.0D0
    decayDataArray(1,4) = 1.0D0
    decayDataArray(1,5) = 1.0D0
    decayDataArray(1,6) = 1.0D0
    decayDataArray(2,1) = 2.0D0
    decayDataArray(2,2) = 0.0D0
    decayDataArray(2,3) = 122.0D0
    decayDataArray(2,4) = 0.9998D0
    decayDataArray(2,5) = 1.0D0
    decayDataArray(2,6) = 1.0D0
    decayDataArray(3,1) = 2.0D0
    decayDataArray(3,2) = 0.0D0
    decayDataArray(3,3) = 1196.0D0
    decayDataArray(3,4) = 1.0D0
    decayDataArray(3,5) = 1.0D0
    decayDataArray(3,6) = 1.0D0
    decayDataArray(4,1) = 2.0D0
    decayDataArray(4,2) = 0.0D0
    decayDataArray(4,3) = 1648.0D0
    decayDataArray(4,4) = 0.9998D0
    decayDataArray(4,5) = 1.0D0
    decayDataArray(4,6) = 1.0D0
    decayDataArray(5,1) = 2.0D0
    decayDataArray(5,2) = 0.0D0
    decayDataArray(5,3) = 7006.0D0
    decayDataArray(5,4) = 1.0D0
    decayDataArray(5,5) = 1.0D0
    decayDataArray(5,6) = 1.0D0
    Allocate(isotopeChange(1:5,1:12))
    print *,""
    isotopeChange = CalcIsotopeAmount(w,decayDataArray,t,1)
    Do n = 1,5
      Print *,n,isotopeChange(n,3),isotopeChange(n,4)
    End Do
    print *,""
    isotopeChange = CalcIsotopeAmount(w,decayDataArray,t,2)
    Do n = 1,5
      Print *,n,isotopeChange(n,3),isotopeChange(n,4)
    End Do
  End Subroutine testCalcIsotopeAmount


  Subroutine testGaverStehfest()
! Tests activity function
    Implicit None ! Force declaration of all variables
! Vars: In/Out
!   None
! Vars: Private
    Real(kind=DoubleReal) :: t, ft, exact
    Real(kind=DoubleReal), Dimension(1:10) :: p
    Print *,"Testing G-S Inversion"
! Init
    t = 0.024D0
    p(1) = 0.731D0
    ft = GaverStehfest(ltExp, t, p, 7)
    exact = exp(p(1)*t)
    print *,"m=7 (double) ",t,exact,ft
    ft = GaverStehfest(ltExp, t, p, 8)
    exact = exp(p(1)*t)
    print *,"m=8 (double) ",t,exact,ft
  End Subroutine testGaverStehfest

  Subroutine testDecayChain()
! Tests activity function
    Implicit None ! Force declaration of all variables
! Vars: In/Out
!   None
! Vars: Private
    Type(decayChainObj) :: decayChain
    Type(activityTimeObj), Dimension(1:100) :: activityTime
    Real(kind=DoubleReal) :: endTime, zeroProductionTime
!-------------------------------------------
    Print *,"Testing decay chain"
! Polonium-218 decay chain
    !decayChain%productionRate = 0.0D0
    decayChain%time = 60.0D0
! Parent
    decayChain%label(1) = "Po-218"
    decayChain%productionRate(1) = 5.0D4
    decayChain%branchFactor(1) = 1.0D0  ! Not used
    decayChain%halfLife(1) = 185.88D0
    decayChain%amountStart(1) = 1.0D6
! Daughter 1
    decayChain%label(2) = "Pb-214"
    decayChain%productionRate(2) = 0.0D0
    decayChain%branchFactor(2) = 0.99981D0   ! bf from 1 to 2
    decayChain%halfLife(2) = 1608.0D0
    decayChain%amountStart(2) = 1.0D5
! Daughter 2
    decayChain%label(3) = "Bi-214"
    decayChain%productionRate(3) = 1.0D4
    decayChain%branchFactor(3) = 1.0D0  ! bf from 2 to 3
    decayChain%halfLife(3) = 1194.0D0
    decayChain%amountStart(3) = 500.0D0
! Daughter 3
    decayChain%label(4) = "Po-214"
    decayChain%productionRate(4) = 0.0D0
    decayChain%branchFactor(4) = 0.99979D0  ! bf from 3 to 4
    decayChain%halfLife(4) = 0.0001637D0
    decayChain%amountStart(4) = 0.0D0
! Daughter 4
    decayChain%label(5) = "Pb-210"
    decayChain%productionRate(5) = 0.0D0
    decayChain%branchFactor(5) = 1.0D0  ! bf from 3 to 4
    decayChain%halfLife(5) = 6.9930D8
    decayChain%amountStart(5) = 0.0D0
! Daughter 5
    decayChain%label(6) = "Bi-210"
    decayChain%productionRate(6) = 0.0D0
    decayChain%branchFactor(6) = 1.0D0  ! bf from 4 to 5
    decayChain%halfLife(6) = 1.600665D-6
    decayChain%amountStart(6) = 0.0D0
! Daughter 6
    decayChain%label(7) = "Tl-206"
    decayChain%productionRate(7) = 0.0D0
    decayChain%branchFactor(7) = 1.0D0  ! bf from 5 to 6
    decayChain%halfLife(7) = 2.5212D2
    decayChain%amountStart(7) = 0.0D0
! Daughter 7
    decayChain%label(8) = "Pb-206"
    decayChain%productionRate(8) = 0.0D0
    decayChain%branchFactor(8) = 1.0D0  ! bf from 6 to 7
    decayChain%halfLife(8) = -1.0D0
    decayChain%amountStart(8) = 0.0D0
! Calculate
    endTime = 10000.0D0
    zeroProductionTime = 1000.0D0
!   print *,endTime,zeroProductionTime
    Call CalcActivities(decayChain,activityTime,endTime,zeroProductionTime)
    Call CalcActivitiesPrint(decayChain,activityTime,.false.)

  End Subroutine testDecayChain



  Subroutine testNeighbourList()
! Uses inverse laplace transform to calculate isotope amounts at time t (after time = 0)
! t time in seconds after t=0
! w production rate of parent isotope
! isotope chain data
    Implicit None ! Force declaration of all variables
! Vars Private
    Type(coordsUnitType), Dimension(1:4) :: coordsUnit
    Type(coordsType), Dimension(1:4) :: coords
    Type(nlType), Dimension(1:4) :: nl
    Integer(kind=StandardInteger) :: i, cKey
! Start test
    print *,"Neighbour List Testing"
    print *,"==========================================================="
    Call initUnitCoords(coordsUnit)
    !Call initCoords(coords)
    coordsUnit%aLat = 4.04D0
    coordsUnit%xCopy = 4
    coordsUnit%yCopy = 4
    coordsUnit%zCopy = 4
    !Call standardCoords("FCC", coordsUnit)
    !Call expandUnitCoords(coordsUnit, coords)
    Do cKey=1,size(coords%length)
      Do i=1,coords(cKey)%length
        print *,coords(cKey)%label(i),&
        coords(cKey)%coords(i,1),coords(cKey)%coords(i,2),coords(cKey)%coords(i,3)
      End Do
    End Do
    print *,coords%length
    Call makeNL(nl, coords, 6.5D0)
  End Subroutine testNeighbourList



  Subroutine testStaticCalc()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Type(coordsUnitType), Allocatable, Dimension(:) :: coordsUnit
    Type(coordsType), Allocatable, Dimension(:) :: coords
    Type(nlType), Allocatable, Dimension(:) :: nl
    Type(potentialType) :: potential
    Type(oBulkProperty), Dimension(1:32) :: bpObj
    Type(potentialSearchType) :: searchObj
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal) :: tempDp
    Character(Len=16) :: inputStr
    Character(Len=2), Dimension(1:4) :: atomLabels
    Integer(kind=StandardInteger) :: i
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(coordsUnit(1:p_confs))
    Allocate(coords(1:p_confs))
    Allocate(nl(1:p_confs))
! Load Vars
    Call loadVars()
    print *,envVars%cwd
    print *,envVars%user
    print *,envVars%home
!    Type(isotopesObj) :: isotopesList
!    Type(isotopeObj) :: isotopeA
! Load isotope data
    !Call loadIsotopes(isotopesList)
    !isotopeA = SearchIsotopes(26,28,isotopesList)
    !Call printIsotopeDetails(isotopeA)
    inputStr = "0.0D0"
    tempDp = StrToDp(inputStr)
    inputStr = "1.3D0"
    tempDp = StrToDp(inputStr)
! Init Page
    Call initPage(mainPage)
! Start test
    Call addLinePage("Static Calculation Testing","T")
! Set up coords
    Call initUnitCoords(coordsUnit)     ! geom.f90
    Call initCoords(coords)             ! geom.f90
    coordsUnit%aLat = 4.04D0
    coordsUnit%xCopy = 4
    coordsUnit%yCopy = 4
    coordsUnit%zCopy = 4
    atomLabels = "AL"
    Call standardCoords("FCC", coordsUnit,1,atomLabels)  ! geom.f90
    Call standardCoords("FCC", coordsUnit,2,atomLabels)  ! geom.f90
    Call expandUnitCoords(coordsUnit, coords) ! geom.f90
! Init potential
    Call addLinePage("Init potential")
    Call initPotential(potential)        ! potentials.f90
! Load potential
    Call addLinePage("Load potential")
    Call loadPotential("test.pot", potential)  ! potentials.f90
! Assign atom IDs
    Call addLinePage("Assign IDs")
    Call atomLabelIDs(potential, coords)
! Print potential summary
    Call printPotentialSummary(potential)      ! potentials.f90
! Print
    Call printAtomLabelIDs(coords)

    Call fitStandardPotentials(potential)
! Output Potential
    !Call outputPotential(potential,envVars%cwd)
! Build neighbour list and make keys
    Call addLinePage("Build neighbour list")
    Call makeNL(nl, coords, 6.5D0)       ! geom.f90
    Call nlPotentialKeys(nl, potential)

    !Call printNLSummary(nl)
    Call calcEFS(nl, potential, 1)
    print *,"Energy EFS ",nl(1)%totalEnergy


    Call calcE(nl, potential, 1)
    print *,"Energy E   ",nl(1)%totalEnergy

    Call calcBP(bpObj,potential)

    Call HeatCoords(coords, 0.05D0, 1)
    Call makeNL(nl, coords, 6.5D0, 1)       ! geom.f90
    Call calcEFS(nl, potential, 1)

    Do i=1,4
      !print *,nl(1)%coordsMD(i,1),&
    !  nl(1)%coordsMD(i,2),&
      !nl(1)%coordsMD(i,3),&
    !  nl(1)%forcesMD(i,1),&
    !  nl(1)%forcesMD(i,2),&
    !  nl(1)%forcesMD(i,3),&
    !  nl(1)%atomEnergy(i,1),&
    !  nl(1)%atomEnergy(i,2),&
    !  nl(1)%atomEnergy(i,3),&
    !  nl(1)%electronDensity(i,1)
      print *, nl(1)%atomEnergy(i,3), nl(1)%electronDensity(i,1)
    End Do

    print *,"--------------"
    !Call printCoords(coords, 1)
    Call optGeom(coords, potential)





!----------------------------------
! Deallocate arrays on HEAP
!----------------------------------
    Deallocate(nl)
    Deallocate(coords)
    Deallocate(coordsUnit)

  End Subroutine testStaticCalc





  Subroutine testCombinations()
    Implicit None ! Force declaration of all variables
  ! Vars Private
    Call combinationsPrint(4,3)
    print *,""
    print *,""
    Call combinationsPrint(5,4)
    print *,""
    print *,""
  End Subroutine testCombinations



  Subroutine testLinearRegression()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Real(kind=DoubleReal), Dimension(1:4,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:4) :: Y
    Real(kind=DoubleReal), Dimension(1:4,1:1) :: X
    Real(kind=DoubleReal), Dimension(1:2) :: parameters

    dataPoints(1,1) = 1.0D0
    dataPoints(1,2) = 13.8D0
    dataPoints(2,1) = 2.0D0
    dataPoints(2,2) = 18.9D0
    dataPoints(3,1) = 3.0D0
    dataPoints(3,2) = 22.6781D0
    dataPoints(4,1) = 4.0D0
    dataPoints(4,2) = 25.2D0

    parameters = BestFitLine(dataPoints)
    print *,"Best fit test: "
    print *,parameters(1),parameters(2)

    dataPoints(1,1) = 1.0D0
    dataPoints(1,2) = 1.58907D0
    dataPoints(2,1) = 2.0D0
    dataPoints(2,2) = 0.78911D0
    dataPoints(3,1) = 3.0D0
    dataPoints(3,2) = 0.39186D0
    dataPoints(4,1) = 8.0D0
    dataPoints(4,2) = 0.011833D0
    parameters = BestFitExp(dataPoints)
    print *,"Best fit exp test: "
    print *,parameters(1),parameters(2)

    dataPoints(1,1) = 1.0D0
    dataPoints(1,2) = 0.537325D0
    dataPoints(2,1) = 2.0D0
    dataPoints(2,2) = 0.10423D0
    dataPoints(3,1) = 3.0D0
    dataPoints(3,2) = 0.0202186D0
    dataPoints(4,1) = 8.0D0
    dataPoints(4,2) = 5.5531D-6
    parameters = BestFitExp(dataPoints)
    print *,"Best fit exp^2 test: "
    print *,parameters(1),parameters(2)


  End Subroutine testLinearRegression

  Subroutine testMorseFit()
! Test static calculation interatomic potentials
    Implicit None ! Force declaration of all variables
! Vars Private
    Real(kind=DoubleReal), Dimension(1:8,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:8,1:2) :: dataPointsDYDX
    Real(kind=DoubleReal), Dimension(1:3) :: parameters
    Real(kind=DoubleReal), Dimension(1:2) :: parametersLJ

    dataPoints(1,1) = 1.0D0
    dataPoints(1,2) = 43.9370D0
    dataPoints(2,1) = 1.6D0
    dataPoints(2,2) = 8.9981D0
    dataPoints(3,1) = 2.0D0
    dataPoints(3,2) = 2.6781D0
    dataPoints(4,1) = 2.6D0
    dataPoints(4,2) = 0.080557D0
    dataPoints(5,1) = 3.0D0
    dataPoints(5,2) = -0.23856D0
    dataPoints(6,1) = 3.4D0
    dataPoints(6,2) = -0.26361D0
    dataPoints(7,1) = 5.0D0
    dataPoints(7,2) = -0.066057D0
    dataPoints(8,1) = 5.6D0
    dataPoints(8,2) = -0.0339987D0

    dataPointsDYDX(1,1) = 1.0D0
    dataPointsDYDX(1,2) = -111.0190D0
    dataPointsDYDX(2,1) = 1.6D0
    dataPointsDYDX(2,2) = -25.2747D0
    dataPointsDYDX(3,1) = 2.0D0
    dataPointsDYDX(3,2) = -8.94672D0
    dataPointsDYDX(4,1) = 2.6D0
    dataPointsDYDX(4,2) = -1.53451D0
    dataPointsDYDX(5,1) = 3.0D0
    dataPointsDYDX(5,2) = -0.28964D0
    dataPointsDYDX(6,1) = 3.4D0
    dataPointsDYDX(6,2) = 0.083474D0
    dataPointsDYDX(7,1) = 5.0D0
    dataPointsDYDX(7,2) = 0.07155D0
    dataPointsDYDX(8,1) = 5.6D0
    dataPointsDYDX(8,2) = 0.038265D0

    parameters = MorseFit(dataPoints)
    print *,parameters(1),parameters(2),parameters(3)


    parameters = MorseFit_Extended(dataPoints, dataPointsDYDX)
    print *,parameters(1),parameters(2),parameters(3)


    parametersLJ = LJFit(dataPoints)
    print *,parametersLJ(1),parametersLJ(2)



  End Subroutine testMorseFit




End Module testMod
