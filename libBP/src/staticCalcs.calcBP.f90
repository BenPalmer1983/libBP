! -------------------------------------------------
!  Include File:   Calc bulk properties
!
! -------------------------------------------------
  Subroutine calcBP(bpObj,potential)
! Calculate bulk properties
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oBulkProperty), Dimension(:) :: bpObj
    Type(potentialType) :: potential
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, cKey
! Vars Private
    Type(coordsUnitType), Allocatable, Dimension(:) :: coordsUnitBP
    Type(coordsType), Allocatable, Dimension(:) :: coordsBP
    Type(nlType), Allocatable, Dimension(:) :: nlBP
    Character(Len=16), Dimension(1:4) :: atomLabels
    Integer(kind=StandardInteger), Dimension(1:4) :: atomIDs
! EOS
    Real(kind=DoubleReal), Dimension(1:5,1:2) :: points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
    Real(kind=DoubleReal), Dimension(1:3) :: polyCoefficients
    Real(kind=DoubleReal) :: optAlat
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortionArray
! Elastic Constants
    Real(kind=DoubleReal) :: sigma
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(coordsUnitBP(1:4))
    Allocate(coordsBP(1:4))
    Allocate(nlBP(1:4))
! Initialise objects
    Call initUnitCoords(coordsUnitBP)     ! geom.f90
    Call initCoords(coordsBP)             ! geom.f90
! Loop through atom IDs
    cKey = 0
    Do i=1,potential%atomID_Count
      atomLabels = potential%atomIDs(i)  ! Set atomLabels array and atomIDs array
      atomIDs = i
!----------------------------------
! FCC Equation of State
!----------------------------------
      cKey = cKey + 1
      coordsUnitBP%aLat = 5.00D0
      coordsUnitBP%xCopy = 3
      coordsUnitBP%yCopy = 3
      coordsUnitBP%zCopy = 3
      Call standardCoords("FCC", coordsUnitBP, cKey, atomLabels, atomIDs)  ! geom.f90
      Call expandUnitCoords(coordsUnitBP, coordsBP)   ! geom.f90
! Find e0 and aLat
      Do j=1,5
        If(j.eq.1)Then
          Call makeNL(nlBP, coordsBP, 6.5D0)       ! geom.f90
          Call nlPotentialKeys(nlBP, potential)    ! staticCalcs.keys.f90
          Call calcE(nlBP, potential, cKey)
        Else
          nlBP(cKey)%aLat = nlBP(cKey)%aLat - 1.0D0
          Call updateNL(nlBP, cKey)       ! geom.f90
          Call calcE(nlBP, potential, cKey)
        End If
        points(6-j,1) = (nlBP(cKey)%aLat)**3/nlBP(cKey)%coordsLength
        points(6-j,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
      End Do
      coefficients = BirchMurnFit(points)
      optAlat = (coefficients(2)*nlBP(cKey)%coordsLength)**(1.0D0/3.0D0)
      coordsUnitBP%aLat = optAlat/3.0D0
      coordsUnitBP%xCopy = 3
      coordsUnitBP%yCopy = 3
      coordsUnitBP%zCopy = 3
      Call standardCoords("FCC", coordsUnitBP, cKey, atomLabels, atomIDs)  ! geom.f90
      Call expandUnitCoords(coordsUnitBP, coordsBP)   ! geom.f90
! Find e0 and aLat (refine)
      Do j=1,5
        If(j.eq.1)Then
          Call makeNL(nlBP, coordsBP, 6.5D0)       ! geom.f90
          Call nlPotentialKeys(nlBP, potential)    ! staticCalcs.keys.f90
          Call calcE(nlBP, potential, cKey)
        Else
          nlBP(cKey)%aLat = optAlat + (3-j) * 0.25D0
          Call updateNL(nlBP, cKey)       ! geom.f90
          Call calcE(nlBP, potential, cKey)
        End If
        points(6-j,1) = (nlBP(cKey)%aLat)**3/nlBP(cKey)%coordsLength
        points(6-j,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
      End Do
      coefficients = BirchMurnFit(points)
      optAlat = (coefficients(2)*nlBP(cKey)%coordsLength)**(1.0D0/3.0D0)
      bpObj(i)%fccAlat = (coefficients(2))**(1.0D0/3.0D0)
      bpObj(i)%fccV0 = (coefficients(2))
      bpObj(i)%fccE0 = coefficients(1)
      bpObj(i)%fccB0 = coefficients(3)
      bpObj(i)%fccB0_GPA = UnitConvert(bpObj(i)%fccB0, "EVAN3", "GPA")   ! Convert to GPA
      bpObj(i)%fccBp0 = coefficients(4)
!----------------------------------
! C11, C12
!----------------------------------
      points(1,1) = 0.0D0
      points(1,2) = bpObj(i)%fccE0
      nlBP(cKey)%aLat = optAlat
      Call makeNL(nlBP, coordsBP, 6.5D0)
      Call nlPotentialKeys(nlBP, potential)    ! staticCalcs.keys.f90
      sigma = 0.025D0
      Do j=1,4
        distortionArray = Tensor_Orthorhombic(j*sigma)
        Call updateNL(nlBP, cKey, distortionArray)       ! geom.f90
        Call calcE(nlBP, potential, cKey)
        points(j+1,1) = j*sigma
        points(j+1,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
      End Do
      polyCoefficients = PolyFit(points,2)
      bpObj(i)%fccC11 = (2.0D0*polyCoefficients(3))/(3.0D0*bpObj(i)%fccV0)+bpObj(i)%fccB0
      bpObj(i)%fccC11_GPA = UnitConvert(bpObj(i)%fccC11, "EVAN3", "GPA")   ! Convert to GPA
      bpObj(i)%fccC12 = (3.0D0*bpObj(i)%fccB0-bpObj(i)%fccC11)/2.0D0
      bpObj(i)%fccC12_GPA = UnitConvert(bpObj(i)%fccC12, "EVAN3", "GPA")   ! Convert to GPA
!----------------------------------
! C44
!----------------------------------
      points(1,1) = 0.0D0
      points(1,2) = bpObj(i)%fccE0
      nlBP(cKey)%aLat = optAlat
      Call makeNL(nlBP, coordsBP, 6.5D0)
      Call nlPotentialKeys(nlBP, potential)    ! staticCalcs.keys.f90
      sigma = 0.025D0
      print *,points(1,1),points(1,2)
      Do j=1,4
        distortionArray = Tensor_Tetragonal(j*sigma)
        Call updateNL(nlBP, cKey, distortionArray)       ! geom.f90
        Call calcE(nlBP, potential, cKey)
        points(j+1,1) = j*sigma
        points(j+1,2) = (nlBP(cKey)%totalEnergy/nlBP(cKey)%coordsLength)
        print *,points(j+1,1),points(j+1,2)
      End Do
      polyCoefficients = PolyFit(points,2)
      bpObj(i)%fccC44 = (2.0D0*polyCoefficients(3))/(bpObj(i)%fccV0)
      bpObj(i)%fccC44_GPA = UnitConvert(bpObj(i)%fccC44, "EVAN3", "GPA")   ! Convert to GPA



! Convert to GPA


      print *,"FCC E0 ",bpObj(i)%fccE0
      print *,"FCC aLat ",bpObj(i)%fccAlat
      print *,"FCC v0 ",bpObj(i)%fccV0
      print *,"FCC B0 ",bpObj(i)%fccB0
      print *,"FCC B0 (GPA) ",bpObj(i)%fccB0_GPA
      print *,"FCC C11 (GPA) ",bpObj(i)%fccC11_GPA
      print *,"FCC C12 (GPA) ",bpObj(i)%fccC12_GPA
      print *,"FCC C44 (GPA) ",bpObj(i)%fccC44_GPA



! BCC
      !cKey = cKey + 1
      !coordsUnitBP%aLat = 2.50D0
      !coordsUnitBP%xCopy = 4
      !coordsUnitBP%yCopy = 4
      !coordsUnitBP%zCopy = 4
      !Call standardCoords("BCC", coordsUnitBP, cKey, atomLabels, atomIDs)  ! geom.f90
      !Call expandUnitCoords(coordsUnitBP, coordsBP)   ! geom.f90
      !Call makeNL(nlBP, coordsBP, 6.5D0, cKey)       ! geom.f90

    End Do




!----------------------------------
! Deallocate arrays on HEAP
!----------------------------------
    Deallocate(coordsUnitBP)
    Deallocate(coordsBP)
    Deallocate(nlBP)





End Subroutine calcBP



Subroutine initBP(bpObj)
! Initialise bulk properties object
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(oBulkProperty) :: bpObj
! Vars:  Private
! Initialise
    bpObj%bccAlat = 0.0D0
    bpObj%bccE0 = 0.0D0

    bpObj%fccAlat = 0.0D0
    bpObj%fccE0 = 0.0D0



End Subroutine initBP
