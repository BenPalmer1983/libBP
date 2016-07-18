! -------------------------------------------------
!  Include File:   Optimise geometry
!
! -------------------------------------------------
  Subroutine optGeom(coordsGeom, potential)
! Optimise Geometry
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Allocatable, Dimension(:) :: coordsGeom
    Type(potentialType) :: potential
! Vars Private
    Type(nlType_Opt), Allocatable, Dimension(:) :: nlGeom
    Integer(kind=StandardInteger) :: i, j, cKey
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(nlGeom(1:4))
! Initialise objects

    cKey = 1
    Call optGeom_BuildMatrix(coordsGeom, potential, cKey)


! Make neighbour list
    !Call makeNL_opt(nlGeom, coordsGeom, 6.5D0, 1)
    !Call nlPotentialKeys_opt(nlGeom, potential)

    !Call calcEF(nlGeom, potential, 1, 1, .true.)
    !print *, nlGeom(1)%atomEnergy(1,3), nlGeom(1)%electronDensity(1,1)

    !Call calcEF(nlGeom, potential, 1, 2, .true.)
    !print *, nlGeom(1)%atomEnergy(2,3), nlGeom(1)%electronDensity(2,1)

    !Call calcEF(nlGeom, potential, 1, 3, .true.)
    !print *, nlGeom(1)%atomEnergy(3,3), nlGeom(1)%electronDensity(3,1)

    !Call calcEF(nlGeom, potential, 1, 4, .true.)
    !print *, nlGeom(1)%atomEnergy(4,3), nlGeom(1)%electronDensity(4,1)




!
!    Do i=1,nlGeom(1)%coordsLength
!      print *,"Neighbours: ",nlGeom(1)%nlOverKey(i,0)
!    End Do



!----------------------------------
! Deallocate arrays on HEAP
!----------------------------------
    Deallocate(nlGeom)





  End Subroutine optGeom





  Subroutine optGeom_BuildMatrix(coordsGeom, potential, cKey)
! Optimise Geometry
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(coordsType), Allocatable, Dimension(:) :: coordsGeom
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars Private
    Type(nlType_Opt), Allocatable, Dimension(:) :: nlGeom
    Integer(kind=StandardInteger) :: i, j
!----------------------------------
! Allocate arrays on HEAP
!----------------------------------
    Allocate(nlGeom(1:1))
! Make neighbour list
    Call makeNL_opt(nlGeom, coordsGeom, 6.5D0, cKey)   ! geom.nl_opt.f90   - add option to set rVerlet, maybe in a data type

! need subroutine to vary atom position and update neighbour list of that individual atom



! Loop through each atom and vary position slightly



!----------------------------------
! Dellocate arrays on HEAP
!----------------------------------
    Deallocate(nlGeom)
  End Subroutine optGeom_BuildMatrix



!
