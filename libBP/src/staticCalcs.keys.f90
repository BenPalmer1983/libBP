! -------------------------------------------------
!  Include File:   Keys
!
! -------------------------------------------------

  Subroutine nlPotentialKeys(nl, potential, cKey_In)
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger), Optional :: cKey_In
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey, cKey_Loop
! Make key arrays
    cKey = 0
    If(Present(cKey_In))Then
      cKey = cKey_In
    End If
    If(cKey.eq.0)Then
      Do cKey_Loop=1,size(nl)
        Call nlPotentialKeys_Action(nl, potential, cKey_Loop)
      End Do
    Else
      Call nlPotentialKeys_Action(nl, potential, cKey)
    End If
  End Subroutine nlPotentialKeys


  Subroutine nlPotentialKeys_Action(nl, potential, cKey)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, fN, pKey, dKey, eKey
    Integer(kind=StandardInteger) :: pMax, deMax, arrayWidth
    Logical :: pairPot
!-------------------------
! Array Management
!-------------------------
    pMax = DoubleKey(nl(cKey)%atomID_Count+1, nl(cKey)%atomID_Count+1)
    deMax = nl(cKey)%atomID_Count+1
    arrayWidth = potential%fCount
! Deallocate
    If(Allocated(nl(cKey)%pairKeyArray))Then
      Deallocate(nl(cKey)%pairKeyArray)
    End If
    If(Allocated(nl(cKey)%densityKeyArray))Then
      Deallocate(nl(cKey)%densityKeyArray)
    End If
    If(Allocated(nl(cKey)%embeddingKeyArray))Then
      Deallocate(nl(cKey)%embeddingKeyArray)
    End If
! Allocate
    Allocate(nl(cKey)%pairKeyArray(1:pMax,0:arrayWidth))
    Allocate(nl(cKey)%densityKeyArray(1:deMax,0:arrayWidth,1:3))
    Allocate(nl(cKey)%embeddingKeyArray(1:deMax,0:arrayWidth,1:3))
!-------------------------
! Pair Potential Keys
!-------------------------
! Fill an array that links the double key of atom A and atom B to the (pair) potential key for quicker searching
    nl(cKey)%pairKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        pairPot = .true.
! If density function or embedding functional
        If(StrMatch(potential%fPotential(fN),"DENS").or.StrMatch(potential%fPotential(fN),"DDEN").or.&
        StrMatch(potential%fPotential(fN),"SDEN"))Then
          pairPot = .false.
        End If
        If(StrMatch(potential%fPotential(fN),"EMBE").or.StrMatch(potential%fPotential(fN),"SEMB").or.&
        StrMatch(potential%fPotential(fN),"DEMB"))Then
          pairPot = .false.
        End If
! If pair pot
        If(pairPot)Then
          pKey = DoubleKey(potential%atomID_A(fN), potential%atomID_B(fN))
          Do j=1,arrayWidth
            If(nl(cKey)%pairKeyArray(pKey,j).eq.0)Then
              nl(cKey)%pairKeyArray(pKey,j) = fN
              nl(cKey)%pairKeyArray(pKey,0) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!-------------------------
! Density Keys  (EAM, 2BMEAM)
!-------------------------
! Fill an array that links the key of each atom to the (density) potential key for quicker searching
    nl(cKey)%densityKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        dKey = potential%atomID_A(fN)
! DENS function
        If(StrMatch(potential%fPotential(fN),"DENS"))Then
          Do j=1,32
            If(nl(cKey)%densityKeyArray(dKey,j,1).eq.0)Then
              nl(cKey)%densityKeyArray(dKey,j,1) = fN
              nl(cKey)%densityKeyArray(dKey,0,1) = j
              Exit
            End If
          End Do
        End If
! DDEN function
        If(StrMatch(potential%fPotential(fN),"DDEN"))Then
          Do j=1,32
            If(nl(cKey)%densityKeyArray(dKey,j,2).eq.0)Then
              nl(cKey)%densityKeyArray(dKey,j,2) = fN
              nl(cKey)%densityKeyArray(dKey,0,2) = j
              Exit
            End If
          End Do
        End If
! SDEN function
        If(StrMatch(potential%fPotential(fN),"SDEN"))Then
          Do j=1,32
            If(nl(cKey)%densityKeyArray(dKey,j,3).eq.0)Then
              nl(cKey)%densityKeyArray(dKey,j,3) = fN
              nl(cKey)%densityKeyArray(dKey,0,3) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!-------------------------
! Embedding Keys  (EAM, 2BMEAM)
!-------------------------
! Fill an array that links the key of each atom to the (density) potential key for quicker searching
    nl(cKey)%embeddingKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        eKey = potential%atomID_A(fN)
! DENS function
        If(StrMatch(potential%fPotential(fN),"EMBE"))Then
          Do j=1,32
            If(nl(cKey)%embeddingKeyArray(eKey,j,1).eq.0)Then
              nl(cKey)%embeddingKeyArray(eKey,j,1) = fN
              nl(cKey)%embeddingKeyArray(eKey,0,1) = j
              Exit
            End If
          End Do
        End If
! DDEN function
        If(StrMatch(potential%fPotential(fN),"DEMB"))Then
          Do j=1,32
            If(nl(cKey)%embeddingKeyArray(eKey,j,2).eq.0)Then
              nl(cKey)%embeddingKeyArray(eKey,j,2) = fN
              nl(cKey)%embeddingKeyArray(eKey,0,2) = j
              Exit
            End If
          End Do
        End If
! SDEN function
        If(StrMatch(potential%fPotential(fN),"SEMB"))Then
          Do j=1,32
            If(nl(cKey)%embeddingKeyArray(eKey,j,3).eq.0)Then
              nl(cKey)%embeddingKeyArray(eKey,j,3) = fN
              nl(cKey)%embeddingKeyArray(eKey,0,3) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!---------------------------------
!
! Store Keys
!---------------------------------
    nl(cKey)%keysSet = .true.
  End Subroutine nlPotentialKeys_Action


  Subroutine nlPotentialKeys_opt(nl, potential, cKey_In)
! Vars:  In/Out
    Type(nlType_Opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger), Optional :: cKey_In
! Vars:  Private
    Integer(kind=StandardInteger) :: cKey, cKey_Loop
! Make key arrays
    cKey = 0
    If(Present(cKey_In))Then
      cKey = cKey_In
    End If
    If(cKey.eq.0)Then
      Do cKey_Loop=1,size(nl)
        Call nlPotentialKeysOpt_Action(nl, potential, cKey_Loop)
      End Do
    Else
      Call nlPotentialKeysOpt_Action(nl, potential, cKey)
    End If
  End Subroutine nlPotentialKeys_opt


  Subroutine nlPotentialKeysOpt_Action(nl, potential, cKey)
! Calculates energy/force/stress of a collection of atoms
!
    Implicit None   ! Force declaration of all variables
! Vars:  In/Out
    Type(nlType_Opt), Dimension(:) :: nl
    Type(potentialType) :: potential
    Integer(kind=StandardInteger) :: cKey
! Vars:  Private
    Integer(kind=StandardInteger) :: i, j, fN, pKey, dKey, eKey
    Integer(kind=StandardInteger) :: pMax, deMax, arrayWidth
    Logical :: pairPot
!-------------------------
! Array Management
!-------------------------
    pMax = DoubleKey(nl(cKey)%atomID_Count+1, nl(cKey)%atomID_Count+1)
    deMax = nl(cKey)%atomID_Count+1
    arrayWidth = potential%fCount
! Deallocate
    If(Allocated(nl(cKey)%pairKeyArray))Then
      Deallocate(nl(cKey)%pairKeyArray)
    End If
    If(Allocated(nl(cKey)%densityKeyArray))Then
      Deallocate(nl(cKey)%densityKeyArray)
    End If
    If(Allocated(nl(cKey)%embeddingKeyArray))Then
      Deallocate(nl(cKey)%embeddingKeyArray)
    End If
! Allocate
    Allocate(nl(cKey)%pairKeyArray(1:pMax,0:arrayWidth))
    Allocate(nl(cKey)%densityKeyArray(1:deMax,0:arrayWidth,1:3))
    Allocate(nl(cKey)%embeddingKeyArray(1:deMax,0:arrayWidth,1:3))
!-------------------------
! Pair Potential Keys
!-------------------------
! Fill an array that links the double key of atom A and atom B to the (pair) potential key for quicker searching
    nl(cKey)%pairKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        pairPot = .true.
! If density function or embedding functional
        If(StrMatch(potential%fPotential(fN),"DENS").or.StrMatch(potential%fPotential(fN),"DDEN").or.&
        StrMatch(potential%fPotential(fN),"SDEN"))Then
          pairPot = .false.
        End If
        If(StrMatch(potential%fPotential(fN),"EMBE").or.StrMatch(potential%fPotential(fN),"SEMB").or.&
        StrMatch(potential%fPotential(fN),"DEMB"))Then
          pairPot = .false.
        End If
! If pair pot
        If(pairPot)Then
          pKey = DoubleKey(potential%atomID_A(fN), potential%atomID_B(fN))
          Do j=1,32
            If(nl(cKey)%pairKeyArray(pKey,j).eq.0)Then
              nl(cKey)%pairKeyArray(pKey,j) = fN
              nl(cKey)%pairKeyArray(pKey,0) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!-------------------------
! Density Keys  (EAM, 2BMEAM)
!-------------------------
! Fill an array that links the key of each atom to the (density) potential key for quicker searching
    nl(cKey)%densityKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        dKey = potential%atomID_A(fN)
! DENS function
        If(StrMatch(potential%fPotential(fN),"DENS"))Then
          Do j=1,32
            If(nl(cKey)%densityKeyArray(dKey,j,1).eq.0)Then
              nl(cKey)%densityKeyArray(dKey,j,1) = fN
              nl(cKey)%densityKeyArray(dKey,0,1) = j
              Exit
            End If
          End Do
        End If
! DDEN function
        If(StrMatch(potential%fPotential(fN),"DDEN"))Then
          Do j=1,32
            If(nl(cKey)%densityKeyArray(dKey,j,2).eq.0)Then
              nl(cKey)%densityKeyArray(dKey,j,2) = fN
              nl(cKey)%densityKeyArray(dKey,0,2) = j
              Exit
            End If
          End Do
        End If
! SDEN function
        If(StrMatch(potential%fPotential(fN),"SDEN"))Then
          Do j=1,32
            If(nl(cKey)%densityKeyArray(dKey,j,3).eq.0)Then
              nl(cKey)%densityKeyArray(dKey,j,3) = fN
              nl(cKey)%densityKeyArray(dKey,0,3) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!-------------------------
! Embedding Keys  (EAM, 2BMEAM)
!-------------------------
! Fill an array that links the key of each atom to the (density) potential key for quicker searching
    nl(cKey)%embeddingKeyArray = 0
! Loop through potentials
    Do fN=1,potential%fCount
      If(potential%switchOn(fN))Then
        eKey = potential%atomID_A(fN)
! DENS function
        If(StrMatch(potential%fPotential(fN),"EMBE"))Then
          Do j=1,32
            If(nl(cKey)%embeddingKeyArray(eKey,j,1).eq.0)Then
              nl(cKey)%embeddingKeyArray(eKey,j,1) = fN
              nl(cKey)%embeddingKeyArray(eKey,0,1) = j
              Exit
            End If
          End Do
        End If
! DDEN function
        If(StrMatch(potential%fPotential(fN),"DEMB"))Then
          Do j=1,32
            If(nl(cKey)%embeddingKeyArray(eKey,j,2).eq.0)Then
              nl(cKey)%embeddingKeyArray(eKey,j,2) = fN
              nl(cKey)%embeddingKeyArray(eKey,0,2) = j
              Exit
            End If
          End Do
        End If
! SDEN function
        If(StrMatch(potential%fPotential(fN),"SEMB"))Then
          Do j=1,32
            If(nl(cKey)%embeddingKeyArray(eKey,j,3).eq.0)Then
              nl(cKey)%embeddingKeyArray(eKey,j,3) = fN
              nl(cKey)%embeddingKeyArray(eKey,0,3) = j
              Exit
            End If
          End Do
        End If
      End If
    End Do
!---------------------------------
!
! Store Keys
!---------------------------------
    nl(cKey)%keysSet = .true.


  End Subroutine nlPotentialKeysOpt_Action



!
