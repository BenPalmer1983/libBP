! --------------------------------------------------------------!
! Keys module
! keysMod
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calculates the forces between atoms
! Calculates the energies of atoms/total energy
! Calculates stresses
!
! ----------------------------------------
! Updated: 15th June 2016
! ----------------------------------------

Module keysMod
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: DoubleKey
  Public :: TripleKey
  Public :: IsotopeKey
! Interfaces
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------!
! Vector Functions
! ------------------------------------------------------------------------!

  Function DoubleKey(keyA, keyB) RESULT (keyU)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
    Implicit None ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger) :: keyA, keyB
! Vars:  Out
    Integer(kind=StandardInteger) :: keyU
! Vars:  Private
    Integer(kind=StandardInteger) :: minKey, maxKey
! Min/Max
    minKey = keyA
    maxKey = keyB
    If(keyA.gt.keyB)Then
      maxKey = keyA
      minKey = keyB
    End If
    keyU = (maxKey*(maxKey-1))/2+minKey
  End Function DoubleKey


  Function TripleKey(keyA, keyB, keyC) RESULT (keyU)
! Calculates unique key for three keys (order of keys NOT important)
! (A,B,C) = (A,C,B) = (B,A,C) = (B,C,A) = (C,A,B) = (C,B,A)
    Implicit None ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger) :: keyA, keyB, keyC
! Vars:  Out
    Integer(kind=StandardInteger) :: keyU
! Vars:  Private
    Integer(kind=StandardInteger) :: minKey, keyAA, keyBB, keyTemp
! Max/Mid/Min
    minKey = keyA
    keyAA = keyB
    keyBB = keyC
    If(keyB.lt.minKey)Then
      minKey = keyB
      keyAA = keyA
      keyBB = keyC
    End If
    If(keyC.lt.minKey)Then
      minKey = keyC
      keyAA = keyA
      keyBB = keyB
    End If
    keyTemp = DoubleKey(keyAA, keyBB)
    keyU = DoubleKey(keyTemp, minKey)
  End Function TripleKey


  Function IsotopeKey(protons,neutrons) RESULT (keyI)
! Calculates unique key for proton/neutron - order important
! maxP  140
! maxN  200
! max key 57830
    Implicit None ! Force declaration of all variables
! Vars:  In
    Integer(kind=StandardInteger) :: protons, neutrons
! Vars:  Out
    Integer(kind=StandardInteger) :: keyI
! Vars:  Private
    Integer(kind=StandardInteger) :: keyP, keyN, maxN
! Calculate key
    maxN = 200
    keyP = maxN + protons
    keyN = neutrons
    keyI = DoubleKey(keyP, keyN)

  End Function IsotopeKey

End Module keysMod

!-----------------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------------------
