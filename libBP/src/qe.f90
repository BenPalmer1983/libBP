! --------------------------------------------------------------!
! QE
! qeTypes
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
!
! ----------------------------------------
! Updated: 4th March 2016
! ----------------------------------------

Module qeTypes
! Setup Modules
  Use kinds

  Type :: qeInputData
! CONTROL
    Character(len=12) ::               calculation="scf         "
    Character(len=16) ::               title="                "
    Character(len=8) ::                verbosity="low     "
    Character(len=12) ::               restart_mode="from_scratch"
    Logical ::                         wf_collect=.false.
    Integer(kind=StandardInteger) ::   nstep=1
    Integer(kind=StandardInteger) ::   iprint=0         ! Don't output
    Logical ::                         tstress=.true.
    Logical ::                         tprnfor=.true.
    Real(kind=DoubleReal) ::           dt=20.D0         ! Don't output
    Character(len=48) ::               outdir="                                                "
    Character(len=48) ::               wfcdir="                                                "
    Character(len=32) ::               prefix="                                "
    Logical ::                         lkpoint_dir=.true.
    Real(kind=DoubleReal) ::           max_seconds=1.0D7
    Real(kind=DoubleReal) ::           etot_conv_thr=1.0D-4
    Real(kind=DoubleReal) ::           forc_conv_thr=1.0D-3
    Character(len=8) ::                disk_io="low     "
    Character(len=48) ::               pseudo_dir="                                                "
    Logical ::                         tefield=.false.
    Logical ::                         dipfield=.false.
    Logical ::                         lelfield=.false.
    Integer(kind=StandardInteger) ::   nberrycyc=1      ! Output if lelfield = .true.
    Logical ::                         lorbm=.false.
    Logical ::                         lberry=.false.
    Integer(kind=StandardInteger) ::   gdir=1
    Integer(kind=StandardInteger) ::   nppstr=1
    Logical ::                         lfcpopt=.false.
! SYSTEM
    Integer(kind=StandardInteger) ::   ibrav=14
    Real(kind=DoubleReal) ::           celldm1 = 1.0D0
    Real(kind=DoubleReal) ::           celldm2 = 1.0D0
    Real(kind=DoubleReal) ::           celldm3 = 1.0D0
    Real(kind=DoubleReal) ::           celldm4 = 0.0D0
    Real(kind=DoubleReal) ::           celldm5 = 0.0D0
    Real(kind=DoubleReal) ::           celldm6 = 0.0D0
    Integer(kind=StandardInteger) ::   nat=0
    Integer(kind=StandardInteger) ::   ntyp=0
    Integer(kind=StandardInteger) ::   nbnd=0
    Real(kind=DoubleReal) ::           tot_charge = 0.0D0
    Real(kind=DoubleReal) ::           tot_magnetization = 0.0D0
    Real(kind=DoubleReal), Dimension(1:2048) :: starting_magnetization = 0.0D0
    Real(kind=DoubleReal) ::           ecutwfc = 0.0D0
    Real(kind=DoubleReal) ::           ecutrho = 0.0D0
    Character(len=12) ::               occupations = "smearing    "
    Character(len=12) ::               smearing = "mv          "
    Real(kind=DoubleReal) ::           degauss = 0.05D0
! ELECTRONS
    Integer(kind=StandardInteger) ::   electron_maxstep = 100
    Real(kind=DoubleReal) ::           conv_thr = 1.0D-6
    Character(len=8) ::                mixing_mode="plain"
    Real(kind=DoubleReal) ::           mixing_beta=0.7D0
    Integer(kind=StandardInteger) ::   mixing_ndim=10
    Character(len=8) ::                diagonalization="david"
! IONS
    Character(len=16) ::               ion_dynamics="bfgs"
! CELL
    Character(len=16) ::               cell_dynamics="bfgs"
    Real(kind=DoubleReal) ::           press=0.0D0
    Real(kind=DoubleReal) ::           cell_factor=1.2D0
! ATOMIC SPECIES
    Character(len=16), Dimension(1:64) :: atomic_species_label
    Real(kind=DoubleReal), Dimension(1:64) :: atomic_species_mass
    Character(len=64), Dimension(1:64) :: atomic_species_pp
! ATOMIC POSITIONS
    Character(len=16) :: atomic_position_type
    Character(len=32), Dimension(1:4096) :: atomic_position_label
    Real(kind=DoubleReal), Dimension(1:4096,1:3) :: atomic_position_coords
! K-Points
    Character(len=16) :: kpoints_type
    Integer(kind=StandardInteger), Dimension(1:3) :: kpoints_grid
    Integer(kind=StandardInteger), Dimension(1:3) :: kpoints_offset

  End Type


End Module qeTypes

Module qe
! Setup Modules
  Use kinds
  Use strings
  Use general
  Use qeTypes
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
!
! Public Subroutines
  Public :: readPwscfInput
  Public :: writePwscfInput



  Contains
! ---------------------------------------------------------------------------------------------------

  Subroutine readPwscfInput(filePath,qeInObj)
! Reset data objects for a plot
    Implicit None  ! Force declaration of all variables
! Input/Output
    Character(*) :: filePath
    Type(qeInputData) :: qeInObj
! Private
    Integer(kind=StandardInteger) :: i, j, dataRows
    Character(len=255), Dimension(1:256) :: pwscfInputData
    Character(len=255) :: tempRow, tempRowU
    Character(len=64) :: bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
    print *,filePath
! Read in pwscf input file
    Call readFile(filePath, pwscfInputData, dataRows)
!    Do i=1,dataRows
    i = 0
    Do While (i.le.dataRows)
      i = i + 1
      tempRow = RemoveTrailing(pwscfInputData(i),",")
      tempRow = StrReplace(tempRow,"="," ")
      tempRow = SingleSpaces(tempRow)
      tempRowU = StrToUpper(tempRow)
!---------------------------------------------------------
! CONTROL
!---------------------------------------------------------
      If(tempRowU(1:11).eq."CALCULATION")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%calculation = trim(bufferB)
      End If
      If(tempRowU(1:5).eq."TITLE")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%title = trim(bufferB)
      End If
      If(tempRowU(1:9).eq."VERBOSITY")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%verbosity = trim(bufferB)
      End If
      If(tempRowU(1:12).eq."RESTART_MODE")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%restart_mode = trim(bufferB)
      End If
      If(tempRowU(1:10).eq."WF_COLLECT")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%wf_collect
      End If
      If(tempRowU(1:5).eq."NSTEP")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%nstep
      End If
      If(tempRowU(1:6).eq."IPRINT")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%iprint
      End If
      If(tempRowU(1:7).eq."TSTRESS")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%tstress
      End If
      If(tempRowU(1:7).eq."TPRNTOR")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%tprnfor
      End If
      If(tempRowU(1:3).eq."DT ")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%dt
      End If
      If(tempRowU(1:6).eq."OUTDIR")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%outdir = trim(bufferB)
      End If
      If(tempRowU(1:6).eq."WFCDIR")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%wfcdir = trim(bufferB)
      End If
      If(tempRowU(1:6).eq."PREFIX")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%prefix = trim(bufferB)
      End If
      If(tempRowU(1:13).eq."ETOT_CONV_THR")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%etot_conv_thr
      End If
      If(tempRowU(1:13).eq."FORC_CONV_THR")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%forc_conv_thr
      End If
      If(tempRowU(1:7).eq."DISK_IO")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%disk_io = trim(bufferB)
      End If
      If(tempRowU(1:10).eq."PSEUDO_DIR")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%pseudo_dir = trim(bufferB)
      End If
!---------------------------------------------------------
! SYSTEM
!---------------------------------------------------------
      If(tempRowU(1:5).eq."IBRAV")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%ibrav
      End If
      If(tempRowU(1:9).eq."CELLDM(1)")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%celldm1
      End If
      If(tempRowU(1:9).eq."CELLDM(2)")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%celldm2
      End If
      If(tempRowU(1:9).eq."CELLDM(3)")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%celldm3
      End If
      If(tempRowU(1:9).eq."CELLDM(4)")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%celldm4
      End If
      If(tempRowU(1:9).eq."CELLDM(5)")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%celldm5
      End If
      If(tempRowU(1:9).eq."CELLDM(6)")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%celldm6
      End If
      If(tempRowU(1:3).eq."NAT")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%nat
      End If
      If(tempRowU(1:4).eq."NTYP")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%ntyp
      End If
      If(tempRowU(1:4).eq."NBND")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%nbnd
      End If
      If(tempRowU(1:10).eq."TOT_CHARGE")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%tot_charge
      End If
      If(tempRowU(1:17).eq."TOT_MAGNETIZATION")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%tot_magnetization
      End If
      If(tempRowU(1:7).eq."ECUTWFC")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%ecutwfc
      End If
      If(tempRowU(1:7).eq."ECUTRHO")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%ecutrho
      End If
      If(tempRowU(1:11).eq."OCCUPATIONS")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%occupations = trim(bufferB)
      End If
      If(tempRowU(1:8).eq."SMEARING")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%smearing = trim(bufferB)
      End If
      If(tempRowU(1:7).eq."DEGAUSS")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%degauss
      End If
!---------------------------------------------------------
! ELECTRONS
!---------------------------------------------------------
      If(tempRowU(1:16).eq."ELECTRON_MAXSTEP")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%electron_maxstep
      End If
      If(tempRowU(1:8).eq."CONV_THR")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%conv_thr
      End If
      If(tempRowU(1:11).eq."MIXING_MODE")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%mixing_mode = trim(bufferB)
      End If
      If(tempRowU(1:11).eq."MIXING_BETA")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%mixing_beta
      End If
      If(tempRowU(1:11).eq."MIXING_NDIM")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%mixing_ndim
      End If
      If(tempRowU(1:15).eq."DIAGONALIZATION")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%diagonalization = trim(bufferB)
      End If
!---------------------------------------------------------
! IONS
!---------------------------------------------------------
      If(tempRowU(1:12).eq."ION_DYNAMICS")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%ion_dynamics = trim(bufferB)
      End If
!---------------------------------------------------------
! CELL
!---------------------------------------------------------
      If(tempRowU(1:13).eq."CELL_DYNAMICS")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%cell_dynamics
      End If
      If(tempRowU(1:5).eq."PRESS")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%press
      End If
      If(tempRowU(1:11).eq."CELL_FACTOR")Then
        Read(tempRow,*) bufferA, bufferB
        Read(bufferB,*) qeInObj%cell_factor
      End If
!---------------------------------------------------------
! ATOMIC SPECIES
!---------------------------------------------------------
      If(tempRowU(1:14).eq."ATOMIC_SPECIES")Then
        Do j = 1,qeInObj%ntyp
          i = i + 1
          tempRow = RemoveTrailing(pwscfInputData(i),",")
          tempRow = StrReplace(tempRow,"="," ")
          tempRow = SingleSpaces(tempRow)
          tempRowU = StrToUpper(tempRow)
          Read(tempRow,*) bufferA, bufferB, bufferC
          qeInObj%atomic_species_label(j) = trim(bufferA)
          Read(bufferB,*) qeInObj%atomic_species_mass(j)
          qeInObj%atomic_species_pp(j) = trim(bufferC)
        End Do
      End If
!---------------------------------------------------------
! ATOMIC POSITIONS
!---------------------------------------------------------
      If(tempRowU(1:16).eq."ATOMIC_POSITIONS")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%atomic_position_type = trim(bufferB)
        Do j = 1,qeInObj%nat
          i = i + 1
          tempRow = RemoveTrailing(pwscfInputData(i),",")
          tempRow = StrReplace(tempRow,"="," ")
          tempRow = SingleSpaces(tempRow)
          tempRowU = StrToUpper(tempRow)
          Read(tempRow,*) bufferA, bufferB, bufferC, bufferD
          qeInObj%atomic_position_label(j) = trim(bufferA)
          Read(bufferB,*) qeInObj%atomic_position_coords(j,1)
          Read(bufferC,*) qeInObj%atomic_position_coords(j,2)
          Read(bufferD,*) qeInObj%atomic_position_coords(j,3)
        End Do
      End If
!---------------------------------------------------------
! K-POINTS
!---------------------------------------------------------
      If(tempRowU(1:8).eq."K_POINTS")Then
        Read(tempRow,*) bufferA, bufferB
        qeInObj%kpoints_type = trim(bufferB)
        bufferB = StrToUpper(bufferB)
        If(bufferB(1:9).eq."AUTOMATIC")Then
          i = i + 1
          tempRow = RemoveTrailing(pwscfInputData(i),",")
          tempRow = StrReplace(tempRow,"="," ")
          tempRow = SingleSpaces(tempRow)
          Read(tempRow,*) bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
          Read(bufferA,*) qeInObj%kpoints_grid(1)
          Read(bufferB,*) qeInObj%kpoints_grid(2)
          Read(bufferC,*) qeInObj%kpoints_grid(3)
          Read(bufferD,*) qeInObj%kpoints_offset(1)
          Read(bufferE,*) qeInObj%kpoints_offset(2)
          Read(bufferF,*) qeInObj%kpoints_offset(3)
        End If
      End If
    End Do
  End Subroutine readPwscfInput




  Subroutine writePwscfInput(filePath,qeInObj)
! Reset data objects for a plot
    Implicit None  ! Force declaration of all variables
! Input/Output
    Character(*) :: filePath
    Type(qeInputData) :: qeInObj
    Character(Len=128) :: tempRow, testRow, buffer
    Integer(kind=StandardInteger) :: i

    Open(UNIT=7777,FILE=Trim(filePath))

!---------------------------------------------------------
! CONTROL
!---------------------------------------------------------
    write(7777,"(A8)") "&CONTROL"
! calculation
    tempRow = BlankString(tempRow)
    tempRow = "calculation = "//char(34)//trim(adjustl(qeInObj%calculation))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! title
    tempRow = BlankString(tempRow)
    testRow = trim(adjustl(qeInObj%title))
    If(testRow(1:1).ne." ")Then
      tempRow = "title = "//char(34)//trim(adjustl(qeInObj%title))//char(34)
      write(7777,"(A)") trim(adjustl(tempRow))
    End If
! verbosity
    tempRow = BlankString(tempRow)
    tempRow = "verbosity = "//char(34)//trim(adjustl(qeInObj%verbosity))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! restart mode
    tempRow = BlankString(tempRow)
    tempRow = "restart_mode = "//char(34)//trim(adjustl(qeInObj%restart_mode))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! wf_collect
    tempRow = BlankString(tempRow)
    If(qeInObj%wf_collect)Then
      tempRow = "wf_collect = .true."
      write(7777,"(A)") trim(adjustl(tempRow))
    End If
! nstep
    tempRow = BlankString(tempRow)
    tempRow = "nstep = "//char(34)//trim(adjustl(IntToStr(qeInObj%nstep)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! tstress
    tempRow = BlankString(tempRow)
    If(qeInObj%tstress)Then
      tempRow = "tstress = .true."
      write(7777,"(A)") trim(adjustl(tempRow))
    End If
! tprnfor
    tempRow = BlankString(tempRow)
    If(qeInObj%tprnfor)Then
      tempRow = "tprnfor = .true."
      write(7777,"(A)") trim(adjustl(tempRow))
    End If
! outdir
    tempRow = BlankString(tempRow)
    testRow = trim(adjustl(qeInObj%outdir))
    If(testRow(1:1).ne." ")Then
      tempRow = "outdir = "//char(34)//trim(adjustl(qeInObj%outdir))//char(34)
      write(7777,"(A)") trim(adjustl(tempRow))
    End If
! wfcdir
    tempRow = BlankString(tempRow)
    testRow = trim(adjustl(qeInObj%wfcdir))
    If(testRow(1:1).ne." ")Then
      tempRow = "wfcdir = "//char(34)//trim(adjustl(qeInObj%wfcdir))//char(34)
      write(7777,"(A)") trim(adjustl(tempRow))
    End If
! prefix
    tempRow = BlankString(tempRow)
    testRow = trim(adjustl(qeInObj%wfcdir))
    If(testRow(1:1).ne." ")Then
      tempRow = "prefix = "//char(34)//trim(adjustl(qeInObj%prefix))//char(34)
      write(7777,"(A)") trim(adjustl(tempRow))
    End If
! etot_conv_thr
    tempRow = BlankString(tempRow)
    tempRow = "etot_conv_thr = "//char(34)//trim(adjustl(DpToStr(qeInObj%etot_conv_thr)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! forc_conv_thr
    tempRow = BlankString(tempRow)
    tempRow = "forc_conv_thr = "//char(34)//trim(adjustl(DpToStr(qeInObj%forc_conv_thr)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! disk_io
    tempRow = BlankString(tempRow)
    tempRow = "disk_io = "//char(34)//trim(adjustl(qeInObj%disk_io))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! pseudo_dir
    tempRow = BlankString(tempRow)
    tempRow = "pseudo_dir = "//char(34)//trim(adjustl(qeInObj%pseudo_dir))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! End of section
    write(7777,"(A1)") "/"
!---------------------------------------------------------
! SYSTEM
!---------------------------------------------------------
    write(7777,"(A7)") "&SYSTEM"
! ibrav
    tempRow = BlankString(tempRow)
    tempRow = "ibrav = "//char(34)//trim(adjustl(IntToStr(qeInObj%ibrav)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! celldm1
    tempRow = BlankString(tempRow)
    tempRow = "celldm(1) = "//char(34)//trim(adjustl(DpToStr(qeInObj%celldm1)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! celldm2
    tempRow = BlankString(tempRow)
    tempRow = "celldm(2) = "//char(34)//trim(adjustl(DpToStr(qeInObj%celldm2)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! celldm3
    tempRow = BlankString(tempRow)
    tempRow = "celldm(3) = "//char(34)//trim(adjustl(DpToStr(qeInObj%celldm3)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! celldm4
    tempRow = BlankString(tempRow)
    tempRow = "celldm(4) = "//char(34)//trim(adjustl(DpToStr(qeInObj%celldm4)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! celldm5
    tempRow = BlankString(tempRow)
    tempRow = "celldm(5) = "//char(34)//trim(adjustl(DpToStr(qeInObj%celldm5)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! celldm6
    tempRow = BlankString(tempRow)
    tempRow = "celldm(6) = "//char(34)//trim(adjustl(DpToStr(qeInObj%celldm6)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! nat
    tempRow = BlankString(tempRow)
    tempRow = "nat = "//char(34)//trim(adjustl(IntToStr(qeInObj%nat)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! ntyp
    tempRow = BlankString(tempRow)
    tempRow = "ntyp = "//char(34)//trim(adjustl(IntToStr(qeInObj%ntyp)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! nbnd
    If(qeInObj%nbnd.gt.0)Then
      tempRow = BlankString(tempRow)
      tempRow = "nbnd = "//char(34)//trim(adjustl(IntToStr(qeInObj%nbnd)))//char(34)
      write(7777,"(A)") trim(adjustl(tempRow))
    End If
! starting_magnetization
! ecutwfc
    tempRow = BlankString(tempRow)
    tempRow = "ecutwfc = "//char(34)//trim(adjustl(DpToStr(qeInObj%ecutwfc)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! ecutrho
    tempRow = BlankString(tempRow)
    tempRow = "ecutrho = "//char(34)//trim(adjustl(DpToStr(qeInObj%ecutrho)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! occupations
    tempRow = BlankString(tempRow)
    tempRow = "occupations = "//char(34)//trim(adjustl(qeInObj%occupations))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! smearing
    tempRow = BlankString(tempRow)
    tempRow = "smearing = "//char(34)//trim(adjustl(qeInObj%smearing))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! degauss
    tempRow = BlankString(tempRow)
    tempRow = "degauss = "//char(34)//trim(adjustl(DpToStr(qeInObj%degauss)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! End of section
    write(7777,"(A1)") "/"
!---------------------------------------------------------
! ELECTRONS
!---------------------------------------------------------
    write(7777,"(A10)") "&ELECTRONS"
! electron_maxstep
    tempRow = BlankString(tempRow)
    tempRow = "electron_maxstep = "//char(34)//trim(adjustl(IntToStr(qeInObj%electron_maxstep)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! conv_thr
    tempRow = BlankString(tempRow)
    tempRow = "conv_thr = "//char(34)//trim(adjustl(DpToStr(qeInObj%conv_thr)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! mixing_mode
    tempRow = BlankString(tempRow)
    tempRow = "mixing_mode = "//char(34)//trim(adjustl(qeInObj%mixing_mode))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! mixing_beta
    tempRow = BlankString(tempRow)
    tempRow = "mixing_beta = "//char(34)//trim(adjustl(DpToStr(qeInObj%mixing_beta)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! mixing_ndim
    tempRow = BlankString(tempRow)
    tempRow = "mixing_ndim = "//char(34)//trim(adjustl(IntToStr(qeInObj%mixing_ndim)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! diagonalization
    tempRow = BlankString(tempRow)
    tempRow = "diagonalization = "//char(34)//trim(adjustl(qeInObj%diagonalization))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! End of section
    write(7777,"(A1)") "/"
!---------------------------------------------------------
! IONS
!---------------------------------------------------------
    write(7777,"(A5)") "&IONS"
! ion_dynamics
    tempRow = BlankString(tempRow)
    tempRow = "ion_dynamics = "//char(34)//trim(adjustl(qeInObj%ion_dynamics))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! End of section
    write(7777,"(A1)") "/"
!---------------------------------------------------------
! CELL
!---------------------------------------------------------
    write(7777,"(A5)") "&CELL"
! cell_dynamics
    tempRow = BlankString(tempRow)
    tempRow = "cell_dynamics = "//char(34)//trim(adjustl(qeInObj%cell_dynamics))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! press
    tempRow = BlankString(tempRow)
    tempRow = "press = "//char(34)//trim(adjustl(DpToStr(qeInObj%press)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! cell_factor
    tempRow = BlankString(tempRow)
    tempRow = "cell_factor = "//char(34)//trim(adjustl(DpToStr(qeInObj%cell_factor)))//char(34)
    write(7777,"(A)") trim(adjustl(tempRow))
! End of section
    write(7777,"(A1)") "/"
!---------------------------------------------------------
! ATOMIC SPECIES
!---------------------------------------------------------
    write(7777,"(A14)") "ATOMIC_SPECIES "
    Do i=1,qeInObj%ntyp
      tempRow = BlankString(tempRow)
      tempRow = trim(adjustl(qeInObj%atomic_species_label(i)))//" "//&
                trim(adjustl(DpToStr(qeInObj%atomic_species_mass(i))))//" "//&
                trim(adjustl(qeInObj%atomic_species_pp(i)))
      write(7777,"(A)") trim(adjustl(tempRow))
    End Do
!---------------------------------------------------------
! ATOMIC POSITIONS
!---------------------------------------------------------
    tempRow = BlankString(tempRow)
    tempRow = "ATOMIC_POSITIONS "//trim(adjustl(qeInObj%atomic_position_type))
    write(7777,"(A)") trim(adjustl(tempRow))
    Do i=1,qeInObj%nat
      tempRow = BlankString(tempRow)
      tempRow = trim(adjustl(qeInObj%atomic_position_label(i)))//" "//&
                trim(adjustl(DpToStr(qeInObj%atomic_position_coords(i,1))))//" "//&
                trim(adjustl(DpToStr(qeInObj%atomic_position_coords(i,2))))//" "//&
                trim(adjustl(DpToStr(qeInObj%atomic_position_coords(i,3))))
      write(7777,"(A)") trim(adjustl(tempRow))
    End Do
!---------------------------------------------------------
! K-POINTS
!---------------------------------------------------------
    tempRow = BlankString(tempRow)
    tempRow = "K_POINTS "//trim(adjustl(qeInObj%kpoints_type))
    write(7777,"(A)") trim(adjustl(tempRow))
    buffer = StrToUpper(qeInObj%kpoints_type)
    If(buffer(1:9).eq."AUTOMATIC")Then
      tempRow = trim(adjustl(IntToStr(qeInObj%kpoints_grid(1))))//" "//&
                trim(adjustl(IntToStr(qeInObj%kpoints_grid(2))))//" "//&
                trim(adjustl(IntToStr(qeInObj%kpoints_grid(3))))//" "//&
                trim(adjustl(IntToStr(qeInObj%kpoints_offset(1))))//" "//&
                trim(adjustl(IntToStr(qeInObj%kpoints_offset(2))))//" "//&
                trim(adjustl(IntToStr(qeInObj%kpoints_offset(3))))
      write(7777,"(A)") trim(adjustl(tempRow))
    End If






    Close(7777)

  End Subroutine writePwscfInput








End Module qe
