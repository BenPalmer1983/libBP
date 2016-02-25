! --------------------------------------------------------------!
! Plot
! plotKinds, plotTypes
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Uses Matplotlib python library to build charts
! Requires matplotlib to be installed
!
! Two modules must be used - plotTypes to load the types, and plot to use the plot functions
! and subroutines.  It must be compiled as a part of this library of functions as there are
! quite a lot of modules it depends on.
!
! There is an option to perform one or multiple fits on the input data.
!
! ----------------------------------------
! Updated: 21st September 2015
! ----------------------------------------

Module plotTypes
! Setup Modules
  Use kinds

  Type :: plotData
    Character(len=128) ::    tempDirectory
    Character(len=128) ::    outputDirectory
    Character(len=64) ::     outputName
    Character(len=32) ::     title
    Character(len=32) ::     xAxis
    Character(len=32) ::     yAxis
    Real(kind=DoubleReal) :: xMin=1.1D99
    Real(kind=DoubleReal) :: xMax=-1.1D99
    Real(kind=DoubleReal) :: yMin=1.1D99
    Real(kind=DoubleReal) :: yMax=-1.1D99
    Integer(kind=StandardInteger) :: dpi=144
    Integer(kind=StandardInteger) :: width=1792
    Integer(kind=StandardInteger) :: height=1008
    Logical ::               cleanPyFile=.true.
    Logical ::               cleanGpFile=.true.
    Logical ::               cleanGpLatexFile=.true.
    Character(Len=32), Dimension(1:100) :: label = "                "                          ! Data set labels
    Integer(kind=StandardInteger), Dimension(1:100,1:2) :: key = -1
    Real(kind=DoubleReal), Dimension(1:10000,1:2) :: dataArr = 0.0D0
    Character(Len=10), Dimension(1:100) :: marker = "."
    Character(Len=10), Dimension(1:100) :: linestyle = "_"
    Character(Len=32), Dimension(1:100) :: gpLinestyle = "points                          "    ! points, lines, lines lt 1 dt 3
    Integer(kind=StandardInteger) :: gpLineColour=0
    Integer(kind=StandardInteger), Dimension(1:100) :: dataSetType=0
    Character(Len=128), Dimension(1:200) :: fittingText = ""  ! POLY2,POLY3,POLY4,POLY5,EXPFIT1,EXPFIT2,EXPFIT3,BM1,BM2,SPLINE,SPLINE1,SPLINE3,SPLINE5,INTERP,INTERP3,INTERP4,INTERP5
    Integer(kind=StandardInteger) :: fittingTextLine = 0
    Logical :: fittingSummaryFile = .false.
    Logical :: dataFile = .false.
    Logical :: gnuPlotFile = .true.
    Logical :: gnuPlotLatexFile = .false.
    Logical :: pyPlotFile = .false.
    Character(Len=512) :: fitListInput = ""
    Character(Len=512) :: csvFileInput = ""
    Character(Len=128) :: gplDir = ""                         ! directory tex and eps files are in relative to LaTeX main file e.g. chapter1/plots/
    Integer(kind=StandardInteger) :: numFitPoints = 500
    Character(Len=3) :: gplOutType = "PNG"
  End Type

! Marker options:     none      No point
!                     point: "."  pixel: ","  circle: "o"  square: "s"  star: "*"
!
!
! Linestyle:          none      No line
!                     -         Solid
!                     --        Dashed
!                     -.        Dash-Dot
!                     :         Dotted
!

End Module plotTypes

Module plot
! Setup Modules
  Use kinds
  Use strings
  Use general
  Use plotTypes
  Use splinesFitting
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
!
! Public Subroutines
  Public :: plotInit
  Public :: plotReadInput
  Public :: plotLoadData
  Public :: plotAdd
  Public :: plotStyle
  Public :: plotMake


!Public :: makePlot

  Contains
! ---------------------------------------------------------------------------------------------------

  Subroutine plotInit(dataObj)
! Reset data objects for a plot
    Implicit None  ! Force declaration of all variables
! Input
    Type(plotData) :: dataObj
! Reset
    dataObj%label = "                "
    dataObj%key=-1
    dataObj%xMin=1.1D99
    dataObj%xMax=-1.1D99
    dataObj%yMin=1.1D99
    dataObj%yMax=-1.1D99
  End Subroutine plotInit

  Subroutine plotReadInput(dataObj, inputFile)
! Reset data objects for a plot
    Implicit None  ! Force declaration of all variables
! In/Out:  Declare variables
    Type(plotData) :: dataObj
    Character(*) :: inputFile
! Private: Declare variables
    Character(Len=256), Dimension(1:100) :: fileArray
    Character(Len=16) ::  bufferA
    Character(Len=128) ::  bufferB
    Integer(kind=StandardInteger) :: rowCount
    Integer(kind=StandardInteger) :: i
! read file into array
    Call readFile(inputFile, fileArray, rowCount)
    Do i=1,rowCount
      Read(fileArray(i),*) bufferA, bufferB
      bufferA = StrToUpper(Trim(Adjustl(bufferA)))
      If(bufferA(1:13).eq."TEMPDIRECTORY")Then
        dataObj%tempDirectory = Trim(Adjustl(bufferB))
      End If
      If(bufferA(1:15).eq."OUTPUTDIRECTORY")Then
        dataObj%outputDirectory = Trim(Adjustl(bufferB))
      End If
      If(bufferA(1:10).eq."OUTPUTNAME")Then
        dataObj%outputName = Trim(Adjustl(bufferB))
      End If
      If(bufferA(1:5).eq."TITLE")Then
        dataObj%title = Trim(Adjustl(bufferB))
      End If
      If(bufferA(1:5).eq."XAXIS")Then
        dataObj%xAxis = Trim(Adjustl(bufferB))
      End If
      If(bufferA(1:5).eq."YAXIS")Then
        dataObj%yAxis = Trim(Adjustl(bufferB))
      End If
      If(bufferA(1:5).eq."WIDTH")Then
        dataObj%width = StrToInt(Trim(Adjustl(bufferB)))
      End If
      If(bufferA(1:6).eq."HEIGHT")Then
        dataObj%height = StrToInt(Trim(Adjustl(bufferB)))
      End If
      If(bufferA(1:8).eq."DATAFILE")Then
        dataObj%dataFile = StrToBool(Trim(Adjustl(bufferB)))
      End If
      If(bufferA(1:11).eq."GNUPLOTFILE")Then
        dataObj%gnuPlotFile = StrToBool(Trim(Adjustl(bufferB)))
      End If
      If(bufferA(1:10).eq."PYPLOTFILE")Then
        dataObj%pyPlotFile = StrToBool(Trim(Adjustl(bufferB)))
      End If
      If(bufferA(1:16).eq."GNUPLOTLATEXFILE")Then
        dataObj%gnuPlotLatexFile = StrToBool(Trim(Adjustl(bufferB)))
      End If
      If(bufferA(1:7).eq."FITLIST")Then
        dataObj%fitListInput = Trim(Adjustl(bufferB))
      End If
      If(bufferA(1:7).eq."CSVFILE")Then
        dataObj%csvFileInput = Trim(Adjustl(bufferB))
      End If
    End Do
  End Subroutine plotReadInput


  Subroutine plotLoadData(dataObj, filePath, fitList)
! Loads data from a csv file and fits each data set
! First row is headers
    Implicit None  ! Force declaration of all variables
! In/Out:  Declare variables
    Type(plotData) :: dataObj
    Character(*) :: filePath, fitList
! Private: Declare variables
    Character(Len=256), Dimension(1:10000) :: fileArray
    Integer(kind=StandardInteger) :: rowCount
    Real(kind=DoubleReal), Dimension(1:10000,1:50) :: csvArray
    Real(kind=DoubleReal), Dimension(1:10000,1:2) :: dataArray
    Character(Len=256), Dimension(1:50) :: dataHeadings
    Character(Len=256), Dimension(1:50) :: colArray
    Integer(kind=StandardInteger) :: cols, dataSets
    Integer(kind=StandardInteger), Dimension(1:50) :: maxRowsArr
    Integer(kind=StandardInteger) :: i, j
! read file into array
    Call readFile(filePath, fileArray, rowCount)
! data headings
    Call explode(fileArray(1), ",", dataHeadings, cols)
    print *,fitList
! read in data
    maxRowsArr = rowCount-1
    Do i=2,rowCount
      Call explode(fileArray(i), ",", colArray, cols)
      Do j=1,cols
        colArray(j) = trim(adjustl(colArray(j)))
        If(colArray(j)(1:1).eq." ")Then
          If(maxRowsArr(j).eq.(rowCount-1))Then
            maxRowsArr(j) = i-2
          End If
        Else
          csvArray(i-1,j) = StrToDp(colArray(j))
        End If
      End Do
    End Do
! Add data
    dataSets = ceiling(cols/2.0D0)
    Do i=1,dataSets
      Do j=1,maxRowsArr(2*(i-1)+1)
        dataArray(j,1) = csvArray(j,2*(i-1)+1)
        dataArray(j,2) = csvArray(j,2*(i-1)+2)
      End Do
      Call plotAdd(dataObj, dataArray, dataHeadings(2*(i-1)+2), fitList, 1, maxRowsArr(2*(i-1)+1) , 1, 2)
    End Do
  End Subroutine plotLoadData


  Recursive Subroutine plotAdd(dataObj, dataArray, labelIn, fitListIn, rowStartIn, rowEndIn, colXIn, colYIn, fitAddIn)
! Add data to the data object
    Implicit None  ! Force declaration of all variables
! Input
    Type(plotData) :: dataObj
    Real(kind=DoubleReal), Dimension(:,:) :: dataArray
    Character(*), Optional :: labelIn
    Integer(kind=StandardInteger), Optional :: rowStartIn, rowEndIn, colXIn, colYIn
    Character(*), Optional :: fitListIn
    Integer(kind=StandardInteger), Optional :: fitAddIn
! Private variables
    Integer(kind=StandardInteger) :: i, n, k, keyFit
    Integer(kind=StandardInteger) :: rowStart, rowEnd, colX, colY, fitAdd
    Character(Len=32) :: label
    Character(Len=64) :: fitList
    Character(Len=12) :: fitType
! Set optional arguments
    label = BlankString(label)
    rowStart = 1
    rowEnd = size(dataArray,1)
    colX = 1
    colY = 2
    fitList = BlankString(fitList)
    fitAdd = 0
    If(Present(labelIn))Then
      label = labelIn
    End If
    If(Present(rowStartIn))Then
      rowStart = rowStartIn
    End If
    If(Present(rowEndIn))Then
      rowEnd = rowEndIn
    End If
    If(Present(colXIn))Then
      colX = colXIn
    End If
    If(Present(colYIn))Then
      colY = colYIn
    End If
    If(Present(fitListIn))Then
      fitList = fitListIn
    End If
    fitList = Trim(Adjustl(StrToUpper(fitList)))
    If(Present(fitAddIn))Then
      fitAdd = fitAddIn
    End If
! Store data - key = k
    Do k=1,size(dataObj%key,1)
      If(dataObj%key(k,1).eq.-1)Then
        If(k.eq.1)Then ! key = 1
! init

! set n
          n = 0
        Else
          n = dataObj%key(k-1,2)
        End If
        dataObj%label(k) = label
        dataObj%key(k,1) = (n + 1)
        Do i=rowStart,rowEnd
          n = n + 1
          dataObj%dataArr(n,1) = dataArray(i,1)
          dataObj%dataArr(n,2) = dataArray(i,2)
        End Do
        dataObj%key(k,2) = n
        dataObj%dataSetType(k) = fitAdd
        Exit
      End If
    End Do
! Fit function to data
    If(fitList(1:1).ne." ")Then
      fitType = BlankString(fitType)
      n = 0
      keyFit = k
      Do i=1,64
        n = n + 1
        If(fitList(i:i).eq.",".or.fitList(i:i).eq." ")Then
          keyFit = keyFit + 1
          Call plotFit(dataObj, dataArray, label, rowStart, rowEnd, colX, colY, fitType, dataObj%numFitPoints)
          Call plotStyle(dataObj,"None","--",keyFit)
          fitType = BlankString(fitType)
          n = 0
          If(fitList(i:i).eq." ")Then
            Exit
          End If
        Else
          fitType(n:n) = fitList(i:i)
        End If
      End Do
    End If
  End Subroutine plotAdd

  Subroutine plotFit(dataObj, dataArray, label, rowStart, rowEnd, colX, colY, fitType, dataPoints)
! Add data to the data object
    Implicit None  ! Force declaration of all variables
! Input
    Type(plotData) :: dataObj
    Real(kind=DoubleReal), Dimension(:,:) :: dataArray
    Character(Len=32) :: label
    Integer(kind=StandardInteger) :: rowStart, rowEnd, colX, colY, dataPoints
    Character(Len=12) :: fitType
! Private
    Real(kind=DoubleReal), Dimension(1:(rowEnd-rowStart+1),1:2) :: inputDataPoints
    Real(kind=DoubleReal), Dimension(1:dataPoints,1:2) :: fitDataPoints
    Integer(kind=StandardInteger) :: i, n
    Character(Len=32) :: labelFit
    !Real(kind=DoubleReal), Dimension(1:(fitPoly+1)) :: coefficients
    !Real(kind=DoubleReal) :: x, xStart, xEnd, xInc, y
! Transfer data
    i = 0
    Do n=rowStart,rowEnd
      i = i + 1
      inputDataPoints(i,1) = dataArray(n,colX)
      inputDataPoints(i,2) = dataArray(n,colY)
    End Do
! Fit data points
    fitDataPoints = FittingPoints(inputDataPoints, fitType, dataPoints)
! Store data
    Do i=1,20
      If(Trim(fittingReport(i)).eq."")Then
        Exit
      Else
        n = dataObj%fittingTextLine
        n = n + 1
        dataObj%fittingText(n) = fittingReport(i)
        dataObj%fittingTextLine = n
      End If
    End Do
! Add data set
    labelFit = trim(adjustl(label))//" "//fitType
    Call plotAdd(dataObj, fitDataPoints, labelFit, "",1,dataPoints,1,2,1)
  End Subroutine plotFit

  Subroutine plotStyle(dataObj, marker, linestyle, dataSetIn)
! Apply style to last dataset added (unless otherwise specified)
    Implicit None  ! Force declaration of all variables
! In:      Declare variables
    Type(plotData) :: dataObj
    Character(*) :: marker
    Character(*) :: linestyle
    Character(Len=10) :: markerUC, linestyleUC
    Character(Len=1) :: lineColour
    Integer(kind=StandardInteger), Optional :: dataSetIn
! Private: Declare variables
    Integer(kind=StandardInteger) :: k, key, keyIn
! Optional
    key = 1
    Do k=1,size(dataObj%key,1)
      If(dataObj%key(k,1).eq.-1)Then
        Exit
      End If
      If(dataObj%dataSetType(k).eq.0)Then
        key = k
      End If
    End Do
    If(Present(dataSetIn))Then
      keyIn = dataSetIn
      If(keyIn.lt.0)Then
        key = key + keyIn
      Else
        key = keyIn
      End If
    End If
! Set marker and linestyle
    dataObj%marker(key) = marker
    dataObj%linestyle(key) = linestyle
    markerUC = StrToUpper(marker)
    linestyleUC = StrToUpper(linestyle)

    dataObj%gpLinestyle(key) = "linespoints"
    If(markerUC.eq."NONE".and.linestyleUC.ne."NONE")Then
      dataObj%gpLinestyle(key) = "lines"
    End If
    If(markerUC.ne."NONE".and.linestyleUC.eq."NONE")Then
      dataObj%gpLinestyle(key) = "points"
    End If

    If(linestyleUC.ne."NONE")Then
      dataObj%gpLineColour = mod(dataObj%gpLineColour,9)+1
      Write(lineColour,"(I1)") dataObj%gpLineColour
      If(linestyleUC(1:2).eq."--")Then
        dataObj%gpLinestyle(key) = trim(dataObj%gpLinestyle(key))//" lt 2 lc "//lineColour
      Elseif(linestyleUC(1:2).eq."-.")Then
        dataObj%gpLinestyle(key) = trim(dataObj%gpLinestyle(key))//" lt 5 lc "//lineColour
      Elseif(linestyleUC(1:2).eq.": ")Then
        dataObj%gpLinestyle(key) = trim(dataObj%gpLinestyle(key))//" lt 4 lc "//lineColour
      Else
        dataObj%gpLinestyle(key) = trim(dataObj%gpLinestyle(key))//" lt 1 lc "//lineColour
      End If
    End If


  End Subroutine plotStyle


  Subroutine plotMake(dataObj)
! Add data to the data object
    Implicit None  ! Force declaration of all variables
! Input
    Type(plotData) :: dataObj
! Private
    Character(len=8) :: fileName
    Integer(kind=StandardInteger) :: i, n, k, a, b
    Integer(kind=StandardInteger) :: iStart, iEnd
    Character(len=128) :: tempDirectory
    Character(len=128) :: outputDirectory
    Character(len=64) :: outputName
    Character(len=256) :: csvFile
    Character(len=14) :: dpStr, dpStrA, dpStrB, intStrA, intStrB
    Character(len=12) :: arrayName, arrayNum, arrayNameX, arrayNameY
    Character(len=512) :: xLine, yLine
    Character(len=512) :: pLine
    Character(len=128) :: cmdLine, gpDataLine
    Character(len=64) :: gpLabel
    Character(len=64) :: outputPy, outputGP
    Integer(kind=StandardInteger) :: termExitStat
    Real(kind=DoubleReal) :: x, y, xMin, xMax, yMin, yMax
    Real(kind=DoubleReal), Dimension(1:10000,1:200) :: csvArray
    Integer(kind=StandardInteger) :: maxRows, maxCols
    Integer(kind=StandardInteger), Dimension(1:200) :: maxRowsArr
! Init vars
    tempDirectory = dataObj%tempDirectory
    outputDirectory = dataObj%outputDirectory
    outputName = dataObj%outputName
    xMin = 0.0D0
    xMax = 0.0D0
    yMin = 0.0D0
    yMax = 0.0D0
! Make directories
    Call makeDir(tempDirectory)
    Call makeDir(outputDirectory)
! File names
    If(dataObj%pyPlotFile.and.dataObj%gnuPlotFile)Then
      outputPy = Trim(Adjustl(outputName))//"_py"
      outputGP =  Trim(Adjustl(outputName))//"_gp"
    Else
      outputPy = Trim(Adjustl(outputName))
      outputGP  = Trim(Adjustl(outputName))
    End If
! --------------------------
! Python File
! --------------------------
! Make temp random name
    fileName = TempFileName()
! Open file
    open(unit=701,file=(trim(tempDirectory)//"/"//fileName//".py"))
! write python headers
    write(701,"(A)") "#!/usr/bin/env python"
    write(701,"(A)") "import numpy as np"
    write(701,"(A)") "import matplotlib"
    write(701,"(A)") "matplotlib.use('Agg')"
    write(701,"(A)") "import matplotlib.pyplot as plt"
! Set figure sizes
    pLine = BlankString(pLine)
    pLine = "plt.figure(figsize=("
    pLine = trim(pLine)//trim(IntToStr(dataObj%width))//"/"//trim(IntToStr(dataObj%dpi))//","
    pLine = trim(pLine)//trim(IntToStr(dataObj%height))//"/"//trim(IntToStr(dataObj%dpi))//"),"
    pLine = trim(pLine)//"dpi="//trim(IntToStr(dataObj%dpi))//")"
    write(701,"(A)") trim(pLine)
    !write(701,"(A)") "plt.figure(figsize=(1792/144, 1008/144), dpi=144)"
! write data arrays
    Do k=1,size(dataObj%key,1)
      If(dataObj%key(k,1).lt.0)Then
        Exit
      End If
! start-end
      iStart = dataObj%key(k,1)
      iEnd = dataObj%key(k,2)
! write x data points to file
      write(arrayNum,"(I8)") k
      arrayName = "arrayX_"//trim(adjustl(arrayNum))
      write(701,"(A)") trim(adjustl(arrayName))//" = ["
      n = 0
      xLine = BlankString(xLine)
      Do i=iStart,iEnd
        n = n + 1
        a = (n-1)*15+1
        b = (n-1)*15+14
        write(dpStr,"(E14.6)") dataObj%dataArr(i,1)
        xLine(a:b) = dpStr
        If(i.lt.iEnd)Then
          xLine(b+1:b+1) = ","
        End If
        If(n.eq.5.or.i.eq.iEnd)Then
          n = 0
          write(701,"(A)") "    "//trim(xLine)
          xLine = BlankString(xLine)
        End If
! find min/may x/y while looping
        If(i.eq.iStart.and.k.eq.1)Then
          xMin = dataObj%dataArr(i,1)
          xMax = dataObj%dataArr(i,1)
          yMin = dataObj%dataArr(i,2)
          yMax = dataObj%dataArr(i,2)
        Else
          If(dataObj%dataArr(i,1).lt.xMin)Then
            xMin = dataObj%dataArr(i,1)
          End If
          If(dataObj%dataArr(i,1).gt.xMax)Then
            xMax = dataObj%dataArr(i,1)
          End If
          If(dataObj%dataArr(i,2).lt.yMin)Then
            yMin = dataObj%dataArr(i,2)
          End If
          If(dataObj%dataArr(i,2).gt.yMax)Then
            yMax = dataObj%dataArr(i,2)
          End If
        End If
      End Do
      write(701,"(A)") "]"
! write x data points to file
      arrayName = "arrayY_"//trim(adjustl(arrayNum))
      write(701,"(A)") trim(adjustl(arrayName))//" = ["
      n = 0
      yLine = BlankString(yLine)
      Do i=iStart,iEnd
        n = n + 1
        a = (n-1)*15+1
        b = (n-1)*15+14
        write(dpStr,"(E14.6)") dataObj%dataArr(i,2)
        yLine(a:b) = dpStr
        If(i.lt.iEnd)Then
          yLine(b+1:b+1) = ","
        End If
        If(n.eq.5.or.i.eq.iEnd)Then
          n = 0
          write(701,"(A)") "    "//trim(yLine)
          yLine = BlankString(yLine)
        End If
      End Do
      write(701,"(A)") "]"
    End Do
! write plot lines
    Do k=1,size(dataObj%key,1)
      If(dataObj%key(k,1).lt.0)Then
        Exit
      End If
! prepare
      dataObj%label(k) = CleanString(dataObj%label(k))
! write
      write(arrayNum,"(I8)") k
      arrayNameX = "arrayX_"//trim(adjustl(arrayNum))
      arrayNameY = "arrayY_"//trim(adjustl(arrayNum))
      pLine = BlankString(pLine)
      pLine = "plt.plot("//trim(adjustl(arrayNameX))//","
      pLine = trim(pLine)//trim(adjustl(arrayNameY))//","
      pLine = trim(pLine)//"marker='"//trim(adjustl(dataObj%marker(k)))//"',"
      pLine = trim(pLine)//"linestyle='"//trim(adjustl(dataObj%linestyle(k)))//"',"
      pLine = trim(pLine)//"label='"//trim(adjustl(dataObj%label(k)))//"'"
      pLine = trim(pLine)//")"
      write(701,"(A)") trim(pLine)
    End Do
! Write Titles
    write(701,"(A)") "plt.title('"//trim(dataObj%title)//"')"
! Axes Labels
    write(701,"(A)") "plt.xlabel('"//trim(dataObj%xAxis)//"')"
    write(701,"(A)") "plt.ylabel('"//trim(dataObj%yAxis)//"')"
! Resize axis
    If(dataObj%xMin.lt.dataObj%xMax)Then
      write(dpStrA,"(E14.6)") dataObj%xMin
      write(dpStrB,"(E14.6)") dataObj%xMax
      write(701,"(A)") "plt.xlim("//dpStrA//","//dpStrB//")"
    Else
      write(dpStrA,"(E14.6)") (xMin-(0.05D0*(xMax-xMin)))
      write(dpStrB,"(E14.6)") (xMax+(0.05D0*(xMax-xMin)))
      write(701,"(A)") "plt.xlim("//dpStrA//","//dpStrB//")"
    End If
    If(dataObj%yMin.lt.dataObj%yMax)Then
      write(dpStrA,"(E14.6)") dataObj%yMin
      write(dpStrB,"(E14.6)") dataObj%yMax
      write(701,"(A)") "plt.ylim("//dpStrA//","//dpStrB//")"
    Else
      write(dpStrA,"(E14.6)") (yMin-(0.05D0*(yMax-yMin)))
      write(dpStrB,"(E14.6)") (yMax+(0.05D0*(yMax-yMin)))
      write(701,"(A)") "plt.ylim("//dpStrA//","//dpStrB//")"
    End If
! Set output file
    write(701,"(A)") "plt.savefig('"//trim(outputDirectory)//"/"//&
    trim(outputPy)//"',dpi="//trim(IntToStr(dataObj%dpi))//")"
! Close file
    close(701)
! Run python and the file to create the chart
    If(dataObj%pyPlotFile)Then
      Call execute_command_line("python "//trim(tempDirectory)//"/"//trim(fileName)//".py",&
      exitstat=termExitStat)
    End If
! Clean python file
    If(dataObj%cleanPyFile)Then
      Call system("rm -f "//(trim(tempDirectory)//"/"//fileName//".py"))
    End If
! --------------------------
! Summary File
! --------------------------
    If(dataObj%fittingSummaryFile)Then
! Make temp random name
      fileName = TempFileName()
      open(unit=702,file=(trim(outputDirectory)//"/summary_"//fileName//".txt"))
      Do i=1,200
        If(trim(dataObj%fittingText(i)).eq."")Then
          Exit
        Else
          write(702,"(A)") trim(dataObj%fittingText(i))
        End If
      End Do
! Close file
      Close(702)
    End If
! --------------------------
! CSV File and GnuPlot
! --------------------------
    If(dataObj%dataFile.or.dataObj%gnuPlotFile)Then
      csvFile = trim(outputDirectory)//"/"//trim(dataObj%outputName)//".csv"
      open(unit=703,file=trim(csvFile))
! Blank array
      csvArray = 0.0D0
      maxRowsArr = 0
      a = 0
      maxRows = 0
      maxCols = 0
! Loop through data sets
      Do k=1,size(dataObj%key,1)
        If(dataObj%key(k,1).lt.0)Then
          maxCols = k-1
          Exit
        End If
! start-end
        iStart = dataObj%key(k,1)
        iEnd = dataObj%key(k,2)
        n = 0
        a = a + 1
        xLine = BlankString(xLine)
! Loop through data points
        Do i=iStart,iEnd
          n = n + 1
          csvArray(n,2*(a-1)+1) = dataObj%dataArr(i,1)
          csvArray(n,2*(a-1)+2) = dataObj%dataArr(i,2)
        End Do
        maxRowsArr(a) = n
        If(maxRows.lt.n)Then
          maxRows = n
        End If
      End Do
! write to file
      Do n=1,maxRows
        pLine = BlankString(pLine)
        If(n.eq.1)Then
          Do k=1,maxCols
            gpLabel = BlankString(gpLabel)
            gpLabel = Trim(Adjustl(dataObj%label(k)))
            If(gpLabel(1:1).eq." ")Then
              write(intStrA,"(I8)") k
              gpLabel = "Set "//Trim(Adjustl(intStrA))
            End If
            If(k.eq.maxCols)Then
              pLine = Trim(pLine)//" -, "//Trim(gpLabel)
            Else
              pLine = Trim(pLine)//" -, "//Trim(gpLabel)//", "
            End If
          End Do
          write(703,"(A)") trim(pLine)
          pLine = BlankString(pLine)
        End If
        Do i=1,maxCols
          If(n.le.maxRowsArr(i))Then
            x = csvArray(n,2*(i-1)+1)
            y = csvArray(n,2*(i-1)+2)
            pLine = Trim(pLine)//Trim(DpToStr(x))//&
            ","//Trim(DpToStr(y))
            If(i.lt.maxCols)Then
              pLine = Trim(pLine)//","
            End If
          Else
            pLine = Trim(pLine)//","
            If(i.lt.maxCols)Then
              pLine = Trim(pLine)//","
            End If
          End If
        End Do
        write(703,"(A)") trim(pLine)
      End Do
! Close file
      Close(703)
! --------------------------
! GNU Plot File (only if CSV file has been made)
! --------------------------
      If(dataObj%gnuPlotFile)Then
! Make temp random name
        fileName = TempFileName()
! GnuPlot file
! ------------------
        open(unit=704,file=(trim(tempDirectory)//"/gp_"//fileName//".gplot"))
        dataObj%gplOutType = StrToUpper(dataObj%gplOutType)
        If(dataObj%gplOutType.eq."PNG")Then
          write(704,"(A)") "set terminal pngcairo size "&
          //trim(IntToStr(dataObj%width))//","&
          //trim(IntToStr(dataObj%height))//" enhanced font 'Verdana,10'"
          write(704,"(A)") "set output "//char(34)//trim(outputDirectory)//"/"&
          //trim(outputGp)//".png"//char(34)
        End If
        If(dataObj%gplOutType.eq."EPS")Then
          write(704,"(A)") "set terminal postscript eps monochrome enhanced blacktext size "&
          //trim(IntToStr(dataObj%width))//","&
          //trim(IntToStr(dataObj%height))//""
          write(704,"(A)") "set output "//char(34)//trim(outputDirectory)//"/"&
          //trim(outputGp)//".eps"//char(34)
        End If
        write(704,"(A)") "set grid xtics mxtics ytics mytics back"
        write(704,"(A)") "set datafile separator "//char(34)//","//char(34)
        write(704,"(A)") "set title "//char(34)//trim(dataObj%title)//char(34)
        write(704,"(A)") "set xlabel "//char(34)//trim(dataObj%xAxis)//char(34)
        write(704,"(A)") "set ylabel "//char(34)//trim(dataObj%yAxis)//char(34)
        If(dataObj%yMin.lt.dataObj%yMax)Then
          write(dpStrA,"(E14.6)") dataObj%yMin
          write(dpStrB,"(E14.6)") dataObj%yMax
          write(704,"(A)") "set yrange ["//dpStrA//":"//dpStrB//"]"
        End If
        write(704,"(A)") "set key autotitle columnheader"
        write(704,"(A)") "plot \"
! Loop through data sets
        Do k=1,size(dataObj%key,1)
          If(dataObj%key(k,1).lt.0)Then
            maxCols = k-1
            Exit
          End If
        End Do
        Do k=1,maxCols
          write(intStrA,"(I8)") (2*(k-1)+1)
          write(intStrB,"(I8)") (2*(k-1)+2)
          gpDataLine = BlankString(gpDataLine)
          If(k.gt.1)Then  ! start with a comma
            gpDataLine = ","
          End If
          gpDataLine = trim(gpDataLine)//"'"//trim(csvFile)//"' using "&
          //trim(adjustl(intStrA))//":"//trim(adjustl(intStrB))//" with "&
          //trim(dataObj%gpLinestyle(k))//" "
          If(k.lt.maxCols)Then  ! end line with new line \ unless last data set
           gpDataLine = trim(gpDataLine)//" \ "
          End If
          write(704,"(A)") trim(gpDataLine)
        End Do
        Close(704)
! Run gnuplot
        cmdLine = "gnuplot "//trim(tempDirectory)//"/gp_"//fileName//".gplot"
        Call execute_command_line(trim(cmdLine),&
        exitstat=termExitStat)
! Clean gp file
        If(dataObj%cleanGpFile)Then
          Call system("rm -f "//trim(tempDirectory)//"/gp_"//fileName//".gplot")
        End If
      End If
! GnuPlot file for latex
! ------------------
      If(dataObj%gnuPlotLatexFile)Then
        open(unit=705,file=(trim(tempDirectory)//"/gpl_"//fileName//".gplot"))
        write(705,"(A)") "cd "//char(34)//trim(outputDirectory)//char(34)
        write(705,"(A)") "set terminal epslatex monochrome size 17cm,9.56cm "
        write(705,"(A)") "set output "//char(34)//trim(outputGp)//".tex"//char(34)
        write(705,"(A)") "set grid xtics mxtics ytics mytics back"
        write(705,"(A)") "set datafile separator "//char(34)//","//char(34)
        write(705,"(A)") "set title "//char(34)//trim(dataObj%title)//char(34)
        write(705,"(A)") "set xlabel "//char(34)//trim(dataObj%xAxis)//char(34)
        write(705,"(A)") "set ylabel "//char(34)//trim(dataObj%yAxis)//char(34)
        write(705,"(A)") "set key autotitle columnheader"
        write(705,"(A)") "plot \"
! Loop through data sets
        Do k=1,size(dataObj%key,1)
          If(dataObj%key(k,1).lt.0)Then
            maxCols = k-1
            Exit
          End If
        End Do
        Do k=1,maxCols
          write(intStrA,"(I8)") (2*(k-1)+1)
          write(intStrB,"(I8)") (2*(k-1)+2)
          gpDataLine = BlankString(gpDataLine)
          If(k.gt.1)Then  ! start with a comma
            gpDataLine = ","
          End If
          gpDataLine = trim(gpDataLine)//"'"//trim(csvFile)//"' using "&
          //trim(adjustl(intStrA))//":"//trim(adjustl(intStrB))//" with "&
          //trim(dataObj%gpLinestyle(k))//" "
          If(k.lt.maxCols)Then  ! end line with new line \ unless last data set
           gpDataLine = trim(gpDataLine)//" \ "
          End If
          write(705,"(A)") trim(gpDataLine)
        End Do
        Close(705)
! Run gnuplot
        cmdLine = "gnuplot "//trim(tempDirectory)//"/gpl_"//fileName//".gplot"
        Call execute_command_line(trim(cmdLine),&
        exitstat=termExitStat)
! Clean gpl file
        If(dataObj%cleanGpLatexFile)Then
          Call system("rm -f "//trim(tempDirectory)//"/gpl_"//fileName//".gplot")
        End If
      End If
! Clean csv
      If(dataObj%dataFile)Then
        ! Keep csv file
      Else
        Call system("rm -f "//trim(csvFile))
      End If
    End If



  End Subroutine plotMake


End Module plot
