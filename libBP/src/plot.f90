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
  Type :: plotSettings
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
  End Type  
  
  Type :: plotData
    Character(Len=32), Dimension(1:100) :: label = "                "
    Integer(kind=StandardInteger), Dimension(1:100,1:2) :: key = -1
    Real(kind=DoubleReal), Dimension(1:10000,1:2) :: dataArr = 0.0D0
    Character(Len=10), Dimension(1:100) :: marker = "."
    Character(Len=10), Dimension(1:100) :: linestyle = "_"
    Integer(kind=StandardInteger), Dimension(1:100) :: dataSetType=0  ! 0 = input, 1 = fit data set
    Character(Len=128), Dimension(1:200) :: fittingText = ""
    Integer(kind=StandardInteger) :: fittingTextLine = 0
    Logical :: fittingSummaryFile = .false.
    Logical :: dataFile = .false.
  End Type  
  
! Marker options:     point: "."  pixel: ","  circle: "o"  square: "s"  star: "*"
! Linestyle:          none: "None"  _ - -- :
!
!  
  
  
End Module plotTypes

Module plot
! Setup Modules
  Use kinds
  Use strings
  Use plotTypes
  Use fitting
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
!    
! Public Subroutines
  Public :: plotInit
  Public :: plotAdd
  Public :: plotStyle
  Public :: plotMake


!Public :: makePlot
  
  Contains
! ---------------------------------------------------------------------------------------------------

  Subroutine plotInit(settings,dataObj) 
! Reset data objects for a plot  
    Implicit None  ! Force declaration of all variables
! Input    
    Type(plotSettings) :: settings
    Type(plotData) :: dataObj
! Reset
    dataObj%label = "                "
    dataObj%key=-1  
    settings%xMin=1.1D99
    settings%xMax=-1.1D99
    settings%yMin=1.1D99
    settings%yMax=-1.1D99  
  End Subroutine plotInit
  
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
        If(k.eq.1)Then
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
          Call plotFit(dataObj, dataArray, label, rowStart, rowEnd, colX, colY, fitType, 200)
          Call plotStyle(dataObj,",","--",keyFit)
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
  End Subroutine plotStyle
  

  Subroutine plotMake(settings, dataObj) 
! Add data to the data object  
    Implicit None  ! Force declaration of all variables
! Input    
    Type(plotSettings) :: settings
    Type(plotData) :: dataObj
! Private
    Character(len=8) :: fileName
    Integer(kind=StandardInteger) :: i, n, k, a, b
    Integer(kind=StandardInteger) :: iStart, iEnd
    Character(len=128) :: tempDirectory 
    Character(len=128) :: outputDirectory 
    Character(len=64) :: outputName     
    Character(len=14) :: dpStr, dpStrA, dpStrB
    Character(len=12) :: arrayName, arrayNum, arrayNameX, arrayNameY
    Character(len=512) :: xLine, yLine
    Character(len=512) :: pLine
    Integer(kind=StandardInteger) :: termExitStat  
    Real(kind=DoubleReal) :: x, y, xMin, xMax, yMin, yMax
    Real(kind=DoubleReal), Dimension(1:10000,1:200) :: csvArray
    Integer(kind=StandardInteger) :: maxRows, maxCols
    Integer(kind=StandardInteger), Dimension(1:200) :: maxRowsArr
! Init vars    
    tempDirectory = settings%tempDirectory
    outputDirectory = settings%outputDirectory 
    outputName = settings%outputName 
    xMin = 0.0D0
    xMax = 0.0D0
    yMin = 0.0D0
    yMax = 0.0D0
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
    pLine = trim(pLine)//trim(IntToStr(settings%width))//"/"//trim(IntToStr(settings%dpi))//","
    pLine = trim(pLine)//trim(IntToStr(settings%height))//"/"//trim(IntToStr(settings%dpi))//"),"
    pLine = trim(pLine)//"dpi="//trim(IntToStr(settings%dpi))//")"
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
    write(701,"(A)") "plt.title('"//trim(settings%title)//"')"
! Axes Labels
    write(701,"(A)") "plt.xlabel('"//trim(settings%xAxis)//"')"
    write(701,"(A)") "plt.ylabel('"//trim(settings%yAxis)//"')"
! Resize axis
    If(settings%xMin.lt.settings%xMax)Then
      write(dpStrA,"(E14.6)") settings%xMin      
      write(dpStrB,"(E14.6)") settings%xMax      
      write(701,"(A)") "plt.xlim("//dpStrA//","//dpStrB//")"     
    Else
      write(dpStrA,"(E14.6)") (xMin-(0.05D0*(xMax-xMin)))
      write(dpStrB,"(E14.6)") (xMax+(0.05D0*(xMax-xMin)))     
      write(701,"(A)") "plt.xlim("//dpStrA//","//dpStrB//")"  
    End If
    If(settings%yMin.lt.settings%yMax)Then
      write(dpStrA,"(E14.6)") settings%yMin      
      write(dpStrB,"(E14.6)") settings%yMax      
      write(701,"(A)") "plt.ylim("//dpStrA//","//dpStrB//")"   
    Else
      write(dpStrA,"(E14.6)") (yMin-(0.05D0*(yMax-yMin)))
      write(dpStrB,"(E14.6)") (yMax+(0.05D0*(yMax-yMin)))     
      write(701,"(A)") "plt.ylim("//dpStrA//","//dpStrB//")"  
    End If
! Set output file
    write(701,"(A)") "plt.savefig('"//trim(outputDirectory)//"/"//&
    trim(outputName)//"',dpi="//trim(IntToStr(settings%dpi))//")"   
! Close file
    close(701)
! Run python and the file to create the chart
    Call execute_command_line("python "//trim(tempDirectory)//"/"//trim(fileName)//".py",&
    exitstat=termExitStat)
! Clean python file    
    If(settings%cleanPyFile)Then
      Call system("rm -f "//(trim(tempDirectory)//"/"//fileName//".py"))
    End If
! --------------------------
! Summary File
! --------------------------    
    If(dataObj%fittingSummaryFile)Then
! Make temp random name
      fileName = TempFileName() 
      open(unit=702,file=(trim(tempDirectory)//"/summary_"//fileName//".txt"))
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
! CSV File
! --------------------------    
    If(dataObj%dataFile)Then 
! Make temp random name      
      fileName = TempFileName()  
      open(unit=703,file=(trim(tempDirectory)//"/data"//fileName//".csv"))
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
    End If
  End Subroutine plotMake
  

End Module plot