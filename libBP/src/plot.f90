! --------------------------------------------------------------!
! Plot
! plotKinds, plotTypes
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Uses Matplotlib python library to build charts
! Requires matplotlib to be installed
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
    Logical ::               cleanPyFile=.true.
  End Type  
  
  Type :: plotData
    Character(Len=32), Dimension(1:100) :: label = "                "
    Integer(kind=StandardInteger), Dimension(1:100,1:2) :: key = -1
    Real(kind=DoubleReal), Dimension(1:10000,1:2) :: dataArr = 0.0D0
    Character(Len=10), Dimension(1:100) :: marker = "."
    Character(Len=10), Dimension(1:100) :: linestyle = "_"
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
  
  Recursive Subroutine plotAdd(dataObj, dataArray, labelIn, fitListIn, rowStartIn, rowEndIn, colXIn, colYIn) 
! Add data to the data object  
    Implicit None  ! Force declaration of all variables
! Input      
    Type(plotData) :: dataObj
    Real(kind=DoubleReal), Dimension(:,:) :: dataArray
    Character(*), Optional :: labelIn
    Integer(kind=StandardInteger), Optional :: rowStartIn, rowEndIn, colXIn, colYIn
    Character(*), Optional :: fitListIn
! Private variables
    Integer(kind=StandardInteger) :: i, n, k, keyFit
    Integer(kind=StandardInteger) :: rowStart, rowEnd, colX, colY
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
! Store data
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
        keyFit = keyFit + 1
        If(fitList(i:i).eq.",".or.fitList(i:i).eq." ")Then
          Call plotFit(dataObj, dataArray, label, rowStart, rowEnd, colX, colY, fitType, 200)
          Call plotStyle(dataObj,",","--")
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
! Add data set
    labelFit = trim(adjustl(label))//" "//fitType
    Call plotAdd(dataObj, fitDataPoints, labelFit, "")
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
      key = k
    End Do
    If(Present(dataSetIn))Then
      keyIn = dataSetIn
      If(keyIn.lt.0)Then
        key = key + keyIn
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
    Integer(kind=StandardInteger) :: termExitStat  
    Real(kind=DoubleReal) :: xMin, xMax, yMin, yMax
! Init vars    
    tempDirectory = settings%tempDirectory
    outputDirectory = settings%outputDirectory 
    outputName = settings%outputName 
    xMin = 0.0D0
    xMax = 0.0D0
    yMin = 0.0D0
    yMax = 0.0D0
! Make temp random name
    fileName = RandName()    
! Open file
    open(unit=701,file=(trim(tempDirectory)//"/"//fileName//".py"))
! write python headers
    write(701,"(A)") "#!/usr/bin/env python"
    write(701,"(A)") "import numpy as np"
    write(701,"(A)") "import matplotlib"
    write(701,"(A)") "matplotlib.use('Agg')"
    write(701,"(A)") "import matplotlib.pyplot as plt"
! Set figure sizes
    write(701,"(A)") "plt.figure(figsize=(1792/144, 1008/144), dpi=144)"    
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
      write(701,"(A)") "plt.plot("//trim(adjustl(arrayNameX))//","&
          //trim(adjustl(arrayNameY))//", "&
          //"marker='"//trim(adjustl(dataObj%marker(k)))//"',"& 
          //"linestyle='"//trim(adjustl(dataObj%linestyle(k)))//"',"& 
          //"label='"//trim(adjustl(dataObj%label(k)))//&
          "')"
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
    write(701,"(A)") "plt.savefig('"//trim(outputDirectory)//"/"//trim(outputName)//"',dpi=144)"   
! Close file
    close(701)
! Run python and the file to create the chart
    Call execute_command_line("python "//trim(tempDirectory)//"/"//trim(fileName)//".py",&
    exitstat=termExitStat)
! Clean python file    
    If(settings%cleanPyFile)Then
      !Call system("rm -f "//(trim(tempDirectory)//"/"//fileName//".py"))
    End If      
  End Subroutine plotMake
  
  

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

  Function RandName() RESULT (fileName)
! force declaration of all variables
    Implicit None  
    Integer(kind=StandardInteger) :: i, characterNum
    Character(len=8) :: fileName
    Real(kind=DoubleReal) :: randNumber
    Character(len=52), Parameter :: alpha = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
! Prepare string
    fileName = "tmp     "
    Do i=4,8
      Call RANDOM_NUMBER(randNumber)
      characterNum = Ceiling(52.0E0*randNumber+1.0E0)
      If(characterNum.lt.1)Then
        characterNum = 1
      End If
      If(characterNum.gt.52)Then
        characterNum = 52
      End If
      fileName(i:i) = alpha(characterNum:characterNum)
    End Do
  End Function RandName

End Module plot