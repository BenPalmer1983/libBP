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
    Character(Len=16), Dimension(1:100) :: label = "                "
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
  Use plotTypes
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
  
  Subroutine plotAdd(dataObj, dataArray, labelIn, rowStartIn, rowEndIn, colXIn, colYIn, fitPolyIn) 
! Add data to the data object  
    Implicit None  ! Force declaration of all variables
! Input      
    Type(plotData) :: dataObj
    Real(kind=DoubleReal), Dimension(:,:) :: dataArray
    Character(*), Optional :: labelIn
    Integer(kind=StandardInteger), Optional :: rowStartIn, rowEndIn, colXIn, colYIn, fitPolyIn
! Private variables
    Integer(kind=StandardInteger) :: i, n, k
    Integer(kind=StandardInteger) :: rowStart, rowEnd, colX, colY, fitPoly
    Character(Len=16) :: label
! Set optional arguments    
    label = "               "
    rowStart = 1
    rowEnd = size(dataArray,1)
    colX = 1
    colY = 2
    fitPoly = 0
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
    If(Present(fitPolyIn))Then
      fitPoly = fitPolyIn
    End If
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
! Fit polynomial    
    If(fitPoly.gt.1)Then    
      Call plotFit(dataObj, dataArray, rowStart, rowEnd, colX, colY, fitPoly, (k+1))
    End If
  End Subroutine plotAdd
    
  Subroutine plotFit(dataObj, dataArray, rowStart, rowEnd, colX, colY, fitPoly, key)
! Add data to the data object  
    Implicit None  ! Force declaration of all variables
! Input    
    Type(plotData) :: dataObj
    Real(kind=DoubleReal), Dimension(:,:) :: dataArray
    Integer(kind=StandardInteger) :: rowStart, rowEnd, colX, colY, fitPoly
! Private
    Real(kind=DoubleReal), Dimension(1:(rowEnd-rowStart+1),1:2) :: fitData
    Integer(kind=StandardInteger) :: i, n, k, kStart, kEnd
    Integer(kind=StandardInteger) :: key, dataPoints
    Real(kind=DoubleReal), Dimension(1:(fitPoly+1)) :: coefficients
    Real(kind=DoubleReal) :: x, xStart, xEnd, xInc, y
! Init    
    dataPoints = 100
    kStart = (dataObj%key(key-1,2) + 1)
    kEnd = kStart+dataPoints-1
! Transfer data    
    i = 0
    Do n=rowStart,rowEnd
      i = i + 1
      fitData(i,1) = dataArray(n,colX)
      fitData(i,2) = dataArray(n,colY)
    End Do  
! Fit polynomial
    coefficients = PolyFit(dataArray,fitPoly)
! save data
    dataObj%key(key,1) = kStart
    dataObj%key(key,2) = kEnd
    dataObj%label(key) = "Fit: "//trim(dataObj%label(key-1))
! make points
    xStart = dataArray(rowStart,colX)
    xEnd = dataArray(rowEnd,colX)
    xInc = (xEnd-xStart)/(1.0D0*(dataPoints-1))     
    x = xStart
    k = kStart
    Do i=1,dataPoints   
      y = CalcPolynomial(coefficients,x)
      dataObj%dataArr(k,1) = x
      dataObj%dataArr(k,2) = y
! Increment      
      x = x + xInc
      k = k + 1
    End Do  
  End Subroutine plotFit
  
  Subroutine plotStyle(dataObj, dataSet, marker, linestyle)
! Add data to the data object  
    Implicit None  ! Force declaration of all variables
! Input    
    Type(plotData) :: dataObj
    Integer(kind=StandardInteger) :: dataSet
    Character(*) :: marker
    Character(*) :: linestyle
! Set marker and linestyle
    dataObj%marker(dataSet) = marker
    dataObj%linestyle(dataSet) = linestyle    
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
          //"linestyle='"//trim(adjustl(dataObj%linestyle(k)))//"')"
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
! Plot maths functions                                                   
! ------------------------------------------------------------------------!
  
  Function PolyFit(points,order) RESULT (coefficients)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: k,col,row,exponentValue
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Integer(kind=StandardInteger) :: order
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficients
    Real(kind=DoubleReal), Dimension(1:(order+1),1:(order+1)) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: yMatrix
! Step 1 - Standard fit with Vandermonde matrix
! Build Least Squares Fitting Vandermonde matrix
    Do row=1,(order+1)
      Do col=1,(order+1)
        exponentValue = row+col-2
        xMatrix(row,col) = 0.0D0
        Do k=1,size(points,1)
          xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*points(k,1)&
          **exponentValue
        End Do
      End Do
    End Do
    Do row=1,(order+1)
      exponentValue = row-1
      yMatrix(row) = 0.0D0
      Do k=1,size(points,1)
        yMatrix(row) = 1.0D0*yMatrix(row)+1.0D0*points(k,2)*&
        points(k,1)**exponentValue
      End Do
    End Do
! invert xMatrix
    xMatrix = InvertMatrix(xMatrix)
! multiply inverse by y to get coefficients
    coefficients = matMul(xMatrix,yMatrix)
  End Function PolyFit
  
   Function InvertMatrix(xMatrix) RESULT (xMatrixInverse)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Integer(kind=StandardInteger) :: row,col,rowb
    Integer(kind=StandardInteger) :: matrixSize
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:2*size(xMatrix,1)) :: xMatrixWorking
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,1)) :: xMatrixInverse
    Real(kind=DoubleReal), Dimension(1:2*size(xMatrix,1)) :: xMatrixRow
! matrix(row,column)
! Initialise variables
    row = 0
    rowb = 0
    col = 0
    matrixSize = size(xMatrix,1)
    xMatrixWorking = 0.0D0
    xMatrixInverse = 0.0D0
    xMatrixRow = 0.0D0
! if a square matrix
    If(size(xMatrix,1).eq.size(xMatrix,2))Then
! Fill working array
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
        End Do
      End Do
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.eq.col)Then
            xMatrixWorking(row,col+matrixSize) = 1.0D0
          End If
        End Do
      End Do
! make lower triangle of zeros
      Do row=1,matrixSize-1
        Do rowb=row+1,matrixSize
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the lower triangle
        Do rowb=row+1,matrixSize
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! re-force zeros in the lower triangle
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.gt.col)Then
            xMatrixWorking(row,col) = 0.0D0
          End If
        End Do
      End Do
! make upper triangle of zeros
      Do row=matrixSize,2,-1
        Do rowb=row-1,1,-1
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the upper triangle
        Do rowb=row-1,1,-1
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! Divide rhs by diagonal on lhs and store in inverse
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixInverse(row,col) = 1.0D0*&
          xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
        End Do
      End Do
    End If
  End Function InvertMatrix
  
  Function CalcPolynomial(polyCoefficientsIn, x, derivIn) RESULT (y)
! Calculates p(x) by default, p'(x) for derivativeIn = 1 etc
! p(x) = polyCoefficientsIn(1) x^0 +  polyCoefficientsIn(2) x^1  + ... + polyCoefficientsIn(n) x^(n-1)
    Implicit None  !Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: polyCoefficientsIn
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger), Optional :: derivIn
! Out
    Real(kind=DoubleReal) :: y
! Private Variables
    Real(kind=DoubleReal), Dimension(1:size(polyCoefficientsIn,1)) :: polyCoefficients
    Integer(kind=StandardInteger) :: deriv
    Integer(kind=StandardInteger) :: i, j
! Optional
    deriv = 0
    If(Present(derivIn))Then
      deriv = derivIn
    End If
! Transfer coeffs to internal/private array
    polyCoefficients = polyCoefficientsIn
! Modify as necessary
    Do i=1,deriv
      Do j=1,size(polyCoefficients,1)
        polyCoefficients(j) = (j-i)*polyCoefficients(j)
      End Do
    End Do
! calculate output
    y = 0.0D0
    Do j=1,size(polyCoefficients,1)
      y = y + polyCoefficients(j)*x**(j-1-deriv)
    End Do
  End Function CalcPolynomial
  

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!
  
  Function BlankString(input) RESULT (output)
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i
    Do i=1,Len(input)
      output(i:i) = " "
    End Do
  End Function BlankString
  
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
  
  Function RemoveSpaces(input)RESULT(output)
! force declaration of all variables
    Implicit None    
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! Local variables
      Integer(kind=StandardInteger) :: i, j
! Copy input string
      outputTemp = input
! Blank output
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
! transfer outputtemp to output without spaces
      j = 0
      Do i = 1, LEN(outputTemp)
        If(outputTemp( i:i ).ne." ")Then
          j = j + 1
          output( j:j ) = outputTemp( i:i )
        End If
      End Do
  End Function RemoveSpaces
End Module plot