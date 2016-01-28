
Module printModTypes
! Setup Modules
  Use kinds
  Type :: tableObj 
    Logical :: printHeaderRow = .false.
    Logical :: printHeaderColumn = .false.
    Logical :: colAutoWidth = .false.
    Character(Len=64) :: headerRowColumn = "  "   ! only used if row and column headers are there
    Character(Len=64), Dimension(1:25) :: headerRow = "  "
    Character(Len=64), Dimension(1:100) :: headerColumn = "  "
    Character(Len=64), Dimension(1:1000,1:25) :: tableData = "  "
    Integer(kind=StandardInteger) :: padU = 0, padD = 0, padR = 0, padL = 0
    Integer(kind=StandardInteger) :: rows, columns
    Integer(kind=StandardInteger) :: lastRow = 0
    Integer(kind=StandardInteger), Dimension(1:25) :: colWidth = 16
    Integer(kind=StandardInteger) :: printHeaderColumnWidth = 16
  End Type  

  
End Module printModTypes


Module printMod

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Read user input file

! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------

! Setup Modules
  Use kinds
  Use strings
  Use general
  Use units
  Use printModTypes
  
     
  
! Force declaration of all variables
  Implicit None
  

! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: printBR
  Public :: printTableInit
  Public :: printTableAddHeadersR
  Public :: printTableAddHeadersC
  Public :: printTableAddHeadersRC
  Public :: printTableAddRow
  Public :: printTableMake
! Interfaces  
  Interface printTableAddRow
    Module Procedure printTableAddRow_DP, printTableAddRow_Char
  End Interface printTableAddRow
  
  Contains
  
! --------------------------------------------------------------------------  
!     Print Basic
! --------------------------------------------------------------------------   
  
  Subroutine printBR(widthIn, brCharIn)
! Add header row (to head each column)    
    Implicit None   ! Force declaration of all variables
! In/out    
    Integer(kind=StandardInteger), Optional :: widthIn
    Character(Len=1), Optional :: brCharIn
! Private variables
    Integer(kind=StandardInteger) :: i, width
    Character(Len=512) :: printLine
    Character(Len=1) :: brChar
! Optional argument
    width = 90
    If(Present(widthIn))Then
      width = widthIn
    End If
    brChar = "-"
    If(Present(brCharIn))Then
      brChar = brCharIn
    End If
    Do i=1,width
      printLine(i:i) = brChar
    End Do  
    Print *,printLine(1:width)
  End Subroutine printBR    

    
! --------------------------------------------------------------------------  
!     Print Table
! --------------------------------------------------------------------------  
  
  Subroutine printTableInit(table)
! Add header row (to head each column)    
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
! Init values    
    table%printHeaderRow = .false.
    table%printHeaderColumn = .false.
    table%colAutoWidth = .false.
    table%headerRowColumn = BlankString(table%headerRowColumn)
    table%headerRow = BlankStringArray(table%headerRow)
    table%headerColumn = BlankStringArray(table%headerColumn)
    table%tableData = BlankStringArray(table%tableData)
    table%padU = 0
    table%padD = 0
    table%padR = 0
    table%padL = 0
    table%rows = 0
    table%columns = 0
    table%lastRow = 0
    table%colWidth = 16
    table%printHeaderColumnWidth = 16
  End Subroutine printTableInit
  
  Subroutine printTableAddHeadersR(table,headerIn)
! Add header row (to head each column)    
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Character(*) :: headerIn
    Character(Len(headerIn)) :: header
    Character(Len=64), Dimension(1:25) :: headerRow 
!  Character(Len=64) :: headerRowBlank 
    Integer(kind=StandardInteger) :: columns
    Character(Len=1) :: fieldSplit
! Blank array
    headerRow = BlankStringArray(headerRow)
! Store input
    header = headerIn    
    table%printHeaderRow = .true.
! explode input    
    fieldSplit = ","    
    Call explode(header,fieldSplit,headerRow,columns)    
    table%headerRow = headerRow
    table%columns = columns
  End Subroutine printTableAddHeadersR
  
  Subroutine printTableAddHeadersC(table,headerIn)
! Add header column (to head each row)    
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Character(*) :: headerIn
    Character(Len(headerIn)) :: header
    Character(Len=64), Dimension(1:100) :: headerColumn 
    Integer(kind=StandardInteger) :: rows
    Character(Len=1) :: fieldSplit
! Blank array
    headerColumn = BlankStringArray(headerColumn)
! Store input
    header = headerIn    
    table%printHeaderColumn = .true.
! explode input    
    fieldSplit = ","    
    Call explode(header,fieldSplit,headerColumn,rows)  
    table%headerColumn = headerColumn
    table%rows = rows
  End Subroutine printTableAddHeadersC
  
  Subroutine printTableAddHeadersRC(table,headerIn)
! Add header column (to head each row)    
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Character(*) :: headerIn
    Character(Len=64) :: headerRowColumn 
! Blank array
    headerRowColumn = BlankString(headerRowColumn)
! Store input
    headerRowColumn = trim(headerIn)
    table%headerRowColumn = headerRowColumn
  End Subroutine printTableAddHeadersRC
  
  Subroutine printTableAddRow_DP(table,rowIn)
! Add row to tableData
! Character(Len=64), Dimension(1:1000,1:100) :: tableData
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Real(kind=DoubleReal), Dimension(:) :: rowIn
    Integer(kind=StandardInteger) :: row, column    
! Increment row    
    row = table%lastRow + 1
! Add row
    table%lastRow = row    
! Blank array
    If(row.eq.1)Then
      table%tableData = BlankStringArray(table%tableData)
    End If
! Store data    
    Do column=1,size(rowIn,1)
      table%tableData(row,column) = DpToStr(rowIn(column))    
    End Do
    If(size(rowIn,1).gt.table%columns)Then
      table%columns = size(rowIn,1)
    End If
    If(row.gt.table%rows)Then
      table%rows = row
    End If
  End Subroutine printTableAddRow_DP
  
  Subroutine printTableAddRow_Char(table,rowIn)
! Add row to tableData
! Character(Len=64), Dimension(1:1000,1:100) :: tableData
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Character(*), Dimension(:) :: rowIn
    Integer(kind=StandardInteger) :: row, column    
! Increment row    
    row = table%lastRow + 1
! Add row
    table%lastRow = row    
! Blank array
    If(row.eq.1)Then
      table%tableData = BlankStringArray(table%tableData)
    End If
! Store data    
    Do column=1,size(rowIn,1)
      table%tableData(row,column) = Trim(rowIn(column))    
    End Do
    If(size(rowIn,1).gt.table%columns)Then
      table%columns = size(rowIn,1)
    End If
    If(row.gt.table%rows)Then
      table%rows = row
    End If
  End Subroutine printTableAddRow_Char 
  
  
  Subroutine printTableMake(table)
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table   
    If(table%printHeaderRow)Then
      Call lineRow(table)  
      Call headerRow(table)
    End If  
    Call lineRow(table)
    Call dataRows(table)     
    Call lineRow(table)
  End Subroutine printTableMake
  
! --------------------------------------------------------------------

  Subroutine lineRow(table)
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Character(Len=1024) :: tableRow
    Integer(kind=StandardInteger) :: i, pixel, column
! Blank string    
    tableRow = BlankString(tableRow)    
    pixel = 0
! header column
    If(table%printHeaderColumn)Then   
      pixel = pixel + 1
      tableRow(pixel:pixel) = "+"    
      Do i=1,table%padL     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = "-"   
      End Do       
      Do i=1,table%printHeaderColumnWidth
        pixel = pixel + 1
        tableRow(pixel:pixel) = "-"  
      End Do
      Do i=1,table%padR     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = "-"   
      End Do
    End If    
! columns       
    pixel = pixel + 1
    tableRow(pixel:pixel) = "+"  
    Do column=1,table%columns
      Do i=1,table%padL     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = "-"   
      End Do
      Do i=1,table%colWidth(column)
        pixel = pixel + 1
        tableRow(pixel:pixel) = "-"  
      End Do
      Do i=1,table%padR     ! Right Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = "-"   
      End Do
    End Do  
    pixel = pixel + 1
    tableRow(pixel:pixel) = "+"         
! Print
    print *,trim(tableRow)
  End Subroutine lineRow
  
  Subroutine headerRow(table)
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Character(Len=1024) :: tableRow
    Character(Len=128) :: headerText
    Integer(kind=StandardInteger) :: i, pixel, column, strSize    
! Blank string    
    tableRow = BlankString(tableRow)    
    pixel = 0
! header column
    If(table%printHeaderColumn)Then   
      pixel = pixel + 1
      tableRow(pixel:pixel) = "|"   
! centre text      
      Call TrimString(table%headerRowColumn, strSize, " ")
      headerText = BlankString(headerText)
      headerText = table%headerRowColumn
      Call StrCenter(headerText,table%printHeaderColumnWidth)  
      Do i=1,table%padL     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = " "   
      End Do        
      Do i=1,table%printHeaderColumnWidth
        pixel = pixel + 1
        tableRow(pixel:pixel) = headerText(i:i)
      End Do
      Do i=1,table%padR     ! Right Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = " "   
      End Do
    End If    
! columns       
    pixel = pixel + 1
    tableRow(pixel:pixel) = "|"  
    Do column=1,table%columns
      Call TrimString(table%headerRow(column), strSize, " ")
      headerText = BlankString(headerText)
      headerText = table%headerRow(column)
      Call StrCenter(headerText,table%colWidth(column))  
      Do i=1,table%padL     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = " "   
      End Do  
      Do i=1,table%colWidth(column)
        pixel = pixel + 1
        tableRow(pixel:pixel) = headerText(i:i)  
      End Do
      Do i=1,table%padR     ! Right Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = " "   
      End Do
    End Do  
    pixel = pixel + 1
    tableRow(pixel:pixel) = "|"  
    print *,trim(tableRow)
  End Subroutine headerRow
  
  Subroutine dataRow(table, row)
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Integer(kind=StandardInteger) :: row
    Character(Len=1024) :: tableRow
    Character(Len=128) :: dataText
    Integer(kind=StandardInteger) :: i, pixel, column, strSize
! Blank string    
    tableRow = BlankString(tableRow)    
    pixel = 0
! header column
    If(table%printHeaderColumn)Then   
      pixel = pixel + 1
      tableRow(pixel:pixel) = "|"      
! center text
      Call TrimString(table%headerColumn(row), strSize, " ")
      dataText = BlankString(dataText)
      dataText = table%headerColumn(row)
      Call StrCenter(dataText,table%printHeaderColumnWidth)   
! write to line  
      Do i=1,table%padL     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = " "   
      End Do      
      Do i=1,table%printHeaderColumnWidth
        pixel = pixel + 1
        tableRow(pixel:pixel) = dataText(i:i)
      End Do
      Do i=1,table%padR     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = " "   
      End Do  
    End If    
! columns       
    pixel = pixel + 1
    tableRow(pixel:pixel) = "|"  
    Do column=1,table%columns
! center text
      Call TrimString(table%tableData(row,column), strSize, " ")
      dataText = BlankString(dataText)
      dataText = table%tableData(row,column)
      Call StrCenter(dataText,table%colWidth(column))   
! write to line      
      Do i=1,table%padL     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = " "   
      End Do  
      Do i=1,table%colWidth(column)
        pixel = pixel + 1
        tableRow(pixel:pixel) = dataText(i:i)  
      End Do
      Do i=1,table%padR     ! Left Padding
        pixel = pixel + 1
        tableRow(pixel:pixel) = " "   
      End Do  
    End Do  
    pixel = pixel + 1
    tableRow(pixel:pixel) = "|"  
    print *,trim(tableRow)
  End Subroutine dataRow
  
  Subroutine dataRows(table)
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(tableObj) :: table
    Integer(kind=StandardInteger) :: row
! Loop through rows
    Do row=1,table%rows
      Call dataRow(table,row)
    End Do
  End Subroutine dataRows
  
  
  
  
End Module printMod