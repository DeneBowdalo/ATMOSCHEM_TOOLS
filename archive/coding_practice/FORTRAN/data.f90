PROGRAM data
	IMPLICIT NONE

	! This program illustrates the main types of data object which 
	! Fortran90 recognises.	

	! Variable declarations
	INTEGER :: min=4,max=25,number,can = 234
	REAL :: x=2.0,y, con=5.667
	COMPLEX :: z
	LOGICAL :: logical_1,logical_2
	CHARACTER(LEN=10) :: string_1,string_2='Physics'

	! evaluation of variable values
	min = min-2
	max=max+2*min
	number=(max-min)**2 + 7
	y=x**min+ 0.00000345 -1.23E-02 + 45.678E-5
	z=CMPLX(x,y) + (1.0E-2,2.0E+12)
	logical_1=.TRUE.
	logical_2=.FALSE.
	string_1='Hi there'

	! output
	PRINT *, 'min =',min,'  max =',max,'   number = ',number
	PRINT *
	PRINT *, 'x =',x,'  y =',y
	PRINT *
	PRINT *, 'z =',z
	PRINT *
	PRINT *, 'logical_1 =',logical_1,'  logical_2 =',logical_2
	PRINT *
	PRINT *, 'string_1 = ',string_1,'   string_2 = ',string_2
	PRINT *, can, con
END PROGRAM data

!! 1.	A Fortran program is a list of instructions called `statements'.
!!	These command the CPU to perform various actions, one at a time
!!	in the sequence as specified.

!! 2.	Fortran statements perform operations on `data objects' such as 
!!	integers, real numbers, complex numbers, logical values or strings 
!!	of characters. If values are stated explicitly they are called 
!!	`constants'. Data objects which are allowed to change their values
!!	are called `variables'. Note, however, that a Fortran variable is a
!!	named memory cell, i.e. a pigeonhole with a fixed NAME but 
!!	in which a variable value may be stored.

!! 3.	Fortran recognises 6 intrinsic data types.  The five which we shall 
!!	use in these sessions are
!!
!!		INTEGER		:for storing signed integers
!!		REAL		:for storing real numbers
!!		COMPLEX		:for storing a complex number (a pair of 
!!				 ordered REALs)
!!		LOGICAL		:for storing the logical values .TRUE. and
!!				 .FALSE.
!!		CHARACTER	:for storing one or more characters
!!
!!	The declaration of variables is illustrated in the program. Note that
!!	the CHARACTER(LEN=10) indicates the longest string of characters 
!!	that can be stored is 10. You may insert what you like here. Also 	
!!	CHARACTER variables must be inserted between apostrophes or quotes.
!!
!!	[Note that f90 allows programmers to create their own data types.]
!!	[Higher precision may be obtained by using the sixth intrinsic data 
!!	type: DOUBLE PRECISION. This is a relic of f77: f90 has alternative  
!!	features which allow improved precision obviating the need to 
!!	use double precision - see textbook for these.]

!! 4.	Examples are:
!!
!!	INTEGER		0, -1, 34512, +34512, -0
!!	REAL		0.0, 0.0E0, 1.234E-08, +123.000E0, 324.896E5 
!!	COMPLEX		(1.2, 4.67), (1.2, 467.0E-02)
!!	LOGICAL		two values only:  .TRUE. or .FALSE.
!!	CHARACTER	'What do you think of it so far', "bye", '28-3-49'

!! 5.	Declaration of variables
!!
!!	The inclusion of the IMPLICIT NONE statement means that ALL variables 
!!	must be declared. The declarations must appear at the beginning of 
!!	the program before any executable statements. Examine the format in 
!!	the above program. The syntax is
!!
!!		TYPE ::	name1, name2, .... , namen
!!
!!	where TYPE may be
!!
!!		INTEGER
!!		REAL
!!		COMPLEX
!!		CHARACTER(LEN=len)
!!		LOGICAL
!!
!!	name1, name2, ..., namen  are variable names (storage cell identifiers).
!!	'len' is the maximum number of characters that can be held in the 
!!	character variable declared. If a string has less than 'len' characters
!!	then unused spaces are filled by blanks. If 'LEN=len' is absent, then 
!!	the default is 'LEN=1'. The values of character variables must be 
!!	enclosed within apostrophes or quotes!

!! 6.	Compile and execute the above program. Examine the output and ensure 
!!	that it matches up to what you expect.
!!
!!	If you wish your executable file to be called other than 'a.out', then
!!	typing 
!!
!!		f90  prog.f90  -o  prog.out
!!
!!	will produce an executable file called 'prog.out'.

!! 7.	PAY ATTENTION TO DETAIL!!!!
!!
!!	A missing comma, apostrophe, etc will ensure that the program will not 
!!	compile. Wrong syntax will be diagnosed by the compiler; usually 
!!	a location (line number) is given and some information about what is 
!!	wrong. For example: edit out a comma from one of the print statements.
!!	Now compile your program and observe what happens. Correct the program.

!! 8.	Initial values

!!	WARNING: Declaring a variable reserves a storage cell in which a 
!!	value will be put. It does NOT set an initial value. Variables should 
!!	therefore be initialised, via READ or assignment statements, prior to 
!!	their first use in the program. Indeed, not doing so is a common 
!!	source of error. 
!!
!!	Initial values are most easily set by following the name in the 
!!	declaration by an equals sign and the initial value (see min, max, x,
!!	string_2 in the program above).

!! 9.	Exercises 2
!!
!!	a) State the type of data of each of the following constants  
!!	and if they are illegal explain why:
!!
!!		-872869		.FALSE.		5.0E5		byebye
!!		'Hello  '	(6.0,7.0)	TRUE		2+3i
!!		32,000		"is this OK"	{0.0,-5.4}	'hcb"
!!
!!
!!	b) Edit the program data.f90 above to include and print out a variable 
!!	of each type of your own choosing.

!!	End of file: data2.f90

