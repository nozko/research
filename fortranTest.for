      program testSample
	! needs 6 spaces before start and end declaration
c	also comments

	! declaration is allowed only at the top of program  haha
	integer :: a(5) = (/ 5, 10, 15, 20, 25 /)

	! same
	write(*,*) 'Hello World'
	print *, 'Hello World'

	x = 0.4
		print *, 'x =', x		! output : 0.400000006  hahahahaha
		print *, '-x =', -X		! not distinguish large and small X=x
	y = .08			! means 0.08
	z = 1			! automatically float type
		print *, 'y =', y, '  z =', z

	! starts from 1 not 0. ※ a(0) -> error
	! use (), not []  hahaha
	print *, a(2)	! a(1):5, a(2):10, a(3):15 ...
	! output from a(1) to a(3)
	print *, (a(i), i=1, 3)

	! empty valiable -> 0.00000000
	print *, 'b =', b

	! function doesn't need 'call'...?
	print *, 'x + y =', plus(x, y)

	! subroutine needs 'call'
	call minus(x, y)

	stop


!!!----- functions -----!!!
	contains

	! function -> return
	function plus(x, y)
		p = x + y
		return		! return p
	end function

	! subroutine -> no return
	subroutine minus(x, y)
		dif = x - y
		print *, 'x - y =', dif
	end subroutine

      end program testSample
