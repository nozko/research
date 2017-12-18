      program testSample
	! needs 6 spaces before start and end declaration
c	also comments

	! 変数定義は最初しかできないポンコツ仕様 hahaha
	integer :: a(5) = (/ 5, 10, 15, 20, 25 /)

	! both available
	write(*,*) 'Hello World'
	print *, 'Hello World!'

	x = 0.4
	print *, x		! なぜか出力は0.400000006
	print *, -x
	y = .08			! means 0.08
	z = 1			! automatically float type
	print *, y, z

	! starts from 1 not 0. ※ a(0) -> error
	print *, a(2)	! a(1):5, a(2):10, a(3):15 ...

	! output from a(1) to a(3)
	print *, (a(i), i=1, 3)

	! empty valiable -> 0
	print *, B

	stop
      end program testSample
