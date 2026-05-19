module polylog_module
implicit none
contains

recursive function Li2(x) result(res)
  implicit none
  real(8), intent(in) :: x
  real(8) :: res, y, pi2
  real(8) :: term
  integer :: k

  pi2 = 1.6449340668482264365d0  ! pi^2/6

  if (x == 1.d0) then
     res = pi2
     return
  endif

  if (x == 0.d0) then
     res = 0.d0
     return
  endif

  if (abs(x) <= 0.5d0) then
     res = 0.d0
     do k = 1, 10000
        term = x**k / dble(k*k)
        res = res + term
        if (abs(term) < 1.d-15) exit
     end do
     return
  endif

  if (x < 1.d0) then
     y = 1.d0 - x
     res = pi2 - log(x)*log(y) - Li2(y)
     return
  endif

  if (x > 1.d0) then
     y = 1.d0/x
     res = pi2 - 0.5d0*log(x)**2 - Li2(y)
     return
  endif

end function Li2

end module polylog_module