subroutine calcrh(temp,dewpt,rhvalue)

  !abstract: calculate relative humidity using temperature,pressure,and dew point
  !use w3lib routine w3fa09 to calculate saturation vapor pressure at temp and dew point

  use constants, only: r100,r1000,tiny_r_kind,zero
  use kinds
  implicit none

  real(4),intent(in)::temp,dewpt
  real(4),intent(out)::rhvalue
  real(4)::vp,satvp,w3fa09

  if (temp.lt.zero.or.dewpt.lt.zero) then
     rhvalue=tiny_r_kind !missing value is tiny_r_kind to ensure it is not max
  else
     vp=w3fa09(dewpt)
     satvp=w3fa09(temp)
     rhvalue=(vp/satvp)*r100
     if (rhvalue.gt.r100) then
        rhvalue=r100
     endif
  endif

end subroutine calcrh
