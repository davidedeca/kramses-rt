subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_d converts mass density from user units into g/cc
  scale_d = units_density
  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3

  ! scale_t converts time from user units into seconds
  scale_t = units_time
  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)

  ! scale_l converts distance from user units into cm
  scale_l = units_length
  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)

  ! scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

end subroutine units

subroutine krome_units(scale_d,krome_scale,inv_krome_scale)

  use amr_parameters, ONLY: dp
  use krome_user, ONLY: krome_nmols,KROME_idx_H,KROME_idx_E,KROME_idx_Hj,&
    &                   KROME_idx_HE,KROME_idx_HEj,KROME_idx_HEjj,KROME_idx_Hk,&
    &                   KROME_idx_H2,KROME_idx_H2j

  implicit none
  real(dp)::scale_d
  real(dp),dimension(1:krome_nmols)::krome_scale,inv_krome_scale

  krome_scale(KROME_idx_H   ) = scale_d/1.67353251819d-24 !H
  krome_scale(KROME_idx_E   ) = scale_d/9.10938188d-28 !E
  krome_scale(KROME_idx_Hj  ) = scale_d/1.67262158d-24 !H+
  krome_scale(KROME_idx_HE  ) = scale_d/6.69206503638d-24 !HE
  krome_scale(KROME_idx_HEj ) = scale_d/6.69115409819d-24 !HE+
  krome_scale(KROME_idx_HEjj) = scale_d/6.69024316d-24 !HE++
  krome_scale(KROME_idx_Hk  ) = scale_d/1.67444345638d-24 !H-
  krome_scale(KROME_idx_H2  ) = scale_d/3.34706503638d-24 !H2
  krome_scale(KROME_idx_H2j ) = scale_d/3.34615409819d-24 !H2+

  inv_krome_scale(1:krome_nmols) = 1.d0/krome_scale(1:krome_nmols)

end subroutine krome_units
