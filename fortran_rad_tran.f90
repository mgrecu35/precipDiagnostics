subroutine calculate_tb_prof(iprof, incAngle, freqs, temperature_2m, u_10m, v_10m, qv, geopotential, &
      temperature, press, tb_all_freq) 
      ! Inputs
      integer, intent(in) :: iprof
      real, intent(in) :: incAngle
      real, intent(in), dimension(:) :: freqs
      real, intent(in), dimension(:) :: temperature_2m, u_10m, v_10m, press
      real, intent(in), dimension(:,:) :: qv, geopotential, temperature
  
      ! Outputs
      real, intent(out), dimension(size(freqs)) :: tb_all_freq
  
      ! Local variables
      integer :: nlyr, ifreq, k, nk, ireturn
      real :: fisot, umu, emis, ebar, w, btemp
      real, dimension(:), allocatable :: hL, lyrtemp, lyrhgt
      real, dimension(29) :: kext1D, asym1D, salb1D
      real :: rho, absair, abswv, tavg, ts
      integer :: nz(2)
      ! Constants
      fisot = 2.7
      umu = cos(incAngle / 180.0 * acos(-1.0))
      nlyr = 29
  
      ! Initialize output
      tb_all_freq = 0.0
  
      ! Loop over frequencies
      allocate(hL(0:nlyr))
      allocate(lyrtemp(0:nlyr))
      allocate(lyrhgt(0:nlyr))
      nz=shape(temperature)
      do ifreq = 1, size(freqs)
        ! Initialize profile-specific variables
        
  
        hL(0) = 0.0
        lyrtemp(0) = temperature_2m(iprof)
  
        do k = 1, nlyr
          nk = nz(2)+1-k
          ireturn = 0
          rho = press(nk) * 100.0 / (287.05 * temperature(iprof,nk))
          call gasabsr98(freqs(ifreq), temperature(iprof, nk), qv(iprof, nk) * rho, press(nk)*100, absair, abswv, ireturn)
          hL(k) = (geopotential(iprof, nk) + geopotential(iprof, nk-1)) / 2.0 / 9.8
          tavg = (temperature(iprof, nk) + temperature(iprof, nk-1)) / 2.0
          lyrtemp(k) = tavg
  
          kext1D(k) = absair + abswv
          asym1D(k) = 0.0
          salb1D(k) = 0.0
          !print*, kext1D(k), asym1D(k), salb1D(k), hL(k), lyrtemp(k), freqs(ifreq)
        end do
  
        ts = temperature_2m(iprof)
        w = sqrt(u_10m(iprof)**2 + v_10m(iprof)**2)
        call emit(freqs(ifreq), 1, ts, w, umu, emis, ebar)
  
        ! Convert heights to km
        
        lyrhgt = hL / 1.0e3
  
        btemp = temperature_2m(iprof)
        call radtran(umu, nlyr, tb_all_freq(ifreq), btemp, lyrtemp, lyrhgt, kext1D, salb1D, asym1D, fisot, emis, ebar, 0 )
  
        ! Deallocate local arrays
        
      end do
      deallocate(hL, lyrtemp, lyrhgt)
end subroutine calculate_tb_prof
  
  