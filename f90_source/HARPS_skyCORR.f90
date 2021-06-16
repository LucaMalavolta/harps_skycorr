program HARPS_skyCORR

use Common
use Fits
use Instrument_HARPS
implicit none

! - lettura e2ds di flusso fibra A
!
! - lettura e2ds di flusso fibra B
! - lettura CCF fibra A
! - lettura CCF fibra B
! - lettura coefficienti fibra A
!
! - determinazione rapporto trasmissione fibra (ordine per ordine)
! - correzione per trasmissione
! - correzione per coefficienti
!  -check: si ottiene la CCF finale della fibra A?
! - determinazione della CCF corretta del cielo
! - sottrazione della CCF del cielo alla CCF dell'oggetto

character (len=nch_file):: archive_dir, output_dir, file_rad_obs, file_rad_cal, &
  file_2d_A, file_2d_B, file_ccf_A, file_ccf_B, file_ccf_O,  &
  file_blaze, file_wave, file_mask

integer :: x_e2ds, y_e2ds, x_ccf, y_ccf
integer :: lun_fits, lun_out, io_stat
integer :: ii

real (PR), dimension(:,:), allocatable :: wave_A, wave_B, flux_A, flux_B, ratio_AB, &
 ccf_A, ccf_B, ccf_S
real (4), dimension(:,:), allocatable :: ccf_O

real(PR), dimension(:), allocatable :: factor_B1, factor_B2
real(PR), dimension(6) :: p_coeff

character (len= 20) :: hdutype
integer :: morekeys , blocksize

character (len= 1) :: same_wave_AB

real (PR), dimension(:,:), allocatable :: wave_B_wa, flux_B_wa, ccf_B_wa !! workaround arrays for missing order in CCF_B
integer,parameter :: order_AB_gap = 45

 call getarg(1,archive_dir)
 call getarg(2,file_rad_obs)
 call getarg(3,file_mask)

 if  (iargc() .ge. 4) then
  call getarg(4,output_dir)
 else
   output_dir = '.'
 end if

 if  (iargc() .ge. 5) then
  call getarg(5,same_wave_AB)
 else
   same_wave_AB = 'F'
 end if

 if  (iargc() .lt. 3) then
   write(*,*)
   write(*,*) 'Calling sequence: '
   write(*,*) ' HARPN_skyCORR.e +'
   write(*,*) ' 1) archive directory, including the day of observation'
   write(*,*) ' 2) root of the HARPS/-N fits file '
   write(*,*) ' 3) mask used to extract the CCF'
   write(*,*) ' 4) optional: output directory'
   write(*,*) ' 5) optional: use wave_A for fiber B as well, T o F allowed'
   stop
 end if

!! file_rad is the star observation. we need the rad of the blaze file

 call get_lun(lun_fits)
 call get_lun(lun_out)

 file_ccf_A = trim(archive_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A.fits'
 file_ccf_B = trim(archive_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_B.fits'

 write(*,*) trim(file_ccf_A)
 write(*,*) trim(file_ccf_B)


 call fits_actv(lun_fits,file_ccf_A)
 call singleheader_BLAZE (lun_fits,file_blaze,io_stat)
 write(*,*) trim(file_blaze)
 call singleheader_WAVE  (lun_fits,file_wave ,io_stat)
 write(*,*) trim(file_wave)
 call singleheader_pCOEFF(lun_fits,p_coeff,io_stat)
 write(*,*) p_coeff
 call fits_close(lun_fits)

 !getting wavelength calibration
 file_rad_cal = file_wave(1:29)

 file_2d_A = trim(archive_dir) // '/' // trim(file_rad_cal) // '_wave_A.fits'
 file_2d_B = trim(archive_dir) // '/' // trim(file_rad_cal) // '_wave_B.fits'

 call fits_size_2d(file_2d_A,x_e2ds,y_e2ds,verbose=6)
 allocate(wave_A(x_e2ds,y_e2ds), wave_B(x_e2ds,y_e2ds))
 allocate(flux_A(x_e2ds,y_e2ds), flux_B(x_e2ds,y_e2ds))

 call fits_open_s2d(lun_fits,file_2d_A,wave_A,verbose=6)
 call fits_close(lun_fits)

 if (same_wave_AB=='T') then
   wave_B = wave_A
 else
   allocate(wave_B_wa(x_e2ds,y_e2ds-1))
   call fits_open_s2d(lun_fits,file_2d_B,wave_B_wa,verbose=6)
   call fits_close(lun_fits)
   ! workaround for the missing order in the CCFs
   wave_B = wave_A
   wave_B(:,1:order_AB_gap) = wave_B_wa(:,1:order_AB_gap)
   wave_B(:,order_AB_gap+2:y_e2ds) = wave_B_wa(:,order_AB_gap+1:y_e2ds-1)
   deallocate(wave_B_wa)
 end if


 file_rad_cal = file_blaze(1:29)
 write(*,*) file_rad_cal
 file_2d_A = trim(archive_dir) // '/' // trim(file_rad_cal) // '_lamp_A.fits'
 file_2d_B = trim(archive_dir) // '/' // trim(file_rad_cal) // '_lamp_B.fits'

 write(*,*) trim(file_blaze)
 write(*,*) trim(file_2d_A)
 write(*,*) trim(file_2d_B)



 call fits_open_s2d(lun_fits,file_2d_A,flux_A,verbose=6)
 call fits_close(lun_fits)

 !! workaround for CCF_B missing order
 allocate(flux_B_wa(x_e2ds,y_e2ds-1))
 call fits_open_s2d(lun_fits,file_2d_B,flux_B_wa,verbose=6)
 call fits_close(lun_fits)
 flux_B = flux_A
 flux_B(:,1:order_AB_gap) = flux_B_wa(:,1:order_AB_gap)
 flux_B(:,order_AB_gap+2:y_e2ds) = flux_B_wa(:,order_AB_gap+1:y_e2ds-1)



 allocate(ratio_AB(x_e2ds,y_e2ds), factor_B1(y_e2ds), factor_B2(y_e2ds) )

 write(*,*) x_e2ds, y_e2ds

 ratio_AB = flux_A/flux_B

 factor_B1 = 1._PR
 factor_B2 = 1._PR

  do ii=1,y_e2ds
    !write(*,*) ' A ',flux_A(2048:2052,y_e2ds)
    !write(*,*) ' B ', flux_B(2048:2052,y_e2ds)
    !write(*,*) ' R ', ratio_AB(2048:2052,y_e2ds)
    factor_B1(ii) = sum(ratio_AB(:, ii), dim=1)/real(x_e2ds)

    factor_B2(ii) = p_coeff(1) + p_coeff(2)*wave_B(x_e2ds/2,ii) +  p_coeff(3)*wave_B(x_e2ds/2,ii)**2 + &
    p_coeff(4)*wave_B(x_e2ds/2,ii)**3 + p_coeff(5)*wave_B(x_e2ds/2,ii)**4 + p_coeff(6)*wave_B(x_e2ds/2,ii)**5

    write(*,*) ii, wave_B(x_e2ds/2,ii), factor_B1(ii), factor_B2(ii), factor_B1(ii)*factor_B2(ii)

  end do



 call fits_size_2d(file_ccf_A,x_ccf,y_ccf,verbose=6)

 allocate(ccf_A(x_ccf,y_ccf),ccf_B(x_ccf,y_ccf),ccf_S(x_ccf,y_ccf),ccf_O(x_ccf,y_ccf))


 call fits_open_s2d(lun_fits,file_ccf_A,ccf_A,verbose=6)
 call fits_close(lun_fits)

 allocate(ccf_B_wa(x_ccf,y_ccf-1))
 call fits_open_s2d(lun_fits,file_ccf_B,ccf_B_wa,verbose=6)
 call fits_close(lun_fits)

 ccf_B = 0._PR
 ccf_B(:,1:order_AB_gap) = ccf_B_wa(:,1:order_AB_gap)
 ccf_B(:,order_AB_gap+2:y_ccf) = ccf_B_wa(:,order_AB_gap+1:y_ccf-1)


 do ii=1,y_e2ds
   ccf_S(:,ii) = ccf_B(:,ii)* (factor_B1(ii)/factor_B2(ii))
 end do

 do ii=1,x_ccf
  ccf_S(ii,y_ccf) = sum(ccf_S(ii,1:y_e2ds),dim=1)
 end do

 ccf_O = ccf_A - ccf_S

 call fits_actv(lun_fits,file_ccf_A,verbose=6)

 file_ccf_O =  trim(output_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A_skyR.fits'
 !call fits_write_2d(lun_out,file_ccf_O,ccf_O,verbose=6)

 blocksize = 1
 call FTINIT(lun_out,file_ccf_O,blocksize,io_stat); io_stat= 0
 call FTCPHD(lun_fits, lun_out, io_stat)
 call FTP2DE(lun_out,0,x_ccf,x_ccf,y_ccf,ccf_O,io_stat); io_stat=0

 call fits_close(lun_fits)
 call fits_close(lun_out)



deallocate(wave_A, wave_B, flux_A, flux_B, ratio_AB)
deallocate(factor_B1, factor_B2)
deallocate(ccf_A, ccf_B, ccf_S, ccf_O)
end program
