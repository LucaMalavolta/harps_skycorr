program HARPN_skyCORR

use Common
use Fits
use Instrument_HARPN
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

character (len=nch_file):: archive_dir, calib_dir, file_list, &
  output_dir, file_rad_obs, file_rad_cal, &
  file_2d_A, file_2d_B, file_ccf_A, file_ccf_B, file_ccf_O,  &
  file_blaze, file_wave, file_mask, file_date

integer :: x_e2ds, y_e2ds, x_ccf, y_ccf
integer :: lun_fits, lun_out, io_stat, lun_input
integer :: ii

real (PR), dimension(:,:), allocatable :: wave_A, wave_B, flux_A, flux_B, ratio_AB, &
 ccf_A, ccf_B, ccf_S
real (4), dimension(:,:), allocatable :: ccf_O

real(PR), dimension(:), allocatable :: factor_B1, factor_B2
real(PR), dimension(6) :: p_coeff

character (len= 20) :: hdutype
integer :: morekeys , blocksize

integer :: date_calib, date_ccf, sts, sts_out, date_output

character (len= 1) :: same_wave_AB
integer :: read_iostat
character (len = nch_file) :: arch_red, opt_in

 date_calib = 1
 date_ccf= 1

 if  (iargc() .lt. 6) then
   write(*,*)
   write(*,*) 'Calling sequence: '
   write(*,*) ' HARPN_skyCORR.e +'
   write(*,*) ' 1) input file list  '
   write(*,*) ' 2) directory with CCF files '
   write(*,*) ' 3) use the date to search for CCF files'
   write(*,*) ' 4) calibration directory'
   write(*,*) ' 5) use the date to search for calibration files'
   write(*,*) ' 6) output directory'
   write(*,*) ' 7) use the date for the output directory'
   write(*,*) ' 8) optional: use wave_A for fiber B as well, T o F allowed'
   stop
 end if

 call getarg(1,file_list)
 call getarg(2,archive_dir)
 call getarg(3,opt_in); read (opt_in,*) date_ccf
 call getarg(4,calib_dir)
 call getarg(5,opt_in); read (opt_in,*) date_calib
 call getarg(6,output_dir)
 call getarg(7,opt_in); read (opt_in,*) date_output

 if  (iargc() .ge. 8) then
  call getarg(8,same_wave_AB)
 else
   same_wave_AB = 'F'
 end if



!! file_rad is the star observation. we need the rad of the blaze file

 call get_lun(lun_fits)
 call get_lun(lun_out)
 call get_lun(lun_input)

 read_iostat = 0
 open(lun_input,file=trim(file_list), status='old')

 do while (read_iostat.eq.0)
   read(lun_input,*,iostat=read_iostat) file_date, file_rad_obs, file_mask
   if (read_iostat.ne.0) exit

   if (date_ccf.eq.1) then
     file_ccf_A = trim(archive_dir) // '/' // trim(file_date)// '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A.fits'
     file_ccf_B = trim(archive_dir) // '/' // trim(file_date)// '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_B.fits'
   else
     file_ccf_A = trim(archive_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A.fits'
     file_ccf_B = trim(archive_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_B.fits'
  end if

  write(*,*)
  write(*,*) trim(file_date), trim(file_rad_obs), trim(file_mask)
   write(*,*) trim(file_ccf_A)
   write(*,*) trim(file_ccf_B)

   call fits_check(file_ccf_B,sts_out)
   write(*,*) sts_out
   if (sts_out.eq.0) then
     if (date_output.eq.1) then
         call system('mkdir -p ' // trim(output_dir) // '/' // trim(file_date) )
         call system('cp ' // trim(file_ccf_A) // ' ' //  trim(output_dir) // '/' // trim(file_date) // '/')
     else
         call system('mkdir -p ' // trim(output_dir))
         call system('cp ' // trim(file_ccf_A) // ' ' //  trim(output_dir) // '/' )
     end if
     write(*,*) 'File copied: ',  trim(file_ccf_A)
     cycle
   end if

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
   if (date_calib.eq.1) then
     file_2d_A = trim(calib_dir) // '/' // trim(file_date)// '/' // trim(file_rad_cal) // '_wave_A.fits'
     file_2d_B = trim(calib_dir) // '/' // trim(file_date)// '/' // trim(file_rad_cal) // '_wave_B.fits'
   else
     file_2d_A = trim(calib_dir) // '/' // trim(file_rad_cal) // '_wave_A.fits'
     file_2d_B = trim(calib_dir) // '/' // trim(file_rad_cal) // '_wave_B.fits'
   end if

   call fits_size_2d(file_2d_A,x_e2ds,y_e2ds,verbose=6)
   allocate(wave_A(x_e2ds,y_e2ds), wave_B(x_e2ds,y_e2ds))
   allocate(flux_A(x_e2ds,y_e2ds), flux_B(x_e2ds,y_e2ds))

   call fits_open_s2d(lun_fits,file_2d_A,wave_A,verbose=6)
   call fits_close(lun_fits)

   if (same_wave_AB=='T') then
     wave_B = wave_A
   else
     call fits_open_s2d(lun_fits,file_2d_B,wave_B,verbose=6)
     call fits_close(lun_fits)
   end if


   file_rad_cal = file_blaze(1:29)

   if (date_calib.eq.1) then
     file_2d_A = trim(calib_dir) // '/' // trim(file_date)// '/' // trim(file_rad_cal) // '_lamp_A.fits'
     file_2d_B = trim(calib_dir) // '/' // trim(file_date)// '/' // trim(file_rad_cal) // '_lamp_A.fits'
   else
     file_2d_A = trim(calib_dir) // '/' // trim(file_rad_cal) // '_lamp_A.fits'
     file_2d_B = trim(calib_dir) // '/' // trim(file_rad_cal) // '_lamp_A.fits'
   end if

   write(*,*) trim(file_blaze)
   write(*,*) trim(file_2d_A)
   write(*,*) trim(file_2d_B)



   call fits_open_s2d(lun_fits,file_2d_A,flux_A,verbose=6)
   call fits_close(lun_fits)

   call fits_open_s2d(lun_fits,file_2d_B,flux_B,verbose=6)
   call fits_close(lun_fits)

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

     write(*,*) ii, factor_B1(ii), factor_B2(ii), factor_B1(ii)*factor_B2(ii)

   end do

  call fits_size_2d(file_ccf_A,x_ccf,y_ccf,verbose=6)

  allocate(ccf_A(x_ccf,y_ccf),ccf_B(x_ccf,y_ccf),ccf_S(x_ccf,y_ccf),ccf_O(x_ccf,y_ccf))


  call fits_open_s2d(lun_fits,file_ccf_A,ccf_A,verbose=6)
  call fits_close(lun_fits)

  call fits_open_s2d(lun_fits,file_ccf_B,ccf_B,verbose=6)
  call fits_close(lun_fits)

  do ii=1,y_e2ds
    ccf_S(:,ii) = ccf_B(:,ii)* (factor_B1(ii)/factor_B2(ii))
  end do

  do ii=1,x_ccf
   ccf_S(ii,y_ccf) = sum(ccf_S(ii,1:y_e2ds),dim=1)
  end do

  ccf_O = ccf_A - ccf_S

  call fits_actv(lun_fits,file_ccf_A,verbose=6)
  if (date_output.eq.1) then
      call system('mkdir -p ' // trim(output_dir) // '/' // trim(file_date) )
      call system('cp ' // trim(file_ccf_A) // ' ' //  trim(output_dir) // '/' // trim(file_date) // '/')
      call system('cp ' // trim(file_ccf_B) // ' ' //  trim(output_dir) // '/' // trim(file_date) // '/')
      file_ccf_O =  trim(output_dir) // '/' // trim(file_date) // '/' // trim(file_rad_obs) // &
                '_ccf_'//trim(file_mask)//'_A_skyR.fits'
  else
      call system('mkdir -p ' // trim(output_dir))
      call system('cp ' // trim(file_ccf_A) // ' ' //  trim(output_dir) // '/' )
      call system('cp ' // trim(file_ccf_B) // ' ' //  trim(output_dir) // '/' )
      file_ccf_O =  trim(output_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A_skyR.fits'
  end if
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

end do

close(lun_input)
end program
