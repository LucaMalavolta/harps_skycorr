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

character (len=nch_file):: archive_dir, output_dir, file_rad_obs, file_rad_cal, &
  file_2d_A, file_2d_B, file_ccf_A, file_ccf_B, file_ccf_O,  &
  file_blaze, file_wave, file_mask

integer :: x_e2ds, y_e2ds, x_ccf, y_ccf
integer :: lun_fits, lun_out, lun_input, io_stat
integer :: ii

real (PR), dimension(:,:), allocatable :: wave_A, wave_B, flux_A, flux_B, ratio_AB, &
 ccf_A, ccf_B, ccf_S
real (4), dimension(:,:), allocatable :: ccf_O

real(PR), dimension(:), allocatable :: factor_B1, factor_B2
real(PR), dimension(6) :: p_coeff

character (len= 20) :: hdutype
integer :: morekeys , blocksize

character (len= 1) :: same_wave_AB




!!! New variables added for chromatic RVs
integer :: chromatic_n_ranges, nw, ic
real(PR) :: wave_start, wave_end
real (PR), dimension(:,:), allocatable :: CCF_A_chromatic, CCF_O_chromatic, chromatic_ranges
character (len=nch_file) :: chromatic_list
real :: crval, crpix, cdelt
integer ::  sts
character (len=64) :: comment
integer, dimension(:,:), allocatable :: chromatic_orders
character (len=nch_file) :: order_list_file
character (len=1) :: xi

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

!!! New optional input for chromatic RVs
if  (iargc() .ge. 6) then
    call getarg(6,chromatic_list)
else
    chromatic_list = 'none'
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
    write(*,*) ' 6) optional: chromatic ranges'
    stop
end if

! file_rad is the star observation. we need the rad of the blaze file

call get_lun(lun_fits)
call get_lun(lun_out)
call get_lun(lun_input)

if (chromatic_list == 'none') then
    chromatic_n_ranges = 3
    allocate(chromatic_ranges(chromatic_n_ranges,2))

    chromatic_ranges(1,1) = 3771.5
    chromatic_ranges(1,2) = 4468.9

    chromatic_ranges(2,1) = 4432.6
    chromatic_ranges(2,2) = 5138.7

    chromatic_ranges(3,1) = 5103.9
    chromatic_ranges(3,2) = 7906.4

else
    io_stat = 0
    open(lun_input,file=trim(chromatic_list), status='old')
    read(lun_input,*) chromatic_n_ranges
    allocate(chromatic_ranges(chromatic_n_ranges,2))

    do while (io_stat.eq.0)
        read(lun_input,*,iostat=io_stat) nw, wave_start, wave_end
        if (io_stat.ne.0) exit
        chromatic_ranges(nw, 1) = wave_start
        chromatic_ranges(nw, 2) = wave_end
    end do

    close(lun_input)

end if


file_ccf_A = trim(archive_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A.fits'
file_ccf_B = trim(archive_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_B.fits'

write(*,*) ' CCF file fiber A:  ', trim(file_ccf_A)
write(*,*) ' CCF file fiber B:  ', trim(file_ccf_B)

call fits_actv(lun_fits,file_ccf_A)
call singleheader_BLAZE (lun_fits,file_blaze,io_stat)
call singleheader_WAVE  (lun_fits,file_wave ,io_stat)
call singleheader_pCOEFF(lun_fits,p_coeff,io_stat)
call fits_close(lun_fits)

!getting wavelength calibration
file_rad_cal = file_wave(1:29)

file_2d_A = trim(archive_dir) // '/' // trim(file_rad_cal) // '_wave_A.fits'
file_2d_B = trim(archive_dir) // '/' // trim(file_rad_cal) // '_wave_B.fits'

write(*,*) ' WAVE file fiber A: ',trim(file_2d_A)
write(*,*) ' WAVE file fiber B: ',trim(file_2d_B)

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
    !!! Empty space replaced by workaround for missing order
    !!! in HARPS_skyCORR_chromatic.f90




end if


file_rad_cal = file_blaze(1:29)

file_2d_A = trim(archive_dir) // '/' // trim(file_rad_cal) // '_lamp_A.fits'
file_2d_B = trim(archive_dir) // '/' // trim(file_rad_cal) // '_lamp_B.fits'

write(*,*) ' LAMP file fiber A: ',trim(file_2d_A)
write(*,*) ' LAMP file fiber B: ',trim(file_2d_B)

call fits_open_s2d(lun_fits,file_2d_A,flux_A,verbose=6)
call fits_close(lun_fits)

call fits_open_s2d(lun_fits,file_2d_B,flux_B,verbose=6)
call fits_close(lun_fits)

!!! Empty space replaced by workaround for missing order
!!! in HARPS_skyCORR_chromatic.f90



allocate(ratio_AB(x_e2ds,y_e2ds), factor_B1(y_e2ds), factor_B2(y_e2ds) )

write(*,*) ' e2ds fits size: ', x_e2ds, y_e2ds

ratio_AB = flux_A/flux_B

factor_B1 = 1._PR
factor_B2 = 1._PR

write(*,*)
write(*,*) 'order, wavelength, lamp ratio, flux correction'

do ii=1,y_e2ds
    !write(*,*) ' A ',flux_A(2048:2052,y_e2ds)
    !write(*,*) ' B ', flux_B(2048:2052,y_e2ds)
    !write(*,*) ' R ', ratio_AB(2048:2052,y_e2ds)
    factor_B1(ii) = sum(ratio_AB(:, ii), dim=1)/real(x_e2ds)

    factor_B2(ii) = p_coeff(1) + p_coeff(2)*wave_B(x_e2ds/2,ii) +  p_coeff(3)*wave_B(x_e2ds/2,ii)**2 + &
    p_coeff(4)*wave_B(x_e2ds/2,ii)**3 + p_coeff(5)*wave_B(x_e2ds/2,ii)**4 + p_coeff(6)*wave_B(x_e2ds/2,ii)**5

    write(*,*) ii, wave_B(x_e2ds/2,ii), factor_B1(ii), factor_B2(ii)

end do
write(*,*)

call fits_size_2d(file_ccf_A,x_ccf,y_ccf,verbose=6)

allocate(ccf_A(x_ccf,y_ccf),ccf_B(x_ccf,y_ccf),ccf_S(x_ccf,y_ccf),ccf_O(x_ccf,y_ccf))

call fits_open_s2d(lun_fits,file_ccf_A,ccf_A,verbose=6)
call fits_close(lun_fits)


call fits_open_s2d(lun_fits,file_ccf_B,ccf_B,verbose=6)
call fits_close(lun_fits)
!!! Empty space replaced by workaround for missing order
!!! in HARPS_skyCORR_chromatic.f90




!!! the CCF profile of each order is corrected for fiber transmission and for flux correction
do ii=1,y_e2ds
    ccf_S(:,ii) = ccf_B(:,ii)* (factor_B1(ii)/factor_B2(ii))
end do

!!! each order is corrected for flux correction and stacked for the final CCF
do ii=1,x_ccf
    ccf_S(ii,y_ccf) = sum(ccf_S(ii,1:y_e2ds),dim=1)
end do

!!! contribute of the sky is removed
ccf_O = ccf_A - ccf_S

call fits_actv(lun_fits,file_ccf_A,verbose=6)

file_ccf_O =  trim(output_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A_skyR.fits'
!call fits_write_2d(lun_out,file_ccf_O,ccf_O,verbose=6)

call system('rm -f '//file_ccf_O)

blocksize = 1
call FTINIT(lun_out,file_ccf_O,blocksize,io_stat); io_stat= 0
call FTCPHD(lun_fits, lun_out, io_stat)
call FTP2DE(lun_out,0,x_ccf,x_ccf,y_ccf,ccf_O,io_stat); io_stat=0

call fits_close(lun_fits)
call fits_close(lun_out)

!!! starting the computation of the chrmatic CCFs
allocate(CCF_A_chromatic(x_ccf,chromatic_n_ranges+1))
allocate(CCF_O_chromatic(x_ccf,chromatic_n_ranges+1))

allocate(chromatic_orders(y_e2ds,chromatic_n_ranges))

do ii=starting_order,y_e2ds
    do ic=1,chromatic_n_ranges

        !if ((wave_A(x_e2ds,ii) .gt. chromatic_ranges(ic,1)) .and. (wave_A(1,ii) .lt. chromatic_ranges(ic,2))) then
        if ((wave_A(x_e2ds/2,ii) .gt. chromatic_ranges(ic,1)) .and. (wave_A(x_e2ds/2,ii) .lt. chromatic_ranges(ic,2))) then

            CCF_A_chromatic(:,ic+1) = CCF_A_chromatic(:,ic+1) + ccf_A(:,ii)
            CCF_O_chromatic(:,ic+1) = CCF_O_chromatic(:,ic+1) + ccf_O(:,ii)
            chromatic_orders(ii,ic) = ii
            ! First line is dedicated to the full CCF

        end if
    end do

   ! First line is dedicated to the full CCF
   CCF_A_chromatic(:,1) = CCF_A_chromatic(:,1) + ccf_A(:,ii)
   CCF_O_chromatic(:,1) = CCF_O_chromatic(:,1) + ccf_O(:,ii)

end do

do ic=1,chromatic_n_ranges
    write(xi, fmt='(i1)') ic

    order_list_file = trim(output_dir)//'/'//trim(file_rad_obs)//'_ccf_'//trim(file_mask)//'_A_chromatic_'// trim(xi) //'.dat'
    call system('rm -f '//order_list_file)

    open(unit=lun_out, file=order_list_file)
    write(lun_out,*) '# Order numbers are already in Python standard (starting from zero)'

    do ii=starting_order,y_e2ds
        if (chromatic_orders(ii,ic) .gt. 0) then
            write(lun_out,*) ii-1, wave_A(x_e2ds/2,ii)
        end if
    end do
    close(lun_out)
end do

call fits_actv(lun_fits,file_ccf_A,verbose=6);
sts = 0
call FTGKYE(lun_fits,'CDELT1',cdelt,comment,sts); call reset_sts(sts,0,0,'CDELT1')
call FTGKYE(lun_fits,'CRVAL1',crval,comment,sts); call reset_sts(sts,0,0,'CRVAL1')
!call FTGKYE(lun_fits,'CRPIX1',crpix,comment,sts); call reset_sts(sts,0,0,'CRPIX1')
crpix = 1
call fits_close(lun_fits)

file_ccf_O =  trim(output_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A_chromatic.fits'
call system('rm -f '//file_ccf_O)
blocksize = 1
call fits_write_2d(lun_out,file_ccf_O,CCF_A_chromatic,verbose=6)
call ccfscale_to_wcs(lun_out, crval, crpix, cdelt)
call fits_close(lun_out)

file_ccf_O =  trim(output_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A_chromatic_skyR.fits'
call system('rm -f '//file_ccf_O)
blocksize = 1
call fits_write_2d(lun_out,file_ccf_O,CCF_O_chromatic,verbose=6)
call ccfscale_to_wcs(lun_out, crval, crpix, cdelt)
call fits_close(lun_out)


deallocate(CCF_A_chromatic, CCF_O_chromatic, chromatic_orders, chromatic_ranges)

deallocate(wave_A, wave_B, flux_A, flux_B, ratio_AB)
deallocate(factor_B1, factor_B2)
deallocate(ccf_A, ccf_B, ccf_S, ccf_O)
end program
