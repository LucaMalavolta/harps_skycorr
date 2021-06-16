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


!!! New variables added for chromatic RVs
integer :: chromatic_n_ranges, nw, ic
real(PR) :: wave_start, wave_end
real (PR), dimension(:,:), allocatable :: CCF_A_chromatic, CCF_O_chromatic, chromatic_ranges
character (len=nch_file) :: chromatic_list
real :: crval, crpix, cdelt
character (len=64) :: comment
integer, dimension(:,:), allocatable :: chromatic_orders
character (len=nch_file) :: order_list_file
character (len=1) :: xi


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
    write(*,*) ' 9) optional: chromatic ranges'
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

!!! New optional input for chromatic RVs
if  (iargc() .ge. 9) then
    call getarg(9,chromatic_list)
else
    chromatic_list = 'none'
end if

!! file_rad is the star observation. we need the rad of the blaze file

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

    write(*,*) ' CCF file fiber A:  ', trim(file_ccf_A)
    write(*,*) ' CCF file fiber B:  ', trim(file_ccf_B)

    call fits_actv(lun_fits,file_ccf_A)
    call singleheader_BLAZE (lun_fits,file_blaze,io_stat)
    call singleheader_WAVE  (lun_fits,file_wave ,io_stat)
    call singleheader_pCOEFF(lun_fits,p_coeff,io_stat)
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
    end if


    file_rad_cal = file_blaze(1:29)

    if (date_calib.eq.1) then
        file_2d_A = trim(calib_dir) // '/' // trim(file_date)// '/' // trim(file_rad_cal) // '_lamp_A.fits'
        file_2d_B = trim(calib_dir) // '/' // trim(file_date)// '/' // trim(file_rad_cal) // '_lamp_B.fits'
    else
        file_2d_A = trim(calib_dir) // '/' // trim(file_rad_cal) // '_lamp_A.fits'
        file_2d_B = trim(calib_dir) // '/' // trim(file_rad_cal) // '_lamp_B.fits'
    end if

    write(*,*) ' LAMP file fiber A: ',trim(file_2d_A)
    write(*,*) ' LAMP file fiber B: ',trim(file_2d_B)

    call fits_open_s2d(lun_fits,file_2d_A,flux_A,verbose=6)
    call fits_close(lun_fits)

    call fits_open_s2d(lun_fits,file_2d_B,flux_B,verbose=6)
    call fits_close(lun_fits)

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

    call fits_size_2d(file_ccf_A,x_ccf,y_ccf,verbose=6)

    allocate(ccf_A(x_ccf,y_ccf),ccf_B(x_ccf,y_ccf),ccf_S(x_ccf,y_ccf),ccf_O(x_ccf,y_ccf))


    call fits_open_s2d(lun_fits,file_ccf_A,ccf_A,verbose=6)
    call fits_close(lun_fits)

    call fits_open_s2d(lun_fits,file_ccf_B,ccf_B,verbose=6)
    call fits_close(lun_fits)

    !!! the CCF profile of each order is corrected for fiber transmission and for flux correction
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


        if (date_output.eq.1) then
            order_list_file =  trim(output_dir) // '/' // trim(file_date) // '/' // trim(file_rad_obs) // &
              '_ccf_'//trim(file_mask)//'_A_chromatic_'// trim(xi) //'.dat'
        else
            order_list_file =  trim(output_dir) // '/' // trim(file_rad_obs) // &
              '_ccf_'//trim(file_mask)//'_A_chromatic_'// trim(xi) //'.dat'
        end if

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

    if (date_output.eq.1) then
        order_list_file =  trim(output_dir) // '/' // trim(file_date) // '/' // trim(file_rad_obs) // &
          '_ccf_'//trim(file_mask)//'_A_chromatic.fits'
    else
        order_list_file =  trim(output_dir) // '/' // trim(file_date) // '/' // trim(file_rad_obs) // &
          '_ccf_'//trim(file_mask)//'_A_chromatic.fits'
    end if

    call system('rm -f '//file_ccf_O)
    blocksize = 1
    call fits_write_2d(lun_out,file_ccf_O,CCF_A_chromatic,verbose=6)
    call ccfscale_to_wcs(lun_out, crval, crpix, cdelt)
    call fits_close(lun_out)

    if (date_output.eq.1) then
        order_list_file =  trim(output_dir) // '/' // trim(file_date) // '/' // trim(file_rad_obs) // &
          '_ccf_'//trim(file_mask)//'_A_chromatic_skyR.fits'
    else
        order_list_file =  trim(output_dir) // '/' // trim(file_date) // '/' // trim(file_rad_obs) // &
          '_ccf_'//trim(file_mask)//'_A_chromatic_skyR.fits'
    end if

    file_ccf_O =  trim(output_dir) // '/' // trim(file_rad_obs) // '_ccf_'//trim(file_mask)//'_A_chromatic_skyR.fits'
    call system('rm -f '//file_ccf_O)
    blocksize = 1
    call fits_write_2d(lun_out,file_ccf_O,CCF_O_chromatic,verbose=6)
    call ccfscale_to_wcs(lun_out, crval, crpix, cdelt)
    call fits_close(lun_out)

    deallocate(CCF_A_chromatic, CCF_O_chromatic, chromatic_orders)

    deallocate(wave_A, wave_B, flux_A, flux_B, ratio_AB)
    deallocate(factor_B1, factor_B2)
    deallocate(ccf_A, ccf_B, ccf_S, ccf_O)

end do

deallocate(chromatic_ranges)

close(lun_input)
end program
