module filepaths
  use common
  implicit none

  character (len=*), parameter :: archive_raw_harpn='/media/TR2/malavolta/HARPN/data/'
  character (len=*), parameter :: archive_red_harpn='/media/TR2/malavolta/HARPN/data/'

  character (len=*), parameter :: archive_raw_harps='/media/TR2/malavolta/HARPS/data/'
  character (len=*), parameter :: archive_red_harps='/media/TR2/malavolta/HARPS/data/'
  
  character (len=*), parameter :: &
       make_atms_m4='/home/malavolta/CODE/PreMoog/modprog_1/makekurandy/makekurucz3_cl_aM4.go', &
       make_atms='/home/malavolta/CODE/PreMoog/modprog_1/makekurandy/makekurucz3_cl.go', &
       moog_call='/home/malavolta/CODE/MOOG/MOOG_1/MOOGSILENT > MOOG_screen_output.dat', &
       moog2013_call='/home/malavolta/CODE/MOOG/MOOG2013_1/MOOGSILENT > MOOG_screen_output.dat'

  character (len=*), parameter :: &
       code_path='/home/malavolta/CODE/' , &
       routinesF90_path='/home/malavolta/CODE/Routines_F90/' , &
       tables_path='/home/malavolta/CODE/Routines_F90/tables/' , &
       masks_path='/home/malavolta/CODE/Routines_F90/masks/' , &
       stellib_path='/media/TR1/malavolta/STELLAR_LIBS/'

  character (len=*), parameter :: GAPS_mastercat= trim(tables_path) // 'gaps-mastercatalog-v7.7.cat'
  
  character (len=*), parameter :: coelho2005_sun= trim(tables_path) // 'Coelho2005_5750_45_p00p00.ms.fits.gz'
  
  character (len=*), parameter :: work_dir_common = '/home/malavolta/moog_workdir/wd_1/'
  character (len=*), parameter :: temp_file1='file_temp1.dat', temp_file2='file_tmp2.dat', temp_file3='file_tmp3.dat'

  character (len=*), parameter :: batch_file='batch_1.par', buffer_file='buffer_1.source'

  integer :: nfiles,nopts
  character (len=nch_file), allocatable :: filelist(:), optslist(:)
  character (len=nch_rad ), allocatable :: prelist(:) , radlist(:)
  
contains
  subroutine int2chr4(int,chr4)
    integer :: int
    character (len=4) :: chr4
    
    chr4 = '----'
    if (int.ge.0    .and.int.lt.10   ) write(chr4,'(A3,I1)') '000', int
    if (int.ge.10   .and.int.lt.100  ) write(chr4,'(A2,I2)') '00' , int
    if (int.ge.100  .and.int.lt.1000 ) write(chr4,'(A1,I3)') '0'  , int
    if (int.ge.1000 .and.int.lt.1000 ) write(chr4,'(A1,I3)') '0'  , int
    if (int.gt.1000                  ) write(chr4,'(   I4)')        int

    if (int.lt.0    .and.int.gt.-10  ) write(chr4,'(A3,I1)') '-00', abs(int)
    if (int.le.-10  .and.int.gt.-100 ) write(chr4,'(A2,I2)') '-0' , abs(int)
    if (int.le.-100 .and.int.gt.-1000) write(chr4,'(A1,I3)') '-'  , abs(int)
    return
  end subroutine int2chr4

  subroutine read_list(nopts_t,fits)
    !     It reads the first $nopts command line arguments simply as strings
    !     and the following $nargs-$nopts arguments as ".fits" filenames.
    !     If the latter arguments do not have a ".fits" suffix it generates
    !     an error. If the arguments are a pure fits filelist, use the
    !     similar readfl routine. The only input parameter is $nopts.

    character (len=nch_file) :: filename
    integer, intent(in) :: nopts_t
    logical, intent(in), optional :: fits
    logical :: fits_t
    integer :: i, n,ind
    integer :: nargs, nslash

    !     nargs = number of arguments
    !     nopts = number of options (at the beginning of the argument list)
    !     nfiles = number of files (at the end of the argument list)
    !     nargs = nopts + nfiles

    !     pre()    = list of file prefixes (only for >$nopts indexes)
    !     rad()    = list of file radicals (only for >$nopts indexes)
    !     prelen() = list of file prefix lengths (only for >$nopts indexes)
    !     radlen() = list of file radicl lengths (only for >$nopts indexes)

    !     NOTE: the optslist() is a character array, so in order to read the
    !     $nopts arguments ad numbers, one must use an internal read. E.g.
    !     read (optslist(1), *) int
    !     for an integer argument. This is also a standard f77 feature.

    fits_t = .true.
    if (present(fits)) fits_t=fits
    nopts=nopts_t
    nargs = iargc()
    nfiles = nargs - nopts

    if (nargs == 0) then
       write(*,*) 'ERROR: no arguments'
       stop
    endif

    if (nargs .lt. (nopts+1)) then
       write(*,*) 'ERROR: too few arguments (<', nopts+1, ')'
       stop
    endif

    if (nopts.gt.0)  then
       allocate(optslist(nopts))
       do i=1, nopts
          call getarg(i,optslist(i))
       enddo
    endif

    allocate(filelist(nfiles),prelist(nfiles),radlist(nfiles))

    do i=1, nfiles
       call getarg(i+nopts,filelist(i))

       nslash = 0
       filename = filelist(i)

       if (filename(1:4).eq.'NULL') then
          radlist(i)='NULL'
          prelist(i)='NULL'
          cycle
       endif

       ind=0
       do n=1, len(filename)
          if (filename(n:n).eq.'/') nslash = n

          if ( (fits_t.eqv. .true.) .and. filename(n:n+4).eq.'.fits') then
             ind=n
             exit
          endif

          if  ((fits_t.eqv. .false.) .and. filename(n:n).eq.'.') ind=n

       enddo

       if (ind == 0) then
          write (*,*) 'ERROR: ', filename
          write (*,*) 'cannot find the .fits extension'
          stop
       endif

       radlist(i) = filename((nslash+1):(ind-1))

       if (nslash .ne. 0) then
          prelist(i)      = filename(1:nslash)
       else
          prelist(i)      = './'
       endif

    enddo

    return
  end subroutine read_list


  subroutine delete_list()

    if (allocated(optslist)) deallocate(optslist)
    if (allocated(prelist )) deallocate(prelist )
    if (allocated(radlist )) deallocate(radlist )
    if (allocated(filelist)) deallocate(filelist)

  end subroutine delete_list


end module filepaths
