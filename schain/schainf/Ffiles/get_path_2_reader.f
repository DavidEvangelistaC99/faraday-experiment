
      fqual_temp = TRIM(ppath)//TRIM(cpath)
      the_len = LEN_TRIM(fqual_temp)

      return
      end

      subroutine get_path_reader(fqual,the_len)
c
c      create table of magnetic field components
c
c      character*47 L

      integer the_len
      character(the_len) :: fqual
      character(1024) :: fqual_temp

      call get_len(the_len,fqual_temp)
      fqual = TRIM(fqual_temp)
      write(*,*) "Done"

      return
      end

