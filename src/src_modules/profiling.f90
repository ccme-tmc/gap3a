module profiling

  type timer
    real(8) :: total_time = 0.0D0
    real(8) :: record = 0.0D0
  end type timer

  contains

    subroutine timer_start(tm)
      type (timer) :: tm
      call cpu_time(tm%record)
    end subroutine timer_start
    
    subroutine timer_end(tm)
      type (timer) :: tm
      real(8) :: now
      call cpu_time(now)
      if (tm%record.ge.1.0D-8) then
        tm%total_time = tm%total_time + now - tm%record
        tm%record = 0.0D0
      endif
    end subroutine timer_end

end module profiling
