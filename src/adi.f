
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine  adi
        include  'header.h'
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       call compute_rhs(0)

       call txinvr

       call x_solve

       call y_solve

       call z_solve

       call add
       
       return
       end
