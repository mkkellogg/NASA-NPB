
         module lhs_kernels
         implicit none
         contains



         attributes(global) subroutine krnl_lhs_init(lhs_x_dv_r,
     >                 lhsp_x_dv_r, lhsm_x_dv_r,  imax, jmax, kmax, 
     >                            imaxp, jmaxp, ni)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r

        integer, value :: ni

        integer j,k

        k = blockIDx%x
        j = blockIDx%y

        if(threadIDx%x .eq. 3) then
          lhs_x_dv_r(0,j,k,threadIDx%x) = 1.0d0
          lhsp_x_dv_r(0,j,k,threadIDx%x) = 1.0d0
          lhsm_x_dv_r(0,j,k,threadIDx%x) = 1.0d0

          lhs_x_dv_r(ni,j,k,threadIDx%x) = 1.0d0
          lhsp_x_dv_r(ni,j,k,threadIDx%x) = 1.0d0
          lhsm_x_dv_r(ni,j,k,threadIDx%x) = 1.0d0
        else
          lhs_x_dv_r(0,j,k,threadIDx%x) = 0.0d0
          lhsp_x_dv_r(0,j,k,threadIDx%x) = 0.0d0
          lhsm_x_dv_r(0,j,k,threadIDx%x) = 0.0d0

          lhs_x_dv_r(ni,j,k,threadIDx%x) = 0.0d0
          lhsp_x_dv_r(ni,j,k,threadIDx%x) = 0.0d0
          lhsm_x_dv_r(ni,j,k,threadIDx%x) = 0.0d0
        endif      

        end subroutine krnl_lhs_init







        attributes(global) subroutine krnl_lhs_initj(
     >       lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r, 
     >         imax, jmax, kmax, imaxp, jmaxp, nj)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r


        integer, value :: nj

        integer i,k 

        k = blockIDx%x
        i = blockIDx%y

        if(threadIDx%x .eq. 3 ) then
       
          lhs_x_dv_r(i,0,k,3) = 1.0d0
          lhsp_x_dv_r(i,0,k,3) = 1.0d0
          lhsm_x_dv_r(i,0,k,3) = 1.0d0
          lhs_x_dv_r(i,nj,k,3) = 1.0d0
          lhsp_x_dv_r(i,nj,k,3) = 1.0d0
          lhsm_x_dv_r(i,nj,k,3) = 1.0d0

        else

             lhs_x_dv_r(i,0,k,threadIDx%x) = 0.0d0
             lhsp_x_dv_r(i,0,k,threadIDx%x) = 0.0d0
             lhsm_x_dv_r(i,0,k,threadIDx%x) = 0.0d0
             lhs_x_dv_r(i,nj,k,threadIDx%x) = 0.0d0
             lhsp_x_dv_r(i,nj,k,threadIDx%x) = 0.0d0
             lhsm_x_dv_r(i,nj,k,threadIDx%x) = 0.0d0
        endif

        end subroutine krnl_lhs_initj













        attributes(global) subroutine krnl_lhs_initk(
     >     lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r,  
     >     imax, jmax, kmax, imaxp, jmaxp, nk)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r

        integer, value :: nk

        integer i,j

        j = blockIDx%x
        i = blockIDx%y

        if(threadIDx%x .eq. 3 ) then
            lhs_x_dv_r(i,j,0,3) = 1.0d0
             lhsp_x_dv_r(i,j,0,3) = 1.0d0
             lhsm_x_dv_r(i,j,0,3) = 1.0d0
             lhs_x_dv_r(i,j,nk,3) = 1.0d0
             lhsp_x_dv_r(i,j,nk,3) = 1.0d0
             lhsm_x_dv_r(i,j,nk,3) = 1.0d0
        else
             lhs_x_dv_r(i,j,0,threadIDx%x) = 0.0d0
             lhsp_x_dv_r(i,j,0,threadIDx%x) = 0.0d0
             lhsm_x_dv_r(i,j,0,threadIDx%x) = 0.0d0
             lhs_x_dv_r(i,j,nk,threadIDx%x) = 0.0d0
             lhsp_x_dv_r(i,j,nk,threadIDx%x) = 0.0d0
             lhsm_x_dv_r(i,j,nk,threadIDx%x) = 0.0d0
        endif
        

        end subroutine krnl_lhs_initk



         end module lhs_kernels


