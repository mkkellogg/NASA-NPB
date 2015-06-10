         module txinvr_kernels
         implicit none
         contains




        attributes(global) subroutine krnl_txinvr(
     >      u_dv_r, us_dv, vs_dv, ws_dv,  qs_dv,rho_i_dv,
     >                     speed_dv, square_dv, rhs_dv_r,
     >                      imax,  jmax, kmax, imaxp, jmaxp,
     >                      c2, bt,iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP, 0:JMAXP,
     >                              0:KMAX-1) :: rho_i_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: us_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: vs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: ws_dv
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1) :: square_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: qs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1) :: speed_dv
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r

         double precision, value :: c2, bt

        integer, value :: iSize, subSections
       
         integer :: i,j,k,  modF, blockSize, m
         double precision t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3, 
     >                  r4, r5, ac2inv 

        double precision, shared :: temp(*) 


        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)

        if(i .le. iSize) then

                ru1 = rho_i_dv(i,j,k)
                uu = us_dv(i,j,k)
                vv = vs_dv(i,j,k)
                ww = ws_dv(i,j,k)
                ac = speed_dv(i,j,k)
                ac2inv = ac*ac

                r1 = rhs_dv_r(i,j,k,1)
                r2 = rhs_dv_r(i,j,k,2)
                r3 = rhs_dv_r(i,j,k,3)
                r4 = rhs_dv_r(i,j,k,4)
                r5 = rhs_dv_r(i,j,k,5)

                t1 = c2 / ac2inv * ( qs_dv(i,j,k)*r1 - uu*r2  - 
     >                  vv*r3 - ww*r4 + r5 )
                t2 = bt * ru1 * ( uu * r1 - r2 )
                t3 = ( bt * ru1 * ac ) * t1

                rhs_dv_r(i,j,k,1) = r1 - t1
                rhs_dv_r(i,j,k,2) = - ru1 * ( ww*r1 - r4 )
                rhs_dv_r(i,j,k,3) =   ru1 * ( vv*r1 - r3 )
                rhs_dv_r(i,j,k,4) = - t2 + t3
                rhs_dv_r(i,j,k,5) =   t2 + t3

        endif


         end subroutine krnl_txinvr




         end module txinvr_kernels
