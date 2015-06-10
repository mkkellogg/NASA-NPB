        module tzetar_kernels
        contains

        attributes(global) subroutine krnl_tzetar(
     >              imax, jmax, kmax, imaxp, jmaxp, rhs_dv_r,
     >                  u_dv_r, us_dv, vs_dv, ws_dv, qs_dv, speed_dv,
     >                  bt, c2iv, iSize, subSections)

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: us_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: vs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: ws_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: qs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1) :: speed_dv


       double precision,value :: bt, c2iv
       integer, value :: iSize, subSections

        integer i,j,k, blockSize

        double precision r1,r2,r3,r4,r5,t1,t2,t3
        double precision xvel, yvel, zvel, ac, ac2u, uzik1, btuz       

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)

          if (i .le. iSize) then


                xvel = us_dv(i,j,k)
                yvel = vs_dv(i,j,k)
                zvel = ws_dv(i,j,k)
                ac   = speed_dv(i,j,k)

                ac2u = ac*ac

                r1 = rhs_dv_r(i,j,k,1)
                r2 = rhs_dv_r(i,j,k,2)
                r3 = rhs_dv_r(i,j,k,3)
                r4 = rhs_dv_r(i,j,k,4)
                r5 = rhs_dv_r(i,j,k,5)      

                uzik1 = u_dv_r(i,j,k,1)
                btuz  = bt * uzik1

                t1 = btuz/ac * (r4 + r5)
                t2 = r3 + t1
                t3 = btuz * (r4 - r5)

                rhs_dv_r(i,j,k,1) = t2
                rhs_dv_r(i,j,k,2) = -uzik1*r2 + xvel*t2
                rhs_dv_r(i,j,k,3) =  uzik1*r1 + yvel*t2
                rhs_dv_r(i,j,k,4) =  zvel*t2  + t3
                rhs_dv_r(i,j,k,5) =  uzik1*(-xvel*r2 + yvel*r1) + 
     >                    qs_dv(i,j,k)*t2 + c2iv*ac2u*t1 + zvel*t3
           

          endif

      
 
        end subroutine krnl_tzetar








        end module tzetar_kernels

