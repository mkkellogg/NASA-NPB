
        module ninvr_kernels
        contains

        attributes(global) subroutine krnl_ninvr(
     >              imax, jmax, kmax, imaxp, jmaxp, rhs_dv_r,
     >                  bt, iSize, subSections)

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r

       double precision,value :: bt
       integer, value :: iSize, subSections

        integer i,j,k,  blockSize, ioff

        double precision :: r1,r2,r3,r4,r5,t1,t2

        double precision, shared :: temp(*)


        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        ioff = (mod(blockIDx%y-1,subSections)*blockSize)
        i = threadIDx%x + ioff 



          if (i .le. iSize) then
          

                r1 = rhs_dv_r(i,j,k,1)
                r2 = rhs_dv_r(i,j,k,2)
                r3 = rhs_dv_r(i,j,k,3)
                r4 = rhs_dv_r(i,j,k,4)
                r5 = rhs_dv_r(i,j,k,5)

                t1 = bt * r3
                t2 = 0.5d0 * ( r4 + r5 )

                rhs_dv_r(i,j,k,1) = -r2
                rhs_dv_r(i,j,k,2) =  r1
                rhs_dv_r(i,j,k,3) = bt * ( r4 - r5 )
                rhs_dv_r(i,j,k,4) = -t1 + t2
                rhs_dv_r(i,j,k,5) =  t1 + t2



          endif

        end subroutine krnl_ninvr


        end module ninvr_kernels

