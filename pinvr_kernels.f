        module pinvr_kernels
        contains

        attributes(global) subroutine krnl_pinvr(
     >              imax, jmax, kmax, imaxp, jmaxp, rhs_dv_r,
     >                  bt, iSize, subSections)

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r

       double precision,value :: bt
       integer, value :: iSize, subSections

        integer i,j,k, blockSize, joff

        double precision r1,r2,r3,r4,r5,t1,t2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        joff = (mod(blockIDx%y-1,subSections)*blockSize)
        i = threadIDx%x + joff 
       
        if (i .le. iSize) then

                r1 = rhs_dv_r(i,j,k,1)
                r2 = rhs_dv_r(i,j,k,2)
                r3 = rhs_dv_r(i,j,k,3)
                r4 = rhs_dv_r(i,j,k,4)
                r5 = rhs_dv_r(i,j,k,5)

                t1 = bt * r1
                t2 = 0.5d0 * ( r4 + r5 )

                rhs_dv_r(i,j,k,1) =  bt * ( r4 - r5 )
                rhs_dv_r(i,j,k,2) = -r3
                rhs_dv_r(i,j,k,3) =  r2
                rhs_dv_r(i,j,k,4) = -t1 + t2
                rhs_dv_r(i,j,k,5) =  t1 + t2
           
        endif
        
        end subroutine krnl_pinvr


        end module pinvr_kernels

