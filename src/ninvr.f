
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine  ninvr
       use glob
       use cudafor
       use ninvr_kernels
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   block-diagonal matrix-vector multiplication              
c---------------------------------------------------------------------

       include 'header.h'

       integer  i, j, k
       double precision r1, r2, r3, r4, r5, t1, t2
       integer kSize,jSize, sharedSize, ni, nj,
     >         blockSize, istat, upperi1, upperi2, fullSize, iSize,
     >         subSections, istart

       if(doGPU .gt. 0) then !BEGIN GPU SWITCH

        

       iSize = nx2
       jSize = ny2
       kSize = nz2       

       blockSize = fixedBlockSize
       subSections = iSize/blockSize
       if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
       blockSize = (blockSize+2)

       dimGrid = dim3(kSize,jSize*subSections,1)
       dimBlock = dim3(blockSize,1,1)


       call krnl_ninvr<<<dimGrid,dimBlock>>>(
     >               imax, jmax, kmax, imaxp, jmaxp, rhs_dv_r,
     >                  bt, iSize, subSections)  

        


       else !BEGIN NON-GPU MODE


         if (timeron) call timer_start(t_ninvr)
         do k = 1, nz2
          do j = 1, ny2
             do i = 1, nx2

                r1 = rhs(1,i,j,k)
                r2 = rhs(2,i,j,k)
                r3 = rhs(3,i,j,k)
                r4 = rhs(4,i,j,k)
                r5 = rhs(5,i,j,k)
               
                t1 = bt * r3
                t2 = 0.5d0 * ( r4 + r5 )

                rhs(1,i,j,k) = -r2
                rhs(2,i,j,k) =  r1
                rhs(3,i,j,k) = bt * ( r4 - r5 )
                rhs(4,i,j,k) = -t1 + t2
                rhs(5,i,j,k) =  t1 + t2
             enddo    
          enddo
         enddo
         if (timeron) call timer_stop(t_ninvr)


       endif !END GPU SWITCH
 

       return
       end
