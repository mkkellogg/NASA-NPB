
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine pinvr
       use glob
       use cudafor
       use pinvr_kernels
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   block-diagonal matrix-vector multiplication                       
c---------------------------------------------------------------------

       include 'header.h'

       integer i, j, k
       double precision r1, r2, r3, r4, r5, t1, t2
       integer kSize,jSize, sharedSize, ni, nj,
     >         blockSize, istat, upperi1, upperi2, iSize,
     >         subSections, istart, subISize
      
       if(doGPU .gt. 0) then !BEGIN GPU SWITCH

         iSize = nx2
         jSize = ny2
         kSize = nz2

         blockSize = fixedBlockSize
         subSections = iSize/blockSize
         if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1

         dimGrid = dim3(kSize,jSize*subSections,1)
         dimBlock = dim3(blockSize,1,1)

         call krnl_pinvr<<<dimGrid,dimBlock>>>(
     >               imax, jmax, kmax, imaxp, jmaxp, rhs_dv_r,
     >                  bt, iSize, subSections)  




       else !BEGIN NON-GPU MODE



        do   k = 1, nz2
          do   j = 1, ny2
             do   i = 1, nx2

                r1 = rhs(1,i,j,k)
                r2 = rhs(2,i,j,k)
                r3 = rhs(3,i,j,k)
                r4 = rhs(4,i,j,k)
                r5 = rhs(5,i,j,k)

                t1 = bt * r1
                t2 = 0.5d0 * ( r4 + r5 )

                rhs(1,i,j,k) =  bt * ( r4 - r5 )
                rhs(2,i,j,k) = -r3
                rhs(3,i,j,k) =  r2
                rhs(4,i,j,k) = -t1 + t2
                rhs(5,i,j,k) =  t1 + t2
             end do
          end do
         end do
     


       endif


     

       return
       end



