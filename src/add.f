
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine  add
       use glob
       use cudafor
       use add_kernels
       use util_kernels
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c addition of update to the vector u
c---------------------------------------------------------------------

       include 'header.h'

       integer i, j, k, m, istat, haredSize, iSize,jSize,kSize,
     >         istart, jstart, kstart, blockSize, numBlocks,
     >         ierr, subISize, subSections, sharedSize
       double precision t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3,
     >                  r4, r5, ac2inv
       real tink


       if(doGPU .gt. 0 )then !BEGIN GPU SWITCH

         if (timeron) call timer_start(t_add)

         iSize = nx2
         jSize = ny2
         kSize = nz2

         blockSize = fixedBlockSize
         subSections = iSize/blockSize
         if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1

         dimGrid = dim3(kSize, jSize*subSections,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(krnlStart,0)
         call krnl_add<<<dimGrid,dimBlock>>>(
     >    u_dv_r, rhs_dv_r,imax,  jmax, kmax, imaxp, jmaxp,
     >                         iSize, subSections)
         istat = cudaEventRecord(krnlEnd,0)
       
        if (timeron) call timer_stop(t_add)


        if(syncForTimers .gt. 0) then
          istat= cudaThreadSynchronize()
          istat = cudaEventElapsedTime(tink, krnlStart, krnlEnd)
          call timer_inc(t_kernel,dble(tink/1000.0))
          call timer_inc(t_add,dble(tink/1000.0))

        endif

         




       else !BEGIN NON-GPU MODE




       if (timeron) call timer_start(t_add)
       do k = 1, nz2
          do j = 1, ny2
             do i = 1, nx2
                do m = 1, 5
                   u(m,i,j,k) = u(m,i,j,k) + rhs(m,i,j,k)
                end do
             end do
          end do
       end do
       if (timeron) call timer_stop(t_add)



       endif !END GPU MODE



       return
       end

