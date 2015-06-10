
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine  txinvr
       use glob
       use cudafor
       use txinvr_kernels
       use util_kernels

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c block-diagonal matrix-vector multiplication                  
c---------------------------------------------------------------------

       include 'header.h'

       integer i, j, k, istat, haredSize, iSize,jSize,kSize,
     >         istart, jstart, kstart, blockSize, numBlocks,
     >         ierr, subISize, subSections, sharedSize
       double precision t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3, 
     >                  r4, r5, ac2inv
       real tink

       if (doGPU .gt. 0) then !BEGIN GPU SWITCH

           if (timeron) call timer_start(t_txinvr)

           istat = cudaEventRecord(krnlStart,0)

           iSize = nx2
           jSize = ny2
           kSize = nz2

           blockSize = fixedBlockSize     
           subSections = iSize/blockSize
           if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
           dimGrid = dim3(kSize,jSize*subSections,1)
           dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(krnlStart,0)
         call krnl_txinvr<<<dimGrid,dimBlock>>>(
     >         u_dv_r, us_dv, vs_dv, ws_dv,  qs_dv, rho_i_dv,
     >                     speed_dv, square_dv, rhs_dv_r,
     >                      imax,  jmax, kmax, imaxp, jmaxp,
     >                      c2, bt, iSize, subSections)
         istat = cudaEventRecord(krnlEnd,0)

        if (timeron) call timer_stop(t_txinvr)


        if(syncForTimers .gt. 0) then
          istat= cudaThreadSynchronize()
          istat = cudaEventElapsedTime(tink, krnlStart, krnlEnd)
          call timer_inc(t_kernel,dble(tink/1000.0))
          call timer_inc(t_txinvr,dble(tink/1000.0))
        endif





       else !BEGIN NON-GPU MODE




        if (timeron) call timer_start(t_txinvr)
       do    k = 1, nz2
          do    j = 1, ny2
             do    i = 1, nx2

                ru1 = rho_i(i,j,k)
                uu = us(i,j,k)
                vv = vs(i,j,k)
                ww = ws(i,j,k)
                ac = speed(i,j,k)
                ac2inv = ac*ac

                r1 = rhs(1,i,j,k)
                r2 = rhs(2,i,j,k)
                r3 = rhs(3,i,j,k)
                r4 = rhs(4,i,j,k)
                r5 = rhs(5,i,j,k)

                t1 = c2 / ac2inv * ( qs(i,j,k)*r1 - uu*r2  - 
     >                  vv*r3 - ww*r4 + r5 )
                t2 = bt * ru1 * ( uu * r1 - r2 )
                t3 = ( bt * ru1 * ac ) * t1

                rhs(1,i,j,k) = r1 - t1
                rhs(2,i,j,k) = - ru1 * ( ww*r1 - r4 )
                rhs(3,i,j,k) =   ru1 * ( vv*r1 - r3 )
                rhs(4,i,j,k) = - t2 + t3
                rhs(5,i,j,k) =   t2 + t3

             end do
          end do
       end do
       if (timeron) call timer_stop(t_txinvr) 
       
       endif



      

       return
       end


