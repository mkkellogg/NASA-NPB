
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine x_solve
       use glob
       use cudafor
       use lhs_kernels
       use x_solve_kernels
       use util_kernels
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c this function performs the solution of the approximate factorization
c step in the x-direction for all five matrix components
c simultaneously. The Thomas algorithm is employed to solve the
c systems for the x-lines. Boundary conditions are non-periodic
c---------------------------------------------------------------------

       include 'header.h'

       integer i, j, k, i1, i2, m, kSize,jSize, sharedSize, ni, nj,
     >         blockSize, istat, upperi1, upperi2, iSize,
     >         subSections, istart, subISize, subJSize, numBlocks,
     >         blockDimSize, iBlocks, jBlocks
       double precision  ru1, fac1, fac2
       real tink

c---------------------------------------------------------------------
c---------------------------------------------------------------------

       if(doGPU .gt. 0) then !BEGIN GPU SWITCH

       if (timeron) call timer_start(t_xsolve)



c---------------------------------------------------------------------
c Computes the left hand side for the three x-factors  
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      first fill the lhs for the u-eigenvalue                   
c---------------------------------------------------------------------
    

       upperi1 = grid_points(1)-1
       upperi2 = nx2
       jSize = ny2
       kSize = nz2
       iSize = nx2

      
       blockSize = fixedBlockSize
       subSections = iSize/blockSize
       if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
       blockSize = (blockSize+2)

       dimGrid = dim3(kSize,jSize*subSections,1)
       dimBlock = dim3(blockSize,1,1)
       sharedSize = blockSize * sizeOfDouble * 2 

       istat = cudaEventRecord(dynamicEvents(1),0)
       call krnl_x_solve_init<<<dimGrid,dimBlock,sharedSize>>>(
     >    lhs_x_dv_r, rho_i_dv, us_dv, imax, jmax, kmax, imaxp, jmaxp, 
     >                         dttx1, dttx2, c2dttx1,
     >                        con43, c1c5, c3c4, dx1, dx2, dx5,dxmax,
     >                      upperi1, upperi2, subSections)  
       istat = cudaEventRecord(dynamicEvents(2),0)




c--------------------------------------------------------------------
c      add fourth order dissipation
c---------------------------------------------------------------------
            

       jSize = ny2
       kSize = nz2
       istart = -1

       blockSize = fixedBlockSize
       subSections = jSize/blockSize
       if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1

       dimGrid = dim3(kSize*subSections,1,1)
       dimBlock = dim3(blockSize,1,1)

       istat = cudaEventRecord(dynamicEvents(3),0)
       call krnl_x_solve_dis1<<<dimGrid,dimBlock>>>(
     >           lhs_x_dv_r, imax, jmax, kmax, imaxp, jmaxp, 
     >             comz1,comz4,comz5,comz6, istart,jSize,  subSections)
       istat = cudaEventRecord(dynamicEvents(4),0)




       iSize = (grid_points(1)-4) - 3 + 1
       jSize = ny2
       kSize = nz2

       blockSize = fixedBlockSize
       subSections= (iSize/blockSize)
       if(mod(iSize,blockSize) .ne. 0)subSections = subSections+1

       dimGrid = dim3(kSize,jSize*subSections,1)
       dimBlock = dim3(blockSize,1,1)

       istat = cudaEventRecord(dynamicEvents(5),0)
       call krnl_x_solve_dis2<<<dimGrid,dimBlock>>>(
     >               imax, jmax, kmax, imaxp, jmaxp, lhs_x_dv_r,
     >                  comz1,comz4,comz5,comz6, iSize, subSections)  
       istat = cudaEventRecord(dynamicEvents(6),0)
      


       jSize = ny2
       kSize = nz2
       istart = grid_points(1)-3

       blockSize = fixedBlockSize
       subSections = jSize/blockSize
       if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1

       dimGrid = dim3(kSize*subSections,1,1)
       dimBlock = dim3(blockSize,1,1)

       istat = cudaEventRecord(dynamicEvents(7),0)
       call krnl_x_solve_dis1<<<dimGrid,dimBlock>>>(
     >      lhs_x_dv_r,  imax, jmax, kmax, imaxp, jmaxp, 
     >           comz1,comz4,comz5,comz6,istart,jSize,  subSections)
       istat = cudaEventRecord(dynamicEvents(8),0)



c---------------------------------------------------------------------
c      subsequently, fill the other factors (u+c), (u-c) by adding to 
c      the first  
c---------------------------------------------------------------------

       iSize = nx2
       jSize = ny2
       kSize = nz2

       blockSize = fixedBlockSize
       subSections = jSize/blockSize
       if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1
       blockSize = blockSize+2

       dimGrid = dim3(kSize,jSize*subSections,1)
       dimBlock = dim3(blockSize,1,1)
       
       istat = cudaEventRecord(dynamicEvents(9),0)
       call krnl_x_solve_postdis<<<dimGrid,dimBlock>>>(
     >        lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r,speed_dv,
     >                     imax, jmax, kmax, imaxp, jmaxp, 
     >                   dttx2, iSize, subSections)  
       istat = cudaEventRecord(dynamicEvents(10),0)


c---------------------------------------------------------------------
c         swap 'i' and 'j' dimensions from rhs_dv_r to rhs_dv_j_r,
c         from lhs_x_dv_r to lhs_x_dv_j_r,from lhsp_x_dv_r to lhsp_x_dv_j_r,
c         and from lhsm_x_dv_r to lhsm_x_dv_j_r
c--------------------------------------------------------------------

        iSize = IMAXP+1
        jSize = JMAXP+1
        kSize = KMAX

        blockDimSize = 31
        iBlocks = iSize/blockDimSize
	if(mod(iSize,blockDimSize) .ne. 0)iBlocks=iBlocks+1
        jBlocks = jSize/blockDimSize
        if(mod(jSize,blockDimSize) .ne. 0)jBlocks=jBlocks+1

        dimGrid = dim3(kSize,iBlocks*jBlocks,1)
        dimBlock = dim3(blockDimSize,blockDimSize,1)
        sharedSize = blockDimSize * blockDimSize * sizeOfDouble 

        istat = cudaEventRecord(dynamicEvents(11),0)
        call krnl_util_swap4D_itoj_c<<<dimGrid,dimBlock,sharedSize>>>(
     >                         lhs_x_dv_r,lhs_x_dv_j_r,
     >                         imax, jmax, kmax, imaxp, imaxp,
     >                         blockDimSize, iBlocks, jBlocks)
        istat = cudaEventRecord(dynamicEvents(12),0)

        istat = cudaEventRecord(dynamicEvents(13),0)
        call krnl_util_swap4D_itoj_c<<<dimGrid,dimBlock,sharedSize>>>(
     >                         lhsp_x_dv_r,lhsp_x_dv_j_r,
     >                         imax, jmax, kmax, imaxp, imaxp,
     >                         blockDimSize, iBlocks, jBlocks)
        istat = cudaEventRecord(dynamicEvents(14),0)

        istat = cudaEventRecord(dynamicEvents(15),0)
        call krnl_util_swap4D_itoj_c<<<dimGrid,dimBlock,sharedSize>>>(
     >                          lhsm_x_dv_r,lhsm_x_dv_j_r,
     >                         imax, jmax, kmax, imaxp, imaxp,
     >                         blockDimSize,iBlocks,jBlocks)
        istat = cudaEventRecord(dynamicEvents(16),0)

        istat = cudaEventRecord(dynamicEvents(17),0)
       call krnl_util_swap4D_itoj_c<<<dimGrid,dimBlock,sharedSize>>>(
     >                         rhs_dv_r,rhs_dv_j_r,
     >                         imax, jmax, kmax, imaxp, jmaxp,
     >                         blockDimSize, iBlocks,jBlocks)   
        istat = cudaEventRecord(dynamicEvents(18),0)

c---------------------------------------------------------------------
c                          FORWARD ELIMINATION  
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      perform the Thomas algorithm; first, FORWARD ELIMINATION     
c---------------------------------------------------------------------


       jSize = ny2
       iSize = grid_points(1)-3+1
       kSize = nz2

       blockSize = fixedBlockSize
       subSections = jSize/blockSize
       if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1
      
       dimGrid = dim3(kSize*subSections,1,1)
       dimBlock = dim3(blockSize,1,1)

       istat = cudaEventRecord(dynamicEvents(19),0)
       
       call  krnl_x_solve_thomas<<<dimGrid,dimBlock>>>(
     >                    lhs_x_dv_j_r, rhs_dv_j_r,
     >                     imax, jmax, kmax, imaxp, jmaxp, iSize,
     >                    jSize, subSections)


       istat = cudaEventRecord(dynamicEvents(20),0)





c---------------------------------------------------------------------
c      The last two rows in this grid block are a bit different, 
c      since they do not have two more rows available for the
c      elimination of off-diagonal entries
c---------------------------------------------------------------------


        jSize = ny2
        kSize = nz2
        istart = grid_points(1)-2

        blockSize = fixedBlockSize
        subSections = jSize/blockSize
        if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1

        dimGrid = dim3(kSize*subSections,1,1)
        dimBlock = dim3(blockSize,1,1)

        istat = cudaEventRecord(dynamicEvents(21),0)
        call  krnl_x_solve_pthomas1<<<dimGrid,dimBlock>>>(
     >                      rhs_dv_j_r,lhs_x_dv_j_r, 
     >                      imax, jmax, kmax, imaxp, jmaxp, 
     >                      jSize, istart, subSections)
        istat = cudaEventRecord(dynamicEvents(22),0)




c---------------------------------------------------------------------
c      do the u+c and the u-c factors                 
c---------------------------------------------------------------------


        jSize = ny2
        iSize = grid_points(1)-3+1
        kSize = nz2

        blockSize = fixedBlockSize
        subSections = jSize/blockSize
        if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1

        dimGrid = dim3(kSize*subSections,1,1)
        dimBlock = dim3(blockSize,1,1)

        istat = cudaEventRecord(dynamicEvents(23),0)
        call  krnl_x_solve_pthomas2<<<dimGrid,dimBlock>>>(
     >                    lhsp_x_dv_j_r,lhsm_x_dv_j_r, rhs_dv_j_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, iSize,
     >                    jSize, subSections)
        istat = cudaEventRecord(dynamicEvents(24),0)


c---------------------------------------------------------------------
c         And again the last two rows separately
c---------------------------------------------------------------------
 

         jSize = ny2
         kSize = nz2
         istart = grid_points(1)-2

        blockSize = fixedBlockSize
        subSections = jSize/blockSize
        if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1

        dimGrid = dim3(kSize*subSections,1,1)
        dimBlock = dim3(blockSize,1,1)


         istat = cudaEventRecord(dynamicEvents(25),0)
         call  krnl_xsolve_pthomas3<<<dimGrid,dimBlock>>>(
     >                      rhs_dv_j_r,lhsm_x_dv_j_r, lhsp_x_dv_j_r,
     >                      imax, jmax, kmax, imaxp, jmaxp, 
     >                      jSize, istart,subSections)
         istat = cudaEventRecord(dynamicEvents(26),0) 


c---------------------------------------------------------------------
c                         BACKSUBSTITUTION 
c---------------------------------------------------------------------


        jSize = ny2
        kSize = nz2
        istart = grid_points(1)-2
       
        blockSize = fixedBlockSize
        subSections = jSize/blockSize
        if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1

        dimGrid = dim3(kSize*subSections,1,1)
        dimBlock = dim3(blockSize,1,1)

        istat = cudaEventRecord(dynamicEvents(27),0)
        call  krnl_xsolve_bfill1<<<dimGrid,dimBlock>>>(
     >        rhs_dv_j_r,lhs_x_dv_j_r, lhsm_x_dv_j_r, lhsp_x_dv_j_r,
     >                      imax, jmax, kmax, imaxp, jmaxp, 
     >                      jSize, istart, subSections)
        istat = cudaEventRecord(dynamicEvents(28),0)


c---------------------------------------------------------------------
c      The first three factors
c---------------------------------------------------------------------

        jSize = ny2
        iSize = grid_points(1)-3+1
        kSize = nz2

        blockSize = fixedBlockSize
        subSections = jSize/blockSize
        if(mod(jSize,blockSize) .ne. 0)subSections=subSections+1

        dimGrid = dim3(kSize*subSections,1,1)
        dimBlock = dim3(blockSize,1,1)

        istat = cudaEventRecord(dynamicEvents(29),0)
        call  krnl_x_solve_bfill2<<<dimGrid,dimBlock>>>(
     >           rhs_dv_j_r,lhs_x_dv_j_r,lhsp_x_dv_j_r,lhsm_x_dv_j_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, iSize,
     >                     jSize, subSections)
        istat = cudaEventRecord(dynamicEvents(30),0)
       



c---------------------------------------------------------------------
c         swap 'j' and 'i' dimensions from rhs_dv_j_r to rhs_dv_r
c--------------------------------------------------------------------

        iSize = IMAXP+1
        jSize = JMAXP+1
        kSize = KMAX

        iBlocks = iSize/blockDimSize
        if(mod(iSize,blockDimSize) .ne. 0)iBlocks=iBlocks+1
        jBlocks = jSize/blockDimSize
        if(mod(jSize,blockDimSize) .ne. 0)jBlocks=jBlocks+1

        dimGrid = dim3(kSize,iBlocks*jBlocks,1)
        dimBlock = dim3(blockDimSize,blockDimSize,1)
        sharedSize = blockDimSize * blockDimSize * sizeOfDouble


	istat = cudaEventRecord(dynamicEvents(31),0)
        call krnl_util_swap4D_jtoi_c<<<dimGrid,dimBlock,sharedSize>>>(
     >                         rhs_dv_j_r,rhs_dv_r,
     >                         imax, jmax, kmax, imaxp, jmaxp,
     >                         blockDimSize, iBlocks, jblocks)  
	istat = cudaEventRecord(dynamicEvents(32),0)

        if (timeron) call timer_stop(t_xsolve)



c---------------------------------------------------------------------
c      Do the block-diagonal inversion          
c---------------------------------------------------------------------

        if (timeron) call timer_start(t_ninvr)       
        istat = cudaEventRecord(dynamicEvents(33),0)
        call ninvr
        istat = cudaEventRecord(dynamicEvents(34),0)
        if (timeron) call timer_stop(t_ninvr)



        if(syncForTimers .gt. 0) then
           istat= cudaThreadSynchronize()
           tink = 0
           do i=1,31,2
              istat = cudaEventElapsedTime(tink, dynamicEvents(i),
     >                                dynamicEvents(i+1))
              call timer_inc(t_kernel,dble(tink/1000.0)) 
              call timer_inc(t_xsolve,dble(tink/1000.0)) 
              !print *, 'x-> ', tink/1000.0
           end do
           istat = cudaEventElapsedTime(tink, dynamicEvents(33),
     >                                dynamicEvents(34))
           call timer_inc(t_kernel,dble(tink/1000.0))
           call timer_inc(t_ninvr,dble(tink/1000.0))
           !print *, 'x-> ', tink/1000.0

        endif



     
       else !BEGIN NON GPU MODE






       if (timeron) call timer_start(t_xsolve)


        if(doGPU .eq. 2) then
          if (timeron) call timer_start(t_gcomm)
             istat = cudaMemcpy(rhs, rhs_dv,
     >                         5*(IMAXP+1)*(JMAXP+1)*KMAX)
          if (timeron) call timer_stop(t_gcomm)
        endif

       do  k = 1, nz2

          call lhsinit(nx2+1, ny2)

c---------------------------------------------------------------------
c Computes the left hand side for the three x-factors  
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      first fill the lhs for the u-eigenvalue                   
c---------------------------------------------------------------------
          do  j = 1, ny2
             do  i = 0, grid_points(1)-1
                ru1 = c3c4*rho_i(i,j,k)
                cv(i) = us(i,j,k)
                rhon(i) = dmax1(dx2+con43*ru1, 
     >                          dx5+c1c5*ru1,
     >                          dxmax+ru1,
     >                          dx1)
             end do

             do  i = 1, nx2
                lhs(1,i,j) =   0.0d0
                lhs(2,i,j) = - dttx2 * cv(i-1) - dttx1 * rhon(i-1)
                lhs(3,i,j) =   1.0d0 + c2dttx1 * rhon(i)
                lhs(4,i,j) =   dttx2 * cv(i+1) - dttx1 * rhon(i+1)
                lhs(5,i,j) =   0.0d0
             end do
          end do

c---------------------------------------------------------------------
c      add fourth order dissipation                             
c---------------------------------------------------------------------

          do  j = 1, ny2
             i = 1
             lhs(3,i,j) = lhs(3,i,j) + comz5
             lhs(4,i,j) = lhs(4,i,j) - comz4
             lhs(5,i,j) = lhs(5,i,j) + comz1
  
             lhs(2,i+1,j) = lhs(2,i+1,j) - comz4
             lhs(3,i+1,j) = lhs(3,i+1,j) + comz6
             lhs(4,i+1,j) = lhs(4,i+1,j) - comz4
             lhs(5,i+1,j) = lhs(5,i+1,j) + comz1
          end do

          do  j = 1, ny2
             do   i=3, grid_points(1)-4
                lhs(1,i,j) = lhs(1,i,j) + comz1
                lhs(2,i,j) = lhs(2,i,j) - comz4
                lhs(3,i,j) = lhs(3,i,j) + comz6
                lhs(4,i,j) = lhs(4,i,j) - comz4
                lhs(5,i,j) = lhs(5,i,j) + comz1
             end do
          end do


          do  j = 1, ny2
             i = grid_points(1)-3
             lhs(1,i,j) = lhs(1,i,j) + comz1
             lhs(2,i,j) = lhs(2,i,j) - comz4
             lhs(3,i,j) = lhs(3,i,j) + comz6
             lhs(4,i,j) = lhs(4,i,j) - comz4

             lhs(1,i+1,j) = lhs(1,i+1,j) + comz1
             lhs(2,i+1,j) = lhs(2,i+1,j) - comz4
             lhs(3,i+1,j) = lhs(3,i+1,j) + comz5
          end do

c---------------------------------------------------------------------
c      subsequently, fill the other factors (u+c), (u-c) by adding to 
c      the first  
c---------------------------------------------------------------------
          do  j = 1, ny2
             do   i = 1, nx2
                lhsp(1,i,j) = lhs(1,i,j)
                lhsp(2,i,j) = lhs(2,i,j) - 
     >                            dttx2 * speed(i-1,j,k)
                lhsp(3,i,j) = lhs(3,i,j)
                lhsp(4,i,j) = lhs(4,i,j) + 
     >                            dttx2 * speed(i+1,j,k)
                lhsp(5,i,j) = lhs(5,i,j)
                lhsm(1,i,j) = lhs(1,i,j)
                lhsm(2,i,j) = lhs(2,i,j) + 
     >                            dttx2 * speed(i-1,j,k)
                lhsm(3,i,j) = lhs(3,i,j)
                lhsm(4,i,j) = lhs(4,i,j) - 
     >                            dttx2 * speed(i+1,j,k)
                lhsm(5,i,j) = lhs(5,i,j)
             end do
          end do

c---------------------------------------------------------------------
c                          FORWARD ELIMINATION  
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      perform the Thomas algorithm; first, FORWARD ELIMINATION     
c---------------------------------------------------------------------

          do  j = 1, ny2
             do    i = 0, grid_points(1)-3
                i1 = i  + 1
                i2 = i  + 2
                fac1      = 1.d0/lhs(3,i,j)
                lhs(4,i,j)  = fac1*lhs(4,i,j)
                lhs(5,i,j)  = fac1*lhs(5,i,j)
                do    m = 1, 3
                   rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
                end do
                lhs(3,i1,j) = lhs(3,i1,j) -
     >                         lhs(2,i1,j)*lhs(4,i,j)
                lhs(4,i1,j) = lhs(4,i1,j) -
     >                         lhs(2,i1,j)*lhs(5,i,j)
                do    m = 1, 3
                   rhs(m,i1,j,k) = rhs(m,i1,j,k) -
     >                         lhs(2,i1,j)*rhs(m,i,j,k)
                end do
                lhs(2,i2,j) = lhs(2,i2,j) -
     >                         lhs(1,i2,j)*lhs(4,i,j)
                lhs(3,i2,j) = lhs(3,i2,j) -
     >                         lhs(1,i2,j)*lhs(5,i,j)
                do    m = 1, 3
                   rhs(m,i2,j,k) = rhs(m,i2,j,k) -
     >                         lhs(1,i2,j)*rhs(m,i,j,k)
                end do
             end do
          end do

c---------------------------------------------------------------------
c      The last two rows in this grid block are a bit different, 
c      since they do not have two more rows available for the
c      elimination of off-diagonal entries
c---------------------------------------------------------------------

          do  j = 1, ny2
             i  = grid_points(1)-2
             i1 = grid_points(1)-1
             fac1      = 1.d0/lhs(3,i,j)
             lhs(4,i,j)  = fac1*lhs(4,i,j)
             lhs(5,i,j)  = fac1*lhs(5,i,j)
             do    m = 1, 3
                rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
             end do
             lhs(3,i1,j) = lhs(3,i1,j) -
     >                      lhs(2,i1,j)*lhs(4,i,j)
             lhs(4,i1,j) = lhs(4,i1,j) -
     >                      lhs(2,i1,j)*lhs(5,i,j)
             do    m = 1, 3
                rhs(m,i1,j,k) = rhs(m,i1,j,k) -
     >                      lhs(2,i1,j)*rhs(m,i,j,k)
             end do
c---------------------------------------------------------------------
c            scale the last row immediately 
c---------------------------------------------------------------------
             fac2             = 1.d0/lhs(3,i1,j)
             do    m = 1, 3
                rhs(m,i1,j,k) = fac2*rhs(m,i1,j,k)
             end do
          end do

c---------------------------------------------------------------------
c      do the u+c and the u-c factors                 
c---------------------------------------------------------------------

          do  j = 1, ny2
             do    i = 0, grid_points(1)-3
                i1 = i  + 1
                i2 = i  + 2
                m = 4
                fac1       = 1.d0/lhsp(3,i,j)
                lhsp(4,i,j)  = fac1*lhsp(4,i,j)
                lhsp(5,i,j)  = fac1*lhsp(5,i,j)
                rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
                lhsp(3,i1,j) = lhsp(3,i1,j) -
     >                        lhsp(2,i1,j)*lhsp(4,i,j)
                lhsp(4,i1,j) = lhsp(4,i1,j) -
     >                        lhsp(2,i1,j)*lhsp(5,i,j)
                rhs(m,i1,j,k) = rhs(m,i1,j,k) -
     >                        lhsp(2,i1,j)*rhs(m,i,j,k)
                lhsp(2,i2,j) = lhsp(2,i2,j) -
     >                        lhsp(1,i2,j)*lhsp(4,i,j)
                lhsp(3,i2,j) = lhsp(3,i2,j) -
     >                        lhsp(1,i2,j)*lhsp(5,i,j)
                rhs(m,i2,j,k) = rhs(m,i2,j,k) -
     >                        lhsp(1,i2,j)*rhs(m,i,j,k)
                m = 5
                fac1       = 1.d0/lhsm(3,i,j)
                lhsm(4,i,j)  = fac1*lhsm(4,i,j)
                lhsm(5,i,j)  = fac1*lhsm(5,i,j)
                rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
                lhsm(3,i1,j) = lhsm(3,i1,j) -
     >                        lhsm(2,i1,j)*lhsm(4,i,j)
                lhsm(4,i1,j) = lhsm(4,i1,j) -
     >                        lhsm(2,i1,j)*lhsm(5,i,j)
                rhs(m,i1,j,k) = rhs(m,i1,j,k) -
     >                        lhsm(2,i1,j)*rhs(m,i,j,k)
                lhsm(2,i2,j) = lhsm(2,i2,j) -
     >                        lhsm(1,i2,j)*lhsm(4,i,j)
                lhsm(3,i2,j) = lhsm(3,i2,j) -
     >                        lhsm(1,i2,j)*lhsm(5,i,j)
                rhs(m,i2,j,k) = rhs(m,i2,j,k) -
     >                        lhsm(1,i2,j)*rhs(m,i,j,k)
             end do
          end do

c---------------------------------------------------------------------
c         And again the last two rows separately
c---------------------------------------------------------------------
          do  j = 1, ny2
             i  = grid_points(1)-2
             i1 = grid_points(1)-1
             m = 4
             fac1       = 1.d0/lhsp(3,i,j)
             lhsp(4,i,j)  = fac1*lhsp(4,i,j)
             lhsp(5,i,j)  = fac1*lhsp(5,i,j)
             rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
             lhsp(3,i1,j) = lhsp(3,i1,j) -
     >                      lhsp(2,i1,j)*lhsp(4,i,j)
             lhsp(4,i1,j) = lhsp(4,i1,j) -
     >                      lhsp(2,i1,j)*lhsp(5,i,j)
             rhs(m,i1,j,k) = rhs(m,i1,j,k) -
     >                      lhsp(2,i1,j)*rhs(m,i,j,k)
             m = 5
             fac1       = 1.d0/lhsm(3,i,j)
             lhsm(4,i,j)  = fac1*lhsm(4,i,j)
             lhsm(5,i,j)  = fac1*lhsm(5,i,j)
             rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
             lhsm(3,i1,j) = lhsm(3,i1,j) -
     >                      lhsm(2,i1,j)*lhsm(4,i,j)
             lhsm(4,i1,j) = lhsm(4,i1,j) -
     >                      lhsm(2,i1,j)*lhsm(5,i,j)
             rhs(m,i1,j,k) = rhs(m,i1,j,k) -
     >                      lhsm(2,i1,j)*rhs(m,i,j,k)
c---------------------------------------------------------------------
c               Scale the last row immediately
c---------------------------------------------------------------------
             rhs(4,i1,j,k) = rhs(4,i1,j,k)/lhsp(3,i1,j)
             rhs(5,i1,j,k) = rhs(5,i1,j,k)/lhsm(3,i1,j)
          end do


c---------------------------------------------------------------------
c                         BACKSUBSTITUTION 
c---------------------------------------------------------------------


          do  j = 1, ny2
             i  = grid_points(1)-2
             i1 = grid_points(1)-1
             do   m = 1, 3
                rhs(m,i,j,k) = rhs(m,i,j,k) -
     >                             lhs(4,i,j)*rhs(m,i1,j,k)
             end do

             rhs(4,i,j,k) = rhs(4,i,j,k) -
     >                          lhsp(4,i,j)*rhs(4,i1,j,k)
             rhs(5,i,j,k) = rhs(5,i,j,k) -
     >                          lhsm(4,i,j)*rhs(5,i1,j,k)
          end do

c---------------------------------------------------------------------
c      The first three factors
c---------------------------------------------------------------------
          do  j = 1, ny2
             do    i = grid_points(1)-3, 0, -1
                i1 = i  + 1
                i2 = i  + 2
                do   m = 1, 3
                   rhs(m,i,j,k) = rhs(m,i,j,k) - 
     >                          lhs(4,i,j)*rhs(m,i1,j,k) -
     >                          lhs(5,i,j)*rhs(m,i2,j,k)
                end do

c---------------------------------------------------------------------
c      And the remaining two
c---------------------------------------------------------------------
                rhs(4,i,j,k) = rhs(4,i,j,k) - 
     >                          lhsp(4,i,j)*rhs(4,i1,j,k) -
     >                          lhsp(5,i,j)*rhs(4,i2,j,k)
                rhs(5,i,j,k) = rhs(5,i,j,k) - 
     >                          lhsm(4,i,j)*rhs(5,i1,j,k) -
     >                          lhsm(5,i,j)*rhs(5,i2,j,k)
             end do
          end do
       end do
       if (timeron) call timer_stop(t_xsolve)

c---------------------------------------------------------------------
c      Do the block-diagonal inversion          
c---------------------------------------------------------------------
       call ninvr


       endif
       return
       end
    






