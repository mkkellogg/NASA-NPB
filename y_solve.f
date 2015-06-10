
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine y_solve
       use glob
       use cudafor
       use lhs_kernels
       use y_solve_kernels
       use util_kernels
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c this function performs the solution of the approximate factorization
c step in the y-direction for all five matrix components
c simultaneously. The Thomas algorithm is employed to solve the
c systems for the y-lines. Boundary conditions are non-periodic
c---------------------------------------------------------------------

       include 'header.h'

       integer i, j, k, j1, j2, m, kSize,jSize, sharedSize, ni, nj,
     >         blockSize, istat, upperj1, upperj2,  iSize,
     >         istart, jstart, subSections, subISize, blockDimSize,
     >         iBlocks, jBlocks, blockWidth,blockHeight, subSectionsW,
     >         subSectionsH
       real tink
       double precision ru1, fac1, fac2


c---------------------------------------------------------------------
c---------------------------------------------------------------------

       if(doGPU .gt. 0) then !BEGIN GPU SWITCH

       if (timeron) call timer_start(t_ysolve)


c---------------------------------------------------------------------
c Computes the left hand side for the three y-factors   
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c         swap 'i' and 'j' dimensions from rho_i_dv to rho_i_dv_j_r
c         and from vs_dv to vs_dv_j_r
c---------------------------------------------------------------------

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

        istat = cudaEventRecord(dynamicEvents(1),0)
        call krnl_util_swap3D_itoj_c<<<dimGrid,dimBlock,sharedSize>>>(
     >                         rho_i_dv,rho_i_dv_j_r,
     >                         imax, jmax, kmax, imaxp, imaxp,
     >                         blockDimSize, iBlocks, jBlocks)
        istat = cudaEventRecord(dynamicEvents(2),0)

        istat = cudaEventRecord(dynamicEvents(3),0)
        call krnl_util_swap3D_itoj_c<<<dimGrid,dimBlock,sharedSize>>>(
     >                         vs_dv,vs_dv_j_r,
     >                         imax, jmax, kmax, imaxp, imaxp,
     >                         blockDimSize, iBlocks, jBlocks)
        istat = cudaEventRecord(dynamicEvents(4),0)

c---------------------------------------------------------------------
c      first fill the lhs for the u-eigenvalue         
c---------------------------------------------------------------------

       upperj1 = grid_points(2)-1
       upperj2 = grid_points(2)-2

       iSize = grid_points(1)-2
       jSize = grid_points(2)-2
       kSize = grid_points(3)-2
       
       blockSize = fixedBlockSize
       subSections = iSize/fixedBlockSize
       if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1
       blockSize = (blockSize+2)

       dimGrid = dim3(kSize,iSize*subSections,1)
       dimBlock = dim3(blockSize,1,1)
       sharedSize = blockSize * sizeOfDouble * 2


       istat = cudaEventRecord(dynamicEvents(5),0)
       call  krnl_y_solve_init<<<dimGrid,dimBlock,sharedSize>>>(
     >      lhs_x_dv_j_r, rho_i_dv_j_r, vs_dv_j_r,  imax, jmax, kmax, 
     >                       imaxp, jmaxp,  dtty1, dtty2, c2dtty1,
     >                  con43, c1c5, c3c4, dy1, dy3, dy5,dymax,
     >                     upperj1,upperj2, subSections)
       istat = cudaEventRecord(dynamicEvents(6),0)

   

c---------------------------------------------------------------------
c         swap 'j' and 'i' dimensions from lhs_x_dv_j_r to lhs_x_dv_r
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

        istat = cudaEventRecord(dynamicEvents(7),0)
        call krnl_util_swap4D_jtoi_c<<<dimGrid,dimBlock,sharedSize>>>(
     >                         lhs_x_dv_j_r,lhs_x_dv_r,
     >                         imax, jmax, kmax, imaxp, jmaxp,

     >                         blockDimSize, iBlocks, jblocks)
        istat = cudaEventRecord(dynamicEvents(8),0)
 


c---------------------------------------------------------------------
c      add fourth order dissipation                             
c---------------------------------------------------------------------

       jstart = -1
       iSize = grid_points(1)-2
       jSize = 1
       kSize = grid_points(3)-2


       blockSize = fixedBlockSize
       subSections = iSize/fixedBlockSize
       if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1

       dimGrid = dim3(kSize*subSections,1,1)
       dimBlock = dim3(blockSize,1,1)


       istat = cudaEventRecord(dynamicEvents(9),0)
       call  krnl_y_solve_dis1<<<dimGrid,dimBlock>>>(
     >           lhs_x_dv_r, imax, jmax, kmax, imaxp, jmaxp,  
     >                    comz1,comz4,comz5,comz6, jstart,
     >                    iSize, subSections)

        istat = cudaEventRecord(dynamicEvents(10),0)




        iSize = grid_points(1)-2
        jSize = (grid_points(2)-4) -3 + 1
        kSize = grid_points(3)-2


       blockSize = fixedBlockSize
       subSections= (iSize/fixedBlockSize)
       if(mod(iSize,fixedBlockSize) .ne. 0)subSections = subSections+1

       dimGrid = dim3(kSize,jSize*subSections,1)
       dimBlock = dim3(blockSize,1,1)


        istat = cudaEventRecord(dynamicEvents(11),0)
        call  krnl_y_solve_dis2<<<dimGrid,dimBlock>>>(
     >           lhs_x_dv_r, imax, jmax, kmax, imaxp, jmaxp,  
     >                    comz1,comz4,comz5,comz6,
     >                    iSize, subSections)
        istat = cudaEventRecord(dynamicEvents(12),0)






        jstart = grid_points(2)-3
        iSize = grid_points(1)-2
        jSize = 1
        kSize = grid_points(3)-2

        blockSize = fixedBlockSize
        subSections = iSize/fixedBlockSize
        if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1

        dimGrid = dim3(kSize*subSections,1,1)
        dimBlock = dim3(blockSize,1,1)

        istat = cudaEventRecord(dynamicEvents(13),0)
        call  krnl_y_solve_dis1<<<dimGrid,dimBlock>>>(
     >           lhs_x_dv_r, imax, jmax, kmax, imaxp, jmaxp,  
     >                    comz1,comz4,comz5,comz6, jstart,
     >                    iSize,subSections)
        istat = cudaEventRecord(dynamicEvents(14),0)






c---------------------------------------------------------------------
c      subsequently, do the other two factors                    
c---------------------------------------------------------------------


        iSize = grid_points(1)-2
        jSize = grid_points(2)-2
        kSize = grid_points(3)-2

        blockSize = fixedBlockSize
        subSections= (iSize/fixedBlockSize)
        if(mod(iSize,fixedBlockSize) .ne. 0)subSections = subSections+1

        dimGrid = dim3(kSize,jSize*subSections,1)
        dimBlock = dim3(blockSize,1,1)

        istat = cudaEventRecord(dynamicEvents(15),0)
        call  krnl_y_solve_postdis<<<dimGrid,dimBlock>>>(
     >          lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r, speed_dv,
     >                    imax, jmax, kmax, imaxp, jmaxp, dtty2,  
     >                    iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(16),0)



c--------------------------------------------------------------------
c                          FORWARD ELIMINATION
c---------------------------------------------------------------------

         jSize =  grid_points(2)-3 + 1
         iSize = grid_points(1)-2
         kSize = grid_points(3)-2

         blockSize = fixedBlockSize
         subSections = iSize/fixedBlockSize
         if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1

         dimGrid = dim3(kSize*subSections,1,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(17),0)
         call  krnl_y_solve_forward1<<<dimGrid,dimBlock>>>(
     >                    lhs_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jSize,
     >                    iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(18),0)




c---------------------------------------------------------------------
c      The last two rows in this grid block are a bit different,
c      since they do not have two more rows available for the
c      elimination of off-diagonal entries
c---------------------------------------------------------------------

         jstart= grid_points(2)-2
         iSize = grid_points(1)-2
         kSize = grid_points(3)-2

         blockSize = fixedBlockSize
         subSections = iSize/fixedBlockSize
         if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1

         dimGrid = dim3(kSize*subSections,1,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(19),0)
         call  krnl_y_solve_forward2<<<dimGrid,dimBlock>>>(
     >                    lhs_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jstart,
     >                    iSize,  subSections)
         istat = cudaEventRecord(dynamicEvents(20),0)






c---------------------------------------------------------------------
c      do the u+c and the u-c factors
c---------------------------------------------------------------------

         jSize = grid_points(2)-3 + 1
         iSize = grid_points(1)-2
         kSize = grid_points(3)-2

         blockSize = fixedBlockSize
         subSections = iSize/fixedBlockSize
         if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1

         dimGrid = dim3(kSize*subSections,1,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(21),0)
         call  krnl_y_solve_forward3<<<dimGrid,dimBlock>>>(
     >                    lhsp_x_dv_r,lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jSize,
     >                    iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(22),0)



         !
         ! Experimenting with 2-Dimensional blocks
         !

         !blockWidth = fixedBlockWidth
         !subSectionsW = iSize/fixedBlockWidth
         !if(mod(iSize,fixedBlockWidth) .ne. 0)
        !>                           subSectionsW=subSectionsW+1

         !blockHeight = fixedBlockHeight
         !subSectionsH = kSize/fixedBlockHeight
         !if(mod(kSize,fixedBlockHeight) .ne. 0)
        !>                           subSectionsH=subSectionsH+1

         !dimGrid = dim3(subSectionsW,subSectionsH,1)
         !dimBlock = dim3(fixedBlockWidth,fixedBlockHeight,1)
         !sharedSize = 32

         !istat = cudaEventRecord(dynamicEvents(21),0)
         !call  krnl_y_solve_forward3B<<<dimGrid,dimBlock,sharedSize>>>(
        !>                    lhsp_x_dv_r,lhsm_x_dv_r, rhs_dv_r,
        !>                    imax, jmax, kmax, imaxp, jmaxp,kSize, jSize,
        !>                    iSize, subSectionsW,subSectionsH)
         !istat = cudaEventRecord(dynamicEvents(22),0)




c---------------------------------------------------------------------
c         And again the last two rows separately
c---------------------------------------------------------------------

         jstart= grid_points(2)-2
         iSize = grid_points(1)-2
         kSize = grid_points(3)-2


         blockSize = fixedBlockSize
         subSections = iSize/fixedBlockSize
         if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1

         dimGrid = dim3(kSize*subSections,1,1)
         dimBlock = dim3(blockSize,1,1)
         sharedSize = 32


         istat = cudaEventRecord(dynamicEvents(23),0)
         call  krnl_y_solve_forward4<<<dimGrid,dimBlock,sharedSize>>>(
     >                    lhsp_x_dv_r, lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jstart,
     >                    iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(24),0)



c---------------------------------------------------------------------
c                         BACKSUBSTITUTION
c---------------------------------------------------------------------

         jstart= grid_points(2)-2
         iSize = grid_points(1)-2
         kSize = grid_points(3)-2


         blockSize = fixedBlockSize
         subSections = iSize/fixedBlockSize
         if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1

         dimGrid = dim3(kSize*subSections,1,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(25),0)
         call  krnl_y_solve_back1<<<dimGrid,dimBlock>>>(
     >            lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jstart,
     >                    iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(26),0)



c---------------------------------------------------------------------
c      The first three factors
c---------------------------------------------------------------------

         jSize = grid_points(2)-3 + 1
         iSize = grid_points(1)-2
         kSize = grid_points(3)-2

         blockSize = fixedBlockSize
         subSections = iSize/fixedBlockSize
         if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1

         dimGrid = dim3(kSize*subSections,1,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(27),0)
         call  krnl_y_solve_back2<<<dimGrid,dimBlock>>>(
     >             lhs_x_dv_r,lhsp_x_dv_r,lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jSize,
     >                    iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(28),0)


         if (timeron) call timer_stop(t_ysolve)


      

         if (timeron) call timer_start(t_pinvr) 
         istat = cudaEventRecord(dynamicEvents(29),0)
         call pinvr
         istat = cudaEventRecord(dynamicEvents(30),0)
         if (timeron) call timer_stop(t_pinvr)

          

         if(syncForTimers .gt. 0) then
            istat= cudaThreadSynchronize()
            tink = 0
            do i=1,27,2
               istat =cudaEventElapsedTime(tink, dynamicEvents(i),
     >                                dynamicEvents(i+1))
               call timer_inc(t_kernel,dble(tink/1000.0))
               call timer_inc(t_ysolve,dble(tink/1000.0))
               !print *, 'y-> ', tink/1000.0
            end do
            istat =cudaEventElapsedTime(tink, dynamicEvents(29),
     >                                dynamicEvents(30))
            call timer_inc(t_pinvr,dble(tink/1000.0))
            call timer_inc(t_kernel,dble(tink/1000.0))
            !print *, 'y-> ', tink/1000.0

         endif





       else ! BEGIN NON GPU MODE






      
       if (timeron) call timer_start(t_ysolve)
       do  k = 1, grid_points(3)-2

          call lhsinitj(ny2+1, nx2)

c---------------------------------------------------------------------
c Computes the left hand side for the three y-factors   
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      first fill the lhs for the u-eigenvalue         
c---------------------------------------------------------------------

          do  i = 1, grid_points(1)-2
             do  j = 0, grid_points(2)-1
                ru1 = c3c4*rho_i(i,j,k)
                cv(j) = vs(i,j,k)
                rhoq(j) = dmax1( dy3 + con43 * ru1,
     >                           dy5 + c1c5*ru1,
     >                           dymax + ru1,
     >                           dy1)
             end do
            
             do  j = 1, grid_points(2)-2
                lhs(1,i,j) =  0.0d0
                lhs(2,i,j) = -dtty2 * cv(j-1) - dtty1 * rhoq(j-1)
                lhs(3,i,j) =  1.0 + c2dtty1 * rhoq(j)
                lhs(4,i,j) =  dtty2 * cv(j+1) - dtty1 * rhoq(j+1)
                lhs(5,i,j) =  0.0d0
             end do
          end do

c---------------------------------------------------------------------
c      add fourth order dissipation                             
c---------------------------------------------------------------------

          do  i = 1, grid_points(1)-2
             j = 1
             lhs(3,i,j) = lhs(3,i,j) + comz5
             lhs(4,i,j) = lhs(4,i,j) - comz4
             lhs(5,i,j) = lhs(5,i,j) + comz1
       
             lhs(2,i,j+1) = lhs(2,i,j+1) - comz4
             lhs(3,i,j+1) = lhs(3,i,j+1) + comz6
             lhs(4,i,j+1) = lhs(4,i,j+1) - comz4
             lhs(5,i,j+1) = lhs(5,i,j+1) + comz1
          end do

          do   j=3, grid_points(2)-4
             do  i = 1, grid_points(1)-2

                lhs(1,i,j) = lhs(1,i,j) + comz1
                lhs(2,i,j) = lhs(2,i,j) - comz4
                lhs(3,i,j) = lhs(3,i,j) + comz6
                lhs(4,i,j) = lhs(4,i,j) - comz4
                lhs(5,i,j) = lhs(5,i,j) + comz1
             end do
          end do

          do  i = 1, grid_points(1)-2
             j = grid_points(2)-3
             lhs(1,i,j) = lhs(1,i,j) + comz1
             lhs(2,i,j) = lhs(2,i,j) - comz4
             lhs(3,i,j) = lhs(3,i,j) + comz6
             lhs(4,i,j) = lhs(4,i,j) - comz4

             lhs(1,i,j+1) = lhs(1,i,j+1) + comz1
             lhs(2,i,j+1) = lhs(2,i,j+1) - comz4
             lhs(3,i,j+1) = lhs(3,i,j+1) + comz5
          end do

c---------------------------------------------------------------------
c      subsequently, do the other two factors                    
c---------------------------------------------------------------------
          do    j = 1, grid_points(2)-2
             do  i = 1, grid_points(1)-2
                lhsp(1,i,j) = lhs(1,i,j)
                lhsp(2,i,j) = lhs(2,i,j) - 
     >                            dtty2 * speed(i,j-1,k)
                lhsp(3,i,j) = lhs(3,i,j)
                lhsp(4,i,j) = lhs(4,i,j) + 
     >                            dtty2 * speed(i,j+1,k)
                lhsp(5,i,j) = lhs(5,i,j)
                lhsm(1,i,j) = lhs(1,i,j)
                lhsm(2,i,j) = lhs(2,i,j) + 
     >                            dtty2 * speed(i,j-1,k)
                lhsm(3,i,j) = lhs(3,i,j)
                lhsm(4,i,j) = lhs(4,i,j) - 
     >                            dtty2 * speed(i,j+1,k)
                lhsm(5,i,j) = lhs(5,i,j)
             end do
          end do


c---------------------------------------------------------------------
c                          FORWARD ELIMINATION  
c---------------------------------------------------------------------

          do    j = 0, grid_points(2)-3
             j1 = j  + 1
             j2 = j  + 2
             do  i = 1, grid_points(1)-2
                fac1      = 1.d0/lhs(3,i,j)
                lhs(4,i,j)  = fac1*lhs(4,i,j)
                lhs(5,i,j)  = fac1*lhs(5,i,j)
                do    m = 1, 3
                   rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
                end do
                lhs(3,i,j1) = lhs(3,i,j1) -
     >                         lhs(2,i,j1)*lhs(4,i,j)
                lhs(4,i,j1) = lhs(4,i,j1) -
     >                         lhs(2,i,j1)*lhs(5,i,j)
                do    m = 1, 3
                   rhs(m,i,j1,k) = rhs(m,i,j1,k) -
     >                         lhs(2,i,j1)*rhs(m,i,j,k)
                end do
                lhs(2,i,j2) = lhs(2,i,j2) -
     >                         lhs(1,i,j2)*lhs(4,i,j)
                lhs(3,i,j2) = lhs(3,i,j2) -
     >                         lhs(1,i,j2)*lhs(5,i,j)
                do    m = 1, 3
                   rhs(m,i,j2,k) = rhs(m,i,j2,k) -
     >                         lhs(1,i,j2)*rhs(m,i,j,k)
                end do
             end do
          end do

c---------------------------------------------------------------------
c      The last two rows in this grid block are a bit different, 
c      since they do not have two more rows available for the
c      elimination of off-diagonal entries
c---------------------------------------------------------------------

          j  = grid_points(2)-2
          j1 = grid_points(2)-1
          do  i = 1, grid_points(1)-2
             fac1      = 1.d0/lhs(3,i,j)
             lhs(4,i,j)  = fac1*lhs(4,i,j)
             lhs(5,i,j)  = fac1*lhs(5,i,j)
             do    m = 1, 3
                rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
             end do
             lhs(3,i,j1) = lhs(3,i,j1) -
     >                      lhs(2,i,j1)*lhs(4,i,j)
             lhs(4,i,j1) = lhs(4,i,j1) -
     >                      lhs(2,i,j1)*lhs(5,i,j)
             do    m = 1, 3
                rhs(m,i,j1,k) = rhs(m,i,j1,k) -
     >                      lhs(2,i,j1)*rhs(m,i,j,k)
             end do
c---------------------------------------------------------------------
c            scale the last row immediately 
c---------------------------------------------------------------------
             fac2      = 1.d0/lhs(3,i,j1)
             do    m = 1, 3
                rhs(m,i,j1,k) = fac2*rhs(m,i,j1,k)
             end do
          end do

c---------------------------------------------------------------------
c      do the u+c and the u-c factors                 
c---------------------------------------------------------------------
          do    j = 0, grid_points(2)-3
             j1 = j  + 1
             j2 = j  + 2
             do  i = 1, grid_points(1)-2
                m = 4
                fac1       = 1.d0/lhsp(3,i,j)
                lhsp(4,i,j)  = fac1*lhsp(4,i,j)
                lhsp(5,i,j)  = fac1*lhsp(5,i,j)
                rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
                lhsp(3,i,j1) = lhsp(3,i,j1) -
     >                       lhsp(2,i,j1)*lhsp(4,i,j)
                lhsp(4,i,j1) = lhsp(4,i,j1) -
     >                       lhsp(2,i,j1)*lhsp(5,i,j)
                rhs(m,i,j1,k) = rhs(m,i,j1,k) -
     >                       lhsp(2,i,j1)*rhs(m,i,j,k)
                lhsp(2,i,j2) = lhsp(2,i,j2) -
     >                       lhsp(1,i,j2)*lhsp(4,i,j)
                lhsp(3,i,j2) = lhsp(3,i,j2) -
     >                       lhsp(1,i,j2)*lhsp(5,i,j)
                rhs(m,i,j2,k) = rhs(m,i,j2,k) -
     >                       lhsp(1,i,j2)*rhs(m,i,j,k)
                m = 5
                fac1       = 1.d0/lhsm(3,i,j)
                lhsm(4,i,j)  = fac1*lhsm(4,i,j)
                lhsm(5,i,j)  = fac1*lhsm(5,i,j)
                rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
                lhsm(3,i,j1) = lhsm(3,i,j1) -
     >                       lhsm(2,i,j1)*lhsm(4,i,j)
                lhsm(4,i,j1) = lhsm(4,i,j1) -
     >                       lhsm(2,i,j1)*lhsm(5,i,j)
                rhs(m,i,j1,k) = rhs(m,i,j1,k) -
     >                       lhsm(2,i,j1)*rhs(m,i,j,k)
                lhsm(2,i,j2) = lhsm(2,i,j2) -
     >                       lhsm(1,i,j2)*lhsm(4,i,j)
                lhsm(3,i,j2) = lhsm(3,i,j2) -
     >                       lhsm(1,i,j2)*lhsm(5,i,j)
                rhs(m,i,j2,k) = rhs(m,i,j2,k) -
     >                       lhsm(1,i,j2)*rhs(m,i,j,k)
             end do
          end do

c---------------------------------------------------------------------
c         And again the last two rows separately
c---------------------------------------------------------------------
          j  = grid_points(2)-2
          j1 = grid_points(2)-1
          do  i = 1, grid_points(1)-2
             m = 4
             fac1       = 1.d0/lhsp(3,i,j)
             lhsp(4,i,j)  = fac1*lhsp(4,i,j)
             lhsp(5,i,j)  = fac1*lhsp(5,i,j)
             rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
             lhsp(3,i,j1) = lhsp(3,i,j1) -
     >                    lhsp(2,i,j1)*lhsp(4,i,j)
             lhsp(4,i,j1) = lhsp(4,i,j1) -
     >                    lhsp(2,i,j1)*lhsp(5,i,j)
             rhs(m,i,j1,k)   = rhs(m,i,j1,k) -
     >                    lhsp(2,i,j1)*rhs(m,i,j,k)
             m = 5
             fac1       = 1.d0/lhsm(3,i,j)
             lhsm(4,i,j)  = fac1*lhsm(4,i,j)
             lhsm(5,i,j)  = fac1*lhsm(5,i,j)
             rhs(m,i,j,k) = fac1*rhs(m,i,j,k)
             lhsm(3,i,j1) = lhsm(3,i,j1) -
     >                    lhsm(2,i,j1)*lhsm(4,i,j)
             lhsm(4,i,j1) = lhsm(4,i,j1) -
     >                    lhsm(2,i,j1)*lhsm(5,i,j)
             rhs(m,i,j1,k)   = rhs(m,i,j1,k) -
     >                    lhsm(2,i,j1)*rhs(m,i,j,k)
c---------------------------------------------------------------------
c               Scale the last row immediately 
c---------------------------------------------------------------------
             rhs(4,i,j1,k)   = rhs(4,i,j1,k)/lhsp(3,i,j1)
             rhs(5,i,j1,k)   = rhs(5,i,j1,k)/lhsm(3,i,j1)
          end do


c---------------------------------------------------------------------
c                         BACKSUBSTITUTION 
c---------------------------------------------------------------------

          j  = grid_points(2)-2
          j1 = grid_points(2)-1
          do  i = 1, grid_points(1)-2
             do   m = 1, 3
                rhs(m,i,j,k) = rhs(m,i,j,k) -
     >                           lhs(4,i,j)*rhs(m,i,j1,k)
             end do

             rhs(4,i,j,k) = rhs(4,i,j,k) -
     >                           lhsp(4,i,j)*rhs(4,i,j1,k)
             rhs(5,i,j,k) = rhs(5,i,j,k) -
     >                           lhsm(4,i,j)*rhs(5,i,j1,k)
          end do

c---------------------------------------------------------------------
c      The first three factors
c---------------------------------------------------------------------
          do   j = grid_points(2)-3, 0, -1
             j1 = j  + 1
             j2 = j  + 2
             do  i = 1, grid_points(1)-2
                do   m = 1, 3
                   rhs(m,i,j,k) = rhs(m,i,j,k) - 
     >                          lhs(4,i,j)*rhs(m,i,j1,k) -
     >                          lhs(5,i,j)*rhs(m,i,j2,k)
                end do

c---------------------------------------------------------------------
c      And the remaining two
c---------------------------------------------------------------------
                rhs(4,i,j,k) = rhs(4,i,j,k) - 
     >                          lhsp(4,i,j)*rhs(4,i,j1,k) -
     >                          lhsp(5,i,j)*rhs(4,i,j2,k)
                rhs(5,i,j,k) = rhs(5,i,j,k) - 
     >                          lhsm(4,i,j)*rhs(5,i,j1,k) -
     >                          lhsm(5,i,j)*rhs(5,i,j2,k)
             end do
          end do
       end do
       if (timeron) call timer_stop(t_ysolve)

       call pinvr
        


       endif !END GPU SWITCH

       return
       end
    






