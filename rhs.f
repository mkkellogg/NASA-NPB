c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine compute_rhs(lastRHS)
       use glob
       use rhs_kernels
       use util_kernels
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       include 'header.h'

       integer i, j, k, m, iSize, jSize, kSize, blockSize, numBlocks,
     >         istat, sharedSize, istart, jstart, kstart, subSections,
     >         subISize, subSectionsW, subSectionsH, blockWidth, 
     >         blockHeight
       double precision aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1,
     >                  wijk, wp1, wm1
       real tink 
       integer :: lastRHS
       
       if(doGPU .gt. 0)  then !BEGIN GPU SWITCH


       if (timeron) call timer_start(t_rhs)

c---------------------------------------------------------------------
c         compute the reciprocal of density, and the kinetic energy,
c         and the speed of sound.
c---------------------------------------------------------------------
                  
       iSize = grid_points(1)
       jSize = grid_points(2)
       kSize = grid_points(3)

 
       blockSize = fixedBlockSize
       subSections = iSize/blockSize
       if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
       dimGrid = dim3(kSize,jSize*subSections,1)
       dimBlock = dim3(blockSize,1,1)
      
       istat = cudaEventRecord(dynamicEvents(1),0)
       call krnl_rhs_init<<<dimGrid,dimBlock>>>(u_dv_r,us_dv,
     >                                   vs_dv, ws_dv,qs_dv,
     >                                 rho_i_dv, speed_dv, square_dv,
     >                                        imax,jmax,kmax, 
     >                    c1c2, imaxp, jmaxp, iSize, subSections) 
        istat = cudaEventRecord(dynamicEvents(2),0)  



c---------------------------------------------------------------------
c copy the exact forcing term to the right hand side;  because 
c this forcing term is known, we can store it on the whole grid
c including the boundary                   
c---------------------------------------------------------------------

        if (timeron) call timer_start(t_gcomm)
        istat = cudaMemcpy(rhs_dv_r,forcing_dv_r,
     >                          5*(IMAXP+1)*(JMAXP+1)*KMAX)   

        if (timeron) call timer_stop(t_gcomm)


c---------------------------------------------------------------------
c         compute xi-direction fluxes
c---------------------------------------------------------------------

        if (timeron) call timer_start(t_rhsx)

        iSize = nx2
        jSize = ny2
        kSize = nz2
     
        blockSize = fixedBlockSize     
        subSections = iSize/blockSize
        if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
        dimGrid = dim3(kSize,jSize*subSections,1)
        dimBlock = dim3(blockSize,1,1)


        istat = cudaEventRecord(dynamicEvents(3),0)

       call krnl_rhs_xflux<<<dimGrid,dimBlock>>>(u_dv_r,
     >                  us_dv,  vs_dv, ws_dv,  qs_dv, rho_i_dv, 
     >                         square_dv,rhs_dv_r,imax, jmax, kmax, 
     >                            imaxp, jmaxp,
     >                              dx1tx1, dx2tx1, dx3tx1,
     >                              dx4tx1,dx5tx1, tx2,
     >                               xxcon2, xxcon3, xxcon4, 
     >           xxcon5, con43, c2, c1, iSize, subSections) 


        istat = cudaEventRecord(dynamicEvents(4),0)




c---------------------------------------------------------------------
c      add fourth order xi-direction dissipation               
c--------------------------------------------------------------------- 

        jSize = ny2
        kSize = nz2
        istart = -1

        blockSize = 5      
        dimGrid = dim3(kSize,jSize,1)
        dimBlock = dim3(blockSize,1,1)
        sharedSize = 40 * sizeOfDouble

        istat = cudaEventRecord(dynamicEvents(5),0)
        call krnl_rhs_xdis1<<<dimGrid,dimBlock,sharedSize>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,istart)
        istat = cudaEventRecord(dynamicEvents(6),0)




        iSize = (nx2-2) - 3 + 1
        jSize = ny2
        kSize = nz2

        blockSize = fixedBlockSize
        subSections = iSize/blockSize
        if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
        dimGrid = dim3(kSize,jSize*subSections,1)
        dimBlock = dim3(blockSize,1,1)


        istat = cudaEventRecord(dynamicEvents(7),0)
        call krnl_rhs_xdis2<<<dimGrid,dimBlock>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp, iSize, subSections)
        istat = cudaEventRecord(dynamicEvents(8),0)




        jSize = ny2
        kSize = nz2
        istart = nx2-1

        blockSize = 5
        dimGrid = dim3(kSize,jSize,1)
        dimBlock = dim3(blockSize,1,1)
        sharedSize = 40 * sizeOfDouble


        istat = cudaEventRecord(dynamicEvents(9),0)
        call krnl_rhs_xdis1<<<dimGrid,dimBlock,sharedSize>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp, istart)
        istat = cudaEventRecord(dynamicEvents(10),0)


        if (timeron) call timer_stop(t_rhsx)






        if (timeron) call timer_start(t_rhsy)       

c---------------------------------------------------------------------
c      compute eta-direction fluxes 
c---------------------------------------------------------------------    

        iSize = nx2
        jSize = ny2
        kSize = nz2
     

        !blockSize = fixedBlockSize     
        !subSections = iSize/fixedBlockSize
        !if(mod(iSize,fixedBlockSize) .ne. 0)subSections=subSections+1
        !dimGrid = dim3(kSize,jSize*subSections,1)
        !dimBlock = dim3(blockSize,1,1)


        !istat = cudaEventRecord(dynamicEvents(11),0)
        !call krnl_rhs_yflux<<<dimGrid,dimBlock>>>(
       !>             u_dv_r,us_dv, vs_dv, ws_dv,  qs_dv, rho_i_dv, 
       !>               square_dv,rhs_dv_r,imax, jmax, kmax, imaxp, jmaxp,
       !>                              dy1ty1, dy2ty1, dy3ty1,
       !>                              dy4ty1,dy5ty1, ty2,
       !>                               yycon2, yycon3, yycon4, 
       !>                     yycon5, con43, c2, c1, iSize, subSections) 
        !istat = cudaEventRecord(dynamicEvents(12),0)



         blockWidth = fixedBlockWidth
         subSectionsW = iSize/fixedBlockWidth
         if(mod(iSize,fixedBlockWidth) .ne. 0)
     >                           subSectionsW=subSectionsW+1

         blockHeight = fixedBlockHeight
         subSectionsH = jSize/fixedBlockHeight
         if(mod(jSize,fixedBlockHeight) .ne. 0)
     >                           subSectionsH=subSectionsH+1

         dimGrid = dim3(kSize,subSectionsW * subSectionsH,1)
         dimBlock = dim3(fixedBlockWidth,fixedBlockHeight,1)
         sharedSize = 32

        istat = cudaEventRecord(dynamicEvents(11),0)
        call krnl_rhs_yfluxB<<<dimGrid,dimBlock,sharedSize>>>(
     >             u_dv_r,us_dv, vs_dv, ws_dv,  qs_dv, rho_i_dv,
     >               square_dv,rhs_dv_r,imax, jmax, kmax, imaxp, jmaxp,
     >                              dy1ty1, dy2ty1, dy3ty1,
     >                              dy4ty1,dy5ty1, ty2,
     >                               yycon2, yycon3, yycon4,
     >              yycon5, con43, c2, c1, iSize,jSize, subSectionsW,
     >                    subSectionsH)
        istat = cudaEventRecord(dynamicEvents(12),0)





c---------------------------------------------------------------------
c      add fourth order eta-direction dissipation         
c---------------------------------------------------------------------


         iSize = nx2
         kSize = nz2
         jstart = -1

         blockSize = fixedBlockSize
         subSections = iSize/blockSize
         if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
         dimGrid = dim3(kSize*subSections,1,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(13),0)
         call krnl_rhs_ydis1<<<dimGrid,dimBlock>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,jstart, iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(14),0)





         iSize = nx2
         jSize = (ny2-2) - 3 + 1
         kSize = nz2
           
         blockSize = fixedBlockSize
         subSections = iSize/blockSize
         if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
         dimGrid = dim3(kSize,jSize*subSections,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(15),0)
         call krnl_rhs_ydis2<<<dimGrid,dimBlock>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,iSize,subSections)
         istat = cudaEventRecord(dynamicEvents(16),0)


		  

         iSize = nx2
         kSize = nz2
         jstart = ny2-1

         blockSize = fixedBlockSize
         subSections = iSize/blockSize
         if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
         dimGrid = dim3(kSize*subSections,1,1)
         dimBlock = dim3(blockSize,1,1)
         

         istat = cudaEventRecord(dynamicEvents(17),0)
         call krnl_rhs_ydis1<<<dimGrid,dimBlock>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp, jstart, iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(18),0)




       



       if (timeron) call timer_stop(t_rhsy)


       if (timeron) call timer_start(t_rhsz)
  
c---------------------------------------------------------------------
c      compute zeta-direction fluxes 
c---------------------------------------------------------------------
        iSize = nx2
        jSize = ny2
        kSize = nz2          
     
           
        blockSize = fixedBlockSize     
        subSections = iSize/blockSize
        if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
        dimGrid = dim3(kSize,jSize*subSections,1)
        dimBlock = dim3(blockSize,1,1)

        istat = cudaEventRecord(dynamicEvents(19),0)
        call krnl_rhs_zflux<<<dimGrid,dimBlock>>>(
     >         u_dv_r,us_dv, vs_dv, ws_dv,  qs_dv, rho_i_dv, 
     >                       square_dv,rhs_dv_r,imax, jmax, kmax, 
     >                            imaxp, jmaxp,
     >                              dz1tz1, dz2tz1, dz3tz1,
     >                              dz4tz1,dz5tz1, tz2,
     >                               zzcon2, zzcon3, zzcon4, 
     >           zzcon5, con43, c2, c1, iSize, subSections) 
        istat = cudaEventRecord(dynamicEvents(20),0)



c---------------------------------------------------------------------
c      add fourth order zeta-direction dissipation                
c---------------------------------------------------------------------

        iSize = nx2
        jSize = ny2
        kstart = -1

        blockSize = fixedBlockSize
        subSections = iSize/blockSize
        if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
        dimGrid = dim3(jSize*subSections,1,1)
        dimBlock = dim3(blockSize,1,1)

        istat = cudaEventRecord(dynamicEvents(21),0)
        call krnl_rhs_zdis1<<<dimGrid,dimBlock>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,kstart, iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(22),0)



         iSize = nx2
         jSize = ny2
         kSize = (nz2-2) - 3 + 1

         blockSize = fixedBlockSize
         subSections = iSize/blockSize
         if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
         dimGrid = dim3(kSize,jSize*subSections,1)
         dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(23),0)
         call krnl_rhs_zdis2<<<dimGrid,dimBlock>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,iSize,subSections)
         istat = cudaEventRecord(dynamicEvents(24),0)



         iSize = nx2
         jSize = ny2
         kstart = nz2-1
         
        blockSize = fixedBlockSize
        subSections = iSize/blockSize
        if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
        dimGrid = dim3(jSize*subSections,1,1)
        dimBlock = dim3(blockSize,1,1)

         istat = cudaEventRecord(dynamicEvents(25),0)
         call krnl_rhs_zdis1<<<dimGrid,dimBlock>>>
     >                (u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp, kstart, iSize, subSections)
         istat = cudaEventRecord(dynamicEvents(26),0)

         if (timeron) call timer_stop(t_rhsz)
         


         iSize = nx2
         jSize = ny2
         kSize = nz2
         
         blockSize =fixedBlockSize     
         subSections = iSize/blockSize
         if(mod(iSize,blockSize) .ne. 0)subSections=subSections+1
         dimGrid = dim3(kSize,jSize*subSections,1)
         dimBlock = dim3(blockSize,1,1)
  
         istat = cudaEventRecord(dynamicEvents(27),0)
         call krnl_rhs_final<<<dimGrid,dimBlock>>>
     >                (rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dt,iSize,subSections)
         istat = cudaEventRecord(dynamicEvents(28),0)
         if (timeron) call timer_stop(t_rhs)


         if(syncForTimers .gt. 0) then
            istat= cudaThreadSynchronize()
            tink = 0
             istat = cudaEventElapsedTime(tink, dynamicEvents(1),
     >                                          dynamicEvents(2))
            call timer_inc(t_kernel,dble(tink/1000.0))
            call timer_inc(t_rhs,dble(tink/1000.0))
            !print *, 'r-> ', tink/1000.0

            do i=3,9,2
               istat = cudaEventElapsedTime(tink, dynamicEvents(i),
     >                                          dynamicEvents(i+1))
               call timer_inc(t_kernel,dble(tink/1000.0))          
               call timer_inc(t_rhsx,dble(tink/1000.0))
               call timer_inc(t_rhs,dble(tink/1000.0))
               !print *, 'r-> ', tink/1000.0
            end do
            do i=11,17,2
               istat = cudaEventElapsedTime(tink, dynamicEvents(i),
     >                                          dynamicEvents(i+1))
               call timer_inc(t_kernel,dble(tink/1000.0))
               call timer_inc(t_rhsy,dble(tink/1000.0))
               call timer_inc(t_rhs,dble(tink/1000.0))
               !print *, 'r-> ', tink/1000.0
            end do
            do i=19,25,2
               istat = cudaEventElapsedTime(tink, dynamicEvents(i),
     >                                          dynamicEvents(i+1))
               call timer_inc(t_kernel,dble(tink/1000.0))
               call timer_inc(t_rhsz,dble(tink/1000.0))
               call timer_inc(t_rhs,dble(tink/1000.0))
               !print *, 'r-> ', tink/1000.0
            end do
            istat = cudaEventElapsedTime(tink, dynamicEvents(27),
     >                                          dynamicEvents(28))
            call timer_inc(t_kernel,dble(tink/1000.0))
            call timer_inc(t_rhs,dble(tink/1000.0))
            !print *, 'r-> ', tink/1000.0
            
         endif






       else !BEGIN NO GPU




 
       if (timeron) call timer_start(t_rhs)
c---------------------------------------------------------------------
c      compute the reciprocal of density, and the kinetic energy, 
c      and the speed of sound. 
c---------------------------------------------------------------------

       do    k = 0, grid_points(3)-1
          do    j = 0, grid_points(2)-1
             do    i = 0, grid_points(1)-1
                rho_inv = 1.0d0/u(1,i,j,k)
                rho_i(i,j,k) = rho_inv
                us(i,j,k) = u(2,i,j,k) * rho_inv
                vs(i,j,k) = u(3,i,j,k) * rho_inv
                ws(i,j,k) = u(4,i,j,k) * rho_inv
                square(i,j,k)     = 0.5d0* (
     >                        u(2,i,j,k)*u(2,i,j,k) + 
     >                        u(3,i,j,k)*u(3,i,j,k) +
     >                        u(4,i,j,k)*u(4,i,j,k) ) * rho_inv
                qs(i,j,k) = square(i,j,k) * rho_inv
c---------------------------------------------------------------------
c               (don't need speed and ainx until the lhs computation)
c---------------------------------------------------------------------
                aux = c1c2*rho_inv* (u(5,i,j,k) - square(i,j,k))
                speed(i,j,k) = dsqrt(aux)
             end do
          end do
       end do

c---------------------------------------------------------------------
c copy the exact forcing term to the right hand side;  because 
c this forcing term is known, we can store it on the whole grid
c including the boundary                   
c---------------------------------------------------------------------

       do    k = 0, grid_points(3)-1
          do    j = 0, grid_points(2)-1
             do    i = 0, grid_points(1)-1
                do    m = 1, 5
                   rhs(m,i,j,k) = forcing(m,i,j,k)
                end do
             end do
          end do
       end do

c---------------------------------------------------------------------
c      compute xi-direction fluxes 
c---------------------------------------------------------------------
       if (timeron) call timer_start(t_rhsx)
       do    k = 1, nz2
          do    j = 1, ny2
             do    i = 1, nx2
                uijk = us(i,j,k)
                up1  = us(i+1,j,k)
                um1  = us(i-1,j,k)

                rhs(1,i,j,k) = rhs(1,i,j,k) + dx1tx1 * 
     >                    (u(1,i+1,j,k) - 2.0d0*u(1,i,j,k) + 
     >                     u(1,i-1,j,k)) -
     >                    tx2 * (u(2,i+1,j,k) - u(2,i-1,j,k))

                rhs(2,i,j,k) = rhs(2,i,j,k) + dx2tx1 * 
     >                    (u(2,i+1,j,k) - 2.0d0*u(2,i,j,k) + 
     >                     u(2,i-1,j,k)) +
     >                    xxcon2*con43 * (up1 - 2.0d0*uijk + um1) -
     >                    tx2 * (u(2,i+1,j,k)*up1 - 
     >                           u(2,i-1,j,k)*um1 +
     >                           (u(5,i+1,j,k)- square(i+1,j,k)-
     >                            u(5,i-1,j,k)+ square(i-1,j,k))*
     >                            c2)

                rhs(3,i,j,k) = rhs(3,i,j,k) + dx3tx1 * 
     >                    (u(3,i+1,j,k) - 2.0d0*u(3,i,j,k) +
     >                     u(3,i-1,j,k)) +
     >                    xxcon2 * (vs(i+1,j,k) - 2.0d0*vs(i,j,k) +
     >                              vs(i-1,j,k)) -
     >                    tx2 * (u(3,i+1,j,k)*up1 - 
     >                           u(3,i-1,j,k)*um1)

                rhs(4,i,j,k) = rhs(4,i,j,k) + dx4tx1 * 
     >                    (u(4,i+1,j,k) - 2.0d0*u(4,i,j,k) +
     >                     u(4,i-1,j,k)) +
     >                    xxcon2 * (ws(i+1,j,k) - 2.0d0*ws(i,j,k) +
     >                              ws(i-1,j,k)) -
     >                    tx2 * (u(4,i+1,j,k)*up1 - 
     >                           u(4,i-1,j,k)*um1)

                rhs(5,i,j,k) = rhs(5,i,j,k) + dx5tx1 * 
     >                    (u(5,i+1,j,k) - 2.0d0*u(5,i,j,k) +
     >                     u(5,i-1,j,k)) +
     >                    xxcon3 * (qs(i+1,j,k) - 2.0d0*qs(i,j,k) +
     >                              qs(i-1,j,k)) +
     >                    xxcon4 * (up1*up1 -       2.0d0*uijk*uijk + 
     >                              um1*um1) +
     >                    xxcon5 * (u(5,i+1,j,k)*rho_i(i+1,j,k) - 
     >                              2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
     >                              u(5,i-1,j,k)*rho_i(i-1,j,k)) -
     >                    tx2 * ( (c1*u(5,i+1,j,k) - 
     >                             c2*square(i+1,j,k))*up1 -
     >                            (c1*u(5,i-1,j,k) - 
     >                             c2*square(i-1,j,k))*um1 )
             end do
          end do

c---------------------------------------------------------------------
c      add fourth order xi-direction dissipation               
c---------------------------------------------------------------------

          do    j = 1, ny2
             i = 1
             do    m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k)- dssp * 
     >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) +
     >                            u(m,i+2,j,k))
             end do

             i = 2
             do    m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (-4.0d0*u(m,i-1,j,k) + 6.0d0*u(m,i,j,k) -
     >                      4.0d0*u(m,i+1,j,k) + u(m,i+2,j,k))
             end do
          end do

          do    j = 1, ny2
             do  i = 3, nx2-2
                do     m = 1, 5
                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (  u(m,i-2,j,k) - 4.0d0*u(m,i-1,j,k) + 
     >                     6.0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) + 
     >                         u(m,i+2,j,k) )
                end do
             end do
          end do

          do    j = 1, ny2
             i = nx2-1
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i-2,j,k) - 4.0d0*u(m,i-1,j,k) + 
     >                      6.0d0*u(m,i,j,k) - 4.0d0*u(m,i+1,j,k) )
             end do

             i = nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i-2,j,k) - 4.d0*u(m,i-1,j,k) +
     >                      5.d0*u(m,i,j,k) )
             end do
          end do
       end do
       if (timeron) call timer_stop(t_rhsx)

c---------------------------------------------------------------------
c      compute eta-direction fluxes 
c---------------------------------------------------------------------
       if (timeron) call timer_start(t_rhsy)
       do     k = 1, nz2
          do     j = 1, ny2
             do     i = 1, nx2
                vijk = vs(i,j,k)
                vp1  = vs(i,j+1,k)
                vm1  = vs(i,j-1,k)
                rhs(1,i,j,k) = rhs(1,i,j,k) + dy1ty1 * 
     >                   (u(1,i,j+1,k) - 2.0d0*u(1,i,j,k) + 
     >                    u(1,i,j-1,k)) -
     >                   ty2 * (u(3,i,j+1,k) - u(3,i,j-1,k))
                rhs(2,i,j,k) = rhs(2,i,j,k) + dy2ty1 * 
     >                   (u(2,i,j+1,k) - 2.0d0*u(2,i,j,k) + 
     >                    u(2,i,j-1,k)) +
     >                   yycon2 * (us(i,j+1,k) - 2.0d0*us(i,j,k) + 
     >                             us(i,j-1,k)) -
     >                   ty2 * (u(2,i,j+1,k)*vp1 - 
     >                          u(2,i,j-1,k)*vm1)
                rhs(3,i,j,k) = rhs(3,i,j,k) + dy3ty1 * 
     >                   (u(3,i,j+1,k) - 2.0d0*u(3,i,j,k) + 
     >                    u(3,i,j-1,k)) +
     >                   yycon2*con43 * (vp1 - 2.0d0*vijk + vm1) -
     >                   ty2 * (u(3,i,j+1,k)*vp1 - 
     >                          u(3,i,j-1,k)*vm1 +
     >                          (u(5,i,j+1,k) - square(i,j+1,k) - 
     >                           u(5,i,j-1,k) + square(i,j-1,k))
     >                          *c2)
                rhs(4,i,j,k) = rhs(4,i,j,k) + dy4ty1 * 
     >                   (u(4,i,j+1,k) - 2.0d0*u(4,i,j,k) + 
     >                    u(4,i,j-1,k)) +
     >                   yycon2 * (ws(i,j+1,k) - 2.0d0*ws(i,j,k) + 
     >                             ws(i,j-1,k)) -
     >                   ty2 * (u(4,i,j+1,k)*vp1 - 
     >                          u(4,i,j-1,k)*vm1)
                rhs(5,i,j,k) = rhs(5,i,j,k) + dy5ty1 * 
     >                   (u(5,i,j+1,k) - 2.0d0*u(5,i,j,k) + 
     >                    u(5,i,j-1,k)) +
     >                   yycon3 * (qs(i,j+1,k) - 2.0d0*qs(i,j,k) + 
     >                             qs(i,j-1,k)) +
     >                   yycon4 * (vp1*vp1       - 2.0d0*vijk*vijk + 
     >                             vm1*vm1) +
     >                   yycon5 * (u(5,i,j+1,k)*rho_i(i,j+1,k) - 
     >                             2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
     >                             u(5,i,j-1,k)*rho_i(i,j-1,k)) -
     >                   ty2 * ((c1*u(5,i,j+1,k) - 
     >                           c2*square(i,j+1,k)) * vp1 -
     >                          (c1*u(5,i,j-1,k) - 
     >                           c2*square(i,j-1,k)) * vm1)
             end do
          end do

c---------------------------------------------------------------------
c      add fourth order eta-direction dissipation         
c---------------------------------------------------------------------

          j = 1
          do     i = 1, nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k)- dssp * 
     >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) +
     >                            u(m,i,j+2,k))
             end do
          end do

          j = 2
          do     i = 1, nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (-4.0d0*u(m,i,j-1,k) + 6.0d0*u(m,i,j,k) -
     >                      4.0d0*u(m,i,j+1,k) + u(m,i,j+2,k))
             end do
          end do

          do    j = 3, ny2-2
             do  i = 1,nx2
                do     m = 1, 5
                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (  u(m,i,j-2,k) - 4.0d0*u(m,i,j-1,k) + 
     >                     6.0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) + 
     >                         u(m,i,j+2,k) )
                end do
             end do
          end do
 
          j = ny2-1
          do     i = 1, nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i,j-2,k) - 4.0d0*u(m,i,j-1,k) + 
     >                      6.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j+1,k) )
             end do
          end do

          j = ny2
          do     i = 1, nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i,j-2,k) - 4.d0*u(m,i,j-1,k) +
     >                      5.d0*u(m,i,j,k) )
             end do
          end do
       end do
       if (timeron) call timer_stop(t_rhsy)

c---------------------------------------------------------------------
c      compute zeta-direction fluxes 
c---------------------------------------------------------------------
       if (timeron) call timer_start(t_rhsz)
       do    k = 1, nz2
          do     j = 1, ny2
             do     i = 1, nx2
                wijk = ws(i,j,k)
                wp1  = ws(i,j,k+1)
                wm1  = ws(i,j,k-1)

                rhs(1,i,j,k) = rhs(1,i,j,k) + dz1tz1 * 
     >                   (u(1,i,j,k+1) - 2.0d0*u(1,i,j,k) + 
     >                    u(1,i,j,k-1)) -
     >                   tz2 * (u(4,i,j,k+1) - u(4,i,j,k-1))
                rhs(2,i,j,k) = rhs(2,i,j,k) + dz2tz1 * 
     >                   (u(2,i,j,k+1) - 2.0d0*u(2,i,j,k) + 
     >                    u(2,i,j,k-1)) +
     >                   zzcon2 * (us(i,j,k+1) - 2.0d0*us(i,j,k) + 
     >                             us(i,j,k-1)) -
     >                   tz2 * (u(2,i,j,k+1)*wp1 - 
     >                          u(2,i,j,k-1)*wm1)
                rhs(3,i,j,k) = rhs(3,i,j,k) + dz3tz1 * 
     >                   (u(3,i,j,k+1) - 2.0d0*u(3,i,j,k) + 
     >                    u(3,i,j,k-1)) +
     >                   zzcon2 * (vs(i,j,k+1) - 2.0d0*vs(i,j,k) + 
     >                             vs(i,j,k-1)) -
     >                   tz2 * (u(3,i,j,k+1)*wp1 - 
     >                          u(3,i,j,k-1)*wm1)
                rhs(4,i,j,k) = rhs(4,i,j,k) + dz4tz1 * 
     >                   (u(4,i,j,k+1) - 2.0d0*u(4,i,j,k) + 
     >                    u(4,i,j,k-1)) +
     >                   zzcon2*con43 * (wp1 - 2.0d0*wijk + wm1) -
     >                   tz2 * (u(4,i,j,k+1)*wp1 - 
     >                          u(4,i,j,k-1)*wm1 +
     >                          (u(5,i,j,k+1) - square(i,j,k+1) - 
     >                           u(5,i,j,k-1) + square(i,j,k-1))
     >                          *c2)
                rhs(5,i,j,k) = rhs(5,i,j,k) + dz5tz1 * 
     >                   (u(5,i,j,k+1) - 2.0d0*u(5,i,j,k) + 
     >                    u(5,i,j,k-1)) +
     >                   zzcon3 * (qs(i,j,k+1) - 2.0d0*qs(i,j,k) + 
     >                             qs(i,j,k-1)) +
     >                   zzcon4 * (wp1*wp1 - 2.0d0*wijk*wijk + 
     >                             wm1*wm1) +
     >                   zzcon5 * (u(5,i,j,k+1)*rho_i(i,j,k+1) - 
     >                             2.0d0*u(5,i,j,k)*rho_i(i,j,k) +
     >                             u(5,i,j,k-1)*rho_i(i,j,k-1)) -
     >                   tz2 * ( (c1*u(5,i,j,k+1) - 
     >                            c2*square(i,j,k+1))*wp1 -
     >                           (c1*u(5,i,j,k-1) - 
     >                            c2*square(i,j,k-1))*wm1)
             end do
          end do
       end do

c---------------------------------------------------------------------
c      add fourth order zeta-direction dissipation                
c---------------------------------------------------------------------

       k = 1
       do     j = 1, ny2
          do     i = 1, nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k)- dssp * 
     >                    ( 5.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) +
     >                            u(m,i,j,k+2))
             end do
          end do
       end do

       k = 2
       do     j = 1, ny2
          do     i = 1, nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (-4.0d0*u(m,i,j,k-1) + 6.0d0*u(m,i,j,k) -
     >                      4.0d0*u(m,i,j,k+1) + u(m,i,j,k+2))
             end do
          end do
       end do

       do     k = 3, nz2-2
          do     j = 1, ny2
             do     i = 1,nx2
                do     m = 1, 5
                   rhs(m,i,j,k) = rhs(m,i,j,k) - dssp * 
     >                    (  u(m,i,j,k-2) - 4.0d0*u(m,i,j,k-1) + 
     >                     6.0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) + 
     >                         u(m,i,j,k+2) )
                end do
             end do
          end do
       end do
 
       k = nz2-1
       do     j = 1, ny2
          do     i = 1, nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i,j,k-2) - 4.0d0*u(m,i,j,k-1) + 
     >                      6.0d0*u(m,i,j,k) - 4.0d0*u(m,i,j,k+1) )
             end do
          end do
       end do

       k = nz2
       do     j = 1, ny2
          do     i = 1, nx2
             do     m = 1, 5
                rhs(m,i,j,k) = rhs(m,i,j,k) - dssp *
     >                    ( u(m,i,j,k-2) - 4.d0*u(m,i,j,k-1) +
     >                      5.d0*u(m,i,j,k) )
             end do
          end do
       end do
       if (timeron) call timer_stop(t_rhsz)

       do    k = 1, nz2
          do    j = 1, ny2
             do    i = 1, nx2
                do    m = 1, 5
                   rhs(m,i,j,k) = rhs(m,i,j,k) * dt
                end do
             end do
          end do
       end do
       if (timeron) call timer_stop(t_rhs)

       endif !END GPU SWITCH
    
       return
       end




