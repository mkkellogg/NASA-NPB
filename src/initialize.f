
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine  initialize
       use glob
       use cudafor
       use util_kernels
       use lhs_kernels
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c This subroutine initializes the field variable u using 
c tri-linear transfinite interpolation of the boundary values     
c---------------------------------------------------------------------

       include 'header.h'
  
       integer i, j, k, m, ix, iy, iz, istat
       double precision  xi, eta, zeta, Pface(5,3,2), Pxi, Peta, 
     >                   Pzeta, temp(5)
       integer iSize, jSize, kSize, sharedSize, blockSize, ni, nj, nk
 
c---------------------------------------------------------------------
c  Later (in compute_rhs) we compute 1/u for every element. A few of 
c  the corner elements are not used, but it convenient (and faster) 
c  to compute the whole thing with a simple loop. Make sure those 
c  values are nonzero by initializing the whole thing here. 
c---------------------------------------------------------------------


      if(doGPU .gt. 0) then !BEGIN GPU MODE


      do k = 0, grid_points(3)-1
         do j = 0, grid_points(2)-1
            do i = 0, grid_points(1)-1
               u_r(i,j,k,1) = 1.0
               u_r(i,j,k,2) = 0.0
               u_r(i,j,k,3) = 0.0
               u_r(i,j,k,4) = 0.0
               u_r(i,j,k,5) = 1.0
            end do
         end do
      end do

c---------------------------------------------------------------------
c first store the "interpolated" values everywhere on the grid    
c---------------------------------------------------------------------
          do  k = 0, grid_points(3)-1
             zeta = dble(k) * dnzm1
             do  j = 0, grid_points(2)-1
                eta = dble(j) * dnym1
                do   i = 0, grid_points(1)-1
                   xi = dble(i) * dnxm1
                  
                   do ix = 1, 2
                      Pxi = dble(ix-1)
                      call exact_solution(Pxi, eta, zeta, 
     >                                    Pface(1,1,ix))
                   end do

                   do    iy = 1, 2
                      Peta = dble(iy-1)
                      call exact_solution(xi, Peta, zeta, 
     >                                    Pface(1,2,iy))
                   end do

                   do    iz = 1, 2
                      Pzeta = dble(iz-1)
                      call exact_solution(xi, eta, Pzeta,   
     >                                    Pface(1,3,iz))
                   end do

                   do   m = 1, 5
                      Pxi   = xi   * Pface(m,1,2) + 
     >                        (1.0d0-xi)   * Pface(m,1,1)
                      Peta  = eta  * Pface(m,2,2) + 
     >                        (1.0d0-eta)  * Pface(m,2,1)
                      Pzeta = zeta * Pface(m,3,2) + 
     >                        (1.0d0-zeta) * Pface(m,3,1)
 
                      u_r(i,j,k,m) = Pxi + Peta + Pzeta - 
     >                          Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + 
     >                          Pxi*Peta*Pzeta

                   end do
                end do
             end do
          end do


c---------------------------------------------------------------------
c now store the exact values on the boundaries        
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c west face                                                  
c---------------------------------------------------------------------

       xi = 0.0d0
       i  = 0
       do  k = 0, grid_points(3)-1
          zeta = dble(k) * dnzm1
          do   j = 0, grid_points(2)-1
             eta = dble(j) * dnym1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u_r(i,j,k,m) = temp(m)
             end do
          end do
       end do

c---------------------------------------------------------------------
c east face                                                      
c---------------------------------------------------------------------

       xi = 1.0d0
       i  = grid_points(1)-1
       do   k = 0, grid_points(3)-1
          zeta = dble(k) * dnzm1
          do   j = 0, grid_points(2)-1
             eta = dble(j) * dnym1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u_r(i,j,k,m) = temp(m)
             end do
          end do
       end do

c---------------------------------------------------------------------
c south face                                                 
c---------------------------------------------------------------------

       eta = 0.0d0
       j   = 0
       do  k = 0, grid_points(3)-1
          zeta = dble(k) * dnzm1
          do   i = 0, grid_points(1)-1
             xi = dble(i) * dnxm1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u_r(i,j,k,m) = temp(m)
             end do
          end do
       end do


c---------------------------------------------------------------------
c north face                                    
c---------------------------------------------------------------------

       eta = 1.0d0
       j   = grid_points(2)-1
       do   k = 0, grid_points(3)-1
          zeta = dble(k) * dnzm1
          do   i = 0, grid_points(1)-1
             xi = dble(i) * dnxm1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u_r(i,j,k,m) = temp(m)
             end do
          end do
       end do

c---------------------------------------------------------------------
c bottom face                                       
c---------------------------------------------------------------------

       zeta = 0.0d0
       k    = 0
       do   j = 0, grid_points(2)-1
          eta = dble(j) * dnym1
          do   i =0, grid_points(1)-1
             xi = dble(i) *dnxm1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u_r(i,j,k,m) = temp(m)
             end do
          end do
       end do

c---------------------------------------------------------------------
c top face     
c---------------------------------------------------------------------

       zeta = 1.0d0
       k    = grid_points(3)-1
       do   j = 0, grid_points(2)-1
          eta = dble(j) * dnym1
          do   i =0, grid_points(1)-1
             xi = dble(i) * dnxm1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u_r(i,j,k,m) = temp(m)
             end do
          end do
       end do

      
       if (timeron) call timer_start(t_gcomm)
       istat = cudaMemcpy(forcing_dv_r,
     >                   forcing_r,5*(IMAXP+1)*(JMAXP+1)*KMAX)
       istat = cudaMemcpy(u_dv_r, u_r,5 * (IMAXP+1) * (JMAXP+1) * KMAX)
       if (timeron) call timer_stop(t_gcomm)
     
c---------------------------------------------------------------------
c         swap 'i' and 'm' dimensions from u_dv_r to u_dv
c--------------------------------------------------------------------

        iSize = IMAXP+1
        jSize = JMAXP+1
        kSize = KMAX

        blockSize = iSize

        dimGrid = dim3(kSize,jSize,1)
        dimBlock = dim3(blockSize,1,1)

       call krnl_util_swap4D_mtofirst<<<dimGrid,dimBlock>>>(
     >                         u_dv_r,u_dv,
     >                         imax, jmax, kmax, imaxp, jmaxp)

c---------------------------------------------------------------------
c         initialize lhs
c--------------------------------------------------------------------

       kSize = nz2
       jSize = ny2
       nj = jSize
       ni = nx2+1

       blockSize = 5
       dimGrid = dim3(kSize,jSize,1)
       dimBlock = dim3(blockSize,1,1)
            
       call krnl_lhs_init<<<dimGrid,dimBlock>>>
     >                (lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r,  
     >                imax, jmax, kmax, imaxp, jmaxp, ni)



c---------------------------------------------------------------------
c         initialize lhsj
c--------------------------------------------------------------------


        nj = ny2+1
        iSize = nx2
        jSize = ny2
        kSize = nz2

        blockSize = 5
        dimGrid = dim3(kSize,iSize,1)
        dimBlock = dim3(blockSize,1,1)

        call krnl_lhs_initj<<<dimGrid,dimBlock>>>(
     >          lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r,
     >          imax, jmax, kmax, imaxp, jmaxp, nj)


c---------------------------------------------------------------------
c         initialize lhsk
c--------------------------------------------------------------------


       iSize = nx2
       jSize = ny2
       kSize = nz2

       nk = nz2+1

        blockSize = 5
        dimGrid = dim3(jSize,iSize,1)
        dimBlock = dim3(blockSize,1,1)

        call krnl_lhs_initk<<<dimGrid,dimBlock>>>(
     >          lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r,
     >          imax, jmax, kmax, imaxp, jmaxp, nk)

c---------------------------------------------------------------------
c         swap 'i' and 'm' dimensions from lhs_x_dv_r to lhs_x_dv
c
c--------------------------------------------------------------------

        iSize = IMAXP+1
        jSize = JMAXP+1
        kSize = KMAX

        blockSize = iSize
        dimGrid = dim3(kSize,jSize,1)
        dimBlock = dim3(blockSize,1,1)


       call krnl_util_swap4D_mtofirst<<<dimGrid,dimBlock>>>(
     >                         lhs_x_dv_r,lhs_x_dv,
     >                         imax, jmax, kmax, imaxp, imaxp)
 
       call krnl_util_swap4D_mtofirst<<<dimGrid,dimBlock>>>(
     >                         lhsp_x_dv_r,lhsp_x_dv,
     >                         imax, jmax, kmax, imaxp, imaxp)

       call krnl_util_swap4D_mtofirst<<<dimGrid,dimBlock>>>(
     >                         lhsm_x_dv_r,lhsm_x_dv,
     >                         imax, jmax, kmax, imaxp, imaxp)




       call krnl_util_swap4D_itoj<<<dimGrid,dimBlock>>>(
     >                         lhs_x_dv_r,lhs_x_dv_j_r,
     >                         imax, jmax, kmax, imaxp, imaxp)

       call krnl_util_swap4D_itoj<<<dimGrid,dimBlock>>>(
     >                         lhsp_x_dv_r,lhsp_x_dv_j_r,
     >                         imax, jmax, kmax, imaxp, imaxp)

       call krnl_util_swap4D_itoj<<<dimGrid,dimBlock>>>(
     >                         lhsm_x_dv_r,lhsm_x_dv_j_r,
     >                         imax, jmax, kmax, imaxp, imaxp)




       call krnl_util_swap4D_itok<<<dimGrid,dimBlock>>>(
     >                         lhs_x_dv_r,lhs_x_dv_k_r,
     >                         imax, jmax, kmax, imaxp, imaxp)






       else !BEGIN NON-GPU MODE








       do k = 0, grid_points(3)-1
         do j = 0, grid_points(2)-1
            do i = 0, grid_points(1)-1
               u(1,i,j,k) = 1.0
               u(2,i,j,k) = 0.0
               u(3,i,j,k) = 0.0
               u(4,i,j,k) = 0.0
               u(5,i,j,k) = 1.0
            end do
         end do
       end do

c---------------------------------------------------------------------
c first store the "interpolated" values everywhere on the grid    
c---------------------------------------------------------------------
          do  k = 0, grid_points(3)-1
             zeta = dble(k) * dnzm1
             do  j = 0, grid_points(2)-1
                eta = dble(j) * dnym1
                do   i = 0, grid_points(1)-1
                   xi = dble(i) * dnxm1
                  
                   do ix = 1, 2
                      Pxi = dble(ix-1)
                      call exact_solution(Pxi, eta, zeta, 
     >                                    Pface(1,1,ix))
                   end do

                   do    iy = 1, 2
                      Peta = dble(iy-1)
                      call exact_solution(xi, Peta, zeta, 
     >                                    Pface(1,2,iy))
                   end do

                   do    iz = 1, 2
                      Pzeta = dble(iz-1)
                      call exact_solution(xi, eta, Pzeta,   
     >                                    Pface(1,3,iz))
                   end do

                   do   m = 1, 5
                      Pxi   = xi   * Pface(m,1,2) + 
     >                        (1.0d0-xi)   * Pface(m,1,1)
                      Peta  = eta  * Pface(m,2,2) + 
     >                        (1.0d0-eta)  * Pface(m,2,1)
                      Pzeta = zeta * Pface(m,3,2) + 
     >                        (1.0d0-zeta) * Pface(m,3,1)
 
                      u(m,i,j,k) = Pxi + Peta + Pzeta - 
     >                          Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + 
     >                          Pxi*Peta*Pzeta

                   end do
                end do
             end do
          end do


c---------------------------------------------------------------------
c now store the exact values on the boundaries        
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c west face                                                  
c---------------------------------------------------------------------

       xi = 0.0d0
       i  = 0
       do  k = 0, grid_points(3)-1
          zeta = dble(k) * dnzm1
          do   j = 0, grid_points(2)-1
             eta = dble(j) * dnym1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u(m,i,j,k) = temp(m)
             end do
          end do
       end do

c---------------------------------------------------------------------
c east face                                                      
c---------------------------------------------------------------------

       xi = 1.0d0
       i  = grid_points(1)-1
       do   k = 0, grid_points(3)-1
          zeta = dble(k) * dnzm1
          do   j = 0, grid_points(2)-1
             eta = dble(j) * dnym1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u(m,i,j,k) = temp(m)
             end do
          end do
       end do

c---------------------------------------------------------------------
c south face                                                 
c---------------------------------------------------------------------

       eta = 0.0d0
       j   = 0
       do  k = 0, grid_points(3)-1
          zeta = dble(k) * dnzm1
          do   i = 0, grid_points(1)-1
             xi = dble(i) * dnxm1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u(m,i,j,k) = temp(m)
             end do
          end do
       end do


c---------------------------------------------------------------------
c north face                                    
c---------------------------------------------------------------------

       eta = 1.0d0
       j   = grid_points(2)-1
       do   k = 0, grid_points(3)-1
          zeta = dble(k) * dnzm1
          do   i = 0, grid_points(1)-1
             xi = dble(i) * dnxm1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u(m,i,j,k) = temp(m)
             end do
          end do
       end do

c---------------------------------------------------------------------
c bottom face                                       
c---------------------------------------------------------------------

       zeta = 0.0d0
       k    = 0
       do   j = 0, grid_points(2)-1
          eta = dble(j) * dnym1
          do   i =0, grid_points(1)-1
             xi = dble(i) *dnxm1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u(m,i,j,k) = temp(m)
             end do
          end do
       end do

c---------------------------------------------------------------------
c top face     
c---------------------------------------------------------------------

       zeta = 1.0d0
       k    = grid_points(3)-1
       do   j = 0, grid_points(2)-1
          eta = dble(j) * dnym1
          do   i =0, grid_points(1)-1
             xi = dble(i) * dnxm1
             call exact_solution(xi, eta, zeta, temp)
             do   m = 1, 5
                u(m,i,j,k) = temp(m)
             end do
          end do
       end do




       endif !END GPU SWITCH
 
 




       return
       end







       subroutine lhsinit(ni, nj)
       use glob 

       include 'header.h'

       integer ni, nj

       integer j, m

c---------------------------------------------------------------------
c     zap the whole left hand side for starters
c     set all diagonal values to 1. This is overkill, but convenient
c---------------------------------------------------------------------
       do j = 1, nj

          do   m = 1, 5
             lhs (m,0,j) = 0.0d0
             lhsp(m,0,j) = 0.0d0
             lhsm(m,0,j) = 0.0d0
          end do
          lhs (3,0,j) = 1.0d0
          lhsp(3,0,j) = 1.0d0
          lhsm(3,0,j) = 1.0d0


          do   m = 1, 5
             lhs (m,ni,j) = 0.0d0
             lhsp(m,ni,j) = 0.0d0
             lhsm(m,ni,j) = 0.0d0
          end do
          lhs (3,ni,j) = 1.0d0
          lhsp(3,ni,j) = 1.0d0
          lhsm(3,ni,j) = 1.0d0


       end do
 
       return
       end


       subroutine lhsinitj(nj, ni)
       use glob

       include 'header.h'

       integer nj, ni

       integer i, m

c---------------------------------------------------------------------
c     zap the whole left hand side for starters
c     set all diagonal values to 1. This is overkill, but convenient
c---------------------------------------------------------------------
       do i = 1, ni
          do   m = 1, 5
             lhs (m,i,0) = 0.0d0
             lhsp(m,i,0) = 0.0d0
             lhsm(m,i,0) = 0.0d0
             lhs (m,i,nj) = 0.0d0
             lhsp(m,i,nj) = 0.0d0
             lhsm(m,i,nj) = 0.0d0
          end do
          lhs (3,i,0) = 1.0d0
          lhsp(3,i,0) = 1.0d0
          lhsm(3,i,0) = 1.0d0
          lhs (3,i,nj) = 1.0d0
          lhsp(3,i,nj) = 1.0d0
          lhsm(3,i,nj) = 1.0d0
       end do
 
       return
       end





