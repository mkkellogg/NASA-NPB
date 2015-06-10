         module z_solve_kernels
         implicit none
         contains



         attributes(global) subroutine krnl_z_solve_init( 
     >   lhs_x_dv_k_r, rho_i_dv_k_r, ws_dv_k_r,  imax, jmax, kmax,
     >                  imaxp, jmaxp,  dttz1, dttz2, c2dttz1,
     >                  con43, c1c5, c3c4, dz1, dz4, dz5,dzmax,
     >                    upperk1, upperk2, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:KMAX-1,0:IMAXP,
     >                          0:IMAXP,5) :: lhs_x_dv_k_r
        double precision, dimension(0:KMAX-1, 0:JMAXP, 
     >                          0:IMAXP) :: rho_i_dv_k_r
        double precision, dimension(0:KMAX-1, 0:JMAXP,
     >                          0:IMAXP)::ws_dv_k_r

        double precision,value :: dttz1,dttz2,c2dttz1, con43, c1c5,
     >                            c3c4, dz1, dz4, dz5, dzmax

        integer, value :: upperk1, upperk2, subSections

        integer i,j,k, modF, blockSize, ioff, subISize

        double precision tempd, ru1, o1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        subISize = blockSize - 2
        j = blockIDx%x
        i = (blockIDx%y-1)/subSections +1
        ioff = (mod(blockIDx%y-1,subSections)*subISize)
        k = threadIDx%x + ioff -1

             if(k .le. upperk1) then

                ru1 = c3c4*rho_i_dv_k_r(k,j,i)
                temp(threadIDx%x) = ws_dv_k_r(k,j,i)
                temp(threadIDx%x+blockSize) = dmax1(dz4 + con43 * ru1,
     >                                        dz5 + c1c5 * ru1,
     >                                        dzmax + ru1,
     >                                        dz1)


             endif

             call syncthreads()

             if(k .le. upperk2 .and. threadIDx%x .gt. 1
     >                         .and. threadIDx%x .lt. blockSize) then

                lhs_x_dv_k_r(k,j,i,1) =   0.0d0
                lhs_x_dv_k_r(k,j,i,2) = -dttz2 * temp(threadIDx%x-1) -
     >                         dttz1 * temp(threadIDx%x+blockSize-1)
                lhs_x_dv_k_r(k,j,i,3) =   1.0d0 + c2dttz1 *
     >                          temp(threadIDx%x+blockSize)
                lhs_x_dv_k_r(k,j,i,4) =   dttz2 * temp(threadIDx%x+1) -
     >                          dttz1 * temp(threadIDx%x+blockSize+1)
                lhs_x_dv_k_r(k,j,i,5) =   0.0d0

             endif

         end subroutine krnl_z_solve_init










         attributes(global) subroutine krnl_z_solve_init_nc(
     >   lhs_x_dv_r, rho_i_dv, ws_dv,  imax, jmax, kmax,
     >                  imaxp, jmaxp,  dttz1, dttz2, c2dttz1,
     >                  con43, c1c5, c3c4, dz1, dz4, dz5,dzmax,
     >                    upperk1, upperk2, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:IMAXP,
     >                          0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP, 0:JMAXP,
     >                          0:KMAX-1) :: rho_i_dv
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                          0:KMAX-1)::ws_dv

        double precision,value :: dttz1,dttz2,c2dttz1, con43, c1c5,
     >                            c3c4, dz1, dz4, dz5, dzmax

        integer, value :: upperk1, upperk2, subSections

        integer i,j,k, modF, blockSize, ioff, subISize

        double precision tempd, ru1, o1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        subISize = blockSize - 2
        j = blockIDx%x
        i = (blockIDx%y-1)/subSections +1
        ioff = (mod(blockIDx%y-1,subSections)*subISize)
        k = threadIDx%x + ioff -1

             if(k .le. upperk1) then

                ru1 = c3c4*rho_i_dv(i,j,k)
                temp(threadIDx%x) = ws_dv(i,j,k)
                temp(threadIDx%x+blockSize) = dmax1(dz4 + con43 * ru1,
     >                                        dz5 + c1c5 * ru1,
     >                                        dzmax + ru1,
     >                                        dz1)


             endif

             call syncthreads()

             if(k .le. upperk2 .and. threadIDx%x .gt. 1
     >                         .and. threadIDx%x .lt. blockSize) then

                lhs_x_dv_r(i,j,k,1) =   0.0d0
                lhs_x_dv_r(i,j,k,2) = -dttz2 * temp(threadIDx%x-1) -
     >                         dttz1 * temp(threadIDx%x+blockSize-1)
                lhs_x_dv_r(i,j,k,3) =   1.0d0 + c2dttz1 *
     >                          temp(threadIDx%x+blockSize)
                lhs_x_dv_r(i,j,k,4) =   dttz2 * temp(threadIDx%x+1) -
     >                          dttz1 * temp(threadIDx%x+blockSize+1)
                lhs_x_dv_r(i,j,k,5) =   0.0d0

             endif

         end subroutine krnl_z_solve_init_nc










        attributes(global) subroutine krnl_z_solve_dis1(
     >        lhs_x_dv_r, imax, jmax, kmax, imaxp, jmaxp,
     >                    comz1,comz4,comz5,comz6, kstart,
     >                    iSize,subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r

        double precision,value :: comz1,comz4,comz5,comz6
        integer, value :: kstart, iSize, subSections

        integer i,j,k, blockSize, ioff

        double precision m

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        j = (blockIDx%x-1)/subSections + 1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff

        if(i .le. iSize) then


           if(kstart .le. 0) then
             k = 1
             lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz5
             lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4
             lhs_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5) + comz1

             k = 2
             lhs_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - comz4
             lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz6
             lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4
             lhs_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5) + comz1
           else
            k = kstart

             lhs_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1) + comz1
             lhs_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - comz4
             lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz6
             lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4

             k = kstart+1
             lhs_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1) + comz1
             lhs_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - comz4
             lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz5

           endif


        endif

        end subroutine krnl_z_solve_dis1






        attributes(global) subroutine krnl_z_solve_dis2(
     >           lhs_x_dv_r,imax, jmax, kmax, imaxp, jmaxp,
     >                    comz1,comz4,comz5,comz6,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r

        double precision,value :: comz1,comz4,comz5,comz6
        integer, value :: iSize, subSections

        integer i,j,k, blockSize, ioff
        double precision m

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = blockIDx%x+2
        j = (blockIDx%y-1)/subSections +1
        ioff = (mod(blockIDx%y-1,subSections)*blockSize)
        i = threadIDx%x + ioff 

        if(i .le. iSize) then

                lhs_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1) + comz1
                lhs_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - comz4
                lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz6
                lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4
                lhs_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5) + comz1

        endif

        end subroutine krnl_z_solve_dis2









       attributes(global) subroutine krnl_z_solve_postdis(
     >               lhs_x_dv_r, lhsm_x_dv_r, lhsp_x_dv_r, speed_dv,
     >                    imax, jmax, kmax, imaxp, jmaxp, dttz2,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1) :: speed_dv


        double precision,value ::dttz2
        integer, value :: iSize, subSections
        double precision m
        double precision, shared :: temp(*)
        integer i,j,k, blockSize, ioff

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        ioff = (mod(blockIDx%y-1,subSections)*blockSize)
        i = threadIDx%x + ioff 


        if(i .le. iSize) then
                lhsp_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1)
                lhsp_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - 
     >                            dttz2 * speed_dv(i,j,k-1)
                lhsp_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3)
                lhsp_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) + 
     >                            dttz2 * speed_dv(i,j,k+1)
                lhsp_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5)
                lhsm_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1)
                lhsm_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) + 
     >                            dttz2 * speed_dv(i,j,k-1)
                lhsm_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3)
                lhsm_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - 
     >                            dttz2 * speed_dv(i,j,k+1)
                lhsm_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5)        

        endif

        end subroutine krnl_z_solve_postdis









        attributes(global) subroutine krnl_z_solve_forward1(
     >                    lhs_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, kSize,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: kSize, iSize, subSections

        integer i,j,k, blockSize, ioff, m, k1, k2
 
        double precision fac1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        j = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff

        do k = 0, kSize-1
             k1 = k  + 1
             k2 = k  + 2

             if(i .le. iSize) then

                fac1      = 1.d0/lhs_x_dv_r(i,j,k,3)
                lhs_x_dv_r(i,j,k,4)  = fac1*lhs_x_dv_r(i,j,k,4)
                lhs_x_dv_r(i,j,k,5)  = fac1*lhs_x_dv_r(i,j,k,5)
                do    m = 1, 3
                   rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
                end do
                lhs_x_dv_r(i,j,k1,3) = lhs_x_dv_r(i,j,k1,3) -
     >                         lhs_x_dv_r(i,j,k1,2)*lhs_x_dv_r(i,j,k,4)
                lhs_x_dv_r(i,j,k1,4) = lhs_x_dv_r(i,j,k1,4) -
     >                         lhs_x_dv_r(i,j,k1,2)*lhs_x_dv_r(i,j,k,5)
                do    m = 1, 3
                   rhs_dv_r(i,j,k1,m) = rhs_dv_r(i,j,k1,m) -
     >                         lhs_x_dv_r(i,j,k1,2)*rhs_dv_r(i,j,k,m)
                end do
                lhs_x_dv_r(i,j,k2,2) = lhs_x_dv_r(i,j,k2,2) -
     >                         lhs_x_dv_r(i,j,k2,1)*lhs_x_dv_r(i,j,k,4)
                lhs_x_dv_r(i,j,k2,3) = lhs_x_dv_r(i,j,k2,3) -
     >                         lhs_x_dv_r(i,j,k2,1)*lhs_x_dv_r(i,j,k,5)
                do    m = 1, 3
                   rhs_dv_r(i,j,k2,m) = rhs_dv_r(i,j,k2,m) -
     >                         lhs_x_dv_r(i,j,k2,1)*rhs_dv_r(i,j,k,m)
       
               end do
             endif

        end do


        end subroutine krnl_z_solve_forward1
             








        attributes(global) subroutine krnl_z_solve_forward2(
     >                    lhs_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, kstart,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r

        integer, value :: kstart, iSize, subSections

        integer i,j,k, blockSize, ioff, m, k1, k2

        double precision fac1, fac2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        j = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff

        k = kstart
        k1 = k+1

        if(i .le. iSize) then

             fac1      = 1.d0/lhs_x_dv_r(i,j,k,3)
             lhs_x_dv_r(i,j,k,4)  = fac1*lhs_x_dv_r(i,j,k,4)
             lhs_x_dv_r(i,j,k,5)  = fac1*lhs_x_dv_r(i,j,k,5)
             do    m = 1, 3
                rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
             end do
             lhs_x_dv_r(i,j,k1,3) = lhs_x_dv_r(i,j,k1,3) -
     >                      lhs_x_dv_r(i,j,k1,2)*lhs_x_dv_r(i,j,k,4)
             lhs_x_dv_r(i,j,k1,4) = lhs_x_dv_r(i,j,k1,4) -
     >                      lhs_x_dv_r(i,j,k1,2)*lhs_x_dv_r(i,j,k,5)
             do    m = 1, 3
                rhs_dv_r(i,j,k1,m) = rhs_dv_r(i,j,k1,m) -
     >                      lhs_x_dv_r(i,j,k1,2)*rhs_dv_r(i,j,k,m)
             end do
c---------------------------------------------------------------------
c               scale the last row immediately
c---------------------------------------------------------------------
             fac2      = 1.d0/lhs_x_dv_r(i,j,k1,3)
             do    m = 1, 3
                rhs_dv_r(i,j,k1,m) = fac2*rhs_dv_r(i,j,k1,m)
             end do

        endif

        end subroutine krnl_z_solve_forward2







       attributes(global) subroutine krnl_z_solve_forward3(
     >                    lhsp_x_dv_r,lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, kSize,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: kSize, iSize,subSections

        integer i,j,k, blockSize, ioff, m, k1, k2

        double precision fac1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        j = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff


        !if(i .le. iSize) then
        do k = 0, kSize-1
             k1 = k  + 1
             k2 = k  + 2

             if(i .le. iSize) then

                m = 4
                fac1       = 1.d0/lhsp_x_dv_r(i,j,k,3)
                lhsp_x_dv_r(i,j,k,4)  = fac1*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j,k,5)  = fac1*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k,m)  = fac1*rhs_dv_r(i,j,k,m)
                lhsp_x_dv_r(i,j,k1,3) = lhsp_x_dv_r(i,j,k1,3) -
     >                       lhsp_x_dv_r(i,j,k1,2)*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j,k1,4) = lhsp_x_dv_r(i,j,k1,4) -
     >                       lhsp_x_dv_r(i,j,k1,2)*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k1,m) = rhs_dv_r(i,j,k1,m) -
     >                       lhsp_x_dv_r(i,j,k1,2)*rhs_dv_r(i,j,k,m)
                lhsp_x_dv_r(i,j,k2,2) = lhsp_x_dv_r(i,j,k2,2) -
     >                       lhsp_x_dv_r(i,j,k2,1)*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j,k2,3) = lhsp_x_dv_r(i,j,k2,3) -
     >                       lhsp_x_dv_r(i,j,k2,1)*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k2,m) = rhs_dv_r(i,j,k2,m) -
     >                       lhsp_x_dv_r(i,j,k2,1)*rhs_dv_r(i,j,k,m)
                m = 5
                fac1       = 1.d0/lhsm_x_dv_r(i,j,k,3)
                lhsm_x_dv_r(i,j,k,4)  = fac1*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j,k,5)  = fac1*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k,m)  = fac1*rhs_dv_r(i,j,k,m)
                lhsm_x_dv_r(i,j,k1,3) = lhsm_x_dv_r(i,j,k1,3) -
     >                       lhsm_x_dv_r(i,j,k1,2)*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j,k1,4) = lhsm_x_dv_r(i,j,k1,4) -
     >                       lhsm_x_dv_r(i,j,k1,2)*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k1,m) = rhs_dv_r(i,j,k1,m) -
     >                       lhsm_x_dv_r(i,j,k1,2)*rhs_dv_r(i,j,k,m)
                lhsm_x_dv_r(i,j,k2,2) = lhsm_x_dv_r(i,j,k2,2) -
     >                       lhsm_x_dv_r(i,j,k2,1)*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j,k2,3) = lhsm_x_dv_r(i,j,k2,3) -
     >                       lhsm_x_dv_r(i,j,k2,1)*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k2,m) = rhs_dv_r(i,j,k2,m) -
     >                       lhsm_x_dv_r(i,j,k2,1)*rhs_dv_r(i,j,k,m)
              
             endif

        end do
        !endif


        end subroutine krnl_z_solve_forward3







        attributes(global) subroutine krnl_z_solve_forward4(
     >                    lhsp_x_dv_r, lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, kstart,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: kstart, iSize, subSections

        integer i,j,k, blockSize, ioff, m, k1, k2

        double precision fac1, fac2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        j = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff

       k = kstart
       k1 = k + 1      

       if(i .le. iSize) then
             m = 4
             fac1       = 1.d0/lhsp_x_dv_r(i,j,k,3)
             lhsp_x_dv_r(i,j,k,4)  = fac1*lhsp_x_dv_r(i,j,k,4)
             lhsp_x_dv_r(i,j,k,5)  = fac1*lhsp_x_dv_r(i,j,k,5)
             rhs_dv_r(i,j,k,m)  = fac1*rhs_dv_r(i,j,k,m)
             lhsp_x_dv_r(i,j,k1,3) = lhsp_x_dv_r(i,j,k1,3) -
     >                    lhsp_x_dv_r(i,j,k1,2)*lhsp_x_dv_r(i,j,k,4)
             lhsp_x_dv_r(i,j,k1,4) = lhsp_x_dv_r(i,j,k1,4) -
     >                    lhsp_x_dv_r(i,j,k1,2)*lhsp_x_dv_r(i,j,k,5)
             rhs_dv_r(i,j,k1,m) = rhs_dv_r(i,j,k1,m) -
     >                    lhsp_x_dv_r(i,j,k1,2)*rhs_dv_r(i,j,k,m)
             m = 5
             fac1       = 1.d0/lhsm_x_dv_r(i,j,k,3)
             lhsm_x_dv_r(i,j,k,4)  = fac1*lhsm_x_dv_r(i,j,k,4)
             lhsm_x_dv_r(i,j,k,5)  = fac1*lhsm_x_dv_r(i,j,k,5)
             rhs_dv_r(i,j,k,m)  = fac1*rhs_dv_r(i,j,k,m)
             lhsm_x_dv_r(i,j,k1,3) = lhsm_x_dv_r(i,j,k1,3) -
     >                    lhsm_x_dv_r(i,j,k1,2)*lhsm_x_dv_r(i,j,k,4)
             lhsm_x_dv_r(i,j,k1,4) = lhsm_x_dv_r(i,j,k1,4) -
     >                    lhsm_x_dv_r(i,j,k1,2)*lhsm_x_dv_r(i,j,k,5)
             rhs_dv_r(i,j,k1,m) = rhs_dv_r(i,j,k1,m) -
     >                    lhsm_x_dv_r(i,j,k1,2)*rhs_dv_r(i,j,k,m)
c---------------------------------------------------------------------
c               Scale the last row immediately (some of this is overkill
c               if this is the last cell)
c---------------------------------------------------------------------
          rhs_dv_r(i,j,k1,4) = rhs_dv_r(i,j,k1,4)/lhsp_x_dv_r(i,j,k1,3)
          rhs_dv_r(i,j,k1,5) = rhs_dv_r(i,j,k1,5)/lhsm_x_dv_r(i,j,k1,3)

       endif


       end subroutine krnl_z_solve_forward4






         attributes(global) subroutine krnl_z_solve_back1(
     >                  lhs_x_dv_r,lhsp_x_dv_r, lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, kstart,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: kstart, iSize, subSections

        integer i,j,k, blockSize, ioff, m, k1, k2

        double precision fac1, fac2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        j = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff

        k = kstart
        k1 = k + 1    

        if(i .le. iSize) then

             do   m = 1, 3
                rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) -
     >              lhs_x_dv_r(i,j,k,4)*rhs_dv_r(i,j,k1, m)
             end do

             rhs_dv_r(i,j,k,4) = rhs_dv_r(i,j,k,4) -
     >            lhsp_x_dv_r(i,j,k,4)*rhs_dv_r(i,j,k1,4)
             rhs_dv_r(i,j,k,5) = rhs_dv_r(i,j,k,5) -
     >             lhsm_x_dv_r(i,j,k,4)*rhs_dv_r(i,j,k1,5)
        endif


        end subroutine krnl_z_solve_back1







       attributes(global) subroutine krnl_z_solve_back2(
     >             lhs_x_dv_r,lhsp_x_dv_r, lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, kSize,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: kSize, iSize, subSections

        integer i,j,k, blockSize, ioff, m, k1, k2

        double precision fac1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        j = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff

      
        do k = kSize-1, 0, -1
             k1 = k  + 1
             k2 = k  + 2

             if(i .le. iSize) then
                do   m = 1, 3
                   rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - 
     >                   lhs_x_dv_r(i,j,k,4)*rhs_dv_r(i,j,k1,m) -
     >                   lhs_x_dv_r(i,j,k,5)*rhs_dv_r(i,j,k2,m)
                end do

c---------------------------------------------------------------------
c      And the remaining two
c---------------------------------------------------------------------
                rhs_dv_r(i,j,k,4) = rhs_dv_r(i,j,k,4) - 
     >                   lhsp_x_dv_r(i,j,k,4)*rhs_dv_r(i,j,k1,4) -
     >                   lhsp_x_dv_r(i,j,k,5)*rhs_dv_r(i,j,k2,4)
                rhs_dv_r(i,j,k,5) = rhs_dv_r(i,j,k,5) - 
     >                   lhsm_x_dv_r(i,j,k,4)*rhs_dv_r(i,j,k1,5) -
     >                   lhsm_x_dv_r(i,j,k,5)*rhs_dv_r(i,j,k2,5)
            
             endif

        end do


        end subroutine krnl_z_solve_back2










        end module z_solve_kernels
