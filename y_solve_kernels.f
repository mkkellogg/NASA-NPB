         module y_solve_kernels
         implicit none
         contains

         attributes(device) subroutine shift_or_resety(a,hi,lo,i)
             integer a
             integer,value :: hi,lo,i

             if(hi .ge .lo) then
               if(a .eq. hi) then
                 a = lo
               else
                 a = a + i
               endif
             else
               if(a .eq. hi) then
                 a = lo
               else
                 a = a - i
               endif
             endif

         end subroutine shift_or_resety





         attributes(global) subroutine krnl_y_solve_init(
     >                       lhs_x_dv_j_r,
     >                  rho_i_dv_j_r, vs_dv_j_r,  imax, jmax, kmax, 
     >                       imaxp, jmaxp,  dtty1, dtty2, c2dtty1,
     >                  con43, c1c5, c3c4, dy1, dy3, dy5,dymax,
     >                      upperj1,upperj2, subSections)


        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision,dimension(0:IMAXP,0:IMAXP,
     >                                  0:KMAX-1,5)::lhs_x_dv_j_r
        double precision,dimension(0:JMAXP,0:IMAXP,
     >                                  0:KMAX-1)::rho_i_dv_j_r
        double precision,dimension(0:JMAXP,0:IMAXP,0:KMAX-1)::vs_dv_j_r

        double precision,value :: dtty1,dtty2,c2dtty1, con43, c1c5,
     >                           c3c4, dy1, dy3, dy5, dymax
        integer, value ::  upperj1, upperj2, subSections

        integer i,j,k, modF, blockSize, ioff, subISize

        double precision tempd, ru1, o1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        subISize = blockSize - 2
        k = blockIDx%x
        i = (blockIDx%y-1)/subSections +1
        ioff = (mod(blockIDx%y-1,subSections)*subISize)
        j = threadIDx%x + ioff -1

             if(j .le. upperj1) then

                ru1 = c3c4*rho_i_dv_j_r(j,i,k)
                temp(threadIDx%x) = vs_dv_j_r(j,i,k)
                temp(threadIDx%x+blockSize) = dmax1(dy3+con43*ru1,
     >                          dy5+c1c5*ru1,
     >                          dymax+ru1,
     >                          dy1)
             endif

             call syncthreads()

             if(j .le. upperj2 .and. threadIDx%x .gt. 1
     >                         .and. threadIDx%x .lt. blockSize) then

                lhs_x_dv_j_r(j,i,k,1) =   0.0d0
                lhs_x_dv_j_r(j,i,k,2) = - dtty2 * temp(threadIDx%x-1) -
     >                          dtty1 * temp(threadIDx%x+blockSize-1)
                lhs_x_dv_j_r(j,i,k,3) =   1.0d0 + c2dtty1 *
     >                          temp(threadIDx%x+blockSize)
                lhs_x_dv_j_r(j,i,k,4) =   dtty2 * temp(threadIDx%x+1) -
     >                          dtty1 * temp(threadIDx%x+blockSize+1)
                lhs_x_dv_j_r(j,i,k,5) =   0.0d0

             endif


         end subroutine krnl_y_solve_init










        attributes(global) subroutine krnl_y_solve_dis1(
     >                    lhs_x_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp,  
     >                    comz1,comz4,comz5,comz6, jstart,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r

        double precision,value :: comz1,comz4,comz5,comz6
        integer, value :: jstart, subISize, subSections, iSize

        integer i,j,k, blockSize, ioff

        double precision m

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections + 1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff


        if(jstart .le. 0) then
           j = 1
        else 
           j = jstart
        endif
        

        if(i .le. iSize) then


         if(jstart .le. 0) then

             j = 1
             lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz5
             lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4
             lhs_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5) + comz1

             lhs_x_dv_r(i,j+1,k,2) = lhs_x_dv_r(i,j+1,k,2) - comz4
             lhs_x_dv_r(i,j+1,k,3) = lhs_x_dv_r(i,j+1,k,3) + comz6
             lhs_x_dv_r(i,j+1,k,4) = lhs_x_dv_r(i,j+1,k,4) - comz4
             lhs_x_dv_r(i,j+1,k,5) = lhs_x_dv_r(i,j+1,k,5) + comz1


         else

             j=jstart
             lhs_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1) + comz1
             lhs_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - comz4
             lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz6
             lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4

             lhs_x_dv_r(i,j+1,k,1) = lhs_x_dv_r(i,j+1,k,1) + comz1
             lhs_x_dv_r(i,j+1,k,2) = lhs_x_dv_r(i,j+1,k,2) - comz4
             lhs_x_dv_r(i,j+1,k,3) = lhs_x_dv_r(i,j+1,k,3) + comz5


         endif

        endif

        end subroutine krnl_y_solve_dis1















        attributes(global) subroutine krnl_y_solve_dis2(
     >                    lhs_x_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp,  
     >                    comz1,comz4,comz5,comz6, 
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r

        double precision,value :: comz1,comz4,comz5,comz6
        integer, value :: subISize, subSections, iSize

        integer i,j,k, blockSize, joff

        double precision m

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +3
        joff = (mod(blockIDx%y-1,subSections)*blockSize)
        i = threadIDx%x + joff 

        if (i .le. iSize) then
        
                lhs_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1) + comz1
                lhs_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - comz4
                lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz6
                lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4
                lhs_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5) + comz1


        endif


        end subroutine krnl_y_solve_dis2
 











        attributes(global) subroutine krnl_y_solve_postdis(
     >           lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r, speed_dv,
     >                    imax, jmax, kmax, imaxp, jmaxp, dtty2,  
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1) :: speed_dv


        double precision,value ::dtty2
        integer, value :: iSize, subSections

        integer i,j,k,  blockSize, joff

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections+1 
        joff = (mod(blockIDx%y-1,subSections)*blockSize)
        i = threadIDx%x + joff



        if(i .le. iSize) then

                lhsp_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1)
                lhsp_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) -
     >                            dtty2 * speed_dv(i,j-1,k)
                lhsp_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3)
                lhsp_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) +
     >                            dtty2 * speed_dv(i,j+1,k)
                lhsp_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5)
                lhsm_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1)
                lhsm_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) +
     >                            dtty2 * speed_dv(i,j-1,k)
                lhsm_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3)
                lhsm_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) -
     >                            dtty2 * speed_dv(i,j+1,k)
                lhsm_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5)


        endif


        end subroutine krnl_y_solve_postdis









        attributes(global) subroutine krnl_y_solve_forward1(
     >                    lhs_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jSize,  
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision,dimension(0:IMAXP,0:IMAXP,
     >                                     0:KMAX-1,5)::lhs_x_dv_r
        double precision,dimension(0:IMAXP,0:JMAXP,
     >                                     0:KMAX-1,5)::rhs_dv_r

        integer, value :: jSize, iSize, subSections

        integer i,j,k, blockSize, ioff,m, j1, j2
        
        double precision fac1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff


        do j = 0, jSize-1
             j1 = j  + 1
             j2 = j  + 2

             if(i .le. iSize) then
               
                fac1      = 1.d0/lhs_x_dv_r(i,j,k,3)
                lhs_x_dv_r(i,j,k,4)  = fac1*lhs_x_dv_r(i,j,k,4)
                lhs_x_dv_r(i,j,k,5)  = fac1*lhs_x_dv_r(i,j,k,5)
                do    m = 1, 3
                   rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
                end do
                lhs_x_dv_r(i,j1,k,3) = lhs_x_dv_r(i,j1,k,3) -
     >                         lhs_x_dv_r(i,j1,k,2)*lhs_x_dv_r(i,j,k,4)
                lhs_x_dv_r(i,j1,k,4) = lhs_x_dv_r(i,j1,k,4) -
     >                         lhs_x_dv_r(i,j1,k,2)*lhs_x_dv_r(i,j,k,5)
                do    m = 1, 3
                   rhs_dv_r(i,j1,k,m) = rhs_dv_r(i,j1,k,m) -
     >                         lhs_x_dv_r(i,j1,k,2)*rhs_dv_r(i,j,k,m)
                end do
                lhs_x_dv_r(i,j2,k,2) = lhs_x_dv_r(i,j2,k,2) -
     >                         lhs_x_dv_r(i,j2,k,1)*lhs_x_dv_r(i,j,k,4)
                lhs_x_dv_r(i,j2,k,3) = lhs_x_dv_r(i,j2,k,3) -
     >                         lhs_x_dv_r(i,j2,k,1)*lhs_x_dv_r(i,j,k,5)
                do    m = 1, 3
                   rhs_dv_r(i,j2,k,m) = rhs_dv_r(i,j2,k,m) -
     >                         lhs_x_dv_r(i,j2,k,1)*rhs_dv_r(i,j,k,m)
                end do


             endif
              
        end do


        end subroutine krnl_y_solve_forward1









        attributes(global) subroutine krnl_y_solve_forward2(
     >                    lhs_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jstart,  
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhs_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: jstart, iSize, subSections

        integer i,j,k, i_idx, modF, blockSize, ioff, m, j1

        double precision fac1, fac2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff


        if(i .le. iSize) then 

          j  = jstart
          j1 = jstart+1
      
             fac1      = 1.d0/lhs_x_dv_r(i,j,k,3)
             lhs_x_dv_r(i,j,k,4)  = fac1*lhs_x_dv_r(i,j,k,4)
             lhs_x_dv_r(i,j,k,5)  = fac1*lhs_x_dv_r(i,j,k,4)
             do    m = 1, 3
                rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
             end do
             lhs_x_dv_r(i,j1,k,3) = lhs_x_dv_r(i,j1,k,3) -
     >                      lhs_x_dv_r(i,j1,k,2)*lhs_x_dv_r(i,j,k,4)
             lhs_x_dv_r(i,j1,k,4) = lhs_x_dv_r(i,j1,k,4) -
     >                      lhs_x_dv_r(i,j1,k,2)*lhs_x_dv_r(i,j,k,5)
             do    m = 1, 3
                rhs_dv_r(i,j1,k,m) = rhs_dv_r(i,j1,k,m) -
     >                      lhs_x_dv_r(i,j1,k,2)*rhs_dv_r(i,j,k,m)
             end do
c---------------------------------------------------------------------
c            scale the last row immediately
c---------------------------------------------------------------------
             fac2      = 1.d0/lhs_x_dv_r(i,j1,k,3)
             do    m = 1, 3
                rhs_dv_r(i,j1,k,m) = fac2*rhs_dv_r(i,j1,k,m)
             end do
               
        endif


        end subroutine krnl_y_solve_forward2










        attributes(global) subroutine krnl_y_solve_forward3(
     >                    lhsp_x_dv_r,lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jSize,  
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: jSize, subISize, subSections, iSize

        integer i,j,k, blockSize, ioff, m, j1, j2

        double precision fac1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff


        do j = 0, jSize-1
             j1 = j  + 1
             j2 = j  + 2

             if(i .le. iSize) then

                m = 4
                fac1       = 1.d0/lhsp_x_dv_r(i,j,k,3)
                lhsp_x_dv_r(i,j,k,4)  = fac1*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j,k,5)  = fac1*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
                lhsp_x_dv_r(i,j1,k,3) = lhsp_x_dv_r(i,j1,k,3) -
     >                       lhsp_x_dv_r(i,j1,k,2)*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j1,k,4) = lhsp_x_dv_r(i,j1,k,4) -
     >                       lhsp_x_dv_r(i,j1,k,2)*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j1,k,m) = rhs_dv_r(i,j1,k,m) -
     >                       lhsp_x_dv_r(i,j1,k,2)*rhs_dv_r(i,j,k,m)
                lhsp_x_dv_r(i,j2,k,2) = lhsp_x_dv_r(i,j2,k,2) -
     >                       lhsp_x_dv_r(i,j2,k,1)*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j2,k,3) = lhsp_x_dv_r(i,j2,k,3) -
     >                       lhsp_x_dv_r(i,j2,k,1)*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j2,k,m) = rhs_dv_r(i,j2,k,m) -
     >                       lhsp_x_dv_r(i,j2,k,1)*rhs_dv_r(i,j,k,m)
                m = 5
                fac1       = 1.d0/lhsm_x_dv_r(i,j,k,3)
                lhsm_x_dv_r(i,j,k,4)  = fac1*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j,k,5)  = fac1*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
                lhsm_x_dv_r(i,j1,k,3) = lhsm_x_dv_r(i,j1,k,3) -
     >                       lhsm_x_dv_r(i,j1,k,2)*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j1,k,4) = lhsm_x_dv_r(i,j1,k,4) -
     >                       lhsm_x_dv_r(i,j1,k,2)*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j1,k,m) = rhs_dv_r(i,j1,k,m) -
     >                       lhsm_x_dv_r(i,j1,k,2)*rhs_dv_r(i,j,k,m)
                lhsm_x_dv_r(i,j2,k,2) = lhsm_x_dv_r(i,j2,k,2) -
     >                       lhsm_x_dv_r(i,j2,k,1)*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j2,k,3) = lhsm_x_dv_r(i,j2,k,3) -
     >                       lhsm_x_dv_r(i,j2,k,1)*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j2,k,m) = rhs_dv_r(i,j2,k,m) -
     >                       lhsm_x_dv_r(i,j2,k,1)*rhs_dv_r(i,j,k,m)


             endif
              
        end do


        end subroutine krnl_y_solve_forward3









        attributes(global) subroutine krnl_y_solve_forward3B(
     >                    lhsp_x_dv_r,lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, kSize,jSize,
     >                    iSize, subSectionsW, subSectionsH)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: kSize, jSize, iSize, subSectionsW, 
     >                    subSectionsH

        integer i,j,k, blockSizeW,blockSizeH, ioff, m, j1, j2

        double precision fac1

        double precision, shared :: temp(*)

        blockSizeW = blockDim%x
        blockSizeH = blockDim%y

        i = ((blockIDx%x-1)*blockSizeW) + threadIDx%x
        k = ((blockIDx%y-1)*blockSizeH) + threadIDx%y


        do j = 0, jSize-1
             j1 = j  + 1
             j2 = j  + 2

             if(i .le. iSize .and. k .le. kSize) then

             m = 4
                fac1       = 1.d0/lhsp_x_dv_r(i,j,k,3)
                lhsp_x_dv_r(i,j,k,4)  = fac1*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j,k,5)  = fac1*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
                lhsp_x_dv_r(i,j1,k,3) = lhsp_x_dv_r(i,j1,k,3) -
     >                       lhsp_x_dv_r(i,j1,k,2)*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j1,k,4) = lhsp_x_dv_r(i,j1,k,4) -
     >                       lhsp_x_dv_r(i,j1,k,2)*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j1,k,m) = rhs_dv_r(i,j1,k,m) -
     >                       lhsp_x_dv_r(i,j1,k,2)*rhs_dv_r(i,j,k,m)
                lhsp_x_dv_r(i,j2,k,2) = lhsp_x_dv_r(i,j2,k,2) -
     >                       lhsp_x_dv_r(i,j2,k,1)*lhsp_x_dv_r(i,j,k,4)
                lhsp_x_dv_r(i,j2,k,3) = lhsp_x_dv_r(i,j2,k,3) -
     >                       lhsp_x_dv_r(i,j2,k,1)*lhsp_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j2,k,m) = rhs_dv_r(i,j2,k,m) -
     >                       lhsp_x_dv_r(i,j2,k,1)*rhs_dv_r(i,j,k,m)
                m = 5
                fac1       = 1.d0/lhsm_x_dv_r(i,j,k,3)
                lhsm_x_dv_r(i,j,k,4)  = fac1*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j,k,5)  = fac1*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
                lhsm_x_dv_r(i,j1,k,3) = lhsm_x_dv_r(i,j1,k,3) -
     >                       lhsm_x_dv_r(i,j1,k,2)*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j1,k,4) = lhsm_x_dv_r(i,j1,k,4) -
     >                       lhsm_x_dv_r(i,j1,k,2)*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j1,k,m) = rhs_dv_r(i,j1,k,m) -
     >                       lhsm_x_dv_r(i,j1,k,2)*rhs_dv_r(i,j,k,m)
                lhsm_x_dv_r(i,j2,k,2) = lhsm_x_dv_r(i,j2,k,2) -
     >                       lhsm_x_dv_r(i,j2,k,1)*lhsm_x_dv_r(i,j,k,4)
                lhsm_x_dv_r(i,j2,k,3) = lhsm_x_dv_r(i,j2,k,3) -
     >                       lhsm_x_dv_r(i,j2,k,1)*lhsm_x_dv_r(i,j,k,5)
                rhs_dv_r(i,j2,k,m) = rhs_dv_r(i,j2,k,m) -
     >                       lhsm_x_dv_r(i,j2,k,1)*rhs_dv_r(i,j,k,m)


             endif

        end do


        end subroutine krnl_y_solve_forward3B








        attributes(global) subroutine krnl_y_solve_forward4(
     >                    lhsp_x_dv_r, lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jstart,  
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsp_x_dv_r
        double precision, dimension(0:IMAXP,
     >                          0:IMAXP,0:KMAX-1,5) :: lhsm_x_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r


        integer, value :: jstart,  subSections, iSize

        integer i,j,k, blockSize, ioff, m, j1, j2

        double precision fac1, fac2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff

       
       j=jstart
       j1 = j + 1

       if(i .le. iSize) then


             m = 4
             fac1       = 1.d0/lhsp_x_dv_r(i,j,k,3)
             lhsp_x_dv_r(i,j,k,4)  = fac1*lhsp_x_dv_r(i,j,k,5)
             lhsp_x_dv_r(i,j,k,5)  = fac1*lhsp_x_dv_r(i,j,k,5)
             rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
             lhsp_x_dv_r(i,j1,k,3) = lhsp_x_dv_r(i,j1,k,3) -
     >                    lhsp_x_dv_r(i,j1,k,2)*lhsp_x_dv_r(i,j,k,4)
             lhsp_x_dv_r(i,j1,k,4) = lhsp_x_dv_r(i,j1,k,4) -
     >                    lhsp_x_dv_r(i,j1,k,2)*lhsp_x_dv_r(i,j,k,5)
             rhs_dv_r(i,j1,k,m)   = rhs_dv_r(i,j1,k,m) -
     >                    lhsp_x_dv_r(i,j1,k,2)*rhs_dv_r(i,j,k,m)
             m = 5
             fac1       = 1.d0/lhsm_x_dv_r(i,j,k,3)
             lhsm_x_dv_r(i,j,k,4)  = fac1*lhsm_x_dv_r(i,j,k,4)
             lhsm_x_dv_r(i,j,k,5)  = fac1*lhsm_x_dv_r(i,j,k,5)
             rhs_dv_r(i,j,k,m) = fac1*rhs_dv_r(i,j,k,m)
             lhsm_x_dv_r(i,j1,k,3) = lhsm_x_dv_r(i,j1,k,3) -
     >                    lhsm_x_dv_r(i,j1,k,2)*lhsm_x_dv_r(i,j,k,4)
             lhsm_x_dv_r(i,j1,k,4) = lhsm_x_dv_r(i,j1,k,4) -
     >                    lhsm_x_dv_r(i,j1,k,2)*lhsm_x_dv_r(i,j,k,5)
             rhs_dv_r(i,j1,k,m)   = rhs_dv_r(i,j1,k,m) -
     >                    lhsm_x_dv_r(i,j1,k,2)*rhs_dv_r(i,j,k,m)
c---------------------------------------------------------------------
c               Scale the last row immediately
c---------------------------------------------------------------------
             rhs_dv_r(i,j1,k,4)   = rhs_dv_r(i,j1,k,4)/
     >                                       lhsp_x_dv_r(i,j1,k,3)
             rhs_dv_r(i,j1,k,5)   = rhs_dv_r(i,j1,k,5)/
     >                                       lhsm_x_dv_r(i,j1,k,3)            

   
        endif

        end subroutine krnl_y_solve_forward4












        attributes(global) subroutine krnl_y_solve_back1(
     >             lhs_x_dv_r,lhsp_x_dv_r, lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jstart,  
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


        integer, value :: jstart, subSections, iSize

        integer i,j,k, blockSize, ioff, m, j1, j2

        double precision fac1, fac2

        double precision, shared :: temp(*)


        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff


        if(i .le. iSize) then

          j  = jstart
          j1 = j+1

             do   m = 1, 3
                rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) -
     >                           lhs_x_dv_r(i,j,k,4)*rhs_dv_r(i,j1,k,m)
             end do

             rhs_dv_r(i,j,k,4) = rhs_dv_r(i,j,k,4) -
     >                           lhsp_x_dv_r(i,j,k,4)*rhs_dv_r(i,j1,k,4)
             rhs_dv_r(i,j,k,5) = rhs_dv_r(i,j,k,5) -
     >                           lhsm_x_dv_r(i,j,k,4)*rhs_dv_r(i,j1,k,5)


               
        endif

        end subroutine krnl_y_solve_back1











        attributes(global) subroutine krnl_y_solve_back2(
     >        lhs_x_dv_r, lhsp_x_dv_r,lhsm_x_dv_r, rhs_dv_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, jSize,  
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


        integer, value :: jSize,  subSections, iSize

        integer i,j,k,  blockSize, ioff, m, j1, j2

        double precision fac1

        double precision, shared :: temp(*)


        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        ioff = (mod(blockIDx%x-1,subSections)*blockSize)
        i = threadIDx%x + ioff


        do j = jSize-1, 0, -1
             j1 = j  + 1
             j2 = j  + 2

             if(i .le. iSize) then

             j1 = j  + 1
             j2 = j  + 2
    
                do   m = 1, 3
                   rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) -
     >                      lhs_x_dv_r(i,j,k,4)*rhs_dv_r(i,j1,k,m) -
     >                      lhs_x_dv_r(i,j,k,5)*rhs_dv_r(i,j2,k,m)
                end do

c---------------------------------------------------------------------
c      And the remaining two
c---------------------------------------------------------------------
                rhs_dv_r(i,j,k,4) = rhs_dv_r(i,j,k,4) -
     >                      lhsp_x_dv_r(i,j,k,4)*rhs_dv_r(i,j1,k,4) -
     >                       lhsp_x_dv_r(i,j,k,5)*rhs_dv_r(i,j2,k,4)
                rhs_dv_r(i,j,k,5) = rhs_dv_r(i,j,k,5) -
     >                     lhsm_x_dv_r(i,j,k,4)*rhs_dv_r(i,j1,k,5) -
     >                        lhsm_x_dv_r(i,j,k,5)*rhs_dv_r(i,j2,k,5)
             

             endif

        end do


        end subroutine krnl_y_solve_back2





         end module y_solve_kernels

