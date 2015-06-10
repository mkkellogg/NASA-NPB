        
         module x_solve_kernels
         implicit none
         contains



        attributes(device) subroutine shift_or_reset(a,hi,lo,i)
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
 
        end subroutine shift_or_reset




        attributes(global) subroutine krnl_x_solve_init(lhs_x_dv_r,
     >               rho_i_dv, us_dv, imax, jmax, kmax, imaxp, jmaxp, 
     >              dttx1, dttx2, c2dttx1,con43, c1c5, c3c4, dx1, dx2,
     >              dx5,dxmax, upperi1, upperi2, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision,dimension(0:IMAXP,0:IMAXP,
     >                                  0:KMAX-1,5)::lhs_x_dv_r
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::rho_i_dv
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::us_dv

        double precision,value :: dttx1,dttx2,c2dttx1, con43, c1c5,
     >                           c3c4, dx1, dx2, dx5, dxmax
        integer, value ::  upperi1, upperi2, subSections

        integer i,j,k, modF, blockSize, ioff, subISize

        double precision tempd, ru1, o1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        subISize = blockSize - 2
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        ioff = (mod(blockIDx%y-1,subSections)*subISize)
        i = threadIDx%x + ioff -1

             if(i .le. upperi1) then
           
                ru1 = c3c4*rho_i_dv(i,j,k)
                temp(threadIDx%x) = us_dv(i,j,k)
                temp(threadIDx%x+blockSize) = dmax1(dx2+con43*ru1, 
     >                          dx5+c1c5*ru1,
     >                          dxmax+ru1,
     >                          dx1)
             endif

             call syncthreads()

             if(i .le. upperi2 .and. threadIDx%x .gt. 1 
     >                         .and. threadIDx%x .lt. blockSize) then

                lhs_x_dv_r(i,j,k,1) =   0.0d0
                lhs_x_dv_r(i,j,k,2) = - dttx2 * temp(threadIDx%x-1) - 
     >                          dttx1 * temp(threadIDx%x+blockSize-1)
                lhs_x_dv_r(i,j,k,3) =   1.0d0 + c2dttx1 * 
     >                          temp(threadIDx%x+blockSize)
                lhs_x_dv_r(i,j,k,4) =   dttx2 * temp(threadIDx%x+1) -
     >                          dttx1 * temp(threadIDx%x+blockSize+1)
                lhs_x_dv_r(i,j,k,5) =   0.0d0

             endif
 
       end subroutine krnl_x_solve_init 









        attributes(global) subroutine krnl_x_solve_dis1(lhs_x_dv_r,
     >                            imax, jmax, kmax, imaxp, jmaxp, 
     >             comz1,comz4,comz5,comz6, istart, jSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision,dimension(0:IMAXP,0:IMAXP,
     >                                0:KMAX-1,5)::lhs_x_dv_r

        double precision,value :: comz1,comz4,comz5,comz6
        integer, value :: istart, jSize, subSections

        integer i,j,k, modF, blockSize, joff

        double precision m1,m2,m3,m4, max4, tempd, ru1, o1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections + 1 
        joff = (mod(blockIDx%x-1,subSections)*blockSize)
        j = threadIDx%x + joff 

        if(j .le. jSize) then
          if(iStart .le. 0) then
             i = 1
             lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz5
             lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4
             lhs_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5) + comz1
  
             lhs_x_dv_r(i+1,j,k,2) = lhs_x_dv_r(i+1,j,k,2) - comz4
             lhs_x_dv_r(i+1,j,k,3) = lhs_x_dv_r(i+1,j,k,3) + comz6
             lhs_x_dv_r(i+1,j,k,4) = lhs_x_dv_r(i+1,j,k,4) - comz4
             lhs_x_dv_r(i+1,j,k,5) = lhs_x_dv_r(i+1,j,k,5) + comz1
          else
             i = istart
             lhs_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1) + comz1
             lhs_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - comz4
             lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz6
             lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4

             lhs_x_dv_r(i+1,j,k,1) = lhs_x_dv_r(i+1,j,k,1) + comz1
             lhs_x_dv_r(i+1,j,k,2) = lhs_x_dv_r(i+1,j,k,2) - comz4
             lhs_x_dv_r(i+1,j,k,3) = lhs_x_dv_r(i+1,j,k,3) + comz5

         endif
        endif

        end subroutine krnl_x_solve_dis1









        attributes(global) subroutine krnl_x_solve_dis2(
     >              imax, jmax, kmax, imaxp, jmaxp, lhs_x_dv_r,
     >                  comz1,comz4,comz5,comz6, iSize,
     >                  subSections)

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision, dimension(0:IMAXP,0:IMAXP,
     >                               0:KMAX-1,5)::lhs_x_dv_r

       double precision,value :: comz1,comz4,comz5,comz6
       integer, value ::  iSize, subSections

        integer i,j,k, modF, blockSize, ioff

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        ioff = (mod(blockIDx%y-1,subSections)*blockSize)
        i = threadIDx%x + ioff +2
        
 
        if (i .le. 3 + iSize-1 ) then                
                lhs_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1) + comz1
                lhs_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - comz4
                lhs_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3) + comz6
                lhs_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - comz4
                lhs_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5) + comz1
        endif

        
        end subroutine krnl_x_solve_dis2











       attributes(global) subroutine krnl_x_solve_postdis(
     >            lhs_x_dv_r, lhsp_x_dv_r, lhsm_x_dv_r,speed_dv,
     >                         imax, jmax, kmax, imaxp, jmaxp, 
     >                     dttx2,iSize,subSections) 

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhs_x_dv_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsp_x_dv_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsm_x_dv_r
       double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: speed_dv

       integer, value :: iSize, subSections
       double precision, value :: dttx2

        integer i,j,k, modF, blockSize, ioff, subISize

        double precision m1,t

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        subISize = blockSize - 2

        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1        
        ioff = (mod(blockIDx%y-1,subSections)*subISize)
        i = threadIDx%x + ioff -1


          if (threadIDx%x .gt. 1 .and. threadIDx%x .lt. blockSize
     >        .and. i .le. iSize) then
                lhsp_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1)
                lhsp_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) - 
     >                            dttx2 * speed_dv(i-1,j,k)
                lhsp_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3)
                lhsp_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) + 
     >                            dttx2 * speed_dv(i+1,j,k)
                lhsp_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5)
                lhsm_x_dv_r(i,j,k,1) = lhs_x_dv_r(i,j,k,1)
                lhsm_x_dv_r(i,j,k,2) = lhs_x_dv_r(i,j,k,2) + 
     >                            dttx2 * speed_dv(i-1,j,k)
                lhsm_x_dv_r(i,j,k,3) = lhs_x_dv_r(i,j,k,3)
                lhsm_x_dv_r(i,j,k,4) = lhs_x_dv_r(i,j,k,4) - 
     >                            dttx2 * speed_dv(i+1,j,k)
                lhsm_x_dv_r(i,j,k,5) = lhs_x_dv_r(i,j,k,5)

          endif

        

        end subroutine krnl_x_solve_postdis











       attributes(global) subroutine krnl_x_solve_thomas(
     >                    lhs_x_dv_j_r, rhs_dv_j_r,
     >                    imax, jmax, kmax, imaxp, jmaxp, iSize,
     >                    jSize, subSections)

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision,dimension(0:IMAXP,0:IMAXP,
     >                            0:KMAX-1,5)::lhs_x_dv_j_r
       double precision,dimension(0:JMAXP,0:IMAXP,
     >                            0:KMAX-1,5)::rhs_dv_j_r

       integer, value :: iSize, jSize, subSections

       integer i,j,k, blockSize, joff, i1,i2,m

       double precision fac1

       double precision, shared :: temp(*)

       blockSize = blockDim%x
       k = (blockIDx%x-1)/subSections +1
       joff = (mod(blockIDx%x-1,subSections)*blockSize)
       j = threadIDx%x + joff 


       do i = 0, iSize-1
             i1 = i  + 1
             i2 = i  + 2

             if(j .le. jSize) then

                i1 = i  + 1
                i2 = i  + 2
                fac1      = 1.d0/lhs_x_dv_j_r(j,i,k,3)
                lhs_x_dv_j_r(j,i,k,4)  = fac1*lhs_x_dv_j_r(j,i,k,4)
                lhs_x_dv_j_r(j,i,k,5)  = fac1*lhs_x_dv_j_r(j,i,k,5)
                do    m = 1, 3
                   rhs_dv_j_r(j,i,k,m) = fac1*rhs_dv_j_r(j,i,k,m)
                end do
                lhs_x_dv_j_r(j,i1,k,3) = lhs_x_dv_j_r(j,i1,k,3) -
     >                 lhs_x_dv_j_r(j,i1,k,2)*lhs_x_dv_j_r(j,i,k,4)
                lhs_x_dv_j_r(j,i1,k,4) = lhs_x_dv_j_r(j,i1,k,4) -
     >                   lhs_x_dv_j_r(j,i1,k,2)*lhs_x_dv_j_r(j,i,k,5)
                do    m = 1, 3
                   rhs_dv_j_r(j,i1,k,m) = rhs_dv_j_r(j,i1,k,m) -
     >           lhs_x_dv_j_r(j,i1,k,2)*rhs_dv_j_r(j,i,k,m)
                end do
                lhs_x_dv_j_r(j,i2,k,2) = lhs_x_dv_j_r(j,i2,k,2) -
     >                 lhs_x_dv_j_r(j,i2,k,1)*lhs_x_dv_j_r(j,i,k,4)
                lhs_x_dv_j_r(j,i2,k,3) = lhs_x_dv_j_r(j,i2,k,3) -
     >                 lhs_x_dv_j_r(j,i2,k,1)*lhs_x_dv_j_r(j,i,k,5)
                do    m = 1, 3
                   rhs_dv_j_r(j,i2,k,m) = rhs_dv_j_r(j,i2,k,m) -
     >                lhs_x_dv_j_r(j,i2,k,1)*rhs_dv_j_r(j,i,k,m)
                end do
             endif
           
        end do

        end subroutine krnl_x_solve_thomas


















       attributes(global) subroutine krnl_x_solve_pthomas1(
     >                      rhs_dv_j_r,lhs_x_dv_j_r, 
     >                      imax, jmax, kmax, imaxp, jmaxp, 
     >                      jSize, istart, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision,dimension(0:IMAXP,0:IMAXP,
     >                            0:KMAX-1,5)::lhs_x_dv_j_r
       double precision,dimension(0:JMAXP,0:IMAXP,
     >                            0:KMAX-1,5)::rhs_dv_j_r

        integer, value :: jSize, istart, subSections

        integer i,j,k, i1, i2,m, modF, joff,  blockSize

        double precision fac1, fac2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        joff = (mod(blockIDx%x-1,subSections)*blockSize)
        j = threadIDx%x + joff 

        if(j .le. jSize) then

             i  = istart
             i1 = istart+1
             fac1      = 1.d0/lhs_x_dv_j_r(j,i,k,3)
             lhs_x_dv_j_r(j,i,k,4)  = fac1*lhs_x_dv_j_r(j,i,k,4)
             lhs_x_dv_j_r(j,i,k,5)  = fac1*lhs_x_dv_j_r(j,i,k,5)
             do    m = 1, 3
                rhs_dv_j_r(j,i,k,m) = fac1*rhs_dv_j_r(j,i,k,m)
             end do
             lhs_x_dv_j_r(j,i1,k,3) = lhs_x_dv_j_r(j,i1,k,3) -
     >                      lhs_x_dv_j_r(j,i1,k,2)*lhs_x_dv_j_r(j,i,k,4)
             lhs_x_dv_j_r(j,i1,k,4) = lhs_x_dv_j_r(j,i1,k,4) -
     >                      lhs_x_dv_j_r(j,i1,k,2)*lhs_x_dv_j_r(j,i,k,5)
             do    m = 1, 3
                rhs_dv_j_r(j,i1,k,m) = rhs_dv_j_r(j,i1,k,m) -
     >                      lhs_x_dv_j_r(j,i1,k,2)*rhs_dv_j_r(j,i,k,m)
             end do
c---------------------------------------------------------------------
c            scale the last row immediately 
c---------------------------------------------------------------------
             fac2             = 1.d0/lhs_x_dv_j_r(j,i1,k,3)
             do    m = 1, 3
                rhs_dv_j_r(j,i1,k,m) = fac2*rhs_dv_j_r(j,i1,k,m)
             end do

        endif

       end subroutine krnl_x_solve_pthomas1










       attributes(global) subroutine krnl_x_solve_pthomas2(
     >                 lhsp_x_dv_j_r,lhsm_x_dv_j_r, rhs_dv_j_r,
     >                      imax, jmax, kmax, imaxp, jmaxp, 
     >                      iSize, jSize, subSections)

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsp_x_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsm_x_dv_j_r
       double precision,dimension(0:JMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::rhs_dv_j_r

       integer, value :: jSize, iSize, subSections

       integer i,j,k, i1, i2,m, joff, blockSize

       double precision fac1

       double precision, shared :: temp(*)

       blockSize = blockDim%x
       k = (blockIDx%x-1)/subSections +1
       joff = (mod(blockIDx%x-1,subSections)*blockSize)
       j = threadIDx%x + joff 

       do i=0,iSize-1            

                i1 = i  + 1
                i2 = i  + 2

                if(j .le. jSize) then
                m = 4
                fac1       = 1.d0/lhsp_x_dv_j_r(j,i,k,3)
                lhsp_x_dv_j_r(j,i,k,4)  = fac1*lhsp_x_dv_j_r(j,i,k,4)
                lhsp_x_dv_j_r(j,i,k,5)  = fac1*lhsp_x_dv_j_r(j,i,k,5)
                rhs_dv_j_r(j,i,k,m) = fac1*rhs_dv_j_r(j,i,k,m)
                lhsp_x_dv_j_r(j,i1,k,3) = lhsp_x_dv_j_r(j,i1,k,3) -
     >                        lhsp_x_dv_j_r(j,i1,k,2)*lhsp_x_dv_j_r(j,i,
     >                        k,4)
                lhsp_x_dv_j_r(j,i1,k,4) = lhsp_x_dv_j_r(j,i1,k,4) -
     >                        lhsp_x_dv_j_r(j,i1,k,2)*lhsp_x_dv_j_r(j,i,
     >                        k,5)
                rhs_dv_j_r(j,i1,k,m) = rhs_dv_j_r(j,i1,k,m) -
     >                        lhsp_x_dv_j_r(j,i1,k,2)*rhs_dv_j_r(j,i,k,
     >                        m)
                lhsp_x_dv_j_r(j,i2,k,2) = lhsp_x_dv_j_r(j,i2,k,2) -
     >                        lhsp_x_dv_j_r(j,i2,k,1)*lhsp_x_dv_j_r(j,i,
     >                        k,4)
                lhsp_x_dv_j_r(j,i2,k,3) = lhsp_x_dv_j_r(j,i2,k,3) -
     >                        lhsp_x_dv_j_r(j,i2,k,1)*lhsp_x_dv_j_r(j,i,
     >                        k,5)
                rhs_dv_j_r(j,i2,k,m) = rhs_dv_j_r(j,i2,k,m) -
     >                        lhsp_x_dv_j_r(j,i2,k,1)*rhs_dv_j_r(j,i,k,
     >                        m)
                m = 5
                fac1       = 1.d0/lhsm_x_dv_j_r(j,i,k,3)
                lhsm_x_dv_j_r(j,i,k,4)  = fac1*lhsm_x_dv_j_r(j,i,k,4)
                lhsm_x_dv_j_r(j,i,k,5)  = fac1*lhsm_x_dv_j_r(j,i,k,5)
                rhs_dv_j_r(j,i,k,m) = fac1*rhs_dv_j_r(j,i,k,m)
                lhsm_x_dv_j_r(j,i1,k,3) = lhsm_x_dv_j_r(j,i1,k,3) -
     >                        lhsm_x_dv_j_r(j,i1,k,2)*lhsm_x_dv_j_r(j,i,
     >                        k,4)
                lhsm_x_dv_j_r(j,i1,k,4) = lhsm_x_dv_j_r(j,i1,k,4) -
     >                        lhsm_x_dv_j_r(j,i1,k,2)*lhsm_x_dv_j_r(j,i,
     >                        k,5)
                rhs_dv_j_r(j,i1,k,m) = rhs_dv_j_r(j,i1,k,m) -
     >                        lhsm_x_dv_j_r(j,i1,k,2)*rhs_dv_j_r(j,i,k,
     >                        m)
                lhsm_x_dv_j_r(j,i2,k,2) = lhsm_x_dv_j_r(j,i2,k,2) -
     >                        lhsm_x_dv_j_r(j,i2,k,1)*lhsm_x_dv_j_r(j,i,
     >                        k,4)
                lhsm_x_dv_j_r(j,i2,k,3) = lhsm_x_dv_j_r(j,i2,k,3) -
     >                        lhsm_x_dv_j_r(j,i2,k,1)*lhsm_x_dv_j_r(j,i,
     >                        k,5)
                rhs_dv_j_r(j,i2,k,m) = rhs_dv_j_r(j,i2,k,m) -
     >                        lhsm_x_dv_j_r(j,i2,k,1)*rhs_dv_j_r(j,i,k,
     >                        m)
                endif
        end do


       end subroutine krnl_x_solve_pthomas2












       attributes(global) subroutine krnl_xsolve_pthomas3(
     >                      rhs_dv_j_r,lhsm_x_dv_j_r,lhsp_x_dv_j_r,
     >                      imax, jmax, kmax, imaxp, jmaxp, 
     >                      jSize, istart, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision,dimension(0:JMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::rhs_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsp_x_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsm_x_dv_j_r


        integer, value :: jSize, istart, subSections

        integer i,j,k, i1, i2,m, joff, blockSize

        double precision fac1, fac2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        joff = (mod(blockIDx%x-1,subSections)*blockSize)
        j = threadIDx%x + joff 


        if(j .le. jSize) then
              i  = istart
             i1 = istart+1
             m = 4
             fac1       = 1.d0/lhsp_x_dv_j_r(j,i,k,3)
             lhsp_x_dv_j_r(j,i,k,4)  = fac1*lhsp_x_dv_j_r(j,i,k,4)
             lhsp_x_dv_j_r(j,i,k,5)  = fac1*lhsp_x_dv_j_r(j,i,k,5)
             rhs_dv_j_r(j,i,k,m) = fac1*rhs_dv_j_r(j,i,k,m)
             lhsp_x_dv_j_r(j,i1,k,3) = lhsp_x_dv_j_r(j,i1,k,3) -
     >                      lhsp_x_dv_j_r(j,i1,k,2)*lhsp_x_dv_j_r(j,i,k,
     >                      4)
             lhsp_x_dv_j_r(j,i1,k,4) = lhsp_x_dv_j_r(j,i1,k,4) -
     >                      lhsp_x_dv_j_r(j,i1,k,2)*lhsp_x_dv_j_r(j,i,k,
     >                      5)
             rhs_dv_j_r(j,i1,k,m) = rhs_dv_j_r(j,i1,k,m) -
     >                      lhsp_x_dv_j_r(j,i1,k,2)*rhs_dv_j_r(j,i,k,m)
             m = 5
             fac1       = 1.d0/lhsm_x_dv_j_r(j,i,k,3)
             lhsm_x_dv_j_r(j,i,k,4)  = fac1*lhsm_x_dv_j_r(j,i,k,4)
             lhsm_x_dv_j_r(j,i,k,5)  = fac1*lhsm_x_dv_j_r(j,i,k,5)
             rhs_dv_j_r(j,i,k,m) = fac1*rhs_dv_j_r(j,i,k,m)
             lhsm_x_dv_j_r(j,i1,k,3) = lhsm_x_dv_j_r(j,i1,k,3) -
     >                      lhsm_x_dv_j_r(j,i1,k,2)*lhsm_x_dv_j_r(j,i,k,
     >                      4)
             lhsm_x_dv_j_r(j,i1,k,4) = lhsm_x_dv_j_r(j,i1,k,4) -
     >                      lhsm_x_dv_j_r(j,i1,k,2)*lhsm_x_dv_j_r(j,i,k,
     >                      5)
             rhs_dv_j_r(j,i1,k,m) = rhs_dv_j_r(j,i1,k,m) -
     >                      lhsm_x_dv_j_r(j,i1,k,2)*rhs_dv_j_r(j,i,k,m)
c---------------------------------------------------------------------
c               Scale the last row immediately
c---------------------------------------------------------------------
             rhs_dv_j_r(j,i1,k,4) = rhs_dv_j_r(j,i1,k,4)/lhsp_x_dv_j_r(
     >       j,i1,k,3)
             rhs_dv_j_r(j,i1,k,5) = rhs_dv_j_r(j,i1,k,5)/lhsm_x_dv_j_r(
     >       j,i1,k,3)
       
             endif

       end subroutine krnl_xsolve_pthomas3












       attributes(global) subroutine krnl_xsolve_bfill1(
     >       rhs_dv_j_r,lhs_x_dv_j_r, lhsm_x_dv_j_r,lhsp_x_dv_j_r,
     >                      imax, jmax, kmax, imaxp, jmaxp, 
     >                      jSize, istart, subSections)

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision,dimension(0:JMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::rhs_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhs_x_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsp_x_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsm_x_dv_j_r

        integer, value :: jSize, istart, subSections

        integer i,j,k, i1, i2,m, joff, blockSize

        double precision fac1, fac2

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        joff = (mod(blockIDx%x-1,subSections)*blockSize)
        j = threadIDx%x + joff 

        if(j .le. jSize) then
             i  = istart
             i1 = istart+1
             do   m = 1, 3
                rhs_dv_j_r(j,i,k,m) = rhs_dv_j_r(j,i,k,m) -
     >                             lhs_x_dv_j_r(j,i,k,4)*rhs_dv_j_r(j,
     >                             i1,k,m)
             end do

             rhs_dv_j_r(j,i,k,4) = rhs_dv_j_r(j,i,k,4) -
     >                          lhsp_x_dv_j_r(j,i,k,4)*rhs_dv_j_r(j,i1,
     >                          k,4)
             rhs_dv_j_r(j,i,k,5) = rhs_dv_j_r(j,i,k,5) -
     >                          lhsm_x_dv_j_r(j,i,k,4)*rhs_dv_j_r(j,i1,
     >                          k,5)
        endif
       
       end subroutine krnl_xsolve_bfill1
















       attributes(global) subroutine krnl_x_solve_bfill2(
     >          rhs_dv_j_r,lhs_x_dv_j_r,lhsp_x_dv_j_r,lhsm_x_dv_j_r, 
     >                      imax, jmax, kmax, imaxp, jmaxp, 
     >                      iSize, jSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision,dimension(0:JMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::rhs_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhs_x_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsp_x_dv_j_r
       double precision,dimension(0:IMAXP,0:IMAXP,
     >                                    0:KMAX-1,5)::lhsm_x_dv_j_r

        integer, value :: iSize, jSize, subSections

        integer i,j,k, i1, i2,m, joff, blockSize

        double precision fac1

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = (blockIDx%x-1)/subSections +1
        joff = (mod(blockIDx%x-1,subSections)*blockSize)
        j = threadIDx%x + joff 
       
       
        do i=iSize -1,0,-1

            if (j .le. jSize) then
                i1 = i  + 1
                i2 = i  + 2
                do   m = 1, 3
                   rhs_dv_j_r(j,i,k,m) = rhs_dv_j_r(j,i,k,m) - 
     >                     lhs_x_dv_j_r(j,i,k,4)*rhs_dv_j_r(j,i1,k,m) -
     >                     lhs_x_dv_j_r(j,i,k,5)*rhs_dv_j_r(j,i2,k,m)
                end do

c---------------------------------------------------------------------
c      And the remaining two
c---------------------------------------------------------------------
                rhs_dv_j_r(j,i,k,4) = rhs_dv_j_r(j,i,k,4) - 
     >                          lhsp_x_dv_j_r(j,i,k,4)*rhs_dv_j_r(j,i1,
     >                          k,4) -
     >                          lhsp_x_dv_j_r(j,i,k,5)*rhs_dv_j_r(j,i2,
     >                          k,4)
                rhs_dv_j_r(j,i,k,5) = rhs_dv_j_r(j,i,k,5) - 
     >                          lhsm_x_dv_j_r(j,i,k,4)*rhs_dv_j_r(j,i1,
     >                          k,5) -
     >                          lhsm_x_dv_j_r(j,i,k,5)*rhs_dv_j_r(j,i2,
     >                          k,5)
            endif
        end do

       end subroutine krnl_x_solve_bfill2








       attributes(global) subroutine krnl_x_solve_bfill2_nc(
     >                    lhs_x_dv, lhsp_x_dv,lhsm_x_dv, rhs_dv,
     >                    imax, jmax, kmax, imaxp, jmaxp, iSize,
     >                    subJSize, subSections, jSize)

       integer, value :: imax, jmax, kmax, imaxp, jmaxp

       double precision,dimension(5,0:IMAXP,0:IMAXP,0:KMAX-1)::lhs_x_dv
       double precision,dimension(5,0:IMAXP,0:IMAXP,0:KMAX-1)::lhsp_x_dv
       double precision,dimension(5,0:IMAXP,0:IMAXP,0:KMAX-1)::lhsm_x_dv
       double precision,dimension(5,0:IMAXP,0:JMAXP,0:KMAX-1) :: rhs_dv

       integer, value :: iSize, subJSize, subSections, jSize

       integer i,j,k, j_idx, doff, blockSize, koff, j_idx_base,
     >             localj, local_j_idx, basej, m, i1, i2
       integer t1,t2,t3,t4,t5,t6, t7, t8, t9

       double precision fac1

       double precision, shared :: temp(*)

       k = (blockIDx%x-1)/subSections +1
       koff = mod((blockIDx%x-1),subSections)
       blockSize = blockDim%x

       j_idx_base = (koff * blockSize)
       basej = (j_idx_base /5) +1
       j_idx = j_idx_base +  threadIDx%x -1
       j = (j_idx / 5) +1
       localj = (threadIDx%x -1)/5 +1
       local_j_idx = threadIDx%x-1
       doff = mod(j_idx,5)

       t1 =0
       t2 = blockSize
       t3 = blockSize * 2
       t4 = blockSize * 3
       t5 = blockSize * 4
       t6 = blockSize * 5

       i = iSize-1

       temp(threadIDx%x + t1) = lhs_x_dv(doff+1,i,j,k)
       temp(threadIDx%x + t2) = lhsp_x_dv(doff+1,i,j,k)
       temp(threadIDx%x + t3) = lhsm_x_dv(doff+1,i,j,k)
       temp(threadIDx%x + t4) = rhs_dv(doff+1,i,j,k)
       temp(threadIDx%x + t5) = rhs_dv(doff+1,i+1,j,k)
       temp(threadIDx%x + t6) = rhs_dv(doff+1,i+2,j,k)

       call syncthreads()

       do i = iSize-1, 0, -1
             i1 = i  + 1
             i2 = i  + 2

             if(doff .eq. 0 .and. j .le. jSize) then

                do   m = 1, 3
                   temp(t4 +((localj-1)*5)+ m) =
     >                        temp(t4 +((localj-1)*5)+ m) -
     >                       temp(t1 +((localj-1)*5)+ 4) *
     >                            temp(t5 +((localj-1)*5)+ m) -
     >                        temp(t1 +((localj-1)*5)+ 5) *
     >                            temp(t6 +((localj-1)*5)+ m)
                end do

c---------------------------------------------------------------------
c      And the remaining two
c---------------------------------------------------------------------
                temp(t4 +((localj-1)*5)+ 4) =
     >                         temp(t4 +((localj-1)*5)+ 4) -
     >                          temp(t2 +((localj-1)*5)+ 4) *
     >                               temp(t5 +((localj-1)*5)+ 4) -
     >                          temp(t2 +((localj-1)*5)+ 5) *
     >                               temp(t6 +((localj-1)*5)+ 4)
                temp(t4 +((localj-1)*5)+ 5) =
     >                          temp(t4 +((localj-1)*5)+ 5) -
     >                          temp(t3 +((localj-1)*5)+ 4) *
     >                              temp(t5 +((localj-1)*5)+ 5) -
     >                          temp(t3 +((localj-1)*5)+ 5) *
     >                              temp(t6 +((localj-1)*5)+ 5)


             endif

             call syncthreads()

             if(j .le. jSize) then
              rhs_dv(doff+1,i+2,j,k) = temp(threadIDx%x+t6)
             endif

              if(i .eq. 0 .and. j .le. jSize) then
                rhs_dv(doff+1,0,j,k) = temp(threadIDx%x+t4)
                rhs_dv(doff+1,1,j,k) = temp(threadIDx%x+t5)
              endif

              call syncthreads()

          call shift_or_reset(t4,blockSize*3,blockSize*5,blockSize)
          call shift_or_reset(t5,blockSize*3,blockSize*5,blockSize)
          call shift_or_reset(t6,blockSize*3,blockSize*5,blockSize)


              if(i .gt. 0) then
                 temp(threadIDx%x+t1) = lhs_x_dv(doff+1,i-1,j,k)
                 temp(threadIDx%x+t2) = lhsp_x_dv(doff+1,i-1,j,k)
                 temp(threadIDx%x+t3) = lhsm_x_dv(doff+1,i-1,j,k)
                 temp(threadIDx%x+t4) = rhs_dv(doff+1,i-1,j,k)
                 call syncthreads()
              endif


        end do


        end subroutine krnl_x_solve_bfill2_nc










        end module x_solve_kernels
