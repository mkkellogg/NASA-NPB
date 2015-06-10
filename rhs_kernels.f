
        module rhs_kernels
         implicit none
         contains






        attributes(global) subroutine krnl_rhs_init(u_dv_r, us_dv,vs_dv,
     >                 ws_dv,  qs_dv, rho_i_dv, speed_dv, square_dv,
     >       IMAX, JMAX, KMAX, c1c2, IMAXP, JMAXP, iSize, subSections)

        integer, value :: IMAX, JMAX, KMAX, IMAXP, JMAXP

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP, 0:JMAXP,0:KMAX-1)::rho_i_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: us_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: vs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: ws_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::square_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: qs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::speed_dv
      
        integer, value :: iSize, subSections

        double precision, value :: c1c2 

        double precision aux, rho_inv,sqr_tmp          
        integer :: i,j,k, blockSize

        double precision, shared :: temp(*) 

        blockSize = blockDim%x
        k = blockIDx%x-1
        j = (blockIDx%y-1)/subSections
        i = threadIDx%x -1 + (mod(blockIDx%y-1,subSections)*blockSize)

        if(i .le. iSize-1 ) then
          rho_inv = 1.0d0/u_dv_r(i,j,k,1)
          rho_i_dv(i,j,k) = rho_inv
          us_dv(i,j,k) = u_dv_r(i,j,k,2) * rho_inv
          vs_dv(i,j,k) = u_dv_r(i,j,k,3) * rho_inv
          ws_dv(i,j,k) = u_dv_r(i,j,k,4) * rho_inv
          square_dv(i,j,k)     = 0.5d0* (
     >                      u_dv_r(i,j,k,2)*u_dv_r(i,j,k,2) + 
     >                      u_dv_r(i,j,k,3)*u_dv_r(i,j,k,3) +
     >                      u_dv_r(i,j,k,4)*u_dv_r(i,j,k,4) ) * rho_inv
          qs_dv(i,j,k) = square_dv(i,j,k) * rho_inv
c---------------------------------------------------------------------
c               (don't need speed and ainx until the lhs computation)
c---------------------------------------------------------------------
          aux = c1c2*rho_inv* (u_dv_r(i,j,k,5) - square_dv(i,j,k))
          speed_dv(i,j,k) = dsqrt(aux)
        endif
     

       end subroutine krnl_rhs_init











       attributes(global) subroutine krnl_rhs_xflux(u_dv_r,us_dv, 
     >         vs_dv, ws_dv,  qs_dv, rho_i_dv, square_dv,rhs_dv_r,
     >           imax, jmax, kmax,  imaxp, jmaxp, dx1tx1, dx2tx1,
     >          dx3tx1, dx4tx1,dx5tx1, tx2, xxcon2, xxcon3, xxcon4, 
     >        xxcon5, con43, c2, c1, iSize, subSections) 

        integer, value :: IMAX, JMAX, KMAX, IMAXP, JMAXP

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::rho_i_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: us_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: vs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: ws_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::square_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: qs_dv
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dx1tx1,dx2tx1,dx3tx1,dx4tx1,dx5tx1, 
     >                tx2,xxcon2, xxcon3, xxcon4, xxcon5,con43, c2, c1  
       
        integer, value :: iSize, subSections

        integer :: i,j,k, blockSize 
        double precision uijk, up1, um1
        double precision, shared :: temp(*) 


        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections + 1
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)


        if(i .le. iSize) then
              uijk = us_dv(i,j,k)
                up1  = us_dv(i+1,j,k)
                um1  = us_dv(i-1,j,k)

                rhs_dv_r(i,j,k,1) = rhs_dv_r(i,j,k,1) + dx1tx1 *
     >                    (u_dv_r(i+1,j,k,1) - 2.0d0*u_dv_r(i,j,k,1) +
     >                     u_dv_r(i-1,j,k,1)) -
     >                    tx2 * (u_dv_r(i+1,j,k,2) - u_dv_r(i-1,j,k,2))

                rhs_dv_r(i,j,k,2) = rhs_dv_r(i,j,k,2) + dx2tx1 *
     >                    (u_dv_r(i+1,j,k,2) - 2.0d0*u_dv_r(i,j,k,2) +
     >                     u_dv_r(i-1,j,k,2)) +
     >                    xxcon2*con43 * (up1 - 2.0d0*uijk + um1) -
     >                    tx2 * (u_dv_r(i+1,j,k,2)*up1 -
     >                           u_dv_r(i-1,j,k,2)*um1 +
     >                         (u_dv_r(i+1,j,k,5)- square_dv(i+1,j,k)-
     >                          u_dv_r(i-1,j,k,5)+ square_dv(i-1,j,k))*
     >                            c2)

                rhs_dv_r(i,j,k,3) = rhs_dv_r(i,j,k,3) + dx3tx1 *
     >                    (u_dv_r(i+1,j,k,3) - 2.0d0*u_dv_r(i,j,k,3) +
     >                     u_dv_r(i-1,j,k,3)) +
     >                  xxcon2 * (vs_dv(i+1,j,k) - 2.0d0*vs_dv(i,j,k) +
     >                              vs_dv(i-1,j,k)) -
     >                    tx2 * (u_dv_r(i+1,j,k,3)*up1 -
     >                           u_dv_r(i-1,j,k,3)*um1)

                rhs_dv_r(i,j,k,4) = rhs_dv_r(i,j,k,4) + dx4tx1 *
     >                    (u_dv_r(i+1,j,k,4) - 2.0d0*u_dv_r(i,j,k,4) +
     >                     u_dv_r(i-1,j,k,4)) +
     >                  xxcon2 * (ws_dv(i+1,j,k) - 2.0d0*ws_dv(i,j,k) +
     >                              ws_dv(i-1,j,k)) -
     >                    tx2 * (u_dv_r(i+1,j,k,4)*up1 -
     >                           u_dv_r(i-1,j,k,4)*um1)

                rhs_dv_r(i,j,k,5) = rhs_dv_r(i,j,k,5) + dx5tx1 *
     >                    (u_dv_r(i+1,j,k,5) - 2.0d0*u_dv_r(i,j,k,5) +
     >                     u_dv_r(i-1,j,k,5)) +
     >                 xxcon3 * (qs_dv(i+1,j,k) - 2.0d0*qs_dv(i,j,k) +
     >                              qs_dv(i-1,j,k)) +
     >                    xxcon4 * (up1*up1 -       2.0d0*uijk*uijk +
     >                              um1*um1) +
     >                   xxcon5 * (u_dv_r(i+1,j,k,5)*rho_i_dv(i+1,j,k) -
     >                          2.0d0*u_dv_r(i,j,k,5)*rho_i_dv(i,j,k) +
     >                            u_dv_r(i-1,j,k,5)*rho_i_dv(i-1,j,k)) -
     >                    tx2 * ( (c1*u_dv_r(i+1,j,k,5) -
     >                             c2*square_dv(i+1,j,k))*up1 -
     >                            (c1*u_dv_r(i-1,j,k,5) -
     >                             c2*square_dv(i-1,j,k))*um1 )
       endif

       end subroutine krnl_rhs_xflux














       attributes(global) subroutine krnl_rhs_xdis1(u_dv_r,rhs_dv_r,
     >                            imax,  jmax, kmax,imaxp, jmaxp,
     >                                 dssp,istart)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dssp

        integer,value ::  istart

        integer i,j,k, blockSize, m,tSize, modF, ti
        integer uidx1, uidx2, uidx3, uidx4, uidx5, uidx6, uidx7

        double precision, shared :: temp(*)
        

        k = blockIDx%x            
        j = blockIDx%y
        tSize = 20
        modF = mod((threadIDx%x-1),5)
        ti = (threadIDx%x-1)/5+1

        if(istart .lt. 0) then
        if(threadIDx%x .le. 4) then
         do m=1,5
          i =1
          temp((m-1)*4+threadIDx%x ) = u_dv_r(threadIDx%x-1+i ,j,k,m)
          temp((m-1)*4+threadIDx%x +tSize)=
     >                                rhs_dv_r(threadIDx%x-1+i,j,k,m)
         end do
        endif
        else
         if(threadIDx%x .le. 4) then
         do m=1,5
          i = istart
          temp((m-1)*4+threadIDx%x ) = u_dv_r(threadIDx%x-1+i-2 ,j,k,m)
          temp((m-1)*4+threadIDx%x +tSize)=
     >                                rhs_dv_r(threadIDx%x-1+i ,j,k,m)
         end do
        endif


        endif

        call syncthreads()
    

        if(istart .lt. 0) then
            
             m = modF+1
             i = 1
                temp((m-1)*4+1+tSize) = temp((m-1)*4+1+tSize)- dssp * 
     >            ( 5.0d0* temp((m-1)*4+1)- 4.0d0*
     >                       temp((m-1)*4+2) +
     >                         temp((m-1)*4+3))

             i = 2
                temp((m-1)*4+2+tSize) = temp((m-1)*4+2+tSize) - dssp * 
     >             (-4.0d0*temp((m-1)*4+1) + 6.0d0*temp((m-1)*4+2) -
     >                 4.0d0*temp((m-1)*4+3) + temp((m-1)*4+4))

        else

             i = istart
             do     m = 1, 5
              temp((m-1)*4+1+tSize) = temp((m-1)*4+1+tSize) - dssp *
     >            ( temp((m-1)*4+1) - 4.0d0*temp((m-1)*4+2) + 
     >             6.0d0*temp((m-1)*4+3) - 4.0d0*temp((m-1)*4+4) )
             end do

             i = istart+1
             do     m = 1, 5
              temp((m-1)*4+2+tSize) = temp((m-1)*4+2+tSize) - dssp *
     >              ( temp((m-1)*4+2) - 4.d0*temp((m-1)*4+3) +
     >                      5.d0*temp((m-1)*4+4) )
             end do

        endif

        call syncthreads()


        if(istart .lt. 0) then
        if(threadIDx%x .le. 4) then
         i =1
         do m=1,5
          rhs_dv_r(threadIDx%x-1+i,j,k,m) =
     >                               temp((m-1)*4+threadIDx%x+tSize) 
         end do
        endif
        else
        if(threadIDx%x .le. 4) then
         i = istart
         do m=1,5
          rhs_dv_r(threadIDx%x-1+i ,j,k,m) =
     >                               temp((m-1)*4+threadIDx%x+tSize) 
         end do
        endif


        endif


 

       end subroutine krnl_rhs_xdis1












        attributes(global) subroutine krnl_rhs_xdis2(u_dv_r,rhs_dv_r,
     >                                  imax,  jmax, kmax,imaxp, jmaxp,
     >                                 dssp, iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                                             0:KMAX-1,5)::rhs_dv_r

        integer, value :: iSize, subSections

        double precision,value :: dssp

        integer :: i,j,k, i_idx, modF, blockSize, m
        
        double precision, shared :: temp(*)


        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/(subSections) + 1
	i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize) + 2

        if(i .le. iSize+2) then
          do     m = 1, 5
           ! m = modF+1
            rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp * 
     >            (  u_dv_r(i-2,j,k,m) - 4.0d0*u_dv_r(i-1,j,k,m) + 
     >               6.0*u_dv_r(i,j,k,m) - 4.0d0*u_dv_r(i+1,j,k,m) + 
     >                   u_dv_r(i+2,j,k,m) )
          end do
        endif

        end subroutine krnl_rhs_xdis2











       attributes(global) subroutine krnl_rhs_yflux(
     >           u_dv_r,us_dv, vs_dv, ws_dv,  qs_dv, rho_i_dv,
     >                      square_dv,rhs_dv_r, imax, jmax, kmax, 
     >                            imaxp, jmaxp,
     >                            dy1ty1, dy2ty1, dy3ty1,
     >         dy4ty1,dy5ty1, ty2, yycon2, yycon3, yycon4,yycon5,
     >         con43, c2, c1, iSize, subSections)
 
        integer, value :: IMAX, JMAX, KMAX, IMAXP, JMAXP

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::rho_i_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: us_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: vs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: ws_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::square_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: qs_dv
        double precision, dimension(0:IMAXP,0:JMAXP, 
     >                             0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dy1ty1, dy2ty1, dy3ty1,dy4ty1,dy5ty1, 
     >                 ty2,yycon2, yycon3, yycon4, yycon5,con43, c2, c1  
       
        integer, value :: iSize, subSections
 
        integer :: i,j,k,  blockSize

        double precision vijk, vp1, vm1
        double precision, shared :: temp(*) 

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections + 1
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)

        !blockSize = blockDim%x
        !k = blockIDx%x
        !j=blockIDx%y
        !i = threadIDx%x
     
        if(i .le. iSize) then 

                vijk = vs_dv(i,j,k)
                vp1  = vs_dv(i,j+1,k)
                vm1  = vs_dv(i,j-1,k)

                rhs_dv_r(i,j,k,1) = rhs_dv_r(i,j,k,1) + dy1ty1 * 
     >                   (u_dv_r(i,j+1,k,1) - 2.0d0*u_dv_r(i,j,k,1) + 
     >                    u_dv_r(i,j-1,k,1)) -
     >                   ty2 * (u_dv_r(i,j+1,k,3) - u_dv_r(i,j-1,k,3))

                rhs_dv_r(i,j,k,2) = rhs_dv_r(i,j,k,2) + dy2ty1 * 
     >                   (u_dv_r(i,j+1,k,2) - 2.0d0*u_dv_r(i,j,k,2) + 
     >                    u_dv_r(i,j-1,k,2)) +
     >                  yycon2 * (us_dv(i,j+1,k) - 2.0d0*us_dv(i,j,k) + 
     >                             us_dv(i,j-1,k)) -
     >                   ty2 * (u_dv_r(i,j+1,k,2)*vp1 - 
     >                          u_dv_r(i,j-1,k,2)*vm1)
  
              rhs_dv_r(i,j,k,3) = rhs_dv_r(i,j,k,3) + dy3ty1 * 
     >                   (u_dv_r(i,j+1,k,3) - 2.0d0*u_dv_r(i,j,k,3) + 
     >                    u_dv_r(i,j-1,k,3)) +
     >                   yycon2*con43 * (vp1 - 2.0d0*vijk + vm1) -
     >                   ty2 * (u_dv_r(i,j+1,k,3)*vp1 - 
     >                          u_dv_r(i,j-1,k,3)*vm1 +
     >                        (u_dv_r(i,j+1,k,5) - square_dv(i,j+1,k) - 
     >                          u_dv_r(i,j-1,k,5) + square_dv(i,j-1,k))
     >                          *c2)

                rhs_dv_r(i,j,k,4) = rhs_dv_r(i,j,k,4) + dy4ty1 * 
     >                   (u_dv_r(i,j+1,k,4) - 2.0d0*u_dv_r(i,j,k,4) + 
     >                    u_dv_r(i,j-1,k,4)) +
     >                  yycon2 * (ws_dv(i,j+1,k) - 2.0d0*ws_dv(i,j,k) + 
     >                             ws_dv(i,j-1,k)) -
     >                   ty2 * (u_dv_r(i,j+1,k,4)*vp1 - 
     >                          u_dv_r(i,j-1,k,4)*vm1)
     
           rhs_dv_r(i,j,k,5) = rhs_dv_r(i,j,k,5) + dy5ty1 * 
     >                   (u_dv_r(i,j+1,k,5) - 2.0d0*u_dv_r(i,j,k,5) + 
     >                    u_dv_r(i,j-1,k,5)) +
     >                  yycon3 * (qs_dv(i,j+1,k) - 2.0d0*qs_dv(i,j,k) + 
     >                             qs_dv(i,j-1,k)) +
     >                   yycon4 * (vp1*vp1       - 2.0d0*vijk*vijk + 
     >                             vm1*vm1) +
     >                  yycon5 * (u_dv_r(i,j+1,k,5)*rho_i_dv(i,j+1,k) - 
     >                          2.0d0*u_dv_r(i,j,k,5)*rho_i_dv(i,j,k) +
     >                           u_dv_r(i,j-1,k,5)*rho_i_dv(i,j-1,k)) -
     >                  ty2 * ((c1*u_dv_r(i,j+1,k,5) - 
     >                           c2*square_dv(i,j+1,k)) * vp1 -
     >                          (c1*u_dv_r(i,j-1,k,5) - 
     >                           c2*square_dv(i,j-1,k)) * vm1)

       endif

       end subroutine krnl_rhs_yflux








        attributes(global) subroutine krnl_rhs_yfluxB(
     >           u_dv_r,us_dv, vs_dv, ws_dv,  qs_dv, rho_i_dv,
     >                      square_dv,rhs_dv_r, imax, jmax, kmax,
     >                            imaxp, jmaxp,
     >                            dy1ty1, dy2ty1, dy3ty1,
     >         dy4ty1,dy5ty1, ty2, yycon2, yycon3, yycon4,yycon5,
     >         con43, c2, c1, iSize, jSize, subSectionsW, 
     >         subSectionsH)

        integer, value :: IMAX, JMAX, KMAX, IMAXP, JMAXP

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::rho_i_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: us_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: vs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: ws_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::square_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: qs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                             0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dy1ty1, dy2ty1, dy3ty1,dy4ty1,dy5ty1,
     >                 ty2,yycon2, yycon3, yycon4, yycon5,con43, c2, c1

        integer, value :: iSize, jSize, subSectionsW, subSectionsH

        integer :: i,j,k,  blockWidth, blockHeight

        double precision vijk, vp1, vm1
        double precision, shared :: temp(*)

        blockWidth = blockDim%x
        blockHeight = blockDim%y

        k = blockIDx%x
        j = ((blockIDx%y-1)/ subSectionsW) * blockHeight + threadIDx%y
        i = mod(blockIDx%y-1, subSectionsW) * blockWidth + threadIDx%x



        if(i .le. iSize .and. j .le. jSize) then
            vijk = vs_dv(i,j,k)
                vp1  = vs_dv(i,j+1,k)
                vm1  = vs_dv(i,j-1,k)

                rhs_dv_r(i,j,k,1) = rhs_dv_r(i,j,k,1) + dy1ty1 *
     >                   (u_dv_r(i,j+1,k,1) - 2.0d0*u_dv_r(i,j,k,1) +
     >                    u_dv_r(i,j-1,k,1)) -
     >                   ty2 * (u_dv_r(i,j+1,k,3) - u_dv_r(i,j-1,k,3))

                rhs_dv_r(i,j,k,2) = rhs_dv_r(i,j,k,2) + dy2ty1 *
     >                   (u_dv_r(i,j+1,k,2) - 2.0d0*u_dv_r(i,j,k,2) +
     >                    u_dv_r(i,j-1,k,2)) +
     >                  yycon2 * (us_dv(i,j+1,k) - 2.0d0*us_dv(i,j,k) +
     >                             us_dv(i,j-1,k)) -
     >                   ty2 * (u_dv_r(i,j+1,k,2)*vp1 -
     >                          u_dv_r(i,j-1,k,2)*vm1)

              rhs_dv_r(i,j,k,3) = rhs_dv_r(i,j,k,3) + dy3ty1 *
     >                   (u_dv_r(i,j+1,k,3) - 2.0d0*u_dv_r(i,j,k,3) +
     >                    u_dv_r(i,j-1,k,3)) +
     >                   yycon2*con43 * (vp1 - 2.0d0*vijk + vm1) -
     >                   ty2 * (u_dv_r(i,j+1,k,3)*vp1 -
     >                          u_dv_r(i,j-1,k,3)*vm1 +
     >                        (u_dv_r(i,j+1,k,5) - square_dv(i,j+1,k) -
     >                          u_dv_r(i,j-1,k,5) + square_dv(i,j-1,k))
     >                          *c2)

                rhs_dv_r(i,j,k,4) = rhs_dv_r(i,j,k,4) + dy4ty1 *
     >                   (u_dv_r(i,j+1,k,4) - 2.0d0*u_dv_r(i,j,k,4) +
     >                    u_dv_r(i,j-1,k,4)) +
     >                  yycon2 * (ws_dv(i,j+1,k) - 2.0d0*ws_dv(i,j,k) +
     >                             ws_dv(i,j-1,k)) -
     >                   ty2 * (u_dv_r(i,j+1,k,4)*vp1 -
     >                          u_dv_r(i,j-1,k,4)*vm1)

           rhs_dv_r(i,j,k,5) = rhs_dv_r(i,j,k,5) + dy5ty1 *
     >                   (u_dv_r(i,j+1,k,5) - 2.0d0*u_dv_r(i,j,k,5) +
     >                    u_dv_r(i,j-1,k,5)) +
     >                  yycon3 * (qs_dv(i,j+1,k) - 2.0d0*qs_dv(i,j,k) +
     >                             qs_dv(i,j-1,k)) +
     >                   yycon4 * (vp1*vp1       - 2.0d0*vijk*vijk +
     >                             vm1*vm1) +
     >                  yycon5 * (u_dv_r(i,j+1,k,5)*rho_i_dv(i,j+1,k) -
     >                          2.0d0*u_dv_r(i,j,k,5)*rho_i_dv(i,j,k) +
     >                           u_dv_r(i,j-1,k,5)*rho_i_dv(i,j-1,k)) -
     >                  ty2 * ((c1*u_dv_r(i,j+1,k,5) -
     >                           c2*square_dv(i,j+1,k)) * vp1 -
     >                          (c1*u_dv_r(i,j-1,k,5) -
     >                           c2*square_dv(i,j-1,k)) * vm1)

       endif

       end subroutine krnl_rhs_yfluxB







      

         attributes(global) subroutine krnl_rhs_ydis1(u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,jstart, iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dssp

        integer,value :: jstart,iSize, subSections 
        integer :: i,j,k, i_idx, blockSize, m

        double precision, shared :: temp(*)

        blockSize = blockDim%x

        k = (blockIDx%x-1)/subSections +1
        i = threadIDx%x + (mod(blockIDx%x-1,subSections)*blockSize)

        if(i .le. iSize) then 

           if(jstart .le. 0) then

              j = 1
            ! do     i = 1, nx2
                do     m = 1, 5
                  rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m)- dssp * 
     >             ( 5.0d0*u_dv_r(i,j,k,m) - 4.0d0*u_dv_r(i,j+1,k,m) +
     >                          u_dv_r(i,j+2,k,m))
                end do
            ! end do

             j = 2
            ! do     i = 1, nx2
                 do     m = 1, 5
                  rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp * 
     >            (-4.0d0*u_dv_r(i,j-1,k,m) + 6.0d0*u_dv_r(i,j,k,m) -
     >             4.0d0*u_dv_r(i,j+1,k,m) + u_dv_r(i,j+2,k,m))
                 end do
            ! end do

           else

          j = jstart
          !do     i = 1, nx2
             do     m = 1, 5
                rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp *
     >              ( u_dv_r(i,j-2,k,m) - 4.0d0*u_dv_r(i,j-1,k,m) + 
     >               6.0d0*u_dv_r(i,j,k,m) - 4.0d0*u_dv_r(i,j+1,k,m) )
             end do
          !end do

          j = jstart+1
          !do     i = 1, nx2
             do     m = 1, 5
                rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp *
     >             ( u_dv_r(i,j-2,k,m) - 4.d0*u_dv_r(i,j-1,k,m) +
     >                 5.d0*u_dv_r(i,j,k,m) )
             end do
          !end do

           endif

        endif
       

        end subroutine krnl_rhs_ydis1













        attributes(global) subroutine krnl_rhs_ydis2(u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,iSize,subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dssp

        integer,value :: iSize, subSections

        integer :: i,j,k, modF, blockSize, m
       
        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/(subSections) +3
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)

        if(i .le. iSize) then 

         do   m = 1, 5
         !m = modF+1
               rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp * 
     >           (  u_dv_r(i,j-2,k,m) - 4.0d0*u_dv_r(i,j-1,k,m) + 
     >            6.0*u_dv_r(i,j,k,m) - 4.0d0*u_dv_r(i,j+1,k,m) + 
     >                u_dv_r(i,j+2,k,m) )
          end do
     
        endif

        end subroutine krnl_rhs_ydis2













       attributes(global) subroutine krnl_rhs_zflux(u_dv_r,us_dv, 
     >         vs_dv, ws_dv,  qs_dv, rho_i_dv,square_dv,rhs_dv_r,
     >         imax, jmax, kmax,imaxp, jmaxp, dz1tz1, dz2tz1, dz3tz1,
     >           dz4tz1,dz5tz1, tz2, zzcon2, zzcon3, zzcon4, 
     >           zzcon5, con43, c2, c1, iSize, subSections) 

        integer, value :: IMAX, JMAX, KMAX, IMAXP, JMAXP

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP, 0:KMAX-1)::rho_i_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: us_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: vs_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: ws_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::square_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1) :: qs_dv
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dz1tz1, dz2tz1, dz3tz1,dz4tz1,dz5tz1,
     >          tz2,zzcon2, zzcon3, zzcon4,zzcon5,con43, c2, c1  
       
        integer, value :: iSize, subSections
 
        integer :: i,j,k, i_idx, modF, blockSize
        double precision wijk, wp1, wm1
        double precision, shared :: temp(*) 

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)
   

        if(i .le. iSize) then 

                wijk = ws_dv(i,j,k)
                wp1  = ws_dv(i,j,k+1)
                wm1  = ws_dv(i,j,k-1)

         rhs_dv_r(i,j,k,1) = rhs_dv_r(i,j,k,1) + dz1tz1 * 
     >                   (u_dv_r(i,j,k+1,1) - 2.0d0*u_dv_r(i,j,k,1) + 
     >                    u_dv_r(i,j,k-1,1)) -
     >                   tz2 * (u_dv_r(i,j,k+1,4) - u_dv_r(i,j,k-1,4))

                rhs_dv_r(i,j,k,2) = rhs_dv_r(i,j,k,2) + dz2tz1 * 
     >                   (u_dv_r(i,j,k+1,2) - 2.0d0*u_dv_r(i,j,k,2) + 
     >                    u_dv_r(i,j,k-1,2)) +
     >                 zzcon2 * (us_dv(i,j,k+1) - 2.0d0*us_dv(i,j,k) + 
     >                             us_dv(i,j,k-1)) -
     >                   tz2 * (u_dv_r(i,j,k+1,2)*wp1 - 
     >                          u_dv_r(i,j,k-1,2)*wm1)

                rhs_dv_r(i,j,k,3) = rhs_dv_r(i,j,k,3) + dz3tz1 * 
     >                   (u_dv_r(i,j,k+1,3) - 2.0d0*u_dv_r(i,j,k,3) + 
     >                    u_dv_r(i,j,k-1,3)) +
     >                  zzcon2 * (vs_dv(i,j,k+1) - 2.0d0*vs_dv(i,j,k) + 
     >                             vs_dv(i,j,k-1)) -
     >                   tz2 * (u_dv_r(i,j,k+1,3)*wp1 - 
     >                          u_dv_r(i,j,k-1,3)*wm1)

                rhs_dv_r(i,j,k,4) = rhs_dv_r(i,j,k,4) + dz4tz1 * 
     >                   (u_dv_r(i,j,k+1,4) - 2.0d0*u_dv_r(i,j,k,4) + 
     >                    u_dv_r(i,j,k-1,4)) +
     >                   zzcon2*con43 * (wp1 - 2.0d0*wijk + wm1) -
     >                   tz2 * (u_dv_r(i,j,k+1,4)*wp1 - 
     >                          u_dv_r(i,j,k-1,4)*wm1 +
     >                        (u_dv_r(i,j,k+1,5) - square_dv(i,j,k+1) - 
     >                          u_dv_r(i,j,k-1,5) + square_dv(i,j,k-1))
     >                          *c2)

                rhs_dv_r(i,j,k,5) = rhs_dv_r(i,j,k,5) + dz5tz1 * 
     >                   (u_dv_r(i,j,k+1,5) - 2.0d0*u_dv_r(i,j,k,5) + 
     >                    u_dv_r(i,j,k-1,5)) +
     >                zzcon3 * (qs_dv(i,j,k+1) - 2.0d0*qs_dv(i,j,k) + 
     >                             qs_dv(i,j,k-1)) +
     >                   zzcon4 * (wp1*wp1 - 2.0d0*wijk*wijk + 
     >                             wm1*wm1) +
     >                  zzcon5 * (u_dv_r(i,j,k+1,5)*rho_i_dv(i,j,k+1) - 
     >                          2.0d0*u_dv_r(i,j,k,5)*rho_i_dv(i,j,k) +
     >                           u_dv_r(i,j,k-1,5)*rho_i_dv(i,j,k-1)) -
     >                 tz2 * ( (c1*u_dv_r(i,j,k+1,5) - 
     >                            c2*square_dv(i,j,k+1))*wp1 -
     >                           (c1*u_dv_r(i,j,k-1,5) - 
     >                            c2*square_dv(i,j,k-1))*wm1)
        endif

       
       end subroutine krnl_rhs_zflux








         attributes(global) subroutine krnl_rhs_zdis1(u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,kstart, iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dssp

        integer,value :: kstart,iSize, subSections

        integer :: i,j,k, blockSize, m

        double precision, shared :: temp(*)

        blockSize =blockDim%x
        j = (blockIDx%x-1)/subSections +1
        i = threadIDx%x + (mod(blockIDx%x-1,subSections)*blockSize)

        if(i .le. iSize) then 
        
          if(kstart .le. 0) then
        
            k = 1
            do     m = 1, 5
                rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m)- dssp * 
     >           ( 5.0d0*u_dv_r(i,j,k,m) - 4.0d0*u_dv_r(i,j,k+1,m) +
     >                            u_dv_r(i,j,k+2,m))
            end do

            k = 2
            do     m = 1, 5
                rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp * 
     >          (-4.0d0*u_dv_r(i,j,k-1,m) + 6.0d0*u_dv_r(i,j,k,m) -
     >             4.0d0*u_dv_r(i,j,k+1,m) + u_dv_r(i,j,k+2,m))
            end do
 
          else

             k = kstart
             do     m = 1, 5
                rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp *
     >              ( u_dv_r(i,j,k-2,m) - 4.0d0*u_dv_r(i,j,k-1,m) + 
     >               6.0d0*u_dv_r(i,j,k,m) - 4.0d0*u_dv_r(i,j,k+1,m) )
             end do

             k = kstart + 1
             do     m = 1, 5
                rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp *
     >              ( u_dv_r(i,j,k-2,m) - 4.d0*u_dv_r(i,j,k-1,m) +
     >                      5.d0*u_dv_r(i,j,k,m) )
             end do

          endif
        endif


        end subroutine krnl_rhs_zdis1










        attributes(global) subroutine krnl_rhs_zdis2(u_dv_r,rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dssp,iSize,subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::rhs_dv_r

        double precision,value :: dssp

        integer,value ::  iSize, subSections

        integer :: i,j,k, modF, blockSize, m

        double precision, shared :: temp(*)

        blockSize = blockDim%x
        k = blockIDx%x +2
        j = (blockIDx%y-1)/(subSections) +1
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)
        if( i .le. iSize) then 

            do     m = 1, 5
              rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) - dssp * 
     >             (  u_dv_r(i,j,k-2,m) - 4.0d0*u_dv_r(i,j,k-1,m) + 
     >               6.0*u_dv_r(i,j,k,m) - 4.0d0*u_dv_r(i,j,k+1,m) + 
     >                         u_dv_r(i,j,k+2,m) )
            end do

        endif

        end subroutine krnl_rhs_zdis2








        attributes(global) subroutine krnl_rhs_final(rhs_dv_r,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                dt,iSize,subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp, maxcells

        double precision,dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::rhs_dv_r
      
        double precision, value :: dt 

        integer,value :: iSize, subSections

        integer :: i,j,k,  modF, blockSize, m
       
        double precision, shared :: temp(*) 

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)

        if(i .le. iSize) then
          do    m = 1, 5
            rhs_dv_r(i,j,k,m) = rhs_dv_r(i,j,k,m) * dt
          end do
        endif

       end subroutine krnl_rhs_final
    




       end module rhs_kernels 
