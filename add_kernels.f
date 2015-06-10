         module add_kernels
         implicit none
         contains




        attributes(global) subroutine krnl_add(
     >     u_dv_r, rhs_dv_r,imax,  jmax, kmax, imaxp, jmaxp,
     >                    iSize, subSections)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::u_dv_r
        double precision, dimension(0:IMAXP,0:JMAXP,
     >                              0:KMAX-1,5) :: rhs_dv_r

        integer, value :: iSize, subSections
       
         integer :: i,j,k, blockSize, m, joff

        double precision, shared :: temp(*) 

        blockSize = blockDim%x
        k = blockIDx%x
        j = (blockIDx%y-1)/subSections +1
        i = threadIDx%x + (mod(blockIDx%y-1,subSections)*blockSize)

         if(i .le. iSize) then
         do m = 1, 5
                   u_dv_r(i,j,k,m) = u_dv_r(i,j,k,m) + rhs_dv_r(i,j,k,m)
         end do
         endif

         end subroutine krnl_add




         end module add_kernels
