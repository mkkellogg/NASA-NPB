         module util_kernels
         implicit none
         contains




         attributes(global) subroutine krnl_util_swap4D_mtolast(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(5,0:IMAXP,0:JMAXP,0:KMAX-1)::src_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::dest_dv

        integer :: i,j,k,m

        double precision :: t
       
        double precision, shared :: temp(*) 

        k = blockIDx%x-1        
        j = blockIDx%y-1        
        i = threadIDx%x-1


        do m =1,5  
          
        t = src_dv(m,i,j,k)
        dest_dv(i,j,k,m) = t

        end do
      

        end subroutine krnl_util_swap4D_mtolast
    






        attributes(global) subroutine krnl_util_swap4D_mtofirst(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::src_dv
        double precision, dimension(5,0:IMAXP,0:JMAXP,0:KMAX-1)::dest_dv

        integer :: i,j,k,m

        double precision :: t
       
        double precision, shared :: temp(*) 

        k = blockIDx%x-1        
        j = blockIDx%y-1        
        i = threadIDx%x-1


        do m =1,5  
          
        t = src_dv(i,j,k,m)
        dest_dv(m,i,j,k) = t

        end do
      

        end subroutine krnl_util_swap4D_mtofirst




    





        attributes(global) subroutine krnl_util_swap4D_itoj_c(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                blockDimSize,iBlocks,jBlocks)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::src_dv
        double precision, dimension(0:JMAXP,0:IMAXP,0:KMAX-1,5)::dest_dv

        integer, value :: blockDimSize, iBlocks, jBlocks

        integer :: i,j,k,basei,basej,baseti,basetj,m, row, col

        double precision t

        double precision, shared :: temp(*)

        k = blockIDx%x-1

        do m = 1,5

        basej = ((blockIDx%y-1)/iBlocks) * blockDimSize
        basei = (mod(blockIDx%y-1,iBlocks)) * blockDimSize

           if(basej+threadIDx%y-1 .le. jmaxp .and.
     >         basei+threadIDx%x-1 .le. imaxp) then
           temp((threadIDx%y-1)*blockDimSize+threadIDx%x) =
     >           src_dv(basei+threadIDx%x-1,basej+threadIDx%y-1,k,m)
           endif

        call syncthreads()

        !basei = ((blockIDx%y-1)/jBlocks) * blockDimSize
        !basej = (mod(blockIDx%y-1,jBlocks)) * blockDimSize

            if(basei+ threadIDx%y-1 .le. imaxp .and. 
     >         basej+threadIDx%x-1 .le. jmaxp) then   
            dest_dv(basej+threadIDx%x-1,basei+threadIDx%y-1,k,m) =        
     >                temp((threadIDx%x-1)*blockDimSize+threadIDx%y)
            endif

        call syncthreads()

        end do
           

        end subroutine krnl_util_swap4D_itoj_c












        attributes(global) subroutine krnl_util_swap4D_jtoi_c(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                blockDimSize,iBlocks,jBlocks)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:JMAXP,0:IMAXP,0:KMAX-1,5)::src_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::dest_dv

        integer, value :: blockDimSize, iBlocks, jBlocks

        integer :: i,j,k,basei,basej,baseti,basetj,m, row, col

        double precision t

        double precision, shared :: temp(*)

        k = blockIDx%x-1

        do m = 1,5

        basei = ((blockIDx%y-1)/jBlocks) * blockDimSize
        basej = (mod(blockIDx%y-1,jBlocks)) * blockDimSize

           if(basej+threadIDx%x-1 .le. jmaxp .and.
     >         basei+threadIDx%y-1 .le. imaxp) then
           temp((threadIDx%y-1)*blockDimSize+threadIDx%x) =
     >           src_dv(basej+threadIDx%x-1,basei+threadIDx%y-1,k,m)
           endif

        call syncthreads()


            if(basei+ threadIDx%x-1 .le. imaxp .and.
     >         basej+threadIDx%y-1 .le. jmaxp) then
            dest_dv(basei+threadIDx%x-1,basej+threadIDx%y-1,k,m) =
     >                temp((threadIDx%x-1)*blockDimSize+threadIDx%y)
            endif

        call syncthreads()

        end do


        end subroutine krnl_util_swap4D_jtoi_c









        attributes(global) subroutine krnl_util_swap4D_itok_c(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                blockDimSize,iBlocks,kBlocks)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::src_dv
        double precision, dimension(0:KMAX-1,0:JMAXP,0:IMAXP,5)::dest_dv

        integer, value :: blockDimSize, iBlocks, kBlocks

        integer :: i,j,k,basei,basek,baseti,basetk,m, row, col

        double precision t

        double precision, shared :: temp(*)

        j = blockIDx%x-1

        do m = 1,5

        basek = ((blockIDx%y-1)/iBlocks) * blockDimSize
        basei = (mod(blockIDx%y-1,iBlocks)) * blockDimSize

           if(basek+threadIDx%y-1 .le. kmax-1 .and.
     >         basei+threadIDx%x-1 .le. imaxp) then
           temp((threadIDx%y-1)*blockDimSize+threadIDx%x) =
     >           src_dv(basei+threadIDx%x-1,j,basek+threadIDx%y-1,m)
           endif

        call syncthreads()

        !basei = ((blockIDx%y-1)/jBlocks) * blockDimSize
        !basej = (mod(blockIDx%y-1,jBlocks)) * blockDimSize

            if(basei+ threadIDx%y-1 .le. imaxp .and. 
     >         basek+threadIDx%x-1 .le. kmax-1) then   
            dest_dv(basek+threadIDx%x-1,j,basei+threadIDx%y-1,m) =        
     >                temp((threadIDx%x-1)*blockDimSize+threadIDx%y)
            endif

        call syncthreads()

        end do
           

        end subroutine krnl_util_swap4D_itok_c









        attributes(global) subroutine krnl_util_swap4D_ktoi_c(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                blockDimSize,iBlocks,kBlocks)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:KMAX-1,0:JMAXP,0:IMAXP,5)::src_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::dest_dv

        integer, value :: blockDimSize, iBlocks, kBlocks

        integer :: i,j,k,basei,basek,baseti,basetk,m, row, col

        double precision t

        double precision, shared :: temp(*)

        j = blockIDx%x-1

        do m = 1,5

        basei = ((blockIDx%y-1)/kBlocks) * blockDimSize
        basek = (mod(blockIDx%y-1,kBlocks)) * blockDimSize

           if(basek+threadIDx%x-1 .le. kmax-1 .and.
     >         basei+threadIDx%y-1 .le. imaxp) then
           temp((threadIDx%y-1)*blockDimSize+threadIDx%x) =
     >           src_dv(basek+threadIDx%x-1,j,basei+threadIDx%y-1,m)
           endif

        call syncthreads()


            if(basei+ threadIDx%x-1 .le. imaxp .and.
     >         basek+threadIDx%y-1 .le. kmax-1) then
            dest_dv(basei+threadIDx%x-1,j,basek+threadIDx%y-1,m) =
     >                temp((threadIDx%x-1)*blockDimSize+threadIDx%y)
            endif

        call syncthreads()

        end do


        end subroutine krnl_util_swap4D_ktoi_c











        attributes(global) subroutine krnl_util_swap4D_itoj(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::src_dv
        double precision, dimension(0:JMAXP,0:IMAXP,0:KMAX-1,5)::dest_dv

        integer :: i,j,k,m

        double precision t
       
        double precision, shared :: temp(*) 

        k = blockIDx%x-1        
        j = blockIDx%y-1        
        i = threadIDx%x-1


        do m =1,5  
          
        t = src_dv(i,j,k,m)
        dest_dv(j,i,k,m) = t

        end do

        end subroutine krnl_util_swap4D_itoj








        attributes(global) subroutine krnl_util_swap4D_jtoi(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:JMAXP,0:IMAXP,0:KMAX-1,5)::src_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::dest_dv

        integer :: i,j,k,m
 
        double precision t

        double precision, shared :: temp(*)

        k = blockIDx%x-1
        j = blockIDx%y-1
        i = threadIDx%x-1


        do m =1,5

        t = src_dv(j,i,k,m)
        dest_dv(i,j,k,m) = t

        end do


        end subroutine krnl_util_swap4D_jtoi









        attributes(global) subroutine krnl_util_swap4D_itok(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::src_dv
        double precision, dimension(0:KMAX-1,0:JMAXP,0:IMAXP,5)::dest_dv

        integer :: i,j,k,m

        double precision t
       
        double precision, shared :: temp(*) 

        k = blockIDx%x-1        
        j = blockIDx%y-1        
        i = threadIDx%x-1


        do m =1,5  
          
        t = src_dv(i,j,k,m)
        dest_dv(k,j,i,m) = t

        end do

        end subroutine krnl_util_swap4D_itok








        attributes(global) subroutine krnl_util_swap4D_ktoi(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:KMAX-1,0:JMAXP,0:IMAXP,5)::src_dv
        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1,5)::dest_dv

        integer :: i,j,k,m
 
        double precision t

        double precision, shared :: temp(*)

        k = blockIDx%x-1
        j = blockIDx%y-1
        i = threadIDx%x-1


        do m =1,5

        t = src_dv(k,j,i,m)
        dest_dv(i,j,k,m) = t

        end do


        end subroutine krnl_util_swap4D_ktoi







        attributes(global) subroutine krnl_util_swap3D_itoj_c(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                blockDimSize,iBlocks,jBlocks)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::src_dv
        double precision, dimension(0:JMAXP,0:IMAXP,0:KMAX-1)::dest_dv

        integer, value :: blockDimSize, iBlocks, jBlocks

        integer :: i,j,k,basei,basej,baseti,basetj,m, row, col

        double precision t

        double precision, shared :: temp(*)

        k = blockIDx%x-1

        basej = ((blockIDx%y-1)/iBlocks) * blockDimSize
        basei = (mod(blockIDx%y-1,iBlocks)) * blockDimSize

           if(basej+threadIDx%y-1 .le. jmaxp .and.
     >         basei+threadIDx%x-1 .le. imaxp) then
           temp((threadIDx%y-1)*blockDimSize+threadIDx%x) =
     >           src_dv(basei+threadIDx%x-1,basej+threadIDx%y-1,k)
           endif

        call syncthreads()

        !basei = ((blockIDx%y-1)/jBlocks) * blockDimSize
        !basej = (mod(blockIDx%y-1,jBlocks)) * blockDimSize

            if(basei+ threadIDx%y-1 .le. imaxp .and.
     >         basej+threadIDx%x-1 .le. jmaxp) then
            dest_dv(basej+threadIDx%x-1,basei+threadIDx%y-1,k) =
     >                temp((threadIDx%x-1)*blockDimSize+threadIDx%y)
            endif

        call syncthreads()


        end subroutine krnl_util_swap3D_itoj_c







        attributes(global) subroutine krnl_util_swap3D_itok_c(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp,
     >                blockDimSize,iBlocks,kBlocks)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::src_dv
        double precision, dimension(0:KMAX-1, 0:JMAXP,0:IMAXP)::dest_dv

        integer, value :: blockDimSize, iBlocks, kBlocks

        integer :: i,j,k,basei,basek,baseti,basetk,m, row, col

        double precision t

        double precision, shared :: temp(*)

        j = blockIDx%x-1

        basek = ((blockIDx%y-1)/iBlocks) * blockDimSize
        basei = (mod(blockIDx%y-1,iBlocks)) * blockDimSize

           if(basek+threadIDx%y-1 .le. kmax-1 .and.
     >         basei+threadIDx%x-1 .le. imaxp) then
           temp((threadIDx%y-1)*blockDimSize+threadIDx%x) =
     >           src_dv(basei+threadIDx%x-1,j,basek+threadIDx%y-1)
           endif

        call syncthreads()

        !basei = ((blockIDx%y-1)/jBlocks) * blockDimSize
        !basej = (mod(blockIDx%y-1,jBlocks)) * blockDimSize

            if(basei+ threadIDx%y-1 .le. imaxp .and.
     >         basek+threadIDx%x-1 .le. kmax-1) then
            dest_dv(basek+threadIDx%x-1,j,basei+threadIDx%y-1) =
     >                temp((threadIDx%x-1)*blockDimSize+threadIDx%y)
            endif

        call syncthreads()


        end subroutine krnl_util_swap3D_itok_c







        attributes(global) subroutine krnl_util_swap3D_itok(
     >                src_dv, dest_dv,
     >                imax,  jmax, kmax,imaxp, jmaxp)

        integer, value :: imax, jmax, kmax, imaxp, jmaxp

        double precision, dimension(0:IMAXP,0:JMAXP,0:KMAX-1)::src_dv
        double precision, dimension(0:KMAX-1,0:JMAXP,0:IMAXP)::dest_dv

        integer :: i,j,k

        double precision t
       
        double precision, shared :: temp(*) 

        k = blockIDx%x-1        
        j = blockIDx%y-1        
        i = threadIDx%x-1

          
        t = src_dv(i,j,k)
        dest_dv(k,j,i) = t


        end subroutine krnl_util_swap3D_itok





        end module util_kernels 
