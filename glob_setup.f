        subroutine glob_setup
         use glob
         use cudafor
         include "header.h"  

         integer :: istat,i
   

        if(doGPU .gt. 0) then
           istat = cudaStreamCreate(stream1)
           istat = cudaStreamCreate(stream2)
           istat = cudaStreamCreate(stream3)

           MaxDynamicEvents = 36
           dynamicEventsInUse = 0

           allocate(dynamicEvents(1:MaxDynamicEvents))
           do i =1,MaxDynamicEvents
             istat = cudaEventCreate(dynamicEvents(i))
           end do

           istat = cudaEventCreate(krnlStart);
           istat = cudaEventCreate(krnlEnd);


           allocate(u_dv(5, 0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(rho_i_dv(0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(us_dv(0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(vs_dv(0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(ws_dv(0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(qs_dv(0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(ainv_dv(0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(speed_dv(0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(square_dv(0:IMAXP, 0:JMAXP, 0:KMAX-1))

           allocate(rho_i_dv_j_r(0:JMAXP, 0:IMAXP, 0:KMAX-1))
           allocate(vs_dv_j_r(0:JMAXP, 0:IMAXP, 0:KMAX-1))

           allocate(rho_i_dv_k_r(0:KMAX-1, 0:JMAXP, 0:IMAXP))
           allocate(ws_dv_k_r(0:KMAX-1, 0:JMAXP, 0:IMAXP))

           allocate(lhs_dv(5, 0:IMAXP, 0:IMAXP))
           allocate(lhsp_dv(5,0:IMAXP,0:IMAXP))
           allocate(lhsm_dv(5,0:IMAXP,0:IMAXP))

           allocate(lhs_x_dv(5, 0:IMAXP, 0:IMAXP,0:KMAX-1))
           allocate(lhsp_x_dv(5,0:IMAXP,0:IMAXP,0:KMAX-1))
           allocate(lhsm_x_dv(5,0:IMAXP,0:IMAXP,0:KMAX-1))

           allocate(rhs_dv(5, 0:IMAXP, 0:JMAXP, 0:KMAX-1))
           allocate(forcing_dv(5, 0:IMAXP, 0:JMAXP, 0:KMAX-1))

           allocate(u_dv_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
           allocate(lhs_x_dv_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
           allocate(lhsp_x_dv_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
           allocate(lhsm_x_dv_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
           allocate(rhs_dv_r(0:IMAXP, 0:JMAXP, 0:KMAX-1,5))
           allocate(forcing_dv_r(0:IMAXP, 0:JMAXP, 0:KMAX-1,5))

           allocate(lhs_x_dv_j_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
           allocate(lhsp_x_dv_j_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
           allocate(lhsm_x_dv_j_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
           allocate(rhs_dv_j_r(0:JMAXP, 0:IMAXP, 0:KMAX-1,5))

           allocate(lhs_x_dv_k_r(0:KMAX-1, 0:IMAXP, 0:IMAXP,5))
        endif

        allocate(lhs(5, 0:IMAXP, 0:IMAXP))
        allocate(lhsp(5,0:IMAXP,0:IMAXP))
        allocate(lhsm(5,0:IMAXP,0:IMAXP))

        allocate(lhs_x(5, 0:IMAXP, 0:IMAXP,0:KMAX-1))
        allocate(lhsp_x(5,0:IMAXP,0:IMAXP,0:KMAX-1))
        allocate(lhsm_x(5,0:IMAXP,0:IMAXP,0:KMAX-1))

        allocate(rhs(5, 0:IMAXP, 0:JMAXP, 0:KMAX-1))

        allocate(u_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
        allocate(forcing_r(0:IMAXP, 0:JMAXP, 0:KMAX-1,5))
        allocate(lhs_x_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
        allocate(lhsp_x_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
        allocate(lhsm_x_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
        allocate(rhs_r(0:IMAXP, 0:JMAXP, 0:KMAX-1,5))

 
        allocate(lhs_x_j_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
        allocate(lhsp_x_j_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
        allocate(lhsm_x_j_r(0:IMAXP, 0:IMAXP,0:KMAX-1,5))
        allocate(rhs_j_r(0:JMAXP, 0:IMAXP, 0:KMAX-1,5))

        !istat = cudaDeviceSetCacheConfig(cudaFuncCachePreferShared)

        end subroutine glob_setup





        real function dynamicTimeElapsed(eventCount)
        use glob
        use cudafor
        include "header.h"
 
        integer::eventCount

        integer :: i, istat

        real tink, total

        total = 0
        do i=1,eventCount,2
          istat = cudaEventElapsedTime(tink, dynamicEvents(i), 
     >                                dynamicEvents(i+1))
          !call timer_inc(t_kernel,dble(tink/1000.0))
           total = total + tink
        end do

        dynamicTimeElapsed = total

        end 





        subroutine glob_cleanup
          use glob 
          use cudafor
      
          include "header.h"

        integer :: istat, i


        if(doGPU .gt. 0) then

          istat = cudaStreamDestroy(stream1)
          istat = cudaStreamDestroy(stream2)
          istat = cudaStreamDestroy(stream3)

          do i =1,MaxDynamicEvents
            istat = cudaEventDestroy(dynamicEvents(i))
          end do
          deallocate(dynamicEvents)


          istat = cudaEventDestroy(krnlStart);
          istat = cudaEventDestroy(krnlEnd);


          deallocate(u_dv)
          deallocate(us_dv)
          deallocate(vs_dv)
          deallocate(ws_dv)
          deallocate(qs_dv)
          deallocate(ainv_dv)
          deallocate(rho_i_dv)
          deallocate(speed_dv)
          deallocate(square_dv)
          deallocate(rhs_dv)
          deallocate(forcing_dv)

          deallocate(rho_i_dv_j_r)
          deallocate(vs_dv_j_r)

          deallocate(rho_i_dv_k_r)
          deallocate(ws_dv_k_r)



          deallocate(lhs_dv)
          deallocate(lhsp_dv)
          deallocate(lhsm_dv)

          deallocate(lhs_x_dv)
          deallocate(lhsp_x_dv)
          deallocate(lhsm_x_dv)

          deallocate(u_dv_r)
          deallocate(lhs_x_dv_r)
          deallocate(forcing_dv_r)
          deallocate(rhs_dv_r)

          deallocate(lhs_x_dv_j_r)
          deallocate(rhs_dv_j_r)

          deallocate(lhs_x_dv_k_r)
        endif

        deallocate(lhs)
        deallocate(lhsp)
        deallocate(lhsm)

        deallocate(lhs_x)
        deallocate(lhsp_x)
        deallocate(lhsm_x)

        deallocate(rhs)

        deallocate(u_r)
        deallocate(rhs_r)
        deallocate(forcing_r)
        deallocate(lhs_x_r)

        deallocate(rhs_j_r)
        deallocate(lhs_x_j_r)

        end subroutine glob_cleanup
