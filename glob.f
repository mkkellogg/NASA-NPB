         module glob
          use cudafor
           IMPLICIT NONE
    
          type(dim3) :: dimGrid, dimBlock
          type(cudaEvent) :: krnlStart, krnlEnd
          
          integer MaxDynamicEvents
          integer dynamicEventsInUse
          type(cudaEvent), allocatable :: dynamicEvents(:) 

          integer stream1, stream2, stream3
     
	!$acc mirror(lhs_1D) 
	!$acc mirror(cv_full)
	!$acc mirror(rhon_full)

        double precision, device, allocatable :: u_dv(:,:,:,:)
        double precision, device, allocatable :: rho_i_dv(:,:,:)
        double precision, device, allocatable :: us_dv(:,:,:)
        double precision, device, allocatable :: vs_dv(:,:,:)
        double precision, device, allocatable :: ws_dv(:,:,:)
        double precision, device, allocatable :: square_dv(:,:,:)
        double precision, device, allocatable :: qs_dv(:,:,:)
        double precision, device, allocatable :: speed_dv(:,:,:)
        double precision, device, allocatable :: ainv_dv(:,:,:)

        double precision, device, allocatable :: rho_i_dv_j_r(:,:,:)
        double precision, device, allocatable :: vs_dv_j_r(:,:,:)

        double precision, device, allocatable :: rho_i_dv_k_r(:,:,:)
        double precision, device, allocatable :: ws_dv_k_r(:,:,:)

        double precision, device, allocatable :: forcing_dv(:,:,:,:)
        double precision, device, allocatable :: rhs_dv(:,:,:,:)
        double precision, device, allocatable :: lhs_dv(:,:,:)
        double precision, device, allocatable :: lhsp_dv(:,:,:)
        double precision, device, allocatable :: lhsm_dv(:,:,:)

        double precision, device, allocatable :: lhs_x_dv(:,:,:,:)
        double precision, device, allocatable :: lhsp_x_dv(:,:,:,:)
        double precision, device, allocatable :: lhsm_x_dv(:,:,:,:)

        double precision, device, allocatable :: in_buffer_dv(:)
        double precision, device, allocatable :: out_buffer_dv(:)

        double precision, device, allocatable :: rhs_dv_r(:,:,:,:)
        double precision, device, allocatable :: forcing_dv_r(:,:,:,:)
        double precision, device, allocatable :: u_dv_r(:,:,:,:)
        double precision, device, allocatable :: lhs_x_dv_r(:,:,:,:)
	double precision, device, allocatable :: lhsp_x_dv_r(:,:,:,:)
	double precision, device, allocatable :: lhsm_x_dv_r(:,:,:,:)
 

        double precision, device, allocatable :: rhs_dv_j_r(:,:,:,:)
        double precision, device, allocatable :: lhs_x_dv_j_r(:,:,:,:)
        double precision, device, allocatable :: lhsp_x_dv_j_r(:,:,:,:)
        double precision, device, allocatable :: lhsm_x_dv_j_r(:,:,:,:)

        double precision, device, allocatable :: lhs_x_dv_k_r(:,:,:,:)

         double precision, allocatable, pinned ::
     >     lhs(:,:,:)
         double precision, allocatable, pinned ::
     >     lhsp(:,:,:)
         double precision, allocatable, pinned ::
     >     lhsm(:,:,:)


         double precision, allocatable, pinned ::
     >     lhs_x(:,:,:,:)
         double precision, allocatable, pinned ::
     >     lhsp_x(:,:,:,:)
         double precision, allocatable, pinned ::
     >     lhsm_x(:,:,:,:)

         double precision, allocatable, pinned ::
     >     rhs(:,:,:,:)



         double precision, allocatable, pinned ::
     >     u_r(:,:,:,:)
         double precision, allocatable, pinned ::
     >     lhs_x_r(:,:,:,:)
         double precision, allocatable, pinned ::
     >     lhsp_x_r(:,:,:,:)
         double precision, allocatable, pinned ::
     >     lhsm_x_r(:,:,:,:)
         double precision, allocatable, pinned ::
     >     forcing_r(:,:,:,:)
         double precision, allocatable, pinned ::
     >     rhs_r(:,:,:,:)

         double precision, allocatable, pinned ::
     >     lhs_x_j_r(:,:,:,:)
         double precision, allocatable, pinned ::
     >     lhsp_x_j_r(:,:,:,:)
         double precision, allocatable, pinned ::
     >     lhsm_x_j_r(:,:,:,:)
         double precision, allocatable, pinned ::
     >     rhs_j_r(:,:,:,:)


         double precision, allocatable, pinned ::
     >     in_buffer(:), out_buffer(:)
        



        !interface
        ! integer function cudaDeviceSetCacheConfig ( cacheconfig )
        !    integer, intent(in) :: cacheconfig
        ! end function
        ! end interface
       end module glob
