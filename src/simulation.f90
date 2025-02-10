!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,            only: WP
   use geometry,             only: cfg
   use hypre_str_class,      only: hypre_str
   use ddadi_class,          only: ddadi
   use tpns_class,           only: tpns
   use vfs_class,            only: vfs
   use timetracker_class,    only: timetracker
   use ensight_class,        only: ensight
   use surfmesh_class,       only: surfmesh
   use event_class,          only: event
   use monitor_class,        only: monitor
   implicit none
   private

   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),      public :: ps
   type(ddadi),          public :: vs
   type(tpns),           public :: fs
   type(vfs),            public :: vf
   type(timetracker),    public :: time

   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,bubblefile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi

   !> Repulsive force
   real(WP), dimension(:,:,:), allocatable :: Frep_x,Frep_y,Frep_z
   ! !> interface cell positions
   ! real(WP), dimension(:), allocatable :: pos_x,pos_y,pos_z
   ! !> interface cell indexes
   ! integer, dimension(:), allocatable :: idx_x,idx_y,idx_z
   ! !> interface cell normal vectors
   ! real(WP), dimension(:), allocatable :: n_x,n_y,n_z
   ! !> interface cells number
   ! integer :: int_cells

   !> Problem definition
   logical :: moving_domain
   real(WP), dimension(3) :: center1,center2,gravity
   real(WP) :: volume,radius,Ycent,Vrise
   real(WP) :: Vin,Vin_old,Vrise_ref,Ycent_ref,G,ti

   real(WP) :: cent1, cent2, ratio1, ratio2

contains


   !> Function that defines a level set function for a rising bubble problem
   function levelset_rising_bubble(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: d1,d2
      real(WP) :: G
      ! Create the bubble
      d1=-radius+sqrt(sum((xyz-center1)**2))
      d2=-radius+sqrt(sum((xyz-center2)**2))
      G = min(d1,d2)
   end function levelset_rising_bubble


   !> Function that localizes the y+ side of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator


   !> Function that localizes the y- side of the domain
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


   !> Routine that computes rise velocity
   subroutine rise_vel()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,ierr
      real(WP) :: myYcent,myVrise,myvol,bubble_vol
      myVrise=0.0_WP
      myvol=0.0_WP
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               myYcent=myYcent+vf%cfg%ym(j)*(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
               myVrise=myVrise+Vi(i,j,k)*(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
               myvol=myvol+(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(myYcent,Ycent     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myVrise,Vrise     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myvol  ,bubble_vol,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      Ycent=Ycent/bubble_vol
      Vrise=Vrise/bubble_vol
   end subroutine


   !> Subroutine that communicate and get repulsive force
   !> It is only for 2 bubbles
   subroutine get_repulsive_force(Krep, index)
      use mpi_f08
      use parallel
      implicit none
      real(WP), intent(in) :: Krep, index
      !> Local variables
      integer :: i,j,k,ierr,n
      integer :: local_int_cells, global_int_cells
      !> local interface cell positions
      real(WP), dimension(:), allocatable :: pos_x_,pos_y_,pos_z_
      !> local interface cell indexes
      integer, dimension(:), allocatable :: idx_x_,idx_y_,idx_z_
      !> local interface cell normal vectors
      real(WP), dimension(:), allocatable :: n_x_,n_y_,n_z_

      real(WP) :: dist, force
      real(WP), dimension(3) :: n1,n2

      !> look over local domain to know the number of interface cells
      local_int_cells=0
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               if (vf%VF(i,j,k).gt.0.01_WP .and. vf%VF(i,j,k).lt.0.99_WP) then
                  local_int_cells=local_int_cells+1
               end if
            end do
         end do
      end do

      !> allocate arrays
      !> local interface cell infomations, to send to other processors
      allocate(pos_x_(local_int_cells),pos_y_(local_int_cells),pos_z_(local_int_cells))
      allocate(idx_x_(local_int_cells),idx_y_(local_int_cells),idx_z_(local_int_cells))
      allocate(n_x_(local_int_cells),n_y_(local_int_cells),n_z_(local_int_cells))

      !> prepare local interface cell informations
      n=1
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               if (vf%VF(i,j,k).gt.0.01_WP .and. vf%VF(i,j,k).lt.0.99_WP) then
                  pos_x_(n)=vf%cfg%xm(i)
                  pos_y_(n)=vf%cfg%ym(j)
                  pos_z_(n)=vf%cfg%zm(k)
                  idx_x_(n)=i
                  idx_y_(n)=j
                  idx_z_(n)=k
                  n_x_(n)=0.5_WP*(sum(fs%divu_x(:,i,j,k)*vf%VF(i-1:i,j,k)) + sum(fs%divu_x(:,i+1,j,k)*vf%VF(i:i+1,j,k)))
                  n_y_(n)=0.5_WP*(sum(fs%divv_y(:,i,j,k)*vf%VF(i,j-1:j,k)) + sum(fs%divv_y(:,i,j+1,k)*vf%VF(i,j:j+1,k)))
                  n_z_(n)=0.5_WP*(sum(fs%divw_z(:,i,j,k)*vf%VF(i,j,k-1:k)) + sum(fs%divw_z(:,i,j,k+1)*vf%VF(i,j,k:k+1)))
                  n=n+1
               end if
            end do
         end do
      end do

      Frep_x=0.0_WP
      Frep_y=0.0_WP
      Frep_z=0.0_WP

      ! normalize the normal vector
      do i=1, local_int_cells
         n1=[n_x_(i),n_y_(i),n_z_(i)]
         n1 = n1/norm2(n1)
         n_x_(i)=n1(1)
         n_y_(i)=n1(2)
         n_z_(i)=n1(3)
      end do

      !> compute the repulsive force
      do i=1, local_int_cells
         n1=[n_x_(i),n_y_(i),n_z_(i)]
         do j=1, local_int_cells
            if (i .eq. j) cycle

            n2=[n_x_(j),n_y_(j),n_z_(j)]

            if (dot_product(n1,n2) .le. 0.0_WP) then

               dist=sqrt((pos_x_(i)-pos_x_(j))**2 + (pos_y_(i)-pos_y_(j))**2 + (pos_z_(i)-pos_z_(j))**2)
               dist = max(dist, 1.0E-5_WP)
               
               if (dist .lt. 0.2_WP*radius) then
                  force = -1.0_WP*Krep/(dist**index)

                  !> apply the force to the interface cell, in normal direction
                  Frep_x(idx_x_(i),idx_y_(i),idx_z_(i)) = Frep_x(idx_x_(i),idx_y_(i),idx_z_(i)) + force*n1(1)
                  Frep_y(idx_x_(i),idx_y_(i),idx_z_(i)) = Frep_y(idx_x_(i),idx_y_(i),idx_z_(i)) + force*n1(2)
                  Frep_z(idx_x_(i),idx_y_(i),idx_z_(i)) = Frep_z(idx_x_(i),idx_y_(i),idx_z_(i)) + force*n1(3)
               end if

            end if

         end do
      end do

      ! !> convert force to pressure
      ! do k=vf%cfg%kmin_,vf%cfg%kmax_
      !    do j=vf%cfg%jmin_,vf%cfg%jmax_
      !       do i=vf%cfg%imin_,vf%cfg%imax_
      !          if (vf%cfg%vol(i,j,k)*vf%SD(i,j,k) .gt. 0.0_WP) then
      !             Frep_x(i,j,k) = Frep_x(i,j,k)/(vf%cfg%vol(i,j,k)*vf%SD(i,j,k))
      !             Frep_y(i,j,k) = Frep_y(i,j,k)/(vf%cfg%vol(i,j,k)*vf%SD(i,j,k))
      !             Frep_z(i,j,k) = Frep_z(i,j,k)/(vf%cfg%vol(i,j,k)*vf%SD(i,j,k))
      !          end if
      !       end do
      !    end do
      ! end do

      deallocate(pos_x_,pos_y_,pos_z_)
      deallocate(idx_x_,idx_y_,idx_z_)
      deallocate(n_x_,n_y_,n_z_)
   end subroutine

   subroutine get_pos_and_shape()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k, ierr
      real(WP) :: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1
      real(WP) :: xmin2, xmax2, ymin2, ymax2, zmin2, zmax2
      real(WP) :: myzvol1, myzvol2, myvol1, myvol2, zvol1, zvol2, vol1, vol2


      myzvol1 = 0.0_WP; myzvol2 = 0.0_WP; myvol1 = 0.0_WP; myvol2 = 0.0_WP
      xmin1 = vf%cfg%imax_ ; xmax1 = vf%cfg%imin_ ; ymin1 = vf%cfg%jmax_ ; ymax1 = vf%cfg%jmin_ ; zmin1 = vf%cfg%kmax_ ; zmax1 = vf%cfg%kmin_
      xmin2 = vf%cfg%imax_ ; xmax2 = vf%cfg%imin_ ; ymin2 = vf%cfg%jmax_ ; ymax2 = vf%cfg%jmin_ ; zmin2 = vf%cfg%kmax_ ; zmax2 = vf%cfg%kmin_

      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               if (vf%VF(i,j,k).gt.0.01_WP .and. vf%VF(i,j,k).lt.0.99_WP) then
                 if (vf%cfg%zm(k) .lt. 0.0_WP) then
                    myzvol1 = myzvol1 + vf%cfg%zm(k)*vf%cfg%vol(i,j,k)
                    myvol1 = myvol1 + vf%cfg%vol(i,j,k)
                 else
                     myzvol2 = myzvol2 + vf%cfg%zm(k)*vf%cfg%vol(i,j,k)
                    myvol2 = myvol2 + vf%cfg%vol(i,j,k)
                 end if
               end if
            end do
         end do
      end do

      call MPI_ALLREDUCE(myzvol1,zvol1     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myzvol2,zvol2     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myvol1,vol1     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myvol2,vol2     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)

      cent1 = zvol1/vol1
      cent2 = zvol2/vol2

   end subroutine

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none


      ! Check first if we use a moving domain
      call param_read('Moving domain',moving_domain,default=.false.)


      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

         allocate(Frep_x(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Frep_y(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Frep_z(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))

      end block allocate_work_arrays


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker


      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: elvira,VFhi,VFlo,remap
         use mathtools, only: Pi
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=elvira,transport_method=remap,name='VOF')
         !vf%cons_correct=.false.
         ! Initialize a bubble
         call param_read('Bubble position1',center1,default=[0.0_WP,0.0_WP,0.0_WP])
         call param_read('Bubble position2',center2,default=[0.0_WP,0.0_WP,0.0_WP])
         call param_read('Bubble radius',radius)
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_rising_bubble,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof


      ! Prepare PID controller if domain is moving
      if (moving_domain) then
         prepare_controller: block
            ! Store target data
            Ycent_ref=center1(2)
            Vrise_ref=0.0_WP
            ! Controller parameters
            G=0.5_WP
            ti=time%dtmax
         end block prepare_controller
      end if


      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: clipped_neumann,dirichlet
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Assign acceleration of gravity
         call param_read('Gravity',gravity); fs%gravity=gravity
         ! Dirichlet inflow at the top and clipped Neumann outflow at the bottom
         if (moving_domain) then
            call fs%add_bcond(name='inflow' ,type=dirichlet      ,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
            call fs%add_bcond(name='outflow',type=clipped_neumann,face='y',dir=-1,canCorrect=.true. ,locator=ym_locator)
         end if
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Zero inflow velocity
         Vin=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='RisingBubble')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         call ens_out%add_scalar("Frep_x", Frep_x)
         call ens_out%add_scalar("Frep_y", Frep_y)
         call ens_out%add_scalar("Frep_z", Frep_z)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call rise_vel()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create bubble monitor
         bubblefile=monitor(fs%cfg%amRoot,'bubble')
         call bubblefile%add_column(time%n,'Timestep number')
         call bubblefile%add_column(time%t,'Time')
         call bubblefile%add_column(Ycent,'Y centroid')
         call bubblefile%add_column(Vrise,'Rise velocity')
         call bubblefile%add_column(Vin,'Inflow velocity')
         call bubblefile%add_column(cent1, 'bubble1 z')
         call bubblefile%add_column(cent2, 'bubble2 z')
      end block create_monitor


   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: harmonic_visc
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Control inflow condition
         if (moving_domain) then
            control_inflow: block
               use tpns_class, only: bcond
               type(bcond), pointer :: mybc
               integer  :: n,i,j,k
               ! Get new inflow velocity
               Vin_old=Vin
               Vin=G*((Vrise_ref-Vrise)+(Ycent_ref-Ycent)/ti)
               ! Apply inflow at top
               call fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%V(i,j,k)=Vin
               end do
            end block control_inflow
         end if

         ! Remember old VOF
         vf%VFold=vf%VF

         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)

         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)

         ! prepare repulsive force
         call get_repulsive_force(1.0E-9_WP, 3.0_WP)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add momentum source terms - adjust gravity if accelerating frame of reference
            if (moving_domain) fs%gravity(2)=gravity(2)+(Vin-Vin_old)/time%dt
            call fs%addsrc_gravity(resU,resV,resW)

            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW

            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)



            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump_frep(dt=time%dt,div=fs%div,vf=vf,Frep_x=Frep_x,Frep_y=Frep_y,Frep_z=Frep_z)
            ! call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         call get_pos_and_shape()

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call rise_vel()
         call mfile%write()
         call cflfile%write()
         call bubblefile%write()

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker

      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)

   end subroutine simulation_final


end module simulation
