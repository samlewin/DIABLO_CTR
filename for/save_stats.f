

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_CHAN(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      CHARACTER*35 FNAME
      CHARACTER*20 GNAME
      LOGICAL FINAL
      integer i,j,k,n
      real*8 uc, ubulk
    
! This variable is used to add up scalar diagnostics
      real*8 thsum(0:NY+1), thsum_prime(0:NY+1), psum(0:NY+1)
! These variables are used to store and write 2D slices and 3D arrays
      real*8 varxy(0:NXM,1:NY),varzy(0:NZP-1,1:NY),varxz(0:NXM,0:NZP-1)
      real*8 varxyz(0:NXM,0:NZP-1,NY)

! These variable are used for HDF5 writing
      real*8 Diag(1:NY)
      real*8 Diag_Scalar
      real*8 DiagX(0:Nxp - 1)
      real*8 DiagZ(0:TNKZ)

      IF (RANK.EQ.0) 
     &     WRITE(6,*) 'Saving flow statistics.'     

      IF (USE_MPI) THEN
        call mpi_barrier(MPI_COMM_WORLD,ierror)
        CALL GHOST_CHAN_MPI
      END IF

C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

      if (FINAL) then
! We are done with the simulation
        
#ifdef HDF5
        FNAME='stats.h5'
        if (USE_MPI) then
          call mpi_barrier(MPI_COMM_WORLD,ierror)
        end if
        IF (RANKZ.EQ.0) THEN
          Diag=GYF(1:NY)
          gname='GYF'
          call WriteStatH5(FNAME,gname,Diag)
          Diag=UBAR(1:NY)
          gname='UBAR'
          call WriteStatH5(FNAME,gname,Diag)
          Diag=VBAR(1:NY)
          gname='VBAR'
          call WriteStatH5(FNAME,gname,Diag)
          Diag=WBAR(1:NY)
          gname='WBAR'
          call WriteStatH5(FNAME,gname,Diag)
          do n=1,N_TH
            Diag=THBAR(1:NY,n)
            gname='THBAR'
     &           //CHAR(MOD(N,100)/10+48)
     &           //CHAR(MOD(N,10)+48)
            call WriteStatH5(FNAME,gname,Diag)
          end do
        END IF

#else
! Here we aren't using HDF5, so save to text files
        IF (RANKZ.EQ.0) THEN
        IF (USE_MPI) THEN
          FNAME='stats'//trim(MPI_IO_NUM)//'.txt'
        ELSE
          FNAME='stats.txt'
        END IF
        open(20,file=FNAME,form='formatted',status='unknown')
        do j=1,NY
          write(20,201) j,GYF(j),UBAR(j),VBAR(j),WBAR(j)
        end do
201     format(I3,',',F16.9,',',F16.9,',',F16.9,',',F16.9)
        do n=1,N_TH
        do j=1,NY
          write(20,202) j,GYF(j),THBAR(j,n)
        end do
        end do
202     format(I3,',',F16.9,',',F16.9)
        close(20)
        END IF
#endif

      else     

! Compute and write out the centerline velocity
      IF (NPROCY.EQ.1) THEN
      if (int(float(NY)/2.) .eq. float(NY)/2.) then
! IF NY is even
        uc=dble(CU1(0,0,int(float(NY)/2.))) 
      else
        uc=0.5*(dble(CU1(0,0,int(float(NY)/2.)-1))
     +         +dble(CU1(0,0,int(float(NY)/2.))))
      end if
      write(*,*) 'Centerline velocity = ', uc 
! Compute and write out bulk velocity
      END IF

! We are in the middle of a run, compile statistics
! First get the number of samples taken so far
      IF (RANK.EQ.0) write(*,*) 'TIME, DELTA_T: ',TIME, DELTA_T
      IF (RANKZ.EQ.0) THEN
         NSAMPLES=NSAMPLES+1
! Get the mean velocity
         do j=1,NY
            UBAR(j)=(1./float(NSAMPLES))*dble(CU1(0,0,j))
     &           +((float(NSAMPLES)-1.)/float(NSAMPLES))*UBAR(j)
            VBAR(j)=(1./float(NSAMPLES))*dble(CU2(0,0,j))
     &           +((float(NSAMPLES)-1.)/float(NSAMPLES))*VBAR(j)
            WBAR(j)=(1./float(NSAMPLES))*dble(CU3(0,0,j))
     &           +((float(NSAMPLES)-1.)/float(NSAMPLES))*WBAR(j)
            do n=1,N_TH
               THBAR(j,n)=(1./float(NSAMPLES))*dble(CTH(0,0,j,n))
     &         +((float(NSAMPLES)-1.)/float(NSAMPLES))*THBAR(j,n)
            end do
         end do

! Integrate the instantaneous mean profile numerically at GY points
         UME=CU1(0,0,:)
      ELSE
         UME=0.d0
      END IF
      CALL INTEGRATE_Y_VAR(UME,UBULK,MPI_COMM_WORLD)
! Write out UBULK
      IF (RANK.EQ.0) write(*,*) 'UBULK: ',UBULK
      IF (RANK.EQ.0) write(*,*) 'Umin: ',UME(2)

! Save CUi
      do k=0,TNKZ
        do i=0,NXP-1 ! NKX
          do j=0,NY+1
            CR1(i,k,j)=CU1(i,k,j)
            CR2(i,k,j)=CU2(i,k,j)
            CR3(i,k,j)=CU3(i,k,j)
          end do
        end do
      end do 

! Get the mean value of the velocities
      IF (RANKZ.EQ.0) THEN
         ume=dble(CU1(0,0,:))
         vme=dble(CU2(0,0,:))
         wme=dble(CU3(0,0,:)) 
         pme=dble(CP(0,0,:))
         DO n=1,N_TH
            thme(:,n)=dble(CTH(0,0,:,n))
         END DO
      END IF
      CALL MPI_BCAST(ume,NY+2,MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_Z,ierror)
      CALL MPI_BCAST(vme,NY+2,MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_Z,ierror)
      CALL MPI_BCAST(wme,NY+2,MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_Z,ierror)
      IF (N_TH.GT.0) CALL MPI_BCAST(thme,(NY+2)*N_TH,
     &     MPI_DOUBLE_PRECISION,0,MPI_COMM_Z,ierror)

! Convert to physical space
      call fft_xz_to_physical(CU1,U1,0,NY+1)
      call fft_xz_to_physical(CU2,U2,0,NY+1)
      call fft_xz_to_physical(CU3,U3,0,NY+1)

! Calculate the dissipation rate
      call tkebudget_chan
      ! call calcbpe
      IF (LES) call tkebudget_chan_les

! Get the total kinetic energy at each level 
      do j=0,NY
        u_ke(j)=0.
        v_ke(j)=0.
        w_ke(j)=0.
        ke_tot(j)=0.
      do k=0,NZP-1
      do i=0,NXM 
        u_ke(j)=u_ke(j)+(U1(i,k,j))**2.
        v_ke(j)=v_ke(j)+(0.5*((U2(i,k,j  )) +
     &                       (U2(i,k,j+1))))**2.
        w_ke(j)=w_ke(j)+(U3(i,k,j))**2.
      end do
      end do
      end do
      

      call mpi_allreduce(mpi_in_place,u_ke,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,v_ke,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,w_ke,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      
      do j=0,NY
        u_ke(j)=u_ke(j)/(float(NZ)*float(NX))
        v_ke(j)=v_ke(j)/(float(NZ)*float(NX))
        w_ke(j)=w_ke(j)/(float(NZ)*float(NX))
        ke_tot(j)=0.5*(u_ke(j) + v_ke(j) + w_ke(j))
      end do 
      
      CALL INTEGRATE_Y_VAR(ke_tot,ke_tot_bulk,MPI_COMM_Y)

! calculate total mean KE
      do j=0, NY
      ke_mean_tot(j)=0
      end do
      do j=0, NY
      ke_mean_tot(j) = ke_mean_tot(j) + 
     &                 0.5*((ume(j))**2. + (vme(j))**2. + (wme(j))**2.)
      end do
      CALL INTEGRATE_Y_VAR(ke_mean_tot,ke_mean_tot_bulk,MPI_COMM_Y)

! Get the turbulent kinetic energy at each level 
      do j=0,NY
        urms(j)=0.
        wrms(j)=0.
      do k=0,NZP-1
      do i=0,NXM 
        urms(j)=urms(j)+(U1(i,k,j)-ume(j))**2.
        wrms(j)=wrms(j)+(U3(i,k,j)-wme(j))**2.
      end do
      end do
      end do

      do j=1,NY
        vrms(j)=0.
      do k=0,NZP-1
      do i=0,NXM 
        vrms(j)=vrms(j)+(0.5*((U2(i,k,j  )-vme(j  )) +
     &                       (U2(i,k,j+1)-vme(j+1)) ))**2
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,urms,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,vrms,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,wrms,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)

      do j=0,NY
        urms(j)=sqrt(urms(j)/(float(NZ)*float(NX)))
        vrms(j)=sqrt(vrms(j)/(float(NZ)*float(NX)))
        wrms(j)=sqrt(wrms(j)/(float(NZ)*float(NX)))
      end do 

      ! Get the bulk rms value
      CALL INTEGRATE_Y_VAR(urms,urms_b,MPI_COMM_Y)
      CALL INTEGRATE_Y_VAR(vrms,vrms_b,MPI_COMM_Y)
      CALL INTEGRATE_Y_VAR(wrms,wrms_b,MPI_COMM_Y)

! Compute the Reynolds stress and mean velocity gradient
! Here, uv and wv are defined on the GY grid
! uw is defined on the GYF grid
      do j=1,NY
        uv(j)=0. 
        uw(j)=0.
        wv(j)=0.
      do k=0,NZP-1
      do i=0,NXM
        uv(j)=uv(j)+
     &     (DYF(j-1)*(U1(i,k,j)-ume(j))+DYF(j)*(U1(i,k,j-1)-ume(j-1)))
     &             /(2.d0*DY(J))
     &     *0.5*(U2(i,k,j)-vme(j)+U2(i,k,j+1)-vme(j+1))
        wv(j)=wv(j)+
     &     (DYF(j-1)*(U3(i,k,j)-wme(j))+DYF(j)*(U3(i,k,j-1)-wme(j-1)))
     &             /(2.d0*DY(J))
     &     *0.5*(U2(i,k,j)-vme(j)+U2(i,k,j+1)-vme(j+1))
        uw(j)=uw(j)+(U1(i,k,j)-ume(j))
     +    *(U3(i,k,j)-wme(j))
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,uv,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,uw,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,wv,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      
      do j=1,NY
        uv(j)=uv(j)/(float(NZ)*float(NX))
        uw(j)=uw(j)/(float(NZ)*float(NX))
        wv(j)=wv(j)/(float(NZ)*float(NX))
      end do
              
! Get the derivatives of the mean velocity at GY points
      do j=1,NY
        dudy(j)=(ume(j)-ume(j-1))/DY(j)
        dwdy(j)=(wme(j)-wme(j-1))/DY(j)
      end do  
! dudx
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 !NKX
        CS1(i,k,j)= CIKX(i)*CR1(i,k,j)
      end do
      end do
      end do

      call fft_xz_to_physical(CS1,S1,0,NY+1)

      do j=1,NY
        dudx(j)=0.d0
      do k=0,NZP-1
      do i=0,NXM
        dudx(j)=dudx(j)+S1(i,k,j)
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,dudx,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      dudx(j)=dudx(j)/(dble(NX)*dble(NZ))
      end do
! dudz
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 !NKX
        CS1(i,k,j)= CIKZ(K)*CR1(i,k,j)
      end do
      end do
      end do

      call fft_xz_to_physical(CS1,S1,0,NY+1)

      do j=1,NY
        dudz(j)=0.d0
        shear(j)=0.d0
      do k=0,NZP-1
      do i=0,NXM
        dudz(j)=dudz(j)+S1(i,k,j)
        shear(j)=shear(j)+S1(i,k,j)**2
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,dudz,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      dudz(j)=dudz(j)/(dble(NX)*dble(NZ))
      end do

! dvdz
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
         CS1(i,k,j)=CIKZ(k)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
      end do
      end do
      end do

      call fft_xz_to_physical(CS1,S1,0,NY+1)

      do j=1,NY
        dvdz(j)=0.d0
      do k=0,NZP-1
      do i=0,NXM
        dvdz(j)=dvdz(j)+S1(i,k,j)
        shear(j)=shear(j)+S1(i,j,k)**2
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,dudz,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,shear,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      dvdz(j)=dvdz(j)/(dble(NX)*dble(NZ))
      shear(j)=shear(j)/(dble(NX)*dble(NZ))
      end do

! dwdx
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 !NKX
        CS1(i,k,j)= CIKX(K)*CR3(i,k,j)
      end do
      end do
      end do

      call fft_xz_to_physical(CS1,S1,0,NY+1)

      do j=1,NY
        dwdx(j)=0.d0
      do k=0,NZP-1
      do i=0,NXM
        dwdx(j)=dwdx(j)+S1(i,k,j)
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,dwdx,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      dwdx(j)=dwdx(j)/(dble(NX)*dble(NZ))
      end do

! dwdx
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 !NKX
        CS1(i,k,j)= CIKZ(K)*CR3(i,k,j)
      end do
      end do
      end do

      call fft_xz_to_physical(CS1,S1,0,NY+1)

      do j=1,NY
        dwdz(j)=0.d0
      do k=0,NZP-1
      do i=0,NXM
        dwdz(j)=dwdz(j)+S1(i,k,j)
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,dwdz,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      dwdz(j)=dwdz(j)/(dble(NX)*dble(NZ))
      end do

! Write out the bulk rms velocity
      if (RANK.eq.0) then
         write(*,*) '<U_rms>: ',urms_b
         write(*,*) '<V_rms>: ',vrms_b
         write(*,*) '<W_rms>: ',wrms_b

      end if

! Get the rms vorticity
! First, get the x-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 !NKX
        CS1(i,k,j)=(CR3(i,k,j+1)-CR3(i,k,j-1))/(2.d0*DYF(j))
     &            -CIKZ(K)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
      do j=1,NY
      omega_x(j)=0.d0
      do k=0,NZP-1
      do i=0,NXM
        omega_x(j)=omega_x(j)+S1(i,k,j)**2.d0
      end do
      end do
      end do
      call mpi_allreduce(mpi_in_place,omega_x,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      omega_x(j)=sqrt(omega_x(j)/(dble(NX)*dble(NZ)))
      end do

! Now, get the y-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 !NKX
        CS1(i,k,j)=CIKZ(k)*CR1(i,k,j)-CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
      do j=1,NY
      omega_y(j)=0.d0
      do k=0,NZP-1
      do i=0,NXM
        omega_y(j)=omega_y(j)+S1(i,k,j)**2.d0
      end do
      end do
      end do
      call mpi_allreduce(mpi_in_place,omega_y,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      omega_y(j)=sqrt(omega_y(j)/(dble(NX)*dble(NZ)))
      end do

! Now, get the z-component in fourier space
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 ! NKX
        CS1(i,k,j)=CIKX(i)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
     &             -(CR1(i,k,j+1)-CR1(i,k,j-1))/(2.d0*DYF(j))
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
! Get the rms value
      do j=1,NY
      omega_z(j)=0.d0
      do k=0,NZP-1
      do i=0,NXM
        omega_z(j)=omega_z(j)+S1(i,k,j)**2.d0
      end do
      end do
      end do
      call mpi_allreduce(mpi_in_place,omega_z,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      omega_z(j)=sqrt(omega_z(j)/(dble(NX)*dble(NZ)))
      end do
      
! Compute the 3D Spanwise Vorticity omega_z_xy
      F6(:,:,:) = 0.d0
! Compute the relevant component in Fourier Space 
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 ! NKX
        CS1(i,k,j)=CIKX(i)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
     &             -(CR1(i,k,j+1)-CR1(i,k,j-1))/(2.d0*DYF(j))

      end do
      end do
      end do

! Convert to physical space
      call fft_xz_to_physical(CS1,F6,0,NY+1) 
! Compute vertical vorticity
      F5(:,:,:) = 0.d0
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1 ! NKX
        CS1(i,k,j)=CIKZ(k)*CR1(i,k,j)-CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
      call fft_xz_to_physical(CS1,F5,0,NY+1) 

! Get 3D KE equation terms 
 
      do i=0,NXM
      do j=0,NY+1
      S2(i,:,j) = 0.d0
      S1(i,:,j) = 0.d0
      R1(i,:,j) = 0.d0
      do k=0,NZP-1   
        S1(i,0,j) = S1(i,0,j) + (U1(i,k,j)-ume(j))/(float(NZ)) !u_kh
        S2(i,0,j) = S2(i,0,j) + U2(i,k,j)/(float(NZ)) !v_kh
        R1(i,0,j) = R1(i,0,j) + (TH(i,k,j,1)-thme(j,1))/(float(NZ))
      end do
      end do
      k=0
      call mpi_allreduce(mpi_in_place,S1(i,k,:),NY+2,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,S2(i,k,:),NY+2,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,R1(i,k,:),NY+2,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      end do

      do j=0,NY+1
      do i=0, NXM
      do k=0, NZP-1
        S1(i,k,j) = S1(i,0,j)
        S2(i,k,j) = S2(i,0,j)
        R3(i,k,j) = S1(i,0,j)
        F1(i,k,j) = U1(i,k,j) - S1(i,0,j) - ume(j)
        F2(i,k,j) = U2(i,k,j) - S2(i,0,j)
        F3(i,k,j) = U3(i,k,j)
        F4(i,k,j) = TH(i,k,j,1) - R1(i,0,j) - thme(j,1)
      end do
      end do
      end do

      do j=1,NY
        u_3d(j) = 0.d0
        v_3d(j) = 0.d0
        w_3d(j) = 0.d0
        u_kh(j) = 0.d0
        v_kh(j) = 0.d0
      do i=0,NXM 
      k=0
        u_kh(j) = u_kh(j) + S1(i,k,j)**2.
        v_kh(j) = v_kh(j) + (0.5*(S2(i,k,j)+S2(i,k,j+1)))**2.
      do k=0,NZP-1
        u_3d(j) = u_3d(j) + F1(i,k,j)**2
        v_3d(j) = v_3d(j) + (0.5*(F2(i,k,j)+F2(i,k,j+1)))**2.
        w_3d(j) = w_3d(j) + F3(i,k,j)**2
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,u_3d,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,v_3d,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,w_3d,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      
      do j=1,NY
        u_kh(j)=u_kh(j)/(float(NX))
        v_kh(j)=v_kh(j)/(float(NX))
        u_3d(j)=u_3d(j)/(float(NZ)*float(NX))
        v_3d(j)=v_3d(j)/(float(NZ)*float(NX))
        w_3d(j)=w_3d(j)/(float(NZ)*float(NX))       
        ke_kh(j)=0.5*(u_kh(j) + v_kh(j))
        ke_3d(j)=0.5*(u_3d(j)+v_3d(j)+w_3d(j))
      end do 

      !compute du_kh/dx
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      do j=1,NY
      do i=0,NXP-1 !NKX
      do k=0, TNKZ
        CS1(i,k,j)= CIKX(i)*CS1(i,k,j)
      end do
      end do
      end do
      CALL FFT_XZ_TO_PHYSICAL(CS1,R1,0,NY+1)     

      !compute dv_kh/dy
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
      R2(i,k,j)=((S2(i,k,j+1))-(S2(i,k,j-1)))
     &            /(GY(j+1)-GY(j-1))
      end do
      end do
      end do

      do j=2,NY
        bckgrd_ext(j)=0. 
        thv_3d(j)=0.
        tilt_ext(j)=0.
      do k=0,NZP-1
      do i=0,NXM
        bckgrd_ext(j) = bckgrd_ext(j)+
     &     (DYF(j-1)*F1(i,k,j) + DYF(j)*F1(i,k,j-1))
     &             /(2.d0*DY(J))
     &     *0.5*(F2(i,k,j)+F2(i,k,j+1))*dudy(j)
        thv_3d(j) = thv_3d(j)+
     &   (DYF(j-1)*F4(i,k,j) + DYF(j)*F4(i,k,j-1))
     &             /(2.d0*DY(J))*
     &   0.5*(F2(i,k,j)+F2(i,k,j+1))
        tilt_ext(j) = tilt_ext(j) + 
     &     (((DYF(j-1)*F1(i,k,j) + DYF(j)*F1(i,k,j-1))
     &     /(2.d0*DY(J)))**2-(0.5*(F2(i,k,j)+F2(i,k,j+1)))**2)*
     &     ((( DYF(j-1)*R1(i,k,j) + DYF(j)*R1(i,k,j-1))
     &             /(2.d0*DY(J)))-R2(i,k,j))
      end do 
      end do  
      end do

      !compute dv_kh/dx
      CALL FFT_XZ_TO_FOURIER(S2,CS1,0,NY+1)
      do j=1,NY
      do i=0,NXP-1 !NKX
      do k=0, TNKZ
        CS1(i,k,j)= CIKX(i)*0.5*(CS1(i,k,j)+CS1(i,k,j+1))
      end do
      end do
      end do
      CALL FFT_XZ_TO_PHYSICAL(CS1,R1,0,NY+1)     

      !compute du_kh/dy
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
      R2(i,k,j)=((R3(i,k,j))
     &           -(R3(i,k,j-1)))
     &            /DY(j)
      end do
      end do
      end do

      do j=2,NY
        kh_ext(j)=0.
      do k=0,NZP-1
      do i=0,NXM
      kh_ext(j) = kh_ext(j) + (R1(i,k,j) + 
     &     R2(i,k,j))*
     &     0.5*(F2(i,k,j)+F2(i,k,j+1))*
     &   (DYF(j-1)*F1(i,k,j) + DYF(j)*F1(i,k,j-1))
     &             /(2.d0*DY(J))
      end do
      end do
      end do

      call mpi_allreduce(mpi_in_place,bckgrd_ext,NY+2,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,thv_3d,NY+2,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place, kh_ext,NY+2,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place, tilt_ext,NY+2,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)

      do j=1,NY
      bckgrd_ext(j)=bckgrd_ext(j)/(float(NZ)*float(NX))
      thv_3d(j)=thv_3d(j)/(float(NZ)*float(NX))
      kh_ext(j)=kh_ext(j)/(float(NZ)*float(NX))
      tilt_ext(j)=0.5*tilt_ext(j)/(float(NZ)*float(NX))
      end do
      
      call tkebudget_chan_3d

      IF (RANK.EQ.0) 
     &     write(*,*) 'done computing velocity stats'


#ifdef HDF5
      FNAME='mean.h5'

      gname='time'
      call WriteHDF5_real(FNAME,gname,TIME)

      gname='ke_tot_bulk'
      call WriteHDF5_real(FNAME, gname, ke_tot_bulk)

      gname='ke_mean_tot_bulk'
      call WriteHDF5_real(FNAME, gname, ke_mean_tot_bulk)

      ! IF (MOD(TIME_STEP,SAVE_STATS_INT*5).EQ.0) THEN
      !  gname='time_red'
      !  call WriteHDF5_real(FNAME,gname,TIME)
      ! END IF
       
      IF (RANKZ.eq.0) then

      gname='gyf'
      Diag=gyf(1:NY)      
      call WriteStatH5(FNAME,gname,Diag)

      gname='ume'
      Diag=ume(1:NY)
      call WriteStatH5(FNAME,gname,Diag)
    
      gname='vme'
      Diag=0.5*(vme(1:NY)+vme(2:NY+1))
      call WriteStatH5(FNAME,gname,Diag)

      gname='wme'
      Diag=wme(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='ke_mean_tot'
      Diag=ke_mean_tot(1:NY)
      call writeStatH5(FNAME,gname,Diag)

      gname='urms'
      Diag=urms(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='vrms'
      Diag=vrms(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='wrms'
      Diag=wrms(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='uv'
      Diag=uv(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='uw'
      Diag=uw(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='wv'
      Diag=wv(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='dudx'
      Diag=dudx(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='dudz'
      Diag=dudz(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='dvdz'
      Diag=dvdz(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='dwdx'
      Diag=dwdx(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='dwdz'
      Diag=dwdz(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='dudy'
      Diag=dudy(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='dwdy'
      Diag=dwdy(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='ke_3d'
      Diag=ke_3d(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='ke_kh'
      Diag=ke_kh(1:NY)
      call WriteStatH5(FNAME, gname, Diag)

      gname='bckgrd_ext'
      Diag=bckgrd_ext(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='thv_3d'
      Diag=thv_3d(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='tilt_ext'
      Diag=tilt_ext(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='kh_ext'
      Diag=kh_ext(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='cp'
      Diag=dble(cp(0,0,1:NY))
      call WriteStatH5(FNAME,gname,Diag)

      gname='pv'
      Diag=dble(pv(1:NY))
      call WriteStatH5(FNAME,gname,Diag)

      gname='shear'
      Diag=shear(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='omega_x'
      Diag=omega_x(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='omega_y'
      Diag=omega_y(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='omega_z'
      Diag=omega_z(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      END IF

      IF (RANK.EQ.0) 
     &     write(*,*) 'done writing mean part 1'

#else
! Here are aren't using HDF5, so write mean statistcs to text files
! Write out the mean statistics at each time
      IF (RANKZ.EQ.0) THEN
      IF (USE_MPI) THEN
        FNAME='mean'//trim(MPI_IO_NUM)//'.txt'
      ELSE
        FNAME='mean.txt'
      END IF
      open(40,file=FNAME,form='formatted',status='unknown')
      write(40,*) TIME_STEP,TIME,DELTA_T
      write(40,*) UBULK
      do j=1,NY
        write(40,401) j,GYF(J),ume(j)
     +      ,0.5*(vme(j+1)+vme(j))
     +      ,wme(j),urms(j),vrms(j),wrms(j)
     +      ,uv(j),uw(j),wv(j),dudy(j),dwdy(j),dble(cp(0,0,j)),shear(j)
     &      ,omega_x(j),omega_y(j),omega_z(j)
      end do
      END IF
401   format(I3,' ',17(F30.20,' '))
#endif


! Do over the number of passive scalars
      do n=1,N_TH

! Save CTH
      do k=0,TNKZ
        do i=0,NXP-1 ! NKX
          do j=0,NY+1
            CRTH(i,k,j,n)=CTH(i,k,j,n)
          end do
        end do
      end do

! Compute the scalar gradient and store in CRi
      do j=1,NY
        do k=0,TNKZ
          do i=0,NXP-1 ! NKX
! Store gradients of TH(:,:,:,n) (if it is used) in CRi
          CR1(i,k,j)=CIKX(i)*CTH(i,k,j,n)
          CR2(i,k,j)=(CTH(i,k,j+1,n)-CTH(i,k,j-1,n))/(GYF(j+1)-GYF(j-1))
          CR3(i,k,j)=CIKZ(k)*CTH(i,k,j,n)
          end do
        end do
      end do
! Convert gradients to physical space
      CALL FFT_XZ_TO_PHYSICAL(CR1,R1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CR2,R2,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CR3,R3,0,NY+1)

! Convert to physical space

      call mpi_barrier(MPI_COMM_WORLD,ierror)

      CS1(:,:,:)=CTH(:,:,:,N)
      CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
      TH(:,:,:,N)=S1(:,:,:)

      do j=1,NY
        thsum(j)=0.
      do k=0,NZP-1
      do i=0,NXM
        thsum(j)=thsum(j)+(abs(TH(i,k,j,n)-thme(j,n)))**2.
      end do
      end do
      end do
      call mpi_allreduce(mpi_in_place,thsum,(NY+2),
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
        thrms(j,n)=sqrt(thsum(j)/(float(NZ)*float(NX)))
      end do
! Compute the Reynolds stress and mean velocity gradient
      do j=1,NY
        thsum(j)=0.
      do k=0,NZP-1
      do i=0,NXM
       thsum(j)=thsum(j)+(TH(i,k,j,n)
     +    *U3(i,k,j))
      end do
      end do
      end do
      call mpi_allreduce(mpi_in_place,thsum,(NY+2),
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      thv(j,n)=thsum(j)/(float(NZ)*float(NX))
      end do

      CS1(:,:,:)=CP(:,:,:)
      CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)

      do j=1,NY
        psum(j)=0.
      do k=0,NZP-1
      do i=0,NXM
       psum(j)=psum(j)+(S1(i,k,j)-pme(j))
     +    *(0.5*(U2(i,k,j)+U2(i,k,j+1))
     &      -0.5*(vme(j)+vme(j+1)))
      end do
      end do
      end do
      call mpi_allreduce(mpi_in_place,psum,(NY+2),
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
      pv(j)=psum(j)/(float(NZ)*float(NX))
      end do

! Get the y-derivative of the mean scalar at GYF points
      do j=1,NY
        dthdy(j,n)=(thme(j+1,n)-thme(j-1,n))/(GYF(j+1)-GYF(j-1))
      end do

! Compute the potential energy dissipation, grad(TH) \cdot grad(TH)
      F2(:,:,:)=0.d0
      F3(:,:,:)=0.d0
      do j=1,NY
        thsum(j)=0.d0
        thsum_prime(j)=0.d0
        do k=0,NZP-1
          do i=0,NXM
            F2(i,k,j)=F2(i,k,j)
     &       + R1(i,k,j)**2.d0+R2(i,k,j)**2.d0+R3(i,k,j)**2.d0
            F3(i,k,j)=F3(i,k,j)
     &           +R1(i,k,j)**2.d0
     &           +(abs(R2(i,k,j)-dthdy(j,n)))**2.d0 
     &           +R3(i,k,j)**2.d0
            thsum(j)=thsum(j)
     &          +R1(i,k,j)**2.d0+R2(i,k,j)**2.d0+R3(i,k,j)**2.d0
            thsum_prime(j)=thsum_prime(j)
     &           +R1(i,k,j)**2.d0
     &           +(R2(i,k,j)-dthdy(j,n))**2.d0 
     &           +R3(i,k,j)**2.d0
          end do
        end do
      end do
      call mpi_allreduce(mpi_in_place,thsum,(NY+2),
     &       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,thsum_prime,(NY+2),
     &       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
      do j=1,NY
        pe_diss(j,n)=thsum(j)/dble(NX*NZ)
        pe_diss_prime(j,n)=thsum_prime(j)/dble(NX*NZ)
      end do

      IF (RANK.EQ.0) 
     &     write(*,*) 'done doing theta stats'

    !   call bpebudget_chan
    !   IF (RANK.EQ.0) 
    !  &     write(*,*) 'done computing bpe'
      ! call bpebudget_chan_2

#ifdef HDF5 
      if (MOVIE) then
         FNAME='thfield.h5'
        if (n.eq.1) then    
        IF (MOD(TIME_STEP,SAVE_STATS_INT*5).EQ.0) THEN
        if (USE_MPI) then
        call mpi_barrier(MPI_COMM_WORLD,ierror)
        end if
        do I=0,NXM
        do J=1,NY
        do K=0,NZP-1
            varxyz(i,k,j)=THMEZ(i,k,j)
        end do
        end do
        end do
        GNAME='th1_xyz'
        ! call writeHDF5_xyz(FNAME,GNAME,varxyz)
         END IF
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=TH(i,NzMovie,j,n)
            end do
            end do
            GNAME='th1_xy'
            FNAME = 'movie.h5'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
        !  END IF
         END IF 
         if (USE_MPI) then
          call mpi_barrier(MPI_COMM_WORLD,ierror)
          end if
          IF (RANKZ.EQ.RANKZMOVIE) THEN
              do I=0,NXM
              do J=1,NY
                varxy(i,j)=F6(i,NzMovie,j)
              end do
              end do
              GNAME='omega_z_xy'
              FNAME='movie.h5'
              call writeHDF5_xyplane(FNAME,GNAME,varxy)
          END IF 
          if (USE_MPI) then
          call mpi_barrier(MPI_COMM_WORLD,ierror)
          end if
          IF (RANKZ.EQ.RANKZMOVIE) THEN
              do I=0,NXM
              do J=1,NY
                varxy(i,j)=F5(i,NzMovie,j)
              end do
              end do
              GNAME='omega_y_xy'
              FNAME='movie.h5'
              call writeHDF5_xyplane(FNAME,GNAME,varxy)
          END IF 
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=F2(i,NzMovie,j)
            end do
            end do
            GNAME='pediss_xy'
            FNAME = 'movie.h5'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
         END IF
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=F3(i,NzMovie,j)
            end do
            end do
            GNAME='pediss_prime_xy'
            FNAME = 'movie.h5'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
         END IF
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
               varxz(i,j)=TH(i,j,NyMovie,n)
            end do
            end do
            GNAME='th1_xz'
            FNAME='movie.h5'
            call writeHDF5_xzplane(FNAME,GNAME,varxz)
         END IF
         if (USE_MPI) then
          call mpi_barrier(MPI_COMM_WORLD,ierror)
          end if
          IF (RANKY.EQ.RANKYMOVIE) THEN
             do I=0,NXM
             do J=0,NZP-1
                varxz(i,j)=F3(i,j,NyMovie)
             end do
             end do
             GNAME='pediss_prime_xz'
             FNAME='movie.h5'
             call writeHDF5_xzplane(FNAME,GNAME,varxz)
          END IF
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
               varxz(i,j)=F6(i,j,NyMovie)
            end do
            end do
            GNAME='omega_z_xz'
            FNAME='movie.h5'
            call writeHDF5_xzplane(FNAME,GNAME,varxz)
         END IF
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
               varxz(i,j)=F6(i,j,NyMovie)
            end do
            end do
            GNAME='omega_y_xz'
            FNAME='movie.h5'
            call writeHDF5_xzplane(FNAME,GNAME,varxz)
         END IF
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=TH(NxMovie,i,j,n)
         end do
         end do
         GNAME='th1_zy'
         FNAME='movie.h5'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)
         if (USE_MPI) then
          call mpi_barrier(MPI_COMM_WORLD,ierror)
          end if
          do I=0,NZP-1
          do J=1,NY
              varzy(i,j)=F3(NxMovie,i,j)
          end do
          end do
          GNAME='pediss_prime_zy'
          FNAME='movie.h5'
          call writeHDF5_zyplane(FNAME,GNAME,varzy)
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=F6(NxMovie,i,j)
         end do
         end do
         GNAME='omega_z_zy'
         FNAME='movie.h5'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=F5(NxMovie,i,j)
         end do
         end do
         GNAME='omega_y_zy'
         FNAME='movie.h5'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)
      else if (n.eq.2) then
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=TH(i,NzMovie,j,n)
            end do
            end do
            GNAME='th2_xy'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
         END IF
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
               varxz(i,j)=TH(i,j,NyMovie,n)
            end do
            end do
            GNAME='th2_xz'
            call writeHDF5_xzplane(FNAME,GNAME,varxz)
         END IF
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=TH(NxMovie,i,j,n)
         end do
         end do
         GNAME='th2_zy'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)
      end if

      END IF

      IF (RANK.EQ.0) 
     &     write(*,*) 'done theta movie bit'
#endif

! Convert back to Fourier space
      S1(:,:,:)=TH(:,:,:,N)
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      CTH(:,:,:,N)=CS1(:,:,:)

! End do over number of passive scalars, n
      end do


#ifdef HDF5
 
      FNAME='mean.h5' 
  
      IF (RANKZ.eq.0) THEN

     
      do n=1,N_TH
 
        Diag=thme(1:NY,n)
        gname='thme'
     &           //CHAR(MOD(N,100)/10+48)
     &           //CHAR(MOD(N,10)+48)
        call WriteStatH5(FNAME,gname,Diag)

        Diag=dthdy(1:NY,n)
        gname='dthdy'
     &           //CHAR(MOD(N,100)/10+48)
     &           //CHAR(MOD(N,10)+48)
        call WriteStatH5(FNAME,gname,Diag)

        Diag=thrms(1:NY,n)
        gname='thrms'
     &           //CHAR(MOD(N,100)/10+48)
     &           //CHAR(MOD(N,10)+48)
        call WriteStatH5(FNAME,gname,Diag)

        Diag=thv(1:NY,n)
        gname='thv'
     &           //CHAR(MOD(N,100)/10+48)
     &           //CHAR(MOD(N,10)+48)
        call WriteStatH5(FNAME,gname,Diag)

        Diag=pe_diss(1:NY,n)
        gname='pe_diss'
     &           //CHAR(MOD(N,100)/10+48)
     &           //CHAR(MOD(N,10)+48)
        call WriteStatH5(FNAME,gname,Diag)

        Diag=pe_diss_prime(1:NY,n)
        gname='pe_diss_prime'
     &           //CHAR(MOD(N,100)/10+48)
     &           //CHAR(MOD(N,10)+48)
        call WriteStatH5(FNAME,gname,Diag)
 
      end do

      END IF

      IF (RANK.EQ.0) 
     &     write(*,*) 'done theta part of mean.h5'

#else
! Here we aren't using HDF5, so write to a text file
! Write out the mean statistics at each time
      IF (RANKZ.EQ.0) THEN
      IF (USE_MPI) THEN
        FNAME='mean_th'//trim(MPI_IO_NUM)//'.txt'
      ELSE
        FNAME='mean_th.txt'
      END IF
      open(41,file=FNAME,form='formatted',status='unknown')
      write(41,*) TIME_STEP,TIME,DELTA_T
      write(41,*) UBULK
      do n=1,N_TH 
      do j=1,NY
        write(41,402) j,GYF(J),thme(j,n)
     +      ,dthdy(j,n),thrms(j,n),thv(j,n),pe_diss(j,n)
      end do
      end do
      END IF
402   format(I3,' ',6(F30.20,' '))
#endif

      IF (RANK.EQ.0) 
     &     write(*,*) 'VERBOSITY: ',VERBOSITY
      if (VERBOSITY.gt.4) then 
      IF (RANK.EQ.0) 
     &        write(*,*) 'Outputting info for gnuplot...'
      open (unit=10, file="solution")
      do i=2,NXM
        do j=2,NYM
          write (10,*) i, j, U1(i,0,j)
        end do
        write (10,*) ""
      end do
      close (10)
      call system ('gnuplot <gnuplot.in') 
      end if

#ifdef HDF5
      if (MOVIE) then
         FNAME='movie.h5'
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=U1(i,NzMovie,j)
            end do
            end do
            GNAME='u_xy'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
         END IF
         IF (RANK.EQ.0) 
     &     write(*,*) 'u_xy'

         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=0.5*(U2(i,NzMovie,j)+U2(i,NzMovie,j+1))
            end do
            end do
            GNAME='v_xy'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
          END IF
          IF (RANK.EQ.0) 
     &     write(*,*) 'v_xy'

         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=U3(i,NzMovie,j)
            end do
            end do
            GNAME='w_xy'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
         END IF
         IF (RANK.EQ.0) 
     &     write(*,*) 'w_xy'

         if (LES) then
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=NU_T(i,NzMovie,j)
            end do
            end do
            GNAME='nu_t_xy'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
         END IF
         end if

         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
               varxz(i,j)=U1(i,j,NyMovie)
            end do
            end do
            GNAME='u_xz'
            call writeHDF5_xzplane(FNAME,GNAME,varxz)
         END IF
         IF (RANK.EQ.0) 
     &     write(*,*) 'u_xz'

         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
               varxz(i,j)=0.5*(U2(i,j,NyMovie)+U2(i,j,NyMovie+1))
            end do
            end do
            GNAME='v_xz'
            call writeHDF5_xzplane(FNAME,GNAME,varxz)
         END IF 
         IF (RANK.EQ.0) 
     &     write(*,*) 'v_xz'

         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
               varxz(i,j)=U3(i,j,NyMovie)
            end do
            end do
            GNAME='w_xz'
            call writeHDF5_xzplane(FNAME,GNAME,varxz)
          END IF
          IF (RANK.EQ.0) 
     &     write(*,*) 'w_xz'

         IF (LES) then
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
              varxz(i,j)=NU_T(i,j,NyMovie)
            end do
            end do
            GNAME='nu_t_xz'
            call writeHDF5_xzplane(FNAME,GNAME,varxz) 
         end if
         END IF
         
         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=U1(NxMovie,i,j)
         end do
         end do
         GNAME='u_zy'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)
         IF (RANK.EQ.0) 
     &     write(*,*) 'u_zy'

         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=0.5*(U2(NxMovie,i,j)+U2(NxMovie,i,j+1))
         end do
         end do
         GNAME='v_zy'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)
         IF (RANK.EQ.0) 
     &     write(*,*) 'v_zy'

         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=U3(NxMovie,i,j)
         end do
         end do
         GNAME='w_zy'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)
         IF (RANK.EQ.0) 
     &     write(*,*) 'w_zy'

         if (USE_MPI) then
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (LES) then
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=NU_T(NxMovie,i,j)
         end do
         end do
         GNAME='nu_t_zy'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)
         END IF ! END IF LES

        !  if (USE_MPI) then
        !  call mpi_barrier(MPI_COMM_WORLD,ierror)
        !  end if
        END IF ! END IF MOVIE

        IF (RANK.EQ.0) 
     &     write(*,*) 'done velocity part of movie.h5'

#endif

C Convert velocity back to Fourier space
      call fft_xz_to_fourier(U1,CU1,0,NY+1)
      call fft_xz_to_fourier(U2,CU2,0,NY+1)
      call fft_xz_to_fourier(U3,CU3,0,NY+1)

      ! save spectra
      cs1 = cu1;
      if (RANKZ == 0) then
        cs1(0, 0, :) = 0.;
      end if
      cs1 = cs1*conjg(cs1);
  
      do j = 0, Ny + 1
        do i = 0, Nxp - 1
          cuu1_yx(j, i) = cs1(i, 0, j)
        end do
      end do
  
      DiagX = 0.d0
      do i = 0, Nxp - 1
        do j = 2, Ny
          DiagX(i) = DiagX(i) + 0.5 * 
     &               (cuu1_yx(j, i) + cuu1_yx(j - 1, i)) * dy(j) ! Integrate cuu1 in y
        end do
      end do
      if (use_mpi) then
        call mpi_allreduce(mpi_in_place, DiagX, Nxp,
     &       mpi_double_precision, mpi_sum, mpi_comm_y, ierror)
      end if
      DiagX = DiagX / Ly
      if(rankY == 0) then
        fname = 'mean.h5'
        gname = 'FTx_uu'
        call WriteStatH5_X(fname, gname, DiagX) 
      end if

      !Now compute z spectra

      cs1 = cu1;
      if (RANKZ == 0) then
        cs1(0, 0, :) = 0.;
      end if
      cs1 = cs1*conjg(cs1);
  
      do j = 0, Ny + 1
        do k = 0, TNKZ
          cuu1_yz(j, k) = cs1(0, k, j)
        end do
      end do
  
      DiagX = 0.d0
      do k = 0, TNKZ
        do j = 2, Ny
          DiagZ(k) = DiagZ(i) + 0.5 * 
     &               (cuu1_yz(j, k) + cuu1_yz(j - 1, k)) * dy(j) ! Integrate cuu1 in y
        end do
      end do
      if (use_mpi) then
        call mpi_allreduce(mpi_in_place, DiagZ, TNKZ,
     &       mpi_double_precision, mpi_sum, mpi_comm_y, ierror)
      end if
      DiagZ = DiagZ / Ly
      if(rankY == 0) then
        fname = 'mean.h5'
        gname = 'FTz_uu'
        call WriteStatH5_Z(fname, gname, DiagZ) 
      end if
      
! END IF FINAL
      end if

      IF (RANK.EQ.0) 
     &     write(*,*) 'done save_stats chan' 

      if (USE_MPI) then
      call mpi_barrier(MPI_COMM_WORLD,ierror)
      end if

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine bpebudget_chan
      !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      ! Compute the Background Potential Energy (BPE) -- Tseng & Ferziger 2001
      ! th(:,:,:,1) already in Physical Space
      include 'header'

      character(len=35) fname
      character(len=20) gname
      integer*8 i, j, k, bin, m
      integer, parameter :: Np_ = 5000
      integer, parameter :: N_edge = 500
      !-------------------------chebyshev----------
      integer, parameter :: Nbin = 3*Np_ + 2*N_edge
      real*8 phi(1:Np_), thb(1:Np_), thbmin(1:Np_)
      real*8 thbmax(1:Np_), thbmid(1:Np_)
      real*8 dTH_, dTH_2
      !----------------------------------------------------
      real*8 bincount(1:Nbin-1)
      real*8 thorpe_scale, total_bins
      real*8 thmin, thmax, BPE
      real*8 PDF(1:Nbin-1)
      real*8 Y_r(1:Nbin)
      real*8 Y_pos(1:Nbin-1)
      real*8 TH_bins(1:Nbin)

      ! Compute bounds of theta

      thmin = minval(th(:, :, :, 1))
      call MPI_ALLREDUCE(mpi_in_place,thmin,1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN,MPI_COMM_WORLD,ierror)
      thmax = maxval(th(:, :, :, 1))
      call MPI_ALLREDUCE(mpi_in_place,thmax,1,MPI_DOUBLE_PRECISION,
     &        MPI_MAX,MPI_COMM_WORLD,ierror)

    ! -------------chebyshev (tanh profile)----------------
      PI=4.D0*ATAN(1.D0)
      phi = 0.d0 
      dTH_ = 1.4d0
      do i = 1, Np_
        phi(i) = PI*(i-1)/Np_ 
        thb(i) = -0.7d0 + dTH_/2*(1.d0 + cos(phi(i)))
      end do
      do i = 1, Np_
        phi(i) = phi(i)/2
      end do 
      do i=1, Np_
        thbmax(i) = 0.7d0 + ((1.d0-0.7d0)*cos(phi(i))) 
        thbmin(i) = -0.7d0-((-0.7d0+1.d0)*cos(phi(Np_- i+1))) 
        thbmid(i) = thb(i) 
      end do     

    !   IF (RANK.EQ.0) 
    !  &     write(*,*) thbmax
    !   IF (RANK.EQ.0) 
    !  &     write(*,*) thbmin
      
      TH_bins = 0.d0

      do i = 1, N_edge
        TH_bins(i) = thmax + (i-1)*(1.d0 - thmax)/N_edge
        TH_bins(Nbin-i+1) = thmin - (i-1)*(thmin + 1.d0)/N_edge
      end do
      do i = 1, Np_
        TH_bins(i+N_edge) = thbmax(i)
        TH_bins(Np_+i+N_edge) = thb(i)
        TH_bins(2*Np_ +i+N_edge) = thbmin(i)
      end do

      bincount = 0.d0
      PDF = 0.d0
      Y_pos = 0.d0
      do i=0,NXM        
        do j=1,NY-1
          do k=0,NZP-1
            do m=1,Nbin-1
              if ((TH(i,k,j,1).GT.TH_bins(m+1)) .and.
     &            (TH(i,k,j,1).LE.TH_bins(m))) THEN
                  bincount(m) = bincount(m) + 1
                  PDF(m) = PDF(m) + DYF(j)
                  Y_pos(m) = Y_pos(m) + GYF(j)
              end if
            end do
          end do
        end do 
      end do

      call mpi_allreduce(mpi_in_place, bincount, Nbin, 
     &     mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(mpi_in_place, Y_pos, Nbin, 
     &     mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      !--------------------------------------

      call mpi_allreduce(mpi_in_place, PDF, Nbin, mpi_double_precision,
     &                   mpi_sum, mpi_comm_world, ierror)


      ! Compute Y_r (at b-mid-points)
      Y_r(1) = 0.d0
      do i = 2, Nbin-1
        Y_r(i) = Y_r(i - 1) + PDF(i-1)/(Nx * Nz) !
      end do
      ! Y_r = Y_r * Ly
      Y_r(Nbin) = Ly

      ! Compute BPE
      BPE = 0.d0
      thorpe_scale=0.d0
      do i = 1, Nbin-1 ! Integrate
        if (bincount(i).GT.0.0) then
          total_bins = total_bins + 1
          Y_pos(i) = -Y_pos(i)/bincount(i) + LY/2
    !       if (Y_r(i).GT.-5.0 .and. Y_r(i).LT.5.0) then
    !         thorpe_scale = thorpe_scale 
    !  &                   +(0.5*(Y_r(i + 1) + Y_r(i)) -
    !  &                     Y_pos(i))**2
          ! end if
        end if
      end do
      do i = 1, Nbin-1 
        BPE = BPE - (0.5 * (TH_bins(i)+TH_bins(i+1)) * 0.5 * 
     +        (Y_r(i + 1) + Y_r(i))) * (Y_r(i + 1) - Y_r(i)) / Ly
      end do
      thorpe_scale = thorpe_scale/total_bins
      
      fname = 'mean.h5'
      gname = 'BPE'
      call WriteHDF5_real(fname, gname, BPE)

      fname = 'mean.h5'
      gname = 'y_b'
      call WriteHDF5_real_long(fname, gname, Y_r, Nbin)

      fname = 'mean.h5'
      gname = 'th_b'
      call WriteHDF5_real_long(fname, gname, TH_bins, Nbin)

      ! fname = 'mean.h5'
      ! gname = 'bincount'
      ! call WriteHDF5_real_long(fname, gname, bincount, Nbin-1)

      end

            !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine bpebudget_chan_2
      !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      ! Compute the Background Potential Energy (BPE) -- Tseng & Ferziger 2001
      ! th(:,:,:,1) already in Physical Space
      include 'header'

      character(len=35) fname
      character(len=20) gname
      integer i, j, k, bin
      integer, parameter :: Nbin = 10000000
      real*8 thmin, thmax, dTH, BPE
      real*8 thorpe_scale, total_bins
      real*8 PDF(Nbin)
      real*8 Y_pos(1:Nbin)
      real*8 Y_r(Nbin+1)
      real*8 bincount(1:Nbin)

      ! Compute bounds of theta
      thmin = minval(th(:, :, :, 1))
      call MPI_ALLREDUCE(mpi_in_place,thmin,1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN,MPI_COMM_WORLD,ierror)
      thmax = maxval(th(:, :, :, 1))
      call MPI_ALLREDUCE(mpi_in_place,thmax,1,MPI_DOUBLE_PRECISION,
     &        MPI_MAX,MPI_COMM_WORLD,ierror)
      IF (RANK.EQ.0) write(*,*) 'rank: ', rankY
      dTH = (thmax - thmin) / Nbin + 1.d-14

      ! Compile the PDF in b
      PDF = 0.d0
      Y_pos = 0.d0
      bincount = 0.d0
      do j = jstart_th(1), jend_th(1) !j = 1, NY-1 
        do k = 0, NZP - 1
          do i = 0, NXM
            bin = int((th(i, k, j, 1) - thmin) / dTH) + 1
            PDF(bin) = PDF(bin) + dyf(j)
            Y_pos(bin) = Y_pos(bin) + GYF(j)
            bincount(bin) = bincount(bin) + 1
          end do
        end do
      end do
    

      ! if (rankY==NPROCY-1) then
      !   write(*,*) 'gyf: ', GYF(NY)
      !   do k = 0, NZP - 1
      !     do i = 0, NXM
      !       bin = int((th(i, k, NY, 1) - thmin) / dTH) + 1
      !       PDF(bin) = PDF(bin) + dyf(NY)
      !       Y_pos(bin) = Y_pos(bin) + GYF(NY)
      !       bincount(bin) = bincount(bin) + 1
      !     end do
      !   end do
      ! end if

      call mpi_allreduce(mpi_in_place, PDF, Nbin, mpi_double_precision,
     &                   mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(mpi_in_place, Y_pos, Nbin, 
     &     mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(mpi_in_place, bincount, Nbin,
     &                   mpi_double_precision,
     &                   mpi_sum, mpi_comm_world, ierror)

      ! Enforce \int_B PDF dB = 1 exactly (small dyf/2 at BCs...)  vs. /(Ly * Nx * Nz)
      PDF = PDF / (sum(PDF) * dTH)

      ! Compute Y_r (at b-mid-points)
      Y_r(1) = 0.d0
      do i = 2, Nbin
        Y_r(i) = Y_r(i - 1) + PDF(i-1) * dTH
      end do
      Y_r = Y_r * Ly
      Y_r(Nbin + 1) = Ly

      ! Compute BPE
      BPE = 0.d0
      ! thorpe_scale = 0.d0
      ! total_bins = 0.d0
      do i = 1, Nbin  ! Integrate
    !     if (bincount(i).GT.0.0) then
    !       total_bins = total_bins + 1
    !       Y_pos(i) = Y_pos(i)/bincount(i) + LY/2
    !       thorpe_scale = thorpe_scale 
    !  &                   +(0.5*(Y_r(i + 1) + Y_r(i)) -
    !  &                     Y_pos(i))**2
    !     end if
        BPE = BPE - (((i - 0.5) * dTH + thmin) * 0.5 * 
     +        (Y_r(i + 1) + Y_r(i))) * (Y_r(i + 1) - Y_r(i)) / Ly
      end do

      thorpe_scale = 0.d0
    !   do j = 1, NY !jstart_th(1), jend_th(1)
    !     do k = 0, NZP - 1
    !       do i = 0, NXM
    !         bin = int((th(i, k, j, 1) - thmin) / dTH) + 1            
    !         thorpe_scale = thorpe_scale + (0.5*(Y_r(bin)+Y_r(bin+1))
    !  &                     -(GYF(j)+LY/2))**2
    !       end do
    !     end do
    !   end do

      call mpi_allreduce(mpi_in_place, thorpe_scale, 1, 
     &                   mpi_double_precision,
     &                   mpi_sum, mpi_comm_world, ierror)

      thorpe_scale = thorpe_scale/(NY*NZP*(NXM+1)*NPROCY*NPROCZ)

      fname = 'mean.h5'
      gname = 'BPE2'
      call WriteHDF5_real(fname, gname, BPE)

      ! fname = 'mean.h5'
      ! gname = 'thorpe_scale2'
      ! call WriteHDF5_real(fname, gname, thorpe_scale)

      ! fname = 'mean.h5'
      ! gname = 'y_pos'
      ! call WriteHDF5_real_long(fname, gname, Y_r, Nbin+1)

      ! fname = 'mean.h5'
      ! gname = 'bincount2'
      ! call WriteHDF5_real_long(fname, gname, bincount, Nbin)
      end

      subroutine tkebudget_chan_les
! Calculate the componet of th SGS dissipation rate 
! only includes the terms timestepped implicitly
      include 'header'
      include 'header_les'

      character*35 FNAME
      CHARACTER*20 GNAME
      real*8 epsilon_sgs(NY)
      real*8 Diag(NY)
      integer i,j,k

! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>

      DO J=1,NY
        DO K=0,NZP-1
          DO I=0,NXM
            TEMP(I,K,J)=U1(I,K,J)*
     &        (  (NU_T(I,K,J+1) * (U1(I,K,J+1) - U1(I,K,J)) / DY(J+1)
     &         -  NU_T(I,K,J) * (U1(I,K,J)   - U1(I,K,J-1)) / DY(J))
     &               /DYF(J)  )
     &           +U3(I,K,J)*
     &        (  (NU_T(I,K,J+1) * (U3(I,K,J+1) - U3(I,K,J)) / DY(J+1)
     &        - NU_T(I,K,J) * (U3(I,K,J)   - U3(I,K,J-1)) / DY(J))
     &              /DYF(J)  )
     &           +U2(I,K,J)*
     &     ((0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*(U2(I,K,J+1)-U2(I,K,J))
     &                                              / DYF(J)
     &    -0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*(U2(I,K,J)-U2(I,K,J-1))
     &                                          / DYF(J-1))   /DY(J)  )
          END DO
        END DO
      END DO
 
! Now calculate the horizontal average
        do j=1,NY
          epsilon_sgs(j)=0.d0
          do i=0,NXM
          do k=0,NZP-1
            epsilon_sgs(j)=epsilon_sgs(j)+TEMP(I,K,J)
          end do
          end do
        end do

      call mpi_allreduce(mpi_in_place,epsilon_sgs,NY+2
     &    ,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)        


#ifdef HDF5
      FNAME='tke.h5'

      IF (RANKZ.eq.0) THEN
        gname='epsilon_sgs'
        Diag=epsilon_sgs(1:NY)
        call WriteStatH5(FNAME,gname,Diag)
      END IF
#else
! Here are aren't using HDF5, so write to text files
      IF (RANKZ.EQ.0) THEN
      IF (USE_MPI) THEN
        FNAME='tke_les'//trim(MPI_IO_NUM)//'.txt'
      ELSE
        FNAME='tke_les.txt'
      END IF
      open(46,file=FNAME,form='formatted',status='unknown')

      write(46,*) TIME_STEP,TIME,DELTA_T
        do j=1,NY
          write(46,460) j,GYF(J),epsilon_sgs(J)
        end do
      END IF
460     format(I3,' ',2(F30.20,' '))
#endif

      END


      subroutine tkebudget_chan
! Calculate the turbulent dissipation rate, epsilon
! Note that this is actually the pseudo-dissipation (see Pope, Turb. Flows)
! for an explanation
      include 'header'

      character*35 FNAME
      character*20 GNAME
      real*8 Diag(1:NY)
      real*8 varxy(0:NXM,1:NY),varzy(0:NZP-1,1:NY),varxz(0:NXM,0:NZP-1)
      integer i,j,k

! Store the 3D dissipation rate in F1, F4
      F1(:,:,:)=0.d0
      F4(:,:,:)=0.d0

! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
! epsilon will be calculated on the GY grid
      do j=1,NY
        epsilon(j)=0.
        epsilon_prime(j)=0.
      end do
! Store du/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKX(i)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
         epsilon(j)=epsilon(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         epsilon_prime(j)=epsilon_prime(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))    
         F1(i,k,j)=F1(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         F4(i,k,j)=F4(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKX(i)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
        epsilon_prime(j)=epsilon_prime(j)+(S1(i,k,j)**2.0)
        F1(i,k,j)=F1(i,k,j)+(S1(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dy at GY gridpoints, note remove mean
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
         S1(i,k,j)=((U1(i,k,j)-ume(j))
     &          -(U1(i,k,j-1)-ume(j-1)))
     &             /DY(j)
         S2(i,k,j)=((U1(i,k,j))
     &          -(U1(i,k,j-1)))
     &             /DY(j)
      end do
      end do
      end do
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S2(i,k,j)**2.0)
        epsilon_prime(j)=epsilon_prime(j)+(S1(i,k,j)**2.0)
        F1(i,k,j)=F1(i,k,j)+(S2(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
         epsilon(j)=epsilon(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         epsilon_prime(j)=epsilon_prime(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         F1(i,k,j)=F1(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         F4(i,k,j)=F4(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
      end do
      end do
      end do
! Compute du/dz at GY gridpoints
! Store du/dz in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKZ(k)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
         epsilon(j)=epsilon(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         epsilon_prime(j)=epsilon_prime(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         F1(i,k,j)=F1(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         F4(i,k,j)=F4(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
      end do
      end do
      end do
! Compute dv/dy at GY gridpoints, note remove mean
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
       S1(i,k,j)=((U2(i,k,j+1)-vme(j+1))
     &        -(U2(i,k,j-1)-vme(j-1)))
     &            /(GY(j+1)-GY(j-1))
       S2(i,k,j)=((U2(i,k,j+1))
     &        -(U2(i,k,j-1)))
     &            /(GY(j+1)-GY(j-1))
      end do
      end do
      end do
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S2(i,k,j)**2.0)
        F1(i,k,j)=F1(i,k,j)+(S2(i,k,j)**2.0)
        epsilon_prime(j)=epsilon_prime(j)+(S1(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dw/dy at GY gridpoints, note remove mean
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
         S1(i,k,j)=((U3(i,k,j)-wme(j))
     &          -(U3(i,k,j-1)-wme(j-1)))
     &             /DY(j)
         S2(i,k,j)=((U3(i,k,j))
     &          -(U3(i,k,j-1)))
     &             /DY(j)
      end do
      end do
      end do
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S2(i,k,j)**2.0)
        F1(i,k,j)=F1(i,k,j)+(S2(i,k,j)**2.0)
        epsilon_prime(j)=epsilon_prime(j)+(S1(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dz in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
         CS1(i,k,j)=CIKZ(k)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
        epsilon_prime(j)=epsilon_prime(j)+(S1(i,k,j)**2.0)
        F1(i,k,j)=F1(i,k,j)+(S1(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKZ(k)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
         epsilon(j)=epsilon(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         epsilon_prime(j)=epsilon_prime(j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         F1(i,k,j)=F1(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
         F4(i,k,j)=F4(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
      end do
      end do
      end do
      do j=1,NY
        epsilon(j)=NU*epsilon(j)/dble(NX*NZ)
        epsilon_prime(j)=NU*epsilon_prime(j)/dble(NX*NZ)
      end do
      F1=NU*F1
      F4=NU*F4
      call mpi_allreduce(mpi_in_place,epsilon,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,epsilon_prime,NY+2,
     &     MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_Z,ierror)


#ifdef HDF5
      FNAME='tke.h5'

      ! gname='time'
      ! call WriteHDF5_real(FNAME,gname,TIME)

      IF (RANKZ.eq.0) THEN
        ! gname='gyf'
        ! Diag=gyf(1:NY)
        ! call WriteStatH5(FNAME,gname,Diag)

        gname='epsilon' 
        Diag=epsilon(1:NY)
        call WriteStatH5(FNAME,gname,Diag)

        gname='epsilon_prime'
        Diag=epsilon_prime(1:NY)
        call WriteStatH5(FNAME,gname,Diag)

      END IF

      if (MOVIE) then
         FNAME='movie.h5'
         if (USE_MPI) then
           call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKZ.EQ.RANKZMOVIE) THEN
            do I=0,NXM
            do J=1,NY
               varxy(i,j)=F1(i,NzMovie,j)
            end do
            end do
            GNAME='epsilon_xy'
            call writeHDF5_xyplane(FNAME,GNAME,varxy)
         END IF
         if (USE_MPI) then
           call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         IF (RANKY.EQ.RANKYMOVIE) THEN
            do I=0,NXM
            do J=0,NZP-1
               varxz(i,j)=F1(i,j,NyMovie)
            end do
            end do
            GNAME='epsilon_xz'
            call writeHDF5_xzplane(FNAME,GNAME,varxz)
         END IF
         if (USE_MPI) then
           call mpi_barrier(MPI_COMM_WORLD,ierror)
         end if
         do I=0,NZP-1
         do J=1,NY
            varzy(i,j)=F1(NxMovie,i,j)
         end do
         end do
         GNAME='epsilon_zy'
         call writeHDF5_zyplane(FNAME,GNAME,varzy)

       END IF

#else
! Here we aren't using HDF5, so write to a text file
      IF (RANKZ.EQ.0) THEN
! Write out the mean statistics at each time
      IF (USE_MPI) THEN
        FNAME='tke'//trim(MPI_IO_NUM)//'.txt'
      ELSE
        FNAME='tke.txt'
      END IF
      open(45,file=FNAME,form='formatted',status='unknown')
      write(45,*) TIME_STEP,TIME,DELTA_T
      do j=2,NYM
        write(45,401) j,GYF(J),epsilon(j)
      end do
401   format(I3,' ',2(F20.9,' '))
      end if
#endif

      return 
      end

      subroutine tkebudget_chan_3d

! Calculate the turbulent dissipation rate, epsilon
! Note that this is actually the pseudo-dissipation (see Pope, Turb. Flows)
! for an explanation
      include 'header'

      character*35 FNAME
      character*20 GNAME
      real*8 Diag(1:NY)
      real*8 varxy(0:NXM,1:NY),varzy(0:NZP-1,1:NY),varxz(0:NXM,0:NZP-1)
      integer i,j,k

! Store the 3D dissipation rate in F4
      CALL FFT_XZ_TO_FOURIER(F1,CR1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(F2,CR2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(F3,CR3,0,NY+1)
      F4(:,:,:)=0.d0

! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
! epsilon will be calculated on the GY grid
      do j=1,NY
        epsilon_3d(j)=0.
      end do
! Store du/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKX(i)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
         epsilon_3d(j)=epsilon_3d(j)+((DYF(j-1)*S1(i,k,j)
     &             +DYF(j)*S1(i,k,j-1))/(2.d0*DY(j)))**2
        ! epsilon_3d(j) = epsilon_3d(j) +S1(i,k,j-1)**2
         F4(i,k,j)=F4(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKX(i)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon_3d(j)=epsilon_3d(j)+(S1(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dy at GY gridpoints, note remove mean
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
          CS1(i,k,j)=CR1(i,k,j)
      end do
      end do
      end do
  ! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
         S2(i,k,j)=((S1(i,k,j))
     &          -(S1(i,k,j-1)))
     &             /DY(j)
      end do
      end do
      end do
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon_3d(j)=epsilon_3d(j)+(S2(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S2(i,k,j)**2.0)
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKX(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
         epsilon_3d(j)=epsilon_3d(j)+((DYF(j-1)*S1(i,k,j)
     &             +DYF(j)*S1(i,k,j-1))/(2.d0*DY(j)))**2
         F4(i,k,j)=F4(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
      end do
      end do
      end do
! ! Compute du/dz at GY gridpoints
! Store du/dz in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKZ(k)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
         epsilon_3d(j)=epsilon_3d(j)+((DYF(j-1)*S1(i,k,j)
     &             +DYF(j)*S1(i,k,j-1))/(2.d0*DY(j)))**2
         F4(i,k,j)=F4(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
      end do
      end do
      end do
! ! Compute dv/dy at GY gridpoints, note remove mean
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
          CS1(i,k,j)=CR2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY-1
      do k=0,NZP-1
      do i=0,NXM
       S2(i,k,j)=((S1(i,k,j+1))
     &        -(S1(i,k,j-1)))
     &            /(GY(j+1)-GY(j-1))
      end do
      end do
      end do
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon_3d(j)=epsilon_3d(j)+(S2(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S2(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dw/dy at GY gridpoints, note remove mean
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
          CS1(i,k,j)=CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
         S2(i,k,j)=((S1(i,k,j))
     &          -(S1(i,k,j-1)))
     &             /DY(j)
      end do
      end do
      end do
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon_3d(j)=epsilon_3d(j)+(S2(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S2(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dz in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
         CS1(i,k,j)=CIKZ(k)*0.5d0*(CR2(i,k,j+1)+CR2(i,k,j))
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=1,NY
      do k=0,NZP-1
      do i=0,NXM
        epsilon_3d(j)=epsilon_3d(j)+(S1(i,k,j)**2.0)
        F4(i,k,j)=F4(i,k,j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=1,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CS1(i,k,j)=CIKZ(k)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xz_to_physical(CS1,S1,0,NY+1)
      do j=2,NY
      do k=0,NZP-1
      do i=0,NXM
         epsilon_3d(j)=epsilon_3d(j)+((DYF(j-1)*S1(i,k,j)
     &             +DYF(j)*S1(i,k,j-1))/(2.d0*DY(j)))**2
         F4(i,k,j)=F4(i,k,j)+(DYF(j-1)*S1(i,k,j)**2.d0
     &             +DYF(j)*S1(i,k,j-1)**2.d0)/(2.d0*DY(j))
      end do
      end do
      end do
      do j=1,NY
        epsilon_3d(j)=NU*epsilon_3d(j)/dble(NX*NZ)
      end do
      F4=NU*F4
      
      call mpi_allreduce(mpi_in_place,epsilon_3d,NY+2,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)


#ifdef HDF5
      FNAME='tke.h5'

      IF (RANKZ.eq.0) THEN
        gname='epsilon_3d'
        Diag=epsilon_3d(1:NY)
        call WriteStatH5(FNAME,gname,Diag)
      END IF

      if (MOVIE) then
      END IF
#endif
      return
       end 

      SUBROUTINE calcbpe

      include 'header'
      REAL*8 LOWER, UPPER
      integer i,j,k,l
      INTEGER BIN
      INTEGER NBINS ! should be a multiple of 5
      REAL RESORT(0:750) !size should be NBINS + 1
      INTEGER*8 BINCOUNT(1:750), sum
      REAL CHEBY, DRHO, DRHO2
      REAL UPPER2, LOWER2

      LOWER=MINVAL(TH(:,:,:,1))
      UPPER=MAXVAL(TH(:,:,:,1))
      if (USE_MPI) then
        call mpi_barrier(MPI_COMM_WORLD,ierror)
      end if
      ! PRINT *, 'local ', LOWER, UPPER
      call mpi_allreduce(mpi_in_place,LOWER,1,MPI_DOUBLE_PRECISION,
     &     MPI_MIN,MPI_COMM_WORLD,ierror)
      call mpi_allreduce(mpi_in_place,UPPER,1,MPI_DOUBLE_PRECISION,
     &     MPI_MAX,MPI_COMM_WORLD,ierror)
      IF (RANK.EQ.0) write(*,*) 'GLOBAL LOWER: ',LOWER
      IF (RANK.EQ.0) write(*,*) 'GLOBAL UPPER: ',UPPER
      NBINS=750
      DRHO=(UPPER-LOWER)/float(NBINS)

      do i=0,INT(NBINS/5)
        CHEBY = PI/2 - i*PI/(2*NBINS/5)
        IF (i.eq.0) THEN 
          RESORT(i) = LOWER
          RESORT(NBINS-i) = UPPER
        END IF
        IF (i.ge.1) THEN
          RESORT(i) = RESORT(i-1) + cos(CHEBY)*DRHO
          RESORT(NBINS-i) = RESORT(NBINS-i+1) - cos(CHEBY)*DRHO
        END IF
      end do
      
      UPPER2=RESORT(INT(NBINS * 4/5))
      LOWER2=RESORT(INT(NBINS/5))
      DRHO2 = (UPPER2 - LOWER2)/float(NBINS * 3/5)
      IF (RANK.EQ.0) write(*,*) 'GLOBAL LOWER2: ',LOWER2
      IF (RANK.EQ.0) write(*,*) 'GLOBAL UPPER2: ',UPPER2

      do i=INT(NBINS/5 +1), INT(NBINS*4/5)
        RESORT(i) = RESORT(i-1) + DRHO2
      end do

      do i=1,750
        BINCOUNT(i)=0
      end do

      do i=0,NXM
      do j=1,NY-1
      do k=0,NZP-1
      do l=1,750
      if (TH(i,k,j,1)>=RESORT(l-1) .and. TH(i,k,j,1)<=RESORT(l)) THEN
        BINCOUNT(l) = BINCOUNT(l) + 1
      end if
      end do
      end do
      end do  
      end do
      call mpi_allreduce(mpi_in_place,BINCOUNT,750,MPI_INTEGER8,
     &     MPI_SUM,MPI_COMM_WORLD,ierror)
      sum=0
      do i=1,750
        sum = sum+BINCOUNT(i)
      end do
      IF (RANK.EQ.0) write(*,*) 'BINCOUNT sum: ', sum
      END SUBROUTINE calcbpe
        

