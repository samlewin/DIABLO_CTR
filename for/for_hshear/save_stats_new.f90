!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine save_stats_chan(movie,final)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  include 'header'

  character(len=35) fname
  character(len=20) gname
  logical movie,final
  integer i, j, k, n
  real(rkind) uc, ubulk

  ! This variable is used to add up scalar diagnostics
  real(rkind) thsum(0:Nyp + 1)
  ! These variables are used to store and write 2D slices
  real(rkind) varxy(0:Nxm1, 1:Nyp), varzy(0:Nzp - 1, 1:Nyp), varxz(0:Nxm1, 0:Nzp - 1)

  ! These variable are used for HDF5 writing
  real(rkind) Diag(1:Nyp)
  real(rkind) DiagX(0:Nxp - 1)

  if (rank == 0) &
    write (*, '("Saving flow statistics for time step " I10)') time_step

  if (use_mpi) then
    call mpi_barrier(mpi_comm_world, ierror)
    call apply_BC_vel_mpi ! Apply BCs FIRST (it screws up ghost cells...)
    call ghost_chan_mpi
    call ghost_chan_mpi_j0 ! Need the j = 0 boundary filled for les output
  else
    call apply_BC_vel_lower
    call apply_BC_vel_upper
  end if


  if (final) then
    ! We are done with the simulation

    fname = 'stats.h5'
    if (use_mpi) then
      call mpi_barrier(mpi_comm_world, ierror)
    end if
    if (rankZ == 0) then
      Diag = gyf(1:Nyp)
      gname = 'gzf'
      call WriteStatH5_Y(fname, gname, Diag)
      Diag = ubar(1:Nyp)
      gname = 'ubar'
      call WriteStatH5_Y(fname, gname, Diag)
      Diag = vbar(1:Nyp)
      gname = 'wbar'
      call WriteStatH5_Y(fname, gname, Diag)
      Diag = wbar(1:Nyp)
      gname = 'vbar'
      call WriteStatH5_Y(fname, gname, Diag)
      do n = 1, N_th
        Diag = thbar(1:Nyp, n)
        write (gname,'("thbar", I0.2)') n
        call WriteStatH5_Y(fname, gname, Diag)
      end do
    end if

  else

    ! Compute and write out the centerline velocity
    if (NprocY == 1) then
      if (int(float(Nyp) / 2.) == float(Nyp) / 2.) then
        ! IF Nyp is even
        uc = dble(cu1(0, 0, int(float(Nyp) / 2.)))
      else
        uc = 0.5 * (dble(cu1(0, 0, int(float(Nyp) / 2.) - 1)) &
                    + dble(cu1(0, 0, int(float(Nyp) / 2.))))
      end if
      write (*, '("Centerline velocity = " ES26.18)') uc
      ! Compute and write out bulk velocity
    end if

    ! We are in the middle of a run, compile statistics
    ! First get the number of samples taken so far
    if (rank == 0) write (*, '("Time    = " ES15.8 "    dt = " ES15.8)') time, dt ! Note: dt is the physical / CFL-constrained time-step
    if (rankZ == 0) then
      Nsamples = Nsamples + 1
      ! Get the mean velocity
      do j = 1, Nyp
        ubar(j) = (1./float(Nsamples)) * dble(cu1(0, 0, j)) &
                  + ((float(Nsamples) - 1.) / float(Nsamples)) * ubar(j)
        vbar(j) = (1./float(Nsamples)) * dble(cu2(0, 0, j)) &
                  + ((float(Nsamples) - 1.) / float(Nsamples)) * vbar(j)
        wbar(j) = (1./float(Nsamples)) * dble(cu3(0, 0, j)) &
                  + ((float(Nsamples) - 1.) / float(Nsamples)) * wbar(j)
        do n = 1, N_th
          thbar(j, n) = (1./float(Nsamples)) * dble(cth(0, 0, j, n)) &
                        + ((float(Nsamples) - 1.) / float(Nsamples)) * thbar(j, n)
        end do
      end do

      ! Integrate the instantaneous mean profile numerically at GY points
      ume = cu1(0, 0, :)
    else
      ume = 0.d0
    end if
    call integrate_y_var(ume, ubulk, mpi_comm_world)
    ! Write out UBULK
    if (rank == 0) write (*,  '("U Bulk  = " ES26.18)') ubulk

    ! Save CUi
    do k = 0, twoNkz
      do i = 0, Nxp - 1 ! Nkx
        do j = 0, Nyp + 1
          cr1(i, k, j) = cu1(i, k, j)
          cr2(i, k, j) = cu2(i, k, j)
          cr3(i, k, j) = cu3(i, k, j)
        end do
      end do
    end do

    ! Get the mean value of the velocities
    if (rankZ == 0) then
      ume = dble(cu1(0, 0, :))
      vme = dble(cu2(0, 0, :))
      wme = dble(cu3(0, 0, :))
      do n = 1, N_th
        thme(:, n) = dble(cth(0, 0, :, n))
      end do
    end if
    call mpi_bcast(ume, Nyp + 2, mpi_double_precision, 0, &
                   mpi_comm_z, ierror)
    call mpi_bcast(vme, Nyp + 2, mpi_double_precision, 0, &
                   mpi_comm_z, ierror)
    call mpi_bcast(wme, Nyp + 2, mpi_double_precision, 0, &
                   mpi_comm_z, ierror)
    if (N_th > 0) call mpi_bcast(thme, (Nyp + 2) * N_th, &
                                    mpi_double_precision, 0, mpi_comm_z, ierror)

    ! Convert to physical space
    call fft_xz_to_physical(cu1, u1, 0, Nyp + 1)
    call fft_xz_to_physical(cu2, u2, 0, Nyp + 1)
    call fft_xz_to_physical(cu3, u3, 0, Nyp + 1)

    ! Calculate the dissipation rate
    call tkebudget_chan(movie)
    if (LES) then
      call ghost_les_mpi ! Share nu_t
      call tkebudget_chan_les
    end if

    ! Get the turbulent kinetic energy at each level
    do j = 1, Nyp
      urms(j) = 0.
      vrms(j) = 0.
      wrms(j) = 0.
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          urms(j) = urms(j) + (u1(i, k, j) - ume(j))**2.
          vrms(j) = vrms(j) + 0.5 * ((u2(i, k, j) - vme(j))**2.+ &
                                     (u2(i, k, j + 1) - vme(j + 1))**2.)
          wrms(j) = wrms(j) + (u3(i, k, j) - wme(j))**2.
        end do
      end do
    end do

    call mpi_allreduce(mpi_in_place, urms, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, vrms, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, wrms, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)

    urms = sqrt(urms / float(Nx * Nz))
    vrms = sqrt(vrms / float(Nx * Nz))
    wrms = sqrt(wrms / float(Nx * Nz))

    ! Get the bulk rms value
    call integrate_y_var(urms, urms_b, mpi_comm_y)
    call integrate_y_var(vrms, vrms_b, mpi_comm_y)
    call integrate_y_var(wrms, wrms_b, mpi_comm_y)

    ! Compute the Reynolds stress and mean velocity gradient
    ! Here, uv and wv are defined on the GY grid
    ! uw is defined on the GYF grid
    do j = 1, Nyp
      uv(j) = 0.
      uw(j) = 0.
      wv(j) = 0.
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          uv(j) = uv(j) + &
                  (dyf(j - 1) * (u1(i, k, j) - ume(j)) + dyf(j) * (u1(i, k, j - 1) - ume(j - 1))) &
                  / (2.d0 * dy(j)) &
                  * (u2(i, k, j) - vme(j))
          wv(j) = wv(j) + &
                  (dyf(j - 1) * (u3(i, k, j) - wme(j)) + dyf(j) * (u3(i, k, j - 1) - wme(j - 1))) &
                  / (2.d0 * dy(j)) &
                  * (u2(i, k, j) - vme(j))
          uw(j) = uw(j) + (u1(i, k, j) - ume(j)) &
                  * (u3(i, k, j) - wme(j))
        end do
      end do
    end do

    call mpi_allreduce(mpi_in_place, uv, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, uw, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, wv, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)

    uv = uv / float(Nx * Nz)
    uw = uw / float(Nx * Nz)
    wv = wv / float(Nx * Nz)

    ! Get the y-derivative of the mean velocity at GY points
    do j = 1, Nyp
      dudy(j) = (ume(j) - ume(j - 1)) / (gyf(j) - gyf(j - 1))
      dwdy(j) = (wme(j) - wme(j - 1)) / (gyf(j) - gyf(j - 1))
    end do

    ! Calculate the mean square shear
    do j = 1, Nyp
      shear(j) = 0.d0
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          shear(j) = shear(j) &
                     + ((u1(i, k, j + 1) - u1(i, k, j - 1)) / (2.d0 * dyf(j)))**2.d0 &
                     + ((u3(i, k, j + 1) - u3(i, k, j - 1)) / (2.d0 * dyf(j)))**2.d0
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, shear, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)

    shear = shear / float(Nx * Nz)

    ! Write out the bulk rms velocity
    if (rank == 0) then
      write (*,  '("<U_rms> = " ES26.18)') urms_b
      write (*,  '("<V_rms> = " ES26.18)') wrms_b
      write (*,  '("<W_rms> = " ES26.18)') vrms_b
    end if

    ! Get the rms vorticity
    ! First, get the x-component in fourier space
    do j = 1, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1 !Nkx
          cs1(i, k, j) = (cr3(i, k, j + 1) - cr3(i, k, j - 1)) / (2.d0 * dyf(j)) &
                         - cikz(k) * 0.5d0 * (cr2(i, k, j + 1) + cr2(i, k, j))
        end do
      end do
    end do
    ! Convert to physical space
    call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
    ! Get the rms value
    do j = 1, Nyp
      omega_x(j) = 0.d0
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          omega_x(j) = omega_x(j) + s1(i, k, j)**2.d0
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, omega_x, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)

    omega_x = sqrt(omega_x / float(Nx * Nz))

    ! Now, get the y-component in fourier space
    do j = 1, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1 !Nkx
          cs1(i, k, j) = cikz(k) * cr1(i, k, j) - cikx(i) * cr3(i, k, j)
        end do
      end do
    end do
    ! Convert to physical space
    call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
    ! Get the rms value
    do j = 1, Nyp
      omega_y(j) = 0.d0
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          omega_y(j) = omega_y(j) + s1(i, k, j)**2.d0
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, omega_y, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)

    omega_y = sqrt(omega_y / float(Nx * Nz))

    ! Now, get the y-component in fourier space
    do j = 1, Nyp
      do k = 0, twoNkz
        do i = 0, Nxp - 1 ! Nkx
          cs1(i, k, j) = cikx(i) * 0.5d0 * (cr2(i, k, j + 1) + cr2(i, k, j)) &
                         - (cr1(i, k, j + 1) - cr1(i, k, j - 1)) / (2.d0 * dyf(j))
        end do
      end do
    end do
    ! Convert to physical space
    call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
    ! Get the rms value
    do j = 1, Nyp
      omega_z(j) = 0.d0
      do k = 0, Nzp - 1
        do i = 0, Nxm1
          omega_z(j) = omega_z(j) + s1(i, k, j)**2.d0
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, omega_z, Nyp + 2, mpi_double_precision, &
                       mpi_sum, mpi_comm_z, ierror)

    omega_z = sqrt(omega_z / float(Nx * Nz))

    fname = 'mean.h5'

    gname = 'time'
    call WriteHDF5_real(fname, gname, time)

    if (rankZ == 0) then

      gname = 'gzf'
      Diag = gyf(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'ume'
      Diag = ume(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'wme'
      Diag = vme(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'vme'
      Diag = wme(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'urms'
      Diag = urms(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'wrms'
      Diag = vrms(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'vrms'
      Diag = wrms(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'uw'
      Diag = uv(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'uv'
      Diag = uw(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'wv'
      Diag = wv(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'dudz'
      Diag = dudy(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'dvdz'
      Diag = dwdy(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'cp'
      Diag = dble(cp(0, 0, 1:Nyp))
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'shear'
      Diag = shear(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'omega_x'
      Diag = omega_x(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'omega_z'
      Diag = omega_y(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'omega_y'
      Diag = omega_z(1:Nyp)
      call WriteStatH5_Y(fname, gname, Diag)

    end if

    ! Do over the number of passive scalars
    do n = 1, N_th

      ! Save CTH
      do k = 0, twoNkz
        do i = 0, Nxp - 1 ! Nkx
          do j = 0, Nyp + 1
            crth(i, k, j, n) = cth(i, k, j, n)
          end do
        end do
      end do

      ! Compute the scalar gradient and store in CRi
      do j = 1, Nyp
        do k = 0, twoNkz
          do i = 0, Nxp - 1 ! Nkx
            ! Store gradients of TH(:,:,:,n) (if it is used) in CRi
            cr1(i, k, j) = cikx(i) * cth(i, k, j, n)
            cr2(i, k, j) = (cth(i, k, j + 1, n) - cth(i, k, j - 1, n)) / (gyf(j + 1) - gyf(j - 1))
            cr3(i, k, j) = cikz(k) * cth(i, k, j, n)
          end do
        end do
      end do
      ! Convert gradients to physical space
      call fft_xz_to_physical(cr1, r1, 0, Nyp + 1)
      call fft_xz_to_physical(cr2, r2, 0, Nyp + 1)
      call fft_xz_to_physical(cr3, r3, 0, Nyp + 1)

      ! Convert to physical space

      call mpi_barrier(mpi_comm_world, ierror)

      cs1(:, :, :) = cth(:, :, :, n)
      call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
      th(:, :, :, n) = s1(:, :, :)

      do j = 1, Nyp
        thsum(j) = 0.
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            thsum(j) = thsum(j) + (abs(th(i, k, j, n) - thme(j, n)))**2.
          end do
        end do
      end do
      call mpi_allreduce(mpi_in_place, thsum, (Nyp + 2), &
                         mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

      thrms(:, n) = sqrt(thsum / float(Nx * Nz))

      ! Compute the Reynolds stress and mean velocity gradient
      do j = 1, Nyp
        thsum(j) = 0.
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            thsum(j) = thsum(j) + (th(i, k, j, n) - thme(j, n)) &
                       * (0.5 * (u2(i, k, j) + u2(i, k, j + 1)) &
                          - 0.5 * (vme(j) + vme(j + 1)))
          end do
        end do
      end do
      call mpi_allreduce(mpi_in_place, thsum, (Nyp + 2), &
                         mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

      thv(:, n) = thsum / float(Nx * Nz)

      ! Get the y-derivative of the mean scalar at GY points
      do j = 1, Nyp
        dthdy(j, n) = (thme(j, n) - thme(j - 1, n)) / (gyf(j) - gyf(j - 1))
      end do

      ! Compute the potential energy dissipation, grad(TH) \cdot grad(TH)
      do j = 1, Nyp
        thsum(j) = 0.d0
        do k = 0, Nzp - 1
          do i = 0, Nxm1
            thsum(j) = thsum(j) &
                       + r1(i, k, j)**2.d0 + r2(i, k, j)**2.d0 + r3(i, k, j)**2.d0
          end do
        end do
      end do
      call mpi_allreduce(mpi_in_place, thsum, (Nyp + 2), &
                         mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

      pe_diss(:, n) = thsum / float(Nx * Nz)


      call bpebudget_chan ! Compute PE budget terms (BPE)


      if (movie) then
        fname = 'movie.h5'
        if (n == 1) then
          if (use_mpi) then
            call mpi_barrier(mpi_comm_world, ierror)
          end if
          if (rankZ == rankzmovie) then
            do j = 1, Nyp
              do i = 0, Nxm1
                varxy(i, j) = th(i, NzMovie, j, n)
              end do
            end do
            gname = 'th1_xz'
            call writeHDF5_xyplane(fname, gname, varxy)
          end if
          if (use_mpi) then
            call mpi_barrier(mpi_comm_world, ierror)
          end if
          if (rankY == rankymovie) then
            do j = 0, Nzp - 1
              do i = 0, Nxm1
                varxz(i, j) = th(i, j, NyMovie, n)
              end do
            end do
            gname = 'th1_xy'
            call writeHDF5_xzplane(fname, gname, varxz)
          end if
          if (use_mpi) then
            call mpi_barrier(mpi_comm_world, ierror)
          end if
          do j = 1, Nyp
            do i = 0, Nzp - 1
              varzy(i, j) = th(NxMovie, i, j, n)
            end do
          end do
          gname = 'th1_yz'
          call writeHDF5_zyplane(fname, gname, varzy)
        else if (n == 2) then
          if (use_mpi) then
            call mpi_barrier(mpi_comm_world, ierror)
          end if
          if (rankZ == rankzmovie) then
            do j = 1, Nyp
              do i = 0, Nxm1
                varxy(i, j) = th(i, NzMovie, j, n)
              end do
            end do
            gname = 'th2_xz'
            call writeHDF5_xyplane(fname, gname, varxy)
          end if
          if (use_mpi) then
            call mpi_barrier(mpi_comm_world, ierror)
          end if
          if (rankY == rankymovie) then
            do j = 0, Nzp - 1
              do i = 0, Nxm1
                varxz(i, j) = th(i, j, NyMovie, n)
              end do
            end do
            gname = 'th2_xy'
            call writeHDF5_xzplane(fname, gname, varxz)
          end if
          if (use_mpi) then
            call mpi_barrier(mpi_comm_world, ierror)
          end if
          do j = 1, Nyp
            do i = 0, Nzp - 1
              varzy(i, j) = th(NxMovie, i, j, n)
            end do
          end do
          gname = 'th2_yz'
          call writeHDF5_zyplane(fname, gname, varzy)
        end if

      end if

      ! Convert back to Fourier space
      s1(:, :, :) = th(:, :, :, n)
      call fft_xz_to_fourier(s1, cs1, 0, Nyp + 1)
      cth(:, :, :, n) = cs1(:, :, :)

    end do ! Over number of passive scalars, n

    fname = 'mean.h5'

    if (rankZ == 0) then

      do n = 1, N_th

        Diag = thme(1:Nyp, n)
        write (gname,'("thme", I0.2)') n
        call WriteStatH5_Y(fname, gname, Diag)

        Diag = dthdy(1:Nyp, n)
        write (gname,'("dthdz", I0.2)') n
        call WriteStatH5_Y(fname, gname, Diag)

        Diag = thrms(1:Nyp, n)
        write (gname,'("thrms", I0.2)') n
        call WriteStatH5_Y(fname, gname, Diag)

        Diag = thv(1:Nyp, n) ! \bar{w'b'}
        write (gname,'("thw", I0.2)') n
        call WriteStatH5_Y(fname, gname, Diag)

        Diag = pe_diss(1:Nyp, n)
        write (gname,'("pe_diss", I0.2)') n
        call WriteStatH5_Y(fname, gname, Diag)

      end do

    end if

    if (movie) then
      fname = 'movie.h5'
      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      if (rankZ == rankzmovie) then
        do j = 1, Nyp
          do i = 0, Nxm1
            varxy(i, j) = u1(i, NzMovie, j)
          end do
        end do
        gname = 'u_xz'
        call writeHDF5_xyplane(fname, gname, varxy)
      end if

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      if (rankZ == rankzmovie) then
        do j = 1, Nyp
          do i = 0, Nxm1
            varxy(i, j) = 0.5 * (u2(i, NzMovie, j) + u2(i, NzMovie, j + 1))
          end do
        end do
        gname = 'w_xz'
        call writeHDF5_xyplane(fname, gname, varxy)
      end if

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      if (rankZ == rankzmovie) then
        do j = 1, Nyp
          do i = 0, Nxm1
            varxy(i, j) = u3(i, NzMovie, j)
          end do
        end do
        gname = 'v_xz'
        call writeHDF5_xyplane(fname, gname, varxy)
      end if

      if (LES) then
        if (use_mpi) then
          call mpi_barrier(mpi_comm_world, ierror)
        end if
        if (rankZ == rankzmovie) then
          do j = 1, Nyp
            do i = 0, Nxm1
              varxy(i, j) = nu_t(i, NzMovie, j)
            end do
          end do
          gname = 'nu_t_xz'
          call writeHDF5_xyplane(fname, gname, varxy)
        end if
      end if

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      if (rankY == rankymovie) then
        do j = 0, Nzp - 1
          do i = 0, Nxm1
            varxz(i, j) = u1(i, j, NyMovie)
          end do
        end do
        gname = 'u_xy'
        call writeHDF5_xzplane(fname, gname, varxz)
      end if

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      if (rankY == rankymovie) then
        do j = 0, Nzp - 1
          do i = 0, Nxm1
            varxz(i, j) = 0.5 * (u2(i, j, NyMovie) + u2(i, j, NyMovie + 1))
          end do
        end do
        gname = 'w_xy'
        call writeHDF5_xzplane(fname, gname, varxz)
      end if

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      if (rankY == rankymovie) then
        do j = 0, Nzp - 1
          do i = 0, Nxm1
            varxz(i, j) = u3(i, j, NyMovie)
          end do
        end do
        gname = 'v_xy'
        call writeHDF5_xzplane(fname, gname, varxz)
      end if

      if (LES) then
        if (use_mpi) then
          call mpi_barrier(mpi_comm_world, ierror)
        end if
        if (rankY == rankymovie) then
          do j = 0, Nzp - 1
            do i = 0, Nxm1
              varxz(i, j) = nu_t(i, j, NyMovie)
            end do
          end do
          gname = 'nu_t_xy'
          call writeHDF5_xzplane(fname, gname, varxz)
        end if
      end if

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      do j = 1, Nyp
        do i = 0, Nzp - 1
          varzy(i, j) = u1(NxMovie, i, j)
        end do
      end do
      gname = 'u_yz'
      call writeHDF5_zyplane(fname, gname, varzy)

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      do j = 1, Nyp
        do i = 0, Nzp - 1
          varzy(i, j) = 0.5 * (u2(NxMovie, i, j) + u2(NxMovie, i, j + 1))
        end do
      end do
      gname = 'w_yz'
      call writeHDF5_zyplane(fname, gname, varzy)

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      do j = 1, Nyp
        do i = 0, Nzp - 1
          varzy(i, j) = u3(NxMovie, i, j)
        end do
      end do
      gname = 'v_yz'
      call writeHDF5_zyplane(fname, gname, varzy)

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
      if (LES) then
        do j = 1, Nyp
          do i = 0, Nzp - 1
            varzy(i, j) = nu_t(NxMovie, i, j)
          end do
        end do
        gname = 'nu_t_yz'
        call writeHDF5_zyplane(fname, gname, varzy)
      end if ! END IF LES

      if (use_mpi) then
        call mpi_barrier(mpi_comm_world, ierror)
      end if
    end if ! END IF MOVIE


    ! Convert velocity back to Fourier space
    call fft_xz_to_fourier(u1, cu1, 0, Nyp + 1)
    call fft_xz_to_fourier(u2, cu2, 0, Nyp + 1)
    call fft_xz_to_fourier(u3, cu3, 0, Nyp + 1)

    ! Save Spectra
    cs1 = cu1;
    if (rankZ == 0) then
      cs1(0, 0, :) = 0.;
    end if
    cs1 = cs1*conjg(cs1);

    do j = 0, Nyp + 1
      do i = 0, Nxp - 1 ! Nkx
        cuu1_yx(j, i) = cs1(i, 0, j)
      end do
    end do

    DiagX = 0.
    do i = 0, Nxp - 1
      do j = 2, Nyp
        DiagX(i) = DiagX(i) + 0.5 * (cuu1_yx(j, i) + cuu1_yx(j - 1, i)) * dy(j) ! Integrate cuu1 in y
      end do
    end do
    if (use_mpi) then
      call mpi_allreduce(mpi_in_place, DiagX, Nxp, &
                         mpi_double_precision, mpi_sum, mpi_comm_y, ierror)
    end if
    DiagX = DiagX / Ly

    if (rankY == 0) then
      fname = 'mean.h5'
      gname = 'FTx_uu'
      call WriteStatH5_X(fname, gname, DiagX)
    end if


    ! END IF FINAL
  end if

  if (use_mpi) then
    call mpi_barrier(mpi_comm_world, ierror)
  end if

  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine bpebudget_chan
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Compute the Background Potential Energy (BPE) -- Tseng & Ferziger 2001
  ! th(:,:,:,1) already in Physical Space
  include 'header'

  character(len=35) fname
  character(len=20) gname
  integer i, j, k, bin
  integer, parameter :: Nbin = 16*Ny
  real(rkind) thmin, thmax, dTH, BPE
  real(rkind) PDF(Nbin)
  real(rkind) Y_r(Nbin)


  ! Compute bounds of theta

  thmin = minval(th(:, :, :, 1))
  call get_minimum_mpi(thmin)
  thmax = maxval(th(:, :, :, 1))
  call get_maximum_mpi(thmax)

  dTH = (thmax - thmin) / Nbin + 1.d-14

  ! Compile the PDF in b
  PDF = 0.d0
  do j = jstart_th(1), jend_th(1)
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        bin = int((th(i, k, j, 1) - thmin) / dTH) + 1
        PDF(bin) = PDF(bin) + dyf(j)
      end do
    end do
  end do

  call mpi_allreduce(mpi_in_place, PDF, Nbin, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, ierror)

  ! Enforce \int_B PDF dB = 1 exactly (small dyf/2 at BCs...)  vs. /(Ly * Nx * Nz)
  PDF = PDF / (sum(PDF) * dTH)

  ! Compute Y_r (at b-mid-points)
  Y_r(1) = 0.d0
  do i = 2, Nbin
    Y_r(i) = Y_r(i - 1) + PDF(i - 1) * dTH
  end do
  Y_r = Y_r * Ly

  ! Compute BPE
  BPE = 0.d0
  do i = 1, Nbin - 1  ! Integrate
    BPE = BPE - (((i - 0.5) * dTH + thmin) * 0.5 * (Y_r(i + 1) + Y_r(i))) * &
                    (Y_r(i + 1) - Y_r(i)) / Ly
  end do

  fname = 'mean.h5'
  gname = 'BPE'
  call WriteHDF5_real(fname, gname, BPE)

end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine tkebudget_chan_les
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Calculate the componet of the SGS dissipation rate
  ! only includes the terms timestepped implicitly
  include 'header'
  include 'header_les'

  character(len=35) fname
  character(len=20) gname
  real(rkind) eps_sgs2(1:Nyp)
  real(rkind) Diag(Nyp)
  integer i, j, k

  ! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
  ! At *consistent* GY points!

  ! Compute the contribution at GYF first. Store in S1.
  do j = 1, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = u1(i, k, j) * &
                        ((nu_t(i, k, j + 1) * (u1(i, k, j + 1) - u1(i, k, j)) / dy(j + 1) &
                               - nu_t(i, k, j) * (u1(i, k, j) - u1(i, k, j - 1)) / dy(j)) &
                           / dyf(j)) &  ! At GYF
                     + u3(i, k, j) * &
                        ((nu_t(i, k, j + 1) * (u3(i, k, j + 1) - u3(i, k, j)) / dy(j + 1) &
                               - nu_t(i, k, j) * (u3(i, k, j) - u3(i, k, j - 1)) / dy(j)) &
                           / dyf(j))  ! At GYF
      end do
    end do
  end do

! Then, interpolate the u1 & u2 contribution onto GY
! so that it conserves the dissipation as in code
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        temp(i, k, j) = 0.5d0 * (s1(i, k, j) + s1(i, k, j - 1)) & ! u1 & u3 at GY
                      + u2(i, k, j) * &
                        ((0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j + 1)) * (u2(i, k, j + 1) - u2(i, k, j)) &
                            / dyf(j) &
                        - 0.5d0 * (nu_t(i, k, j) + nu_t(i, k, j - 1)) * (u2(i, k, j) - u2(i, k, j - 1)) &
                            / dyf(j - 1)) / dy(j)) ! At GY
      end do
    end do
  end do

  ! Now calculate the horizontal average
  do j = 1, Nyp
    eps_sgs2(j) = 0.d0
    do i = 0, Nxm1
      do k = 0, Nzp - 1
        eps_sgs2(j) = eps_sgs2(j) + temp(i, k, j)
      end do
    end do
  end do

  call mpi_allreduce(mpi_in_place, eps_sgs2, Nyp &
                     , mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  fname = 'mean.h5'

  if (rankZ == 0) then
    gname = 'eps_sgs2'
    Diag = eps_sgs2(1:Nyp) / float(Nx * Nz)
    call WriteStatH5_Y(fname, gname, Diag)
  end if

end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine save_stats_mean_xy_chan
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Calculates the z-line averaged quantities to save as xy planes.
  include 'header'

  character(len=35) fname
  character(len=20) gname
  logical final
  integer i, j, k, n

  real(rkind) var_xy(0:Nxm1, 1:Nyp + 1)
  real(rkind) uvar_xy(0:Nxm1, 1:Nyp), vvar_xy(0:Nxm1, 1:Nyp), wvar_xy(0:Nxm1, 1:Nyp)
  real(rkind) thvar_xy(0:Nxm1, 1:Nyp)
  real(rkind) ume_xy(0:Nxm1, 1:Nyp), vme_xy(0:Nxm1, 1:Nyp), wme_xy(0:Nxm1, 1:Nyp)
  real(rkind) thme_xy(0:Nxm1, 1:Nyp)

  if (rank == 0) write (*, '("Saving XZ mean flow statistics.")')

  fname = 'mean_xz.h5'

  ! Ghost cells & mean quantities in ume, etc already computed in save_stats_chan


  ! Calculate the potential vorticity (at GYF points) and store in CS1
  ! Use CFi for working arrays, so
  ! NOTE: ** This should only be done after full RK steps **
  ! Note adding background velocity and buoyancy gradients
  vvar_xy = 0.d0
  ! Calculate dudz in XY plane (dv/dz-dw/dy)*dth/dx
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf1(i, k, j) = 0.5d0 * cikz(k) * (cu2(i, k, j) + cu2(i, k, j + 1)) &
                     - 0.5d0 * (cu3(i, k, j + 1) - cu3(i, k, j - 1)) / dyf(j)
        cf2(i, k, j) = cikx(i) * cth(i, k, j, 1)
      end do
    end do
  end do

  call fft_xz_to_physical(cf1, f1, 0, Nyp + 1)
  call fft_xz_to_physical(cf2, f2, 0, Nyp + 1)
  do j = 1, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        vvar_xy(i, j) = vvar_xy(i, j) + (f1(i, k, j) - dTHdx(1) / (Ro_inv / delta) ) & ! Include Thermal Wind
                                      * (f2(i, k, j) + dTHdx(1)) ! Include Background Gradient
      end do
    end do
  end do

  ! (du/dy - dv/dx)*(dth/dz)
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf1(i, k, j) = 0.5d0 * (cu1(i, k, j + 1) - cu1(i, k, j - 1)) / dyf(j) &
                     - 0.5d0 * cikx(i) * (cu2(i, k, j + 1) + cu2(i, k, j))
        cf2(i, k, j) = cikz(k) * cth(i, k, j, 1)
      end do
    end do
  end do
  call fft_xz_to_physical(cf1, f1, 0, Nyp + 1)
  call fft_xz_to_physical(cf2, f2, 0, Nyp + 1)
  do j = 1, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        vvar_xy(i, j) = vvar_xy(i, j) + (f1(i, k, j) - dTHdz(1) / (Ro_inv / delta)) &
                                      * (f2(i, k, j) + dTHdz(1))
      end do
    end do
  end do

  ! (dw/dx - du/dz)*dth/dy
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cf1(i, k, j) = cikx(i) * cu3(i, k, j) &
                     - cikz(k) * cu1(i, k, j)
        cf2(i, k, j) = 0.5d0 * (cth(i, k, j + 1, 1) - cth(i, k, j - 1, 1)) / dyf(j)
      end do
    end do
  end do
  call fft_xz_to_physical(cf1, f1, 0, Nyp + 1)
  call fft_xz_to_physical(cf2, f2, 0, Nyp + 1)
  do j = 1, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        ! Include planetary vorticity
        vvar_xy(i, j) = vvar_xy(i, j) + (f1(i, k, j) + Ro_inv / delta) * f2(i, k, j)
      end do
    end do
  end do

  vvar_xy = vvar_xy / float(Nz)

  if (use_mpi) then
    ! Need to sum across the NprocZ
    if (rankZ == 0) then
      call mpi_reduce(mpi_in_place, vvar_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    else ! Don't use mpi_in_place for other processes, except for allreduce...
      call mpi_reduce(vvar_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    end if
    call mpi_barrier(mpi_comm_world, ierror)
  end if
  if (rankZ == 0) then
    gname = 'q_xz'
    call writeHDF5_xyplane(fname, gname, vvar_xy)
  end if


  ! Convert to physical space
  call fft_xz_to_physical(cu1, u1, 0, Nyp + 1)
  call fft_xz_to_physical(cu2, u2, 0, Nyp + 1)
  call fft_xz_to_physical(cu3, u3, 0, Nyp + 1)
  call fft_xz_to_physical(cp, p, 0, Nyp + 1)
  do n = 1, N_th
     cs1(:, :, :) = cth(:, :, :, N)
     call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
     th(:, :, :, N) = s1(:, :, :)
  end do


  ! Get the z-mean XY plane velocities and buoyancy
  do j = 1, Nyp
    do i = 0, Nxm1
      ume_xy(i, j) = 0.d0
      vme_xy(i, j) = 0.d0
      wme_xy(i, j) = 0.d0
      thme_xy(i, j) = 0.d0
      do k = 0, Nzp - 1
        ume_xy(i, j) = ume_xy(i, j) + u1(i, k, j)
        vme_xy(i, j) = vme_xy(i, j) + 0.5d0 * (u2(i, k, j) + u2(i, k, j + 1))
        wme_xy(i, j) = wme_xy(i, j) + u3(i, k, j)
        thme_xy(i, j) = thme_xy(i, j) + th(i, k, j, 1)
      end do
    end do
  end do

  var_xy = 0.d0
  do j = 1, Nyp + 1
    do i = 0, Nxm1
      do k = 0, Nzp - 1
        var_xy(i, j) = var_xy(i, j) + u2(i, k, j) ! Mean on GY
      end do
    end do
  end do

  ume_xy = ume_xy / float(Nz)
  vme_xy = vme_xy / float(Nz)
  var_xy = var_xy / float(Nz)
  wme_xy = wme_xy / float(Nz)
  thme_xy = thme_xy / float(Nz)
  if (use_mpi) then
    ! Need to sum across the NprocZ, and communicate back to everyone
    if (rankZ == 0) then
      call mpi_reduce(mpi_in_place, vme_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    else ! Don't use mpi_in_place for other processes, except for allreduce...
      call mpi_reduce(vme_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    end if
    call mpi_allreduce(mpi_in_place, ume_xy, Nx * Nyp, &
                       mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, var_xy, Nx * (Nyp + 1), &
                       mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, wme_xy, Nx * Nyp, &
                       mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(mpi_in_place, thme_xy, Nx * Nyp, &
                    mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
    call mpi_barrier(mpi_comm_world, ierror)
  end if
  if (rankZ == 0) then
    gname = 'ume_xz'
    call writeHDF5_xyplane(fname, gname, ume_xy)
    gname = 'wme_xz'
    call writeHDF5_xyplane(fname, gname, vme_xy)
    gname = 'vme_xz'
    call writeHDF5_xyplane(fname, gname, wme_xy)
    gname = 'thme_xz'
    call writeHDF5_xyplane(fname, gname, thme_xy)
  end if


  ! Get the TKE and Buoyancy variance in each XY plane
  do j = 1, Nyp
    do i = 0, Nxm1
      uvar_xy(i, j) = 0.d0
      vvar_xy(i, j) = 0.d0
      wvar_xy(i, j) = 0.d0
      thvar_xy(i, j) = 0.d0
      do k = 0, Nzp - 1
        uvar_xy(i, j) = uvar_xy(i, j) + (u1(i, k, j) - ume_xy(i, j))**2.d0
        vvar_xy(i, j) = vvar_xy(i, j) + &
                          0.5d0 * ((u2(i, k, j) - var_xy(i, j))**2.d0 &
                                 + (u2(i, k, j + 1) - var_xy(i, j + 1))**2.d0)
        wvar_xy(i, j) = wvar_xy(i, j) + (u3(i, k, j) - wme_xy(i, j))**2.d0
        thvar_xy(i, j) = thvar_xy(i, j) + (th(i, k, j, 1) - thme_xy(i, j))**2.d0
      end do
    end do
  end do

  if (use_mpi) then
    ! Need to sum across the NprocZ, and communicate back to everyone
    if (rankZ == 0) then
      call mpi_reduce(mpi_in_place, uvar_xy, Nx * Nyp, &
                         mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(mpi_in_place, vvar_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(mpi_in_place, wvar_xy, Nx * Nyp, &
                         mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(mpi_in_place, thvar_xy, Nx * Nyp, &
                         mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    else ! Don't use mpi_in_place for other processes, except for allreduce...
      call mpi_reduce(uvar_xy, 0, Nx * Nyp, &
                         mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(vvar_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(wvar_xy, 0, Nx * Nyp, &
                         mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(thvar_xy, 0, Nx * Nyp, &
                         mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    end if
    call mpi_barrier(mpi_comm_world, ierror)
  end if
  if (rankZ == 0) then
    uvar_xy = sqrt(uvar_xy / float(Nz))
    vvar_xy = sqrt(vvar_xy / float(Nz))
    wvar_xy = sqrt(wvar_xy / float(Nz))
    thvar_xy = sqrt(thvar_xy / float(Nz))

    gname = 'urms_xz'
    call writeHDF5_xyplane(fname, gname, uvar_xy)
    gname = 'wrms_xz'
    call writeHDF5_xyplane(fname, gname, vvar_xy)
    gname = 'vrms_xz'
    call writeHDF5_xyplane(fname, gname, wvar_xy)
    gname = 'thrms_xz'
    call writeHDF5_xyplane(fname, gname, thvar_xy)
  end if


  ! Compute the momentum fluxes
  do j = 1, Nyp
    do i = 0, Nxm1
      uvar_xy(i, j) = 0.d0 ! uv
      vvar_xy(i, j) = 0.d0 ! wv
      wvar_xy(i, j) = 0.d0 ! uw
      do k = 0, Nzp - 1
        uvar_xy(i, j) = uvar_xy(i, j) + (u1(i, k, j) - ume_xy(i, j)) &
                                  * 0.5d0 * ((u2(i, k, j) - var_xy(i, j)) &
                                           + (u2(i, k, j + 1) - var_xy(i, j + 1)))
        vvar_xy(i, j) = vvar_xy(i, j) + (u3(i, k, j) - wme_xy(i, j)) &
                                  * 0.5d0 * ((u2(i, k, j) - var_xy(i, j)) &
                                           + (u2(i, k, j + 1) - var_xy(i, j + 1)))
        wvar_xy(i, j) = wvar_xy(i, j) + (u1(i, k, j) - ume_xy(i, j)) &
                                      * (u3(i, k, j) - wme_xy(i, j))
      end do
    end do
  end do

  uvar_xy = uvar_xy / float(Nz)
  vvar_xy = vvar_xy / float(Nz)
  wvar_xy = wvar_xy / float(Nz)
  if (use_mpi) then
    ! Need to sum across the NprocZ
    if (rankZ == 0) then
      call mpi_reduce(mpi_in_place, uvar_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(mpi_in_place, vvar_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(mpi_in_place, wvar_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    else ! Don't use mpi_in_place for other processes, except for allreduce...
      call mpi_reduce(uvar_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(vvar_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(wvar_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    end if
    call mpi_barrier(mpi_comm_world, ierror)
  end if
  if (rankZ == 0) then
    gname = 'uw_xz'
    call writeHDF5_xyplane(fname, gname, uvar_xy)
    gname = 'wv_xz'
    call writeHDF5_xyplane(fname, gname, vvar_xy)
    gname = 'uv_xz'
    call writeHDF5_xyplane(fname, gname, wvar_xy)
  end if


  ! Compute the buoyancy fluxes
  do j = 1, Nyp
    do i = 0, Nxm1
      uvar_xy(i, j) = 0.d0 ! thu
      vvar_xy(i, j) = 0.d0 ! thv
      wvar_xy(i, j) = 0.d0 ! thw
      do k = 0, Nzp - 1
        uvar_xy(i, j) = uvar_xy(i, j) + (th(i, k, j, 1) - thme_xy(i, j)) &
                                      * (u1(i, k, j) - ume_xy(i, j))
        vvar_xy(i, j) = vvar_xy(i, j) + (th(i, k, j, 1) - thme_xy(i, j)) &
                                  * 0.5d0 * ((u2(i, k, j) - var_xy(i, j)) &
                                           + (u2(i, k, j + 1) - var_xy(i, j + 1)))
        wvar_xy(i, j) = wvar_xy(i, j) + (th(i, k, j, 1) - thme_xy(i, j)) &
                                      * (u3(i, k, j) - wme_xy(i, j))
      end do
    end do
  end do

  uvar_xy = uvar_xy / float(Nz)
  vvar_xy = vvar_xy / float(Nz)
  wvar_xy = wvar_xy / float(Nz)
  if (use_mpi) then
    ! Need to sum across the NprocZ
    if (rankZ == 0) then
      call mpi_reduce(mpi_in_place, uvar_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(mpi_in_place, vvar_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(mpi_in_place, wvar_xy, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    else ! Don't use mpi_in_place for other processes, except for allreduce...
      call mpi_reduce(uvar_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(vvar_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
      call mpi_reduce(wvar_xy, 0, Nx * Nyp, &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_z, ierror)
    end if
    call mpi_barrier(mpi_comm_world, ierror)
  end if
  if (rankZ == 0) then
    gname = 'thu_xz'
    call writeHDF5_xyplane(fname, gname, uvar_xy)
    gname = 'thw_xz'
    call writeHDF5_xyplane(fname, gname, vvar_xy)
    gname = 'thv_xz'
    call writeHDF5_xyplane(fname, gname, wvar_xy)
  end if


  if (verbosity > 2 .and. rank == 0) write (*, '("Done saving XZ mean flow statistics.")')

  ! Convert velocity back to Fourier space
  call fft_xz_to_fourier(u1, cu1, 0, Nyp + 1)
  call fft_xz_to_fourier(u2, cu2, 0, Nyp + 1)
  call fft_xz_to_fourier(u3, cu3, 0, Nyp + 1)
  call fft_xz_to_fourier(p, cp, 0, Nyp + 1)
  do n = 1, N_th
     s1(:, :, :) = th(:, :, :, N)
     call fft_xz_to_fourier(s1, cs1, 0, Nyp + 1)
     cth(:, :, :, N) = cs1(:, :, :)
  end do

  return
end




!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine tkebudget_chan(movie)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  ! Calculate the turbulent dissipation rate, epsilon
  ! Note that this is actually the pseudo-dissipation (see Pope, Turb. Flows)
  ! for an explanation
  include 'header'

  character(len=35) fname
  character(len=20) gname
  logical movie
  real(rkind) Diag(1:Nyp)
  real(rkind) varxy(0:Nxm1, 1:Nyp), varzy(0:Nzp - 1, 1:Nyp), varxz(0:Nxm1, 0:Nzp - 1)
  integer i, j, k

  ! Store the 3D dissipation rate in F1
  f1(:, :, :) = 0.d0

  ! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
  ! epsilon will be calculated on the GY grid (2:Nyp)
  do j = 1, Nyp
    epsilon(j) = 0.
  end do
  ! Store du/dx in CS1
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikx(i) * cr1(i, k, j)
      end do
    end do
  end do
  ! Convert to physical space
  call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                   + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
        f1(i, k, j) = f1(i, k, j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                     + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
      end do
    end do
  end do
  ! Store dv/dx in CS1
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikx(i) * cr2(i, k, j)
      end do
    end do
  end do
  ! Convert to physical space
  call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do
  ! Compute du/dy at GY gridpoints, note remove mean
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = ((u1(i, k, j) - ume(j)) &
                       - (u1(i, k, j - 1) - ume(j - 1))) &
                      / dy(j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do
  ! Store dw/dx in CS1
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikx(i) * cr3(i, k, j)
      end do
    end do
  end do
  ! Convert to physical space
  call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                   + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
        f1(i, k, j) = f1(i, k, j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                     + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
      end do
    end do
  end do
  ! Compute du/dz at GY gridpoints
  ! Store du/dz in CS1
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikz(k) * cr1(i, k, j)
      end do
    end do
  end do
  ! Convert to physical space
  call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                   + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
        f1(i, k, j) = f1(i, k, j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                     + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
      end do
    end do
  end do
  ! Compute dv/dy at GY gridpoints, note remove mean
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = ((u2(i, k, j + 1) - vme(j + 1)) &
                       - (u2(i, k, j - 1) - vme(j - 1))) &
                      / (gy(j + 1) - gy(j - 1))
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do
  ! Compute dw/dy at GY gridpoints, note remove mean
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        s1(i, k, j) = ((u3(i, k, j) - wme(j)) &
                       - (u3(i, k, j - 1) - wme(j - 1))) &
                      / dy(j)
      end do
    end do
  end do
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do
  ! Store dv/dz in CS1
  do j = 2, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikz(k) * cr2(i, k, j)
      end do
    end do
  end do
  ! Convert to physical space
  call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (s1(i, k, j)**2.0)
        f1(i, k, j) = f1(i, k, j) + (s1(i, k, j)**2.0)
      end do
    end do
  end do
  ! Store dw/dz in CS1
  do j = 1, Nyp
    do k = 0, twoNkz
      do i = 0, Nxp - 1
        cs1(i, k, j) = cikz(k) * cr3(i, k, j)
      end do
    end do
  end do
  ! Convert to physical space
  call fft_xz_to_physical(cs1, s1, 0, Nyp + 1)
  do j = 2, Nyp
    do k = 0, Nzp - 1
      do i = 0, Nxm1
        epsilon(j) = epsilon(j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                   + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
        f1(i, k, j) = f1(i, k, j) + (dyf(j - 1) * s1(i, k, j)**2.d0 &
                                     + dyf(j) * s1(i, k, j - 1)**2.d0) / (2.d0 * dy(j))
      end do
    end do
  end do
  epsilon = nu * epsilon / float(Nx * Nz)
  f1 = nu * f1
  call mpi_allreduce(mpi_in_place, epsilon, Nyp + 2, mpi_double_precision, &
                     mpi_sum, mpi_comm_z, ierror)

  fname = 'mean.h5'

  if (rankZ == 0) then
    gname = 'epsilon'
    Diag = epsilon(1:Nyp)
    call WriteStatH5_Y(fname, gname, Diag)
  end if

  if (movie) then
    fname = 'movie.h5'
    if (use_mpi) then
      call mpi_barrier(mpi_comm_world, ierror)
    end if
    if (rankZ == rankzmovie) then
      do i = 0, Nxm1
        do j = 1, Nyp
          varxy(i, j) = f1(i, NzMovie, j)
        end do
      end do
      gname = 'epsilon_xz'
      call writeHDF5_xyplane(fname, gname, varxy)
    end if
    if (use_mpi) then
      call mpi_barrier(mpi_comm_world, ierror)
    end if
    if (rankY == rankymovie) then
      do i = 0, Nxm1
        do j = 0, Nzp - 1
          varxz(i, j) = f1(i, j, NyMovie)
        end do
      end do
      gname = 'epsilon_xy'
      call writeHDF5_xzplane(fname, gname, varxz)
    end if
    if (use_mpi) then
      call mpi_barrier(mpi_comm_world, ierror)
    end if
    do i = 0, Nzp - 1
      do j = 1, Nyp
        varzy(i, j) = f1(NxMovie, i, j)
      end do
    end do
    gname = 'epsilon_yz'
    call writeHDF5_zyplane(fname, gname, varzy)

  end if

  return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine save_stats_LES_OOL(blank)
  !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  include 'header'

  character(len=35) fname
  character(len=20) gname
  logical blank
  real(rkind) :: Diag(1:Nyp)

  if (blank) then
    fname = 'mean.h5'

    if (rankZ == 0) then
      Diag = 0.d0
      gname = 'nu_sgs'
      call WriteStatH5_Y(fname, gname, Diag)

      gname = 'eps_sgs1'
      call WriteStatH5_Y(fname, gname, Diag)

    end if
  else
    ! Needed to write out LES Statistics without timestepping...
    ! DON'T run this except for when stopping the simulation!

    rk_step = 1
    flag_save_LES = .true.

    call les_chan

    call fft_xz_to_fourier(u1, cu1, 0, Nyp + 1)
    call fft_xz_to_fourier(u2, cu2, 0, Nyp + 1)
    call fft_xz_to_fourier(u3, cu3, 0, Nyp + 1)

  end if

  return

end
