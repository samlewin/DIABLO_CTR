       SUBROUTINE USER_RHS_CHAN_PHYSICAL
       include 'header'
! Here, you can add terms to the right hand side
! of the momentum and scalar equations.  
! The right hand side forcing arrays, CF2, CF2, CF3, CFTH
! are in Fourier space.  The velocity and scalars are available 
! in physical space.
! S1 is available as a working variable

       integer i,j,k,n

       real*8 alpha

! For example, to add a linear damping term (e.g. -alpha*U) to the RHS:
!       alpha=-0.1d0
!       DO J=JSTART,JEND
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*U1(I,K,J)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=JSTART,JEND
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO

! For U2 do this...
! Note that the only thing that changes are the bounds of the J index
!       DO J=2,NY
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*U2(I,K,J)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=2,NY
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO

! For scalars, do this...
! Loop over all scalars
!       DO N=1,N_TH
!       DO J=JSTART,JEND
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*TH(I,K,J,N)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=JSTART,JEND
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CFTH(I,K,J,N)=CFTH(I,K,J,N)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO
!      END DO 

      RETURN 
      END


       SUBROUTINE USER_RHS_CHAN_FOURIER
       include 'header'
! Here, you can add terms to the right hand side
! of the momentum and scalar equations.  
! The right hand side forcing arrays, CF2, CF2, CF3, CFTH
! are in Fourier space.  The velocity and scalars are available 
! in Fourier space.
! S1 is available as a working variable

       integer i,j,k,n

       real*8 alpha

! For example, to add a linear damping term (e.g. -alpha*U) to the RHS:
!       alpha=1.d0
!       DO J=JSTART,JEND
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CF1(I,K,J)=CF1(I,K,J)-alpha*CU1(I,K,J)
!           END DO
!         END DO
!       END DO

! ! For U2 do this...
! ! Note that the only thing that changes are the bounds of the J index
!       DO J=2,NY
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CF2(I,K,J)=CF2(I,K,J)-alpha*CU2(I,K,J)
!           END DO
!         END DO
!       END DO

! ! For scalars, do this...
!       DO J=JSTART,JEND
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CFTH(I,K,J,N)=CFTH(I,K,J,N)-alpha*CTH(I,K,J,N)
!           END DO
!         END DO
!       END DO
      DO N=1,N_TH
        CALL SPONGE_TH(N)
      END DO
      CALL SPONGE_VEL

      RETURN 
      END


      SUBROUTINE USER_RHS_PER_PHYSICAL
      include 'header'
C Optionally, add forcing terms to the right hand side
C of the momentum and scalar evolution equations
C Here, CF1, CF2, CF3, and CFTH are the forcing arrays in Fourier space
C U1, U2, U3, and TH are the velocity and scalar variables in Physical space

      integer i,j,k,n

C For forced homogeneous turbulence:
C Add some forcing to the system to keep the Batchelor scale fixed
      EK=0.d0
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
              EK=EK+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
          END DO
        END DO
      END DO
C Note, that each cell has the same volume, so we can just average over all points
      EK=EK/dble(NX*NY*NZ)
! Scale EK by an amount to compensate for dissipation from 2/3 de-aliasing:
      EK=0.8d0*EK
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U1(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U2(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF3(I,K,J)=CF3(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END

      SUBROUTINE USER_RHS_DUCT_PHYSICAL
      include 'header'
      RETURN
      END

      SUBROUTINE USER_RHS_CAVITY_PHYSICAL
      include 'header'
      RETURN
      END 

            !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine sponge_th(n)
      !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      ! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
      ! specified background state for the temperature field
      ! The intention is to allow an open boundary
      include 'header'
      integer i, j, k, n
      real*8 L_sponge, L_top
      real*8 sponge_amp

      ! The following variables will store the background state
      real*8 th_0(-1:NY + 1)

      ! real*8 ri_b(0:NY + 1)

      ! This variable will hold the forcing rate
      real*8 sponge_sigma(0:NY + 1)

      ! Set the amplitude of the sponge
      sponge_amp = 1.d0
      ! Set the top of the sponge layer in physical units
      L_sponge = 17.5d0
      ! Set the bottom of the computational domain in physical units
      L_top = LY/2
      do j = 0, NY + 1
        ! Quadratic damping at lower wall
        if (gyf(j) .GT. L_sponge) then
          sponge_sigma(j) = sponge_amp * ((L_sponge - gyf(j)) 
     &                       / (L_sponge - L_top))**2.d0
        else if (gyf(j) .LT. -L_sponge) then
          sponge_sigma(j) = sponge_amp * ((L_sponge + gyf(j)) 
     &                       / (L_sponge - L_top))**2.d0
        else
          sponge_sigma(j) = 0.d0
        end if
      end do

      ! Set the profile for relaxing the mean TH
      do j = 0, NY + 1
        th_0(j) = gyf(j)
      end do

      ! Add damping to R-K terms
      ! Damp the perturbations towards 0
      do k = 0, TNKZ
        do i = 0, NXP-1
          if ((rankZ.NE.0) .or. (i.NE.0) .or. (k.NE.0)) then
            do j = jstart_th(n), jend_th(n)
              cfth(i, k, j, n) = cfth(i, k, j, n) 
     &                       - sponge_sigma(j) * (cth(i, k, j, n) - 0.)
            end do
          end if
        end do
      end do
      ! Damp the mean gradient towards TH_0
      do j = jstart_th(n), jend_th(n)
        cfth(0, 0, j, n) = cfth(0, 0, j, n) - sponge_sigma(j) 
     &                       * (cth(0, 0, j, n) - th_0(j))
      end do
      return
      end
      !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine sponge_vel
      !----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      ! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
      ! specified background state
      ! The intention is to allow an open boundary
      include 'header'
      integer i, j, k

      real*8 L_sponge, L_top
      real*8 sponge_amp

      ! The following variables will store the background state
      real*8 u1_0(-1:NY + 1), u2_0(0:NY + 1), u3_0(-1:NY + 1)

      ! This variable will hold the forcing rate
      real*8 sponge_sigma(0:NY + 1)

      ! Set the amplitude of the sponge
      sponge_amp = 1.d0
      ! Set the top of the sponge layer in physical units
      L_sponge = 17.5d0
      ! Set the bottom of the computational domain in physical units
      L_top = LY/2
      do j = 0, NY + 1
        ! Quadratic damping at lower wall
        if (gyf(j) > L_sponge) then
          sponge_sigma(j) = sponge_amp * ((L_sponge - gyf(j)) 
     &                           / (L_sponge - L_top))**2.d0
        else if ((gyf(j) < -L_sponge)) then
          sponge_sigma(j) = sponge_amp * ((L_sponge + gyf(j)) 
     &                           / (L_sponge - L_top))**2.d0
        else
          sponge_sigma(j) = 0.d0
        end if
      end do

      ! Set the background state 
      ! Here, set the background to be geostrophic, with a linear temperature profile
      do j = 0, NY + 1
        u1_0(j) = tanh(gyf(j))
        u3_0(j) = 0.d0
      end do
      do j = 0, NY + 1
        u2_0(j) = 0.d0
      end do

      ! Add damping function to explicit R-K
      do k = 0, TNKZ
        do i = 0, NXP-1
          if ((RANKZ.NE.0).or.(i.ne.0).or.(k.ne.0)) then
            do j = jstart, jend
              cf1(i, k, j) = cf1(i, k, j) 
     &         - sponge_sigma(j) * (cu1(i, k, j) - 0.d0)
              cf3(i, k, j) = cf3(i, k, j) 
     &         - sponge_sigma(j) * (cu3(i, k, j) - 0.d0)
            end do
            do j = 2, NY
              cf2(i, k, j) = cf2(i, k, j) - 
     &                         0.5 * (sponge_sigma(j) 
     &         + sponge_sigma(j + 1)) * (cu2(i, k, j) - 0.d0)
            end do
          end if
        end do
      end do
      ! Damp mean flow
      if (rankz.eq.0) then
      do j = jstart, jend
        cf1(0, 0, j) = cf1(0, 0, j) - sponge_sigma(j)
     &    * (cu1(0, 0, j) - u1_0(j))
        cf3(0, 0, j) = cf3(0, 0, j) - sponge_sigma(j)
     &     * (cu3(0, 0, j) - u3_0(j))
      end do
      do j = 2, NY
        cf2(0, 0, j) = cf2(0, 0, j) - sponge_sigma(j)
     &      * (cu2(0, 0, j) - u2_0(j))
      end do
      end if 

      return
      end





