Module Kernels

! NEW METHOD
! Use Checkpoints from CSS to do seismology
! instead of worrying about each timestep, just do every hour or so
! also don't worry about parallelization. see how fast it is serial.



! START OLD NOTES

! This module handles generating results that are comparable to ring-diagram
!  seismology results. This involves convolving flow fields with sensitivity
!  kernels.

! Tiles are tracked for some amount of time (~25hr) and along longitude. This
!  means that the same tile will average over the velocity field not only over
!  an extended physical region but an extended time. These tiles need to be 
!  continuously 'tracked' within CSS so that the appropriate flow field can be
!  used.
! Tiles are created at some time during the run and updated for some duration.

! Kernels will need to be interpolated on to the simulation grid
! Since the r/z-spacing does not change, might as well re-interpolate the 
!  kernels on to the correct z-levels as a pre-computation. The x and y
!  positions change by arbitrary amounts during the run, so we can't 
!  pre-compute those.

! NOTE TO SELF: write this assuming no mpi, then figure out parallelization
!  once you know exactly how the domain is split among processors

! more notes:
! r_ind seems to be the radial position (in cm) of each radial point (r_ind(nr))
! theta(nth) and phi(nphi) give positions of those too in radians

! sim_time is the current simulation time in seconds

! mynr is the number of radial grid points stored locally
! myr1 is the first index of the local r-grid
! myr2 is the last index of the local r-grid
! I assume mynth,myth1,myth2, mynphi,myphi1,myphi2 are similar

Real*8 :: r_surface = 6.96D10

! Data structure for a single kernel
Type :: Kernel
	Real*8, Allocatable :: sens(:,:,:) ! 3d sensitivity
	Real*8, Allocatable :: sens_interp(:,:,:) ! interpolated sensitivity
	Integer :: k, n ! integer identifiers
	Integer :: kersize ! size in degrees
	Character(256) :: filename ! filename of the kernel in question
End Type Kernel

! Data structure for a single tile
Type :: Tile
	Real*8 :: clat, clon ! central position
	Real*8 :: trackrate ! in some units..
	Real*8 :: duration, ctime ! tracking duration and central time in seconds
	Integer :: tilesize ! size in degrees, should match kers(kerind)%kersize
	Real*8, Allocatable :: vx(:), vy(:) ! one for each kernel, each timestep
	Real*8, Allocatable :: weightx(:), weighty(:) ! for keeping track of weighting
End Type Tile

! Kernel description, applies to all kernels
Integer :: ker_dim(3) ! initial kernel dimensions
Integer :: numkers, numkersizes
Integer, Allocatable :: numkers_ts(:), kersize_deg(:)
Real*8, Allocatable :: ker_depth(:) ! initial z-positions of kernels
Type(Kernel), Allocatable :: kers(:) ! main array of kernels. includes all sizes

! Tile set
Integer :: numtiles
Type(Tile), Allocatable :: tiles(:)

Contains

	Subroutine Step_Tiles
		Implicit None
		Integer :: ii, ix, iy, iz
		Real*8 :: currphi, currtht, tsrad

		! main loop over every tile
		Do ii=1,numtiles

		Enddo


		! for each tile, add to the running total of vx,vy
		!  for each kernel. also increment the weighting

		! since the current processor doesn't keep the whole grid in local
		! memory, try to only loop over tiles that overlap the local grid.
		! weighting should reflect this.
	End Subroutine Step_Tiles

	! wrap a value between 0 and 2pi
	REAL*8 FUNCTION wrapphi (x)
		IMPLICIT NONE
		REAL*8 :: x
		wrapphi = x
		IF (x .GT. 6.28318530718D0) wrapphi = x - 6.28318530718D0
		IF (x .LT. 0D0) wrapphi = x + 6.28318530718D0
	END FUNCTION wrapphi
	
	! Load a specialized kernel set
	Subroutine Load_Kernels (nr, r_ind)
		Implicit None
		Integer :: ii, ij, c
		Character(256) :: fname_kerset, fname_ker
		Real*8 :: d1, d2, r1, r2
		Integer :: nr
		Real*8 :: r_ind(nr)
		Real*8, Allocatable :: slice(:,:)

		r1 = r_ind(1)
		r2 = r_ind(nr)
		! Start by reading a kernel reference file. Format:
		! ker_dim(1) ker_dim(2) ker_dim(3)
		! numkers_total
		! numtilesizes
		! numkers_ts(1)
		! numkers_ts(2)
		! ...
		!  [ start ts(1) ]
		! nind kind fname
		! nind kind fname
		! ...
		!  [ start ts(2) ]
		! ...

		Write(fname_kerset, '(A,A,A)') "idl", "/", "kernelset_128.dat"

		Open(55,file=fname_kerset) ! TODO: un-hardcode this filename
		! Read dimensions
		Read(55,'(3I8)') ker_dim(1), ker_dim(2), ker_dim(3)
		! Read depths
		Allocate(ker_depth(ker_dim(3)))
		Do ii=1,ker_dim(3)
			Read(55,'(F14.6)') ker_depth(ii)
		Enddo
		! Read number of kernels and sizes
		Read(55,'(I8)') numkers
		Read(55,'(I8)') numkersizes
		! Read info on each size
		Allocate(numkers_ts(numkersizes))
		Allocate(kersize_deg(numkersizes))
		Do ii=1,numkersizes
			Read(55,'(I8)') kersize_deg(ii)
			Read(55,'(I8)') numkers_ts(ii)
		Enddo
		! Begin reading info on each kernel
		Allocate(kers(numkers))
		c = 1
		Do ij=1,numkersizes
			Do ii=1,numkers_ts(ij)
				Read(55,*) kers(c)%n, kers(c)%k, kers(c)%filename
				kers(c)%kersize = kersize_deg(ij)
				c = c + 1
			Enddo
		Enddo
		Close(55)

		! Once the main list is read in, start loading the kernels individually
		Do ii=1,numkers
			! create space for the kernel
			Allocate(kers(ii)%sens(ker_dim(1), ker_dim(2), ker_dim(3)))
			fname_ker = "idl/" // TRIM(kers(ii)%filename)
			Open(56,file=fname_ker,form='unformatted',access='stream')
			Read(56) kers(ii)%sens ! fairly certain this works
			Close(56)
		Enddo


		! Then, interpolate the kernels on to the simulation z-grid
		
!		print*, r_ind
!		print*, nr, r1, r2

		Allocate(slice(ker_dim(1), ker_dim(2)))
		DO ii=1,numkers
!			DO ij=1,ker_dim(3)
!				PRINT*, ii, 1, ker_depth(ij), SUM(kers(ii)%sens(:,:,ij))
!			ENDDO
			! make space for interpolated kernel
			ALLOCATE(kers(ii)%sens_interp(ker_dim(1), ker_dim(2), nr))
			DO ij=1,nr
				! r_ind starts at bottom and gets bigger
				! ker%sens starts at surface

				! determine where to integrate between
				! this could be precomputed....
				IF (ij .EQ. 1) THEN ! bottom of sim
					d2 = MAX(ker_depth(ker_dim(3)), (r_surface-r1)/1D8)
				ELSE
					d2 = (r_surface - 0.5D0*(r_ind(ij) + r_ind(ij-1)))/1D8
				ENDIF

				IF (ij .EQ. nr) THEN ! top of sim
					d1 = MIN(ker_depth(1), (r_surface-r2)/1D8)
				ELSE
					d1 = (r_surface - 0.5D0*(r_ind(ij) + r_ind(ij+1)))/1D8
				ENDIF

!				PRINT*, ii, ij, d1, d2
				Call Integrate_Kernel(kers(ii), d1, d2, slice)

				kers(ii)%sens_interp(:,:,ij) = slice
				
!				PRINT*, ii, 2, (r_surface-r_ind(ij))/1D8, SUM(slice)
			ENDDO
			! get rid of non-interpolated sensitivity
			DEALLOCATE(kers(ii)%sens)
		ENDDO
		Deallocate(slice)
	End Subroutine Load_Kernels

	! create a slice of ker that is the integral from depth1 to depth2
	Subroutine Integrate_Kernel(ker, depth1, depth2, value)
		Implicit None
		Type(Kernel) :: ker
		Real*8 :: depth1, depth2
		Real*8 :: value(ker_dim(1),ker_dim(2))
		Real*8 :: leftval(ker_dim(1),ker_dim(2)), rightval(ker_dim(1),ker_dim(2))
		Real*8 :: leftpos, rightpos
		Integer :: ind1, ind2, ii, ij
		Logical :: done

		value(:,:) = 0D0

		! depth1 and depth2 can be outside the range of the kernel

		IF (depth2 .LT. ker_depth(1)) THEN ! entire range is too shallow
			value(:,:) = 0D0
		ELSE IF (depth1 .GT. ker_depth(ker_dim(3))) THEN ! entire range is too deep
			value(:,:) = 0D0
		ELSE
			! here, we are guaranteed that either depth1 or 2 are within the kernel
			done = .False.
			leftpos = depth1
			rightpos = MIN(ker_depth(next_higher(depth1)), depth2)
			DO WHILE (.NOT. done)
				CALL ker_interp(ker, leftpos, leftval)
				CALL ker_interp(ker, rightpos, rightval)
				! area of trapezoid
				value(:,:) = value(:,:) + 0.5D0*(rightval(:,:)+leftval(:,:))*(rightpos-leftpos)
!				PRINT*, 'int', leftpos, SUM(leftval), SUM(value)
				leftpos = rightpos
				IF (leftpos .EQ. depth2 .OR. leftpos .EQ. ker_depth(ker_dim(3))) THEN
					done = .TRUE.
				ELSE
					rightpos = MIN(depth2, ker_depth(next_higher(rightpos)))
				ENDIF
			ENDDO
		ENDIF

	End Subroutine Integrate_Kernel

	! return a slice of ker interpolated at depth x
	! NOTE: x may be exactly on a grid point, handle this case
	SUBROUTINE ker_interp (ker, x, ret)
		IMPLICIT NONE
		TYPE(kernel) :: ker
		REAL*8 :: x, x1, x2
		REAL*8 :: y1(ker_dim(1),ker_dim(2)), y2(ker_dim(1),ker_dim(2))
		REAL*8 :: ret(ker_dim(1),ker_dim(2))
		INTEGER :: ind1, ind2, ii

		IF (x .LT. ker_depth(1) .OR. x .GT. ker_depth(ker_dim(3))) THEN
			ret(:,:) = 0D0
		ELSE
			! not super elegant..
			ind1 = 1
			DO ii=1,ker_dim(3)
				IF (x .GT. ker_depth(ind1)) THEN
					ind1 = ind1 + 1
				ENDIF
			ENDDO
			ind2 = ker_dim(3)
			DO ii=ker_dim(3),1,-1
				IF (x .LT. ker_depth(ind2)) THEN
					ind2 = ind2 - 1
				ENDIF
			ENDDO
			IF (ind1 .LT. ind2) THEN
				x1 = ker_depth(ind1)
				x2 = ker_depth(ind2)
				y1(:,:) = ker%sens(:,:,ind1)
				y2(:,:) = ker%sens(:,:,ind2)
				ret(:,:) = y1(:,:) + (x-x1)*(y2(:,:)-y1(:,:))/(x2-x1)
			ELSE
				ret(:,:) = ker%sens(:,:,ind1)
			ENDIF
		ENDIF
!		PRINT*, x, ind1, ind2, x1, x2, SUM(ret)
	END SUBROUTINE ker_interp

	INTEGER FUNCTION next_lower (x)
		IMPLICIT NONE
		REAL*8 :: x
		INTEGER :: ii
		ii = ker_dim(3)
		DO WHILE (ker_depth(ii) .GE. x)
			ii = ii - 1
		ENDDO
		next_lower = ii
	END FUNCTION next_lower

	INTEGER FUNCTION next_higher (x)
		IMPLICIT NONE
		REAL*8 :: x
		INTEGER :: ii
		ii = 1
		DO WHILE (ker_depth(ii) .LE. x)
			ii = ii + 1
		ENDDO
		next_higher = ii
	END FUNCTION next_higher

	! turn depth in Mm into radius in cm
	REAL FUNCTION Translate_Depth (x)
		IMPLICIT NONE
		REAL :: x
		Translate_depth = 6.955D10 - (x * 1D8)
	END FUNCTION Translate_Depth


	! Take each tile and move it in longitude by an amount determined by the
	!  tracking rate and the timestep dt
	Subroutine Nudge_Tiles(dt)
		Implicit None
		Real*8, Intent(IN) :: dt
		Integer :: ii

		Do ii=1,numtiles
			tiles(ii)%clon = tiles(ii)%clon + dt * tiles(ii)%trackrate
		Enddo

	End Subroutine Nudge_Tiles



	Subroutine Unload_Kernels
		! lol, memory management
	End Subroutine Unload_Kernels

	! Create a set of tiles
	! Parameters needed: central lat/lon, range of lat/lon, densepack factor,
	! tile sizes to use, 
	Subroutine Create_Tiles(numts, tilesizes, clon, clat, lonrn, latrn, densepack)
		Implicit None
		Real*8, Intent(IN) :: clon, clat, lonrn, latrn, densepack
		Integer, Intent(IN) :: numts
		Integer, Intent(IN) :: tilesizes(numts)
		Real*8 :: space, apode, latmin, lonmin
		Integer :: ii, ij, ix, iy, kertsind, cnt, nx, ny, ntiles(numts)

		apode = 0.9375D0

		! count the number of tiles needed for each tilesize
		Do ii=1,numts
			space = tilesizes(ii) * densepack * apode
			nx = INT(lonrn/space+1)
			ny = INT(latrn/space+1)
			ntiles(ii) = nx*ny
		Enddo

		! allocate space for all of the tiles
		numtiles = SUM(ntiles)
		Allocate(tiles(numtiles))

		! go through and load data into each tile
		cnt = 1
		Do ii=1,numts
			space = tilesizes(ii) * densepack * apode
			nx = INT(lonrn/space+1)
			ny = INT(latrn/space+1)
			lonmin = clon - 0.5*lonrn
			latmin = clat - 0.5*latrn
			kertsind = 1
			Do While (kersize_deg(kertsind) .NE. tilesizes(ii))
				kertsind = kertsind + 1
			Enddo
			If (kertsind .GT. numkersizes) Then
				PRINT*, "ALERT: could not find a matching kernel size. Expect issues."
			Endif

			Do ix=1,nx
				Do iy=1,ny
					! set tilesize
					tiles(cnt)%tilesize = tilesizes(ii)
					! set central position
					tiles(cnt)%clon = lonmin + (ix-1)*space
					tiles(cnt)%clat = latmin + (iy-1)*space
					! wrap longitude
					If (tiles(cnt)%clon .GT. 360D0) &
						tiles(cnt)%clon = tiles(cnt)%clon - 360D0
					If (tiles(cnt)%clon .LT. 0D0) &
						tiles(cnt)%clon = tiles(cnt)%clon + 360D0
					! set other attributes
					tiles(cnt)%trackrate = 1D-4 ! in deg per sec
					tiles(cnt)%duration = 2048D0*45D0 ! in seconds
					tiles(cnt)%ctime = 1024D0*45D0 ! in seconds
					! there are numkers_ts(kertsind) kernels for this tile
					Allocate(tiles(cnt)%vx(numkers_ts(kertsind)))
					Allocate(tiles(cnt)%vy(numkers_ts(kertsind)))
					Allocate(tiles(cnt)%weightx(numkers_ts(kertsind)))
					Allocate(tiles(cnt)%weighty(numkers_ts(kertsind)))
					! make damn sure they start at zero
					Do ij=1,numkers_ts(kertsind)
						tiles(cnt)%vx(ij) = 0D0
						tiles(cnt)%vy(ij) = 0D0
						tiles(cnt)%weightx(ij) = 0D0
						tiles(cnt)%weighty(ij) = 0D0
					Enddo
					cnt = cnt + 1
				Enddo
			Enddo
		Enddo

	End Subroutine Create_Tiles

	! Have processor 0 gather weights and output final results
	Subroutine Output_Tileset
		! vars?(mynth, mynphi, mynr, nv)
	End Subroutine Output_Tileset

	! cubic interpolation kernel
	! useful range of x is -2 to +2
	REAL FUNCTION CUBICKERNEL (x)
		IMPLICIT NONE
		REAL :: x, res
		x = ABS(x)
		IF (x .le. 1.0) THEN
			res = 1.0 + x*x*(-2.5 + 1.5*x)
		ELSE IF (x .lt. 2.0) THEN
			res = 2.0 + x*(-4.0 + x*(2.5 - 0.5*x))
		ELSE
			res = 0.0
		ENDIF
		CUBICKERNEL = res
	END FUNCTION CUBICKERNEL

End Module Kernels
