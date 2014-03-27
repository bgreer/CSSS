
! each checkpoint has the same dimensions and the same bounds
! don't try to track the tiles, that requires advection and stuff
! since each tile stays still, might as well pre-compute which 
! simulation pixels each tile interacts with. that way we don't need
! to search through the sim domain each time to figure out which pixels
! are relevant.
! when 'integrating' amode of a tile for a single time step, just go 
! through the list of relevant horizontal points at every depth. at
! each depth, interpolate the kernel on to the sim grid and stat summing
! estimated speed up? instead of searching through all grid points
! (hundreds of thousands of points) and computing the great-circle distance
! for each (6 trig operations), know relevant points instantly. seems
! worthwhile.

! want to also compute some statistics. the timeaveraged flow field will
! have different statistics than the individual checkpoints, so 
! collect stats for both.

! TODO
! fix kernel integration
! compute flow stats
! parallelize
! modify CSS to print radial grid in checkpoint header
! fix massive memory leaks
! modularize code, separate io, stats, kernels, etc

PROGRAM Seismology
	USE Kernels

	IMPLICIT NONE

	CHARACTER(256) :: fname_base, workdir
	INTEGER :: nr, nth, nphi, nv, istep
	REAL*8 :: deltat, sim_time, r1, r2, th1, th2, phi1, phi2
	REAL*8, ALLOCATABLE :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
	REAL*4, ALLOCATABLE :: vx_avg(:,:,:), vy_avg(:,:,:), vz_avg(:,:,:)
	REAL*8, ALLOCATABLE :: lons(:), lats(:)
	REAL*8, ALLOCATABLE :: r_ind(:)
	INTEGER*4 :: reclen
	INTEGER :: ii, ij, junk, ix, iy, iz, ik
	REAL*8 :: d3r, sens
	REAL*8 :: PI, TWOPI, RADTODEG, DEGTORAD
	REAL*8 :: dist

	! tile stuff
	INTEGER :: numts, ind
	INTEGER, ALLOCATABLE :: ts(:)
	REAL*8 :: clon, clat, lonrn, latrn, densepack

	! list of checkpoints
	INTEGER :: checkstart, checkend, checkskip
	INTEGER :: numchecks, currcheck

	checkstart = 420000
	checkend = 420000
	checkskip = 1000

	! computed, not set:
	numchecks = (checkend-checkstart)/checkskip + 1
	PRINT*, 'Checkpoints to Load: ', numchecks

	! set up pies using machine precision
	PI = 4D0*atan(1D0)
	TWOPI = 2D0*PI
	RADTODEG = 180D0/PI
	DEGTORAD = PI/180D0
	
	! format the filenames
	workdir = "." ! read this one in, eventually
	
	WRITE(fname_base, '(2A,I0.7,A)') TRIM(workdir), "/Checkpoints/", checkstart, '/'

	! read header
	PRINT*, 'Opening header file: ', TRIM(fname_base)//"header"
	OPEN(30, FILE=TRIM(fname_base)//"header")
	! read resolution
	READ(30, FMT='(5i9)') nr, nth, nphi, nv, istep
	! read physical boundaries
	READ(30, FMT='(8e14.7)') deltat, sim_time, r1, r2, th1, th2, phi1, phi2
	CLOSE(30)

	PRINT*, 'Simulation dimensions (r,th,phi): ', nr, nth, nphi

	! read non-uniform grid
	ALLOCATE(r_ind(nr))
	OPEN(50, FILE="grid")
	DO ii=1,nr
		READ(50, *) junk, r_ind(ii)
	ENDDO
	CLOSE(50)

	! create uniform grid for lat/lon
	ALLOCATE(lons(nphi))
	ALLOCATE(lats(nth))
	DO ix=1,nphi
		lons(ix) = RADTODEG*(phi1 + (phi2-phi1)*FLOAT(ix-1)/FLOAT(nphi-1))
	ENDDO
	DO iy=1,nth
		lats(iy) = 90D0-RADTODEG*(th1 + (th2-th1)*FLOAT(iy-1)/FLOAT(nth-1))
	ENDDO

	! load kernels
	CALL Load_Kernels (nr, r_ind)
	PRINT*, "Kernels Loaded: ", numkers

	! create a tile set
	! Init seismology stuff
	numts = 1 ! number of tile sizes
	Allocate(ts(numts))
	ts(1) = 16 ! tile size in heliographic degrees
	! ANINT rounds to nearest integer, returns real
	clon = ANINT(0.5D0*RADTODEG*(phi1 + phi2))
	clat = ANINT(90D0-0.5D0*RADTODEG*(th1 + th2))
	lonrn = ANINT(RADTODEG*(phi2-phi1))
	latrn = ANINT(RADTODEG*(th2-th1))
	densepack = 0.5D0 ! 1-this is overlap fraction. silly..
	CALL Create_Tiles(numts, ts, clon, clat, lonrn, latrn, densepack)
	PRINT*, "Tiles Created: ", numtiles


	! time to pre-compute which grid points are relevant for each tile!
	! scan through each tile counting how many grid points are relevant
	PRINT*, 'Precomputing Tile Grid Points..'
	DO ii=1,numtiles
		tiles(ii)%numgridpts = 0
		! print grid points that fall in this tile
		DO ix=1,nphi
			DO iy=1,nth
				! compute great circle distance between this sim pixel
				! and the tile center
				dist  = sin(DEGTORAD*lats(iy))*sin(DEGTORAD*tiles(ii)%clat) + &
					cos(DEGTORAD*lats(iy))*cos(DEGTORAD*tiles(ii)%clat)* &
					cos(DEGTORAD*(lons(ix) - tiles(ii)%clon))
				dist = RADTODEG*acos(dist)
				
				IF (dist .LT. FLOAT(tiles(ii)%tilesize)*apode*0.5D0) THEN
					tiles(ii)%numgridpts = tiles(ii)%numgridpts + 1
				ENDIF
			ENDDO
		ENDDO
		! allocate space
		ALLOCATE(tiles(ii)%gridptx(tiles(ii)%numgridpts))
		ALLOCATE(tiles(ii)%gridpty(tiles(ii)%numgridpts))
		! go back through and fill
		ind = 1
		DO ix=1,nphi
			DO iy=1,nth
				! compute great circle distance between this sim pixel
				! and the tile center
				dist  = sin(DEGTORAD*lats(iy))*sin(DEGTORAD*tiles(ii)%clat) + &
					cos(DEGTORAD*lats(iy))*cos(DEGTORAD*tiles(ii)%clat)* &
					cos(DEGTORAD*(lons(ix) - tiles(ii)%clon))
				dist = RADTODEG*acos(dist)
				IF (dist .LT. FLOAT(tiles(ii)%tilesize)*apode*0.5D0) THEN
					tiles(ii)%gridptx(ind) = ix
					tiles(ii)%gridpty(ind) = iy
					ind = ind + 1
				ENDIF
			ENDDO
		ENDDO
	ENDDO

	! allocate space for averaged flow field
	ALLOCATE(vx_avg(nth,nphi,nr))
	ALLOCATE(vy_avg(nth,nphi,nr))
	ALLOCATE(vz_avg(nth,nphi,nr))
	vx_avg = 0D0
	vy_avg = 0D0
	vz_avg = 0D0

	! allocate space for a single checkpoint
	ALLOCATE(vx(nth,nphi,nr))
	ALLOCATE(vy(nth,nphi,nr))
	ALLOCATE(vz(nth,nphi,nr))

	! loop over each checkpoint
	DO currcheck=checkstart,checkend,checkskip
		! open checkpoint file
		PRINT*, 'Loading velocity files..'
		reclen = nth*nphi*nr*2 ! double precision
		WRITE(fname_base, '(2A,I0.7,A)') TRIM(workdir), "/Checkpoints/", currcheck, '/'
		! u = vr
		! v = vy
		! w = vx
		OPEN(40, FILE=TRIM(fname_base)//"checkpoint.w",FORM='unformatted',ACCESS='direct',STATUS='unknown',RECL=reclen)
		READ(40, REC=1) vx
		CLOSE(40)
		OPEN(40, FILE=TRIM(fname_base)//"checkpoint.v",FORM='unformatted',ACCESS='direct',STATUS='unknown',RECL=reclen)
		READ(40, REC=1) vy
		CLOSE(40)
		OPEN(40, FILE=TRIM(fname_base)//"checkpoint.u",FORM='unformatted',ACCESS='direct',STATUS='unknown',RECL=reclen)
		READ(40, REC=1) vz
		CLOSE(40)

		! add to running average of flow field
		vx_avg = vx_avg + vx / FLOAT(numchecks)
		vy_avg = vy_avg + vy / FLOAT(numchecks)
		vz_avg = vz_avg + vz / FLOAT(numchecks)

		PRINT*, 'Convolving Flow Field..'
		!$OMP PARALLEL DO PRIVATE(ii,iz,d3r,ij,ix,iy,ik,sens)
		DO ii=1,1!numtiles
			! go through each depth
			DO iz=1,nr-1
				d3r = 1.1377777778D0 * (r_ind(iz+1)-r_ind(iz))*1D-8 / FLOAT(tiles(ii)%numgridpts)
				! sum horizontally only within tile
				DO ij=1,tiles(ii)%numgridpts
					ix = tiles(ii)%gridptx(ij)
					iy = tiles(ii)%gridpty(ij)
					! compute for each kernel
					DO ik=1,numkers
						IF (iz .LE. kers(ik)%lastdepth) THEN
!							sens = SUM(kers(ik)%sens_interp(:,:,iz))
							sens = 1D0
							tiles(ii)%vx(ik) = tiles(ii)%vx(ik) + vx(iy,ix,iz)*sens*d3r
							tiles(ii)%vy(ik) = tiles(ii)%vy(ik) + vy(iy,ix,iz)*sens*d3r
							tiles(ii)%weightx(ik) = tiles(ii)%weightx(ik) + sens*d3r
							tiles(ii)%weighty(ik) = tiles(ii)%weighty(ik) + sens*d3r
						ENDIF
					ENDDO
				ENDDO
			ENDDO
			! normalize?
			DO ik=1,numkers
				tiles(ii)%vx(ik) = tiles(ii)%vx(ik) / tiles(ii)%weightx(ik)
				tiles(ii)%vy(ik) = tiles(ii)%vy(ik) / tiles(ii)%weighty(ik)
			ENDDO
			PRINT*, ii
		ENDDO
		!$OMP END PARALLEL DO

	ENDDO ! end loop over checkpoints

	OPEN(50, FILE="invdata")
	DO ii=1,numtiles
		DO ik=1,numkers
			WRITE(50,'(2E15.6,2I5,4E15.6)') tiles(ii)%clon, tiles(ii)%clat, kers(ik)%k, kers(ik)%n, tiles(ii)%vx(ik), 0D0, tiles(ii)%vy(ik), 0D0
		ENDDO
	ENDDO
	CLOSE(50)

	! write out averaged flow field
	PRINT*, 'Writing Averaged Flow Field..'
	reclen = nth*nphi*nr ! single precision
	OPEN(60, FILE="flowfield.x",FORM='unformatted',ACCESS='direct',STATUS='unknown',RECL=reclen)
	WRITE(60, REC=1) vx_avg
	CLOSE(60)
	OPEN(60, FILE="flowfield.y",FORM='unformatted',ACCESS='direct',STATUS='unknown',RECL=reclen)
	WRITE(60, REC=1) vy_avg
	CLOSE(60)
	OPEN(60, FILE="flowfield.z",FORM='unformatted',ACCESS='direct',STATUS='unknown',RECL=reclen)
	WRITE(60, REC=1) vz_avg
	CLOSE(60)

	! deallocate space
	CALL Unload_Kernels()
	CALL Unload_Tiles()
	DEALLOCATE(vx)
	DEALLOCATE(vy)
	DEALLOCATE(vz)
	DEALLOCATE(lons)
	DEALLOCATE(lats)
	DEALLOCATE(vx_avg)
	DEALLOCATE(vy_avg)
	DEALLOCATE(vz_avg)
	
	DEALLOCATE(r_ind)
	DEALLOCATE(ts)

END PROGRAM Seismology
