
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

PROGRAM Seismology
	USE Kernels

	IMPLICIT NONE

	CHARACTER(256) :: fname_base, workdir
	INTEGER :: nr, nth, nphi, nv, istep
	INTEGER :: checknum
	REAL*8 :: deltat, sim_time, r1, r2, th1, th2, phi1, phi2
	REAL*8, ALLOCATABLE :: vx(:,:,:), vy(:,:,:)
	REAL*8, ALLOCATABLE :: lons(:), lats(:)
	REAL*8, ALLOCATABLE :: r_ind(:)
	INTEGER*4 :: reclen
	INTEGER :: ii, ij, junk, ix, iy
	REAL*8 :: ux, uyi
	REAL*8 :: PI, TWOPI, RADTODEG, DEGTORAD
	REAL*8 :: dist

	! tile stuff
	INTEGER :: numts
	INTEGER, ALLOCATABLE :: ts(:)
	REAL*8 :: clon, clat, lonrn, latrn, densepack

	! set up pies using machine precision
	PI = 4D0*atan(1D0)
	TWOPI = 2D0*PI
	RADTODEG = 180D0/PI
	DEGTORAD = PI/180D0
	
	! format the filenames
	workdir = "." ! read this one in, eventually
	checknum = 420000 ! and this one
	
	WRITE(fname_base, '(2A,I0.7,A)') TRIM(workdir), "/Checkpoints/", checknum, '/'

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

	! allocate space for a single checkpoint
	ALLOCATE(vx(nth,nphi,nr))
	ALLOCATE(vy(nth,nphi,nr))

	! open checkpoint file
	PRINT*, 'Loading velocity files..'
	reclen = nth*nphi*nr*2 ! double precision

	! u = vr
	! v = vy
	! w = vx
	OPEN(40, FILE=TRIM(fname_base)//"checkpoint.w",FORM='unformatted',ACCESS='direct',STATUS='unknown',RECL=reclen)
	READ(40, REC=1) vx
	CLOSE(40)
	OPEN(40, FILE=TRIM(fname_base)//"checkpoint.v",FORM='unformatted',ACCESS='direct',STATUS='unknown',RECL=reclen)
	READ(40, REC=1) vy
	CLOSE(40)

	! time to pre-compute which grid points are relevant for each tile!
	! scan through each tile counting how many grid points are relevant
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
		PRINT*, tiles(ii)%clon, tiles(ii)%clat, tiles(ii)%numgridpts
		! allocate..
	ENDDO

	! deallocate space
	DEALLOCATE(vx)
	DEALLOCATE(vy)
	DEALLOCATE(lons)
	DEALLOCATE(lats)
	
	DEALLOCATE(r_ind)
	DEALLOCATE(ts)

END PROGRAM Seismology
