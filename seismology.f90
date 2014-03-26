
PROGRAM Seismology
	USE Kernels

	IMPLICIT NONE

	CHARACTER(256) :: fname_base, workdir
	INTEGER :: nr, nth, nphi, nv, istep
	INTEGER :: checknum
	REAL*8 :: deltat, sim_time, r1, r2, th1, th2, phi1, phi2
	REAL*8, ALLOCATABLE :: vx(:,:,:), vy(:,:,:)
	REAL*8, ALLOCATABLE :: r_ind(:)
	INTEGER*4 :: reclen
	INTEGER :: ii, ij, junk
	REAL*8 :: ux, uy

	! tile stuff
	INTEGER :: numts
	INTEGER, ALLOCATABLE :: ts(:)
	REAL*8 :: clon, clat, lonrn, latrn, densepack

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

	! load kernels
	CALL Load_Kernels (nr, r_ind)
	PRINT*, "Kernels Loaded: ", numkers

	! create a tile set
	! Init seismology stuff
	numts = 1 ! number of tile sizes
	Allocate(ts(numts))
	ts(1) = 16 ! tile size in heliographic degrees
	! ANINT rounds to nearest integer, returns real
	clon = ANINT(0.5D0*57.2957795131D0*(phi1 + phi2))
	clat = ANINT(90D0-0.5D0*57.2957795131D0*(th1 + th2))
	lonrn = ANINT(57.2957795131D0*(phi2-phi1))
	latrn = ANINT(57.2957795131D0*(th2-th1))
	densepack = 5D-1
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

	! do stuff

	! deallocate space
	DEALLOCATE(meanflowx)
	DEALLOCATE(meanflowy)
	DEALLOCATE(vx)
	DEALLOCATE(vy)
	
	DEALLOCATE(r_ind)
	DEALLOCATE(ts)

END PROGRAM Seismology
