c--------------------------------------------------------------------
c  ALE3D STATE VARIABLE INITIALIZATION SUBROUTINE
c--------------------------------------------------------------------

      subroutine sdvinia(statev,nstatv,seed,orient,ngrains)
      implicit none
c     input
      integer nstatv, orient, ngrains
      double precision seed
c     output
      double precision statev(nstatv)
c     static params
      integer nstatvsingle, nstatvpoly, nmax, bigint
      parameter(nstatvsingle=25,nstatvpoly=31,nmax=1000)
      parameter(bigint=1790989824) !arbitrary
c     util
      integer i, itemp, iseed, i8_uniform
      double precision phi1, theta, phi2, rmat(3,3)
      double precision angs(nmax, 3)

c     Sanity checks for number of state variables declared
      if (ngrains.eq.1) then
        if(nstatv.ne.nstatvsingle) then
          print*, '=== WARNING IN SDVINIA === '
          print*,'nstatv is: ', nstatv, ' but should be ', 
     &      nstatvsingle, ' for a single crystal'
        end if
      else if (ngrains.gt.1) then
        if(mod(nstatv,nstatvpoly).gt.0) then
          print*, '=== WARNING IN SDVINIA === '
          print*,'nstatv is: ', nstatv, ' and is not a multiple of ', 
     &      nstatvpoly
        else if (ngrains*nstatvpoly.ne.nstatv) then
          print*, '=== WARNING IN SDVINIA === '
          print*, 'simulation has ', ngrains, ' grains and nstatv is',
     &      nstatv
          print*, 'nstatv should be set to: ', nstatvpoly*ngrains
        end if
      else
        print*, '=== WARNING IN SDVINIA === '
        print*, 'ngrains should be a positive integer, but it is', 
     &    ngrains
        print*, 'ngrains is controlled by initializing statev3'
      end if

C     Load angles for orientation from fortran in hardcode.f
      if (orient.le.3) then !load ang from fortran
        iseed = int(bigint*seed)
        if(orient.eq.1) then !az31 rolled data
          call azrolledtex(angs)
        else if (orient.eq.2) then !az31 ecae data
          call azecaetex(angs)
        else if (orient.eq.3) then !amx602 data
          call amxtex(angs)
        else
          print*, '=== INVALID TEXTURE FILE CHOICE, EXITING === '
          print*, '-- Change value of ORIENTATION in input file'
          stop  
        end if
      end if

c     single xtal      
      if (ngrains.eq.1) then
        if (orient.ge.4) then
          call single_xtal_orients(orient, phi1, theta, phi2)
        else
          itemp = i8_uniform(iseed, nmax)
          phi1 = angs(itemp, 1)
          theta = angs(itemp, 2)
          phi2 = angs(itemp, 3)     
        end if
        call create_rmat_bunge(phi1, theta, phi2, rmat)     
        statev(1) = rmat(1,1)
        statev(2) = rmat(1,2)
        statev(3) = rmat(1,3)
        statev(4) = rmat(2,1)
        statev(5) = rmat(2,2)
        statev(6) = rmat(2,3)
        statev(7) = rmat(3,1)
        statev(8) = rmat(3,2)
        statev(9) = rmat(3,3)
c     polyxtal        
      else
        do i=1,ngrains      
          itemp = i8_uniform(iseed, nmax)
          phi1 = angs(itemp, 1)
          theta = angs(itemp, 2)
          phi2 = angs(itemp, 3)
          call create_rmat_bunge(phi1, theta, phi2, rmat)
          statev((i-1)*nstatvpoly+1) =  rmat(1,1)
          statev((i-1)*nstatvpoly+2) =  rmat(1,2)
          statev((i-1)*nstatvpoly+3) =  rmat(1,3)
          statev((i-1)*nstatvpoly+4) =  rmat(2,1)
          statev((i-1)*nstatvpoly+5) =  rmat(2,2)
          statev((i-1)*nstatvpoly+6) =  rmat(2,3)
          statev((i-1)*nstatvpoly+7) =  rmat(3,1)
          statev((i-1)*nstatvpoly+8) =  rmat(3,2)
          statev((i-1)*nstatvpoly+9) =  rmat(3,3)
        end do
      end if

      return
      end


c--------------------------------------------------------------------
c  ABAQUS STATE VARIABLE INITIALIZATION SUBROUTINE
c--------------------------------------------------------------------

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)
C
      IMPLICIT NONE
      INTEGER NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT
      DOUBLE PRECISION STATEV, COORDS
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
C

      double precision phi1, theta, phi2 !Bunge Euler Angle (rad)
      double precision rmat(3,3)

      !TODO - CHANGE SO READS IN FROM ELSEWHERE

      if (NCRDS.le.3) then !load an angle file
        print*, '=== INVALID ORIENTATION FOR SINGLE XTAL, EXITING === '
        print*, '-- Note 1-3 are reserved for polyxtal deformation'
        print*, '-- Change value of ORIENTATION in input file'
        stop  
      else
        call single_xtal_orients(NCRDS, phi1, theta, phi2)
      end if

      call create_rmat_bunge(phi1, theta, phi2, rmat)

      STATEV(1) = rmat(1,1)
      STATEV(2) = rmat(1,2)
      STATEV(3) = rmat(1,3)
      STATEV(4) = rmat(2,1)
      STATEV(5) = rmat(2,2)
      STATEV(6) = rmat(2,3)
      STATEV(7) = rmat(3,1)
      STATEV(8) = rmat(3,2)
      STATEV(9) = rmat(3,3)
      STATEV(10) = 0.0d+0 !gamma basal
      STATEV(11) = 0.0d+0 !gamma twin 1
      STATEV(12) = 0.0d+0 !gamma twin 2
      STATEV(13) = 0.0d+0 !gamma twin 3
      STATEV(14) = 0.0d+0 !gamma twin 4
      STATEV(15) = 0.0d+0 !gamma twin 5
      STATEV(16) = 0.0d+0 !gamma twin 6
      STATEV(17) = 0.0d+0 !epsilon non-basal slip
      STATEV(18) = 0.0d+0 !epsilondot non-basal slip
      STATEV(19) = 0.0d+0 !energy per unit reference volume
      STATEV(20) = 0.0d+0 !temperature stored in isv      
      STATEV(21) = 0.0d+0 !epeff due to gambas
      STATEV(22) = 0.0d+0 !epeff due to gamtw
      STATEV(23) = 0.0d+0 !vf of all twinning
      STATEV(24) = 0.0d+0 !total effective plastic strain
      STATEV(25) = 0.0d+0 !total effective plastic strain rate      
      
C
      RETURN
      END 

c--------------------------------------------------------------------
c  ABAQUS STATE VARIABLE INITIALIZATION SUBROUTINE - POLYCRYSTAL
c--------------------------------------------------------------------

      SUBROUTINE SDVINIP(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT) 
      IMPLICIT NONE
      
C     Variables required by abaqus
      INTEGER NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT
      DOUBLE PRECISION STATEV, COORDS
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
      
C     Variables used for iteration
      double precision phi1, theta, phi2 !Bunge Euler Angle (rad)
      double precision rmat(3,3)
      integer I, NSTATVPG, NGRAINS
      
C     Variables associated with sampling random grains
      integer u, norients, nmax, iseed, iseedinit, itemp, i8_uniform
      parameter (u=20,nmax=10000, iseedinit=1111)
      double precision angs(nmax, 3)

C     Load angles for orientation files if necessary
      if (NCRDS.le.3) then !load an angle file
        iseed = iseedinit
        if(NCRDS.eq.1) then !az31 rolled data
          open (u, FILE='orients/az31roll.txt', STATUS='OLD')
        else if (NCRDS.eq.2) then !az31 ecae data
          open (u, FILE='orients/az31ecae4bc.txt', STATUS='OLD')
        else if (NCRDS.eq.3) then !amx602 data
          open (u, FILE='orients/amx602extr.txt', STATUS='OLD')
        else
          print*, '=== INVALID TEXTURE FILE CHOICE, EXITING === '
          print*, '-- Change value of ORIENTATION in input file'
          stop  
        end if
        read(u,*) norients
        call inpdat(u,norients,nmax,angs)
        close(u)      
      else
        call single_xtal_orients(NCRDS, phi1, theta, phi2)
      end if

      NSTATVPG = 31
      NGRAINS = NSTATV / NSTATVPG

      DO I=1,NGRAINS
      
        if (NCRDS.ne.0) then
          itemp = i8_uniform(iseed, norients)
          phi1 = angs(itemp, 1)
          theta = angs(itemp, 2)
          phi2 = angs(itemp, 3)
          call create_rmat_bunge(phi1, theta, phi2, rmat)
        end if

        call create_rmat_bunge(phi1, theta, phi2, rmat)

        STATEV((I-1)*NSTATVPG+1) = rmat(1,1)
        STATEV((I-1)*NSTATVPG+2) = rmat(1,2)
        STATEV((I-1)*NSTATVPG+3) = rmat(1,3)
        STATEV((I-1)*NSTATVPG+4) = rmat(2,1)
        STATEV((I-1)*NSTATVPG+5) = rmat(2,2)
        STATEV((I-1)*NSTATVPG+6) = rmat(2,3)
        STATEV((I-1)*NSTATVPG+7) = rmat(3,1)
        STATEV((I-1)*NSTATVPG+8) = rmat(3,2)
        STATEV((I-1)*NSTATVPG+9) = rmat(3,3)
        STATEV((I-1)*NSTATVPG+10) = 0.0d+0 !gamma basal
        STATEV((I-1)*NSTATVPG+11) = 0.0d+0 !gamma twin 1
        STATEV((I-1)*NSTATVPG+12) = 0.0d+0 !gamma twin 2
        STATEV((I-1)*NSTATVPG+13) = 0.0d+0 !gamma twin 3
        STATEV((I-1)*NSTATVPG+14) = 0.0d+0 !gamma twin 4
        STATEV((I-1)*NSTATVPG+15) = 0.0d+0 !gamma twin 5
        STATEV((I-1)*NSTATVPG+16) = 0.0d+0 !gamma twin 6
        STATEV((I-1)*NSTATVPG+17) = 0.0d+0 
        STATEV((I-1)*NSTATVPG+18) = 0.0d+0 
        STATEV((I-1)*NSTATVPG+19) = 0.0d+0 
        STATEV((I-1)*NSTATVPG+20) = 0.0d+0 
        STATEV((I-1)*NSTATVPG+21) = 0.0d+0 
        STATEV((I-1)*NSTATVPG+22) = 0.0d+0         
        STATEV((I-1)*NSTATVPG+23) = 0.0d+0 
        STATEV((I-1)*NSTATVPG+24) = 0.0d+0 
        STATEV((I-1)*NSTATVPG+25) = 0.0d+0  
        STATEV((I-1)*NSTATVPG+26) = 0.0d+0  !stress comp 1
        STATEV((I-1)*NSTATVPG+27) = 0.0d+0  !stress comp 2
        STATEV((I-1)*NSTATVPG+28) = 0.0d+0  !stress comp 3
        STATEV((I-1)*NSTATVPG+29) = 0.0d+0  !stress comp 4
        STATEV((I-1)*NSTATVPG+30) = 0.0d+0  !stress comp 5
        STATEV((I-1)*NSTATVPG+31) = 0.0d+0  !stress comp 6
                
      END DO

      RETURN
      END 

c====================================================================
c====================================================================
c     Determines the bunge angles based on an integer input NCRDS
c--------------------------------------------------------------------

      subroutine single_xtal_orients(NCRDS, phi1, theta, phi2)
      implicit none

c     input      
      integer NCRDS

c     output
      double precision phi1, theta, phi2

c     util
      double precision pi
      
      pi = 3.141592653589793d+0

      if (NCRDS.eq.4) then
c     c axis
c     HOSFORD ---  a (P13), e(32), f(P31)
        phi1 = 0.0d+0
        theta = 0.0d+0
        phi2 = 0.0d+0    
      else if (NCRDS.eq.5) then
c     HOSFORD ---  b (P13)
        phi1 = pi/6.0d+0 !30 deg
        theta = 0.0d+0
        phi2 = 0.0d+0    
      else if (NCRDS.eq.6) then
c     HOSFORD ---  c (P13)
        phi1 = 0.0d+0    
        theta = pi/2.0d+0 !90 deg
        phi2 = 0.0d+0    
      else if (NCRDS.eq.7) then
c     HOSFORD ---  d (P13)
        phi1 = 0.0d+0    
        theta = pi/2.0d+0 !90 deg
        phi2 = pi/6.0d+0 !30 deg
      else if (NCRDS.eq.8) then
c     HOSFORD ---  g (P13)
        phi1 = pi/2.0d+0 !90 deg
        theta = pi/4.0d+0  !45 deg
        phi2 = -pi/2.0d+0  !-90 deg
      else if (NCRDS.eq.9) then
c     slightly off c axis
        phi1 = 0.0d+0
        theta = 0.1d+0
        phi2 = 0.0d+0      
      else if (NCRDS.eq.10) then
c     random single xtal
        phi1 = 1.4d+0
        theta = 0.7d+0
        phi2 = 0.2d+0 
      else if (NCRDS.eq.11) then
        phi1 = 0.0d+0    
        theta = pi/4.0d+0 !45 deg
        phi2 = pi/6.0d+0 !30 deg
      else if (NCRDS.eq.12) then
c       ecae idealized texture 
c       143, 285, 0 or (37, 75, 0)
        phi1 = 2.4958d+0
        theta = 4.9742d+0
        phi2 = 0.0d+0    
      else
        print*, '=== INVALID TEXTURE FILE CHOICE, EXITING === '
        print*, '-- Change value of ORIENTATION in input file'
        stop  
      end if
      
      return
      end



      subroutine pumat( stress,  statev,  ddsdde,  sse,     spd,
     &                  scd,     rpl,     ddsddt,  drplde,  drpldt,
     &                  strain,  dstrain, time,    dtime,   temp,
     &                  dtemp,   predef,  dpred,   cmname,  ndi,
     &                  nshr,    ntens,   nstatv,  props,   nprops,
     &                  coords,  drot,    pnewdt,  celent,  dfgrd0,
     &                  dfgrd1,  noel,    npt,     layer,   kspt,
     &                  kstep,   kinc ) 
     
      implicit none

c     -- UMAT DELCARATIONS --       
      ! Dimension variables passed into the UMAT sub (not all are used)
      integer ndi      ! Number of direct stress components
      integer nshr     ! Number of shear stress components
      integer ntens    ! Size of stess or stran array (ndi + nshr)
      integer nstatv   ! Number of SDVs
      integer nprops   ! Number of material constants
      integer noel     ! Element number
      integer layer    ! Layer number (for composites)
      integer kspt     ! Section point number within layer
      integer kstep    ! Step number
      integer kinc     ! Increment number
      integer npt      ! Integration point number
      character*7 cmname ! Material name
      double precision sse ! Specific elastic stain energy
      
      double precision
     & celent,         ! Characteristic element length
     & dtime,          ! Time increment
     & temp,           ! Temperature at start of increment
     & dtemp,          ! Temperature increment
     & pnewdt,         ! Ratio of new time increment to time
                       ! increment being used
     & spd,            ! Specific plastic dissipation
     & scd,            ! Specific creep dissipation
     & rpl,            ! Volumetic heat generation per unit time
     & drpldt,         ! Varation of rpl with temperature
     & coords(3),      ! Coordinates of Gauss pt. being evaluated
     & ddsdde(ntens,ntens), ! Tangent Stiffness Matrix
     & ddsddt(ntens),  ! Change in stress per change in temperature
     & dfgrd1(3,3),    ! Deformation gradient at end of step
     & dfgrd0(3,3),    ! Deformation gradient at beginning of step
     & dpred(1),       ! Change in predefined state variables
     & drplde(ntens),  ! Change in heat generation per change in strain
     & drot(3,3),      ! Rotation matrix
     & dstrain(ntens), ! Strain increment tensor stored in vector form
     & predef(1),      ! Predefined state vars dependent on field
                       ! variables
     & props(nprops),  ! Material properties passed in
     & statev(nstatv), ! State Variables
     & strain(ntens),  ! Strain tensor stored in vector form
     & stress(ntens),  ! Cauchy stress tensor stored in vector form
     & time(2)         ! Step Time and Total Time

c     -- UTIL --
      integer i,j, ii, jj, nstatvpg, nstatvreg, ngrains
      parameter (nstatvreg=25, nstatvpg=31)
      double precision stressg(6), statevg(nstatvreg), ddsddeg(6,6)
      double precision temp0, tempg

      ngrains = nstatv / nstatvpg

c     initialize vars that will be averaged over
      temp0 = temp
      temp = 0.0d+0
      do i=1,6
        stress(i) = 0.0d+0
        do j=1,6
          ddsdde(i,j) = 0.0d+0
        end do
      end do

c     loop over grains
      do i=1,ngrains
c       initialize state variables for each grain
        do j=1,nstatvreg
          statevg(j) = statev((i-1)*nstatvpg+j)
        end do
        
c       access stored stress components        
        stressg(1) = statev((i-1)*nstatvpg+(nstatvreg+1))
        stressg(2) = statev((i-1)*nstatvpg+(nstatvreg+2))
        stressg(3) = statev((i-1)*nstatvpg+(nstatvreg+3))
        stressg(4) = statev((i-1)*nstatvpg+(nstatvreg+4))
        stressg(5) = statev((i-1)*nstatvpg+(nstatvreg+5))
        stressg(6) = statev((i-1)*nstatvpg+(nstatvreg+6))

c       initialize temperature to avg init temperature
        tempg = temp0
        
c       call umat
        call cumat(stressg,  statevg,  ddsddeg,  sse,     spd,
     &   scd,     rpl,     ddsddt,  drplde,  drpldt,
     &   strain,  dstrain, time,    dtime,   tempg,
     &   dtemp,   predef,  dpred,   cmname,  ndi,
     &   nshr,    ntens,   nstatv,  props,   nprops,
     &   coords,  drot,    pnewdt,  celent,  dfgrd0,
     &   dfgrd1,  noel,    npt,     layer,   kspt,
     &   kstep,   kinc)

c       record state variables from each grain
        do j=1,nstatvreg
          statev((i-1)*nstatvpg+j) = statevg(j) 
        end do
c       store calculated stress components        
        statev((i-1)*nstatvpg+nstatvreg+1) = stressg(1) 
        statev((i-1)*nstatvpg+nstatvreg+2) = stressg(2) 
        statev((i-1)*nstatvpg+nstatvreg+3) = stressg(3) 
        statev((i-1)*nstatvpg+nstatvreg+4) = stressg(4) 
        statev((i-1)*nstatvpg+nstatvreg+5) = stressg(5) 
        statev((i-1)*nstatvpg+nstatvreg+6) = stressg(6) 

c       add up variables that are averaged over
        do ii=1,6
          stress(ii) = stress(ii) + stressg(ii) 
          do jj=1,6
            ddsdde(ii,jj) = ddsdde(ii,jj) + ddsddeg(ii,jj) 
          end do
        end do
        temp = temp + tempg
      end do

c     divide averaged quantities by number of grains      
      do ii=1,6
        stress(ii) = stress(ii) / dble(ngrains)
        do jj=1,6
          ddsdde(ii,jj) = ddsdde(ii,jj)  / dble(ngrains)
        end do
      end do      
      temp = temp / dble(ngrains)
      
      return
      end     

c====================================================================
c====================================================================
c   Reads entries open file u with n vars. Nmax is the max # of entries
c     possible to read, and x is the matrix data is stored in.
c   
c   -- This is for reading in data files with 3 columns
c--------------------------------------------------------------------
      
      subroutine inpdat(u,n,nmax,x)
      implicit none
      
      integer i,n,u,nmax
      double precision x(nmax,3)

      if (n.GT.nmax) then
         write(*,*) 'Error: n = ', n, 'is larger than nmax =', nmax
         goto 9999
      endif
      
      do i= 1, n
         read(u,100) x(i,1), x(i,2), x(i,3)
      end do
  100 format (3(F10.4))

      return
 9999 stop      
      end

c====================================================================
c====================================================================
c   Given a seed, return a random integer n from 1 to max
c--------------------------------------------------------------------
      function i8_uniform( seed, n)

      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer i8_uniform
      integer seed, n

      if (seed .eq. 0 ) then
        print*, ' ----------------------'
        print*, 'i8_uniform - Fatal error!'
        print*, '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      i8_uniform = nint((dble(seed)*4.656612875D-10)*(n-1.0d+0))+1

      return
      end
