! Author: Edoardo ascani
!
! Description: 
! Molecular dynamics program using Velocity Verlet integration and Lennard-Jones potential.

module utilities

    public :: lecture, get_distances, get_forces, get_kinetic, get_potential
  
  contains
  
    subroutine lecture(filename, positions, atom, velocities, nstep, mass, timestep, sig, eps, natom)
      implicit none
      integer :: i, nstep, natom
      double precision :: timestep, mass, eps, sig
      double precision, allocatable, dimension(:,:) :: positions, velocities
      character(len=100) :: filename
      character(len=2), allocatable, dimension(:) :: atom
  
      open(1, file=filename, status="old")
      read(1,*) 
      read(1,*) nstep
      read(1,*) timestep
      read(1,*) natom
  
      allocate(positions(natom,3))
      allocate(atom(natom))
      allocate(velocities(natom,3))
  
      read(1,*) mass
      do i = 1, natom
        read(1,*) atom(i), positions(i,:)
      end do
      read(1,*) sig
      read(1,*) eps
      do i = 1, natom
        read(1,*) velocities(i,:)
      end do
  
      close(1)
    end subroutine lecture
  
    subroutine get_distances(positions, distances, natom)
      implicit none
      integer, intent(in) :: natom
      double precision, intent(in) :: positions(natom,3)
      double precision, intent(out) :: distances(natom,natom)
      integer :: i, j
  
      distances = 0.0d0
      do i = 1, natom
        do j = 1, natom
          distances(i,j) = dsqrt((positions(i,1) - positions(j,1))**2 + &
                                 (positions(i,2) - positions(j,2))**2 + &
                                 (positions(i,3) - positions(j,3))**2)
        end do
      end do
    end subroutine get_distances
  
    subroutine get_forces(forces, distances, positions, sig, eps, natom)
      implicit none
      integer, intent(in) :: natom
      double precision, intent(in) :: eps, sig
      double precision, intent(inout) :: forces(natom,3)
      double precision, intent(in) :: positions(natom,3), distances(natom,natom)
      integer :: i, j, k
  
      forces = 0.d0
      do i = 1, natom
        do k = 1, 3
          do j = 1, natom
            if (i == j) cycle
            forces(i,k) = forces(i,k) + 48*eps*(positions(i,k)-positions(j,k)) * &
                          (sig**12/distances(i,j)**14 - 0.5d0*sig**6/distances(i,j)**8)
          end do
        end do
      end do
    end subroutine get_forces
  
    subroutine get_kinetic(velocities, natom, kin_ener, mass)
      implicit none
      double precision, intent(in) :: velocities(natom,3), mass
      double precision, intent(out) :: kin_ener
      integer, intent(in) :: natom
      integer :: i, k
  
      kin_ener = 0.0d0
      do i = 1, natom
        do k = 1, 3
          kin_ener = kin_ener + 0.5d0 * mass * velocities(i,k)**2
        end do
      end do
    end subroutine get_kinetic
  
    subroutine get_potential(distances, natom, sig, eps, pot_ener)
      implicit none
      double precision, intent(in) :: distances(natom,natom), sig, eps
      double precision, intent(out) :: pot_ener
      integer, intent(in) :: natom
      integer :: i, j
  
      pot_ener = 0.0d0
      do i = 1, natom-1
        do j = i+1, natom
          pot_ener = pot_ener + 4.d0 * eps * ((sig/distances(i,j))**12 - (sig/distances(i,j))**6)
        end do
      end do
    end subroutine get_potential
  
  end module utilities
  
  program md
    use utilities
    implicit none
  
    double precision, parameter  :: kcal_to_Hartree = 0.0015936015d0
    double precision, parameter  :: ev_to_Hartree = 0.0367492929d0
    double precision, parameter  :: ang_to_bohr = 1.8897261d0
    double precision, parameter  :: fs_to_timeau = 41.34137314d0
    double precision, parameter  :: amu_to_au = 1822.8884850d0
    double precision, parameter  :: timeau_to_fs = 0.02418884254d0
  
    integer :: i, j, step, nstep, natom
    double precision :: timestep, mass, eps, sig, kin_ener, pot_ener, energy, time
    double precision, allocatable :: positions(:,:), velocities(:,:), distances(:,:), forces(:,:)
    character(len=100) :: filename = "input"
    character(len=2), allocatable :: atom(:)
  
    open(10, file="md.xyz")
    open(11, file="energy.dat")
  
    call lecture(filename, positions, atom, velocities, nstep, mass, timestep, sig, eps, natom)
  
    allocate(forces(natom,3), distances(natom,natom))
  
    ! Conversion to atomic units
    mass = mass * amu_to_au
    timestep = timestep * fs_to_timeau
    sig = sig * ang_to_bohr
    eps = eps * ev_to_Hartree
    positions = positions * ang_to_bohr
    velocities = velocities * ang_to_bohr / fs_to_timeau
  
    call get_distances(positions, distances, natom)
    call get_forces(forces, distances, positions, sig, eps, natom)
  
    call get_kinetic(velocities, natom, kin_ener, mass)
    call get_potential(distances, natom, sig, eps, pot_ener)
    energy = kin_ener + pot_ener
  
    step = 0
    time = step * timestep * timeau_to_fs
    write(10,"(I3)") natom
    write(10,"(A,I4,A,F6.3,3(A,E15.6))") "# step: ", step, ",  t(fs): ", time, &
                                       & ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener
    do i = 1, natom
       write(10,*) atom(i), positions(i,:)
    end do
    write(11,*) step, time, energy, kin_ener, pot_ener
  
    do step = 1, nstep
      do i = 1, natom
        do j = 1, 3
          velocities(i,j) = velocities(i,j) + 0.5d0 * timestep * forces(i,j) / mass
        end do
      end do
  
      do i = 1, natom
        do j = 1, 3
          positions(i,j) = positions(i,j) + timestep * velocities(i,j)
        end do
      end do
  
      call get_distances(positions, distances, natom)
      call get_forces(forces, distances, positions, sig, eps, natom)
  
      do i = 1, natom
        do j = 1, 3
          velocities(i,j) = velocities(i,j) + 0.5d0 * timestep * forces(i,j) / mass
        end do
      end do
  
      call get_kinetic(velocities, natom, kin_ener, mass)
      call get_potential(distances, natom, sig, eps, pot_ener)
      energy = kin_ener + pot_ener
  
      time = step * timestep * timeau_to_fs
      write(10,"(I3)") natom                                
      write(10,"(A,I4,A,F9.3,3(A,E15.6))") "# step: ", step, ",  t(fs): ", time, &
                                           & ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener
      do i = 1, natom
         write(10,*) atom(i), positions(i,:)
      end do
      write(11,*) step, time, energy, kin_ener, pot_ener
    end do
  
    close(10)
    close(11)
  
  end program md
  