program pfasst_rttddft
  use pf_mod_mpi
  integer :: rank, err

  call mpi_init(err)
  if (err .ne. 0) stop "ERROR: Can't initialize MPI."
  call mpi_comm_rank(MPI_COMM_WORLD, rank, err)

  call rttddft()
  call mpi_finalize(err)

contains
  subroutine rttddft()
      use pfasst     !< library module defining highest-level pfasst object
      use pf_mod_mpi !< library module for mpi-related business
      use pf_mod_zndarray    !< library module containing solution type information

      use probin     !< prog-specified module for reading/parsing problem parameters
      use sweeper    !< prog-specified module containing sweeper information
      use hooks      !< prog-specified module containing program hooks

      implicit none

      type(pf_pfasst_t)  :: pf          !< pfasst data structure
      type(pf_comm_t)    :: comm        !< mpi communicator

      type(zndarray) :: dmat_t0, dmat_tfinal

      character(256) :: probin_fname       !<  file name for input
      integer    :: err, l
      real(pfdp) :: start, finish

      probin_fname = "probin.nml"
      if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
      call probin_init(probin_fname)

      call pf_mpi_create(comm, MPI_COMM_WORLD)
      call pf_pfasst_create(pf, comm, fname=probin_fname)


      !---- Create the levels -------------------------------------------------------
      do l = 1, pf%nlevels
          allocate(pf%levels(l)%shape(2))
          pf%levels(l)%shape = nparticles
          pf%levels(l)%mpibuflen = nparticles * nparticles * 2

          allocate(magpicard_context::pf%levels(l)%ulevel)
          allocate(zndarray_factory::pf%levels(l)%ulevel%factory)
          allocate(magpicard_sweeper_t::pf%levels(l)%ulevel%sweeper)

          call initialize_magpicard_sweeper(pf%levels(l)%ulevel%sweeper, &
               l, pf%qtype, pf%debug)

          if (pf%qtype == 5) then
            pf%levels(l)%nnodes = nnodes(l)+2
          else
            pf%levels(l)%nnodes = nnodes(l)
          endif

      end do

      print *,'Initializing mpi and pfasst...'
      call pf_pfasst_setup(pf)

      call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_residual)
      if (save_solutions) call pf_add_hook(pf, -1, PF_POST_CONVERGENCE, save_solution)

      call zndarray_build(dmat_t0, [nparticles,nparticles])
      call initial(dmat_t0)  
      call zndarray_build(dmat_tfinal,[nparticles,nparticles])


      call mpi_barrier(MPI_COMM_WORLD, err)

      start = MPI_Wtime()
      print*, 'Running pfasst...'
      call pf_pfasst_run(pf, dmat_t0, dt, 0.0_pfdp, nsteps, dmat_tfinal)


      finish = MPI_Wtime()
      print '("processor ", i2, " returned from pfasst CPU time: ", f16.10, " seconds")', &
           pf%rank, finish-start

      call mpi_barrier(MPI_COMM_WORLD, err)

      if(pf%rank == comm%nproc-1) then
         call dmat_tfinal%write_to_disk('sol_final') !necessary for pfasst.py
         print *,'solution at end of run'
         if (pf%debug) call dmat_tfinal%eprint() !only for debug purpose
          call dmat_tfinal%eprint() !only for debug purpose         

      endif

      call zndarray_destroy(dmat_t0)
      call zndarray_destroy(dmat_tfinal)

      call pf_pfasst_destroy(pf)

    end subroutine rttddft

  end program pfasst_rttddft
