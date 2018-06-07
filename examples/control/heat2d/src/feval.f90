module feval
  use pf_mod_dtype
  use pf_mod_ndarray_oc 
  use pf_mod_restrict
  use pf_mod_imexQ_oc
  use probin
  use solutions
  
  implicit none
  
  include 'fftw3.f03'

!   real(pfdp), parameter :: pi = 3.141592653589793_pfdp    !defined in probin already
!   real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp
  
  type, extends(pf_user_level_t) :: ad_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type ad_level_t

  type, extends(pf_imexQ_oc_t) :: ad_sweeper_t  !generalize so we can use misdc as well?
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:,:)              ! work space
     complex(pfdp), allocatable ::  lap(:,:) ! operators
     real(pfdp), pointer    :: u(:,:,:)
     real(pfdp), pointer    :: ydesired(:,:,:)
     real(pfdp)             :: alpha
   contains

     procedure :: f_eval
     procedure :: f_comp
!     final :: destroy0, destroy1

  end type ad_sweeper_t

contains
  
  function as_ad_sweeper(sweeper) result(r)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(ad_sweeper_t), pointer :: r
    select type(sweeper)
    type is (ad_sweeper_t)
       r => sweeper
    class default
       stop
    end select
  end function as_ad_sweeper
  
      ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece, flags, idx, step)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece
    integer, intent(in)              :: flags
    integer, intent(in), optional    :: idx       ! index of quadrature node
    integer, intent(in), optional    :: step   ! time step for sequential version
    
    real(pfdp),    pointer :: yvec(:,:), fvec(:,:), p(:,:)
    complex(pfdp), pointer :: wk(:,:)
    integer :: l, loopstart, loopend,nx,ny

!     print *, "f_eval, piece = ", piece, "flags = ", flags, "idx = ", idx, "step = ", step
      
    select case (piece)
    case (1)  ! Explicit piece
      select case (flags)
      case (0)
        ! first component: y
        fvec => get_array2d_oc(f, 1)
!         yvec  => get_array2d_oc(y, 1)
        if (do_imex .eq. 1) then
          fvec = this%u(idx, :,:)
        else
          fvec = this%u(idx, :,:)
       endif
       ! second component: p
       yvec  => get_array2d_oc(y, 1)
       fvec => get_array2d_oc(f, 2)      
       if (do_imex .eq. 1) then
!          p  => get_array2d_oc(y, 2)
         fvec = (yvec - this%ydesired(idx,:,:))
       else
         fvec = (yvec - this%ydesired(idx,:,:))
       end if 
      case (1)
         fvec => get_array2d_oc(f, 1)
         if (do_imex .eq. 1) then
!             yvec  => get_array2d_oc(y, 1)
            fvec = this%u(idx, :,:) 
         else
            fvec = this%u(idx, :,:)
         endif
       case (2)
         ! evaluate y-y_d
         fvec => get_array2d_oc(f, 2)
         yvec  => get_array2d_oc(y, 1)
         if (do_imex .eq. 1) then
!            p  => get_array2d_oc(y, 2)
           fvec = (yvec -this%ydesired(idx,:,:)) 
         else
           fvec = (yvec -this%ydesired(idx,:,:)) 
         end if
       case default
         stop "ERROR in f_eval: only 0, 1, 2 allowed as flags"
       end select    
       
    case (2)  ! Implicit piece
       select case (flags)
       case (0)
         loopstart = 1
         loopend = 2
       case (1)
         loopstart = 1
         loopend = 1
       case (2)
         loopstart = 2
         loopend = 2
       case default
         stop "ERROR in f_eval: only 0, 1, 2 allowed as flags"
       end select
       do l = loopstart, loopend	
        yvec  => get_array2d_oc(y, l)
        fvec  => get_array2d_oc(f, l)
        nx=size(fvec,1)
        ny=size(fvec,2)
        wk => this%wk
        wk = yvec
        call fftw_execute_dft(this%ffft, wk, wk)
        wk = this%lap * wk / dble(nx*ny)
        call fftw_execute_dft(this%ifft, wk, wk)
        fvec = real(wk)
       end do
      
    case DEFAULT
      print *,'Bad case for piece in f_eval ', piece
      call exit(0)
    end select
  end subroutine f_eval

  
  ! Solve for y and return f2 also.
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f, piece, flags)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
    integer,             intent(in   ) :: flags

    real(pfdp),      pointer :: s(:,:), rhsvec(:,:), fvec(:,:)
    complex(pfdp),   pointer :: wk(:,:)
    integer :: l, loopstart, loopend,nx,ny
    
!     print *, "f_comp, piece = ", piece, "flags = ", flags, "dtq = ", dtq
    
    if (piece == 2) then
      select case (flags)
      case (0)
        loopstart = 1
        loopend = 2
      case (1)
        loopstart = 1
        loopend = 1
      case (2)
        loopstart = 2
        loopend = 2
      case default
        stop "ERROR in f2comp: only 0, 1, 2 allowed as flags"
      end select
    
      do l = loopstart, loopend	   
        s  => get_array2d_oc(y, l)
        rhsvec => get_array2d_oc(rhs, l)
        fvec => get_array2d_oc(f, l)
        
        nx = size(s,1)
        ny = size(s,2)
        
        wk => this%wk
      
        wk = rhsvec
        call fftw_execute_dft(this%ffft, wk, wk)
        wk = wk / (1.0_pfdp - dtq*this%lap) /dble(nx*ny)
        call fftw_execute_dft(this%ifft, wk, wk)
        s = real(wk)
        fvec = (s - rhsvec) / dtq
      end do
    else
      print *,'Bad piece in f_comp ',piece
      call exit(0)
    end if
  end subroutine f_comp


  
  subroutine setup(sweeper, shape, nnodes)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,     intent(in)  :: shape(2), nnodes

    class(ad_sweeper_t), pointer :: this
    integer     :: i,j,nx,ny
    type(c_ptr) :: wk
    real(pfdp)  :: kxi,kxj

    this => as_ad_sweeper(sweeper)

    nx = shape(1)
    ny = shape(2)
    
 ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(nx*ny, c_size_t))
    call c_f_pointer(wk, this%wk, [nx,ny])

    this%ffft = fftw_plan_dft_2d(nx,ny, &
         this%wk, this%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    this%ifft = fftw_plan_dft_2d(nx,ny, &
         this%wk, this%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! create operators
    allocate(this%lap(nx,ny))
    do j = 1, ny
       if (j <= ny/2+1) then
          kxj = two_pi / Ly * dble(j-1)
       else
          kxj = two_pi / Ly * dble(-ny + j - 1)
       end if
       do i = 1, nx
       if (i <= nx/2+1) then
          kxi = two_pi / Lx * dble(i-1)
       else
          kxi = two_pi / Lx * dble(-nx + i - 1)
       end if

       if (kxi**2+kxj**2 < 1e-13) then
          this%lap(i,j) = 0.0_pfdp
       else
          this%lap(i,j) = -(kxi**2+kxj**2)
       end if
    end do
    end do

!     print *, this%lap
    
    ! allocate control and desired state
    allocate(this%u(nnodes,nx,ny))
    allocate(this%ydesired(nnodes,nx,ny))
  end subroutine setup

  
  subroutine destroy(this, lev)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_level_t), intent(inout)   :: lev

    deallocate(this%wk)
    deallocate(this%lap)
    call fftw_destroy_plan(this%ffft)
    call fftw_destroy_plan(this%ifft)
    
    deallocate(this%u)
    deallocate(this%ydesired)
    
    call this%imexQ_oc_destroy(lev)

  end subroutine destroy

  
  subroutine compute_exp_lap_dt_times_vec(s, dt, sol)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),          intent(in)    :: dt
    class(pf_encap_t),   intent(inout) :: sol

    integer :: nx, ny, i, j    
    real(pfdp),          pointer :: vec(:,:), expdtL(:,:)
    complex(pfdp),       pointer :: wk(:,:)
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)    
    vec     => get_array2d_oc(sol, 2) ! for adjoint 
    nx = size(vec, 1)
    ny = size(vec, 2)
    
    allocate(expdtL(nx,ny))
    do i=1, nx
       do j=1, ny
          expdtL(i,j) = exp(dt*(sweeper%lap(i,j)))
       end do
    end do
    
    wk => sweeper%wk
    wk = vec
    call fftw_execute_dft(sweeper%ffft, wk, wk)
    wk = expdtL * wk / dble(nx*ny) 
    call fftw_execute_dft(sweeper%ifft, wk, wk)
    vec = real(wk)
    
    deallocate(expdtL)
  end subroutine compute_exp_lap_dt_times_vec

  
  
  subroutine initialize_oc(s, t0, dt, nodes, shape)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp), intent(in)    :: nodes(:), t0, dt
    integer,             intent(in)    :: shape(2)

    integer :: nnodes, m
    real(pfdp) :: t
    
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nnodes = size(nodes)
    
    do m = 1, nnodes
       t = t0 + dt*nodes(m)
!        print *, t
       sweeper%u(m,:,:) = 0.0_pfdp
!        call exact_u(sweeper%u(m,:,:), shape, t)
       call y_desired(sweeper%ydesired(m,:,:), shape, t)
    end do
    sweeper%alpha = alpha
  end subroutine initialize_oc
  
  
  subroutine dump_control(s, pf, fbase)
    class(pf_sweeper_t), intent(inout) :: s
    type(pf_pfasst_t),  intent(inout) :: pf
    character(len = *), intent(in   ) :: fbase

    integer :: nnodes, m
    character(len=256)     :: fname
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nnodes = pf%levels(pf%nlevels)%nnodes
    do m = 1, nnodes
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase), pf%rank, m

      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), sweeper%u(m,:,:))
    end do
  end subroutine dump_control
  
  subroutine dump_ydesired(s, pf, t0, dt, fbase, fbase2, fbase3, fbase4, fbase5)
    class(pf_sweeper_t), intent(inout) :: s
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),  intent(in) :: t0, dt
    character(len = *), intent(in   ) :: fbase, fbase2, fbase3, fbase4, fbase5

    integer :: nnodes, i
    character(len=256)     :: fname
    real(pfdp), pointer :: y(:,:), ydiff(:,:), ymyd(:,:)
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nnodes = pf%levels(pf%nlevels)%nnodes
    do i = 1, nnodes
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase), pf%rank, i

      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), sweeper%ydesired(i,:,:))
           
      y => get_array2d_oc(pf%levels(pf%nlevels)%Q(i), 1)
      allocate(ydiff(size(y,1),size(y,2)))
      allocate(ymyd(size(y,1),size(y,2)))
      ymyd(:,:) = y(:,:)-sweeper%ydesired(i,:,:)
      
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase2), pf%rank, i
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), ymyd(:,:))
      
      call ymyd_exact(ydiff, shape(ydiff), t0+dt*pf%levels(pf%nlevels)%nodes(i))
      ydiff = ydiff - ymyd
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase3), pf%rank, i
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), ydiff(:,:))
      
      y => get_array2d_oc(pf%levels(pf%nlevels)%F(i,1), 2)
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase4), pf%rank, i
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), y(:,:))
           
!       call lap_p_tilde(ydiff, shape(ydiff), t0+dt*pf%levels(pf%nlevels)%nodes(i), t0+dt)
      call ymyd_exact(ydiff, shape(ydiff), t0+dt*pf%levels(pf%nlevels)%nodes(i))
      ydiff = ydiff - y
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase5), pf%rank, i
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), ydiff(:,:))
      
      deallocate(ymyd)
      deallocate(ydiff)
      
      
    end do
  end subroutine dump_ydesired
  

  
    subroutine dump_exact_adjoint(s, pf, t0, dt, fbase, fbase2, fbase3, L2error, LinfError, L2exact, LinfExact, L2norm, LinfNorm)
    class(pf_sweeper_t), intent(inout) :: s
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp), intent(in)    ::  t0, dt
    character(len = *), intent(in   ) :: fbase, fbase2, fbase3
    real(pfdp),           intent(  out) :: L2error, LinfError, L2exact, LinfExact, L2norm, LinfNorm

    real(pfdp), pointer :: pexact(:,:), pdiff(:,:), ex(:), diff(:), p(:,:), nodes(:), nrm(:)
    
    
    integer :: nnodes, m, i, j, nx, ny
    character(len=256)     :: fname
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nx=pf%levels(pf%nlevels)%shape(1)
    ny=pf%levels(pf%nlevels)%shape(2)
    nnodes = pf%levels(pf%nlevels)%nnodes
    
    allocate(pexact(nx,ny))   
    allocate(pdiff(nx,ny)) 
    allocate(ex(nnodes))
    allocate(diff(nnodes)) 
    allocate(nrm(nnodes)) 
    
    allocate(nodes(nnodes))    
    nodes = pf%levels(pf%nlevels)%nodes

    do m = 1, nnodes
      pexact = 0.0_pfdp
      call exact_p(pexact, shape(pexact), t0+dt*pf%levels(pf%nlevels)%nodes(m))
!       call p_tilde(pexact, shape(pexact), t0+dt*pf%levels(pf%nlevels)%nodes(m), t0+dt)
!       print *, pf%rank, m, t0+dt*pf%levels(pf%nlevels)%nodes(m), t0+dt
      
      
      p => get_array2d_oc(pf%levels(pf%nlevels)%Q(m), 2)
      
      pdiff(:,:) = pexact(:,:) - p(:,:)
      
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase), pf%rank, m

      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), pexact(:,:))
                                 
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase2), pf%rank, m
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), p(:,:))
           
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase3), pf%rank, m
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), pdiff(:,:))
      ! compute norm
      
      ex(m) = 0.0_pfdp
      diff(m) = 0.0_pfdp
      nrm(m) = 0.0_pfdp
      do j = 1, ny
      do i = 1, nx
         ex(m) = ex(m) + pexact(i,j)**2.0
         diff(m) = diff(m) + pdiff(i,j)**2.0
         nrm(m) = p(i,j)**2
      end do
      end do
      ex(m) = ex(m)*Lx*Ly/dble(nx*ny)
      diff(m) = diff(m)*Lx*Ly/dble(nx*ny)
      nrm(m) = nrm(m)*Lx*Ly/dble(nx*ny)
      
      LinfExact =  max(maxval(abs(pexact(:,:))), LinfExact)
      LinfError =  max(maxval(abs(pdiff(:,:))), LinfError)
      LinfNorm =  max(maxval(abs(p(:,:))), LinfNorm)
    end do

    L2error = 0.0
    L2exact = 0.0
    L2norm = 0.0
    do m=1, nnodes-1
       L2error = L2error + (diff(m)+diff(m+1))*(nodes(m+1)-nodes(m))*dt
       L2exact = L2exact + (ex(m)+ex(m+1))*(nodes(m+1)-nodes(m))*dt
       L2norm = L2norm + (nrm(m)+nrm(m+1))*(nodes(m+1)-nodes(m))*dt
    end do

    L2error = 0.5*L2error
    L2exact = 0.5*L2exact    
    L2norm  = 0.5*L2norm
    
    deallocate(nodes)
    deallocate(nrm)
    deallocate(diff)
    deallocate(ex)
    deallocate(pdiff)
    deallocate(pexact)
  end subroutine dump_exact_adjoint

  
  
  subroutine dump_exact_state(s, pf, t0, dt, fbase, fbase2, fbase3, L2error, LinfError, L2exact, LinfExact)
    class(pf_sweeper_t), intent(inout) :: s
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp), intent(in)    ::  t0, dt
    character(len = *), intent(in   ) :: fbase, fbase2, fbase3
    real(pfdp),           intent(  out) :: L2error, LinfError, L2exact, LinfExact

    real(pfdp), pointer :: yexact(:,:), ydiff(:,:), ex(:), diff(:), y(:,:), nodes(:)
    
    
    integer :: nnodes, m, i, j, nx, ny
    character(len=256)     :: fname
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nx=pf%levels(pf%nlevels)%shape(1)
    ny=pf%levels(pf%nlevels)%shape(2)
    nnodes = pf%levels(pf%nlevels)%nnodes
    
    allocate(yexact(nx,ny))   
    allocate(ydiff(nx,ny)) 
    allocate(ex(nnodes))
    allocate(diff(nnodes)) 

    allocate(nodes(nnodes))    
    nodes = pf%levels(pf%nlevels)%nodes


    do m = 1, nnodes
      yexact = 0.0_pfdp
      call exact_y(yexact, shape(yexact), t0+dt*pf%levels(pf%nlevels)%nodes(m))
      
      y => get_array2d_oc(pf%levels(pf%nlevels)%Q(m), 1)
      
      ydiff(:,:) = yexact(:,:) - y(:,:)
      
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase), pf%rank, m

      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), yexact(:,:))
           
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase2), pf%rank, m
                 
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), y(:,:))
           
      write(fname, "(A,'r',i0.2,'m',i0.2,'.npy')") trim(fbase3), pf%rank, m
      
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), ydiff(:,:))
      ! compute norm
      
      ex(m) = 0.0_pfdp
      diff(m) = 0.0_pfdp
      do j = 1, ny
      do i = 1, nx
         ex(m) = ex(m) + yexact(i,j)**2.0
         diff(m) = diff(m) + ydiff(i,j)**2.0
      end do
      end do
      ex(m) = ex(m)*Lx*Ly/dble(nx*ny)
      diff(m) = diff(m)*Lx*Ly/dble(nx*ny)
      
      LinfExact =  max(maxval(abs(yexact(:,:))), LinfExact)
      LinfError =  max(maxval(abs(ydiff(:,:))), LinfError)
    end do

    L2error = 0.0
    L2exact = 0.0
    do m=1, nnodes-1
       L2error = L2error + (diff(m)+diff(m+1))*(nodes(m+1)-nodes(m))*dt
       L2exact = L2exact + (ex(m)+ex(m+1))*(nodes(m+1)-nodes(m))*dt
    end do

    L2error = 0.5*L2error
    L2exact = 0.5*L2exact    
    
    deallocate(nodes)
    deallocate(diff)
    deallocate(ex)
    deallocate(ydiff)
    deallocate(yexact)
  end subroutine dump_exact_state
  
  
  
  subroutine dump_exact_control(s, pf, t0, dt, nodes, L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfExactCtrl)
    type(pf_pfasst_t),    intent(inout) :: pf
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),           intent(in   ) :: nodes(:), t0, dt
    real(pfdp),           intent(  out) :: L2errorCtrl, LinfErrorCtrl, L2exactCtrl, LinfExactCtrl

    integer                :: nnodes, i,j,m,nx,ny
    character(len=256)     :: fname
    real(pfdp), pointer    :: uexact(:,:), udiff(:,:), ex(:), diff(:)
    complex(pfdp), pointer :: wk(:,:)

    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nx=pf%levels(pf%nlevels)%shape(1)
    ny=pf%levels(pf%nlevels)%shape(2)

    nnodes = size(nodes,1)
    
    allocate(uexact(nx,ny))   
    allocate(udiff(nx,ny))   
    allocate(ex(nnodes))
    allocate(diff(nnodes))   
    
    LinfErrorCtrl = 0.0_pfdp
    LinfExactCtrl = 0.0_pfdp
    
    do m = 1, nnodes
      call exact_u(uexact, shape(uexact), t0+dt*nodes(m))

      write(fname, "('uexact','r',i0.2,'m',i0.2,'.npy')") pf%rank, m
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), uexact)
           
      udiff(:,:) = uexact(:,:) - sweeper%u(m,:,:)
      write(fname, "('udiff','r',i0.2,'m',i0.2,'.npy')") pf%rank, m
      call ndarray_dump_numpy(trim(pf%outdir)//c_null_char,trim(fname)//c_null_char, '<f8'//c_null_char//c_null_char, &
           ndim, pf%levels(pf%nlevels)%shape, product(pf%levels(pf%nlevels)%shape), udiff(:,:))
           
      ! compute norm
      ex(m) = 0.0_pfdp
      diff(m) = 0.0_pfdp
      do j = 1, ny
      do i = 1, nx
         ex(m) = ex(m) + uexact(i,j)**2.0
         diff(m) = diff(m) + udiff(i,j)**2.0
      end do
      end do
      ex(m) = ex(m)*Lx*Ly/dble(nx*ny)
      diff(m) = diff(m)*Lx*Ly/dble(nx*ny)
      
      LinfExactCtrl =  max(maxval(abs(uexact(:,:))), LinfExactCtrl)
      LinfErrorCtrl =  max(maxval(abs(udiff(:,:))), LinfErrorCtrl)
    end do

    L2errorCtrl = 0.0
    L2exactCtrl = 0.0
    do m=1, nnodes-1
       L2errorCtrl = L2errorCtrl + (diff(m)+diff(m+1))*(nodes(m+1)-nodes(m))*dt
       L2exactCtrl = L2exactCtrl + (ex(m)+ex(m+1))*(nodes(m+1)-nodes(m))*dt
    end do

    L2errorCtrl = 0.5*L2errorCtrl
    L2exactCtrl = 0.5*L2exactCtrl

    deallocate(ex)
    deallocate(diff)          
    deallocate(uexact)
    deallocate(udiff)          
  end subroutine dump_exact_control
  
  
  subroutine construct_gradient(s, grad, nodes, LinftyNormGrad, L2NormGradSq)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),          intent(inout) :: grad(:,:,:)
    real(pfdp),          intent(in)    :: nodes(:)
    real(pfdp),          intent(out)   :: LinftyNormGrad, L2NormGradSq
    
    integer              :: m, i,j, nnodes, nx,ny
    real(pfdp),  pointer :: obj(:)
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)

    nnodes = size(grad,1)
    nx  = size(grad,2)
    ny  = size(grad,3)
    allocate(obj(nnodes))

    LinftyNormGrad =  0
    L2NormGradSq = 0
    do m = 1, nnodes
       grad(m,:,:) = grad(m,:,:) + sweeper%alpha * sweeper%u(m,:,:)
       LinftyNormGrad = max(maxval(abs(grad(m,:,:))), LinftyNormGrad)  
       !obj(m) = sum((grad(m,:)**2.0))*Lx/nvars
       obj(m) = 0.0_pfdp
       do j = 1, nx
          do i = 1, ny
             obj(m) = obj(m) + grad(m,i,j)**2.0
          end do
       end do
       obj(m) = obj(m)*Lx*Ly/dble(nx*ny)
     end do

    L2NormGradSq = 0.0
    do m=1, nnodes-1
       L2normGradSq = L2normGradSq + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
    end do

    L2NormGradSq = 0.5*L2NormGradSq !0.5 for trapezoidal rule

    deallocate(obj)
  end subroutine construct_gradient
  

  subroutine update_control(s, delta, step)
    class(pf_sweeper_t), intent(inout) :: s
    real(pfdp),          intent(in)    :: delta(:,:,:), step
    
    integer              :: m, nnodes   
    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nnodes = size(delta,1)

    do m = 1, nnodes
       sweeper%u(m,:,:) = sweeper%u(m,:,:) + step * delta(m,:,:)
    end do
 end subroutine update_control




  subroutine control_L2Q(s, dt, nodes, shape, L2normSq)
    class(pf_sweeper_t), intent(inout) :: s
    integer,             intent(in)    :: shape(2)
    real(pfdp),          intent(in)    :: nodes(:), dt
    real(pfdp),          intent(out)   :: L2normSq
    
    real(pfdp),  pointer   :: obj(:)
    integer                :: m, i,j,nx,ny, nnodes
 
    class(ad_sweeper_t), pointer :: sweeper;
    sweeper => as_ad_sweeper(s)
    nx=shape(1)
    ny=shape(2)
    
    nnodes = size(nodes)
    allocate(obj(nnodes))

    do m=1, nnodes
       !print *, 'sum u/nvars = ', sum(work%u(m,:,:))*Lx/nvars, 'sum u2/nvars = ', sum(work%u(m,:,:)**2.0)*Lx/nvars, &
       !          maxval(work%u(m,:,:)), minval(work%u(m,:,:))
       obj(m) = 0.0_pfdp
       !obj(m) = sum(work%u(m,:,:)**2.0)*Lx/nvars !rectangle rule
       do j = 1, ny
          do i = 1, nx
             obj(m) = obj(m) + sweeper%u(m,i,j)**2.0
          end do
       end do
       obj(m) = obj(m)*Lx*Ly/dble(nx*ny)
    end do

    L2normSq = 0.0
    do m=1, nnodes-1
       L2normSq = L2normSq + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
    end do

    L2normSq = 0.5*L2normSq !0.5 for trapezoidal rule

    deallocate(obj)
  end subroutine control_L2Q

  
  subroutine objective_function(s, sol, shape, m, objective) 
  ! actually just the tracking part of the objective
    class(pf_sweeper_t), intent(inout) :: s
    class(pf_encap_t), intent(in   )   :: sol
    integer,             intent(in)    :: shape(2), m
    real(pfdp),          intent(out)   :: objective
    
    real(pfdp),  pointer   :: y(:,:), f(:,:) !, obj(:)
    integer                :: nx,ny,i,j !, nnodes
    complex(pfdp), pointer :: wk(:,:)

    class(ad_sweeper_t), pointer :: sweeper
    sweeper => as_ad_sweeper(s)
    
    nx = shape(1)
    ny = shape(2)
    allocate(f(nx,ny))

       y => get_array2d_oc(sol, 1)
       f = (y -sweeper%ydesired(m,:,:))   
       objective = 0.0_pfdp
       do i = 1, nx
          do j = 1, ny
             objective = objective + f(i,j)**2
          end do
       end do
       objective = objective * Lx*Ly / dble(nx*ny)

    deallocate(f)
  end subroutine objective_function
  
  
  ! this shouldbe in a separate module
  function compute_scalar_prod(f, g, nodes, dt) result(r)
    real(pfdp), intent(in)  :: f(:,:,:), g(:,:,:), nodes(:), dt
    real(pfdp)              :: r
    real(pfdp), pointer     :: obj(:)
    integer                 :: nnodes,  m, i,j,nx,ny
    r = 0.0_pfdp
    nnodes = size(f,1)
    nx  = size(f,2)
    ny  = size(f,3)
    allocate(obj(nnodes))
    do m = 1, nnodes
      obj(m) = 0.0_pfdp
      do i = 1, nx
        do j = 1, ny
        obj(m) = obj(m) + f(m,i,j)*g(m,i,j) 
      end do
      end do
      obj(m) = obj(m)*Lx*Ly/dble(nx*ny)
    end do
    do m=1, nnodes-1
      r = r + (obj(m)+obj(m+1))*(nodes(m+1)-nodes(m))*dt
    end do
    r = r * 0.5
    deallocate(obj)
  end function compute_scalar_prod

  
  subroutine restrict_control(sG, sF)
    class(pf_sweeper_t), intent(inout) :: sG, sF
        
    real(pfdp), pointer  :: uF(:,:,:), uG(:,:,:)
    integer :: nvarF, nvarG, xrat, nnodesF, nnodesG, trat, m
    
    class(ad_sweeper_t), pointer :: sweeperF, sweeperG
    sweeperF => as_ad_sweeper(sF)
    sweeperG => as_ad_sweeper(sG)

    uF => sweeperF%u
    uG => sweeperG%u

    nnodesF = size(uF,1)
    nnodesG = size(uG,1)
    nvarF = size(uF,2)
    nvarG = size(uG,2)

    xrat  = nvarF / nvarG
    trat  = ceiling(real(nnodesF) / real(nnodesG))
!     print *, 'restrict u', xrat, trat

    !do m=1,nnodesG
       uG(:,:,:) = uF(::trat,::xrat,::xrat)
    !end do
  end subroutine restrict_control

  
  subroutine restrict_ydesired(sG, sF)
    class(pf_sweeper_t), intent(inout) :: sG, sF
    
    real(pfdp), pointer  :: ydesiredF(:,:,:), ydesiredG(:,:,:)
    integer :: nvarF, nvarG, xrat, nnodesF, nnodesG, trat, m
    
    class(ad_sweeper_t), pointer :: sweeperF, sweeperG
    sweeperF => as_ad_sweeper(sF)
    sweeperG => as_ad_sweeper(sG)

    ydesiredF => sweeperF%ydesired
    ydesiredG => sweeperG%ydesired

    nnodesF = size(ydesiredF,1)
    nnodesG = size(ydesiredG,1)
    nvarF = size(ydesiredF,2)
    nvarG = size(ydesiredG,2)

    xrat  = nvarF / nvarG
    trat  = ceiling(real(nnodesF) / real(nnodesG))
!     print *, 'restrict ydesired', xrat, trat

    !do m=1,nnodesG
       ydesiredG(:,:,:) = ydesiredF(::trat,::xrat,::xrat)
    !end do
  end subroutine restrict_ydesired
  
  
  subroutine restrict_for_adjoint(pf, which)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: which
    real(pfdp), pointer :: tF(:), tG(:) !zF(:), zG(:)
    integer :: l, m !, nnodesF, nnodesG, nvarsF, nvarsG
    
    do l = pf%nlevels, 2, -1
      allocate(tF(pf%levels(l)%nnodes))
      allocate(tG(pf%levels(l-1)%nnodes))
      tF = pf%state%t0 + pf%state%dt*pf%levels(l)%nodes
      tG = pf%state%t0 + pf%state%dt*pf%levels(l-1)%nodes
      call restrict_sdc(pf%levels(l), pf%levels(l-1), pf%levels(l)%Q, pf%levels(l-1)%Q, .false. ,tF, which)

        call pf%levels(l-1)%ulevel%sweeper%evaluate_all(pf%levels(l-1), tG, which)
      deallocate(tF)
      deallocate(tG)
    end do
  end subroutine restrict_for_adjoint
  
    
  subroutine interp2(qF, qG, adF, adG)
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),  intent(inout) :: qF(:,:), qG(:,:)

    complex(pfdp), pointer :: wkF(:,:), wkG(:,:)
    integer      :: nxF, nxG, nyF, nyG, xrat,yrat,i,j,ii,jj
    nxF=size(qF,1)
    nyF=size(qF,2)
    nxG=size(qG,1)
    nyG=size(qG,2)
    xrat  = nxF/nxG
    yrat  = nyF/nyG
    

    if (xrat == 1 .and. yrat==1) then
       qF = qG
       return
    endif

    
    wkF => adF%wk
    wkG => adG%wk
       
    wkG = qG
    call fftw_execute_dft(adG%ffft, wkG, wkG)
    wkG = wkG / (nxG*nyG)
       
    wkF = 0.0d0
    do j = 1, nyG
       if (j <= nyG/2) then
          jj = j
       else if (j > nyG/2+1) then
          jj = nyF - nyG + j
       else
          cycle
       end if
       
       do i = 1, nxG
          if (i <= nxG/2) then
             ii = i
          else if (i > nxG/2+1) then
             ii = nxF - nxG + i
          else
             cycle
          end if
          
          wkF(ii, jj) = wkG(i, j)
       end do
    end do

    call fftw_execute_dft(adF%ifft, wkF, wkF)
!     wkG = qG
!     call fftw_execute_dft(adG%ffft, wkG, wkG)
!     wkG = wkG / (nxG*nyG)
! 
!     wkF = 0
!     do j = 1, nyG/2
!        do i = 1, nxG/2
!           wkF(i,j)       = wkG(i,j)
!        end do
!        do i = 0, nxG/2-2
!           wkF(nxF-i,j) = wkG(nxG-i,j)
!        end do
!     end do
!     do j = 0, nyG/2-2
!        do i = 1, nxG/2
!           wkF(i,nyF-j)       = wkG(i,nyG-j)
!        end do
!        do i = 0, nxG/2-2
!           wkF(nxF-i,nyF-j) = wkG(nxG-i,nyG-j)
!        end do
!     end do
! 
!     call fftw_execute_dft(adF%ifft, wkF, wkF)
       
    qF = real(wkF)    
  end subroutine interp2
    
    
  subroutine interpolate(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags
    
    integer :: which
    class(ad_sweeper_t), pointer :: adF, adG
    real(pfdp),          pointer :: f(:,:), g(:,:)
             
    which = 0
    if(present(flags)) which = flags
    if(.not.present(flags)) print *, "interpolate without flags"

    adG => as_ad_sweeper(levelG%ulevel%sweeper)
    adF => as_ad_sweeper(levelF%ulevel%sweeper)
  
    if ((which .eq. 0) .or. (which .eq. 1)) then
      f => get_array2d_oc(qF,1)
      g => get_array2d_oc(qG,1)
      call interp2(f, g, adF, adG)  
    end if
    if ((which .eq. 0) .or. (which .eq. 2)) then
      f => get_array2d_oc(qF,2)
      g => get_array2d_oc(qG,2)
      call interp2(f, g, adF, adG)
    end if
  end subroutine interpolate
  

  subroutine restrict(this, levelF, levelG, qF, qG, t, flags)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    integer, intent(in), optional :: flags

    real(pfdp), pointer :: f(:,:), g(:,:)

    integer :: NxF, NxG,NyF, NyG, xrat, yrat,which
    which = 1
    if(present(flags)) which = flags
    if(.not.present(flags)) print *, "restrict without flags"
    
    if ((which .eq. 0) .or. (which .eq. 1)) then
      f => get_array2d_oc(qF,1)
      g => get_array2d_oc(qG,1)
      NxF=size(f,1)
      NyF=size(f,2)
      NxG=size(g,1)
      NyG=size(g,2)
      xrat  = NxF/NxG
      yrat  = NyF/NyG
      g = f(::xrat,::yrat)
    end if
    if ((which .eq. 0) .or. (which .eq. 2)) then
      f => get_array2d_oc(qF,2)
      g => get_array2d_oc(qG,2)
      NxF=size(f,1)
      NyF=size(f,2)
      NxG=size(g,1)
      NyG=size(g,2)
      xrat  = NxF/NxG
      yrat  = NyF/NyG
      g = f(::xrat,::yrat)
    end if    
  end subroutine restrict
 
end module feval