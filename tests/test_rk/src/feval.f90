! Test feval.f90 for the RK stepper

! RHS routines for advection/diffusion example.

module feval
  use pf_mod_dtype
  use pf_mod_ndarray
  use pf_mod_imexQ
  use pf_mod_rkstepper
  implicit none

  real(pfdp), parameter :: &
       a = 1.0_pfdp, & ! "advection"
       d = 2.0_pfdp    ! "diffusion"

  ! Define the derived user_lever type
  type, extends(pf_user_level_t) :: ad_level_t

   contains

     procedure :: restrict    => restrict
     procedure :: interpolate => interpolate

  end type ad_level_t

  ! Define the derived sweeper type
  type, extends(pf_imexQ_t) :: ad_sweeper_t

   contains

     procedure :: f_eval => ad_sweeper_f_eval
     procedure :: f_comp => ad_sweeper_f_comp

  end type ad_sweeper_t

  ! Define the derived stepper type
  type, extends(pf_ark_t) :: ad_stepper_t
     
   contains
       
     procedure :: f_eval => ad_stepper_f_eval
     procedure :: f_comp => ad_stepper_f_comp
     
  end type ad_stepper_t

contains

  ! Helper functions

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

  subroutine setup(sweeper, nvars)
    class(pf_sweeper_t), intent(inout) :: sweeper
    integer,             intent(in   ) :: nvars

    class(ad_sweeper_t), pointer :: this
    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    this => as_ad_sweeper(sweeper)

    ! to be implemented

  end subroutine setup

  subroutine destroy(this, lev)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_level_t), intent(inout)   :: lev
    
    call this%imexQ_destroy(lev)

  end subroutine destroy

  ! Set initial condition.
  subroutine initial(q0)
    type(ndarray), intent(inout) :: q0

    ! to be implemented

  end subroutine initial

  ! Sweeper functions

  ! Evaluate the explicit function at y, t.
  subroutine ad_sweeper_f_eval(this, y, t, level_index, f, piece)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece    

    ! to be implemented

  end subroutine ad_sweeper_f_eval

  ! Solve for y and return f2 also.
  subroutine ad_sweeper_f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
    
    ! to be implemented

  end subroutine ad_sweeper_f_comp

  ! RK stepper functions

  ! Evaluate the explicit function at y, t.
  subroutine ad_stepper_f_eval(this, y, t, level_index, f, piece)
    class(ad_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece    

    ! to be implemented

  end subroutine ad_stepper_f_eval

  ! Solve for y and return f2 also.
  subroutine ad_stepper_f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    class(ad_stepper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dtq
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level_index
    class(pf_encap_t),   intent(inout) :: f
    integer,             intent(in   ) :: piece
    
    ! to be implemented

  end subroutine ad_stepper_f_comp

  ! User_level functions

  subroutine interpolate(this, levelF, levelG, qF, qG, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t
    
    ! to be implemented 

  end subroutine interpolate

  subroutine restrict(this, levelF, levelG, qF, qG, t)
    class(ad_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout) :: levelF
    class(pf_level_t), intent(inout) :: levelG
    class(pf_encap_t), intent(inout) :: qF, qG
    real(pfdp),        intent(in   ) :: t

    ! to be implemented

  end subroutine restrict

end module feval

