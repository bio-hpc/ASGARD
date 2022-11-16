#ifdef NDEBUG
#  define ASSERT(e)
#else
#  define ASSERT(e) if(.not.(e))call Aass('e',__FILE__,__LINE__)
#endif

!   DEBUG( integer known_atm(natom) )
!   DEBUG( known_atm = 0 )
!   DEBUG_DO( do j=1,nneedatm ; known_atm(need_atmlist(j))=1 ; end do )

!   transforms into:

!   integer known_atm(natom)
!   known_atm = 0
!   do j=1,nneedatm ; known_atm(need_atmlist(j))=1 ; end do

!   The 2 form is necessary since the comma in the do loop is recognized by
!   the preprocessor as a second argument.

#ifdef NDEBUG
#  define DEBUG(e)
#  define DEBUG_DO(e,f)
#else
#  define DEBUG(e) e
#  define DEBUG_DO(e,f) e,f
#endif

#define REQUIRE(e) if(.not.(e))call Aass('e',__FILE__,__LINE__)
