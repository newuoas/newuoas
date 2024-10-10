#include "fintrf.h"
C
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C [f, info] = evalfun(fun, x)

C     Declarations
      implicit none

C     mexFunction arguments:
      integer(kind=4), intent(in) :: nlhs, nrhs
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(inout) :: plhs(nlhs)

      mwSize :: mxGetM, mxGetN 
      mwPointer :: mxGetPr, mxCreateDoubleScalar
      mwPointer :: mxCreateNumericMatrix
      integer(kind=4) :: mxGetString!, mexPrintf
      integer(kind=4) :: mxClassIDFromClassName

      integer(kind=4) :: n, info, allocate_status, str_status 
      character(len=:), allocatable :: fun
      real(kind=8), allocatable :: x(:) 
      real(kind=8) :: f

      if (allocated(fun)) deallocate(fun)
      allocate(character(len=mxGetN(prhs(1)))::fun)
      str_status = mxGetString(prhs(1), fun, mxGetN(prhs(1)))
C     Check if mxGetString is successful.
      if (str_status .ne. 0) then 
         call mexErrMsgIdAndTxt ('MATLAB:evalfun:ReadError',
     +                           'Error reading Fun.')
      endif

      n = mxGetM(prhs(2))
      if (allocated(x)) deallocate(x)
      allocate(x(n), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgTxt('evalfun:InsufficientMemory',
     +'allocate(x) failed.')
      endif
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), x, n)  

      f = 1.0D308
      info = 0 

      call evalfun(fun, n, x, f, info)

      plhs(1) = mxCreateDoubleScalar(f)
      plhs(2) = mxCreateNumericMatrix(1, 1,
     +          mxClassIDFromClassName('int32'), 0)
      call mxCopyInteger4ToPtr(info, mxGetPr(plhs(2)), 1)

      if (info .eq. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:evalfun:UnknownFun', 
     +        'Unknown function name')
      endif

      if (info .eq. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:evalfun:WrongDimension', 
     +        'The test function is not defined at the given dimension')
      endif

      deallocate(fun)
      deallocate(x)

      return
      endsubroutine mexFunction
