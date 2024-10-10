#include "fintrf.h"
C
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C [x, rhobeg, fopt, info] = setuptest(fun, n)

C     Declarations
      implicit none

C     mexFunction arguments:
      integer(kind=4), intent(in) :: nlhs, nrhs
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(inout) :: plhs(nlhs)

      mwSize :: mxGetM, mxGetN 
      mwPointer :: mxGetPr, mxCreateDoubleScalar 
      mwPointer :: mxCreateDoubleMatrix, mxCreateNumericMatrix
      integer(kind=4) :: mxGetString, mexPrintf
      integer(kind=4) :: mxClassIDFromClassName
      logical :: mxIsClass

      integer(kind=4) :: n, info, str_status 
      character(len=:), allocatable :: fun
      real(kind=8), allocatable :: x(:) 
      real(kind=8) :: rhobeg, fopt, n_r
      character(10000) :: STR

      if (allocated(fun)) deallocate(fun)
      allocate(character(len=mxGetN(prhs(1))) :: fun)
      str_status = mxGetString(prhs(1), fun, mxGetN(prhs(1)))
C     Check if mxGetString is successful.
      if (str_status .ne. 0) then 
         call mexErrMsgIdAndTxt ('MATLAB:setuptest:ReadError',
     +                           'Error reading Fun.')
      endif

      if (mxIsClass(prhs(2), 'double') .or. mxIsClass(prhs(2),'single'))
     + then 
         call mxCopyPtrToReal8(mxGetPr(prhs(2)), n_r, 1)
         n = int(n_r)
      else
         call mxCopyPtrToInteger4(mxGetPr(prhs(2)), n, 1)
      endif
      
      if (allocated(x)) deallocate(x)
      allocate(x(n))

      x = 1.0D308
      rhobeg = 1.0D0
      fopt = 1.0D308
      info = 0 

      call setuptest(fun, n, x, rhobeg, fopt, info)
      
      plhs(1) = mxCreateDoubleMatrix(n, 1, 0)
      call mxCopyReal8ToPtr(x, mxGetPr(plhs(1)), n)
      plhs(2) = mxCreateDoubleScalar(rhobeg)
      plhs(3) = mxCreateDoubleScalar(fopt)
      plhs(4) = mxCreateNumericMatrix(1, 1,
     +          mxClassIDFromClassName('int32'), 0)
      call mxCopyInteger4ToPtr(info, mxGetPr(plhs(4)), 1)

      if (info .eq. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:setuptest:UnknownFun', 
     +        'Unknown function name')
      endif

      if (info .eq. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:setuptest:WrongDimension', 
     +        'The test function is not defined at the given dimension')
      endif

      deallocate(fun)
      deallocate(x)

      return
      endsubroutine mexFunction
