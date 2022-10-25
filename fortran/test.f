      PROGRAM main
        REAL(kind = 8) :: A(3,3), LAM(3), V(3,3)
        A(1,1) = 1 
        A(1,2) = 1 
        A(1,3) = 0
        A(2,1) = 10 
        A(2,2) = 1
        A(2,3) = 0
        A(3,1) = 1
        A(3,2) = 3
        A(3,3) = 1
        CALL EIG(A, 3, 3, LAM, V)
        write(*,*) A
        write(*,*) LAM 
        write(*,*) V
      END
