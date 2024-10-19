!----------------------------------------------------------
! This is a subset of the module matrix used in music. 
! Used here only to calculate the inverse of a 3x3 matrix.
!----------------------------------------------------------
Module matrix

  Use defaults, Only: RDbl
  Use matrixops, Only: matrixops_ludcmp, matrixops_lubksb

  Implicit None
  Save

  Private
  Public :: MatrixType, matrix_inverse

  Integer, Parameter   :: nrows = 3
  Integer, Parameter   :: ncols = 3
  
  Type MatrixType
    Real(kind=RDbl), Dimension(nrows, ncols)  :: comp
  End Type MatrixType


Contains



  !----------------------------------------------------------------------------
  ! Inverts a 3x3 matrix numerically
  !----------------------------------------------------------------------------
  Type(MatrixType) Function matrix_inverse(m1)
    Implicit None
    Type(MatrixType), Intent(in) :: m1

    Integer                      :: i,j
    !variables for the NumRec subroutines
    Integer                      :: np,n
    Integer, Dimension(3)        :: indx
    Real(kind=RDbl)              :: a(3,3),y(3,3),d

    !** Invert the matrix using techniques from Numerical Recipes
    np = 3
    n = 3
    d = 1.0d0
    
    !setup an identity matrix
    Do i = 1,n
       Do j = 1,n
          y(i,j) = 0
       EndDo
       y(i,i) = 1
    EndDo
    
    !** Perform LU decomposition 
    a = m1%comp
    Call matrixops_ludcmp(a,n,np,indx,d)
    
    !** find the inverse by solving an equation for each column
    Do j = 1,n
       Call matrixops_lubksb(a,n,np,indx,y(1:n,j))
    EndDo
    
    matrix_inverse%comp = y

  End Function matrix_inverse

  
End Module matrix
