# why you are still using fortran?

```
!!! display matrix
subroutine print_double_matrix(BBB)
double precision,dimension(:,:)::BBB
integer::ii
integer,dimension(:),allocatable::nn
nn=shape(BBB)
do ii=1,nn(1)
write(*,"(20f12.6)")BBB(ii,:)
end do
end subroutine

!!! inverse matrix AA^(-1)
!!! AA: matrix. nn: dimension of matrix
function inv(AA,nn) result(AAinv)
integer,intent(in)::nn
double precision,intent(in) :: AA(nn,nn)
double precision            :: AAinv(nn,nn)
double precision            :: work(nn)            ! work array for LAPACK
integer         :: info,ipiv(nn)     ! pivot indices
AAinv = AA
call DGETRF(nn,nn,AAinv,nn,ipiv,info)
call DGETRI(nn,AAinv,nn,ipiv,work,nn,info)
end function inv

!!! solve linear system AA^(-1)*bb
!!! AA: matrix. BB: vector. nn: dimension of matrix
function solve(AA,bb,nn) result(bbcopy)
integer,intent(in)::nn
double precision,intent(in)::AA(nn,nn),bb(nn)
double precision::AAcopy(nn,nn),bbcopy(nn)
integer::iipp(nn),info
AAcopy=AA
bbcopy=bb
call dgesv(nn,1,AAcopy,nn,iipp,bbcopy,nn,info)
end function solve

!!! Calculate the square-root matrix: AA^(1/2)
!!! AA: matrix. nn: dimension of matrix
function sqrtmat(AA,nn) result(sqrtmatAA)
integer,intent(in)::nn
double precision,intent(in)::AA(nn,nn)
double precision::eigen_vectors(nn,nn),sqrtmatAA(nn,nn),eigen_values(nn),work(3*nn),diag_matrix(nn)
integer::info,lwork,ii,jj,kk
lwork=3*nn
eigen_vectors=AA
call dsyev('V','U',nn,eigen_vectors,nn,eigen_values,work,lwork,info)
diag_matrix=sqrt(eigen_values)
sqrtmatAA=0d0
do ii=1,nn
 do jj=1,nn
  do kk=1,nn
   sqrtmatAA(ii,jj)=sqrtmatAA(ii,jj)+eigen_vectors(ii,kk)*diag_matrix(kk)*eigen_vectors(jj,kk)
  end do
 end do
end do
end function sqrtmat

!!! Calculate the inverse of square-root matrix: AA^(-1/2)
!!! AA: matrix. nn: dimension of matrix
function invsqrtmat(AA,nn) result(invsqrtmatAA)
integer,intent(in)::nn
double precision,intent(in)::AA(nn,nn)
double precision::eigen_vectors(nn,nn),invsqrtmatAA(nn,nn),eigen_values(nn),work(3*nn),diag_matrix(nn)
integer::info,lwork,ii,jj,kk
lwork=3*nn
eigen_vectors=AA
call dsyev('V','U',nn,eigen_vectors,nn,eigen_values,work,lwork,info)
diag_matrix=1d0/sqrt(eigen_values)
invsqrtmatAA=0d0
do ii=1,nn
 do jj=1,nn
  do kk=1,nn
   invsqrtmatAA(ii,jj)=invsqrtmatAA(ii,jj)+eigen_vectors(ii,kk)*diag_matrix(kk)*eigen_vectors(jj,kk)
  end do
 end do
end do
end function invsqrtmat

```
