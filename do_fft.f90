subroutine do_FFT(fftbox,Nfft,sign)

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  include 'fftw_f77.i'

  integer, dimension(3),INTENT(IN) :: Nfft
  integer, INTENT(IN) :: sign
  complex(DPC), dimension(Nfft(1),Nfft(2),Nfft(3)),INTENT(INOUT) :: fftbox
  integer*8, save :: plus_plan, minus_plan
    
  call fftwnd_f77_create_plan(plus_plan,3,Nfft,FFTW_BACKWARD, &
       FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
  call fftwnd_f77_create_plan(minus_plan,3,Nfft,FFTW_FORWARD, &
       FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)

   if (sign == 1) then
    call fftwnd_f77_one(plus_plan,fftbox,0)
  else if (sign == -1) then
    call fftwnd_f77_one(minus_plan,fftbox,0)
  else
    print *,"sign is not 1 or -1 in do_FFT"
    stop
  endif

  return
end subroutine do_FFT

!**************************************************************
subroutine cI(fft_in,Nfft)

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  integer, dimension(3),INTENT(IN) :: Nfft
  complex(DPC), dimension(Nfft(1)*Nfft(2)*Nfft(3)),INTENT(INOUT) :: fft_in
  complex(DPC), dimension(Nfft(1),Nfft(2),Nfft(3)) :: fftbox
  integer :: sign

  fftbox = RESHAPE( fft_in, (/ Nfft(1), Nfft(2), Nfft(3) /))
  sign = 1
  call do_FFT(fftbox,Nfft,sign)
  fft_in = RESHAPE( fftbox, (/ Nfft(1)*Nfft(2)*Nfft(3) /))
  
  return
end subroutine cI

!*******************************************************************
subroutine cJ(fft_in,Nfft)

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  integer, dimension(3),INTENT(IN) :: Nfft
  complex(DPC), dimension(Nfft(1)*Nfft(2)*Nfft(3)),INTENT(INOUT) :: fft_in
  complex(DPC), dimension(Nfft(1),Nfft(2),Nfft(3)) :: fftbox
  integer :: sign,fftsize
  
  fftsize = Nfft(1)*Nfft(2)*Nfft(3)
  fftbox = RESHAPE( fft_in, (/ Nfft(1), Nfft(2), Nfft(3) /))
  sign = -1
  call do_FFT(fftbox,Nfft,sign)
  fft_in = RESHAPE( fftbox, (/ Nfft(1)*Nfft(2)*Nfft(3) /))
  fft_in = fft_in / REAL(fftsize)

  return
end subroutine cJ

!*********************************************************************
subroutine dmult( vec, mat_in, mat_out, nline, ncol)
  
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  
  INTEGER :: nline,ncol,i
  REAL(DP), dimension(nline), INTENT(IN) :: vec
  COMPLEX(DPC), dimension(nline,ncol), INTENT(IN) :: mat_in
  COMPLEX(DPC), dimension(nline,ncol), INTENT(OUT) :: mat_out
  
  DO i=1,ncol
     mat_out(:,i) = vec * mat_in(:,i)
  END DO

  return
end subroutine dmult
!********************************************************************
subroutine setup(inverseR,S,G,G2)

  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  real(DP), dimension(3,3), intent(in) :: inverseR
  integer, dimension(3), intent(in) :: S
  integer :: i,lenS,error
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ms,m1,m2,m3,n1,n2,n3
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: M,N
  real(DP), dimension(S(1)*S(2)*S(3),3), intent(out) :: G
  real(DP), dimension(S(1)*S(2)*S(3)), intent(out) :: G2
  real(DP) :: pi

  lenS = S(1)*S(2)*S(3)
  pi = 2.*ACOS(0.0)
  !**************allocate stuff***************************
  ALLOCATE( ms(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for ms"
     STOP
  END IF

  ALLOCATE( m3(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for m3"
     STOP
  END IF

  ALLOCATE( m2(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for m2"
     STOP
  END IF

  ALLOCATE( m1(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for m1"
     STOP
  END IF

  ALLOCATE( n1(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for n1"
     STOP
  END IF

  ALLOCATE( n2(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for n2"
     STOP
  END IF

  ALLOCATE( n3(lenS) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for n3"
     STOP
  END IF

  ALLOCATE( M(lenS,3) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for M"
     STOP
  END IF

  ALLOCATE( N(lenS,3) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for N"
     STOP
  END IF
  !**************done allocating******************

  DO i=1,lenS
     ms(i) = i-1
  END DO

  m3 = MOD(ms,S(3))
  m2 = MOD(FLOOR( REAL(ms)/S(3)),S(2))
  m1 = MOD(FLOOR( REAL(ms)/(S(3)*S(2))),S(1))

  M(:,1) = m1; M(:,2) = m2; M(:,3) = m3

  DO i=1,lenS

     IF( m1(i) > REAL(S(1))/2 ) THEN
        n1(i) = 1
     ELSE
        n1(i) = 0
     END IF

     IF( m2(i) > REAL(S(2))/2 ) THEN
        n2(i) = 1
     ELSE
        n2(i) = 0
     END IF

     IF( m3(i) > REAL(S(3))/2 ) THEN
        n3(i) = 1
     ELSE
        n3(i) = 0
     END IF

  END DO

  n1 = m1-n1*S(1)
  n2 = m2-n2*S(2)
  n3 = m3-n3*S(3)
  N(:,1) = n1; N(:,2) = n2; N(:,3) = n3

  G  = 2.*pi*MATMUL(N,inverseR)
  G2 = SUM(G**2,2)

  return
end subroutine setup
!*************************************************
subroutine get_localarray_sizes(matrow,matcol,blocksize,nprow,short_row,short_col,nblocks_row,nblocks_col,row_list,col_list)

  integer, intent(in) :: matrow,matcol,blocksize,nprow
  integer, intent(out):: short_row,short_col,nblocks_row,nblocks_col
  integer, dimension(nprow),intent(out):: row_list,col_list
  integer :: i,j

  short_row = mod(matrow,blocksize)
  short_col = mod(matcol,blocksize)
  nblocks_row  = (matrow-short_row)/blocksize
  nblocks_col  = (matcol-short_col)/blocksize

  if( short_row .ne. 0) then
     nblocks_row = nblocks_row + 1
  else
     short_row = blocksize
  end if

  if( short_col .ne. 0) then
     nblocks_col = nblocks_col + 1
  else
     short_col = blocksize
  end if

  !**********find row_list***************************
  row_list = 0
  do i = 1,nprow

     j=i
     do while(j <= nblocks_row)

        if( j == nblocks_row) then
           row_list(i) = row_list(i) + short_row
        else
           row_list(i) = row_list(i) + blocksize
        end if
        j = j + nprow
     end do
  end do

  !**********find col_list***************************
  col_list = 0
  do i = 1,nprow

     j=i
     do while(j <= nblocks_col)

        if( j == nblocks_col) then
           col_list(i) = col_list(i) + short_col
        else
           col_list(i) = col_list(i) + blocksize
        end if
        j = j + nprow
     end do
  end do


  return
end subroutine get_localarray_sizes
!******************************************************
subroutine get_ownerlist(nprow,nblocks_row,nblocks_col,owner_row,owner_col,block_x,block_y)

  integer, intent(in) :: nprow,nblocks_row,nblocks_col
  integer, dimension(nblocks_row,nblocks_col), intent(out) :: owner_row,owner_col,block_x,block_y  
  integer :: rowind,colind,blockrow,blockcol,i,j

  rowind = 0
  blockrow = 0
  blockcol = 0
  do i = 1,nblocks_row

     colind   = 0
     blockcol = 0
     do j = 1,nblocks_col

        owner_row(i,j) = rowind
        owner_col(i,j) = colind
        block_x(i,j)   = blockrow
        block_y(i,j)   = blockcol

        colind = colind + 1
        if(colind == nprow) then
           colind   = 0
           blockcol = blockcol + 1
        end if
     end do

     rowind = rowind + 1
     if(rowind == nprow) then
        rowind   = 0
        blockrow = blockrow + 1
     end if
  end do

  return
end subroutine get_ownerlist

