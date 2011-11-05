program diel

  implicit none
  include "mpif.h"
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  !variables for dielectric part
  INTEGER :: lenS,nbands,nvd,ncd,nvu,ncu,number_partial_occup,tag,req,tnv,tis
  INTEGER :: nprow,npcol,reorder,ndim,grid_comm,me,iam_blacs,ictxt,myrow,mycol
  INTEGER :: lrow,lcol,shortsize,mb,li,lj,rsrc,csrc,nrank,nrow,ncol,info
  INTEGER :: ns,error,i,j,lenAct,index,Ind,is,iv,ic,Ngk,npair,ok,start_ind
  INTEGER :: nv_start,lnv,nprocs,nprocsUp,nprocsDown,iam,master,ierr,nv_len,blocksize
  integer :: Phi_short_row,Phi_short_col,Phi_nblocks_row,Phi_nblocks_col,jb,sizerow,sizecol
  integer :: vecs_short_row,vecs_short_col,vecs_nblocks_row,vecs_nblocks_col,vecslrow,vecslcol
  integer :: nfound,ik,jk,LWORK,LRWORK,LIWORK
  INTEGER, DIMENSION(3) :: S
  INTEGER, DIMENSION(2) :: nv,nc,dims,ncoords
  INTEGER, dimension(9) :: desca,descz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: active,indx,spin_list,nv_list,start_list,IWORK
  integer, dimension(MPI_STATUS_SIZE) :: status
  LOGICAL, dimension(2) :: period
  integer, allocatable, dimension(:,:) :: Phi_owner_row,Phi_owner_col,Phi_block_x,Phi_block_y
  integer, allocatable, dimension(:,:) :: vecs_owner_row,vecs_owner_col,vecs_block_x,vecs_block_y
  integer, allocatable, dimension(:) :: Phi_row_list,Phi_col_list,vecs_row_list,vecs_col_list,ranksvec
  integer :: descvecs(9),descAllPhi(9),descPhiKPhi(9)
  integer :: world_comm,comm_worker,group_world,group_worker,nprocsNew,wfblocksize

  REAL(DP) :: G2max, Vcell,dV,pi,absq,matfac
  REAL(DP), DIMENSION(3) :: absqvec
  REAL(DP), DIMENSION(3,3) :: R,inverseR
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: G2,Gact2,eigs1,eigs2,de,w,K,KX,RWORK
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: G,Gact,eigs,wf1,wf2,wfq1,wfq2

  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:) :: wfv,wfc,Phi,wfvq,qRem,WORK,lqRem,tempqRem
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:) :: AllPhi,KAllPhi,wfu,wfd,vecs,lAllPhi,tempAllPhi,la,lz,na,nz
  complex(dpc), allocatable, dimension(:,:) :: lKAllPhi,nAllPhi,nKAllPhi,lvecs,lMcorr
  COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: wf,wfq
  character(len=20) :: fil_indx,fil_evc1,fil_evc2,fil_evcq1,fil_evcq2,fil_e1,fil_e2

  !variables for sigma part
  integer :: IV_new, IS_new, nfreq,ib,ind_start,lenind,len
  integer, dimension(2) :: nc_new
  integer, allocatable, dimension(:) :: index_list,nindex_list,len_list
  real(DP) :: eta, SigX,stR,corr,ws_start,ws_step,elda,VMfac,corrtmp
  real(DP), allocatable, dimension(:) :: ws,ReSig,ImSig
  real(DP), allocatable, dimension(:,:) :: eigs_new

  complex(DPC) :: Ieta
  complex(DPC), allocatable, dimension(:) :: WFs,coh,sex,wfn,PhiX,Phi_new
  complex(DPC), allocatable, dimension(:) :: PhiKPhi,V,V2,Wp2,wp_2,wp,VMvec,lVMvec,lPhiKPhi,nVMvec,nPhiKPhi,lV,mV
  complex(DPC), allocatable, dimension(:,:) :: VM,VMC,vrem,lVM,lVMC,lvrem,nvrem
  complex(DPC), allocatable, dimension(:,:,:) :: wf_new  

  master = 0
  tag    = 0
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,iam,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
  
  S(1) = 90; S(2) = 90; S(3) = 90
  lenS = S(1)*S(2)*S(3)

  R(1,1) = 13.5;  R(1,2) =  0.0;  R(1,3) =  0.0
  R(2,1) =  0.0;  R(2,2) = 13.5;  R(2,3) =  0.0
  R(3,1) =  0.0;  R(2,3) =  0.0;  R(3,3) = 13.5 

  inverseR(1,1) = 1./R(1,1); inverseR(1,2) = 0.0;       inverseR(1,3) = 0.0
  inverseR(2,1) = 0.0;       inverseR(2,2) = 1./R(2,2); inverseR(2,3) = 0.0
  inverseR(3,1) = 0.0;       inverseR(3,2) = 0.0;       inverseR(3,3) = 1./R(3,3)   

  !WARNING: only works for special case of rectangular cell
  Vcell = R(1,1)*R(2,2)*R(3,3)
  dV    = Vcell/lenS

  nprow = 12; npcol = 12
  G2max = 12.0;
  wfblocksize = 4
  blocksize = 64
  nbands = 500
  nvd =  132
  ncd =  200
  number_partial_occup = 2
  Ngk = 41755
  nvu = nvd + number_partial_occup
  ncu = ncd - number_partial_occup
  npair = nvu*ncu + nvd*ncd
  nv(1) = nvu
  nv(2) = nvd
  nc(1) = ncu
  nc(2) = ncd
  ns  =   2
  if(nvd+ncd>nbands) then 
     print *,"bad input: not enough bands"
     stop
  end if

  fil_indx = "PBEdat/indx"
  fil_evc1 = "PBEdat/evc1"
  fil_evc2 = "PBEdat/evc2"
  fil_evcq1= "PBEdat/evcq1"
  fil_evcq2= "PBEdat/evcq2"
  fil_e1   = "PBEdat/eigenval1"
  fil_e2   = "PBEdat/eigenval2"

  pi = 2.*ACOS(0.0)
  if(iam == master) then
     print *, "nvu=",nvu,"  ncu=",ncu
     print *, "nvd=",nvd,"  ncd=",ncd
     print *, "npair=",npair
  end if
  !****************************************************************
  ! ALLOCATE VARIABLES for setup
  !****************************************************************
  ALLOCATE( G(lenS,3) , STAT=error)
  ALLOCATE( G2(lenS) , STAT=error)
  call setup(inverseR,S,G,G2)

  lenAct = 0
  DO i=1,lenS
     IF( G2(i) < G2max) THEN
        lenAct = lenAct + 1
     END IF
  END DO
  if(iam==master) print *,"lenAct=", lenAct

  ALLOCATE( active(lenAct) , STAT=error)
  ALLOCATE( Gact(lenAct,3) , STAT=error)
  ALLOCATE( Gact2(lenAct) , STAT=error)
  ALLOCATE( K(lenAct) , STAT=error)
  ALLOCATE( KX(lenS) , STAT=error)
  
  index = 1
  DO i=1,lenS
     IF( G2(i) < G2max) THEN
        active(index) = i
        index = index + 1
     END IF
  END DO

  Gact  = G(active,:)
  Gact2 = SUM(Gact**2,2) 
  
  K     = 4.*pi / Gact2;
  K(1)  = 0.0;
  KX    = 4.*pi/ G2;
  KX(1) = (48.0* Vcell**2 /pi)**(1.0/3.0);

  absqvec(1) = 0.0
  absqvec(2) = 0.0
  absqvec(3) = 0.001
  absqvec = 2.*pi*MATMUL( absqvec,inverseR)
  absq = SQRT( absqvec(1)**2 + absqvec(2)**2 + absqvec(3)**2 )

  deallocate(G,G2,Gact,Gact2)
  !*******************************************************************
  !***************END SETUP*******************************************
  !**********ALLOCATE STUFF FOR DIELECTRIC MATRIX*******************
  !*****************************************************************
  ALLOCATE( wfv(lenS) , STAT=error)
  ALLOCATE( wfvq(lenS) , STAT=error)
  ALLOCATE( wfc(lenS) , STAT=error)
  ALLOCATE( Phi(lenS) , STAT=error)
  ALLOCATE( indx(Ngk) , STAT=error)
  ALLOCATE( w(npair) , STAT=error)
  ALLOCATE( qRem(npair) , STAT=error)
  ALLOCATE( de(npair) , STAT=error)
  ALLOCATE( eigs(nv(1)+nc(1),ns) , STAT=error)
  
  ALLOCATE( wf(Ngk,nv(1)+nc(1),2) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for wf"; STOP
  END IF
  
  ALLOCATE( wfq(Ngk,nv(1),2) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for wfq"; STOP
  END IF
  
  ALLOCATE( AllPhi(lenAct,npair) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for AllPhi"; STOP
  END IF

  ALLOCATE( KAllPhi(lenAct,npair) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for KAllPhi"; STOP
  END IF
  !*****************************************************
  !*********DONE ALLOCATING*****************************
  !*********START READING IN DATA***********************
  if(iam == master) then
        
     ALLOCATE( wfu(Ngk,nbands) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wfu"; STOP
     END IF
     
     ALLOCATE( wfd(Ngk,nbands) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wfd"; STOP
     END IF
     
     ALLOCATE( wf1(nbands*Ngk,2) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wf1"; STOP
     END IF

     ALLOCATE( wf2(nbands*Ngk,2) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wfq2"; STOP
     END IF

     ALLOCATE( wfq1(nbands*Ngk,2) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wfq1"; STOP
     END IF

     ALLOCATE( wfq2(nbands*Ngk,2) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wfq2"; STOP
     END IF   

     ALLOCATE( eigs1(nbands) , STAT=error)
     ALLOCATE( eigs2(nbands) , STAT=error)
     !*****************************************************
     print *,"start reading"
     OPEN(10,FILE=fil_indx,form='formatted')
     DO i=1,Ngk
        READ(10,*) indx(i)
     END DO
     CLOSE(10)
     
     OPEN(11,FILE=fil_e1,form='formatted')
     DO i=1,nbands
        READ(11,*) eigs1(i)
     END DO
     CLOSE(11)
     
     OPEN(12,FILE=fil_e2,form='formatted')
     DO i=1,nbands
        READ(12,*) eigs2(i)
     END DO
     CLOSE(12)

     print *,"reading evc1"
     OPEN(13,FILE=fil_evc1,form='formatted')
     DO i=1,nbands*Ngk
        READ(13,*) wf1(i,1),wf1(i,2)
     END DO
     CLOSE(13)
     
     print *,"reading evc2"
     OPEN(14,FILE=fil_evc2,form='formatted')
     DO i=1,nbands*Ngk
        READ(14,*) wf2(i,1),wf2(i,2)
     END DO
     CLOSE(14)
     
     print *,"reading evcq1"
     OPEN(15,FILE=fil_evcq1,form='formatted')
     DO i=1,nbands*Ngk
        READ(15,*) wfq1(i,1),wfq1(i,2)
     END DO
     CLOSE(15)

     print *,"reading evcq2"
     OPEN(16,FILE=fil_evcq2,form='formatted')
     DO i=1,nbands*Ngk
        READ(16,*) wfq2(i,1),wfq2(i,2)
     END DO
     CLOSE(16)
     
     wfu = RESHAPE( cmplx( wf1(:,1), wf1(:,2) ), (/ Ngk,nbands /) )
     wf(:,:,1) = wfu(:,1:nv(1)+nc(1)) /SQRT(Vcell)
     wfd = RESHAPE( cmplx( wf2(:,1), wf2(:,2) ), (/ Ngk,nbands /) )
     wf(:,:,2) = wfd(:,1:nv(1)+nc(1)) /SQRT(Vcell)

     wfu = RESHAPE( cmplx( wfq1(:,1), wfq1(:,2) ), (/ Ngk,nbands /) )
     wfq(:,:,1) = wfu(:,1:nv(1)) /SQRT(Vcell)
     wfd = RESHAPE( cmplx( wfq2(:,1), wfq2(:,2) ), (/ Ngk,nbands /) )
     wfq(:,:,2) = wfd(:,1:nv(1)) /SQRT(Vcell)

     eigs(:,1) = eigs1(1:nv(1)+nc(1))
     eigs(:,2) = eigs2(1:nv(1)+nc(1))

     Ind = 1
     do is = 1,ns
        do iv = 1,nv(is)
           do ic = 1,nc(is)
              de(Ind) = eigs(nv(is)+ic,is) - eigs(iv,is)
              Ind = Ind + 1
           end do
        end do
     end do
     
     deallocate(eigs2,eigs1,wfu,wfd,wf1,wf2,wfq1,wfq2)
  end if
  
  call MPI_BCAST(indx, Ngk, MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
  call MPI_BCAST(de,2*npair,MPI_REAL,master,MPI_COMM_WORLD,ierr) !is factor of 2 correct?
  call MPI_BCAST(eigs,2*ns*(nv(1)+nc(1)),MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(wf, 4*(nv(1)+nc(1))*Ngk, MPI_COMPLEX,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(wfq, 4*nv(1)*Ngk, MPI_COMPLEX,master,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  !***************************************************************************
  !***************start calculating dielectric matrix*************************
  !***************************************************************************
  nprocsUp   = ceiling( real(nprocs)/2.0)
  nprocsDown = nprocs - nprocsUp

  if( wfblocksize*nprocsUp < nv(1) .and. wfblocksize*nprocsDown < nv(2)) then
     print *,"not enough procs for wfblocksize"
     stop
  end if
  
  nprocsUp   = ceiling( real(nv(1)) / real(wfblocksize))
  nprocsDown = ceiling( real(nv(2)) / real(wfblocksize))   
  nprocsNew  = nprocsUp + nprocsDown 
  
  if(iam == master) then
     print *,"nprocsNew",nprocsNew
     print *,"nprocsup=",nprocsUp
     print *,"nprocsdown=",nprocsDown
  end if
  
  allocate(spin_list(nprocsNew))
  allocate(nv_list(nprocsNew))
  allocate(start_list(nprocsNew))

  do i = 0,nprocsNew-1
     if(i < nprocsUp) then
        is = 1
        lnv = floor( real(nv(1)) / real(nprocsUp-1) )
        nv_start  = i*lnv + 1
        nv_len    = lnv
        if(i == nprocsUp-1) nv_len = nv(1) - (nprocsUp-1)*lnv
     else
        is = 2
        lnv = floor( real(nv(2)) / real(nprocsDown-1) )
        nv_start  = (i-nprocsUp)*lnv + 1
        nv_len    = lnv
        if(i == nprocsNew-1) nv_len = nv(2) - (nprocsDown-1)*lnv
     end if
     
     spin_list(i+1) = is
     nv_list(i+1)   = nv_len
     start_list(i+1)= nv_start
  end do
  
  if(iam==master) then
     do i=1,nprocsNew
        print *,i,spin_list(i),nv_list(i),start_list(i)
     end do
  end if

  call MPI_COMM_DUP(MPI_COMM_WORLD,world_comm,ierr)
  call MPI_Comm_group(world_comm, group_world, ierr)
  allocate(ranksvec(nprocsNew))
  do i=1,nprocsNew
     ranksvec(i) = i-1
  end do
  call MPI_Group_incl(group_world, nprocsNew, ranksvec, group_worker, ierr)
  call MPI_Comm_create(world_comm, group_worker, comm_worker, ierr)

  if(iam < nprocsNew) then
     
     is       = spin_list(iam + 1)
     nv_len   =   nv_list(iam + 1)
     nv_start = start_list(iam + 1)

     ALLOCATE( lAllPhi(lenAct,nv_len*nc(is)) , STAT=error)
     ALLOCATE( lqRem(nv_len*nc(is)) , STAT=error)

     Ind = 1
     DO iv = nv_start, nv_start+nv_len-1 
     
        wfv = (0.0, 0.0)
        wfv(indx) = wf(:,iv,is)
        call cI(wfv,S)
        PRINT *,"is=",is," state=",iv," norm=",dV*SUM( ABS(wfv)**2 )
     
        wfvq = (0.0, 0.0)
        wfvq(indx) = wfq(:,iv,is)
        call cI(wfvq,S)
        !PRINT *,"norm of shifted-q wf=",dV*SUM( ABS(wfvq)**2 )
     
        DO ic=1,nc(is)
        
           wfc = (0.0, 0.0)
           wfc(indx) = wf(:,nv(is)+ic,is)
           call cI(wfc,S)        
           Phi = wfv * CONJG(wfc) 
           call cJ( Phi,S)
           lAllPhi(:,Ind) = Phi(active)
        
           lqRem(Ind) =  dV* SUM( CONJG(wfc) * wfvq ) /absq ;
           Ind = Ind + 1
        END DO
     END DO

     call MPI_Barrier(comm_worker,ierr)

     !*****************sending AllPhi to master********************************
     call MPI_Isend(lAllPhi, 2*lenAct*nv_len*nc(is), MPI_COMPLEX, master, tag, comm_worker, req, ierr)
  
     if(iam  == master) then
        do i = 0,nprocsNew-1

           tnv = nv_list(i+1)
           tis = spin_list(i+1)

           allocate(tempAllPhi(lenAct, tnv*nc(tis) ) )
           call MPI_Recv(tempAllPhi, 2*lenAct*tnv*nc(tis), MPI_COMPLEX, i, tag, comm_worker, status, ierr)
           start_ind = (tis-1)*nv(1)*nc(1) + (start_list(i+1) - 1)*nc(tis) + 1
           AllPhi(:, start_ind : start_ind + tnv*nc(tis) - 1 ) = tempAllPhi
        end do
        deallocate(tempAllPhi)
     end if
     call MPI_Barrier(comm_worker,ierr)
     !*********************sending qRem to master********************************************
     call MPI_Isend(lqRem, 2*nv_len*nc(is), MPI_COMPLEX, master, tag, comm_worker, req, ierr)
     if(iam  == master) then
        do i = 0,nprocsNew-1
           
           tnv = nv_list(i+1)
           tis = spin_list(i+1)

           allocate(tempqRem(tnv*nc(tis) ) )
           call MPI_Recv(tempqRem, 2*tnv*nc(tis), MPI_COMPLEX, i, tag, comm_worker, status, ierr)
           start_ind = (tis-1)*nv(1)*nc(1) + (start_list(i+1) - 1)*nc(tis) + 1
           qRem(start_ind : start_ind + tnv*nc(tis) - 1 ) = tempqRem
        end do
     end if
  
  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_BCAST(qRem, 2*npair, MPI_COMPLEX, master, MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  !***************compute polarizability************************************  
  ndim = 2
  dims(1) = nprow; dims(2) = npcol
  period(1) = .false.; period(2) = .false.
  reorder = 1

  call MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,period,reorder,grid_comm,ierr)
  call MPI_Comm_rank(grid_comm, me, ierr)
  
  call blacs_pinfo(iam_blacs,nprocs)
  call blacs_get(-1,0,ictxt)
  call blacs_gridinit(ictxt,'r',nprow,npcol)
  call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)

  if(iam==master) then
               
     print *,"doing dmult"
     call dmult( K, AllPhi, KAllPhi, lenAct, npair )
  end if

  !*********distribute AllPhi and KAllPhi*********************************
  !note: row_list(prow_row) contains number of rows of local array on proc
  !note: owner_row(block_row,block_list) gives row of proc where block lives
  !note: block_x(block_row,block_list) gives position of blocks in the local array

  allocate(Phi_row_list(nprow))
  allocate(Phi_col_list(nprow))
  call get_localarray_sizes(lenAct,npair,blocksize,nprow,Phi_short_row,Phi_short_col, & 
       Phi_nblocks_row,Phi_nblocks_col,Phi_row_list,Phi_col_list)
  
  allocate(Phi_owner_row( Phi_nblocks_row, Phi_nblocks_col))
  allocate(Phi_owner_col( Phi_nblocks_row, Phi_nblocks_col))
  allocate(Phi_block_x( Phi_nblocks_row, Phi_nblocks_col))
  allocate(Phi_block_y( Phi_nblocks_row,Phi_nblocks_col))
  
  call get_ownerlist(nprow,Phi_nblocks_row,Phi_nblocks_col,& 
       Phi_owner_row,Phi_owner_col,Phi_block_x,Phi_block_y)
  
  lrow = Phi_row_list(myrow+1)
  lcol = Phi_col_list(mycol+1)
  allocate(lAllPhi(lrow,lcol))
  allocate(lKAllPhi(lrow,lcol))
  call MPI_BCAST(AllPhi, 2*npair*lenAct, MPI_COMPLEX, master, MPI_COMM_WORLD,ierr)
  !***********distributing AllPhi to processor grid*************************
  do ib = 1,Phi_nblocks_row
     do jb = 1,Phi_nblocks_col

        if(Phi_owner_row(ib,jb) == myrow .and. Phi_owner_col(ib,jb) == mycol) then
           
           sizerow = blocksize
           sizecol = blocksize
           if(ib == Phi_nblocks_row) sizerow = Phi_short_row
           if(jb == Phi_nblocks_col) sizecol = Phi_short_col

           do li= 1, sizerow
              do lj=1, sizecol
                 lAllPhi(li + blocksize*Phi_block_x(ib,jb), lj + blocksize*Phi_block_y(ib,jb)) = &
                      AllPhi(li + blocksize*(ib-1) , lj + blocksize*(jb-1) )
              end do
           end do
        end if
     end do
  end do
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  deallocate(AllPhi)
  if(iam==master) print *,"done distributing AllPhi"
  !*************distributing KAllPhi to processor grid************************************
  call MPI_BCAST(KAllPhi, 2*npair*lenAct, MPI_COMPLEX, master, MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done bcast"

  do ib = 1,Phi_nblocks_row
     do jb = 1,Phi_nblocks_col

        if(Phi_owner_row(ib,jb) == myrow .and. Phi_owner_col(ib,jb) == mycol) then

           sizerow = blocksize
           sizecol = blocksize
           if(ib == Phi_nblocks_row) sizerow = Phi_short_row
           if(jb == Phi_nblocks_col) sizecol = Phi_short_col
           
           do li= 1, sizerow
              do lj=1, sizecol
                 lKAllPhi(li + blocksize*Phi_block_x(ib,jb), lj + blocksize*Phi_block_y(ib,jb)) = &
                      KAllPhi(li + blocksize*(ib-1) , lj + blocksize*(jb-1) )
              end do
           end do
        end if
     end do
  end do
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done with KAllPhi"
  
  allocate(vecs_row_list(nprow))
  allocate(vecs_col_list(nprow))
  call get_localarray_sizes(npair,npair,blocksize,nprow,vecs_short_row,vecs_short_col,&
       vecs_nblocks_row,vecs_nblocks_col,vecs_row_list,vecs_col_list)
  
  vecslrow = vecs_row_list(myrow+1)
  vecslcol = vecs_col_list(mycol+1)
  allocate(lvecs(vecslrow,vecslcol))
  allocate(lMcorr(vecslrow,vecslcol))

  allocate(vecs_owner_row( vecs_nblocks_row, vecs_nblocks_col))
  allocate(vecs_owner_col( vecs_nblocks_row, vecs_nblocks_col))
  allocate(vecs_block_x( vecs_nblocks_row, vecs_nblocks_col))
  allocate(vecs_block_y( vecs_nblocks_row,vecs_nblocks_col))
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"start ownerlist",vecs_nblocks_row,vecs_nblocks_col
  
  call get_ownerlist(nprow,vecs_nblocks_row,vecs_nblocks_col,& 
       vecs_owner_row,vecs_owner_col,vecs_block_x,vecs_block_y)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done with vecslists"
  !****************making Mcorr = conjg(qRem) * qRem*******************************
  do ib = 1,vecs_nblocks_row
     do jb = 1,vecs_nblocks_col

        if(vecs_owner_row(ib,jb) == myrow .and. vecs_owner_col(ib,jb) == mycol) then

           sizerow = blocksize
           sizecol = blocksize
           if(ib == vecs_nblocks_row) sizerow = vecs_short_row
           if(jb == vecs_nblocks_col) sizecol = vecs_short_col
           
           do li= 1, sizerow
              do lj=1, sizecol
                 lMcorr(li + blocksize*vecs_block_x(ib,jb), lj + blocksize*vecs_block_y(ib,jb)) = &
                      conjg(qRem(li + blocksize*(ib-1))) * qRem(lj + blocksize*(jb-1) )
                 
              end do
           end do
        end if
     end do
  end do
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done with lMcorr"

  lvecs = lMcorr
  matfac = 4.0*pi/Vcell
  
  rsrc=0; csrc=0
  call descinit(descAllPhi,lenAct,npair,blocksize,blocksize,rsrc,csrc,ictxt,lrow,info)
  call descinit(descvecs,npair,npair,blocksize,blocksize,rsrc,csrc,ictxt,vecslrow,info)
  !**********compute: vecs = Vcell * AllPhi' * KAllPhi + matfac*Mcorr********************************************
  call pzgemm('C','N',npair,npair,lenAct,cmplx(Vcell,0.0,kind=dpc),lKAllPhi,1,1,descAllPhi,lAllPhi,1,1,descAllPhi,&
       cmplx(matfac,0.0,kind=dpc),lvecs,1,1,descvecs)

  if(iam==master) print *,"done with AllPhi'*KAllPhi"

  !**************doing vecs = 2*sqrt(de)*vecs*sqrt(de) + de^2 **************************
  do ib = 1,vecs_nblocks_row
     do jb = 1,vecs_nblocks_col

        if(vecs_owner_row(ib,jb) == myrow .and. vecs_owner_col(ib,jb) == mycol) then

           sizerow = blocksize
           sizecol = blocksize
           if(ib == vecs_nblocks_row) sizerow = vecs_short_row
           if(jb == vecs_nblocks_col) sizecol = vecs_short_col
           
           do li= 1, sizerow
              do lj=1, sizecol
                 lvecs(li + blocksize*vecs_block_x(ib,jb), lj + blocksize*vecs_block_y(ib,jb)) = &
                      2*sqrt(de(li + blocksize*(ib-1)))* &
                      lvecs(li + blocksize*vecs_block_x(ib,jb), lj + blocksize*vecs_block_y(ib,jb)) * &
                      sqrt(de(lj + blocksize*(jb-1)))
                 
                 if( li + blocksize*(ib-1) ==  lj + blocksize*(jb-1) ) then
                    lvecs(li + blocksize*vecs_block_x(ib,jb), lj + blocksize*vecs_block_y(ib,jb)) = &
                         lvecs(li + blocksize*vecs_block_x(ib,jb), lj + blocksize*vecs_block_y(ib,jb)) +&
                         de(li + blocksize*(ib-1))**2
                 end if

              end do
           end do
        end if
        
     end do
  end do
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !*************parallel inversion**********************************
  allocate(lz(vecslrow,vecslcol))
  rsrc=0; csrc=0
  allocate(WORK(1))
  allocate(RWORK(1))
  allocate(IWORK(1))
  call pzheevd('V', 'U', npair, lvecs,1,1,descvecs, w,lz,1,1,descvecs, WORK, -1, RWORK,1,IWORK,1, ok)
  LWORK = int(real(WORK(1)))
  LRWORK = int(RWORK(1))
  LIWORK = IWORK(1)
  allocate(WORK(LWORK))
  allocate(RWORK(LRWORK))
  allocate(IWORK(LIWORK))

  if(iam==master) print *,"doing inversion of size", npair
  call pzheevd('V', 'U', npair, lvecs,1,1,descvecs, w,lz,1,1,descvecs, WORK,LWORK, RWORK,LRWORK,IWORK,LIWORK, ok)
  if(iam==master) print *,"done with inversion"
  w = SQRT(w)
  lvecs = lz
  deallocate(WORK,RWORK,IWORK,lz)
  
  !***************doing vecs = sqrt(de)*vecs/sqrt(w) ************************
  do ib = 1,vecs_nblocks_row
     do jb = 1,vecs_nblocks_col

        if(vecs_owner_row(ib,jb) == myrow .and. vecs_owner_col(ib,jb) == mycol) then
           
           sizerow = blocksize
           sizecol = blocksize
           if(ib == vecs_nblocks_row) sizerow = vecs_short_row
           if(jb == vecs_nblocks_col) sizecol = vecs_short_col
           
           do li= 1, sizerow
              do lj=1, sizecol
                 lvecs(li + blocksize*vecs_block_x(ib,jb), lj + blocksize*vecs_block_y(ib,jb)) = &
                      sqrt(de(li + blocksize*(ib-1))) *&
                      lvecs(li + blocksize*vecs_block_x(ib,jb), lj + blocksize*vecs_block_y(ib,jb)) /&
                      sqrt(w(lj + blocksize*(jb-1)))
                 
              end do
           end do
        end if
        
     end do
  end do
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam == master) print *,"done with sqrt(de)*vecs/sqrt(w)"
  !***********************************************************************
  !***************done with dielectric part ******************************
  !***************starting sigma part*************************************
  IV_new = 135
  IS_new = 2
  Ieta  = cmplx( 0.0, 0.005/27.21)
  nfreq = 600
  ws_start = 11.0 !eV
  ws_step  = 0.02 !eV
  nbands = 500
  nc_new(1) = ncu
  nc_new(2) = ncd
  fil_indx = "PBEdat/indx"
  fil_evc1 = "PBEdat/evc1"
  fil_evc2 = "PBEdat/evc2"
  fil_e1   = "PBEdat/eigenval1"
  fil_e2   = "PBEdat/eigenval2"
  
  ALLOCATE( ws(nfreq) , STAT=error)
  ALLOCATE( ReSig(nfreq) , STAT=error)
  ALLOCATE( ImSig(nfreq) , STAT=error)
  ALLOCATE( eigs_new(nv(1)+nc(1),2) , STAT=error)
  ALLOCATE( PhiX(lenS) , STAT=error)
  ALLOCATE( Phi_new(lenAct) , STAT=error)
  ALLOCATE( WFs(lenS) , STAT=error)
  ALLOCATE( wfn(lenS) , STAT=error)
  ALLOCATE( coh(nfreq) , STAT=error)
  ALLOCATE( sex(nfreq) , STAT=error)
  ALLOCATE( PhiKPhi(npair) , STAT=error)
  ALLOCATE( V(npair) , STAT=error)
  ALLOCATE( V2(npair) , STAT=error)
  ALLOCATE( Wp2(lenS) , STAT=error)
  ALLOCATE( wp(lenS) , STAT=error)
  ALLOCATE( VMvec(npair) , STAT=error)   
  ALLOCATE( wp_2(lenS) , STAT=error)
  ALLOCATE( lPhiKPhi(vecslrow))
  ALLOCATE( lV(vecslrow))
  ALLOCATE( lVM(vecslrow,vecslcol) )  
  ALLOCATE( lVMC(vecslrow,vecslcol) )

  ALLOCATE( wf_new(Ngk, nv(1)+nc(1), 2) , STAT=error)
  IF (error /= 0) THEN
     PRINT *, "could not allocate space for wf_new"; STOP
  END IF

    indx = 0
  !*********START READING IN DATA***********************
  if(iam == master) then
        
     ALLOCATE( wfu(Ngk,nbands) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wfu"; STOP
     END IF

     ALLOCATE( wfd(Ngk,nbands) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wfd"; STOP
     END IF
     
     ALLOCATE( wf1(nbands*Ngk,2) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wf1"; STOP
     END IF

     ALLOCATE( wf2(nbands*Ngk,2) , STAT=error)
     IF (error /= 0) THEN
        PRINT *, "could not allocate space for wfq2"; STOP
     END IF

     ALLOCATE( eigs1(nbands) , STAT=error)
     ALLOCATE( eigs2(nbands) , STAT=error)
     !*****************************************************
     print *,"start reading"
     OPEN(10,FILE=fil_indx,form='formatted')
     DO i=1,Ngk
        READ(10,*) indx(i)
     END DO
     CLOSE(10)
     
     print *,"reading e1"
     OPEN(11,FILE=fil_e1,form='formatted')
     DO i=1,nbands
        READ(11,*) eigs1(i)
     END DO
     CLOSE(11)
     
     OPEN(12,FILE=fil_e2,form='formatted')
     DO i=1,nbands
        READ(12,*) eigs2(i)
     END DO
     CLOSE(12)

     print *,"reading evc1"
     OPEN(13,FILE=fil_evc1,form='formatted')
     DO i=1,nbands*Ngk
        READ(13,*) wf1(i,1),wf1(i,2)
     END DO
     CLOSE(13)
     
     print *,"reading evc2"
     OPEN(14,FILE=fil_evc2,form='formatted')
     DO i=1,nbands*Ngk
        READ(14,*) wf2(i,1),wf2(i,2)
     END DO
     CLOSE(14)
     
     wfu = RESHAPE( cmplx( wf1(:,1), wf1(:,2) ), (/ Ngk,nbands /) )
     wf_new(:,:,1) = wfu(:,1:nv(1)+nc_new(1)) /SQRT(Vcell)
     wfd = RESHAPE( cmplx( wf2(:,1), wf2(:,2) ), (/ Ngk,nbands /) )
     wf_new(:,:,2) = wfd(:,1:nv(1)+nc_new(1)) /SQRT(Vcell)

     eigs_new(:,1) = eigs1(1:nv(1)+nc_new(1))
     eigs_new(:,2) = eigs2(1:nv(1)+nc_new(1))
     
     deallocate(eigs2,eigs1,wfu,wfd,wf1,wf2)
  end if
  
  call MPI_BCAST(indx, Ngk, MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
  call MPI_BCAST(eigs_new,2*ns*(nv(1)+nc_new(1)),MPI_REAL,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(wf_new, 4*(nv(1)+nc_new(1))*Ngk, MPI_COMPLEX,master,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !wf_new = wf
  deallocate(wf,wfq)
  !eigs_new = eigs
  
  !*******compute q=0 correction*********************************************
  VMfac = 4.0*pi/Vcell**2 * (48./pi*Vcell**2)**(1./3.)
  lVM = VMfac *lMcorr 
  call descinit(descPhiKPhi,npair,1,blocksize,1,rsrc,csrc,ictxt,vecslrow,info)
  !*************compute: VM = vecs' * VM * vecs *****************************
  call pzgemm('N','N',npair,npair,npair,cmplx(1.0,0.0,kind=dpc),lVM,1,1,descvecs,lvecs,1,1,descvecs,&
       cmplx(0.0,0.0,kind=dpc),lVMC,1,1,descvecs)
  call pzgemm('C','N',npair,npair,npair,cmplx(1.0,0.0,kind=dpc),lvecs,1,1,descvecs,lVMC,1,1,descvecs,&
       cmplx(0.0,0.0,kind=dpc),lVM,1,1,descvecs)

  allocate(lVMvec(vecslrow))
  allocate(index_list(vecslrow))
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam == master) print *,"done with vecs*VM*vecs"

  nfound = 1
  index_list = 0
  !*****************find diagonal elements of VM************************
  do ib = 1,vecs_nblocks_row
     do jb = 1,vecs_nblocks_col

        if(vecs_owner_row(ib,jb) == myrow .and. vecs_owner_col(ib,jb) == mycol) then
           sizerow = blocksize
           sizecol = blocksize
           if(ib == vecs_nblocks_row) sizerow = vecs_short_row
           if(jb == vecs_nblocks_col) sizecol = vecs_short_col

           do li= 1, sizerow
              do lj=1, sizecol

                 if( li + blocksize*(ib-1) == lj + blocksize*(jb-1)) then
                    index_list( nfound ) = li + blocksize*(ib-1)
                    lVMvec(nfound ) = &
                         lVM(li + blocksize*vecs_block_x(ib,jb), lj + blocksize*vecs_block_y(ib,jb))
                    nfound = nfound + 1
                 end if
              end do
           end do
        end if
     end do
  end do
           
  call MPI_Isend(lVMvec, 2*vecslrow, MPI_COMPLEX, master, tag, grid_comm, req, ierr)
  call MPI_Isend(index_list, vecslrow, MPI_INTEGER, master, tag, grid_comm, req, ierr)
  !***********collect data and put it into VMvec*********************************
  if(iam == master) then

     do i = 0,nprow-1
        do j = 0,nprow-1
           
           nrow = vecs_row_list(i+1)
           allocate(nVMvec(nrow))
           allocate(nindex_list(nrow))
           ncoords(1) = i; ncoords(2) = j
           call MPI_Cart_rank(grid_comm,ncoords,nrank,ierr)
           call MPI_Recv(nVMvec,2*nrow, MPI_COMPLEX, nrank, tag, grid_comm, status, ierr)
           call MPI_Recv(nindex_list,nrow, MPI_INTEGER, nrank, tag, grid_comm, status, ierr)

           li = 1
           do while( nindex_list(li) .ne. 0 .and. li <= nrow)
              VMvec( nindex_list(li)) = nVMvec(li)
              li = li + 1
           end do
        end do
     end do
  end if
  
  deallocate(lVM,lVMC)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam == master) print *,"done with VMvec"
  call MPI_BCAST(VMvec,2*npair,MPI_COMPLEX,master,MPI_COMM_WORLD,ierr)
  !**************setup stuff for self-energy*************************
  ws(1) = ws_start
  DO i =2,nfreq
     ws(i) = ws(i-1) + ws_step 
  END DO
  ws = ws/27.21

  WFs = (0.0, 0.0)
  WFs(indx) = wf_new(:,IV_new,IS_new) 
  call cI( WFs, S);
  
  coh = (0.0, 0.0)
  sex = (0.0, 0.0)
  SigX = 0.0
  stR  = 0.0   
  
  DO ib = 1, nv(IS_new)+nc_new(IS_new)

     lPhiKPhi = cmplx(0.0, 0.0,kind=dpc)

     if(iam == master) then
        wfn = (0.0, 0.0)
        wfn(indx) = wf_new(:,ib,IS_new)
        call cI(wfn, S);
        print *,"norm of state ",ib,"=",dV*sum(abs(wfn)**2)
        
        PhiX = WFs * conjg( wfn)
        call cJ(PhiX, S)
        Phi_new  = PhiX(active)
        
        PhiKPhi = Vcell * matmul( conjg(Phi_new) , KAllPhi)
  
        !**********distribute PhiKPhi*********************
        do i=0, nprow-1
           
           nrow = vecs_row_list(i+1)
           allocate(nPhiKPhi(nrow))

           do ik = 1,vecs_nblocks_row
              if(vecs_owner_row(ik,1) == i ) then
                 
                 sizerow = blocksize
                 if(ik == vecs_nblocks_row) sizerow = vecs_short_row
                 
                 do li= 1, sizerow
                    nPhiKPhi(li + blocksize*vecs_block_x(ik,1) ) = PhiKPhi(li + blocksize*(ik-1) )
                 end do
              end if 
           end do

           ncoords(1) = i; ncoords(2) = 0
           call MPI_Cart_rank(grid_comm,ncoords,nrank,ierr)
           call MPI_Isend(nPhiKPhi, 2*nrow, MPI_COMPLEX, nrank, tag, grid_comm, req, ierr)
        end do
     end if
        
     if(mycol == 0) then
        call MPI_Recv(lPhiKPhi,2*vecslrow, MPI_COMPLEX, master, tag, grid_comm, &
             status, ierr)
     end if
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     !***************compute V = PhiKPhi * vecs by doing: V^T = vecs^T * PhiKPhi^T ***************
     call pzgemv('T',npair,npair,cmplx(1.0,0.0,kind=dpc),lvecs,1,1,descvecs,lPhiKPhi,1,1,descPhiKPhi,1,& 
          cmplx(0.0,0.0,kind=dpc),lV,1,1,descPhiKPhi,1)
     
     if(mycol == 0) then
        call MPI_Isend(lV, 2*vecslrow, MPI_COMPLEX, master, tag, grid_comm, req, ierr)
     end if
     
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     if(iam == master) then
        !*****************collect data*****************************************
        do i = 0,nprow-1
           nrow = vecs_row_list(i+1)
           allocate(mV(nrow))
           ncoords(1) = i; ncoords(2) = 0
           call MPI_Cart_rank(grid_comm,ncoords,nrank,ierr)
           call MPI_Recv(mV,2*nrow, MPI_COMPLEX, nrank, tag, grid_comm, status, ierr)

           do ik = 1,vecs_nblocks_row
              if(vecs_owner_row(ik,1) == i ) then

                 sizerow = blocksize
                 if(ik == vecs_nblocks_row) sizerow = vecs_short_row

                 do li= 1, sizerow
                    V(li + blocksize*(ik-1)) = mV( li + blocksize*vecs_block_x(ik,1))
                 end do
              end if

           end do
        end do
     
        V2 = abs(V)**2
        
        if(ib == IV_new) then
           DO i=1,npair
              V2(i) = V2(i) + VMvec(i)
           END DO
        end if
        
        DO Ind = 1,npair
           
           if(ib <= nv(IS_new)) then
              sex = sex- V2(Ind) /( ws - eigs_new(ib,IS_new) - w(Ind) +Ieta )
              sex = sex+ V2(Ind) /( ws - eigs_new(ib,IS_new) + w(Ind) -Ieta )
              coh = coh+ V2(Ind) /( ws - eigs_new(ib,IS_new) - w(Ind) +Ieta )
              stR = stR+ V2(Ind) / w(Ind)
           else
              coh = coh+ V2(Ind) /( ws - eigs_new(ib,IS_new) - w(Ind) +Ieta)
              stR = stR+ V2(Ind) / w(Ind)
           end if
        END DO
        
        if(ib <= nv(IS_new)) then
           SigX = SigX -Vcell* real( sum( PhiX * KX * conjg( PhiX )))
        end if
        
     end if
  END DO
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam == master) print *,"computing static remainder correction"
  !**************compute: vrem = KAllPhi * vecs  *************
  lrow = Phi_row_list(myrow+1)
  lcol = Phi_col_list(mycol+1)
  allocate( lvrem( lrow, lcol) , STAT=error)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam == master) print *,"start pzgemm"

  call pzgemm('N','N',lenAct,npair,npair,cmplx(1.0,0.0,kind=dpc),lKAllPhi,1,1,descAllPhi,lvecs,1,1,descvecs,&
       cmplx(0.0,0.0,kind=dpc),lvrem,1,1,descAllPhi)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam == master) print *,"done pzgemm"

  call MPI_Isend(lvrem, 2*lrow*lcol, MPI_COMPLEX, master, tag, grid_comm, req, ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  deallocate(KAllPhi,lvecs,lKAllPhi,wf_new,lAllPhi)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done deallocating"
  
  allocate(vrem(lenAct,npair))
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done allocating vrem"
  !************master collects data for vrem****************************************
  if(iam == master) then
     do i = 0,nprow-1
        do j =0,nprow-1
           
           nrow = Phi_row_list(i+1)
           ncol = Phi_col_list(j+1)
           allocate(nvrem(nrow,ncol))
           ncoords(1) = i; ncoords(2) = j
           call MPI_Cart_rank(grid_comm,ncoords,nrank,ierr)
           call MPI_Recv(nvrem,2*nrow*ncol, MPI_COMPLEX, nrank, tag, grid_comm, status, ierr)

           do ik = 1,Phi_nblocks_row
              do jk = 1,Phi_nblocks_col
                 if(Phi_owner_row(ik,jk) == i .and. Phi_owner_col(ik,jk) == j ) then
              
                    sizerow = blocksize
                    sizecol = blocksize
                    if(ik == Phi_nblocks_row) sizerow = Phi_short_row
                    if(jk == Phi_nblocks_col) sizecol = Phi_short_col
                    
                    do li= 1, sizerow
                       do lj = 1, sizecol
                          vrem(li + blocksize*(ik-1), lj + blocksize*(jk-1)) = &
                               nvrem( li + blocksize*Phi_block_x(ik,jk), lj + blocksize*Phi_block_y(ik,jk))
                       end do
                    end do
                 end if
              end do
           end do
        end do
     end do
     print *,"done making vrem"
  end if
  
  call MPI_BCAST(vrem,2*npair*lenAct,MPI_COMPLEX,master,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done bcast"

  allocate(len_list(nprocs))
  allocate(start_list(nprocs))
  do i = 0,nprocs-1
     len = floor( real(npair) / real(nprocs-1) )
     ind_start = i*len + 1
     lenind    = len
     if(i == nprocs-1) then
        lenind = npair - (nprocs-1)*len
     endif
     
     len_list(i+1)  = lenind
     start_list(i+1)= ind_start
  end do

  lenind   =   len_list(iam + 1)
  ind_start= start_list(iam + 1)

  Wp2  = (0.0, 0.0)
   
  do Ind = ind_start, ind_start+lenind-1
     wp = (0.0, 0.0)
     wp(active) = vrem(:,Ind)
     call cI(wp, S)
     Wp2  = Wp2 + (abs(wp)**2 + VMvec(Ind) )/w(Ind) 
  end do
     
  corr = -dV * sum( abs(WFs)**2 * Wp2)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done with corr"
  call MPI_REDUCE(corr,corrtmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(iam==master) print *,"done with Reduce",corrtmp,corr

  if(iam == master) then
     
     corr = corrtmp
     stR  = ( stR + corr) /2.0
     elda = eigs_new(IV_new,IS_new)
     
     print *,"corr=",corr
     print *,"state ",IV_new," with spin ",IS_new
     print *,"G2max =", G2max
     print *,"nvu=",nvu,"  ncu=",ncu
     print *,"nvd=",nvd,"  ncd=",ncd
     print *,"npair=",npair
     print *,"stR  in eV: ",stR*27.21
     print *,"SigX in eV: ",SigX*27.21
     print *,"Ieta in eV: ",Ieta*27.21
     print *,"elda in eV: ",elda*27.21
     
     ReSig = real(coh+sex+SigX+stR)
     ImSig = aimag(coh+sex+SigX+stR)
     OPEN(13,FILE='sig_out',form='formatted')
     DO i=1,nfreq
        WRITE(13,*) ws(i)*27.21, ReSig(i)*27.21, ImSig(i)*27.21
     END DO
     CLOSE(13)
  end if

     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_Finalize(ierr)
  print *,"done"
end program diel
