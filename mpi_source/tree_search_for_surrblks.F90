      subroutine tree_search_for_surrblks ()

      use local_tree_common
      use physicaldata
      use tree
      use paramesh_dimensions
      use paramesh_interfaces
      use paramesh_comm_data
      use mpi_morton, only : lperiodicx, lperiodicy, lperiodicz
      use constants

      implicit none

      include 'mpif.h'

      real :: neigh_coord(3), neigh_coord2(3)
      integer :: ierr, nprocs, mype
      integer :: neigh_lb,neigh_proc, neigh_nodetype
      integer :: i,j,k,ii,jj,kk,iboun,lb,ir,ib
      integer :: neigh_list(2,mfaces,maxblocks),surr_blks_temp(2,3,3,maxblocks)
      integer :: nrecv, reqr(maxblocks), statr(MPI_STATUS_SIZE,maxblocks)
      integer :: temp(3)

      logical :: found

      real :: time_exe, time_max
      real :: eps,accuracy

      accuracy = 100./10.**precision(accuracy)
      if (accuracy > 1.0e-10) then
         eps = 1.e-10
      else
         eps = accuracy
      end if

      time_exe = MPI_WTIME()

      call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

! FIND SURROUNDING BLOCKS
!$OMP PARALLEL DO &
!$OMP PRIVATE (kk, k, jj, j, ii, i, neigh_coord, neigh_coord2, neigh_lb, &
!$OMP          neigh_proc, neigh_nodetype, found, iboun)
      do lb = 1, lnblocks

         kk = -k3d
         do k = 1,1+2*k3d

         neigh_coord(3) = coord(3,lb) + kk*bsize(3,lb)
         if (lperiodicz.and.neigh_coord(3).lt.grid_zmin)               &
          neigh_coord(3) = neigh_coord(3) + (grid_zmax-grid_zmin)
         if (lperiodicz.and.neigh_coord(3).gt.grid_zmax)               &
          neigh_coord(3) = neigh_coord(3) - (grid_zmax-grid_zmin)

         jj = -k2d
         do j = 1,1+2*k2d

         neigh_coord(2) = coord(2,lb) + jj*bsize(2,lb)
         if (lperiodicy.and.neigh_coord(2).lt.grid_ymin)               &
          neigh_coord(2) = neigh_coord(2) + (grid_ymax-grid_ymin)
         if (lperiodicy.and.neigh_coord(2).gt.grid_ymax)               &
          neigh_coord(2) = neigh_coord(2) - (grid_ymax-grid_ymin)

         ii = -1
         do i = 1,3

         neigh_coord(1) = coord(1,lb) + ii*bsize(1,lb)
         if (lperiodicx.and.neigh_coord(1).lt.grid_xmin)               &
          neigh_coord(1) = neigh_coord(1) + (grid_xmax-grid_xmin)
         if (lperiodicx.and.neigh_coord(1).gt.grid_xmax)               &
          neigh_coord(1) = neigh_coord(1) - (grid_xmax-grid_xmin)

         neigh_coord2(:) = neigh_coord(:)

!--------Reset coordinates of neighbor in spherical coordinates
         If (spherical_pm) Then
            If ((((jj == -1).and.(abs(bnd_box(1,2,lb)) < eps)) .or.    &  
                 ((jj ==  1).and.(abs(bnd_box(2,2,lb)-pi) < eps)))     & 
                 .and. lsingular_line ) Then
               neigh_coord2(2) = coord(2,lb)
               If (neigh_coord2(3) < pi) Then
                  neigh_coord2(3) = neigh_coord2(3) + pi
               Elseif(neigh_coord2(3) > pi)  Then
                  neigh_coord2(3) = neigh_coord2(3) - pi
               End If
            End If
         End If

         neigh_lb = -1
         neigh_proc = -1
         neigh_nodetype = -1
         found = .false.

!--------First check boundaries

         do iboun = 1, nboundaries
            if (ndim == 1) then
            if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
                neigh_coord2(1) < boundary_box(2,1,iboun)) then
               found = .true.
               neigh_lb = boundary_index(iboun)
               neigh_proc = boundary_index(iboun)
            end if
            elseif (ndim == 2) then
            if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
                neigh_coord2(1) < boundary_box(2,1,iboun) .and.        &
                neigh_coord2(2) > boundary_box(1,1+k2d,iboun) .and.    &
                neigh_coord2(2) < boundary_box(2,1+k2d,iboun)) then
               found = .true.
               neigh_lb = boundary_index(iboun)
               neigh_proc = boundary_index(iboun)
            end if
            elseif (ndim == 3) then
            if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
                neigh_coord2(1) < boundary_box(2,1,iboun) .and.        &
                neigh_coord2(2) > boundary_box(1,1+k2d,iboun) .and.    &
                neigh_coord2(2) < boundary_box(2,1+k2d,iboun) .and.    &
                neigh_coord2(3) > boundary_box(1,1+2*k3d,iboun) .and.  &
                neigh_coord2(3) < boundary_box(2,1+2*k3d,iboun)) then
               found = .true.
               neigh_lb = boundary_index(iboun)
               neigh_proc = boundary_index(iboun)
            end if
            end if
         end do

         surr_blks(:,i,j,k,lb) = -1
         if (found) then
          surr_blks(1,i,j,k,lb) = neigh_lb
          surr_blks(2,i,j,k,lb) = neigh_proc
          surr_blks(3,i,j,k,lb) = neigh_nodetype
!         else
!          call search_sub_tree(local_tree,neigh_coord2,lrefine(lb),     &
!                               neigh_lb,neigh_proc,neigh_nodetype,found)
!          if (nodetype(lb) == 1) then
!           surr_blks(1,i,j,k,lb) = neigh_lb
!           surr_blks(2,i,j,k,lb) = neigh_proc
!           surr_blks(3,i,j,k,lb) = neigh_nodetype
!          else
! TESTING
!           surr_blks(1,i,j,k,lb) = -100000
!           surr_blks(2,i,j,k,lb) = -100000
!           surr_blks(3,i,j,k,lb) = neigh_nodetype
!          end if
         end if

         ii = ii + 1
      end do
      jj = jj + k2d
      end do
      kk = kk + k3d
      end do

    end do
!$OMP END PARALLEL DO

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (mype == 0) print *,' STARTING search_around !!!!!!!!!!!!!'

      call search_around()

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (mype == 0) print *,' DONE setting surr_blks !!!!!!!!!!!!!'

      end subroutine tree_search_for_surrblks

      subroutine search_around()

      use physicaldata
      use tree
      use paramesh_dimensions
      use paramesh_interfaces
      use paramesh_comm_data
      use mpi_morton, only : lperiodicx, lperiodicy, lperiodicz

      include 'mpif.h'

         integer :: proc, tag, ierr, nprocs, mype, ireq, bid
         integer :: bid_from, bid_sent, proc_from, proc_sent
         logical :: flag
         integer :: nsend, send_req, send_status(MPI_STATUS_SIZE)
         integer :: nrecv, statr(MPI_STATUS_SIZE,maxblocks_tr), reqr(maxblocks_tr)
         integer :: npot(maxblocks_tr), npot_temp(maxblocks_tr), npotr(6,maxblocks_tr)
 !        real :: block_data(8,maxblocks) ! contains coord and size of blocks
         real :: block_datar(10,1000,maxblocks) ! contains coord and size of blocks
         real :: block_datas(10,1000,maxblocks) ! contains coord and size of blocks
!         integer :: tesurr_blks(3,3,1+2*k2d,1+2*k3d,maxblocks)
         real :: dx, dy, dz, pot_neighs(3)
         
!      do lb = 1, lnblocks
!        do k = 1, 1+k3d*2
!        do j = 1, 1+k2d*2
!        do i = 1, 3
!         if (surr_blks(1,i,j,k,lb) <= -1) then
!         tesurr_blks(:,i,j,k,lb) = surr_blks(:,i,j,k,lb)
!         else
!         tesurr_blks(:,i,j,k,lb) = -1
!         end if
!        end do
!        end do
!        end do
!      end do

         call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
         call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

         npot(:) = 0

         do lb = 1, lnblocks
            block_datas(1:3,1,lb) = coord(1:3,lb)
            block_datas(4:6,1,lb) = bsize(1:3,lb)
            block_datas(7,1,lb) = lb
            block_datas(8,1,lb) = mype
            block_datas(9,1,lb) = nodetype(lb)
            block_datas(10,1,lb) = lrefine(lb)
            surr_blks(1,2,k2d+1,k3d+1,lb) = lb
            surr_blks(2,2,k2d+1,k3d+1,lb) = mype
            surr_blks(3,2,k2d+1,k3d+1,lb) = nodetype(lb)
!            tesurr_blks(1,2,k2d+1,k3d+1,lb) = lb
!            tesurr_blks(2,2,k2d+1,k3d+1,lb) = mype
!            tesurr_blks(3,2,k2d+1,k3d+1,lb) = nodetype(lb)
            npot(lb) = 1
         end do

         do lb = 1, lnblocks
            block_datar(:,1:npot(lb),lb) = block_datas(:,1:npot(lb),lb)
            npot_temp(lb) = npot(lb)
         end do

! send from children to parent
         nrecv = 0
         do lb = 1, lnblocks
           do ichild = 1, nchild
              if (child(1,ichild,lb) > 0) then
              npot(lb) = npot(lb) + 1
              if (child(2,ichild,lb) .ne. mype) then               
                 nrecv = nrecv + 1 ! per lb
                 call MPI_IRECV(block_datar(1,npot(lb),lb),  &
                                10,                          &
                                amr_mpi_real,                &
                                child(2,ichild,lb),          &
                                child(1,ichild,lb),          &
                                MPI_COMM_WORLD,              &
                                reqr(nrecv),                 &
                                ierr)
              else
                 block_datar(:,npot(lb),lb) = block_datas(:,1,child(1,ichild,lb))
              end if
              end if
           end do
         end do

         do lb = 1, lnblocks
              ! send block_datas(1,nsend,lb) to parent
            if (parent(1,lb) > 0) then
            if (parent(2,lb) .ne. mype) then
              CALL MPI_SSEND(block_datas(1,1,lb),   &
                       10,                          &
                       amr_mpi_real,                &
                       parent(2,lb),                & 
                       lb,                          & 
                       MPI_COMM_WORLD,              &
                       ierr)
            end if
            end if
         end do

         if (nrecv > 0) then
             Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         print *,' HERE 1 '

! send to neighbors in x,y,z dirs
         do nstep = 1, ndim

         if (nstep == 1) then
            is = 1
            ie = 2
         else if (nstep == 2) then
            is = 3
            ie = 4
         else if (nstep == 3) then
            is = 5
            ie = 6
         end if

         do lb = 1, lnblocks
            block_datas(:,1:npot(lb),lb) = block_datar(:,1:npot(lb),lb)
            npot_temp(lb) = npot(lb)
         end do

         nrecv = 0
         do lb = 1, lnblocks
           npotr(:,lb) = 0
           do ineigh = is, ie
              if (neigh(1,ineigh,lb) > 0) then
              if (neigh(2,ineigh,lb) .ne. mype) then
                 nrecv = nrecv + 1 ! per lb
                 call MPI_IRECV(npotr(ineigh,lb),            &
                                1,                           &  ! need to know size of message !
                                MPI_INTEGER,                 &
                                neigh(2,ineigh,lb),          &
                                neigh(1,ineigh,lb),          &
                                MPI_COMM_WORLD,              &
                                reqr(nrecv),                 &
                                ierr)
              else
                 npotr(ineigh,lb) = npot(neigh(1,ineigh,lb))
              end if
              end if
           end do
         end do

         do lb = 1, lnblocks
            do ineigh = is, ie
            if (neigh(1,ineigh,lb) > 0) then
            if (neigh(2,ineigh,lb) .ne. mype) then
              CALL MPI_SSEND(npot(lb),               &
                       1,                           &
                       MPI_INTEGER,                 &
                       neigh(2,ineigh,lb),          & 
                       lb,                          & 
                       MPI_COMM_WORLD,              &
                       ierr)
            end if
            end if
            end do
         end do

         if (nrecv > 0) then
             Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         nrecv = 0
         do lb = 1, lnblocks
           do ineigh = is, ie
              if (neigh(1,ineigh,lb) > 0) then
              if (neigh(2,ineigh,lb) .ne. mype) then
                 nrecv = nrecv + 1 ! per lb
                 call MPI_IRECV(block_datar(1,npot(lb)+1,lb),  &
                                10*npotr(ineigh,lb),         &  ! need to know size of message !
                                amr_mpi_real,                &
                                neigh(2,ineigh,lb),          &
                                neigh(1,ineigh,lb),          &
                                MPI_COMM_WORLD,              &
                                reqr(nrecv),                 &
                                ierr)
               else
                  block_datar(:,npot(lb)+1:npot(lb)+npotr(ineigh,lb),lb) =  &
                    block_datas(:,1:npotr(ineigh,lb),neigh(1,ineigh,lb))
               end if
               npot(lb) = npot(lb) + npotr(ineigh,lb)
              end if
           end do
         end do

         do lb = 1, lnblocks
            do ineigh = is, ie
            if (neigh(1,ineigh,lb) > 0) then
            if (neigh(2,ineigh,lb) .ne. mype) then
              CALL MPI_SSEND(block_datas(1,1,lb),          &
                       10*npot_temp(lb),                   &
                       amr_mpi_real,                       &
                       neigh(2,ineigh,lb),                 & 
                       lb,                                 & 
                       MPI_COMM_WORLD,                     &
                       ierr)
            end if
            end if
            end do
         end do

         if (nrecv > 0) then
             Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         end do

! NOW SEND FROM PARENTS TO CHILDREN

         do lb = 1, lnblocks
            block_datas(:,1:npot(lb),lb) = block_datar(:,1:npot(lb),lb)
            npot_temp(lb) = npot(lb)
         end do

         nrecv = 0
         do lb = 1, lnblocks
              if (parent(1,lb) > 0) then
              npotr(:,lb) = 0
              if (parent(2,lb) .ne. mype) then
              nrecv = nrecv + 1 ! per lb
              call MPI_IRECV(npotr(1,lb),            &
                             1,                      &  
                             MPI_INTEGER,            &
                             parent(2,lb),           &
                             parent(1,lb),           &
                             MPI_COMM_WORLD,         &
                             reqr(nrecv),            &
                             ierr)
              else
              npotr(1,lb) = npot(parent(1,lb))
              end if
              end if
         end do

         do lb = 1, lnblocks
            do ichild = 1, nchild
            if (child(1,ichild,lb) > 0) then
            if (child(2,ichild,lb) .ne. mype) then
              CALL MPI_SSEND(npot(lb),              &
                             1,                     &
                             MPI_INTEGER,           &
                             child(2,ichild,lb),    & 
                             lb,                    & 
                             MPI_COMM_WORLD,        &
                             ierr)
            end if
            end if
            end do
         end do

         if (nrecv > 0) then
             Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         nrecv = 0
         do lb = 1, lnblocks
            if (parent(1,lb) > 0) then
            if (parent(2,lb) .ne. mype) then
               nrecv = nrecv + 1 ! per lb
               call MPI_IRECV(block_datar(1,npot(lb)+1,lb),  &
                              10*npotr(1,lb),                &  ! need to know size of message !
                              amr_mpi_real,                  &
                              parent(2,lb),                  &
                              parent(1,lb),                  &
                              MPI_COMM_WORLD,                &
                              reqr(nrecv),                   &
                              ierr)
             else
             block_datar(:,npot(lb)+1:npot(lb)+npotr(1,lb),lb) =  &
                block_datas(:,1:npotr(1,lb),parent(1,lb))
             end if
             npot(lb) = npot(lb) + npotr(1,lb)

            end if
         end do

         do lb = 1, lnblocks
            do ichild = 1, nchild
            if (child(1,ichild,lb) > 0) then
            if (child(2,ichild,lb) .ne. mype) then
              CALL MPI_SSEND(block_datas(1,1,lb),          &
                       10*npot_temp(lb),                   &
                       amr_mpi_real,                       &
                       child(2,ichild,lb),                 & 
                       lb,                                 & 
                       MPI_COMM_WORLD,                     &
                       ierr)
            end if
            end if
            end do
         end do

         if (nrecv > 0) then
             Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! PERIODIC

         do lb = 1, lnblocks

!            do ipot = 1, npot(lb) ! loop over potential neighbors of blocks
!	    if (lb == 66 .and. mype == 0) then
!	    print *,' pot neigh = ',block_datar(7:9,ipot,lb),ipot
!	    end if
!	    end do

         do k = 1, k3d*2+1
         do j = 1, k2d*2+1
         do i = 1,3

            pot_neighs(1) = coord(1,lb) + (i-2)*bsize(1,lb)
            if (lperiodicx) then
            if (i == 1 .and. pot_neighs(1) < grid_xmin) then
                pot_neighs(1) = coord(1,lb) + (grid_xmax-grid_xmin)
                pot_neighs(1) = pot_neighs(1) - bsize(1,lb)
            end if
            if (i == 3 .and. pot_neighs(1) > grid_xmax) then
                pot_neighs(1) = coord(1,lb) - (grid_xmax-grid_xmin)
                pot_neighs(1) = pot_neighs(1) + bsize(1,lb)
            end if
            end if

            if (ndim >= 2) then
            pot_neighs(2) = coord(2,lb) + (j-2)*bsize(2,lb)
            if (lperiodicy) then
            if (j == 1 .and. pot_neighs(2) < grid_ymin) then
                pot_neighs(2) = coord(2,lb) + (grid_ymax-grid_ymin)
                pot_neighs(2) = pot_neighs(2) - bsize(2,lb)
            end if
            if (j == 3 .and. pot_neighs(2) > grid_ymax) then
                pot_neighs(2) = coord(2,lb) - (grid_ymax-grid_ymin)
                pot_neighs(2) = pot_neighs(2) + bsize(2,lb)
            end if
            end if
            end if

            if (ndim == 3) then
            pot_neighs(3) = coord(3,lb) + (k-2)*bsize(3,lb)
            if (lperiodicz) then
            if (k == 1 .and. pot_neighs(3) < grid_zmin) then
                pot_neighs(3) = coord(3,lb) + (grid_zmax-grid_zmin)
                pot_neighs(3) = pot_neighs(3) - bsize(3,lb)
            end if
            if (k == 3 .and. pot_neighs(3) > grid_zmax) then
                pot_neighs(3) = coord(3,lb) - (grid_zmax-grid_zmin)
                pot_neighs(3) = pot_neighs(3) + bsize(3,lb)
            end if
            end if
            end if
            do ipot = 1, npot(lb) ! loop over potential neighbors of blocks

               if (bsize(1,lb) == block_datar(4,ipot,lb)) then

                   ! compute distance to block = dx
                   dx = block_datar(1,ipot,lb) - pot_neighs(1)
                   dy = 0.
                   dz = 0.
                   if (ndim >= 2) then
                     dy = block_datar(2,ipot,lb) - pot_neighs(2)
                     if (ndim == 3) then
                       dz = block_datar(3,ipot,lb) - pot_neighs(3)
                     end if
                   end if
                   ! if dx <= size of block, then its a neighbor and goes into
                   ! surr_blks for "this" block = lb
                   ! if dx > 0, ix = 3, dx == 0, ix = 2, dx < 0, ix = 1
                   if (abs(dx) < 0.5*bsize(1,lb) .and. abs(dy) < 0.5*bsize(2,lb)  &
                       .and. abs(dz) < 0.5*bsize(3,lb)) then
                     surr_blks(1,i,j,k,lb) = int(block_datar(7,ipot,lb))
                     surr_blks(2,i,j,k,lb) = int(block_datar(8,ipot,lb))
                     surr_blks(3,i,j,k,lb) = int(block_datar(9,ipot,lb))

!                     tesurr_blks(1,i,j,k,lb) = int(block_datar(7,ipot,lb))
!                     tesurr_blks(2,i,j,k,lb) = int(block_datar(8,ipot,lb))
!                     tesurr_blks(3,i,j,k,lb) = int(block_datar(9,ipot,lb))

                     exit

                   end if
               end if
            end do

         end do
         end do
         end do
       
       end do

#ifdef NOT_NOW
         do lb = 1, lnblocks
         do ipot= 1, npot(lb) ! loop over potential neighbors of blocks

               if (bsize(1,lb) == block_datar(4,ipot,lb)) then
                   ! compute distance to block = dx
                   iy = 1
                   iz = 1
                   dx = block_datar(1,ipot,lb) - coord(1,lb)
                   ix = int(dx/bsize(1,lb)) + 2
                   dy = 0.
                   dz = 0.
                   if (ndim >= 2) then
                     dy = block_datar(2,ipot,lb) - coord(2,lb)
                     iy = int(dy/bsize(2,lb)) + 2
                     if (ndim == 3) then
                       dz = block_datar(3,ipot,lb) - coord(3,lb)
                       iz = int(dz/bsize(3,lb)) + 2
                     end if
                   end if
                   ! if dx <= size of block, then its a neighbor and goes into
                   ! surr_blks for "this" block = lb
                   ! if dx > 0, ix = 3, dx == 0, ix = 2, dx < 0, ix = 1
                   if (abs(dx) <= bsize(1,lb) .and. abs(dy) <= bsize(2,lb)  &
                       .and. abs(dz) <= bsize(3,lb)) then
                     surr_blks(1,ix,iy,iz,lb) = int(block_datar(7,ipot,lb))
                     surr_blks(2,ix,iy,iz,lb) = int(block_datar(8,ipot,lb))
                     surr_blks(3,ix,iy,iz,lb) = int(block_datar(9,ipot,lb))
!                     tesurr_blks(1,ix,iy,iz,lb) = int(block_datar(7,ipot,lb))
!                     tesurr_blks(2,ix,iy,iz,lb) = int(block_datar(8,ipot,lb))
!                     tesurr_blks(3,ix,iy,iz,lb) = int(block_datar(9,ipot,lb))
                   end if
               end if
         end do
      end do
#endif

! TESTING FOR CORRECTNESS: REMOVE LATER
!      print *,' '
!      do lb = 1, lnblocks
!        do k = 1, 1+k3d*2
!        do j = 1, 1+k2d*2
!        do i = 1, 3
!!!          surr_blks(:,i,j,k,lb) = tesurr_blks(:,i,j,k,lb)
!         if (any(surr_blks(:,i,j,k,lb) .ne. tesurr_blks(:,i,j,k,lb))) then
!          print *,' ERROR: in tesurr_blk at ',i,j,k,lb,mype
!          print *,surr_blks(1:3,i,j,k,lb)
!          print *,tesurr_blks(1:3,i,j,k,lb)
!          print *,' '

!          do kk = 1, 1+k3d*2
!          do jj = 1, 1+k2d*2
!          do ii = 1, 3
!             print *,' i,j,k = ',ii,jj,kk
!             print *,' surr_blks = ',surr_blks(1:3,ii,jj,kk,lb)
!          end do
!          end do
!          end do
     
!          do ii = 1, 6
!             print *,' i = ',ii
!             print *,' neigh = ',neigh(1:2,ii,lb)
!          end do
     
!          stop
!         end if
!        end do
!        end do
!        end do
!      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      end subroutine search_around
