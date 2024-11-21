!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_save_soln
!! NAME
!!
!!   amr_1blk_save_soln
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_save_soln()
!!
!! ARGUMENTS
!!
!!   None
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!
!! CALLS
!!
!!   No other paramesh routines called.
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!!   This routine saves a global solution update into the time 
!!   synchronized global solution arrays, as is required when 
!!   using NO_PERMANENT_GUARDCELLS and the amr_1blk_guardcell routines.
!!
!! AUTHORS
!!
!! Written :     Peter MacNeice          May 1999
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_1blk_save_soln

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree

      Implicit None

      Integer :: lb

!-----Begin Executable Code

      If (no_permanent_guardcells) Then

!$OMP PARALLEL DO
        Do lb = 1, maxblocks

        If (nvar > 0) Then
          unk(:,:,:,:,lb) = gt_unk(:,:,:,:,lb)
        End If
        If (nbndvar > 0) Then
          facevarx(:,:,:,:,lb) = gt_facevarx(:,:,:,:,lb)
          If (ndim >= 2) Then
            facevary(:,:,:,:,lb) = gt_facevary(:,:,:,:,lb)
          End If
          If (ndim == 3) Then
            facevarz(:,:,:,:,lb) = gt_facevarz(:,:,:,:,lb)
          End If
        End If ! End If (nbndvar > 0)

        If (ndim > 1) Then
        If (nvaredge > 0) Then
          unk_e_x(:,:,:,:,lb) = gt_unk_e_x(:,:,:,:,lb)
          unk_e_y(:,:,:,:,lb) = gt_unk_e_y(:,:,:,:,lb)
          If (ndim == 3) Then
          unk_e_z(:,:,:,:,lb) = gt_unk_e_z(:,:,:,:,lb)
          End If
        End If  ! End If (nvaredge > 0)
        End If  ! End If (ndim > 1)

        If (nvarcorn > 0) Then
          unk_n(:,:,:,:,lb) = gt_unk_n(:,:,:,:,lb)
        End If

      End Do ! End Do lb
!$OMP END PARALLEL DO

      End If  ! End If (no_permanent_guardcells)

      Return
      End Subroutine amr_1blk_save_soln
