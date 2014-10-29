!requires lapack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Modules !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module inputData

        integer deltaStep

endmodule inputData

module openmp

        integer np !number of threads

endmodule openmp

module atomData

        integer nAtoms
        real, allocatable :: atomPos(:,:,:)
        real, allocatable :: atomCharges(:)
        real, allocatable :: atomMasses(:)

        integer nRes                               ! number of residues in psf file
        integer, allocatable :: atomResNumber(:)   ! residue numbers directly from psf file
        integer, allocatable :: resNumber(:)       ! residue numbers renumbered to be sequential
        integer, allocatable :: resAtomStart(:)    ! first atom number in each residue.  To be used for boundary atoms

        character*80 atomPsfFile
        character*80 atomDcdFile

        integer nSteps

endmodule atomData


module cgData

        integer nCg
        real, allocatable :: cgPos(:,:,:)
        real, allocatable :: cgCharges(:)
        character*80 cgDcdFile

endmodule cgData

module output

        character*80 outputFile

endmodule output

module timing

        real (kind=8) ti,t1,t2
        real (kind=8) readTime,computeTime,accumTime, fitTime

endmodule timing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Main Program !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cg_charge_fit
        use atomData
        use cgData
        use output
        use timing
        use openmp
        implicit none
        integer, parameter :: maxIter=10
        real, allocatable :: A(:,:)
        real, allocatable :: B(:,:)
        real, allocatable :: C(:,:)
        integer, allocatable :: boundaryRes(:)
        integer, allocatable :: minBoundaryRes(:)
        real minChi2
        real chi2
        real (kind=8) omp_get_wtime
        integer iter
!        integer, allocatable :: gen_rand_ordered_seq(:)

        ti = omp_get_wtime()

        call parse_command_line(atomPsfFile,atomDcdFile,cgDcdFile,outputFile,nCg)

        call read_psf_file

        allocate(A(nCg,nCg),B(nCg,nAtoms),C(nAtoms,nAtoms),cgCharges(nCg))
        allocate(boundaryRes(nCg-1),minBoundaryRes(nCg-1))!,gen_rand_ordered_seq(nCg-1))

        call read_atom_trajectory(C)

        allocate(cgPos(nCg,3,nSteps))

        do iter=1,maxIter

                !generate boundary atoms
                call gen_rand_ordered_seq(nCg-1,1,nRes,boundaryRes)

                !make CG trajectory
                call create_CG_traj(atomPos,atomMasses,nAtoms,resAtomStart,nRes,cgPos,nCg,nSteps,boundaryRes)

                !Accumulate A and B
                call compute_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B,nSteps)

                !Fit charges
                call fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCg,chi2)

        enddo

        t2 = omp_get_wtime()
        write(*,'("Total time elapsed:",f8.3)') t2-ti

endprogram cg_charge_fit



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine create_CG_traj(atomPos,atomMasses,nAtoms,resAtomStart,nRes,cgPos,nCg,nSteps,boundaryRes)
        implicit none
        integer nAtoms
        integer nSteps
        integer nCg
        integer nRes
        real atomPos(nAtoms,3,nSteps)
        real atomMasses(nAtoms)
        integer resAtomStart(nRes)
        real cgPos(nCg,3,nSteps)
        integer boundaryRes(nCg-1)
        integer startAtom, stopAtom
        integer cg, atom, step
        real temp(3)
        real mass


        do step=1,nSteps
                startAtom=1
                do cg=1,nCg
                
                        temp=0
                        mass=0
                        if (cg<nCg) then
                                stopAtom = resAtomStart(boundaryRes(cg))-1
                        else
                                stopAtom = nAtoms
                        endif
                        do atom=startAtom,stopAtom
                                temp = atomMasses(atom)*atomPos(atom,:,step)
                                mass = mass+ atomMasses(atom)
                        enddo
                        cgPos(cg,:,step) = temp/mass
                        startAtom = resAtomStart(boundaryRes(cg))
                enddo
        enddo


endsubroutine create_CG_traj


subroutine gen_rand_ordered_seq(numSeq,minSeq,maxSeq,seq)
        implicit none
        integer numSeq
        integer minSeq
        integer maxSeq
        integer seq(numSeq)
        real temp
        integer tempInt
        integer i, j, k
        logical diff

        do i=1,numSeq
                diff = .false.
                do while (diff.eqv..false.) 
                        call random_number(temp)
                        tempInt = int(temp*(maxSeq-minSeq))+minSeq
                        diff = .true.
                        if (i.gt.1) then
                                do j=1,i-1
                                        if (tempInt==seq(j)) then
                                                diff = .false.
                                                exit 
                                        endif
                                enddo
                        endif
                enddo
                seq(i) = tempInt
                if (i.gt.1) then
                        do j=1,i-1
                                if (seq(i)<seq(j)) then
                                        tempInt=seq(i)
                                        do k=i,j+1  !downwards
                                                seq(k) = seq(k-1)
                                        enddo
                                        seq(j) = tempInt
                                        exit
                                endif
                        enddo
                endif
        enddo
       
endsubroutine gen_rand_ordered_seq

subroutine read_atom_trajectory(C)
        use inputData
        use cgData
        use atomData
        use openmp
        use timing
        implicit none
        integer step
        integer atom1, atom2, j
        real C(nAtoms,nAtoms)
        real dist, temp

        C=0

        ! Read atom dcd header and grab nAtoms, nSteps
        call read_dcd_header(atomDcdFile,nAtoms,nSteps,20)

        allocate(atomPos(nAtoms,3,nSteps))

        do step=1,nSteps
                call read_dcd_step(atomPos(:,:,step),nAtoms,20)
                !Atom-atom distance for later
                !$omp parallel private(atom1,atom2,j,dist,temp) SHARED(atomPos,C) num_threads(np)
                !$omp do schedule(dynamic,10)
                do atom1 = 1, nAtoms-1
                        do atom2 = atom1+1,nAtoms
                                dist = 0
                                do j=1,3
                                        temp = atomPos(atom1,j,step)-atomPos(atom2,j,step)
                                        dist = dist + temp*temp
                                enddo
                                dist = sqrt(dist)
                                C(atom1,atom2) = C(atom1,atom2)-dist
                                C(atom2,atom1) = C(atom2,atom1)-dist
                        enddo
                enddo
                !$omp end do
                !$omp end parallel 

        enddo

        close(20)

endsubroutine read_atom_trajectory


!read command line information
subroutine parse_command_line(atomPsfFile,atomDcdFile,cgDcdFile,outputFile,nCg)
        use inputData
        use openmp
        implicit none
        integer i
        character*30 arg
        character*80 atomPsfFile
        character*80 atomDcdFile
        character*80 cgDcdFile
        character*80 outputFile
        logical deltaStepFlag
        logical atomPsfFlag
        logical atomDcdFlag
        logical cgDcdFlag
        logical nCgFlag
        logical npFlag
        logical outputFlag
        integer nCg

        atomPsfFlag = .false.
        atomDcdFlag = .false.
        cgDcdFlag = .false.
        nCgFlag = .false.
        npFlag = .false.
        deltaStepFlag = .false.
        outputFlag=.false.
        i=1
        do 
   
                call get_command_argument(i, arg) 
   
                select case (arg) 
   
                case ('-ncg')
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') nCg
                        nCgFlag=.true.
                case ('-np')
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') np
                        npFlag=.true.
                case ('-psf')
                        i = i+1
                        call get_command_argument(i,atomPsfFile)
                        atomPsfFlag=.true.
                case ('-o')
                        i = i+1
                        call get_command_argument(i,outputFile)
                        outputFlag=.true.
                case ('-adcd')
                        i = i+1
                        call get_command_argument(i,atomDcdFile)
                        atomDcdFlag=.true.
                case ('-cgdcd')
                        i = i+1
                        call get_command_argument(i,cgDcdFile)
                        cgDcdFlag=.true.
                case ('-step') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') deltaStep
                        deltaStepFlag = .true.
                case default 
                        print '(a,a,/)', 'Unrecognized command-line option: ', arg 
                        print*, 'Usage: fit_charges.x -ncg [# cg sites] -psf [atom psf file] -adcd [atom dcd file] -cgdcd [cg dcd file] -step [delta step size] -o [outputfile name] -np [# of threads]'
                        stop 
                end select 
                i = i+1
                if (i.ge.command_argument_count()) exit
        enddo

        if (nCgFlag.eqv..false.) then
                write(*,'("Must provide number of CG sites using command line argument -ncg [number of CG sites]")')
                stop
        endif
        if (atomPsfFlag.eqv..false.) then
                write(*,'("Must provide a psf file using command line argument -psf [psf file name]")')
                stop
        endif
        if (atomDcdFlag.eqv..false.) then
                write(*,'("Must provide a atom dcd file using command line argument -adcd [atom dcd file name]")')
                stop
        endif
        if (cgDcdFlag.eqv..false.) then
                write(*,'("Must provide a CG dcd file using command line argument -cgdcd [CG dcd file name]")')
                stop
        endif
        if (outputFlag.eqv..false.) then
                write(*,'("No output file name provided.  Please provide an output file name with the following command line argument: -o [output file name]")')
                stop
        endif
        if (deltaStepFlag.eqv..false.) then
                write(*,'("Using default step size of 1.  Change this with command line argument -step [delta step size]")')
                deltaStep=1
        endif
        if (npFlag.eqv..false.) then
                write(*,'("Running on 1 thread.  Change this with command line argument -np [number of threads]")')
                np = 1
        else
                write(*,'("Running on",i3," thread(s)")') np
        endif

endsubroutine parse_command_line


!Read atomic charges from some file
subroutine read_psf_file
        use atomData
        implicit none
        integer atom, j
        character*6 check           !character to check if NATOM is in the line
        character*8 numChar         !character to read number of atoms.  must be converted to an integer
        character*4 atomCheck
        character*24 posChar
        integer ios
        integer resCount

        !open the psf file
        open(10,file=atomPsfFile)

        !run through the header of the file and look for the number of atoms
        do 

                read(10,'(a8,2x,a6)') numChar, check

                !if we read the number of atoms exit this do loop
                if (check.eq.'NATOM ') then
                        !this converts the character natoms_char to the integer natoms
                        read(numChar,*) nAtoms
                        write(*,*) "Number of atoms=", nAtoms
                        !Now that we know the number of atoms we must allocate the arrays
                        allocate(atomCharges(nAtoms),atomResNumber(nAtoms),resNumber(nAtoms),atomMasses(nAtoms))
                        !Now we loop through the number of atoms and read the pertinent information
                        resCount=1
                        do atom=1,nAtoms

                                read(10,'(14x,i4,16x,f10.6,4x,f10.4)') atomResNumber(atom),atomCharges(atom),atomMasses(atom)
                                ! make residue numbers sequential
                                if (atomResNumber(atom) .ne. atomResNumber(atom-1)) then
                                        resCount = resCount + 1
                                endif
                                resNumber(atom) = resCount

                        enddo
                elseif (check.eq.'NBOND:') then
                        exit
                endif

        enddo

        close(10)

        nRes = resCount
        !Create array of residue start atom numbers
        allocate(resAtomStart(nRes))
        resCount=1
        do atom=1,nAtoms

                if (atom==1 .or. resNumber(atom) .ne. resNumber(atom-1)) then
                        resAtomStart(resCount) = atom
                        resCount = resCount +1
                endif

        enddo

        write(*,*) "Total Charge of the system = ", sum(atomCharges)

endsubroutine read_psf_file


!requires lapack
subroutine compute_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B,nSteps)
        use openmp
        implicit none
        integer nAtoms
        integer nCg
        integer nSteps
        real atomPos(nAtoms,3,nSteps)
        real cgPos(nCg,3,nSteps)
        real A(nCg,nCg)
        real B(nCg,nAtoms)
        real dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1, atom2
        !open MP stuff
        integer omp_get_thread_num, omp_get_num_procs, omp_get_num_threads, omp_get_max_threads
        integer nProcs, nThreads, threadID
        integer step

        A=0
        B=0

        !Compute the distance between the CG sites
        !$omp parallel private(step,cgSite1,cgSite2,atom1,j,dist,temp) SHARED(cgPos,A,B,atomPos) num_threads(np)
        !$omp do schedule(dynamic)
        do step=1,nSteps
                do cgSite1 = 1, nCg-1
                        do cgSite2 = cgSite1+1,nCg
                                dist = 0
                                do j=1,3
                                        temp = cgPos(cgSite1,j,step)-cgPos(cgSite2,j,step)
                                        dist = dist + temp*temp
                                enddo
                                dist = sqrt(dist)
                                !$omp atomic
                                A(cgSite1,cgSite2) = A(cgSite1,cgSite2)-dist
                                !$omp atomic
                                A(cgSite2,cgSite1) = A(cgSite2,cgSite1)-dist
                             !   print*, "Distance between site", cgSite1, "and site", cgSite2,":",dist
                        enddo
                enddo

                !Compute the distance between cg sites and atoms
                do cgSite1 = 1, nCg
                        do atom1 = 1,nAtoms
                                dist = 0
                                do j=1,3
                                        temp = cgPos(cgSite1,j,step)-atomPos(atom1,j,step)
                                        dist = dist + temp*temp
                                enddo
                                !$omp atomic
                                B(cgSite1,atom1) = B(cgSite1,atom1)-sqrt(dist)
                        enddo
                enddo
        enddo
        !$omp end do 
        !$omp end parallel 

endsubroutine compute_A_B_matrices


subroutine fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCg, rss)
        implicit none
        integer nAtoms
        integer nCg
        real A(nCg,nCg)
        real ATemp(nCg+1,nCg)
        real D(nCg-1,nCg)
        real B(nCg,nAtoms)
        real C(nAtoms,nAtoms)
        real BTemp(nCg+1,nAtoms)
        real cgCharges(nCg,1)
        real atomCharges(nAtoms)
        real atomChargesM(nAtoms,1)
        real newB(nCg+1)
        real temp(1,1)
        real rss
        !lapack routine variables
        real (kind=8) workQuery(1)
        real (kind=8), allocatable :: work(:)
        integer info
        integer lwork
        integer j, i, k

        !First we need to modify A and B to have the correct matrix properties 
        !create D matrices
        D=0
        do i=1,nCg-1
                D(i,i)=1
                D(i,i+1)=-1
        enddo 
        !multiply A by D0 and B by D1 giving new matrices A and B the correct behavior
        ATemp(1:nCg,:) = matmul(D,A)
        BTemp(1:nCg,:) = matmul(D,B)
        !generate new matrices with last line having 1s meaing charge is conserved
        ATemp(nCg+1,:) = 1.0
        BTemp(nCg+1,:) = 1.0


        !Now use lapack routine to solve least squares problem A*cgCharges=B*atomCharges
        newB = matmul(BTemp,atomCharges)
        call dgels("N",nCg+1,nCg,1,ATemp,nCg+1,newB,nCg+1,workQuery,-1,info)
        lwork = int(workQuery(1))
        allocate(work(lwork))
        call dgels("N",nCg+1,nCg,1,ATemp,nCg+1,newB,nCg+1,work,lwork,info)
        deallocate(work)
        
        cgCharges(:,1) = newB(1:nCg)
        atomChargesM(:,1) = atomCharges

        temp = matmul(transpose(atomChargesM),matmul(C,atomChargesM))+matmul(transpose(cgCharges),matmul(A,cgCharges))-2*matmul(transpose(cgCharges),matmul(B,atomChargesM))
        rss = temp(1,1)

!        write(*,'(a17,f8.3)') "Total CG Charge:",sum(cgCharges)
!        write(*,'("Residual Sum of Squares:",f8.3,f8.3)') newB(nCg+1), rss(1,1)


endsubroutine fit_charges


