
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
        real (kind=8), allocatable :: atomPos(:,:)
        real (kind=8), allocatable :: atomCharges(:)
        character*80 atomPsfFile
        character*80 atomDcdFile

endmodule atomData


module cgData

        integer nCg
        real (kind=8), allocatable :: cgPos(:,:)
        real (kind=8), allocatable :: cgCharges(:)
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
        real (kind=8), allocatable :: AtA(:,:)
        real (kind=8), allocatable :: AtB(:,:)
        real (kind=8) omp_get_wtime
        
        ti = omp_get_wtime()
        readTime=0
        computeTime=0
        accumTime=0
        fitTime=0

        call parse_command_line(atomPsfFile,atomDcdFile,cgDcdFile,outputFile,nCg)

        call read_psf_file

        allocate(AtA(nCg,nCg),AtB(nCg,nAtoms),cgCharges(nCg))

        call read_trajectories(AtA,AtB)

        print*, "fitting charges"
        t1 = omp_get_wtime()
        call fit_charges(AtA, AtB, atomCharges, cgCharges, nAtoms, nCg,outputFile)
        t2 = omp_get_wtime()
        fitTime = t2-t1

        !Correct times for multiple cores
!        accumTime = accumTime/real(np)
!        computeTime = computeTime/real(np)

        write(*,'("Total time elapsed:",f8.3,f8.3)') readTime+accumTime+computeTime+fitTime,t2-ti
        write(*,'("Time to read dcd  :",f8.3)') readTime
        write(*,'("Time to accumulate:",f8.3)') accumTime
        write(*,'("Time to compute   :",f8.3)') computeTime
        write(*,'("Time to fit       :",f8.3)') fitTime

endprogram cg_charge_fit



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_trajectories(AtA,AtB)
        use inputData
        use cgData
        use atomData
        use openmp
        use timing
        implicit none
        integer atom1, cg1, cg2
        integer k
        integer nSteps
        integer step
        integer omp_get_thread_num, omp_get_num_procs, omp_get_num_threads, omp_get_max_threads
        real (kind=8) omp_get_wtime
        integer nProcs, nThreads, threadID
        real (kind=8) temp
        real (kind=8) A(nCg,nCg)
        real (kind=8) AtA(nCg,nCg)
        real (kind=8) B(nCg,nAtoms)
        real (kind=8) AtB(nCg,nAtoms)

        !Zero the arrays we will be accumulating over the trajectory
        do cg1=1,nCg
                do cg2=1,nCg
                        AtA(cg1,cg2)=0
                        A(cg1,cg2) = 0
                enddo
                do atom1=1,nAtoms
                        AtB(cg1,atom1) = 0
                enddo
        enddo

        !blah
        call read_dcd_header(atomDcdFile,nAtoms,nSteps,20)
        call read_dcd_header(cgDcdFile,nCg,nSteps,30)

        write(*,'("Starting step loop for ",i10," number of steps")') nSteps
        allocate(atomPos(nAtoms,3),cgPos(nCg,3))

        do step=1,nSteps
                t1 = omp_get_wtime()
                call read_dcd_step(atomPos,nAtoms,20)
                call read_dcd_step(cgPos,nCg,30)
                t2 = omp_get_wtime()
                readTime = readTime+(t2-t1)

                if (mod(step,deltaStep)==0) then
                        write(*,'("Working on step ",i10," of ",i10)') step, nSteps
                        !compute A and B for this step (distances)
                        t1 = omp_get_wtime()
                        call compute_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B)
                        t2 = omp_get_wtime()
                        computeTime = computeTime + (t2-t1)

                        call omp_set_num_threads(np)
                        t1 = omp_get_wtime()
                        !$OMP  PARALLEL SHARED(AtA,AtB,A,B,nCg,nAtoms) PRIVATE(cg1,cg2,k,atom1)
                        threadID = omp_get_thread_num()
                        !accumulate AtA and AtB
                        !$OMP DO
                        do cg1=1,nCg
                                !AtA
                                do cg2=1,nCg
                                        do k=1,nCg
                                                AtA(cg1,cg2) = AtA(cg1,cg2) + A(k,cg1)*A(k,cg2)
                                        enddo 
                                enddo
                                !AtB 
                                do atom1=1,nAtoms
                                        do k=1,nCg
                                                AtB(cg1,atom1) = AtB(cg1,atom1) + A(k,cg1)*B(k,atom1)
                                        enddo
                                enddo
                                                
                        enddo
                        !$OMP END DO NOWAIT
                        !$OMP END PARALLEL
                        t2 = omp_get_wtime()
                        accumTime = accumTime + (t2-t1)
                endif
        enddo
        !

        close(20)
        close(30)

endsubroutine read_trajectories

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
                        allocate(atomCharges(nAtoms))
                        !Now we loop through the number of atoms and read the pertinent information
                        do atom=1,nAtoms

                                read(10,'(34x,f10.6)') atomCharges(atom)
                                !add to total charge

                        enddo
                elseif (check.eq.'NBOND:') then
                        exit
                endif

        enddo

        close(10)

        write(*,*) "Total Charge of the system = ", sum(atomCharges)

endsubroutine read_psf_file


!requires lapack
subroutine compute_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B)
        use openmp
        implicit none
        integer nAtoms
        integer nCg
        real (kind=8) atomPos(nAtoms,3)
        real (kind=8) cgPos(nCg,3)
        real (kind=8) A(nCg,nCg)
        real (kind=8) ATemp(nCg,nCg)
        real (kind=8) AcolAvg(nCg)
        real (kind=8) ArowAvg(nCg)
        real (kind=8) Aavg
        real (kind=8) B(nCg,nAtoms)
        real (kind=8) BcolAvg(nAtoms)
        real (kind=8) BrowAvg(nCg)
        real (kind=8) Bavg
        real (kind=8) dist
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1, atom2
        !open MP stuff
        integer omp_get_thread_num, omp_get_num_procs, omp_get_num_threads, omp_get_max_threads
        integer nProcs, nThreads, threadID

        AcolAvg=0
        ArowAvg=0
        Aavg=0

        !Compute the distance between the CG sites
!        !$omp parallel private(cgSite1,cgSite2,j,dist) SHARED(cgPos,A,AcolAvg,ArowAvg,Aavg) num_threads(np)
!        !$omp do
        do cgSite1 = 1, nCg-1
                do cgSite2 = cgSite1+1,nCg
                        dist = 0
                        do j=1,3
                                dist = dist + (cgPos(cgSite1,j)-cgPos(cgSite2,j))**2
                        enddo
                        A(cgSite1,cgSite2) = -sqrt(dist)
                        A(cgSite2,cgSite1) = A(cgSite1,cgSite2)
                        AcolAvg(cgSite1) = AcolAvg(cgSite1) + A(cgSite2,cgSite1)
                        AcolAvg(cgSite2) = AcolAvg(cgSite2) + A(cgSite1,cgSite2)
                        ArowAvg(cgSite1) = ArowAvg(cgSite1) + A(cgSite1,cgSite2)
                        ArowAvg(cgSite2) = ArowAvg(cgSite2) + A(cgSite2,cgSite1)
!                        !$omp critical
                        Aavg = Aavg + 2*A(cgSite1,cgSite2)
!                        !$omp end critical
                enddo
        enddo
!        !$omp end do nowait
!        !$omp end parallel
        !finish averages
        Aavg = Aavg/dble(nCg**2)
        AcolAvg = AcolAvg/dble(nCg)
        ArowAvg = ArowAvg/dble(nCg)
        !complete the matrix
        do cgSite1 = 1, nCg
                A(cgSite1,cgSite1) = 0
        enddo


       !Subtract means and such 
        !$omp parallel private(cgSite1,cgSite2) SHARED(A,AcolAvg,ArowAvg,Aavg) num_threads(np)
        !$omp do
        do cgSite1 = 1, nCg
                do cgSite2 = 1,nCg
                       A(cgSite1,cgSite2) = A(cgSite1,cgSite2) - (AcolAvg(cgSite2)+ArowAvg(cgSite1)) + Aavg + 1 
                enddo
        enddo
        !$omp end do nowait
        !$omp end parallel


        BcolAvg=0
        BrowAvg=0
        Bavg=0
        !Compute the distance between cg sites and atoms
!        !$omp parallel private(cgSite1,atom1,j,dist) SHARED(cgPos,atomPos,B,BcolAvg,BrowAvg,Bavg) num_threads(np)
!        !$omp do
        do cgSite1 = 1, nCg
                do atom1 = 1,nAtoms
                        dist = 0
                        do j=1,3
                                dist = dist + (cgPos(cgSite1,j)-atomPos(atom1,j))**2
                        enddo
                        B(cgSite1,atom1) = -sqrt(dist)
                        BcolAvg(atom1) = BcolAvg(atom1) + B(cgSite1,atom1)
                        BrowAvg(cgSite1) = BrowAvg(cgSite1) + B(cgSite1,atom1)
!                        !$omp critical
                        Bavg = Bavg + B(cgSite1,atom1)
!                        !$omp end critical
                enddo
        enddo
!        !$omp end do nowait
!        !$omp end parallel 
        !finish averages
        Bavg = Bavg/dble(nAtoms*nCg)
        BcolAvg = BcolAvg/dble(nCg)
        BrowAvg = BrowAvg/dble(nAtoms)

        !$omp parallel private(cgSite1,atom1) SHARED(B,BcolAvg,BrowAvg,Bavg) num_threads(np)
        !$omp do
        !Subtract means and such 
        do cgSite1 = 1, nCg
                do atom1 = 1,nAtoms
                       B(cgSite1,atom1) = B(cgSite1,atom1) - (BcolAvg(atom1)+BrowAvg(cgSite1)) + Bavg + 1 
                enddo
        enddo
        !$omp end do nowait
        !$omp end parallel 


endsubroutine compute_A_B_matrices



subroutine fit_charges(AtA, AtB, atomCharges, cgCharges, nAtoms, nCg, outputFile)
        implicit none
        integer nAtoms
        integer nCg
        real (kind=8) AtA(nCg,nCg)
        real (kind=8) ATemp(nCg,nCg)
        real (kind=8) AtB(nCg,nAtoms)
        real (kind=8) cgCharges(nCg)
        real (kind=8) atomCharges(nAtoms)
        real (kind=8) newB(nCg)
        !lapack routine variables
        real (kind=8) workQuery(1)
        real (kind=8), allocatable :: work(:)
        character*80 outputFile
        integer info
        integer lwork
        integer j

        !Now use lapack routine to solve least squares problem A*cgCharges=B*atomCharges

        newB = matmul(AtB,atomCharges)
        ATemp=AtA
        call dgels("N",nCg,nCg,1,AtA,nCg,newB,nCg,workQuery,-1,info)
        lwork = int(workQuery(1))
        allocate(work(lwork))
        call dgels("N",nCg,nCg,1,AtA,nCg,newB,nCg,work,lwork,info)
        deallocate(work)
        
        cgCharges = newB

        open(50,file=outputFile)
        do j=1,nCg

                write(50,'(f8.3)') cgCharges(j)

        enddo
        close(50)

        write(*,'(a17,f8.3)') "Total CG Charge:",sum(cgCharges)


endsubroutine fit_charges


