      program tdCIS 
!
!     This program perfoms a real-time time-dependent truncated CI calculation.
!
!     L. M. Thompson, 2021
!
      use mqc_gaussian
      use iso_fortran_env, only: int32, int64, real64
!
!****x* Main/TDCIS
!*    NAME
!*      Real-Time Truncated Configuration Interaction Time Evolution Calculator
!*
!*    SYNOPSIS
!*      Propagates a truncated configuration interaction wavefunction in time
!*      using specified initial conditions.
!
      implicit none
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      character(len=:),allocatable::command,fileName,help_path
      character(len=256),dimension(:),allocatable::fileList
      character(len=256)::vecString,root_string='',sub_string='',field_string='',&
        pulseShape='rectangle',tcf_file='tcf',file_tmp,ci_string='orthogonal',gauss_exe='$g16root/g16/g16',&
        mem='8GB',ncpu='2'
      integer(kind=int64)::iOut=6,iPrint=1,iUnit,i,j,k,maxsteps,flag,stat_num,numFile=0,nullSize
      integer(kind=int64),dimension(:),allocatable::isubs
      real(kind=real64)::delta_t=0.1,field_size=0.005,simTime=100.0,t0=0.0,sigma=0.5,omega=10.0,&
        beta=3.0
      real(kind=real64),allocatable::tcf_start
      complex(kind=real64)::imag=(0.0,1.0)
      logical::UHF,file_exists,doDirect=.false.,found,doProcMem=.false.,keep_intermediate_files=.false.
      type(mqc_pscf_wavefunction)::wavefunction
      type(mqc_molecule_data)::moleculeInfo
      type(mqc_twoERIs),dimension(:),allocatable::eris
      type(mqc_twoERIs)::mo_ERIs
      type(mqc_scalar)::Vnn,final_energy,nIJ,pnIJ,hij,dij
      type(mqc_determinant)::determinants
      type(mqc_scf_integral)::mo_core_ham,density
      type(mqc_scf_integral),dimension(2)::rho
      type(mqc_scf_integral),dimension(3)::dipole,dipoleMO
      type(mqc_scf_integral),dimension(:),allocatable::mo_list
      type(mqc_matrix)::CI_Hamiltonian,exp_CI_Hamiltonian,CI_Overlap,iden,sh2AtMp,shlTyp,nPrmSh,prmExp,&
        conCoef,conCoTwo,shCoor
      type(mqc_matrix),dimension(3)::CI_Dipole,dipole_eigvecs,Xmat,invXmat
      type(mqc_vector)::subs,nuclear_dipole,total_dipole,state_coeffs,td_ci_coeffs,field_vector,&
        td_field_vector,tcf_ci_epsilon,tcf
      type(mqc_vector),dimension(3)::dipole_eigvals
      type(mqc_vector),dimension(2)::nRoot
      real(kind=real64),parameter::zero_thresh=1.00E-8
!
!*    USAGE
!*      TDCIS [-f <matrix_file>] [--print-level <print_level>] [--sub-levels substitutions] 
!*        [--ci-type type_string] [--direct] [--gauss-exe gaussian_executable_string] 
!*        [--do-proc-mem] [--mem mem] [--ncpu ncpu] [--keep-inters] [--intial-state weight_vector] 
!*        [--pulse-shape pulse] [--field-vector field_vector] [--field-size magnitude] [--t0 time] 
!*        [--omega frequency] [--sigma width] [--beta shift] [--tcf-start time] 
!*        [--time-step time_step] [--simulation-time time] [--help]
!*
!*    OPTIONS
!
!
!     Print program information.
!
      write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',73),NEW_LINE('a'), &
          '  Real-Time Truncated Configuration Interaction Time Evolution Calculator',NEW_LINE('a'), &
          ' ',repeat('*',73),NEW_LINE('a'), &
          NEW_LINE('a'),repeat(' ',30),'Version 21.10.1',NEW_LINE('a'),NEW_LINE('a'),&
          ' L. M. Thompson, Louisville KY, 2021.',NEW_LINE('a')
!
!     Parse input options.
!
      j = 1
      do i=1,command_argument_count()
        if(i.ne.j) cycle
        call mqc_get_command_argument(i,command)
        if(command.eq.'-f') then
!
!*      -f matrix_file                   Input matrix file with initial set of molecular orbitals.
!*
          call mqc_get_command_argument(i+1,fileName)
          j = i+2
        elseif(command.eq.'--print-level') then
!
!*      --print-level print_level        Verbosity of output. Default print level is 1. Options
!*                                       0-4.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I1)') iPrint
          j = i + 2
        elseIf(command.eq.'--sub-levels') then
!
!*      --sub-levels substitutions       Substitution levels permitted in truncated CI calculation. 
!*
!*                                       Example: [1,2] specifies single and double substitutions.
!*
          call mqc_get_command_argument(i+1,command)
          sub_string = command
          j = i+2
        elseIf(command.eq.'--ci-type') then
!
!*      --ci-type type_string            Specifies the type of configuration interaction. Options
!*                                       are:
!*                                       1) orthogonal (default)
!*                                          Perform orthogonal configuration interaction with 
!*                                          determinant expansion specified by sub-levels option.
!*                                          Only one matrix file should be specified in the input 
!*                                          file is expected and additional inputs will result in 
!*                                          an error.
!*                                       2) nonorthogonal
!*                                          Perform nonorthogonal configuration interaction with 
!*                                          determinant expansion specified by matrix files 
!*                                          listed in input file.
!*                                       3) hybrid (not yet implemented)
!*                                          Perform the orthogonal expansion specified by the sub-
!*                                          levels option on each matrix file in the input file
!*                                          which are generally nonorthogonal.
!*
          call mqc_get_command_argument(i+1,command)
          ci_string = command
          j = i+2
        elseif(command.eq.'--direct') then
 !
 !*      --direct                         Avoid the use of 2ERIs in nonorthogonal code. Gaussian 
 !*                                       is required if requested.
 !*
          doDirect=.true.
          j = i + 1
        elseIf(command.eq.'--gauss-exe') then
!
!*      --gauss-exe gaussian_executable  Path, executable and command line flags for Gaussian
!*                                       executable used when running Gaussian jobs from within
!*                                       the program.
!*
          call mqc_get_command_argument(i+1,command)
          gauss_exe = command
          j = i+2
        elseIf(command.eq.'--do-proc-mem') then
!
!*      --do-proc-mem                    Request %mem and %nproc lines be added to the output
!*                                       Gaussian input files if doing direct matrix elements.
!*
          doProcMem = .true.
          j = i+1
        elseIf(command.eq.'--mem') then
!
!*      --mem mem                        Amount of memory to be included in %mem output to
!*                                       Gaussian input files. --do-proc-mem must be set.
!*                                       Default is 8GB. Note input is not checked for validity.
!*
          call mqc_get_command_argument(i+1,command)
          mem = command
          j = i+2
        elseIf(command.eq.'--nproc') then
!
!*      --ncpu ncpu                      Number of processors to be included in %nproc output
!*                                       to Gaussian input files. --do-proc-mem must be set.
!*                                       Default is 2. Note input is not checked for validity.
!*
          call mqc_get_command_argument(i+1,command)
          ncpu = command
          j = i+2
        elseif(command.eq.'--keep-inters') then
!
!*      --keep-inters                    Do not delete intemediate files used in the calculation.
!*                                       Default is to delete intermediate files.
!*
          keep_intermediate_files=.true.
          j = i + 1
        elseIf(command.eq.'--initial-state') then
!
!*      --initial-state weight_vector    Initial state for simulation. For a pure state, the  
!*                                       desired root can be input as a single integer. For a 
!*                                       mixed state, each non-zero weighted root should be 
!*                                       included along wih the specified weight (Default is
!*                                       a ground state population of 1.0). Note that the input 
!*                                       vector is normalized regardles of input values.
!*                                   
!*                                       Example: [1,0.5:2,0.5] specifies an initial state with 
!*                                       50% weight on the lowest root and 50% weight on the 
!*                                       fourth root.
!*
          call mqc_get_command_argument(i+1,command)
          root_string = command
          j=i+2
        elseIf(command.eq.'--pulse-shape') then
!
!*      --pulse-shape pulse              Pulse shape used in simulation. Options are:
!*                                       1) rectangle
!*                                          E(t) = E(0)*H(t-t0)
!*                                       2) delta
!*                                          E(t) = E(0)*d(t-t0) 
!*                                       3) continuous 
!*                                          E(t) = E(0)*sin(w(t-t0))*H(t-t0)
!*                                       4) transform limited
!*                                          E(t) = E(0)*exp(-(t-t0)^2/2*sigma^2)*sin(omega*(t-t0))
!*                                       5) chirped pulse
!*                                          E(t) = E(0)*exp(-(t-t0)^2/2*sigma^2)*
!*                                                   sin((omega+beta(t-t0))*(t-t0))
!*
          call mqc_get_command_argument(i+1,command)
          pulseShape = command
          j = i+2
        elseIf(command.eq.'--field-vector') then
!
!*      --field-vector field_vector      Field polarization vector. Should be vector of length
!*                                       three. The default is to orient on the z-axis. Note that
!*                                       the initial vector is normalized regardless of input 
!*                                       values.
!*
!*                                       Example: [0.5,0,-0.5] specifies a field orientation of 
!*                                       1/sqrt(2) along the x axis and -1/sqrt(2) along the z axis.
!*
          call mqc_get_command_argument(i+1,command)
          field_string = command
          j = i+2
        elseif(command.eq.'--field-size') then
!
!*      --field-size field_magnitude     Maximum magnitude of field vector in simulation (au)
!*                                       (default is 0.005 au = 25.7 MV/cm = 1.752 W/m^2 in free
!*                                       space conditions).
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F12.6)') field_size
          j = i + 2
        elseif(command.eq.'--t0') then
!
!*      --t0 time                        Time for either pulse onset or pulse maximum depending on
!*                                       pulse type (default is 0.0 au).
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F12.6)') t0
          j = i + 2
        elseif(command.eq.'--omega') then
!
!*      --omega frequency                Frequency for pulse shapes that use it (default is 10.0 au). 
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F12.6)') omega
          j = i + 2
        elseif(command.eq.'--sigma') then
!
!*      --sigma width                    Pulse width for pulse shapes that use it (default is 0.5 au). 
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F12.6)') sigma
          j = i + 2
        elseif(command.eq.'--beta') then
!
!*      --beta shift                     Chirp parameter for pulse shapes that use it (default is 3.0 au). 
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F12.6)') sigma
          j = i + 2
        elseif(command.eq.'--tcf-start') then
!
!*      --tcf-start time                 Compute time correlation function starting from input time. The 
!*                                       Fourier transform of this function when the start time is the time 
!*                                       step after the delta pulse provides the absorption spectrum.
!*
          call mqc_get_command_argument(i+1,command)
          allocate(tcf_start)
          read(command,'(F12.6)') tcf_start
          j = i + 2
        elseif(command.eq.'--time-step') then
!
!*      --time-step time_step            Time step for simulation in au (default 0.1 au = 2.42 as).
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F12.6)') delta_t
          j = i + 2
        elseif(command.eq.'--simulation-time') then
!
!*      --simulation-time time           Total simulation time in au (default is 100.0 au = 2.42 fs).
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F20.8)') simTime
          j = i + 2
        elseIf(command.eq.'--help') then
!
!*      --help                           Output help documentation to terminal.
!*
          if(command_argument_count().gt.1) call mqc_error_I('Help output requested with multiple arguments',6, &
            'command_argument_count()',command_argument_count())
          call mqc_get_command_argument(0,help_path)
          help_path = 'less ' // trim(help_path(1:scan(help_path,'/',.true.))) // 'doc/TDCIS.txt'
          call execute_command_line(help_path,exitstat=flag)
          if(flag.ne.0) call mqc_error('Help output command failed')
          stop
        else
          call mqc_error_A('Unrecognised input flag',6,'command',command)
        endIf
      endDo
!
!     Parse input file and extract required data from matrix files.
!
      if(.not.allocated(fileName)) call mqc_error('No input file provided',iOut)
      open(newunit=iUnit,file=fileName,status='old',iostat=stat_num)
      if(stat_num/=0) call mqc_error_a('Error opening file',iOut,'fileName',fileName)
      read(unit=iUnit,fmt='(i20)',iostat=stat_num) numFile
      if(stat_num/=0) call mqc_error('Error reading file number',iOut)
      allocate(fileList(numFile))
      do i = 1, numFile
        read(unit=iUnit,fmt='(A)',iostat=stat_num) fileList(i)
        if((stat_num<0).and.(i<=numFile)) call mqc_error('File EOF reached early',iOut)
      endDo
      close(unit=iUnit)
!
      if(ci_string.eq.'orthogonal'.and.numFile.ne.1) then
        call mqc_error('Multiple matrix files input to requested orthogonal CI expansion.&
          & Use sub-levels option to provide expansion',iOut)
      elseIf(ci_string.eq.'orthogonal'.or.ci_string.eq.'nonorthogonal') then
        allocate(mo_list(numFile))
        do i = 1, numFile
          call fileInfo%getESTObj('mo coefficients',est_integral=mo_list(i),filename=fileList(i))
          if(iPrint.ge.4) call mo_list(i)%print(iOut,'MO coefficients from matrix file '//trim(num2char(i)))
        endDo
      elseIf(ci_string.eq.'hybrid') then
        call mqc_error('Hybrid ci type is not yet implemented',iOut)
      else
        call mqc_error_a('Unrecognized CI type string provided',iOut,'ci_string',ci_string)
      endIf
!
      call fileInfo%load(fileList(1))
      call fileInfo%getMolData(moleculeInfo)
      call fileInfo%getESTObj('wavefunction',wavefunction)
      if(iPrint.ge.4) call wavefunction%print(iOut,'all')
      call fileInfo%getESTObj('dipole x',est_integral=dipole(1))
      if(iPrint.ge.4) call dipole(1)%print(6,'AO dipole x integrals')
      call fileInfo%getESTObj('dipole y',est_integral=dipole(2))
      if(iPrint.ge.4) call dipole(2)%print(6,'AO dipole y integrals')
      call fileInfo%getESTObj('dipole z',est_integral=dipole(3))
      if(iPrint.ge.4) call dipole(3)%print(6,'AO dipole z integrals')
!
      if(.not.doDirect) then
        if(ci_string.eq.'orthogonal') then
          allocate(eris(1))
        else
          allocate(eris(3))
        endIf
        call fileInfo%get2ERIs('regular',eris(1),foundERI=found)
        if(found.and.ci_string.ne.'orthogonal') then
          eris(2) = eris(1)
          eris(3) = eris(1)
          call mqc_twoeris_transform(eris(1),'raffenetti1')
          call mqc_twoeris_transform(eris(2),'raffenetti2')
          call mqc_twoeris_transform(eris(3),'raffenetti3')
        elseIf(ci_string.ne.'orthogonal') then
          call fileInfo%get2ERIs('raffenetti1',eris(1),foundERI=found)
          if(found) call fileInfo%get2ERIs('raffenetti2',eris(2),foundERI=found)
          if(found) call fileInfo%get2ERIs('raffenetti3',eris(3),foundERI=found)
        endIf
        if(.not.found) call mqc_error_l('Integrals were not found on matrix file',6, &
          'found',found)
        if(iPrint.ge.4) then
          do i = 1, size(eris) 
            call eris(i)%print(iOut,'AO 2ERIs set '//trim(num2char(i)))
          endDo
        endIf
      endIf
!
      if(ci_string.ne.'nonorthogonal') then
        if(wavefunction%wf_type.eq.'U') then
          UHF = .true.
          if(iPrint.ge.3) write(iOut,'(1X,A)') 'Found UHF wavefunction'//NEW_LINE('A')
        elseIf(wavefunction%wf_type.eq.'R') then
          UHF = .False.
          if(iPrint.ge.3) write(iOut,'(1X,A)') 'Found RHF wavefunction'//NEW_LINE('A')
        else
          call mqc_error_A('Unsupported wavefunction type',iOut, &
            'Wavefunction%wf_type',wavefunction%wf_type)
        endIf 
        if (wavefunction%wf_complex) call mqc_error('Complex wavefunctions unsupported')
      endIf
!
      if(doDirect) then
        call fileInfo%getArray('SHELL TO ATOM MAP',sh2AtMp)
        call fileInfo%getArray('SHELL TYPES',shlTyp)
        call fileInfo%getArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
        call fileInfo%getArray('PRIMITIVE EXPONENTS',prmExp)
        call fileInfo%getArray('CONTRACTION COEFFICIENTS',conCoef)
        call fileInfo%getArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
        call fileInfo%getArray('COORDINATES OF EACH SHELL',shCoor)
        if(ci_string.eq.'orthogonal') call mqc_error('Direct matrix elements not possible with orthogonal CI')
      endIf
!
!
!     Compute the nuclear-nuclear repulsion energy.
!
      call moleculeInfo%print(iOut)
      Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
      call Vnn%print(iOut,'Nuclear Repulsion Energy (au)',Blank_At_Bottom=.true.) 
!
!     Generate Slater determinants wavefunction expansion if orthogonal expansion requested. 
!
      if(ci_string.ne.'nonorthogonal') then
        call substitution_builder(sub_string,subs)
        if(iPrint.ge.1) then
          write(iOut,'(1X,A)') 'Building Determinant Strings'
          call subs%print(6,'Permitted substitution levels',Blank_At_Bottom=.true.)
        endIf
        isubs = [(i, i=1,int(maxval(subs)))]
        call trci_dets_string(iOut,iPrint,wavefunction%nBasis,wavefunction%nAlpha, &
          Wavefunction%nBeta,isubs,determinants)
      endIf
!
!     Transform one and two-electron integrals to MO basis (only if performing orthogonal CI).
!
      if(ci_string.ne.'nonorthogonal') then
        if(iPrint.ge.1) write(iOut,'(1X,A)') 'Transforming MO integrals'//NEW_LINE('A')
        mo_core_ham = matmul(transpose(wavefunction%MO_Coefficients),matmul(wavefunction%core_Hamiltonian, &
            Wavefunction%MO_Coefficients))
        if(IPrint.ge.4) call mo_core_ham%print(iOut,'MO Basis Core Hamiltonian') 
        call twoERI_trans(iOut,iPrint,wavefunction%MO_Coefficients,ERIs(1),mo_ERIs)
      endIf
!
!     Compute nuclear dipole moment.
!
      nuclear_dipole = matmul(transpose(moleculeInfo%Nuclear_Charges),&
        transpose(moleculeInfo%Cartesian_Coordinates))
      if(iPrint.ge.1) call nuclear_dipole%print(6,'Nuclear dipole',Blank_At_Bottom=.true.)
      call total_dipole%init(3)
!
!     Generate field-free Hamiltonian matrix and dipole matrices (in length form) in CI basis.
!
      if(ci_string.eq.'orthogonal') then
        if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building orthogonal CI Hamiltonian matrix'
        call subs%unshift(0)
        isubs = subs
        call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis,determinants, &
          mo_core_ham,mo_ERIs,UHF,CI_Hamiltonian,isubs)
        if(iPrint.ge.4) call CI_Hamiltonian%print(6,'CI Hamiltonian',Blank_At_Bottom=.true.)
        if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building orthogonal CI dipole matrices'
        do i = 1, 3
          dipoleMO(i) = matmul(dagger(wavefunction%MO_Coefficients),&
            matmul(dipole(i),Wavefunction%MO_Coefficients))
          if(iprint.ge.4) call dipoleMO(i)%print(6,'MO dipole integrals axis '//trim(num2char(i)),&
            Blank_At_Bottom=.true.)
          call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis,determinants, &
            dipoleMO(i),UHF=UHF,CI_Hamiltonian=CI_Dipole(i),subs=isubs)
          if(iprint.ge.4) call CI_Dipole(i)%print(6,'SD dipole integrals axis '//trim(num2char(i)),&
            Blank_At_Bottom=.true.)
          call iden%identity(size(CI_Dipole(i),1),size(CI_Dipole(i),2),nuclear_dipole%at(i))
          CI_Dipole(i) = iden - CI_Dipole(i)
          if(iprint.ge.4) call CI_Dipole(i)%print(6,'SD dipole including nuclear term axis '//&
            trim(num2char(i)),Blank_At_Bottom=.true.)
        endDo
      elseIf(ci_string.eq.'nonorthogonal') then
        if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building nonorthogonal CI Hamiltonian, CI overlap and CI dipole matrices'
!       loop over pairs of Slater determinants and build transition density matrices. Get the H, N and Mu CI matrices
        call CI_Hamiltonian%init(numFile,numFile,storage='StorHerm') 
        call CI_Overlap%init(numFile,numFile,storage='StorHerm') 
        do i = 1, 3
          call CI_Dipole(i)%init(numFile,numFile,storage='StorHerm') 
        endDo
        do i = 1, numFile
          do j = 1, i
            call get_rhos(rho,nIJ,pnIJ,nullSize,mo_list(i),mo_list(j),wavefunction%overlap_matrix,wavefunction%nBasis,&
              wavefunction%nAlpha,wavefunction%nBeta)
            call CI_Overlap%put(nIJ,i,j,'hermitian')
            hij = get_hij(pnij,nullSize,rho,wavefunction%core_Hamiltonian,eris,doDirect,fileInfo,fileName,&
              symIndexHash(i,j),sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor,gauss_exe,doProcMem,&
              mem,ncpu,keep_intermediate_files)
            call CI_Hamiltonian%put(hIJ,i,j,'hermitian')
            do k = 1, 3
              dij = get_dij(pnij,rho(1),dipole(k))
              call CI_Dipole(k)%put(dij,i,j,'hermitian')
            endDo
          endDo
        endDo
        if(iPrint.ge.4) call CI_Overlap%print(6,'CI Overlap',Blank_At_Bottom=.true.)
        if(iPrint.ge.4) call CI_Hamiltonian%print(6,'CI Hamiltonian',Blank_At_Bottom=.true.)
        do i = 1, 3
          call iden%identity(size(CI_Dipole(i),1),size(CI_Dipole(i),2),nuclear_dipole%at(i))
          CI_Dipole(i) = iden - CI_Dipole(i)
          if(iPrint.ge.4) call CI_Dipole(i)%print(6,'SD dipole including nuclear term axis '//&
            trim(num2char(i)),Blank_At_Bottom=.true.)
        endDo
      endIf
!
!     Diagonalize Hamiltonian
!
      if(iPrint.ge.1) write(iOut,'(1X,A)') 'Diagonalizing CI Hamiltonian'//NEW_LINE('A')
      if(ci_string.eq.'orthogonal') then
        call CI_Hamiltonian%diag(wavefunction%pscf_energies,wavefunction%pscf_amplitudes)
      elseIf(ci_string.eq.'nonorthogonal') then
        call CI_Hamiltonian%eigensys(CI_Overlap,wavefunction%pscf_energies,wavefunction%pscf_amplitudes)
      endIf

      if(iPrint.ge.4) call wavefunction%pscf_amplitudes%print(iOut,'CI Eigenvectors',Blank_At_Bottom=.true.)
      if(iPrint.ge.3) call wavefunction%pscf_energies%print(iOut,'CI Eigenvalues',Blank_At_Bottom=.true.)
!
!     Diagonalize the dipole moment matrix to get transformation matrix U
!
      do i = 1, 3
        if(iPrint.eq.1) write(iOut,'(1X,A)') 'Diagonalizing SD dipole axis '//trim(num2char(i))//NEW_LINE('A')
        if (ci_string.eq.'orthogonal') then
          call CI_Dipole(i)%diag(dipole_eigvals(i),dipole_eigvecs(i))
        elseIf(ci_string.eq.'nonorthogonal') then
          call CI_Dipole(i)%eigensys(CI_Overlap,dipole_eigvals(i),dipole_eigvecs(i))
        endIf
        if(iPrint.ge.4) call dipole_eigvecs(i)%print(iOut,'CI Dipole Eigenvectors axis '//trim(num2char(i)),&
          Blank_At_Bottom=.true.)
        if(iprint.ge.3) call dipole_eigvals(i)%print(iOut,'CI Dipole Eigenvalues axis '//trim(num2char(i)),&
          Blank_At_Bottom=.true.)
        if (ci_string.eq.'orthogonal') then
          Xmat(i) = matmul(dagger(dipole_eigvecs(i)),wavefunction%pscf_amplitudes)
          invXmat(i) = dagger(Xmat(i))
        elseIf(ci_string.eq.'nonorthogonal') then
          Xmat(i) = matmul(mqc_matrix_inverse(dipole_eigvecs(i)),wavefunction%pscf_amplitudes)
          invXmat(i) = mqc_matrix_inverse(Xmat(i))
        endIf
      endDo
!
!     Build the electric field pulse vector
!
      call field_vector_builder(field_string,field_vector)
      field_vector = field_vector*field_size
!
!     Determine initial conditions in the basis of CI states.
!
      call root_vector_builder(root_string,size(wavefunction%pscf_energies),nRoot)
      call state_coeffs%init(size(wavefunction%pscf_energies),0.0)
      maxsteps = simTime/delta_t 
      do i = 1, size(nroot(1))
        call state_coeffs%put(nRoot(2)%at(i),nRoot(1)%at(i))
      endDo
      state_coeffs = state_coeffs/state_coeffs%norm()
      vecString = '# 4) State weights:'//repeat(' ',24)//'#'//NEW_LINE('a')
      do i = 1, size(nRoot(1))
        vecString = trim(vecString)//' #'//repeat(' ',26)//trim(num2char(nRoot(1)%at(i),'I3'))//' = '//&
          trim(num2char(state_coeffs%at(nRoot(1)%at(i)),'F7.4'))//repeat(' ',6)//'#'
        if(i.ne.size(nRoot(1))) vecString = trim(vecString)//NEW_LINE('a')
      endDo
      write(iOut,'(1X,A)')               '############################################'
      write(iOut,'(1X,A)')               '#                                          #'
      write(iOut,'(1X,A)')               '#            INITIAL CONDITIONS            #'
      write(iOut,'(1X,A)')               '#            ------------------            #'
      write(iOut,'(1X,A)')               '#                                          #'
      write(iOut,'(1X,A,1X,F15.6,A)')      '# 1) Simulation time:',simTime,' au   #'
      write(iOut,'(1X,A,1X,I15,A)')      '# 2) Number of steps:',maxsteps,'      #'
      write(iOut,'(1X,A,1X,F21.6,1X,A)') '# 3) Time step:',delta_t,'au   #'
      write(iOut,'(1X,A)')               trim(vecString)
      write(iOut,'(1X,A)')               '# 5) Field parameters:                     #' 
      write(iOut,'(1X,A,1X,A19,A)')      '#    Pulse shape:',trim(pulseShape),'      #'
      write(iOut,'(1X,A,1X,F21.6,A)')      '#    Magnitude:',field_size,' au   #'
      write(iOut,'(1X,A,1X,F23.6,A)')      '#    Field X:',float(field_vector%at(1)),'      #'
      write(iOut,'(1X,A,1X,F23.6,A)')      '#    Field Y:',float(field_vector%at(2)),'      #'
      write(iOut,'(1X,A,1X,F23.6,A)')      '#    Field Z:',float(field_vector%at(3)),'      #'
      write(iOut,'(1X,A,1X,F28.6,A)')      '#    t0:',t0,'      #'
      write(iOut,'(1X,A,1X,F25.6,A)')      '#    omega:',omega,'      #'
      write(iOut,'(1X,A,1X,F25.6,A)')      '#    sigma:',sigma,'      #'
      write(iOut,'(1X,A,1X,F26.6,A)')      '#    beta:',beta,'      #'
      write(iOut,'(1X,A)')               '#                                          #'
      write(iOut,'(1X,A)')               '############################################'
      write(iOut,'(1X,A)')               ''
!
!     Do loops over time-steps and propagate CI coefficients in each time step using eq. 16 of 
!     Krause et al. J. Phys. Chem., 2007, 127, 034107.
!
      write(iOut,'(1X,A)') 'STARTING TIME PROPAGATION'//NEW_LINE('A')
      do i = 1, maxsteps
        if(iPrint.ge.1) write(iOut,'(1X,A)') repeat('=',44)//NEW_LINE('A')
        if(i.eq.maxsteps) write(iOut,'(1X,A)') repeat(' ',16)//'FINAL TIME STEP'//NEW_LINE('A')//NEW_LINE('A')//&
        ' '//repeat('=',44)//NEW_LINE('A')
        if(iPrint.ge.1) write(iOut,'(1X,A,1X,F12.6,1X,A)') 'Time step:',delta_t*i,'au'//NEW_LINE('A')

        td_field_vector = get_field_vector(delta_t,i,field_vector,pulseShape,t0,omega,sigma,beta)
        if(iPrint.ge.1) call td_field_vector%print(6,'Applied field vector',Blank_At_Bottom=.true.)

        state_coeffs = exp((-1)*imag*delta_t*wavefunction%pscf_energies).ewp.state_coeffs
        do j = 3, 1, -1
          state_coeffs = matmul(Xmat(j),state_coeffs)
          state_coeffs = exp(imag*delta_t*td_field_vector%at(j)*dipole_eigvals(j)).ewp.state_coeffs
          state_coeffs = matmul(invXmat(j),state_coeffs)
        endDo
        call print_coeffs_and_pops(iOut,iPrint,1,state_coeffs,'TD State')

        call td_ci_coeffs%init(size(wavefunction%pscf_amplitudes,1))
        do j = 1, size(state_coeffs)
          td_ci_coeffs = td_ci_coeffs + state_coeffs%at(j)*wavefunction%pscf_amplitudes%vat([0],[j])
        endDo
        call print_coeffs_and_pops(iOut,iPrint,2,td_ci_coeffs,'TD CI')

        if(allocated(tcf_start)) then
          if(delta_t*i.ge.tcf_start.and.delta_t*i.lt.tcf_start+delta_t) tcf_ci_epsilon = td_ci_coeffs
          if(delta_t*i.ge.tcf_start) call tcf%push(dot_product(dagger(tcf_ci_epsilon),td_ci_coeffs))
        endIf

        final_energy = get_CI_Energy(CI_Hamiltonian,td_ci_coeffs) 
        if(iPrint.ge.1.or.i.eq.maxsteps) call final_energy%print(6,'Energy (au)',Blank_At_Bottom=.true.,&
          FormatStr='F14.8')

        density = get_one_gamma_matrix(iOut,iPrint,wavefunction%nBasis,determinants,td_ci_coeffs,UHF,subs=isubs)
        if(iPrint.ge.1.or.i.eq.maxsteps) call density%print(6,'MO Density matrix',Blank_At_Bottom=.true.)
        density = 0.5*matmul(matmul(wavefunction%mo_coefficients,density),dagger(wavefunction%mo_coefficients))
        if(iPrint.ge.1.or.i.eq.maxsteps) call density%print(6,'AO Density matrix',Blank_At_Bottom=.true.)
        do j = 1, 3
          call total_dipole%put((-1)*contraction(density,dipole(j)) + nuclear_dipole%at(j),j)
        endDo
        call total_dipole%print(6,'Total dipole',Blank_At_Bottom=.true.)
      endDo

      if(allocated(tcf_start)) then
        i = 1
        file_exists = .true.
        file_tmp = trim(tcf_file)
        do while (file_exists) 
          inquire(file=trim(file_tmp)//'.dat',exist=file_exists)
          if(file_exists) then
            i = i+1
            file_tmp = trim(tcf_file)
            call build_string_add_int(i,file_tmp,20)
          else
            open(newunit=iunit,file=trim(file_tmp)//'.dat',status='new',iostat=stat_num)
            if(stat_num.ne.0) call mqc_error_a('Could not open file',6,'file name',trim(file_tmp)//'.dat')
            call tcf%print(iUnit,'# time ('//trim(num2char(delta_t))//'as)'//repeat(' ',10)//'time correlation function')
            write(iOut,'(A)') 'Saving data to file name '//trim(file_tmp)//'.dat'
          endIf
        endDo
      endIf

!
      contains
!
!     
!     PROCEDURE get_CI_energy
!
!     get_CI_energy is a function that returns the energy given the Hamiltonian and CI vectors.
!
      function get_CI_Energy(CI_Hamiltonian,CI_vectors) result(energy)

      implicit none
      type(mqc_matrix),intent(in)::CI_Hamiltonian
      type(mqc_vector),intent(in)::CI_Vectors
      type(mqc_scalar)::energy

      energy = dot_product(matmul(dagger(CI_vectors),CI_Hamiltonian),CI_vectors)

      end function get_CI_energy
!
!     
!     PROCEDURE get_field_vector
!
!     get_field_vector is a function that returns the field vector at time dt*step. 
!
      function get_field_vector(dt,step,static_field,pulseShape,t0,omega,sigma,beta) result(td_field)

      implicit none
      real(kind=real64),intent(in)::dt,t0,omega,sigma,beta 
      integer(kind=int64),intent(in)::step
      type(mqc_vector),intent(in)::static_field
      type(mqc_vector)::td_field
      character(len=*)::pulseShape

      select case(pulseShape)
      case('rectangle')
        if(dt*step.ge.t0.and.dt*step.lt.t0+sigma) then
          td_field = static_field
        else
          td_field = [0.0,0.0,0.0]
        endIf
      case('delta')
        if(dt*step.ge.t0.and.dt*step.lt.t0+dt) then
          td_field = static_field
        else
          td_field = [0.0,0.0,0.0]
        endIf
      case('continuous')
        if(dt*step.ge.t0.and.dt*step.lt.t0+sigma) then
          td_field = static_field*sin(omega*(dt*step-t0))
        else
          td_field = [0.0,0.0,0.0]
        endIf
      case('transform limited')
        td_field = static_field*exp((-(dt*step-t0)**2)/(2*sigma**2))*&
          sin(omega*(dt*step-t0))
      case('chirped pulse')
        td_field = static_field*exp((-(dt*step-t0)**2)/(2*sigma**2))*&
          sin((omega+beta*(dt*step-t0))*(dt*step-t0))
      case default
        call mqc_error_a('Unrecognized pulse shape requested',6,'pulseShape',pulseShape)
      end select

      end function get_field_vector
!
!
!     PROCEDURE substitution_builder
!
!     substitution_builder is a subroutine that builds the input wavefunction substituion levels.
!
      subroutine substitution_builder(subs_in,subs_out)

      implicit none

!     input/output variables
      character(len=*),intent(in)::subs_in
      type(mqc_vector),intent(inOut)::subs_out

!     text parsing variables
      character(len=:),allocatable::subsString
      character(len=10)::val
      integer::i,ival
      logical::newNum=.false.

      subsString = trim(subs_in)
      if(len(subsString).eq.0) then
        write(iOut,'(1X,A)') 'Defaulting to CIS.'//NEW_LINE('A')
        subs_out = [1]
      else
        do i = 1,len(subsString)
          select case (subsString(i:i))
          case('[','\(')
            if(i.eq.1) then
              newNum = .true.
              cycle
            else
              call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
            endIf
          case(']','\)')
            if(i.eq.len(subsString).and..not.newNum.and.i.ne.1) then
              read(val,'(I10)') ival
              call subs_out%push(ival)
              newNum = .true.
              cycle
            else
              call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
            endIf
          case(',',' ')
            if(i.eq.1.or.i.eq.len(subsString).or.i.eq.2.or.i.eq.len(subsString)-1.or.newNum) then
              call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
            else
              read(val,'(I10)') ival
              call subs_out%push(ival)
              newNum = .true.
              cycle
            endIf
          case('0':'9')
            if(i.eq.1.or.i.eq.len(subsString)) then
              call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
            else
              if(newNum) then
                val = subsString(i:i)
              else
                val = trim(val)//subsString(i:i)
              endIf
              newNum = .false.
              cycle
            endIf
          case default
            call mqc_error_A('Substitution level input format incorrect',6,'subsString',subsString)
          end select
        endDo
      endIf

      end subroutine substitution_builder
!    
!
!     PROCEDURE field_vector_builder
!
!     field_vector_builder is a subroutine that builds the normalized input t=0 field vector.
!
      subroutine field_vector_builder(field_in,field_out)

      implicit none

!     input/output variables
      character(len=*),intent(in)::field_in
      type(mqc_vector),intent(inOut)::field_out

!     text parsing variables
      character(len=:),allocatable::fieldString
      character(len=12)::val
      integer::i
      real(kind=real64)::rval
      logical::newNum=.false.

      fieldString = trim(field_in)
      if(len(fieldString).eq.0) then
        write(iOut,'(1X,A)') 'Defaulting to t=0 field aligned along z-axis.'//NEW_LINE('A')
        field_out = [0,0,1]
      else
        do i = 1,len(fieldString)
          select case (fieldString(i:i))
          case('[','\(')
            if(i.eq.1) then
              newNum = .true.
              cycle
            else
              call mqc_error_A('Field (t=0) vector input format incorrect',6,'fieldString',fieldString)
            endIf
          case(']','\)')
            if(i.eq.len(fieldString).and..not.newNum.and.i.ne.1) then
              read(val,'(F12.6)') rval
              call field_out%push(rval)
              newNum = .true.
              cycle
            else
              call mqc_error_A('Field (t=0) vector input format incorrect',6,'fieldString',fieldString)
            endIf
          case(',',' ')
            if(i.eq.1.or.i.eq.len(fieldString).or.i.eq.2.or.i.eq.len(fieldString)-1.or.newNum) then
              call mqc_error_A('Field (t=0) vector input format incorrect',6,'fieldString',fieldString)
            else
              read(val,'(F12.6)') rval
              call field_out%push(rval)
              newNum = .true.
              cycle
            endIf
          case('0':'9','.')
            if(i.eq.1.or.i.eq.len(fieldString)) then
              call mqc_error_A('Field (t=0) vector input format incorrect',6,'fieldString',fieldString)
            else
              if(newNum) then
                val = fieldString(i:i)
              else
                val = trim(val)//fieldString(i:i)
              endIf
              newNum = .false.
              cycle
            endIf
          case('-')
            if(i.eq.1.or.i.eq.len(fieldString).or..not.newNum) then
              call mqc_error_A('Field (t=0) vector input format incorrect',6,'fieldString',fieldString)
            else
              val = fieldString(i:i)
              newNum = .false.
              cycle
            endIf
          case default
            call mqc_error_A('Field (t=0) vector input format incorrect',6,'fieldString',fieldString)
          end select
        endDo
      endIf

      if(size(field_out).ne.3) call mqc_error_I('input field (t=0) vector should be length 3',6,'size(field_out)',&
        size(field_out))
      field_out = field_out/field_out%norm()

      end subroutine field_vector_builder
!    
!
!     PROCEDURE root_vector_builder
!     
!     root_vector_builder is a subroutine builds the vector of initial state weights.
!
      subroutine root_vector_builder(root_string_in,detLen,root_vector_out)

      implicit none

!     input/output variables
      character(len=*),intent(in)::root_string_in
      integer,intent(in)::detLen
      type(mqc_vector),dimension(2),intent(inOut)::root_vector_out

!     text parsing variables
      character(len=:),allocatable::rootString
      character(len=80)::processString
      integer::i,root
      real(kind=real64)::weight
      logical::rootNum,weightNum,rootSet,weightSet,parenStr

      rootString = trim(root_string_in)
      if(len(rootString).eq.0) then
        call root_vector_out(1)%push(1)
        call root_vector_out(2)%push(1.0)
      else
        rootNum = .false.
        rootSet = .false.
        weightNum = .false.
        weightSet = .false.
        parenStr = .false.
        processString = ''
        do i = 1, len(rootString)
          select case(rootString(i:i))
          case('[')
            if(i.eq.1) then
              rootNum = .true.
              parenStr = .true.
              cycle
            else
              call mqc_error("Malformed root weight string.  Expected format:  [#,#:#,#:...]")
            end if
          case("0":"9")
            processString = trim(processString)//rootString(i:i)
          case(",")
            if(.not.rootNum.or..not.parenStr) then
              call mqc_error("Malformed root weight string.  Expected format:  [#,#:#,#:...]")
            else
              read(processString,'(I3)') root
              processString = ''
              rootSet = .true.
              rootNum = .false.
              weightNum = .true.
            end if
          case(".")
            if(.not.weightNum.or..not.parenStr) then
              call mqc_error("Malformed root weight string.  Expected format:  [#,#:#,#:...]")
            else
              processString = trim(processString)//rootString(i:i)
            end if
          case(":")
            if(rootSet.and.weightNum.and.(.not.weightSet).and.parenStr) then
              read(processString,'(F10.3)') weight
              processString = ''
              weightSet = .true.
              rootNum = .true.
              weightNum = .false.
            else
              call mqc_error("Malformed root weight string.  Expected format:  [#,#:#,#:...]")
            end if
          case("]")
            if(.not.weightNum.or..not.parenStr) then
              call mqc_error("Malformed root weight string.  Expected format:  [#,#:#,#:...]")
            else
              read(processString,'(F10.3)') weight
              processString = ''
              weightSet = .true.
              weightNum = .false.
            end if
          case default
            call mqc_error("Unrecognized character in root weight input:  "//rootString(i:i))
          end select
          if(i.eq.len(rootString).and..not.parenStr) then
            read(processString,'(I3)') root
            weight = 1.0
            rootSet = .true.
            weightSet = .true.
          endIf

          !  Algorithm has identified a valid pair of numbers for processing.
          if(rootSet.and.weightSet) then
            if(root.gt.detLen) call mqc_error("Root requested exceeds dimension of Hamiltonian matrix.") 
            if(weight.gt.1.0e0) call mqc_error("Root weight is greater than one.") 
!            call root_vector_out%put(weight,root)
            call root_vector_out(1)%push(root)
            call root_vector_out(2)%push(weight)
            rootSet = .false.
            weightSet = .false.
            rootNum = .true.
            weightNum = .false.
          endIf
        endDo
      endIf

      end subroutine root_vector_builder
!    
!
!     PROCEDURE print_coeffs_and_pops
!     
!     print_coeffs_and_pops is a subroutine that prints coefficients and 
!     populations of a vector, using the input string to build the label
!
      subroutine print_coeffs_and_pops(iOut,iPrint,printThresh,vector,titleString)

      implicit none

!     input/output variables
      integer(kind=int64),intent(in)::iOut,iPrint,printThresh
      type(mqc_vector),intent(in)::vector
      character(len=*),intent(in)::titleString

!     subroutine variables
      type(mqc_vector)::tmpVec
      integer(kind=int64)::i,elems
      character(len=256)::vecString


      if(iPrint.ge.printThresh+1) then
        elems = min(3**iPrint,size(vector))
        tmpVec = vector
        write(6,'(1x,A)') 'Largest '//trim(num2char(elems))//' '//trim(titleString)//' coefficients'
        vecString = ''
        do i = 1, elems
          vecString = trim(vecString)//' '//trim(num2char(maxLoc(abs(tmpVec)),'I3'))//' = '//&
            trim(num2char(tmpVec%at(maxLoc(abs(tmpVec))),'F10.6'))
          if(i.ne.elems) vecString = trim(vecString)//', '
          call tmpVec%put(0.0,maxLoc(abs(tmpVec)))
          if(i.eq.10.or.i.eq.elems) then
            write(6,'(1x,A)') trim(vecString)
            vecString = ''
          endIf
        endDo
        write(iOut,'(A)') ''
      endIf
      if(iPrint.ge.printThresh) then
        elems = min(3**iPrint,size(vector))
        tmpVec = real(conjg(vector).ewp.vector)
        write(6,'(1x,A)') 'Largest '//trim(num2char(elems))//' '//trim(titleString)//' populations'
        vecString = ''
        do i = 1, elems
          vecString = trim(vecString)//' '//trim(num2char(maxLoc(abs(tmpVec)),'I3'))//' = '//&
            trim(num2char(tmpVec%at(maxLoc(abs(tmpVec))),'F10.6'))
          if(i.ne.elems) vecString = trim(vecString)//', '
          call tmpVec%put(0.0,maxLoc(abs(tmpVec)))
          if(i.eq.10.or.i.eq.elems) then
            write(6,'(1x,A)') trim(vecString)
            vecString = ''
          endIf
        endDo
        write(iOut,'(A)') ''
      endIf
!
      end subroutine print_coeffs_and_pops
!    
!
!     PROCEDURE get_rhos
!     
!     get_rhos is a subroutine that returns the atomic orbital transition density 
!     matrices of two (nonorthogonal) Slater determinants, as well as the dimension
!     of the overlap null space, the overlap and psuedo-overlap matrix elements. The
!     sign of nullSize gives the multiple of matrix elements accounting for antisymmetry
!     due to permutation of orbitals in the SVD.
!
      subroutine get_rhos(rho,nIJ,pnIJ,nullSize,mo_I,mo_J,overlap,nBasis,nAlpha,nBeta)

      implicit none

!     input/output variables
      type(mqc_scf_integral),dimension(2),intent(inOut)::rho
      type(mqc_scalar),intent(inOut)::nIJ,pnIJ
      integer(kind=int64)::nullSize
      type(mqc_scf_integral),intent(in)::mo_I,mo_J,overlap
      type(mqc_scalar),intent(in)::nBasis,nAlpha,nBeta

!     subroutine variables
      type(mqc_matrix)::mo_I_occ,mo_J_occ,mIJ,uMat,vMat,tmoI,tmoJ,rhoMat,tMat1,tMat2,&
        tMat3,tMat4
      type(mqc_vector)::sigmaMat
      logical::orthflag
      integer(kind=int64)::i

!
      mo_I_occ = mqc_integral_output_block(mo_I%orbitals('occupied',[int(nAlpha)],[int(nBeta)]),'full') 
      mo_J_occ = mqc_integral_output_block(mo_J%orbitals('occupied',[int(nAlpha)],[int(nBeta)]),'full') 
      
      mIJ = matmul(matmul(dagger(mo_I_occ),overlap%getBlock('full')),mo_J_occ)
      nIJ = mIJ%det()

      orthflag = .false.
      if((nIJ%abs()).lt.zero_thresh) then
        call mIJ%svd(EVals=sigmaMat,EUVecs=uMat,EVVecs=vMat)
        if(minval(abs(sigmaMat)).lt.zero_thresh) orthflag = .true.
      endIf

      if(orthflag) then
        pnIJ = 1.0
        nullSize = 0
        do i = 1,size(sigmaMat)
          if(sigmaMat%at(i).gt.zero_thresh) then
            pnIJ = pnIJ*sigmaMat%at(i)
          else
            nullSize = nullSize + 1
          endIf
        endDo
        nullSize = sign(1.0,real(uMat%det()))*sign(1.0,real(vMat%det()))*nullSize
      else
        pnIJ = nIJ
        nullSize = 0
      endIf

      if(orthflag) then
        call signCheckSVD(uMat,sigmaMat,vMat,int(nAlpha+nBeta)) 
        tmoI = matmul(dagger(uMat),dagger(mo_I_occ))
        tmoJ = matmul(mo_J_occ,dagger(vMat))
        call rhoMat%init(int(nBasis)*2,int(nBasis)*2)
        call pairDensity(sigmaMat,tmoI,tmoJ,rhoMat,int(nBasis),1)
      else
        rhoMat = matmul(matmul(mo_J_occ,MIJ%inv()),dagger(mo_I_occ))
      endif

      tMat1 = rhoMat%mat([1,int(nBasis)],[1,int(nBasis)])
      tMat2 = rhoMat%mat([int(nBasis)+1,int(nBasis)*2],[int(nBasis)+1,int(nBasis)*2])
      tMat3 = rhoMat%mat([int(nBasis)+1,int(nBasis)*2],[1,int(nBasis)])
      tMat4 = rhoMat%mat([1,int(nBasis)],[int(nBasis)+1,int(nBasis)*2])

      if(MQC_Matrix_Norm((tMat1-tMat2)).lt.1.0e-14.and.&
        tMat3%norm().lt.1.0e-14.and.tMat4%norm().lt.1.0e-14) then
        call mqc_integral_allocate(rho(1),'','space',tMat1)
      elseIf(tMat3%norm().lt.1.0e-14.and.tMat4%norm().lt.1.0e-14) then
        call mqc_integral_allocate(rho(1),'','spin',tMat1,tMat2)
      else
        call mqc_integral_allocate(rho(1),'','general',tMat1,tMat2,tMat3,tMat4)
      endIf

      if(orthflag) then
        call rhoMat%init(int(nbasis)*2,int(nbasis)*2)
        call pairDensity(sigmaMat,tmoI,tmoJ,rhoMat,int(nBasis),2)
        tmat1 = rhoMat%mat([1,int(nBasis)],[1,int(nBasis)])
        tmat2 = rhoMat%mat([int(nBasis)+1,int(nBasis)*2],[int(nBasis+1),int(nBasis)*2])
        tmat3 = rhoMat%mat([int(nBasis)+1,int(nBasis)*2],[1,int(nBasis)])
        tmat4 = rhoMat%mat([1,int(nBasis)],[int(nBasis)+1,int(nBasis)*2])
      endIf
      if(MQC_Matrix_Norm((tMat1-tMat2)).lt.1.0e-14.and.&
        tMat3%norm().lt.1.0e-14.and.tMat4%norm().lt.1.0e-14) then
        call mqc_integral_allocate(rho(2),'','space',tMat1)
      elseIf(tMat3%norm().lt.1.0e-14.and.tMat4%norm().lt.1.0e-14) then
        call mqc_integral_allocate(rho(2),'','spin',tMat1,tMat2)
      else
        call mqc_integral_allocate(rho(2),'','general',tMat1,tMat2,tMat4,tMat3)
      endIf
!
      end subroutine get_rhos
!
!
!     PROCEDURE signCheckSVD
!     
!     signCheckSVD is a subroutine that returns the phases of left and right
!     singular vectors such that they are located in the upper left quadrant of
!     the argand diagram. Optional inputs offset and dimen provide the ability
!     to update the phases of a subset of singular vectors.
!    
      subroutine signCheckSVD(U,S,Vdag,dimen,offset)
      
      implicit none
      type(mqc_matrix)::U
      type(mqc_matrix),optional::Vdag
      type(mqc_vector),optional::S
      type(mqc_scalar)::zero,det1,theta1,exp1
      integer(kind=int64)::dimen,start
      integer(kind=int64),optional::offset
      
      zero = 1.0e-13
      
      if(present(offset)) then
        start = offset
      else
        start = 1
      endIf
      
      do k = start, start+dimen-1
        det1 = U%at(mqc_vector_maxloc(abs(U%vat([0],[k]))),k)
        if(abs(real(det1)).le.zero.and.abs(aimag(det1)).gt.zero) then
          theta1 = (-1)*asin(aimag(det1)/(sqrt(real(det1)**2+aimag(det1)**2)))
        elseIf(abs(aimag(det1)).le.zero.and.abs(real(det1)).gt.zero) then
          theta1 = (-1)*acos(real(det1)/(sqrt(real(det1)**2+aimag(det1)**2)))
        elseIf(abs(real(det1)).gt.zero.and.abs(aimag(det1)).gt.zero) then
          theta1 = (-1)*atan2(det1)
        else
          theta1 = zero
        endIf
        exp1 = cmplx(cos(theta1),sin(theta1))
        if(abs(aimag(exp1)).le.zero) exp1 = real(exp1)
        call U%vput(exp1*U%vat([0],[k]),[0],[k])
        if(present(S).and.present(Vdag)) then
          if(S%at(k).gt.zero_thresh) then
            call Vdag%vput(dagger(exp1*dagger(Vdag%vat([k],[0]))),[k],[0])
          else
            det1 = conjg(Vdag%at(k,mqc_vector_maxloc(abs(Vdag%vat([k],[0])))))
            if(abs(real(det1)).le.zero.and.abs(aimag(det1)).gt.zero) then
              theta1 = (-1)*asin(aimag(det1)/(sqrt(real(det1)**2+aimag(det1)**2)))
            elseIf(abs(aimag(det1)).le.zero.and.abs(real(det1)).gt.zero) then
              theta1 = (-1)*acos(real(det1)/(sqrt(real(det1)**2+aimag(det1)**2)))
            elseIf(abs(real(det1)).gt.zero.and.abs(aimag(det1)).gt.zero) then
              theta1 = (-1)*atan2(det1)
            else
              theta1 = zero
            endIf
            exp1 = cmplx(cos(theta1),sin(theta1))
            if(abs(aimag(exp1)).le.zero) exp1 = real(exp1)
            call Vdag%vput(dagger(exp1*dagger(Vdag%vat([k],[0]))),[k],[0])
          endIf
        endIf
      endDo
      
      end subroutine signCheckSVD
!
!
!     PROCEDURE pairDens
!
!     pairDensity is a subroutine to calculate transition density matrix from MOs
!     given in the diagonal overlap basis between bra and ket orbitals, i.e. 
!     from transformed orbitals MOs:
!     C_J' = C_J.V
!     C_I' = U*.C_I*
!     Incoming rho must be zeroed.
!
      subroutine pairDensity(Svals,tMOI,tMOJ,total_rho,numBasis,DenNum)

      implicit none
      type(mqc_scalar)::hold
      type(mqc_vector)::Svals
      type(mqc_matrix)::tMOI,tMOJ,rho_pair,total_rho,i_pair,j_pair,t_mat
      real(kind=real64)::temp,det_m = 1.0
      integer::i,nvals,numBasis,nZero,DenNum,Flag
   
      nvals = Svals%size()
      nZero = 0
      Flag = 0
      do i = 1, nvals
        if(MQC_Scalar_Get_Intrinsic_Real(Svals%at(i)).lt.zero_thresh) nZero = nZero + 1
      endDo
   
      if(nZero.gt.2) return
   
      do i=1, nvals
        temp = Svals%at(i)
   
        i_pair = tMOI%mat([i,i],[1,numBasis*2])
        j_pair = tMOJ%mat([1,numBasis*2],[i,i])
   
        if(nZero.eq.0) then
          if(temp.gt.zero_thresh) then
            cycle
          else
            rho_pair = matmul(j_pair,i_pair)
          end if
        elseIf(nZero.eq.1) then
          if(DenNum.eq.1) then
            if(temp.lt.zero_thresh) then
              cycle
            else
              rho_pair = matmul(j_pair,i_pair)/Svals%at(i)
            end if
          elseIf(DenNum.eq.2) then
            if(temp.ge.zero_thresh) then
              cycle
            else
              rho_pair = matmul(j_pair,i_pair)
            end if
          else
            call mqc_error_i('Density number not recognised',6,'DenNum',DenNum)
          endIf
        elseIf(nZero.eq.2) then
          if(DenNum.eq.1) then
            if(temp.ge.zero_thresh) then
              cycle
            elseIf(flag.eq.0) then
              rho_pair = matmul(j_pair,i_pair)
              flag = 1
            elseIf(flag.eq.1) then
              cycle
            end if
          elseIf(DenNum.eq.2) then
            if(temp.ge.zero_thresh) then
              cycle
            elseIf(flag.eq.0) then
              flag = 1
              cycle
            elseIf(flag.eq.1) then
              rho_pair = matmul(j_pair,i_pair)
            end if
          else
            call mqc_error_i('Density number not recognised',6,'DenNum',DenNum)
          endIf
        endIf
        total_rho = total_rho + rho_pair
      end do
   
      end subroutine pairDensity
!
!
!     PROCEDURE get_hij
!
!     get_hij is a function that returns the value of a Hamiltonian matrix element between
!     two (nonorthogonal) determinants, given the transition density matrices and required
!     integrals.
!
      function get_hij(pnij,nullSize,rho,coreHam,eris,doDirectIn,fileInfo,fileName,&
          loopNumber,sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor,gauss_exe,doProcMem,&
          mem,ncpu,keep_intermediate_files) result(hIJ)

      implicit none
      type(mqc_scalar),intent(in)::pnij
      integer(kind=int64),intent(in)::nullSize
      type(mqc_scf_integral),dimension(2),intent(in)::rho
      type(mqc_scf_integral),intent(in)::coreHam
      type(mqc_twoeris),dimension(:),allocatable,intent(in)::eris
      logical,optional,intent(in)::doDirectIn
      type(mqc_gaussian_unformatted_matrix_file),optional,intent(inOut)::fileInfo
      character(len=*),optional,intent(in)::fileName
      integer(kind=int64),optional,intent(in)::loopNumber
      type(mqc_matrix),optional,intent(in)::sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor
      character(len=*),optional,intent(in)::gauss_exe,mem,ncpu
      logical,optional,intent(in)::doProcMem,keep_intermediate_files
      type(mqc_scalar)::hIJ

      logical::doDirect
      type(mqc_scf_integral)::fMat
      type(mqc_scalar)::core_con,gmat_con
      real(kind=real64)::zero=0.0,half=0.5

      if(present(doDirectIn)) then
        doDirect = doDirectIn
        if(.not.present(fileInfo).or.&
          .not.present(fileName).or.&
          .not.present(loopNumber).or.&
          .not.present(sh2AtMp).or.&
          .not.present(shlTyp).or.&
          .not.present(nPrmSh).or.&
          .not.present(prmExp).or.&
          .not.present(conCoef).or.&
          .not.present(conCoTwo).or.&
          .not.present(shCoor).or.&
          .not.present(gauss_exe).or.&
          .not.present(mem).or.&
          .not.present(ncpu).or.&
          .not.present(doProcMem).or.&
          .not.present(keep_intermediate_files)) &
          call mqc_error('Missing data required for direct Fock matrix computation')
      else
        doDirect = .false.
      endIf

      if(doDirect) then
        fMat = do_external_fock_build(fileInfo,fileName,loopNumber,rho(1),sh2AtMp,shlTyp,nPrmSh,prmExp,&
          conCoef,conCoTwo,shCoor,gauss_exe,doProcMem,mem,ncpu,keep_intermediate_files)
      else
        fMat = mqc_eri_integral_contraction(eris,rho(1))
        fMat = coreHam + fMat
      endIf
      core_con = mqc_scf_integral_contraction(rho(2),coreham)
      gmat_con = mqc_scf_integral_contraction(rho(2),fMat-coreHam)

      if(abs(nullSize).gt.0) then
        if(abs(nullSize).eq.1) then
          hij = pnij*(core_con+gmat_con)
        elseIf(abs(nullSize).eq.2) then
          hij = pnij*(gmat_con)
        elseIf(abs(nullSize).gt.2) then
          hij = zero
        else
          call MQC_error_I('Number of zero-overlap orbital pairs is zero, but then &
            &NIJ should not be zero',6,'nullSize',nullSize)
        endIf
        ! take care of antisymmetry
        hij = sign(1.0,nullSize)*hij
      else
        hij = pnij*(core_con+(half*gmat_con))
      endif

      end function get_hij
!
!
!     PROCEDURE get_dij
!
!     get_dij is a function that returns the value of a one-electron operator matrix element
!     between two (nonorthogonal) determinants, given the transition density matrices and 
!     required integrals.
!
      function get_dij(pnij,rho,oneElInt) result(dIJ)

      implicit none
      type(mqc_scalar),intent(in)::pnij
      type(mqc_scf_integral),intent(in)::rho
      type(mqc_scf_integral),intent(in)::oneElInt
      type(mqc_scalar)::dIJ

      dIJ = pnij*contraction(rho,oneElInt)

      end function get_dij
!
!
!     PROCEDURE do_external_fock_build
!
!     do_external_fock_build is a function that returns the Fock matrix between two 
!     (nonorthogonal) determinants, given a transition density matrix. The code requires
!     a modified version of Gaussian 16.
!
      function do_external_fock_build(temp_file,fileName,loopNumber,rho,sh2AtMp,shlTyp,nPrmSh,prmExp,&
          conCoef,conCoTwo,shCoor,gauss_exe,doProcMem,mem,ncpu,keep_intermediate_files) result(fMat)

      implicit none
      type(mqc_gaussian_unformatted_matrix_file),intent(inOut)::temp_file
      character(len=*),intent(in)::fileName
      integer(kind=int64),intent(in)::loopNumber
      type(mqc_scf_integral),intent(in)::rho
      type(mqc_matrix),intent(in)::sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor
      character(len=*),intent(in)::gauss_exe,mem,ncpu
      logical,intent(in)::doProcMem,keep_intermediate_files
      type(mqc_scf_integral)::fMat
      
      character(len=256)::newFileName

      newFileName = trim(fileName)//'-matrix-'
      call build_string_add_int(loopNumber,newFileName,20)
      newFileName = trim(newFileName) // '.mat'

      temp_file%icgu = 221
      call temp_file%create(newFileName)
      call temp_file%writeArray('SHELL TO ATOM MAP',sh2AtMp)
      call temp_file%writeArray('SHELL TYPES',shlTyp)
      call temp_file%writeArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
      call temp_file%writeArray('PRIMITIVE EXPONENTS',prmExp)
      call temp_file%writeArray('CONTRACTION COEFFICIENTS',conCoef)
      call temp_file%writeArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
      call temp_file%writeArray('COORDINATES OF EACH SHELL',shCoor)
      call temp_file%writeESTObj('density',est_integral=rho,override='general')
      call Close_MatF(temp_file%UnitNumber)

      call write_gau_fock_file(newFileName,mqc_matrix_test_symmetric(rho%getblock(),'hermitian'), &
        mqc_matrix_test_symmetric(rho%getblock(),'antihermitian'),doProcMem,mem,ncpu)
      write(6,'(A,A)') 'Processing input file ',newFileName
      call execute_command_line(gauss_exe//(trim(newFileName(1:(len(trim(newFileName))-4)))//'.com'))

      call temp_file%load(trim(newFileName(1:(len(trim(newFileName))-4)))//'-out.mat')
      call temp_file%getESTObj('fock',est_integral=Fmat)

      if(.not.keep_intermediate_files) then
        call execute_command_line('rm '//(trim(newFileName(1:(len(trim(newFileName))-4)))//'.com'))
        call execute_command_line('rm '//(trim(newFileName(1:(len(trim(newFileName))-4)))//'.log'))
        call execute_command_line('rm '//(trim(newFileName(1:(len(trim(newFileName))-4)))//'-out.mat'))
        call execute_command_line('rm '//newFileName)
      endIf

      end function do_external_fock_build
!
!
!     PROCEDURE write_gau_fock_file
!
      subroutine write_gau_fock_file(matrixFile,diag,ahrm,doProcMem,mem,ncpu)

      implicit none
      character(len=80)::matrixFile,tempFile,outFile,mem,ncpu
      logical::diag,ahrm,doProcMem
      integer::io_stat_number,unitno
    
      tempFile = matrixFile(1:(len(trim(matrixFile))-4))
      tempFile = trim(tempFile) // '.com'
      outFile = matrixFile(1:(len(trim(matrixFile))-4))
      outFile = trim(outFile) // '-out.mat'
    
      open(newunit=unitno,file=tempFile,status='replace',iostat=io_stat_number)
      write(unitno,"(A11,A80)")'%oldmatrix=',adjustl(matrixFile)
      if(doProcMem) write(unitno,"(A13,A6)")'%nprocshared=',adjustl(ncpu)
      if(doProcMem) write(unitno,"(A5,A6)")'%mem=',adjustl(mem)
      write(unitno,"(A9)")'#P nonstd'
!
!     nonstandard route section; uses iop 4/5=19,4/200=4
!
      write(unitno,*) '1/29=7,38=1,172=1/1;'
      write(unitno,*) '2/12=2,15=1,40=1/2;'
      write(unitno,*) '3/5=7,6=2,11=9,14=-4,25=1,30=1,67=1,116=7/1,2,3;'
      if(diag) then
        write(unitno,*) '4/5=19,6=2,200=5,33=2/1;'
        write(unitno,*) '6/33=3,200=1/20;'
      elseIf(ahrm) then
        write(unitno,*) '4/5=19,6=2,200=6,33=2/1;'
        write(unitno,*) '6/33=3,200=11/20;'
      else
        write(unitno,*) '4/5=19,6=2,200=4,33=2/1;'
        write(unitno,*) '6/33=3,200=11/20;'
      endIf
      write(unitno,*) '99/5=1,6=100000,9=1/99;'
      write(unitno,*) ''
      write(unitno,"(A)") outFile
      write(unitno,*) ''
    
      close(unit=unitno)
    
      end subroutine write_gau_fock_file
!
!
!*    NOTES
!*      Compilation of this program requires the MQC library (https://github.com/MQCPack/mqcPack)
!*      and the gauopen utility (http://gaussian.com/g16/gauopen.zip) and compilation with the
!*      f08 standard.
!*
!*      Compilation tested using: gfortran 9.2.0
!*
!*      Note that subroutine Wr_LCBuf needs modifying in gauopen/qcmatrix.F as follows:
!*        line 58:       LenBX = (LenBuf/(2*abs(NR)))*abs(NR)
!*        line 60:       Call Wr_CBuf(IU,NTot*abs(NR),LenBX,X)
!*
!*      Documentation generated with robodoc. To update documentation edit robodoc.rc to
!*      determine documentation output type and then run robodoc at the command line in the
!*      main directory.
!*
!*    AUTHORS
!*      Lee M. Thompson, University of Louisville, lee.thompson.1@lousiville.edu
!*
!*    COPYRIGHT
!*      (c) 2021-2020 by Lee M. Thompson distributed under terms of the MIT license.
!*
!****
!
 999  End Program tdCIS
