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
      character(len=:),allocatable::command,fileName,help_path,nucFile
      character(len=256),dimension(:),allocatable::fileList
      character(len=256)::vecString,root_string='',sub_string='',field_string='',&
        pulseShape='rectangle',tcf_file='tcf',file_tmp,ci_string='oci',gauss_exe='$g16root/g16/g16',&
        mem='8GB',ncpu='2',alterations='',activeSpace='',command_message,nucUpdate='sudden'
      character(len=4),dimension(:),allocatable::atomList
      integer(kind=int64)::iOut=6,iPrint=1,iUnit,i,j,k,maxsteps,flag,stat_num,numFile=0,nullSize,&
        saveDen=0,nCore=0,nVirt=0,nNucSteps=0,nucStep,nAtoms,io_stat_number,exit_stat_number,&
        cmd_stat_number
      integer(kind=int64),dimension(:),allocatable::isubs,inactiveList,activeList,alphaList,betaList
      integer(kind=int64),allocatable::maxNucSteps
      real(kind=real64)::delta_t=0.1,field_size=0.005,simTime=100.0,t0=0.0,sigma=0.5,omega=10.0,&
        beta=3.0,newCharge,xCoord,yCoord,zCoord,current_time,nucStart=0.0
      real(kind=real64),allocatable::tcf_start,nucTime
      real(kind=real64),parameter::zero_thresh=1.00E-8
      complex(kind=real64)::imag=(0.0,1.0)
      logical::UHF,file_exists,doVelDip=.false.,doDirect=.false.,found,doProcMem=.false.,&
        keep_intermediate_files=.false.
      type(mqc_pscf_wavefunction)::wavefunction
      type(mqc_molecule_data)::moleculeInfo
      type(mqc_twoERIs),dimension(:),allocatable::eris
      type(mqc_twoERIs)::mo_ERIs
      type(mqc_scalar)::Vnn,final_energy,nIJ,pnIJ,hij,dij
      type(mqc_determinant)::determinants
      type(mqc_scf_integral)::mo_core_ham,density
      type(mqc_scf_integral),dimension(2)::rho
      type(mqc_scf_integral),dimension(3)::dipole,dipoleMO,veldipole
      type(mqc_scf_integral),dimension(:),allocatable::mo_list
      type(mqc_matrix)::CI_Hamiltonian,exp_CI_Hamiltonian,CI_Overlap,iden,sh2AtMp,shlTyp,nPrmSh,prmExp,&
        conCoef,conCoTwo,shCoor,saved_pscf_amplitudes,newGeom
      type(mqc_matrix),dimension(3)::CI_Dipole,dipole_eigvecs,Xmat,invXmat
      type(mqc_vector)::subs,nuclear_dipole,total_dipole,state_coeffs,td_ci_coeffs,field_vector,&
        td_field_vector,tcf_ci_epsilon,tcf,total_veldip,new_state_coeffs,chargeList
      type(mqc_vector),dimension(3)::dipole_eigvals
      type(mqc_vector),dimension(2)::nRoot
!
!*    USAGE
!*      TDCIS [-f <matrix_file>] [--print-level <print_level>] [--print-veldip] 
!*        [--sub-levels <substitutions>] [--core-orbitals <core-orbitals>] 
!*        [--virt-orbitals <virt-orbitals>] [--active-space <active_space>] [--alter <alter_list>] 
!*        [--ci-type <type_string>] [--direct] [--gauss-exe <gaussian_executable_string>] 
!*        [--do-proc-mem] [--mem <mem>] [--ncpu <ncpu>] [--keep-inters] 
!*        [--intial-state <weight_vector>] [--pulse-shape <pulse>] [--field-vector <field_vector>] 
!*        [--field-size <magnitude>] [--t0 <time>] [--omega <frequency>] [--sigma <width>] 
!*        [--beta <shift>] [--tcf-start <time>] [--time-step <time_step>] [--simulation-time <time>] 
!*        [--save <steps>] [--nuclear-step <update_file>] [--nuclear-time <time>] 
!*        [--nuclear_start <time>] [--maxNucSteps <steps>] [--nuclear-update <method>] [--help]
!*
!*    OPTIONS
!* 
!
!     Print program information.
!
      write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',73),NEW_LINE('a'), &
          '  Real-Time Truncated Configuration Interaction Time Evolution Calculator',NEW_LINE('a'), &
          ' ',repeat('*',73),NEW_LINE('a'), &
          NEW_LINE('a'),repeat(' ',30),'Version 22.1.1',NEW_LINE('a'),NEW_LINE('a'),&
          ' L. M. Thompson, A. M. Kinyua, Louisville KY, 2022.',NEW_LINE('a')
!
!     Parse input options.
!
!*   1. Input/output
!*
      j = 1
      do i=1,command_argument_count()
        if(i.ne.j) cycle
        call mqc_get_command_argument(i,command)
        if(command.eq.'-f') then
!
!*      -f matrix_file                   Input matrix file with initial set of molecular orbitals. 
!*                                       The first line contains the number of matrix files in the 
!*                                       input, and then on each line is a separate matrix file.
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
        elseif(command.eq.'--print-veldip') then
!
!*      --print-veldip                   Compute and print the velocity gauge dipole along with the
!*                                       length gauge dipole.
!*
          doVelDip=.true.
          j = i + 1
!
!*   2. Determinant expansion 
!*
        elseIf(command.eq.'--sub-levels') then
!
!*      --sub-levels substitutions       Substitution levels permitted in truncated CI calculation. 
!*                                       The default is all single substitutions.
!*
!*                                       Example: [1,2] specifies single and double substitutions.
!*
          call mqc_get_command_argument(i+1,command)
          sub_string = command
          j = i+2
        elseIf(command.eq.'--core-orbitals') then
!
!*      --core-orbitals core-orbitals    Number of occupied orbitals to exclude from the truncated 
!*                                       CI determinant expansion.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I3)') nCore
          j = i+2
        elseIf(command.eq.'--virt-orbitals') then
!
!*      --virt-orbitals virt-orbitals    Number of virtual orbitals to exclude from the truncated 
!*                                       CI determinant expansion.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I3)') nVirt
          j = i+2
        elseIf(command.eq.'--active-space') then
!
!*      --active-space active_space      Defines orbital space in which to construct initial basis
!*                                       determinant expansion through orbital swaps. For each input
!*                                       solution the number of electrons and orbitals in the active
!*                                       space should be specified.
!*
!*                                       Example: [4,4:3,6] specifies an expansion of four electrons
!*                                       in four orbitals in the first input orbitals and an
!*                                       expansion of three electrons in six orbitals in the second
!*                                       input orbitals.
!*
          call mqc_get_command_argument(i+1,command)
          activeSpace = command
          j = i+2
        elseIf(command.eq.'--alter') then
!
!*      --alter alter_list               Changes the order of orbitals in the input molecular
!*                                       orbitals. Alterations to each input orbitals should be
!*                                       colon separated and the two orbital numbers to be swapped
!*                                       should be separated by a if two alpha orbitals will be
!*                                       swapped, or b if two beta orbitals will be swapped.
!*                                       Different orbital swaps should be comma separated.
!*
!*                                       Example: [3a4,2b4:5a6] swaps alpha orbitals 3 and 4 and
!*                                       beta orbitals 2 and 4 in input orbitals 1, and swaps alpha
!*                                       orbitals 5 and 6 in input orbitals 2.
!*
          call mqc_get_command_argument(i+1,command)
          alterations = command
          j = i+2
        elseIf(command.eq.'--ci-type') then
!
!*      --ci-type type_string            Specifies the type of configuration interaction. Options
!*                                       are:
!*                                       1) oci (default)
!*                                          Perform orthogonal configuration interaction with 
!*                                          determinant expansion specified by sub-levels option.
!*                                          Only one matrix file should be specified in the input 
!*                                          file is expected and additional inputs will result in 
!*                                          an error.
!*                                       2) ocas
!*                                          Perform orthogonal complete active space determinant 
!*                                          expansion specified by active-space option. Only one 
!*                                          matrix file should be specified in the input file is 
!*                                          expected and additional inputs will result in an 
!*                                          error.
!*                                       3) noci
!*                                          Perform nonorthogonal configuration interaction with 
!*                                          determinant expansion specified by matrix files 
!*                                          listed in input file.
!*                                       4) hci (not yet implemented)
!*                                          Perform a hybrid configuration determinant expansion.
!*                                          That is, an orthogonal truncated configuration 
!*                                          interaction determinant expansion on molecular 
!*                                          orbitals in each matrix file in the input file, which 
!*                                          generally leads to mixed orthogonal and nonorthogonal 
!*                                          determinants. Substitution levels are specified by 
!*                                          sub-levels option.
!*                                       5) hcas (not yet implemented)
!*                                          Perform a hybrid complete active space determinant 
!*                                          expansion. That is, an orthogonal complete active 
!*                                          space determinant expansion on molecular orbitals in 
!*                                          each matrix file in the input file, which generally 
!*                                          leads to mixed orthogonal and nonorthogonal 
!*                                          determinants. Active spaces are specified by active-
!*                                          space option.
!*
          call mqc_get_command_argument(i+1,command)
          ci_string = command
          j = i+2
!
!*   3. Hamiltonian matrix construction
!*
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
!
!*   4. Initial conditions
!*
        elseIf(command.eq.'--initial-state') then
!
!*      --initial-state weight_vector    Initial state for simulation. For a pure state, the  
!*                                       desired root can be input as a single integer. For a 
!*                                       mixed state, each non-zero weighted root should be 
!*                                       included along wih the specified weight (Default is
!*                                       a ground state population of 1.0). Note that the input 
!*                                       vector is normalized regardless of input values.
!*                                   
!*                                       Example: [1,0.5:2,0.5] specifies an initial state with 
!*                                       50% weight on the lowest root and 50% weight on the 
!*                                       fourth root.
!*
          call mqc_get_command_argument(i+1,command)
          root_string = command
          j=i+2
!
!*   5. Applied field 
!*
        elseIf(command.eq.'--pulse-shape') then
!
!*      --pulse-shape pulse              Pulse shape used in simulation. Options are:
!*                                       1) rectangle (default)
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
!*                                       6) cos squared
!*                                          E(t) = E(0)*(cos((pi/(2*sigma))*(sigma-t+t0)))**2*
!*                                                   sin(omega*(t-t0))*H(t-t0)*H(sigma-t)
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
!
!*   6. Simulation parameters
!* 
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
        elseif(command.eq.'--save') then
!
!*      --save steps                     Save the density to a matrix file each time the specified interval
!*                                       of steps is reached, starting with the first step. A value of zero
!*                                       (default) means that no densities are saved.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I5)') saveDen
          j = i + 2
!
!*   7. Nuclear potential 
!*
        elseIf(command.eq.'--nuclear-step') then
!
!*      --nuclear-step file_name         File giving increment to nuclear coordinates and nuclear charge for
!*                                       each step in the nuclear potential. 
!*
!*                                       WARNING! The current implementation overwrites the input matrix files.
!*
          call mqc_get_command_argument(i+1,nucFile)
          j = i + 2
        elseif(command.eq.'--nuclear-time') then
!
!*      --nuclear-time nuclear_time      Time step for each change in nuclear geometry in au (default is equal
!*                                       to simulation time so no nuclear step is taken).
!*
          call mqc_get_command_argument(i+1,command)
          allocate(nucTime)
          read(command,'(F12.6)') nucTime
          j = i + 2
        elseif(command.eq.'--nuclear-start') then
!
!*      --nuclear-start nuclear_start    Time at which to start the nuclear steps so the first step is taken
!*                                       at nuclear_start+nuclear_time. Default is 0.0 au.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(F12.6)') nucStart
          j = i + 2
        elseif(command.eq.'--maxNucSteps') then
!
!*      --maxNucSteps nuclear_steps      Maximum number of nuclear steps. Can be used to specify the number of
!*                                       times the nuclear increment is applied.
!*
          call mqc_get_command_argument(i+1,command)
          allocate(maxNucSteps)
          read(command,'(I3)') maxNucSteps
          j = i + 2
        elseIf(command.eq.'--nuclear-update') then
!
!*      --nuclear-update nuclear_update  Method for updating the state weights after a nuclear step. Options 
!*                                       are:
!*                                       1) sudden (default)
!*                                          Reweight eigenstates so that wavefunction is the same after the 
!*                                          nuclear step, i.e. the wavefunction does not have time to respond
!*                                          to the perturbation.
!*                                       2) adiabatic
!*                                          Keep eigenstate weights the same after the nuclear step, i.e. the
!*                                          perturbation is slow enough that the wavefunction adapts to stay 
!*                                          in the same state.
!*
          call mqc_get_command_argument(i+1,command)
          nucUpdate = command
          call string_change_case(nucUpdate,'l')
          j = i + 2
!
!*   8. Help 
!*
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
!     Build the electric field pulse vector
!
      call field_vector_builder(field_string,field_vector)
      field_vector = field_vector*field_size
!
!     Output simulation conditions
!
      call root_vector_builder(root_string,nRoot)
      vecString = '# 4) State weights:'//repeat(' ',24)//'#'//NEW_LINE('a')
      do i = 1, size(nRoot(1))
        vecString = trim(vecString)//' #'//repeat(' ',26)//trim(num2char(nRoot(1)%at(i),'I3'))//' = '//&
          trim(num2char(nRoot(2)%at(i),'F7.4'))//repeat(' ',6)//'#'
        if(i.ne.size(nRoot(1))) vecString = trim(vecString)//NEW_LINE('a')
      endDo
      write(iOut,'(1X,A)')               '############################################'
      write(iOut,'(1X,A)')               '#                                          #'
      write(iOut,'(1X,A)')               '#            INITIAL CONDITIONS            #'
      write(iOut,'(1X,A)')               '#            ------------------            #'
      write(iOut,'(1X,A)')               '#                                          #'
      write(iOut,'(1X,A,1X,F15.6,A)')      '# 1) Simulation time:',simTime,' au   #'
      write(iOut,'(1X,A,1X,I15,A)')      '# 2) Number of steps:',int(simTime/delta_t),'      #'
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
!     Parse input file.
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
!     Update nuclear potential if necessary.
!
      if(allocated(nucTime)) then
        nNucSteps = (simTime-nucStart)/nucTime
        if(nNucSteps.lt.0) call mqc_error('Number of nuclear steps less than zero. Check &
          &nuclear step input parameters.')
      else
        allocate(nucTime)
        nucTime = simTime
      endIf
      
      if(allocated(maxNucSteps)) then
        if(maxNucSteps.lt.nNucSteps) nNucSteps = maxNucSteps
      endIf

      if(nNucSteps.ge.1) then
        if(.not.allocated(nucFile)) call mqc_error('No nuclear update file provided',6)
        open(newunit=iUnit,file=nucFile,status='old',iostat=io_stat_number)
        if(io_stat_number/=0) then
          call mqc_error('Error opening nuclear update file',6)
        endIf
        read(unit=iUnit, fmt='(i2)', iostat=io_stat_number) nAtoms
        if(nAtoms/=fileInfo%getVal('natoms',filename=fileList(1))) &
          call mqc_error_i('Import of nuclear change failed',6,'atoms on file',&
          nAtoms,'atoms on disk',MQC_Scalar_Get_Intrinsic_Integer(moleculeInfo%getNumAtoms()))
        call newGeom%init(3,nAtoms)
        call chargeList%init(nAtoms)
        read(unit=iUnit, fmt=*)
        do j = 1,nAtoms
          read(unit=iUnit, fmt=*, iostat=io_stat_number) newCharge, xCoord, yCoord, zCoord
          call chargeList%put(newCharge,j)
          call newGeom%put(xCoord/angPBohr,1,j)
          call newGeom%put(yCoord/angPBohr,2,j)
          call newGeom%put(zCoord/angPBohr,3,j)
        endDo
      endIf
      close(unit=iUnit)

      nucloop: do nucStep = 0, nNucSteps
        
        if(nucStep.ge.1) then
          write(iOut,'(1X,A)') repeat('*',50)//NEW_LINE('A')
          write(iOut,'(1X,A)') repeat(' ',12)//'UPDATING NUCLEAR POTENTIAL'//NEW_LINE('A')//NEW_LINE('A')//&
          ' '//repeat('*',50)//NEW_LINE('A')
          if(allocated(eris)) deallocate(eris)
          atomlist = MQC_Get_Nuclear_Symbols(moleculeInfo)
          do i = 1,size(fileList)
            call write_GauIn_file(iPrint,'nuclear_update_temp.com',fileList(i),i.gt.1,doProcMem,ncpu,mem,'chkbas',&
              .false.,(i.eq.1.and..not.doDirect),atomList,angPBohr*(moleculeInfo%Cartesian_Coordinates+NewGeom),&
              int(wavefunction%charge),int(wavefunction%multiplicity),moleculeInfo%Nuclear_Charges+chargeList)
          endDo
          call EXECUTE_COMMAND_LINE(gauss_exe//' nuclear_update_temp.com',exitstat=exit_stat_number,&
            cmdstat=cmd_stat_number,cmdmsg=command_message)
          if(exit_stat_number/=0) then
            call mqc_error_i('Error executing Gaussian calculation',6,'exit flag',exit_stat_number)
          elseIf(cmd_stat_number/=0) then
            call mqc_error_a('Could not execute Gaussian calculation',6,'cmdnsg',command_message)
          endIf
        endIf
!
!       Extract required data from matrix files.
!
        if((ci_string.eq.'oci'.or.ci_string.eq.'ocas').and.numFile.ne.1) then
          call mqc_error('Multiple matrix files input to requested orthogonal CI expansion.&
            & Use sub-levels option to provide expansion',iOut)
        elseIf(ci_string.eq.'oci'.or.ci_string.eq.'ocas'.or.ci_string.eq.'noci') then
          if(.not.allocated(mo_list)) allocate(mo_list(numFile))
          do i = 1, numFile
            call fileInfo%getESTObj('mo coefficients',est_integral=mo_list(i),filename=fileList(i))
            if(iPrint.ge.4) call mo_list(i)%print(iOut,'MO coefficients from matrix file '//trim(num2char(i)))
          endDo
        elseIf(ci_string.eq.'hci'.or.ci_string.eq.'hcas') then
          call mqc_error('Hybrid ci type is not yet implemented',iOut)
        else
          call mqc_error_a('Unrecognized CI type string provided',iOut,'ci_string',ci_string)
        endIf
!
        if(len_trim(alterations).ne.0) then
          call orbital_swapper(alterations,mo_list)
          if(iPrint.ge.4) then
            do i = 1, numFile
              call mo_list(i)%print(iOut,'Altered MO coefficients '//trim(num2char(i)))
            endDo
          endIf
        endIf
!
!        call MQC_Gaussian_SetDEBUG(.true.)
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
        if(doVelDip) then
          call fileInfo%getESTObj('vel dipole x',est_integral=veldipole(1))
          if(iPrint.ge.4) call veldipole(1)%print(6,'AO velocity dipole x integrals')
          call fileInfo%getESTObj('vel dipole y',est_integral=veldipole(2))
          if(iPrint.ge.4) call veldipole(2)%print(6,'AO velocity dipole y integrals')
          call fileInfo%getESTObj('vel dipole z',est_integral=veldipole(3))
          if(iPrint.ge.4) call veldipole(3)%print(6,'AO velocity dipole z integrals')
        endIf
!        call MQC_Gaussian_SetDEBUG(.false.)
!
        if(.not.doDirect) then
          if(ci_string.eq.'oci'.or.ci_string.eq.'ocas') then
            allocate(eris(1))
          else
            allocate(eris(3))
          endIf
          call fileInfo%get2ERIs('regular',eris(1),foundERI=found)
          if(found.and.ci_string.ne.'oci'.and.ci_string.ne.'ocas') then
            eris(2) = eris(1)
            eris(3) = eris(1)
            call mqc_twoeris_transform(eris(1),'raffenetti1')
            call mqc_twoeris_transform(eris(2),'raffenetti2')
            call mqc_twoeris_transform(eris(3),'raffenetti3')
          elseIf(ci_string.ne.'oci'.and.ci_string.ne.'ocas') then
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
        if(ci_string.eq.'oci'.or.ci_string.eq.'ocas') then
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
        if(doDirect.or.abs(saveDen).gt.0) then
          call fileInfo%getArray('SHELL TO ATOM MAP',sh2AtMp)
          call fileInfo%getArray('SHELL TYPES',shlTyp)
          call fileInfo%getArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
          call fileInfo%getArray('PRIMITIVE EXPONENTS',prmExp)
          call fileInfo%getArray('CONTRACTION COEFFICIENTS',conCoef)
          call fileInfo%getArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
          call fileInfo%getArray('COORDINATES OF EACH SHELL',shCoor)
        endIf
        if(doDirect.and.(ci_string.eq.'oci'.or.ci_string.eq.'ocas')) &
          call mqc_error('Direct matrix elements not possible with orthogonal CI')
!
!
!       Compute the nuclear-nuclear repulsion energy.
!
        call moleculeInfo%print(iOut)
        Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
        call Vnn%print(iOut,'Nuclear Repulsion Energy (au)',Blank_At_Bottom=.true.) 
!
!       Generate Slater determinants wavefunction expansion if orthogonal expansion requested. 
!
        if(ci_string.eq.'oci') then
          call substitution_builder(sub_string,subs)
          if(iPrint.ge.1) then
            write(iOut,'(1X,A)') 'Building Determinant Strings'
            call subs%print(6,'Permitted substitution levels',Blank_At_Bottom=.true.)
            write(iOut,'(1X,A,1X,I3,1X,A)') 'Excluding',nCore,'core orbitals'
            write(iOut,'(1X,A,1X,I3,1X,A)') 'Excluding',nVirt,'virtual orbitals'
          endIf
          if(nCore.gt.min(int(wavefunction%nAlpha),int(Wavefunction%nBeta))) &
            call mqc_error_i('Impossible number of core orbitals requested in truncated CI expansion',&
            6,'nCore',nCore,'wavefunction%nAlpha',int(wavefunction%nAlpha),'wavefunction%nBeta',&
            int(wavefunction%nBeta))
          if(nVirt.gt.min(int(wavefunction%nBasis-wavefunction%nAlpha),int(wavefunction%nBasis-Wavefunction%nBeta))) &
            call mqc_error_i('Impossible number of virtual orbitals excluded from truncated CI expansion',&
            6,'nVirt',nVirt,'wwavefunction%nBasis-avefunction%nAlpha',int(wavefunction%nBasis-wavefunction%nAlpha),&
            'wavefunction%nBasis-wavefunction%nBeta',int(wavefunction%nBasis-wavefunction%nBeta))
          isubs = [(i, i=1,int(maxval(subs)))]
          call trci_dets_string(iOut,iPrint,wavefunction%nBasis-nCore-nVirt,wavefunction%nAlpha-nCore, &
            Wavefunction%nBeta-nCore,isubs,determinants,nCore)
        elseIf(ci_string.eq.'ocas') then
          call parse_active_space(activeSpace,numFile,wavefunction%nBasis,wavefunction%nAlpha,&
            Wavefunction%nBeta,Wavefunction%nElectrons,activeList,inactiveList,alphaList,betaList)
          if(iPrint.ge.1) then
            write(iOut,'(1X,A)') 'Building Determinant Strings'
            call mqc_print(6,activeList,'Active orbitals')
            call mqc_print(6,inactiveList,'Core orbitals')
            call mqc_print(6,alphaList,'Active alpha electrons')
            call mqc_print(6,betaList,'Active beta electrons')
          endIf
          call gen_det_str(iOut,iPrint,activeList(1),alphaList(1),betaList(1),determinants,inactiveList(1))
          nCore = inactiveList(1)
          nVirt = wavefunction%nBasis-nCore-activeList(1)
        endIf
!
!       Transform one and two-electron integrals to MO basis (only if performing orthogonal CI).
!
        if(ci_string.ne.'noci') then
          if(iPrint.ge.1) write(iOut,'(1X,A)') 'Transforming MO integrals'//NEW_LINE('A')
          mo_core_ham = matmul(dagger(wavefunction%MO_Coefficients),matmul(wavefunction%core_Hamiltonian, &
              Wavefunction%MO_Coefficients))
          if(IPrint.ge.4) call mo_core_ham%print(iOut,'MO Basis Core Hamiltonian') 
          call twoERI_trans(iOut,iPrint,wavefunction%MO_Coefficients,ERIs(1),mo_ERIs)
        endIf
!
!       Compute nuclear dipole moment.
!
        nuclear_dipole = matmul(transpose(moleculeInfo%Nuclear_Charges),&
          transpose(moleculeInfo%Cartesian_Coordinates))
        if(iPrint.ge.1) call nuclear_dipole%print(6,'Nuclear dipole',Blank_At_Bottom=.true.)
        call total_dipole%init(3)
        if(doVelDip) call total_veldip%init(3)
!
!       Generate field-free Hamiltonian matrix and dipole matrices (in length form) in CI basis.
!
        if(ci_string.eq.'oci') then
          if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building orthogonal CI Hamiltonian matrix'
          call subs%unshift(0)
          isubs = subs
          call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis-nVirt,determinants, &
            mo_core_ham,mo_ERIs,UHF,CI_Hamiltonian,isubs)
          if(iPrint.ge.4) call CI_Hamiltonian%print(6,'CI Hamiltonian',Blank_At_Bottom=.true.)
          if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building orthogonal CI dipole matrices'
          do i = 1, 3
            dipoleMO(i) = matmul(dagger(wavefunction%MO_Coefficients),&
              matmul(dipole(i),Wavefunction%MO_Coefficients))
            if(iprint.ge.4) call dipoleMO(i)%print(6,'MO dipole integrals axis '//trim(num2char(i)),&
              Blank_At_Bottom=.true.)
            call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis-nVirt,determinants, &
              dipoleMO(i),UHF=UHF,CI_Hamiltonian=CI_Dipole(i),subs=isubs)
            if(iprint.ge.4) call CI_Dipole(i)%print(6,'SD dipole integrals axis '//trim(num2char(i)),&
              Blank_At_Bottom=.true.)
            call iden%identity(size(CI_Dipole(i),1),size(CI_Dipole(i),2),nuclear_dipole%at(i))
            CI_Dipole(i) = iden - CI_Dipole(i)
            if(iprint.ge.4) call CI_Dipole(i)%print(6,'SD dipole including nuclear term axis '//&
              trim(num2char(i)),Blank_At_Bottom=.true.)
          endDo
        elseIf(ci_string.eq.'ocas') then
          if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building orthogonal CI Hamiltonian matrix'
          call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis-nVirt,determinants, &
            mo_core_ham,mo_ERIs,UHF,CI_Hamiltonian)
          if(iPrint.ge.4) call CI_Hamiltonian%print(6,'CI Hamiltonian',Blank_At_Bottom=.true.)
          if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building orthogonal CI dipole matrices'
          do i = 1, 3
            dipoleMO(i) = matmul(dagger(wavefunction%MO_Coefficients),&
              matmul(dipole(i),Wavefunction%MO_Coefficients))
            if(iprint.ge.4) call dipoleMO(i)%print(6,'MO dipole integrals axis '//trim(num2char(i)),&
              Blank_At_Bottom=.true.)
            call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis-nVirt,determinants, &
              dipoleMO(i),UHF=UHF,CI_Hamiltonian=CI_Dipole(i))
            if(iprint.ge.4) call CI_Dipole(i)%print(6,'SD dipole integrals axis '//trim(num2char(i)),&
              Blank_At_Bottom=.true.)
            call iden%identity(size(CI_Dipole(i),1),size(CI_Dipole(i),2),nuclear_dipole%at(i))
            CI_Dipole(i) = iden - CI_Dipole(i)
            if(iprint.ge.4) call CI_Dipole(i)%print(6,'SD dipole including nuclear term axis '//&
              trim(num2char(i)),Blank_At_Bottom=.true.)
          endDo
        elseIf(ci_string.eq.'noci') then
          if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building nonorthogonal CI Hamiltonian, CI overlap and CI dipole matrices'
!         loop over pairs of Slater determinants and build transition density matrices. Get the H, N and Mu CI matrices
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
                dij = get_dij(pnij,nullSize,rho(2),dipole(k))
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
!       Diagonalize Hamiltonian
!
        if(iPrint.ge.1) write(iOut,'(1X,A)') 'Diagonalizing CI Hamiltonian'//NEW_LINE('A')
        if(ci_string.eq.'oci'.or.ci_string.eq.'ocas') then
          call CI_Hamiltonian%diag(wavefunction%pscf_energies,wavefunction%pscf_amplitudes)
        elseIf(ci_string.eq.'noci') then
          call CI_Hamiltonian%eigensys(CI_Overlap,wavefunction%pscf_energies,wavefunction%pscf_amplitudes)
        endIf

        if(iPrint.ge.4) call wavefunction%pscf_amplitudes%print(iOut,'CI Eigenvectors',Blank_At_Bottom=.true.)
        if(iPrint.ge.3) call wavefunction%pscf_energies%print(iOut,'CI Eigenvalues',Blank_At_Bottom=.true.)
!
!       Update state coefficients if nuclear potential has changed.
!
        if(nucStep.eq.0) then
          call state_coeffs%init(size(wavefunction%pscf_energies),0.0)
          do i = 1, size(nroot(1))
            call state_coeffs%put(nRoot(2)%at(i),nRoot(1)%at(i))
          endDo
        else
          select case (nucUpdate)
          case('sudden')
            call new_state_coeffs%init(size(state_coeffs),0.0)
            do i = 1, size(wavefunction%pscf_amplitudes, 1)
              do j = 1, size(wavefunction%pscf_amplitudes, 1)
                call new_state_coeffs%put(new_state_coeffs%at(j)+dot_product(dagger(saved_pscf_amplitudes%vat([0],[i])),&
                  wavefunction%pscf_amplitudes%vat([0],[j]))*state_coeffs%at(i),j)
              endDo
            endDo
            state_coeffs = new_state_coeffs
          case('adiabatic')
            continue
          case default 
            call mqc_error_a('Unrecognized nuclear step update argument given',6,'nucUpdate',nucUpdate)
          end select
        endIf
        state_coeffs = state_coeffs/state_coeffs%norm()
        if(nNucSteps.gt.0) saved_pscf_amplitudes = wavefunction%pscf_amplitudes
!
!       Diagonalize the dipole moment matrix to get transformation matrix U
!
        do i = 1, 3
          if(iPrint.eq.1) write(iOut,'(1X,A)') 'Diagonalizing SD dipole axis '//trim(num2char(i))//NEW_LINE('A')
          if (ci_string.eq.'oci'.or.ci_string.eq.'ocas') then
            call CI_Dipole(i)%diag(dipole_eigvals(i),dipole_eigvecs(i))
          elseIf(ci_string.eq.'noci') then
            call CI_Dipole(i)%eigensys(CI_Overlap,dipole_eigvals(i),dipole_eigvecs(i))
          endIf
          if(iPrint.ge.4) call dipole_eigvecs(i)%print(iOut,'CI Dipole Eigenvectors axis '//trim(num2char(i)),&
            Blank_At_Bottom=.true.)
          if(iprint.ge.3) call dipole_eigvals(i)%print(iOut,'CI Dipole Eigenvalues axis '//trim(num2char(i)),&
            Blank_At_Bottom=.true.)
          if (ci_string.eq.'oci'.or.ci_string.eq.'ocas') then
            Xmat(i) = matmul(dagger(dipole_eigvecs(i)),wavefunction%pscf_amplitudes)
            invXmat(i) = dagger(Xmat(i))
            if(iPrint.ge.2) call mqc_print(matmul(matmul(dagger(wavefunction%pscf_amplitudes),CI_Dipole(i)),&
              wavefunction%pscf_amplitudes),6,'State basis dipoles and transition dipoles on axis '//trim(num2char(i)),&
              Blank_At_Bottom=.true.)
          elseIf(ci_string.eq.'noci') then
            Xmat(i) = matmul(mqc_matrix_inverse(dipole_eigvecs(i)),wavefunction%pscf_amplitudes)
            invXmat(i) = mqc_matrix_inverse(Xmat(i))
            if(iPrint.ge.2) call mqc_print(matmul(matmul(wavefunction%pscf_amplitudes%inv(),CI_Dipole(i)),&
              wavefunction%pscf_amplitudes),6,'State basis dipoles and transition dipoles on axis '//trim(num2char(i)),&
              Blank_At_Bottom=.true.)
          endIf
        endDo
!
!       Determine initial conditions in the basis of CI states.
!
        if(nucStep.eq.nNucSteps.and.nucStep.ne.0) then
          maxsteps = nint((simtime-nucStart-nucTime*nNucSteps)/delta_t)
        elseIf(nucStep.eq.0) then
          maxsteps = (nucStart+nucTime)/delta_t
        else
          maxsteps = nucTime/delta_t
        endIf
!
!       Do loops over time-steps and propagate CI coefficients in each time step using eq. 16 of 
!       Krause et al. J. Phys. Chem., 2007, 127, 034107.
!
        write(iOut,'(1X,A)') 'STARTING TIME PROPAGATION'//NEW_LINE('A')
        do i = 0, maxsteps
          if(nucStep.ge.1) then
            current_time = nucStart
          else
            current_time = 0.0
          endIf
          current_time = current_time+(nucTime*nucStep)+(delta_t*i)
          write(iOut,'(1X,A)') repeat('=',44)//NEW_LINE('A')
          if(i.eq.maxsteps) write(iOut,'(1X,A)') repeat(' ',16)//'FINAL TIME STEP'//NEW_LINE('A')//NEW_LINE('A')//&
          ' '//repeat('=',44)//NEW_LINE('A')
          write(iOut,'(1X,A,1X,F12.6,1X,A)') 'Time step:',current_time,'au'//NEW_LINE('A')

          td_field_vector = get_field_vector(delta_t,(maxsteps*nucStep)+nucStep+i,field_vector,pulseShape,t0,omega,sigma,beta)
          call td_field_vector%print(6,'Applied field vector',Blank_At_Bottom=.true.,FormatStr='F16.10')

          if(i.gt.0) then
            state_coeffs = exp((-1)*imag*delta_t*wavefunction%pscf_energies).ewp.state_coeffs
            do j = 3, 1, -1
              state_coeffs = matmul(Xmat(j),state_coeffs)
              state_coeffs = exp(imag*delta_t*td_field_vector%at(j)*dipole_eigvals(j)).ewp.state_coeffs
              state_coeffs = matmul(invXmat(j),state_coeffs)
            endDo
          endIf
          call print_coeffs_and_pops(iOut,iPrint,1,state_coeffs,'TD State')

          call td_ci_coeffs%init(size(wavefunction%pscf_amplitudes,1))
          do j = 1, size(state_coeffs)
            td_ci_coeffs = td_ci_coeffs + state_coeffs%at(j)*wavefunction%pscf_amplitudes%vat([0],[j])
          endDo
          call print_coeffs_and_pops(iOut,iPrint,2,td_ci_coeffs/td_ci_coeffs%norm(),'TD CI')

          if(allocated(tcf_start)) then
            if(current_time.ge.tcf_start.and.current_time.lt.tcf_start+delta_t) tcf_ci_epsilon = td_ci_coeffs/td_ci_coeffs%norm()
            if(current_time.ge.tcf_start) call tcf%push(dot_product(dagger(tcf_ci_epsilon),td_ci_coeffs/td_ci_coeffs%norm()))
          endIf

          final_energy = get_CI_Energy(CI_Hamiltonian,td_ci_coeffs) 
          if(iPrint.ge.1.or.i.eq.maxsteps) call final_energy%print(6,'Energy (au)',Blank_At_Bottom=.true.,&
            FormatStr='F14.8')

          if(iPrint.ge.1.or.i.eq.maxsteps) call mqc_print(state_coeffs.outer.transpose(state_coeffs),6,'State Density matrix', &
            Blank_At_Bottom=.true.)
          if(ci_string.eq.'oci'.or.ci_string.eq.'ocas') then
            density = get_one_gamma_matrix(iOut,iPrint,wavefunction%nBasis-nVirt,determinants,td_ci_coeffs,UHF,&
              nOrbsIn=int(wavefunction%nBasis),subs=isubs)
            if(iPrint.ge.1.or.i.eq.maxsteps) call density%print(6,'MO Density matrix',Blank_At_Bottom=.true.)
            density = matmul(matmul(wavefunction%mo_coefficients,density),dagger(wavefunction%mo_coefficients))
            if(iPrint.ge.1.or.i.eq.maxsteps) call density%print(6,'AO Density matrix',Blank_At_Bottom=.true.)
          else
            density = get_noci_density(td_ci_coeffs,mo_list,wavefunction%overlap_matrix,wavefunction%nBasis,&
              wavefunction%nAlpha,wavefunction%nBeta) 
            if(iPrint.ge.1.or.i.eq.maxsteps) call mqc_print(matmul(matmul(matmul(dagger(mo_list(1)),&
              wavefunction%overlap_matrix),density),matmul(wavefunction%overlap_matrix,mo_list(1))),&
              6,'MO Density matrix (projected on MOs 1)',Blank_At_Bottom=.true.)
            if(iPrint.ge.1.or.i.eq.maxsteps) call density%print(6,'AO Density matrix',Blank_At_Bottom=.true.)
          endIf
          do j = 1, 3
            call total_dipole%put((-1)*contraction(density,dipole(j)) + nuclear_dipole%at(j),j)
            if(doVelDip) call total_veldip%put((-1)*contraction(density,veldipole(j)) + nuclear_dipole%at(j),j)
          endDo
          call total_dipole%print(6,'Total dipole',Blank_At_Bottom=.true.,FormatStr='F16.10')
          if(doVelDip) call total_veldip%print(6,'Total velocity gauge dipole',Blank_At_Bottom=.true.)

          if(abs(saveDen).gt.0) then
            if(mod(i,abs(saveDen)).eq.0) call outputCheckFile(fileInfo,&
              'MO-matrix-'//trim(num2char(nucStep))//'-'//trim(num2char(i)),&
              density,sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor,wavefunction)
          endIf

        endDo

        if(allocated(tcf_start)) then
          if(tcf_start.gt.simTime) then
            write(iOut,'(A)') 'Time correlation function requested to start after simulation time' 
          else
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
                call tcf%print(iUnit,'# time (time step: '//trim(num2char(delta_t))//'as, tcf start: '//&
                  trim(num2char(tcf_start))//')'//repeat(' ',10)//'time correlation function')
                write(iOut,'(A)') 'Saving data to file name '//trim(file_tmp)//'.dat'
              endIf
            endDo
          endIf
        endIf
      endDo nucLoop
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
      case('cos squared')
        if(dt*step.ge.t0.and.dt*step.lt.t0+2*sigma) then
          td_field = static_field*(cos((pi/(2*sigma))*(sigma-dt*step+t0)))**2*&
            sin(omega*(dt*step-t0))
        else
          td_field = [0.0,0.0,0.0]
        endIf
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
      type(mqc_vector),intent(out)::subs_out

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
      subroutine root_vector_builder(root_string_in,root_vector_out)

      implicit none

!     input/output variables
      character(len=*),intent(in)::root_string_in
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
      function get_dij(pnij,nullSize,rho,oneElInt) result(dIJ)

      implicit none
      type(mqc_scalar),intent(in)::pnij
      integer(kind=int64),intent(in)::nullSize
      type(mqc_scf_integral),intent(in)::rho
      type(mqc_scf_integral),intent(in)::oneElInt
      type(mqc_scalar)::dIJ

      real(kind=real64)::zero=0.0

      if(abs(nullSize).lt.2) then
        dIJ = pnij*contraction(rho,oneElInt)
      else
        dIJ = zero
      endIf

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
!     Write out a Gaussian input file with non-standard route used for computing Fock matrix 
!     with direct matrix elements. Code requires additional development Gaussian code.
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
!     PROCEDURE get_noci_density
!
!     get_noci_density is a function that returns the NOCI one-PDM.
!
      function get_noci_density(eigenvecs,mo_list,overlap,nBasis,nAlpha,nBeta) result(oneDM)

      implicit none

!     input/output variables
      type(mqc_vector),intent(in)::eigenvecs
      type(mqc_scf_integral),dimension(:),allocatable,intent(in)::mo_list
      type(mqc_scf_integral),intent(in)::overlap
      type(mqc_scalar),intent(in)::nBasis,nAlpha,nBeta
      type(mqc_scf_integral)::oneDM

!     subroutine variables
      type(mqc_scf_integral),dimension(2)::rho
      type(mqc_scalar)::nIJ,pnIJ,weight
      integer(kind=int64)::nullSize
      type(mqc_matrix)::oneDMmat,tMat1,tMat2,tMat3,tMat4
      integer::i,j
!
      call oneDMmat%init(int(nBasis)*2,int(nBasis)*2)
      do i = 1, numFile
        do j = 1, i
          call get_rhos(rho,nIJ,pnIJ,nullSize,mo_list(i),mo_list(j),overlap,nBasis,nAlpha,nBeta)
          if(nullSize.ge.2) cycle
          weight = eigenvecs%at(i)*pnIJ*conjg(eigenvecs%at(j))
          oneDMmat = oneDMmat + weight*rho(2)%getBlock('full')
          if(i.ne.j) oneDMmat = oneDMmat + conjg(weight)*dagger(rho(2)%getBlock('full'))
        endDo
      endDo
      tmat1 = oneDMmat%mat([1,int(nBasis)],[1,int(nBasis)])
      tmat2 = oneDMmat%mat([int(nBasis)+1,int(nBasis)*2],[int(nBasis+1),int(nBasis)*2])
      tmat3 = oneDMmat%mat([int(nBasis)+1,int(nBasis)*2],[1,int(nBasis)])
      tmat4 = oneDMmat%mat([1,int(nBasis)],[int(nBasis)+1,int(nBasis)*2])
      if(MQC_Matrix_Norm((tMat1-tMat2)).lt.1.0e-14.and.&
        tMat3%norm().lt.1.0e-14.and.tMat4%norm().lt.1.0e-14) then
        call mqc_integral_allocate(oneDM,'','space',tMat1)
      elseIf(tMat3%norm().lt.1.0e-14.and.tMat4%norm().lt.1.0e-14) then
        call mqc_integral_allocate(oneDM,'','spin',tMat1,tMat2)
      else
        call mqc_integral_allocate(oneDM,'','general',tMat1,tMat2,tMat4,tMat3)
      endIf
!
      end function get_noci_density
!
!
!     PROCEDURE outputCheckFile
!
!     outputMatFile is a subroutine that outputs a checkpoint file that can then be read into 
!     Gaussian or Gaussview. The purpose of this routine is to enable visualization of
!     the density matrix time evolution in Gaussview.
!
      subroutine outputCheckFile(temp_file,fileName,density,sh2AtMp,shlTyp,nPrmSh,prmExp,&
          conCoef,conCoTwo,shCoor,wavefunction) 

      implicit none

!     Input/output variables
      type(mqc_gaussian_unformatted_matrix_file),intent(inOut)::temp_file
      character(len=*),intent(in)::fileName
      type(mqc_scf_integral),intent(in)::density
      type(mqc_matrix),intent(in)::sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor
      class(mqc_wavefunction),intent(in)::wavefunction
!
!     subroutine variables
      character(len=80)::spinSymStr,compSpinStr
      
      temp_file%icgu = 111
      compSpinStr = 'real'
      if(density%type().eq.'space') then
        spinSymStr = 'space'
        if(MQC_Matrix_HaveComplex(density%getBlock('full'))) then
          temp_file%icgu = temp_file%icgu + 10
          compSpinStr = 'complex'
        endIf
      elseIf(density%type().eq.'spin') then
        spinSymStr = 'spin'
        temp_file%icgu = temp_file%icgu + 1
        if(MQC_Matrix_HaveComplex(density%getBlock('full'))) then
          temp_file%icgu = temp_file%icgu + 10
          compSpinStr = 'complex'
        endIf
      else
        spinSymStr = 'general'
        temp_file%icgu = temp_file%icgu + 100
        temp_file%icgu = temp_file%icgu + 10
        compSpinStr = 'complex'
      endIf
      if(temp_file%isOpen()) call Close_MatF(temp_file%UnitNumber)
      call temp_file%create(trim(fileName)//'.mat')
      call temp_file%writeArray('SHELL TO ATOM MAP',sh2AtMp)
      call temp_file%writeArray('SHELL TYPES',shlTyp)
      call temp_file%writeArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
      call temp_file%writeArray('PRIMITIVE EXPONENTS',prmExp)
      call temp_file%writeArray('CONTRACTION COEFFICIENTS',conCoef)
      call temp_file%writeArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
      call temp_file%writeArray('COORDINATES OF EACH SHELL',shCoor)
      call temp_file%writeESTObj('mo energies',est_eigenvalues=wavefunction%mo_energies,override=spinSymStr)
      call temp_file%writeESTObj('mo coefficients',est_integral=wavefunction%mo_coefficients,&
        override=spinSymStr,imagORide=compSpinStr)
      call temp_file%writeESTObj('density',est_integral=density,override=spinSymStr,imagORide=compSpinStr)
      call temp_file%writeESTObj('scf density',est_integral=wavefunction%scf_density_matrix,override=spinSymStr,&
        imagORide=compSpinStr)
      call Close_MatF(temp_file%UnitNumber)

      call execute_command_line("unfchk -matrix "//trim(fileName)//".mat "//trim(fileName)//".chk")
      call execute_command_line("formchk "//trim(fileName)//".chk")
      call execute_command_line("cubegen 1 density=scf "//trim(fileName)//".fchk "//trim(fileName)//".cube -3 h")

      end subroutine outputCheckFile 
!
!
!     PROCEDURE orbital_swapper
!
!     orbital_swapper is a subroutine that parses input and swaps molecular orbitals on the input matrix files.
!
      subroutine orbital_swapper(alterations,mo_list)

      implicit none

!     input/output variables
      type(mqc_scf_integral),dimension(:),allocatable,intent(inOut)::mo_list
      character(len=*),intent(in)::alterations

!     text parsing variables
      integer::i,leftNum,rightNum,solutionNum
      character(len=:),allocatable::alterString
      logical::lNum,rNum,lSet,rSet,betaPair,newSoln,hasNum
      character(len=80)::processString

      alterString = trim(alterations)
      if(len(alterString).ne.0) then
        !  Initialize orbital swap variables
        leftNum = 0
        rightNum = 0
        !  Flags for determining current state of the string processing.  lNum and rNum
        !  indicate whether the orbital pair is being assembled, lSet and rSet inidicate
        !  whether the left and right orbital numbers have been processed and stored.
        !  processString contains the characters that are being examined.
        lNum = .false.
        rNum = .false.
        lSet = .false.
        rSet = .false.
        betaPair = .false.
        newSoln = .false.
        hasNum = .false.
        processString = ''
        solutionNum = 1

        do i = 1, len(alterString)
          select case(alterString(i:i))
          case('[')
            if(i.eq.1) then
              lNum = .true.
              cycle
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            endif
          case("0":"9")
            if(i.lt.len(alterString)) then
              processString = trim(processString)//alterString(i:i)
              if(.not.hasNum) hasNum = .true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            endif
          case("a")
            if(lNum.and.hasNum) then
              read(processString,'(I3)') leftNum
              processString = ''
              betaPair = .false.
              hasNum = .false.
              lSet = .true.
              lNum = .false.
              rNum = .true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case("b")
            if(lNum.and.hasNum) then
              read(processString,'(I3)') leftNum
              processString = ''
              betaPair = .true.
              hasNum = .false.
              lSet = .true.
              lNum = .false.
              rNum=.true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case(",")
            if((lSet.and.rNum).and.hasNum) then
              read(processString,'(I3)') rightNum
              processString = ''
              hasNum = .false.
              rNum = .false.
              rSet = .true.
              lNum = .true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case(":")
            if((lSet.and.rNum).and.hasNum) then
              read(processString,'(I3)') rightNum
              processString = ''
              hasNum = .false.
              lNum = .true.
              rNum = .false.
              rSet = .true.
              newSoln = .true.
            else if((.not.lSet.and..not.rSet.and.lNum).and..not.(hasNum)) then
              newSoln = .true.
            else
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case("]")
            if(rNum.and.hasNum) then
              read(processString,'(I3)') rightNum
              rNum = .false.
              rSet = .true.
            else if(.not.(alterString(i-1:i-1).eq.':')) then
              call mqc_error("Malformed alteration string.  Expected format:  [#a#,#b#:#a#,...]")
            end if
          case default
            call mqc_error("Unrecognized character in alteration string input:  "//alterString(i:i))
          end select

          if(lSet.and.rSet) then
            if(.not.betaPair) then
              mo_list(solutionNum) = mo_list(solutionNum)%swap([leftNum,rightNum])
            else
              mo_list(solutionNum) = mo_list(solutionNum)%swap(betaOrbsIn=[leftNum,rightNum])
            endIf
            rSet = .false.
            lSet = .false.
          end if
          if(newSoln) then
            solutionNum = solutionNum + 1
            if(solutionNum.gt.size(mo_list)) call mqc_error_i('Orbital alterations requested on &
              &more determinants than available',6,'SolutionNum',SolutionNum,'size(mo_list)',&
              size(mo_list))
            newSoln = .false.
          end if
        end do
      end if

      end subroutine orbital_swapper
!    
!
!     PROCEDURE parse_active_space
!
!     parse_active_space is a subroutine that parses input for complete active space determinant 
!     expansions.
!
      subroutine parse_active_space(active_space,numFile,nBasis,nAlpha,nBeta,nElec,activeList,&
        inactiveList,alphaList,betaList)

      implicit none

!     input/output variables
      character(len=*),intent(in)::active_space
      integer,intent(in)::numFile
      type(mqc_scalar),intent(in)::nBasis,nElec,nAlpha,nBeta
      integer,dimension(:),allocatable,intent(out)::inactiveList,activeList,alphaList,betaList

!     text parsing variables
      integer::i,occNum,elecNum,solutionNum,maxActive
      character(len=:),allocatable::activeString
      logical::occSet,elecSet,oNum,eNum
      character(len=80)::processString
!
!     parse active space information
!
      allocate(inactiveList(numFile))
      allocate(activeList(numFile))
      allocate(alphaList(numFile))
      allocate(betaList(numFile))
      activeString = trim(active_space)
      if(len(activeString).eq.0) then
        do i = 1, numFile
          activeList(i) = nBasis
          inactiveList(i) = 0
          alphaList(i) = nAlpha
          betaList(i) = nBeta
          write(6,'(A)') ' Defaulting to substitutions over all orbitals in solution '//trim(num2char(i))
        endDo
      else
        !  Initialize occupied orbitals and electron count to zero
        occNum = 0
        elecNum = 0
        !  Flags for determining current state of the string processing.  oNum
        !  and eNum indicate whether the occupied orbital or active electron
        !  number for solution 'n' is currently being assembled.  occSet and
        !  elecSet indicate whether the occupied orbital number or active
        !  electron number have been fully assembled and stored.  processString
        !  contains the characters that correspond to occNum or elecNum.
        oNum = .false.
        eNum = .false.
        occSet = .false.
        elecSet = .false.
        processString = ''
        solutionNum = 1
        do i = 1, len(activeString)
          select case(activeString(i:i))
          case('[')
            if(i.eq.1) then
              oNum = .true.
              cycle
            else
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            end if
          case("0":"9")
            if(i.lt.len(activeString)) then
              processString = trim(processString)//activeString(i:i)
            else
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            endif
          case(",")
            if(.not.oNum) then
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            else
              read(processString,'(I3)') occNum
              processString = ''
              occSet = .true.
              oNum = .false.
              eNum = .true.
            end if
          case("]")
            if(.not.eNum) then
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            else
              read(processString,'(I3)') elecNum
              elecSet = .true.
              eNum = .false.
            end if
          case(":")
            if(occSet.and.eNum.and.(.not.elecSet)) then
              read(processString,'(I3)') elecNum
              elecSet = .true.
              oNum = .true.
              eNum = .false.
            else
              call mqc_error("Malformed active space string.  Expected format:  [#,#:#,#:...]")
            end if
          case default
            call mqc_error("Unrecognized character in active space input:  "//activeString(i:i))
          end select

          !  Algorithm has identified a valid pair of numbers for processing.
          if(occSet.and.elecSet) then

            maxActive = nBasis - (nElec/2) - mod(int(nElec),2) + (elecNum / 2) + mod(elecNum,2)

            !  Check users requested orbitals or electrons are sensible
            !  actually exist in the SCF results. I would hope that's not something
            !  the user would forget to check, but I'll err on the side of caution
            if(elecNum.gt.nElec) then
              call mqc_error("User has requested more active electrons than exist for solution "//&
                num2char(solutionNum))
            else if(occNum.gt.nBasis) then
              call mqc_error("User has requested more active orbitals than exist for solution "//&
                num2char(solutionNum))
            else if(occNum.gt.maxActive) then
              call mqc_error("User has requested too many active orbitals for solution "//&
                num2char(solutionNum))
            end if
   
            activeList(solutionNum) = occNum
            inactiveList(solutionNum) = (nelec-elecNum)/2 + mod(int(nelec-elecNum),2)
            alphaList(solutionNum) = nAlpha - inactiveList(solutionNum)
            betaList(solutionNum) = nBeta - inactiveList(solutionNum)
   
            occSet = .false.
            elecSet = .false.
            oNum = .true.
            eNum = .false.
            processString = ''
            solutionNum = solutionNum + 1
   
          end if
          if((solutionNum-1).gt.numFile) call mqc_error_i('Active space specifies more input solutions than&
            & present',6,'solutionNum',solutionNum-1,'numFile',numFile)
        end do
        if((solutionNum-1).lt.numFile) call mqc_error_i('Active space specifies fewer input solutions than&
          & present',6,'solutionNum',solutionNum-1,'numFile',numFile)
      end if

      end subroutine parse_active_space
!
!
!      PROCEDURE write_GauIn_file
!
!      write_GauIn_file is a subroutine that writes a Gaussian input file.
!
       subroutine write_GauIn_file(iPrint,filename,matFile,doAppend,doProcMem,nProc,mem,route_addition,saveMat,&
         doTwoERIs,atomList,cartesians,charge,multiplicity,nucCharge)
!
       implicit none
       integer,intent(in)::iPrint
       character(len=*),intent(in)::filename,matFile,nProc,mem
       character(len=*),optional,intent(in)::route_addition
       logical,intent(in)::doAppend,doProcMem,doTwoERIs
       logical,intent(in),optional::saveMat
       character(len=*),dimension(:),allocatable,intent(in),optional::atomList
       type(mqc_matrix),intent(in),optional::cartesians
       integer,intent(in),optional::charge,multiplicity
       type(mqc_vector),intent(in),optional::nucCharge
       integer::unitNumber,i
       character(len=256)::chkFile,matFileSave,geomcmd
!
       if(.not.doAppend) then
         open(newunit=unitNumber,file=filename,status='UNKNOWN')
       else
         open(newunit=unitNumber,file=filename,position='APPEND',status='OLD')
         write(unitNumber,'(A9)') '--link1--'
       endIf

       if(doProcMem) then
         write(unitnumber,'(A7,A)') '%nproc=',trim(nProc)
         write(unitnumber,'(A5,A)') '%mem=',trim(mem)
       endIf
       write(unitnumber,'(A11,A)') '%oldmatrix=',trim(matFile)
       chkFile = matFile(1:(len(trim(matFile))-4))
       chkFile = trim(chkFile) // '.chk'
       write(unitnumber,'(A5,A)') '%chk=',trim(chkFile)
       if(present(cartesians)) then
         geomcmd = ''
       else
         geomcmd = 'geom=allcheck'
       endIf
       if(doTwoERIs) then
         write(unitnumber,'(A,A,A)') '#P hf guess=read '//trim(geomcmd)//' scf=(conven,skip) int=noraf nosymm output=matrix'
       else
         write(unitnumber,'(A,A,A)') '#P hf guess=read '//trim(geomcmd)//' scf=(skip) int=noraf nosymm output=matrix'
       endIf
       if(present(route_addition)) write(unitnumber,'(A2,A)') '# ',trim(route_addition)
       write(unitnumber,'(A)') ''
       if(present(atomList).and.present(cartesians).and.present(charge).and.present(multiplicity)) then
         write(unitnumber,'(A)') 'Tile Card Required'
         write(unitnumber,'(A)') ''
         write(unitnumber,'(I3,1x,I2)') charge, multiplicity
         if(size(atomList).ne.size(cartesians,2)) call mqc_error_i('atom name and coordinate lists are not the&
           & same size in write_GauIn_file',6,'size(atomList)',size(atomList),'size(cartesians,2)',size(cartesians,2))
         do i = 1, size(atomList)
           if(present(nucCharge)) then
             write(unitnumber,'(1x,A,A,G0,A,1x,F15.8,1x,F15.8,1x,F15.8)') trim(atomList(i)),'(znuc=',&
               MQC_Scalar_Get_Intrinsic_Real(nucCharge%at(i)),')',MQC_Scalar_Get_Intrinsic_Real(cartesians%at(1,i)),&
               MQC_Scalar_Get_Intrinsic_Real(cartesians%at(2,i)),MQC_Scalar_Get_Intrinsic_Real(cartesians%at(3,i))
           else
             write(unitnumber,'(1x,A,F15.8,1x,F15.8,1x,F15.8)') atomList(i),MQC_Scalar_Get_Intrinsic_Real(cartesians%at(1,i)),&
               MQC_Scalar_Get_Intrinsic_Real(cartesians%at(2,i)),MQC_Scalar_Get_Intrinsic_Real(cartesians%at(3,i))
           endIf
         endDo
         write(unitnumber,'(A)') ''
       endIf
       if(present(saveMat).and.saveMat) then
         matFileSave = matFile(1:(len(trim(matFile))-4))
         matFileSave = trim(matFileSave) // '-save.mat'
         write(unitnumber,'(A)') matFileSave
       else
         write(unitnumber,'(A)') matFile
       endIf
       close(unit=unitNumber)

       end subroutine write_GauIn_file
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
!*      Adam M. Kinyua, University of Louisville, adam.kinyua@lousiville.edu
!*
!*    COPYRIGHT
!*      (c) 2021 by Lee M. Thompson distributed under terms of the MIT license.
!*
!****
!
 999  End Program tdCIS
