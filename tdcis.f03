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
      character(len=:),allocatable::command,fileName,help_path
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      integer(kind=int64)::iOut=6,iPrint=1,iUnit,i,j,maxsteps,flag,stat_num
      real(kind=real64)::delta_t=0.1,field_size=0.005,simTime=100.0,t0=0.0,sigma=0.5,omega=10.0,&
        beta=3.0
      real(kind=real64),allocatable::tcf_start
      complex(kind=real64)::imag=(0.0,1.0)
      logical::UHF,file_exists
      type(mqc_pscf_wavefunction)::wavefunction
      type(mqc_molecule_data)::moleculeInfo
      type(mqc_twoERIs)::eris,mo_ERIs
      type(mqc_scalar)::Vnn,final_energy
      type(mqc_determinant)::determinants
      type(mqc_scf_integral)::mo_core_ham,density
      type(mqc_scf_integral),dimension(3)::dipole,dipoleMO
      type(mqc_matrix)::CI_Hamiltonian,exp_CI_Hamiltonian,iden
      type(mqc_matrix),dimension(3)::CI_Dipole,dipole_eigvecs
      type(mqc_vector)::subs,nuclear_dipole,total_dipole,state_coeffs,td_ci_coeffs,field_vector,&
        td_field_vector,tcf_ci_epsilon,tcf
      type(mqc_vector),dimension(3)::dipole_eigvals
      integer(kind=int64),dimension(:),allocatable::isubs
      type(mqc_vector),dimension(2)::nRoot
      character(len=256)::vecString,root_string='',sub_string='',field_string='',&
        pulseShape='rectangle',tcf_file='tcf',file_tmp
!
!*    USAGE
!*      TDCIS [-f <matrix_file>] [--print-level <print_level>] [--sub-levels substitutions] 
!*        [--intial-state weight_vector] [--pulse-shape pulse] [--field-vector field_vector] 
!*        [--field-size magnitude] [--time-step time_step] [--simulation-time time] [--help]
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
!     Get required data from matrix file.
!
      call fileInfo%load(filename)
      call fileInfo%getMolData(moleculeInfo)
      call fileInfo%getESTObj('wavefunction',wavefunction)
      if(iPrint.ge.4) call wavefunction%print(iOut,'all')
      call fileInfo%getESTObj('dipole x',est_integral=dipole(1))
      if(iPrint.ge.4) call dipole(1)%print(6,'AO dipole x integrals')
      call fileInfo%getESTObj('dipole y',est_integral=dipole(2))
      if(iPrint.ge.4) call dipole(2)%print(6,'AO dipole y integrals')
      call fileInfo%getESTObj('dipole z',est_integral=dipole(3))
      if(iPrint.ge.4) call dipole(3)%print(6,'AO dipole z integrals')
      call fileInfo%get2ERIs('regular',eris)
      if(iPrint.ge.4) call eris%print(iOut,'AO 2ERIs')
!
      if(wavefunction%wf_type.eq.'U') then
        UHF = .true.
        if(iPrint.ge.3) write(iOut,'(1X,A)') 'Found UHF wavefunction'//NEW_LINE('A')
      elseIf(wavefunction%wf_type.eq.'R') then
        UHF = .False.
        if(iPrint.ge.3) write(iOut,'(1X,A)') 'Found RHF wavefunction'//NEW_LINE('A')
      else
        call mqc_error_A('Unsupported wavefunction type in fullci',iOut, &
          'Wavefunction%wf_type',wavefunction%wf_type)
      endIf 
!
      if (wavefunction%wf_complex) call mqc_error('Complex wavefunctions unsupported in fullci')
!
!     Compute the nuclear-nuclear repulsion energy.
!
      call moleculeInfo%print(iOut)
      Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
      call Vnn%print(iOut,'Nuclear Repulsion Energy (au)',Blank_At_Bottom=.true.) 
!
!     Generate Slater determinants wavefunction expansion      
!
      call substitution_builder(sub_string,subs)
      if(iPrint.ge.1) then
        write(iOut,'(1X,A)') 'Building Determinant Strings'
        call subs%print(6,'Permitted substitution levels',Blank_At_Bottom=.true.)
      endIf
      isubs = [(i, i=1,int(maxval(subs)))]
      call trci_dets_string(iOut,iPrint,wavefunction%nBasis,wavefunction%nAlpha, &
        Wavefunction%nBeta,isubs,determinants)
!
!     Transform one and two-electron integrals to MO basis
!
      if(iPrint.ge.1) write(iOut,'(1X,A)') 'Transforming MO integrals'//NEW_LINE('A')
      mo_core_ham = matmul(transpose(wavefunction%MO_Coefficients),matmul(wavefunction%core_Hamiltonian, &
          Wavefunction%MO_Coefficients))
      if(IPrint.ge.4) call mo_core_ham%print(iOut,'MO Basis Core Hamiltonian') 
      call twoERI_trans(iOut,iPrint,wavefunction%MO_Coefficients,ERIs,mo_ERIs)
!
!     Generate field-free Hamiltonian matrix
!
      if(iPrint.ge.1) write(iOut,'(1X,A)') 'Building CI Hamiltonian'
      call subs%unshift(0)
      isubs = subs
      call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis,determinants, &
        mo_core_ham,mo_ERIs,UHF,CI_Hamiltonian,isubs)
      if(iPrint.ge.4) call CI_Hamiltonian%print(6,'CI Hamiltonian')
!
!     Diagonalize Hamiltonian
!
      if(iPrint.ge.1) write(iOut,'(1X,A)') 'Diagonalizing CI Hamiltonian'//NEW_LINE('A')
      call CI_Hamiltonian%diag(wavefunction%pscf_energies,wavefunction%pscf_amplitudes)
      if(iPrint.ge.4) call wavefunction%pscf_amplitudes%print(iOut,'CI Eigenvectors')
      if(iPrint.ge.3) call wavefunction%pscf_energies%print(iOut,'CI Eigenvalues')
!
!     Transform dipole integrals to the length form in the CI basis.
!
      nuclear_dipole = matmul(transpose(moleculeInfo%Nuclear_Charges),&
        transpose(moleculeInfo%Cartesian_Coordinates))
      if(iPrint.ge.1) call nuclear_dipole%print(6,'Nuclear dipole',Blank_At_Bottom=.true.)
      call total_dipole%init(3)

      do i = 1, 3
        dipoleMO(i) = matmul(dagger(wavefunction%MO_Coefficients),&
          matmul(dipole(i),Wavefunction%MO_Coefficients))
        if(iprint.ge.4) call dipoleMO(i)%print(6,'MO dipole integrals axis '//trim(num2char(i)))
        call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis,determinants, &
          dipoleMO(i),UHF=UHF,CI_Hamiltonian=CI_Dipole(i),subs=isubs)
        if(iprint.ge.4) call CI_Dipole(i)%print(6,'SD dipole integrals axis '//trim(num2char(i)))
        call iden%identity(size(CI_Dipole(i),1),size(CI_Dipole(i),2),nuclear_dipole%at(i))
        CI_Dipole(i) = iden - CI_Dipole(i)
        if(iprint.ge.4) call CI_Dipole(i)%print(6,'SD dipole including nuclear term axis '//trim(num2char(i)))
      endDo
!
!     Diagonalize the dipole moment matrix to get transformation matrix U
!
      do i = 1, 3
        if(iPrint.eq.1) write(iOut,'(1X,A)') 'Diagonalizing SD dipole axis '//trim(num2char(i))//NEW_LINE('A')
        call CI_Dipole(i)%diag(dipole_eigvals(i),dipole_eigvecs(i))
        if(iPrint.ge.4) call dipole_eigvecs(i)%print(iOut,'CI Dipole Eigenvectors axis'//trim(num2char(i)))
        if(iprint.ge.3) call dipole_eigvals(i)%print(iOut,'CI Dipole Eigenvalues axis '//trim(num2char(i)))
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
          state_coeffs = matmul(matmul(dagger(dipole_eigvecs(j)),wavefunction%pscf_amplitudes),state_coeffs)
          state_coeffs = exp(imag*delta_t*td_field_vector%at(j)*dipole_eigvals(j)).ewp.state_coeffs
          state_coeffs = matmul(matmul(dagger(wavefunction%pscf_amplitudes),dipole_eigvecs(j)),state_coeffs)
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
