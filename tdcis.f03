      program tdCIS 
!
!     This program perfoms a real-time time-dependent truncated CI calculation.
!
!     L. M. Thompson, 2021
!
      use mqc_gaussian
      use iso_fortran_env
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
      character(len=:),allocatable::command,fileName,pLevel,subs_in,root_in,delta_t_in,maxsteps_in,&
        help_path
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      integer(kind=int64)::iOut=6,iPrint=1,ival,i,j,maxsteps=1000,elems,flag
      real(kind=real64)::delta_t=0.1  
      complex(kind=real64)::imag=(0.0,1.0)
      logical::UHF,newNum=.false.
      type(mqc_pscf_wavefunction)::wavefunction
      type(mqc_molecule_data)::moleculeInfo
      type(mqc_twoERIs)::eris,mo_ERIs
      type(mqc_scalar)::Vnn,final_energy
      type(mqc_determinant)::determinants
      type(mqc_scf_integral)::mo_core_ham,dipolex,dipoley,dipolez
      type(mqc_matrix)::CI_Hamiltonian,CI_Dipole_X,CI_Dipole_Y,CI_Dipole_Z,dipolex_eigvecs,&
        dipoley_eigvecs,dipolez_eigvecs,exp_CI_Hamiltonian,density
      type(mqc_vector)::subs,dipolex_eigvals,dipoley_eigvals,dipolez_eigvals,&
        nuclear_dipole,state_coeffs,td_ci_coeffs,Dmatrix
      character(len=10)::val
      integer(kind=int64),dimension(:),allocatable::isubs
      type(mqc_vector),dimension(2)::nRoot
      character(len=256)::vecString,root_string=''
!
!*    USAGE
!*      TDCIS [-f <matrix_file>] [--print-level <print_level>] [--sub-levels substitutions] 
!*        [--intial-cond weight_vector] [--time-step time_step] [--steps steps] [--help]
!*
!*    OPTIONS
!
      write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',73),NEW_LINE('a'), &
          '  Real-Time Truncated Configuration Interaction Time Evolution Calculator',NEW_LINE('a'), &
          ' ',repeat('*',73),NEW_LINE('a'), &
          NEW_LINE('a'),repeat(' ',30),'Version 21.10.1',NEW_LINE('a'),NEW_LINE('a'),&
          ' L. M. Thompson, Louisville KY, 2021.',NEW_LINE('a')
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
          call mqc_get_command_argument(i+1,pLevel)
          read(pLevel,'(I1)') iPrint
          j = i + 2
        elseIf(command.eq.'-sub-levels') then
!
!*      --sub-levels substitutions       Substitution levels permitted in truncated CI calculation. 
!*
          call mqc_get_command_argument(i+1,subs_in)
          j = i+2
        elseIf(command.eq.'--initial-cond') then
!
!*      --initial-cond weight_vector     Initial state for simulation. For a pure state, the  
!*                                       desired root can be input as a single integer. For a 
!*                                       mixed state, each non-zero weighted root should be 
!*                                       included along wih the specified weight (Default is
!*                                       a ground state population of 1.0).  
!*                                   
!*                                       Example: [1,0.5:2,0.5] specifies an initial state with 
!*                                       50% weight on the lowest root and 50% weight on the 
!*                                       fourth root.
!*
          call mqc_get_command_argument(i+1,root_in)
          root_string = root_in
          j=i+2
        elseif(command.eq.'--time-step') then
!
!*      --time-step time_step            Time step for simulation (default 0.1 au).
!*
          call mqc_get_command_argument(i+1,delta_t_in)
          read(delta_t_in,'(F12.6)') delta_t
          j = i + 2
        elseif(command.eq.'--steps') then
!
!*      --steps steps                    Number of steps to include in the calculation (default is 1000). 
!*
          call mqc_get_command_argument(i+1,maxsteps_in)
          read(maxsteps_in,'(I20)') maxsteps
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
        deallocate(command)
      endDo
!
!     Get determinant expansion substitution levels
!
      if(.not.allocated(subs_in)) then
        write(iOut,'(A)') 'Defaulting to CIS.'
        subs = [1]
      else
        do i = 1,len(subs_in)
          select case (subs_in(i:i))
          case('[','\(')
            if(i.eq.1) then
              newNum = .true.
              cycle
            else
              call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
            endIf
          case(']','\)')
            if(i.eq.len(subs_in).and..not.newNum.and.i.ne.1) then
              read(val,'(I10)') ival
              call subs%push(ival)
              newNum = .true.
              cycle
            else
              call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
            endIf
          case(',',' ')
            if(i.eq.1.or.i.eq.len(subs_in).or.i.eq.2.or.i.eq.len(subs_in)-1.or.newNum) then
              call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
            else
              read(val,'(I10)') ival
              call subs%push(ival)
              newNum = .true.
              cycle
            endIf
          case('0':'9')
            if(i.eq.1.or.i.eq.len(subs_in)) then
              call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
            else
              if(newNum) then
                val = subs_in(i:i)
              else
                val = trim(val)//subs_in(i:i)
              endIf
              newNum = .false.
              cycle
            endIf
          case default
            call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
          end select
        endDo
      endIf
!
      call fileInfo%load(filename)
      call fileInfo%getMolData(moleculeInfo)
      call fileInfo%getESTObj('wavefunction',wavefunction)
      if(iPrint.ge.4) call wavefunction%print(iOut,'all')
      call fileInfo%getESTObj('dipole x',est_integral=dipolex)
      if(iPrint.ge.4) call dipolex%print(6,'dipolex')
      call fileInfo%getESTObj('dipole y',est_integral=dipoley)
      if(iPrint.ge.4) call dipoley%print(6,'dipoley')
      call fileInfo%getESTObj('dipole z',est_integral=dipolez)
      if(iPrint.ge.4) call dipolez%print(6,'dipolez')
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
!     Generate Slater Determinants       
!
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
!     Transform dipole integrals to the length form in the CI basis
!     Need to check if nuclear dipole should be included here.
!
!      nuclear_dipole = matmul(transpose(moleculeInfo%Nuclear_Charges),&
!        transpose(moleculeInfo%Cartesian_Coordinates))
!      call nuclear_dipole%print(6,'nuclear_dipole')
!
!      dipolex = matmul(transpose(wavefunction%MO_Coefficients),&
!        matmul(dipolex,Wavefunction%MO_Coefficients))
!      call dipolex%print(6,'Dipole x in MO basis')
!      call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis,determinants, &
!        dipolex,UHF=UHF,CI_Hamiltonian=CI_Dipole_X,subs=isubs)
!      call CI_Dipole_X%print(6,'CI Dipole X')
!          
!      dipoley = matmul(transpose(wavefunction%MO_Coefficients),&
!        matmul(dipoley,Wavefunction%MO_Coefficients))
!      call dipoley%print(6,'Dipole y in MO basis')
!      call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis,determinants, &
!        dipoley,UHF=UHF,CI_Hamiltonian=CI_Dipole_Y,subs=isubs)
!      call CI_Dipole_Y%print(6,'CI Dipole Y')
!
!      dipolez = matmul(transpose(wavefunction%MO_Coefficients),&
!        matmul(dipolez,Wavefunction%MO_Coefficients))
!      call dipolez%print(6,'Dipole z in MO basis')
!      call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis,determinants, &
!        dipolez,UHF=UHF,CI_Hamiltonian=CI_Dipole_Z,subs=isubs)
!      call CI_Dipole_Z%print(6,'CI Dipole Z')
!
!     Diagonalize the dipole moment matrix to get transformation matrix U
!
!      if(iPrint.eq.1) write(iOut,*) 'Diagonalizing CI Dipole X'
!      call CI_Dipole_X%diag(dipolex_eigvals,dipolex_eigvecs)
!      if(iPrint.ge.1) then 
!        call dipolex_eigvecs%print(iOut,'CI Dipole X Eigenvectors')
!        call dipolex_eigvals%print(iOut,'CI Dipole X Eigenvalues')
!      endIf
!
!      if(iPrint.eq.1) write(iOut,*) 'Diagonalizing CI Dipole Y'
!      call CI_Dipole_Y%diag(dipoley_eigvals,dipoley_eigvecs)
!      if(iPrint.ge.1) then 
!        call dipoley_eigvecs%print(iOut,'CI Dipole Y Eigenvectors')
!        call dipoley_eigvals%print(iOut,'CI Dipole Y Eigenvalues')
!      endIf
!
!      if(iPrint.eq.1) write(iOut,*) 'Diagonalizing CI Dipole Z'
!      call CI_Dipole_Z%diag(dipolez_eigvals,dipolez_eigvecs)
!      if(iPrint.ge.1) then 
!        call dipolez_eigvecs%print(iOut,'CI Dipole Z Eigenvectors')
!        call dipolez_eigvals%print(iOut,'CI Dipole Z Eigenvalues')
!      endIf
!
!     Determine initial conditions in the basis of CI states.
!
      call root_vector_builder(root_string,size(wavefunction%pscf_energies),nRoot)
      call state_coeffs%init(size(wavefunction%pscf_energies),0.0)
      do i = 1, size(nroot(1))
        call state_coeffs%put(nRoot(2)%at(i),nRoot(1)%at(i))
      endDo
      state_coeffs = state_coeffs/state_coeffs%norm()
      vecString = '# 3) State weights:                        #'//NEW_LINE('a')
      do i = 1, size(nRoot(1))
        vecString = trim(vecString)//' #                     '//trim(num2char(nRoot(1)%at(i),'I3'))//' = '//&
          trim(num2char(state_coeffs%at(nRoot(1)%at(i)),'F7.4'))//'           #'
        if(i.ne.size(nRoot(1))) vecString = trim(vecString)//NEW_LINE('a')
      endDo
      write(iOut,'(1X,A)')               '############################################'
      write(iOut,'(1X,A)')               '#                                          #'
      write(iOut,'(1X,A)')               '#            INITIAL CONDITIONS            #'
      write(iOut,'(1X,A)')               '#            ------------------            #'
      write(iOut,'(1X,A)')               '#                                          #'
      write(iOut,'(1X,A,1X,I10,A)')      '# 1) Number of steps:',maxsteps,'           #'
      write(iOut,'(1X,A,5X,F12.6,1X,A)') '# 2) Time step:',delta_t,'au        #'
      write(iOut,'(1X,A)')               trim(vecString)
      write(iOut,'(1X,A)')               '#                                          #'
      write(iOut,'(1X,A)')               '############################################'
      write(iOut,'(1X,A)')               ''
!
!     Do loops over time-steps and propagate CI coefficients in each time step using eq. 16 of 
!     Krause et al. J. Phys. Chem., 2007, 127, 034107.
!
      write(iOut,'(1X,A)') 'STARTING TIME PROPAGATION'
      write(iOut,'(1X,A)')               ''
      do i = 1, maxsteps
        if(iPrint.ge.1) write(iOut,'(1X,A)') repeat('=',44)//NEW_LINE('A')
        if(i.eq.maxsteps) write(iOut,'(1X,A)') repeat(' ',18)//'FINAL RESULT'//NEW_LINE('A')//NEW_LINE('A')//&
        ' '//repeat('=',44)//NEW_LINE('A')
        if(iPrint.ge.1) write(iOut,'(1X,A,1X,F12.6,1X,A)') 'Time step:',delta_t*i,'au'//NEW_LINE('A')

        state_coeffs = exp((-1)*imag*delta_t*wavefunction%pscf_energies).ewp.state_coeffs
        if(iPrint.ge.2) call state_coeffs%print(6,'State coefficients',Blank_At_Bottom=.true.)

        call td_ci_coeffs%init(size(wavefunction%pscf_amplitudes,1))
        do j = 1, size(state_coeffs)
          td_ci_coeffs = td_ci_coeffs + state_coeffs%at(j)*wavefunction%pscf_amplitudes%vat([0],[j])
        endDo

        if(iPrint.ge.2) then
          elems = min(3**iPrint,size(td_ci_coeffs))
          Dmatrix = td_ci_coeffs
          write(6,'(1x,A)') 'Largest '//trim(num2char(elems))//' TD CI coefficients'
          vecString = ''
          do j = 1, elems
            vecString = trim(vecString)//' '//trim(num2char(maxLoc(abs(Dmatrix)),'I3'))//' = '//&
              trim(num2char(Dmatrix%at(maxLoc(abs(Dmatrix))),'F10.6'))
            if(j.ne.elems) vecString = trim(vecString)//', '
            call Dmatrix%put(0.0,maxLoc(abs(Dmatrix)))
            if(j.eq.10.or.j.eq.elems) then
              write(6,'(1x,A)') trim(vecString)
              vecString = ''
            endIf
          endDo
          write(iOut,'(A)') ''
        endIf

        final_energy = get_CI_Energy(CI_Hamiltonian,td_ci_coeffs) 
        if(iPrint.ge.1.or.i.eq.maxsteps) call final_energy%print(6,'Energy (au)',Blank_At_Bottom=.true.,&
          FormatStr='F14.8')

        density = get_one_gamma_matrix(iOut,iPrint,wavefunction%nBasis,determinants,td_ci_coeffs,subs=isubs)
        if(iPrint.ge.1.or.i.eq.maxsteps) call density%print(6,'MO Density matrix',Blank_At_Bottom=.true.)
      endDo
!
      contains

      function get_CI_Energy(CI_Hamiltonian,CI_vectors) result(energy)

      implicit none
      type(mqc_matrix),intent(in)::CI_Hamiltonian
      type(mqc_vector),intent(in)::CI_Vectors
      type(mqc_scalar)::energy

      energy = dot_product(matmul(dagger(CI_vectors),CI_Hamiltonian),CI_vectors)

      end function get_CI_energy
!
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

      end subroutine
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
