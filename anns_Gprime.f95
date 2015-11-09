! Bhadeshia et al., ISIJ International, 1998, 38, 495
! http://www.msm.cam.ac.uk/map/nickel/programs/lattice.html
subroutine anns_LP_GPrime(T,XAl,XCr,alpha)

  implicit none

!  real(kind=8):: X_comp(1)
  double precision:: T,XAL,XCr,alpha

  Integer(kind=8), parameter:: MaxHN=20
  Integer(kind=8):: no_inputs, no_parameters, no_committee
  Integer(kind=8):: i,j,k,l,faulse
  real(kind=8), allocatable:: weigh(:,:,:),bias(:,:),extrema(:,:)
! for weigh(#input,#neuron,#model in the committee)
! for bias(#input, #model in the committee)

  real(kind=8), allocatable:: inputs(:)
  integer, allocatable:: hidnodes(:)
  real(kind=8):: A,B
  real(kind=8):: totalpha, parameters(100000)
  real(kind=8):: time1, time2
  real(kind=8):: w_main
  integer, parameter:: infile=19

  open(unit=infile,file='P_GPrime.dat')

  read(infile,*) no_inputs
  read(infile,*) no_committee

  allocate(inputs(no_inputs))
  allocate(extrema(no_inputs+1,2))
  allocate(hidnodes(no_committee))

  w_main=1d0
!  do i=1,1
!    if(i/=maj_ele) w_main=w_main-X_comp(i)
!  end do

  w_main=w_main-XAl-XCr

!======== databse ========
! 1 AL  2 CO  3 CR  4 MO  5 NI  6 TA  7 TI  8 W
!========         ========

! 1) Temperature Kelvin (it will be convert to deg. C)
  inputs(1)=T-273.15d0
! 2) Nickel  (at%)
  inputs(2)=w_main
! 3) Cobalt (at%)
  inputs(3)=0d0   !X_comp(2)
! 4) Chromium (at%)
  inputs(4)=XCr
! 5) Molybdenum (at%)
  inputs(5)=0d0   !X_comp(4)
! 6) Tungsten (at%)
  inputs(6)=0d0   !X_comp(8)
! 7) Aluminium (at%)
  inputs(7)=XAl
! 8) Titantium (at%)
  inputs(8)=0d0   !X_comp(7)
! 9) Niobium (at%)
  inputs(9)=0d0
!10) Tantalum (at%)
  inputs(10)=0d0  !X_comp(6)
!11) Hafnium (at%)
  inputs(11)=0d0
!12) Rhenium (at%)
  inputs(12)=0d0
!13) Vanadium (at%)
  inputs(13)=0d0
!14) Iron (at%)
  inputs(14)=0d0
!15) Gallium (at%)
  inputs(15)=0d0

! # of weighting is ( input + 1 (output) ) * hidden_nodes
  allocate(weigh(no_inputs+1,MaxHN,no_committee))
! # of the bias is hidden_nodes + 1 (output)
  allocate(bias(MaxHN+1,no_committee))

! Read the boundary of each input and one output
! extrema(i,1): minimum
! extrema(i,2): Maximum

  do i=1,no_inputs+1
    read(infile,*) extrema(i,1),extrema(i,2)
  end do

  do l=1, no_committee
    read(infile,*) hidnodes(l)

! ====================================== Read the weighting parameters
    no_parameters=(no_inputs+2)*hidnodes(l)+1

    parameters=0.0d0
    read(infile,*) (parameters(i),i=1,no_parameters)

!=== The format for Singh 1998

    k=1
    do j=1,hidnodes(l)
      bias(j,l)=parameters(k)
      k=k+1
      do i=1,no_inputs
        weigh(i,j,l)=parameters(k)
        k=k+1
      end do
    end do

    j=hidnodes(l)+1
    bias(j,l)=parameters(k)
    k=k+1

    i=no_inputs+1
    do j=1,hidnodes(l)
      weigh(i,j,l)=parameters(k)
      k=k+1
    end do

  end do

  do i=1,no_inputs
    inputs(i)=(inputs(i)-extrema(i,1))/(extrema(i,2)-extrema(i,1))-0.5d0
  end do

  totalpha=0.0d0
  do l=1, no_committee
    B=0.0d0
    do j=1,hidnodes(l)
      A=0.0d0
      do i=1,no_inputs
        A=A+inputs(i)*weigh(i,j,l)
      end do

      i=no_inputs+1
      B=B+tanh(A+bias(j,l))*weigh(i,j,l)
    end do

    i=no_inputs+1
    j=hidnodes(l)+1
    totalpha=totalpha+(B+bias(j,l)+0.5d0)*(extrema(i,2)-extrema(i,1))+extrema(i,1)
  end do

  alpha=totalpha/no_committee

  close(infile)

  return
end subroutine anns_LP_GPrime

