!  Initialization of variables for the canonical-isokinetic algorithm:

   if(ntt.eq.9) then

      boltz=2.d0*boltz2
      erlxt=1.d0/gammai
      erlxt2=erlxt**2
      erlixt2=1.d0/erlxt2
      sinsh2=dtx/2.d0/nkija
      sinsh4=sinsh2/2.d0
      sinsh8=sinsh2/4.d0
      tktk=temp0*boltz
      davalev=0.d0
      clfs=0.d0

!     Read or generate the atomic velocities and chain-thermostat variables:

      if(irest.eq.1) then
         INQUIRE(FILE='vfreez.rst',EXIST=exstf1)
         INQUIRE(FILE='tfreez.rst',EXIST=exstf2)
      end if

!     Reading the atomic velocities and chain-thermostat variables:

#ifdef MPI
      if(master.and.irest.eq.1.and.exstf1.and.exstf2) then
         write(6,576) 'Using the existing files vfreez.rst and tfreez.rst', &
                'to read initial velocities and thermostat variables'
  576       format(/,a,/,a,/)

         open(578,file='vfreez.rst')
         iii = 0
         do j=1,natom
            read(578,*) svs(iii+1),svs(iii+2),svs(iii+3)
            iii = iii+3
         end do
         close(578)

         open(579,file='tfreez.rst')
         do j=1,natom
            read(579,*) s1s(j),s2s(j),s3s(j)
         end do
         close(579)

      end if

      if(irest.eq.1.and.exstf1.and.exstf2) then

         v(1:3*natom)=0.d0

         call mpi_allreduce(svs,v,3*natom,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)

         s1=0.d0
         s2=0.d0
         s3=0.d0

         call mpi_allreduce(s1s,s1,natom,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
         call mpi_allreduce(s2s,s2,natom,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
         call mpi_allreduce(s3s,s3,natom,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)

      end if

#else

      if(irest.eq.1.and.exstf1.and.exstf2) then
         write(6,576) 'Using the existing files vfreez.rst and tfreez.rst', &
                'to read initial velocities and thermostat variables'
  576       format(/,a,/,a,/)

         open(578,file='vfreez.rst')
         iii = 0
         do j=1,natom
            read(578,*) v(iii+1),v(iii+2),v(iii+3)
            iii = iii+3
         end do
         close(578)

         open(579,file='tfreez.rst')
         do j=1,natom
            read(579,*) s1(j),s2(j),s3(j)
         end do
         close(579)

      end if

#endif

      if(irest.eq.0.or..not.exstf1.or..not.exstf2) then

         if(master) write(6,577) 'Using random generation of initial &
                                 &velocities and thermostat variables'
  577       format(/,a,/)

!        At the very beginning of the simulations the velocities of atoms
!        are generated randomly to satisfy the kinetic energy constraint

         iseed(1)=0
         iseed(2)=1
         iseed(3)=2
         iseed(4)=3

         iii=0

         do j=1,natom

!           Random generation within the allowed interval

            call dlarnv(2,iseed,3,rndbim)

            v(iii+1)=rndbim(1)*dsqrt(tktk/amass(j))
            v(iii+2)=rndbim(2)*dsqrt(tktk/amass(j))
            v(iii+3)=rndbim(3)*dsqrt(tktk/amass(j))

!           Scaling to exactly satisfy the required initial temperature:

            if(tempi.ge.0.d0) then

               tkik=(amass(j)*(v(iii+1)**2+v(iii+2)**2+v(iii+3)**2))/3.d0/boltz*4.d0/3.d0

               v(iii+1)=v(iii+1)*dsqrt(tempi/tkik)
               v(iii+2)=v(iii+2)*dsqrt(tempi/tkik)
               v(iii+3)=v(iii+3)*dsqrt(tempi/tkik)
 
               tkik=(amass(j)*(v(iii+1)**2+v(iii+2)**2+v(iii+3)**2))/3.d0/boltz*4.d0/3.d0

            end if

!           Satisfying the kinetic energy constraint:

            hhin=amass(j)*(v(iii+1)**2+v(iii+2)**2+v(iii+3)**2)
            s1(j)=dsqrt(4.d0*(1.d0-hhin/(3.d0*tktk)))/erlxt

!           while the chain-thermostat variables are put to zero:

            s2(j)=0.d0
            s3(j)=0.d0
            iii=iii+3
         end do
      end if

!     Initialization for local conservation laws:

      do j=1,natom
         glcl(j)=0.d0
         etlc(j)=(s2(j)**2+s3(j)**2)*erlxt2/2.d0+glcl(j)
      end do
      iconfs=0

!        Initialization for distribution of dynamical variables:

      if(idistr.ne.0) then

         iconfd=0
         dargvs=5.d0/lep
         do jl=1,lep
            argvel(jl)=0.5d0*dargvs+(jl-1)*dargvs
            argvel(-jl)=-argvel(jl)
            distrs(jl)=0.d0
            distrs(-jl)=0.d0
            distrh1(jl)=0.d0
            distrh2(jl)=0.d0
            distrh2(-jl)=0.d0
            distrh3(jl)=0.d0
            distrh3(-jl)=0.d0
         end do

      end if

   end if
