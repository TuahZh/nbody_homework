program nbody
!use omp_lib
implicit none
integer num,i,k,thread_num,thread_count,z
integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
!real(16) time,x_diff,y_diff,z_diff,dist,dt,time_max
!real(16) pos(10,3),vel(10,3),force(10,3)
!real(16) mass(10),force_3(3)     !!变量声明位置固定？
real(16) time,x_diff,y_diff,z_diff,dist,dt,time_max,loo,er,cj,dt_max,mass_s, loo1
real(16) pos(99,3),vel(99,3),r(99,6),forces(99,3),pos_r(99,3),vel_r(99,3),pos_bc(3),vel_bc(3)
real(16) force(5,99,3),f(14,99,6)
!real,external :: fof
real(16) mass(99),force_3(3)     

! unit: (m,kg,s)
real*8, parameter :: G0=6.67D-11
real*8, parameter :: R0=1.5D11
real*8, parameter :: M0=2.D30
real*8, parameter :: F0=G0*M0*M0/R0/R0
real*8, parameter :: V0=sqrt(G0*M0/R0)
real*8, parameter :: T0=sqrt(R0*R0*R0/G0/M0)
open(unit=10,file='./2.dat')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print*,'number of stars(less than 99)'
!read*,num
!!print*,'put in the timescale and the timestep'
!!read*,time_max,dt
!i=1
!do while (i<=num)
!print*,'put in the position,velocity and mass of star number ',i
!read*,pos(i,1),pos(i,2),pos(i,3),vel(i,1),vel(i,2),vel(i,3),mass(i)
!i=i+1
!end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
time=0
dt=1.D-6
dt_max=1.
time_max=3.
num=3
mass(3)=1.D3
mass(2)=1.D0
mass(1)=1.D-3
!mass(4)=0.001
pos(3,1)=0.0D0
pos(3,2)=0.0D0
pos(3,3)=0.0D0
pos(2,1)=1.D0
pos(2,2)=0.0
pos(2,3)=0.0
pos(1,1)=pos(2,1)+1.D-3
pos(1,2)=0.0
pos(1,3)=0.0
!pos(4,1)=-84.0
!pos(4,2)=0.0
!pos(4,3)=0.0
vel(3,1)=0.0
vel(3,2)=0.0
vel(3,3)=0.0
vel(2,1)=0.0
vel(2,2)=3.1623D1
vel(2,3)=0.0
vel(1,1)=0.0
vel(1,2)=vel(2,2)+3.1623D1
vel(3,3)=0.0
!vel(4,1)=
!vel(4,3)=0.0
er=1.E-3
!!!!!calc barycentre!!!!!
pos_bc=0
vel_bc=0
do k=1,3
do i=1,num
pos_bc(k)=pos_bc(k)+mass(i)*pos(i,k)
vel_bc(k)=vel_bc(k)+vel(i,k)*mass(i)
end do
end do
mass_s=0
do i=1,num
mass_s=mass_s+mass(i)
end do
do k=1,3
pos_bc(k)=pos_bc(k)/mass_s
vel_bc(k)=vel_bc(k)/mass_s
end do
!!!!pos'=pos_bc+time*vel_bc
!num=99
!do i=1,num/2
!pos(i,1)=mod(mod(i,16),4)*10.0-0.13000000 
!pos(i,2)=(mod(i,16))/4*10.0-0.1300000000
!pos(i,3)=i/16*10.0-0.1300000000
!vel(i,1)=-(mod(i,16))/4-0.13000000
!vel(i,2)=i/16-0.13000000
!vel(i,3)=mod(mod(i,16),4)-0.13000000
!mass(i)=1    !!!const G equals 1!!!
!end do
!do i=num/2+1,num
!pos(i,1)=pos(i-num/2,1)+400.0000000
!pos(i,2)=pos(i-num/2,2)
!pos(i,3)=pos(i-num/2,3)
!vel(i,1)=vel(i-num/2,1)-500.00000000
!vel(i,2)=vel(i-num/2,2)
!vel(i,3)=vel(i-num/2,3)
!mass(i)=1    !!!const G equals 1!!!
!end do
do while (time<time_max) !!final time=1
!!!!!calc forces!!!!!
!计算k(n)
    do    i=1,num
!    write(10,'(g18.8,1x,g18.8,1x,g18.8,1x,g18.10)',advance='NO') pos(i,1)-pos_s(1),pos(i,2)-pos_s(2),&
!        &pos(i,3)-pos_s(3)
    write(10,'(i5,1x,g18.8,1x,g18.8,1x,g18.8,1x,g18.10)') i, pos(i,1)-pos_bc(1),pos(i,2)-pos_bc(2),&
        &pos(i,3)-pos_bc(3), time
    end do
f=0
	call fof(f,1,num,pos,vel,mass) !k(0)
  !  print*,f(1,1,1)
	do i=1,num
		do k=1,3
		pos_r(i,k)=2.0/27*dt*f(1,i,k)+pos(i,k)
		end do
  !  print*,f(1,1,1),f(1,1,4),pos_r(1,1)
		do k=4,6
		vel_r(i,k-3)=2.0/27*dt*f(1,i,k)+vel(i,k-3)
		end do
	!k(0)计算完毕，准备k(1)
	end do

	call fof(f,2,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(1.0/36*f(1,i,k)+1.0/12*f(2,i,k))+pos(i,k)
		end do
		do k=4,6
		!vel_r(i,k-3)=2.0/27*dt*f(1,i,k)+vel(i,k-3)
		vel_r(i,k-3)=dt*(1.0/36*f(1,i,k)+1.0/12*f(2,i,k))+vel(i,k-3)
		end do
	!准备k(2)
	end do
	
	call fof(f,3,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(1.0/24*f(1,i,k)+1.0/8*f(3,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(1.0/24*f(1,i,k)+1.0/8*f(2,i,k))+vel(i,k-3)
		end do
	!准备k(3)
	end do

	call fof(f,4,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(5.0/12*f(1,i,k)-25.0/16*f(3,i,k)+25.0/16*f(4,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(5.0/12*f(1,i,k)-25.0/16*f(3,i,k)+25.0/16*f(4,i,k))+vel(i,k-3)
		end do
	!准备k(4)
	end do

	call fof(f,5,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(1.0/20*f(1,i,k)+1.0/5*f(5,i,k)+1.0/4*f(4,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(1.0/20*f(1,i,k)+1.0/5*f(5,i,k)+1.0/4*f(4,i,k))+vel(i,k-3)
		end do
	!准备k(5)
	end do

	call fof(f,6,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(-25.0/108*f(1,i,k)-65.0/27*f(5,i,k)+125.0/108*f(4,i,k)+125.0/54*f(6,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(-25.0/108*f(1,i,k)-65.0/27*f(5,i,k)+125.0/108*f(4,i,k)+125.0/54*f(6,i,k))+vel(i,k-3)
		end do
	!准备k(6)
	end do

	call fof(f,7,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(31.0/300*f(1,i,k)+61.0/225*f(5,i,k)-2.0/9*f(6,i,k)+13.0/900*f(7,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(31.0/300*f(1,i,k)+61.0/225*f(5,i,k)-2.0/9*f(6,i,k)+13.0/900*f(7,i,k))+vel(i,k-3)
		end do
	!准备k(7)
	end do

	call fof(f,8,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(2*f(1,i,k)-53.0/6*f(4,i,k)+704.0/45*f(5,i,k)-107.0/9*f(6,i,k)+67.0/90*f(7,i,k)+3.0*f(8,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(2*f(1,i,k)-53.0/6*f(4,i,k)+704.0/45*f(5,i,k)-107.0/9*f(6,i,k)+67.0/90*f(7,i,k)+3.0*f(8,i,k))+vel(i,k-3)
		end do
	!准备k(8)
	end do

	call fof(f,9,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(-91.0/108*f(1,i,k)+23.0/108*f(4,i,k)-976.0/135*f(5,i,k)+311.0/54*f(6,i,k)&
            &-19.0/60*f(7,i,k)+17.0/6*f(8,i,k)-1.0/12*f(9,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(-91.0/108*f(1,i,k)+23.0/108*f(4,i,k)-976.0/135*f(5,i,k)+311.0/54*f(6,i,k)&
            &-19.0/60*f(7,i,k)+17.0/6*f(8,i,k)-1.0/12*f(9,i,k))+vel(i,k-3)
		end do
	!准备k(9)
	end do

	call fof(f,10,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(2383.0/4100*f(1,i,k)-341.0/164*f(4,i,k)+4496.0/1025*f(5,i,k)-301.0/82*f(6,i,k)&
            &+2133.0/4100*f(7,i,k)+45.0/82*f(8,i,k)+45.0/164*f(9,i,k)+18.0/41*f(10,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(2383.0/4100*f(1,i,k)-341.0/164*f(4,i,k)+4496.0/1025*f(5,i,k)-301.0/82*f(6,i,k)&
            &+2133.0/4100*f(7,i,k)+45.0/82*f(8,i,k)+45.0/164*f(9,i,k)+18.0/41*f(10,i,k))+vel(i,k-3)
		end do
	!准备k(10)
	end do

	call fof(f,11,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(3.0/205*f(1,i,k)-6.0/41*f(6,i,k)-3.0/205*f(7,i,k)-3.0/41*f(8,i,k)+3.0/41*f(9,i,k)+6.0/41*f(10,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(3.0/205*f(1,i,k)-6.0/41*f(6,i,k)-3.0/205*f(7,i,k)-3.0/41*f(8,i,k)+3.0/41*f(9,i,k)+6.0/41*f(10,i,k))+vel(i,k-3)
		end do
	!准备k(11)
	end do

	call fof(f,12,num,pos_r,vel_r,mass) 
	do i=1,num
		do k=1,3
		pos_r(i,k)=dt*(-1777.0/4100*f(1,i,k)-341.0/164*f(4,i,k)+4496.0/1025*f(5,i,k)-289.0/82*f(6,i,k)&
            &+2193.0/4100*f(7,i,k)+51.0/82*f(8,i,k)+33.0/164*f(9,i,k)+12.0/41*f(10,i,k)+f(12,i,k))+pos(i,k)
		end do
		do k=4,6
		vel_r(i,k-3)=dt*(-1777.0/4100*f(1,i,k)-341.0/164*f(4,i,k)+4496.0/1025*f(5,i,k)-289.0/82*f(6,i,k)&
            &+2193.0/4100*f(7,i,k)+51.0/82*f(8,i,k)+33.0/164*f(9,i,k)+12.0/41*f(10,i,k)+f(12,i,k))+vel(i,k-3)
		end do
	!准备k(12)
	end do
	
	call fof(f,13,num,pos_r,vel_r,mass)
	
	loo=1E-6
	do i=1,num
		do k=1,3
			loo1=abs(41.0/810*dt*(f(1,i,k)+f(11,i,k)-f(12,i,k)-f(13,i,k)))
			if (loo<loo1) then
				loo = loo1
			end if
		end do 
	end do
!	print *, f(1,3,1)
	if (loo > er ) then
	dt=dt/2.0
    print*,'loser!'
	time=time-dt
	else

    do i=1,num
        do k=i+1,num
            x_diff=pos_r(i,1)-pos_r(k,1)
            y_diff=pos_r(i,2)-pos_r(k,2)
            z_diff=pos_r(i,3)-pos_r(k,3)
!			if (abs(x_diff)<1.D-7 .or. abs(y_diff)<1.D-7 .or. abs(z_diff)<1.D-7) then
!			if (abs(x_diff)<1.D-30) then
!				print *, 'close', i, x_diff, y_diff, z_diff
!			end if
		end do
	end do

	dist=sqrt(x_diff**2+y_diff**2+z_diff**2)
	if (abs(dist)<1.D-4) then
		print *, 'close', dist
		time=time-dt
		dt=dt/2.0
		CYCLE
	end if

	do i=1,num
		do k=1,3
		pos(i,k)=dt*(34.0/105*f(6,i,k)+9.0/35*f(7,i,k)+9.0/35*f(8,i,k)+9.0/280*f(9,i,k)&
            &+9.0/280*f(10,i,k)+41.0/840*f(12,i,k)+41.0/840*f(13,i,k))+pos(i,k)
		end do
		do k=4,6
		vel(i,k-3)=dt*(34.0/105*f(6,i,k)+9.0/35*f(7,i,k)+9.0/35*f(8,i,k)+9.0/280*f(9,i,k)&
            &+9.0/280*f(10,i,k)+41.0/840*f(12,i,k)+41.0/840*f(13,i,k))+vel(i,k-3)
		end do
	end do
!     if (( loo < er*dt ).and.(dt<dt_max)) then
!                   dt=dt*2
!    end if
!   	write(10,'(f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f8.5)')&
   ! pos(19,1),pos(19,2),pos(19,3),pos(89,1),pos(89,2),pos(89,3),pos(29,1),pos(29,2),pos(29,3),time
	!write(10,'(f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f18.8,1x,f8.5)') pos(10,1),pos(10,2),pos(10,3),pos(2,1),pos(2,2),pos(2,3),time
    do i=1,3
    pos_bc(i)=pos_bc(i)+dt*vel_bc(i)
    end do
!    do    i=1,num-1
!    write(10,'(f28.8,1x,f28.8,1x,f28.8,1x)',advance='NO') pos(i,1)-pos_bc(1),pos(i,2)-pos_bc(2),&
!        &pos(i,3)-pos_bc(3)
!    end do
!    write(10,'(f28.8,1x,f28.8,1x,f28.8,1x,f28.10)') pos(num,1)-pos_bc(1),pos(num,2)-pos_bc(2),&
!        &pos(num,3)-pos_bc(3),time
!print*,pos(1,2)
	end if
time=time+dt
end do
!i=1
!do while (i<=num)
!print*,pos(i,1),pos(i,2),pos(i,3),vel(i,1),vel(i,2),vel(i,3),mass(i)
!i=i+1
!end do
end program


subroutine fof(f,rang,num,pos,vel,mass)
implicit none
integer i,k,thread_num,thread_count
integer, intent(in) :: rang, num  !!rang for k(rang)
integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
real(16) x_diff,y_diff,z_diff,dist
real(16), intent(in) :: pos(99,3),vel(99,3)
real(16)  forces(99,3)
real(16) force(5,99,3)
real(16), intent(inout) :: f(14,99,6)
real(16) mass(99),force_3(3)     
!print*,'num',num,'mass',mass(2)
!print*,'pos',pos(2,1),pos(2,2),pos(2,3)
!print*,'vel',vel(1,1),vel(1,2),vel(1,3)
!real fof  !10 for number of stars,14 for rkf range
!real(16),dimension(14,10,6) fof  !10 for number of stars,14 for rkf range
 !$omp parallel private(x_diff,y_diff,z_diff,dist,i,force_3),shared(pos,num,mass,vel,force,forces)
    thread_num=OMP_GET_THREAD_NUM()+1
    thread_count=OMP_GET_NUM_THREADS()
    !$omp do schedule(dynamic)
    do i=1,num
    do k=1,thread_count
    force(k,i,1)=0
    force(k,i,2)=0
    force(k,i,3)=0
    end do
    end do
    !$omp end do
    !$omp do schedule(guided)
    do i=1,num
        do k=i+1,num
            x_diff=pos(i,1)-pos(k,1)
            y_diff=pos(i,2)-pos(k,2)
            z_diff=pos(i,3)-pos(k,3)
!			if (abs(x_diff)<1.D-7 .or. abs(y_diff)<1.D-7 .or. abs(z_diff)<1.D-7) then
!			if (abs(x_diff)<1.D-30) then
!				print *, 'close', i, x_diff, y_diff, z_diff
!			end if
            dist=sqrt(x_diff**2+y_diff**2+z_diff**2)
            force_3(1)=-mass(i)*mass(k)/dist**3*x_diff
            force_3(2)=-mass(i)*mass(k)/dist**3*y_diff
            force_3(3)=-mass(i)*mass(k)/dist**3*z_diff
            force(thread_num,i,1)=force(thread_num,i,1)+force_3(1)
            force(thread_num,i,2)=force(thread_num,i,2)+force_3(2)
            force(thread_num,i,3)=force(thread_num,i,3)+force_3(3)
            force(thread_num,k,1)=force(thread_num,k,1)-force_3(1)
            force(thread_num,k,2)=force(thread_num,k,2)-force_3(2)
            force(thread_num,k,3)=force(thread_num,k,3)-force_3(3)
        end do
    end do
    !$omp end do
    !$omp do schedule(guided)
    do i=1,num
    forces(i,1)=0
    forces(i,2)=0
    forces(i,3)=0
        do k=1,thread_count  
        forces(i,1)=force(k,i,1)+forces(i,1)
        forces(i,2)=force(k,i,2)+forces(i,2)
        forces(i,3)=force(k,i,3)+forces(i,3)
        end do
    end do
    !$omp end do
    !$omp end parallel
	do i=1,num
	f(rang,i,1)=vel(i,1)
	f(rang,i,2)=vel(i,2)
	f(rang,i,3)=vel(i,3)
	f(rang,i,4)=forces(i,1)/mass(i)
	f(rang,i,5)=forces(i,2)/mass(i)
	f(rang,i,6)=forces(i,3)/mass(i)
	end do
!   print*,forces(1,1)
return 
end 


!!!FAILED!!!rkf方法需要知道每个步长内12个给定点的相应值，必须用到相互作用。每个h步长都要考虑N-1个多体的力函数。精度从h提高到h^7相应计算量增加20倍起。需要大量改动！
