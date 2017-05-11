      subroutine ODE(t,Nx,Ny,Nc,beta,gamma,ro,cells)
      
      implicit none
      double precision t
      integer Nx, Ny,Nc,i,j

      double precision beta(Nx,Ny),gamma(Nx,Ny),ro(Nx,Ny)
      double precision beta0(Nx,Ny),gamma0(Nx,Ny),ro0(Nx,Ny)
      double precision betak1(Nx,Ny),gammak1(Nx,Ny),rok1(Nx,Ny)
      double precision betak2(Nx,Ny),gammak2(Nx,Ny),rok2(Nx,Ny)
      double precision betak3(Nx,Ny),gammak3(Nx,Ny),rok3(Nx,Ny)
      double precision betak4(Nx,Ny),gammak4(Nx,Ny),rok4(Nx,Ny)
      double precision betak5(Nx,Ny),gammak5(Nx,Ny),rok5(Nx,Ny)
      double precision b1(Nx,Ny),g1(Nx,Ny),r1(Nx,Ny)
      double precision vdx(Nx,Ny), vdy(Nx,Ny)
      double precision cells(Nc,2)
      integer grid(Nx,Ny)

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision tau,h,err
      integer iteration, index

      tau=0.d0
      h=dt
      call flow(t,Nx,Ny,vdx,vdy)

 13   call FromCellToGrid(Nx,Ny,Nc,cells,grid)

      do j=1,Ny
       do i=1,Nx
        beta0(i,j)=beta(i,j)
        gamma0(i,j)=gamma(i,j)
        ro0(i,j)=ro(i,j)
       enddo
      enddo
      iteration=0

      call rs(t,Nx,Ny,beta0,gamma0,ro0,betak1,gammak1,rok1,grid,vdx)
!     Runge-Kutta-Merson Method

 16   do j=1,Ny
       do i=1,Nx
        beta(i,j)=beta0(i,j)+h*betak1(i,j)/3
        gamma(i,j)=gamma0(i,j)+h*gammak1(i,j)/3
        ro(i,j)=ro0(i,j)+h*rok1(i,j)/3
       enddo
      enddo

      call rs(t+h/3,Nx,Ny,beta,gamma,ro,betak2,gammak2,rok2,grid,vdx)

      do j=1,Ny
       do i=1,Nx
        beta(i,j)=beta0(i,j)+h*(betak1(i,j)+betak2(i,j))/6
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+gammak2(i,j))/6
        ro(i,j)=ro0(i,j)+h*(rok1(i,j)+rok2(i,j))/6
       enddo
      enddo

      call rs(t+h/3,Nx,Ny,beta,gamma,ro,betak3,gammak3,rok3,grid,vdx)

      do j=1,Ny
       do i=1,Nx
        beta(i,j)=beta0(i,j)+h*(betak1(i,j)+3*betak3(i,j))/8
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+3*gammak3(i,j))/8
        ro(i,j)=ro0(i,j)+h*(rok1(i,j)+3*rok3(i,j))/8
       enddo
      enddo

      call rs(t+h/2,Nx,Ny,beta,gamma,ro,betak4,gammak4,rok4,grid,vdx)

       do j=1,Ny
       do i=1,Nx
        beta(i,j)=beta0(i,j)+h*(betak1(i,j)-3*betak3(i,j)
     .   +4*betak4(i,j))/2
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)-3*gammak3(i,j)
     .   +4*gammak4(i,j))/2
        ro(i,j)=ro0(i,j)+h*(rok1(i,j)-3*rok3(i,j)
     .   +4*rok4(i,j))/2
       enddo
      enddo

      call rs(t+h,Nx,Ny,beta,gamma,ro,betak5,gammak5,rok5,grid,vdx)

      do j=1,Ny
       do i=1,Nx
        beta(i,j)=beta0(i,j)+h*(betak1(i,j)+4*betak4(i,j)
     .   +betak5(i,j))/6
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+4*gammak4(i,j)
     .   +gammak5(i,j))/6
        ro(i,j)=ro0(i,j)+h*(rok1(i,j)+4*rok4(i,j)
     .   +rok5(i,j))/6
       enddo
      enddo

      do j=1,Ny
       do i=1,Nx
        b1(i,j)=beta(i,j)-h*(betak1(i,j)+betak2(i,j)+betak3(i,j)
     .   +betak4(i,j)+betak5(i,j))/5-beta0(i,j)
        g1(i,j)=gamma(i,j)-h*(gammak1(i,j)+gammak2(i,j)+gammak3(i,j)
     .   +gammak4(i,j)+gammak5(i,j))/5-gamma0(i,j)
        r1(i,j)=ro(i,j)-h*(rok1(i,j)+rok2(i,j)+rok2(i,j)+rok3(i,j)
     .   +rok4(i,j)+rok5(i,j))/5-ro0(i,j)
       enddo
      enddo

      err=0.d0
      index=0
      do j=1,Ny
       do i=1,Nx
       err=max(abs(b1(i,j)),abs(g1(i,j)),abs(r1(i,j)),err)
      if (beta(i,j) .lt. 0 .or. gamma(i,j) .lt. 0 .or. ro(i,j) .lt. 0)
     . then
      index=1

      endif
       enddo
      enddo

      if (err .gt. tol .or. index .eq. 1) then
      h=h/2
      iteration=iteration+1

        if (iteration .gt. 2) then
            write(6,*) 't =',t,' index =',index, 'iteration=',iteration
        endif
        if (iteration .gt. 40) then
            write(6,*) 'Emergency Exit'
            call EXIT(0)
        endif

      go to 16
      endif


      t=t+h
      tau=tau+h

      call UpdateDiscreteCells(h,Nx,Ny,Nc, gamma, ro, cells)

      h=dt

      if (tau + h .le. tout+tol*dt) then


       go to 13
      elseif(tau .lt. tout-tol*dt)then
         h = tout - tau

         go to 13
      endif


      return
      end




