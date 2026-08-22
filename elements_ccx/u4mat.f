!
!     SPAX: mixed displacement/pressure tetrahedron for near-incompressible
!     phases -- the element CalculiX does not have and Abaqus calls C3D4H.
!
!     FORMULATION: MINI, P1(+bubble)/P1.
!
!     Why not the literal P1/P0 that Abaqus documents for C3D4H: on a linear
!     tet the strain is constant over the element, so an element-constant
!     pressure represents the divergence exactly, static condensation of that
!     pressure is algebraically the mean-dilatation B-bar operator, and B-bar
!     on a one-point tet is the identity. A faithful P1/P0 element in ccx is
!     therefore bit-for-bit plain C3D4 -- measured, in
!     patches/0002-bbar-mean-dilatation.patch. It cannot close any gap.
!
!     P1/P1 without enrichment violates the inf-sup condition (checkerboard
!     pressure). MINI restores it parameter-free by enriching the displacement
!     with a cubic bubble, which is element-internal and so condenses away
!     locally, leaving 4 nodes x 4 dof.
!
!     DOF LAYOUT: node-major, 4 dof per node -- 1,2,3 displacement, 4 pressure.
!     Declare in the deck with
!
!         *USER ELEMENT,TYPE=U4,NODES=4,INTEGRATIONPOINTS=15,MAXDOF=4
!
!     MAXDOF=4 raises mi(2) in allocation.f, which is what gives mastruct.c
!     room for the fourth nodal dof.
!
!     THE SYSTEM IS SYMMETRIC INDEFINITE (saddle point). Incomplete-Cholesky
!     PCG cannot be used on it; solve with SOLVER=SPOOLES.
!
!
!     u4mat: the whole element operator in one place. e_c3d_u4 turns it into a
!     stiffness matrix and resultsmech_u4 turns it into internal forces, so the
!     two cannot drift apart -- the failure mode that leaves
!     0002-bbar-mean-dilatation.patch with an equilibrium gap of 1.0.
!
!     Returns the bubble-condensed blocks aull (12x12), bl (4x12), cc (4x4),
!     and the pieces abi/alb/bb needed to reconstruct the bubble amplitude.
!
      subroutine u4mat(xl,um,xk,aull,bl,cc,abi,alb,bb,nelem)
!
      implicit none
!
      integer i,j,c,d,ii,jj,m,kk,nelem
      real*8 xl(3,4),um,xk,aull(12,12),bl(4,12),cc(4,4),abi(3,3),
     &     alb(12,3),bb(4,3),abb(3,3),shp(4,4),xsj,xi,et,ze,weight,w,
     &     gl(3,4),lv(4),bv,dbv(3),tmp(3,12),tmp2(3,4),det,gh,fac
!
      include "gauss.f"
!
      do i=1,12
        do j=1,12
          aull(i,j)=0.d0
        enddo
        do j=1,3
          alb(i,j)=0.d0
        enddo
      enddo
      do i=1,3
        do j=1,3
          abb(i,j)=0.d0
        enddo
      enddo
      do i=1,4
        do j=1,12
          bl(i,j)=0.d0
        enddo
        do j=1,3
          bb(i,j)=0.d0
        enddo
        do j=1,4
          cc(i,j)=0.d0
        enddo
      enddo
!
!     15-point rule: the bubble gradient is quadratic, so bubble-bubble terms
!     are quartic and the 1- and 4-point tet rules do not integrate them.
!
      do kk=1,15
        xi=gauss3d6(1,kk)
        et=gauss3d6(2,kk)
        ze=gauss3d6(3,kk)
        weight=weight3d6(kk)
!
        call shape4tet(xi,et,ze,xl,xsj,shp,3)
        w=weight*xsj
!
        do i=1,4
          lv(i)=shp(4,i)
          do j=1,3
            gl(j,i)=shp(j,i)
          enddo
        enddo
!
        bv=256.d0*lv(1)*lv(2)*lv(3)*lv(4)
        do m=1,3
          dbv(m)=256.d0*(gl(m,1)*lv(2)*lv(3)*lv(4)
     &                  +lv(1)*gl(m,2)*lv(3)*lv(4)
     &                  +lv(1)*lv(2)*gl(m,3)*lv(4)
     &                  +lv(1)*lv(2)*lv(3)*gl(m,4))
        enddo
!
!       deviatoric stiffness. For scalar shape functions with gradients g and
!       h and components c,d:
!         eps:eps = (delta_cd (g.h) + g_d h_c)/2 ,  div div = g_c h_d
!       so 2*mu*dev gives mu*(delta_cd (g.h) + g_d h_c) - 2*mu/3*g_c*h_d
!
        do i=1,4
          do c=1,3
            ii=3*(i-1)+c
            do j=1,4
              do d=1,3
                jj=3*(j-1)+d
                gh=gl(1,i)*gl(1,j)+gl(2,i)*gl(2,j)+gl(3,i)*gl(3,j)
                fac=0.d0
                if(c.eq.d) fac=gh
                aull(ii,jj)=aull(ii,jj)+w*(um*(fac+gl(d,i)*gl(c,j))
     &               -2.d0*um/3.d0*gl(c,i)*gl(d,j))
              enddo
            enddo
            do d=1,3
              gh=gl(1,i)*dbv(1)+gl(2,i)*dbv(2)+gl(3,i)*dbv(3)
              fac=0.d0
              if(c.eq.d) fac=gh
              alb(ii,d)=alb(ii,d)+w*(um*(fac+dbv(d)*gl(c,i))
     &             -2.d0*um/3.d0*gl(c,i)*dbv(d))
            enddo
          enddo
        enddo
        do c=1,3
          do d=1,3
            gh=dbv(1)*dbv(1)+dbv(2)*dbv(2)+dbv(3)*dbv(3)
            fac=0.d0
            if(c.eq.d) fac=gh
            abb(c,d)=abb(c,d)+w*(um*(fac+dbv(d)*dbv(c))
     &           -2.d0*um/3.d0*dbv(c)*dbv(d))
          enddo
        enddo
!
        do i=1,4
          do j=1,4
            do d=1,3
              bl(i,3*(j-1)+d)=bl(i,3*(j-1)+d)+w*lv(i)*gl(d,j)
            enddo
          enddo
          do d=1,3
            bb(i,d)=bb(i,d)+w*lv(i)*dbv(d)
          enddo
          do j=1,4
            cc(i,j)=cc(i,j)+w*lv(i)*lv(j)/xk
          enddo
        enddo
      enddo
!
!     static condensation of the bubble
!
      det=abb(1,1)*(abb(2,2)*abb(3,3)-abb(2,3)*abb(3,2))
     &   -abb(1,2)*(abb(2,1)*abb(3,3)-abb(2,3)*abb(3,1))
     &   +abb(1,3)*(abb(2,1)*abb(3,2)-abb(2,2)*abb(3,1))
      if(dabs(det).lt.1.d-30) then
        write(*,*) '*ERROR in u4mat: singular bubble block, element',
     &       nelem
        write(*,*) '       mu,K =',um,xk
        write(*,*) '       node1=',xl(1,1),xl(2,1),xl(3,1)
        write(*,*) '       node2=',xl(1,2),xl(2,2),xl(3,2)
        write(*,*) '       node3=',xl(1,3),xl(2,3),xl(3,3)
        write(*,*) '       node4=',xl(1,4),xl(2,4),xl(3,4)
        write(*,*) '       abb  =',abb(1,1),abb(2,2),abb(3,3)
        call exit(201)
      endif
      abi(1,1)=(abb(2,2)*abb(3,3)-abb(2,3)*abb(3,2))/det
      abi(1,2)=(abb(1,3)*abb(3,2)-abb(1,2)*abb(3,3))/det
      abi(1,3)=(abb(1,2)*abb(2,3)-abb(1,3)*abb(2,2))/det
      abi(2,1)=(abb(2,3)*abb(3,1)-abb(2,1)*abb(3,3))/det
      abi(2,2)=(abb(1,1)*abb(3,3)-abb(1,3)*abb(3,1))/det
      abi(2,3)=(abb(1,3)*abb(2,1)-abb(1,1)*abb(2,3))/det
      abi(3,1)=(abb(2,1)*abb(3,2)-abb(2,2)*abb(3,1))/det
      abi(3,2)=(abb(1,2)*abb(3,1)-abb(1,1)*abb(3,2))/det
      abi(3,3)=(abb(1,1)*abb(2,2)-abb(1,2)*abb(2,1))/det
!
      do i=1,3
        do j=1,12
          tmp(i,j)=abi(i,1)*alb(j,1)+abi(i,2)*alb(j,2)+abi(i,3)*alb(j,3)
        enddo
        do j=1,4
          tmp2(i,j)=abi(i,1)*bb(j,1)+abi(i,2)*bb(j,2)+abi(i,3)*bb(j,3)
        enddo
      enddo
!
      do i=1,12
        do j=1,12
          aull(i,j)=aull(i,j)-(alb(i,1)*tmp(1,j)+alb(i,2)*tmp(2,j)
     &         +alb(i,3)*tmp(3,j))
        enddo
      enddo
      do i=1,4
        do j=1,12
          bl(i,j)=bl(i,j)-(bb(i,1)*tmp(1,j)+bb(i,2)*tmp(2,j)
     &         +bb(i,3)*tmp(3,j))
        enddo
        do j=1,4
          cc(i,j)=cc(i,j)+(bb(i,1)*tmp2(1,j)+bb(i,2)*tmp2(2,j)
     &         +bb(i,3)*tmp2(3,j))
        enddo
      enddo
!
      return
      end
