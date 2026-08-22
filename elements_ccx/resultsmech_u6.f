!
!     SPAX: internal forces for the U6 volumetric patch.
!
!     f = K V_a b (b . u), the same rank-one operator e_c3d_u6 assembles, and
!     obtained through the same u6patch call so the two cannot disagree. The
!     patch carries no material volume and no stress: stx is left alone and U6
!     must stay out of every *EL PRINT set.
!
      subroutine resultsmech_u6(co,kon,ipkon,lakon,ne,v,elcon,nelcon,
     &     ielmat,mi,ncmat_,ntmat_,fn,calcul_fn,nal,nelem)
!
      implicit none
!
      character*8 lakon(*)
      integer mi(*)
      integer kon(*),ipkon(*),ielmat(mi(3),*),ncmat_,ntmat_,
     &     nelcon(2,*),nelem,ne,calcul_fn,nal,i,c,nope,indexe,acen,
     &     imat,konl(255)
      real*8 co(3,*),v(0:mi(2),*),elcon(0:ncmat_,ntmat_,*),
     &     fn(0:mi(2),*),e,un,xk,va,bb(3,255),theta,fac
!
      nope=ichar(lakon(nelem)(8:8))
      indexe=ipkon(nelem)
      do i=1,nope
        konl(i)=kon(indexe+i)
      enddo
      acen=konl(1)
!
      imat=ielmat(1,nelem)
      e=elcon(1,1,imat)
      un=elcon(2,1,imat)
      xk=e/(3.d0*(1.d0-2.d0*un))
!
      call u6patch(co,kon,ipkon,lakon,ne,konl,nope,acen,va,bb,nelem)
      if(va.le.0.d0) return
!
!     theta = b . u
!
      theta=0.d0
      do i=1,nope
        do c=1,3
          theta=theta+bb(c,i)*v(c,konl(i))
        enddo
      enddo
      fac=xk*va*theta
!
      if(calcul_fn.eq.1) then
        do i=1,nope
          do c=1,3
            fn(c,konl(i))=fn(c,konl(i))+fac*bb(c,i)
          enddo
        enddo
      endif
!
      nal=nal+nope
!
      return
      end
