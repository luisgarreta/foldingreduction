c        1         2         3         4         5         6         7 !
c 3456789012345678901234567890123456789012345678901234567890123456789012345678

      subroutine getTMscore (proteinModel, proteinNative, resTMscore)
      character*500 proteinModel, proteinNative
      double precision  resTMscore

      PARAMETER(nmax=5000)
      
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)

      character*500 fnam,pdb(2)
      character*3 aa(-1:20),seqA(nmax),seqB(nmax)
      character*500 s,du
      character*1 chA(nmax),chB(nmax)
      character*1 chain1(90000),chain2(90000)
      character seq1A(nmax),seq1B(nmax),ali(nmax),du1
      character sequenceA(nmax),sequenceB(nmax),sequenceM(nmax)
      character ins1(nmax),ins2(nmax),ains1(90000),ains2(90000)
      
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score
      dimension xa(nmax),ya(nmax),za(nmax)
      
      character*10 aa1,ra1,aa2,ra2
      dimension ia1(90000),aa1(90000),ra1(90000),ir1(90000)
      dimension xa1(90000),ya1(90000),za1(90000)
      dimension ia2(90000),aa2(90000),ra2(90000),ir2(90000)
      dimension xa2(90000),ya2(90000),za2(90000)

      dimension nres1(nmax,32:122),nres2(nmax,32:122) !number of atoms
      character*500 atom1(nmax,30),atom2(nmax,30) !atom name
c     here atom1(i,j) should be atom2(i,j,k) but will beyond dimension

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

      data aa/ 'BCK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP','CYX'/
      character*1 slc(-1:20)
      data slc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W','C'/

******* options ----------->
      m_out=-1
      narg=iargc()
      m_complex=0

      pdb(1)=proteinModel
      pdb(2)=proteinNative
      
ccccccccc read data from first CA file:
      
      open(unit=10,file=pdb(1),status='old')
      i=0
      na1=0
      ntt=0
 101  read(10,104,end=102) s
      if(m_complex.eq.0)then
         if(s(1:3).eq.'TER') goto 102
      endif
      if(s(1:4).eq.'ATOM')then
***   
         if(s(27:27).eq.' ')then
            du1=' '
         else
            read(s(27:27),*)du1 !Residue_Insertion_tag
         endif
         mk=1
         if(s(17:17).ne.' ')then !ignor alternative residues/atoms
            read(s(13:16),*)du
            read(s(23:26),*)i8
            do i1=1,nres1(i8,ichar(du1))
               if(du.eq.atom1(i8,i1))then !atom is a ALTloc
                  mk=-1
               endif
            enddo
         endif
         if(mk.eq.1)then
            read(s(13:16),*)du
            read(s(23:26),*)i8
            if(nres1(i8,ichar(du1)).lt.30)then
               nres1(i8,ichar(du1))=nres1(i8,ichar(du1))+1
               atom1(i8,nres1(i8,ichar(du1)))=du
            endif
***   
            ntt=ntt+1
            if(ntt.ge.90000)goto 102
            na1=na1+1
            read(s,8999)du,ia1(na1),du,aa1(na1),du,ra1(na1),du,
     &           chain1(na1),ir1(na1),du,xa1(na1),ya1(na1),za1(na1)
            ains1(na1)=du1
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &           eq.'  CA')then
               i=i+1
               read(s,103)du,seqA(i),du,chA(i),nresA(i),du,
     &              xa(i),ya(i),za(i)
               ins1(i)=du1
               do j=-1,20
                  if(seqA(i).eq.aa(j))then
                     seq1A(i)=slc(j)
                     goto 21
                  endif
               enddo
               seq1A(i)=slc(-1)
 21            continue
               if(i.ge.nmax)goto 102
            endif
         endif
      endif
      goto 101
 102  continue
 103  format(A17,A3,A1,A1,i4,A4,3F8.3)
 104  format(A100)
 8999 format(a6,I5,a1,A4,a1,A3,a1,a1,I4,a4,3F8.3)
      close(10)
      nseqA=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ccccccccc read data from first CA file:
      open(unit=10,file=pdb(2),status='old')
      i=0
      na2=0
      ntt=0
 201  read(10,204,end=202) s
      if(m_complex.eq.0)then
         if(s(1:3).eq.'TER') goto 202
      endif
      if(s(1:4).eq.'ATOM')then
***   
         if(s(27:27).eq.' ')then
            du1=' '
         else
            read(s(27:27),*)du1
         endif
         mk=1
         if(s(17:17).ne.' ')then
            read(s(13:16),*)du
            read(s(23:26),*)i8
            do i1=1,nres2(i8,ichar(du1))
               if(du.eq.atom2(i8,i1))then
                  mk=-1
               endif
            enddo
         endif
         if(mk.eq.1)then
            read(s(13:16),*)du
            read(s(23:26),*)i8
            if(nres2(i8,ichar(du1)).lt.30)then
               nres2(i8,ichar(du1))=nres2(i8,ichar(du1))+1
               atom2(i8,nres2(i8,ichar(du1)))=du
            endif
***   
            ntt=ntt+1
            if(ntt.ge.90000)goto 202
            na2=na2+1
            read(s,8999)du,ia2(na2),du,aa2(na2),du,ra2(na2),du,
     &           chain2(na2),ir2(na2),du,xa2(na2),ya2(na2),za2(na2)
            ains2(na2)=du1
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &           eq.'  CA')then
               i=i+1
               read(s,203)du,seqB(i),du,chB(i),nresB(i),du,
     &              xb(i),yb(i),zb(i)
               ins2(i)=du1
               do j=-1,20
                  if(seqB(i).eq.aa(j))then
                     seq1B(i)=slc(j)
                     goto 22
                  endif
               enddo
               seq1B(i)=slc(-1)
 22            continue
               if(i.ge.nmax)goto 202
            endif
         endif
      endif
      goto 201
 202  continue
 203  format(A17,A3,A1,A1,i4,A4,3F8.3)
 204  format(A100)
      close(10)
      nseqB=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
******************************************************************
*     pickup the aligned residues:
******************************************************************
      if(m_complex.eq.0)then    !monomer
         k=0
         do i=1,nseqA
            do j=1,nseqB
               if(nresA(i).eq.nresB(j))then
                  if(ins1(i).eq.ins2(j))then
                     k=k+1
                     iA(k)=i
                     iB(k)=j
                     goto 205
                  endif
               endif
            enddo
 205        continue
         enddo
      else                      !complex
         k=0
         do i=1,nseqA
            do j=1,nseqB
               if(nresA(i).eq.nresB(j).and.chA(i).eq.chB(j))then
                  if(ins1(i).eq.ins2(j))then !Residue_insertion_code
                     k=k+1
                     iA(k)=i
                     iB(k)=j
                     goto 206
                  endif
               endif
            enddo
 206        continue
         enddo
      endif
      n_ali=k                   !number of aligned residues
      if(n_ali.lt.1)then
         write(*,*)'There is no common residues in the input structures'
         goto 9999
      endif
      
******************************************************************
*     check the residue serial numbers ------------->
      if(m_complex.eq.0)then
         nA_repeat=0
         do i=1,nseqA
            do j=i+1,nseqA
               if(nresA(i).eq.nresA(j))then
                  if(ins1(i).eq.ins1(j))then
                     nA_repeat=nA_repeat+1
                  endif
               endif
            enddo
         enddo
         if(nA_repeat.gt.0)then
            write(*,380)nA_repeat
         endif
 380     format('Warning: TMscore calculation can be wrong because ',
     &        'there are ',I3,' residues with the serial number same ',
     &        'as others in first Chain. Please modify the PDB file ',
     &        'and rerun the program!!')
         nB_repeat=0
         do i=1,nseqB
            do j=i+1,nseqB
               if(nresB(i).eq.nresB(j))then
                  if(ins2(i).eq.ins2(j))then
                     nB_repeat=nB_repeat+1
                  endif
               endif
            enddo
         enddo
         if(nB_repeat.gt.0)then
            write(*,381)nB_repeat
         endif
 381     format('Warning: TMscore calculation can be wrong because ',
     &        'there are ',I3,' residues with the serial number same ',
     &        'as others in second Chain. Please modify the PDB file ',
     &        'and rerun the program!!')
      endif
      
************/////
*     parameters:
*****************
***   d0------------->
      d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      if(d0.lt.0.5)d0=0.5
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue
      
******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              ka=ka+1
              k_ali(ka)=k
              LL=LL+1
           enddo
           if(i_init.eq.1)then  !global superposition
              call u3b(w,r_1,r_2,LL,2,rms,u,t,ier) !0:rmsd; 1:u,t; 2:rmsd,u,t
              armsd=dsqrt(rms/LL)
              rmsd_ali=armsd
           else
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M
      
******************************************************************
*     Output
******************************************************************
***   output TM-scale ---------------------------->
      resTMscore=score_max
      return
 9999 END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun
      PARAMETER(nmax=5000)

      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score

      d_tmp=d
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xt(i)-xb(j))**2+(yt(i)-yb(j))**2+(zt(i)-zb(j))**2)
***   for iteration:
         if(dis.lt.d_tmp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
***   for TM-score:
         score_sum=score_sum+1/(1+(dis/d0)**2)
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         goto 21
      endif
      score=score_sum/float(nseqB) !TM-score
      
      return
      end

cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m                    (given)
c  x    - x(i,m) are coordinates of atom c m in set x          (given)
c  y    - y(i,m) are coordinates of atom c m in set y          (given)
c  n    - n is number of atom pairs                            (given)
c  mode  - 0:calculate rms     only                            (given,short)
c          1:calculate     u,t only                            (given,medium)
c          2:calculate rms,u,t                                 (given,longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs            (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      double precision w(*), x(3,*), y(3,*)
      integer n, mode
      
      double precision rms, u(3,3), t(3)
      integer ier
      
      integer i, j, k, l, m1, m
      integer ip(9), ip2312(4)
      double precision r(3,3), xc(3), yc(3), wc
      double precision a(3,3), b(3,3), e(3), rr(6), ss(6)
      double precision e0, d, spur, det, cof, h, g
      double precision cth, sth, sqrth, p, sigma
      double precision c1x, c1y, c1z, c2x, c2y, c2z
      double precision s1x, s1y, s1z, s2x, s2y, s2z
      double precision sxx, sxy, sxz, syx, syy, syz, szx, szy, szz
      
      double precision sqrt3, tol, zero
      
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
      
      wc  = zero
      rms = zero
      e0  = zero
      s1x = zero
      s1y = zero
      s1z = zero
      s2x = zero
      s2y = zero
      s2z = zero
      sxx = zero
      sxy = zero
      sxz = zero
      syx = zero
      syy = zero
      syz = zero
      szx = zero 
      szy = zero
      szz = zero
      
      do i=1, 3
         xc(i) = zero
         yc(i) = zero
         t(i) = zero
         do j=1, 3
            r(i,j) = zero
            u(i,j) = zero
            a(i,j) = zero
            if( i .eq. j ) then
               u(i,j) = 1.0
               a(i,j) = 1.0
            end if
         end do
      end do
      
      ier = -1
      if( n .lt. 1 ) return
      ier = -2
      
      do m=1, n
         c1x=x(1, m)
         c1y=x(2, m)
         c1z=x(3, m)
         
         c2x=y(1, m)
         c2y=y(2, m)
         c2z=y(3, m)
         
         s1x = s1x + c1x
         s1y = s1y + c1y;
         s1z = s1z + c1z;
         
         s2x = s2x + c2x;
         s2y = s2y + c2y;
         s2z = s2z + c2z;
         
         sxx = sxx + c1x*c2x; 
         sxy = sxy + c1x*c2y; 
         sxz = sxz + c1x*c2z; 
         
         syx = syx + c1y*c2x; 
         syy = syy + c1y*c2y; 
         syz = syz + c1y*c2z;
         
         szx = szx + c1z*c2x; 
         szy = szy + c1z*c2y; 
         szz = szz + c1z*c2z;
      end do
      
      xc(1) = s1x/n;
      xc(2) = s1y/n; 	
      xc(3) = s1z/n;
      
      yc(1) = s2x/n;
      yc(2) = s2y/n; 	
      yc(3) = s2z/n;
      if(mode.eq.2.or.mode.eq.0) then ! need rmsd                     		
         do m=1, n		
            do i=1, 3
               e0 = e0+ (x(i, m)-xc(i))**2 + (y(i, m)-yc(i))**2			
            end do				
         end do
      endif
      
      r(1, 1) = sxx-s1x*s2x/n;
      r(2, 1) = sxy-s1x*s2y/n;
      r(3, 1) = sxz-s1x*s2z/n;
      r(1, 2) = syx-s1y*s2x/n;
      r(2, 2) = syy-s1y*s2y/n;
      r(3, 2) = syz-s1y*s2z/n;
      r(1, 3) = szx-s1z*s2x/n;
      r(2, 3) = szy-s1z*s2y/n;
      r(3, 3) = szz-s1z*s2z/n;
      
      det = r(1,1) * ( (r(2,2)*r(3,3)) - (r(2,3)*r(3,2)) )
     &     - r(1,2) * ( (r(2,1)*r(3,3)) - (r(2,3)*r(3,1)) )
     &     + r(1,3) * ( (r(2,1)*r(3,2)) - (r(2,2)*r(3,1)) )
      
      sigma = det
      
      m = 0
      do j=1, 3
         do i=1, j
            m = m+1
            rr(m) = r(1,i)*r(1,j) + r(2,i)*r(2,j) + r(3,i)*r(3,j)
         end do
      end do
      
      spur = (rr(1)+rr(3)+rr(6)) / 3.0
      cof = (((((rr(3)*rr(6) - rr(5)*rr(5)) + rr(1)*rr(6))
     &     - rr(4)*rr(4)) + rr(1)*rr(3)) - rr(2)*rr(2)) / 3.0
      det = det*det
      
      do i=1, 3
         e(i) = spur
      end do
      if( spur .le. zero ) goto 40
      d = spur*spur
      h = d - cof
      g = (spur*cof - det)/2.0 - spur*h
      if( h .le. zero ) then
         if( mode .eq. 0 ) then
            goto 50
         else
            goto 30
         end if
      end if
      sqrth = dsqrt(h)
      d = h*h*h - g*g
      if( d .lt. zero ) d = zero
      d = datan2( dsqrt(d), -g ) / 3.0
      cth = sqrth * dcos(d)
      sth = sqrth*sqrt3*dsin(d)
      e(1) = (spur + cth) + cth
      e(2) = (spur - cth) + sth
      e(3) = (spur - cth) - sth
      
      if( mode .eq. 0 ) then
         goto 50
      end if
      
      do l=1, 3, 2
         d = e(l)
         ss(1) = (d-rr(3)) * (d-rr(6))  - rr(5)*rr(5)
         ss(2) = (d-rr(6)) * rr(2)      + rr(4)*rr(5)
         ss(3) = (d-rr(1)) * (d-rr(6))  - rr(4)*rr(4)
         ss(4) = (d-rr(3)) * rr(4)      + rr(2)*rr(5)
         ss(5) = (d-rr(1)) * rr(5)      + rr(2)*rr(4)
         ss(6) = (d-rr(1)) * (d-rr(3))  - rr(2)*rr(2)
         
         if( dabs(ss(1)) .ge. dabs(ss(3)) ) then
            j=1
            if( dabs(ss(1)) .lt. dabs(ss(6)) ) j = 3
         else if( dabs(ss(3)) .ge. dabs(ss(6)) ) then
            j = 2
         else
            j = 3
         end if
         
         d = zero
         j = 3 * (j - 1)
         
         do i=1, 3
            k = ip(i+j)
            a(i,l) = ss(k)
            d = d + ss(k)*ss(k)
         end do
         if( d .gt. zero ) d = 1.0 / dsqrt(d)
         do i=1, 3
            a(i,l) = a(i,l) * d
         end do
      end do
      
      d = a(1,1)*a(1,3) + a(2,1)*a(2,3) + a(3,1)*a(3,3)
      if ((e(1) - e(2)) .gt. (e(2) - e(3))) then
         m1 = 3
         m = 1
      else
         m1 = 1
         m = 3
      endif
      
      p = zero
      do i=1, 3
         a(i,m1) = a(i,m1) - d*a(i,m)
         p = p + a(i,m1)**2
      end do
      if( p .le. tol ) then
         p = 1.0
         do 21 i=1, 3
            if (p .lt. dabs(a(i,m))) goto 21
            p = dabs( a(i,m) )
            j = i
 21      continue
         k = ip2312(j)
         l = ip2312(j+1)
         p = dsqrt( a(k,m)**2 + a(l,m)**2 )
         if( p .le. tol ) goto 40
         a(j,m1) = zero
         a(k,m1) = -a(l,m)/p
         a(l,m1) =  a(k,m)/p
      else
         p = 1.0 / dsqrt(p)
         do i=1, 3
            a(i,m1) = a(i,m1)*p
         end do
      end if
      
      a(1,2) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
      a(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
      a(3,2) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
      
 30   do l=1, 2
         d = zero
         do i=1, 3
            b(i,l) = r(i,1)*a(1,l) + r(i,2)*a(2,l) + r(i,3)*a(3,l)
            d = d + b(i,l)**2
         end do
         if( d .gt. zero ) d = 1.0 / dsqrt(d)
         do i=1, 3
            b(i,l) = b(i,l)*d
         end do
      end do
      d = b(1,1)*b(1,2) + b(2,1)*b(2,2) + b(3,1)*b(3,2)
      p = zero
      
      do i=1, 3
         b(i,2) = b(i,2) - d*b(i,1)
         p = p + b(i,2)**2
      end do
      if( p .le. tol ) then
         p = 1.0
         do 22 i=1, 3
            if(p.lt.dabs(b(i,1)))goto 22
            p = dabs( b(i,1) )
            j = i
 22      continue
         k = ip2312(j)
         l = ip2312(j+1)
         p = dsqrt( b(k,1)**2 + b(l,1)**2 )
         if( p .le. tol ) goto 40
         b(j,2) = zero
         b(k,2) = -b(l,1)/p
         b(l,2) =  b(k,1)/p
      else
         p = 1.0 / dsqrt(p)
         do i=1, 3
            b(i,2) = b(i,2)*p
         end do
      end if
      
      b(1,3) = b(2,1)*b(3,2) - b(2,2)*b(3,1)
      b(2,3) = b(3,1)*b(1,2) - b(3,2)*b(1,1)
      b(3,3) = b(1,1)*b(2,2) - b(1,2)*b(2,1)
      
      do i=1, 3
         do j=1, 3
            u(i,j) = b(i,1)*a(j,1) + b(i,2)*a(j,2) + b(i,3)*a(j,3)
         end do
      end do
      
 40   do i=1, 3
         t(i) = ((yc(i) - u(i,1)*xc(1)) - u(i,2)*xc(2)) - u(i,3)*xc(3)
      end do
 50   do i=1, 3
         if( e(i) .lt. zero ) e(i) = zero
         e(i) = dsqrt( e(i) )
      end do
      
      ier = 0
      if( e(2) .le. (e(1) * 1.0d-05) ) ier = -1
      
      d = e(3)
      if( sigma .lt. 0.0 ) then
         d = - d
         if( (e(2) - e(3)) .le. (e(1) * 1.0d-05) ) ier = -1
      end if
      d = (d + e(2)) + e(1)
      
      if(mode .eq. 2.or.mode.eq.0) then ! need rmsd                     		
         rms = (e0 - d) - d
         if( rms .lt. 0.0 ) rms = 0.0
      endif
      
      return
      end
      

