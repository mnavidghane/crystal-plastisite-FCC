      subroutine vumat(
     1 nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, dt, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     5 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
     7 stressNew, stateNew, enerInternNew, enerInelasNew)

      include 'vaba_param.inc'

c----- Arguments -----------------------------------------------------
      character*80 cmname
      integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
      double precision stepTime, totalTime, dt
      double precision coordMp(3,nblock), charLength(nblock)
      double precision props(nprops), density(nblock)
      double precision strainInc(ndir+nshr,nblock), relSpinInc(nshr,nblock)
      double precision tempOld(nblock), tempNew(nblock)
      double precision stretchOld(3,3,nblock), stretchNew(3,3,nblock)
      double precision defgradOld(3,3,nblock), defgradNew(3,3,nblock)
      double precision fieldOld(nfieldv,nblock), fieldNew(nfieldv,nblock)
      double precision stressOld(ndir+nshr,nblock), stressNew(ndir+nshr,nblock)
      double precision stateOld(nstatev,nblock), stateNew(nstatev,nblock)
      double precision enerInternOld(nblock), enerInternNew(nblock)
      double precision enerInelasOld(nblock), enerInelasNew(nblock)

c----- Parameters / material -----------------------------------------
      integer NSLIP
      parameter (NSLIP = 12)

c  props(1)=E, props(2)=nu, props(3)=gam0, props(4)=mexp
c  props(5)=tau0, props(6)=taus, props(7)=h0  (Pa)
      double precision E, nu, gam0, mexp, tau0, taus, h0
      double precision lambda, mu

      double precision D(6,6)             ! isotropic elastic stiffness
      double precision Ms(3,3,NSLIP)      ! Schmid tensors

      double precision Tau(NSLIP)
      double precision Gamma_old(NSLIP), Gamma_new(NSLIP)
      double precision TauCR_old(NSLIP), TauCR_new(NSLIP)
      double precision Dgamma(NSLIP)

      double precision de(6), ds(6)
      double precision Einc(3,3), Ep(3,3), Eel(3,3)
      double precision SigmaOld(3,3), SigmaTrial(3,3)

      integer ib, i, j, islip, ntens
      double precision one, two, zero, small
      double precision at, r, signTau, gdot
      double precision arg, fac, gamma_eq

      parameter (one=1.d0, two=2.d0, zero=0.d0, small=1.d-12)

c----- Read material parameters --------------------------------------
      E    = props(1)
      nu   = props(2)
      gam0 = props(3)
      mexp = props(4)
      tau0 = props(5)
      taus = props(6)
      h0   = props(7)

c  Lame constants
      mu     = E / (two*(one+nu))
      lambda = E*nu / ((one+nu)*(one-two*nu))

c----- Isotropic elastic stiffness D(6x6) ----------------------------
      do i = 1,6
        do j = 1,6
          D(i,j) = 0.d0
        end do
      end do

      D(1,1) = lambda + 2.d0*mu
      D(2,2) = lambda + 2.d0*mu
      D(3,3) = lambda + 2.d0*mu

      D(1,2) = lambda
      D(1,3) = lambda
      D(2,1) = lambda
      D(2,3) = lambda
      D(3,1) = lambda
      D(3,2) = lambda

      D(4,4) = 2.d0*mu
      D(5,5) = 2.d0*mu
      D(6,6) = 2.d0*mu

c----- FCC slip systems ----------------------------------------------
      call init_fcc_schmid(Ms)

      ntens = ndir + nshr

c=====================================================================
c  Loop over blocks
c=====================================================================
      do ib = 1, nblock

c  temperature / fields: copy
        tempNew(ib) = tempOld(ib)
        do i = 1, nfieldv
          fieldNew(i,ib) = fieldOld(i,ib)
        end do

c  strain increment (Voigt)
        do i = 1,6
          if (i .le. ntens) then
            de(i) = strainInc(i,ib)
          else
            de(i) = 0.d0
          end if
        end do

c=====================================================================
c  CASE 1: stepTime = 0  → purely elastic check for dilatational modulus
c=====================================================================
        if (stepTime .eq. 0.d0) then

c         elastic increment: ds = D * de
          do i = 1,6
            ds(i) = 0.d0
            do j = 1,6
              ds(i) = ds(i) + D(i,j)*de(j)
            end do
          end do

c         update stress
          do i = 1, ntens
            stressNew(i,ib) = stressOld(i,ib) + ds(i)
          end do

c         just copy state variables (or init CRSS if صفر است)
          do islip = 1, NSLIP
            Gamma_old(islip) = stateOld(islip,ib)
            TauCR_old(islip) = stateOld(NSLIP+islip,ib)
            if (TauCR_old(islip) .le. 0.d0) TauCR_old(islip) = tau0

            stateNew(islip,ib)       = Gamma_old(islip)
            stateNew(NSLIP+islip,ib) = TauCR_old(islip)
          end do

c         energies
          enerInternNew(ib) = enerInternOld(ib)
          enerInelasNew(ib) = enerInelasOld(ib)

          cycle   ! برو سراغ ib بعدی، پلاستیسیته فعلاً خاموش
        end if

c=====================================================================
c  CASE 2: stepTime > 0  → crystal plasticity active
c=====================================================================

c  init state at very first increment (totalTime=0)
        if (totalTime .eq. 0.d0) then
          do islip = 1, NSLIP
            stateNew(islip,ib)       = 0.d0
            stateNew(NSLIP+islip,ib) = tau0
          end do
        end if

c  read old state
        do islip = 1, NSLIP
          Gamma_old(islip) = stateOld(islip,ib)
          TauCR_old(islip) = stateOld(NSLIP+islip,ib)
          if (TauCR_old(islip) .le. 0.d0) TauCR_old(islip) = tau0
        end do

c  total strain tensor
        call mat_from_vec(de, Einc)

c  trial stress tensor (elastic increment only)
        call mat_from_vec(ds, SigmaTrial)

c  add current (old) stress state -> total trial stress
        do i = 1,3
          do j = 1,3
            SigmaTrial(i,j) = SigmaTrial(i,j) + SigmaOld(i,j)
          end do
        end do
        
c  resolved shear stresses
        do islip = 1, NSLIP
          Tau(islip) = 0.d0
          do i = 1,3
            do j = 1,3
              Tau(islip) = Tau(islip) + Ms(i,j,islip)*SigmaTrial(i,j)
            end do
          end do
        end do

c  slip increments
        do islip = 1, NSLIP
          at = dabs(Tau(islip))
          if (at .le. small) then
            Dgamma(islip) = 0.d0
          else
            if (TauCR_old(islip) .le. small) TauCR_old(islip) = tau0
            r = at / TauCR_old(islip)
            if (r .lt. small) r = small
            signTau = Tau(islip) / at
            gdot = gam0 * (r**(one/mexp)) * signTau
            Dgamma(islip) = gdot * dt
          end if
        end do

c  Ep = sum_s Ms * Dgamma_s
        call zero_tensor(Ep)
        do islip = 1, NSLIP
          do i = 1,3
            do j = 1,3
              Ep(i,j) = Ep(i,j) + Ms(i,j,islip)*Dgamma(islip)
            end do
          end do
        end do

c  elastic strain increment
        do i = 1,3
          do j = 1,3
            Eel(i,j) = Einc(i,j) - Ep(i,j)
          end do
        end do

c  back to Voigt و ds = D * de_el
        call vec_from_mat(Eel, de)
        do i = 1,6
          ds(i) = 0.d0
          do j = 1,6
            ds(i) = ds(i) + D(i,j)*de(j)
          end do
        end do

c  new stress
        do i = 1, ntens
          stressNew(i,ib) = stressOld(i,ib) + ds(i)
        end do

c  update Gamma و CRSS (Voce)
        do islip = 1, NSLIP
          Gamma_new(islip) = Gamma_old(islip) + dabs(Dgamma(islip))

          if (taus .gt. tau0) then
            gamma_eq = Gamma_new(islip)
            arg = -h0 * gamma_eq / (taus - tau0)
            fac = 1.d0 - dexp(arg)
            TauCR_new(islip) = tau0 + (taus - tau0)*fac
          else
            TauCR_new(islip) = tau0 + h0*Gamma_new(islip)
          end if

          stateNew(islip,ib)       = Gamma_new(islip)
          stateNew(NSLIP+islip,ib) = TauCR_new(islip)
        end do

c  energies
        enerInternNew(ib) = enerInternOld(ib)
        enerInelasNew(ib) = enerInelasOld(ib)

      end do

      return
      end

c=====================================================================
c   FCC slip systems {111}<110>
c=====================================================================
      subroutine init_fcc_schmid(Ms)
      include 'vaba_param.inc'
      integer NSLIP
      parameter (NSLIP = 12)
      double precision Ms(3,3,NSLIP)
      double precision n(3), s(3)
      double precision sqrt2, sqrt3
      integer i,j,islip

      sqrt2 = dsqrt(2.d0)
      sqrt3 = dsqrt(3.d0)

      do islip = 1, NSLIP
        do i = 1,3
          do j = 1,3
            Ms(i,j,islip) = 0.d0
          end do
        end do
      end do

c  Plane (1 1 1)
      n(1)= 1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)= 1.d0/sqrt3

c  Sys1: [0 -1 1]
      s(1)= 0.d0
      s(2)=-1.d0/sqrt2
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,1)

c  Sys2: [1 0 -1]
      s(1)= 1.d0/sqrt2
      s(2)= 0.d0
      s(3)=-1.d0/sqrt2
      call build_schmid(n,s,Ms,2)

c  Sys3: [-1 1 0]
      s(1)=-1.d0/sqrt2
      s(2)= 1.d0/sqrt2
      s(3)= 0.d0
      call build_schmid(n,s,Ms,3)

c  Plane (1 -1 1)
      n(1)= 1.d0/sqrt3
      n(2)=-1.d0/sqrt3
      n(3)= 1.d0/sqrt3

c  Sys4: [0 1 1]
      s(1)= 0.d0
      s(2)= 1.d0/sqrt2
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,4)

c  Sys5: [1 0 1]
      s(1)= 1.d0/sqrt2
      s(2)= 0.d0
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,5)

c  Sys6: [-1 -1 0]
      s(1)=-1.d0/sqrt2
      s(2)=-1.d0/sqrt2
      s(3)= 0.d0
      call build_schmid(n,s,Ms,6)

c  Plane (-1 1 1)
      n(1)=-1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)= 1.d0/sqrt3

c  Sys7: [0 1 -1]
      s(1)= 0.d0
      s(2)= 1.d0/sqrt2
      s(3)=-1.d0/sqrt2
      call build_schmid(n,s,Ms,7)

c  Sys8: [1 0 1]
      s(1)= 1.d0/sqrt2
      s(2)= 0.d0
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,8)

c  Sys9: [-1 1 0]
      s(1)=-1.d0/sqrt2
      s(2)= 1.d0/sqrt2
      s(3)= 0.d0
      call build_schmid(n,s,Ms,9)

c  Plane (1 1 -1)
      n(1)= 1.d0/sqrt3
      n(2)= 1.d0/sqrt3
      n(3)=-1.d0/sqrt3

c  Sys10: [0 1 1]
      s(1)= 0.d0
      s(2)= 1.d0/sqrt2
      s(3)= 1.d0/sqrt2
      call build_schmid(n,s,Ms,10)

c  Sys11: [1 0 -1]
      s(1)= 1.d0/sqrt2
      s(2)= 0.d0
      s(3)=-1.d0/sqrt2
      call build_schmid(n,s,Ms,11)

c  Sys12: [-1 -1 0]
      s(1)=-1.d0/sqrt2
      s(2)=-1.d0/sqrt2
      s(3)= 0.d0
      call build_schmid(n,s,Ms,12)

      return
      end

c=====================================================================
      subroutine build_schmid(n,s,Ms,idx)
      include 'vaba_param.inc'
      double precision n(3), s(3)
      double precision Ms(3,3,*)
      integer idx, i, j

      do i = 1,3
        do j = 1,3
          Ms(i,j,idx) = 0.5d0*(s(i)*n(j) + n(i)*s(j))
        end do
      end do

      return
      end

c=====================================================================
      subroutine mat_from_vec(v, A)
      include 'vaba_param.inc'
      double precision v(6), A(3,3)

      A(1,1)=v(1)
      A(2,2)=v(2)
      A(3,3)=v(3)
      A(1,2)=v(4)
      A(2,1)=v(4)
      A(1,3)=v(5)
      A(3,1)=v(5)
      A(2,3)=v(6)
      A(3,2)=v(6)
      return
      end

c=====================================================================
      subroutine vec_from_mat(A, v)
      include 'vaba_param.inc'
      double precision A(3,3), v(6)

      v(1)=A(1,1)
      v(2)=A(2,2)
      v(3)=A(3,3)
      v(4)=A(1,2)
      v(5)=A(1,3)
      v(6)=A(2,3)
      return
      end

c=====================================================================
      subroutine zero_tensor(T)
      include 'vaba_param.inc'
      double precision T(3,3)
      integer i,j

      do i=1,3
        do j=1,3
          T(i,j)=0.d0
        end do
      end do
      return
      end
