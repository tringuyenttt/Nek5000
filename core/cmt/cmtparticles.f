c----------------------------------------------------------------------
c routine called in case of particle calls only in .usr file (i.e.,
c  with nek5000 particles and not cmt-nek particles. must use 
c  bdf/ext time integration. Otherwise, cmt-nek will not call this fxn.
      subroutine stokes_particles
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      nstage_part = 3
      if (abs(time_integ) .eq. 2) nstage_part = 1

      if (istep.eq.0) then
         call usr_particles_init(0)
      else

         call set_tstep_coef_part(dt) ! in nek5000 with rk3
         do stage=1,nstage_part
            call usr_particles_solver
         enddo
         call compute_phig_qtl(dt,usrdiv) ! nek5000 (see Zwick 2018)
      endif

      if(mod(istep,iostep).eq.0.or. istep.eq.1) then
         call usr_particles_io
      endif

      return
      end
c----------------------------------------------------------------------
c     setup routines
c----------------------------------------------------------------------
      subroutine usr_particles_init(idum)
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'CMTDATA'
      include 'CMTPART'

      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload

      icmtp = idum

      call set_part_pointers
      if (llpart .ne. 1) call read_particle_input
      call set_bounds_box
      call set_part_params ! n initialized here
      call place_particles
      call update_particle_location   ! move outlier particles
      call set_check_spl_params ! in case spl/collisions are set
      call move_particles_inproc          ! initialize fp & cr comm handles
      if (red_part .le. 2) call init_interpolation ! barycentric weights for interpolation
      if (two_way.gt.1) then
         call compute_neighbor_el_proc    ! compute list of neigh. el. ranks 
         call create_extra_particles
         call send_ghost_particles
c        call point_to_grid_corr_init    ! for gamma correction integrat
         call spread_props_grid           ! put particle props on grid

c        do i = 1,nitspl
c           call interp_props_part_location ! interpolate
c           call correct_spl
c           call create_extra_particles
c           call send_ghost_particles
c           call spread_props_grid           ! put particle props on grid
c           if (nid.eq.0) write(6,*) i,'Pre-SPL iteration'
c        enddo
         
         call set_check_spl_params        !  in case spl has changed!
         call create_extra_particles
         call send_ghost_particles
         call spread_props_grid           ! put particle props on grid
      endif
      call interp_props_part_location ! interpolate again for two-way

      if (time_integ .lt. 0) call pre_sim_collisions ! e.g., settling p

      resetFindpts = 0
c     call computeRatio
c     call reinitialize
c     call printVerify

      ntmp  = iglsum(n,1)
      if (nid.eq.0) write(6,*) 'Passed usr_particles_init'

      return
      end
c----------------------------------------------------------------------
      subroutine place_particles
c
c     Place particles in this routine, also called for injection
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      integer icalld
      save    icalld
      data    icalld  /-1/

      icalld = icalld + 1

      if (ipart_restartr .eq. 0) then

c        correct nwe if discrepancy on rank 0
         nwe         = int(nw/np)                ! num. part per proc
         nw_tmp      = iglsum(nwe,1)
         ndef        = nw - nw_tmp
         if (nid .lt. ndef) nwe = nwe + 1
c        if ((nw_tmp .ne. nw) .and. (nid.eq.0)) nwe = nwe +(nw - nw_tmp)

c        main loop to distribute particles
         do i_pt_part = 1,nwe
            n = n + 1

            call place_particles_user

            do j=0,2
               rpart(jx +j,n) = x_part(j)
               rpart(jx1+j,n) = x_part(j)
               rpart(jx2+j,n) = x_part(j)
               rpart(jx3+j,n) = x_part(j)

               rpart(jv0+j,n) = v_part(j)
               rpart(jv1+j,n) = v_part(j)
               rpart(jv2+j,n) = v_part(j)
               rpart(jv3+j,n) = v_part(j)
            enddo
         
c           set some rpart values for later use
            rpart(jdp,n)   = d_part                              ! particle diameter
            rpart(jtaup,n) = rpart(jdp,n)**2*rho_p/18.0d+0/mu_0  ! particle time scale
            rpart(jrhop,n) = rho_p                               ! particle density 
            rpart(jvol,n)  = pi*rpart(jdp,n)**3/6.               ! particle volume
            rpart(jspl,n)  = rspl                                ! super particle loading

            rpart(jrpe,n) = rpart(jspl,n)**(1./3.)*rpart(jdp,n)/2.
            rpart(jvol,n) = rpart(jspl,n)*rpart(jvol,n)
         
            rpart(jtemp,n)  = tp_0                               ! intial particle temp
            rpart(jtempf,n) = tp_0                               ! intial fluid temp (overwritten)
            rpart(jrho,n)   = param(1)                           ! initial fluid density (overwritten interp)
         
c           set global particle id (3 part tag)
            ipart(jpid1,n) = nid 
            ipart(jpid2,n) = i_pt_part
            ipart(jpid3,n) = icalld
         enddo

      else
         ! read in data
         nread_part = 1
         do j=1,nread_part
            call read_parallel_restart_part
c           do i=1,n
c           ry = rpart(jy,i) + rpart(jrpe,i)
c           if (ry .gt. -0.165) then
c               rpart(jy,i) = 1E8
c           endif
c           enddo
            call update_particle_location   ! move outlier particles
            call move_particles_inproc
         enddo

      endif

      ! Error checking
      if (n.gt.llpart)then 
         if (nid.eq.0)
     >      write(6,*)'Not enough space to store more particles'
         call exitt
      endif

      if (.not. if3d) then
         do i=1,n
            if (abs(rpart(jz,i)-1.0) .gt. 1E-16) then
               if (nid.eq.0)
     >            write(6,*)'Particle zstart is not right for 2d case'
               call exitt
            endif
         enddo
      endif

      if (nid.eq.0)
     >            write(6,*)'Finished placing/reading particles'

      return
      end
c----------------------------------------------------------------------
      subroutine set_bounds_box
c
c     set domain and element bounds for a box geometry. Notice that
c     this ONLY works with non curved elements.
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'CMTTIMERS'
      include 'CMTPART'

c     added by keke
      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload
c     end added by keke

      rdum  = 1E8
      rleng = 1E8

c     if(istep.eq.0.or.istep.eq.1)then
      if((istep.eq.0) .or. (istep.eq.1).or.(resetFindpts.eq.1)) then !resetFindpts .eq. 0 added by keke
        call domain_size(xdrange(1,1),xdrange(2,1),xdrange(1,2)
     $                  ,xdrange(2,2),xdrange(1,3),xdrange(2,3))
        ntot = lx1*ly1*lz1*nelt
        nxyz = lx1*ly1*lz1
        do ie = 1,nelt
           xerange(1,1,ie) = vlmin(xm1(1,1,1,ie),nxyz)
           xerange(2,1,ie) = vlmax(xm1(1,1,1,ie),nxyz)
           xerange(1,2,ie) = vlmin(ym1(1,1,1,ie),nxyz)
           xerange(2,2,ie) = vlmax(ym1(1,1,1,ie),nxyz)
           xerange(1,3,ie) = vlmin(zm1(1,1,1,ie),nxyz)
           xerange(2,3,ie) = vlmax(zm1(1,1,1,ie),nxyz)

           rdum1 = min(xerange(2,1,ie) - xerange(1,1,ie),
     >                 xerange(2,2,ie) - xerange(1,2,ie),
     >                 xerange(2,3,ie) - xerange(1,3,ie))
           if (rdum1 .lt. rdum) rdum = rdum1
        enddo  

        rdum1 = glmin(rdum,1)
        if (rdum1 .lt. rleng) then
           rleng = rdum1
        endif
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine set_part_params
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*132 deathmessage

      if (icalld .lt. 0) then
         rdum   = ran2(-nrandseed*np-nid-1) ! initialize random number generator
         icalld = icalld + 1
      endif

c     setup items
      nlxyze = nx1*ny1*nz1*nelt
      call rzero(rpart,lr*llpart)
      call izero(ipart,li*llpart)
      call rzero(ptw,nlxyze*8)
      call rzero(ptdum,iptlen)
      call rzero(pttime,iptlen)

c     filter width setup (note deltax is implicit in expressions b4 def)
      rtmp_rle = 0.0 ! dummy number, max value it can be

      ! do nothing, no spreading
      if (npro_method .eq. 0) then

      ! box filter in this element, still no spreading
      elseif (npro_method .eq. 1) then

      ! gaussian set by user input parameters
      elseif (npro_method .eq. 2) then

         rtmp_rle = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))

         if (rtmp_rle .gt. 0.5) then

            if (rtmp_rle .lt. 1.0) then
            if (nrect_assume .eq. 2) then
               goto 123
            else
               deathmessage = 'Resetting to full projection for filter'
               if (nid.eq. 0) write(6,*) deathmessage
               nrect_assume = 2
               goto 123
            endif
            endif

            deathmessage = 'WARNING filter width is too large'
            if (nid.eq. 0)write(6,*) deathmessage,rtmp_rle,icalld
            call exittr(deathmessage,rtmp_rle,icalld)
  123 continue
         endif

      endif

      d2chk(1) = rtmp_rle*rleng
      d2chk(2) = d2chk(1)
      d2chk(3) = d2chk(1)

      rsig     = dfilt*rleng/(2.*sqrt(2.*log(2.))) ! gaussian filter std.

      return
      end
c----------------------------------------------------------------------
      subroutine set_check_spl_params
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      character*132 deathmessage

      ! now, check this filter width against collision width
      if (two_way .gt. 2) then

         rdeff_max = dp(2)
         do i = 1,n
            if (rpart(jrpe,i) .gt. rdeff_max) rdeff_max=rpart(jrpe,i)
         enddo
         rdeff_max = glmax(rdeff_max,1)
         rdeff_max = rdeff_max*2. ! to get diameter

         rtmp_rle2 = d2chk(1)/rleng
         rtmp_rle_col = rdeff_max*1.00/rleng

         if ( abs(npro_method) .gt. 1) then
         if ( rtmp_rle_col .gt. rtmp_rle2) then
            if (rtmp_rle_col .gt. 0.5) then
               if (nrect_assume .eq. 2) then
                  if (rtmp_rle_col .gt. 1) then
                     goto 1234
                  else
                     deathmessage =  
     >                  'Collision > filter width-Resetting width'
                     if (nid.eq. 0)write(6,*) deathmessage
                  endif
               elseif (nrect_assume .eq. 1) then
                  deathmessage =  
     >              'Collision > filter width-Using full project now'
                  if (nid.eq. 0)write(6,*) deathmessage
                  nrect_assume = 2
               endif
             endif
         endif
            rtmp_rle2 = max(rtmp_rle_col,rtmp_rle2)
         else
            ! no filter dependent spreading, so use collision width
            rtmp_rle2 = rtmp_rle_col ! note r/Le > 0.5 is already caught
         endif
      endif
      goto 1237

 1234 continue
            deathmessage =  
     >        'Collision/filter width is too large!'
            if (nid.eq. 0)write(6,*) deathmessage,rdeff_max,icalld
            call exittr(deathmessage,rdeff_max,icalld)
 1237 continue

      d2chk(1) = max(d2chk(1),rtmp_rle2*rleng)
      d2chk(2) = d2chk(2)
      d2chk(3) = rtmp_rle_col*rleng

      return
      end
c----------------------------------------------------------------------
      subroutine set_dt_particles(rdt_part)
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real dt_dum,dt_col,cflp,cflt

      common /save_dt_part/ rdpe_max, rdpe_min

      integer icalld
      save    icalld
      data    icalld  /-1/
      
      if (llpart .eq. 1) goto 1234

      icalld = icalld + 1
      if (icalld .eq. 0) then
         rdpe_max = 0.
         rdpe_min = 100.
         do i=1,n
            if (rpart(jrpe,i).lt.rdpe_min) rdpe_min = rpart(jrpe,i)
            if (rpart(jrpe,i).gt.rdpe_max) rdpe_max = rpart(jrpe,i)
         enddo
         rdpe_min = glmin(rdpe_min,1)
         rdpe_max = glmax(rdpe_max,1)
         rdpe_min = 2.*rdpe_min ! to get diameter
         rdpe_max = 2.*rdpe_max ! to get diameter
      endif

c     ! particle cfl, particles cant move due to velocity
      dt_dum = rdt_part
      dt_part = 1000.
      cflp = 0.10
      rvmag_max = 0.
      rhmax = -1E8 
      do i=1,n
         rvmag  = sqrt(rpart(jv0,i)**2 + rpart(jv0+1,i)**2 
     >                 + rpart(jv0+2,i)**2)
         if (rpart(jy,i) + rpart(jrpe,i) .gt. rhmax) 
     >                   rhmax = rpart(jy,i) + rpart(jrpe,i)

         cflt = dt_dum*rvmag/rdpe_min
         if (cflt .lt. cflp) dt_part = dt_dum
         if (cflt .ge. cflp) dt_part = cflp*rdpe_min/rvmag ! make sure smallest small overlap

         if (rvmag .gt. rvmag_max) rvmag_max = rvmag
      enddo
      rvmag_max = glmax(rvmag_max,1)
      dt_part  = glmin(dt_part,1)
      rhmax    = glmax(rhmax,1)

      ! resolving collisions
      rm1      = rho_p*pi/6.*rdpe_min**3 ! max
      rm2      = rho_p*pi/6.*rdpe_max**3 ! max
      rm12     = 1./(1./rm1 + 1./rm2)
      n_resolve= 10
      dt_col   = sqrt(rm12/ksp*(log(e_rest)**2+pi**2))/n_resolve

      if (two_way .gt. 2) then
         rdt_part = min(dt_part,dt_col)
      else
         rdt_part = dt_part ! don't set if no collisions!
      endif

c     if (nid.eq.0) write(6,*) 'PART DT:', 
c    >      rdt_part,dt_part,dt_col,rvmag_max,rhmax

 1234 continue

      return
      end
c----------------------------------------------------------------------
      subroutine set_part_pointers
      include 'SIZE'
      include 'CMTPART'

      nr   = lr     ! Mandatory for proper striding
      ni   = li     ! Mandatory
      nrgp = lrgp
      nigp = ligp
      nrf  = lrf
      nif  = lif
      nw   = 0
      n    = 0

c     ipart pointers ------------------------------------------------
      jrc   = 1 ! Pointer to findpts return code
      jpt   = 2 ! Pointer to findpts return processor id
      je0   = 3 ! Pointer to findpts return element id
      jps   = 4 ! Pointer to proc id for data swap
      jpid1 = 5 ! initial proc number
      jpid2 = 6 ! initial local particle id
      jpid3 = 7 ! initial time step introduced
      jicx  = 8 ! initial time step introduced
      jicy  = 9 ! initial time step introduced
      jicz  = 10 ! initial time step introduced
      jai   = 11 ! Pointer to auxiliary integers

      nai = ni - (jai-1)  ! Number of auxiliary integers
      if (nai.le.0) call exitti('Error in nai:$',ni)

c     rpart pointers ------------------------------------------------
      jr  = 1         ! Pointer to findpts return rst variables
      jd  = jr + 3    ! Pointer to findpts return distance
      jx  = jd + 1    ! Pointer to findpts input x value
      jy  = jx + 1    ! Pointer to findpts input y value
      jz  = jy + 1    ! Pointer to findpts input z value
      jv0 = jz + 1    ! particle velocity at this timestep
      ju0 = jv0 + 3   ! fluid velocity at this time step
      jf0 = ju0 + 3   ! particle total force at this timestep
      jq0 = jf0 + 3   ! temperature forcing
      jg0 = jq0 + 1   ! work done by forces
      jquu= jg0 + 1   ! undisturbed unsteady temp. forcing
      jqqs= jquu+ 1   ! quasi-steady temp. forcing

c     forcing
      ii  = jqqs + 1
c     if (part_force(1).ne.0) then ! user specified force
         jfusr = ii
         ii    = ii + 3
c     endif
c     if (part_force(2).ne.0) then ! quasi-steady force
         jfqs  = ii
         ii    = ii + 3
c     endif
c     if (part_force(3).ne.0) then ! undisturbed force
         jfun  = ii
         ii    = ii + 3
c     endif
c     if (part_force(4).ne.0) then ! inviscid unsteady force
         jfiu  = ii
         ii    = ii + 3
c     endif

      jfcol  = ii  ! collisional force

c     other parameters (some may not be used; all at part. location)
      jtaup   = jfcol   + 3 ! particle time scale
      jcd     = jtaup   + 1 ! drag coeff
      jdrhodt = jcd     + 3 ! density material time derivative
      jre     = jdrhodt + 1 ! Relative Reynolds number
      jDuDt   = jre     + 1 ! fluid velocity time derivative
      jtemp   = jDuDt   + 3 ! part. temperature 
      jtempf  = jtemp   + 1 ! fluid temperature at part. loc.
      jrho    = jtempf  + 1 ! fluid denisty 
      jrhop   = jrho    + 1 ! particle material density
      ja      = jrhop   + 1 ! fluid mach number
      jvol    = ja      + 1 ! particle volume 
      jvol1   = jvol    + 1 ! particle volume fraction at part. loc.
      jdp     = jvol1   + 1 ! particle diameter
      jrpe    = jdp     + 1 ! particle effective diameter spl
      jgam    = jrpe    + 1 ! spread to grid correction
      jspl    = jgam    + 1 ! super particle loading
      jcmiu   = jspl    + 1 ! added mass coefficient

c     bdf/ext integration
      jx1 = jcmiu+1 ! Pointer to xyz at t^{n-1}
      jx2 = jx1 +3 ! Pointer to xyz at t^{n-1}
      jx3 = jx2 +3 ! Pointer to xyz at t^{n-1}

      jv1 = jx3+ 3 ! Pointer to particle velocity at t^{n-1}
      jv2 = jv1+ 3 ! Pointer to particle velocity at t^{n-2}
      jv3 = jv2+ 3 ! Pointer to particle velocity at t^{n-3}

      ju1 = jv3+ 3 ! Pointer to fluid velocity at t^{n-1}
      ju2 = ju1+ 3 ! Pointer to fluid velocity at t^{n-2}
      ju3 = ju2+ 3 ! Pointer to fluid velocity at t^{n-3}

      jar = ju3+ 3 ! Pointer to auxiliary reals

      nar = nr - (jar-1)  ! Number of auxiliary reals
      if (nar.le.0) call exitti('Error in nar:$',nr)

c     ghost particle integer pointers -------------------------------
      jgpiic  = 1 ! if gp used in collisions or not
      jgpps   = 2 ! Pointer to proc id for data swap
      jgppt   = 3 ! findpts return processor id
      jgpes   = 4 ! Destination element to be sent to

c     ghost particle real pointers ----------------------------------
      jgpx    = 1 ! ghost particle xloc
      jgpy    = 2 ! ghost particle yloc
      jgpz    = 3 ! ghost particle zloc
      jgpfh   = 4 ! ghost particle hydrodynamic xforce (i+1 > y, i+2 > z)
      jgpvol  = jgpfh+3  ! ghost particle volume
      jgprpe  = jgpvol+1  ! ghost particle effective radius
      jgpspl  = jgprpe+1 ! spreading correction (if used)
      jgpg0   = jgpspl+1 ! 
      jgpq0   = jgpg0 +1 ! 
      jgpv0   = jgpq0 +1 ! velocity (3 components)

      return
      end
c----------------------------------------------------------------------
c     particle force routines
c----------------------------------------------------------------------
      subroutine usr_particles_solver
c
c     call routines in ordered way - main solver structure
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      logical ifinject
      integer icalld
      save    icalld
      data    icalld  /-1/

      if (icalld .eq. -1) then
         pttime(1) = 0.
      else
         pttime(1) = pttime(1) + dnekclock() - ptdum(1)
      endif

      icalld = icalld + 1

c     should we inject particles at this time step?
      ifinject = .false.
      if (inject_rate .gt. 0) then
      if ((mod(istep,inject_rate).eq.0)) then 
         ifinject = .true. 
      endif
      endif

      if (istep .gt. time_delay) then

      if (stage.eq.1) then
         ! Update coordinates if particle moves outside boundary
         ptdum(2) = dnekclock()
            call update_particle_location  
         pttime(2) = pttime(2) + dnekclock() - ptdum(2)

         ! Inject particles if needed
         if (ifinject) call place_particles

         ! Update where particle is stored at
         ptdum(3) = dnekclock()
            call move_particles_inproc
         pttime(3) = pttime(3) + dnekclock() - ptdum(3)

         if (two_way.gt.1) then
            ! Create ghost/wall particles
            ptdum(4) = dnekclock()
               call create_extra_particles
               call sort_local_particles
            pttime(4) = pttime(4) + dnekclock() - ptdum(4)

            ! Send ghost particles
            ptdum(5) = dnekclock()
               call send_ghost_particles
            pttime(5) = pttime(5) + dnekclock() - ptdum(5)
   
            ! Projection to Eulerian grid
            ptdum(6) = dnekclock()
               call spread_props_grid
            pttime(6) = pttime(6) + dnekclock() - ptdum(6)
   
         endif
      endif

      ! Interpolate Eulerian properties to particle location
      ptdum(7) = dnekclock()
         call interp_props_part_location
      pttime(7) = pttime(7) + dnekclock() - ptdum(7)

      ! Evaluate particle force models
      ptdum(8) = dnekclock()
         call usr_particles_forcing  
      pttime(8) = pttime(8) + dnekclock() - ptdum(8)

      ! Integrate in time
      ptdum(9) = dnekclock()
         if (abs(time_integ) .eq. 1) call rk3_integrate
         if (abs(time_integ) .eq. 2) call bdf_integrate
      pttime(9) = pttime(9) + dnekclock() - ptdum(9)

      ! Update forces
      ptdum(10) = dnekclock()
         call compute_forcing_post_part
      pttime(10) = pttime(10) + dnekclock() - ptdum(10)

      endif ! time_delay

      ptdum(1) = dnekclock()

      return
      end
c----------------------------------------------------------------------
c     subroutine correct_spl
c
c     correct initial super particle loading
c
c     include 'SIZE'
c     include 'INPUT'
c     include 'GEOM'
c     include 'SOLN'
c     include 'CMTDATA'
c     include 'MASS'
c     include 'CMTPART'

c     real rdumvol(llpart,2*3)

c     ! this routine needs updating!

c     do i=1,n
c        rdumvol(i,1) = rpart(jvol,i)  ! particle volume
c        rdumvol(i,2) = rpart(jvol1,i) ! interp vol frac @ part loc
c        rdumvol(i,3) = rpart(jspl,i)  ! super part. loading
c     enddo

c     call usr_particles_io

c     begin diagnostics ----
c
c     eulerian volume frac 
c     nxyze = nx1*ny1*nz1*nelt
c     rmu1  = glsc2(bm1,ptw(1,1,1,1,4),nxyze)
c     rmu1  = rmu1/vol_distrib

c
c     lagrangian volume frac
c     rmu2  = glsum(rdumvol(1,2),n)
c     rmu2  = rmu2/nw
c     rmin2 = glmin(rdumvol(1,2),n)
c     rmax2 = glmax(rdumvol(1,2),n)

c
c     what spl mean should be
c     rdumt   = glsum(rdumvol(1,1),n)
c     rsplavg = phi_desire*vol_distrib/rdumt

c
c     what spl mean actually is
c     rmu3  = glsum(rdumvol(1,3),n)
c     rmu3  = rmu3/nw
c     rmin3 = glmin(rdumvol(1,3),n)
c     rmax3 = glmax(rdumvol(1,3),n)

c
c     variance and skew stuff
c     do i=1,n
c        rdumvol(i,4) = (rpart(jvol1,i) - rmu2)**2
c        rdumvol(i,5) = (rpart(jspl,i) - rmu3)**2
c     enddo

c     rvar2 = glsum(rdumvol(1,4),n)
c     rvar2 = rvar2/nw

c     rvar3 = glsum(rdumvol(1,5),n)
c     rvar3 = rvar3/nw

c     if (nid.eq.0) write(6,*) '-DZ- Md,bd,Me --'
c     if (nid.eq.0) write(6,*) phi_desire,rsplavg,rmu1
c     if (nid.eq.0) write(6,*) '-DZ-- Ml,Sl,Minl,Maxl'
c     if (nid.eq.0) write(6,*) rmu2,sqrt(rvar2),rmin2,rmax2
c     if (nid.eq.0) write(6,*) '-DZ--- Mb,Sb,Minb,Maxb'
c     if (nid.eq.0) write(6,*) rmu3,sqrt(rvar3),rmin3,rmax3
c     end diagnostics ----


c     do ip=1,n
c        linear roll off near bed edges
c        rthresh = 2.*rleng
c        if (rpart(jx,ip) .lt. rxbo(1,1)+rthresh) then
c           rxmin = rxbo(1,1)
c           rxmax = rxbo(1,1)+rthresh
c           rm    = (0.-phi_desire)/(rxmin-rxmax)
c           phi_val = 0. + rm*(rpart(jx,ip)-rxmin)
c        elseif (rpart(jx,ip) .gt. rxbo(2,1)-rthresh) then
c           rxmin = rxbo(2,1)-rthresh
c           rxmax = rxbo(2,1)
c           rm    = (phi_desire-0.)/(rxmin-rxmax)
c           phi_val = phi_desire + rm*(rpart(jx,ip)-rxmin)
c        else
c           phi_val = phi_desire
c        endif

c        tanh roll off near bed edges ! do this one if bed edge!!
c        rthresh = 4.*rleng
c        rxboavg = (rxbo(2,1)+rxbo(1,1))/2.
c        if (rpart(jx,ip) .lt. rxboavg) then
c           ra = rxbo(1,1)
c           rb = rxbo(1,1)+rthresh
c           rap=0.
c           rbp=phi_desire
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           phi_val = rap + rd*(0.5 + 0.5*tanh(200.*(rpart(jx,ip)-rc)))
c        elseif (rpart(jx,ip) .gt. rxboavg) then
c           ra = rxbo(2,1)-rthresh
c           rb = rxbo(2,1)
c           rap=phi_desire
c           rbp=0.
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           phi_val = rap + rd*(0.5 + 0.5*tanh(200.*(rpart(jx,ip)-rc)))
c        else
c           phi_val = phi_desire
c        endif

c        exponential roll off
c        ryyl = 0.26
c        ryyr = 0.28
c        ry = rpart(jy,ip)
c        if (ry .le. ryyl) then
c           phi_val = phi_desire
c        else
c           rssig = sqrt(-(ryyr - ryyl)**2/2./log(0.01))
c           phi_val = phi_desire*exp(-(ry-ryyl)**2/2./rssig**2)
c        endif

c        phi_val = phi_desire ! comment out if bed and uncomment above
c        rtmp = 0.30*rpart(jspl,ip)
c        rxi = rtmp*(1. - rpart(jvol1,ip)/phi_val)
c        rpart(jspl,ip)=rpart(jspl,ip) + rxi
c        if (rpart(jspl,ip).lt.0.) rpart(jspl,ip) = 0.

c        rthresh = 2.*rleng
c        if (rpart(jx,ip) .lt. rxbo(1,1)+rthresh) then
c           rxmin = rxbo(1,1)
c           rxmax = rxbo(1,1)+rthresh
c           rm    = (1.-rsplavg)/(rxmin-rxmax)
c           rpart(jspl,ip) = 1. + rm*(rpart(jx,ip)-rxmin)

c           ra = rxbo(1,1)
c           rb = rxbo(1,1)+rthresh
c           rap=1.
c           rbp=rsplavg
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           rpart(jspl,ip) = rbp + rd*(0.5 + 
c    >                 0.5*tanh(100.*(rb-ra)/rc*(rpart(jx,ip)-rc)))
c        elseif (rpart(jx,ip) .gt. rxbo(2,1)-rthresh) then
c           rxmin = rxbo(2,1)-rthresh
c           rxmax = rxbo(2,1)
c           rm    = (rsplavg-1.)/(rxmin-rxmax)
c           rpart(jspl,ip) = rsplavg + rm*(rpart(jx,ip)-rxmin)

c           ra = rxbo(2,1)-rthresh
c           rb = rxbo(2,1)
c           rap=rsplavg
c           rbp=1.
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           rpart(jspl,ip) = rbp + rd*(0.5 + 
c    >                 0.5*tanh(100.*(rb-ra)/rc*(rpart(jx,ip)-rc)))
c        else
c           rtmp = 0.30*rpart(jspl,ip)
c           rxi = rtmp*(1. - rpart(jvol1,ip)/phi_desire)
c           rpart(jspl,ip)=rpart(jspl,ip) + rxi
c        endif
c        if (rpart(jspl,ip).lt.0.) rpart(jspl,ip) = 0.


c        linear roll off near bed edges
c        rthresh = 1.*rleng
c        if (rpart(jx,ip) .lt. rxbo(1,1)+rthresh) then
c           rxmin = rxbo(1,1)
c           rxmax = rxbo(1,1)+rthresh
c           rm    = (0.-phi_desire)/(rxmin-rxmax)
c           phi_val = 0. + rm*(rpart(jx,ip)-rxmin)
c        elseif (rpart(jx,ip) .gt. rxbo(2,1)-rthresh) then
c           rxmin = rxbo(2,1)-rthresh
c           rxmax = rxbo(2,1)
c           rm    = (phi_desire-0.)/(rxmin-rxmax)
c           phi_val = phi_desire + rm*(rpart(jx,ip)-rxmin)
c        else
c           phi_val = phi_desire
c        endif

c        tanh roll off near bed edges
c        rthresh = 2.*rleng
c        rxboavg = (rxbo(2,1)+rxbo(1,1))/2.
c        if (rpart(jx,ip) .lt. rxboavg) then
c           ra = rxbo(1,1)
c           rb = rxbo(1,1)+rthresh
c           rap=0.
c           rbp=phi_desire
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           phi_val = rbp + rd*(0.5 + 
c    >                 0.5*tanh(100.*(rb-ra)/rc*(rpart(jx,ip)-rc)))
c        elseif (rpart(jx,ip) .gt. rxboavg) then
c           ra = rxbo(2,1)-rthresh
c           rb = rxbo(2,1)
c           rap=phi_desire
c           rbp=0.
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           phi_val = rbp + rd*(0.5 + 
c    >                 0.5*tanh(100.*(rb-ra)/rc*(rpart(jx,ip)-rc)))
c        else
c           phi_val = phi_desire
c        endif

c1511 continue
c     enddo

c     return
c     end
c----------------------------------------------------------------------
      subroutine spread_props_grid
c
c     spread particle properties at fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'MASS'
      include 'TSTEP'
      include 'CMTPART'

      real    xx,yy,zz,vol,pfx,pfy,pfz,pmass,pmassf,vcell,multfc
     >       ,qgqf,rvx,rvy,rvz,rcountv(8,nelt)
      integer e

      nlxyze = lx1*ly1*lz1*lelt
      nxyze  = nx1*ny1*nz1*nelt
      call rzero(ptw,nlxyze*8)

      ! do nothing
      if (abs(npro_method) .eq. 0) then

      ! finite volume-like spreading
      elseif (abs(npro_method) .eq. 1) then
      call rzero(rcountv,8*nelt)

c     ! for grid spreading line in finite volume.. no ghost particles
c     local mpi rank effects
      do ip=1,n
         e   = ipart(je0,ip) + 1
         xx  = rpart(jx,ip)
         yy  = rpart(jy,ip)
         zz  = rpart(jz,ip)

         pfx = -rpart(jf0,ip)
         pfy = -rpart(jf0+1,ip)
         pfz = -rpart(jf0+2,ip)
         vol = rpart(jvol,ip)
         qgqf= -(rpart(jg0,ip) + rpart(jq0,ip))
         rvx = rpart(jv0  ,ip)*vol
         rvy = rpart(jv0+1,ip)*vol
         rvz = rpart(jv0+2,ip)*vol

         rcountv(1,e) = rcountv(1,e) + pfx
         rcountv(2,e) = rcountv(2,e) + pfy
         rcountv(3,e) = rcountv(3,e) + pfz
         rcountv(4,e) = rcountv(4,e) + vol
         rcountv(5,e) = rcountv(5,e) + qgqf
         rcountv(6,e) = rcountv(6,e) + rvx
         rcountv(7,e) = rcountv(7,e) + rvy
         rcountv(8,e) = rcountv(8,e) + rvz
      enddo

         call local_part_to_grid_fv(ptw(1,1,1,1,1),ptw(1,1,1,1,2),
     >                              ptw(1,1,1,1,3),ptw(1,1,1,1,4),
     >                              ptw(1,1,1,1,5),ptw(1,1,1,1,6),
     >                              ptw(1,1,1,1,7),ptw(1,1,1,1,8),
     >                              rcountv)


      ! gaussian spreading
      elseif (abs(npro_method).eq.2) then


      pi       = 4.0d+0*atan(1.0d+0)
      rbexpi   = 1./(-2.*rsig**2)

      ralph    = dfilt/2.*sqrt(-log(ralphdecay)/log(2.))*rleng
      multfc   = 1./(sqrt(2.*pi)**3 * rsig**3) ! exponential
      ralph2   = ralph**2

c     local mpi rank effects
      do ip=1,n
         e   = ipart(je0,ip) + 1
         xx  = rpart(jx,ip)
         yy  = rpart(jy,ip)
         zz  = rpart(jz,ip)

         pfx = -rpart(jf0,ip)*multfc
         pfy = -rpart(jf0+1,ip)*multfc
         pfz = -rpart(jf0+2,ip)*multfc
         vol = rpart(jvol,ip)*multfc
         qgqf= -(rpart(jg0,ip) + rpart(jq0,ip))*multfc
         rvx = rpart(jv0  ,ip)*vol
         rvy = rpart(jv0+1,ip)*vol
         rvz = rpart(jv0+2,ip)*vol

         call local_part_to_grid(ptw(1,1,1,1,1),ptw(1,1,1,1,2),
     >                           ptw(1,1,1,1,3),ptw(1,1,1,1,4),
     >                           ptw(1,1,1,1,5),ptw(1,1,1,1,6),
     >                           ptw(1,1,1,1,7),ptw(1,1,1,1,8),
     >                           pfx,pfy,pfz,vol,qgqf,rvx,rvy,rvz,
     >                           xx,yy,zz,rbexpi,
     >                           ralph,ralph2,e)
      enddo

c     ntmp1 = iglmax(n,1)
c     ntmp2 = iglmax(nfptsgp,1)
c     if (nid.eq.0) write(6,*) 'Passed local spreading to grid'
c    >                                                     ,ntmp1,ntmp2

c     remote mpi rank effects
      do ip=1,nfptsgp
         e   = iptsgp(jgpes,ip) + 1
         xx  = rptsgp(jgpx,ip)
         yy  = rptsgp(jgpy,ip)
         zz  = rptsgp(jgpz,ip)

         pfx = -rptsgp(jgpfh,ip)*multfc
         pfy = -rptsgp(jgpfh+1,ip)*multfc
         pfz = -rptsgp(jgpfh+2,ip)*multfc
         vol = rptsgp(jgpvol,ip)*multfc
         qgqf= -(rptsgp(jgpg0,ip) + rptsgp(jgpq0,ip))*multfc
         rvx = rptsgp(jgpv0  ,ip)*vol
         rvy = rptsgp(jgpv0+1,ip)*vol
         rvz = rptsgp(jgpv0+2,ip)*vol

         if (iptsgp(jgpiic,ip) .ne. 0) then
            call remote_part_to_grid(ptw(1,1,1,1,1),ptw(1,1,1,1,2),
     >                               ptw(1,1,1,1,3),ptw(1,1,1,1,4),
     >                               ptw(1,1,1,1,5),ptw(1,1,1,1,6),
     >                               ptw(1,1,1,1,7),ptw(1,1,1,1,8),
     >                               pfx,pfy,pfz,vol,qgqf,rvx,rvy,rvz,
     >                               xx,yy,zz,rbexpi,
     >                               ralph,ralph2,e)
         endif
      enddo

      endif

c     ntmp = iglsum(n,1)
c     if (nid.eq.0) write(6,*) 'Passed remote spreading to grid'

c     wght = 1.0
c     ncut = 1
c     do i=1,8
c        call filter_s0(ptw(1,1,1,1,i),wght,ncut,'ptw') 
c     enddo

      wght = 1.0
      ncut = 1
      call filter_s0(ptw(1,1,1,1,4),wght,ncut,'phip') 

c     rvfmax = 0.7
      rvfmin = 0.0
      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
c        if (ptw(i,j,k,ie,4) .gt. rvfmax) ptw(i,j,k,ie,4) = rvfmax
         if (ptw(i,j,k,ie,4) .lt. rvfmin) ptw(i,j,k,ie,4) = rvfmin
         phig(i,j,k,ie) = 1. - ptw(i,j,k,ie,4)
c        rhs_fluidp(i,j,k,ie,8) = phig(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo


c     do ie=1,nelt
c     do k=1,nz1
c     do j=1,ny1
c     do i=1,nx1
c        if (rhs_fluidp(i,j,k,ie,8) .gt. 1.0) rhs_fluidp(i,j,k,ie,8)=1.0
c     enddo
c     enddo
c     enddo
c     enddo

      
      return
      end
c----------------------------------------------------------------------
      subroutine local_part_to_grid_fv(fvalgx,fvalgy,fvalgz,fvalgv,
     >                              fvalgg,fvalv1,fvalv2,fvalv3,
     >                              rcountv)
c
c     spread a local particle property to local fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      integer e,er
      real    fvalgx(nx1,ny1,nz1,nelt),fvalgy(nx1,ny1,nz1,nelt),
     >        fvalgz(nx1,ny1,nz1,nelt),fvalgv(nx1,ny1,nz1,nelt),
     >        fvalgg(nx1,ny1,nz1,nelt),fvalv1(nx1,ny1,nz1,nelt),
     >        fvalv2(nx1,ny1,nz1,nelt),fvalv3(nx1,ny1,nz1,nelt),
     >        rcountv(8,nelt)

         do ie=1,nelt

            rvole = (xerange(2,1,ie) - xerange(1,1,ie))*
     >              (xerange(2,2,ie) - xerange(1,2,ie))*
     >              (xerange(2,3,ie) - xerange(1,3,ie))
            rvolei=1./rvole
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            fvalgx(i,j,k,ie) = rcountv(1,ie)*rvolei
            fvalgy(i,j,k,ie) = rcountv(2,ie)*rvolei
            fvalgz(i,j,k,ie) = rcountv(3,ie)*rvolei
            fvalgv(i,j,k,ie) = rcountv(4,ie)*rvolei
            fvalgg(i,j,k,ie) = rcountv(5,ie)*rvolei
            fvalv1(i,j,k,ie) = rcountv(6,ie)*rvolei
            fvalv2(i,j,k,ie) = rcountv(7,ie)*rvolei
            fvalv3(i,j,k,ie) = rcountv(8,ie)*rvolei
         enddo
         enddo
         enddo
         enddo

      return
      end
c----------------------------------------------------------------------
      subroutine local_part_to_grid(fvalgx,fvalgy,fvalgz,fvalgv,fvalgg,
     >                              fvalv1,fvalv2,fvalv3,
     >                              pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                              ppv1,ppv2,ppv3,
     >                              xx,yy,zz,rbexpi,ralph,ralph2,e)
c
c     spread a local particle property to local fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      common /gpfix/ ilgp_f(lelt,6),ilgp_e(lelt,12),ilgp_c(lelt,8)

      integer e,er
      real    fvalgx(lx1,ly1,lz1,lelt),fvalgy(lx1,ly1,lz1,lelt),
     >        fvalgz(lx1,ly1,lz1,lelt),fvalgv(lx1,ly1,lz1,lelt),
     >        fvalgg(lx1,ly1,lz1,lelt),fvalv1(lx1,ly1,lz1,lelt),
     >        fvalv2(lx1,ly1,lz1,lelt),fvalv3(lx1,ly1,lz1,lelt),
     >        pvalpx,pvalpy,pvalpz,pvalpv,xx,yy,zz,ppg,ppv1,ppv2,ppv3

c     this element
      call point_to_grid(fvalgx(1,1,1,e),fvalgy(1,1,1,e),
     >                   fvalgz(1,1,1,e),fvalgv(1,1,1,e),
     >                   fvalgg(1,1,1,e),fvalv1(1,1,1,e),
     >                   fvalv2(1,1,1,e),fvalv3(1,1,1,e),
     >                   xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)

c     faces
      do ii=1,nfacegp
         er=el_face_el_map(e,ii) + 1
         impi=el_face_proc_map(e,ii)
         if (ilgp_f(e,ii) .eq. 0) then
         if (impi .eq. nid) then
            call point_to_grid(fvalgx(1,1,1,er),fvalgy(1,1,1,er),
     >                   fvalgz(1,1,1,er),fvalgv(1,1,1,er),
     >                   fvalgg(1,1,1,er),fvalv1(1,1,1,er),
     >                   fvalv2(1,1,1,er),fvalv3(1,1,1,er),
     >                   xm1(1,1,1,er),ym1(1,1,1,er),zm1(1,1,1,er),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)
         endif
         endif
      enddo

c     edges
      do ii=1,nedgegp
         er=el_edge_el_map(e,ii) + 1
         impi=el_edge_proc_map(e,ii)
         if (ilgp_e(e,ii) .eq. 0) then
         if (impi .eq. nid) then
            call point_to_grid(fvalgx(1,1,1,er),fvalgy(1,1,1,er),
     >                   fvalgz(1,1,1,er),fvalgv(1,1,1,er),
     >                   fvalgg(1,1,1,er),fvalv1(1,1,1,er),
     >                   fvalv2(1,1,1,er),fvalv3(1,1,1,er),
     >                   xm1(1,1,1,er),ym1(1,1,1,er),zm1(1,1,1,er),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)
         endif
         endif
      enddo

c     corners
      do ii=1,ncornergp
         er=el_corner_el_map(e,ii) + 1
         impi=el_corner_proc_map(e,ii)
         if (ilgp_c(e,ii) .eq. 0) then
         if (impi .eq. nid)  then
            call point_to_grid(fvalgx(1,1,1,er),fvalgy(1,1,1,er),
     >                   fvalgz(1,1,1,er),fvalgv(1,1,1,er),
     >                   fvalgg(1,1,1,er),fvalv1(1,1,1,er),
     >                   fvalv2(1,1,1,er),fvalv3(1,1,1,er),
     >                   xm1(1,1,1,er),ym1(1,1,1,er),zm1(1,1,1,er),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)
         endif
         endif
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine remote_part_to_grid(fvalgx,fvalgy,fvalgz,fvalgv,fvalgg,
     >                                fvalv1,fvalv2,fvalv3,
     >                                pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                                ppv1,ppv2,ppv3,
     >                                xx,yy,zz,rbexpi,ralph,ralph2,e)
c
c     spread a remote particle property to local fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      integer e
      real    fvalgx(lx1,ly1,lz1,lelt),fvalgy(lx1,ly1,lz1,lelt),
     >        fvalgz(lx1,ly1,lz1,lelt),fvalgv(lx1,ly1,lz1,lelt),
     >        fvalgg(lx1,ly1,lz1,lelt),fvalv1(lx1,ly1,lz1,lelt),
     >        fvalv2(lx1,ly1,lz1,lelt),fvalv3(lx1,ly1,lz1,lelt),
     >        pvalpx,pvalpy,pvalpz,pvalpv,xx,yy,zz,ppg,
     >        ppv1,ppv2,ppv3

      call point_to_grid(fvalgx(1,1,1,e),fvalgy(1,1,1,e),
     >                   fvalgz(1,1,1,e),fvalgv(1,1,1,e),
     >                   fvalgg(1,1,1,e),fvalv1(1,1,1,e),
     >                   fvalv2(1,1,1,e),fvalv3(1,1,1,e),
     >                   xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)

      return
      end
c----------------------------------------------------------------------
      subroutine point_to_grid(gval1,gval2,gval3,gval4,gval5,
     >                      gval6,gval7,gval8,
     >                      xgd,ygd,zgd,
     >                      pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                      ppv1,ppv2,ppv3,
     >                      xx,yy,zz,rbexpi,
     >                      ralph,ralph2)
c
c     spreads point onto grid in element e
c       gval: grid value
c       pval: particle value
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      integer i,j,k,e,ip
      real    gval1(lx1,ly1,lz1),gval2(lx1,ly1,lz1),gval3(lx1,ly1,lz1),
     >        gval4(lx1,ly1,lz1),gval5(lx1,ly1,lz1),gval6(lx1,ly1,lz1),
     >        gval7(lx1,ly1,lz1),gval8(lx1,ly1,lz1),
     >        xgd(lx1,ly1,lz1),ygd(lx1,ly1,lz1),zgd(lx1,ly1,lz1),
     >        pvalpx,pvalpy,pvalpz,pvalpv,pi,
     >        distx,disty,distz,xx,yy,zz,distx2,disty2,distz2,multfc,
     >        ppg,ppv1,ppv2,ppv3

c     optimized code ------------------------------------------------
c     can we skip this entire element?
c     rxl      = abs(xx - xgd(1,1,1))
c     rxr      = abs(xx - xgd(nx1,1,1))
c     if (rxl.gt.ralph .and. rxr.gt.ralph) goto 1514
c     rxl      = abs(yy - ygd(1,1,1))
c     rxr      = abs(yy - ygd(1,ny1,1))
c     if (rxl.gt.ralph .and. rxr.gt.ralph) goto 1514
c     rxl      = abs(zz - zgd(1,1,1))
c     rxr      = abs(zz - zgd(1,1,nz1))
c     if (rxl.gt.ralph .and. rxr.gt.ralph) goto 1514

      do k=1,nz1
         rtmp = 0.
         distz = zz - zgd(1,1,k) ! even element spacing only!
         distz2 = distz**2
! DZ
         if (distz2 .gt. ralph2) goto 1513
      do j=1,ny1
         disty = yy - ygd(1,j,1) ! even element spacing only!
         disty2 = disty**2
         rtmp1   =  distz2 + disty2
! DZ
         if (rtmp1 .gt. ralph2) goto 1512
      do i=1,nx1
         distx = xx - xgd(i,1,1) ! even element spacing only!
         distx2 = distx**2
         rtmp2 = rtmp1 + distx2
! DZ
         if (rtmp2 .gt. ralph2) goto 1511

         rdum = rtmp2*rbexpi
         rexp = exp(rdum)

         gval1(i,j,k) = gval1(i,j,k) + pvalpx*rexp
         gval2(i,j,k) = gval2(i,j,k) + pvalpy*rexp
         gval3(i,j,k) = gval3(i,j,k) + pvalpz*rexp
         gval4(i,j,k) = gval4(i,j,k) + pvalpv*rexp
         gval5(i,j,k) = gval5(i,j,k) + ppg   *rexp
         gval6(i,j,k) = gval6(i,j,k) + ppv1  *rexp
         gval7(i,j,k) = gval7(i,j,k) + ppv2  *rexp
         gval8(i,j,k) = gval8(i,j,k) + ppv3  *rexp
 1511 continue
      enddo
 1512 continue
      enddo
 1513 continue
      enddo
 1514 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine rk3_integrate
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      common /PARTRK3/ kv_stage_p
      real   kv_stage_p(llpart,7)

      integer fdim
      real    pmass

      jx0 = jx

c     rk3 stage one items ---------------------------------------------
      if (stage.eq.1) then
c        used for time derivative of v in iu force
         call get_bdf_ext_coefs(beta,alpha,times)

c        move data to previous positions
         do j=0,ndim-1
         do i=1,n
            rpart(ju3+j,i)=rpart(ju2+j,i)
            rpart(ju2+j,i)=rpart(ju1+j,i)
            rpart(ju1+j,i)=rpart(ju0+j,i)
            rpart(jv3+j,i)=rpart(jv2+j,i)
            rpart(jv2+j,i)=rpart(jv1+j,i)
            rpart(jv1+j,i)=rpart(jv0+j,i)
            rpart(jx3+j,i)=rpart(jx2+j,i)
            rpart(jx2+j,i)=rpart(jx1+j,i)
            rpart(jx1+j,i)=rpart(jx0+j,i)
         enddo
         enddo

         do i=1,n
            kv_stage_p(i,1) = rpart(jx0  ,i)
            kv_stage_p(i,2) = rpart(jx0+1,i)
            kv_stage_p(i,3) = rpart(jx0+2,i)
            kv_stage_p(i,4) = rpart(jv0  ,i)
            kv_stage_p(i,5) = rpart(jv0+1,i)
            kv_stage_p(i,6) = rpart(jv0+2,i)
            kv_stage_p(i,7) = rpart(jtemp,i)
         enddo
      endif

c     all rk3 stages items --------------------------------------------
      do i=1,n
         rpart(jv0  ,i) = tcoef(1,stage)*kv_stage_p(i,4) +
     >                    tcoef(2,stage)*rpart(jv0  ,i)  +
     >                    tcoef(3,stage)*rpart(jf0  ,i)
         rpart(jv0+1,i) = tcoef(1,stage)*kv_stage_p(i,5) +
     >                    tcoef(2,stage)*rpart(jv0+1,i)  +
     >                    tcoef(3,stage)*rpart(jf0+1,i)
         rpart(jv0+2,i) = tcoef(1,stage)*kv_stage_p(i,6) +
     >                    tcoef(2,stage)*rpart(jv0+2,i)  +
     >                    tcoef(3,stage)*rpart(jf0+2,i)
         rpart(jx0  ,i) = tcoef(1,stage)*kv_stage_p(i,1) +
     >                    tcoef(2,stage)*rpart(jx0  ,i)  +
     >                    tcoef(3,stage)*rpart(jv0  ,i)
         rpart(jx0+1,i) = tcoef(1,stage)*kv_stage_p(i,2) +
     >                    tcoef(2,stage)*rpart(jx0+1,i)  +
     >                    tcoef(3,stage)*rpart(jv0+1,i)
         rpart(jx0+2,i) = tcoef(1,stage)*kv_stage_p(i,3) +
     >                    tcoef(2,stage)*rpart(jx0+2,i)  +
     >                    tcoef(3,stage)*rpart(jv0+2,i)
         rpart(jtemp,i) = tcoef(1,stage)*kv_stage_p(i,7) +
     >                    tcoef(2,stage)*rpart(jtemp,i)  +
     >                    tcoef(3,stage)*rpart(jq0  ,i)

c     if (i .eq. 1) then
c     if (nid.eq.700) then
c        write(6,'(A,10F15.10)')'uu',rpart(jy,i)
c    >                             ,rpart(jf0+1,i)
c    >                             ,rpart(ju0+1,i)
c    >                             ,rpart(jv0+1,i)
c    >                             ,rpart(jfqs+1,i)
c    >                             ,rpart(jfcol+1,i)
c    >                             ,rpart(jfusr+1,i)
c    >                             ,rpart(jtaup,i)
c    >                             ,rpart(jvol,i)
c    >                             ,rpart(jdp,i)

c     endif
c     endif

c     write(6,*) 'forcing',rpart(jx,i),rpart(jv0,i),rpart(jf0,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_forcing
c
c     calculate the rhs of particle equation
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'
      include 'PERFECTGAS'

      real vel_diff,pmass,pmassf

      common /PARTRK3/ kv_stage_p
      real kv_stage_p(llpart,7)

      real MixtPerf_C_GRT_part
      external MixtPerf_C_GRT_part

      pi  = 4.0d+0*atan(1.0d+0)

      if ((part_force(2) .ne. 0) .or. (part_force(3) .ne. 0))
     >   call calc_substantial_derivative

      do i=1,n
         rpart(jfcol  ,i) = 0.
         rpart(jfcol+1,i) = 0.
         rpart(jfcol+2,i) = 0.
      enddo

      do i=1,n
c        setup values ------------------------------------------------
         pmass = rpart(jvol,i)*rpart(jrhop,i)
         pmassf= rpart(jvol,i)*rpart(jrho,i)
         if(part_force(3).ne.0) pmass = pmass + rpart(jcmiu,i)*pmassf ! am

         vel_diff = sqrt((rpart(ju0  ,i)-rpart(jv0  ,i))**2+
     >                   (rpart(ju0+1,i)-rpart(jv0+1,i))**2+
     >                   (rpart(ju0+2,i)-rpart(jv0+2,i))**2)
         rpart(ja,i)  = MixtPerf_C_GRT_part(gmaref,rgasref,
     >                           rpart(jtempf,i),icmtp)
         rpart(ja,i)  = vel_diff/rpart(ja,i) ! relative mach number
         rpart(jre,i) = rpart(jrho,i)*rpart(jdp,i)*vel_diff/mu_0 ! Re

c        momentum rhs ------------------------------------------------

         call usr_particles_f_col(i)     ! colision force all at once

         call usr_particles_f_user(i)    ! user/body force
            rpart(jfusr+0,i) = f_part(0)
            rpart(jfusr+1,i) = f_part(1)
            rpart(jfusr+2,i) = f_part(2)

         do j=0,ndim-1
            call usr_particles_f_qs(i,j)
            call usr_particles_f_un(i,j)
            call usr_particles_f_iu(i,j)

            rdum = 0.
            rdum = rdum + rpart(jfusr+j,i)
            rdum = rdum + rpart(jfqs+j,i)
            rdum = rdum + rpart(jfun+j,i)
            rdum = rdum + rpart(jfiu+j,i)
            rdum = rdum + rpart(jfcol+j,i)

            rpart(jf0+j,i) = rdum/pmass ! mass weighted force
         enddo

c        energy rhs --------------------------------------------------
         call usr_particles_q_uu(i)
         call usr_particles_q_qs(i)

         rdum = 0. 
         rdum = rdum + rpart(jquu,i)
         rdum = rdum + rpart(jqqs,i)

         pmass = rpart(jvol,i)*rpart(jrhop,i)
         rpart(jq0,i) = rdum/(pmass*cp_p)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_forcing_post_part
c
c     post calculate forces due to factoring of equations
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real uvel(0:2), vvel(0:2), pmass, pmassf,S_qs

      common /PARTRK3/ kv_stage_p
      real kv_stage_p(llpart,7)

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)

      pi  = 4.0d+0*atan(1.0d+0)

      do i=1,n
         pmass = rpart(jvol,i)*rpart(jrhop,i)
         pmassf= rpart(jvol,i)*rpart(jrho,i)
         if (part_force(3).ne.0) pmass =pmass + rpart(jcmiu,i)*pmassf

c        momentum forcing to fluid
         do j=0,ndim-1
            rdum = 0.

            rdvdt = rpart(jf0+j,i) ! note already divided by Mp + am
            ram_s = rdvdt*rpart(jcmiu,i)*pmassf
            rpart(jfiu+j,i) = rpart(jfiu+j,i) - ram_s

c           note that no coupled f_un in this formulation
c           rdum = rdum + rpart(jfun+j,i)
            rdum = rdum + rpart(jfiu+j,i)
            rdum = rdum + rpart(jfqs+j,i)

            rpart(jf0+j,i) = rdum ! now the force to couple with gas
         enddo

c        energy forcing to fluid (quasi-steady)
         rpart(jg0,i) = rpart(jv0  ,i)*rpart(jfqs  ,i) + !force work
     >                  rpart(jv0+1,i)*rpart(jfqs+1,i) +
     >                  rpart(jv0+2,i)*rpart(jfqs+2,i)
         rpart(jg0,i) = rpart(jg0,i) + 
     >                  rpart(ju0  ,i)*rpart(jfiu  ,i) + !iu
     >                  rpart(ju0+1,i)*rpart(jfiu+1,i) +
     >                  rpart(ju0+2,i)*rpart(jfiu+2,i)

         rdum = 0.
         rdum = rdum + rpart(jqqs,i)

         rpart(jq0,i) = rdum

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_substantial_derivative
c 
c     calculate rhs of NS, which is just pressure gradient. used for
c     undisturbed and invisicid unsteady force
c
c     no forcing included...should it be?
c
c     rhs_fluidp(i,j,k,e,l)
c   
c        l = 1     >    dP/dx
c        l = 2     >    dP/dy
c        l = 3     >    dP/dz
c        l = 4     >    -P* div(phi_p v), v is eulerian particle vel
c        l = 5     >    P * d phi_g/dx
c        l = 6     >    P * d phi_g/dy
c        l = 7     >    P * d phi_g/dz
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      integer e
      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1),
     >         pm1(lx1,ly1,lz1,lelt,3),udum(lx1,ly1,lz1)

      nxyz=nx1*ny1*nz1
      nlxyze = lx1*ly1*lz1*lelt

      call rzero(rhs_fluidp,nlxyze*7)

      ! if pressure solved on different mesh, map to vel mesh
      if (lx2 .ne. lx1) then
         call mappr(pm1(1,1,1,1,1),pr,pm1(1,1,1,1,2),pm1(1,1,1,1,3))
      else
         call copy(pm1(1,1,1,1,1),pr(1,1,1,1),nlxyze)
      endif

c     compute grad pr
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        pm1(1,1,1,e,1),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,1) = 1.0d+0/JACM1(i,1,1,e)* !d/dx
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,2) = 1.0d+0/JACM1(i,1,1,e)* !d/dy
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,3) = 1.0d+0/JACM1(i,1,1,e)* !d/dz
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
        endif ! end 3d
      enddo

      ! div (phi_p * v)
      do e=1,nelt
        do i=1,nxyz
           udum(i,1,1) = pm1(i,1,1,e,1)*ptw(1,1,1,e,6)
        enddo
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! x dir
     >                                        udum(1,1,1),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = 1.0d+0/JACM1(i,1,1,e)* !d/dx
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo
        endif ! end 3d

        do i=1,nxyz
           udum(i,1,1) = pm1(i,1,1,e,1)*ptw(1,1,1,e,7)
        enddo
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! y dir
     >                                        udum(1,1,1),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = rhs_fluidp(i,1,1,e,4) +
     >             1.0d+0/JACM1(i,1,1,e)* !d/dy
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo
         endif

        do i=1,nxyz
           udum(i,1,1) = pm1(i,1,1,e,1)*ptw(1,1,1,e,8)
        enddo
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1), ! z dir
     >                                        udum(1,1,1),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,4) = rhs_fluidp(i,1,1,e,4) +
     >             1.0d+0/JACM1(i,1,1,e)* !d/dz
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
        endif
      enddo

      do e=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rhs_fluidp(i,j,k,e,4) = -rhs_fluidp(i,j,k,e,4)
      enddo
      enddo
      enddo
      enddo

c     compute grad phi_g
      do e=1,nelt
        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        phig(1,1,1,e),lx1,if3d)
        if(if3d) then ! 3d
            do i=1,nxyz
              rhs_fluidp(i,1,1,e,5) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,6) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo

            do i=1,nxyz
              rhs_fluidp(i,1,1,e,7) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
        endif ! end 3d
      enddo

      do e=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         rhs_fluidp(i,j,k,e,5) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,5)
         rhs_fluidp(i,j,k,e,6) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,6)
         rhs_fluidp(i,j,k,e,7) = pm1(i,j,k,e,1)*rhs_fluidp(i,j,k,e,7)
      enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_q_uu(ii)
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      if (part_force(5) .eq. 0) then
         rpart(jquu,ii) = rpart(jvol,ii)*0. ! none for now (del . q_f) = 0
      endif

      ! end timer

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_q_qs(ii)
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real nu_g,S_qs,rpra,kappa_g

      pi      = 4.0d+0*atan(1.0d+0)

      kappa_g = abs(param(8))

      nu_g    = mu_0/rpart(jrho,ii) ! kinematic viscosity
      rpra    = nu_g/kappa_g
      S_qs    = 2.*pi*kappa_g*rpart(jdp,ii)

      if (part_force(4) .gt. 0) then
         S_qs    = S_qs*(1. + 0.3*sqrt(rpart(jre,ii))*rpra**(2./3.)) !re < 500
         rpart(jqqs,ii) = S_qs*(rpart(jtempf,ii) - rpart(jtemp,ii))
      elseif (part_force(4) .eq. 0) then
         S_qs    = 0.
         rpart(jqqs,ii) = S_qs
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_un(ii,jj)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      if (part_force(2) .ne. 0) then
         rpart(jfun+jj,ii) = -rpart(jvol,ii)*rpart(jDuDt+jj,ii)
      else
         rpart(jfun+jj,ii) = 0.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_iu(ii,jj)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real S_qs

      S_qs = 0.5
      if (part_force(3) .gt. 0) then
         S_qs = S_qs*(1. + 1.8*rpart(ja,ii)**2 + 7.6*rpart(ja,ii)**4)
         rphip = rpart(jvol1,ii)
         if (rphip .gt. 0.3) rphip = 0.3
         S_qs = S_qs*(1.+2.*rphip)/(1.-rphip)

         rpart(jcmiu,ii) = S_qs
         
c        need to fix, fake, no D(rho_f u)/Dt contribution
         rpart(jfiu+jj,ii) = -1.*0.*S_qs*rpart(jvol,ii)*
     >                       rpart(jDuDt+jj,ii)

      elseif (part_force(3) .eq. 0) then
         S_qs = 0.

         rpart(jcmiu,ii) = S_qs

         rpart(jfiu+jj,ii) = 0.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_qs(i,j)
c     calculate quasi steady force with drag corrections
      include 'SIZE'
      include 'TOTAL'
      include 'PERFECTGAS'
      include 'CMTDATA'
      include 'CMTPART'

      real cd,S_qs

c     cd_std = 24/re_p already taken into account below
      S_qs = rpart(jvol,i)*rpart(jrhop,i)/rpart(jtaup,i)
      ! note rpart(jvol,i) already has super particle loading!

      if (part_force(1).eq.1) then
         rrep = rpart(jre,i)
         rrma = rpart(ja,i)
         rmacr= 0.6
         rcd_std = 1.+0.15*rrep**(0.687) + 
     >               0.42*(1.+42500./rrep**(1.16))**(-1)
         rcd_mcr = 1.+0.15*rrep**(0.684) + 
     >               0.513*(1. + 483./rrep**(0.669))**(-1)
         rcd1 = rcd_std + (rcd_mcr - rcd_std)*rrma/rmacr
         
         S_qs = S_qs*rcd1

         rphip = rpart(jvol1,i)
         if (rphip .gt. 0.3) rphip = 0.3
         rcd2 = (1. - 2.*rphip)/(1. - rphip)**3

         S_qs = S_qs*rcd2

         rpart(jfqs+j,i) = S_qs*(rpart(ju0+j,i) - rpart(jv0+j,i))
      elseif (part_force(1).eq.2) then
         rrep = rpart(jre,i)
         rcd_std = 1.+0.15*rrep**(0.687) + 
     >               0.42*(1.+42500./rrep**(1.16))**(-1)
         S_qs = S_qs*rcd_std

         rphip = rpart(jvol1,i)
         if (rphip .gt. 0.3) rphip = 0.3
         rcd2 = (1. - 2.*rphip)/(1. - rphip)**3

         S_qs = S_qs*rcd2

         rpart(jfqs+j,i) = S_qs*(rpart(ju0+j,i) - rpart(jv0+j,i))
      elseif (part_force(1).eq.3) then
         vel_diff = sqrt((rpart(ju0  ,i)-rpart(jv0  ,i))**2+
     >                   (rpart(ju0+1,i)-rpart(jv0+1,i))**2+
     >                   (rpart(ju0+2,i)-rpart(jv0+2,i))**2)
         
         rphip = rpart(jvol1,i)
         rphig = 1. - rpart(jvol1,i)
         rrep = rphig*rpart(jre,i)

         ! ergun
         rbeta1 = 150.*rphip*rphip*mu_0/rphig/rpart(jdp,i)**2 +
     >            1.75*rphip*rpart(jrho,i)*vel_diff/rpart(jdp,i)

         ! wen-yu
         if (rrep .lt. 1000) then
            rcd = 24./rrep*(1. + 0.15*rrep**(0.687))
         else
            rcd = 0.44
         endif

         rbeta2 = 0.75*rcd*rphig**(-2.65)*
     >             rphip*rphig*rpart(jrho,i)*vel_diff/rpart(jdp,i)

         ! stiching
         rs = 0.2
         rpp = atan(150 * 1.75*(rphip - rs))/pi + 0.5
         rbeta = rpp*rbeta1 + (1.-rpp)*rbeta2


         rpart(jfqs+j,i) = rpart(jvol,i)*rbeta/rphip    
     >              *(rpart(ju0+j,i) - rpart(jv0+j,i))


      elseif (part_force(1).eq.0) then
         S_qs = 0.

         rpart(jfqs+j,i) = S_qs
      endif

      rpart(jcd+j,i)  = 0.

      return
      end
c-----------------------------------------------------------------------
      subroutine usr_particles_f_col(i)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real mcfac, rdum, rdum3(3)

      real rrp1,rrp2,rvol1,rvol2,rrho1,rrho2,rx1(3),rx2(3),rv1(3),rv2(3)
      common /tcol_b/ rrp1,rvol1,rrho1,rx1,rv1

      ptdum(11) = dnekclock()

      mcfac  = 2.*sqrt(ksp)*log(e_rest)/sqrt(log(e_rest)**2+pi**2)

      if (two_way .gt. 2) then

      rrp1   = rpart(jrpe ,i)
      rvol1  = rpart(jvol ,i)
      rrho1  = rpart(jrhop,i)
      rx1(1) = rpart(jx   ,i)
      rx1(2) = rpart(jx+1 ,i)
      rx1(3) = rpart(jx+2 ,i)
      rv1(1) = rpart(jv0  ,i)
      rv1(2) = rpart(jv0+1,i)
      rv1(3) = rpart(jv0+2,i)

      icx1 = ipart(jicx,i)
      icy1 = ipart(jicy,i)
      icz1 = ipart(jicz,i)

      icxm = icx1 -1
      icxp = icx1 +1
      icym = icy1 -1
      icyp = icy1 +1
      iczm = icz1 -1
      iczp = icz1 +1

c     let every particle search for itself
c        particles in local elements
         do j = i+1,n
            if (ipart(je0,i) .eq. ipart(je0,j)) then
c           if (i .ne. j) then

               ! finally excude on basis of sub element mesh
               icx2 = ipart(jicx,j)
               icy2 = ipart(jicy,j)
               icz2 = ipart(jicz,j)
               if ((icx2.ge.icxm).and.(icx2.le.icxp)) then
               if ((icy2.ge.icym).and.(icy2.le.icyp)) then
               if ((icz2.ge.iczm).and.(icz2.le.iczp)) then

               rrp2   = rpart(jrpe ,j)
               rvol2  = rpart(jvol ,j)
               rrho2  = rpart(jrhop,j)
               rx2(1) = rpart(jx   ,j)
               rx2(2) = rpart(jx+1 ,j)
               rx2(3) = rpart(jx+2 ,j)
               rv2(1) = rpart(jv0  ,j)
               rv2(2) = rpart(jv0+1,j)
               rv2(3) = rpart(jv0+2,j)

               idum = 1
               call compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
     >                              rpart(jfcol,i),rpart(jfcol,j),idum)

               endif
               endif
               endif
               
c           endif
            endif
         enddo

c        search list of ghost particles
         do j = 1,nfptsgp
            if (ipart(je0,i) .eq. iptsgp(jgpes,j)) then
            ! exclude if not meant for collisions
            if (iptsgp(jgpiic,j) .eq. 2) goto 1235
               rrp2   = rptsgp(jgprpe ,j)
               rvol2  = rptsgp(jgpvol ,j)
               rrho2  = rho_p            ! assume same density. Need2fix
               rx2(1) = rptsgp(jgpx   ,j)
               rx2(2) = rptsgp(jgpx+1 ,j)
               rx2(3) = rptsgp(jgpx+2 ,j)
               rv2(1) = rptsgp(jgpv0  ,j)
               rv2(2) = rptsgp(jgpv0+1,j)
               rv2(3) = rptsgp(jgpv0+2,j)

               idum = 0
               call compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
     >                              rpart(jfcol,i),rdum3,idum)

        endif
         enddo
 1235 continue

c        collision with 6 walls, but only when specified by .inp file 
         do j = 1,6
            if (bc_part(j) .eq. -1) then
               nj1 = mod(j,2)
               if (nj1.ne.0) nj1 = 1
               if (nj1.eq.0) nj1 = 2
               nj2 = int((j-1)/2) + 1

               rrp2   = 0.
               rvol2  = 1.
               rrho2  = 1E8
               rx2(1) = rpart(jx  ,i)
               rx2(2) = rpart(jx+1,i)
               rx2(3) = rpart(jx+2,i)
               rv2(1) = 0.
               rv2(2) = 0.
               rv2(3) = 0.
               rx2(nj2) = xdrange(nj1,nj2)

               idum = 0
               call compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
     >                              rpart(jfcol,i),rdum3,idum)
            endif
         enddo

!        collision with cylinders? put in rxbo(j,1) ...
c        do j=1,2
c           ! at this particles_location
c           rtheta = atan2(rpart(jx+1,i),rpart(jx+0,i))
c           rrad   = rxbo(j,1) ! assumes radius stored here

c           rrp2   = 0.
c           rvol2  = 1.
c           rrho2  = 1E8
c           rx2(1) = rrad*cos(rtheta)
c           rx2(2) = rrad*sin(rtheta)
c           rx2(3) = rpart(jx+2,i)
c           rv2(1) = 0.
c           rv2(2) = 0.
c           rv2(3) = 0.

c           idum = 0
c           call compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,
c    >                              rpart(jfcol,i),rdum3,idum)
c        enddo


      endif

      pttime(11) = pttime(11) + dnekclock() - ptdum(11)

      return
      end
c----------------------------------------------------------------------
      subroutine compute_collide(mcfac,rrp2,rvol2,rrho2,rx2,rv2,fcf1,
     >                                                      fcf2,iflg)
c
c     extra body forces
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real fcf1(3),fcf2(3),er,eta,ere,mcfac
      integer iflg

      real rrp1,rrp2,rvol1,rvol2,rrho1,rrho2,rx1(3),rx2(3),rv1(3),rv2(3)
      common /tcol_b/ rrp1,rvol1,rrho1,rx1,rv1

      rthresh  = rrp1 + rrp2
      rthresh2 = rthresh**2

      rxdiff  = rx2(1) - rx1(1)
      rsum2 = rxdiff**2
      if (rsum2 .gt. rthresh2) goto 1511

      rydiff = rx2(2) - rx1(2)
      rydiff2 = rydiff**2
      rsum2 = rsum2 + rydiff2
      if (rsum2 .gt. rthresh2) goto 1511

      rzdiff = rx2(3) - rx1(3)
      rzdiff2 = rzdiff**2
      rsum2 = rsum2 + rzdiff2
      if (rsum2 .gt. rthresh2) goto 1511

      rdiff = sqrt(rsum2)
      rm1   = rrho1*rvol1
      rm2   = rrho2*rvol2

c     rm12 = rm1*rm2/(rm1 + rm2)
c     eta  = 2.*sqrt(ksp*rm12)*log(e_rest)/sqrt(log(e_rest)**2+pi**2)
c     rtop  = rm1*rm2
c     rbot  = rm1 + rm2
c     rmult = rtop/rbot
c     rmult = rmult**(0.5)
      rmult = 1./sqrt(1./rm1 + 1./rm2)
      eta  = mcfac*rmult

      ! first, handle normal collision part
      rbot     = 1./rdiff
      rn_12x   = rxdiff*rbot
      rn_12y   = rydiff*rbot
      rn_12z   = rzdiff*rbot

      rdelta12 = rthresh - rdiff

      rv12_mag = (rv2(1) - rv1(1))*rn_12x +
     >           (rv2(2) - rv1(2))*rn_12y +
     >           (rv2(3) - rv1(3))*rn_12z

      rv12_mage = rv12_mag*eta

      rksp_max = ksp*rdelta12

      rnmag = -rksp_max - rv12_mage

      rfn1 = rnmag*rn_12x
      rfn2 = rnmag*rn_12y
      rfn3 = rnmag*rn_12z

      fcf1(1) = fcf1(1) + rfn1
      fcf1(2) = fcf1(2) + rfn2
      fcf1(3) = fcf1(3) + rfn3

      if (iflg .eq. 1) then
          fcf2(1) = fcf2(1) - rfn1
          fcf2(2) = fcf2(2) - rfn2
          fcf2(3) = fcf2(3) - rfn3
      endif

 1511 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine send_ghost_particles
c
c     send only ghost particles
c
c     bc_part = -1,1  => non-periodic search
c     bc_part = 0  => periodic search
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      logical partl         ! dummy used in c_t_t()

c     send ghost particles
      call fgslib_crystal_tuple_transfer(i_cr_hndl,nfptsgp,llpart_gp
     $           , iptsgp,nigp,partl,0,rptsgp,nrgp,jgpps) ! jgpps is overwri

c     sort ghost particles by jgpiic for quick discard in collision
c     algorithm. Note that jgpiic loc in iptsgp has values of 0 (used
c     in collisions) or 1 (not used in collisions, but for projection)
      if (two_way .gt. 2) then
      call fgslib_crystal_tuple_sort    (i_cr_hndl,nfptsgp
     $              , iptsgp,nigp,partl,0,rptsgp,nrgp,jgpiic,1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine create_extra_particles
c
c     create ghost and wall particles
c
c     bc_part = -1,1  => non-periodic search
c     bc_part = 0  => periodic search
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      character*132 deathmessage

c     create ghost particles
      call create_ghost_particles_rect
c     if (nrect_assume .gt. 0) call create_wall_particles_image2

      nmax   = iglmax(n,1)
      ngpmax = iglmax(nfptsgp,1)

      if (nmax .gt. llpart) then
         if (nid.eq.0) write(6,1) nmax, llpart, nid
         call exitt
      elseif (ngpmax .gt. llpart_gp) then
         if (nid.eq.0) write(6,2) ngpmax, llpart_gp, nid
         call exitt
      endif
    1 format('Max number of real particles:',
     >   i9,'. Not moving because llpart =',i9, ' on nid = ', i9)
    2 format('Max number of ghost particles:',
     >   i9,'. Not moving because llpart_gp =',i9, ' on nid = ', i9)

      return
      end
c----------------------------------------------------------------------
      subroutine create_ghost_particles_rect
c
c     this routine will create ghost particles by checking if particle
c     is within d2chk of element faces
c
c     ghost particle x,y,z list will be in rptsgp(jgpx,j),rptsgp(jgpy,j),
c     rptsgp(jgpz,j), while processor and local element id are in
c     iptsgp(jgppt,j) and iptsgp(jgpes,j)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      nfptsgp = 0
      do i = 1,n
         ie = ipart(je0,i)+1

         ! projection
         iip1 = 0
         iip2 = 0
         iip3 = 0
         rxdum1 = abs(rpart(jx,i) - xerange(1,1,ie))
         rxdum2 = abs(rpart(jx,i) - xerange(2,1,ie))
         rxdum = min(rxdum1,rxdum2)
         if (rxdum1 .lt. rxdum2) iip1 = -1
         if (rxdum2 .lt. rxdum1) iip1 =  1
         rxdum1 = abs(rpart(jy,i) - xerange(1,2,ie))
         rxdum2 = abs(rpart(jy,i) - xerange(2,2,ie))
         rydum = min(rxdum1,rxdum2)
         if (rxdum1 .lt. rxdum2) iip2 = -1
         if (rxdum2 .lt. rxdum1) iip2 =  1
         rxdum1 = abs(rpart(jz,i) - xerange(1,3,ie))
         rxdum2 = abs(rpart(jz,i) - xerange(2,3,ie))
         rzdum = min(rxdum1,rxdum2)
         if (rxdum1 .lt. rxdum2) iip3 = -1
         if (rxdum2 .lt. rxdum1) iip3 =  1

c        if (abs(rpart(jx,i) - xerange(1,1,ie)).lt.d2chk(1)) iip1 = -1
c        if (abs(rpart(jx,i) - xerange(2,1,ie)).lt.d2chk(1)) iip1 = 1
c        if (abs(rpart(jy,i) - xerange(1,2,ie)).lt.d2chk(1)) iip2 = -1
c        if (abs(rpart(jy,i) - xerange(2,2,ie)).lt.d2chk(1)) iip2 = 1
c        if (abs(rpart(jz,i) - xerange(1,3,ie)).lt.d2chk(1)) iip3 = -1
c        if (abs(rpart(jz,i) - xerange(2,3,ie)).lt.d2chk(1)) iip3 = 1

         ! collisions
         iic1 = 0
         iic2 = 0
         iic3 = 0
c        if (abs(rpart(jx,i) - xerange(1,1,ie)).lt.d2chk(3)) iic1 = -1
c        if (abs(rpart(jx,i) - xerange(2,1,ie)).lt.d2chk(3)) iic1 = 1
c        if (abs(rpart(jy,i) - xerange(1,2,ie)).lt.d2chk(3)) iic2 = -1
c        if (abs(rpart(jy,i) - xerange(2,2,ie)).lt.d2chk(3)) iic2 = 1
c        if (abs(rpart(jz,i) - xerange(1,3,ie)).lt.d2chk(3)) iic3 = -1
c        if (abs(rpart(jz,i) - xerange(2,3,ie)).lt.d2chk(3)) iic3 = 1

         if (rxdum .lt. d2chk(3)) iic1 = iip1
         if (rydum .lt. d2chk(3)) iic2 = iip2
         if (rzdum .lt. d2chk(3)) iic3 = iip3

         do j=1,3*nfacegp-2,3   ! faces
            ii = el_face_num(j) 
            jj = el_face_num(j+1) 
            kk = el_face_num(j+2) 

            iip = 0
            iic = 2
            if (ii .ne. 0) then
               if (iip1 .eq. ii) iip = 1
               if (iic1 .eq. ii) iic = 1
            endif
            if (jj .ne. 0) then
               if (iip2 .eq. jj) iip = 1
               if (iic2 .eq. jj) iic = 1
            endif
            if (kk .ne. 0) then
               if (iip3 .eq. kk) iip = 1
               if (iic3 .eq. kk) iic = 1
            endif

            if (nrect_assume .eq. 2) iip = 1
            if (iip .eq. 1) then
            call gp_create(ii,jj,kk,iic,i,
     >       nfacegp,el_face_num,el_face_proc_map,el_face_el_map)
            endif
         enddo

         do j=1,3*nedgegp-2,3   ! edges
            ii = el_edge_num(j) 
            jj = el_edge_num(j+1) 
            kk = el_edge_num(j+2) 

            iip = 0
            iic = 2
            if (ii .ne. 0) then
               if (jj .ne. 0) then
               ! ii jj
                  if (iip1 .eq. ii) then
                  if (iip2 .eq. jj) then
                     iip = 1
                  endif
                  endif
                  if (iic1 .eq. ii) then
                  if (iic2 .eq. jj) then
                     iic = 1
                  endif
                  endif

               elseif (kk .ne. 0) then
               ! ii kk
                  if (iip1 .eq. ii) then
                  if (iip3 .eq. kk) then
                     iip = 1
                  endif
                  endif
                  if (iic1 .eq. ii) then
                  if (iic3 .eq. kk) then
                     iic = 1
                  endif
                  endif

               endif

            elseif (jj .ne. 0) then
               if (kk .ne. 0) then
               ! jj kk
                  if (iip2 .eq. jj) then
                  if (iip3 .eq. kk) then
                     iip = 1
                  endif
                  endif
                  if (iic2 .eq. jj) then
                  if (iic3 .eq. kk) then
                     iic = 1
                  endif
                  endif

               endif

            endif

            if (nrect_assume .eq. 2) iip = 1

            if (iip .eq. 1) then
            call gp_create(ii,jj,kk,iic,i,
     >       nedgegp,el_edge_num,el_edge_proc_map,el_edge_el_map)
            endif
         enddo

         do j=1,3*ncornergp-2,3   ! corners
            ii = el_corner_num(j) 
            jj = el_corner_num(j+1) 
            kk = el_corner_num(j+2) 

            iip = 0
            iic = 2
            if (iic1 .eq. ii) then
            if (iic2 .eq. jj) then
            if (iic3 .eq. kk) then
               iic = 1
            endif
            endif
            endif
            if (iip1 .eq. ii) then
            if (iip2 .eq. jj) then
            if (iip3 .eq. kk) then
               iip = 1
            endif
            endif
            endif

            if (nrect_assume .eq. 2) iip = 1
            if (iip .eq. 1) then
            call gp_create(ii,jj,kk,iic,i,
     >     ncornergp,el_corner_num,el_corner_proc_map,el_corner_el_map)
            endif
         enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gp_create(ii,jj,kk,iic,i,
     >             nnl,el_tmp_num,el_tmp_proc_map,el_tmp_el_map)
c
c     this routine will create a ghost particle and append its position
c     to rptsgp and its processor and element to iptsgp. nfptsgp will then
c     be incremented. Note that ghost particles will not be created if 
c     they are to be created on the same processor. In the near future, 
c     this might not be true if periodic conditions are needed.
c
c     el_tmp_num holds vector coordinates of tmp=face,edge, or corners
c     el_tmp_proc_map holds MPI rank of neighbor elements in el_tmp_num
c                     order
c     el_tmp_el_map holds local element number of neighbor elements
c
c     ii,jj,kk are vectors that tell what element a ghost particle
c     should be sent to
c
c     i is which particle is creating the ghost particle from rpart,etc
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      integer el_tmp_proc_map(lelt,12)  ,el_tmp_el_map(lelt,12),
     >        el_tmp_num(36)

      real rdumpos(3)

      xdlen = xdrange(2,1) - xdrange(1,1)
      ydlen = xdrange(2,2) - xdrange(1,2)
      zdlen = xdrange(2,3) - xdrange(1,3)

            ie = ipart(je0,i)+1

      xedlen = xerange(2,1,ie) - xerange(1,1,ie)
      yedlen = xerange(2,2,ie) - xerange(1,2,ie)
      zedlen = xerange(2,3,ie) - xerange(1,3,ie)


      ic = 0
      do j=1,3*nnl-2,3
         ic = ic + 1
         if (el_tmp_num(j)  .eq.ii) then
         if (el_tmp_num(j+1).eq.jj) then
         if (el_tmp_num(j+2).eq.kk) then

c           if (nid .eq. el_tmp_proc_map(ie,ic)) then 
c           if (ie .eq.  el_tmp_el_map(ie,ic) + 1) then 
c              goto 1511
c           endif
c           endif


            nfptsgp = nfptsgp + 1
            iitmp1 = 0
            iitmp2 = 0
            iitmp3 = 0

            ! note that altering locs is for bc in periodic ..
            xloc = rpart(jx,i)
            if (xloc+xedlen*ii .gt. xdrange(2,1)) then
                 xloc = rpart(jx,i) - xdlen
                 iitmp1 = 1
                 goto 123
            endif
            if (xloc+xedlen*ii .lt. xdrange(1,1))then
                 xloc = rpart(jx,i) + xdlen
                 iitmp1 = 1
                 goto 123
            endif
  123 continue
            yloc = rpart(jy,i)
            if (yloc+yedlen*jj .gt. xdrange(2,2))then
                 yloc = rpart(jy,i) - ydlen
                 iitmp2 = 1
                 goto 124
            endif
            if (yloc+yedlen*jj .lt. xdrange(1,2))then
                 yloc = rpart(jy,i) + ydlen
                 iitmp2 = 1
                 goto 124
            endif
  124 continue
            zloc = rpart(jz,i)
            if (zloc+zedlen*kk .gt. xdrange(2,3))then
                 zloc = rpart(jz,i) - zdlen
                 iitmp3 = 1
                 goto 125
            endif
            if (zloc+zedlen*kk .lt. xdrange(1,3))then
                 zloc = rpart(jz,i) + zdlen
                 iitmp3 = 1
                 goto 125
            endif
  125 continue

            rptsgp(jgpx,nfptsgp)    = xloc           ! x loc
            rptsgp(jgpy,nfptsgp)    = yloc           ! y log
            rptsgp(jgpz,nfptsgp)    = zloc           ! z log
            rptsgp(jgpfh,nfptsgp)   = rpart(jf0,i)   ! hyd. force x
            rptsgp(jgpfh+1,nfptsgp) = rpart(jf0+1,i) ! hyd. force y
            rptsgp(jgpfh+2,nfptsgp) = rpart(jf0+2,i) ! hyd. force z
            rptsgp(jgpvol,nfptsgp)  = rpart(jvol,i)  ! particle volum
            rptsgp(jgprpe,nfptsgp)  = rpart(jrpe,i)  ! particle rp eff
            rptsgp(jgpspl,nfptsgp)  = rpart(jspl,i)  ! spl
            rptsgp(jgpg0,nfptsgp)   = rpart(jg0,i)   ! work done by forc
            rptsgp(jgpq0,nfptsgp)   = rpart(jq0,i)   ! heating from part 
            rptsgp(jgpv0,nfptsgp)   = rpart(jv0,i)   ! particle velocity
            rptsgp(jgpv0+1,nfptsgp) = rpart(jv0+1,i) ! particle velocity
            rptsgp(jgpv0+2,nfptsgp) = rpart(jv0+2,i) ! particle velocity

            iptsgp(jgpiic,nfptsgp)  = iic            ! use in collisions

               ipdum  = el_tmp_proc_map(ie,ic)
               iedum  = el_tmp_el_map(ie,ic)

! DZ
               if (ipdum .lt. 0 .or. iedum .lt.0) then
                  nfptsgp=nfptsgp-1
                  goto 1511
               endif


            iptsgp(jgpps,nfptsgp)   = ipdum  ! overwritten mpi
            iptsgp(jgppt,nfptsgp)   = ipdum  ! dest. mpi rank
            iptsgp(jgpes,nfptsgp)   = iedum    ! dest. elment

c           check if extra particles have been created on the same mpi
c           rank and also take care of boundary particles
            ibctype = abs(bc_part(1))+abs(bc_part(3))+abs(bc_part(5))

c           if (nid.eq. ipdum) then
c           if (iic .ne. 0) then
c              iptsgp(jgpiic,nfptsgp) = 0
c              goto 1511
c           endif
c           endif

! DZ
c           take care of periodic stuff first
            if (nid.eq.iptsgp(jgppt,nfptsgp)) then ! dont create gp on own rank 
                                                   ! unless moved and periodic
            if (ibctype .eq. 0) then            ! all three sides periodic
               if (iitmp1+iitmp2+iitmp3 .eq.0) then
c                 nfptsgp=nfptsgp-1
                  iptsgp(jgpiic,nfptsgp) = 0
                  goto 1511
               endif
            elseif (ibctype .eq. 1) then        ! only two sides periodic
               if (abs(bc_part(1)) .eq. 1) then
                  if (iitmp2+iitmp3 .eq. 0) then
c                    nfptsgp=nfptsgp-1
                     iptsgp(jgpiic,nfptsgp) = 0
                     goto 1511
                  endif
               elseif (abs(bc_part(3)) .eq. 1) then
                  if (iitmp1+iitmp3 .eq. 0) then
c                    nfptsgp=nfptsgp-1
                     iptsgp(jgpiic,nfptsgp) = 0
                     goto 1511
                  endif
               elseif (abs(bc_part(5)) .eq. 1) then
                  if (iitmp1+iitmp2 .eq. 0) then
c                    nfptsgp=nfptsgp-1
                     iptsgp(jgpiic,nfptsgp) = 0
                     goto 1511
                  endif
               endif
            elseif (ibctype .eq. 2) then        ! only one side periodic
               if (bc_part(1) .eq. 0) then
                  if (iitmp1 .eq. 0) then
c                    nfptsgp=nfptsgp-1
                     iptsgp(jgpiic,nfptsgp) = 0
                     goto 1511
                  endif
               elseif (bc_part(3) .eq. 0) then
                  if (iitmp2 .eq. 0) then
c                    nfptsgp=nfptsgp-1
                     iptsgp(jgpiic,nfptsgp) = 0
                     goto 1511
                  endif
               elseif (bc_part(5) .eq. 0) then
                  if (iitmp3 .eq. 0) then
c                    nfptsgp=nfptsgp-1
                     iptsgp(jgpiic,nfptsgp) = 0
                     goto 1511
                  endif
               endif
            elseif (ibctype .eq. 3) then        ! no sides periodic 
c              nfptsgp=nfptsgp-1
               iptsgp(jgpiic,nfptsgp) = 0
               goto 1511
            endif
            endif ! end if(nid.eq. ...)


c           take care of non-periodic stuff second
c           if (ibctype .gt. 0) then
c              if (ibctype .eq. 3) then         ! no sides periodic
c                 if (iitmp1+iitmp2+iitmp3 .gt.0) then
c                    nfptsgp=nfptsgp-1
c                    goto 1511
c                 endif
c              elseif (ibctype .eq.1) then      ! two sides periodic
c                 if (abs(bc_part(1)) .eq. 1) then
c                    if (iitmp1 .gt. 0) then
c                       nfptsgp=nfptsgp-1
c                       goto 1511
c                    endif
c                 elseif (abs(bc_part(3)) .eq. 1) then
c                    if (iitmp2 .gt. 0) then
c                       nfptsgp=nfptsgp-1
c                       goto 1511
c                    endif
c                 elseif (abs(bc_part(5)) .eq. 1) then
c                    if (iitmp3 .gt. 0) then
c                       nfptsgp=nfptsgp-1
c                       goto 1511
c                    endif
c                 endif
c              elseif (ibctype .eq.2) then      ! one side periodic
c                 if (bc_part(1) .eq. 0) then
c                    if (iitmp2+iitmp3.gt.0) then
c                       nfptsgp=nfptsgp-1
c                       goto 1511
c                    endif
c                 elseif (bc_part(3) .eq. 0) then
c                    if (iitmp1+iitmp3.gt.0) then
c                       nfptsgp=nfptsgp-1
c                       goto 1511
c                    endif
c                 elseif (bc_part(5) .eq. 0) then
c                    if (iitmp1+iitmp2.gt.0) then
c                       nfptsgp=nfptsgp-1
c                       goto 1511
c                    endif
c                 endif
c              endif
c           endif

            goto 1511
         endif
         endif
         endif
      enddo
 1511 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_neighbor_el_proc
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'
c
c     This routine is called once at the beginning of the particle
c     simulation. At the end of this routine, the common blocks
c     /neighbor_proc/ & /neighbor_el_number/ are set. The idea behind
c     this routine is to know what processor owns neighboring spectral
c     elements and what local element number the neighboring element is.
c
c     el_*_proc_map holds: *(face,edge,corner) neighboring element 
c                           MPI rank number
c     el_*_el_map holds:   *(face,edge,corner) neighboring element
c                           local numbers
c
c     The ordering of faces, edges, and corners are given in el_*_num
c
c     el_*_proc_map(i,j) and el_*_el_map(i,j) are ordered by elements 
c     1 <= i <= nelt, and 1 <= j <= 26, where j=1,nfacegp are element
c     faces, j=nfacegp+1,nfacegp+nedgegp are element edges, and 
c     j = nfacegp+nedgegp+1,nfacegp+nedgegp+ncornergp are corners

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      common /gpfix/ ilgp_f(lelt,6),ilgp_e(lelt,12),ilgp_c(lelt,8)

      real  rimp(7,lelt*26)
      integer iimp(4,lelt*26)

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/
     >                 0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                 0,-1,1,  1,0,1,  0,1,1,  -1,0,1,
     >                 -1,-1,0, 1,-1,0, 1,1,0,  -1,1,0
     >              /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1
     >              /)
      nfacegp   = 6  ! number of faces
      nedgegp   = 12 ! number of edges
      ncornergp = 8  ! number of corners
      idum      = 0  ! dummy arguement

      rtmult = 1.5

      icount = 0 
      do i=1,nelt
         do j=1,3*nfacegp-2,3   ! faces
            icount = icount + 1

            xlen = xerange(2,1,i)-xerange(1,1,i)
            ylen = xerange(2,2,i)-xerange(1,2,i)
            zlen = xerange(2,3,i)-xerange(1,3,i)

            xmid = (xerange(2,1,i)+xerange(1,1,i))/2.0
            ymid = (xerange(2,2,i)+xerange(1,2,i))/2.0
            zmid = (xerange(2,3,i)+xerange(1,3,i))/2.0

            isignxx = el_face_num(j) 
            isignyy = el_face_num(j+1) 
            isignzz = el_face_num(j+2) 

            xloc = xmid + (xlen*rtmult)*isignxx/2.0
            yloc = ymid + (ylen*rtmult)*isignyy/2.0
            zloc = zmid + (zlen*rtmult)*isignzz/2.0
            iimp(4,icount) = 0
            idum = 1
            call bounds_p_check(xloc,xdrange(1,1),xdrange(2,1),idum)
            if (idum .ne. 0) iimp(4,icount) = -1
            idum = 2
            call bounds_p_check(yloc,xdrange(1,2),xdrange(2,2),idum)
            if (idum .ne. 0) iimp(4,icount) = -1
            idum = 3
            call bounds_p_check(zloc,xdrange(1,3),xdrange(2,3),idum)
            if (idum .ne. 0) iimp(4,icount) = -1

            rimp(1,icount) = xloc
            rimp(2,icount) = yloc
            rimp(3,icount) = zloc

            iimp(1,icount) = nid
            iimp(2,icount) = 0
            iimp(3,icount) = i-1
         enddo
         do j=1,3*nedgegp-2,3    ! edges
            icount = icount + 1

            xlen = xerange(2,1,i)-xerange(1,1,i)
            ylen = xerange(2,2,i)-xerange(1,2,i)
            zlen = xerange(2,3,i)-xerange(1,3,i)

            xmid = (xerange(2,1,i)+xerange(1,1,i))/2.0
            ymid = (xerange(2,2,i)+xerange(1,2,i))/2.0
            zmid = (xerange(2,3,i)+xerange(1,3,i))/2.0

            isignxx = el_edge_num(j) 
            isignyy = el_edge_num(j+1) 
            isignzz = el_edge_num(j+2) 

            xloc = xmid + (xlen*rtmult)*isignxx/2.0
            yloc = ymid + (ylen*rtmult)*isignyy/2.0
            zloc = zmid + (zlen*rtmult)*isignzz/2.0
            iimp(4,icount) = 0
            idum = 1
            call bounds_p_check(xloc,xdrange(1,1),xdrange(2,1),idum)
            if (idum .ne. 0) iimp(4,icount) = -1
            idum = 2
            call bounds_p_check(yloc,xdrange(1,2),xdrange(2,2),idum)
            if (idum .ne. 0) iimp(4,icount) = -1
            idum = 3
            call bounds_p_check(zloc,xdrange(1,3),xdrange(2,3),idum)
            if (idum .ne. 0) iimp(4,icount) = -1

            rimp(1,icount) = xloc
            rimp(2,icount) = yloc
            rimp(3,icount) = zloc

            iimp(1,icount) = nid
            iimp(2,icount) = 0
            iimp(3,icount) = i-1
         enddo
         do j=1,3*ncornergp-2,3   ! corners
            icount = icount + 1

            xlen = xerange(2,1,i)-xerange(1,1,i)
            ylen = xerange(2,2,i)-xerange(1,2,i)
            zlen = xerange(2,3,i)-xerange(1,3,i)

            xmid = (xerange(2,1,i)+xerange(1,1,i))/2.0
            ymid = (xerange(2,2,i)+xerange(1,2,i))/2.0
            zmid = (xerange(2,3,i)+xerange(1,3,i))/2.0

            isignxx = el_corner_num(j) 
            isignyy = el_corner_num(j+1) 
            isignzz = el_corner_num(j+2) 

            xloc = xmid + (xlen*rtmult)*isignxx/2.0
            yloc = ymid + (ylen*rtmult)*isignyy/2.0
            zloc = zmid + (zlen*rtmult)*isignzz/2.0
            iimp(4,icount) = 0
            idum = 1
            call bounds_p_check(xloc,xdrange(1,1),xdrange(2,1),idum)
            if (idum .ne. 0) iimp(4,icount) = -1
            idum = 2
            call bounds_p_check(yloc,xdrange(1,2),xdrange(2,2),idum)
            if (idum .ne. 0) iimp(4,icount) = -1
            idum = 3
            call bounds_p_check(zloc,xdrange(1,3),xdrange(2,3),idum)
            if (idum .ne. 0) iimp(4,icount) = -1

            rimp(1,icount) = xloc
            rimp(2,icount) = yloc
            rimp(3,icount) = zloc

            iimp(1,icount) = nid
            iimp(2,icount) = 0
            iimp(3,icount) = i-1
         enddo
      enddo

c     get processor and local element number of neighboring elemetns
      call fgslib_findpts(i_fp_hndl !  stride     !   call fgslib_findpts( ihndl,
     $           , iimp(2,1),4        !   $             rcode,1,
     $           , iimp(1,1),4        !   &             proc,1,
     $           , iimp(3,1),4        !   &             elid,1,
     $           , rimp(5,1),7        !   &             rst,ndim,
     $           , rimp(4,1),7        !   &             dist,1,
     $           , rimp(1,1),7        !   &             pts(    1),1,
     $           , rimp(2,1),7        !   &             pts(  n+1),1,
     $           , rimp(3,1),7 ,icount)    !   &             pts(2*n+1),1,n)

c     set common block values to be used later
      do i = 1,nelt
         nstride = (i-1)*(nfacegp+nedgegp+ncornergp)
         do j = 1,nfacegp
            ijloc = nstride + j
            if (iimp(2,ijloc) .eq. 0) then
            el_face_proc_map(i,j) = iimp(1,ijloc)
            el_face_el_map(i,j) = iimp(3,ijloc)
            else
            el_face_proc_map(i,j) = -1
            el_face_el_map(i,j) = -1
            endif
            ilgp_f(i,j) = iimp(4,ijloc)
         enddo
         nstride = nstride + nfacegp 
         do j = 1,nedgegp
            ijloc = nstride + j
            if (iimp(2,ijloc) .eq. 0) then
            el_edge_proc_map(i,j) = iimp(1,ijloc)
            el_edge_el_map(i,j) = iimp(3,ijloc)
            else
            el_edge_proc_map(i,j) = -1
            el_edge_el_map(i,j) = -1
            endif
            ilgp_e(i,j) = iimp(4,ijloc)
         enddo
         nstride = nstride + nedgegp 
         do j = 1,ncornergp
            ijloc = nstride + j
            if (iimp(2,ijloc) .eq. 0) then
            el_corner_proc_map(i,j) = iimp(1,ijloc)
            el_corner_el_map(i,j) = iimp(3,ijloc)
            else
            el_corner_proc_map(i,j) = -1
            el_corner_el_map(i,j) = -1
            endif
            ilgp_c(i,j) = iimp(4,ijloc)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine bounds_p_check(xx,xl,xr,ifmove)
      include 'SIZE'
      include 'CMTPART'
c     
c     check if xx is between domain bounds of left (xl) and right (xr)
c
c     if it is outside of bounds, move periodically to other domain side
c     and set ifmove to 1 so we know if xx has been changed
c

      ! use fact that both sides have to be set as periodic...
      rdum_val = 2.*xr         ! dummy value definitley not in domain
      if (ifmove .eq. 1) then
         if (bc_part(1) .ne. 0) then
            if (xx .gt. xr .or. xx .lt. xl) then
               ifmove = -1
               goto 1511
            endif
         endif
      elseif (ifmove .eq. 2) then
         if (bc_part(3) .ne. 0) then
            if (xx .gt. xr .or. xx .lt. xl) then
               ifmove = -1
               goto 1511
            endif
         endif
      elseif (ifmove .eq. 3) then
         if (bc_part(5) .ne. 0) then
            if (xx .gt. xr .or. xx .lt. xl) then
               ifmove = -1
               goto 1511
            endif
         endif
      endif

      ifmove = 0
      if (xx .gt. xr) then
         xx = xx - (xr -xl)
c        xx = abs(xx - xr) + xl
c        xx = xx - xr
         ifmove = 1
         goto 1511
      endif
      if (xx .lt. xl) then
         xx = xx + (xr -xl)
c        xx = xr - abs(xx - xl) 
c        xx = xx + xl
         ifmove = 1
         goto 1511
      endif

 1511 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine get_bdf_ext_coefs(beta,alpha,times)
      include 'SIZE'
      include 'TOTAL'
      include 'CMTPART'

      real beta(0:3),alpha(0:3),times(0:3)
      real c(0:8)

      integer ilast,ncoef
      save    ilast,ncoef
      data    ilast,ncoef / -9 , 0 /

      do i=3,1,-1
         times(i)=times(i-1)
      enddo
      times(0) = time

      call rzero(beta ,4)
      call rzero(alpha,4)
      if (istep.ne.ilast) then
         ilast = istep
         ncoef = ncoef + 1
         ncoef = min(ncoef,3) ! Maximum 3rd order in time
      endif
      ncoefm1 = ncoef - 1

      call fd_weights_full(times(0),times(1),ncoefm1,0,alpha(1))
      call fd_weights_full(times(0),times(0),ncoef,1,c)
      do j=0,ncoef
         beta(j) = c(ncoef+1+j)
      enddo
      do j=1,ncoef
         beta(j) = -beta(j)  ! Change sign, for convenience
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine update_particle_location
c     check if particles are outside domain
c     > if bc_part = 0 then it is periodic
c     > if bc_part = -1,1 then particles are killed (outflow)
      include 'SIZE'
      include 'CMTDATA'
      include 'CMTPART'

      integer in_part(llpart)

      jx0 = jx

      do i=1,n
         in_part(i) = 0
         do j=0,ndim-1
            if (rpart(jx0+j,i).lt.xdrange(1,j+1))then
c           if (bc_part(1).eq.0) then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx0+j,i) = xdrange(2,j+1) - 
     &                             abs(xdrange(1,j+1) - rpart(jx0+j,i))
                  rpart(jx1+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx1+j,i))
                  rpart(jx2+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx2+j,i))
                  rpart(jx3+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx3+j,i))
                  goto 1512
c              elseif (bc_part(1).eq. 1) then
               elseif (((bc_part(1).ne.0) .and. (j.eq.0)) .or. ! outflow
     >                 ((bc_part(3).ne.0) .and. (j.eq.1)) .or.     
     >                 ((bc_part(5).ne.0) .and. (j.eq.2)) ) then
                  in_part(i) = -1
                  goto 1511
               endif
            endif
            if (rpart(jx0+j,i).gt.xdrange(2,j+1))then
c              if (bc_part(1).eq. 0) then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx0+j,i) = xdrange(1,j+1) +
     &                             abs(rpart(jx0+j,i) - xdrange(2,j+1))
                  rpart(jx1+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx1+j,i) - xdrange(2,j+1))
                  rpart(jx2+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx2+j,i) - xdrange(2,j+1))
                  rpart(jx3+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx3+j,i) - xdrange(2,j+1))
                  goto 1512
c              elseif (bc_part(1).eq. 1) then
               elseif (((bc_part(1).ne.0) .and. (j.eq.0)) .or. ! outflow
     >                 ((bc_part(3).ne.0) .and. (j.eq.1)) .or.     
     >                 ((bc_part(5).ne.0) .and. (j.eq.2)) ) then
                  in_part(i) = -1
                  goto 1511
               endif
            endif
 1512 continue
         enddo
 1511 continue
      enddo

      nbc_sum = abs(bc_part(1)) + abs(bc_part(2)) + 
     >          abs(bc_part(3)) + abs(bc_part(4)) +
     >          abs(bc_part(5)) + abs(bc_part(6)) ! all periodic, don't search
      if (nbc_sum .gt. 0) then
      ic = 0
      do i=1,n
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               call copy(rpart(1,ic),rpart(1,i),nr)
               call icopy(ipart(1,ic),ipart(1,i),ni)
            endif
         endif
      enddo
      n = ic
      endif

      return
      end
c-----------------------------------------------------------------------
c     interpolation routines
c-----------------------------------------------------------------------
      subroutine init_interpolation
      include 'SIZE' 
      include 'INPUT' 
      include 'CMTPART' 
c
c     calculates the barycentric lagrange weights
c

c     get gll points in all directions
      call zwgll(xgll,wxgll,lx1)
      call zwgll(ygll,wygll,ly1)
      call rone(zgll,lz1)
      if(if3d) call zwgll(zgll,wzgll,lz1)
c     set all weights to ones first
      call rone(wxgll,lx1)
      call rone(wygll,ly1)
      call rone(wzgll,lz1)
c
c     copy for reduced interpolation
      nx1r = nx1
      if (red_interp.eq.1) then
         nx1r = (nx1+1)/2
         ic = 0
         do j=1,lx1,2
            ic = ic + 1
            xgll(ic) = xgll(j)
            ygll(ic) = ygll(j)
            zgll(ic) = zgll(j)
         enddo
      endif

c     calc x bary weights
      do j=1,nx1r
         do k=1,nx1r
            if (j .NE. k) then
               wxgll(j) = wxgll(j)/(xgll(j) - xgll(k))
            endif
         enddo
      enddo
c     calc y bary weights
      do j=1,nx1r
         do k=1,nx1r
            if (j .NE. k) then
               wygll(j) = wygll(j)/(ygll(j) - ygll(k))
            endif
         enddo
      enddo
c     calc z bary weights
      do j=1,nx1r
         do k=1,nx1r
            if (j .NE. k) then
               wzgll(j) = wzgll(j)/(zgll(j) - zgll(k))
            endif
         enddo
      enddo

      return 
      end
c-----------------------------------------------------------------------
      subroutine init_baryinterp(x,y,z,nxyz)
c     used for 3d interpolation only
      include 'SIZE'
      include 'CMTPART'
      common /BARRYREP/ rep, bot
      real              rep(lx1,ly1,lz1), bot

      real x, y, z, repy, repz,repx,diff
      real bwgtx(lx1),bwgty(ly1),bwgtz(lz1)

      bot= 0.00
      do k=1,nx1r
         diff = z - zgll(k)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgtz(k) = wzgll(k)/diff
      enddo
      do i=1,nx1r
         diff = x - xgll(i)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgtx(i) = wxgll(i)/diff
      enddo 
      do j=1,nx1r
         diff = y-ygll(j)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgty(j) = wygll(j)/diff
      enddo

      do k=1,nx1r
      do j=1,nx1r
         repdum = bwgty(j)*bwgtz(k)
      do i=1,nx1r
         rep(i,j,k) =  repdum* bwgtx(i)
         bot        =  bot + rep(i,j,k)
      enddo
      enddo
      enddo 

      do k=1,nx1r
      do j=1,nx1r
      do i=1,nx1r
         rep(i,j,k) =  rep(i,j,k)/bot
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine baryinterp(field,pofx,nxyz)
c     used for 3d interpolation only
      include 'SIZE'
      include 'CMTPART'
      common /BARRYREP/ rep, bot
      real              rep(lx1,ly1,lz1), bot
      real field(1),pofx,top

      pofx = 0.00
      if (nx1r.eq.lx1) then ! full interpolation
      do i=1,nxyz
         pofx =  pofx + rep(i,1,1)*field(i)
      enddo
      else                  ! reduced interpolation

      kk = 0 
      do k=1,nx1,2
         kk = kk + 1
         jj = 0
         ijk3 = (k-1)*nx1**2
      do j=1,nx1,2
         jj = jj + 1
         ii = 0
         ijk2 = ijk3+(j-1)*nx1
      do i=1,nx1,2
         ii   = ii + 1
         ijk1 = ijk2 + i
         pofx =  pofx + rep(ii,jj,kk)*field(ijk1)
      enddo
      enddo
      enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine triinterp(xf,yf,zf,field,x,y,z,r,s,t,ie,pval)
c     
c     used for 3d trilinear interpolation
c
      include 'SIZE'
      include 'CMTPART'

      real field(nx1,ny1,nz1),xf(nx1,ny1,nz1),yf(nx1,ny1,nz1),
     >                        zf(nx1,ny1,nz1)
      real x,y,z,pval,c00,c01,c10,c11,c0,c1_0,c1_1,r,s,t

      rdelta = 2./(nx1-1.)
      sdelta = 2./(ny1-1.)
      tdelta = 2./(nz1-1.)

      mxx = floor((1.+r)/rdelta)+1
      myy = floor((1.+s)/sdelta)+1
      mzz = floor((1.+t)/tdelta)+1

      xd = (x - xf(mxx,myy,mzz))/(xf(mxx+1,myy,mzz)-xf(mxx,myy,mzz))
      yd = (y - yf(mxx,myy,mzz))/(yf(mxx,myy+1,mzz)-yf(mxx,myy,mzz))
      zd = (z - zf(mxx,myy,mzz))/(zf(mxx,myy,mzz+1)-zf(mxx,myy,mzz))

      c00=field(mxx,myy,mzz)*(1.-xd)+field(mxx+1,myy,mzz)*xd
      c01=field(mxx,myy,mzz+1)*(1.-xd)+field(mxx+1,myy,mzz+1)*xd
      c10=field(mxx,myy+1,mzz)*(1.-xd)+field(mxx+1,myy+1,mzz)*xd
      c11=field(mxx,myy+1,mzz+1)*(1.-xd)+field(mxx+1,myy+1,mzz+1)*xd

      c1_0 = c00*(1.-yd) + c10*yd
      c1_1 = c01*(1.-yd) + c11*yd

      pval = c1_0*(1.-zd) + c1_1*zd

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_props_part_location
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'
      include 'GEOM'
      common /BARRYREP/ rep, bot
      real              rep(lx1,ly1,lz1), bot

      common /point2gridc/ p2gc
      real   p2gc(lx1,ly1,lz1,lelt,4)

      nxyz = nx1*ny1*nz1
        do i=1,n
           rrdum = 1.0
           if(if3d) rrdum = rpart(jr+2,i)
           ie  =  ipart(je0,i) + 1

           ! Barycentric or reduced barycentric lagrange interpolation
           if (red_interp .le. 2) then 

           call init_baryinterp(rpart(jr,i),rpart(jr+1,i),rrdum,nxyz)

           call baryinterp(vx(1,1,1,ie),rpart(ju0,i),nxyz)   !fluid uvel
           call baryinterp(vy(1,1,1,ie),rpart(ju0+1,i),nxyz) !fluid vvel
           call baryinterp(vz(1,1,1,ie),rpart(ju0+2,i),nxyz) !fluid wvel

           call baryinterp(t(1,1,1,ie,1),rpart(jtempf,i),nxyz)    !fluid temp
           call baryinterp(vtrans(1,1,1,ie,1),rpart(jrho,i),nxyz) !fluid dens
           call baryinterp(rhs_fluidp(1,1,1,ie,1),rpart(jDuDt,i),nxyz)  !dp/dx
           call baryinterp(rhs_fluidp(1,1,1,ie,2),rpart(jDuDt+1,i),nxyz)!dp/dy
           call baryinterp(rhs_fluidp(1,1,1,ie,3),rpart(jDuDt+2,i),nxyz)!dp/dz
           call baryinterp(ptw(1,1,1,ie,4),rpart(jvol1,i),nxyz)   !phi_p

c          call baryinterp(p2gc(1,1,1,ie,4),rpart(jgam,i),nxyz)   !gam

           ! trilinear interpolation between grid points
           else
              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),vx(1,1,1,ie),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(ju0,i))
              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),vy(1,1,1,ie),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(ju0+1,i))
              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),vz(1,1,1,ie),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(ju0+2,i))

              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),t(1,1,1,ie,1),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(jtempf,i))
              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),vtrans(1,1,1,ie,1),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(jrho,i))
              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),rhs_fluidp(1,1,1,ie,1),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(jDuDt,i))
              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),rhs_fluidp(1,1,1,ie,2),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(jDuDt+1,i))
              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),rhs_fluidp(1,1,1,ie,3),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(jDuDt+2,i))
              call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
     >                       zm1(1,1,1,ie),ptw(1,1,1,ie,4),
     >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
     >                       rpart(jr,i),rpart(jr+1,i),rrdum,
     >                       ie,rpart(jvol1,i))

c             call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >                       zm1(1,1,1,ie),p2gc(1,1,1,ie,4),
c    >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
c    >                       rpart(jr,i),rpart(jr+1,i),rrdum,
c    >                       ie,rpart(jgam,i))
           endif

        enddo

      return
      end
c----------------------------------------------------------------------
c     particle input/output/ restart routines
c----------------------------------------------------------------------
      subroutine usr_particles_io
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'
      include 'mpif.h'

      rdumt = dnekclock()

c     always output fields if 2 or 4 way coupled 
      if (two_way .gt. 1) then
         call output_two_way_io
      endif

c     output diagnostics to logfile
      if (npio_method .lt. 0) then
         call output_particle_timers 
         call output_particle_diagnostics
      endif

c     output particle  information
      if     (abs(npio_method) .eq. 1) then
         call output_parallel_lagrangian_parts

      elseif (abs(npio_method) .eq. 2) then
c        call output_parallel_lagrangian_parts
         call output_along_line_avg_x
         call output_along_line_avg_y
         call output_along_line_avg_z

      elseif (abs(npio_method) .eq. 3) then
c        call output_parallel_lagrangian_parts
         call output_along_line_avg_x

      elseif (abs(npio_method) .eq. 4) then
c        call output_parallel_lagrangian_parts
         call output_along_line_avg_y

      elseif (abs(npio_method) .eq. 5) then
c        call output_parallel_lagrangian_parts
         call output_along_line_avg_z

      elseif (abs(npio_method) .eq. 6) then
         call output_along_line_avg_x
         call output_along_line_avg_y
         call output_along_line_avg_z

      elseif (abs(npio_method) .eq. 7) then
         call output_along_line_avg_x

      elseif (abs(npio_method) .eq. 8) then
         call output_along_line_avg_y

      elseif (abs(npio_method) .eq. 9) then
         call output_along_line_avg_z

      endif

      ! Always output restart files
c     output restart information if needed
c     if (ipart_restarto .gt. 0) then
c        if (mod(nistep,ipart_restarto) .eq. 0) then
            call output_parallel_restart_part
c        endif
c     endif

      pttime(1) = pttime(1) - (dnekclock() - rdumt)

      return
      end
c----------------------------------------------------------------------
      subroutine output_parallel_restart_part
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'
      include 'mpif.h'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      integer icalld
      save    icalld
      data    icalld  /0/

      character*16 locstring, datastring
      integer*8    disp, stride_len 
      integer      status_mpi(MPI_STATUS_SIZE)
      integer      prevs(0:np-1),npt_total,e,oldfile
      integer      color,particle_io_comm,key

      logical partl         ! This is a dummy placeholder, used in cr()

! ----------------------------------------
! Setup file names to write to mpi
! ----------------------------------------
      icalld = icalld+1
      write(locstring,'(A8,I5.5,A3)') 'rpartxyz', icalld, '.3D' 
      write(datastring,'(A9,I5.5)')   'rpartdata', icalld

      nptot = iglsum(n,1)
      if (nid.eq.0) then
      open(364, file=datastring, action="write")
         write(364,*) nptot
      close(364)
      endif

! ----------------------------------------
! Calculate how many processors to be used
! ----------------------------------------
      npmax_set = int(nptot/llpart) + 1 ! use as few procs as possible
      npmax = np
      if (npmax .gt. npmax_set) npmax = npmax_set

! ----------------------------------------
! Create communicator with that many ranks
! ----------------------------------------
      color = MPI_UNDEFINED
      if (nid.lt.npmax) color = 0
      key=nid
      call MPI_Comm_split(nekcomm,color,key,particle_io_comm,ierr)

! ------------------------------------
! Get a global prefix sum of partilces
! ------------------------------------
      call MPI_Send(n, 1, MPI_INTEGER,0,0, nekcomm, ierr)
      
      ! keep track of how many particles are on previous procs
      if (nid.eq. 0) then
          prevs(0) = n
          do i=0,np-1
             call MPI_Recv(prevs(i),1,MPI_INTEGER,i,
     >                     0,nekcomm,status_mpi,ierr)
          enddo
      endif
      call MPI_BCAST(prevs,np,MPI_INTEGER,0,nekcomm,ierr)

! ------------------------------------------------------
! Each rank counts how many particles in ranks before it
! ------------------------------------------------------
      stride_len = 0
      if (nid .ne. 0) then
      do i=1,nid
         stride_len = stride_len + prevs(i-1)
      enddo
      endif

! -------------------------------------
! Set new rank to send to and then send
! -------------------------------------
      do i = 1,n
         ipart(jps,i) = int((stride_len + i)/llpart)
      enddo
      nl = 0
      call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
     >                  , ipart,ni,partl,nl,rpart,nr,jps)

! -----------------------------------------
! Only use npmax ranks to read in particles
! -----------------------------------------
      if (nid .lt. npmax) then

         do i = 1,n
            rfpts(1,i) = rpart(jx,i)
            rfpts(2,i) = rpart(jy,i)
            rfpts(3,i) = rpart(jz,i)
            rfpts(4,i) = rpart(jx1+0,i)
            rfpts(5,i) = rpart(jx1+1,i)
            rfpts(6,i) = rpart(jx1+2,i)
            rfpts(7,i) = rpart(jx2+0,i)
            rfpts(8,i) = rpart(jx2+1,i)
            rfpts(9,i) = rpart(jx2+2,i)
            rfpts(10,i) = rpart(jx3+0,i)
            rfpts(11,i) = rpart(jx3+1,i)
            rfpts(12,i) = rpart(jx3+2,i)
         
            rfpts(13,i) = rpart(jv0,i)
            rfpts(14,i) = rpart(jv0+1,i)
            rfpts(15,i) = rpart(jv0+2,i)
            rfpts(16,i) = rpart(jv1+0,i)
            rfpts(17,i) = rpart(jv1+1,i)
            rfpts(18,i) = rpart(jv1+2,i)
            rfpts(19,i) = rpart(jv2+0,i)
            rfpts(20,i) = rpart(jv2+1,i)
            rfpts(21,i) = rpart(jv2+2,i)
            rfpts(22,i) = rpart(jv3+0,i)
            rfpts(23,i) = rpart(jv3+1,i)
            rfpts(24,i) = rpart(jv3+2,i)
         
            rfpts(25,i) = rpart(ju0,i)
            rfpts(26,i) = rpart(ju0+1,i)
            rfpts(27,i) = rpart(ju0+2,i)
            rfpts(28,i) = rpart(ju1+0,i)
            rfpts(29,i) = rpart(ju1+1,i)
            rfpts(30,i) = rpart(ju1+2,i)
            rfpts(31,i) = rpart(ju2+0,i)
            rfpts(32,i) = rpart(ju2+1,i)
            rfpts(33,i) = rpart(ju2+2,i)
            rfpts(34,i) = rpart(ju3+0,i)
            rfpts(35,i) = rpart(ju3+1,i)
            rfpts(36,i) = rpart(ju3+2,i)
         
            rfpts(37,i) = rpart(jdp,i)
            rfpts(38,i) = rpart(jspl,i)
            rfpts(39,i) = rpart(jtemp,i)
            rfpts(40,i) = real(ipart(jpid1,i))
            rfpts(41,i) = real(ipart(jpid2,i))
            rfpts(42,i) = real(ipart(jpid3,i))
         enddo

! -----------------------------------------------
! Each rank sends how many particles it will read
! -----------------------------------------------
         call MPI_Send(n, 1, MPI_INTEGER,0,0, particle_io_comm, ierr)
         
         ! keep track of how many particles are on previous procs
         if (nid.eq. 0) then
             prevs(0) = n
             do i=0,npmax-1
                call MPI_Recv(prevs(i),1,MPI_INTEGER,i,
     >                        0,particle_io_comm,status_mpi,ierr)
             enddo
         endif
         call MPI_BCAST(prevs,npmax,MPI_INTEGER,0,particle_io_comm,ierr)

! ------------------------------------------------------
! Each rank counts how many particles in ranks before it
! ------------------------------------------------------
         stride_len = 0
         if (nid .ne. 0) then
         do i=1,nid
            stride_len = stride_len + prevs(i-1)
         enddo
         endif

! -------------------------
! Parallel MPI file read in
! -------------------------
         call MPI_FILE_OPEN(particle_io_comm, locstring,
     >                   MPI_MODE_CREATE + MPI_MODE_WRONLY, 
     >                   MPI_INFO_NULL, oldfile, ierr) 
   
         disp = stride_len*lrf*8  ! lrf properties each with 8 bytes
         call MPI_FILE_SET_VIEW(oldfile, disp, MPI_DOUBLE_PRECISION,
     >                       MPI_DOUBLE_PRECISION, "native", 
     >                       MPI_INFO_NULL, ierr) 
         call MPI_FILE_WRITE(oldfile, rfpts(1,1), n*lrf,
     >                  MPI_DOUBLE_PRECISION,
     >                  MPI_STATUS_IGNORE, ierr) 

         call MPI_FILE_CLOSE(oldfile, ierr) 

      endif

! ----------------------------------------
! Move particles back to where they belong
! ----------------------------------------
      call move_particles_inproc

      return
      end
c----------------------------------------------------------------------
      subroutine read_parallel_restart_part
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'
      include 'mpif.h'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      character*16 locstring, datastring
      integer*8    disp, stride_len 
      integer      status_mpi(MPI_STATUS_SIZE)
      integer      prevs(0:np-1),npt_total,e,oldfile
      integer      color,particle_io_comm,key

      rpi    = 4.0d+0*atan(1.d+0) ! pi

c     setup files to write to mpi 
      write(locstring,'(A8,I5.5,A3)') 'rpartxyz', ipart_restartr, '.3D' 
      write(datastring,'(A9,I5.5)')   'rpartdata', ipart_restartr

      open(364, file=datastring, action="read")
         read(364,*) nw
      close(364)

! ----------------------------------------
! Calculate how many processors to be used
! ----------------------------------------
      npmax_set = int(nw/llpart) + 1 ! use as few procs as possible
      npmax = np
      if (npmax .gt. npmax_set) npmax = npmax_set

! ----------------------------------------
! Create communicator with that many ranks
! ----------------------------------------
      color = MPI_UNDEFINED
      if (nid.lt.npmax) color = 0
      key=nid
      call MPI_Comm_split(nekcomm,color,key,particle_io_comm,ierr)

! -----------------------------------------
! Only use npmax ranks to read in particles
! -----------------------------------------
      if (nid .lt. npmax) then

! -------------------------------------------------------
! Each rank computes the number of particles it will read
! -------------------------------------------------------
         nnp       = int(nw/npmax) 
         if (nid.gt.npmax-1) nnp = 0
         nw_tmp    = nnp*npmax
         ndef      = nw - nw_tmp
         if (nid .lt. ndef) nnp = nnp + 1
         
! -----------------------------------------------
! Each rank sends how many particles it will read
! -----------------------------------------------
         call MPI_Send(nnp, 1, MPI_INTEGER,0,0, particle_io_comm, ierr)
         
         ! keep track of how many particles are on previous procs
         if (nid.eq. 0) then
             prevs(0) = nnp
             do i=0,npmax-1
                call MPI_Recv(prevs(i),1,MPI_INTEGER,i,
     >                        0,particle_io_comm,status_mpi,ierr)
             enddo
         endif
         call MPI_BCAST(prevs,npmax,MPI_INTEGER,0,particle_io_comm,ierr)

! ------------------------------------------------------
! Each rank counts how many particles in ranks before it
! ------------------------------------------------------
         stride_len = 0
         if (nid .ne. 0) then
         do i=1,nid
            stride_len = stride_len + prevs(i-1)
         enddo
         endif

! -------------------------
! Parallel MPI file read in
! -------------------------
         call MPI_FILE_OPEN(particle_io_comm, locstring,
     >                   MPI_MODE_RDONLY, 
     >                   MPI_INFO_NULL, oldfile, ierr) 
   
         disp = stride_len*lrf*8  ! lrf properties each with 8 bytes
         call MPI_FILE_SET_VIEW(oldfile, disp, MPI_DOUBLE_PRECISION,
     >                       MPI_DOUBLE_PRECISION, "native", 
     >                       MPI_INFO_NULL, ierr) 
         call MPI_FILE_READ(oldfile, rfpts(1,1), nnp*lrf,
     >                  MPI_DOUBLE_PRECISION,
     >                  MPI_STATUS_IGNORE, ierr) 

         call MPI_FILE_CLOSE(oldfile, ierr) 

! ----------------------------------------
! Assign values read in to rpart and ipart
! ----------------------------------------
         i = n ! if there are previous particles
         do ii = 1,nnp
            i = n + ii
            rpart(jx,i)    =  rfpts(1,ii) 
            rpart(jy,i)    =  rfpts(2,ii)
            rpart(jz,i)    =  rfpts(3,ii)  
            rpart(jx1+0,i) =  rfpts(4,ii)  
            rpart(jx1+1,i) =  rfpts(5,ii)
            rpart(jx1+2,i) =  rfpts(6,ii)  
            rpart(jx2+0,i) =  rfpts(7,ii)  
            rpart(jx2+1,i) =  rfpts(8,ii)
            rpart(jx2+2,i) =  rfpts(9,ii)  
            rpart(jx3+0,i) =  rfpts(10,ii)
            rpart(jx3+1,i) =  rfpts(11,ii)
            rpart(jx3+2,i) =  rfpts(12,ii)
                                           
            rpart(jv0,i)   =  rfpts(13,ii)
            rpart(jv0+1,i) =  rfpts(14,ii)
            rpart(jv0+2,i) =  rfpts(15,ii)
            rpart(jv1+0,i) =  rfpts(16,ii)
            rpart(jv1+1,i) =  rfpts(17,ii)
            rpart(jv1+2,i) =  rfpts(18,ii)
            rpart(jv2+0,i) =  rfpts(19,ii)
            rpart(jv2+1,i) =  rfpts(20,ii)
            rpart(jv2+2,i) =  rfpts(21,ii)
            rpart(jv3+0,i) =  rfpts(22,ii)
            rpart(jv3+1,i) =  rfpts(23,ii)
            rpart(jv3+2,i) =  rfpts(24,ii)
                                           
            rpart(ju0,i)   =  rfpts(25,ii)
            rpart(ju0+1,i) =  rfpts(26,ii)
            rpart(ju0+2,i) =  rfpts(27,ii)
            rpart(ju1+0,i) =  rfpts(28,ii)
            rpart(ju1+1,i) =  rfpts(29,ii)
            rpart(ju1+2,i) =  rfpts(30,ii)
            rpart(ju2+0,i) =  rfpts(31,ii)
            rpart(ju2+1,i) =  rfpts(32,ii)
            rpart(ju2+2,i) =  rfpts(33,ii)
            rpart(ju3+0,i) =  rfpts(34,ii)
            rpart(ju3+1,i) =  rfpts(35,ii)
            rpart(ju3+2,i) =  rfpts(36,ii)
                                           
            rpart(jdp,i)   =  rfpts(37,ii)
            rpart(jspl,i)  =  rfpts(38,ii)
            rpart(jtemp,i) =  rfpts(39,ii)
            ipart(jpid1,i) =  nint(rfpts(40,ii))
            ipart(jpid2,i) =  nint(rfpts(41,ii))
            ipart(jpid3,i) =  nint(rfpts(42,ii))
        
            ! extra stuff
            rpart(jtaup,i) = rpart(jdp,i)**2*rho_p/18.0d+0/mu_0
            rpart(jrhop,i) = rho_p      ! material density of particle
            rpart(jvol,i)  = rpart(jspl,i)*rpi*rpart(jdp,i)**3/6.d+0! particle volume
            rpart(jgam,i)  = 1.          ! initial integration correction
            rpart(jrpe,i)  = rpart(jspl,i)**(1./3.)*rpart(jdp,i)/2.
         enddo
         n = i
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine output_parallel_lagrangian_parts
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'
      include 'mpif.h'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*15 locstring, datastring
      integer*8    disp, stride_len 
      integer      status_mpi(MPI_STATUS_SIZE)
      integer      prevs(0:np-1),npt_total,e,oldfile
      real         one_t

      rpi    = 4.0*atan(1.) ! pi

c     setup files to write to mpi 
      icalld = icalld+1
      write(locstring,'(A7,I5.5,A3)') 'partxyz', icalld, '.3D' 
      write(datastring,'(A8,I5.5)')   'partdata', icalld

      ! these are the values that will be output in .3D binary files
      ! kind of a dumb way to do it, but good for memory  and needed for
      ! striding...
      icount = 0
      do i = 1,n
         icount = icount + 1
         rfpts(icount,1) = rpart(jx,i)
         icount = icount + 1
         rfpts(icount,1) = rpart(jy,i)
         icount = icount + 1
         rfpts(icount,1) = rpart(jz,i)
         icount = icount + 1
         rfpts(icount,1) = 2.*rpart(jrpe,i) ! output diameter
      enddo
     
      call MPI_Send(n, 1, MPI_INTEGER, 0, 0, nekcomm, ierr)
      npt_total = iglsum(n,1)

      ! keep track of how many particles are on previous procs
      if (nid.eq. 0) then
c         output data so files can be easily converted to binary
          open(364, file=datastring, action="write")
             write(364,*) npt_total
          close(364)

          prevs(0) = n
          do i=1,np-1
             call MPI_Recv(prevs(i),1,MPI_INTEGER,i,
     >                        0,nekcomm,status_mpi,ierr)
          enddo
      endif
      call MPI_BCAST(prevs,np, MPI_INTEGER,0,nekcomm,ierr) 

      stride_len = 0
      if (nid .ne. 0) then
      do i=1,nid
         stride_len = stride_len + prevs(i-1)
      enddo
      endif

      call MPI_FILE_OPEN(nekcomm, locstring,
     >                   MPI_MODE_CREATE + MPI_MODE_WRONLY, 
     >                   MPI_INFO_NULL, oldfile, ierr) 
   
      disp = stride_len*4*8  ! 4 properties each with 8 bytes
      call MPI_FILE_SET_VIEW(oldfile, disp, MPI_DOUBLE_PRECISION,
     >                       MPI_DOUBLE_PRECISION, "native", 
     >                       MPI_INFO_NULL, ierr) 
      call MPI_FILE_WRITE(oldfile, rfpts(1,1), n*4,
     >                  MPI_DOUBLE_PRECISION,
     >                  MPI_STATUS_IGNORE, ierr) 

      call MPI_FILE_CLOSE(oldfile, ierr) 

      return
      end
c----------------------------------------------------------------------
      subroutine output_two_way_io
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'
      include 'mpif.h'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real rfpfluid(3),rfpfluidl(3),msum,msum_tot(3,2)

      itmp = 1
      call outpost2(ptw(1,1,1,1,1),         ! fhyd_x
     >              ptw(1,1,1,1,2),         ! fhyd_y
     >              ptw(1,1,1,1,3),         ! fhyd_z
     >              ptw(1,1,1,1,4),         ! phi_p (but not if lx1!=lx2
     >              ptw(1,1,1,1,4),         ! phi_p
     >              itmp          ,        
     >              'ptw')

c     eulerian integrations -----------------------------------------
c     fluid momentum 
      msum_tot(1,1) = glsc3(bm1,vtrans,vx,nx1*ny1*nz1*nelv)
      msum_tot(2,1) = glsc3(bm1,vtrans,vy,nx1*ny1*nz1*nelv)
      msum_tot(3,1) = glsc3(bm1,vtrans,vz,nx1*ny1*nz1*nelv)
c     particle volume fraction
      vf_part_e     = glsc2(bm1,ptw(1,1,1,1,4),nx1*ny1*nz1*nelt)
c     particle forces on fluid
      rfpfluid(1)   = glsc2(bm1,ptw(1,1,1,1,1),nx1*ny1*nz1*nelt)
      rfpfluid(2)   = glsc2(bm1,ptw(1,1,1,1,2),nx1*ny1*nz1*nelt)
      rfpfluid(3)   = glsc2(bm1,ptw(1,1,1,1,3),nx1*ny1*nz1*nelt)
c     lagrangian integrations ---------------------------------------
c     particle momentum
      do ieq=0,2
         msum = 0.0
         rsum = 0.0
         do i=1,n
           msum = msum + 
     >       rpart(jv0+ieq,i)*rpart(jrhop,i)*rpart(jvol,i)
           rsum = rsum + rpart(jf0+ieq,i)
         enddo
         msum_tot(ieq+1,2) = glsum(msum,1)
         rfpfluidl(1+ieq)  = glsum(rsum,1)
      enddo
c     particle volume fraction
      msum = 0.0
      do i=1,n
         msum = msum + rpart(jvol,i)
      enddo
      vf_part_l = glsum(msum,1)

c     print to files ------------------------------------------------
c     print properties to logfile
      if (nid.eq.0) write(6,500) "--- Eulerian Properties ------"
      if (nid.eq.0) write(6,500) "Fluid Momentum :              ", 
     >                  msum_tot(1,1),msum_tot(2,1),msum_tot(3,1)
      if (nid.eq.0) write(6,500) "Particle forces:              ", 
     >                  rfpfluid(1),rfpfluid(2),rfpfluid(3)         
      if (nid.eq.0) write(6,500) "Particle Volume:              ", 
     >                  vf_part_e
      if (nid.eq.0) write(6,500) "--- Lagrangian Properties --- "
      if (nid.eq.0) write(6,500) "Particle Momentum :           ", 
     >                  msum_tot(1,2),msum_tot(2,2),msum_tot(3,2)
      if (nid.eq.0) write(6,500) "Particle forces:              ", 
     >                  rfpfluidl(1),rfpfluidl(2),rfpfluidl(3)         
      if (nid.eq.0) write(6,500) "Particle Volume:              ", 
     >                  vf_part_l

  500 FORMAT(A30,9ES20.10)

      return
      end
c----------------------------------------------------------------------
      subroutine output_along_line_avg_x
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*15 outstring
      integer*8 disp, stride_len 
      real send_vals(1+50)
      real rxgls(lx1),uf(lx1,ly1,lz1,lelt,50),rcount(50),rdum(50),
     >     rtmp(50)
      integer nlfl

      common /running_avgs/ rec_vals
      real rec_vals(1+50,8000*15) !8000 elements, by nx1=15 max

      nlfl = 21 

      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         uf(i,j,k,ie,1) = ptw(i,j,k,ie,1)
         uf(i,j,k,ie,2) = ptw(i,j,k,ie,2)
         uf(i,j,k,ie,3) = ptw(i,j,k,ie,3)
         uf(i,j,k,ie,4) = ptw(i,j,k,ie,4)
         uf(i,j,k,ie,5) = ptw(i,j,k,ie,5)
         uf(i,j,k,ie,6) = ptw(i,j,k,ie,6)
         uf(i,j,k,ie,7) = ptw(i,j,k,ie,7) 
         uf(i,j,k,ie,8) = ptw(i,j,k,ie,8) 
         uf(i,j,k,ie,9) = rhs_fluidp(i,j,k,ie,1)
         uf(i,j,k,ie,10)= rhs_fluidp(i,j,k,ie,2)
         uf(i,j,k,ie,11)= rhs_fluidp(i,j,k,ie,3)
         uf(i,j,k,ie,12)= rhs_fluidp(i,j,k,ie,4)
         uf(i,j,k,ie,13)= rhs_fluidp(i,j,k,ie,5)
         uf(i,j,k,ie,14)= rhs_fluidp(i,j,k,ie,6)
         uf(i,j,k,ie,15)= rhs_fluidp(i,j,k,ie,7)
         uf(i,j,k,ie,16)= vx(i,j,k,ie)
         uf(i,j,k,ie,17)= vy(i,j,k,ie)
         uf(i,j,k,ie,18)= vz(i,j,k,ie)
         uf(i,j,k,ie,19)= pr(i,j,k,ie) ! not set up for lx2 mesh!!
         uf(i,j,k,ie,20)= vtrans(i,j,k,ie,1)
         uf(i,j,k,ie,21)= t(i,j,k,ie,1)
      enddo
      enddo
      enddo
      enddo

      icalld = icalld+1
      write(outstring,'(A9,I5.5)') 'avgsdatax', icalld

      rlengx = xm1(nx1,1,1,1) - xm1(1,1,1,1)
      do i=1,nx1 
         rxgls(i) = (xgll(i) + 1.)*rlengx/2.
      enddo
      
      rthresh = rlengx/100.
      rxs = xdrange(1,1)
      rxt = rxs
      icm = 1

      do while (rxs .le. xdrange(2,1))

         do i=1,nx1
            rxs = rxt + rxgls(i)

            call rzero(rdum,nlfl)
            call rzero(rcount,nlfl)

            do ie=1,nelt
            do ik=1,nz1
            do ij=1,ny1
            do ii=1,nx1
               rxv = xm1(ii,ij,ik,ie)
               if (abs(rxv - rxs) .lt. rthresh) then
                  do j = 1,nlfl
                     rdum(j) = rdum(j) + uf(ii,ij,ik,ie,j)
                     rcount(j) = rcount(j) + 1.
                  enddo
               endif
            enddo
            enddo
            enddo
            enddo

            isz = 1
            rec_vals(1,icm) = rxs

            do j = 1,nlfl
               rec_vals(j+1,icm) =  glsum(rdum(j),1)
               rtmp(j) = glsum(rcount(j),1)
               rec_vals(j+1,icm) = rec_vals(j+1,icm)/rtmp(j)
            enddo

            icm = icm + 1
         enddo

         rxt = rxs
      enddo

      if (nid.eq. 0) then
          open(364, file=outstring, action="write",position="append")
          do i =1,icm-1 ! last point has issues
             write(364,600) rec_vals(1,i),
     >                      rec_vals(2,i),
     >                      rec_vals(3,i),
     >                      rec_vals(4,i),
     >                      rec_vals(5,i),
     >                      rec_vals(6,i),
     >                      rec_vals(7,i),
     >                      rec_vals(8,i),
     >                      rec_vals(9,i),
     >                      rec_vals(10,i),
     >                      rec_vals(11,i),
     >                      rec_vals(12,i),
     >                      rec_vals(13,i),
     >                      rec_vals(14,i),
     >                      rec_vals(15,i),
     >                      rec_vals(16,i),
     >                      rec_vals(17,i),
     >                      rec_vals(18,i),
     >                      rec_vals(19,i),
     >                      rec_vals(20,i),
     >                      rec_vals(21,i),
     >                      rec_vals(22,i) ! add one
          enddo
          close(364)
      endif
      
  600 FORMAT(22ES20.10)
      return
      end
c----------------------------------------------------------------------
      subroutine output_along_line_avg_y
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*15 outstring
      integer*8 disp, stride_len 
      real send_vals(1+50)
      real rygls(ly1),uf(lx1,ly1,lz1,lelt,50),rcount(50),rdum(50),
     >     rtmp(50)
      integer nlfl

      real*8 rlengy

      common /running_avgs/ rec_vals
      real rec_vals(1+50,8000*15) !8000 elements, by nx1=15 max

      nlfl = 21 

      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         uf(i,j,k,ie,1) = ptw(i,j,k,ie,1)
         uf(i,j,k,ie,2) = ptw(i,j,k,ie,2)
         uf(i,j,k,ie,3) = ptw(i,j,k,ie,3)
         uf(i,j,k,ie,4) = ptw(i,j,k,ie,4)
         uf(i,j,k,ie,5) = ptw(i,j,k,ie,5)
         uf(i,j,k,ie,6) = ptw(i,j,k,ie,6)
         uf(i,j,k,ie,7) = ptw(i,j,k,ie,7) 
         uf(i,j,k,ie,8) = ptw(i,j,k,ie,8) 
         uf(i,j,k,ie,9) = rhs_fluidp(i,j,k,ie,1)
         uf(i,j,k,ie,10)= rhs_fluidp(i,j,k,ie,2)
         uf(i,j,k,ie,11)= rhs_fluidp(i,j,k,ie,3)
         uf(i,j,k,ie,12)= rhs_fluidp(i,j,k,ie,4)
         uf(i,j,k,ie,13)= rhs_fluidp(i,j,k,ie,5)
         uf(i,j,k,ie,14)= rhs_fluidp(i,j,k,ie,6)
         uf(i,j,k,ie,15)= rhs_fluidp(i,j,k,ie,7)
         uf(i,j,k,ie,16)= vx(i,j,k,ie)
         uf(i,j,k,ie,17)= vy(i,j,k,ie)
         uf(i,j,k,ie,18)= vz(i,j,k,ie)
         uf(i,j,k,ie,19)= pr(i,j,k,ie) ! not set up for lx2 mesh!!
         uf(i,j,k,ie,20)= vtrans(i,j,k,ie,1)
         uf(i,j,k,ie,21)= t(i,j,k,ie,1)
      enddo
      enddo
      enddo
      enddo

      icalld = icalld+1
      write(outstring,'(A9,I5.5)') 'avgsdatay', icalld

      rlengy = ym1(1,ny1,1,1) - ym1(1,1,1,1)
      rlengy = glmin(rlengy,1)
      do i=1,ny1 
         rygls(i) = (ygll(i) + 1.)*rlengy/2.
      enddo

      rthresh = rlengy/100.
      rys = xdrange(1,2)
      ryt = rys
      icm = 1

      ! DZ FAKE
c     do while (rys .le. xdrange(2,2))
      do while (rys .le. 0.23)

         do i=1,ny1
            rys = ryt + rygls(i)

            call rzero(rdum,nlfl)
            call rzero(rcount,nlfl)

            do ie=1,nelt
            do ik=1,nz1
            do ij=1,ny1
            do ii=1,nx1
               ryv = ym1(ii,ij,ik,ie)
               if (abs(ryv - rys) .lt. rthresh) then
                  do j = 1,nlfl
                     rdum(j) = rdum(j) + uf(ii,ij,ik,ie,j)
                     rcount(j) = rcount(j) + 1.
                  enddo
               endif
            enddo
            enddo
            enddo
            enddo

            isz = 1
            rec_vals(1,icm) = rys

            do j = 1,nlfl
               rec_vals(j+1,icm) =  glsum(rdum(j),1)
               rtmp(j) = glsum(rcount(j),1)
               rec_vals(j+1,icm) = rec_vals(j+1,icm)/rtmp(j)
            enddo

            icm = icm + 1
         enddo

         ryt = rys
      enddo

      if (nid.eq. 0) then
          open(365, file=outstring, action="write",position="append")
          do i =1,icm-1 ! last point has issues
             write(365,600) rec_vals(1,i),
     >                      rec_vals(2,i),
     >                      rec_vals(3,i),
     >                      rec_vals(4,i),
     >                      rec_vals(5,i),
     >                      rec_vals(6,i),
     >                      rec_vals(7,i),
     >                      rec_vals(8,i),
     >                      rec_vals(9,i),
     >                      rec_vals(10,i),
     >                      rec_vals(11,i),
     >                      rec_vals(12,i),
     >                      rec_vals(13,i),
     >                      rec_vals(14,i),
     >                      rec_vals(15,i),
     >                      rec_vals(16,i),
     >                      rec_vals(17,i),
     >                      rec_vals(18,i),
     >                      rec_vals(19,i),
     >                      rec_vals(20,i),
     >                      rec_vals(21,i),
     >                      rec_vals(22,i) ! add one
          enddo
          close(365)
      endif
      
  600 FORMAT(22ES20.10)
      return
      end
c----------------------------------------------------------------------
      subroutine output_along_line_avg_z
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'CMTDATA'
      include 'CMTPART'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      integer icalld
      save    icalld
      data    icalld  /-1/

      character*15 outstring
      integer*8 disp, stride_len 
      real send_vals(1+50)
      real rzgls(lz1),uf(lx1,ly1,lz1,lelt,50),rcount(50),rdum(50),
     >     rtmp(50)
      integer nlfl

      common /running_avgs/ rec_vals
      real rec_vals(1+50,8000*15) !8000 elements, by nx1=15 max

      nlfl = 21 

      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         uf(i,j,k,ie,1) = ptw(i,j,k,ie,1)
         uf(i,j,k,ie,2) = ptw(i,j,k,ie,2)
         uf(i,j,k,ie,3) = ptw(i,j,k,ie,3)
         uf(i,j,k,ie,4) = ptw(i,j,k,ie,4)
         uf(i,j,k,ie,5) = ptw(i,j,k,ie,5)
         uf(i,j,k,ie,6) = ptw(i,j,k,ie,6)
         uf(i,j,k,ie,7) = ptw(i,j,k,ie,7) 
         uf(i,j,k,ie,8) = ptw(i,j,k,ie,8) 
         uf(i,j,k,ie,9) = rhs_fluidp(i,j,k,ie,1)
         uf(i,j,k,ie,10)= rhs_fluidp(i,j,k,ie,2)
         uf(i,j,k,ie,11)= rhs_fluidp(i,j,k,ie,3)
         uf(i,j,k,ie,12)= rhs_fluidp(i,j,k,ie,4)
         uf(i,j,k,ie,13)= rhs_fluidp(i,j,k,ie,5)
         uf(i,j,k,ie,14)= rhs_fluidp(i,j,k,ie,6)
         uf(i,j,k,ie,15)= rhs_fluidp(i,j,k,ie,7)
         uf(i,j,k,ie,16)= vx(i,j,k,ie)
         uf(i,j,k,ie,17)= vy(i,j,k,ie)
         uf(i,j,k,ie,18)= vz(i,j,k,ie)
         uf(i,j,k,ie,19)= pr(i,j,k,ie) ! not set up for lx2 mesh!!
         uf(i,j,k,ie,20)= vtrans(i,j,k,ie,1)
         uf(i,j,k,ie,21)= t(i,j,k,ie,1)
      enddo
      enddo
      enddo
      enddo

      icalld = icalld+1
      write(outstring,'(A9,I5.5)') 'avgsdataz', icalld

      rlengz = zm1(1,1,nz1,1) - zm1(1,1,1,1)
      do i=1,nz1 
         rzgls(i) = (zgll(i) + 1.)*rlengz/2.
      enddo

      rthresh = rlengz/100.
      rzs = xdrange(1,3)
      rzt = rzs
      icm = 1

      do while (rzs .le. xdrange(2,3))

         do i=1,nz1
            rzs = rzt + rzgls(i)

            call rzero(rdum,nlfl)
            call rzero(rcount,nlfl)

            do ie=1,nelt
            do ik=1,nz1
            do ij=1,ny1
            do ii=1,nx1
               rzv = zm1(ii,ij,ik,ie)
               if (abs(rzv - rzs) .lt. rthresh) then
                  do j = 1,nlfl
                     rdum(j) = rdum(j) + uf(ii,ij,ik,ie,j)
                     rcount(j) = rcount(j) + 1.
                  enddo
               endif
            enddo
            enddo
            enddo
            enddo

            isz = 1
            rec_vals(1,icm) = rzs

            do j = 1,nlfl
               rec_vals(j+1,icm) =  glsum(rdum(j),1)
               rtmp(j) = glsum(rcount(j),1)
               rec_vals(j+1,icm) = rec_vals(j+1,icm)/rtmp(j)
            enddo

            icm = icm + 1
         enddo

         rzt = rzs
      enddo

      if (nid.eq. 0) then
          open(366, file=outstring, action="write",position="append")
          do i =1,icm-1 ! last point has issues
             write(366,600) rec_vals(1,i),
     >                      rec_vals(2,i),
     >                      rec_vals(3,i),
     >                      rec_vals(4,i),
     >                      rec_vals(5,i),
     >                      rec_vals(6,i),
     >                      rec_vals(7,i),
     >                      rec_vals(8,i),
     >                      rec_vals(9,i),
     >                      rec_vals(10,i),
     >                      rec_vals(11,i),
     >                      rec_vals(12,i),
     >                      rec_vals(13,i),
     >                      rec_vals(14,i),
     >                      rec_vals(15,i),
     >                      rec_vals(16,i),
     >                      rec_vals(17,i),
     >                      rec_vals(18,i),
     >                      rec_vals(19,i),
     >                      rec_vals(20,i),
     >                      rec_vals(21,i),
     >                      rec_vals(22,i) ! add one
          enddo
          close(366)
      endif
      
  600 FORMAT(22ES20.10)
      return
      end
c----------------------------------------------------------------------
      subroutine read_particle_input
      include 'SIZE'
      include 'INPUT'
      include 'CMTTIMERS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'CMTPART'

      character*72 dum_str

      open(unit=81,file="particles.inp",form="formatted")

c - - - - PARTICLE PHYSICAL PROPERTIES - - - - - - - - - - - - - - - - 
      read(81,*) dum_str
      read(81,*) nw 
      read(81,*) dp(1), dp(2)
      read(81,*) tp_0
      read(81,*) rho_p
      read(81,*) cp_p
c - - - - PARTICLE FORCE MODELS - - - - - - - - - - - - - - - - - - - -
      read(81,*) dum_str
      read(81,*) part_force(1)
      read(81,*) part_force(2)
      read(81,*) part_force(3)
      read(81,*) part_force(4)
      read(81,*) part_force(5)
c - - - - USER OPTIONS - - - - - - - - - - - - - - - - - - - - - - - - 
      read(81,*) dum_str
      read(81,*) time_integ
      read(81,*) two_way
      read(81,*) red_interp
      read(81,*) npio_method
      read(81,*) inject_rate
      read(81,*) time_delay
      read(81,*) nrandseed
c - - - - PROJECTION OPTIONS - - - - - - - - - - - - - - - - - - - - - 
      read(81,*) dum_str
      read(81,*) npro_method
      read(81,*) nitspl
      read(81,*) rspl
      read(81,*) dfilt
      read(81,*) ralphdecay
c - - - - BOUNDARY TREATMENT - - - - - - - - - - - - - - - - - - - - - 
      read(81,*) dum_str
      read(81,*) nrect_assume
      read(81,*) bc_part(1), bc_part(2)
      read(81,*) bc_part(3), bc_part(4)
      read(81,*) bc_part(5), bc_part(6)
c - - - - RESTARTING - - - - - - - - - - - - - - - - - - - - - - - - - 
      read(81,*) dum_str
      read(81,*) ipart_restarto
      read(81,*) ipart_restartr
c - - - - DEM - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      read(81,*) dum_str
      read(81,*) ksp
      read(81,*) e_rest

      close(81)

      mu_0   = abs(param(2))

      return
      end
c----------------------------------------------------------------------
      subroutine output_particle_diagnostics
      include 'SIZE'
      include 'CMTTIMERS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'CMTPART'

      real rdiags_part(3,2), rvels(4)

      do i=1,3
         rdiags_part(i,1) =  1E8
         rdiags_part(i,2) = -1E8
      enddo
      do i=1,4
         rvels(i) = 0.0
      enddo
      
      do i=1,n
         rvx = rpart(jv0,i)
         rvy = rpart(jv0+1,i)
         rvz = rpart(jv0+2,i)
         rvmag = sqrt(rvx**2 + rvy**2 + rvz**2)

         rxx = rpart(jx,i)
         ryy = rpart(jy,i)
         rzz = rpart(jz,i)
         rrp = rpart(jrpe,i)

         ! maxes
         if (rxx + rrp .gt. rdiags_part(1,2)) rdiags_part(1,2) = rxx+rrp
         if (ryy + rrp .gt. rdiags_part(2,2)) rdiags_part(2,2) = ryy+rrp
         if (rzz + rrp .gt. rdiags_part(3,2)) rdiags_part(3,2) = rzz+rrp

         ! mins
         if (rxx - rrp .lt. rdiags_part(1,1)) rdiags_part(1,1) = rxx-rrp
         if (ryy - rrp .lt. rdiags_part(2,1)) rdiags_part(2,1) = ryy-rrp
         if (rzz - rrp .lt. rdiags_part(3,1)) rdiags_part(3,1) = rzz-rrp

         ! velocities
         if ( abs(rvx)   .gt. abs(rvels(1)) )  rvels(1) = rvx
         if ( abs(rvy)   .gt. abs(rvels(2)) )  rvels(2) = rvy
         if ( abs(rvz)   .gt. abs(rvels(3)) )  rvels(3) = rvz
         if ( abs(rvmag) .gt. abs(rvels(4)) )  rvels(4) = rvmag
      enddo

      ! compute globally now
      do i=1,3
         rdum = rdiags_part(i,1)
         rdiags_part(i,1) =  glmin(rdum,1)
         rdum = rdiags_part(i,2)
         rdiags_part(i,2) =  glmax(rdum,1)
      enddo
      do i=1,4
         rdum  = rvels(i)
         rdum1 = glmin(rdum,1)
         rdum  = rvels(i)
         rdum2 = glmax(rdum,1)
         rvels(i) = rdum1
         if (abs(rdum1) .lt. abs(rdum2)) rvels(i) = rdum2
      enddo

      if (nid .eq. 0) then
       write(6,*)'----- START PARTICLE DIAGNOSTICS: -----'
       write(6,*)'XMIN,XMAX       :',rdiags_part(1,1),rdiags_part(1,2)
       write(6,*)'YMIN,YMAX       :',rdiags_part(2,1),rdiags_part(2,2)
       write(6,*)'ZMIN,ZMAX       :',rdiags_part(3,1),rdiags_part(3,2)
       write(6,*)'MAX(VX,VY,VZ,|V|:',rvels(1),rvels(2),rvels(3),rvels(4)
       write(6,*)'----- END PARTICLE DIAGNOSTICS: -----'
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine output_particle_timers
      include 'SIZE'
      include 'CMTTIMERS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'CMTPART'

      rftime_t = 0.
      rptime_t = 0.

      if(nid.eq.0) then
         write(6,*) 'TIMER H: ', istep
      endif

      do i = 1,iptlen
         rdum  = pttime(i)/istep
         dtime = glsum(rdum,1)
         rtime = dtime/np
         if(nid.eq.0)  write(6,*) 'TIMER #:',i,rtime

         ! fluid and particle total time: note i == iptlen is f_col
         if (i .eq. 1) rftime_t = rtime
         if ((i .gt. 1) .and. (i.ne.iptlen)) rptime_t = rptime_t +
     >                                                  rtime
      enddo

      
      if (nid.eq.0) then
         write (6,*) 'TOTAL F:', rftime_t
         write (6,*) 'TOTAL P:', rptime_t
      endif


      return
      end
c----------------------------------------------------------------------
c     effeciently move particles between processors routines
c----------------------------------------------------------------------
      subroutine move_particles_inproc
c     Interpolate fluid velocity at current xyz points and move
c     data to the processor that owns the points.
c     Input:    n = number of points on this processor
c     Output:   n = number of points on this processor after the move
c     Code checks for n > llpart and will not move data if there
c     is insufficient room.
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'CMTPART'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      logical partl         ! This is a dummy placeholder, used in cr()

      nl = 0                ! No logicals exchanged

c     if (icalld1.eq.0) then
      if (icalld1.eq.0 .or. (resetFindpts .eq. 1)) then
         tolin = 1.e-12
         if (wdsize.eq.4) tolin = 1.e-6
         call intpts_setup  (tolin,i_fp_hndl)
         call fgslib_crystal_setup (i_cr_hndl,nekcomm,np)

         if (resetFindpts .eq. 1) icalld = 0
      endif

      icalld1 = icalld1 + 1

      if (icalld .le. 2 .or. (resetFindpts .eq. 1)) then
         call particles_in_nid
         call fgslib_findpts(i_fp_hndl !  stride     !   call fgslib_findpts( ihndl,
     $           , ifpts(jrc,1),lif        !   $             rcode,1,
     $           , ifpts(jpt,1),lif        !   &             proc,1,
     $           , ifpts(je0,1),lif        !   &             elid,1,
     $           , rfpts(jr ,1),lrf        !   &             rst,ndim,
     $           , rfpts(jd ,1),lrf        !   &             dist,1,
     $           , rfpts(jx ,1),lrf        !   &             pts(    1),1,
     $           , rfpts(jy ,1),lrf        !   &             pts(  n+1),1,
     $           , rfpts(jz ,1),lrf ,nfpts)    !   &             pts(2*n+1),1,n)
         call update_findpts_info
      else
         call findpts_box
      endif

      nmax = iglmax(n,1)
      if (nmax.gt.llpart) then
         if (nid.eq.0) write(6,1) nmax,llpart
    1    format('WARNING: Max number of particles:',
     $   i9,'.  Not moving because llpart =',i9,'.')
      else
c        Move particle info to the processor that owns each particle
c        using crystal router in log P time:

         jps = jpid1-1     ! Pointer to temporary proc id for swapping
         do i=1,n        ! Can't use jpt because it messes up particle info
            ipart(jps,i) = ipart(jpt,i)
         enddo
         call fgslib_crystal_tuple_transfer(i_cr_hndl,n,llpart
     $              , ipart,ni,partl,nl,rpart,nr,jps)
c        Sort by element number - for improved local-eval performance
         call fgslib_crystal_tuple_sort    (i_cr_hndl,n 
     $              , ipart,ni,partl,nl,rpart,nr,je0,1)
      endif
         call reset_rst_part

      return
      end
c-----------------------------------------------------------------------
      subroutine reset_rst_part
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      do ip=1,n
         xloc = rpart(jx,ip)
         yloc = rpart(jy,ip)
         zloc = rpart(jz,ip)
         ie = ipart(je0,ip) + 1
         rloc = -1.0 + 2.0*(xloc - xerange(1,1,ie))/
     $          (xerange(2,1,ie)-xerange(1,1,ie))
         sloc = -1.0 + 2.0*(yloc - xerange(1,2,ie))/
     $          (xerange(2,2,ie)-xerange(1,2,ie))
         tloc = -1.0 + 2.0*(zloc - xerange(1,3,ie))/
     $          (xerange(2,3,ie)-xerange(1,3,ie))
         rpart(jr  ,ip) = rloc
         rpart(jr+1,ip) = sloc
         rpart(jr+2,ip) = tloc

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine findpts_box
c
c     this routine will create ghost particles by checking if particle
c     is within d2chk of element faces
c
c     ghost particle x,y,z list will be in rptsgp(jgpx,j),rptsgp(jgpy,j),
c     rptsgp(jgpz,j), while processor and local element id are in
c     iptsgp(jgppt,j) and iptsgp(jgpes,j)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      do i=1,n
         ie = ipart(je0,i) + 1

         iip1 = 0
         iip2 = 0
         iip3 = 0

         rxdum1 = rpart(jx,i) - xerange(1,1,ie)
         rxdum2 = rpart(jx,i) - xerange(2,1,ie)
         if ((rpart(jx,i) .lt. xerange(2,1,ie)) .and.
     >       (rpart(jx,i) .gt. xerange(1,1,ie))) goto 123
         if (abs(rxdum1).lt. abs(rxdum2)) iip1 = -1
         if (abs(rxdum2).lt. abs(rxdum1)) iip1 =  1
 123  continue

         rxdum1 = rpart(jy,i) - xerange(1,2,ie)
         rxdum2 = rpart(jy,i) - xerange(2,2,ie)
         if ((rpart(jy,i) .lt. xerange(2,2,ie)) .and.
     >       (rpart(jy,i) .gt. xerange(1,2,ie))) goto 124
         if (abs(rxdum1).lt. abs(rxdum2)) iip2 = -1
         if (abs(rxdum2).lt. abs(rxdum1)) iip2 =  1
 124  continue

         rxdum1 = rpart(jz,i) - xerange(1,3,ie)
         rxdum2 = rpart(jz,i) - xerange(2,3,ie)
         if ((rpart(jz,i) .lt. xerange(2,3,ie)) .and.
     >       (rpart(jz,i) .gt. xerange(1,3,ie))) goto 125
         if (abs(rxdum1).lt. abs(rxdum2)) iip3 = -1
         if (abs(rxdum2).lt. abs(rxdum1)) iip3 =  1
 125  continue

         itype = abs(iip1) + abs(iip2) + abs(iip3)
         if (itype .eq. 0) then
            ier  = ipart(je0,i)
            impi = ipart(jpt,i)
            goto 1511
         elseif (itype .eq. 1) then
            do j=1,3*nfacegp-2,3   ! faces
               ii = el_face_num(j) 
               jj = el_face_num(j+1) 
               kk = el_face_num(j+2) 

               if (ii .eq. iip1) then
               if (jj .eq. iip2) then
               if (kk .eq. iip3) then
                  jdum = j/3+1
                  ier=el_face_el_map(ie,jdum)
                  impi=el_face_proc_map(ie,jdum)
                  goto 1511
               endif
               endif
               endif
            enddo
         elseif (itype .eq. 2) then
            do j=1,3*nedgegp-2,3   ! edges
               ii = el_edge_num(j) 
               jj = el_edge_num(j+1) 
               kk = el_edge_num(j+2) 

               if (ii .eq. iip1) then
               if (jj .eq. iip2) then
               if (kk .eq. iip3) then
                  jdum = j/3+1
                  ier=el_edge_el_map(ie,jdum)
                  impi=el_edge_proc_map(ie,jdum)
                  goto 1511
               endif
               endif
               endif
            enddo
         elseif (itype .eq. 3) then
            do j=1,3*ncornergp-2,3   ! corners
               ii = el_corner_num(j) 
               jj = el_corner_num(j+1) 
               kk = el_corner_num(j+2) 

               if (ii .eq. iip1) then
               if (jj .eq. iip2) then
               if (kk .eq. iip3) then
                  jdum = j/3+1
                  ier=el_corner_el_map(ie,jdum)
                  impi=el_corner_proc_map(ie,jdum)
                  goto 1511
               endif
               endif
               endif
            enddo
         endif
         
 1511 continue
         ipart(jpt,i) = impi
         ipart(je0,i) = ier
         ipart(jrc,i) = 0
         rpart(jd,i)  = 1.0

         if (impi .lt. 0) rpart(jx,i) = 1E12 ! kill the particle
         if (ier .lt. 0) rpart(jx,i) = 1E12 ! kill the particle

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine particles_in_nid
      include 'SIZE'
      include 'CMTPART'

      integer icalld
      save    icalld
      data    icalld  /-1/

      icalld = icalld + 1

      nfpts = 0
      do ip = 1,n
         xloc = rpart(jx,ip)
         yloc = rpart(jy,ip)
         zloc = rpart(jz,ip)
         itest = 0
         do ie=1,nelt
            if (xloc.ge.xerange(1,1,ie).and.xloc.le.xerange(2,1,ie))then
            if (yloc.ge.xerange(1,2,ie).and.yloc.le.xerange(2,2,ie))then
            if (zloc.ge.xerange(1,3,ie).and.zloc.le.xerange(2,3,ie))then
                ipart(je0 ,ip) = ie-1
                ipart(jrc ,ip) = 0
                ipart(jpt ,ip) = nid
                rpart(jd  ,ip) = 1.0 
                rloc = -1.0 + 2.0*(xloc - xerange(1,1,ie))/
     $                 (xerange(2,1,ie)-xerange(1,1,ie))
                sloc = -1.0 + 2.0*(yloc - xerange(1,2,ie))/
     $                 (xerange(2,2,ie)-xerange(1,2,ie))
                tloc = -1.0 + 2.0*(zloc - xerange(1,3,ie))/
     $                 (xerange(2,3,ie)-xerange(1,3,ie))
                rpart(jr  ,ip) = rloc
                rpart(jr+1,ip) = sloc
                rpart(jr+2,ip) = tloc
                itest = 1
                goto 123
            endif
            endif
            endif
         enddo
         if (itest.eq.0)then
            nfpts = nfpts + 1
            call copy (rfpts(1,nfpts),rpart(1,ip),7) ! only copy 1st 7 up to jz
            call icopy(ifpts(1,nfpts),ipart(1,ip),3) ! only copy 1st 3
            ifpts(4,nfpts) = ip                      ! 4th is map
            if(nfpts.gt.llpart)then
               write(6,*)'Too many points crossing over ',
     $                      nfpts,llpart,nid
               call exitt
            endif
         endif
123      continue
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine update_findpts_info
      include 'SIZE'
      include 'CMTPART'

      do ifp = 1,nfpts
         call copy(rpart(1,ifpts(4,ifp)),rfpts(1,ifp),7) ! only copy 1st 7 up to jz
         call icopy(ipart(1,ifpts(4,ifp)),ifpts(1,ifp),3) ! only copy 1st 3
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine intpts_setup(tolin,ih)
c
c setup routine for interpolation tool
c tolin ... stop point seach interation if 1-norm of the step in (r,s,t) 
c           is smaller than tolin 
c
      include 'SIZE'
      include 'GEOM'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      tol = tolin
      if (tolin.lt.0) tol = 1e-13 ! default tolerance 

      n       = lx1*ly1*lz1*lelt 
      npt_max = 256
      nxf     = 2*nx1 ! fine mesh for bb-test
      nyf     = 2*ny1
      nzf     = 2*nz1
      bb_t    = 0.1 ! relative size to expand bounding boxes by
c
      if(nid.eq.0) write(6,*) 'initializing intpts(), tol=', tol
      call fgslib_findpts_setup(ih,nekcomm,npp,ndim,
     &                     xm1,ym1,zm1,nx1,ny1,nz1,
     &                     nelt,nxf,nyf,nzf,bb_t,n,n,
     &                     npt_max,tol)
c       
      return
      end
c----------------------------------------------------------------------
c     routines not used currently
c----------------------------------------------------------------------
      subroutine point_to_grid_corr_init
c
c
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'CMTDATA'
      include 'CMTPART'

      integer i,j,k,e,ip
      real    xx,yy,zz,msum

      common /point2gridc/ p2gc
      real   p2gc(lx1,ly1,lz1,lelt,4)

      if (nid.eq.0) write(6,*) 'Starting point_to_grid_corr_init'

c     local mpi rank effects
      do ie=1,nelt
         xs = xerange(1,1,ie)
         xe = xerange(2,1,ie)
         ys = xerange(1,2,ie)
         ye = xerange(2,2,ie)
         zs = xerange(1,3,ie)
         ze = xerange(2,3,ie)
         xdelta = (xe-xs)/(nx1-1)
         ydelta = (ye-ys)/(ny1-1)
         zdelta = (ze-zs)/(nz1-1)
         do k=1,nz1
            zz = zs + (k-1)*zdelta
         do j=1,ny1
            yy = ys + (j-1)*ydelta
         do i=1,nx1
            xx = xs + (i-1)*xdelta
c           p2gc(i,j,k,ie,1) = xx
c           p2gc(i,j,k,ie,2) = yy 
c           p2gc(i,j,k,ie,3) = zz 

            p2gc(i,j,k,ie,1) = xm1(i,j,k,ie)
            p2gc(i,j,k,ie,2) = ym1(i,j,k,ie) 
            p2gc(i,j,k,ie,3) = zm1(i,j,k,ie) 
         enddo
         enddo
         enddo
      enddo

      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         xx = p2gc(i,j,k,ie,1) 
         yy = p2gc(i,j,k,ie,2) 
         zz = p2gc(i,j,k,ie,3) 
         call compute_gamma_grid(ie,xx,yy,zz,p2gc(i,j,k,ie,4))
      enddo
      enddo
      enddo
      enddo

      if (nid.eq.0) write(6,*) 'Ending point_to_grid_corr_init'

      return
      end
c----------------------------------------------------------------------
      subroutine compute_gamma_grid(ie,xx,yy,zz,gam_val)
c
c
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'CMTDATA'
      include 'CMTPART'

      integer e,er
      real    msum,msum_total,pi,multfc
      real    xx,yy,zz,Lx,Ly,Lz,rdum(lx1,ly1,lz1)
      real    mesharound(81),gam_val,dumval(lx1,ly1,lz1,27)
      real    xgd(lx1,ly1,lz1),ygd(lx1,ly1,lz1),zgd(lx1,ly1,lz1)

      xs = xerange(1,1,ie)
      xe = xerange(2,1,ie)
      ys = xerange(1,2,ie)
      ye = xerange(2,2,ie)
      zs = xerange(1,3,ie)
      ze = xerange(2,3,ie)

      Lx = xe - xs
      Ly = ye - ys
      Lz = ze - zs

      mesharound = (/ 
     >                0. ,0. ,0. , ! 1
     >                -Lx,-Ly,0. , ! 2
     >                0. ,-Ly,0. , ! 3
     >                Lx ,-Ly,0. , ! 4
     >                -Lx,0. ,0. , ! 5
     >                Lx ,0. ,0. , ! 6
     >                -Lx,Ly ,0. , ! 7
     >                0. ,Ly ,0. , ! 8
     >                Lx ,Ly ,0. , ! 9
     >                0. ,0. ,-Lz, ! 10
     >                -Lx,-Ly,-Lz, ! 11
     >                0. ,-Ly,-Lz, ! 12
     >                Lx ,-Ly,-Lz, ! 13
     >                -Lx,0. ,-Lz, ! 14
     >                Lx ,0. ,-Lz, ! 15
     >                -Lx,Ly ,-Lz, ! 16
     >                0. ,Ly ,-Lz, ! 17
     >                Lx ,Ly ,-Lz, ! 18
     >                0. ,0. ,Lz , ! 19
     >                -Lx,-Ly,Lz , ! 20
     >                0. ,-Ly,Lz , ! 21
     >                Lx ,-Ly,Lz , ! 22
     >                -Lx,0. ,Lz , ! 23
     >                Lx ,0. ,Lz , ! 24
     >                -Lx,Ly ,Lz , ! 25
     >                0. ,Ly ,Lz , ! 26
     >                Lx ,Ly ,Lz   ! 27
     >                            /)

      call rzero(dumval,lx1*ly1*lz1*27)
      call rzero(rdum,lx1*ly1*lz1)

      pi       = 4.0d+0*atan(1.0d+0)
      multfc   = 1./(sqrt(2.*pi)**3 * rsig**3)
      rbexpi   = 1./(-2.*rsig**2)

c     ralphd   = 1E10   ! dummy so it will spread everywhere
c     ralphd2   = 1E10  ! dummy so it will spread everywhere

      ralphd    = d2chk(2)     ! assume all directions same!
      ralphd2   = ralphd**2

      do iie=1,27         
         ioff = (iie-1)*3
         do k=1,ny1
         do j=1,ny1
         do i=1,nx1
            xgd(i,j,k) = xm1(i,j,k,ie) + mesharound(ioff+1)
            ygd(i,j,k) = ym1(i,j,k,ie) + mesharound(ioff+2)
            zgd(i,j,k) = zm1(i,j,k,ie) + mesharound(ioff+3)
         enddo
         enddo
         enddo
         
         rpass = 1.*multfc

c        call point_to_grid(dumval(1,1,1,iie),1.,xx,yy,zz,1.,
c    >             xgd,ygd,zgd)

         call point_to_grid(rdum(1,1,1),rdum(1,1,1),
     >                  rdum(1,1,1),dumval(1,1,1,iie),
     >                  rdum(1,1,1),rdum(1,1,1),rdum(1,1,1),rdum(1,1,1),
     >                  xgd,ygd,zgd,
     >                  1.,1.,1.,rpass,1.,1.,1.,1.,xx,yy,zz,rbexpi,
     >                  ralphd,ralphd2)

c        call point_to_grid(fvalgx(1,1,1,e),fvalgy(1,1,1,e),
c    >                   fvalgz(1,1,1,e),fvalgv(1,1,1,e),
c    >                   xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),
c    >                   pvalpx,pvalpy,pvalpz,pvalpv,xx,yy,zz,rbexpi,
c    >                   ralph,ralph2)
c        do k=1,ny1
c        do j=1,ny1
c        do i=1,nx1
c           print *, i,j,k,xgd(i,j,k),ygd(i,j,k),zgd(i,j,k),
c    >               dumval(i,j,k,iie)
c        enddo
c        enddo
c        enddo
      enddo

      msum = 0.
      do iie=1,27
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            msum = msum + dumval(i,j,k,iie)*bm1(i,j,k,ie)
         enddo
         enddo
         enddo
      enddo
c     msum_total = glsum(msum,1)
      if (abs(msum) .lt. 1E-16) then
         gam_val = 1.
      else
         gam_val = 1./msum
      endif

      return
      end
c----------------------------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     $        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c Long period (> 2 ! 1018 ) random number generator of L’Ecuyer with 
c Bays-Durham shuffle and added safeguards. Returns a uniform random deviate 
c between 0.0 and 1.0 (exclusive of the endpoint values). 
c Call with idum a negative integer to initialize; thereafter, do not alter 
c idum between successive deviates in a sequence. RNMX should approximate the 
c largest floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then 
         idum1=max(-idum,1) 
         idum2=idum1
         do j=NTAB+8,1,-1
            k=idum1/IQ1
            idum1=IA1*(idum1-k*IQ1)-k*IR1 
            if (idum1.lt.0) idum1=idum1+IM1 
            if (j.le.NTAB) iv(j)=idum1
         enddo
         iy=iv(1) 
      endif
      k=idum1/IQ1 
      idum1=IA1*(idum1-k*IQ1)-k*IR1
      if (idum1.lt.0) idum1=idum1+IM1 
      k=idum2/IQ2 
      idum2=IA2*(idum2-k*IQ2)-k*IR2 
      if (idum2.lt.0) idum2=idum2+IM2 
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum1 
      if(iy.lt.1)iy=iy+IMM1 
      ran2=min(AM*iy,RNMX)
      return
      END
c----------------------------------------------------------------------
      function unif_random(rxl,rxr)
c
c     must initialize ran2 first
c
      real xl,xr,unif_random

      rdum       = ran2(2)
      rlen       = rxr - rxl
      unif_random= rxl + rdum*rlen

      return
      end
c-----------------------------------------------------------------------
      function unif_random_norm(rxl,rxr,rstd)
c
c     must initialize ran2 first
c
      real xl,xr,unif_random_norm,rstd,rxfne(1000),rcdf(1000)

      rmu = (rxr + rxl)/2.
      nxfn  = 1000
      rxlf  = rxl
      rxrf  = rxr
      rdxf  = (rxrf-rxlf)/(nxfn-1.)

      do i=1,nxfn
         rxfne(i) = rxlf + (i-1.)*rdxf
         rcdf(i)  = 0.5*(1. + erf((rxfne(i)-rmu)/(rstd*sqrt(2.))))
      enddo

      rdum = unif_random(0.,1.)

!     find lower min value for inverse sampling
      idum = 0
      rmin = 100.
      do i=1,nxfn
         if (abs(rdum - rcdf(i)) .lt. rmin) then
            rmin = abs(rdum -rcdf(i))
            idum = i
         endif
      enddo
      ml = idum
      if (rdum .lt. rcdf(idum)) ml = ml + 1

      if (rdum .gt. rcdf(nxfn)) then
         unif_random_norm = rxrf 
      elseif (rdum .lt. rcdf(1)) then
         unif_random_norm = rxlf
      else
         rm = (rxfne(ml+1) - rxfne(ml))/(rcdf(ml+1) - rcdf(ml))
         unif_random_norm = rxfne(ml) + rm*(rdum - rcdf(ml))
      endif

      return
      end
c-----------------------------------------------------------------------
      function unif_random_cyl(rxl,rxr)
c
c     must initialize ran2 first
c
      parameter (nxfn = 10000)
      real xl,xr,unif_random_cyl,rstd,rxfne(nxfn),rcdf(nxfn)

      real    unif_random
      external unif_random

      rxlf  = rxl
      rxrf  = rxr
      rdxf  = (rxrf-rxlf)/(nxfn-1.)
      rdxf  = (1. - (rxlf/rxrf)**2)/(nxfn-1.)

      rnormalized_cdf = (rxr)**2

      do i=1,nxfn
         rxfne(i) = (rxlf/rxrf)**2 + (i-1.)*rdxf
         rcdf(i)  = (rxfne(i))**2 / rnormalized_cdf
      enddo

      rdum = unif_random((rxl/rxr)**2,rxr/rxr)

      unif_random_cyl = rxr*sqrt(rdum)

      return
      end
c----------------------------------------------------------------------
C> Compute coefficients for Runge-Kutta stages \cite{TVDRK}
      subroutine set_tstep_coef_part(dt_in)

      real tcoef(3,3),dt_cmt,time_cmt
      COMMON /TIMESTEPCOEF/ tcoef,dt_cmt,time_cmt

      real dt_in

      dt_cmt = dt_in

      tcoef(1,1) = 0.0
      tcoef(2,1) = 1.0 
      tcoef(3,1) = dt_cmt
      tcoef(1,2) = 3.0/4.0
      tcoef(2,2) = 1.0/4.0 
      tcoef(3,2) = dt_cmt/4.0 
      tcoef(1,3) = 1.0/3.0
      tcoef(2,3) = 2.0/3.0 
      tcoef(3,3) = dt_cmt*2.0/3.0 

      return
      end
c----------------------------------------------------------------------
      subroutine compute_phig_qtl(rdt_in,div)
c
c     Computes modified divergence constraint for multiphase dense
c     compressible flow
c
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      real div(lx2,ly2,lz2,lelv)

      common /phig_qtl_blk/ phig_last,phig_qtl
      real phig_last(lx1,ly1,lz1,lelt,4),phig_qtl(lx1,ly1,lz1,lelt)

      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1)

      integer icalld
      save    icalld
      data    icalld  /-1/

      icalld = icalld + 1

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         phig_last(ix,iy,iz,ie,1) = 1. - ptw(ix,iy,iz,ie,4)
         phig_qtl(ix,iy,iz,ie) = 0.
      enddo
      enddo
      enddo
      enddo

      if (icalld .eq. 0) then
         do ie=1,nelt
         do iz=1,nz1
         do iy=1,ny1
         do ix=1,nx1
            phig_last(ix,iy,iz,ie,2) = 1. - ptw(ix,iy,iz,ie,4)
            phig_last(ix,iy,iz,ie,3) = 1. - ptw(ix,iy,iz,ie,4)
            phig_last(ix,iy,iz,ie,4) = 1. - ptw(ix,iy,iz,ie,4)
         enddo
         enddo
         enddo
         enddo

         goto 123
      endif
      
      if (icalld .lt. 5) goto 123

      do ie=1,nelt
         call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                   phig_last(1,1,1,ie,1),lx1,if3d)
         
         do iz=1,nz1
         do iy=1,ny1
         do ix=1,nx1
            phig_qtl(ix,iy,iz,ie) = phig_last(ix,iy,iz,ie,1) -
     >                              phig_last(ix,iy,iz,ie,2)
            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie)/rdt_in

c           phig_qtl(ix,iy,iz,ie) =  (11./6.)*phig_last(ix,iy,iz,ie,1) 
c    >                              -(3.    )*phig_last(ix,iy,iz,ie,2) 
c    >                              +(3./2. )*phig_last(ix,iy,iz,ie,3) 
c    >                              -(1./3. )*phig_last(ix,iy,iz,ie,4) 
c           phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie)/rdt_in


            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vx(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dx
     >             (ur(ix,iy,iz)*RXM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SXM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TXM1(ix,iy,iz,ie))

            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vy(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dy
     >             (ur(ix,iy,iz)*RYM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SYM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TYM1(ix,iy,iz,ie))

            phig_qtl(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie) + 
     >          vz(ix,iy,iz,ie)/JACM1(ix,iy,iz,ie)* !d/dz
     >             (ur(ix,iy,iz)*RZM1(ix,iy,iz,ie) +
     >              us(ix,iy,iz)*SZM1(ix,iy,iz,ie) +
     >              ut(ix,iy,iz)*TZM1(ix,iy,iz,ie))


            phig_qtl(ix,iy,iz,ie) = -1.0*phig_qtl(ix,iy,iz,ie)/
     >                              phig_last(ix,iy,iz,ie,1)
         enddo
         enddo
         enddo
      enddo

  123 continue

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         phig_last(ix,iy,iz,ie,4) = phig_last(ix,iy,iz,ie,3)
         phig_last(ix,iy,iz,ie,3) = phig_last(ix,iy,iz,ie,2)
         phig_last(ix,iy,iz,ie,2) = phig_last(ix,iy,iz,ie,1)
      enddo
      enddo
      enddo
      enddo

      do ie=1,nelt
      do iz=1,nz1
      do iy=1,ny1
      do ix=1,nx1
         div(ix,iy,iz,ie) = phig_qtl(ix,iy,iz,ie)
      enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pre_sim_collisions
c
c     time stepping routine for pre-simulation collisions/settling
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'CMTDATA'
      include 'CMTPART'

      common /elementload/ gfirst, inoassignd, resetFindpts, pload(lelg)
      integer gfirst, inoassignd, resetFindpts, pload

      nmax_step = nsteps  ! number of pre-iteration steps
      ninj_step = 3000

      nstage_part = 3
      if (abs(time_integ) .eq. 2) nstage_part = 1

      ! pre simulation iteration for packed bed
      do istep=0,nmax_step

c        if (istep .eq. 0) then
c        if (ipart_restartr .ne. 0) then
c            do i=1,n
c            do j=0,2
c               rpart(jv0+j,i) = 0.
c               rpart(jv1+j,i) = 0.
c               rpart(jv2+j,i) = 0.
c               rpart(jv3+j,i) = 0.
c            enddo
c            enddo
c        endif
c        endif

         if (istep.eq.0) then
            time = 0
            pttime(1) = 0.
         else
            pttime(1) = pttime(1) + dnekclock() - ptdum(1)
         endif

         if (nid.eq. 0) write(6,*) 'pre-sim_io time',istep,time,dt_cmt
         if(mod(istep,iostep).eq.0) then
            call usr_particles_io
         endif

         do stage=1,nstage_part

            if (stage .eq. 1) then
               rdt_part = abs(param(12))
               call set_dt_particles(rdt_part)
               dt_cmt = rdt_part
               dt     = dt_cmt
               time   = time + dt_cmt
               if (abs(time_integ).eq.1)
     >            call set_tstep_coef_part(rdt_part)

               ! Update coordinates if particle moves outside boundary
               ptdum(2) = dnekclock()
                  call update_particle_location  
               pttime(2) = pttime(2) + dnekclock() - ptdum(2)

c              resetFindpts = 0
c              ilbstep = param(78)
c              if ((mod(istep,ilbstep).eq.0)) then
c                 ! Load balance if applicable
c                 resetFindpts = 0
c                 call computeRatio
c                 call reinitialize
c                 !call printVerify
c              else
c                 ! Update where particle is stored at
                  ptdum(3) = dnekclock()
                     call move_particles_inproc
                  pttime(3) = pttime(3) + dnekclock() - ptdum(3)
c              endif

               if (two_way.gt.1) then
                  ! Create ghost/wall particles
                  ptdum(4) = dnekclock()
                     call create_extra_particles
                     call sort_local_particles
                  pttime(4) = pttime(4) + dnekclock() - ptdum(4)
                  
                  ! Send ghost particles
                  ptdum(5) = dnekclock()
                     call send_ghost_particles
                  pttime(5) = pttime(5) + dnekclock() - ptdum(5)
                  
                  ! Projection to Eulerian grid
                  ptdum(6) = dnekclock()
                     call spread_props_grid
                  pttime(6) = pttime(6) + dnekclock() - ptdum(6)
               endif
            endif

            ! Evaluate particle force models
            ptdum(8) = dnekclock()
               call usr_particles_forcing  
            pttime(8) = pttime(8) + dnekclock() - ptdum(8)
   
            ! Integrate in time
            ptdum(9) = dnekclock()
               if (abs(time_integ) .eq. 1) call rk3_integrate
               if (abs(time_integ) .eq. 2) call bdf_integrate
            pttime(9) = pttime(9) + dnekclock() - ptdum(9)
   
            ! Update forces
            ptdum(10) = dnekclock()
               call compute_forcing_post_part
            pttime(10) = pttime(10) + dnekclock() - ptdum(10)

         enddo

         ptdum(1) = dnekclock()
      enddo


      return
      end
c----------------------------------------------------------------------
      subroutine sort_local_particles
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'CMTPART'
      include 'CMTDATA'

      common /save_dt_part/ rdpe_max, rdpe_min

      rat  = rdpe_max/rleng
      ratm = 2.
      if (rat .gt. ratm) then

        nrat = floor(rat)

        do i=1,n
           ie = ipart(je0,i) + 1
           rdx  = (xerange(2,1,ie) - xerange(1,1,ie))/nrat
           rdy  = (xerange(2,2,ie) - xerange(1,2,ie))/nrat
           rdz  = (xerange(2,3,ie) - xerange(1,3,ie))/nrat

           rpx  = rpart(jx,i) - xerange(1,1,ie)
           rpy  = rpart(jy,i) - xerange(1,2,ie)
           rpz  = rpart(jz,i) - xerange(1,3,ie)

           ipart(jicx,i) = floor(rpx/rdx)
           ipart(jicy,i) = floor(rpy/rdy)
           ipart(jicz,i) = floor(rpz/rdz)
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine userf_particles(ix,iy,iz,e,ffxp,ffyp,ffzp,qvolp)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'CMTPART'
      include 'CMTDATA'

      integer ix,iy,iz,e

      real ffxp,ffyp,ffzp,qvolp

      ! particle forcing
      if (two_way .ge.2) then
         if (istep .gt. time_delay) then
            ffxp =  ptw(ix,iy,iz,e,1)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))
            ffyp =  ptw(ix,iy,iz,e,2)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))
            ffzp =  ptw(ix,iy,iz,e,3)/vtrans(ix,iy,iz,e,1)
     >                              /(1.-ptw(ix,iy,iz,e,4))
            ! energy coupling for cmt-nek
            if (icmtp .eq. 1) then
               qvolp= ptw(ix,iy,iz,e,5) + rhs_fluidp(ix,iy,iz,e,4)
            else
               qvolp=0.
            endif
         else
            ffxp = 0.0
            ffyp = 0.0
            ffzp = 0.0
         endif
      else
         ffxp = 0.0
         ffyp = 0.0
         ffzp = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bdf_integrate
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'

      common /myparts/ times(0:3),alpha(0:3),beta(0:3)
      
      real s,pmass

      call get_bdf_ext_coefs(beta,alpha,times)

      jx0 = jx
c     move data to previous positions
      do j=0,ndim-1
      do i=1,n
         rpart(ju3+j,i)=rpart(ju2+j,i)
         rpart(ju2+j,i)=rpart(ju1+j,i)
         rpart(ju1+j,i)=rpart(ju0+j,i)
         rpart(jv3+j,i)=rpart(jv2+j,i)
         rpart(jv2+j,i)=rpart(jv1+j,i)
         rpart(jv1+j,i)=rpart(jv0+j,i)
         rpart(jx3+j,i)=rpart(jx2+j,i)
         rpart(jx2+j,i)=rpart(jx1+j,i)
         rpart(jx1+j,i)=rpart(jx0+j,i)
      enddo
      enddo

c     Solve for velocity at time t^n
      do i=1,n
        do j=0,ndim-1
          rhs = rpart(jf0+j,i)
     $        +     beta (1)*rpart(jv1+j,i)
     $        +     beta (2)*rpart(jv2+j,i)
     $        +     beta (3)*rpart(jv3+j,i)
          rpart(jv0+j,i) = rhs / beta(0)
          rhx = beta (1)*rpart(jx1+j,i)
     $        + beta (2)*rpart(jx2+j,i)
     $        + beta (3)*rpart(jx3+j,i) + rpart(jv0+j,i)
          rpart(jx0+j,i) = rhx / beta(0)     ! Implicit solve for x
        enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      FUNCTION MixtPerf_C_GRT_part(G,R,T,icmt)
      REAL G,R,T! INTENT(IN) 
      REAL MixtPerf_C_GRT_part
      integer icmt
      if (icmt .eq. 0) then
         MixtPerf_C_GRT_part = 1.
      else
         MixtPerf_C_GRT_part = SQRT(G*R*T)
      endif

      END
c-----------------------------------------------------------------------
