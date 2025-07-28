C==========================================================================
      program rel_profile_wave
C==========================================================================
      

C--------------------------------------------------------------------------

      PARAMETER (MAXSTEP=400,MAXWAVE=10000)
      
      REAL PCNU(MAXWAVE),WAVE(MAXWAVE),FLUX(MAXWAVE)
      
      character*100 parfile,wavefile,outfile
      character*5 normalization,fluxunits
      character*1 wavescale,relativistic,version

C-- Call the subroutine PARAMETERS to get the model parameter and other
C--   specifications from a parameter file.

      parfile='parfile.dat'
      
      call parameters(
     &     PARFILE,MAXSTEP,
     &     OLAMBDA,WAVEMIN,WAVEMAX,WAVESCALE,WAVESTEP,WAVEFILE,
     &     NSTEP,RELATIVISTIC,NORMALIZATION,FLUXUNITS,OUTFILE,
     &     Q1,Q2,XIB,ANGI,XI1,XI2,BROAD,
     &     T0,ETA,ANGLAM,VERSION,
     &     AMP,NARMS,AOBS,PITCH,WIDTH,XISPIN,XISPOUT
     &     )      

C-- Call subroutine GETWAVESCALE to construct the wavelength scale in
C-- one of two ways. Either make a linear wavelength scale between the
C-- wavelength limits specified by the user or pull it out of the file
C-- specified by the user.
      
      call getwavescale
     &     (WAVEMIN,WAVEMAX,WAVESCALE,WAVESTEP,WAVEFILE,
     &     wave,npix)

C--   Call subroutione profile to compute the requested line profile on
C--   the wavelength scale supplied above.

      call profile(maxstep,                            ! grid array dimension
     &     xi1,xi2,broad,q1,q2,xib,angi,               ! disk parameters
     &     anglam,t0,eta,version,                      ! wind parameters
     &     amp,narms,aobs,pitch,width,xispin,xispout,  ! spiral parameters
     &     nstep,relativistic,olambda,                 ! computation settings
     &     npix,wave,pcnu)                             ! spectrum arrays

C-- Convert to the requested units (f-nu/f-lambda) and normalize the
C-- profile according to the scheme selected by the user
C-- (max/flux/none).

      call setflux
     &     (fluxunits,normalization,pcnu,
     &     npix,wave,flux)

C-- Write the input parameters and resulting model in an output file

      call writefile
     &     (parfile,outfile,npix,wave,flux,fluxunits)
      
C-- End of story

      write(*,*)
      STOP '>*< End of program rel_profile_wave!'
      END


      
C=============================================================================
      subroutine parameters(
     &     PARFILE,MAXSTEP,
     &     OLAMBDA,WAVEMIN,WAVEMAX,WAVESCALE,WAVESTEP,WAVEFILE,
     &     NSTEP,RELATIVISTIC,NORMALIZATION,FLUXUNITS,OUTFILE,
     &     Q1,Q2,XIB,ANGI,XI1,XI2,BROAD,
     &     T0,ETA,ANGLAM,VERSION,
     &     AMP,NARMS,AOBS,PITCH,WIDTH,XISPIN,XISPOUT
     &     )      
C=============================================================================
C Fetch the parameters needed to sonstruct a disk profile model from a
C parameter file. The parameter file contains quite a few parameters and
C other specifications, not all of which will be needed by any given
C family of models. But this subroutine is meant to be useful for many
C families of models hence it is written so as to fetch all parameters.
C
C A list of all the parameters is given below in the "data parlabel"
C statement. The parameters are identified by unique kewords, which are
C also read in. The program identifies the keywords and assigns the
C parameter values that it reads accordingly. As such, the parameters
C need not be listed in particular order in the parameter file as long
C as they are associated with the correct keywords.
C
C Each row of the parameter file begins with the parameter value and is
C followed by the relevant keword (unbroken string, in capitals) and
C then an explanation. These keyword values and names must befollowed by
C tabs or spaces. The explanation includes instructions on how to
C specify the parameter values (allowed range, constraints, syntax,
C etc). All but one of the parameters are numbers or single, short
C strings and they are read with unformatted read statements (which is
C why they must be followed by tabs or spaces). The only exception is
C WAVEFILE, which contains the path to a file so the string value must
C be enclosed in single quotes.
C
C After reading the parameter file, the program carries out a series of
C checks to make sure that all the necessary parameter have been
C supplied and that they have been specified correctly. Some potential
C problems are fatal, in which case the program issues an error message
C and terminates with extreme prejudice. Other problems can be "fixed"
C in which case the program issues a warning to describe what hapenned
C and continues. Finally some parameter can be specified in a way that
C triggers the program to reset them to other, default values (this is
C just for convenience). In such cases, the program issues a note to
C inform the user that the parameters have been reset.
C-----------------------------------------------------------------------------
      
      parameter (maxl=29)
      
      character*100 parfile,wavefile,outfile
      character*5 normalization,fluxunits
      character*1 wavescale,relativistic,version
      
      character line*100,value*100,keyword*20

      character*15 parlabel(maxl)

      logical debug
      debug=.false.
      
C--   Assign values to the parameters label array 

      data parlabel/
c 7      
     & 'Q1',            ! inner emissivity powerlaw index
     & 'Q2',            ! outer emissivity powerlaw index
     & 'XIB',	        ! power-law index break radius
     & 'ANGI',          ! disk inclination angle (degrees)
     & 'XI1',           ! inner disk radius (GM/c^2)
     & 'XI2',           ! outer disk radius (GM/c^2)
     & 'BROAD',         ! broadening parameter (km/s)
c 5
     & 'NSTEP',         ! integration steps (integer)
     & 'RELATIVISTIC',  ! include relativistic effects? (y/n)
     & 'NORMALIZATION', ! normalization scheme for profile (max/flux/none) [none]
     & 'FLUXUNITS',     ! flux density units (fnu/flam) [fnu]
     & 'OUTFILE',       ! output file where final model is written
c 6
     & 'OLAMBDA',       ! nominal wavelength of the line (Angstrom)
     & 'WAVEMIN',       ! minimum wavelength (Angstrom)
     & 'WAVEMAX',       ! maximum wavelength (angstrom)
     & 'WAVESCALE',     ! definition of the wavelength scale (f=file, s=step)
     & 'WAVESTEP',      ! size of wavelength step (Angstrom)
     & 'WAVEFILE',      ! name of file with wavelength scale (string)
c 4
     & 'T0',            ! optical depth normalization (0=no wind)
     & 'ETA',           ! optical depth power-law index
     & 'ANGLAM',        ! wind opening angle (degrees)
     & 'VERSION',       ! formula for escape probability (f=Flohic, m=Murray)
c 7
     & 'AMP',           ! contrast of spiral arm (0=no arms)
     & 'NARMS',         ! number of arms (integer)
     & 'AOBS',          ! orientation angle of spiral (deg, + inner, - outer)
     & 'PITCH',         ! pitch angle of spiral pattern (+ leading, - trailing)
     & 'WIDTH',         ! angular width of arm (degrees)
     & 'XISPIN',        ! inner spiral arm radius radius (GM/c^2, 0=XI1)
     & 'XISPOUT'        ! outer spiral arm radius (GM/c^2, 0=XI2)
c
     &     /

c      do i=1,maxl
c         write(*,*) i,parlabel(i)
c      enddo
      
C--   Open the parameter file and go through it. Read one line at a
C--   time, separate the keyord name from its value, and go through a
C--   msater IF-block that identifies the keywords and assigns values to
C--   interbal valiable accordingly. 
      
      open (unit=31,file=parfile,status='old')
      n=0
      m=0
      
 100  read(31,'(a)',end=200) line
      if (line(1:1).eq.'#') goto 100
      read(line,*) value,keyword
      n=n+1
c      if (debug) write(*,*) 
c      if (debug) write(*,*) n,value,keyword
      
      if (keyword.eq.'OLAMBDA') then ! ----- wavelength scale 
         read(value,*) olambda
         do i=1,maxl
            if (parlabel(i).eq.'OLAMBDA') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,olambda
      elseif (keyword.eq.'WAVEMIN') then
         read(value,*) wavemin
         do i=1,maxl
            if (parlabel(i).eq.'WAVEMIN') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,wavemin
      elseif (keyword.eq.'WAVEMAX') then
         read(value,*) wavemax
         do i=1,maxl
            if (parlabel(i).eq.'WAVEMAX') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,wavemax
      elseif (keyword.eq.'WAVESCALE') then
         read(value,*) wavescale
         do i=1,maxl
            if (parlabel(i).eq.'WAVESCALE') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,wavescale
      elseif (keyword.eq.'WAVESTEP') then
         read(value,*) wavestep
         do i=1,maxl
            if (parlabel(i).eq.'WAVESTEP') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,wavestep
      elseif (keyword.eq.'WAVEFILE') then
         wavefile=value
         do i=1,maxl
            if (parlabel(i).eq.'WAVEFILE') parlabel(i)='check'
         enddo         
         m=m+1
         if (debug) write(*,*) m,keyword,wavefile
      elseif (keyword.eq.'NSTEP') then ! ----- integration and normalization
         read(value,*) nstep
         do i=1,maxl
            if (parlabel(i).eq.'NSTEP') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,nstep
      elseif (keyword.eq.'RELATIVISTIC') then
         read(value,*) relativistic
         do i=1,maxl
            if (parlabel(i).eq.'RELATIVISTIC') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,relativistic
      elseif (keyword.eq.'NORMALIZATION') then
         read(value,*) normalization
         do i=1,maxl
            if (parlabel(i).eq.'NORMALIZATION') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,normalization
      elseif (keyword.eq.'FLUXUNITS') then
         read(value,*) fluxunits
         do i=1,maxl
            if (parlabel(i).eq.'FLUXUNITS') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,fluxunits
      elseif (keyword.eq.'OUTFILE') then
         outfile=value
         do i=1,maxl
            if (parlabel(i).eq.'OUTFILE') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,outfile
      elseif (keyword.eq.'Q1') then ! ----- circular disk
         read(value,*) q1
         do i=1,maxl
            if (parlabel(i).eq.'Q1') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,q1
      elseif (keyword.eq.'Q2') then
         read(value,*) q2
         do i=1,maxl
            if (parlabel(i).eq.'Q2') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,q2
      elseif (keyword.eq.'XIB') then
         read(value,*) xib
         do i=1,maxl
            if (parlabel(i).eq.'XIB') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,xib
      elseif (keyword.eq.'ANGI') then
         read(value,*) angi
         do i=1,maxl
            if (parlabel(i).eq.'ANGI') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,angi
      elseif (keyword.eq.'XI1') then
         read(value,*) xi1
         do i=1,maxl
            if (parlabel(i).eq.'XI1') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,xi1
      elseif (keyword.eq.'XI2') then
         read(value,*) xi2
         do i=1,maxl
            if (parlabel(i).eq.'XI2') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,xi2
      elseif (keyword.eq.'BROAD') then
         read(value,*) broad
         do i=1,maxl
            if (parlabel(i).eq.'BROAD') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,broad
      elseif (keyword.eq.'T0') then ! ----- wind properties
         read(value,*) t0
         do i=1,maxl
            if (parlabel(i).eq.'T0') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,t0
      elseif (keyword.eq.'ETA') then
         read(value,*) eta
         do i=1,maxl
            if (parlabel(i).eq.'ETA') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,eta
      elseif (keyword.eq.'ANGLAM') then
         read(value,*) anglam
         do i=1,maxl
            if (parlabel(i).eq.'ANGLAM') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,anglam
      elseif (keyword.eq.'VERSION') then
         read(value,*) version
         do i=1,maxl
            if (parlabel(i).eq.'VERSION') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,version
      elseif (keyword.eq.'AMP') then ! ----- spiral arms
         read(value,*) amp
         do i=1,maxl
            if (parlabel(i).eq.'AMP') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,amp
      elseif (keyword.eq.'NARMS') then
         read(value,*) narms
         do i=1,maxl
            if (parlabel(i).eq.'NARMS') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,narms
      elseif (keyword.eq.'AOBS') then
         read(value,*) aobs
         do i=1,maxl
            if (parlabel(i).eq.'AOBS') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,aobs
      elseif (keyword.eq.'PITCH') then
         read(value,*) pitch
         do i=1,maxl
            if (parlabel(i).eq.'PITCH') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,pitch
      elseif (keyword.eq.'WIDTH') then
         read(value,*) width
         do i=1,maxl
            if (parlabel(i).eq.'WIDTH') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,width
      elseif (keyword.eq.'XISPIN') then
         read(value,*) xispin
         do i=1,maxl
            if (parlabel(i).eq.'XISPIN') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,xispin
      elseif (keyword.eq.'XISPOUT') then
c         write (*,*) 'found offending keyword'
         read(value,*) xispout
c         write (*,*) 'read offending keyword value'
         do i=1,maxl
            if (parlabel(i).eq.'XISPOUT') parlabel(i)='check'
         enddo
         m=m+1
         if (debug) write(*,*) m,keyword,xispout
      endif
         
      goto 100
      
 200  close(31)
      
      if (debug) write(*,*)
      if (debug) write(*,*) 'parfile rows=',n
      if (debug) write(*,*) 'matched pars=',m

C--   Check that all the expected parameters were present in the
C--   parameter file. If there are missing parameters, issue an error
C--   message and stop.

      nmiss=0
      do i=1,maxl
         if (parlabel(i).ne.'check') then
            nmiss=nmiss+1
         endif
      enddo

      if (nmiss.eq.0) then
         write(*,500) maxl,parfile
 500     format(/,'All ',i2,' expected parameters ',
     &        'retrieved from ',a35)
      else
         write(*,510) nmiss,maxl
 510     format(/,'ERROR: The following ',i2,' of the expected ',i2,
     &        ' parameters are missing.')
         do i=1,maxl
            if (parlabel(i).ne.'check') then
               write(*,'(7x,a15)') parlabel(i)
            endif
         enddo
         write(*,*)
         stop '> FATAL ERROR: MISSING PARAMETERS, SEE LIST ABOVE'
      endif
      
C--   Check that the parameter values are sensible and issue error
C--   messages or warnings, if they are not. The first block below looks
C--   for fatal errors, i.e., parameters whose values will break the
C--   code and cannot be fixed. Any one of these errors causes the
C--   code to stop.
      
      nerr=0                    ! Tally of fatal errors
      
      if (wavemin.le.0.or.wavemax.le.0.or.wavemin.ge.wavemax) then
         write(*,*) 
         write(*,*) 'ERROR: Check WAVEMIN,WAVEMAX. '
         write(*,*)  '      Need WAVEMAX > WAVEMIN > 0.'
         nerr=nerr+1
      endif
      
      if (wavestep.eq.0.and.wavescale.eq.'s') then
         write(*,*) 
         write(*,*) 'ERROR: When WAVESCALE = s, ',
     &        'need to have WAVESTEP > 0.'
         nerr=nerr+1
      endif
      
      if (wavescale.eq.'s'.and.wavestep.le.0) then
         write(*,*) 
         write(*,*) 'ERROR: When WAVESCALE = s, ',
     &        'need WAVESTEP > 0.'
         nerr=nerr+1
      endif
      
      if ((wavescale.eq.'d'.or.wavescale.eq.'2'.or.wavescale.eq.'3')
     &     .and.wavefile.eq.'none') then
         write(*,*) 
         write(*,*) 'ERROR: When WAVESCALE = d/2/3, ',
     &        'need to specify WAVEFILE.'
         nerr=nerr+1
      endif
      
      if (broad.eq.0.or.angi.eq.0.
     &        or.xi1.eq.0.or.xi2.eq.0.or.xi1.gt.xi2) then
         write(*,*) 
         write(*,*) 'ERROR: Check circular disk parameters.'
         write(*,*) '       Need ANGI, BROAD > 0 and XI2 > XI1 > 0.'
         nerr=nerr+1
      endif
      
      if (xispin.gt.xi2.or.
     &     (xispin.ge.xispout.and.xispin.ne.0.and.xispout.ne.0)) then
         write(*,*) 
         write(*,*) 'ERROR: XISPIN not specified properly. '
         write(*,*) '       Need XI2 >/= XISPOUT > XISPIN >/= XI1'
         write(*,*) '       or XISPIN = 0.' 
         nerr=nerr+1
      endif
      
      if ((xispout.lt.xi1.and.xispout.ne.0).or.
     &     (xispin.ge.xispout.and.xispout.ne.0.)) then
         write(*,*) 
         write(*,*) 'ERROR: XISPOUT not specified properly. '
         write(*,*) '       Need XI2 >/= XISPOUT > XISPIN >/= XI1'
         write(*,*) '       or XISPOUT = 0.' 
         nerr=nerr+1
      endif

      if (nerr.gt.0) then
         write(*,*) 
         stop '> SEE FATAL PARAMETER ERROR(S) ABOVE'
      endif

C--   This block looks for parameters that were not specified correctly
C--   but they can be assigned default values nonetheless. Any such
C--   problem leads to a warning but the code continues. 
      
      if (nstep.gt.maxstep) then
         write(*,*) nstep
         write(*,520) NSTEP,MAXSTEP 
 520     format(/,' WARNING: NSTEP = ',I4,' which exceeds MAXSTEP.'
     &        'Reseting NSTEP = MAXSTEP = ',i4,'.')
      endif
      if (relativistic.ne.'y'.and.relativistic.ne.'n') then
         write(*,*) 
         write(*,*) 'WARNING: RELATIVISTIC not specified properly. ',
     &        'Defaulting to "y".'
         relativistic='y'
      endif
      
      if (normalization.ne.'max'.and.normalization.ne.'flux'.
     &        and.normalization.ne.'none') then
         write(*,*) 
         write(*,*) 'WARNING: NORMALIZATION not specified properly. ',
     &        'Defaulting to "none".'
         normalization='none'
      endif
      
      if (fluxunits.ne.'fnu'.and.fluxunits.ne.'flam') then
         write(*,*) 
         write(*,*) 'WARNING: FLUXUNITS not specified properly. ',
     &        'Defaulting to "fnu".'
         fluxunits='fnu'
      endif
      
      if (t0.eq.0) then
         write(*,*) 
         write(*,*) 'NOTE: T0 set to 0. ',
     &        'Reseting to T0=1e-6. No radiative transfer effects.'
         t0=1.e-6
      endif
      
      if (version.ne.'f'.and.fluxunits.ne.'m') then
         write(*,*) 
         write(*,*) 'WARNING: VERSION not specified properly. ',
     &        'Defaulting to axisymmetric escape probability.'
      endif
      
      if (xib.gt.xi2) then
         write(*,*) 
         write(*,*) 'NOTE: XIB > XI2. ',
     &        'Will default to q = Q1 throughout the disk.'
      endif
      
      if (xib.lt.xi1) then
         write(*,*) 
         write(*,*) 'NOTE: XIB < XI1. ',
     &        'Will default to q = Q2 throughout the disk.'
      endif
      
      if (xispin.eq.0.or.xispin.lt.xi1) then
         write(*,*) 
         write(*,*) 'NOTE: XISPIN = 0 or < XI1. ',
     &        'Reseting XISPIN = XI1.'
         xispin=xi1
      endif
      
      if (xispout.eq.0.or.xixpout.gt.xi2) then
         write(*,*) 
         write(*,*) 'NOTE: XISPOUT = 0 or > XI2. ',
     &        'Reseting XISPOUT = XI2.'
         xispout=xi2
      endif

      return
      end

      
C=============================================================================
      subroutine getwavescale(
     &     WAVEMIN,WAVEMAX,WAVESCALE,WAVESTEP,WAVEFILE, ! parameters (input)
     &     wave,npix                   ! wavelength array and pixels (output)
     &     )      
C=============================================================================
C This soubroutine retrieves the wavelength scale for the model, i.e.,
C the array of wavelengths at which the model flux density wiull be
C calculated. The parameter WAVESCALE (s/f) specifies which of the two
C available options to adopt. These options are
C
C (a) construct the wavelength scale given the minimum and maximum
C wavelengths and wavelength step, and
C
C (b) read the wavelength array from a file and keep only the
C wavelengths that fall between the minimum and maximum wavelength. The
C external file is assumed to be an ASCII file that is organized as
C follows
C - the wavelengths must be in the first column
C - any lines starting with one of the following comment characters are
C   skipped exclamation, semicolon, or pound, i.e. ! or ; or #)      
C - the wavelength scale file can have header lines that do not start
C   with a comment character, which have to be skipped      
C   if WAVESCALE=s, do not skip any lines (assume no header is present)
C   if WAVESCALE=2, skip the first 7 lines (assumed to be the 2cu header)
C   if WAVESCALE=3, skip the first 8 lines (assumed to be the 3cu header)
C----------------------------------------------------------------------------
      
      real wave(*)
      character wavescale*1,wavefile*100,line*100

      if (wavescale.eq.'s') then
         npix=1+int((wavemax-wavemin)/wavestep)
         do i=1,npix
            wave(i)=wavemin+wavestep*float(i-1)
         enddo
         write(*,500) wave(1),wave(npix),wavestep,npix,
     &        wavemin,wavemax
 500     format(
     &        /,'Wavelength scale constructed from input parameters',
     &        /,'Final range: ',f7.1,' - ',f7.1,' A, Step: ',f5.1,
     &        ' A, ',i5,' pixels'
     &        /,'[requested range: ',f7.1,' - ',f7.1,' A]')
      elseif (wavescale.eq.'d'.or.wavescale.eq.'2'.
     &        or.wavescale.eq.'3') then 
         open (unit=51,file=wavefile,status='old',err=890)
         if (wavescale.eq.'2') then
            do i=1,7
               read(51,'(a)') line
            enddo
         elseif (wavescale.eq.'3') then
            do i=1,8
               read(51,'(a)') line
            enddo
         endif
         npix=0
         nw=0
 100     read(51,'(a)',end=110) line
         if (line(1:1).ne.'!'.and.line(1:1).ne.';'.
     &        and.line(1:1).ne.'#') then
            read (line,*) w
            nw=nw+1
            if (nw.eq.1) wfirst=w
            wlast=w
            if (w.ge.wavemin.and.w.le.wavemax) then
               npix=npix+1
               wave(npix)=w
            endif
         endif
         goto 100
 110     close (51)
         if (wfirst.ge.wlast.or.wave(1).ge.wave(npix)) then
            write(*,550)
 550        format(/,'ERROR: Wavelengths in WAVEFILE not ',
     &           'in increasing order.')
            goto 900
         endif
         dw=(wave(npix)-wave(1))/float(npix)
         write(*,555) wavefile,
     &        wave(1),wave(npix),dw,npix,
     &        wfirst,wlast
 555     format(
     &        /,'Wavelength scale from file: ',a50,
     &        /,'Range: ',f7.1,'-',f7.1,' A, Avg step: ',f5.1,' A, ',
     &        i5,' pixels'
     &        /,'[full range in file: ',f7.1,'-',f7.1,']')
      else
         write(*,560) 
 560     format(/,'ERROR: Illegal WAVESCALE parameter in GETWAVESCALE')
         goto 900
      endif

C-- Below is the error message that will be issued if there is a problem
C-- opening the wavelength scale file. Then there is the STOP statement
C-- that the program will jump to if any of the checks above fail. Note
C-- the GOTO statement that precedes the block, which will redirect the
C-- flow if the program gets here smoothly (through its normal,
C-- error-free execution).
      
      goto 990
      
 890  write(*,570)
 570  format(/,'ERROR: Problem opening WAVEFILE')
      
 900  write(*,*)
      stop '> SEE FATAL ERROR(S) ABOVE'
      
 990  return
      end
      
C============================================================================
      subroutine profile(maxstep,                      ! array dimension
     &     xi1,xi2,broad,q1,q2,xib,angi,               ! disk parameters
     &     anglam,t0,eta,version,                      ! wind parameters
     &     amp,narms,aobs,pitch,width,xispin,xispout,  ! spiral parameters
     &     nstep,relativistic,olambda,                 ! computation settings
     &     npix,wave,flux)                             ! spectrum arrays
C----------------------------------------------------------------------------
C Calculates the profile of an emission lines from a circular disk with
C a pattern in the brightness distribution and considering radiative
C transfer through the base of an outflowing wind.
C============================================================================

c      parameter (maxstep=400)
      
      REAL*8 SINCOS(MAXSTEP),SINSIN(MAXSTEP)
      REAL*8 CNU,DOPPLER,BOOST,BEND,XX,ELEMENT

      real*8 sini,cosi,coti,tani,sinlam,cotlam
      real*8 omega,tau,direscprob,arg,radial

      REAL PHIP(MAXSTEP),WAVE(*),flux(*)
      character*1 version,relativistic

c      write(*,*) maxstep,                      ! array dimension
c     &     xi1,xi2,broad,q,angi,                       ! disk parameters
c     &     anglam,t0,eta,version,                      ! wind parameters
c     &     amp,narms,aobs,pitch,width,xispin,xispout,  ! spiral parameters
c     &     nstep,relativistic,olambda,                 ! computation settings
c     &     npix
      
C--   Set the values of useful constants

      cappa=1./4.7 ! ratio of Keplerian speed to terminal speed

      PI=3.14159
      clight=2.9979e5
      deg_per_rad=180./pi

C--   Convert the input parameters of the model into the proper units

C     Normalization constant for outer power law
      if (xib.le.x1.or.xib.gt.xi2) then
         pownorm=1.
      else
         pownorm=xib**(q2-q1)
      endif
      
C     From i(degrees) to sin i, cos i
      SINI = SIN(ANGI/deg_per_rad)
      COSI = COS(ANGI/deg_per_rad)

C     Broadening parameter from km/s to v/c
      BROAD=BROAD/clight   

C     From wind lambda(degrees) to lambda(radians), sin l, cos l, cot l
      ANGLAM=ANGLAM/deg_per_rad
      sinlam=sin(anglam)
      coslam=cos(anglam)
      cotlam=coslam/sinlam

C     Angles for spiral brightness patterm
      IF (AOBS.GE.0.) THEN
         XIREF=XI2
         AOBS=AOBS/deg_per_rad
      ELSE
         XIREF=XI1
         AOBS=-AOBS/deg_per_rad
      ENDIF
      TANPITCH=TAN(PITCH/deg_per_rad)
      SIG=WIDTH/SQRT(8.*LOG(2.))/deg_per_rad

C--   Construct the arrays of trigonometric functions. These functions
C--   depend only on the azimuth and will be very useful later in the
C--   CPU intensive loops.

      XIDEL = ALOG10(XI2/XI1)/NSTEP/2. ! log steps in radius
      XIDEL = 10.**XIDEL               ! 1/2 radial step size
      XIDIFF = (XIDEL-1./XIDEL)
      PHISTEP = 2.*PI/NSTEP       ! phi step size

      DO I=1,NSTEP
         PHIP(I)=0.5*PHISTEP*(2*I-1)
         SINCOS(I)=SINI*COS(PHIP(I))
         SINSIN(I)=SINI*SIN(PHIP(I))
      ENDDO
      
C--> This is the heart of the program. Three nested loops compute the line
C    profile as a function of wavelength, by integrating over the surface of 
C    the disk.

      DO K=1,NPIX
         XX=(OLAMBDA/WAVE(K))-1.
         CNU = 0.
         DO J=1,NSTEP
            XI=XI1*XIDEL**(2*J-1)
            XISTEP=XI*XIDIFF
            ALPHA=SQRT(1.-(3./XI))
c            BETA=SQRT(1.-(2./XI))
            vinf=(1./sqrt(xi))/cappa
            radial=t0*((xi/xi1)**eta)/(vinf/xi)
            if (xi.lt.xib) then
               XIPOW=XI**(1-Q1)
            else
               XIPOW=POWNORM*XI**(1-Q2)
            endif
            sinlam=sin(anglam*xi1/xi)
            coslam=cos(anglam*xi1/xi)
            cotlam=coslam/sinlam
            PSI0=AOBS+LOG10(XI/XIREF)/TANPITCH
            DO I=1,NSTEP
               ARMS=0.
               IF (AMP.NE.0.AND.XI.GE.XISPIN.AND.XI.LE.XISPOUT) THEN
                  DO N=1,NARMS
                     PSI=PSI0+2.*PI*FLOAT(N-1)/FLOAT(NARMS)
                     DPSI=ABS(PHIP(I)-PSI)
                     ARMS=ARMS+EXP(-DPSI*DPSI/2./SIG/SIG)
                     DPSI=2.*PI-ABS(PHIP(I)-PSI)
                     ARMS=ARMS+EXP(-DPSI*DPSI/2./SIG/SIG)
                  ENDDO
               ENDIF
               if (version.eq.'m') then !-> MC97 version
                  omega=sincos(i)*(sincos(i)+1.5*cappa*sinsin(i))
     &                 -cosi*(sincos(i)*cotlam+cosi
     &                        +0.5*cappa*sinsin(i)/sinlam)
               elseif (version.eq.'f') then !-> FEB12 version
                  omega=sincos(i)*(sincos(i)+1.5*cappa*sinsin(i))
     &                 +cosi*(sincos(i)/sinlam+cosi
     &                        +0.5*cappa*sinsin(i)/sinlam)
               else 
                  omega=1.
               endif
               tau=radial/abs(omega)
               direscprob=(1.-exp(-tau))/tau
               DOPPLER=1.+(SINSIN(I)/SQRT(XI))
               if (relativistic.eq.'y') then 
                  DOPPLER=ALPHA/DOPPLER ! relativistic
                  BOOST=DOPPLER*DOPPLER*DOPPLER
                  BEND=(1+((1-SINCOS(I))/(1+SINCOS(I)))/XI)
               else                  
                  DOPPLER=(1.-(SINSIN(I)/SQRT(XI)))/DOPPLER ! non-relativistic
                  DOPPLER=SQRT(DOPPLER) ! non-relativistic
                  BOOST=1.
                  BEND=1.
               endif
               EXPON=(1.+XX-DOPPLER)/DOPPLER/BROAD
               EXPON=EXPON*EXPON/2.
               ARG=
     &              BOOST*             ! Doppler boosting
     &              XIPOW*             ! radial brightness profile
     &              direscprob*        ! wind opacity
     &              (1.+0.5*AMP*ARMS)* ! spiral arm emissivity
     &              BEND*              ! light bending 
     &              EXP(-EXPON)        ! intrinsic line profile
               ELEMENT=ARG*XISTEP*PHISTEP
               CNU=CNU+ELEMENT
            ENDDO
         ENDDO
         flux(k)=cnu
c         write(*,*)  cnu
      ENDDO

      return
      end

      
C======================================================================
      subroutine setflux
     &     (fluxunits,normalization,pcnu,
     &     npix,wave,flux)
C======================================================================
C Put the flux density in the units requested by the user and normalize
C the line profile. The flux units are specified by the parameter
C FLUXUNITS (fnu/flam) and the normalization is specified by
C NORMALIZATION (max/flux/none). The array PCNU carries the model
C profile computed in the main program (in units of f-nu). The array
C FLUX carries the final profile after unit conversion and
C normalization.
C----------------------------------------------------------------------

      character*5 normalization,fluxunits
      real wave(*),flux(*),pcnu(*)

      clightarg=2.9979e-8       ! speed of light in 1e5 km/s

C-- Put the flux density in the units requested by the user. The PCNU
C-- array carries the computed profile in f-nu. The FLUX array will
C-- carry the final profile in the units requested by the user. The
C-- default is f-nu.

      if (fluxunits.eq.'flam') then
         do i=1,npix
            flux(i)=clightarg*pcnu(i)/(wave(i)*wave(i))
         enddo
      else
         do i=1,npix
            flux(i)=pcnu(i)
         enddo
      endif

C-- Normalize the line profile. If NORMALIZATION=max, find the maximum
C-- of the FLUX array (which is already in the units requested by the
C-- user) and divide FLUX array by that. If NORMALIZATION=flux, find the
C-- integrated flux of the PCNU array (which is always in f-nu) and
C-- divide the FLUX array by that. Otherwise do nothing so that the
C-- profile keeps its original normalization (the last option is the
C-- default).
      
      if (normalization.eq.'max') then
         fmax=flux(1)
         do i=1,npix
            fmax=max(fmax,flux(i))
         enddo
         do i=1,npix
            flux(i)=flux(i)/fmax
         enddo
      elseif (normalization.eq.'flux') then
         ftot=0.
         do i=1,npix
            if (i.eq.1) then
               dw=wave(2)-wave(1)
            elseif (i.eq.npix) then
               dw=wave(npix)-wave(npix-1)
            else
               dw=0.5*(wave(i+1)-wave(i-1))
            endif
            ftot=ftot+clightarg*dw*pcnu(i)/(wave(i)*wave(i))
         enddo
         do i=1,npix
            flux(i)=flux(i)/ftot
         enddo
      endif

      return
      end
      
C======================================================================
      subroutine writefile
     &     (parfile,outfile,npix,wave,flux,fluxunits)
C======================================================================
C Write out the final version of the model in a simple ASCII file. As a
C pre-amble, include the entire parameter file so that the user can have
C a complete record of how the model was specified. The lines of the
C pre-amble start with a '#' so that they are treated as comments when
C the file is read it later by another program. Then write out the model
C in two columns (wavelength, flux density). The FLUXUNITS parameter is
C also passed down so that the flude density column in the output file
C can be labelled appropriately.
C----------------------------------------------------------------------
      
      parameter (maxl=100)
      
      real wave(*),flux(*)
      integer lfinal(maxl)
      character*100 parfile,outfile,line(maxl),pline
      character ylabel*35,flabel*15,fluxunits*5
      
C-- Open the parameter file and read teh lines that it contains. Then go
C-- through each line and figure out the position of the last non-blank
C-- character.
      
      open (unit=21,file=parfile,status='old')
      nl=0
 100  read(21,'(a)',end=110) line(nl+1)
      nl=nl+1
      goto 100
 110  close(21)
c      write(*,*) nl,' lines from ',parfile

      do i=1,nl
         do k=100,2,-1
            if (line(i)(k:k).eq.' '.and.
     &           line(i)(k-1:k-1).ne.' ') then
               lfinal(i)=k-1
               goto 200
            endif
         enddo
 200     continue
      enddo

C-- Open the output file and write the lines of the parameters
C-- file. Start each line with a '#' which serves as a comment
C-- character. Obviously, some lines will end up starting with a
C-- '##'. Then, write PGPLOT labels and column headings (also as
C-- comments) the wavelengths and, finally, flux density arrays of the
C-- model profile.

      open (unit=41,file=outfile,status='unknown')

      do i=1,nl
         pline(1:2)='# '
         pline(3:lfinal(1)+2)=line(i)(1:lfinal(i))
         write(41,'(a)') pline(1:lfinal(1)+2)
      enddo

      if (fluxunits.eq.'flam') then
         flabel='f-lambda'
         ylabel='\(2156)\d\gl\u (normalized)'
      else
         flabel='f-nu'
         ylabel='\(2156)\d\gn\u (normalized)'
      endif
      
      write(41,500) ylabel,flabel
 500  format(
     &     '#',70('='),/,
     &     '#-x Wavelength (\A)',/,
     &     '#-y ',A33,/,
     &     '#',/,
     &     '# Wavelength (A)',5x,a11,/,
     &     '# --------------   --------------'
     &     )
      
      do i=1,npix
         write(41,550) wave(i),flux(i)
      enddo
 550  format (2xf10.2,6x,1pe15.8)
      
      close(41)

      return
      end
      
