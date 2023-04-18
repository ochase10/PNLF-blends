c  Program to fit an observed Planetary Nebula Luminosity Function to 
c         an arbitary function with blends
c  For an explanation of the key variables, see the README file
c

	PARAMETER	(NARRAY = 1000, NLF = 20000, NGAUSS = 70000)
	PARAMETER	(N_PN = 1000)
	PARAMETER       (R_V = 3.1, R_5007 = 3.47)
	PARAMETER	(SUN_BOL = 4.74, BOL_CORR = -0.85)

	common /phterr/	errmag(100), err(100), earr(4,100), nerr
	common /param/	amu0, amu1, amures, alpha0, alpha1, alphares,
     1                            gaussres, amagres, znsigma
	common /pnlf/	phimag(NLF), phi1(NLF), phi2(NLF), phi0(NLF), 
     1                             cum(NLF), jcum

	common /pn_data/ pnmag(N_PN), background(N_PN), 
     1                     blend_frac(N_PN), vdisp(N_PN)

	character*80	fname
	character*1	carriage, ans
	real*8		prob(NARRAY,NARRAY)
	real*8		phimag, phi1, phi2, phi0, cum
	real*8		phi(0:NLF), cmag(0:NLF), gauss(0:NGAUSS)
	real*8		dscratch(NARRAY*NARRAY)
	real*8		fnorm, fnorm_25, summation, tot_lum
	real*8		Big_T, Big_T0, Big_T1
	real*8		flimit, flimit_abs, prob_pn, dfrac, vd
	real*8		poly_dnterp, poisson_prob, weight, weight_min
	real*8		alambda, alambda_fact, p1, p2, ptot, pmin
	real*8		alpha_scale(NARRAY), alpha_val(NARRAY)
	real*8		alpha_prob(NARRAY), alpha_temp(NARRAY)

	real*4		pscratch(NARRAY)
	real*4		temp_store(NARRAY,NARRAY)

	logical		blend, verbose

c  Smallest number for Poissonian probability
	weight_min = 1.0D-30

c  pmin is stand-in for log(0.)
	pmin = -1.0D38

c  The command line can either be blank or have an argument
c  If the argument is 'V' or 'v', then verbose mode is turned on
c  Otherwise, verbose is suppressed

	m = iargc()
	verbose = .true.
	if (m .gt. 0) then
	    call getarg (1, ans)
	    call up_case (ans)
	    if (ans .eq. 'Q') verbose = .false.
	end if

	carriage = achar(13)

c  Initial the 2-D likelihood array

	do i = 1, NARRAY
	    do j = 1, NARRAY
		prob(i,j) = 0.
	    end do
	end do

c  Subroutine entre reads all the program's parameters and data

	call entre (nphi, kpn, flimit, zlimit, tot_vmag, ebv, R_V, 
     1              fnorm_25, fname, blend, verbose)

c  This section of code prepares for the probability calculation

c  Create the convoluted PNLF, phi2  (if blending is turned on)

	if (blend) call blender (nphi)
	    
	ncum = nphi + jcum - 1

c  The observed luminosity function will be the true luminosity function
c     convolved with the photometric errors of the observations
c  Make sure the arrays are large enough for this convolution, 
c     compute the Gaussian kernal for convolution, and
c     fit a spline to photometric error versus magnitude data points
c  This is all skipped if no error versus magnitude file is read

	if (nerr .ne. 0) then
	    ngaussize = znsigma / gaussres
	    if (ngaussize .gt. NGAUSS) then
		write (6,302) ngaussize, NGAUSS
302		format ('Error:  N-sigma/convolution resolution (', i6,
     1          ') larger than dimensioned array NGAUSS (', i6, ')' )
		stop
	    end if

	    call gaussian (gauss, ngaussize, gaussres)
	    call spline (errmag, err, nerr, 2, earr)

	end if

c  "kloop" will contain the number of distance moduli analyzed
c  "mloop" will contain the number "alpha" values analyzed

	kloop = nint((amu1 - amu0) / amures)
	mloop = nint((alpha1 - alpha0) / alphares)

c  Make sure the array sizes are big enough for "kloop" and "mloop"

	if (kloop .gt. NARRAY) then
	    write (6,303) kloop, ' distance ', NARRAY
303	    format ('Error: ', i5, a10, 'steps requested, but ', 
     1                     'NARRAY dimension is only', i5)
	    stop
	end if

	if (mloop .gt. NARRAY) then
	    write (6,303) mloop, ' alpha ', NARRAY
	    stop
	end if

c  Calculate the alpha value for every step in the loop ahead of time

	do 4 m = 1, mloop
	alpha_val(m) = alpha0 + (float(m-1) * alphares)
4	continue

c
c  Now the real work begins: loop over all requested distances and alphas
c  Distance loop variable is i; alpha loop variable is m

	do 11 i = 1, kloop

	dmod = amu0 + (float(i) * amures)
	dmod_5007 = dmod + (R_5007 * ebv)
c
c  Check to see if the brightest PN is brighter than possible at the distance
c
	brightest_M = pnmag(1) - dmod_5007
	if (nerr .eq. 0) then
	    zmag = phimag(1)
	else
	    value = phimag(1) + dmod_5007
	    zmag = phimag(1)-
     1            (znsigma*splterp(errmag,err,earr,nerr,value))
	end if
 
c  Set (log) likelihood equal to zero if brightest PN is brighter than the
c     luminosity function allows
 
	if (brightest_M .lt. zmag) then
	    do j = 1, mloop
		prob(i,j) = pmin
	    end do
	    go to 11
	end if

	dist_10mpc = 10.**((dmod - 30.) / 5.)
	flimit_abs = flimit - dmod_5007
	fnorm = poly_dnterp (phimag(jcum),cum(jcum),ncum,flimit_abs,2)

	fnorm = fnorm / fnorm_25
	tot_lum = (SUN_BOL - (tot_vmag + BOL_CORR - dmod_5007)) / 2.5
	tot_lum = 10.**tot_lum
	Big_T0 = tot_lum * alpha0 * fnorm
	Big_T1 = tot_lum * alpha1 * fnorm

	write (6,300,advance='no') carriage, i, dmod
300     format (a1, i5, f10.3, i10, f10.2, f10.2)

c  Loop over all alpha values

	do 9 m = 1, mloop

	Big_T =  tot_lum * alpha_val(m) * fnorm

	weight = poisson_prob (kpn, Big_T)

	if (weight .lt. weight_min) then
	    alpha_prob(m) = pmin
	    go to 9
	end if

c  Calculation for lambda without the term that depends on surface brightness
c
	alambda_fact = 0.017146 * (alpha_val(m)/20D-9) * 
     1                                   fnorm * (dist_10mpc**2)

c  Total PN population (luminosity times alpha) brighter than flim

	summation = 0.

c  Loop over each PN, since each PN has a different LF

	k = 0
	alpha_scale(m) = 0.
	do while (k .lt. kpn)
	    k = k + 1

c  compute with weights for the 1-object and 2-object LFs, and create the LF

	    alambda = alambda_fact * background(k) * blend_frac(k)

		if (k .eq. 1) temp_store(i,m) = alambda
	    p1 = 1.
	    p2 = 0.
	    if (blend) then
	        p1 = alambda * exp(-alambda)
		p2 = alambda * p1 / 2.
	        ptot = p1 + p2
	        p1 = p1 / ptot
	        p2 = p2 / ptot
	    end if
	    alpha_scale(m) = alpha_scale(m) + p1 + (2. * p2)

	    do j = 1, nphi
		phi0(j) = p1 * phi1(j)
		if (blend) phi0(j) = phi0(j) + (p2 * phi2(j))
	    end do

c  Convolve input luminosity function with photometric errors
c
	    call obsfunc (zmag, dmod_5007, flimit, zlimit, nphi,
     1           gauss, ngaussize, cmag, phi, maxpnt, tphi)

	    vd = pnmag(k) - dmod_5007
	    prob_pn = poly_dnterp (cmag, phi, maxpnt, vd, 1)

	    vd = dlog(prob_pn/tphi)
	    summation = summation + vd
	end do
	alpha_prob(m) = summation + dlog10(weight)
	alpha_scale(m) = alpha_scale(m) / float(kpn)
	alpha_scale(m) = alpha_val(m) * alpha_scale(m)

9	continue

  	if (blend) call modify_alpha (alpha_scale, alpha_val, 
     1                     alpha_prob, alpha_temp, mloop, pmin)

	do 10 m = 1, mloop
	prob(i,m) = alpha_prob(m)
10	continue
	
11	continue

c
c  Write output info
c
	write (6,301) kpn
301	format (/'Number of PN analyzed = ', i8)
	call out2d (prob, NARRAY, kloop, mloop, amu0, amures, ebv,
     1       alpha0, alphares, dscratch, pscratch, fname, imax, jmax)

	write (6,*) 'max lambda = ', temp_store(imax,jmax)

	stop
	end
c
c
c
c  Subroutine to read in data and information
c
	subroutine entre (nphi, kpn, flimit, zlimit, tot_vmag, ebv, 
     1           R_V, fnorm_25, fname, blend, verbose)

c  Objects with separations less than SEEING_FRAC * seeing_FWHM (from the
c     PN data file) are considered to be unresolved
c  Objects with velocity separations less than DELTA_V are considered to be
c     unresolved
c  A_5007 = R_5007 * E(B-V)

	PARAMETER	(NLF = 20000)
	PARAMETER	(N_PN = 1000)
	PARAMETER	(SEEING_FRAC = 0.5, DELTA_V = 200.)
  
	common /phterr/	errmag(100), err(100), earr(4,100), nerr
	common /param/	amu0, amu1, amures, alpha0, alpha1, alphares,
     1                            gaussres, amagres, znsigma
	common /pnlf/	phimag(NLF), phi1(NLF), phi2(NLF), phi0(NLF), 
     1                         cum(NLF), jcum

	common /pn_data/ pnmag(N_PN), background(N_PN), 
     1                     blend_frac(N_PN), vdisp(N_PN)

	character*80	fname
	character*1	ans
	real*8		phimag, phi1, phi2, phi0, cum, fnorm_25
	real*8		vsum, flimit, x, blend_wave, blend_spatial
	logical		read_pndata, readxy, readxy_s, logic, blend
	logical		verbose

c
c  Read the empirical PNLF (absolute mag versus number of PN)
c
	if (.not. readxy ('True PNLF filename:  ',
     1              2, phimag, phi1, NLF, nphi, verbose)) stop


c  Read the PN data (5007 mag, background mag, seeing, velocity dispersion)

	if (.not. read_pndata ('Observed PN data filename:  ', 2, 
     1                             npn, verbose)) stop

c  Sorry the PN data by 

	call sort4 (pnmag, background, blend_frac, vdisp, npn)

c  Read the file containing photometric error versus magnitude

	logic = readxy_s ('Gaussian error vs. magnitude filename:  ',
     2            2, errmag, err, 100, nerr, verbose)

c  Enter general information

	if (verbose) write (6,100)
100	format ('Limiting magnitude for complete sample:  ', $)
	read (5,200) flimit
200	format (f20.0)
	if (flimit.lt.0.) stop

	if (verbose) write (6,101)
101	format ('Limiting magnitude for zero detections:  ', $)
	read (5,200) zlimit
	if (zlimit.lt.0.) stop

	if (verbose) write (6,102)
102	format ('E(B-V) in direction of galaxy:  ', $)
	read (5,200) ebv
	if (ebv.lt.0.) stop

	if (verbose) write (6,103)
103	format ('Minimum and Maximum distance moduli to compute:  ', $)
	read (5,*) amu0, amu1
	if (amu0 .lt. 0. .or. amu1 .lt. 0. .or. amu0 .eq. amu1) stop
	if (amu0 .gt. amu1) then
	    v1 = amu0
	    amu0 = amu1
	    amu1 = v1
	end if

	if (verbose) write (6,106)
106	format ('Distance modulus magnitude increment:  ', $)
	read (5,200) amures
	if (amures .le. 0.) stop

	if (verbose) write (6,104)
104	format ('Minimum and maximum PN/L (top 2.5 mag) to test:  ', $)
	read (5,*) alpha0, alpha1
	if (alpha0 .lt. 0. .or. alpha1 .lt. 0.) stop

	if (verbose) write (6,108)
108	format ('PN/L increment:  ', $)
	read (5,200) alphares
	if (alphares .le. 0.) stop

c  Calculation to scale the depth of the survey to alpha_2.5

	alphalim = phimag(1) + 2.5

	vsum = 0.
	k = 0
	do 1 i = 1, nphi
	vsum = vsum + phi1(i)
	cum(i) = vsum
	if (phimag(i) .le. alphalim) k = k + 1
1	continue
	vsum = (alphalim - phimag(k)) / (phimag(k+1) - phimag(k))
	fnorm_25 = cum(k) + (vsum * (phi1(k+1) - phi1(k)))
	jcum = 1

	if (verbose) write (6,105)
105	format ('Approximate magnitude resolution for convolution:  ', $)
	read (5,200) amagres
	if (amagres .le. 0.) stop

	if (verbose) write (6,107)
107	format ('Gaussian convolution resolution and N-sigma:  ', $)
	read (5,*) gaussres, znsigma

	if (verbose) write (6,109)
109	format ('Apparent V-magnitude of underlying population:  ', $)
	read (5,200) tot_vmag

	if (verbose) write (6,110)
110	format ('Line-of-sight velocity dispersion of galaxy:  ',$)
	read (5,200) sigma_gal

	if (sigma_gal .gt. 0.) then
	    x = DELTA_V / sigma_gal / sqrt(2.0d0)
	    blend_wave = erf(x)
	end if
c
c  Translate surface brightness of the galaxy (from PN data file) to
c    linear units and correct for fraction of blends which could be 
c    resolved via the PN's differing velocities.  Also count the
c    number of PN brighter than the limiting magnitude
c
	kpn = 0
	do i = 1, npn
	    background(i) = background(i) - (R_V * ebv) 
	    background(i) = 10.**((22.-background(i)) / 2.5)
	    blend_spatial = (SEEING_FRAC * blend_frac(i))**2
	    if (sigma_gal .lt. 0.) then
		x = DELTA_V / vdisp(i) / sqrt(2.0d0)
		blend_wave = erf(x)
	    end if
	    blend_frac(i) = blend_wave * blend_spatial
	    if (pnmag(i) .le. flimit) kpn = i
	end do
	if (kpn .eq. 0) stop

	if (verbose) write (6,111)
111	format ('Include blends in the analysis:  ',$)
	read (5,201) ans
	call up_case (ans)
	if (ans .eq. 'Y') then
	    blend = .true.
	else if (ans .eq. 'N') then
	    blend = .false.
	else
	    stop
	end if

	if (verbose) write (6,114)
114	format ('Filename for output:  ', $)
	read (5,201) fname
201	format (a)
	if (lnblnk(fname) .eq. 0) stop
c
	if (verbose) write (6,*) ' '
	return
	end
c
c
c
	subroutine get_columns (line, keyword, max_col)

	character*(*)	line
	character*40	word
	character*8	labels(8)
	integer*4	keyword(8), nchar(8)

	data labels / 'PN-MAG', 'SEEING', 'MAG-GAL', 'FLUX-GAL', 
     1                'LAMBDA', 'APER', 'V-DISP', 'KILL' /
	data nchar / 6, 6, 7, 8, 6, 4, 6, 4 /

	do 1 i = 1, 8
	keyword(i) = 0
1	continue
	max_col = 0

	n = lnblnk(line)
	call up_case (line(1:n))

	k = 1

	j = 0
2	j = j + 1
	call wordget (k, line, n, word, m)
	if (m .eq. 0) go to 4

	i = 0
3	i = i + 1
	if (i .eq. 9) go to 2

	if (word(1:m) .ne. labels(i)(1:nchar(i))) go to 3
	keyword(i) = j
	if (j .gt. max_col) max_col = j
	go to 2

4	if (keyword(1) .eq. 0) then
	    write (6,*) 'No keyword found for the PN 5007 magnitude'
	    stop
	else if (keyword(2) .eq. 0) then
	    write (6,*) 'No keyword found for the seeing'
	    stop
	else if (keyword(3) .eq. 0 .and. keyword(4) .eq. 0) then
	    write (6,*) 'No keyword for galaxy surface brightness'
	    stop 
	else if (keyword(4) .ne. 0 .and. (keyword(5) .eq. 0 .or.
     1           keyword(6) .eq. 0)) then
	    write (6,*) 'Aperture or wavelength missing for FLUX-GAL'
	    stop
	end if

	return
	end
c
c
c
c  Subroutine to convolve the input PNLF with photometric errors
c
	subroutine obsfunc (zmag, dmod, flimit, zlimit, nphi,
     1       gauss, ngaussize, cmag, phi, maxpnt, tphi)

	PARAMETER	(NLF = 20000)

	common /phterr/	errmag(100), err(100), earr(4,100), nerr
	common /param/	amu0, amu1, amures, alpha0, alpha1, alphares,
     1                            gaussres, amagres, znsigma
	common /pnlf/	phimag(NLF), phi1(NLF), phi2(NLF), phi0(NLF), 
     1                             cum(NLF), jcum

	real*8		phimag, phi0, phi1, phi2, cum
	real*8		dmag, delta, flimit, total, poly_dnterp
	real*8		cmag(0:NLF), phi(0:NLF)
	real*8		gauss(0:ngaussize), factor(0:NLF/2)
	logical		done

	nloop = 0

	factor(0) = 1.
	truelimit =  flimit - dmod
	truncation = zlimit - dmod
	maxbin = int((truelimit - zmag) / amagres) + 1
	if (mod(maxbin, 2) .ne. 0) maxbin = maxbin + 1
	delta = (truelimit - zmag) / float(maxbin)

	do j = 0, maxbin
	    phi(j) = 0.
	end do

	j = -1
	done = .false.

	do while (.not. done)
	    j = j + 1
	    dmag = zmag + (float(j) * delta)
	    if (j .le. maxbin) cmag(j) = dmag

	    if (dmag .ge. phimag(1)) then

	        phivalue = poly_dnterp (phimag, phi0, nphi, dmag, 2)

		if (nerr .eq. 0) then
		    nloop = 0
		    total = 1.0d0
		else
		    obmag = dmag + dmod
		    sigma = splterp (errmag, err, earr, nerr, obmag)
 		    sg = delta / (sigma * gaussres)
		    nloop = int(znsigma * sigma / delta)
		    if (nloop .gt. 10000) nloop = 9999

		    total = 0.0d0
		    do k = 1, nloop
			index = nint(float(k) * sg)
			if (index .gt. ngaussize) index = index - 1
			factor(k) = gauss(index)
			total = total + factor(k)
		    end do
		    total = (2. * total) + 1.
		end if

		do m = -nloop, nloop
		    n = j + m
		    if (n .ge. 0 .and. n .le. maxbin) phi(n) = 
     1                   phi(n) + (phivalue * factor(iabs(m)) / total)
		end do

	    end if

	    check_mag = zmag + (float(j-nloop) * delta)

	    done = dmag .gt. truncation
	end do

	tphi = -phi(0) - phi(maxbin)
	iweight = 2
	do j = 0, maxbin
	    tphi = tphi + (float(iweight) * phi(j))
	    iweight = 6 - iweight
	end do

	tphi = tphi * delta / 3.
	maxpnt = maxbin + 1

	return
	end
c
c
c
c  Subroutine to output the 2-D probabilities
c
c
	subroutine out2d (prob, NARRAY, kloop, mloop, amu0, amures, ebv,
     1       alpha0, alphares, dscratch, pscratch, fname, imax, jmax)
c
	PARAMETER (NLEVEL = 10)

	character*80	fname
	real*8		prob(NARRAY,NARRAY), dsum, confidence(NLEVEL)
	real*8		dscratch(NARRAY*NARRAY), spmax, pmax
	real*4		zlevel(NLEVEL), pscratch(NARRAY)
	data	confidence/0.38292d0, 0.68269d0, 0.86638d0, 0.95449d0,
     1       0.98758d0, 0.9973002d0, 0.99953474d0, 0.999936656d0,
     2       0.9999932043d0, 0.99999942657/

c
c  First, we're going to re-normalize the likelihoods;
c  To do this, first, find the maximum (log) likelihood
c
	pmax = prob(1,1)
	do i = 1, kloop
	    do j = 1, mloop
		if (prob(i,j) .gt. pmax) pmax = prob(i,j)
	    end do
	end do

	write (6,101) pmax
101	format (/'Maximum unnormalized 2-D log probability = ', 
     1             f20.4)

c
c  Now scale all the other (log) likelihoods so that the most
c     likely log value is 0
c  Then, take the antilog, to get rid of the log
c  And, while you're at it, find the sum of all likelihoods
c
	dsum = 0.0d0
	spmax = 0.0d0
	do i = 1, kloop
	    dscratch(i) = 0.
	    do j = 1, mloop
		prob(i,j) = prob(i,j) - pmax
		if (prob(i,j) .lt. -87.) then
		    prob(i,j) = 0.
		else
		    prob(i,j) = exp(prob(i,j))
		    dsum = dsum + prob(i,j)
		    dscratch(i) =  dscratch(i) + prob(i,j)
		end if
	    end do
	    if (dscratch(i) .gt. spmax) spmax = dscratch(i)
	end do
	spmax = dlog(spmax) + pmax
	write (6,102) spmax
102	format ('Maximum unnormalized 1-D log probability = ', 
     1             f20.4)

c
c  Note where the maximum likelihood is, and
c  Re-normalize all the likelihoods so that they equal 1.0
c  If the search was sufficiently all encompassing, this turns the
c  likelihoods into probabilities
c
	pmax = prob(1,1)
	do i = 1, kloop
	    do j = 1, mloop
		if (prob(i,j) .gt. pmax) then
		    pmax = prob(i,j)
		    dmod = amu0 + (float(i) * amures)
		    imax = i
		    jmax = j
		    alpha = (alpha0 + (float(j) * alphares))
		end if
		prob(i,j) = prob(i,j) / dsum
	    end do
	end do

	alpha = alpha * 1.0e9
	write (6,100) dmod, alpha
100	format (' Best fit (m - M)    =',f8.3,' (dereddened)'/
     1          ' Best fit alpha(2.5) =',f8.3,' x 1.0E-9'/)

c
c  Now, we're going to find the contours which contain
c  0.5, 1.0, 1.5, 2.0, 2.5, and 3.0 sigma probabilities
c
	k=0
	do i = 1, kloop
	    do j = 1, mloop
		k = k + 1
		dscratch(k)=prob(i,j)
	    end do
	end do

	call sort1 (dscratch, k)

	do j = 1, NLEVEL
	    zlevel(j) = 0.
	end do

	dsum = 0.0d0
	j = 1
	i = k + 1
	do while (j .le. NLEVEL .and. i .gt. 0)
	    i = i - 1
	    dsum = dsum + dscratch(i)
	    if (dsum .gt. confidence(j)) then
	      frac = (confidence(j) - (dsum-dscratch(i))) / dscratch(i)
		zlevel(j) = dscratch(i+1) +
     1                           (frac * (dscratch(i) - dscratch(i+1)))
		j = j + 1
	    end if
	end do
c
c  Write out all the probabilities in an unformatted file
c

	open (unit=4, file=fname, status='new', form='unformatted')
	write (4) amu0, amures, alpha0, alphares, ebv, mloop, 
     1                        kloop, NLEVEL, (zlevel(j),j=1,NLEVEL)
	do j = 1, mloop
	    do i = 1, kloop
		pscratch(i) = prob(i,j)
	    end do
	    write (4) (pscratch(i), i = 1, kloop)
	end do
	close (unit=4)
	return
	end
c
c
c
c  Subroutine to output 1-D probabilities into an ascii file
c
c
c  Computation of (unnormalized) Gaussian
c
	subroutine gaussian (gauss, ngaussize, gaussres)
c
	real*8		gauss(0:ngaussize), arg, gss
c
	gss = gaussres

	do i = 0, ngaussize
	    arg = dfloat(i) * gss
	    gauss(i) = exp(-0.5 * (arg**2))
	end do
	end
c
c  Subroutine to perform a spline interpolation
c
c
	function splterp (x, y, a, n, value)

	real*4 x(n), y(n), a(4,n)
c
	if (value.lt.x(1)) then
	    splterp=y(1)
	    return
	end if
c
	if (value.gt.x(n)) then
	    splterp=y(n)
	    return
	end if

	if (n.lt.4) then
	    splterp = poly_interp (x, y, n, value, 1)
	    return
	end if

	k=2
	do while (value.gt.x(k).and.k.lt.n) 
	    k=k+1
	end do
	t=value-x(k-1)

	splterp = 0.
	do l=1,4
	    splterp=(splterp*t)+a(l,k-1)
	end do
	return
	end

c
c  Subroutine to read (x,y) data from a file (in double precision)
c
	logical function readxy (Prompt, lun, x, y, max, nval, verbose)
	character*(*)	Prompt
	real*8		x(max), y(max)
	integer*4	lun, max, nval
	logical		openrfile, verbose

	nval = 0
	readxy = openrfile (Prompt, lun, verbose)
	if (.not. readxy) return
	if (max .gt. 0) then
	    imax = max
	else
	    imax = 32767
	end if
	i = 0
	do while (i .lt. imax)
	    read (lun, *, end=1) xp, yp
	    i = i + 1
	    x(i) = xp
	    y(i) = yp
	end do
1	nval = i
	close (unit=lun)
	return
	end
c
c
c
c
c  Subroutine to read (x,y) data from a file (in single precision)
c
	logical function readxy_s (Prompt, lun, x, y, max, nval, 
     1                        verbose)
	character*(*)	Prompt
	real*4		x(max), y(max)
	integer*4	lun, max, nval
	logical		openrfile, verbose

	nval = 0
	readxy_s = openrfile (Prompt, lun, verbose)
	if (.not. readxy_s) return
	if (max .gt. 0) then
	    imax = max
	else
	    imax = 32767
	end if
	i = 0
	do while (i .lt. imax)
	    read (lun, *, end=1) xp, yp
	    i = i + 1
	    x(i) = xp
	    y(i) = yp
	end do
1	nval = i
	close (unit=lun)
	return
	end
c
c
c
	subroutine sort1 (a,n)
c
c  Binary sort routine based on magic
c
	real*8 a(n), t1
	int=2
10	int=2*int
	if (int.lt.n) go to 10
	int=min0(n,(3*int)/4-1)
20	int=int/2
	ifin=n-int
	do 70 ii=1,ifin
	i=ii
	j=i+int
	if (a(i).le.a(j)) go to 70
	t1=a(j)
40	a(j)=a(i)
	j=i
	i=i-int
	if (i.le.0) go to 60
	if (a(i).gt.t1) go to 40
60	a(j)=t1
70	continue
	if (int.gt.1) go to 20
	return
	end
c
c
c
	subroutine sort4 (a,b,c,d,n)
c
c  Binary sort routine based on magic
c
	real*4 a(n),b(n),c(n),d(n)
	int=2
10	int=2*int
	if (int.lt.n) go to 10
	int=min0(n,(3*int)/4-1)
20	int=int/2
	ifin=n-int
	do 70 ii=1,ifin
	i=ii
	j=i+int
	if (a(i).le.a(j)) go to 70
	t1=a(j)
	t2=b(j)
	t3=c(j)
	t4=d(j)
40	a(j)=a(i)
	b(j)=b(i)
	c(j)=c(i)
	d(j)=d(i)
	j=i
	i=i-int
	if (i.le.0) go to 60
	if (a(i).gt.t1) go to 40
60	a(j)=t1
	b(j)=t2
	c(j)=t3
	d(j)=t4
70	continue
	if (int.gt.1) go to 20
	return
	end
c
c  Subroutine to compute spline coefficients
c
c
	SUBROUTINE SPLINE (X,Y,N,IXTRAP,A)
	REAL*4 X(N), Y(N), A(4,N)
	DX1 = X(2) - X(1)
	DY1 =(Y(2)-Y(1))/DX1*6.
	DO 10 I = 1, N-2
	DX2 = X(I+2) - X(I+1)
	DY2 =(Y(I+2)-Y(I+1))/DX2*6.
	A(1,I) = DX1
	A(2,I) = 2. * (DX1 + DX2)
	A(3,I) = DX2
	A(4,I) = DY2 -DY1
	DX1 = DX2
	DY1 = DY2
10	CONTINUE
	GO TO (21,22,23),IXTRAP
C  PARABOLIC ENDS
22	A(2,1) = A(2,1) + X(2) - X(1)
	A(2,N-2) = A(2,N-2) + X(N) - X(N-1)
C  CUBIC ENDS
	GO TO 21
23	DX1=X(2)-X(1)
	DX2=X(3)-X(2)
	A(2,1)=(DX1+DX2)*(DX1+2.*DX2)/DX2
	A(3,1)=(DX2*DX2-DX1*DX1)/DX2
	DXN2=X(N-1)-X(N-2)
	DXN1=X(N)-X(N-1)
	A(1,N-2)=(DXN2*DXN2-DXN1*DXN1)/DXN2
	A(2,N-2)=(DXN1+DXN2)*(DXN1+2.*DXN2)/DXN2
21	DO 11 I = 2, N-2
	A(2,I) = A(2,I) - A(1,I)/A(2,I-1) * A(3,I-1)
	A(4,I) = A(4,I) - A(1,I)/A(2,I-1) * A(4,I-1)
11	CONTINUE
	A(4,N-2) = A(4,N-2) / A(2,N-2)
	DO 12 I = 2, N-2
	J = N - I - 1
        A(4,J)=(A(4,J)-A(3,J)*A(4,J+1))/A(2,J)
12	CONTINUE
	DO 14 I = 1, N-2
	J = N - I
	A(4,J) = A(4,J-1)
14	CONTINUE
	GO TO (24,25,26),IXTRAP
C  LINEARY BOUNDARY CONDITION
24	A(4,1)=0.
	A(4,N)=0.
	GO TO 28
C  PARABOLIC BOUNDARY CONDITION
25	A(4,1) = A(4,2)
	A(4,N) = A(4,N-1)
	GO TO 28
C  CUBIC END CONDITION
26	A(4,1)=((DX1+DX2)*A(4,2)+DX1*A(4,3))/DX2
	A(4,N)=((DXN2+DXN1)*A(4,N-1)-DXN1*A(4,N))/DXN2
28	DO 15 J = 1, N-1
	TEM = X(J+1) - X(J)
	A(1,J) = (A(4,J+1) - A(4,J)) / (6.*TEM)
	A(2,J) = A(4,J) / 2.
	A(3,J) = ((Y(J+1) - Y(J)) / TEM) -
     *   (TEM*(((2.*A(4,J)) + A(4,J+1))) / 6.)
	A(4,J) = Y(J)
15	CONTINUE
	RETURN
	END
c
c  Open a file
c
	logical function openrfile (Prompt, lun, verbose)
	character*(*)	Prompt
	character*80	fname
	integer*4	lun
	logical		verbose

	openrfile = .false.
	if (verbose) write (6,100) Prompt
100	format (a, $)
	read (5,101) fname
101	format (a)
	if (lnblnk(fname) .eq. 0) return
	open (unit=lun, file=fname, status='old')
	openrfile = .true.
	return
	end
c
c
c
	logical function read_pndata (Prompt, lun, npn, verbose)

	PARAMETER	(N_PN = 1000)

	common /pn_data/ pnmag(N_PN), background(N_PN), 
     1                     blend_frac(N_PN), vdisp(N_PN)
	character*(*)	Prompt
	character*256	line
	character*20	word
	integer*4	keyword(8)
	logical		openrfile, verbose
c
c  arcsec/pixel of MUSE
c
	scale_muse = 0.2

c
c  Open the file
c
	npn = 0
	read_pndata = openrfile (Prompt, lun, verbose)
	if (.not. read_pndata) return
c
c  Ignore the comment line
c
	read (lun,100) line
100	format (a)

c
c  Read the keyword line and figure out what columns are interesting
c
	read (lun,100) line
	call get_columns (line, keyword, max_col)

c  Read the PN data
c
	i = 0
1	i = i + 1
2	read (lun,100,end=4) line
	n = lnblnk(line)
	k = 1
c
c  Read each word in the line, pick out the key numbers, and store them
c    in their arrays
c  Array vdisp is the temporary storage for the velocity dispersion
c
	do 3 j = 1, max_col
	vdisp(i) = 0.
	call wordget (k, line, n, word, m)

c  PN Magnitudes
	if (j .eq. keyword(1)) then
	    read (word(1:m),*) pnmag(i)
c  Seeing
	else if (j .eq. keyword(2)) then
	    read (word(1:m),*) blend_frac(i)

c  V Surface Brightness
	else if (j .eq. keyword(3)) then
	    read (word(1:m),*) background(i)

c  Gal flux in aperture
	else if (j .eq. keyword(4)) then
	    read (word(1:m),*) d3

c  Wavelength for surface brightness calculation
	else if (j .eq. keyword(5)) then
	    read (word(1:m),*) wave

c  Aperture size for surface brightness
	else if (j .eq. keyword(6)) then 
	    read (word(1:m),*) d1

c  Velocity dispersion of underlying stars
	else if (j .eq. keyword(7)) then
	    if (m .ne. 0) read (word(1:m),*) vdisp(i)
   
c  Ignore entry
	else if (j .eq. keyword(8)) then
	    if (word(1:1) .eq. 'x' .or. word(1:1) .eq. 'X') go to 2

	end if
3	continue

	if (keyword(3) .eq. 0.) then
c
c  If V surface brightness not given, translate flux within an aperture
c     to a V-mag surface brightness.  If the flux is zero or negative, a
c     magnitude of m=28 is arbitrarily assigned
c
	    if (d3 .gt. 0.) then
		area = 3.141592654 * ((d1 * scale_muse)**2)
		background(i) = d3 * (wave**2) / 3E18 / area
		background(i) = -2.5 * alog10(background(i)) - 48.6
	    else
		background(i) = 28.
	    end if

	end if

	if (i .lt. N_PN) go to 1
	i = i + 1

4	close (unit=lun)
	npn = i - 1

	return
	end
c
c
c
	subroutine blender (nphi)
c
c  Subroutine to convolve the PNLF with itself, thereby giving the
c     luminosity function of blended sources
c

	implicit real*8 (a-h, o-z)

	PARAMETER	(NLF = 20000)
	common /pnlf/	phimag(NLF), phi1(NLF), phi2(NLF), phi0(NLF), 
     1                             cum(NLF), jcum

	idmag = nint(phimag(2) * 1.E4) - nint(phimag(1) * 1.E4)
	dmag = dfloat(idmag) / 1.d4

	istart = 1 
	xstar = phimag(istart) - dmag

	icount = nint(0.75257 / dmag) 

	if (nphi + icount .gt. NLF) then
	    write (6,*) 'Dimension error on phi values'
	    stop
	end if

	do 2 i = nphi, istart, -1
	phimag(i+icount) = phimag(i)
	phi1(i+icount) = phi1(i)
	cum(i+icount) = cum(i)
2	continue

	do 3 i = 1, icount
	phimag(icount - i + 1) = phimag(istart) - (float(i) * dmag)
	phi1(icount - i + 1) = 0.0
	cum(icount - i + 1) = 0.0

3	continue
	jcum = icount + 1
	nphi = nphi + icount

	istart = icount + 1

	tsum = 0.
	do 4 i = 1, nphi
	tsum = tsum + phi1(i)
	phi2(i) = 0.
4	continue
	pstar = phimag(1) - (dmag / 2.)

	do 7 i = istart, nphi
	if (phi1(i) .eq. 0.) go to 7
	alum1 = 10.**((xstar - phimag(i)) / 2.5)

	do 6 j = istart, nphi
	alum2 = 10.**((xstar - phimag(j)) / 2.5)
	alum = alum1 + alum2
	tmag = xstar - (2.5 * dlog10(alum)) 

	ibin = int((tmag - pstar) / dmag) + 1

	if (ibin .le. 0 .or. ibin .gt. NLF) then
		write (6,*) 'Problem with ibin in blending', ibin
		stop
	end if
	phi2(ibin) = phi2(ibin) + (phi1(i) * (phi1(j) / tsum))
6	continue
7	continue

	return
	end
c
c
c
	real*8 function poly_dnterp (xp, yp, n, x, iflag)
c
c  Polynomial interpolation.  
c    flag = 0  means use the interpolation blindly
c    flag .ne. 0 means make sure behavior is monotonic between the points
c

	implicit real*8 (a-h, o-z)
	real*8  xp(n), yp(n)

	if (n .lt. 2) then
	    write (6,*) 'too few points in interpolation'
	    stop
	end if

	if (n .eq. 2) then
	    ndegree = 2
	    idx = 1
	    poly_dnterp = d_flagrange (xp(idx), yp(idx), ndegree, x)
	    return
	end if

	i = 1
1	i = i + 1
	if (i .eq. n) go to 2
	if (x .gt. xp(i)) go to 1
2	ilow = i - 1
	ihigh = i

	if (n .eq. 3) then
	    ndegree = 3
	    idx = 1
	else
	    ndegree = 4
	    idx = ilow - 1
	    if (idx .le. 0) idx = 1
	    if (ihigh .ge. n) idx = n - 3
	end if

	poly_dnterp = d_flagrange (xp(idx), yp(idx), ndegree, x)

	if (iflag .eq. 0) return

	if (yp(ilow) .eq. yp(ihigh)) then
	    poly_dnterp = yp(ilow)
	    return
	end if

	if (yp(ihigh) .gt. yp(ilow)) then
	     if (poly_dnterp .ge. yp(ilow) .and. 
     1           poly_dnterp .le. yp(ihigh)) return
	else
	     if (poly_dnterp .le. yp(ilow) .and. 
     1           poly_dnterp .ge. yp(ihigh)) return
	end if

	poly_dnterp = d_flagrange (xp(ilow), yp(ilow), 2, x)

        return
        end
c
c
c
	real*8 function d_flagrange (xp, yp, ndegree, x)

	implicit real*8 (a-h, o-z)
	real*8  xp(ndegree), yp(ndegree)

	y = 0.
	do 2 i = 1, ndegree
	top = 1.
	bottom = 1.
	do 1 j = 1, ndegree
	if (j .ne. i) then
	    top = top * (x - xp(j))
	    bottom = bottom * (xp(i) - xp(j))
	end if
1	continue
	y = y + (top * yp(i) / bottom)
2	continue
	d_flagrange = y

	return
	end
c
c
c
	real*4 function poly_interp (xp, yp, n, x, iflag)

        implicit real*4 (a-h, o-z)
        real*4  xp(n), yp(n)

        if (n .lt. 2) then
            write (6,*) 'too few points in interpolation'
            stop
        end if

        if (x .lt. xp(1) .or. x .gt. xp(n)) then
            write (6,*) 'x-coordinate out of range of interpolation'
            stop
        end if

	if (n .eq. 2) then
	    ndegree = 2
	    idx = 1
            poly_interp = flagrange (xp(idx), yp(idx), ndegree, x)
	    return
	end if

	i = 0
1	i = i + 1
	if (i .eq. n) go to 2
	if (x .gt. xp(i)) go to 1
2	ilow = i - 1
	ihigh = i

	if (n .eq. 3) then
            ndegree = 3
            idx = 1
        else
            ndegree = 4
	    idx = ilow - 1
	    if (idx .le. 0) idx = 1
	    if (ihigh .ge. n) idx = n - 3
	end if

        poly_interp = flagrange (xp(idx), yp(idx), ndegree, x)

	if (iflag .eq. 0) return

	if (yp(ilow) .eq. yp(ihigh)) then
	    poly_interp = yp(ilow)
	    return
	end if

	if (yp(ihigh) .gt. yp(ilow)) then
	    if (poly_interp .ge. yp(ilow) .and. 
     1           poly_interp .le. yp(ihigh)) return

        else
             if (poly_interp .le. yp(ilow) .and.
     1           poly_interp .ge. yp(ihigh)) return

        end if

	poly_interp = flagrange (xp(ilow), yp(ilow), 2, x)

	return
        end
c
c
c
        real*4 function flagrange (xp, yp, ndegree, x)

        implicit real*4 (a-h, o-z)
        real*4  xp(ndegree), yp(ndegree)

        y = 0.
        do 2 i = 1, ndegree
        top = 1.
        bottom = 1.
        do 1 j = 1, 4
        if (j .ne. i) then
            top = top * (x - xp(j))
            bottom = bottom * (xp(i) - xp(j))
        end if
1       continue
        y = y + (top * yp(i) / bottom)
2       continue
        flagrange = y

        return
        end
c
c
c
	real*8 function poisson_prob (n, alambda)
	implicit real*8  (a-h, o-z)

	x = float(n)
	xg = float(n+1)
	z = x * dlog(alambda) - alambda - dlgama(xg)
	if (z .gt. -150.) then
	    poisson_prob = exp(z)
	else
	    poisson_prob = 0.
	end if

	return
	end
c
c
c
	subroutine modify_alpha (scale, val, prob, temp, m, pmin)
	real*8	scale(m), val(m), prob(m), temp(m), poly_dnterp, pmin

	i0 = 0
1	i0 = i0 + 1
	if (i0 .gt. m) return
	if (prob(i0) .le. pmin) go to 1

	iend = i0 
2	iend = iend + 1
	if (iend .gt. m) go to 3
	if (prob(iend) .gt. pmin) go to 2

3	iend = iend - 1
	if (i0 .eq. iend) then
	    prob(i0) = pmin
	    return
	end if

	np = iend - i0 + 1

	do 4 i = 1, m

	if (val(i) .le. 0.) go to 4

	if (val(i) .lt. scale(i0)) then
	    temp(i) = pmin
	else if (val(i) .gt. scale(iend)) then
	    temp(i) = pmin
	else
	    temp(i) = poly_dnterp (scale(i0), prob(i0), np, val(i), 2)
	end if

4	continue

	do 5 i = 1, m
	prob(i) = temp(i)
5	continue

	return
	end
c
c
c
	subroutine wordget (k, line, n, outword, m)

c
c  Subroutine to get the next "word" in the line
c  Blanks for TABS denote the end of the word
c  Arguments "k" is the current position on the line
c            "n" is the last non-blank character in the line
c            "outword" is the output string
c            "m" is the number of character is "outword"
c

	character*(*) line, outword

1	if (k .gt. n) then
	    m = 0
	    return
	end if
c
c  read until the first non-blank character is reached
c
	if (line(k:k) .eq. ' ' .or. line(k:k) .eq. '\t') then
	    k = k + 1
	    go to 1
	end if
c
c  now read the word, until the next non-blank character is reached
c
	i = k
2	i = i + 1
	if (i .gt. n) go to 3
	if (line(i:i) .ne. ' ' .and. line(i:i) .ne. '\t') go to 2
c
c  copy the word, and position the pointer at the next non-blank character
c
3	i = i - 1
	outword = line(k:i)
	m = lnblnk(outword)
	k = i
4	k = k + 1
	if (k.le.n .and. (line(k:k).eq.' ' .or. line(k:k).eq.'\t')) 
     1                             go to 4
	return
	end
c
c
c
	subroutine up_case (string)

	character*(*)   string

	n = lnblnk(string)
	do i = 1, n
	    k = ichar(string(i:i))
	    if (k .ge. 97 .and. k .le. 122) then
		k = k - 32
		string(i:i) = char(k)
	    end if
	end do
	return
	end
