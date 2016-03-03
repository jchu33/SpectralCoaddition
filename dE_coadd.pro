PRO dE_coadd, infits
; Executes the spectral co-addition with inverse variance weighting technique.
;read in the *colorBinning files
infits = 'dE_cleanSample_0.65_rebinnedSample.fits'
s = mrdfits(infits,1,/SILENT)
;creating a fits file to store these coadds in
coaddstar = {lbin:FLTARR(6153),spec:FLTARR(6153), ivar:FLTARR(6153),rms:FLTARR(4096),nstars:0d,bimean:0d, bistd:0d}
binnum = 0 
master = replicate(coaddstar,max(binnum)+1)

   print, 'binnum: 0'
   binstars = s
   nstars = n_elements(s)
   ;bicolor = median(binstars.F475MAG-binstars.F814MAG)
   ;bistd = stddev(binstars.F475MAG - binstars.F814MAG)

   lbin = binstars[0].lbin
   print, n_elements(lbin)

;At each lambda value, check if the ivar is a valid float. If it is,
;add it to the ivar sum array (coaddivar)

   print, 'starting ivar'
   coaddivar = fltarr(6153) ;CHANGE THIS NUMBER
   for j = 0, 4095 do begin
      a = fltarr(nstars)
      for k = 0, nstars - 1 do begin
         if binstars[k].ivar[j] gt -99999 and binstars[k].ivar[j] lt 99999 $
         then a[k] = binstars[k].ivar[j] $
         else a[k] = 0
      endfor

      coaddivar[j] = total(a,/PRESERVE_TYPE)
   endfor

;at each lambda value, if ivar != 0, calculate the sum of
;spec*ivar. If it is a value, divide it by the sum of the ivars

   print, 'starting spec'
   coaddspec = fltarr(6153) ;CHANGE THIS NUMBER
   for j = 0, 4095 do begin
      if coaddivar[j] eq 0 then begin
         coaddspec[j] = !VALUES.F_NAN
      endif else begin
         a = fltarr(nstars)
         for k = 0, nstars-1 do begin
            if binstars[k].spec[j] gt -99999 and binstars[k].spec[j] lt 99999 $
            then a[k] = binstars[k].spec[j]*binstars[k].ivar[j] $
            else a[k] = 0
         endfor
         
         b = a[where(a eq a)] ;silly trick for getting the number of finite elements in a
         coaddspec[j] = total(b)/coaddivar[j]
      endelse
   endfor

  ; """print,'starting rms'
  ; coaddrms = fltarr(13333) ;CHANGE THIS NUMBER
  ; for j = 0, 13332 do begin
  ;    if coaddivar[j] eq 0 then begin
  ;       coaddrms[j] = !VALUES.F_NAN
  ;    endif else begin
  ;       a = fltarr(nstars)
  ;       for k = 0, nstars-1 do begin
  ;          if binstars[k].spec[j] gt -99999 and binstars[k].spec[j] lt 99999 $
  ;          then a[k] = binstars[k].spec[j]*binstars[k].spec[j]*binstars[k].ivar[j] $
   ;         else a[k] = 0
   ;      endfor
;
 ;        b = a[where(a eq a)]
 ;        coaddrms[j] = sqrt(total(b)/coaddivar[j] - (coaddspec[j]*coaddspec[j]))/(n_elements(b)-1)
  ;    endelse
  ; endfor"""

   master.lbin = lbin
   master.ivar = coaddivar
   master.spec = coaddspec
  ; master.rms = coaddrms
   master.nstars = nstars
  ; master.bimean = bicolor
  ; master.bistd = bistd

strreplace,infits,'rebinnedSample','coadd'
mwrfits,master,infits,/CREATE
print, "IM DONEEEEEEEEEE"   
END
    
