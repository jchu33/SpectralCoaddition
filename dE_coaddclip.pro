PRO coaddition_sclip,infits
;read in the *colorBinning files
infits = 'dE_cleanSample_0.65_rebinnedSample.fits'
s = mrdfits(infits,1,/SILENT)
binnum = 0
nlambda = n_elements(s[0].lbin)
print, nlambda
coaddstar = {lbin:FLTARR(nlambda),spec:FLTARR(nlambda), ivar:FLTARR(nlambda),nstars:0d,bimean:0d, bistd:0d}
master = replicate(coaddstar,max(binnum)+1)


   print, 'binnum: 0'
   binstars = s
   nstars = n_elements(s)
   ;bicolor = median(binstars.F475MAG-binstars.F814MAG)
   ;bistd = stddev(binstars.F475MAG - binstars.F814MAG)

   lbin = binstars[0].lbin
   print, lbin
   medianspec = fltarr(nlambda)
   for j = 0, nlambda-1 do begin
      a = fltarr(nstars)
      for k = 0, nstars-1 do begin
         if binstars[k].spec[j] gt -999 and binstars[k].ivar[j] gt 0.01 $
            then a[k] = binstars[k].spec[j]
      endfor
      medianspec[j] = median(a)
   endfor

   coaddivar = fltarr(nlambda)
   for j = 0, nlambda-1 do begin
      a = fltarr(nstars)
      for k = 0, nstars-1 do begin
         if binstars[k].ivar[j] gt -9999 $ 
            then a[k] = binstars[k].ivar[j] else a[k] = 0
      endfor
      coaddivar[j] = total(a,/PRESERVE_TYPE)
   endfor

   clipivar = fltarr(nlambda)
   clipspec = fltarr(nlambda)
   for j = 0, nlambda-1 do begin
         
      if coaddivar[j] eq 0 then begin
         clipivar[j] = 0
         clipspec[j] = !VALUES.F_NAN
      endif else begin
         spec = fltarr(nstars)
         ivar = fltarr(nstars)
         n_sigma = fltarr(nstars)
         for k = 0, nstars-1 do begin
            spec[k] = binstars[k].spec[j]
            ivar[k] = binstars[k].ivar[j]
         endfor
         for k = 0, nstars-1 do begin
            if spec[k] ne spec[k] then n_sigma[k] = 100.0 $
            else n_sigma[k] = abs(spec[k] - medianspec[k])*sqrt(ivar[k])
         endfor
         	
         sig_ind = where(n_sigma lt 3.5)
         a_spec = spec[sig_ind]*ivar[sig_ind]
         a_ivar = ivar[sig_ind]

         clipivar[j] = total(a_ivar)
         clipspec[j] = total(a_spec)
      endelse
   endfor

   b = clipspec
   clipspec = fltarr(n_elements(b))
   for j = 0, n_elements(b)-1 do begin
      if clipivar[j] eq 0 then clipspec[j] = !VALUES.F_NAN $
      else clipspec[j] = b[j]/clipivar[j]
   endfor

   ;c = clipspec
   ;clipspec = fltarr(n_elements(c))
   ;for j = 0, n_elements(c) - 1 do begin
   ;   if size(c[j],/TYPE) eq 4 then clipspec[j] = j+(binnum*0.4)) $
   ;   else clipspec[j] = j
   ;endfor

   master.lbin = lbin
   master.ivar = clipivar
   master.spec = clipspec
;   master[i].rms = coaddrms
   master.nstars = nstars
  ; master[i].bimean = bicolor
  ; master[i].bistd = bistd


strreplace,infits,'rebinnedSample','coaddclip'
mwrfits,master,infits,/CREATE

END
