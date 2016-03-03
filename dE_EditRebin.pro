PRO dE_EditRebin, infits, lint
;This is step 2 in creating the fits file that can be used to make coadds
;This code is a modification of newrebin.pro, it normalizes, converts
;velocity to the rest frame, and stitches blue and red


lambda_lo = 7500.0
lambda_hi = 7600.0

lmin = 4650D
lmax = 10000D
;lint = 0.65  because of 600 line grating with the wavelength range
lint =0.65
lbin_int = lint
nlbin = FIX((lmax-lmin)/lbin_int)


;read in fits file that contains red and 	blue values for a set of mastertra
master = mrdfits('dE_cleanSample.fits',1,/SILENT)

;add tags of the appropriate size
tags = ['lbin','spec','ivar']
value = 'fltarr('+strtrim(string(nlbin))+')'
values = [value,value,value]
ADD_TAGS,master,tags,values,masterFull

masterFull.lbin = lmin+lbin_int*DINDGEN(nlbin)


n = n_elements(masterFull)
unitflux = replicate(1,nlbin)

for i = 0, n_elements(masterFull)-1 do begin

;Adjust for doppler shift
   lrshift = masterFull[i].lbinr/(1.+masterFull[i].z)
   lbshift = masterFull[i].lbinb/(1.+masterFull[i].z)

;Rebin red spectra
   varR = 1.0/masterFull[i].ivarr
   x_specrebin,lrshift,masterFull[i].specr,masterFull[i].lbin, newRflux, var = varR, nwvar = newvarR, /SILENT
   x_specrebin, lrshift, unitflux[masterFull.specr], masterFull[i].lbin, unitbin, /SILENT
   
   specbinR = newRflux/unitbin
   newvarR[where(newvarR EQ 0.0)] = -10.0
   ivarbinR = 1.0/(newvarR*unitbin)
   ivarbinR[where(ivarbinr LE 0.0)] = 0.0

;Rebin blue spectra
   varB = 1.0/masterFull[i].ivarb
   x_specrebin,lbshift, masterFull[i].specb, masterFull[i].lbin, newBflux, var = varB, nwvar = newvarB,/SILENT
   x_specrebin, lbshift, unitflux[masterFull.specb], masterFull[i].lbin, unitbin,/SILENT

   specbinB = newBflux/unitbin
   newvarB[where(newvarB EQ 0.0)] = -10.0
   ivarbinB = 1.0/(newvarB*unitbin)
   ivarbinB[where(ivarbinb LE 0.0)] = 0.0

;combine red/blue spectra into a single array
   
   for j=0,nlbin-1 do begin
      if specbinR[j] GE -99999 AND specbinR[j] LE 99999 then begin
         masterFull[i].spec[j] = specbinR[j]
         masterFull[i].ivar[j] = ivarbinR[j]

      endif else begin
         masterFull[i].spec[j] = specbinB[j]
         masterFull[i].ivar[j] = ivarbinB[j]

      endelse
   endfor

   a = masterFull[i].spec
   region = a[where(masterFull[i].lbin ge lambda_lo and masterFull[i].lbin le lambda_hi)]
   med = median(region,/double)
   masterFull[i].spec = masterFull[i].spec/med
   masterFull[i].ivar = masterFull[i].ivar*med*med

  ; print, 'did loop '+string(i)+' out of '+string(n)
endfor

masterFull.ivar(where(masterFull.ivar le -9999 or masterFull.ivar ge 9999)) = 0.0

fitsparts = strsplit("dE_cleanSample.fits",'.',/EXTRACT)
outfits = fitsparts[0]+'_'+strtrim(string(lint,format = '(f4.2)'))+'_rebinnedSample.fits'
mwrfits,masterFull,outfits,/CREATE

END