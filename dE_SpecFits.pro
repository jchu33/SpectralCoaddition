PRO dE_SpecFits, infile
; This is step 1 in creating the fits file that can be used to make
; coadds. It reads in a list of spectra that getCleanSample.pro deemed
; uncrowded and unreddened, and creates one giant fits file that
; contains the relevant extensions to be fed into rebinSpecFits.pro

catalog = mrdfits("VDGC_catalog.fits",1,hdr,/SILENT)
ngvsphot = mrdfits("VDGC_ngvsphot.fits",1,hdr,/SILENT)
obj_array = catalog.OBJTYPE
dE_list = LIST()
numbersRemove = [10, 364, 366, 371, 453]
for i=0, n_elements(numbersRemove) - 1 do begin
    obj_array[numbersRemove[i]] = "GC"
endfor
slitfiles = ["vdgc1/", "vdgc2/", "vdgc3/", "vdgc4B/", "vdgc5/", "vdgc6/", "vdgc7/", "vdgc8/", "vdgc9B/"]
slit_assign = ["vdgc1", "vdgc2", "vdgc3", "vdgc4B", "vdgc5", "vdgc6", "vdgc7", "vdgc8", "vdgc9B"]
for i = 0, n_elements(catalog.OBJTYPE) - 1 do begin
	if catalog.OBJTYPE[i] eq "dENUC" then begin
		dE_list.Add, i
        endif
endfor
print, "DE LIST"	
print, dE_list
;creating the fits file that will get spit out
    specname = LIST()
    vdgc_list = LIST()
    umag = LIST()
    gmag = LIST()
    imag = LIST()
    final_index = LIST()
    spec_list = LIST()
for j=0, n_elements(slit_assign)-1 do begin
      spec_list = where(STRMID(catalog.SPEC1DNAME,7,5) eq slit_assign[j] or STRMID(catalog.SPEC1DNAME,7,6) eq slit_assign[j])
      spec_temp = LIST()
      for nbr = 0, n_elements(spec_list) -1 do begin 
          for nbr2 = 0, n_elements(dE_list) - 1 do begin
              if spec_list[nbr] eq dE_list[nbr2] then begin
                   spec_temp.Add, spec_list[nbr]
                 
              endif    
          endfor  
      endfor               
      for nbr = 0, n_elements(spec_temp)-1 do begin
          specname.Add, STRING(catalog.SPEC1DNAME[spec_temp[nbr]])
          vdgc_list.Add, STRING(slitfiles[j])
          ;umag.Add, ngvsphot.umag[spec_temp[nbr]]
          ;gmag.Add, ngvsphot.gmag[spec_temp[nbr]]
          ;imag.Add, ngvsphot.imag[spec_temp[nbr]]
          final_index.Add, spec_temp[nbr]
      endfor      
endfor 
print, dE_list
    stars = {z:0d,zquality:0.0,maskname:'',slitname:'',objname:'',lbinR:FLTARR(4096),specR:FLTARR(4096),ivarR:FLTARR(4096),lbinB:FLTARR(4096),specB:FLTARR(4096),ivarB:FLTARR(4096),UMAG:0d,GMAG:0d,IMAG:0d, JMAG:0d, HMAG:0d}
    
    goodstars = replicate(stars,n_elements(specname)) ; returns an array with the given dimensions, filled with the scalar value specified as the first parameter.
    
    ; fileparts has multiple indices where each index was separated by a period. E.g ;spec1d.vdgc9B.040.gcN171.fits.gz would have
    ; vdgc9B to be mask, 040 to be slit, and gcN171 to be id
    for i=0,n_elements(specname)-1 do begin
       
       fileparts = strsplit(specname[i],'.',/EXTRACT)
       mask = fileparts[1]
       slit = fileparts[2]
       id = fileparts[3]
    
       print, mask
       print, slit
       print, id
       ;objid = catalog.OBJTYPE 
       maskname = catalog.MASKNAME
       slitname = catalog.SLITNUM
       z = catalog.ZHEL
    
       goodstars[i].maskname = mask
       goodstars[i].slitname = slit
       goodstars[i].objname = id
       goodstars[i].z = catalog.ZHEL[final_index[i]]
    
    ;reading in specvl data from the spec1dfile -> change to VDGC catalog? name  
       specFile =specpath +vdgc_list[i] + specname[i]
       specFile2 = STRCOMPRESS(specFile, /REMOVE_All)
      ; specFile = vdgc_list[i] +  specname[i]
       test = file_test(STRCOMPRESS(specFile2, /REMOVE_All))
       print, test
       print, specFile2
       if test eq 0 then CONTINUE
           fits_info,specFile2,N_EXT = next,/SILENT
    
       print, "I am here"
       sB = mrdfits(specFile2,M1,hdr,/SILENT) ;what are sB and sR? sB is blue spectrum, sR is red spectrum
       sR = mrdfits(specFile2,2,hdr,/SILENT)
       ;print, sR.LAMBDA
       ;print, sB.LAMBDA
    
       if min(sR.LAMBDA) lt max(sB.LAMBDA) then sR = mrdfits(specFile2,3,hdr);,/SILENT)
       print, "HERE"
       goodstars[i].lbinr = sR.LAMBDA
       goodstars[i].specr = sR.SPEC
       goodstars[i].ivarr = sR.IVAR 
       goodstars[i].lbinb = sB.LAMBDA
       goodstars[i].specb = sB.SPEC
       goodstars[i].ivarb = sB.IVAR
    ; plugging in photometric data from the input file to the fits file
       ;goodstars[i].UMAG = umag[i]
       ;goodstars[i].GMAG = gmag[i]
       ;goodstars[i].IMAG = imag[i]
    
       print, 'finished '+specname[i]
    endfor
;print, goodstars
outfile = 'cleanSample_dE.fits'
mwrfits,goodstars,outfile,/CREATE
;makes fitsfile called outfile with highz stars for data

print, "Done!"  
END
