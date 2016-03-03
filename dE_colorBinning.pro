PRO colorBinning, infits
infits = 'dE_cleanSample2_0.65_rebinnedSample2.fits'
s = mrdfits(infits,1,/SILENT)

gi_color = s.GMAG-s.IMAG

a = {binnum:-1, gi:0.0}
if max(gi_color) gt 30. then print, 'at least one spectra has a wonky color'
;bibins = REPLICATE(a, N_ELEMENTS(where(bi_color LE 30.)))
;bibins.bi = bi_color[where(gi_color LE 30.0)]
gibins = REPLICATE(a,N_ELEMENTS(gi_color))
gibins.gi = gi_color

gimin = MIN(gibins.gi)
gimax = MAX(gibins.gi)
range = 0.1

max = gimin
print, gi_color
print, "Fix number"
print, FIX((gimax-gimin/0.1 +0.5))

for i=0, 1 do begin
   print, 'binnum', i
   min = max
   max = min
   remaining = n_elements(gibins(where(gibins.binnum EQ -1)))

   if remaining gt 1 then begin
      print, 'remaining: ', remaining
      while n_elements(gibins(where(gibins.binnum eq i))) le 10 do begin
         max += 0.001
         gibins(where(gibins.gi ge min AND gibins.gi lt max)).binnum = i
         remaining = n_elements(gibins(where(gibins.binnum eq -1)))
         if remaining le 1 then begin
            gibins(where(gibins.binnum eq -1)).binnum = i
            print, 'remaining: ', remaining
            BREAK
         endif

         print, 'remaining: ', remaining
      endwhile
   endif
endfor


tags = ['binnum']
values = ['0']

strreplace,infits,'rebinnedSample2','colorBinned'
ADD_TAGS,s,tags,values,s2
s2.binnum = gibins.binnum
mwrfits,s2,infits,/CREATE
END
