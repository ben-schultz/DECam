;---------------------------------------------------------------------------
;fovsources2.pro
;author: ben schultz
;
; Just like fovsources.pro, but takes input from a .csv file
;  instead of an IDL data structure. Use with outputs of
;  find_sources.pro. Indices have also been adjusted to deal with
;  RA/DEC columns in these files
;
;INPUTS:
; - img: filename of FITS image to examine
; - data: .csv output of find_sources.pro corresponding
;    to image in question
; - CircSize, WindSize = optional values for 
;      crosshair and display window widths in pixels
;      defaults are 9 and 250 px, respectively
;
;EXAMPLE CALL:
;  fovsources, 'testimg', 'testimg_data.csv', CircSize = 5, WindSize = 500
; 
;
;LOG:
; 3/9/16 - modified to make program compatible with cut-data files
;          (where source N /= row # in file). for testing, use the
;          files r0.fits.fz and r0_find.4.1.c1.csv in the /scripts
;          directory. (also changed default wsize = 800). Make sure
;          fpack and funpack are enabled in shell before doing so!
;
;
;------------------------------------------------------------------------------
function data_table, datafile 

READCOL, datafile $
         , n, x, y, d, a, extn, flux, sharp, round $
         , skipline=1, /quick, /silent
; datafiles have extra junk line after header
; but .c#. files don't!

z = transpose([ [n], [x], [y], [d], [a], [extn], [flux] $
                , [sharp], [round] ])
return,z
end

;--------------------------------------------------------
function source_array, z, n, w ;works
; inputs are: data array, n of target, window size 

;id target source
t      = z[*, where(z[0,*] eq n)] ;edit
t_x    = t[1]
t_y    = t[2]
t_extn = t[5]

;id same ext sources
all_extn     = z[5,*]
same_extn    = where(all_extn eq t_extn, n_same)
             ;n_same is # of same extension sources
z_sub1       = z[*, same_extn]

;id fov sources
z_x          = z_sub1[1,*]
z_y          = z_sub1[2,*]

z_xs         = abs(z_x - t_x)
z_ys         = abs(z_y - t_y)

z_xs[where(z_xs gt w, /null)] = -1.0
z_ys[where(z_ys gt w, /null)] = -1.0

con_coords   = [z_xs, z_ys]
sum_coords   = total(con_coords,1)

inds         = intarr(n_same)
for i = 0, n_same-1 do begin
   if (sum_coords[i] ge max(con_coords[*,i])) then $
      inds[i] = 1
endfor

good_inds     = where(inds eq 1)
z_sub2        = z_sub1[*, good_inds]
return, z_sub2

end
;--------------------------------------------------------
function tar_coords, in_vec ;works
; in_vec = [ wsize,
;           (x,y) of target in FITS coords
;           (xmax, ymax) of FITS image]

x1 = in_vec[1] - in_vec[0]
x2 = in_vec[1] + in_vec[0]
y1 = in_vec[2] - in_vec[0]
y2 = in_vec[2] + in_vec[0]

if (x1 lt 0) then x1 = 0
if (y1 lt 0) then y1 = 0
if (x2 gt in_vec[3]) then x2 = in_vec[3] - 1
if (y2 gt in_vec[4]) then y2 = in_vec[4] - 1

xp = in_vec[1] - x1
yp = in_vec[2] - y1
out_vec = [xp, yp, x1, y1, x2, y2] 

return, out_vec
end
;--------------------------------------------------------
function src_coords, in_vec
; in_vec = [(x,y) of subwindow origin in FITS coords,
;           (x,y) of source in FITS coords

xs = in_vec[2] - in_vec[0]
ys = in_vec[3] - in_vec[1]
out_vec = [xs, ys]

return, out_vec
end
;--------------------------------------------------------
pro print_data, v ;v is row of source_array() output

n_char = strn(long(v[0]))
print,"       "
print,"Relevant FIND values for source " + n_char + ' are:'
print,"EXTN:  " + strn(v[5])
print,"X:     " + strn(v[1])
print,"Y:     " + strn(v[2])
print,"RA:    " + strn(v[4])
print,"DEC:   " + strn(v[3])        
print,"FLUX:  " + strn(v[6])
print,"SHARP: " + strn(v[7])
print,"ROUND: " + strn(v[8])
print,"       "
end
;--------------------------------------------------------
pro main_body, img_name, z, n, CircSize=csize, Windsize=wsize

;keyword defaults
if N_Elements(csize) eq 0 then csize=9
if N_Elements(wsize) eq 0 then wsize=800 ; new default
w = ceil(0.5 * wsize) ;rounds up to make integer

sources    = source_array(z, n, w)
sources_s  = size(sources)
ind        = where(sources[0,*] eq n)
tar        = sources[*,ind]
extn       = tar[5]
;tar = [n, x, y, d, a, extn, flux, sharp, round] of target

;read image
;user must supply image file extension!
img       = READFITS(img_name, h, ext=extn, /silent)
img_s     = size(img)

c1        = [w, tar[1], tar[2], img_s[1], img_s[2]] ;old coords
c2        = tar_coords(c1)              ; new window coords

GC, -1
subimg = img[c2[2]:c2[4], c2[3]:c2[5]]
TVS, subimg, med=1000

for i=0, sources_s[2] - 1 do begin
   s        =  sources[*,i]
   in       =  [ c2[2], c2[3], s[1], s[2]]
   xy       =  src_coords(in)
   if (s[0] eq tar[0]) then color = 'red' else color = 'green'
   tvcircle, csize, xy[0], xy[1], color, Device='TVS'
   text = strn(long(ceil( s[0])))
   xtext = xy[0] + csize
   ytext = xy[1] - csize
   ; print source designations
   xyouts, xtext, ytext, text, charsize=1.0, /DEVICE

endfor

print_data, tar

end
;--------------------------------------------------------
;loop wrapper - called in shell
pro fovsources2, img, data, CircSize=csize, Windsize=wsize

z = data_table(data)
s = size(z)
l = s[2]
inds = z[0,*]
;forprint,inds[0:9] ;debug

while (1 eq 1) do begin
   print, "Enter the integer index of the source"
   read, n, prompt="(Negative to quit): "
   if (n lt 0) then break

   ; test if N among sources in file
   n_test = where(inds eq n)
   if (n_test eq -1) then begin
      print, "Source is not present in file (invalid index)."
      print, "     "
      continue
   endif

   main_body, img, z, n, CircSize = csize, WindSize = wsize
   ;main body function does three things, given table data
   ; & source number:
   ; - print info about target source to screen
   ; - show + center + highlight target source in graphics window
   ; - show + highlight other sources in FOV

endwhile
end
