pro tmp 

	coords = "mlt"
	ps_open, './plots/convection_with_wic_image/tmpfig.ps'

	restore, "./from_harald_frey/IMFHWIC_2002_0318_160426i.sav"
	image = imageinfo.Mlt_img
	lats = imageinfo.Mlat
	lons = imageinfo.Mlon
	mlts = imageinfo.Mlt
	image_size = SIZE(image, /DIMENSIONS)
	height = image_size[0]
	width = image_size[1]

	no_plot = 0

	; get color preferences
	foreground  = get_foreground()
	ncolors     = get_ncolors()
	bottom      = get_bottom()
	;scale = [2e3, 3e4]
	;scale = [0., 1.]
	scale = [0., 1000.]

	; scale image
	;cimage = bytscl(image, min=scale[0], max=scale[1], top=(ncolors - bottom - 1), $
	;       /nan) + bottom
	cimage = bytscl(image)
	cgLoadCT, 0
;	rad_load_colortable, 1, bw=bw, whitered=whitered, bluewhitered=bluewhitered, $
;	leicester=leicester, themis=themis, aj=aj, brewer=brewer, default=default, $
;	website=website

	lat_lim = 30.
	lon_lim = -0.5
	cgplot, [-50, 50], [-50, 50], PSYM=1, Aspect=1.0
	for h=0, height-2 do begin
	   for w=0, width-2 do begin
;		   if ~finite(lats[w, h]) || ~finite(lons[w, h]) || $
;			   ~finite(lats[w+1, h]) || ~finite(lons[w+1, h]) || $
;			   ~finite(lats[w+1, h+1]) || ~finite(lons[w+1, h+1]) || $
;			   ~finite(lats[w, h+1]) || ~finite(lons[w, h+1]) then $
;			   continue
		   if lats[h, w] lt lat_lim then continue
		   if lons[h, w] lt lon_lim then continue
		   p1 = calc_stereo_coords(lats[w, h], lons[w, h], mlt=(coords eq 'mlt'))
		   p2 = calc_stereo_coords(lats[w+1, h], lons[w+1, h], mlt=(coords eq 'mlt'))
		   p3 = calc_stereo_coords(lats[w+1, h+1], lons[w+1, h+1], mlt=(coords eq 'mlt'))
		   p4 = calc_stereo_coords(lats[w, h+1], lons[w, h+1], mlt=(coords eq 'mlt'))
		   xx = [p1[0],p2[0],p3[0],p4[0]]
		   yy = [p1[1],p2[1],p3[1],p4[1]]
		   if n_elements(rotate) ne 0 then begin
			   _x1 = cos(rotate*!dtor)*xx - sin(rotate*!dtor)*yy
			   _y1 = sin(rotate*!dtor)*xx + cos(rotate*!dtor)*yy
			   xx = _x1
			   yy = _y1
		   endif
		   tt = where(no_plot eq cimage[w,h], cc)
		   if cc eq 0 then $
			   ;polyfill, xx, yy, color=cimage[w,h], noclip=0
			   cgColorFill, xx, yy, color=cimage[w,h], noclip=0
			   ;polyfill, xx, yy, color=get_gray(), noclip=0
	   endfor
	endfor

        p = [0.02, 0.99, 0.98, 0.98]
        cgColorbar, Position=[p[0], p[1]-0.1, p[2], p[1]-0.05], MAXRANGE=1000, $
        /VERTICAL, /RIGHT


	ps_close, /no_filename

END
