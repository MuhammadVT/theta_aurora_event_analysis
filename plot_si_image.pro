pro plot_si_image, time
        ; This procedure plots IMAGE SI image data

	coords = "mlt"
	ps_open, './plots/convection_with_si_image/tmpfig_si.ps'

	fname_si = "./from_harald_frey/IMF12LSI_2002_0318_" + STRCOMPRESS(STRING(time), /REMOVE_ALL) + "*.sav"
	restore, fname_si

	image = imageinfo.Mlt_img
	lats = imageinfo.Mlat
	lons = imageinfo.Mlon
	mlts = imageinfo.Mlt
	image_size = SIZE(image, /DIMENSIONS)
	height = image_size[0]
	width = image_size[1]

	; scale image
	cimage = bytscl(image)
	;loadct, 1, BOTTOM=0
	cgLoadCT, 0, NColors=256

	thisPostion = [0.1, 0.1, 0.9, 0.9]
	cgIMAGE, image, POSITION=thisPosition , /KEEP_ASPECT_RATIO
	p = [0.02, 0.99, 0.98, 0.98]
	cgColorbar, Position=[p[0], p[1]-0.1, p[2], p[1]-0.05], MAXRANGE=1000, $
	/VERTICAL, /RIGHT
	

	ps_close, /no_filename

END
