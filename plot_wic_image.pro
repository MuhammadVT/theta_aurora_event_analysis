pro plot_wic_image, time

        ; This procedure plots IMAGE WIC image data
	coords = "mlt"
	ps_open, './plots/convection_with_wic_image/tmpfig.ps'

	;restore, "./from_harald_frey/IMFHWIC_2002_0318_160426i.sav"
	fname_wic = "./from_harald_frey/IMFHWIC_2002_0318_" + STRCOMPRESS(STRING(time), /REMOVE_ALL) + "*.sav"
	restore, fname_wic

	image = imageinfo.Mlt_img
	lats = imageinfo.Mlat
	lons = imageinfo.Mlon
	mlts = imageinfo.Mlt
	image_size = SIZE(image, /DIMENSIONS)
	height = image_size[0]
	width = image_size[1]

	no_plot = 0

	; scale image
	;cimage = bytscl(image, min=scale[0], max=scale[1], top=(ncolors - bottom - 1), $
	;       /nan) + bottom
	cimage = bytscl(image)
	;cimage = image
	;loadct, 1, BOTTOM=0
	cgLoadCT, 0, NColors=256
;	rad_load_colortable, 1, bw=bw, whitered=whitered, bluewhitered=bluewhitered, $
;	leicester=leicester, themis=themis, aj=aj, brewer=brewer, default=default, $
;	website=website



;	filename = FILEPATH(SUBDIR=['./'], 'worldelv.dat')
;	OPENR, lun, filename, /GET_LUN
;	READU, lun, image
;	FREE_LUN, lun

	thisPostion = [0.1, 0.1, 0.9, 0.9]
	cgIMAGE, image, POSITION=thisPosition , /KEEP_ASPECT_RATIO
	p = [0.02, 0.99, 0.98, 0.98]
	cgColorbar, Position=[p[0], p[1]-0.1, p[2], p[1]-0.05], MAXRANGE=1000, $
	/VERTICAL, /RIGHT
	

	ps_close, /no_filename

END
