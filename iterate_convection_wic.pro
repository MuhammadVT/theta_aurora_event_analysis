pro iterate_convection_wic, date, stime, etime, coords

	FOR i = 0, (etime-stime)/2 DO BEGIN
		tm = stime + i * 2
		convection_plot_with_wic_image, date, tm, coords
	ENDFOR

END
