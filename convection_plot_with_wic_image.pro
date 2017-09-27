pro convection_plot_with_wic_image, date, time, coords
; This procedure overlays WIC image data in 
; mlt coordinates onto a convection map

    SETENV, 'RAD_GRD_DATA_PATH=./data/%YEAR%/%GRD%/%HEMISPHERE%'
    SETENV, 'RAD_MAP_DATA_PATH=./data/%YEAR%/%MAP%/%HEMISPHERE%'

    int_hemi=0

    common rad_data_blk
    rad_map_read, date, time=time, /north
;    rad_map_read, date, time=time, /north, filemapex="./20060306.north.mapex.bz2"

    ; calculate index from date and time
    stime = time
    sfjul, date, stime, sjul
    dd = min( abs( (*rad_map_data[int_hemi[0]]).mjuls-sjul[0] ), sindex)

    ; check if time distance is not too big
    if dd*1440.d gt 2. then $
        prinfo, 'data found, but '+string(dd*1440.d,format='(I4)')+' minutes off chosen time.'

	indx = sindex
    jul_tmp = (*rad_map_data[int_hemi[0]]).mjuls[indx]
    sfjul, date_tmp, time_tmp, jul_tmp, /jul_to_date
    ;        ps_open, './plots/rad_map_plot/' + STRTRIM(STRING(date_tmp), 1) + '/pot_plot_' + STRTRIM(STRING(date_tmp), 1) + $
    ;                 '.' + STRTRIM(STRING(time_tmp),1) +'.ps'
    
    ps_open, './plots/convection_with_wic_image/' + $
    		 'convection_with_wic_' + STRTRIM(STRING(date_tmp), 1) + $
    		 '.' + STRTRIM(STRING(time_tmp),1) +'.ps'
    
    
    rad_map_plot_for_wic, date=date_tmp, time=time_tmp, coords=coords, scale=[0, 1000], $
    	  xrange = [-30, 30], yrange=[-30, 30], $
	  /no_fill, /north, /grad, /vectors, /isotropic, $
	  /no_plot_imf_ind, /no_fov_names, /noOverlayHMB
;    my_rad_map_plot, date=date_tmp, time=time_tmp, coords=coords, scale=[0, 1000], $
;                   fixed_points=fixed_points, extra_points=extra_points, $
;                   xrange = [-30, 30], yrange=[0, 30], $ 
;                   ;/coast, /no_fill, /north, /grad, /isotropic, $
;                   ;/no_plot_imf_ind, /website, $
;                   ;/no_fill, /north, /grad, /isotropic, /contour_fill, /potential, $
;                   /no_fill, /north, /grad, /isotropic, /vectors, $
;                   ;/contour_fill, /potential, $
;                   ;/no_fill, /north, /grad, /isotropic, /potential, $
;                   /no_plot_imf_ind, /no_fov_names, /noOverlayHMB, /website



    ps_close, /no_filename

END


; e.g.
;	 convection_plot_with_wic_image, 20020318, 1604, "mlt"


