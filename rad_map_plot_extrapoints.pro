pro rad_map_plot_extrapoints, date, time, coords , extra_points=extra_points, fixed_points=fixed_points
; This procedure overlays fitted velocities at random points in 
; mlt coordinates (e.g., mlat and mltlon. 180 mltlon is mlt noon)
; on to the convection maps

    SETENV, 'RAD_GRD_DATA_PATH=/home/muhammad/Dropbox/mypapers/coauthor_papers/eriksson_2017/revision/plots/rad_vel_plot/data/%YEAR%/%GRD%/%HEMISPHERE%'
    SETENV, 'RAD_MAP_DATA_PATH=/home/muhammad/Dropbox/mypapers/coauthor_papers/eriksson_2017/revision/plots/rad_vel_plot/data/%YEAR%/%MAP%/%HEMISPHERE%'


; Note, you need to do <.r /davit/lib/vt/idl/misc/mlt.pro> before running thi procedure.
    ;date = 20030214
    ;time = [1932, 2300]
    int_hemi=0
    ;coords = "magn"
    ;coords = "mlt"

    ;function rad_map_plot_extrapoints, date, time, coords

    ;ff = FINDFILE("./custom_idl_codes")
    ;FOR k=0, n_elements(ff) DO BEGIN
    ;    .run ff[k]
    ;ENDFOR

    ;.r /davit/lib/vt/idl/misc/mlt.pro

    common rad_data_blk
    rad_map_read, date, time=time, /north
;    rad_map_read, date, time=time, /north, filemapex="./20060306.north.mapex.bz2"

    ; calculate index from date and time
    stime = time[0]
    etime = time[1]
    sfjul, date, stime, sjul
    dd = min( abs( (*rad_map_data[int_hemi[0]]).mjuls-sjul[0] ), sindex)
    ; check if time distance is not too big
    if dd*1440.d gt 2. then $
        prinfo, 'data found, but '+string(dd*1440.d,format='(I4)')+' minutes off chosen time.'

    sfjul, date, etime, ejul
    dd = min( abs( (*rad_map_data[int_hemi[0]]).mjuls-ejul[0] ), eindex)
    ; check if time distance is not too big
    if dd*1440.d gt 2. then $
        prinfo, 'data found, but '+string(dd*1440.d,format='(I4)')+' minutes off chosen time.'

    FOR indx=sindex, eindex DO BEGIN
        jul_tmp = (*rad_map_data[int_hemi[0]]).mjuls[indx]
        sfjul, date_tmp, time_tmp, jul_tmp, /jul_to_date
;        ps_open, './plots/rad_map_plot/' + STRTRIM(STRING(date_tmp), 1) + '/pot_plot_' + STRTRIM(STRING(date_tmp), 1) + $
;                 '.' + STRTRIM(STRING(time_tmp),1) +'.ps'

        ps_open, '/home/muhammad/Dropbox/mypapers/coauthor_papers/eriksson_2017/revision/plots/rad_map_plot/' + $
                 ;'fitted_vel_plot_limited_range_' + STRTRIM(STRING(date_tmp), 1) + $
                 'tmpfig_' + STRTRIM(STRING(date_tmp), 1) + $
                 '.' + STRTRIM(STRING(time_tmp),1) +'.ps'

;        ps_open, '/home/muhammad/Dropbox/mypapers/coauthor_papers/eriksson_2017/' + $
;                 'revision/plots/rad_map_plot/2min_nav/' + $
;                 'Fitted_vel_2min_nav_' + STRTRIM(STRING(date_tmp), 1) + $
;                 ;'LOS_2min_nav_' + STRTRIM(STRING(date_tmp), 1) + $
;                 '.' + STRTRIM(STRING(time_tmp),1) +'.ps'

;        rad_map_plot, date=date_tmp, time=time_tmp, coords=coords, scale=[0, 1000], $
;                      /coast, /no_fill, /north, /grad, /isotropic, $
;                      /no_plot_imf_ind, /website

        ;print, time_tmp
        ;time_tmp = [1832, 1850]
        ; plot the fitted velocity
        my_rad_map_plot, date=date_tmp, time=time_tmp, coords=coords, scale=[0, 1000], $
                       fixed_points=fixed_points, extra_points=extra_points, $
                       xrange = [-30, 30], yrange=[0, 30], $ 
                       ;/coast, /no_fill, /north, /grad, /isotropic, $
                       ;/no_plot_imf_ind, /website, $
                       ;/no_fill, /north, /grad, /isotropic, /contour_fill, /potential, $
                       /no_fill, /north, /grad, /isotropic, /vectors, $
                       ;/contour_fill, /potential, $
                       ;/no_fill, /north, /grad, /isotropic, /potential, $
                       /no_plot_imf_ind, /no_fov_names, /noOverlayHMB, /website

;        ; plot the zoomed in velocity
;        my_rad_map_plot, date=date_tmp, time=time_tmp, coords=coords, scale=[0, 1000], $
;                       fixed_points=fixed_points, extra_points=extra_points, $
;                       xrange = [-30, 30], yrange=[0, 30], $ 
;                       zoom_xr=[-15, 15], zoom_yr=[0, 20],$
;                       ;/coast, /no_fill, /north, /grad, /isotropic, $
;                       ;/no_plot_imf_ind, /website, $
;                       ;/no_fill, /north, /grad, /isotropic, /contour_fill, /potential, $
;                       /no_fill, /north, /los, /isotropic, /vectors, $
;                       /contour_fill, /potential, $
;                       /no_plot_imf_ind, /no_fov_names, /noOverlayHMB, /website , /zoom_map

;        ; plot the LOS velocity with rad fov names
;        my_rad_map_plot, date=date_tmp, time=time_tmp, coords=coords, scale=[0, 1000], $
;                       fixed_points=fixed_points, extra_points=extra_points, $
;                       ;xrange = [-35, 35], yrange=[-35, 35], $  ; for 1x1 plot
;                       ;xrange = [-15, 15], yrange=[0, 20], $ 
;                       xrange = [-10, 10], yrange=[0, 10], $ 
;                       ;/coast, /no_fill, /north, /grad, /isotropic, $
;                       ;/no_plot_imf_ind, /website, $
;                       ;/no_fill, /north, /los, /isotropic, /contour_fill, /potential, $
;                       /no_fill, /north, /los, /isotropic, /vectors, $
;                       ;/contour_fill, /potential, $
;                       /no_plot_imf_ind,  /noOverlayHMB, /website; , /zoom_map



;        rad_map_plot, time=${time}, xrange=[${xrange}], yrange=[${yrange}], $
;                 merge=('${model}' eq 'merge'), true=('${model}' eq 'true'), $
;                 los=('${model}' eq 'los'), grad=('${model}' eq 'grad'), $
;                 scale=[${scale}], coords='${coords}', hemisphere=${hemi}, $
;                 /coast, /no_fill, orig_fan=${orig_fan}, fan_radar_ids=[${st_ids}],/isotropic, $
;                 no_fov_names = ${no_fov_names}, no_plot_imf_ind = ${no_plot_imf_ind}, contour_fill = ${rem_cntr_fill}, $
;                 dms_ssies_overlay = ${dms_ssies_overlay}, dms_ssj_overlay = ${dms_ssj_overlay}, poes_flux_overlay = ${poes_flux_overlay}, $
;                 overlay_r1_oval = ${overlay_r1_oval}, overlay_r2_oval = ${overlay_r2_oval}, rotate=${rotate_map}, $
;                 zoom_map = ${zoom_map}, zoom_xr = [${xrange_zoom}], zoom_yr = [${yrange_zoom}], website=1

        ;rad_map_plot, date=date, time=time, coords=coords, /north, /grad
        ps_close, /no_filename
    ENDFOR

END


; e.g.
;        rad_map_plot_extrapoints, 20040325, [0846, 0848], "mlt" 


