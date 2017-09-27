function rad_map_calc_fitvel_vecs, date, time, coords, npoints, fixed_points=fixed_points
; This ifunction is used to extract fitted velocities at
;random points in mlt coordinates (e.g., mlat and mltlon. 180 mltlon is mlt noon). 
; The points are specified by lats and lons variables. 

common rad_data_blk
rad_map_read, date, time=time, /north
int_hemi = 0    ; 0 for north, 1 for south

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

parse_date, date, year, month, day

;npoints = 1    ; number of (lat, lon) points where the fitted velocities are caculated
IF keyword_set(fixed_points) THEN BEGIN
    ; set fixed points where we do the velocity fitting
    lats = 0.5*FINDGEN(npoints) + 78.  
    ;lats = 0.5*FINDGEN(npoints) + 74.  
    ;lats = 0.5*FINDGEN(npoints) + 84.  
    lons = REPLICATE(195, npoints, 1)    ; mltlon value in degree. 180 degree is 12 mlt.
    extra_pos = transpose([[lats], [lons]])

    ; create column names
    lats_only_str = STRING(TRANSPOSE(extra_pos[0,*]), FORMAT="(F10.1)")
    lons_only_str = STRING(TRANSPOSE(extra_pos[1,*]), FORMAT="(F10.1)")
    lats_lons_str = lats_only_str
    FOR i=0, n_elements(lats_only_str)-1 DO BEGIN
        lats_lons_str[i] = "mlat"+STRTRIM(lats_only_str[i],1) + "_mltlon"+STRTRIM(lons_only_str[i],1)
    ENDFOR
ENDIF ELSE BEGIN
    ; create column names
    lats_lons_str = REPLICATE("mlat_mltlon", npoints)
    lats_only_str = REPLICATE("mlat", npoints)
    lons_only_str = REPLICATE("mltlon", npoints)
    FOR i=0, n_elements(lats_only_str)-1 DO BEGIN
        lats_only_str[i] = STRTRIM(lats_only_str[i],1) + "_" + STRTRIM(STRING(i),1)
        lons_only_str[i] = STRTRIM(lons_only_str[i],1) + "_" + STRTRIM(STRING(i),1)
    ENDFOR
ENDELSE

; open a file to write
line_width = 150
OPENW, U, "./data/"+STRTRIM(STRING(date),1)+"_fitvel_angle.txt", /GET_LUN, WIDTH=line_width
OPENW, UU, "./data/"+STRTRIM(STRING(date),1)+"_fitvel_mag.txt", /GET_LUN, WIDTH=line_width
OPENW, UUU, "./data/"+STRTRIM(STRING(date),1)+"_fitvel_mlatloc.txt", /GET_LUN, WIDTH=line_width
OPENW, UUUU, "./data/"+STRTRIM(STRING(date),1)+"_fitvel_mltlonloc.txt", /GET_LUN, WIDTH=line_width
;PRINTF, U, "jul", transpose(extra_pos[0,*]), FORMAT="(" + REPLICATE("F7.2, 1X, ", npoints-1, 1) +  "F7.2, 1X)"
PRINTF, U, "datetime", lats_lons_str
PRINTF, UU, "datetime", lats_lons_str
PRINTF, UUU, "datetime", lats_only_str
PRINTF, UUUU, "datetime", lons_only_str
;PRINTF, U, "datetime", string(transpose(extra_pos[0,*]), FORMAT="(F10.1)")
;PRINTF, UU, "datetime", string(transpose(extra_pos[0,*]), FORMAT="(F10.1)")


FOR index=sindex, eindex DO BEGIN


    ; get the order of the fit
    order = (*rad_map_data[int_hemi]).fit_order[index]

    ; get latitude shift
    lat_shft = (*rad_map_data[int_hemi]).lat_shft[index]

    ; get longitude shift
    lon_shft = (*rad_map_data[int_hemi]).lon_shft[index]

    ; get minimum latitude
    latmin = (*rad_map_data[int_hemi]).latmin[index]
    tmax      = (90.0 - latmin)*!dtor

    ; get coefficients of the expansion
    coeffs = (*(*rad_map_data[int_hemi]).coeffs[index])


    ; calculate lon_shft to convert magnetic longitude into mlt coordinates
    jul = (*rad_map_data[int_hemi[0]]).mjuls[index]
    utsec = (jul - julday(1, 1, year, 0, 0))*86400.d
    lon_shft_tmp = MLTConvertYrSec(year, utsec, 0.)*15.

    IF keyword_set(fixed_points) THEN BEGIN
        extra_pos = transpose([[lats], [lons]])
        ;if coords eq 'mlt' then $
            ;lon_shft += MLTConvertYrSec(year, utsec, 0.)*15.
            ;extra_pos[1,*] = extra_pos[1,*] - lon_shft
        extra_pos[1,*] = extra_pos[1,*] - lon_shft_tmp
        ;print, lon_shft_tmp
        ;print, [lon_shft, lat_shft]

    ENDIF ELSE BEGIN
        ; get "real" velocity data
        IF (N_ELEMENTS( (*(*rad_map_data[int_hemi]).gvecs[index]) ) EQ 0) THEN BEGIN
            prinfo, 'No vectors found.'
        ENDIF ELSE BEGIN
            ;lat_range = [60., 80.]
            lat_range = [82., 90.]
            ;lat_range = [82., 90.]
            ;lon_range = [6., 18.] * 15      ; in mltlon 
            lon_range = [10., 14.] * 15      ; in mltlon 
            ;lon_range = [12.5, 14.] * 15      ; in mltlon (18:16 - 18:18) 
            ;lon_range = [12.0, 13.2] * 15      ; in mltlon (18:19 - 18:20)
            ;lon_range = [13.5, 14.5] * 15      ; in mltlon (18:34 - 18:36)
            ;lon_range = [13.4, 14.5] * 15      ; in mltlon (18:36 - 18:38)
            vecs = (*(*rad_map_data[int_hemi]).gvecs[index])
            vecs_pos = transpose([[vecs[*].mlat], [vecs[*].mlon]])

            ; calculate fitted velocity vectors
            rvecs = rad_map_eval_grad_vecs(vecs_pos, coeffs[*,2], latmin, order)
            rmag = reform(sqrt(rvecs[0,*]^2 + rvecs[1,*]^2))

            ; convert magnetic longitude into mltlon
            vecs_pos[1,*] = vecs_pos[1,*] + lon_shft_tmp

            ; do modulo operation
            vecs_pos[1,*] = vecs_pos[1,*] + 360
            vecs_pos[1,*] = vecs_pos[1,*] MOD 360

            ; filter the vecs by positions in lat_range and lon_range in mlt coordinates
            selected_pos_rowindex = WHERE(LOGICAL_AND( $ 
                LOGICAL_AND(vecs_pos[0,*] GE lat_range[0], vecs_pos[0,*] LE lat_range[1]), $
                LOGICAL_AND(vecs_pos[1,*] GE lon_range[0], vecs_pos[1,*] LE lon_range[1])), count)
            IF count GT 0 THEN BEGIN
                selected_vecs_pos = vecs_pos[*, selected_pos_rowindex]
                selected_vecs_mag = rmag[selected_pos_rowindex]
                selected_nvecs = N_ELEMENTS(selected_vecs_pos[0,*])
                vels_mag_tmp = reform(selected_vecs_mag)
                sorted_vecs_rowindex = SORT(vels_mag_tmp) 
                IF selected_nvecs GT npoints THEN BEGIN
                    max_vecs_rowindex = sorted_vecs_rowindex[selected_nvecs-npoints:selected_nvecs-1]
                ENDIF ELSE BEGIN
                    max_vecs_rowindex = sorted_vecs_rowindex
                ENDELSE
                extra_pos = selected_vecs_pos[*, max_vecs_rowindex]
                ; back to magnetic longitude from mltlon 
                extra_pos[1,*] = extra_pos[1,*] - lon_shft_tmp
            ENDIF ELSE BEGIN
                extra_pos = 0 
            ENDELSE

            ;STOP 

        ENDELSE
    ENDELSE

    ; get "real" velocity data
    ;if n_elements( (*(*rad_map_data[int_hemi]).gvecs[index]) ) eq 0 then begin
    ;        prinfo, 'No vectors found.'
    ;endif
    ;real_vec = (*(*rad_map_data[int_hemi]).gvecs[index])
    ;real_nvec = (*rad_map_data[int_hemi]).vcnum[index]
    ;real_pos = transpose([[real_vec[*].mlat], [real_vec[*].mlon]])

;    ; First shift coordinates into 'model' reference (pole shifted 4 deg nightwards)
;    if lat_shft ne 0. then begin
;      npnts = N_ELEMENTS(extra_pos)/2L
;      kaz = FLTARR(npnts)
;      crd_shft, lon_shft, lat_shft, npnts, extra_pos, kaz
;    endif

    ; Calculate vectors
    IF (N_ELEMENTS(extra_pos) EQ 1) THEN BEGIN
        raz = REPLICATE("NaN", npoints)
        rmag = REPLICATE("NaN", npoints)
        mlat_pos = REPLICATE("NaN", npoints)
        mltlon_pos = REPLICATE("NaN", npoints)
        extra_pos = TRANSPOSE([[mlat_pos], [mltlon_pos]])

    ENDIF ELSE BEGIN
        rvecs = rad_map_eval_grad_vecs(extra_pos, coeffs[*,2], latmin, order)
        rmag = reform(sqrt(rvecs[0,*]^2 + rvecs[1,*]^2))

        q = where (rmag ne 0, qc)
        if (qc eq 0) then begin
          prinfo, 'All "real" vectors have 0 length.'
        endif

        raz    = fltarr(n_elements(rmag))
        if (int_hemi eq 1) then begin
              raz[q] = atan(rvecs[1,q],rvecs[0,q])*!radeg
        endif else begin
              raz[q] = atan(rvecs[1,q],-rvecs[0,q])*!radeg
        endelse

    ;    ; Now shift back into 'real word'
    ;    if (lat_shft ne 0.) then begin
    ;          xat_shft = -lat_shft
    ;          npnts    =  n_elements(rmag)
    ;          crd_shft, lon_shft, xat_shft, npnts, extra_pos, raz
    ;    endif

        ; convert magnetic longitude into mltlon
        extra_pos[1,*] = extra_pos[1,*] + lon_shft_tmp
    ENDELSE

    ; write data into the opened file
    sfjul, date_tmp, time_tmp, jul, /jul_to_date
    PRINTF, U, STRTRIM(STRING(date_tmp) + "-" + STRTRIM(STRING(time_tmp),1), 1), raz
    PRINTF, UU, STRTRIM(STRING(date_tmp) + "-" + STRTRIM(STRING(time_tmp),1), 1), rmag

    PRINTF, UUU, STRTRIM(STRING(date_tmp) + "-" + STRTRIM(STRING(time_tmp),1), 1), TRANSPOSE(extra_pos[0,*]) 
    PRINTF, UUUU, STRTRIM(STRING(date_tmp) + "-" + STRTRIM(STRING(time_tmp),1), 1),TRANSPOSE(extra_pos[1,*]) 
    ;print, [raz]
ENDFOR

extra_pos[1, *] = extra_pos[1, *] / 15.
vel_data  = {pos:extra_pos,  vectors:transpose([[rmag], [raz+180.]])}

;return, vel_data
print, time_tmp
print, vel_data

; close the file
CLOSE, U
CLOSE, UU
CLOSE, UUU
CLOSE, UUUU

RETURN, 1

END

;Note: argument 'coords' is just a place holder. this code only workds for mlt 
;e.g.
        ;.run rad_map_calc_fitvel
        ;real = rad_map_calc_fitvel_vecs(20030214, [1800, 2300], "mlt", 5, /fixed_points)

