;+
; NAME: 
; RAD_MAP_PLOT_FOR_WIC
;
; PURPOSE: 
; This procedure plots a map with potential contours, convection boundary, velocity vectors
; and some scales, colorbar and a title.
; 
; CATEGORY: 
; Graphics
; 
; CALLING SEQUENCE:
; RAD_MAP_PLOT
;
; KEYWORD PARAMETERS:
; DATE: A scalar or 2-element vector giving the time range to plot, 
; in YYYYMMDD or MMMYYYY format.
;
; TIME: A 2-element vector giving the time range to plot, in HHII format.
;
; LONG: Set this keyword to indicate that the TIME value is in HHIISS
; format rather than HHII format.
;
; INDEX: The index number of the map to plot. If not set, that map with the timestamp closest
; to the gien input date will be plotted.
;
; NORTH: Set this keyword to plot the convection pattern for the northern hemisphere.
;
; SOUTH: Set this keyword to plot the convection pattern for the southern hemisphere.
;
; HEMISPHERE; Set this to 1 for the northern and to -1 for the southern hemisphere.
;
; COORDS: Set this to a string containing the name of the coordinate system to plot the map in.
; Allowable systems are geographic ('geog'), magnetic ('magn') and magnetic local time ('mlt').
;
; COAST: Set this keyword to plot coast lines.
;
; NO_FILL: Set this keyword to surpress filling of the coastal lines.
;
; CROSS: Set this keyword to plot a coordinate cross rather than a box.
;
; SCALE: Set this keyword to a 2-element vector containing the minimum and maximum velocity 
; used for coloring the vectors.
;
; MODEL: Set this keyword to include velocity vectors added by the model.
;
; MERGE: Set this keyword to plot velocity vectors
;
; TRUE: Set this keyword to plot velocity vectors
;
; LOS: Set this keyword to plot velocity vectors
;
; GRAD: Set this keyword to plot velocity vectors calculated from the ExB drift using the coefficients
; of the potential.
;
; NEW_PAGE: Set this keyword to plot multiple maps each on a separate page.
;
; COMMON BLOCKS:
; RAD_DATA_BLK: The common block holding map data.
;
; EXAMPLE:
;
; COPYRIGHT:
; Non-Commercial Purpose License
; Copyright © November 14, 2006 by Virginia Polytechnic Institute and State University
; All rights reserved.
; Virginia Polytechnic Institute and State University (Virginia Tech) owns the DaViT
; software and its associated documentation (“Software”). You should carefully read the
; following terms and conditions before using this software. Your use of this Software
; indicates your acceptance of this license agreement and all terms and conditions.
; You are hereby licensed to use the Software for Non-Commercial Purpose only. Non-
; Commercial Purpose means the use of the Software solely for research. Non-
; Commercial Purpose excludes, without limitation, any use of the Software, as part of, or
; in any way in connection with a product or service which is sold, offered for sale,
; licensed, leased, loaned, or rented. Permission to use, copy, modify, and distribute this
; compilation for Non-Commercial Purpose is hereby granted without fee, subject to the
; following terms of this license.
; Copies and Modifications
; You must include the above copyright notice and this license on any copy or modification
; of this compilation. Each time you redistribute this Software, the recipient automatically
; receives a license to copy, distribute or modify the Software subject to these terms and
; conditions. You may not impose any further restrictions on this Software or any
; derivative works beyond those restrictions herein.
; You agree to use your best efforts to provide Virginia Polytechnic Institute and State
; University (Virginia Tech) with any modifications containing improvements or
; extensions and hereby grant Virginia Tech a perpetual, royalty-free license to use and
; distribute such modifications under the terms of this license. You agree to notify
; Virginia Tech of any inquiries you have for commercial use of the Software and/or its
; modifications and further agree to negotiate in good faith with Virginia Tech to license
; your modifications for commercial purposes. Notices, modifications, and questions may
; be directed by e-mail to Stephen Cammer at cammer@vbi.vt.edu.
; Commercial Use
; If you desire to use the software for profit-making or commercial purposes, you agree to
; negotiate in good faith a license with Virginia Tech prior to such profit-making or
; commercial use. Virginia Tech shall have no obligation to grant such license to you, and
; may grant exclusive or non-exclusive licenses to others. You may contact Stephen
; Cammer at email address cammer@vbi.vt.edu to discuss commercial use.
; Governing Law
; This agreement shall be governed by the laws of the Commonwealth of Virginia.
; Disclaimer of Warranty
; Because this software is licensed free of charge, there is no warranty for the program.
; Virginia Tech makes no warranty or representation that the operation of the software in
; this compilation will be error-free, and Virginia Tech is under no obligation to provide
; any services, by way of maintenance, update, or otherwise.
; THIS SOFTWARE AND THE ACCOMPANYING FILES ARE LICENSED “AS IS”
; AND WITHOUT WARRANTIES AS TO PERFORMANCE OR
; MERCHANTABILITY OR ANY OTHER WARRANTIES WHETHER EXPRESSED
; OR IMPLIED. NO WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE IS
; OFFERED. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF
; THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE,
; YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR
; CORRECTION.
; Limitation of Liability
; IN NO EVENT WILL VIRGINIA TECH, OR ANY OTHER PARTY WHO MAY
; MODIFY AND/OR REDISTRIBUTE THE PRORAM AS PERMITTED ABOVE, BE
; LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL,
; INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR
; INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS
; OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED
; BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE
; WITH ANY OTHER PROGRAMS), EVEN IF VIRGINIA TECH OR OTHER PARTY
; HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
; Use of Name
; Users will not use the name of the Virginia Polytechnic Institute and State University nor
; any adaptation thereof in any publicity or advertising, without the prior written consent
; from Virginia Tech in each case.
; Export License
; Export of this software from the United States may require a specific license from the
; United States Government. It is the responsibility of any person or organization
; contemplating export to obtain such a license before exporting.
;
; MODIFICATION HISTORY: 
; Written by Lasse Clausen, Dec, 11 2009
; Modified by Evan Thomas, Apr, 09 2014
;-
pro rad_map_plot_for_wic, date=date, time=time, long=long, $
	coords=coords, index=index, scale=scale, new_page=new_page, $
	north=north, south=south, hemisphere=hemisphere, $
	xrange=xrange, yrange=yrange, rotate=rotate, $
	cross=cross, coast=coast, no_fill=no_fill, orig_fan=orig_fan, $
	model=model, merge=merge, true=true, los=los, grad=grad, $
	vec_radar_ids=vec_radar_ids, fan_radar_ids=fan_radar_ids, $
	vectors=vectors, potential=potential, efield=efield, $
	comp_east_west=comp_east_west, comp_north_south=comp_north_south, $
	isotropic=isotrpoic,dt_diff_newplot=dt_diff_newplot, contour_fill = contour_fill, $
	rad_sct_flg_val = rad_sct_flg_val, no_fov_names = no_fov_names, no_plot_imf_ind = no_plot_imf_ind, $
	dms_ssj_overlay = dms_ssj_overlay, dms_ssies_overlay = dms_ssies_overlay, poes_flux_overlay = poes_flux_overlay, $
	overlay_r1_oval = overlay_r1_oval, overlay_r2_oval = overlay_r2_oval, noOverlayHMB = noOverlayHMB, $
	zoom_map = zoom_map, zoom_xr = zoom_xr, zoom_yr = zoom_yr, website = website

common rad_data_blk
common omn_data_blk
common aur_data_blk
common kpi_data_blk
common pos_data_blk
common r12oval_data_blk


if ~keyword_set(vectors) and ~keyword_set(potential) and ~keyword_set(efield) then begin
	if ~keyword_set(silent) then $
		prinfo, 'Nothing selected to plot, choosing vectors and potential.'
	vectors = 1
	potential = 1
	efield = 0
endif

; check hemisphere and north and south
if ~keyword_set(hemisphere) then begin
	if keyword_set(north) then $
		hemisphere = 1. $
	else if keyword_set(south) then $
		hemisphere = -1. $
	else $
		hemisphere = 1.
endif

hemi_orig_chosen=hemisphere

if ~keyword_set(rad_sct_flg_val) then $
	rad_sct_flg_val = 0

if (~keyword_set(contour_fill)) then $
	contour_fill_final = 1	$
else $
	contour_fill_final = 0

if ~keyword_set(dt_diff_newplot) then $
	dt_diff_newplot = 2.


if ~keyword_set(scale) then $
	scale = [0,2000]

;; set the plotting range for the map	
if ( ~keyword_set(zoom_map) ) then begin
	if n_elements(yrange) ne 2 then $
		yrange = [-46,46]

	if n_elements(xrange) ne 2 then $
		xrange = [-46,46]
endif else begin
	xrange = zoom_xr
	yrange = zoom_yr
endelse

if ~keyword_set(no_fov_names) then $
	no_fov_names = 0	


aspect = float(xrange[1]-xrange[0])/float(yrange[1]-yrange[0])

; this makes int_hemi 0 for north and 1 for south
int_hemi = (hemisphere lt 0)

if rad_map_info[int_hemi].nrecs eq 0L then begin
	if ~keyword_set(silent) then $
		prinfo, 'No data loaded.'
	return
endif

if n_elements(index) ne 0 then $
	sfjul, date, time, (*rad_map_data[int_hemi]).mjuls[index], /jul

if ~keyword_set(date) then begin
	if ~keyword_set(silent) then $
		prinfo, 'No DATE given, trying for scan date.'
	caldat, (*rad_map_data[int_hemi]).sjuls[0], month, day, year
	date = year*10000L + month*100L + day
endif

if n_elements(time) lt 1 then $
	time = 1200

sfjul, date, time, sjul, fjul

; sample time of maps
; in minutes
if rad_map_info[int_hemi].nrecs ge 3L then $
	dt = mean(deriv((*rad_map_data[int_hemi]).mjuls*1440.d)) $
else $
	dt = dt_diff_newplot

; account for sjul being before the
; date/time of the first map
sjul = ( sjul > (*rad_map_data[int_hemi]).sjuls[0] )

if n_elements(time) eq 2 then begin
	npanels = round((fjul-sjul)*1440.d/dt_diff_newplot)
endif else begin
	npanels = 1
endelse
; calculate number of panels per page
if npanels eq 1 then begin
	xmaps = 1
	ymaps = 1
endif else if npanels eq 2 then begin
	xmaps = 2
	ymaps = 1
endif else if npanels le 4 then begin
	xmaps = 2
	ymaps = 2
endif else if npanels le 6 then begin
	xmaps = 3
	ymaps = 2
endif else begin
	xmaps = floor(sqrt(npanels)) > 1
	ymaps = ceil(npanels/float(xmaps)) > 1
endelse

; take into account format of page
; if landscape, make xmaps > ymaps
fmt = get_format(landscape=ls, sardines=sd)
if ls then begin
	if ymaps gt xmaps then begin
		tt = xmaps
		xmaps = ymaps
		ymaps = tt
	endif
; if portrait, make ymaps > xmaps
endif else begin
	if xmaps gt ymaps then begin
		tt = ymaps
		ymaps = xmaps
		xmaps = tt
	endif
endelse

;clear_page
set_format, /sardi

if keyword_set(website) then begin
    colortable = get_colortable()
    rad_load_colortable,/website
endif

mpos = define_panel(xmaps, 1, xmaps-1, 0, aspect=aspect, /bar) - [.06, .075, .06, .075]
;mpos = define_panel(xmaps, 1, xmaps-1, 0, aspect=aspect) - [.09, +.025, .09, +.075]  ; used for 1x1 plot
if (~keyword_set(new_page)) then begin
	cb_pos = define_cb_position(mpos, height=55, gap=.1*(mpos[2]-mpos[0]))
	if keyword_set(vectors) then $
		plot_colorbar, /square, scale=scale, parameter='velocity', position=cb_pos,/keep_first_last_label, $
			/no_rotate, charsize = 0.65
	if keyword_set(orig_fan) then begin
		cb_pos = define_cb_position(mpos, height=35, gap=.065*(mpos[2]-mpos[0]))
		plot_colorbar, /square, scale=.5*[-scale[1],scale[1]], parameter='velocity',/keep_first_last_label, $
			/left, position=cb_pos, legend=' ', charsize = 0.5
	endif
endif

if keyword_set(website) then begin
    rad_load_colortable,colortable
endif


; loop through panels
for b=0, npanels-1 do begin

	asjul = sjul + double(b)*dt_diff_newplot/1440.d
	sfjul, date, time, asjul, /jul_to_date
	; print,dt,dt_diff_newplot,time,date
	; calculate index from date and time
	if n_elements(index) eq 0 then begin
		dd = min( abs( (*rad_map_data[int_hemi]).mjuls-asjul ), _index)
		; check if time ditance is not too big
		if dd*1440.d gt 60. then $
			prinfo, 'Map found, but '+string(dd*1440.d,format='(I4)')+' minutes off chosen time.'
	endif else begin
		asjul = (*rad_map_data[int_hemi]).sjuls[index]
		sfjul, date, time, (*rad_map_data[int_hemi]).sjuls[index], /jul_to_date
		_index = index
	endelse
	amjul = (*rad_map_data[int_hemi]).mjuls[_index]

	if keyword_set(new_page) then begin
		clear_page
		xmaps = 1
		ymaps = 1
		xmap = 0
		ymap = 0
	endif else begin
		xmap = b mod xmaps
		ymap = b/xmaps
	endelse

	if b eq 0 or keyword_set(new_page) then begin
		opos = define_panel(1, 1, 0, 0, aspect=aspect, /bar) - [.06, .005, .06, .005]

		if (n_elements(tm_orig_set) eq 1) then begin
			orange = [(sjul+fjul)/2.d + [-1.d,1.d]*30.d/1440.d]
		endif else begin
			orange = [sjul-(30.d/1440.d),fjul+(30.d/1440.d)]
		endelse
		
; 		orange = [amjul + [-1.d,1.d]*30.d/1440.d]
		
		sfjul, odate, otime, orange, /jul_to_date
        ; Reads data from omni (not omniex) files.
        omn_read, odate, time=otime, /force

		;;check if omni data is present if not print 'NO-DATA'
		jinds_omn_exist=where(finite(omn_data.bz_gsm),check_omn)

		if (check_omn ne 0) then begin
		
			jinds_omn_by_finite=where(finite(omn_data.by_gsm))
			jinds_omn_bz_finite=where(finite(omn_data.bz_gsm))
			jinds_omn_pd_finite=where(finite(omn_data.pd))

			yrange_omn_min_val=min([min(omn_data.by_gsm[jinds_omn_by_finite]),min(omn_data.bz_gsm[jinds_omn_bz_finite])])
			yrange_omn_max_val=max([10,max(omn_data.by_gsm[jinds_omn_by_finite]),max(omn_data.bz_gsm[jinds_omn_bz_finite])])
			yrange_omn=[5*(round(yrange_omn_min_val/5)-1),5*(round(yrange_omn_max_val/5)+1)]
			yrange_imf_bybz_plot=[5*(round(yrange_omn_min_val/5)-1),5*(round(yrange_omn_max_val/5)+1)]


			yrange_pd_min_val=min([0,min(omn_data.pd[jinds_omn_pd_finite])])
			yrange_pd_max_val=max([5,max(omn_data.pd[jinds_omn_pd_finite])])
			yrange_imf_pd_plot=[max([5*(round(yrange_pd_min_val/5)-1),0]),5*(round(yrange_pd_max_val/5))]
		
			yticks_pd_omn = fix(((yrange_imf_pd_plot[1])-(yrange_imf_pd_plot[0]))/5)
		
			if (fix(((yrange_imf_bybz_plot[1])-(yrange_imf_bybz_plot[0]))/5) lt 6) then begin
				yticks_imf_omn = fix(((yrange_imf_bybz_plot[1])-(yrange_imf_bybz_plot[0]))/5) 
			endif else begin
				if (fix(((yrange_imf_bybz_plot[1])-(yrange_imf_bybz_plot[0]))/10) lt 6) then begin
					yticks_imf_omn = fix(((yrange_imf_bybz_plot[1])-(yrange_imf_bybz_plot[0]))/10)
				endif else begin
					yticks_imf_omn = fix(((yrange_imf_bybz_plot[1])-(yrange_imf_bybz_plot[0]))/15)	
				endelse
			endelse
			

			oopos = [opos[0], opos[3]+.11, opos[2], opos[3]+.2]

			;Plot the IMF panel
			if (~keyword_set(no_plot_imf_ind)) then begin
				omn_plot_panel, date=odate, time=otime, position=oopos, yrange=yrange_imf_bybz_plot, ystyle = 1, $
					param='by_gsm', yticks=yticks_imf_omn, charsize='0.5', xstyle=1, /first, linecolor=get_blue(), ytitle='OMNI-IMF[nT]', linethick=2
				omn_plot_panel, date=odate, time=otime, position=oopos, yrange=yrange_imf_bybz_plot, ystyle=5, $
					param='bz_gsm', charsize='0.5', xstyle=5, /first, linecolor=253, ytitle='[nT]', linethick=2
				omn_plot_panel, date=odate, time=otime, position=oopos, yrange=yrange_imf_pd_plot, ystyle=5, $
					param='pd', charsize='0.5', xstyle=5, /first, linecolor=get_black(), ytitle='Pd[nPa]', linethick=2
				axis,yaxis=1,ytitle='Pd[nPa]',ystyle = 1,color=get_black(),yticks=yticks_imf_omn,charsize='0.5',yrange=yrange_imf_pd_plot
				xyouts, amjul - 28.d/1440.d, !y.crange[0]+.13*(!y.crange[1]-!y.crange[0]), $
					'By', color=get_blue(), charsize='0.5', charthick=2
       			xyouts, amjul - 28.d/1440.d, !y.crange[1]-.18*(!y.crange[1]-!y.crange[0]), $
					'Bz', color=253, charsize='0.5', charthick=2
				xyouts, amjul + 25.d/1440.d, !y.crange[1]-.28*(!y.crange[1]-!y.crange[0]), $
					'Pd', color=get_black(), charsize='0.5', charthick=2
				oplot, !x.crange, replicate(.5*!y.crange[0], 2), linestyle=1, color=get_gray()
				oplot, !x.crange, replicate(.5*!y.crange[1], 2), linestyle=1, color=get_gray()
				oplot, replicate(amjul,2), !y.crange, linestyle=2, color=252
				axis, /xaxis, /xstyle, xrange=orange, xticks='6', $;get_xticks(orange[0], orange[1]), $
					charsize='0.6', xtickformat='label_date',xtitle='UT'
			endif
		
        endif
		
;;;; Plot Sym-H, Ae and Kp indices

;;;; Plot Aur/Sym indices

		oopos = [opos[0], 0.766909, opos[2], 0.846909]

		if (~keyword_set(no_plot_imf_ind)) then begin
			kpi_read,odate, time=[0000,2400], /force
			kpi_file_check = kpi_check_loaded( date )

			if ( kpi_file_check eq 1 ) then begin

				;print,'ooo',otime,orange
				kpi_plot_panel, date=odate, time=otime, position=oopos,xtickname=REPLICATE(' ', 40),ytickname=REPLICATE(' ', 40), $
					ytitle=' ', xstyle=1, ystyle = 1,charsize=0.0001, xticks=12, xminor=6, $
					yminor=2,linecolor=get_green(),yticks=1

				for kk=0,n_elements(kpi_data.juls)-1 do begin	
; 					oplot,kpi_data.juls,kpi_data.kp_index,color=get_green(),linestyle=1,thick=2
					oplot,[kpi_data.juls[kk],kpi_data.juls[kk]+3.d/24.d],[kpi_data.kp_index[kk],kpi_data.kp_index[kk]],color=get_green(),thick=2
				endfor

				axis,orange[1]+0.1*(orange[1]-orange[0]),yaxis=1, ytitle='Kp-Index',color=get_green(),charsize='0.5',ystyle=1,yminor=3,yticks=3,yrange=[0,9],ticklen=-0.005
;				axis, yaxis=1, ystyle=1, yrange=[0,9], yticks=3, $
;					ytickname=strtrim(string([0,3,6,9],format=level_format),2), ytitle='Kp-Index', charthick='1', charsize='0.5'

				aur_read, odate, time=otime, /force
				jul_asy_ae_arr=dblarr(n_elements(aur_data.ae_index))
				for aeas=0,n_elements(aur_data.ae_index)-1 do begin
					jul_asy_ae_arr[aeas]=orange[0]+1.d*aeas/1440.d
				endfor

; 				plot,jul_asy_ae_arr,aur_data.ae_index,position=oopos, yrange=[round(min([aur_data.ae_index]))-100,round(max([aur_data.ae_index]))+100], /ystyle, $
;  					yticks=5,xrange=[orange[0],orange[1]],xticks='6',xtickname=REPLICATE(' ', 40), charsize='0.5', xstyle=1, ytitle='AE-Ind[nT]',thick='4'

				plot,jul_asy_ae_arr,aur_data.ae_index,position=oopos, yrange=[0,1500], ystyle = 1, $
					yminor=2,yticks=3,xrange=[orange[0],orange[1]],xticks='6',xtickname=REPLICATE(' ', 40),ytickname=REPLICATE(' ', 40), charsize='0.5', xstyle=1,thick='4',color=get_blue()
				plot,jul_asy_ae_arr,aur_data.sym_h,position=oopos, yrange=[-120,60], /ystyle, $
					yminor=2,yticks=3,xrange=[orange[0],orange[1]],xticks='6',xtickname=REPLICATE(' ', 40), charsize='0.5', xstyle=1, ytitle='Sym-H[nT]',thick='4'
				axis,yaxis=1, ytitle='AE-Ind[nT]',color=get_blue(),yticks=3,charsize='0.5',yrange=[0,1500],yminor=2
			endif
		endif

;;;;Plot Npts and Phi-PC

		oopos = [opos[0], opos[3]+.03, opos[2], opos[3]+.08]
; 		print,odate,otime,b,npanels

        if (~keyword_set(no_plot_imf_ind)) then begin
			rad_map_plot_npoints_panel, date=odate, time=otime, position=oopos, yrange=[1e1,1e3], ystyle=5, $
				charsize='0.5', xstyle=5, /ylog, linethick=2, hemisphere=hemisphere
			rad_map_plot_potential_panel, date=odate, time=otime, position=oopos, yrange=[30,130], ystyle=9, $
				charsize='0.5', xstyle=9, /first, linecolor=200, linethick=2, hemisphere=hemisphere
			xyouts, (sjul+fjul)/2.d - 30.d/1440.d + 2.d/1440.d, !y.crange[0]+.13*(!y.crange[1]-!y.crange[0]), $
				textoidl('\Phi_{pc}'), color=200, charsize=get_charsize(1,2), charthick=2
			xyouts, (sjul+fjul)/2.d - 30.d/1440.d + 2.d/1440.d, !y.crange[1]-.18*(!y.crange[1]-!y.crange[0]), $
				'Npts', charsize='0.5', charthick=2
			axis, /yaxis, ystyle=1, yrange=[1e1,1e3], yticks=2, /ylog, charsize='0.5', ytitle='Npts'
			oplot, replicate(amjul,2), !y.crange, linestyle=2, color=252
            if (check_omn eq 0) then begin
                ; If there's no OMNI data, at least plot the time axis across the top here.
                axis, /xaxis, /xstyle, xrange=orange, xticks='6', $
                    charsize='0.6', xtickformat='label_date',xtitle='UT'
            endif
        endif
	endif

        if ~keyword_set(position) then $
		_position = define_panel(xmaps, ymaps, xmap, ymap, aspect=aspect, /bar) - [.06, .1, .06, .1] $ ;[.06, .075, .06, .075] 
	else $
		_position = position

;	rad_map_plot_imf, xmaps, ymaps, xmap, ymap, gap=.04*get_charsize(xmaps,ymaps), $
;		index=_index, size=.125/sqrt(xmaps > ymaps)*(_position[2]-_position[0]), $
;		int_hemisphere=int_hemi, panel_position=_position

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ; overlay WIC image
        ;fname_wic = "./from_harald_frey/IMFHWIC_2002_0318_160426i.sav"
        IF time < 1631 THEN BEGIN
            fname_wic = "./from_harald_frey/IMFHWIC_2002_0318_" + STRCOMPRESS(STRING(time-1), /REMOVE_ALL) + "*.sav" ; for 1600 to 1630
        ENDIF ELSE BEGIN
            fname_wic = "./from_harald_frey/IMFHWIC_2002_0318_" + STRCOMPRESS(STRING(time), /REMOVE_ALL) + "*.sav"  ; for 1632 to 1658
        ENDELSE
        restore, fname_wic
        image = imageinfo.Mlt_img
        lats = imageinfo.Mlat
        lons = imageinfo.Mlon
        mlts = imageinfo.Mlt
        image_size = SIZE(image, /DIMENSIONS)
        height = image_size[0]
        width = image_size[1]
        cgLoadCT, 0, /REVERSE
        cgIMAGE, image,  POSITION=_position, MAXVALUE=1000
        ;p = [0.02, 0.99, 0.98, 0.98]
        p = cb_pos
        ;cgColorbar, Position=[p[0]+0.1, p[1], p[2]+0.1, p[1]-0.05], MAXRANGE=1000, $
        ;cgColorbar, MAXRANGE=1000, $
        ;/VERTICAL, /RIGHT

        ;if keyword_set(website) then begin
        rad_load_colortable, /website
        ;endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


	;rad_map_plot_vector_scale, xmaps, ymaps, xmap, ymap, gap=.05*get_charsize(xmaps,ymaps), $
	;	scale=scale, xrange=xrange, factor=factor*480., panel_position=_position

	if (keyword_set(new_page)) then begin
		cb_pos = define_cb_position(_position, height=35, gap=.01*(_position[2]-_position[0]))
		if keyword_set(vectors) then $
			plot_colorbar, /square, scale=scale, parameter='velocity',/keep_first_last_label, position=cb_pos, $
				/no_rotate, charsize = 0.65
		if keyword_set(orig_fan) then begin
			cb_pos = define_cb_position(mpos, height=35, gap=.065*(mpos[2]-mpos[0]))
			plot_colorbar, /square, scale=.5*[-scale[1],scale[1]], parameter='velocity',/keep_first_last_label, $
				/left, position=cb_pos, legend=' ', charsize = 0.5
		endif
	endif

	my_rad_map_plot_panel, xmaps, ymaps, xmap, ymap, $
		date=date, time=time, long=long, rotate=rotate, $
		north=north, south=south, hemisphere=hemisphere, $
		coords=coords, index=_index, scale=scale, $
		no_fill=no_fill, cross=cross, coast=coast, $
		model=model, merge=merge, true=true, los=los, grad=grad, $
		xrange=xrange, yrange=yrange, factor=factor, orig_fan=orig_fan, $
		vec_radar_ids=vec_radar_ids, fan_radar_ids=fan_radar_ids, $
		position=_position, fill = contour_fill_final , noOverlayHMB = noOverlayHMB, $
		vectors=vectors, potential=potential, efield=efield, website=website, $
		comp_east_west=comp_east_west, comp_north_south=comp_north_south,isotropic=isotrpoic, rad_sct_flg_val=rad_sct_flg_val, no_fov_names = no_fov_names
	; rad_map_plot_panel somehow recasts hemisphere as a 1 element array, which sucks! ;; Also if plotting southern hemisphere it changes back to +1 after one plot...

	hemisphere = hemi_orig_chosen

	if ( keyword_set(overlay_r1_oval) )  then begin
	
	; Read the R1 oval data
		r12_oval_read, date, time = time, hemisphere = hemisphere, /r1_oval, /silent
		if ( r12oval_data.juls[0] ne 0.d ) then begin ;; Check if the data is there
			r12_oval_map_overlay, date , time = time, linecolor = get_black(), linethick = 5, map_rotate=rotate
		endif
		
	endif

	if ( keyword_set(overlay_r2_oval) )  then begin
		; Read the R1 oval data
		r12_oval_read, date, time = time, hemisphere = hemisphere, /r2_oval, /silent
		if ( r12oval_data.juls[0] ne 0.d ) then begin ;; Check if the data is there
			r12_oval_map_overlay, date , time = time, linecolor = get_black(), linethick = 5, linestyle = 2, map_rotate=rotate
		endif
		
	endif

    if keyword_set(website) then $
        rad_load_colortable,/website
	
	if ( keyword_set(poes_flux_overlay) )  then begin
		rad_map_overlay_poes, date, time, coords = coords, $
			hemisphere = hemisphere, map_rotate=rotate
			
		rad_map_overlay_poes_bnd, date, time, hemisphere = hemisphere, coords = coords, $
			fitline_color = get_red(), fitline_style = 3, $
			fitline_thick = 5, map_rotate=rotate


; Moved after DMSP overlays so SSIES data isn't plotted inside colorbar - EGT, 20140423			
;		if ( ~keyword_set(dms_ssj_overlay) ) then begin
;		
;			plot_colorbar, 1, 2.55, 0.63, 1.03, /square, parameter='power',/keep_first_last_label, $
;				/left, charsize = 0.5, legend = TeXtoIDL('Total. Energy Flux [ergs cm^{-2} s^{-1}]'), scale = [-3, 1 ], level_format = '(f6.2)'
;		
;		endif
	
	endif	

	if ( keyword_set(dms_ssj_overlay) or keyword_set(dms_ssies_overlay) )  then begin
		rad_map_overlay_dmsp, date, time, coords = coords, scale_dmsp_ssies_vel = scale, $
			hemisphere = hemisphere,ssj = dms_ssj_overlay, ssies = dms_ssies_overlay, map_rotate=rotate
		if ( keyword_set(dms_ssj_overlay) ) then begin
			plot_colorbar, 1, 2.55, 0.63, 1.03, /square, parameter='power',/keep_first_last_label, $
				/left, charsize = 0.5, legend = TeXtoIDL('Log Elec. Energy Flux [keV/cm^{2} s sr]');'Log Elec. Energy Flux [Kev/cm^2 s sr]'
		endif
		
	
	endif

	if ( keyword_set(poes_flux_overlay) and ~keyword_set(dms_ssj_overlay) ) then begin
		plot_colorbar, 1, 2.55, 0.63, 1.03, /square, parameter='power',/keep_first_last_label, $
			/left, charsize = 0.5, legend = TeXtoIDL('Total Energy Flux [ergs cm^{-2} s^{-1}]'), scale = [-3, 1], level_format = '(f6.2)'
	endif
	
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120307_0450_0452_sphi1_sdpr2_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120307_0900_0902_sphi1_sdpr2_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120309_0614_0616_sphi1_sdpr2_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120309_0852_0854_sphi1_sdpr2_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120311_1520_1522_sphi1_sdpr2_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120312_1210_1212_ipp350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120312_1620_1622_ipp350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120315_1450_1452_sphi1_sdpr2_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120315_1930_1932_ipp350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120316_2220_2222_sphi1_sdpr2_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120307_0800_0802_sphi1+sdpr2_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120307_0800_0802_sphi1+sdpr2_ele30_SHipx350km.txt'

;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/120312_0920_0922_sphi1_sdpr2_ele30_sine_ipx350km.txt'

;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/150317_0910_0915_sphi1_sdpr2mm_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/150317_1830_1835_sphi1_sdpr2mm_ele30_sine_ipx350km.txt'
    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/150317_1325_1330_sphi1_sdpr2mm_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/150317_1840_1845_sphi1_sdpr2mm_ele30_sine_ipx350km.txt'
;    fipp_igs = '/davit/lib/vt/idl/tec/paul/new/150317_1900_1905_sphi1_sdpr2mm_ele30_sine_ipx350km.txt'

    ipp_data_igs = igsipp_paul(fipp_igs)

	Name_igs_ipp = ipp_data_igs[*,0]
	mlat_igs_ipp = float( ipp_data_igs[*,1] )
	mlt_igs_ipp = float( ipp_data_igs[*,2] )

    ; Overlay pierce points from igs stations
	load_usersym, /circle
	for jii = 0, n_elements( Name_igs_ipp )-1 do begin
;        oplot, [ mlat_igs_ipp[jii] ], [ mlt_igs_ipp[jii] ], psym = 8, symsize = 0.4, thick = 2.5;, color = color;, color = get_red()
	endfor

;    glat = 65.13
;    glon = -147.47
;    tmp = cnvcoorddavit(glat,glon,1.)
;    mlat = reform(tmp[0,*])
;    mlon = reform(tmp[1,*])
;    sfjul,date,time,asdfjul
;    caldat, asdfjul, mm, dd, year
;    yrsec = (asdfjul-julday(1,1,year,0,0,0))*86400.d
;    lon = mltdavit(year,yrsec,mlon[0])
;    in_mlt=!true
;    tmp = calc_stereo_coords(mlat, lon, mlt=in_mlt, rotate=rotate)
;    plots, tmp[0], tmp[1], psym=8,symsize=0.8, thick=2.5

	if keyword_set(website) then $
		rad_load_colortable,colortable

	my_rad_map_plot_title, position=_position, index=_index, coords=coords, $
		charsize=get_charsize(1,2), int_hemisphere=int_hemi, /silent, /no_disclaimer_str;, /no_earth_str

endfor

end
