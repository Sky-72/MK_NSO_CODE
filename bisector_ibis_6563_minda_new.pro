PRO bisector_ibis_6563_minda_new

;;----------------------------------------------------------
;;---------------------
;; pathes and filenames
;;--------------------- 
     iscan = 0
     line = '6563'

     inpath  = '/workspace/mkatilius/IDL/mData/20141025/prefiltcalib/6563_nbn/'
     outpath = '/workspace/mkatilius/IDL/mData/20141025/bisec/'

     file    = line + '_nbn*.sav'
     outfile = '6563_bisec_' + string(iscan)+ '.sav'
     posfile = 'reference_pos_' + line + '.sav'
     posfile1 = 'reference_pos_region_' + line + '.sav'
 ;    fpfilter = '/home/ali/idl/IBIS/pro/prefilter.profile.5896.Jun2006.sav'

;;--------------------------------------------------

     files = file_search(inpath + file, count=nfiles)

;;--------------------------------------------------
;;---------
;; aperture
;;---------
;;;----------
;; variables
;;----------

    nbisec = 11
    ibisec_start = 0.02
    ibisec_steps = 0.08

    nstart = 0
    nend = nfiles-1

    kernel = 9
    level = 2800.

    xa = 570 & xb = 810
    ya = 690 & yb = 810
;;-------------------------------------
;; Define intensity levels for bisector
;;-------------------------------------

    ibisec    = fltarr(nbisec)
    ibisec[0] = ibisec_start
    FOR i = 1, nbisec-1 DO ibisec[i] = ibisec[i-1] + ibisec_steps
;;---------------------------
;; start bisector calculation
;;---------------------------

     for iscan = nstart, nend do begin 

         restore, files[iscan], /verb
         nbdwc = nb_data_crrt
         s = size(nbdwc)
         nz = s[3]
         index = where(nbdwc eq 0)
         if index[0] ne -1 then nbdwc[index]=nbdwc[5,5,0]

;     aperture = reform(mwld[*,*,0]) ne 0. and reform(nbdwc[*,*,0]) ne 0.
;     mask = aperture eq 0.
;     index_out = where(aperture eq 0.)
;     index_in = where(aperture eq 1.)

;;--------------------------------------------------
;;-----------------------------
;; equidistant wavelength grid
;;----------------------------

         wtmp = info_nb.grid
         wlam = findgen(nz)/(nz-1) * (max(wtmp) - min(wtmp)) + min(wtmp)

;;--------------------------------------------------
;;-----------------------------
;; prefilter correction factors
;;-----------------------------

;     RESTORE, fpfilter, /verb
;     pfilter1 = interpol(prefilter_5896_lamp,scan_waves_5896, wtmp1)
;     pfilter1 = interpol(pfilter1, wtmp1, wtmp1-0.1)
;     pfilter = interpol(prefilter_5896_lamp,scan_waves_5896, wtmp)
;     pfilter = interpol(pfilter, wtmp, wtmp-0.1)

;;----------------------------------------------------------


;;----------------------------------------------------------


;;----------------------------------------------------------
;;------------------------------------------
;; arrays for reference wavelength positions
;;------------------------------------------

        wref     = FLTARR(nbisec, nfiles)
        wref_fft = FLTARR(nfiles)
        iqsavg   = fltarr(nz, nfiles)
        time     = fltarr(nz, nfiles)

        wrefr     = FLTARR(nbisec, nfiles)
        wrefr_fft = FLTARR(nfiles)
        iravg     = fltarr(nz, nfiles)



       ;;--------------------------------------------------
       ;; determine QS region from broadband for zero level
       ;;--------------------------------------------------

         mskout = shift(dist(s[1]), s[1]/2, s[2]/2) lt 455
         index_out = where(mskout eq 0.)
         for i=0,nz-1 do begin 
             tmp = reform(nbdwc[*,*,i])
             tmp[index_out] = reform(nbdwc[5,5,i])
             nbdwc[*,*,i] = tmp
             tmp = 0
             ;tmp = reform(mwld[*,*,i])
            ; tmp[index_out] = reform(mwld[5,5,i])
            ; mwld[*,*,i] = tmp
         endfor

         qsmask   = mskout eq 1. ;;;;; (smoothe(tmp, kernel) gt level) * mskout
         index_qs = where(qsmask eq 1.)

;         window,1
;         tvscl, qsmask ;* tmp

       ;;-----------------------------
       ;; determine QS average profile
       ;;-----------------------------

         pavg = fltarr(nz)
         rpavg = fltarr(nz)

         for k = 0, nz-1 do rpavg[k] = avg(nbdwc[xa:xb,ya:yb,k])
         for k = 0, nz-1 do pavg[k] = avg((nbdwc[*,*,k])[index_qs])
;         for k = 0, nz-1 do upavg[k] = avg((nbdwc[*,*,k])[index_u])

         pavg = pavg ;;;; /pfilter
         pavg = interpol(pavg,wtmp,wlam)
         rpavg = rpavg ;;;;; /pfilter
         rpavg = interpol(rpavg,wtmp,wlam)

       ;;--------------------------------------------------
       ;; determine bisector of average qs profile
       ;;--------------------------------------------------
          
          p    = pavg
          a    = max(p,imin)
          dixc = [0-imin]; ,nz-1-imin] 
          wixc = [0]

          BISEC_ABS_LIN, p, dixc, wixc, xl_avg, xr_avg, ycont_avg, ymin_avg, xmin_avg, ibisec = ibisec ;;, /plot
                 
          bavg = (xr_avg + xl_avg)/2. 
          wref[*,iscan] = bavg
          lpff, p, posavg_lpff  
          wref_fft[iscan] = posavg_lpff

          p    = rpavg
          a    = min(p,imin)
          dixc = [0-imin];,nz-1-imin] 
          wixc = [0]

          BISEC_ABS_LIN, p, dixc, wixc, xl_avg, xr_avg, ycont_avg, ymin_avg, xmin_avg, ibisec = ibisec ;;, /plot
                 
          bavg = (xr_avg + xl_avg)/2. 
          wrefr[*,iscan] = bavg
          lpff, p, posavg_lpff  
          wrefr_fft[iscan] = posavg_lpff

          iqsavg[*,iscan] = pavg
          iravg[*,iscan] = rpavg

       ;;-----------------------------
       ;; define arrays
       ;;-----------------------------

            BISEC    = FLTARR(s[1],s[2], nbisec)
            INT_CONT = FLTARR(s[1],s[2])
            INT_CORE = FLTARR(s[1],s[2])
            POS_CORE = FLTARR(s[1],s[2])
            POS_FFT  = FLTARR(s[1],s[2])
            EQVW     = FLTARR(s[1],s[2])
            FWHM     = FLTARR(s[1],s[2])

       ;;----------------------------------
       ;; do it: inner loop over all pixels
       ;;----------------------------------

         FOR cx = 0,s[1]-1 DO BEGIN

             FOR cy = 0,s[2]-1 DO BEGIN

                 p = reform(nbdwc[cx,cy,*]) ;;;; /pfilter
                 p = interpol(p,wtmp,wlam)
                 a = min(p,imin)
                 ave_p = mean(p)
                 
                 if (p[5] lt p[4]) or (p[6] lt p[5]) or (p[7] lt p[6]) and (p[6] lt 2200) and (p[5] lt 2150) $
                 and (p[7] lt 2200) and (p[13] gt 750) and (p[0] gt 900) then begin
                 
                   dixc = [0-imin];,nz-1-imin] 
                   wixc = [0]

                   BISEC_ABS_LIN, p, dixc, wixc, xl, xr, ycont, ymin, xmin,  $
                                ibisec = ibisec ;,/plot

                    BISEC[cx,cy,*] = (xr + xl)/2. 
                    INT_CONT[cx,cy]  = ycont[2]
                    INT_CORE[cx,cy]  = ymin[1]
                    POS_CORE[cx,cy]  = xmin[1] 
                  ;  EQVW[cx,cy]      = eqvw_i

                    FWHM[cx,cy]      = avg(xr[5:7] - xl[5:7])

                    index1 = 4  ;;;; bavg[0]-15.
                    index2 = 10 ;;;; bavg[0]+15.

                    lpff, p[index1:index2], pos_lpff  
                    POS_FFT[cx,cy]  = pos_lpff + index1
                    
                  endif else print, 'No curve'
                  
             ENDFOR ;cy

         ENDFOR ;cx

         print,' ... end of inner loop'

       ;;----------------------------------
       ;; end inner loop over all pixels
       ;;----------------------------------

       ;  irms = fltarr(nz)
      ;   for i = 0, nz-1 do irms[i] = stdev((mwld[*,*,i])[index_qs])/avg((mwld[*,*,i])[index_qs])
      ;   posmax = where(irms eq max(irms))
       ;  bestimage = mwld[*,*,posmax[0]]

;;         bestimage = nbdwc[*,*,0]

        ; tmp = nb_expos.time_string
         time[*,iscan] = time_conv_new(strmid(tmp,12))

         s_file = files[iscan]

         file = STRMID(s_file, (STRSPLIT(s_file, '/'))[N_ELEMENTS(STRSPLIT(s_file, '/'))-1])

         print,'save:', outpath + outfile 

         save,filename = outpath + outfile , $
              BISEC, INT_CONT, INT_CORE, POS_CORE, POS_FFT, EQVW, fwhm, $
              BESTIMAGE, ibisec, dixc, wixc, nb_expos, posavg_lpff, bavg, /verb

         window,2, xs = 500+250,ys = 500
         tvscl,rebin(-pos_core,250,250),0,250
         tvscl,rebin(-pos_fft,250,250),250,250
        ; tvscl,rebin(bestimage,250,250),0,0
         tvscl,rebin(int_core,250,250),250,0
         tvscl,rebin(fwhm,250,250),2*250,0

       ;;---------------------
         endfor ; end loop_out
       ;;---------------------

     print,'save:', outpath + posfile
         
     save, wref, wref_fft, iqsavg, time, filename = outpath + posfile

     print,'save:', outpath + posfile1
         
     save, wrefr, wrefr_fft, iravg, time, filename = outpath + posfile1

;;-----------------------------------------------------------------
;; end outer loop
;;-----------------------------------------------------------------

;;----------------------------------------------------------
;stop
end
