;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;; 
;; Main procedure to calibrate the IBIS observations for darks
;; gain, blueshift and image distortions (destretch). The procedure
;; uses the output of ibis_combine.pro as input data.
;;
;; ref:      reference image(s)
;;    	    if equal to a scalar > 0, then a running mean is used
;; data:     data cube do be destreched
;; ks:       vector with kernel sizes, in the order the alignment is to be done
;;           rigid alignment corresponds to a zero, needs to be the first index
;; sfp:      2x2 array that describes lower left and upper right positions
;;           within which the 'rigid' alignment is to be done
;;           e.g. sfp = [ [xll, yll], [xur, yur] ]
;; dp:       detrend parameter
;; 
;; shifts:   return shifts on finest grid, if desired
;; grid:     return finest grid, if desired
;; 
;;
;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION destretch_ibis, ref, data, ks, sfp, dp, SHIFTS = shifts, GRID = grid

;;***************************************************************
;; copy data
;;***************************************************************

     cube = data

;;***************************************************************
;; dimension of ref
;;***************************************************************

     dimr = SIZE(ref)

     IF (dimr[0] LT 2) THEN BEGIN

        dimr = [dimr[0], 1, 1, 1]
        tr = ref

     ENDIF

     IF (dimr[0] EQ 2) THEN dimr = [dimr[0:2], 1]

;;***************************************************************
;; dimension of cube
;;***************************************************************

     dim = SIZE(cube)
     IF (dim[0] LT 2) THEN message, '2D or 3D data required as input.'
     IF (dim[0] EQ 2) THEN dim = [dim[0:2], 1]
     IF (dim[0] GT 3) THEN message, '2D or 3D data required as input.'

;;***************************************************************
;; computations for multi-reference destretch
;;***************************************************************

     nim1 = dim[3]/dimr[3]	; even number
     nim2 = dim[3] MOD nim1	; remainder

;;***************************************************************
;; flag for combination of shifts
;;***************************************************************

     flag = 0

;;***************************************************************
;; rigid alignment, if wanted
;;***************************************************************

     IF (ks[0] eq 0) THEN BEGIN

        sft = FLTARR(2, 1, 1, dim[3])

;;***************************************************************
;; do the first chunk with the even numbers
;;***************************************************************

        FOR j = 0, dimr[3]-1 DO BEGIN

            FOR i = j*nim1, (j+1)*nim1-1 DO BEGIN

                IF (dimr[0] LT 2) THEN IF (i GT 0) THEN ref = avg(cube[*,*,(i-tr)>0:i], 2) ELSE ref = cube[*,*,i]

                sft[*,0,0,i] = shc( ref[sfp[0,0]:sfp[0,1], sfp[1,0]:sfp[1,1], j], $
                                    cube[sfp[0,0]:sfp[0,1], sfp[1,0]:sfp[1,1], i], /interp )
                cube[*,*,i] = shift_bicub(cube[*,*,i], sft[0,0,0,i], sft[1,0,0,i])

           ENDFOR

        ENDFOR

;;***************************************************************
;; do the remainder, if neccessary
;;***************************************************************

        IF (nim2 GT 0) THEN BEGIN

           FOR i = dimr[3]*nim1, dim[3]-1 DO BEGIN

               IF (dimr[0] LT 2) THEN IF (i GT 0) THEN ref = avg(cube[*,*,(i-tr)>0:i], 2) ELSE ref = cube[*,*,i]
               sft[*,0,0,i] = shc( ref[sfp[0,0]:sfp[0,1], sfp[1,0]:sfp[1,1], dimr[3]-1], $
                                   cube[sfp[0,0]:sfp[0,1], sfp[1,0]:sfp[1,1], i], /interp )
               cube[*,*,i] = shift_bicub(cube[*,*,i], sft[0,0,0,i], sft[1,0,0,i])

           ENDFOR

        ENDIF
    
        flag = 1
        rt = -1
        IF (N_ELEMENTS(ks) eq 1) THEN BEGIN
            IF arg_present(shifts) THEN shifts = sft
            RETURN, cube 
        ENDIF ELSE ks = ks[1:*]

     ENDIF

;;***************************************************************
;; dimension of (remaining) kernel size vector
;;***************************************************************

     run = (SIZE(ks,/dimensions))[0]

;;***************************************************************
; destretch
;;***************************************************************

     FOR j = 0, run-1 DO BEGIN

         ; get shifts
   
         kernel = BYTARR(ks[j], ks[j])
         r      = mkcps(ref, kernel)

         IF (dimr[0] LT 2) THEN BEGIN

            shft = FLTARR((SIZE(r))[1], (SIZE(r))[2], (SIZE(r))[3], dim[3])

            FOR i = 0, dim[3]-1 DO BEGIN
         
                IF (i GT 0) THEN ref = avg(cube[*,*,(i-tr)>0:i], 2) ELSE ref = cube[*,*,i]

                shft[*,*,*,i] = cps(cube[*,*,i], ref, kernel)

            ENDFOR

         ENDIF ELSE IF (dimr[3] gt 1) THEN BEGIN

                       shft = fltarr((SIZE(r))[1], (SIZE(r))[2], (SIZE(r))[3], dim[3])

                       ; normal
 
                       FOR i = 0, dimr[3]-2 DO BEGIN
                           shft[*,*,*,i*nim1:(i+1)*nim1] = cps(REFORM(cube[*,*,i*nim1:(i+1)*nim1]), $ 
                                                               REFORM(ref[*,*,i]), kernel)
                       endfor

                       ; remainder

                       shft[*,*,*,(dimr[3]-1)*nim1:*] = cps(REFORM(cube[*,*,(dimr[3]-1)*nim1:*]), $ 
                                                            REFORM(ref[*,*,dimr[3]-1]), kernel)
         
                   ENDIF ELSE BEGIN

                           shft = cps(cube, ref, kernel)

                           ENDELSE

        ; detrend shifts

        IF (dp NE 0) THEN shft = detrendshifts(shft, r, dp)

        ; register

        FOR i = 0, dim[3]-1 DO cube[*,*,i] = doreg(cube[*,*,i], r, shft[*,*,*,i])

        ; combine shifts

        IF (flag EQ 1) THEN BEGIN
 
           sft = combineshifts(dim[1], dim[2], rt, r, sft, shft)
           rt  = r

        ENDIF

        flag = 1

     ENDFOR

     IF arg_present(shifts) THEN shifts = sft
     IF arg_present(grid)   THEN grid   = r

;;***************************************************************
;; Done
;;***************************************************************

     RETURN, cube

END
