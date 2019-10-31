PRO EHANDLER, ev
  widget_control, /DESTROY, ev.TOP
END

PRO anim, movie, TRACK=track, CYCLE=cycle, GLOBAL_SCALE=global_scale, $
          _EXTRA=extra_keywords

  movieSize = size(movie)
  IF (movieSize[0] NE 3) THEN BEGIN
    print, 'Not a three-dimensional array'
    return
  ENDIF
  xwindowSize = movieSize[1]
  ywindowSize = movieSize[2]
  Ntime = movieSize[3]

  baseWidget    = widget_base(TITLE='Animate plot')
  animateWidget = cw_animate(baseWidget, xwindowSize, ywindowSize, Ntime, $
			     TRACK=keyword_set(TRACK), $
			     CYCLE=keyword_set(CYCLE))
  widget_control, /REALIZE, baseWidget

  IF (keyword_set(GLOBAL_SCALE)) THEN BEGIN
    bmovie = bytscl(movie)
    FOR n_iter=0, Ntime-1 DO $
     cw_animate_load, animateWidget, FRAME=n_iter, IMAGE=bmovie[*, *, n_iter]
  ENDIF ELSE BEGIN
    FOR n_iter=0, Ntime-1 DO $
     cw_animate_load, animateWidget, FRAME=n_iter, $
     IMAGE=bytscl(movie[*, *, n_iter])
  ENDELSE

  cw_animate_run, animateWidget
  xmanager, "Xanim", baseWidget, EVENT_HANDLER="EHANDLER"
END
