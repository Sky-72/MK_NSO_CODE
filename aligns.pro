Function aligns, datapath, wave

filenames = file_search(datapath + wave + '_bb*.sav', count = nfiles)
restore, filenames[0]
refimage = bb_data_dstr[*,*,0]


for i = 0, nfiles - 1 do begin
  restore, filenames[i]
  bbdata = bb_data_dstr
   
  shifts = fltarr(2,10)
  bicubs = fltarr(1000,1000,10)
  
  for a = 0 , 9 do begin
    imreg = bbdata[*,*,a]
    
    initshifts = shc(refimage, imreg, /interpolate)
    shifts[0,a] = initshifts[0]
    shifts[1,a] = initshifts[1]
    
    bicub = shift_bicub(imreg, initshifts[0], initshifts[1])
    
    for x = 0, 999 do begin
      for y = 0, 999 do begin
        
         bicubs[x,y,a] = bicub[x,y]
        
      endfor
    endfor
   
    
  endfor
  save, filename =  '/home/mkatilius/IDL/shifts/6173/' + wave+ 'shifts' + string(i) + '.sav', shifts, bicubs, /verb
endfor


end

