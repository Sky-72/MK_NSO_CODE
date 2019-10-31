Function alignall, datapathfe, datapathha

  filenamefe = file_search(datapathfe + '6173_bb*.sav', count = nfilesfe)
  filenameha = file_search(datapathha + '6563_bb*.sav', count = nfilesha)

  restore, filenamefe[0]
  refimage = bb_data_dstr[*,*,0]


  for i = 0, nfilesha - 1 do begin
    restore, filenameha[i]
    bbdata = bb_data_dstr

    shifts = fltarr(2,14)
    bicubs = fltarr(1000,1000,14)

    for a = 0 , 13 do begin
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
    save, filename =  '/home/mkatilius/IDL/shifts/6563_6173_CO/' + 'Ha_FeI_shifts' + string(i) + '.sav', shifts, bicubs, /verb
  endfor


end

