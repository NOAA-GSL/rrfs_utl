load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  lfile="LightningInGSI.dat"
  lfile="LightningInGSI_bufr.dat"
  numobs = fbinrecread(lfile,3,(/1/),"integer")
  lon = fbinrecread(lfile,4,(/numobs/),"float")
  lat = fbinrecread(lfile,5,(/numobs/),"float")

  print(numobs)
  print(lon(1))
  print(lat(1))

  colors=new(numobs,string)

  wks = gsn_open_wks("pdf","lightning_obs_glm")

  res                             = True
  res@gsnFrame                    = False   ; So we can draw markers
  res@gsnMaximize                 = True    ; affects ps, eps, pdf only
  res@gsnPaperOrientation         = "portrait"  ; force portrait

  res@tiMainString                = "GLM Observation Distribution"

  res@mpFillOn                    = False
  res@mpOutlineDrawOrder          = "PostDraw"
  res@mpFillDrawOrder             = "PreDraw"
  res@mpOutlineBoundarySets       = "GeophysicalAndUSStates"
  res@mpUSStateLineColor          = "Gray10"
  res@mpUSStateLineDashPattern    = 2

  res@pmTickMarkDisplayMode       = "Always"
 
  res@mpCenterLonF                = 270
  res@mpMaxLatF                   = 60           ; choose subregion
  res@mpMinLatF                   = -30
  res@mpMaxLonF                   = 355
  res@mpMinLonF                   = 200

; Draw the map (frame won't get advanced because gsnFrame was set to False).
  map = gsn_csm_map(wks,res)


; Draw markers on the plot in the lat/lon locations.
  mkres               = True
  mkres@gsMarkerIndex = 17     ; Filled circle
  mkres@gsMarkerSizeF = 0.001
  mkres@gsMarkerColor = (/"blue"/)

;  colors = (/"limegreen","orange","green","red","yellow","purple","blue",\
;             "red","brown","red2","skyblue"/)

  do i = 0, numobs-1
     colors(i)="blue"
     mkres@gsMarkerColor = colors(i)
     gsn_polymarker(wks,map,lon(i),lat(i),mkres)
  end do

  frame(wks)   ; Now advance the frame.

  res@mpCenterLonF                = 270
  res@mpMaxLatF                   = 50           ; choose subregion
  res@mpMinLatF                   = 15
  res@mpMaxLonF                   = 300
  res@mpMinLonF                   = 230
  map = gsn_csm_map(wks,res)

  mkres@gsMarkerIndex = 17     ; Filled circle
  mkres@gsMarkerSizeF = 0.001
  mkres@gsMarkerColor = (/"blue"/)
  do i = 0, numobs-1
     colors(i)="blue"
     mkres@gsMarkerColor = colors(i)
     gsn_polymarker(wks,map,lon(i),lat(i),mkres)
  end do
  frame(wks)   ; Now advance the frame.
end
