load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  bk=addfile("geo_em.d01.nc.cdf","r")
  rlon = bk->XLONG_M(0,:,:)
  rlat = bk->XLAT_M(0,:,:)

  nlon=dimsizes(rlon(1,:))
  nlat=dimsizes(rlon(:,1))
  print(nlon)
  print(nlat)

  lfile="LightningInGSI.dat"
  lfile="LightningInGSI_bufr.dat"
  dims= fbinrecread(lfile,0,(/6/),"integer")
  print(dims)
  numobs=dims(3)
  lgtbufr= fbinrecread(lfile,1,(/numobs,3/),"float")
;  print(lgtbufr)
;  exit
  lgt= fbinrecread(lfile,2,(/nlat,nlon/),"float")

  xwks = gsn_open_wks("pdf","HRRR_lightning_glm")
  gsn_define_colormap(xwks,"cosam12")
 
  resources                  =  True              ; plot mods desired
  resources@sfXArray         = rlon
  resources@sfYArray         = rlat


  resources@cnMonoLineColor    = False
  resources@tiMainString       ="HRRR_lightning_glm"

  resources@mpProjection       = "LambertConformal" ; choose projection
  resources@mpLimitMode        = "Corners"          ; choose range of map
  resources@mpLeftCornerLatF      = 20.0
  resources@mpLeftCornerLonF      = -125.
  resources@mpRightCornerLatF     = 46.0
  resources@mpRightCornerLonF     = -50.0

  resources@mpLambertParallel1F   = 38.5        ; ncl adds from grib file
  resources@mpLambertParallel2F   = 38.5       ; ncl adds from grib file
  resources@mpLambertMeridianF    = -97.0      ; ncl adds from grib file
  resources@pmTickMarkDisplayMode = "Always"
  resources@mpFillOn              =  False             ; turn off map fill
  resources@mpOutlineDrawOrder    = "PostDraw"         ; continent outline last
  resources@mpOutlineBoundarySets = "GeophysicalAndUSStates" ; state boundaries
  resources@cnFillOn              = True               ; turn on color
  resources@mpGridLineDashPattern = 2                  ; lat/lon lines as dashed
  resources@cnLinesOn             = False              ; turn off contour lines
  resources@cnLineLabelsOn        = False   
  resources@cnConstFLabelPerimOn  = False               ; 
  resources@mpPerimOn = True
  
  resources@pmLabelBarDisplayMode = "Always" ; Turn on label bar.
  resources@lbPerimOn             = False    ; Turn off perimeter on label bar.

  resources@cnFillColors          = (/-1,3,4,5,6,7,8,9,10,11,12/)
  resources@cnFillDrawOrder       = "PreDraw"         ; make sure fill map on top
; 
  resources@cnLevelSelectionMode = "ExplicitLevels"     ;set explicit contour levels

  resources@vpXF       = 0.1
  resources@vpYF       = 0.8
  resources@vpWidthF   = 0.6
  resources@vpHeightF  = 0.6

  resources@cnLevels    = (/1.0,2.0,3.0,4.0,5.0,8.0,10.0,20.0/)	; set levels
  plot = gsn_contour_map(xwks,lgt, resources)

;  frame(xwks)   ;
  delete(plot)
 
  resources@gsnFrame    = False   ; So we can draw markers
  map = gsn_csm_map(xwks,resources)

  colors=new(numobs,string)

; Draw markers on the plot in the lat/lon locations.
  mkres               = True
  mkres@gsMarkerIndex = 17     ; Filled circle
  mkres@gsMarkerSizeF = 0.001
  mkres@gsMarkerColor = (/"blue"/)
  do i = 0, numobs-1
     ii=toint(lgtbufr(i,0))
     jj=toint(lgtbufr(i,1))
     lgt=lgtbufr(i,2)
     llon=rlon(jj,ii)
     llat=rlat(jj,ii)
     colors(i)="blue"
     mkres@gsMarkerColor = colors(i)
     gsn_polymarker(xwks,map,llon,llat,mkres)
  end do
;  frame(xwks)   ; Now advance the frame.

  delete(resources)
end

