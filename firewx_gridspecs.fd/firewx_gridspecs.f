         program latlon2llc
C
      INTEGER       GRD(18,2)
      CHARACTER*6   DOMAIN
      INTEGER       IINDEX(4)
      INTEGER       JINDEX(4)
      REAL          ALAT(4)
      REAL          ELON(4)
C
      INTEGER       GRIBTYPE
C
      DATA  GRD /0, 255, 3,4289,2753, 20192, -121554,  8,  -95000,
     &   1270, 1270, 0, 64, 0, 25000, 25000, 0, 0, 
     &           0, 255, 5,3297,2209, 40530, -178571,  8,  -150000,
     &   1488, 1488, 0, 64, 0, 0, 0, 0, 0/
C
       READ(45,100) DOMAIN
100    FORMAT(A6)
       READ(45,99) GRIBTYPE
 99    FORMAT(I1)
       READ(11,101) LATSW,LSW,LATSE,LSE,LATNW,LNW,LATNE,LNE
101    format(4(1X,I5,1x,I6,1x))
102    format(4(1X,I5,1x,I7,1x))
       LONSW=LSW-360000
       LONSE=LSE-360000
       LONNW=LNW-360000
       LONNE=LNE-360000
       WRITE(6,102) LATSW,LONSW,LATSE,LONSE,LATNW,LONNW,LATNE,LONNE
C
C        K=1: CONUS GRID; K=2: ALASKA GRID

       IF (DOMAIN .EQ. 'CONUS') THEN
         K=1
       ELSE
         K=2
       ENDIF
C
       IF(K.EQ.1) THEN 

C   fire weather domain is over CONUS
C
       PRINT *,'Fire Weather region= ',DOMAIN

         ignum=90

         alat1 = real(grd(6,k))/1000.
         elon1 = 360.0 + (real(grd(7,k))/1000.)
         dx = real(grd(10,k))
         dx = real(grd(11,k))
         elonv = 360. + (real(grd(9,k))/1000.)
         alatan =  real(grd(15,k))/1000.
         print *,alat1,elon1,dx,elonv,alatan

         alat(1)=real(LATSW)/1000.
         elon(1)=360.0 + real(LONSW)/1000.
         call W3FB11(ALAT(1),ELON(1),ALAT1,ELON1,DX,ELONV,ALATAN,
     1      XISW,XJSW)
         isw=int(xisw)
         jsw=int(xjsw)
         iindex(1)=isw
         jindex(1)=jsw
         print *,'sw pt lat=',alat(1),' lon=',elon(1),' at i,j= ',
     1      xisw,xjsw,isw,jsw
         alat(2)=real(LATSE)/1000.
         elon(2)=360.0 + real(LONSE)/1000.
         call W3FB11(ALAT(2),ELON(2),ALAT1,ELON1,DX,ELONV,ALATAN,
     1      XISE,XJSE)
         ise=int(xise)
         jse=int(xjse)
         iindex(2)=ise
         jindex(2)=jse
         print *,'se pt lat=',alat(2),' lon=',elon(2),' at i,j= ',
     1      xise,xjse,ise,jse
         alat(3)=real(LATNW)/1000.
         elon(3)=360.0 + real(LONNW)/1000.
         call W3FB11(ALAT(3),ELON(3),ALAT1,ELON1,DX,ELONV,ALATAN,
     1      XINW,XJNW)
         inw=int(xinw)
         jnw=int(xjnw)
         iindex(3)=inw
         jindex(3)=jnw
         print *,'nw pt lat=',alat(3),' lon=',elon(3),' at i,j= ',
     1      xinw,xjnw,inw,jnw
         alat(4)=real(LATNE)/1000.
         elon(4)=360.0 + real(LONNE)/1000.
         call W3FB11(ALAT(4),ELON(4),ALAT1,ELON1,DX,ELONV,ALATAN,
     1      XINE,XJNE)
         ine=int(xine)
         jne=int(xjne)
         iindex(4)=ine
         jindex(4)=jne
         print *,'ne pt lat=',alat(4),' lon=',elon(4),' at i,j= ',
     1      xine,xjne,ine,jne
   
           xlatsw=99999.
           xlonsw=99999.
           xlatne=-99999.
           xlonne=-99999.
           do j=1,4
             xlatsw=min(alat(j),xlatsw)
             xlonsw=min(elon(j),xlonsw)
             xlatne=max(alat(j),xlatne)
             xlonne=max(elon(j),xlonne)
           enddo
           print *,xlatsw,xlonsw,xlatne,xlonne
           call W3FB11(xlatsw,xlonsw,ALAT1,ELON1,DX,ELONV,ALATAN,
     1        XSW,YSW)
           call W3FB11(xlatne,xlonne,ALAT1,ELON1,DX,ELONV,ALATAN,
     1        XNE,YNE)
           print *,XSW,YSW,XNE,YNE 
           isw=int(xsw)
           jsw=int(ysw)  
           ine=int(xne)
           jne=int(yne)
           idim=ine-isw+1
           jdim=jne-jsw+1
           print *,isw,jsw,ine,jne,idim,jdim
C
C    get lat/lon on sw corner on master CONUS grid
C
         fxsw=real(isw)
         fysw=real(jsw)
         print *,fxsw,fysw

         CALL W3FB12(FXSW,FYSW,ALAT1,ELON1,DX,ELONV,ALATAN,
     1      ALATOUT,ELONOUT,IERR)
         print *,'lat/lon of SW corner on master grid= ',ALATOUT,ELONOUT
C
C  write out block for copygb to create fire weather output grid
C
         latstartfw=nint(alatout*1000.)
         lonstartfw=nint(elonout*1000.)
         lonv=360000+int(grd(9,k))
         idx=grd(10,k)
         idy=grd(11,k)
         ilatan=grd(15,k)
C
        if(gribtype.eq.1) then 
         open(110,file='copygb_gridnavfw.txt',form='formatted'
     1       ,status='unknown')
         write(110,1000)idim,jdim,latstartfw,lonstartfw,lonv,
     1      idx,idy,ilatan,ilatan
1000     format('255 3 ',2(I4,1x),I6,1x,I7,1x,'8 ',I7,1x,
     1         2(I6,1x),'0 64',2(1x,I6))
         close(110)
        else
         open(110,file='copygb_gridnavfw.txt',form='formatted'
     1       ,status='unknown')
         open(120,file='gridnavfw.ijdims.txt',form='formatted'
     1       ,status='unknown')
         lonv2=lonv/1000
         ilatan2=ilatan/1000
         write(110,1100)lonv2,ilatan2,ilatan2,elonout,idim,
     1         idx,alatout,jdim,idy
1100     format('lambert:',i3,':',i2,':',i2,1x,f7.3,':',i3,':',i4,
     1        1x,f6.3,':',i3,':',i4)
         write(120,1200)idim,jdim
1200     format(2(i4,1x))
         close(110)
         close(120)
        endif
       ELSE

c   fire weather domain is over Alaska

         ignum=92

         PRINT *,'Fire Weather region= ',DOMAIN, K

         alat1 = real(grd(6,k))/1000.
         elon1 = 360.0 + (real(grd(7,k))/1000.)
         dx = real(grd(10,k))
         alonv = 360.0 + (real(grd(9,k))/1000.)
         print *,alat1,elon1,dx,alonv

         alat=real(LATSW)/1000.
         elon=360.0+(real(LONSW)/1000.)
         call W3FB06(ALAT,ELON,ALAT1,ELON1,DX,ALONV,XISW,XJSW)
         isw=int(xisw)
         jsw=int(xjsw)
         iindex(1)=isw
         jindex(1)=jsw
         print *,'sw grid ',ignum,'lat=',alat,' lon=',elon,' at i,j= ',
     1      xisw,xjsw,isw,jsw

         alat=real(LATSE)/1000.
         elon=360.0+(real(LONSE)/1000.)
         call W3FB06(ALAT,ELON,ALAT1,ELON1,DX,ALONV,XISE,XJSE)
         ise=int(xise)
         jse=int(xjse)
         iindex(2)=ise
         jindex(2)=jse
         print *,'se grid ',ignum,'lat=',alat,' lon=',elon,' at i,j= ',
     1      xise,xjse,ise,jse

         alat=real(LATNW)/1000.
         elon=360.0+(real(LONNW)/1000.)
         call W3FB06(ALAT,ELON,ALAT1,ELON1,DX,ALONV,XINW,XJNW)
         inw=int(xinw)
         jnw=int(xjnw)
         iindex(3)=inw
         jindex(3)=jnw
         print *,'nw grid ',ignum,'lat=',alat,' lon=',elon,' at i,j= ',
     1      xinw,xjnw,inw,jnw

         alat=real(LATNE)/1000.
         elon=360.0+(real(LONNE)/1000.)
         call W3FB06(ALAT,ELON,ALAT1,ELON1,DX,ALONV,XINE,XJNE)
         ine=int(xine)
         jne=int(xjne)
         iindex(4)=ine
         jindex(4)=jne
         print *,'ne grid ',ignum,'lat=',alat,' lon=',elon,' at i,j= ',
     1      xine,xjne,ine,jne

         iptmin=99999
         jptmin=99999
         iptmax=-99999
         jptmax=-99999
         do l=1,4
           iptmin=min(iindex(l),iptmin)
           iptmax=max(iindex(l),iptmax)
           jptmin=min(jindex(l),jptmin)
           jptmax=max(jindex(l),jptmax)
         enddo
         print *,iptmin,iptmax,jptmin,jptmax
C
         iptsw=iptmin-100
         jptsw=jptmin-100
         iptne=iptmax+100
         jptne=jptmax+100
         idim=iptne-iptsw+1
         jdim=jptne-jptsw+1
         print *,'corners on target grid=', iptsw,jptsw,iptne,jptne
         print *,'fire weather grid dimensions= ',idim,jdim
C
C    get lat/lon on sw corner on master CONUS grid
C
         xiptsw=real(iptsw)
         yjptsw=real(jptsw)
         xiptne=real(iptne)
         yjptne=real(jptne)

         print *,xiptsw,yjptsw,alat1,elon1,dx,alonv
         call W3FB07(xiptsw,yjptsw,alat1,elon1,dx,alonv,alatout,elonout)
         print *,'lat/lon of SW corner on master grid= ',ALATOUT,ELONOUT
         call W3FB07(xiptne,yjptne,alat1,elon1,dx,alonv,alatoutne,
     1       elonoutne)
         print *,'lat/lon of NE corner on master grid= ',ALATOUTNE,
     1       ELONOUTNE

         latstartfw=nint(alatout*1000.)
         lonstartfw=nint(elonout*1000.)
         lonv=360000+int(grd(9,k))
         idx=grd(10,k)
         idy=grd(11,k)
         ilatan=grd(15,k)
C
        if(gribtype.eq.1) then
         open(110,file='copygb_gridnavfw.txt',form='formatted'
     1       ,status='unknown')
         write(110,1001)idim,jdim,latstartfw,lonstartfw,lonv,
     1      idx,idy
1001     format('255 5 ',2(I4,1x),I6,1x,I7,1x,'8 ',I7,1x,
     1         2(I6,1x),'0 64 0 0')
         close(110)
        else
         open(110,file='copygb_gridnavfw.txt',form='formatted'
     1       ,status='unknown')
         open(120,file='gridnavfw.ijdims.txt',form='formatted'
     1       ,status='unknown')
         lonv2=lonv/1000
         ilatan2=60
         write(110,1101)lonv2,ilatan2,elonout,idim,idx,alatout,
     1      jdim,idy
1101     format('nps:',i3,':',i2,1x,f7.3,':',i3,':',i4,1x,f6.3,':',
     1         i3,':',i4)
         write(120,1201)idim,jdim
1201     format(2(i4,1x))
         close(110)
         close(120)
        endif

       ENDIF
C
         stop
         end
