function polarwind(args)
*************************************************************************
*************************************************************************
*************************************************************************
***
*** GrADS Script Name: "polarwind"
***                     by Kyle Ahern
***                     on 2015 November 1 (Last Mod: 2015 Dec 7)
*** Updates/fixes:
***
***          -'15 Dec 4: Increasing distortion factor discovered with
***            distance from Equator, result of using scaled lat-lon 
***            projection to establish lateral and cross-lateral
***            angular distances from a given center point. 
***          +'15 Dec 7: Script updated to transform lat-lon range
***            to yield equivalent coordinates for an azimuthal
***            equidistant projection centered at a specified point
***            (r=0). Distortion is no longer increasing with absolute
***            latitude, but the polar-based information becomes
***            less meaningful with increasing radial distance. This
***            should be fine when analyzing sub-planetary scale 
***            phenomena.
***
*** Purpose: Facilitate describing motion in polar coordinates
***          assuming a given center (radius=0 @ given X, Y).
***
*** Arguments:
***
***          Four arguments can be passed to the function, none of
***          which are required for the script to run:
***          + Argument 1: longitude of r=0, OR 
***                        'q' to query the position of a mouse click
***                        that coincides with r=0, OR 
***                        any other appropriate string (not 'q') 
***                        to pass center longitude of current domain.
***          + Argument 2: latitude of r=0, OR 
***                        any appropriate string to pass center 
***                        latitude of current domain. 
***                        (not used if Argument 1 is 'q')
***          + Argument 3: name of a defined zonal flow field
***          + Argument 4: name of a defined meridional flow field
***          
***          If Arguments 1/2 are not passed, r=0 is assumed to be at 
***          the midpoint of the line from (XMIN,YMIN) to (XMAX,YMAX)
***          of the immediate domain. If Arguments 3/4 are not passed,
***          the assumed names of the zonal and meridional flow fields
***          are 'u' and 'v', respectively.
***
*** Programming Assumptions:
***
***          ~ A file has already been opened
***          ~ If passing 'q' in Argument 1, something is displayed.
***          ~ X and Y are not fixed.  
***
*** Products:
***          
***          For the dimensions set (X,Y,Z,T) when the function is called, 
***          + radius: distance(meters) from the given center at (X,Y)
***          + radang: azimuth(math) of radial unit vector at (X,Y)
***          + aziang: azimuth(math) of azimuthal unit vector at (X,Y)
***          + wndang: azimuth(math) of horizontal wind at (X,Y,Z,T)
***          + radmag: magnitude of radial wind component at (X,Y,Z,T)
***          + azimag: magnitude of azimuthal wind component at (X,Y,Z,T)
***          NOTE: passed flow fields share units with radmag and azimag.
***
*** Side-Effects:
***          
***          - The fields defined in the script take up memory,
***            so running this script with a large 4-D domain might
***            result in some undesirable effects...
***
*************************************************************************
*************************************************************************
*************************************************************************

**************
*****
* Query current set dimensional bounds
*****
**************

'q dims'
rec = sublin(result,2)
xmin = subwrd(rec,11)
xmax = subwrd(rec,13)

rec = sublin(result,3)
ymin = subwrd(rec,11)
ymax = subwrd(rec,13)

rec = sublin(result,4)
if (subwrd(rec,3) = 'fixed')
 zmin = subwrd(rec,9)
 zmax = zmin
else
 zmin = subwrd(rec,11)
 zmax = subwrd(rec,13)
endif

rec = sublin(result,5)
if (subwrd(rec,3) = 'fixed')
 tmin = subwrd(rec,9)
 tmax = tmin
else
 tmin = subwrd(rec,11)
 tmax = subwrd(rec,13)
endif

'set Z 1'
'set T 1'
'set X 'xmin
'set Y 'ymin
'define lonmin = lon'
'define latmin = lat'
'set X 'xmax
'set Y 'ymax
'define lonmax = lon'
'define latmax = lat'

**************
***** 
* Process passed args (if any); find (midlon,midlat) where r=0
*****
**************

mlon = ''
mlat = ''
if (args != '')
 if (subwrd(args,1) = 'q')
   'q pos'
   rec = sublin(result,1)
   cx = subwrd(rec,3)
   cy = subwrd(rec,4)
   'q xy2w 'cx' 'cy' '
   rec = sublin(result,1)
   mlon = subwrd(rec,3)
   mlat = subwrd(rec,6)
 else
   rc=subwrd(args,1)
   if (valnum(rc) != 0)
     mlon = subwrd(args,1)
   endif
   rc=subwrd(args,2)
   if (valnum(rc) != 0)
     mlat = subwrd(args,2)
   endif
 endif
endif

if (mlon != '')
 'define midlon = 'mlon
else
 'define midlon = (lonmin+lonmax)/2'
endif

if (mlat != '')
 'define midlat = 'mlat
else
 'define midlat = (latmin+latmax)/2'
endif

'undefine lonmin'
'undefine latmin'
'undefine lonmax'
'undefine latmax'

if (subwrd(args,3) != '')
 ufield = subwrd(args,3)
else
 ufield = 'u'
endif 

if (subwrd(args,4) != '')
 vfield = subwrd(args,4)
else
 vfield = 'v'
endif

**************
*****
* Define required constants
*****
**************

'define pi = 3.141592'
'define A = 6371393'

**************
*****
* Define 2D arrays describing X (cross-lateral) 
* and Y (lateral) angular distances relative to the center point
* throughout the domain set (xrange and yrange).
* In other words, we define xrange and yrange, which are 
* the angular distances from the origin of a flat plane orthogonal to the
* Earth's surface at the given center point
*****
**************

'set X 'xmin' 'xmax
'set Y 'ymin' 'ymax

'define adist = acos(sin(midlat*pi/180)*sin(lat*pi/180) + cos(midlat*pi/180)*cos(lat*pi/180)*cos((lon-midlon)*pi/180))'
'define kprime = adist/sin(adist)'

'define xrange = kprime*cos(lat*pi/180)*sin((lon-midlon)*pi/180)'
'define yrange = kprime*(cos(midlat*pi/180)*sin(lat*pi/180) - sin(midlat*pi/180)*cos(lat*pi/180)*cos((lon-midlon)*pi/180))'
'define radius = adist * A'

'undefine adist'
*'undefine kprime'
*'undefine midlon'
*'undefine midlat'

**************
*****
* Define arrays of azimuths for radial, tangential, and full horizontal
* wind (streamline) unit vectors. Angles are mathematical with 0
* degrees for eastward wind and increasing counter-clockwise. 
* Units are in degrees.
*****
**************

'define radang = 180 - atan2(yrange,-xrange)*180/pi'
'define aziang = 180 - atan2(xrange, yrange)*180/pi'

'set Z 'zmin' 'zmax
'set T 'tmin' 'tmax
'define wndang = 180 - atan2('vfield',-'ufield')*180/pi'

**************
*****
* TESTING
*****
**************
test=1
if (test = 1)
*'define delta = 0.00000001'
*'define xpert = delta * cos(wndang * pi / 180)'
*'define ypert = delta * sin(wndang * pi / 180)'

'define delta = 1000'
'define xpert = (delta * cos(wndang * pi / 180)) / (A * cos(lat * pi / 180))'
'define ypert = (delta * sin(wndang * pi / 180)) / A'
'define lonpert = lon + xpert'
'define latpert = lat + ypert'

'define Padist = acos(sin(midlat*pi/180)*sin(latpert*pi/180) + cos(midlat*pi/180)*cos(latpert*pi/180)*cos((lonpert-midlon)*pi/180))'
'define Pkprime = Padist/sin(Padist)'

'define Pxrange = Pkprime*cos(latpert*pi/180)*sin((lonpert-midlon)*pi/180)'
'define Pyrange = Pkprime*(cos(midlat*pi/180)*sin(latpert*pi/180) - sin(midlat*pi/180)*cos(latpert*pi/180)*cos((lonpert-midlon)*pi/180))'
'define Pradius = Padist * A'
'define projang = atan2(Pyrange-yrange, Pxrange-xrange)*180/pi'

*'tmpa=const(maskout(projang,projang),0,-u)'
*'tmpb=const(maskout(projang+360,-projang),0,-u)'
*'projang=tmpa+tmpb'
*'undefine tmpa'
*'undefine tmpb'
'proju=mag(u,v)*cos(projang*3.141592/180)'
'projv=mag(u,v)*sin(projang*3.141592/180)'
'difu=(proju/mag(proju,projv))-(u/mag(u,v))'
'difv=(projv/mag(proju,projv))-(v/mag(u,v))'

endif

**************
*****
* With azimuths of the real wind and radial/tangential unit 
* vectors, define magnitudes of the radial/tangential wind
* components via vector multiplication (note that unit vector
* magnitudes are always 1). Units are the same as u and v.
*****
**************

*'define azimag = mag('ufield','vfield')*cos((wndang - aziang)*pi/180)'
*'define radmag = mag('ufield','vfield')*cos((wndang - radang)*pi/180)'
'define azimag = mag('ufield','vfield')*cos((projang - aziang)*pi/180)'
'define radmag = mag('ufield','vfield')*cos((projang - radang)*pi/180)'
'define aziali = cos((projang - aziang)*pi/180)'
'define radali = cos((projang - radang)*pi/180)'

'undefine A'
'undefine pi'
*'undefine xrange'
*'undefine yrange'
