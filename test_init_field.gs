* This GrADS script reads in a strictly formatted namelist
* block for multi TC testing in FV3 forecasts, and produces
* the initial surface fields associated with that namelist.
*
* Important fields produced:
*   - psfc : surface pressure
*   - u : zonal wind
*   - v : meridional wind

'reinit'

* IMPORTANT NOTE:
* Change the line below so that SOME valid
* file is opened for GrADS to work with.
* The resolution of the fields produced by this
* test will be the same as that of the file opened.
_filename=""
if (_filename = "")
  say "  Error: No valid data file has been specified."
  say "         Please check the script's top comments,"
  say "         and specify a data file to open."
  return
endif
'open '_filename

_PI = 3.1415926

* open/parse namelist
* Assumptions: The namelist has the following structure:
*
* &test_case_nml
*       num_vortex = _N
*       lon = val.0, val.1, val.2, val._(N-1)
*       lat = val.0, val.1, val.2, val._(N-1)
*       dp = val.0, val.1, val.2, val._(N-1)
*       rsize = val.0, val.1, val.2, val._(N-1)
*       q_0 = val
*       T_0 = val
*       p_0 = val
*/
*
* ... with lists of values delimited with a comma and space.
* Spaces should be between variable names and the equals sign, and
* each variable should appear in the above order. All variables must
* be specified.

file="namelist.nml"
rc=read(file)

rc  = read(file);
lin = sublin(rc,2);
_N  = subwrd(lin,3);

rc = read(file); lon_lin    = sublin(rc,2);
rc = read(file); lat_lin    = sublin(rc,2);
rc = read(file); dp_lin     = sublin(rc,2);
rc = read(file); rsize_lin  = sublin(rc,2);
i = 0; wrdoffset = 3;
while (i < _N)
  _lon.i   = subwrd(lon_lin,   wrdoffset + i);
  _lat.i   = subwrd(lat_lin,   wrdoffset + i);
  _dp.i    = subwrd(dp_lin,    wrdoffset + i);
  _rsize.i = subwrd(rsize_lin, wrdoffset + i);
  if (i != _N - 1)
    len = strlen(_lon.i);   _lon.i   = substr(_lon.i,1,len-1);
    len = strlen(_lat.i);   _lat.i   = substr(_lat.i,1,len-1);
    len = strlen(_dp.i);    _dp.i    = substr(_dp.i,1,len-1);
    len = strlen(_rsize.i); _rsize.i = substr(_rsize.i,1,len-1);
  endif
  i=i+1
endwhile

rc = read(file); lin = sublin(rc,2); _q0 = subwrd(lin,3);
rc = read(file); lin = sublin(rc,2); _T0 = subwrd(lin,3);
rc = read(file); lin = sublin(rc,2); _p0 = subwrd(lin,3);

********************************************************

'u = lev * 0'
'v = lev * 0'

* i=0
*while (i<_N)
*  'q w2xy '_lon.i' '_lat.i
*  vx = subwrd(result,3)
*  vy = subwrd(result,6)
*  'draw mark 6 'vx' 'vy' 0.08'
*  i = i + 1
*endwhile

'c'
'denom = 0*lev'

i=0; while (i<_N)
  'run polarwind.gs '_lon.i' '_lat.i' u v'
  'r'i' = radius/1000.0'
  'denom = denom + (1.0 / (r'i' * r'i'))'
  'w'i' = 1.0 / (r'i' * r'i')'
  i=i+1
endwhile
i=0; while (i<_N)
  'w'i' = w'i' / denom'
  _rsize.i = _rsize.i / 1000.0
  i=i+1
endwhile

'psfc = lev*0'
i=0; while (i<_N)
  'p'i' = '_p0' - '_dp.i' * EXP(-POW(r'i'/'_rsize.i',3/2))'
  'psfc = psfc + (w'i' * p'i')'
  i=i+1
endwhile

* temperature/humidity have no horizontal structure at sfc
'Tsfc = '_T0'*lev/lev'
'qsfc = '_q0'*lev/lev'

'u = lev*0'
'v = lev*0'
i=0; while (i<_N)
  'f = 2 * 7.292 * pow(10,-5) * sin('_lat.i'*'_PI'/180.0)'
  'terma = - f * r'i' * 1000.0 / 2.0'
  'termb = f * f * (r'i' * 1000.0) * (r'i' * 1000.0) / 4.0'
  'termc = (3/2) * Tsfc * 287.15 * POW(r'i'/'_rsize.i',3/2)'
  'termd = 1.0 - ( ('_p0'/'_dp.i') * EXP(POW(r'i'/'_rsize.i',3/2)) )'
  'vgrad'i' = terma + SQRT(termb - (termc / termd))'

* convert vgrad to u, v components using DCMIP-2016's method
  'd1 = sin('_lat.i'*'_PI'/180.0)*cos(lat*'_PI'/180.0) - cos('_lat.i'*'_PI'/180.0)*sin(lat*'_PI'/180.0)*cos((lon-'_lon.i')*('_PI'/180.0))'
  'd2 = cos('_lat.i'*'_PI'/180.0)*sin((lon-'_lon.i')*('_PI'/180.0))'
  'd  = sqrt((d1*d1) + (d2*d2))'
  'dtmp = const(maskout(d,d-pow(10,-25)),0,-u)'

  'u'i' = vgrad'i' * d1 / dtmp'
  'v'i' = vgrad'i' * d2 / dtmp'

  'u = u + (w'i' * u'i')'
  'v = v + (w'i' * v'i')'

  i=i+1
endwhile

'u = const(maskout(u,abs(u)-pow(10,-12)),0,-u)'
'v = const(maskout(v,abs(v)-pow(10,-12)),0,-u)'

return
