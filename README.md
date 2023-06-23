# rrfs_utl
Utilities for RRFS applications

# Content:

This repository includes many utilites for RRFS application:
 - cloud analysis for FV3
 - radar reflectivity to temperature tendency
 - preprocess NSSL radar reflectivity mosaic
 - preprocess lightning data (bufr and netcdf)
 - preprocess METAR cloud observations
 - preprocess NASA LaRC cloud products

## Prerequisites

This package requires the following third party libraries:
- [Jasper](http://www.ece.uvic.ca/~mdadams/jasper/)
- [libpng](http://www.libpng.org/pub/png/libpng.html)
- [zlib](http://www.zlib.net/)

This package requires the folling NCEPLIBS libraries:
- [NCEPLIBS-bacio](https://github.com/NOAA-EMC/NCEPLIBS-bacio)
- [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp)
- [NCEPLIBS-ip](https://github.com/NOAA-EMC/NCEPLIBS-ip)
- [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2)
- [NCEPLIBS-g2tmpl](https://github.com/NOAA-EMC/NCEPLIBS-g2tmpl)
- [NCEPLIBS-w3nco](https://github.com/NOAA-EMC/NCEPLIBS-w3nco)
- [NCEPLIBS-w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc)
- [NCEPLIBS-bufr](https://github.com/NOAA-EMC/NCEPLIBS-bufr)
- [NCEPLIBS-wrf_io](https://github.com/NOAA-EMC/NCEPLIBS-wrf_io)

This package requires the folling GSI packages:
- [GSI-ncdiag](https://github.com/NOAA-EMC/GSI-ncdiag)
- [GSI](https://github.com/NOAA-EMC/GSI)

## Installing

```
mkdir build
cd build
cmake PATH2ROOT
make
```

## Disclaimer

```
The United States Department of Commerce (DOC) GitHub project code is
provided on an "as is" basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

