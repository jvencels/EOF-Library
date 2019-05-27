# EOF-Library release notes
## Release file modified:
27-May-2019

## Current release:
0.4.0

### X.Y.Z
**X** - major releases, not back-compatible (0 = beta)
**Y** - new features, back-compatible within the same release
**Z** - patches and bug fixes

## Current issues:
* **OpenFOAM v5.0** - There are many ways to install OpenFOAM v5.0 - Ubuntu repositories, source from git development branch or source pack from OpenFOAM.org. All these versions have different changes in code responsible for MPI communication, therefore user may incounter issues compiling EOF-Library. **Fix:** try compiling OpenFOAM from https://openfoam.org/download/5-0-source/

## Changes for 0.4.0:
* Added support for OpenFOAM's multiRegion

## Changes for 0.3.0:
* Fixed Travis automatic tests.
* Travis pushes Docker images to https://hub.docker.com/u/eoflibrary

## Changes for 0.2.0:
* Added script for dynamically setting user and group ID's inside Docker.

## Changes for 0.1.0:
* First release.
