# constbnd

### Code and data for IAU constellation boundaries

In 1930,  [the International Astronomical Union (IAU) came up with an official set of 88 constellations](https://www.newyorker.com/magazine/2017/05/01/the-irish-constellation),  dividing the sky into areas along RA/dec lines.  This code takes the dividing lines and emits a C array that simplifies the problem of figuring out which constellation corresponds to a particular RA/dec.

I really only needed this code to run once,  to create the array now stored in [`conbound.c` in the `lunar` repository](https://github.com/Bill-Gray/lunar/blob/master/conbound.c).  So the code is of historical/documentary interest only.  The list of boundary data,  of course,  remains of interest.

`constbnd.txt` contains further documentation of the history,  format,  and use of the constellation boundary data.
