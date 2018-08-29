#!/bin/csh
pscoast -R130/230/-25/60 -JM6i -B20/20 -K -W0.5p >! test.ps
pscontour pac.tzavg.dat -R -JM -K -O -Cpac.cpt -I >> test.ps
pscoast -R -JM -O -G140 >> test.ps
#psxy -R -JM -K -O -Sx0.10 -G0 data.hawaii.bp -: >> test.ps
#psxy -R -JM -K -O -St0.15 -G0 kip >> test.ps
#psxy -R -JM -O -Si0.15 -G0 events >> test.ps
#psscale -Csaw.cpt -D0.3i/1.0i/1.75i/1.5i -O -B0.5::/:"TZ Avg dVs (%)":
 >> test.ps
#pscoast -R -JM -W0.5 -O  >> test.ps
#pscontour 2degmelt.xyz -W0.5p -C../CPT/melt.cpt -R -JM -O -: >> test.ps
