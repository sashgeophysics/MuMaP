#!/bin/bash
rm coral_sea.eps
pscoast -R149.9/166.9/-29.34/-11.57 -JM6i -B2/2 -K -W0.5p >> coral_sea.eps
#pscontour b18T1800.dat  -R -JM -K -O -Cpac.cpt -I >> b18t1800.eps
#psxy -R -JM -K -O -Sc0.10 -G0 b18T1800.dat >> b18t1800.eps
