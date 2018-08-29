type fit
   real(sp)::p1,p2,p3,p4,p5
end type fit


function temperature2velocity(reg)
  type(mantle),intent(in)::reg
  type(fit),dimension(3)::f
  if(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.024698741765895_sp
     f(1)%p2 = 0.0131335_sp
     f(1)%p3 = -0.1338732_sp
     f(1)%p4 = 0.067342423105036_sp
     f(1)%p5 = -0.147207542559905_sp
     f(1)%p6 = 4.62002695638021_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.040544552_sp
     f(2)%p2 = 0.020999869_sp
     f(2)%p3 = -0.216600713_sp
     f(2)%p4 = 0.09565541_sp
     f(2)%p5 = -0.164837217_sp
     f(2)%p6 = 8.553836816_sp

     ! Enter coefficients for polynomial fit for rho                          
     f(3)%p1 = 0.011475026_sp
     f(3)%p2 = 0.004669748837907_sp
     f(3)%p3 = 0.058832244_sp
     f(3)%p4 = 0.025607662224470_sp
     f(3)%p5 = 0.019420457648454_sp
     f(3)%p6 = 0.517247507_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.05

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.033569124758877_sp
     f(1)%p2 = 0.008705668391063_sp
     f(1)%p3 = -0.157966508146254_sp
     f(1)%p4 = 0.073347975005621_sp
     f(1)%p5 = -0.121804037718849_sp
     f(1)%p6 = 4.62650502278646_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.052163222_sp
     f(2)%p2 = 0.013556168_sp
     f(2)%p3 = -0.246302606_sp
     f(2)%p4 = 0.107638253_sp
     f(2)%p5 = -0.134896583_sp
     f(2)%p6 = 8.573073743_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.013821107_sp
     f(3)%p2 = 0.002755044_sp
     f(3)%p3 = -0.064526981_sp
     f(3)%p4 = 0.028471266_sp
     f(3)%p5 = -0.01516631_sp
     f(3)%p6 = 3.529530537_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.10

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.03930151619814_sp
     f(1)%p2 = 0.000733747938222_sp
     f(1)%p3 = -0.172109249674277_sp
     f(1)%p4 = 0.089147281545208_sp
     f(1)%p5 = -0.103497251265887_sp
     f(1)%p6 = 4.62765343098958_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.060134205_sp
     f(2)%p2 = 0.002456594_sp
     f(2)%p3 = -0.264999403_sp
     f(2)%p4 = 0.128898322_sp
     f(2)%p5 = -0.112109564_sp
     f(2)%p6 = 8.58734111_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.01554027_sp
     f(3)%p2 = 0.000598703_sp
     f(3)%p3 = -0.068267948_sp
     f(3)%p4 = 0.031778066_sp
     f(3)%p5 = -0.012176441_sp
     f(3)%p6 = 3.541306475_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.15

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.041504524257226_sp
     f(1)%p2 = -0.003784023091918_sp
     f(1)%p3 = -0.174075016461201_sp
     f(1)%p4 = 0.092591525528565_sp
     f(1)%p5 = -0.09649452055014_sp
     f(1)%p6 = 4.6299030859375_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.062904446_sp
     f(2)%p2 = -0.003196849_sp
     f(2)%p3 = -0.265522935_sp
     f(2)%p4 = 0.130509996_sp
     f(2)%p5 = -0.10707011_sp
     f(2)%p6 = 8.603064199_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.015918749_sp
     f(3)%p2 = -0.000385164_sp
     f(3)%p3 = -0.067206599_sp
     f(3)%p4 = 0.031094597_sp
     f(3)%p5 = -0.013336514_sp
     f(3)%p6 = 3.553276846_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.18

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.03701197717091_sp
     f(1)%p2 = -0.004985774478932_sp
     f(1)%p3 = -0.152249270780995_sp
     f(1)%p4 = 0.088063432860998_sp
     f(1)%p5 = -0.107475130879295_sp
     f(1)%p6 = 4.64967819010416_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.057696289_sp
     f(2)%p2 = -0.005732567_sp
     f(2)%p3 = -0.239714163_sp
     f(2)%p4 = 0.126204053_sp
     f(2)%p5 = -0.112272617_sp
     f(2)%p6 = 8.637712256_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.01516792_sp
     f(3)%p2 = -0.001184303_sp
     f(3)%p3 = -0.063388024_sp
     f(3)%p4 = 0.030633568_sp
     f(3)%p5 = -0.010989683_sp
     f(3)%p6 = 3.571851087_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.20

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.042266652130644_sp
     f(1)%p2 = -0.007738160893687_sp
     f(1)%p3 = -0.170863780380063_sp
     f(1)%p4 = 0.095185486673414_sp
     f(1)%p5 = -0.093701649784071_sp
     f(1)%p6 = 4.63214535807291_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.063007713_sp
     f(2)%p2 = -0.009165374_sp
     f(2)%p3 = -0.256370948_sp
     f(2)%p4 = 0.134439251_sp
     f(2)%p5 = -0.109796588_sp
     f(2)%p6 = 8.61650792_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.015604386_sp
     f(3)%p2 = -0.001806935_sp
     f(3)%p3 = -0.063643435_sp
     f(3)%p4 = 0.032010412_sp
     f(3)%p5 = -0.016413248_sp
     f(3)%p6 = 3.564093919_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.25

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.042605700609891_sp
     f(1)%p2 = -0.010975750862198_sp
     f(1)%p3 = -0.166442746006029_sp
     f(1)%p4 = 0.09578245346941_sp
     f(1)%p5 = -0.091548686331383_sp
     f(1)%p6 = 4.63598026367187_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.063153832_sp
     f(2)%p2 = -0.014878907_sp
     f(2)%p3 = -0.24768418_sp
     f(2)%p4 = 0.137851388_sp
     f(2)%p5 = -0.111710574_sp
     f(2)%p6 = 8.630769515_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.015510524_sp
     f(3)%p2 = -0.003407902_sp
     f(3)%p3 = -0.060888157_sp
     f(3)%p4 = 0.033516738_sp
     f(3)%p5 = -0.018758124_sp
     f(3)%p6 = 3.574602526_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.30

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.042292854144701_sp
     f(1)%p2 = -0.013817279577148_sp
     f(1)%p3 = -0.159933474_sp
     f(1)%p4 = 0.095241912_sp
     f(1)%p5 = -0.090852791_sp
     f(1)%p6 = 4.641026403_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.063174852_sp
     f(2)%p2 = -0.020310125_sp
     f(2)%p3 = -0.238827058_sp
     f(2)%p4 = 0.140614542_sp
     f(2)%p5 = -0.113439841_sp
     f(2)%p6 = 8.645885238_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.015625988_sp
     f(3)%p2 = -0.005036731_sp
     f(3)%p3 = -0.05892288_sp
     f(3)%p4 = 0.035171653_sp
     f(3)%p5 = -0.020385918_sp
     f(3)%p6 = 3.584999427_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.35

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.041823706_sp
     f(1)%p2 = -0.016224804_sp
     f(1)%p3 = -0.152532655_sp
     f(1)%p4 = 0.093843508_sp
     f(1)%p5 = -0.091142201_sp
     f(1)%p6 = 4.646993145_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.063699171_sp
     f(2)%p2 = -0.02481963_sp
     f(2)%p3 = -0.230992297_sp
     f(2)%p4 = 0.141766303_sp
     f(2)%p5 = -0.114752384_sp
     f(2)%p6 = 8.661962962_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.016203864_sp
     f(3)%p2 = -0.006240088_sp
     f(3)%p3 = -0.058152775_sp
     f(3)%p4 = 0.036135159_sp
     f(3)%p5 = -0.021261536_sp
     f(3)%p6 = 3.595496185_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.40

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.039052482_sp
     f(1)%p2 = -0.017000918_sp
     f(1)%p3 = -0.138160051_sp
     f(1)%p4 = 0.08840004_sp
     f(1)%p5 = -0.094660502_sp
     f(1)%p6 = 4.654375313_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.060833504_sp
     f(2)%p2 = -0.026517515_sp
     f(2)%p3 = -0.212888474_sp
     f(2)%p4 = 0.136257738_sp
     f(2)%p5 = -0.120284653_sp
     f(2)%p6 = 8.679960553_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.01599011_sp
     f(3)%p2 = -0.006598754_sp
     f(3)%p3 = -0.054844046_sp
     f(3)%p4 = 0.035298378_sp
     f(3)%p5 = -0.023079493_sp
     f(3)%p6 = 3.606348245_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.45

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.032990047_sp
     f(1)%p2 = -0.01636712_sp
     f(1)%p3 = -0.113890055_sp
     f(1)%p4 = 0.079072426_sp
     f(1)%p5 = -0.103116581_sp
     f(1)%p6 = 4.663306257_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.052347668_sp
     f(2)%p2 = -0.025767389_sp
     f(2)%p3 = -0.178078371_sp
     f(2)%p4 = 0.124206449_sp
     f(2)%p5 = -0.13346647_sp
     f(2)%p6 = 8.700166868_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.014048297_sp
     f(3)%p2 = -0.006376013_sp
     f(3)%p3 = -0.046411589_sp
     f(3)%p4 = 0.032982718_sp
     f(3)%p5 = -0.027084709_sp
     f(3)%p6 = 3.617598636_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.50

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.027202412_sp
     f(1)%p2 = -0.016070079_sp
     f(1)%p3 = -0.090614349_sp
     f(1)%p4 = 0.070780745_sp
     f(1)%p5 = -0.111141216_sp
     f(1)%p6 = 4.672457861_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.04358884_sp
     f(2)%p2 = -0.025580195_sp
     f(2)%p3 = -0.142748807_sp
     f(2)%p4 = 0.113819649_sp
     f(2)%p5 = -0.146635043_sp
     f(2)%p6 = 8.720566146_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.011688453_sp
     f(3)%p2 = -0.006423624_sp
     f(3)%p3 = -0.0368329_sp
     f(3)%p4 = 0.031398812_sp
     f(3)%p5 = -0.031425453_sp
     f(3)%p6 = 3.628757852_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.55

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.021946754_sp
     f(1)%p2 = -0.016081429_sp
     f(1)%p3 = -0.069129207_sp
     f(1)%p4 = 0.063463503_sp
     f(1)%p5 = -0.118348314_sp
     f(1)%p6 = 4.681843763_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.035078793_sp
     f(2)%p2 = -0.025990608_sp
     f(2)%p3 = -0.108678265_sp
     f(2)%p4 = 0.10509154_sp
     f(2)%p5 = -0.158616193_sp
     f(2)%p6 = 8.741274014_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.009175135_sp
     f(3)%p2 = -0.006795698_sp
     f(3)%p3 = -0.02708548_sp
     f(3)%p4 = 0.030589567_sp
     f(3)%p5 = -0.035329749_sp
     f(3)%p6 = 3.639916449_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.60

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.016558326_sp
     f(1)%p2 = -0.014049376_sp
     f(1)%p3 = -0.049013716_sp
     f(1)%p4 = 0.051687846_sp
     f(1)%p5 = -0.122302592_sp
     f(1)%p6 = 4.691804854_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.026136988_sp
     f(2)%p2 = -0.023248959_sp
     f(2)%p3 = -0.076155718_sp
     f(2)%p4 = 0.089234995_sp
     f(2)%p5 = -0.165243333_sp
     f(2)%p6 = 8.762930456_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.00654304_sp
     f(3)%p2 = -0.006417016_sp
     f(3)%p3 = -0.017797747_sp
     f(3)%p4 = 0.027979874_sp
     f(3)%p5 = -0.03755133_sp
     f(3)%p6 = 3.65131458_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.65

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.0109134_sp
     f(1)%p2 = -0.011387177_sp
     f(1)%p3 = -0.02969415_sp
     f(1)%p4 = 0.038489494_sp
     f(1)%p5 = -0.124407932_sp
     f(1)%p6 = 4.702029899_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.016253319_sp
     f(2)%p2 = -0.019715201_sp
     f(2)%p3 = -0.043538947_sp
     f(2)%p4 = 0.071091559_sp
     f(2)%p5 = -0.168978461_sp
     f(2)%p6 = 8.785049268_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.003379053_sp
     f(3)%p2 = -0.006102296_sp
     f(3)%p3 = -0.007927441_sp
     f(3)%p4 = 0.025076372_sp
     f(3)%p5 = -0.038955894_sp
     f(3)%p6 = 3.662798197_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.70

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.005375976_sp
     f(1)%p2 = -0.009419329_sp
     f(1)%p3 = -0.011912365_sp
     f(1)%p4 = 0.026005524_sp
     f(1)%p5 = -0.126077485_sp
     f(1)%p6 = 4.712550452_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.007067711_sp
     f(2)%p2 = -0.01768926_sp
     f(2)%p3 = -0.016037651_sp
     f(2)%p4 = 0.05359961_sp
     f(2)%p5 = -0.169839897_sp
     f(2)%p6 = 8.807975876_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.000762383_sp
     f(3)%p2 = -0.00652743_sp
     f(3)%p3 = -0.001414936_sp
     f(3)%p4 = 0.02203551_sp
     f(3)%p5 = -0.038164108_sp
     f(3)%p6 = 3.674702438_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.75

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = 0.001976462_sp
     f(1)%p2 = -0.004990076_sp
     f(1)%p3 = -0.002519194_sp
     f(1)%p4 = 0.009521643_sp
     f(1)%p5 = -0.120845481_sp
     f(1)%p6 = 4.722646012_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.003920241_sp
     f(2)%p2 = -0.010083831_sp
     f(2)%p3 = -0.008891679_sp
     f(2)%p4 = 0.02419484_sp
     f(2)%p5 = -0.157652421_sp
     f(2)%p6 = 8.831121234_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.001400649_sp
     f(3)%p2 = -0.004017158_sp
     f(3)%p3 = -0.004648058_sp
     f(3)%p4 = 0.011682748_sp
     f(3)%p5 = -0.032889549_sp
     f(3)%p6 = 3.687296338_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.80

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = -3.76464E-05_sp
     f(1)%p2 = -0.00072531_sp
     f(1)%p3 = 0.001035542_sp
     f(1)%p4 = -0.005488626_sp
     f(1)%p5 = -0.109460599_sp
     f(1)%p6 = 4.733227051_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = 0.000590765_sp
     f(2)%p2 = -0.001347994_sp
     f(2)%p3 = -0.00117746_sp
     f(2)%p4 = -0.005797057_sp
     f(2)%p5 = -0.145773726_sp
     f(2)%p6 = 8.852928659_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = 0.000487814_sp
     f(3)%p2 = -0.00050459_sp
     f(3)%p3 = -0.001606874_sp
     f(3)%p4 = -4.22179E-05_sp
     f(3)%p5 = -0.03369994_sp
     f(3)%p6 = 3.698255557_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.85

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = -0.005309834_sp
     f(1)%p2 = -0.004908326_sp
     f(1)%p3 = 0.011505968_sp
     f(1)%p4 = -0.000598491_sp
     f(1)%p5 = -0.114153777_sp
     f(1)%p6 = 4.772411234_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = -0.006190071_sp
     f(2)%p2 = -0.005686021_sp
     f(2)%p3 = 0.01333081_sp
     f(2)%p4 = -0.001933834_sp
     f(2)%p5 = -0.151365034_sp
     f(2)%p6 = 8.906918271_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = -0.002948922_sp
     f(3)%p2 = -0.002923342_sp
     f(3)%p3 = 0.006331245_sp
     f(3)%p4 = 0.002866888_sp
     f(3)%p5 = -0.037087786_sp
     f(3)%p6 = 3.71529793_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.90

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = -0.011845216_sp
     f(1)%p2 = -0.011460222_sp
     f(1)%p3 = 0.025291559_sp
     f(1)%p4 = 0.010388758_sp
     f(1)%p5 = -0.125923348_sp
     f(1)%p6 = 4.822881794_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = -0.013833547_sp
     f(2)%p2 = -0.013312537_sp
     f(2)%p3 = 0.029430101_sp
     f(2)%p4 = 0.010847734_sp
     f(2)%p5 = -0.163561582_sp
     f(2)%p6 = 8.974855977_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = -0.006650117_sp
     f(3)%p2 = -0.006608626_sp
     f(3)%p3 = 0.014137745_sp
     f(3)%p4 = 0.008823015_sp
     f(3)%p5 = -0.040456308_sp
     f(3)%p6 = 3.735330059_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 0.95

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = -0.018407236_sp
     f(1)%p2 = -0.018157002_sp
     f(1)%p3 = 0.039278526_sp
     f(1)%p4 = 0.021593614_sp
     f(1)%p5 = -0.138490345_sp
     f(1)%p6 = 4.877084486_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = -0.021537012_sp
     f(2)%p2 = -0.021163343_sp
     f(2)%p3 = 0.045832963_sp
     f(2)%p4 = 0.024002496_sp
     f(2)%p5 = -0.176979886_sp
     f(2)%p6 = 9.047342415_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = -0.010348388_sp
     f(3)%p2 = -0.010329093_sp
     f(3)%p3 = 0.021961238_sp
     f(3)%p4 = 0.014868924_sp
     f(3)%p5 = -0.044034441_sp
     f(3)%p6 = 3.755985264_sp

  elseif(reg%comp%basalt_fraction>0.and.reg%comp%basalt_fraction<up) then
     ! Coefficients for basalt fraction of 1.00

     ! Enter coefficients for polynomial fit for vs                            
     f(1)%p1 = -0.025037086_sp
     f(1)%p2 = -0.024994142_sp
     f(1)%p3 = 0.053467106_sp
     f(1)%p4 = 0.033037613_sp
     f(1)%p5 = -0.151265978_sp
     f(1)%p6 = 4.933536999_sp

     ! Enter coefficients for polynomial fit for vp                            
     f(2)%p1 = -0.029349011_sp
     f(2)%p2 = -0.029212374_sp
     f(2)%p3 = 0.062564025_sp
     f(2)%p4 = 0.037485344_sp
     f(2)%p5 = -0.190914672_sp
     f(2)%p6 = 9.122506357_sp

     ! Enter coefficients for polynomial fit for rho                           
     f(3)%p1 = -0.014053394_sp
     f(3)%p2 = -0.014068838_sp
     f(3)%p3 = 0.029818145_sp
     f(3)%p4 = 0.020957802_sp
     f(3)%p5 = -0.047677764_sp
     f(3)%p6 = 3.776954899_sp

  end if
end function temperature2velocity
