$PROB  Linear growth with first order decay
# dT/dt=Kgrow - Kdecay*T

$PARAM @annotated
TumorB: 300  : Tumor Size Baseline (mm)
Kgrow: 0.03 : Tumor Growth Rate (mm/hour)
Kdecay:0.00002 : Tumor nature shrink rate(1/hour)

$OMEGA
0.3 0.5 0.5

$CMT  @annotated
Tumor  : Tumor size(mm)

$MAIN  

//double Tumor_0 = TumorB;
//double KG = Kgrow;
//double KD = Kdecay;

double Tumor_0 = TumorB* exp(ETA(1));
double KG = Kgrow * exp(ETA(2));
double KD = Kdecay * exp(ETA(3));

$ODE
dxdt_Tumor = KG - KD * Tumor;

$ERROR

$CAPTURE KG KD  Tumor_0