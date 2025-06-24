$PROB Cabozantinib  ADVAN4 TRANS4 
# Model: `pk2cmt`
- Two-compartment PK model
- first-order absorption represent by KA
- There is a lag time on drug absorption


$PARAM @annotated
KA  :  0.02  : Absorption rate constant 1 (1/hour)
CL   :  50  : Clearance (L/hour)
V2   : 347  : Central volume (L)
Q    :  10  : Inter-compartmental clearance (L/hour)
V3   : 70  : Peripheral volume of distribution (L)
//ALAG1: 1.5  : Lag time for KA1 (hour)
TR   : 0.02: Logit transformation factor for F1

$CMT  @annotated
GUT    : Dosing compartment (mg)
CENT   : Central compartment (mg)
PERIPH : Peripheral compartment (mg) 
AUC    : AUC (mg*hour/ml)

$GLOBAL  // visit one time

#define S2 (V2/1000)
#define CP (CENT/S2)
 
$MAIN  // visit every time prior to advancing the systerm

double K = CL/V2;
double K23 = Q/V2;
double K32 = Q/V3;
double F1 = 1 - (TR/(1+TR));

//ALAG_GUT = ALAG1;

$ODE

double CENT2 = std::max(0.0, CENT);

dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - K23*CENT2 + K32*PERIPH - K*CENT2;
dxdt_PERIPH = K23*CENT2 - K32*PERIPH;
dxdt_AUC = CP;

$CAPTURE  @annotated
CP : Plasma concentration (mg/ml)
