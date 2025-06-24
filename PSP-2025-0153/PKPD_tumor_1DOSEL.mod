$PROB Drug indirectly increase Kdecay with an Emax effect 
# PKModel: `pk2cmt`
- Two-compartment PK model
- first-order absorption represent by KA
- CL is not an input parameter

# PDModel:'Tumor size'
- dT/dt=Kgrow - Kdecay(1+INH*exp(-Ktol*t))*T
- INH = EMAX*CP/(CP+EC50)

$PARAM @annotated
KA  :  0.5  : Absorption rate constant 1 (1/hour)
cl   :  10  : Clearance (L/hour)
V2   : 253  : Central volume (L)
Q    :  9.783  : Inter-compartmental clearance (L/hour)
V3   : 70  : Peripheral volume of distribution (L)
TR   : 0.02: Logit transformation factor for F1

TumorB: 300  : Tumor Size Baseline (mm)
Kgrow: 0.03 : Tumor Growth Rate (mm/hour)
Kdecay:0.00002 : Tumor nature shrink rate(1/hour)

emax : 70  : Maximum induction
ec50 : 75  : Concentration for 50% of max induction(ng/ml)
n    : 0.8 : Emax model sigmoidicity  
ktol :0.00065: tolerant rate (1/hour)

$OMEGA
0.6 0.5 0.5 0.15

$CMT  @annotated
GUT    : Dosing compartment (mg)
CENT   : Central compartment (mg)
PERIPH : Peripheral compartment (mg) 
AUC    : AUC (mg*hour/ml)
Tumor  : Tumor size compartment(mm)

$GLOBAL  // visit one time
#define S2 (V2/1000)
#define CP (CENT/S2)
 
$MAIN  

double CL = cl * exp(ETA(1));
double K = CL/V2;
double K23 = Q/V2;
double K32 = Q/V3;
double F1 = 1 - (TR/(1+TR));

double KG = Kgrow;
double KD = Kdecay;
double Tumor_0 = TumorB;

double EMAX = emax * exp(ETA(2));
double EC50 = ec50 * exp(ETA(3));
double Ktol = ktol * exp(ETA(4));

$ODE

double GUT2 = std::max(0.0, GUT);
double CENT2 = std::max(0.0, CENT);

dxdt_GUT = -KA*GUT2;
dxdt_CENT = KA*GUT2 - K23*CENT2 + K32*PERIPH - K*CENT2;
dxdt_PERIPH = K23*CENT2 - K32*PERIPH;
dxdt_AUC = CP;

#define t SOLVERTIME

double Tumor_constrained = std::max(Tumor, 0.0);
dxdt_Tumor = KG - KD *(1+(INH*exp(-Ktol*t))) * Tumor_constrained;

$ERROR
double INH = (pow(EMAX,n)*pow(CP,n))/(pow(CP,n)+pow(EC50,n));

$CAPTURE  @annotated
CP : Plasma concentration (mg/ml)
CL : Drug clearance (L/hour)

$CAPTURE Tumor_0 KG KD EMAX EC50 Ktol KA V2 Q V3 TR