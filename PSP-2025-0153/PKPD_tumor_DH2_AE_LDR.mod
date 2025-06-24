$PROB 
# PKModel: pk2cmt 
- Two-compartment PK model
- first-order absorption represent by KA
- CL is not an input parameter

# PDModel: Tumor size
- dT/dt=Kgrow - Kdecay(1+INH*exp(-Ktol*t))*T
- INH = EMAX*CP/(CP+EC50)

$PARAM @annotated
KA  :  0.5  : Absorption rate constant 1 (1/hour)
CL   :  10  : Clearance (L/hour)
V2   : 253  : Central volume (L)
Q    :  9.783  : Inter-compartmental clearance (L/hour)
V3   : 70  : Peripheral volume of distribution (L)
TR   : 0.02: Logit transformation factor for F1

TumorB: 300  : Tumor Size Baseline (mm)
Kgrow: 0.03 : Tumor Growth Rate (mm/hour)
Kdecay:0.00002 : Tumor nature shrink rate(1/hour)

EMAX : 70  : Maximum effect
EC50 : 75  : Concentration for 50% of max effect(ng/ml)
n    : 0.8 : Emax model sigmoidicity  
Ktol :0.00065: tolerant rate (1/hour)

B : -9 : the intercept of the logit model
SLP: 4 : the slope of the logit model

$PARAM
DOSEI=200, DOSEII=100, DOSEIII=60, INTERVAL=24, UNTIL=24*600, WHERE=1 // Initial dose and dosing interval

$PLUGIN evtools Rcpp

$GLOBAL
evt::regimen reg;
#define S2 (V2/1000)
#define CP (CENT/S2)

$OMEGA
0.1 0.1

$CMT  @annotated
GUT    : Dosing compartment (mg)
CENT   : Central compartment (mg)
PERIPH : Peripheral compartment (mg) 
AUC    : AUC (mg*hour/ml)
Tumor  : Tumor size compartment(mm)

$MAIN  
double K = CL/V2;
double K23 = Q/V2;
double K32 = Q/V3;
double F1 = 1 - (TR/(1+TR));

double KG = Kgrow;
double KD = Kdecay;
double Tumor_0 = TumorB;

double LPBase = B + ETA(1);
double LPSlope = SLP * exp(ETA(2));

if (NEWIND <= 1) {
double condition_met = 0;
double reduction_first = 0;
double reduction_second = 0;
double reduction_third = 0;
double interruption_days = 0;
double current_dose = DOSEI;
reg.init(self);
reg.amt(DOSEI);
reg.cmt(WHERE);
reg.ii(INTERVAL);
reg.until(UNTIL);
}

$ODE
double GUT2 = std::max(0.0000000001, GUT);
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

double LOGIT = LPBase + LPSlope * CP/1200;
double C1=exp(LOGIT);
double PR=C1/(C1+1);
double AE = R::rbinom(1, PR); // binomial distribution to model random events

if (fmod(TIME, 24) == 0) {
  if (interruption_days > 0) {
    interruption_days--;
    reg.amt(0);
  } else {
    if (AE == 1) {
      interruption_days = 6;
      reg.amt(0);
      if (condition_met == 0) {
        condition_met = 1;
      } else if (reduction_first == 1 && reduction_second == 0) {
        reduction_second = 1;
      } else if (reduction_first == 1 && reduction_second == 1 && reduction_third == 1) {
        // Stop treatment
        reg.amt(0);
        current_dose = 0;
      }
    } else {
      if (condition_met == 1 && reduction_first == 0) {
        current_dose = DOSEII;
        reduction_first = 1;
      } else if (reduction_second == 1 && reduction_third == 0) {
        current_dose = DOSEIII;
        reduction_third = 1;
      }
      reg.amt(current_dose);
    }
  }
}

capture Dose = reg.amt(); // Capture the current dose
reg.execute(); // Execute the regimen updates

$CAPTURE  @annotated
CP : Plasma concentration (mg/ml)
CL : Drug clearance (L/hour)
PR : Probability of having adverse event
AE : Adverse event to have dose interruption

$CAPTURE Tumor_0 KG KD EMAX EC50 Ktol  KA V2 Q V3 TR LPBase LPSlope