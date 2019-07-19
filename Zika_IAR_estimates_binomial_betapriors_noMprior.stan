data{
  int<lower=0> x; //Number of areas
  int<lower=0> xB;
  int N[x];       //Population size in each area
  int Npreg[x];   //Number of pregnant women in each area
  int Nbirths[x]; //Number of births in each area
  int C[x];       //Confirmed cases
  int S[x];   //Suspected cases
  int Cpreg[x];       //Confirmed cases in pregnant women
  int Spreg[x];   //Suspected cases in pregnant women
  int M[x];       //Microcephaly cases
  int G[x];       //GBS cases
  int Call;       //country total cases
  int Mall;       //country total Microcephaly cases
  int Gall;       //country total GBS cases
  int Sall;       //country total suspected cases
  int Cpregall;       //country total cases
  int Spregall;       //country total suspected cases
}

parameters{
  real<lower=0> alphaC;
  real<lower=0> alphaS;
  real<lower=0> alphaCpreg;
  real<lower=0> alphaSpreg;
  real<lower=0> betaC;
  real<lower=0> betaS;
  real<lower=0> betaCpreg;
  real<lower=0> betaSpreg;
  real<lower=0,upper=1> pC[x];
  real<lower=0,upper=1> pS[x];
  real<lower=0,upper=1> pCpreg[x];
  real<lower=0,upper=1> pSpreg[x];
  real<lower=0,upper=1> pM;
  real<lower=0,upper=1> pG;
  real<lower=0,upper=1> pSlat;
  real<lower=0,upper=1> varI[x]; //infection rate in each location

}

transformed parameters{

  vector[xB] varA;
  real<lower=0,upper=1> pCall;
  real<lower=0,upper=1> pSall;
  real<lower=0,upper=1> pCpregall;
  real<lower=0,upper=1> pSpregall;


  varA = rep_vector(0,xB);
  pCall = 0; //rep_vector(0,xB);
  pSall = 0; //rep_vector(0,xB);
  pCpregall = 0; //rep_vector(0,xB);
  pSpregall = 0; //rep_vector(0,xB);
      
  for (i in 1:x){
    varA = varA + varI[i]*(N[i]/(1.0*sum(N)));
    pCall = pCall +pC[i]*(N[i]/(1.0*sum(N)));
    pSall = pSall +pS[i]*(N[i]/(1.0*sum(N)));
    pCpregall = pCpregall +pCpreg[i]*(N[i]/(1.0*sum(N)));
    pSpregall = pSpregall +pSpreg[i]*(N[i]/(1.0*sum(N)));
  }
}

model{
  //Priors
  for(i in 1:x){
    varI[i] ~ beta(1,2);
  }
  //pM ~ beta(5,114);
 
  pSlat ~ beta(3.880856,5.344922);
  
  alphaC ~ cauchy(0,25);
  alphaS ~ cauchy(0,25);
  alphaCpreg ~ cauchy(0,25);
  alphaSpreg ~ cauchy(0,25);
  betaC ~ cauchy(0,25);
  betaS ~ cauchy(0,25);
  betaCpreg ~ cauchy(0,25);
  betaSpreg ~ cauchy(0,25);

  for (i in 1:x){
    pC[i] ~ beta(alphaC,betaC);
    pS[i] ~ beta(alphaS,betaS);
    pCpreg[i] ~ beta(alphaCpreg,betaCpreg);
    pSpreg[i] ~ beta(alphaSpreg,betaSpreg);
    C[i] ~ binomial(N[i],varI[i]*pSlat*pC[i]);
    S[i] ~ binomial(N[i],varI[i]*pSlat*pS[i]);
    Cpreg[i] ~ binomial(Npreg[i],varI[i]*pSlat*pCpreg[i]);
    Spreg[i] ~ binomial(Npreg[i],varI[i]*pSlat*pSpreg[i]);
    M[i] ~ binomial(Nbirths[i],pM*varI[i]);
    G[i] ~ binomial(N[i],pG*pSlat*varI[i]);
  }
  Mall ~ binomial(sum(Nbirths),pM*varA);
  Gall ~ binomial(sum(N),varA*pSlat*pG); 
  Call ~ binomial(sum(N),varA*pSlat*pCall);  
  Sall ~ binomial(sum(N),varA*pSlat*pSall);  
  Cpregall ~ binomial(sum(Npreg),varA*pSlat*pCpregall);  
  Spregall ~ binomial(sum(Npreg),varA*pSlat*pSpregall); 
}

