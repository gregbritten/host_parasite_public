data {
  int<lower=0> N; //number of data points
  vector[N]    I;    //read new data
  vector[N]    P;
  real         dt;
  int          time_points[N];
  real         sigma_I[N];
  real         sigma_P[N]; 
  int          maxt;
  real         K;
}
parameters {
	real<lower=0.0,upper=2.0*12080.0> H0; //initial conditions are parameters
	real<lower=0.0> I0;
	real<lower=0.0> P0;
	real<lower=0.0> r;  //biological parameters
	real<lower=0.0> h1;
	real<lower=0.0> h2;
	real<lower=0.0> a;
	real<lower=0.0> m;
	real<lower=0.0> e;
}
transformed parameters{	
	real x[maxt,3];    //matrix to store variables while looping
	
	real growH[maxt];       //declare intermediaries
	real lossH[maxt];
	real dHdt;
	real hand[maxt];
	real dIdt;
	real growP[maxt];
	real lossP[maxt];
	real dPdt;
	
	x[1,1] = H0;   //specify initial conditions in x
	x[1,2] = I0;
	x[1,3] = P0;
     
  	growH[1] = r * x[1,1] * (1-x[1,1]/K);
  	lossH[1] = (a * x[1,1]/(1+(a * h1 * x[1,1]))) * x[1,3];
  	hand[1]  = x[1,2]/h2;
  	growP[1] = e * hand[1];
 	lossP[1] = m * x[1,3];
	 
    for(j in 2:maxt){         //time step the model
	  growH[j] = r * x[j-1,1] * (1-x[j-1,1]/K);
	  lossH[j] = (a * x[j-1,1]/(1+(a * h1 * x[j-1,1]))) * x[j-1,3];
	  dHdt  = growH[j] - lossH[j];
    
	  hand[j] = x[j-1,2]/h2;
	  dIdt = lossH[j] - hand[j];
    
	  growP[j] = e * hand[j];
	  lossP[j] = m * x[j-1,3];
	  dPdt  = growP[j] - lossP[j];

  	  x[j,1] = fmax(x[j-1,1] + dHdt*dt, 1E-15); //store variable in x; fmax disallows negative values
  	  x[j,2] = fmax(x[j-1,2] + dIdt*dt, 1E-15);
  	  x[j,3] = fmax(x[j-1,3] + dPdt*dt, 1E-15);
    }  
}
model {
	
	//--PRIORS-------------------------------------//
	r  ~ normal(0.2666667,    0.07371115);	
	a  ~ normal(1.697333e-07, 2.350159e-07);
	h1 ~ normal(2.736667,    0.2542309);
    h2 ~ normal(0.5,          0.5);
	m  ~ normal(0.5866667,    0.3910669);
	e  ~ normal(396.6667,     0.5*309.8925);

	I0 ~ uniform(0, 10.0*I[1]);
	P0 ~ uniform(0, 10.0*P[1]);

	//--LIKELIHOOD-----------------------------------//
	I ~ normal(x[time_points,2], sigma_I);
	P ~ normal(x[time_points,3], sigma_P);
}

