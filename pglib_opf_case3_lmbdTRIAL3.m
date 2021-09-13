function mpc = pglib_opf_case3_lmbdTRIAL3
mpc.version = '4';
mpc.baseMVA = 1.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 .79  .5	 0.0	 0.0	 1	    1.00000	    0.00000	 1.0	 1	    .95     .86;
	2	 1	 0	 0	 0.0	 0.0	 1	    1.00000	    0.00000	 1.0	 1	    1.05     .97;
	3	 1	 0 	 0	 0.0	 0.0	 1	    1.00000	    0.00000	 1.0	 1      1.21     .80;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
mpc.gen = [
	1	 10.0	 0.0	 0.9	     0	 1.0	 1.0	 1	 2	 0.0;
	2	 10.0	 0.0	 0.2	 0	 1.0	 1.0	 1	 1.2	 0.0;
	3	 0.0	 0.0	 2.0  	 0	 1.0	 1.0	 1	 2   0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	   1.0	   10	   0.000000;
	2	 0.0	 0.0	 3	   1.0	   10	   0.000000;
	2	 0.0	 0.0	 3	   1.0	   10	   0.000000;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 3	 0.01	 0.01	 0	 1.0	 1.0	 1.0	 0.0	 0.0	 1	 -30.0	 30.0;
	3	 2	 0.01	 0.01	 0	 1.0	 1.0	 1.0	 0.0	 0.0	 1	 -30.0	 30.0;
	1	 2	 0.01	 0.01	 0	 1.0	 1.0	 1.0	 0.0	 0.0	 1	 -30.0	 30.0;
];