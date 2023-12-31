# synchronisation engineering nonpairwise interactions
This repository contains the code accompanying my thesis Designing the Dynamics of Coupled Oscillators at the University of Exeter under supervision of Kyle Wedgwood and Christian Bick.
It can be used to design control oscillators using a nonpairwise target interaction function.

Running the ```main.m``` script will produce a ```.txt``` file in ```/data/``` with order parameters calculated after simulations of coupled Brusselators.
Running the ```load_data_and_plot.m``` script loads data from ```/data/[name].txt``` and produces (bifurcation) plots.

Some of the parameters that can be changed within ```main.m```:
- ```xis``` and ```chis``` (parameters determining the target phase model)
     (in particular: ```xi3=1``` for nonpairwise model from thesis, and ```xi3=0``` for the pairwise model from the thesis)
- ```Dtaus``` (a range of overall delay parameters)
- ```which_oscillator``` (the type of oscillator)
- ```M``` (the number of oscillators)
- ```startat``` (initial configuration of the oscillators, e.g. splay configuration)
- ```K``` ("epsilon", overall coupling strength)

Running the ```load_data_and_plot.m``` script loads the data produced by ```main.m``` (```.txt``` files containing order 
parameters for a range of delays Dtaus), and creates plots including the bifurcation plots in the thesis that
plot ```Dtau``` versus ```log(1-R_k)```.


Below is a list of the scripts and function it calls/data it loads.

	main.m
	---> load_isochrones.m
	---> 	BrusselatorFourierAverages.dat
		FHNFourierAverages.dat
		find_prc.m	
	optimisation.m
	---> 	get_h_functions.m
		vector_to_taus_and_ks.m
	simulations.m
	---> 	oscillators.m
		--->  	brus.m
			fitzhughnagumo.m
	find_order_parameter.m
	
other scripts that can be called seperately after running (part of) ```main.m```:

	plot_resulting_feedback.m  	(after running optimisation)
	load_data_and_plot.m		(after running all of main)
	---> 	[name].txt file in /data/
		getcolours.m
