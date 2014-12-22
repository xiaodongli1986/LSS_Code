
Next Steps:
	# model-independent parameters alpha_parallel, alpha_perpendicular 

August 5:
	# Calculation of beta skeleton! All follow-up statistics!!!
	
	Usage: ./BSK_calc -om omegam -w w -input intpufilename -output outputfilname -beta beta -printinfo printinfo -numNNB numNNB -minimalrcut minimalr_cut -maximalrcut maximalr_cut -nglcrosscheck do_nglcrosscheck -ngldir ngldir -numdrop numdrop -dropstep dropstep
	
	Example: ./BSK_calc -om 0.26 -w -1 -input test.xyz.txt -output testoutput -beta 5.0 -printinfo True -numNNB 30 -nglcrosscheck False -numdrop 2 -dropstep 0.5

	 * Fast enough; for beta=5, 10 million points, time consumed is 10 minutes.
	 * Although it can be paralleled easily I did not do so because it's fast enough.
 	 * Now it only supports beta>1.0 but it can be easily extend to also support beta<1.0 (just define the centers, radius of the circles in ap_smooth.f90/BSK_Connected). 
 	 * Outputs:
 	 	* om??????_w???????.BSKIndex: 
 	 		index of points in connections, start from 0
 	 	* om??????_w???????.nglBSK: 
 	 		result of ngl as cross-check
 	 	* om??????_w???????.BSKinfo: 
 	 		information of connections. fmt: distance (middle point), redshift (middle point), length, mu
 	 	* om??????_w???????_xyz.txt: 
 	 		x,y,z of inputfile, converted to that om/w cosmology
 	 	* om??????_w???????_BSKmuinfo_idrop?.***.txt: 
 	 		statistical results of mu
 	 			"jkf" means error bars estimated from jackknife; 
 	 			'noabs' means do not take the absolute of mu 
 
	-input intpufilename 
		Inputfile. Positions of x,y,z in Om=0.26/w=-1 cosmology. 
		
	-om omegam -w w	
		Input x,y,z will be automatically converted into xyz in another cosmology, and do calculation. 
	
	-output outputfilname 
		Name of the outputfile. In fact it will creat a directory and put everything inside it.
	
	-beta beta
		Value of beta
		
	-printinfo printinfo
		True or False.
	
	-numNNB numNNB
		I limit the searching region to numNNB nearest neighbors, making the computing time scale as N. For beta = 10/5/3/1, numNNB = 15/20/25/100 will be enough.
		
	-numdrop numdrop 
		Integer. If numdrop=1 then will use the whole BSK. If numdrop>1 then will output statistical results after dropping a ratio of short connections. 
		
	-dropstep dropstep
		Step size of dropping. For example, numdrop=5, dropstep=0.1 corresponding to results of dropping 0%, 10%, 20%, 30%, 40% connections.
		
	-nglcrosscheck do_nglcrosscheck -ngldir ngldir
		Cross check with ngl. do_nglcrosscheck: True of False; ngldir: directory of NGL.
		
	-minimalrcut minimalr_cut -maximalrcut maximalr_cut
		To make it easy to select out a shell, I add these options. Set them as values like -1.0e30 and 1.0e30 if you do not want to use them.
	
July 22:
	
	Gradient Field Version 4.
	
	# Updated ra/dec random subroutines -- taking into consideration the distorted shape of the circle when a sphere is projected to 2D ra/dec plane
	# information of four BOSS catalogues: detailed information, including recommended settings, are added to settings_init.f90
	# jackknife 

	
July 9:

	Gradient Field Version 3.

	Supported data types:
		# real data
		# constant nbar(r) mock; not constant nbar(r), RSD, noRSD mock
		# vlos added mock
			
