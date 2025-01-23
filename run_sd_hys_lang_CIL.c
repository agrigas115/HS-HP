/* Run MD based on setup_MD.py */
/* Alex Grigas               */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h> 

#define pi 3.14159265358979323846

// Compute non-bonded interaction
// ''Spring - LJ''
// Repulsive linear spring
// Short range attracting spring, set by a,b,c
double compute_slj(double coords[][2], int nonbonded_array1[], int nonbonded_array2[], 
					double total_force[][2], int num_neighs, int vlist[], double Lx, double Ly, 
					double vstress[][2], double sigma_ij_array[], int num_atoms, int adhesion_array[], int* z){
	double distance, mag;
	double V = 0;
	int i, index;
	double x1, x2, y1, y2, delta_x, delta_y;
	double sigma_ij;
	double avg_delta = 0;
	int delta_count = 0;
	*z = 0;
	for (i = 0; i < num_neighs; i++){
		index = vlist[i];
		x1 = coords[nonbonded_array1[index]][0];
		y1 = coords[nonbonded_array1[index]][1];
		x2 = coords[nonbonded_array2[index]][0];
		y2 = coords[nonbonded_array2[index]][1];

		delta_x = x1-x2;
		delta_x -= Lx * nearbyint(delta_x / Lx);
		
		delta_y = y1-y2;
		delta_y -= Ly * nearbyint(delta_y / Ly);

		distance = sqrt((delta_x*delta_x) + (delta_y*delta_y));

		sigma_ij = sigma_ij_array[index];
		
		if (adhesion_array[index] == 1){
			mag = (sigma_ij - distance);
			*z += 1;
			if (distance > sigma_ij){
				avg_delta += -mag;
				delta_count += 1;
			}
		}
		else if (adhesion_array[index] == 0){
		
			if (distance <= sigma_ij){
				mag = sigma_ij - distance;
				adhesion_array[index] = 1;
			}
			else{
				mag = 0.0;
			}
		}

		total_force[nonbonded_array1[index]][0] += (delta_x/distance) * mag;
		total_force[nonbonded_array1[index]][1] += (delta_y/distance) * mag;
		total_force[nonbonded_array2[index]][0] += -(delta_x/distance) * mag;
		total_force[nonbonded_array2[index]][1] += -(delta_y/distance) * mag;

		// Compute Vstress
		vstress[0][0] += ((delta_x/distance) * mag)*delta_x;
		vstress[1][1] += ((delta_y/distance) * mag)*delta_y;

		vstress[0][1] += ((delta_x/distance) * mag)*delta_y;
		vstress[1][0] += ((delta_y/distance) * mag)*delta_x;

	}
	return avg_delta / delta_count;
}

void compute_adhesion_repulsion(double coords[][2], int nonbonded_array1[], int nonbonded_array2[], int num_nonbonded, 
								double total_force[][2], int adhesion_array[], double v0, double Lx, double Ly, 
								double ar_force[][2], int num_atoms, double vstress[][2]){

	int i;
	double x1, x2, y1, y2, delta_x, delta_y, distance;

	for (i=0;i<num_nonbonded;i++){

		if (adhesion_array[i] == 1){

			x1 = coords[nonbonded_array1[i]][0];
			y1 = coords[nonbonded_array1[i]][1];
			x2 = coords[nonbonded_array2[i]][0];
			y2 = coords[nonbonded_array2[i]][1];

			delta_x = x1-x2;
			delta_x -= Lx * nearbyint(delta_x / Lx);
			
			delta_y = y1-y2;
			delta_y -= Ly * nearbyint(delta_y / Ly);

			distance = sqrt((delta_x*delta_x) + (delta_y*delta_y));

			ar_force[nonbonded_array1[i]][0] += (delta_x/distance);
			ar_force[nonbonded_array1[i]][1] += (delta_y/distance);
			ar_force[nonbonded_array2[i]][0] += -(delta_x/distance);
			ar_force[nonbonded_array2[i]][1] += -(delta_y/distance);


		}

	}

	
	for (i=0;i<num_atoms;i++){

		x1 = ar_force[i][0];
		y1 = ar_force[i][1];

		distance = sqrt((x1*x1) + (y1*y1));

		if (distance > 1e-14){
			total_force[i][0] += x1 * v0 / distance;
			total_force[i][1] += y1 * v0 / distance;

			// Compute Vstress
			vstress[0][0] += (x1 * v0 / distance)*(coords[i][0] - Lx*nearbyint(coords[i][0] / Lx));
			vstress[1][1] += (y1 * v0 / distance)*(coords[i][1] - Ly*nearbyint(coords[i][1] / Ly));

			vstress[0][1] += (x1 * v0 / distance)*(coords[i][1] - Ly*nearbyint(coords[i][1] / Ly));
			vstress[1][0] += (y1 * v0 / distance)*(coords[i][0] - Lx*nearbyint(coords[i][0] / Lx));

		}
		
	}

}

void update_bub(double coords[][2], int adhesion_array[], double p_on, double p_off, int num_nonbonded){
	int i;
	double rand_unif;
	for (i = 0; i < num_nonbonded; i++){
		if (adhesion_array[i] == 1){
			rand_unif = ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
			if (rand_unif < p_off){
				adhesion_array[i] = 0;
			}
		}
	
	}
	return;
}

int generate_vlist(double coords[][2], int nonbonded_array1[], int nonbonded_array2[],
					int num_nonbonded, int vlist[], double r_l, double Lx, double Ly, 
					double sigma_ij_array[], int adhesion_array[], int num_atoms){
	int i;
	int num_neighs = 0;
	double x1, x2, y1, y2, delta_x, delta_y;
	double distance, rc, rl;
	double cutoff;
	for (i=0; i<num_nonbonded; i++){
		x1 = coords[nonbonded_array1[i]][0];
		y1 = coords[nonbonded_array1[i]][1];
		x2 = coords[nonbonded_array2[i]][0];
		y2 = coords[nonbonded_array2[i]][1];

		cutoff = sigma_ij_array[i] * r_l;

		delta_x = x1-x2;
		delta_x -= Lx * nearbyint(delta_x / Lx);
		
		delta_y = y1-y2;
		delta_y -= Ly * nearbyint(delta_y / Ly);

		distance = sqrt((delta_x*delta_x) + (delta_y*delta_y));

		if (distance < cutoff || adhesion_array[i] == 1){
			vlist[num_neighs] = i;
			num_neighs += 1;
		}
	}
	return num_neighs;
}


void generate_gaussian(int half_sample_size, double gaussian_array[][2]){
	int i, dim;
	double rand_unif1, rand_unif2;
	double rand_gauss1, rand_gauss2;
	double twopi = 6.28318530718;
	double A, B;
	for (dim = 0; dim < 2; dim++){
		for (i = 0; i < half_sample_size; i++){
			rand_unif1 = ((double) rand()+1.) / ((double) RAND_MAX+2.);
			rand_unif2 = ((double) rand()+1.) / ((double) RAND_MAX+2.);

			A = sqrt(-2.*log(rand_unif1));
			B = twopi * rand_unif2;
			rand_gauss1 = A*cos(B);
			rand_gauss2 = A*sin(B);

			gaussian_array[i*2][dim] = rand_gauss1;
			gaussian_array[(i*2)+1][dim] = rand_gauss2;
		}
	}
}

int main(int argc, char *argv[]) {
	// argv[1]: int for sequence
	// argv[2]: run number
	srand(atoi(argv[1])); // seed RNG
	char buf[0x100]; // buffering strings

	//printf("Welcome!\n");
	double total_potential_energy = 0, total_kinetic_energy = 0;
	double lj_energy, wall_energy; // Energy terms
	double v;
	double avg_x, avg_y; // centering
	clock_t time_1, time_2; // timing
	int i, j, k; // indices
	int num_neighs; // verlet-list
	double scaling; // velocity rescaling
	FILE *fp; // opening files
	double r_l = 1.1; // Multiplier for verlet-list skin
	double dt = 0.1; // Time step
	double dt_sqr = dt*dt;
	double dt_lang = dt;

	double target_temp = atof(argv[5]);

	int z;
	double avg_delta;

	double vstress[2][2] = {0};

	double initial_temp, temp;

	int num_atoms = 256;  
    int half_num_atoms = num_atoms / 2;

	// Initialize Coords //
	double coords[num_atoms][2], previous_coords[num_atoms][2], half_coords[num_atoms][2], displacement_array[num_atoms][2];
	double total_force[num_atoms][2], prev_total_force[num_atoms][2], ar_force[num_atoms][2];
	double velocs[num_atoms][2], half_velocs[num_atoms][2], v_prime[num_atoms][2];
	double sigma_i_array[num_atoms];

	for (i=0; i<num_atoms/2; i++){
		sigma_i_array[i] = 0.5;
		sigma_i_array[i+(num_atoms/2)] = 0.7;
	}

	double p_on = atof(argv[2]);
	double p_off = atof(argv[2]);
	double v0 = atof(argv[3]);

	double target_phi = atof(argv[4]);

	// Find initial Box size //
	double initial_phi = 0.1;
	double vol = 0;
	for (i=0; i<num_atoms; i++){
		vol += pi*sigma_i_array[i]*sigma_i_array[i];
	}
	double box_size = sqrt(vol / initial_phi);
	double prev_box_size = box_size;

	// Assume unit mass
	double mass_array[num_atoms];
	for (i=0; i<num_atoms; i++){
		mass_array[i] = 1.0;
	}

//////////

	// Loading Nonbonded Interactions //
	// Count number of nonboned pairs
	int num_nonbonded = 0; 
	for (i=0; i<num_atoms; i++){
		for(j=i+1; j<num_atoms; j++){
			num_nonbonded += 1;
		}
	}
	// Setting Non-bonded parameters
	double* sigma_ij_array; 
	sigma_ij_array = (double*)malloc(num_nonbonded * sizeof(double));

	int count = 0;
	double sigma_ij;
	for (i=0; i<num_atoms; i++){
		for(j=i+1; j<num_atoms; j++){
			sigma_ij = sigma_i_array[i] + sigma_i_array[j];
			sigma_ij_array[count] = sigma_ij;
			count += 1;
		}
	}


	// Loading Nonbonded Index Pairs 
	int* adhesion_array; 
	adhesion_array = (int*)malloc(num_nonbonded * sizeof(int));
	int* nonbonded_array1; 
	nonbonded_array1 = (int*)malloc(num_nonbonded * sizeof(int));
	int* nonbonded_array2; 
	nonbonded_array2 = (int*)malloc(num_nonbonded * sizeof(int));
	count = 0;
	for (i=0; i<num_atoms; i++){
		for(j=i+1; j<num_atoms; j++){
			nonbonded_array1[count] = i;
			nonbonded_array2[count] = j;
			adhesion_array[count] = 0;
			count += 1;
		}
	}

	// Initialize vlist //
	int* vlist; 
	vlist = (int*)malloc(num_nonbonded * sizeof(int));

	int restart;
	double Lx, Ly, prev_Lx, prev_Ly;
	double prev_coords[num_atoms][2], gaussian_array[num_atoms][2];
	int NVT_count = 0;

	// Check if a checkpoint file exists
	snprintf(buf, sizeof(buf), "chk/coords_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
	fp = fopen(buf, "r");
	if (fp != NULL){
		restart = 1;
		// Load Step num, box size, and coords
		if (fscanf(fp," %d",&NVT_count)==0){
	    	printf("Problem\n");
	    }
	    if (fscanf(fp," %lf",&Lx)==0){
	    	printf("Problem\n");
	    }
	    if (fscanf(fp," %lf",&Ly)==0){
	    	printf("Problem\n");
	    }
		for (i=0; i<num_atoms; i++){
			for (j=0; j<2; j++){
				if (fscanf(fp," %lf",&coords[i][j])==0){
					printf("Problem\n");
				}
			}
		}
		fclose(fp);

		// Load adhesion array
		snprintf(buf, sizeof(buf), "chk/coords_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
		fp = fopen(buf, "r");
		for (i=0; i<num_nonbonded; i++){
			if (fscanf(fp," %d",&adhesion_array[i])==0){
				printf("Problem\n");
			}
		}
	}
	else{
		restart = 0;
		// If not from check point, delete old data files
		snprintf(buf, sizeof(buf), "vstress_data/vstress_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
		remove(buf);
		snprintf(buf, sizeof(buf), "bonded_data/bonded_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
		remove(buf);
		snprintf(buf, sizeof(buf), "traj_data/traj_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
		remove(buf);

		// Initialize random dilute packing
		box_size = sqrt(vol / initial_phi);
		prev_box_size = box_size;
		Lx = box_size;
		Ly = box_size;

		/* Getting Initial Coords - RW */
		
		memset(coords, 0, sizeof(coords));
		memset(gaussian_array, 0, sizeof(gaussian_array));
		double delta_x, delta_y;
		//generate_gaussian(half_sample_size, gaussian_array);
		double random_vec_x = 0, random_vec_y = 0;
		int overlap_check = 0, overlap_count = 0;
		double distance, x1, x2, y1, y2;
		for (i=0; i<num_atoms; i++){
			overlap_check = 0;
			// Normalize and scale
			while (overlap_check == 0){
				overlap_count = 0;
				//generate_gaussian(half_sample_size, gaussian_array);
				random_vec_x = (double) (rand()) / (double) (RAND_MAX);
				random_vec_y = (double) (rand()) / (double) (RAND_MAX);
				//if (i == 0){
					//coords[i][0] = 0.;
					//coords[i][1] = 0.;
					//overlap_check = 1;
				//}
				//else{

				coords[i][0] = (box_size * random_vec_x) - (box_size/2.);
				coords[i][1] = (box_size * random_vec_y) - (box_size/2.);

				//Check for overlaps
				for (j=0; j<i; j++){
					
					x1 = coords[i][0];
					y1 = coords[i][1];
					x2 = coords[j][0];
					y2 = coords[j][1];

					delta_x = x1-x2;
					delta_x -= box_size * nearbyint(delta_x / box_size);
					
					delta_y = y1-y2;
					delta_y -= box_size * nearbyint(delta_y / box_size);

					sigma_ij = sigma_i_array[i] + sigma_i_array[j];

					distance = sqrt(delta_x*delta_x + delta_y*delta_y);
					if (distance/sigma_ij < 1){
						overlap_check = 0;
						overlap_count += 1;
					}
					
				}
				if (overlap_count == 0){
					overlap_check = 1;
				}
				//}
			}	
		}
	}

	generate_gaussian(half_num_atoms, gaussian_array);
	for (i=0; i<num_atoms; i++){
		for (j=0; j<2; j++){
			total_force[i][j] = 0.;
			prev_total_force[i][j] = 0.;
			displacement_array[i][j] = 0.;
			velocs[i][j] = 0.;
		}
	}

	memcpy(prev_coords, coords, sizeof(prev_coords));
	prev_box_size = box_size;

	// Begin Velocity Verlet //
	memset(displacement_array, 0, sizeof(displacement_array));
	memset(total_force, 0, sizeof(total_force));
	memset(vstress, 0, sizeof(vstress));

	/*
	total_kinetic_energy = 0;
	for (j = 0; j < num_atoms; j++){
		v = sqrt( (velocs[j][0]*velocs[j][0]) + (velocs[j][1]*velocs[j][1]) );
		total_kinetic_energy += 0.5 * v * v;
	}
	temp = total_kinetic_energy / num_atoms;

	scaling = sqrt(target_temp/temp);
	for (j = 0; j < num_atoms; j++){
		for (k = 0; k < 2; k++){
			velocs[j][k] *= scaling;
		}
	}
	*/

	// NVE //
	double disp = 0, disp_max1 = 0, disp_max2 = 0;
	double disp_x, disp_y;
	double disp_list[num_atoms];
	double min_rc = 0.5;
	int split = 1000;
	num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms);
	avg_delta = compute_slj(coords, nonbonded_array1, nonbonded_array2, total_force, num_neighs, vlist, Lx, Ly, vstress, sigma_ij_array, num_atoms, adhesion_array, &z);
	memcpy(prev_total_force, total_force, sizeof(prev_total_force));
	i = 0;
	

	// UNCHANGING FIRE VARIABLES
	int npPos = 0;
    int npNeg = 0;
    int npPMin = 0;
    
	double vnorm = 0;
	double fnorm = 0;
	double P = 0;
	int disp_threshold_count = 0;

	double alpha0 = 0.25;
	double finc = 1.01;
	double fdec = 0.25;
	double falpha = 0.99;

	double NNEGMAX = 2000;
	double NDELAY = 1000;

	double fireitmax = 1e8;
	double fireitmin = 50;

    // FIRE Parameters #
    double FIRE_alpha = alpha0;
    double FIRE_dt =  dt;
    
    double dtmax = 0.1;
    double dtmin = 1e-2*FIRE_dt;

	double Ftol = 1e-14;
	double Ptol = 1e-7;
	double fcheck = 10*Ftol;

	/* Beginning Jamming */
	bool jammed = 0;
	bool overcompressed, undercompressed;
	double pcheck = 0;
	int fireit = 0;
	int it = 0;

	double delta_phi = 1e-3;
	double phi_sign = 1;
	double Lscale, phi, delta_L;

	fcheck = 0;
	for (j = 0; j < num_atoms; j++){
		fcheck += sqrt(total_force[j][0]*total_force[j][0] + total_force[j][1]*total_force[j][1]);
	}
	fcheck = fcheck/num_atoms;
	phi = vol / (Lx*Ly);
	int cycle = 0;
	if (restart == 0){
		while (cycle < 2){
			//printf("phi = %lf\n",phi);
			if (fabs(pcheck) > 1e-12 && cycle == 0){
	    		delta_phi = 1e-3;
	    	}

	    	if (phi >= 0.82){
	    		delta_phi = -1e-3;
	    		cycle += 1;
	    	}
	    	if (phi < target_phi && cycle != 0){
	    		delta_phi = 1e-3;
	    		cycle += 1;
	    	}

		    it += 1;

		    Lscale = sqrt(phi/(phi+delta_phi));
		    //printf("Lscale = %lf \n", Lscale);
		    delta_L = box_size - (box_size * Lscale);
		    //printf("delta_L = %lf \n", delta_L);
		    box_size *= Lscale;
		    Lx *= Lscale;
		    Ly *= Lscale;
		    for (i=0; i<num_atoms; i++){
		    	for (j=0; j<2; j++){
		    		coords[i][j] *= 1-(delta_L/box_size);
		    	}
		    }

		    phi = vol / (Lx*Ly);

			num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms);
			memset(displacement_array, 0, sizeof(displacement_array));

			
			fireit = 0;
			FIRE_dt = dt;
			fcheck = 10*Ftol;
		    while (fcheck > Ftol){
		        fireit += 1;
		        //printf("fireit = %d \n", fireit);
		        // Velocity-verlet step
		        for (i = 0; i<num_atoms; i++){
		        	velocs[i][0] += 0.5 * FIRE_dt * total_force[i][0];
		        	velocs[i][1] += 0.5 * FIRE_dt * total_force[i][1];
		        }   
		        
		        // 1. Calculate P, fnorm and vnorm
		        vnorm = 0;
		        for (i = 0; i < num_atoms; i++){
		        	vnorm += sqrt(velocs[i][0]*velocs[i][0] + velocs[i][1]*velocs[i][1]);
		        }
		        fnorm = 0;
		        for (i = 0; i < num_atoms; i++){
		        	fnorm += sqrt(total_force[i][0]*total_force[i][0] + total_force[i][1]*total_force[i][1]);
		        }
		        P = 0;
		        for (i=0; i < num_atoms; i++){
		        	P += velocs[i][0]*total_force[i][0] + velocs[i][1]*total_force[i][1];
		        }
		        
		        
		        // 2. Adjust simulation based on net motion
		        if (P > 0){
		        	npPos += 1;
		            
		            npNeg = 0;
		            
		            //alphat = alpha;
		           
		            // alter simulation if enough positive steps have been taken
		            if (npPos > fireitmin){
		            	// change time step
		                if (FIRE_dt*finc < dtmax){
		                    FIRE_dt = FIRE_dt*finc;
		                }
		            }
		            // decrease alpha
		            FIRE_alpha = FIRE_alpha*falpha;
		        }
		        else {
		        	// reset positive counter
		            npPos = 0;

		            // increase negative counter
		            npNeg = npNeg + 1;

		            // check for stuck simulation
		            if (npNeg > NNEGMAX){
		                printf("Simulation negative for too long, ending program here \n");
		            }

		            // decrease time step if past initial delay
		            if (fireit > fireitmin){
		                // decrease time step
		                if (FIRE_dt*fdec > dtmin){
		                    FIRE_dt = FIRE_dt*fdec;
		                }

		                // change alpha
		                FIRE_alpha = alpha0;
		            }

		            // take a half step backwards
		            for (i = 0; i<num_atoms; i++){
		            	coords[i][0] -= 0.5*FIRE_dt*velocs[i][0];
		            	coords[i][1] -= 0.5*FIRE_dt*velocs[i][1];
		            }
		            

		            // reset velocities to 0
		            memset(velocs, 0, sizeof(velocs));
		        }
		            
		        
		        // update velocities if forces are acting
		        if (fnorm > 0){
		        	for (i = 0; i<num_atoms; i++){
		        		velocs[i][0] = (1 - FIRE_alpha) * velocs[i][0] + FIRE_alpha*(total_force[i][0]/fnorm)*vnorm;
		        		velocs[i][1] = (1 - FIRE_alpha) * velocs[i][1] + FIRE_alpha*(total_force[i][1]/fnorm)*vnorm;
		        	}
		        }

		        // do first verlet update for vertices (assume unit mass)
		        for (i = 0; i<num_atoms; i++){
		        	coords[i][0] += 0.5*FIRE_dt*velocs[i][0];
		        	coords[i][1] += 0.5*FIRE_dt*velocs[i][1];
		        }


		        // Compute new Force
		        memset(total_force, 0, sizeof(total_force));
		        memset(vstress, 0, sizeof(vstress));

		        // Checking Disp for update on Verlet List 
				disp_threshold_count = 0;
				for (j=0; j<num_atoms; j++){
					disp_x = displacement_array[j][0]-coords[j][0];
					disp_y = displacement_array[j][1]-coords[j][1];
					disp = sqrt(disp_x*disp_x + disp_y*disp_y);
					disp_list[j] = disp;
				}

				// Find first max 
				disp_max1 = 0;
				for (j=0; j<num_atoms; j++){
					if (disp_list[j] > disp_max1){
						disp_max1 = disp_list[j];
					}
				}

				// Find second max 
				disp_max2 = 0;
				for (j=0; j<num_atoms; j++){
					if (disp_list[j] > disp_max2 && disp_list[j] < disp_max1){
						disp_max2 = disp_list[j];
					}
				}

		        if (disp_max1 + disp_max2 > (min_rc*r_l) - min_rc){
					memcpy(displacement_array, coords, sizeof(displacement_array));
					num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms);
				}

				avg_delta = compute_slj(coords, nonbonded_array1, nonbonded_array2, total_force, num_neighs, vlist, Lx, Ly, vstress, sigma_ij_array, num_atoms, adhesion_array, &z);

		        // Do second verlet update for velocities
		        for (i = 0; i<num_atoms; i++){
		        	velocs[i][0] += 0.5 * FIRE_dt * total_force[i][0];
		        	velocs[i][1] += 0.5 * FIRE_dt * total_force[i][1];
		        }   
		        // update Kcheck and Fcheck
		        // Calc avg norm of total_force
				fcheck = 0;
				for (j = 0; j < num_atoms; j++){
					fcheck += sqrt(total_force[j][0]*total_force[j][0] + total_force[j][1]*total_force[j][1]);
				}
		    	fcheck = fcheck/num_atoms;

		        if (fcheck < Ftol){
		            npPMin = npPMin + 1;
		        }
		        else{
		            npPMin = 0;
		        }

			}

			pcheck = (vstress[0][0] + vstress[1][1]) / (2*Lx*Ly);
			phi = vol / (Lx*Ly);
			double rho = sqrt(num_atoms / (Lx*Ly));

		    
			// Save coords //
			memcpy(prev_coords, coords, sizeof(prev_coords));
		    prev_Lx = Lx; 
		    prev_Ly = Ly;

		    /*
		    snprintf(buf, sizeof(buf), "bonded.txt");
			FILE *out_file_bonded = fopen(buf, "a+");
			for (i = 0; i < num_nonbonded; i++){
				if (adhesion_array[i] == 1){
					fprintf(out_file_bonded, "%d\t%d\t",nonbonded_array1[i], nonbonded_array2[i]);
				}
			}
			fprintf(out_file_bonded, "\n");
			fclose(out_file_bonded);


			snprintf(buf, sizeof(buf), "traj.txt");
			FILE* out_file_coords_final = fopen(buf, "a+"); // write only 
			fprintf(out_file_coords_final, "%.17g\n", box_size);
			for (j = 0; j < num_atoms; j++){
				fprintf(out_file_coords_final, "%.17g\t%.17g\t%.17g\n", coords[j][0], coords[j][1], sigma_i_array[j]); // write to file
			}
			fclose(out_file_coords_final);
			*/
		
		}
	}

	double gamma = atof(argv[6]);
	double exp_gamma_dt =  exp(-gamma * dt_lang);
	double fluc_coeff = sqrt(1. - exp(-2. * gamma * dt_lang)) * sqrt(target_temp);
	i = -1;
	int cycles = 0;

	double phi_center = phi;
	double phi_0 = 0.1;
	double phi_target = phi;

	t = 0;
	NVT_count = 0
	while (NVT_count < 1e9){

		delta_phi = phi_target - phi; 

		Lscale = sqrt(phi/(phi+delta_phi));
	    //printf("Lscale = %lf \n", Lscale);
	    delta_L = box_size - (box_size * Lscale);
	    //printf("delta_L = %lf \n", delta_L);
	    box_size *= Lscale;
	    Lx *= Lscale;
	    Ly *= Lscale;
	    for (i=0; i<num_atoms; i++){
	    	for (j=0; j<2; j++){
	    		coords[i][j] *= 1-(delta_L/box_size);
	    	}
	    }

	    phi = vol / (Lx*Ly);
	    printf("phi = %lf\n",phi);

	    phi_target = phi_0*sin(omega*t) + phi_center;

		NVT_count += 1;

		// Update Half Velocities //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				half_velocs[j][k] = velocs[j][k] + (0.5 * dt_lang * total_force[j][k]);
			}
		}

		// Update Half Coordinates //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				half_coords[j][k] = coords[j][k] + 0.5 * dt_lang * half_velocs[j][k];
			}
		}

		// Calc v_prime //
		generate_gaussian(num_atoms/2, gaussian_array);
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				v_prime[j][k] = exp_gamma_dt * half_velocs[j][k] + fluc_coeff * gaussian_array[j][k];
			}
		}

		// Update Coordinates //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				coords[j][k] = half_coords[j][k] + 0.5 * dt_lang * v_prime[j][k];
			}
		}
		// Compute new Force
        memset(total_force, 0, sizeof(total_force));
        memset(vstress, 0, sizeof(vstress));
        memset(ar_force, 0, sizeof(ar_force));


        // Checking Disp for update on Verlet List 
		disp_threshold_count = 0;
		for (j=0; j<num_atoms; j++){
			disp_x = displacement_array[j][0]-coords[j][0];
			disp_y = displacement_array[j][1]-coords[j][1];
			disp = sqrt(disp_x*disp_x + disp_y*disp_y);
			disp_list[j] = disp;
		}

		// Find first max 
		disp_max1 = 0;
		for (j=0; j<num_atoms; j++){
			if (disp_list[j] > disp_max1){
				disp_max1 = disp_list[j];
			}
		}

		// Find second max 
		disp_max2 = 0;
		for (j=0; j<num_atoms; j++){
			if (disp_list[j] > disp_max2 && disp_list[j] < disp_max1){
				disp_max2 = disp_list[j];
			}
		}

        if (disp_max1 + disp_max2 > (min_rc*r_l) - min_rc){
			memcpy(displacement_array, coords, sizeof(displacement_array));
			num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, vlist, r_l, Lx, Ly, sigma_ij_array, adhesion_array, num_atoms);
		}

		avg_delta = compute_slj(coords, nonbonded_array1, nonbonded_array2, total_force, num_neighs, vlist, Lx, Ly, vstress, sigma_ij_array, num_atoms, adhesion_array, &z);
		compute_adhesion_repulsion(coords, nonbonded_array1, nonbonded_array2, num_nonbonded, total_force, adhesion_array, v0, Lx, Ly, ar_force, num_atoms, vstress);

		// Update Velocities //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 2; k++){
				velocs[j][k] = v_prime[j][k] + 0.5 * dt_lang * total_force[j][k];
			}
		} 

		// Average Velocs //
		avg_x = 0;
		avg_y = 0;
		for (j=0; j<num_atoms; j++){
			avg_x += velocs[j][0];
			avg_y += velocs[j][1];
		}
		avg_x /= num_atoms;
		avg_y /= num_atoms;

		// Add Veloc comp. to vstress //
		for (j = 0; j < num_atoms; j++){

			vstress[0][0] += (velocs[j][0] - avg_x) * (velocs[j][0] - avg_x);
			vstress[1][1] += (velocs[j][1] - avg_y) * (velocs[j][1] - avg_y);

			vstress[0][1] += (velocs[j][0] - avg_y) * (velocs[j][1] - avg_x);
			vstress[1][0] += (velocs[j][1] - avg_x) * (velocs[j][0] - avg_y);

		}
		
		// Normalized by system volume //
		vstress[0][0] /= (Lx*Ly);
		vstress[1][1] /= (Lx*Ly);
		vstress[0][1] /= (Lx*Ly);
		vstress[1][0] /= (Lx*Ly);

		t += dt_lang;

		if (NVT_count % split == 0){
			update_bub(coords, adhesion_array, p_on, p_off, num_nonbonded);

			snprintf(buf, sizeof(buf), "vstress_data/vstress_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
			FILE *out_file_vstress = fopen(buf, "a+");
			fprintf(out_file_vstress, "%.16g\t%.16g\t%.16g\t%.16g\t%.16g\n",avg_delta,vstress[0][0],vstress[1][1],vstress[0][1],phi);
			fclose(out_file_vstress);

			/*
			snprintf(buf, sizeof(buf), "chk/coords_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
			FILE* out_file_coords_chk = fopen(buf, "w"); // write only 
			fprintf(out_file_coords_chk, "%d\t%.10g\t%.10g\n", NVT_count, Lx, Ly);
			for (j = 0; j < num_atoms; j++){
				fprintf(out_file_coords_chk, "%.10g\t%.10g\n", coords[j][0], coords[j][1]); // write to file
			}
			fclose(out_file_coords_chk);

			snprintf(buf, sizeof(buf), "chk/adhesion_array_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.chk",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
			FILE* out_file_adhesion_chk = fopen(buf, "w"); // write only
			for (j = 0; j < num_nonbonded; j++){
				fprintf(out_file_adhesion_chk, "%d\n", adhesion_array[j]);
			}
			fclose(out_file_adhesion_chk);
			*/
		}

		if (NVT_count % 100000 == 0){
			snprintf(buf, sizeof(buf), "bonded_data/bonded_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
			FILE *out_file_bonded = fopen(buf, "a+");
			for (j = 0; j < num_nonbonded; j++){
				if (adhesion_array[j] == 1){
					fprintf(out_file_bonded, "%d\t%d\t",nonbonded_array1[j], nonbonded_array2[j]);
				}
			}
			fprintf(out_file_bonded, "\n");
			fclose(out_file_bonded);

			snprintf(buf, sizeof(buf), "traj_data/traj_%s_p%s_v0_%s_phi%s_T_%s_gamma_%s.txt",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
			FILE* out_file_coords_final = fopen(buf, "a+"); // write only 
			fprintf(out_file_coords_final, "%.10g\t%.10g\n", Lx, Ly);
			for (j = 0; j < num_atoms; j++){
				fprintf(out_file_coords_final, "%.10g\t%.10g\n", coords[j][0], coords[j][1]); // write to file
			}
			fclose(out_file_coords_final);
		}
	}

}

