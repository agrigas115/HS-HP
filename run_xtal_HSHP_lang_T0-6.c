/* Run MD based on setup_MD.py */
/* Alex Grigas               */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h> 

#define pi 3.14159265358979323846

double pow_int(double x, int n){
	int i;
	double out = 1;
	for (i = 0; i < n; i++){
		out *= x;
	}
	return out;
}

// Non-bonded interaction
double compute_urlj(double coords[][3], int nonbonded_array1[], int nonbonded_array2[], 
					double sigma_array[], double a_array[], double k_array[], 
					double rb_array[], double total_force[][3], 
					int num_neighs, int vlist[], int* z){
	double distance;
	double mag;
	double u_rlj_energy = 0;
	int i, j, index;
	double x1, x2, y1, y2, z1, z2, delta_x, delta_y, delta_z;
	double sigma_ij, a, k, rb;
	*z = 0;
	for (i = 0; i < num_neighs; i++){
		index = vlist[i];
		x1 = coords[nonbonded_array1[index]][0];
		y1 = coords[nonbonded_array1[index]][1];
		z1 = coords[nonbonded_array1[index]][2];
		x2 = coords[nonbonded_array2[index]][0];
		y2 = coords[nonbonded_array2[index]][1];
		z2 = coords[nonbonded_array2[index]][2];

		// Measure distance between atoms
		delta_x = x1-x2;
		//delta_x -= box_size * nearbyint(delta_x / box_size);
		
		delta_y = y1-y2;
		//delta_y -=  box_size * nearbyint(delta_y / box_size);

		delta_z = z1-z2;
		//delta_z -=  box_size * nearbyint(delta_z / box_size);
		sigma_ij = sigma_array[index];
		distance = sqrt((delta_x*delta_x) + (delta_y*delta_y) +(delta_z*delta_z));

		a = a_array[index];
		k = k_array[index];
		rb = rb_array[index];

		// Calculate magnitude of interaction based on distance
		if (distance < rb){
			mag = (1. - (distance/sigma_ij))/sigma_ij;
			if (distance < sigma_ij){
				*z += 1;
				u_rlj_energy += 0.5 * (1. - (distance/sigma_ij)) * (1. - (distance/sigma_ij));
			}
		} 
		else if (distance < a){
			//mag = k*(a - distance) - b;
			mag = k * ((distance/a) - 1.) / a;
			//U += 0.5 * k * (a - distance)*(a - distance) + (b*distance) - 0.5*(b*a + b*sigma_ij) - c;
		}
		else{
			mag = 0.0;
		}

		// Update Force
		// Magnitude of force times the direction of the force
		total_force[nonbonded_array1[index]][0] += (delta_x/distance) * mag;
		total_force[nonbonded_array1[index]][1] += (delta_y/distance) * mag;
		total_force[nonbonded_array1[index]][2] += (delta_z/distance) * mag;
		total_force[nonbonded_array2[index]][0] += -(delta_x/distance) * mag;
		total_force[nonbonded_array2[index]][1] += -(delta_y/distance) * mag;
		total_force[nonbonded_array2[index]][2] += -(delta_z/distance) * mag;

	}
	return u_rlj_energy;
}

// Bond length interaction
double compute_bonded(double coords[][3], int bonded_array1[], int bonded_array2[],
                    double bonded_avg_array[], int num_bonded, double total_force[][3]){
	double distance;
	double mag;
	double bonded_energy = 0;
	int i;
	double x1, x2, y1, y2, z1, z2, delta_x, delta_y, delta_z;
	double std_bonded;
	for (i = 0; i < num_bonded; i++){
		x1 = coords[bonded_array1[i]][0];
		y1 = coords[bonded_array1[i]][1];
		z1 = coords[bonded_array1[i]][2];
		x2 = coords[bonded_array2[i]][0];
		y2 = coords[bonded_array2[i]][1];
		z2 = coords[bonded_array2[i]][2];

		delta_x = x1-x2;
		delta_y = y1-y2;
		delta_z = z1-z2;

		distance = sqrt((delta_x*delta_x) + (delta_y*delta_y) +(delta_z*delta_z));
		
		mag = (distance - bonded_avg_array[i]);
		bonded_energy = 0.5*mag*mag;
		total_force[bonded_array1[i]][0] += -(delta_x/distance) * mag;
		total_force[bonded_array1[i]][1] += -(delta_y/distance) * mag;
		total_force[bonded_array1[i]][2] += -(delta_z/distance) * mag;
		total_force[bonded_array2[i]][0] += (delta_x/distance) * mag;
		total_force[bonded_array2[i]][1] += (delta_y/distance) * mag;
		total_force[bonded_array2[i]][2] += (delta_z/distance) * mag;
		
	}
	return bonded_energy;
}

//Bond angle interaction
double compute_angle(double coords[][3], int angle_array1[], int angle_array2[], int angle_array3[],
                     double angle_avg_array[], int num_angles, double total_force[][3]){
	double va_x, va_y, va_z, vb_x, vb_y, vb_z, vc_x, vc_y, vc_z;
	double vba_x, vba_y, vba_z;
	double vab_x, vab_y, vab_z;
	double vbc_x, vbc_y, vbc_z;
	double vcb_x, vcb_y, vcb_z;
	double vca_x, vca_y, vca_z;
	double lba, lbc;
	double measured_angle;
	double A, B, C, D, H, AB;
	double fa_x, fa_y, fa_z, fb_x, fb_y, fb_z, fc_x, fc_y, fc_z;
	double angle_energy = 0;
	int i;
	double Pzz_angle = 0;
	double std_angle;
	for (i = 0; i < num_angles; i++){
		va_x = coords[angle_array1[i]][0];
		va_y = coords[angle_array1[i]][1];
		va_z = coords[angle_array1[i]][2];
		vb_x = coords[angle_array2[i]][0];
		vb_y = coords[angle_array2[i]][1];
		vb_z = coords[angle_array2[i]][2];
		vc_x = coords[angle_array3[i]][0];
		vc_y = coords[angle_array3[i]][1];
		vc_z = coords[angle_array3[i]][2];

		vba_x = vb_x - va_x;
		vba_y = vb_y - va_y;
		vba_z = vb_z - va_z;

		vbc_x = vb_x - vc_x;
		vbc_y = vb_y - vc_y;
		vbc_z = vb_z - vc_z;

		vab_x = -vba_x;
		vab_y = -vba_y;
		vab_z = -vba_z;

		vcb_x = -vbc_x;
		vcb_y = -vbc_y;
		vcb_z = -vbc_z;

		lba = sqrt(vba_x*vba_x + vba_y*vba_y + vba_z*vba_z);
		lbc = sqrt(vbc_x*vbc_x + vbc_y*vbc_y + vbc_z*vbc_z);

		measured_angle = acos((vba_x*vbc_x + vba_y*vbc_y + vba_z*vbc_z) / (lba * lbc));
			
		A = (measured_angle - angle_avg_array[i]);
	
		angle_energy += 0.5 * A * A;

		B = lba * lbc * fabs(sin(measured_angle));

		C = vab_x*vcb_x + vab_y*vcb_y + vab_z*vcb_z;

		D = pow_int(lba,2);

		H = pow_int(lbc,2);
		
		AB = A / B;

		fa_x = AB * (vcb_x - (C / D)*vab_x);
		fa_y = AB * (vcb_y - (C / D)*vab_y);
		fa_z = AB * (vcb_z - (C / D)*vab_z);

		fc_x = AB * (vab_x - (C / H)*vcb_x);
		fc_y = AB * (vab_y - (C / H)*vcb_y);
		fc_z = AB * (vab_z - (C / H)*vcb_z);

		fb_x = -fa_x-fc_x;
		fb_y = -fa_y-fc_y;
		fb_z = -fa_z-fc_z;

		total_force[angle_array1[i]][0] += fa_x;
		total_force[angle_array1[i]][1] += fa_y;
		total_force[angle_array1[i]][2] += fa_z;
		total_force[angle_array2[i]][0] += fb_x;
		total_force[angle_array2[i]][1] += fb_y;
		total_force[angle_array2[i]][2] += fb_z;
		total_force[angle_array3[i]][0] += fc_x;
		total_force[angle_array3[i]][1] += fc_y;
		total_force[angle_array3[i]][2] += fc_z;
		
	}
	return angle_energy;
}

// Dihedral angle interaction
double compute_omega(double coords[][3], int CA1_omega_array[], int C_omega_array[],
                     int N_omega_array[], int CA2_omega_array[], double omega_array[], 
                     int num_omega, double total_force[][3]){
	double omega_energy = 0;
	double omega;
	double CA1_x, CA1_y, CA1_z;
	double C_x, C_y, C_z;
	double N_x, N_y, N_z;
	double CA2_x, CA2_y, CA2_z;
	double r_ij_x, r_ij_y, r_ij_z;
	double r_jk_x, r_jk_y, r_jk_z;
	double v_3_x, v_3_y, v_3_z;
	double cross_x, cross_y, cross_z, len_cross;
	double norm_1_x, norm_1_y, norm_1_z;
	double norm_2_x, norm_2_y, norm_2_z;
	double v_2_unit_x, v_2_unit_y, v_2_unit_z;
	double len_r_jk;
	double ortho_u_v_x, ortho_u_v_y, ortho_u_v_z;
	double sin_theta, cos_theta;

	//Get rid of me
	double r_kl_x, r_kl_y, r_kl_z;
	//

	double r_mj_x, r_mj_y, r_mj_z;
	double r_nk_x, r_nk_y, r_nk_z;
	double r_kj_len;
	double A, B, C, D;
	double f_CA1_x, f_CA1_y, f_CA1_z;
	double f_CA2_x, f_CA2_y, f_CA2_z;
	double f_C_x, f_C_y, f_C_z;
	double f_N_x, f_N_y, f_N_z;
	double mag, len;
	double Pzz_omega = 0;
	int i;

	for (i = 0; i < num_omega; i++){
		CA1_x = coords[CA1_omega_array[i]][0];
		CA1_y = coords[CA1_omega_array[i]][1];
		CA1_z = coords[CA1_omega_array[i]][2];
	    C_x = coords[C_omega_array[i]][0];
	    C_y = coords[C_omega_array[i]][1];
	    C_z = coords[C_omega_array[i]][2];
	    N_x = coords[N_omega_array[i]][0];
	    N_y = coords[N_omega_array[i]][1];
	    N_z = coords[N_omega_array[i]][2];
	    CA2_x = coords[CA2_omega_array[i]][0];
	    CA2_y = coords[CA2_omega_array[i]][1];
	    CA2_z = coords[CA2_omega_array[i]][2];

	    r_ij_x = CA1_x - C_x;
	    r_ij_y = CA1_y - C_y;
	    r_ij_z = CA1_z - C_z;

	    r_jk_x = C_x - N_x;
	    r_jk_y = C_y - N_y;
	    r_jk_z = C_z - N_z;

	    v_3_x = -(N_x - CA2_x);
	    v_3_y = -(N_y - CA2_y);
	    v_3_z = -(N_z - CA2_z);

	    //r_ij X r_jk
	    cross_x = r_ij_y * r_jk_z - r_jk_y * r_ij_z;
	    cross_y = r_jk_x * r_ij_z - r_ij_x * r_jk_z;
	    cross_z = r_ij_x * r_jk_y - r_jk_x * r_ij_y;

	    len_cross = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);

	    norm_1_x = cross_x / len_cross;
	    norm_1_y = cross_y / len_cross;
	    norm_1_z = cross_z / len_cross;

	    //r_jk X v_3
	    cross_x = r_jk_y * v_3_z - v_3_y * r_jk_z;
	    cross_y = v_3_x * r_jk_z - r_jk_x * v_3_z;
	    cross_z = r_jk_x * v_3_y - v_3_x * r_jk_y;

	    len_cross = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);

	    norm_2_x = cross_x / len_cross;
	    norm_2_y = cross_y / len_cross;
	    norm_2_z = cross_z / len_cross;

	    len_r_jk = sqrt(r_jk_x*r_jk_x + r_jk_y*r_jk_y + r_jk_z*r_jk_z);
	    v_2_unit_x = r_jk_x / len_r_jk;
	    v_2_unit_y = r_jk_y / len_r_jk;
	    v_2_unit_z = r_jk_z / len_r_jk;

	    //v_2_unit X norm_2
	    ortho_u_v_x = v_2_unit_y * norm_2_z - norm_2_y * v_2_unit_z;
	    ortho_u_v_y = norm_2_x * v_2_unit_z - v_2_unit_x * norm_2_z;
	    ortho_u_v_z = v_2_unit_x * norm_2_y - norm_2_x * v_2_unit_y;
	    
	    cos_theta = norm_1_x*norm_2_x + norm_1_y*norm_2_y + norm_1_z*norm_2_z;
	    sin_theta = norm_1_x*ortho_u_v_x + norm_1_y*ortho_u_v_y + norm_1_z*ortho_u_v_z;

	    omega = atan2(sin_theta, cos_theta);
	    if (omega > pi / 2){
	    	omega -= pi;
	    }
	    else if (omega < - pi / 2){
	    	omega += pi;
	    }
	    omega -= omega_array[i];

    	omega_energy += 0.5 * (omega)*(omega);

    	r_kl_x = -v_3_x;
    	r_kl_y = -v_3_y;
    	r_kl_z = -v_3_z;

    	//r_ij X r_kj
	    r_mj_x = r_ij_y * -r_jk_z - -r_jk_y * r_ij_z;
	    r_mj_y = -r_jk_x * r_ij_z - r_ij_x * -r_jk_z;
	    r_mj_z = r_ij_x * -r_jk_y - -r_jk_x * r_ij_y;

	    //r_kj X r_kl
	    r_nk_x = -r_jk_y * r_kl_z - r_kl_y * -r_jk_z;
	    r_nk_y = r_kl_x * -r_jk_z - -r_jk_x * r_kl_z;
	    r_nk_z = -r_jk_x * r_kl_y - r_kl_x * -r_jk_y;

	    r_kj_len = sqrt(r_jk_x*r_jk_x + r_jk_y*r_jk_y + r_jk_z*r_jk_z);

	    A = r_kj_len / (r_mj_x*r_mj_x + r_mj_y*r_mj_y + r_mj_z*r_mj_z);

	    f_CA1_x = A * r_mj_x;
	    f_CA1_y = A * r_mj_y;
	    f_CA1_z = A * r_mj_z;

	    B = r_kj_len / (r_nk_x*r_nk_x + r_nk_y*r_nk_y + r_nk_z*r_nk_z);

	    f_CA2_x = -B * r_nk_x;
	    f_CA2_y = -B * r_nk_y;
	    f_CA2_z = -B * r_nk_z;

	    C = -(r_ij_x*r_jk_x + r_ij_y*r_jk_y + r_ij_z*r_jk_z) / (r_kj_len*r_kj_len);

	    D = -(r_kl_x*r_jk_x + r_kl_y*r_jk_y + r_kl_z*r_jk_z) / (r_kj_len*r_kj_len);

	    f_C_x = ((C - 1) * f_CA1_x) - (D * f_CA2_x);
	    f_C_y = ((C - 1) * f_CA1_y) - (D * f_CA2_y);
	    f_C_z = ((C - 1) * f_CA1_z) - (D * f_CA2_z);

	    f_N_x = ((D - 1) * f_CA2_x) - (C * f_CA1_x);
	    f_N_y = ((D - 1) * f_CA2_y) - (C * f_CA1_y);
	    f_N_z = ((D - 1) * f_CA2_z) - (C * f_CA1_z);


	    mag = (omega);


	    // F_CA1
		f_CA1_x *= -mag;
		f_CA1_y *= -mag;
		f_CA1_z *= -mag;

		// F_C
		f_C_x *= -mag;
		f_C_y *= -mag;
		f_C_z *= -mag;


		// F_N
		f_N_x *= -mag;
		f_N_y *= -mag;
		f_N_z *= -mag;


		// F_CA2
		f_CA2_x *= -mag;
		f_CA2_y *= -mag;
		f_CA2_z *= -mag;

	    total_force[CA1_omega_array[i]][0] += f_CA1_x;
	    total_force[CA1_omega_array[i]][1] += f_CA1_y;
	    total_force[CA1_omega_array[i]][2] += f_CA1_z;

	    total_force[C_omega_array[i]][0] += f_C_x;
	    total_force[C_omega_array[i]][1] += f_C_y;
	    total_force[C_omega_array[i]][2] += f_C_z;

	    total_force[N_omega_array[i]][0] += f_N_x;
	    total_force[N_omega_array[i]][1] += f_N_y;
	    total_force[N_omega_array[i]][2] += f_N_z;

	    total_force[CA2_omega_array[i]][0] += f_CA2_x;
	    total_force[CA2_omega_array[i]][1] += f_CA2_y;
	    total_force[CA2_omega_array[i]][2] += f_CA2_z;

	}

	return omega_energy;
}

// Function for energy minimized via damping (linear drag)
void compute_damping(double velocs[][3], double total_force[][3], double prev_total_force[][3], int num_atoms, double dt, double damping){

	int i, j;
	double before = 0, after = 0;
	double force_dot_velocs = 0;
	for (i=0; i<num_atoms; i++){
		for (j=0; j<3; j++){
			total_force[i][j] -= damping*velocs[i][j] - (0.5*damping*prev_total_force[i][j]*dt);
			total_force[i][j] /=  (1+0.5*damping*dt);
		}
	}
	return;
}

// Draw gaussian numbers
void generate_gaussian(int half_sample_size, double gaussian_array[][3]){
	int i, dim;
	double rand_unif1, rand_unif2;
	double rand_gauss1, rand_gauss2;
	double twopi = 6.28318530718;
	double A, B;
	for (dim = 0; dim < 3; dim++){
		for (i = 0; i < half_sample_size; i++){
			rand_unif1 = ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
			rand_unif2 = ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);

			A = sqrt(-2.*log(rand_unif1));
			B = twopi * rand_unif2;
			rand_gauss1 = A*cos(B);
			rand_gauss2 = A*sin(B);

			gaussian_array[i*2][dim] = rand_gauss1;
			gaussian_array[(i*2)+1][dim] = rand_gauss2;
		}
	}
}

// Create the Verlet list to know who to check for non-bonded interaction
int generate_vlist(double coords[][3], int nonbonded_array1[], int nonbonded_array2[], double sigma_array[],
					int num_nonbonded, int vlist[], double r_l, double a_array[]){
	int i;
	int num_neighs = 0;
	double x1, x2, y1, y2, z1, z2, delta_x, delta_y, delta_z;
	double distance, rc, rl;
	
	for (i=0; i<num_nonbonded; i++){
		x1 = coords[nonbonded_array1[i]][0];
		y1 = coords[nonbonded_array1[i]][1];
		z1 = coords[nonbonded_array1[i]][2];
		x2 = coords[nonbonded_array2[i]][0];
		y2 = coords[nonbonded_array2[i]][1];
		z2 = coords[nonbonded_array2[i]][2];

		delta_x = x1-x2;
		//delta_x -= box_size * nearbyint(delta_x / box_size);
		
		delta_y = y1-y2;
		//delta_y -=  box_size * nearbyint(delta_y / box_size);

		delta_z = z1-z2;
		//delta_z -=  box_size * nearbyint(delta_z / box_size);

		distance = sqrt((delta_x*delta_x) + (delta_y*delta_y) +(delta_z*delta_z));
		
		if (distance < (a_array[i]) * r_l){
			vlist[num_neighs] = i;
			num_neighs += 1;
		}
	}
	return num_neighs;
}

// Set the geometry restraints to be centered at the current value
void reset_bonds(double coords[][3], int bonded_array1[], int bonded_array2[],
                    double bonded_avg_array[], int num_bonded){
	double distance;
	int i;
	double x1, x2, y1, y2, z1, z2, delta_x, delta_y, delta_z;
	for (i = 0; i < num_bonded; i++){
		x1 = coords[bonded_array1[i]][0];
		y1 = coords[bonded_array1[i]][1];
		z1 = coords[bonded_array1[i]][2];
		x2 = coords[bonded_array2[i]][0];
		y2 = coords[bonded_array2[i]][1];
		z2 = coords[bonded_array2[i]][2];

		delta_x = x1-x2;
		delta_y = y1-y2;
		delta_z = z1-z2;

		distance = sqrt((delta_x*delta_x) + (delta_y*delta_y) +(delta_z*delta_z));
		bonded_avg_array[i] = distance;
	}
}

void reset_angles(double coords[][3], int angle_array1[], int angle_array2[], int angle_array3[],
                     double angle_avg_array[], int num_angles){
	double va_x, va_y, va_z, vb_x, vb_y, vb_z, vc_x, vc_y, vc_z;
	double vba_x, vba_y, vba_z;
	double vab_x, vab_y, vab_z;
	double vbc_x, vbc_y, vbc_z;
	double vcb_x, vcb_y, vcb_z;
	double vca_x, vca_y, vca_z;
	double lba, lbc;
	double measured_angle;
	int i;
	for (i = 0; i < num_angles; i++){
		va_x = coords[angle_array1[i]][0];
		va_y = coords[angle_array1[i]][1];
		va_z = coords[angle_array1[i]][2];
		vb_x = coords[angle_array2[i]][0];
		vb_y = coords[angle_array2[i]][1];
		vb_z = coords[angle_array2[i]][2];
		vc_x = coords[angle_array3[i]][0];
		vc_y = coords[angle_array3[i]][1];
		vc_z = coords[angle_array3[i]][2];

		vba_x = vb_x - va_x;
		vba_y = vb_y - va_y;
		vba_z = vb_z - va_z;

		vbc_x = vb_x - vc_x;
		vbc_y = vb_y - vc_y;
		vbc_z = vb_z - vc_z;

		vab_x = -vba_x;
		vab_y = -vba_y;
		vab_z = -vba_z;

		vcb_x = -vbc_x;
		vcb_y = -vbc_y;
		vcb_z = -vbc_z;

		lba = sqrt(vba_x*vba_x + vba_y*vba_y + vba_z*vba_z);
		lbc = sqrt(vbc_x*vbc_x + vbc_y*vbc_y + vbc_z*vbc_z);

		measured_angle = acos((vba_x*vbc_x + vba_y*vbc_y + vba_z*vbc_z) / (lba * lbc));
		angle_avg_array[i] = measured_angle;
	}
}
void reset_omega(double coords[][3], int CA1_omega_array[], int C_omega_array[],
                     int N_omega_array[], int CA2_omega_array[], double omega_array[], int num_omega){
	double omega_energy = 0;
	double omega;
	double CA1_x, CA1_y, CA1_z;
	double C_x, C_y, C_z;
	double N_x, N_y, N_z;
	double CA2_x, CA2_y, CA2_z;
	double r_ij_x, r_ij_y, r_ij_z;
	double r_jk_x, r_jk_y, r_jk_z;
	double v_3_x, v_3_y, v_3_z;
	double cross_x, cross_y, cross_z, len_cross;
	double norm_1_x, norm_1_y, norm_1_z;
	double norm_2_x, norm_2_y, norm_2_z;
	double v_2_unit_x, v_2_unit_y, v_2_unit_z;
	double len_r_jk;
	double ortho_u_v_x, ortho_u_v_y, ortho_u_v_z;
	double sin_theta, cos_theta;

	//Get rid of me
	double r_kl_x, r_kl_y, r_kl_z;
	//

	double r_mj_x, r_mj_y, r_mj_z;
	double r_nk_x, r_nk_y, r_nk_z;
	double r_kj_len;
	double A, B, C, D;
	double f_CA1_x, f_CA1_y, f_CA1_z;
	double f_CA2_x, f_CA2_y, f_CA2_z;
	double f_C_x, f_C_y, f_C_z;
	double f_N_x, f_N_y, f_N_z;
	double mag, len;
	double Pzz_omega = 0;
	int i;

	double std_omega = 0.0;//2.0*sqrt(0.01 / k_omega);
	for (i = 0; i < num_omega; i++){
		CA1_x = coords[CA1_omega_array[i]][0];
		CA1_y = coords[CA1_omega_array[i]][1];
		CA1_z = coords[CA1_omega_array[i]][2];
	    C_x = coords[C_omega_array[i]][0];
	    C_y = coords[C_omega_array[i]][1];
	    C_z = coords[C_omega_array[i]][2];
	    N_x = coords[N_omega_array[i]][0];
	    N_y = coords[N_omega_array[i]][1];
	    N_z = coords[N_omega_array[i]][2];
	    CA2_x = coords[CA2_omega_array[i]][0];
	    CA2_y = coords[CA2_omega_array[i]][1];
	    CA2_z = coords[CA2_omega_array[i]][2];

	    r_ij_x = CA1_x - C_x;
	    r_ij_y = CA1_y - C_y;
	    r_ij_z = CA1_z - C_z;

	    r_jk_x = C_x - N_x;
	    r_jk_y = C_y - N_y;
	    r_jk_z = C_z - N_z;

	    v_3_x = -(N_x - CA2_x);
	    v_3_y = -(N_y - CA2_y);
	    v_3_z = -(N_z - CA2_z);

	    //r_ij X r_jk
	    cross_x = r_ij_y * r_jk_z - r_jk_y * r_ij_z;
	    cross_y = r_jk_x * r_ij_z - r_ij_x * r_jk_z;
	    cross_z = r_ij_x * r_jk_y - r_jk_x * r_ij_y;

	    len_cross = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);

	    norm_1_x = cross_x / len_cross;
	    norm_1_y = cross_y / len_cross;
	    norm_1_z = cross_z / len_cross;

	    //r_jk X v_3
	    cross_x = r_jk_y * v_3_z - v_3_y * r_jk_z;
	    cross_y = v_3_x * r_jk_z - r_jk_x * v_3_z;
	    cross_z = r_jk_x * v_3_y - v_3_x * r_jk_y;

	    len_cross = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);

	    norm_2_x = cross_x / len_cross;
	    norm_2_y = cross_y / len_cross;
	    norm_2_z = cross_z / len_cross;

	    len_r_jk = sqrt(r_jk_x*r_jk_x + r_jk_y*r_jk_y + r_jk_z*r_jk_z);
	    v_2_unit_x = r_jk_x / len_r_jk;
	    v_2_unit_y = r_jk_y / len_r_jk;
	    v_2_unit_z = r_jk_z / len_r_jk;

	    //v_2_unit X norm_2
	    ortho_u_v_x = v_2_unit_y * norm_2_z - norm_2_y * v_2_unit_z;
	    ortho_u_v_y = norm_2_x * v_2_unit_z - v_2_unit_x * norm_2_z;
	    ortho_u_v_z = v_2_unit_x * norm_2_y - norm_2_x * v_2_unit_y;
	    
	    cos_theta = norm_1_x*norm_2_x + norm_1_y*norm_2_y + norm_1_z*norm_2_z;
	    sin_theta = norm_1_x*ortho_u_v_x + norm_1_y*ortho_u_v_y + norm_1_z*ortho_u_v_z;

	    omega = atan2(sin_theta, cos_theta);
	    if (omega > pi / 2){
	    	omega -= pi;
	    }
	    else if (omega < - pi / 2){
	    	omega += pi;
	    }
	    omega_array[i] = omega;

	}    
}

// Main
// General organization:
// -Load the coordinates and all the interaction parameters
// -Energy minimize slightly
// -Run constant temperature dynamics via Langevin
int main(int argc, char *argv[]) {
	
	remove("traj_lang.xyz");
	remove("traj_comp.xyz");
	remove("fcheck.txt");
	
	printf("Welcome!\n");
	srand(time(0));
	double total_potential_energy = 0;
	double lj_energy, bonded_energy, angle_energy, omega_energy, planar_energy, improper_energy;
	double total_kinetic_energy = 0;
	double temp, rescale_factor;
	double target_temp = 1e-6;
	double v;
	//double avg_x, avg_y, avg_z, x_box, y_box, z_box;
	clock_t time_1, time_2;
	int i, j, k;
	int num_neighs;
	double vstress[3] = {0};
	double x_box, y_box, z_box;
	double rl, rc;
	double damping = 0.1;
	double damping_final = atof(argv[2]);
	double dt = 0.05;
	double dt_lang = 0.05;
	double dt_sqr = dt*dt;
	double r_l = 1.1;

	int z = 0;

	double box_size = 0;

	// Counting number of atoms
	FILE *fp;
	int num_atoms = 0;  
    char c;
    char buf[0x100];
	//snprintf(buf, sizeof(buf), "params/%s_params/coords.txt",argv[1]);
	snprintf(buf, sizeof(buf), "em_xtal_coords/coords_%s.txt",argv[1]);
    fp = fopen(buf, "r"); 
    // Check if file exists 
    if (fp == NULL) 
    { 
        printf("Could not open file %s", "coords.txt"); 
        return 0; 
    } 
    // Extract characters from file and store in character c 
    for (c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n') 
            num_atoms = num_atoms + 1; 
    // Close the file 
    fclose(fp); 
    printf("Num_atoms = %d\n",num_atoms);
    int half_num_atoms = num_atoms / 2;
	/* Loading Initial Coords */
	fp = fopen(buf, "r"); 
	double coords[num_atoms][3], previous_coords[num_atoms][3], half_coords[num_atoms][3], displacement_array[num_atoms][3];
	double total_force[num_atoms][3], lj_force[num_atoms][3], bonded_force[num_atoms][3], angle_force[num_atoms][3], omega_force[num_atoms][3], planar_force[num_atoms][3], improper_force[num_atoms][3];
	double prev_total_force[num_atoms][3];
	double velocs[num_atoms][3], half_velocs[num_atoms][3], v_prime[num_atoms][3];
	double gaussian_array[num_atoms][3];
	double sasa_atom_array[num_atoms], voro_atom_array[num_atoms], vol_atom_array[num_atoms];
	int resid_array[num_atoms];
	double trash_double;

	for (i=0; i<num_atoms; i++){
		for (j=0; j<3; j++){
			fscanf(fp," %lf",&coords[i][j]);
		}
	}
	fclose(fp);
	printf("Coords loaded!\n");

	//Centering Coords
	double avg_x = 0;
	double avg_y = 0;
	double avg_z = 0;
	for (i=0; i<num_atoms; i++){
		avg_x += coords[i][0];
		avg_y += coords[i][1];
		avg_z += coords[i][2];
	}
	avg_x /= num_atoms;
	avg_y /= num_atoms;
	avg_z /= num_atoms;
	for (i=0; i<num_atoms; i++){
		coords[i][0] -= avg_x;
		coords[i][1] -= avg_y;
		coords[i][2] -= avg_z;
	}
	/* Loading Initial Velocities */
	for (i=0; i<num_atoms; i++){
		for (j=0; j<3; j++){
			velocs[i][j] = 1e-15;
		}
	}


	snprintf(buf, sizeof(buf), "params/%s_params/sigma_i_array.txt",argv[1]);
	fp = fopen(buf, "r"); 
	double sigma_i_array[num_atoms];
	for (i=0; i<num_atoms; i++){
		fscanf(fp," %lf",&sigma_i_array[i]);
	}
	fclose(fp); 


//////////

	// Loading atomic masses //
	//snprintf(buf, sizeof(buf), "params/%s_params/mass_array.txt",argv[1]);
	//fp = fopen(buf, "r"); 
	double mass_array[num_atoms];
	for (i=0; i<num_atoms; i++){
		//fscanf(fp," %lf",&mass_array[i]);
		mass_array[i] = 1.0;
	}
	//fclose(fp); 

//////////
	
	// Loading Atom Types in 2D Array // 
	char str;
	char atom_array[num_atoms][10];
	//A 2D Array with a buffered col 
	int atom_len_array[num_atoms];
	//Saving the number of characters in each atom type
	int atom_len = 0;
	int count = 0;
	int col = 0;
	char item[2];
	snprintf(buf, sizeof(buf), "params/%s_params/atom_array.txt",argv[1]);
	fp = fopen(buf, "r"); 
	i = 0;
	for (str = getc(fp); str != EOF; str = getc(fp))
        if (str != '\n'){
        	atom_array[i][col] = str;
			count++;
			col++;
			atom_len++;
        } 
		else{
			atom_len_array[i] = atom_len;
			col = 0;
			count++;
			i++;
			atom_len = 0;
		}
	printf("atom names loaded!\n");
//////////

	// Loading Nonbonded Sigma Pairs //
	// Count number of nonboned pairs
	int num_nonbonded = 0; 
	snprintf(buf, sizeof(buf), "params/%s_params/sigma_array.txt",argv[1]);
	fp = fopen(buf, "r"); 
    for (c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n')
            num_nonbonded = num_nonbonded + 1; 
    fclose(fp);
    printf("Num_nonbonded = %d\n",num_nonbonded);
	// Loading //
	snprintf(buf, sizeof(buf), "params/%s_params/sigma_array.txt",argv[1]);
	fp = fopen(buf, "r"); 
	double* sigma_array;
	sigma_array = (double*)malloc(num_nonbonded * sizeof(double));
	for (i=0; i<num_nonbonded; i++){
		fscanf(fp," %lf",&sigma_array[i]);
	}
	fclose(fp);

	// Loading //
	snprintf(buf, sizeof(buf), "params/%s_params/epsilon_ij_array.txt",argv[1]);
	fp = fopen(buf, "r"); 
	double* epsilon_ij_array;
	epsilon_ij_array = (double*)malloc(num_nonbonded * sizeof(double));
	for (i=0; i<num_nonbonded; i++){
		fscanf(fp," %lf",&epsilon_ij_array[i]);
	}
	fclose(fp);


/////////

	printf("Loading non-bonded...\n");
	// Loading Nonbonded Index Pairs //
	snprintf(buf, sizeof(buf), "params/%s_params/nonbonded_array1.txt",argv[1]);
    fp = fopen(buf, "r"); 

	int* nonbonded_array1;
	nonbonded_array1 = (int*)malloc(num_nonbonded * sizeof(int));
	for (i=0; i<num_nonbonded; i++){
		fscanf(fp," %d",&nonbonded_array1[i]);
	}
	fclose(fp); 
	// Loading //
	snprintf(buf, sizeof(buf), "params/%s_params/nonbonded_array2.txt",argv[1]);
    fp = fopen(buf, "r"); 
	int* nonbonded_array2;
	nonbonded_array2 = (int*)malloc(num_nonbonded * sizeof(int));
	for (i=0; i<num_nonbonded; i++){
		fscanf(fp," %d",&nonbonded_array2[i]);
	}
	fclose(fp); 

	// Initialize vlist //
	int* vlist;
	vlist = (int*)malloc(num_nonbonded * sizeof(int));

	printf("Done!\n");
	// Setting Non-bonded parameters
	int* inter_intra_array;
	double* a_array;
	double alpha = atof(argv[2]);
	double beta = atof(argv[3]);
	double* k_array;
	double* c_array;
	double* rb_array;

	a_array = (double*)malloc(num_nonbonded * sizeof(double));
	k_array = (double*)malloc(num_nonbonded * sizeof(double));
	rb_array = (double*)malloc(num_nonbonded * sizeof(double));
	double epsilon_ij;
	for(i=0; i<num_nonbonded; i++){
		epsilon_ij = epsilon_ij_array[i];
		a_array[i] = sigma_array[i] + alpha*sigma_array[i];
		rb_array[i] = sigma_array[i]*(1. + (sigma_array[i]*beta*epsilon_ij));
		k_array[i] = -(a_array[i]*beta*epsilon_ij) / ((rb_array[i]/a_array[i]) - 1.); 
	}
	printf("Nonbonded loaded!\n");

/////////

	printf("Loading bonded...\n");
	// Loading Bonded Index Pairs //
	// Counting number of Bonded Pairs
	int num_bonded = 0;
	snprintf(buf, sizeof(buf), "params/%s_params/bonded_array1.txt",argv[1]);
    fp = fopen(buf, "r"); 
    for (c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n') 
            num_bonded = num_bonded + 1; 
    fclose(fp);

	snprintf(buf, sizeof(buf), "params/%s_params/bonded_array1.txt",argv[1]);
    fp = fopen(buf, "r"); 
	int bonded_array1[num_bonded];
	for (i=0; i<num_bonded; i++){
		fscanf(fp," %d",&bonded_array1[i]);
	}
	fclose(fp); 

	snprintf(buf, sizeof(buf), "params/%s_params/bonded_array2.txt",argv[1]);
    fp = fopen(buf, "r"); 
	int bonded_array2[num_bonded];
	for (i=0; i<num_bonded; i++){
		fscanf(fp," %d",&bonded_array2[i]);
	}
	fclose(fp); 

/////////

	// Loading Bonded stats //
	double bonded_avg_array[num_bonded];

/////////

	// Loading Angle Index Triples //
	//Counting the number of angles
	int num_angles = 0;  // Line counter (result) 
	snprintf(buf, sizeof(buf), "params/%s_params/angle_array1.txt",argv[1]);
    fp = fopen(buf, "r"); 
    for (c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n')
            num_angles = num_angles + 1; 
    fclose(fp); 
	
	snprintf(buf, sizeof(buf), "params/%s_params/angle_array1.txt",argv[1]);
    fp = fopen(buf, "r");
	int angle_array1[num_angles];
	for (i=0; i<num_angles; i++){
		fscanf(fp," %d",&angle_array1[i]);
	}
	fclose(fp); 
	
	snprintf(buf, sizeof(buf), "params/%s_params/angle_array2.txt",argv[1]);
    fp = fopen(buf, "r");
	int angle_array2[num_angles];
	for (i=0; i<num_angles; i++){
		fscanf(fp," %d",&angle_array2[i]);
	}
	fclose(fp); 
	
	snprintf(buf, sizeof(buf), "params/%s_params/angle_array3.txt",argv[1]);
    fp = fopen(buf, "r");
	int angle_array3[num_angles];
	for (i=0; i<num_angles; i++){
		fscanf(fp," %d",&angle_array3[i]);
	}

/////////

	// Loading Angle Stats //
	double angle_avg_array[num_angles];	

/////////

	// Loading Dihedral Omega Indexs //
	// Counting number of omega dihedrals
	snprintf(buf, sizeof(buf), "params/%s_params/CA1_omega_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int num_omega = 0;
    for (c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n') 
            num_omega = num_omega + 1; 
    fclose(fp); 
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/CA1_omega_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int CA1_omega_array[num_omega];
	for (i=0; i<num_omega; i++){
		fscanf(fp," %d",&CA1_omega_array[i]);
	}
	fclose(fp);
	// Loading //
	snprintf(buf, sizeof(buf), "params/%s_params/C_omega_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int C_omega_array[num_omega];
	for (i=0; i<num_omega; i++){
		fscanf(fp," %d",&C_omega_array[i]);
	}
	fclose(fp);
	// Loading //
	snprintf(buf, sizeof(buf), "params/%s_params/N_omega_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int N_omega_array[num_omega];
	for (i=0; i<num_omega; i++){
		fscanf(fp," %d",&N_omega_array[i]);
	}
	fclose(fp);
	// Loading //
	snprintf(buf, sizeof(buf), "params/%s_params/CA2_omega_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int CA2_omega_array[num_omega];
	for (i=0; i<num_omega; i++){
		fscanf(fp," %d",&CA2_omega_array[i]);
	}
	fclose(fp);

/////////

	// Loading Dihedral stats //
	double omega_array[num_omega];

/////////

	// Loading Planar Dihedral Indexs //
	// Counting number of planar dihedrals
	snprintf(buf, sizeof(buf), "params/%s_params/A_planar_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int num_planar = 0;
    for (c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n') 
            num_planar = num_planar + 1; 
    fclose(fp); 
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/A_planar_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int A_planar_array[num_planar];
	for (i=0; i<num_planar; i++){
		fscanf(fp," %d",&A_planar_array[i]);
	}
	fclose(fp);
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/B_planar_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int B_planar_array[num_planar];
	for (i=0; i<num_planar; i++){
		fscanf(fp," %d",&B_planar_array[i]);
	}
	fclose(fp);
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/C_planar_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int C_planar_array[num_planar];
	for (i=0; i<num_planar; i++){
		fscanf(fp," %d",&C_planar_array[i]);
	}
	fclose(fp);
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/D_planar_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int D_planar_array[num_planar];
	for (i=0; i<num_planar; i++){
		fscanf(fp," %d",&D_planar_array[i]);
	}
	fclose(fp);

	double planar_array[num_planar];
	for (i=0; i<num_planar; i++){
		planar_array[i] = 0.0;
	}

/////////

	// Loading Improper Dihedral Indexs //
	// Counting number of improper dihedrals
	snprintf(buf, sizeof(buf), "params/%s_params/A_improper_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int num_improper = 0;
    for (c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n') 
            num_improper = num_improper + 1; 
    fclose(fp); 
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/A_improper_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int A_improper_array[num_improper];
	for (i=0; i<num_improper; i++){
		fscanf(fp," %d",&A_improper_array[i]);
	}
	fclose(fp);
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/B_improper_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int B_improper_array[num_improper];
	for (i=0; i<num_improper; i++){
		fscanf(fp," %d",&B_improper_array[i]);
	}
	fclose(fp);
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/C_improper_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int C_improper_array[num_improper];
	for (i=0; i<num_improper; i++){
		fscanf(fp," %d",&C_improper_array[i]);
	}
	fclose(fp);
	// Loading //
    snprintf(buf, sizeof(buf), "params/%s_params/D_improper_array.txt",argv[1]);
    fp = fopen(buf, "r");
	int D_improper_array[num_improper];
	for (i=0; i<num_improper; i++){
		fscanf(fp," %d",&D_improper_array[i]);
	}
	fclose(fp);

	double improper_array[num_improper];
	for (i=0; i<num_improper; i++){
		improper_array[i] = 0.0;
	}

	printf("Params loaded!\n");

	double Ftol = 1e-5;    // force tolerance

	double fcheck = 10*Ftol;

    int fireit = -1;

    int npPos = 0;
    int npNeg = 0;
    int npPMin = 0;
    
	double vnorm = 0;
	double fnorm = 0;
	double P = 0;
	double disp;
	double disp_x, disp_y, disp_z;
	int disp_threshold_count = 0;


	memset(total_force, 0, sizeof(total_force));

	reset_bonds(coords, bonded_array1, bonded_array2, bonded_avg_array, num_bonded);
	reset_angles(coords, angle_array1, angle_array2, angle_array3, angle_avg_array, num_angles);
	reset_omega(coords, CA1_omega_array, C_omega_array, N_omega_array, CA2_omega_array, omega_array, num_omega);

	// Compute new Force

	double disp_list[num_atoms];
	double disp_max1, disp_max2;
	double min_rc = 1.0+alpha;

	memset(total_force, 0, sizeof(total_force));
    memcpy(displacement_array, coords, sizeof(displacement_array));
	num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, sigma_array, num_nonbonded, vlist, r_l, a_array);
	lj_energy = compute_urlj(coords, nonbonded_array1, nonbonded_array2, sigma_array, a_array, k_array, rb_array, total_force, num_neighs, vlist, &z);
	bonded_energy = compute_bonded(coords, bonded_array1, bonded_array2, bonded_avg_array, num_bonded, total_force);
	angle_energy = compute_angle(coords, angle_array1, angle_array2, angle_array3, angle_avg_array, num_angles, total_force);
	omega_energy = compute_omega(coords, CA1_omega_array, C_omega_array, N_omega_array, CA2_omega_array, omega_array, num_omega, total_force);
	planar_energy = compute_omega(coords, A_planar_array, B_planar_array, D_planar_array, C_planar_array, planar_array, num_planar, total_force);
	improper_energy = compute_omega(coords, B_improper_array, A_improper_array, C_improper_array, D_improper_array, improper_array, num_improper, total_force);
	compute_damping(velocs, total_force, prev_total_force, num_atoms, dt, damping);

	memcpy(prev_total_force, total_force, sizeof(prev_total_force));

	fcheck = 0.0;
	for (j = 0; j < num_atoms; j++){
		fcheck += sqrt(total_force[j][0]*total_force[j][0] + total_force[j][1]*total_force[j][1] + total_force[j][2]*total_force[j][2]);
	}
	fcheck = fcheck/num_atoms;
	printf("initial fcheck = %.8e \n", fcheck);

    fcheck = 10*Ftol;
    fireit = -1;

    
    printf("Initial EM...\n");
    while ((fcheck > Ftol)) {
        fireit += 1;

        // Update Coordinates 
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 3; k++){
				coords[j][k] += dt * velocs[j][k] + 0.5*dt_sqr * (prev_total_force[j][k]);
			}
		}

        // Compute new Force
        memset(total_force, 0, sizeof(total_force));

        // Checking Disp for update on Verlet List 
		disp_threshold_count = 0;
		for (j=0; j<num_atoms; j++){
			disp_x = displacement_array[j][0]-coords[j][0];
			disp_y = displacement_array[j][1]-coords[j][1];
			disp_z = displacement_array[j][2]-coords[j][2];
			disp = sqrt(disp_x*disp_x + disp_y*disp_y + disp_z*disp_z);
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
			num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, sigma_array, num_nonbonded, vlist, r_l, a_array);
		}

		// Now call force functions //
		lj_energy = compute_urlj(coords, nonbonded_array1, nonbonded_array2, sigma_array, a_array, k_array, rb_array, total_force, num_neighs, vlist, &z);
		bonded_energy = compute_bonded(coords, bonded_array1, bonded_array2, bonded_avg_array, num_bonded, total_force);
		angle_energy = compute_angle(coords, angle_array1, angle_array2, angle_array3, angle_avg_array, num_angles, total_force);
		omega_energy = compute_omega(coords, CA1_omega_array, C_omega_array, N_omega_array, CA2_omega_array, omega_array, num_omega, total_force);
		planar_energy = compute_omega(coords, A_planar_array, B_planar_array, D_planar_array, C_planar_array, planar_array, num_planar, total_force);
		improper_energy = compute_omega(coords, B_improper_array, A_improper_array, C_improper_array, D_improper_array, improper_array, num_improper, total_force);

		fcheck = 0.0;
		for (j = 0; j < num_atoms; j++){
			fcheck += sqrt(total_force[j][0]*total_force[j][0] + total_force[j][1]*total_force[j][1] + total_force[j][2]*total_force[j][2]);
		}
    	fcheck = fcheck/num_atoms;

    	compute_damping(velocs, total_force, prev_total_force, num_atoms, dt, damping);

        // Update Velocities 
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 3; k++){
				velocs[j][k] += (0.5 * dt * (total_force[j][k]+prev_total_force[j][k]));
			}
		}

		memcpy(prev_total_force, total_force, sizeof(prev_total_force));

		
        if (fireit % 10000 == 0){

			printf("Fireit = %d \n", fireit);
            printf("\t fcheck = %.8e \n", fcheck);

            snprintf(buf, sizeof(buf), "traj_em.xyz");
			FILE *out_file = fopen(buf, "a+"); // write only 
			// test for files not existing. 
			if (out_file == NULL){   
			    printf("Error! Could not open file\n"); 
			    exit(-1); // must include stdlib.h 
			}
			fprintf(out_file, "%d \n Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" Origin=\"%lf %lf %lf\"  \n", num_atoms, box_size, box_size, box_size, -box_size/2, -box_size/2, -box_size/2);
			avg_x = 0;
			avg_y = 0;
			avg_z = 0;
			for (j = 0; j < num_atoms; j++){
				avg_x += coords[j][0];// - box_size*nearbyint(coords[j][0]/box_size);
				avg_y += coords[j][1];// - box_size*nearbyint(coords[j][1]/box_size);
				avg_z += coords[j][2];// - box_size*nearbyint(coords[j][2]/box_size);
			}
			avg_x = avg_x / num_atoms;
			avg_y = avg_y / num_atoms;
			avg_z = avg_z / num_atoms;

			for (j = 0; j < num_atoms; j++){
				//loop to print variable length atom names
				for (k = 0; k < atom_len_array[j]; k++){
					fprintf(out_file, "%c", atom_array[j][k]);
				}
				x_box = coords[j][0] - avg_x; 
				y_box = coords[j][1] - avg_y; 
				z_box = coords[j][2] - avg_z; 
				fprintf(out_file, "\t%lf\t%lf\t%lf\n", x_box, y_box, z_box); // write to file
			}
			fclose(out_file);
			
        }
        
	
	}
	

	// Langevin NVT //

	int split = 100;
	int tot_steps = 1e7;
	int split_count = 0;

    fcheck = 10*Ftol;
    fireit = 0;
    count = -1;
    temp = 1;
    double gamma = 0.001;
    double min_x,max_x,min_y,max_y,min_z,max_z;
    double exp_gamma_dt =  exp(-gamma * dt);
	double fluc_coeff = sqrt(1. - exp(-2. * gamma * dt)) * sqrt(target_temp);
    while (count < 1e7){
		count += 1;
		// Update Half Velocities //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 3; k++){
				half_velocs[j][k] = velocs[j][k] + (0.5 * dt_lang * total_force[j][k]);
			}
		}

		// Update Half Coordinates //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 3; k++){
				half_coords[j][k] = coords[j][k] + 0.5 * dt_lang * half_velocs[j][k];
			}
		}

		// Calc v_prime //
		generate_gaussian(half_num_atoms, gaussian_array);
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 3; k++){
				v_prime[j][k] = exp_gamma_dt * half_velocs[j][k] + fluc_coeff * gaussian_array[j][k];
			}
		}

		// Update Coordinates //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 3; k++){
				coords[j][k] = half_coords[j][k] + 0.5 * dt_lang * v_prime[j][k];
			}
		}

		// Compute Force //
        memset(total_force, 0, sizeof(total_force));

        // Checking Disp for update on Verlet List 
		disp_threshold_count = 0;
		for (j=0; j<num_atoms; j++){
			disp_x = displacement_array[j][0]-coords[j][0];
			disp_y = displacement_array[j][1]-coords[j][1];
			disp_z = displacement_array[j][2]-coords[j][2];
			disp = sqrt(disp_x*disp_x + disp_y*disp_y + disp_z*disp_z);
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
			num_neighs = generate_vlist(coords, nonbonded_array1, nonbonded_array2, sigma_array, num_nonbonded, vlist, r_l, a_array);
		}

		// Now call force functions //
		lj_energy = compute_urlj(coords, nonbonded_array1, nonbonded_array2, sigma_array, a_array, k_array, rb_array, total_force, num_neighs, vlist, &z);
		bonded_energy = compute_bonded(coords, bonded_array1, bonded_array2, bonded_avg_array, num_bonded, total_force);
		angle_energy = compute_angle(coords, angle_array1, angle_array2, angle_array3, angle_avg_array, num_angles, total_force);
		omega_energy = compute_omega(coords, CA1_omega_array, C_omega_array, N_omega_array, CA2_omega_array, omega_array, num_omega, total_force);
		planar_energy = compute_omega(coords, A_planar_array, B_planar_array, D_planar_array, C_planar_array, planar_array, num_planar, total_force);
		improper_energy = compute_omega(coords, B_improper_array, A_improper_array, C_improper_array, D_improper_array, improper_array, num_improper, total_force);

		// Update Velocities //
		for (j = 0; j < num_atoms; j++){
			for (k = 0; k < 3; k++){
				velocs[j][k] = v_prime[j][k] + 0.5 * dt_lang * total_force[j][k];
			}
		} 

        if (count % 10000 == 0){
        	
        	printf("Count = %d\n",count);

            snprintf(buf, sizeof(buf), "traj_data/traj_%s_alpha_%s_beta_%s_T-6.xyz",argv[1],argv[2],argv[3]);
			FILE *out_file = fopen(buf, "a+"); // write only 
			// test for files not existing. 
			if (out_file == NULL){   
			    printf("Error! Could not open file\n"); 
			    exit(-1); // must include stdlib.h 
			}
			fprintf(out_file, "%d \n Lattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" Origin=\"%lf %lf %lf\"  \n", num_atoms, box_size, box_size, box_size, -box_size/2, -box_size/2, -box_size/2);

			for (j = 0; j < num_atoms; j++){
				//loop to print variable length atom names
				for (k = 0; k < atom_len_array[j]; k++){
					fprintf(out_file, "%c", atom_array[j][k]);
				}
				x_box = coords[j][0]; 
				y_box = coords[j][1]; 
				z_box = coords[j][2]; 
				fprintf(out_file, "\t%lf\t%lf\t%lf\n", x_box, y_box, z_box); // write to file
			}
			fclose(out_file);
			
			

			
        }		



	}



}
