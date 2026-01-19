/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_EM_interaction.h"

void Total_Stop_Pow(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_stop_pow){

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_mass, 64);

  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_stop_pow, 64);

  int i;

  ALIGNED_(64) VAR_COMPUTE v_scaled_T[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_scaled_T[i] = hadron->v_T[i]/hadron->v_mass[i];
  }


  ALIGNED_(64) int v_index[VLENGTH];
//  v_index[vALL] = (int)(v_scaled_T[vALL] / (UMeV*PSTAR_BIN));
//  v_index[vALL] = (int)floor(v_scaled_T[vALL] / (UMeV*PSTAR_BIN));
  ALIGNED_(64) VAR_COMPUTE v_scaled_T2[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_scaled_T2[i] = v_scaled_T[i] / (UMeV*PSTAR_BIN);
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_index[i] = (int)floor(v_scaled_T2[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_data_Energy1[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_data_Energy1[i] = v_index[i] * UMeV * PSTAR_BIN;
  }

  ALIGNED_(64) VAR_COMPUTE v_data_Energy2[VLENGTH];
//  v_data_Energy2[vALL] = v_data_Energy1[vALL] + UMeV * PSTAR_BIN;
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_data_Energy2[i] = (v_index[i]+1) * UMeV * PSTAR_BIN;
  }


  ALIGNED_(64) VAR_COMPUTE v_Stop_Pow1[VLENGTH];
  for(i=0; i<VLENGTH; i++){
    v_Stop_Pow1[i] = (VAR_COMPUTE)material[v_material_label[i]].Stop_Pow[v_index[i]];
  }

  ALIGNED_(64) VAR_COMPUTE v_Stop_Pow2[VLENGTH];
  for(i=0; i<VLENGTH; i++){
    v_Stop_Pow2[i] = (VAR_COMPUTE)material[v_material_label[i]].Stop_Pow[v_index[i]+1];
  }

  vec_Linear_Interpolation(v_scaled_T, v_data_Energy1, v_data_Energy2, v_Stop_Pow1, v_Stop_Pow2, v_stop_pow);

  return;
}


void Total_Hard_Cross_Section(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, VAR_COMPUTE Te_min, VAR_COMPUTE *v_dE_max, DATA_config *config, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_T, 64);

  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_dE_max, 64);
  __assume_aligned(v_result, 64);

  ALIGNED_(64) VAR_COMPUTE v_cross_section[VLENGTH];
  cross_section_ionization(hadron, v_N_el, Te_min, v_cross_section);

  ALIGNED_(64) VAR_COMPUTE v_tmp_result[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(config->Simulate_Nuclear_Interactions == 1){
        total_Nuclear_cross_section(hadron, material, v_material_label, v_density, v_tmp_result);
        v_cross_section[i] += v_tmp_result[i];
      }
  }


  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      tmp.v_T[i] = hadron->v_T[i] - v_dE_max[i];
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(tmp.v_T[i] <= 0) tmp.v_T[i] = hadron->v_T[i];
  }
  Update_Hadron(&tmp);
    
  ALIGNED_(64) VAR_COMPUTE v_cross_section2[VLENGTH];
  cross_section_ionization(&tmp, v_N_el, Te_min, v_cross_section2);

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(config->Simulate_Nuclear_Interactions == 1){
        total_Nuclear_cross_section(&tmp, material, v_material_label, v_density, v_tmp_result);
        v_cross_section2[i] += v_tmp_result[i];
      }
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_result[i] = fmax(v_cross_section[i], v_cross_section2[i]);
  }

  return;
}


void get_interaction_type(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, VAR_COMPUTE Te_min, VAR_COMPUTE *v_dE_max, VAR_COMPUTE *v_tot_section, VSLStreamStatePtr RNG_Stream, DATA_config *config, int *v_result){

  __assume_aligned(&hadron->v_T, 64);

  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_dE_max, 64);
  __assume_aligned(v_tot_section, 64);
  __assume_aligned(v_result, 64);


    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_result[i] = 0;
  }

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(RNG_Stream, v_rnd);

  ALIGNED_(64) VAR_COMPUTE v_ionization_section[VLENGTH];
  cross_section_ionization(hadron, v_N_el, Te_min, v_ionization_section);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_ionization_section[i] = v_ionization_section[i] / v_tot_section[i];
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_rnd[i] <= v_ionization_section[i]){
        v_result[i] = 1;
      }
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(config->Simulate_Nuclear_Interactions == 1){
        ALIGNED_(64) VAR_COMPUTE v_nuclear_section[VLENGTH];
        total_Nuclear_cross_section(hadron, material, v_material_label, v_density, v_nuclear_section);
        v_nuclear_section[i] = (v_nuclear_section[i] / v_tot_section[i]) + v_ionization_section[i];

        if(v_rnd[i] <= v_ionization_section[i]){
          v_result[i] = 1;
        }
        else if(v_rnd[i] <= v_nuclear_section[i]){
          v_result[i] = 2;
        }
      }

      else{
        if(v_rnd[i] <= v_ionization_section[i]){
          v_result[i] = 1;
        }
      }
  }

  return;
}

void cross_section_ionization(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE Te_min, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);

  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);

  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_result, 64);



  ALIGNED_(64) VAR_COMPUTE v_log_result[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_log_result[i] = hadron->v_Te_max[i]/Te_min;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_log_result[i] = log(v_log_result[i]);
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_result[i] =  	2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[i] * hadron->v_charge[i]*hadron->v_charge[i] 
			* (	((1.0/Te_min) - (1.0/hadron->v_Te_max[i])) 
				- (hadron->v_beta2[i]/hadron->v_Te_max[i]) * v_log_result[i] 
				+ (hadron->v_Te_max[i]-Te_min) / (2*hadron->v_E[i]*hadron->v_E[i])
			  )
			/ (hadron->v_beta2[i]);
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_Te_max[i] <= Te_min) v_result[i] = 0.0;
  }

  return;
}




void Compute_L(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, Materials *material, VAR_COMPUTE Te_min, int *v_material_label, VAR_COMPUTE *v_result){

  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);
  
  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);


  ALIGNED_(64) VAR_COMPUTE v_log_result[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_log_result[i] = hadron->v_Te_max[i]/Te_min;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_log_result[i] = log(v_log_result[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_M[VLENGTH];

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_M[i] =	(2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[i] * hadron->v_charge[i]*hadron->v_charge[i] / hadron->v_beta2[i]) 
		* ( 	v_log_result[i]
			- (hadron->v_Te_max[i] - Te_min) * hadron->v_beta2[i] / hadron->v_Te_max[i]
			+ (hadron->v_Te_max[i]*hadron->v_Te_max[i] - Te_min*Te_min) / (4*hadron->v_E[i]*hadron->v_E[i])
		);
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_Te_max[i] <= Te_min) v_M[i] = 0;
  }

  Total_Stop_Pow(hadron, material, v_material_label, v_result);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_result[i] = v_density[i] * hadron->v_charge[i]*hadron->v_charge[i] * v_result[i] - v_M[i];
  }	// Pouvoir d'arrêt restreint en eV / cm

  return;
}


void Compute_dE2(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, Materials *material, VAR_COMPUTE Te_min, int *v_material_label, VAR_COMPUTE *v_s, VAR_COMPUTE *v_result){

  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_s, 64);
  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);
  
  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);


  // Valeur précalculées

  ALIGNED_(64) VAR_COMPUTE v_L[VLENGTH];
  Compute_L(hadron, v_N_el, v_density, material, Te_min, v_material_label, v_L);

  ALIGNED_(64) VAR_COMPUTE v_dE1[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_dE1[i] = v_L[i] * v_s[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_tau1[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_tau1[i] = hadron->v_T[i] / MC2_PRO;
  }

  ALIGNED_(64) VAR_COMPUTE v_e1[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_e1[i] = v_dE1[i] / hadron->v_T[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_C[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_C[i] = v_L[i] * hadron->v_beta2[i];
  }


  // Calcul numérique de la dérivée de C(E)

  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      tmp.v_T[i] = hadron->v_T[i] * CONST_DERIV;
  }
  Update_Hadron(&tmp);

  ALIGNED_(64) VAR_COMPUTE v_L2[VLENGTH];
  Compute_L(&tmp, v_N_el, v_density, material, Te_min, v_material_label, v_L2);

  ALIGNED_(64) VAR_COMPUTE v_C2[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_C2[i] = v_L2[i] * tmp.v_beta2[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_deriv_C[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_deriv_C[i] = (v_C2[i] - v_C[i]) / (tmp.v_T[i] - hadron->v_T[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_b[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_b[i] = hadron->v_T[i] * v_deriv_C[i] / v_C[i];
  }


  // Calcul de dE2

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_result[i] = v_dE1[i] * (	1 
					+ (v_e1[i] / ((1+v_tau1[i]) * (2+v_tau1[i]))) 
					+ (	v_e1[i]*v_e1[i] 
						* (2+2*v_tau1[i]+v_tau1[i]*v_tau1[i]) 
						/ ((1+v_tau1[i])*(1+v_tau1[i])*(2+v_tau1[i])*(2+v_tau1[i]))
					  ) 
					- (v_b[i] * v_e1[i] * (0.5 + 2*v_e1[i]/(3*(1+v_tau1[i])*(2+v_tau1[i])) + (1-v_b[i]) * v_e1[i]/6)) 
			 	);
  }

  return;
}


void Compute_Energy_straggling(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE Te_min, VAR_COMPUTE *v_s, VAR_COMPUTE *v_result){

  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_s, 64);
  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);
  
  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_result[i] = 	2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[i] * hadron->v_charge[i]*hadron->v_charge[i] * v_s[i] 
			* fmin(Te_min, hadron->v_Te_max[i]) 
			* (1 - 0.5*hadron->v_beta2[i]) / hadron->v_beta2[i];
  }

  return;
}


void Compute_MS_Fippel(Hadron *hadron, VAR_COMPUTE *v_s, VAR_COMPUTE *v_X0, VAR_COMPUTE *v_result){

  __assume_aligned(v_s, 64);
  __assume_aligned(v_X0, 64);
  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);
  
  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_result[i] = (CONST_MS_Fippel*UMeV * hadron->v_charge[i] / (hadron->v_beta2[i]*hadron->v_gamma[i]*MC2_PRO)) * sqrt(v_s[i]/v_X0[i]);
  }

  return;
}


void Compute_Ionization_Energy(Hadron *hadron, VAR_COMPUTE Te_min, VSLStreamStatePtr RNG_Stream, VAR_COMPUTE *v_result){

  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_Te[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_g[VLENGTH];

  ALIGNED_(64) VAR_COMPUTE v_mask[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_mask[i] = 1.0;
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_Te_max[i] < Te_min){
        v_result[i] = 0.0;
        v_mask[i] = 0.0;
      }
  }

    int run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_mask[i];
  }

  while(run != 0.0){
    rand_uniform(RNG_Stream, v_rnd);
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        v_Te[i] = ( Te_min * hadron->v_Te_max[i]) / ((1-v_rnd[i]) * hadron->v_Te_max[i] + v_rnd[i] * Te_min);
    }
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        v_g[i] = 1.0 - hadron->v_beta2[i] * (v_Te[i]/hadron->v_Te_max[i]) + v_Te[i]*v_Te[i]/(2*hadron->v_E[i]*hadron->v_E[i]);
    }
    rand_uniform(RNG_Stream, v_rnd);

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_rnd[i] <= v_g[i] && v_mask[i] == 1.0){
          v_result[i] = v_Te[i];
          v_mask[i] = 0.0;
        }
    }

    run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_mask[i];
  }
  }

  return;
}



