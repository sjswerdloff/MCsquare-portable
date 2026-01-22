/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_semi_infinite_slab.h"

void SemiInfiniteSlab_step(Hadron *hadron, Materials *material, Hadron_buffer *hadron_list, ControlPoint_parameters **layer_data, field_parameters **field_data, int *Hadron_ID, int *Nbr_hadrons, VAR_COMPUTE *RS_exit_position, DATA_config *config, machine_parameters *machine, VSLStreamStatePtr RNG_Stream){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64);

  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);

  __assume_aligned(&hadron->v_type, 64);

  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);


  Update_Hadron(hadron);


  int i,j,r;

  // Compute physical quantities
  ALIGNED_(64) int v_material_label[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_init_density[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_N_el[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_X0[VLENGTH];
  for(i=0; i<VLENGTH; i++){
    if(hadron->v_type[i] == Unknown){
      v_material_label[i] = 0;
      v_init_density[i] = 1;
      v_N_el[i] = 1;
      v_X0[i] = 1;
    }
    else{
      r = field_data[Hadron_ID[i]]->RS_num;
      // Bounds check for RS_num (matches compute_range_shifter.c validation)
      if (r < 0 || r >= machine->RS_number) {
          r = 0;  // Clamp to valid index
      }
      v_material_label[i] = machine->RS_Material[r];
      v_init_density[i] = machine->RS_Density[r];
      v_N_el[i] = material[machine->RS_Material[r]].N_el * v_init_density[i];
      v_X0[i] = material[machine->RS_Material[r]].X0 / v_init_density[i];
    }
  }

  // Compute total cross section
  ALIGNED_(64) VAR_COMPUTE v_Dist_Interface[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_Dist_Interface[i] = hadron->v_z[i] - RS_exit_position[i] + 1e-4;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_Dist_Interface[i] < 0) v_Dist_Interface[i] = 0;
  }

  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_stop_pow[i] = v_init_density[i] * hadron->v_charge[i]*hadron->v_charge[i] * v_stop_pow[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_step_max[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_step_max[i] = fmin(fmin(v_Dist_Interface[i], config->D_Max), (config->Epsilon_Max * hadron->v_T[i] / v_stop_pow[i]));
  }

  ALIGNED_(64) VAR_COMPUTE v_dE_max[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_dE_max[i] = v_step_max[i] * v_stop_pow[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_section[VLENGTH];
  Total_Hard_Cross_Section(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, config, v_section); 
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_section[i] += 1e-10;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_section[i] *= 1.017;
  }


  // Compute step length
  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(RNG_Stream, v_rnd);

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_step[i] = -log(v_rnd[i])/v_section[i];
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_step[i] > v_step_max[i]) v_step[i] = v_step_max[i];
  }  // stop at step_max


  // Compute CSDA + MS
  ALIGNED_(64) VAR_COMPUTE v_mean_dE[VLENGTH];
  Compute_dE2(hadron, v_N_el, v_init_density, material, (config->Te_Min*UMeV), v_material_label, v_step, v_mean_dE);  

  ALIGNED_(64) VAR_COMPUTE v_straggling[VLENGTH];
  Compute_Energy_straggling(hadron, v_N_el, (config->Te_Min*UMeV), v_step, v_straggling);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_straggling[i] = sqrt(v_straggling[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_dE[VLENGTH];
  rand_normal(RNG_Stream, v_dE, v_mean_dE, v_straggling);

  ALIGNED_(64) VAR_COMPUTE v_MS[VLENGTH];
  Compute_MS_Fippel(hadron, v_step, v_X0, v_MS);

  ALIGNED_(64) VAR_COMPUTE v_theta[VLENGTH];
  rand_normal_zero(RNG_Stream, v_theta, v_MS);

  ALIGNED_(64) VAR_COMPUTE v_phi[VLENGTH];
  rand_uniform(RNG_Stream, v_phi);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_phi[i] = 2*M_PI*v_phi[i];
  }


  // Lose energy
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_dE[i] = 0;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      hadron->v_T[i] = hadron->v_T[i] - v_dE[i];
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] != Unknown && hadron->v_T[i] <= (config->Ecut_Pro * UMeV)){
        hadron->v_type[i] = Unknown;
        hadron_list[Hadron_ID[i]].type = Unknown;
      }
  }
  Update_Hadron(hadron);

  // Update position and direction
  Update_position(hadron, v_step);
  Update_direction(hadron, v_theta, v_phi);


  // Hard interaction

  ALIGNED_(64) int v_interaction_type[VLENGTH];

  get_interaction_type(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, v_section, RNG_Stream, config, v_interaction_type);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_step[i] == v_step_max[i]) v_interaction_type[i] = 0;
  } // force ficitious interaction if step >= step_max

  // Ionization
  ALIGNED_(64) VAR_COMPUTE v_dE_hard[VLENGTH];
  Compute_Ionization_Energy(hadron, (config->Te_Min*UMeV), RNG_Stream, v_dE_hard);
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_dE_hard[i] = 0;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_interaction_type[i] == 1) hadron->v_T[i] = hadron->v_T[i] - v_dE_hard[i];
  }

  // Nuclear interaction
  DATA_Scoring tmp;
  int previous_Nbr_hadrons;
  for(i=0; i<VLENGTH; i++){
    if(hadron->v_type[i] != Unknown && v_interaction_type[i] == 2){
      previous_Nbr_hadrons = *Nbr_hadrons;
      Compute_Nuclear_interaction(i, hadron, material, v_material_label[i], hadron_list, Nbr_hadrons, &tmp, RNG_Stream, config);
      if(hadron->v_type[i] == Unknown) hadron_list[Hadron_ID[i]].type = Unknown;
      for(j=previous_Nbr_hadrons; j<*Nbr_hadrons; j++){
         layer_data[j] = layer_data[Hadron_ID[i]];
         field_data[j] = field_data[Hadron_ID[i]];
      }
    }
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] != Unknown && hadron->v_T[i] <= (config->Ecut_Pro * UMeV)){
        hadron->v_type[i] = Unknown;
        hadron_list[Hadron_ID[i]].type = Unknown;
      }
  }

}
