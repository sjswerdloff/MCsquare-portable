/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_geometry.h"

void verif_position(Hadron *hadron, DATA_CT *ct){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64);  

  __assume_aligned(&hadron->v_type, 64);

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(	hadron->v_x[i] < 0 || hadron->v_y[i] < 0 || hadron->v_z[i] < 0 || 
    	hadron->v_x[i] >= ct->Length[0] || hadron->v_y[i] >= ct->Length[1] || hadron->v_z[i] >= ct->Length[2]) hadron->v_type[i] = Unknown;
  }

  return;
}


void get_CT_Offset(Hadron *hadron, DATA_CT *ct, int *v_index){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64);

  __assume_aligned(&hadron->v_type, 64);  

  __assume_aligned(v_index, 64);

  // Calcul de l'offset : Offset = x + ct->Nx * y + ct->Nx * ct->Ny * z
  // L'axe x du repère simulation ne correspond pas à l'axe x du repère CT : x_simu = -x_ct + Lx
  // Conversion position -> index CT : index = floor(x/dx)

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_index[i] = 	(int)floor( (-hadron->v_x[i] + ct->Length[0]) / ct->VoxelLength[0] ) 
			+ ct->GridSize[0] * (int)floor( hadron->v_y[i] / ct->VoxelLength[1] ) 
			+ ct->GridSize[0] * ct->GridSize[1] * (int)floor( hadron->v_z[i] / ct->VoxelLength[2] );
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(	hadron->v_x[i] < 0 || hadron->v_y[i] < 0 || hadron->v_z[i] < 0 || 
    	hadron->v_x[i] >= ct->Length[0] || hadron->v_y[i] >= ct->Length[1] || hadron->v_z[i] >= ct->Length[2]) hadron->v_type[i] = Unknown;
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_index[i] < 0 || v_index[i] >= ct->Nbr_voxels) hadron->v_type[i] = Unknown;
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_index[i] = 0;
  }

  return;
}


void Dist_To_Material_Interface(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE dist_max, int *v_init_index, VAR_COMPUTE *v_init_density, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64); 
  __assume_aligned(&hadron->v_type, 64);  

  __assume_aligned(v_init_index, 64);
  __assume_aligned(v_init_density, 64);
  __assume_aligned(v_result, 64);


  ALIGNED_(64) VAR_COMPUTE v_mass_distance[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_mass_distance[i] = 0;
  }

  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);

  ALIGNED_(64) int v_index[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_index[i] = v_init_index[i];
  }

  ALIGNED_(64) int v_index2[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_index2[i] = v_init_index[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_run[i] = 1.0;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_run[i] = 0.0;
  }

    int run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_run[i];
  }

  while(run != 0){

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if((v_mass_distance[i] / v_init_density[i]) < dist_max && ct->material[v_index[i]] == ct->material[v_index2[i]]){

          Dist_To_Interface(&tmp, ct, v_step);
          Update_position(&tmp, v_step);

          v_index[i] = v_index2[i];
          get_CT_Offset(&tmp, ct, v_index2);

          v_mass_distance[i] += v_step[i] * ct->density[v_index[i]];
        }
        else{
          v_run[i] = 0.0;
        }
    }

    run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_run[i];
  }
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_result[i] = v_mass_distance[i] / v_init_density[i];
  }

  return;
}


void Dist_To_Interface(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_result){
  
  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64); 

  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64); 

  __assume_aligned(v_result, 64); 


  ALIGNED_(64) VAR_COMPUTE v_DistX[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_DistX[i] = fabs(((floor(hadron->v_x[i]/ct->VoxelLength[0]) + (hadron->v_u[i] > 0)) * ct->VoxelLength[0] - hadron->v_x[i])/hadron->v_u[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_DistY[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_DistY[i] = fabs(((floor(hadron->v_y[i]/ct->VoxelLength[1]) + (hadron->v_v[i] > 0)) * ct->VoxelLength[1] - hadron->v_y[i])/hadron->v_v[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_DistZ[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_DistZ[i] = fabs(((floor(hadron->v_z[i]/ct->VoxelLength[2]) + (hadron->v_w[i] > 0)) * ct->VoxelLength[2] - hadron->v_z[i])/hadron->v_w[i]);
  }

  // Add safety increment to compensate for rounding errors and to be sure to pass the interface (2e-4 for float, 1.5e-8 for double);
  #if VAR_COMPUTE_PRECISION==1
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        v_result[i] = fmin(v_DistX[i], fmin(v_DistY[i], v_DistZ[i]));
    } 
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_result[i] < 1e-3) v_result[i] += 2e-4;
        else v_result[i] += 5e-5;
    }
  #else
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        v_result[i] = fmin(v_DistX[i], fmin(v_DistY[i], v_DistZ[i])) + 1.5e-8;
    } 
  #endif
  


  return;
}


void Update_position(Hadron *hadron, VAR_COMPUTE *v_step){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64); 

  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64); 

  __assume_aligned(v_step, 64); 

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      hadron->v_x[i] += v_step[i] * hadron->v_u[i];
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      hadron->v_y[i] += v_step[i] * hadron->v_v[i];
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      hadron->v_z[i] += v_step[i] * hadron->v_w[i];
  }

  return;
}


void CT_Transport(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density){

  __assume_aligned(v_s, 64);
  __assume_aligned(v_tau, 64);
  __assume_aligned(v_init_index, 64);
  __assume_aligned(v_init_density, 64);
  __assume_aligned(v_hinge_index, 64);
  __assume_aligned(&hadron->v_type, 64); 

  ALIGNED_(64) VAR_COMPUTE v_MassDistance[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_MassDistance[i] = v_s[i] * v_init_density[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  Dist_To_Interface(hadron, ct, v_step);

  ALIGNED_(64) int v_index[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_index[i] = v_init_index[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_density[i] = ct->density[v_index[i]];
  }

  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_run[i] = 1.0;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_run[i] = 0.0;
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_hinge_index[i] = -1;
  }

  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_HingeDistance[i] = v_MassDistance[i] - v_tau[i] * v_init_density[i];
  }

    int run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_run[i];
  }

  while(run != 0){

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_MassDistance[i] > v_step[i] * v_density[i]){ 
          v_MassDistance[i] -= v_step[i] * v_density[i];
        }
        else{
          v_run[i] = 0.0;
          v_step[i] = 0.0;
        }
    }

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_MassDistance[i] < v_HingeDistance[i] && v_hinge_index[i] == -1){
          v_hinge_index[i] = v_index[i];
        }
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_run[i] = 0.0;
    }
    get_CT_Offset(hadron, ct, v_index);
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        v_density[i] = ct->density[v_index[i]];
    }
    Dist_To_Interface(hadron, ct, v_step);

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_index[i] >= ct->Nbr_voxels || v_index[i] < 0){
          v_run[i] = 0.0;
          v_hinge_index[i] = 0;
        }
    }

    run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_run[i];
  }
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_step[i] = v_MassDistance[i] / v_density[i];
  }
  Update_position(hadron, v_step);

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_hinge_index[i] == -1){
        v_hinge_index[i] = v_index[i];
      }
  }



  return;
}


void CT_Transport_SPR(Hadron *hadron, DATA_CT *ct, Materials *material, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density){

  __assume_aligned(v_s, 64);
  __assume_aligned(v_tau, 64);
  __assume_aligned(v_init_index, 64);
  __assume_aligned(v_init_density, 64);
  __assume_aligned(v_hinge_index, 64);
  __assume_aligned(&hadron->v_type, 64); 


  ALIGNED_(64) int v_index[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_index[i] = v_init_index[i];
  }

  ALIGNED_(64) int v_material_label[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_material_label[i] = ct->material[v_index[i]];
  }

  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_density[i] = v_init_density[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  Dist_To_Interface(hadron, ct, v_step);

  ALIGNED_(64) int v_data_index[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_data_index[i] = (int)ceil( hadron->v_T[i] / (UMeV*PSTAR_BIN * hadron->v_mass[i]));
  }

  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  int i;
  for(i=0; i<VLENGTH; i++){
    v_stop_pow[i] = (VAR_COMPUTE)material[v_material_label[i]].Stop_Pow[v_data_index[i]];
  }
/*
  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);
*/


  ALIGNED_(64) VAR_COMPUTE v_MassDistance[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_MassDistance[i] = v_s[i] * v_init_density[i] * v_stop_pow[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_HingeDistance[i] = v_MassDistance[i] - v_tau[i] * v_init_density[i] * v_stop_pow[i];
  }


  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_run[i] = 1.0;
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_run[i] = 0.0;
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_hinge_index[i] = -1;
  }

    int run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_run[i];
  }

  while(run != 0){

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_MassDistance[i] > v_step[i] * v_density[i] * v_stop_pow[i]){ 
          v_MassDistance[i] -= v_step[i] * v_density[i] * v_stop_pow[i];
        }
        else{
          v_run[i] = 0.0;
          v_step[i] = 0.0;
        }
    }

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_MassDistance[i] < v_HingeDistance[i] && v_hinge_index[i] == -1){
          v_hinge_index[i] = v_index[i];
        }
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_run[i] = 0.0;
    }
    get_CT_Offset(hadron, ct, v_index);
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        v_density[i] = ct->density[v_index[i]];
    }
    Dist_To_Interface(hadron, ct, v_step);

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_material_label[i] != ct->material[v_index[i]]){
    	v_material_label[i] = ct->material[v_index[i]];
    //	Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);
      	for(i=0; i<VLENGTH; i++){
      	  v_stop_pow[i] = (VAR_COMPUTE)material[v_material_label[i]].Stop_Pow[v_data_index[i]];
      	}
        }
    }

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_index[i] >= ct->Nbr_voxels || v_index[i] < 0){
          v_run[i] = 0.0;
          v_hinge_index[i] = 0;
        }
    }

    run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_run[i];
  }
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_step[i] = v_MassDistance[i] / (v_density[i] * v_stop_pow[i]);
  }
  Update_position(hadron, v_step);

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_hinge_index[i] == -1){
        v_hinge_index[i] = v_index[i];
      }
  }



  return;
}


void CT_Transport_Random_Hinge(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density, VAR_COMPUTE *v_mask){

  __assume_aligned(v_s, 64);
  __assume_aligned(v_tau, 64);
  __assume_aligned(v_init_index, 64);
  __assume_aligned(v_init_density, 64);
  __assume_aligned(v_hinge_index, 64);
  __assume_aligned(v_mask, 64);
  __assume_aligned(&hadron->v_type, 64); 

  ALIGNED_(64) VAR_COMPUTE v_MassDistance[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_MassDistance[i] = v_s[i] * v_init_density[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  Dist_To_Interface(hadron, ct, v_step);

  ALIGNED_(64) int v_index[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_index[i] = v_init_index[i];
  }

  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_density[i] = ct->density[v_index[i]];
  }

  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_run[i] = v_mask[i];
  }
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_run[i] = 0.0;
  }

  ALIGNED_(64) unsigned short int v_material[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_material[i] = ct->material[v_index[i]];
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_hinge_index[i] = -1;
  }

  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_HingeDistance[i] = v_MassDistance[i] - v_tau[i] * v_init_density[i];
  }

    int run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_run[i];
  }

  while(run != 0){

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_MassDistance[i] > v_step[i] * v_density[i] && v_run[i] != 0.0){
          v_MassDistance[i] -= v_step[i] * v_density[i];
        }
        else{
          v_run[i] = 0.0;
          v_step[i] = 0.0;
        }
    }

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_MassDistance[i] < v_HingeDistance[i] && v_hinge_index[i] == -1){
          v_hinge_index[i] = v_index[i];
        }
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_type[i] == Unknown) v_run[i] = 0.0;
    }
    get_CT_Offset(hadron, ct, v_index);
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        v_density[i] = ct->density[v_index[i]];
    }
    Dist_To_Interface(hadron, ct, v_step);

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
    if(v_index[i] >= ct->Nbr_voxels || v_index[i] < 0) v_run[i] = 0.0;

        else if(v_material[i] != ct->material[v_index[i]]){
          v_run[i] = 0.0;
          v_mask[i] = 0.0;
        }
    }

    run = 0;
  #pragma omp simd reduction(+:run)
  for (int i = 0; i < VLENGTH; i++) {
      run += v_run[i];
  }
  }  

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_mask[i] == 0.0){
        v_MassDistance[i] = 0.0;
      }
  }

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_step[i] = v_MassDistance[i] / v_density[i];
  }
  Update_position(hadron, v_step);

    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(v_hinge_index[i] == -1 && v_mask[i] != 0.0){
        v_hinge_index[i] = v_index[i];
      }
  }

  return;
}


void Update_direction(Hadron *hadron, VAR_COMPUTE *v_theta, VAR_COMPUTE *v_phi){

  __assume_aligned(v_theta, 64); 
  __assume_aligned(v_phi, 64); 

  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64); 

  ALIGNED_(64) VAR_COMPUTE v_cosT[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_cosT[i] = cos(v_theta[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_sinT[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_sinT[i] = sin(v_theta[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_cosP[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_cosP[i] = cos(v_phi[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_sinP[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_sinP[i] = sin(v_phi[i]);
  }

  ALIGNED_(64) VAR_COMPUTE v_prev_u[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_prev_u[i] = hadron->v_u[i];
  }


//  if(fabs(1.0 - hadron->v_w[vALL]) > 1e-10){	// Si direction non parallèle à l'axe z
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
    if(hadron->v_w[i] < 0.999999 && hadron->v_w[i] > -0.999999){	// Si direction non parallèle à l'axe z

        hadron->v_u[i] = (v_prev_u[i]*v_cosT[i] 
    			+ (v_sinT[i] / (sqrt(1.0-hadron->v_w[i]*hadron->v_w[i]))) * (v_prev_u[i]*hadron->v_w[i]*v_cosP[i] - hadron->v_v[i]*v_sinP[i]));

        hadron->v_v[i] = (hadron->v_v[i]*v_cosT[i] 
    			+ (v_sinT[i] / (sqrt(1.0-hadron->v_w[i]*hadron->v_w[i]))) * (hadron->v_v[i]*hadron->v_w[i]*v_cosP[i] + v_prev_u[i]*v_sinP[i]));

        hadron->v_w[i] = (hadron->v_w[i]*v_cosT[i] - sqrt(1.0-hadron->v_w[i]*hadron->v_w[i])*v_sinT[i]*v_cosP[i]);

      }
      else{
        hadron->v_v[i] = v_sinT[i] * v_sinP[i];

        if(hadron->v_w[i] > 0){		// Si direction parallère à l'axe z
          hadron->v_u[i] = v_sinT[i] * v_cosP[i];
          hadron->v_w[i] = v_cosT[i];
        }
        else{				// Si direction antiparallère à l'axe z
          hadron->v_u[i] = -v_sinT[i] * v_cosP[i];
          hadron->v_w[i] = -v_cosT[i];
        }
      }
  }



  // Si la norme dévie trop de 1, on renormalise
  ALIGNED_(64) VAR_COMPUTE v_norme[VLENGTH];
    #pragma omp simd
  for (int i = 0; i < VLENGTH; i++) {
      v_norme[i] = sqrt(hadron->v_u[i]*hadron->v_u[i] + hadron->v_v[i]*hadron->v_v[i] + hadron->v_w[i]*hadron->v_w[i]);
  }

//  if(fabs(v_norme[vALL]-1) > 1e-10){	

        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        hadron->v_u[i] = hadron->v_u[i] / v_norme[i];
    }
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        hadron->v_v[i] = hadron->v_v[i] / v_norme[i];
    }
        #pragma omp simd
    for (int i = 0; i < VLENGTH; i++) {
        hadron->v_w[i] = hadron->v_w[i] / v_norme[i];
    }
//  }  


  return;
}


void Update_buffer_direction(Hadron_buffer *secondary_hadron, VAR_COMPUTE theta, VAR_COMPUTE phi){
  VAR_COMPUTE cosT = cos(theta);
  VAR_COMPUTE sinT = sin(theta);
  VAR_COMPUTE cosP = cos(phi);
  VAR_COMPUTE sinP = sin(phi);

  if(fabs(secondary_hadron->w) < 0.999999){	// Si direction non parallèle à l'axe z

    VAR_COMPUTE Prev_u = secondary_hadron->u;	// valeur initiale de u
					// nécessaire car particule->u est mis à jour
    secondary_hadron->u = Prev_u*cosT + (sinT/(sqrt(1-secondary_hadron->w*secondary_hadron->w)))*(Prev_u*secondary_hadron->w*cosP - secondary_hadron->v*sinP);
    secondary_hadron->v = secondary_hadron->v*cosT + (sinT/(sqrt(1-secondary_hadron->w*secondary_hadron->w)))*(secondary_hadron->v*secondary_hadron->w*cosP + Prev_u*sinP);
    secondary_hadron->w = secondary_hadron->w*cosT - sqrt(1-secondary_hadron->w*secondary_hadron->w)*sinT*cosP;
  }
  else{
    secondary_hadron->v = sinT*sinP;

    if(secondary_hadron->w > 0){		// Si direction parallère à l'axe z
      secondary_hadron->u = sinT*cosP;
      secondary_hadron->w = cosT;
    }
    else{				// Si direction antiparallère à l'axe z
      secondary_hadron->u = -sinT*cosP;
      secondary_hadron->w = -cosT;
    }
  }

  // Si la norme dévie trop de 1, on renormalise
  VAR_COMPUTE norme = sqrt(secondary_hadron->u*secondary_hadron->u + secondary_hadron->v*secondary_hadron->v + secondary_hadron->w*secondary_hadron->w);
//  if(fabs(norme-1) > 1e-14){	

    secondary_hadron->u = secondary_hadron->u / norme;
    secondary_hadron->v = secondary_hadron->v / norme;
    secondary_hadron->w = secondary_hadron->w / norme;
//  }
}
