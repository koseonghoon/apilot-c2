#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_8130422103780136687);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6233895866051428142);
void car_H_mod_fun(double *state, double *out_7070445437383799799);
void car_f_fun(double *state, double dt, double *out_4955099949312172907);
void car_F_fun(double *state, double dt, double *out_476684378344736844);
void car_h_25(double *state, double *unused, double *out_486402922224284601);
void car_H_25(double *state, double *unused, double *out_4515941117287146063);
void car_h_24(double *state, double *unused, double *out_7781227329859495183);
void car_H_24(double *state, double *unused, double *out_6693155540894296036);
void car_h_30(double *state, double *unused, double *out_5798718511680226446);
void car_H_30(double *state, double *unused, double *out_7034274075794394690);
void car_h_26(double *state, double *unused, double *out_4318790407105879486);
void car_H_26(double *state, double *unused, double *out_774437798413089839);
void car_h_27(double *state, double *unused, double *out_3902676211873354617);
void car_H_27(double *state, double *unused, double *out_4859510763993969779);
void car_h_29(double *state, double *unused, double *out_7646352294368479525);
void car_H_29(double *state, double *unused, double *out_7544505420108786874);
void car_h_28(double *state, double *unused, double *out_8145026344422347853);
void car_H_28(double *state, double *unused, double *out_2462106403039256300);
void car_h_31(double *state, double *unused, double *out_1886074921141698207);
void car_H_31(double *state, double *unused, double *out_148229696179738363);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}