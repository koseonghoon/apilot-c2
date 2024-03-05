#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_6934754482349461241);
void live_err_fun(double *nom_x, double *delta_x, double *out_7134518254598619458);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_3168509789052392768);
void live_H_mod_fun(double *state, double *out_8812197276460084856);
void live_f_fun(double *state, double dt, double *out_5323663110479347867);
void live_F_fun(double *state, double dt, double *out_7291716453244524992);
void live_h_4(double *state, double *unused, double *out_8149121993453717660);
void live_H_4(double *state, double *unused, double *out_2533638448163623402);
void live_h_9(double *state, double *unused, double *out_1352873043042257376);
void live_H_9(double *state, double *unused, double *out_4753580487100824068);
void live_h_10(double *state, double *unused, double *out_2559315150682307729);
void live_H_10(double *state, double *unused, double *out_1646686746551854404);
void live_h_12(double *state, double *unused, double *out_861564556841939270);
void live_H_12(double *state, double *unused, double *out_8914896825206356398);
void live_h_35(double *state, double *unused, double *out_1791863090781450798);
void live_H_35(double *state, double *unused, double *out_6169333792881342689);
void live_h_32(double *state, double *unused, double *out_5578284502534678464);
void live_H_32(double *state, double *unused, double *out_368029420673730131);
void live_h_13(double *state, double *unused, double *out_8648356105277536640);
void live_H_13(double *state, double *unused, double *out_7277191461913882294);
void live_h_14(double *state, double *unused, double *out_1352873043042257376);
void live_H_14(double *state, double *unused, double *out_4753580487100824068);
void live_h_33(double *state, double *unused, double *out_5283889615969827563);
void live_H_33(double *state, double *unused, double *out_3018776788242485085);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}