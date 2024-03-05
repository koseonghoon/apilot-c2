#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8130422103780136687) {
   out_8130422103780136687[0] = delta_x[0] + nom_x[0];
   out_8130422103780136687[1] = delta_x[1] + nom_x[1];
   out_8130422103780136687[2] = delta_x[2] + nom_x[2];
   out_8130422103780136687[3] = delta_x[3] + nom_x[3];
   out_8130422103780136687[4] = delta_x[4] + nom_x[4];
   out_8130422103780136687[5] = delta_x[5] + nom_x[5];
   out_8130422103780136687[6] = delta_x[6] + nom_x[6];
   out_8130422103780136687[7] = delta_x[7] + nom_x[7];
   out_8130422103780136687[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6233895866051428142) {
   out_6233895866051428142[0] = -nom_x[0] + true_x[0];
   out_6233895866051428142[1] = -nom_x[1] + true_x[1];
   out_6233895866051428142[2] = -nom_x[2] + true_x[2];
   out_6233895866051428142[3] = -nom_x[3] + true_x[3];
   out_6233895866051428142[4] = -nom_x[4] + true_x[4];
   out_6233895866051428142[5] = -nom_x[5] + true_x[5];
   out_6233895866051428142[6] = -nom_x[6] + true_x[6];
   out_6233895866051428142[7] = -nom_x[7] + true_x[7];
   out_6233895866051428142[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_7070445437383799799) {
   out_7070445437383799799[0] = 1.0;
   out_7070445437383799799[1] = 0;
   out_7070445437383799799[2] = 0;
   out_7070445437383799799[3] = 0;
   out_7070445437383799799[4] = 0;
   out_7070445437383799799[5] = 0;
   out_7070445437383799799[6] = 0;
   out_7070445437383799799[7] = 0;
   out_7070445437383799799[8] = 0;
   out_7070445437383799799[9] = 0;
   out_7070445437383799799[10] = 1.0;
   out_7070445437383799799[11] = 0;
   out_7070445437383799799[12] = 0;
   out_7070445437383799799[13] = 0;
   out_7070445437383799799[14] = 0;
   out_7070445437383799799[15] = 0;
   out_7070445437383799799[16] = 0;
   out_7070445437383799799[17] = 0;
   out_7070445437383799799[18] = 0;
   out_7070445437383799799[19] = 0;
   out_7070445437383799799[20] = 1.0;
   out_7070445437383799799[21] = 0;
   out_7070445437383799799[22] = 0;
   out_7070445437383799799[23] = 0;
   out_7070445437383799799[24] = 0;
   out_7070445437383799799[25] = 0;
   out_7070445437383799799[26] = 0;
   out_7070445437383799799[27] = 0;
   out_7070445437383799799[28] = 0;
   out_7070445437383799799[29] = 0;
   out_7070445437383799799[30] = 1.0;
   out_7070445437383799799[31] = 0;
   out_7070445437383799799[32] = 0;
   out_7070445437383799799[33] = 0;
   out_7070445437383799799[34] = 0;
   out_7070445437383799799[35] = 0;
   out_7070445437383799799[36] = 0;
   out_7070445437383799799[37] = 0;
   out_7070445437383799799[38] = 0;
   out_7070445437383799799[39] = 0;
   out_7070445437383799799[40] = 1.0;
   out_7070445437383799799[41] = 0;
   out_7070445437383799799[42] = 0;
   out_7070445437383799799[43] = 0;
   out_7070445437383799799[44] = 0;
   out_7070445437383799799[45] = 0;
   out_7070445437383799799[46] = 0;
   out_7070445437383799799[47] = 0;
   out_7070445437383799799[48] = 0;
   out_7070445437383799799[49] = 0;
   out_7070445437383799799[50] = 1.0;
   out_7070445437383799799[51] = 0;
   out_7070445437383799799[52] = 0;
   out_7070445437383799799[53] = 0;
   out_7070445437383799799[54] = 0;
   out_7070445437383799799[55] = 0;
   out_7070445437383799799[56] = 0;
   out_7070445437383799799[57] = 0;
   out_7070445437383799799[58] = 0;
   out_7070445437383799799[59] = 0;
   out_7070445437383799799[60] = 1.0;
   out_7070445437383799799[61] = 0;
   out_7070445437383799799[62] = 0;
   out_7070445437383799799[63] = 0;
   out_7070445437383799799[64] = 0;
   out_7070445437383799799[65] = 0;
   out_7070445437383799799[66] = 0;
   out_7070445437383799799[67] = 0;
   out_7070445437383799799[68] = 0;
   out_7070445437383799799[69] = 0;
   out_7070445437383799799[70] = 1.0;
   out_7070445437383799799[71] = 0;
   out_7070445437383799799[72] = 0;
   out_7070445437383799799[73] = 0;
   out_7070445437383799799[74] = 0;
   out_7070445437383799799[75] = 0;
   out_7070445437383799799[76] = 0;
   out_7070445437383799799[77] = 0;
   out_7070445437383799799[78] = 0;
   out_7070445437383799799[79] = 0;
   out_7070445437383799799[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_4955099949312172907) {
   out_4955099949312172907[0] = state[0];
   out_4955099949312172907[1] = state[1];
   out_4955099949312172907[2] = state[2];
   out_4955099949312172907[3] = state[3];
   out_4955099949312172907[4] = state[4];
   out_4955099949312172907[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4955099949312172907[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4955099949312172907[7] = state[7];
   out_4955099949312172907[8] = state[8];
}
void F_fun(double *state, double dt, double *out_476684378344736844) {
   out_476684378344736844[0] = 1;
   out_476684378344736844[1] = 0;
   out_476684378344736844[2] = 0;
   out_476684378344736844[3] = 0;
   out_476684378344736844[4] = 0;
   out_476684378344736844[5] = 0;
   out_476684378344736844[6] = 0;
   out_476684378344736844[7] = 0;
   out_476684378344736844[8] = 0;
   out_476684378344736844[9] = 0;
   out_476684378344736844[10] = 1;
   out_476684378344736844[11] = 0;
   out_476684378344736844[12] = 0;
   out_476684378344736844[13] = 0;
   out_476684378344736844[14] = 0;
   out_476684378344736844[15] = 0;
   out_476684378344736844[16] = 0;
   out_476684378344736844[17] = 0;
   out_476684378344736844[18] = 0;
   out_476684378344736844[19] = 0;
   out_476684378344736844[20] = 1;
   out_476684378344736844[21] = 0;
   out_476684378344736844[22] = 0;
   out_476684378344736844[23] = 0;
   out_476684378344736844[24] = 0;
   out_476684378344736844[25] = 0;
   out_476684378344736844[26] = 0;
   out_476684378344736844[27] = 0;
   out_476684378344736844[28] = 0;
   out_476684378344736844[29] = 0;
   out_476684378344736844[30] = 1;
   out_476684378344736844[31] = 0;
   out_476684378344736844[32] = 0;
   out_476684378344736844[33] = 0;
   out_476684378344736844[34] = 0;
   out_476684378344736844[35] = 0;
   out_476684378344736844[36] = 0;
   out_476684378344736844[37] = 0;
   out_476684378344736844[38] = 0;
   out_476684378344736844[39] = 0;
   out_476684378344736844[40] = 1;
   out_476684378344736844[41] = 0;
   out_476684378344736844[42] = 0;
   out_476684378344736844[43] = 0;
   out_476684378344736844[44] = 0;
   out_476684378344736844[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_476684378344736844[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_476684378344736844[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_476684378344736844[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_476684378344736844[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_476684378344736844[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_476684378344736844[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_476684378344736844[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_476684378344736844[53] = -9.8000000000000007*dt;
   out_476684378344736844[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_476684378344736844[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_476684378344736844[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_476684378344736844[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_476684378344736844[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_476684378344736844[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_476684378344736844[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_476684378344736844[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_476684378344736844[62] = 0;
   out_476684378344736844[63] = 0;
   out_476684378344736844[64] = 0;
   out_476684378344736844[65] = 0;
   out_476684378344736844[66] = 0;
   out_476684378344736844[67] = 0;
   out_476684378344736844[68] = 0;
   out_476684378344736844[69] = 0;
   out_476684378344736844[70] = 1;
   out_476684378344736844[71] = 0;
   out_476684378344736844[72] = 0;
   out_476684378344736844[73] = 0;
   out_476684378344736844[74] = 0;
   out_476684378344736844[75] = 0;
   out_476684378344736844[76] = 0;
   out_476684378344736844[77] = 0;
   out_476684378344736844[78] = 0;
   out_476684378344736844[79] = 0;
   out_476684378344736844[80] = 1;
}
void h_25(double *state, double *unused, double *out_486402922224284601) {
   out_486402922224284601[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4515941117287146063) {
   out_4515941117287146063[0] = 0;
   out_4515941117287146063[1] = 0;
   out_4515941117287146063[2] = 0;
   out_4515941117287146063[3] = 0;
   out_4515941117287146063[4] = 0;
   out_4515941117287146063[5] = 0;
   out_4515941117287146063[6] = 1;
   out_4515941117287146063[7] = 0;
   out_4515941117287146063[8] = 0;
}
void h_24(double *state, double *unused, double *out_7781227329859495183) {
   out_7781227329859495183[0] = state[4];
   out_7781227329859495183[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6693155540894296036) {
   out_6693155540894296036[0] = 0;
   out_6693155540894296036[1] = 0;
   out_6693155540894296036[2] = 0;
   out_6693155540894296036[3] = 0;
   out_6693155540894296036[4] = 1;
   out_6693155540894296036[5] = 0;
   out_6693155540894296036[6] = 0;
   out_6693155540894296036[7] = 0;
   out_6693155540894296036[8] = 0;
   out_6693155540894296036[9] = 0;
   out_6693155540894296036[10] = 0;
   out_6693155540894296036[11] = 0;
   out_6693155540894296036[12] = 0;
   out_6693155540894296036[13] = 0;
   out_6693155540894296036[14] = 1;
   out_6693155540894296036[15] = 0;
   out_6693155540894296036[16] = 0;
   out_6693155540894296036[17] = 0;
}
void h_30(double *state, double *unused, double *out_5798718511680226446) {
   out_5798718511680226446[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7034274075794394690) {
   out_7034274075794394690[0] = 0;
   out_7034274075794394690[1] = 0;
   out_7034274075794394690[2] = 0;
   out_7034274075794394690[3] = 0;
   out_7034274075794394690[4] = 1;
   out_7034274075794394690[5] = 0;
   out_7034274075794394690[6] = 0;
   out_7034274075794394690[7] = 0;
   out_7034274075794394690[8] = 0;
}
void h_26(double *state, double *unused, double *out_4318790407105879486) {
   out_4318790407105879486[0] = state[7];
}
void H_26(double *state, double *unused, double *out_774437798413089839) {
   out_774437798413089839[0] = 0;
   out_774437798413089839[1] = 0;
   out_774437798413089839[2] = 0;
   out_774437798413089839[3] = 0;
   out_774437798413089839[4] = 0;
   out_774437798413089839[5] = 0;
   out_774437798413089839[6] = 0;
   out_774437798413089839[7] = 1;
   out_774437798413089839[8] = 0;
}
void h_27(double *state, double *unused, double *out_3902676211873354617) {
   out_3902676211873354617[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4859510763993969779) {
   out_4859510763993969779[0] = 0;
   out_4859510763993969779[1] = 0;
   out_4859510763993969779[2] = 0;
   out_4859510763993969779[3] = 1;
   out_4859510763993969779[4] = 0;
   out_4859510763993969779[5] = 0;
   out_4859510763993969779[6] = 0;
   out_4859510763993969779[7] = 0;
   out_4859510763993969779[8] = 0;
}
void h_29(double *state, double *unused, double *out_7646352294368479525) {
   out_7646352294368479525[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7544505420108786874) {
   out_7544505420108786874[0] = 0;
   out_7544505420108786874[1] = 1;
   out_7544505420108786874[2] = 0;
   out_7544505420108786874[3] = 0;
   out_7544505420108786874[4] = 0;
   out_7544505420108786874[5] = 0;
   out_7544505420108786874[6] = 0;
   out_7544505420108786874[7] = 0;
   out_7544505420108786874[8] = 0;
}
void h_28(double *state, double *unused, double *out_8145026344422347853) {
   out_8145026344422347853[0] = state[0];
}
void H_28(double *state, double *unused, double *out_2462106403039256300) {
   out_2462106403039256300[0] = 1;
   out_2462106403039256300[1] = 0;
   out_2462106403039256300[2] = 0;
   out_2462106403039256300[3] = 0;
   out_2462106403039256300[4] = 0;
   out_2462106403039256300[5] = 0;
   out_2462106403039256300[6] = 0;
   out_2462106403039256300[7] = 0;
   out_2462106403039256300[8] = 0;
}
void h_31(double *state, double *unused, double *out_1886074921141698207) {
   out_1886074921141698207[0] = state[8];
}
void H_31(double *state, double *unused, double *out_148229696179738363) {
   out_148229696179738363[0] = 0;
   out_148229696179738363[1] = 0;
   out_148229696179738363[2] = 0;
   out_148229696179738363[3] = 0;
   out_148229696179738363[4] = 0;
   out_148229696179738363[5] = 0;
   out_148229696179738363[6] = 0;
   out_148229696179738363[7] = 0;
   out_148229696179738363[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_8130422103780136687) {
  err_fun(nom_x, delta_x, out_8130422103780136687);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6233895866051428142) {
  inv_err_fun(nom_x, true_x, out_6233895866051428142);
}
void car_H_mod_fun(double *state, double *out_7070445437383799799) {
  H_mod_fun(state, out_7070445437383799799);
}
void car_f_fun(double *state, double dt, double *out_4955099949312172907) {
  f_fun(state,  dt, out_4955099949312172907);
}
void car_F_fun(double *state, double dt, double *out_476684378344736844) {
  F_fun(state,  dt, out_476684378344736844);
}
void car_h_25(double *state, double *unused, double *out_486402922224284601) {
  h_25(state, unused, out_486402922224284601);
}
void car_H_25(double *state, double *unused, double *out_4515941117287146063) {
  H_25(state, unused, out_4515941117287146063);
}
void car_h_24(double *state, double *unused, double *out_7781227329859495183) {
  h_24(state, unused, out_7781227329859495183);
}
void car_H_24(double *state, double *unused, double *out_6693155540894296036) {
  H_24(state, unused, out_6693155540894296036);
}
void car_h_30(double *state, double *unused, double *out_5798718511680226446) {
  h_30(state, unused, out_5798718511680226446);
}
void car_H_30(double *state, double *unused, double *out_7034274075794394690) {
  H_30(state, unused, out_7034274075794394690);
}
void car_h_26(double *state, double *unused, double *out_4318790407105879486) {
  h_26(state, unused, out_4318790407105879486);
}
void car_H_26(double *state, double *unused, double *out_774437798413089839) {
  H_26(state, unused, out_774437798413089839);
}
void car_h_27(double *state, double *unused, double *out_3902676211873354617) {
  h_27(state, unused, out_3902676211873354617);
}
void car_H_27(double *state, double *unused, double *out_4859510763993969779) {
  H_27(state, unused, out_4859510763993969779);
}
void car_h_29(double *state, double *unused, double *out_7646352294368479525) {
  h_29(state, unused, out_7646352294368479525);
}
void car_H_29(double *state, double *unused, double *out_7544505420108786874) {
  H_29(state, unused, out_7544505420108786874);
}
void car_h_28(double *state, double *unused, double *out_8145026344422347853) {
  h_28(state, unused, out_8145026344422347853);
}
void car_H_28(double *state, double *unused, double *out_2462106403039256300) {
  H_28(state, unused, out_2462106403039256300);
}
void car_h_31(double *state, double *unused, double *out_1886074921141698207) {
  h_31(state, unused, out_1886074921141698207);
}
void car_H_31(double *state, double *unused, double *out_148229696179738363) {
  H_31(state, unused, out_148229696179738363);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
