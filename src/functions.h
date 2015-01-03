#ifndef theoretical_functions_h
#define theoretical_functions_h

class inverse;

// theoretical functions ********************************************
double opening_angle(double tilt_angle, double gamma);
double opening_angle_over_min(double tilt_angle, double gamma);
double min_opening_angle(double gamma);
double tilt_angle_prob(double tilt_angle, int helicity);
double tilt_angle_cos_prob(double tilt_angle_cos, int helicity);
double tilt_angle_abs_cos_prob(double tilt_angle_cos, int helicity);
double opening_angle_prob(double angle, int helicity, const inverse& tilt_from_op);

double tilt_angle(double opening_angle, double gamma);
double tilt_angle(double angle, double gamma, double& derivative);
double opening_angle_prob(double angle, int helicity, double gamma);

// ******************************************************************
double cut_angle(double gamma, int num_tilt_op, int num_opt);
double cut_angle(double gamma);

#endif
