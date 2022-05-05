#pragma once
#include<iostream>
#include<vector>
#include<algorithm>
#include<random>
#include"PROBDIST.H"

class Hypergeom_distr
{
private:
	int a = 10, b = 10, k = 10;
	int h_a = 10, h_b = 10, h_k = 10;
	int d_f = 0, n = 100;
	int samples_nmb = 100, sample_sz = 100000;
	int type = 1;
	double chi_sq = 0;
	double p_value = 0;
	double power = 0;
	double alpha = 0.05;
	int start_sz = 100, steps_nmb = 10, step_sz = 100, power_n_sample_sz = 1000;

	int* mod_distr = nullptr;
	double* th_distr = nullptr; //prob
	double* power_n = nullptr;
	double p_lvl_distr[101];
public:
	Hypergeom_distr(){}
	Hypergeom_distr(int a_, int b_, int k_) :a{ a_ }, b{ b_ }, k{ k_ }, h_a{ a_ }, h_b{ b_ }, h_k{ k_ }{}
	~Hypergeom_distr();

	void set_param(int a, int b, int k, int n);

	void set_hypothesis(int a, int b, int k, int sample_sz, int samples_nmb, double alpha = 0.05);

	void set_power_n_param(int start_sz, int steps_nmb, int steps_sz, int sample_sz);

	int get_k() { return k; }
	double get_chi_sq() { return chi_sq; }
	int get_d_f() { return d_f; }
	double get_p_val() { return p_value; }
	double get_power() { return power; }
	double get_alpha() { return alpha; }
	int get_init_pos() { return start_sz; }
	int get_steps() { return steps_nmb; }
	int get_step() { return step_sz; }
	const int* get_gen_distr() { return mod_distr; }
	const double* get_th_distr() { return th_distr; }
	const double* get_p_levels() { return p_lvl_distr; }
	const double* get_power_n() { return power_n; }



	void calc_chi_sq();
	void gen_distr(int n);
	virtual int* gen_distr(int a, int b, int k, int n)=0;
	void calc_distr();
	virtual void gen_p_levels();

	double* calc_distr(int a, int b, int k);
	void calc_p_value(const int* mod_distr, const double* th_distr, int k, int h_k, double& p_value, int n);
	void calc_p_value() { calc_p_value(mod_distr, th_distr, k, k, p_value, n); }
	void power_n_dependence();
};

class Hypergeom_distr_inv : public Hypergeom_distr
{
public:
	Hypergeom_distr_inv(int _a, int _b, int _k) :Hypergeom_distr(_a, _b, _k) {}
	~Hypergeom_distr_inv(){}
	virtual int* gen_distr(int a, int b, int k, int test_nmb);
};

class Hypergeom_distr_bern : public Hypergeom_distr
{
public:
	Hypergeom_distr_bern(int _a, int _b, int _k):Hypergeom_distr(_a,_b,_k){}
	~Hypergeom_distr_bern() {};
	virtual int* gen_distr(int a, int b, int k, int test_nmb);
};
