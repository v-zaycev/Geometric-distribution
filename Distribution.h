#pragma once
#include<iostream>
#include<vector>
#include<algorithm>
#include<random>
#include"PROBDIST.H"

class Hypergeom_distr_0
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
	Hypergeom_distr_0(){}
	Hypergeom_distr_0(int a_, int b_, int k_) :a{ a_ }, b{ b_ }, k{ k_ }, h_a{ a_ }, h_b{ b_ }, h_k{ k_ }{}
	~Hypergeom_distr_0();

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

class Hypergeom_distr_inv : public Hypergeom_distr_0
{
public:
	Hypergeom_distr_inv(int _a, int _b, int _k) :Hypergeom_distr_0(_a, _b, _k) {}
	~Hypergeom_distr_inv(){}
	virtual int* gen_distr(int a, int b, int k, int test_nmb);
};

class Hypergeom_distr_bern : public Hypergeom_distr_0
{
public:
	Hypergeom_distr_bern(int _a, int _b, int _k):Hypergeom_distr_0(_a,_b,_k){}
	~Hypergeom_distr_bern() {};
	virtual int* gen_distr(int a, int b, int k, int test_nmb);
};


class Hypergeom_distr
{
private:
	int a, b, k;
	double* probs;
public:
	Hypergeom_distr() :a(10), b(10), k(10), probs(nullptr) {};
	Hypergeom_distr(int _a, int _b, int _k) :a(_a), b(_b), k(_k), probs(nullptr) {};
	~Hypergeom_distr() { delete[] probs; };
	void set_param(int _a, int _b, int _k);
	void calc_probs();
	int get_a() { return a; }
	int get_b() { return b; }
	int get_k() { return k; }
	const double* get_probs() { return probs; };
};

class Hypergeom_sample
{
protected:
	int* sim_freq;
	int sample_sz;
	int a, b, k;
public:
	Hypergeom_sample() :sample_sz(1000),a(10),b(10), k(10), sim_freq(nullptr) {};
	Hypergeom_sample(const Hypergeom_sample& smpl);
	Hypergeom_sample(int _a, int _b, int _k, int _n = 100) :sample_sz(_n), a(_a), b(_b), k(_k), sim_freq(nullptr) {};
	virtual ~Hypergeom_sample() { delete[] sim_freq; };
	virtual void gen_sample()=0;
	virtual Hypergeom_sample* copy()=0;
	void set_all(int _a, int _b, int _k, int _n);
	void set_param(int _a, int _b, int _k);
	void set_sz(int sz) { sample_sz = sz; };
	int get_n() { return sample_sz; }
	int get_a() { return a; }
	int get_b() { return b; }
	int get_k() { return k; }
	virtual int get_type() = 0;
	const int* get_sim_freq() { return sim_freq; };
};

class Hypergeom_inv : public Hypergeom_sample
{
public:
	Hypergeom_inv() :Hypergeom_sample() {};
	Hypergeom_inv(const Hypergeom_inv& a) :Hypergeom_sample(a) {};
	Hypergeom_inv(int a, int b, int k, int n = 100) :Hypergeom_sample(a, b, k, n) {};
	virtual void gen_sample();
	virtual int get_type() { return 0; };
	virtual Hypergeom_sample* copy() { return new Hypergeom_inv(*this); };
	virtual ~Hypergeom_inv(){};
};

class Hypergeom_bern : public Hypergeom_sample
{
public:
	Hypergeom_bern() :Hypergeom_sample() {};
	Hypergeom_bern(const Hypergeom_inv& a) :Hypergeom_sample(a) {};
	Hypergeom_bern(int a, int b, int k, int n = 100) :Hypergeom_sample(a, b, k, n) {};
	virtual void gen_sample();
	virtual int get_type() { return 1; };
	virtual Hypergeom_sample* copy() { return new Hypergeom_bern(*this); };
	virtual ~Hypergeom_bern() {};
};

class Chi_sq
{
private:
	double chi_sq;
	double p_val;
	double power;
	double alpha;
	int d_f;
	int sample_sz;
	int samples_nmb;
	double p_lvls[101];
	double* power_n_dep;
	int start_sz, steps_nmb, step_sz, power_n_samples_nmb;

public:
	Chi_sq();
	~Chi_sq();

	void set_p_lvls(int _sample_sz, int _samples_nmb, double _alpha);
	void set_power_n(int _start_sz, int _steps_nmb, int _step_sz, int _power_n_samples_nmb, double _alpha);

	double calc_p_val(Hypergeom_sample *smpl, Hypergeom_distr& h0);
	double calc_power(double* vals, int sz);
	void gen_p_levels(Hypergeom_sample* smpl, Hypergeom_distr& h0);
	void power_n_dependence(Hypergeom_sample* smpl, Hypergeom_distr& h0);

	double get_chi_sq() { return chi_sq; }
	int get_d_f() { return d_f; }
	double get_p_val() { return p_val; }
	double get_power() { return power; }
	double get_smpl_sz() { return sample_sz; }
	double get_smpls_nmb() { return samples_nmb; }
	double get_alpha() { return alpha; }
	int get_start_pos() { return start_sz; }
	int get_steps_nmb() { return steps_nmb; }
	int get_power_n_samples_nmb() { return power_n_samples_nmb; }
	int get_step_sz() { return step_sz; }
	const double* get_p_levels() { return p_lvls; }
	const double* get_power_n() { return power_n_dep; }

};



//class Chi_sq
//{
//private:
//	Hypergeom_sample* sample_generator;
//	Hypergeom_distr* h0;
//	double chi_sq;
//	double p_val;
//	double power;
//	double alpha;
//	int d_f;
//	int sample_sz;
//	int samples_nmb;
//	double p_lvls[101];
//	double* power_n_dep;
//	int start_sz, steps_nmb, step_sz, power_n_samples_nmb;
//
//public:
//	Chi_sq();
//	Chi_sq(Hypergeom_sample* smpl, Hypergeom_distr* _h0);
//	~Chi_sq();
//
//	void set_new_generator(Hypergeom_sample* smpl) {
//		delete sample_generator;
//		sample_generator = smpl;
//	}
//	void set_generator_param(int a, int b, int k) { sample_generator->set_param(a, b, k); }
//	void set_generator_all(int a, int b, int k, int sz)
//	{
//		sample_generator->set_param(a, b, k);
//		sample_generator->set_sz(sz);
//
//	}
//	void set_new_H0(Hypergeom_distr* _h0) {
//		delete h0;
//		h0 = _h0;
//	}
//	void set_H0_param(int a, int b, int k) { h0->set_param(a, b, k;) }
//	void set_p_lvls(int _sample_sz, int _samples_nmb, double _alpha);
//	void set_power_n(int _start_sz, int _steps_nmb, int _step_sz, int _power_n_samples_nmb, double _alpha);
//
//	double calc_p_val(Hypergeom_sample* smpl);
//	double calc_power(double* vals, int sz);
//	void gen_p_levels();
//	void power_n_dependence();
//
//	int get_k() { return sample_generator->get_k(); }
//	double get_chi_sq() { return chi_sq; }
//	int get_d_f() { return d_f; }
//	double get_p_val() { return p_val; }
//	double get_power() { return power; }
//	double get_alpha() { return alpha; }
//	int get_init_pos() { return start_sz; }
//	int get_steps() { return steps_nmb; }
//	int get_step() { return step_sz; }
//	const int* get_gen_freq() { return sample_generator->get_sim_freq(); }
//	const double* get_distr() { return h0->get_probs(); }
//	const double* get_p_levels() { return p_lvls; }
//	const double* get_power_n() { return power_n_dep; }
//
//};