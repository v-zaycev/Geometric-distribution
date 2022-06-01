#pragma once
#include<iostream>
#include<vector>
#include<algorithm>
#include<random>
#include"PROBDIST.H"

/// �����, ���������� ���������� � ���������� �������������������� ������������� � ����������� ��������� ����������� ���������.
class Hypergeom_distr
{
private:
	int a, b, k;
	int last_calc;
	double* probs;
public:
	/// �����������, ��������� ��������� � ��������� ����������� �������������
	Hypergeom_distr(int _a = 10, int _b = 10, int _k = 10) :a(_a), b(_b), k(_k), last_calc(-1), probs(new double[_k + 1]) {};
	/// ����������� �����������
	Hypergeom_distr(Hypergeom_distr& distr);
	~Hypergeom_distr();
	/// �������, �������� ��� ��������� �������������
	void set_param(int _a, int _b, int _k);
	/// �������, ����������� ��� ����������� ��� ������� �������������
	void calc_all_probs();
	/// �������, ����������� ����������� ��� ������� ���������, ��� �������� ��� ��� �� ���������� 
	void calc_next_prob();
	/// �������, ����������� ��� ����������� ��� ��������� � 0 �� n-�� � ������������ ��������� ����������� �����������
	double calc_probs_to(int n);
	int get_a() { return a; }
	int get_b() { return b; }
	int get_k() { return k; }
	const double* get_probs() { return probs; };
};

/// ���������� �����, ���������� ���������� ������� ������� � �������������, ��� �������� ��� ����� ��������������.
class Hypergeom_sample
{
protected:
	int* sim_freq;
	Hypergeom_distr& h1;
	int sample_sz;
public:
	/// �����������, ��������� ������ ��� �������� ������������� � ������� �������
	Hypergeom_sample(Hypergeom_distr& _h1, int _n = 1000) :sample_sz(_n), sim_freq(nullptr), h1(_h1) {};
	/// ����������� �����������
	Hypergeom_sample(const Hypergeom_sample& smpl);
	virtual ~Hypergeom_sample() { delete[] sim_freq; };
	/// �������, ������������ �������
	virtual void gen_sample()=0;
	/// �������, ����������� ��������� ����� ���������� ��� ������������ ������
	virtual Hypergeom_sample* copy()=0;	
	void set_sz(int sz) { sample_sz = sz; };
	int get_a() { return h1.get_a(); }
	int get_b() { return h1.get_b(); }
	int get_k() { return h1.get_k(); }
	int get_n() { return sample_sz; }
	/// ������� ������������ �����, ��������������� ������ �������������, ������� ���������� ����� 
	virtual int get_type() = 0;
	const int* get_sim_freq() { return sim_freq; };
};

/// ��������� ������ Hypergeom_sample, ����������� ����� �������� �������.
class Hypergeom_inv : public Hypergeom_sample
{
public:
	Hypergeom_inv(const Hypergeom_inv& a) :Hypergeom_sample(a) {};
	Hypergeom_inv(Hypergeom_distr& h1, int n = 1000) :Hypergeom_sample(h1, n) {};
	/// �������, ������������ ������� ����������� ������� �������� �������
	virtual void gen_sample();
	virtual int get_type() { return 0; };
	virtual Hypergeom_sample* copy() { return new Hypergeom_inv(*this); };
	virtual ~Hypergeom_inv(){};
};

/// ��������� ������ Hypergeom_sample, ����������� ��������� ����� �������� �������.
class Hypergeom_table : public Hypergeom_sample
{
public:
	Hypergeom_table(const Hypergeom_table& a) :Hypergeom_sample(a) {};
	Hypergeom_table(Hypergeom_distr &h1, int n = 1000) :Hypergeom_sample(h1, n) {};
	/// �������, ������������ ������� ��������� ������� �������� �������
	virtual void gen_sample();
	virtual int get_type() { return 2; };
	virtual Hypergeom_sample* copy() { return new Hypergeom_table(*this); };
	virtual ~Hypergeom_table() {};
};

/// ��������� ������ Hypergeom_sample, ������������ ��������� ��������.
class Hypergeom_bern : public Hypergeom_sample
{
public:
	Hypergeom_bern(const Hypergeom_inv& a) :Hypergeom_sample(a) {};
	Hypergeom_bern(Hypergeom_distr& h1, int n = 1000) :Hypergeom_sample(h1, n) {};
	/// �������, ������������ ������� ����������� ������� ��������� ��������
	virtual void gen_sample();
	virtual int get_type() { return 1; };
	virtual Hypergeom_sample* copy() { return new Hypergeom_bern(*this); };
	virtual ~Hypergeom_bern() {};
};

/*!

	\brief
	�����, ����������� ������� ��������� �������� ��� �������� H0 � H1 � �������� ���������, ����������� ��� ������ ����������.

*/
class Chi_sq
{
private:
	double chi_sq; /// �������� ���������� ��-�������
	double p_val; /// �������� p-value 
	double power; /// �������� 
	double alpha;
	int d_f; /// ����� �������� �������
	int sample_sz; /// ������ �������, ��� ��������� ������� p-value
	int samples_nmb; /// ����� �������, ��� ��������� ������� p-value
	double p_lvls[101]; /// �������� ������� ������������� p-value �� ������� [0,1] � ����� 0.01
	double* power_n_dep; /// �������� �������� ��� �������� ������������������ �������� ������� 
	int start_sz, steps_nmb, step_sz, power_n_samples_nmb;

public:
	/// �����������, �������� ������� ���������
	Chi_sq();
	~Chi_sq();

	/// �������, �������� ��������� ����������� ��� ��������� ������� p-value.
	void set_p_lvls(int _sample_sz, int _samples_nmb, double _alpha);
	/// �������, �������� ���������, ����������� ��� ��������� ����������� �������� �� ������� �������.
	void set_power_n(int _start_sz, int _steps_nmb, int _step_sz, int _power_n_samples_nmb, double _alpha); 
	/// �������, ��������� �������� p-value ��� �������� ������� � H0.
	double calc_p_val(Hypergeom_sample *smpl, Hypergeom_distr& h0);
	/// �������, ��������� �������� �������� ��� �������� ������� �� p-value.
	double calc_power(double* vals, int sz);
	/// �������, ������������ ������� �� p-value ��� ��������� H0 � H1.
	void gen_p_levels(Hypergeom_sample* smpl, Hypergeom_distr& h0);
	/// �������, ���������� �������� �������� ��� �������� H0 � H1.
	void power_n_dependance(Hypergeom_sample* smpl, Hypergeom_distr& h0);

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
