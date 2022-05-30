#include "pch.h"
#include "Distribution.h"
using namespace std;

double combination(int n, int k)
{
	if (k<0 || k>n || n < 1)
		return 1;
	else
	{
		int i_0 = max(k, n - k);
		double p = 1;
		for (int j = i_0 + 1; j < n + 1; ++j)
			p *= double(j) / (j - i_0);
		return p;
	}
}

void Hypergeom_distr_0::set_param(int a, int b, int k, int n)
{
	this->a = a;
	this->b = b;
	this->k = k;
	this->n = n;
	chi_sq = 0;
	p_value = 0;
	d_f = 0;
	delete[] mod_distr;
	mod_distr = nullptr;
	delete[] th_distr;
	th_distr = nullptr;

	return;
}

void Hypergeom_distr_0::set_hypothesis(int a, int b, int k, int sample_sz, int samples_nmb, double alpha)
{
	h_a = a;
	h_b = b;
	h_k = k;
	this->samples_nmb = samples_nmb;
	this->sample_sz = sample_sz;
	this->alpha = alpha;
	fill(p_lvl_distr, p_lvl_distr + 101, 0);

	return;
}

void Hypergeom_distr_0::set_power_n_param(int start_sz, int steps_nmb, int steps_sz, int sample_sz)
{
	this->step_sz = steps_sz;
	this->steps_nmb = steps_nmb;
	this->start_sz = start_sz;
	this->power_n_sample_sz = sample_sz;
	return;
}

Hypergeom_distr_0::~Hypergeom_distr_0()
{
	if (mod_distr!=nullptr)
		delete[] mod_distr;
	if (th_distr != nullptr)
		delete[] th_distr;
	if (power_n != nullptr)
		delete[] power_n;
}

void Hypergeom_distr_0::calc_distr()
{
	if (th_distr != nullptr)
		delete[] th_distr;

	th_distr = new double[k + 1];
	for (int i = 0; i < k + 1; ++i)
		th_distr[i] = combination(a, i) * combination(b, k - i) / combination(a + b, k);

	return;

}

double* Hypergeom_distr_0::calc_distr(int a, int b, int k)
{
	double* th_distr = new double[k + 1];
	for (int i = 0; i < k + 1; ++i)
		th_distr[i] = combination(a, i) * combination(b, k - i) / combination(a + b, k);

	return th_distr;
}

void Hypergeom_distr_0::calc_chi_sq()
{
	if (mod_distr == nullptr || th_distr == nullptr)
		return;

	double th = 0, mod = 0, sum = 0;

	for (int i = 0; i < k + 1; ++i)
	{
		sum += th_distr[i];
		th += th_distr[i];
		mod += mod_distr[i];
		if (n * th > 5 && n - n * sum > 5)
		{
			chi_sq += (n * th - mod) * (n * th - mod) / th / n;
			++d_f;
			th = 0;
			mod = 0;
		}
	}
	if (th > 0)
		chi_sq += (n * th - mod) * (n * th - mod) / th / n;
	else
		--d_f;

	p_value = 1 - pChi(chi_sq, d_f);

	return;
}

void Hypergeom_distr_0::calc_p_value(const int* mod_distr, const double* th_distr, int k, int h_k, double& p_value, int n)
{
	double th = 0, mod = 0, sum = 0;
	chi_sq = 0;
	d_f = 0;
	int r_k = count_if(mod_distr, mod_distr + k, [](int x) {return x > 0; });

	if (mod_distr == nullptr || th_distr == nullptr || r_k > h_k)
		return;

	for (int i = 0; i < max(r_k, h_k) + 1; ++i)
	{
		sum += th_distr[i];
		th += th_distr[i];
		if (i <= k)
			mod += mod_distr[i];
		if (n * th > 5 && n - n * sum > 5)
		{
			chi_sq += (n * th - mod) * (n * th - mod) / th / n;
			++d_f;
			th = 0;
			mod = 0;
		}
	}
	if (th > 0)
		chi_sq += (n * th - mod) * (n * th - mod) / th / n;

	p_value = 1 - pChi(chi_sq, d_f);

	return;
}

void Hypergeom_distr_0::gen_p_levels()
{
	p_lvl_distr[0] = 0;

	double* p_val = new double[sample_sz];

	int* mod_distr = nullptr;
	double* th_distr = nullptr;
	th_distr = calc_distr(h_a, h_b, h_k);

	for (int i = 0; i < sample_sz; ++i)
	{
		mod_distr = gen_distr(a, b, k, samples_nmb);


		calc_p_value(mod_distr, th_distr, k, h_k, p_val[i], samples_nmb);
		delete[] mod_distr;
	}
	delete[] th_distr;

	std::sort(p_val, p_val + sample_sz);

	power = 0;
	if (h_a == a && h_b == b && h_k == k)
	{
		for (int i = 0; i < sample_sz; ++i)
			if (p_val[i] < alpha)
				++power;
		power /= sample_sz;
	}
	else
	{
		double* th0_distr = calc_distr(a, b, k);
		double* th1_distr = calc_distr(h_a, h_b, h_k);


		double th = 0, mod = 0, sum = 0;
		int df_0 = 0, d_f = 0;

		for (int i = 0; i < k + 1; ++i)
		{
			sum += th0_distr[i];
			th += th0_distr[i];
			if (samples_nmb * th > 5 && samples_nmb - samples_nmb * sum > 5)
			{
				++df_0;
				th = 0;
				mod = 0;
			}
		}
		if (th == 0)
			--df_0;

		th = 0;
		sum = 0;
		for (int j = 0; j < h_k + 1; ++j)
		{
			sum += th1_distr[j];
			th += th1_distr[j];
			if (samples_nmb * th > 5 && samples_nmb - samples_nmb * sum > 5)
			{
				++d_f;
				th = 0;
			}
		}
		if (th == 0)
			--d_f;



		delete[] th0_distr;
		delete[] th1_distr;

		for (int i = 0; i < sample_sz; ++i)
			if (p_val[i] > pChi(xChi(alpha, df_0), d_f))
				++power;
		power /= sample_sz;
		power = 1 - power;
	}

	int j = 0, i = 1;
	while (j < sample_sz - 1)
		if (p_val[j + 1] < double(i) / 100)
			j++;
		else
		{
			p_lvl_distr[i] = double(j) / sample_sz;
			++i;
		}
	delete[] p_val;
	for (; i < 101; ++i)
		p_lvl_distr[i] = 1;

	return;
}

void Hypergeom_distr_0::power_n_dependence()
{
	if (power_n != nullptr)
	{
		delete[] power_n;
		power_n = nullptr;
	}

	power_n = new double[steps_nmb + 1];
	double* th0_distr = calc_distr(a, b, k);
	double* th1_distr = calc_distr(h_a, h_b, h_k);
	int* mod_distr = nullptr;
	double* p_val = nullptr;
	int current_sz = 0;
	for (int i = 0; i < steps_nmb + 1; ++i)
	{
		current_sz = start_sz + step_sz * i;
		p_val = new double[power_n_sample_sz];
		for (int j = 0; j < power_n_sample_sz; ++j)
		{
			mod_distr = gen_distr(a, b, k, current_sz);
			calc_p_value(mod_distr, th1_distr, k, h_k, p_val[j], current_sz);
			delete[] mod_distr;
		}
		std::sort(p_val, p_val + power_n_sample_sz);


		power_n[i] = 0;
		if (h_a == a && h_b == b && h_k == k)
		{
			for (int j = 0; j < power_n_sample_sz; ++j)
				if (p_val[j] < alpha)
					++power_n[i];
			power_n[i] /= power_n_sample_sz;
		}
		else
		{
			double th = 0, sum = 0;
			int df_0 = 0, d_f = 0;

			for (int j = 0; j < k + 1; ++j)
			{
				sum += th0_distr[j];
				th += th0_distr[j];
				if (current_sz * th > 5 && current_sz - current_sz * sum > 5)
				{
					++df_0;
					th = 0;
				}
			}
			if (th == 0)
				--df_0;
			th = 0;
			sum = 0;
			for (int j = 0; j < h_k + 1; ++j)
			{
				sum += th1_distr[j];
				th += th1_distr[j];
				if (current_sz * th > 5 && current_sz - current_sz * sum > 5)
				{
					++d_f;
					th = 0;
				}
			}
			if (th == 0)
				--d_f;

			for (int j = 0; j < power_n_sample_sz; ++j)
				if (p_val[j] > pChi(xChi(alpha, df_0), d_f))
					++power_n[i];
			power_n[i] /= power_n_sample_sz;
			power_n[i] = 1 - power_n[i];
		}
		delete[] p_val;
	}


	delete[] th0_distr;
	delete[] th1_distr;

	return;
}

void Hypergeom_distr_0::gen_distr(int n)
{
	this->n = n;
	if (mod_distr)
	{
		delete[] mod_distr;
		mod_distr = nullptr;
	}
	this->mod_distr=gen_distr(a, b, k, n);
	return;
}

void Hypergeom_distr::set_param(int _a, int _b, int _k)
{
	a = _a;
	b = _b;
	k = _k;
	if (probs)
		delete[] probs;
	probs = nullptr;
}

void Hypergeom_distr::calc_probs()
{
	if (probs != nullptr)
	{
		delete[] probs;
		probs = nullptr;
	}

	if (a < 1 || b < 1 || k < 1 || k > a || k > b)
		return;

	probs = new double[k + 1];
	for (int i = 0; i < k + 1; ++i)
		probs[i] = combination(a, i) * combination(b, k - i) / combination(a + b, k);

	return;

}

void Hypergeom_bern::gen_sample()
{
	if (sim_freq)
	{
		delete[] sim_freq;
		sim_freq = nullptr;
	}

	if (a < 1 || b < 1 || k < 1 || k > a || k > b)
		return;

	sim_freq = new int[k + 1];
	fill(sim_freq, sim_freq + k + 1, 0);
	int n = a + b, j = 0;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> rndm(0, 1);
	for (int i = 0; i < sample_sz; ++i)
	{
		j = 0;
		for (int i = 0; i < k; ++i)
			if (rndm(gen) < (a - j) / double(n - i))
				++j;
		++sim_freq[j];
	}
	return;
}


void Hypergeom_inv::gen_sample()
{
	if (sim_freq)
	{
		delete[] sim_freq;
		sim_freq = nullptr;
	}

	if (a < 1 || b < 1 || k < 1 || k > a || k > b)
		return;


	sim_freq = new int[k + 1];
	fill(sim_freq, sim_freq + k + 1, 0);
	int n = a + b, i = 0;
	double p_0 = 1, p, l, x;
	bool flg_1 = false, flg_2 = false;

	for (int j = 1; j < 2 * n - a - k + 1; ++j)
	{
		if (!flg_1)
			p_0 *= j;
		else
			p_0 *= j - b;

		if (!flg_2)
			p_0 /= j;
		else
			p_0 /= j - b + k;

		if (j == b)
			flg_1 = true;
		if (j == b - k)
			flg_2 = true;
	}

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> rndm(0, 1);
	for (int j = 0; j < sample_sz; ++j)
	{
		l = p = p_0;
		i = 0;
		x = rndm(gen);
		while (x > l)
		{
			p *= double(a - i) * (k - i) / (i + 1) / (b - k + i + 1);
			l += p;
			++i;
		}
		++sim_freq[i];
	}

	return;
}

Chi_sq::Chi_sq()
	: power_n_dep(nullptr)
	, chi_sq(0)
	, p_val(0)
	, power(0)
	, alpha(0.05)
	, d_f(0)
	, sample_sz(100)
	, samples_nmb(10'000)
	, start_sz(100) 
	, steps_nmb(10)
	, step_sz(100)
	, power_n_samples_nmb(10000)
{ 
	std::fill(p_lvls, p_lvls + 101, 0); 
}

//Chi_sq::Chi_sq(Hypergeom_sample* smpl, Hypergeom_distr* _h0)
//	:sample_generator(smpl)
//	, h0(_h0)
//	, chi_sq(0)
//	, p_val(0)
//	, power(0)
//	, alpha(0.05)
//	, d_f(0)
//	, sample_sz(100)
//	, samples_nmb(10'000)
//	, start_sz(100)
//	, steps_nmb(10)
//	, step_sz(100)
//	, power_n_samples_nmb(10000)
//{ 
//	std::fill(p_lvls, p_lvls + 101, 0); 
//}

Chi_sq::~Chi_sq()
{
	delete power_n_dep;
}

double Chi_sq::calc_p_val(Hypergeom_sample* smpl, Hypergeom_distr &h0)
{
	if (!smpl)
		return -1;

	const int* mod_distr = smpl->get_sim_freq();
	const double* th_distr = h0.get_probs();

	int n = smpl->get_n();
	int k = smpl->get_k();
	int h_k = h0.get_k();


	double th = 0, mod = 0, sum = 0;
	chi_sq = 0;
	d_f = 0;

	if (mod_distr == nullptr || th_distr == nullptr)
		return -1;

	for (int i = 0; i < max(k, h_k); ++i)
	{
		if (i <= k)
			mod += mod_distr[i];
		if (i <= h_k)
		{
			sum += th_distr[i];
			th += th_distr[i];
		}
		if (n * th > 5 && n - n * sum > 5)
		{
			chi_sq += (n * th - mod) * (n * th - mod) / th / n;
			++d_f;
			th = 0;
			mod = 0;
		}
	}
	if (th > 0)
		chi_sq += (n * th - mod) * (n * th - mod) / th / n;

	p_val = 1 - pChi(chi_sq, d_f);

	return p_val;
}

void Chi_sq::set_p_lvls(int _sample_sz, int _samples_nmb, double _alpha)
{
	sample_sz = _sample_sz;
	samples_nmb = _samples_nmb;
	alpha = _alpha;
	std::fill(p_lvls, p_lvls + 101, 0);
}

void Chi_sq::set_power_n(int _start_sz, int _steps_nmb, int _step_sz, int _power_n_samples_nmb, double _alpha)
{
	start_sz = _start_sz;
	steps_nmb = _steps_nmb;
	step_sz = _step_sz;
	power_n_samples_nmb = _power_n_samples_nmb;
	alpha = _alpha;
}

double Chi_sq::calc_power(double* vals, int sz)
{
	int k = 0;

	for (int i = 0; i < sz; ++i)
		if (vals[i] < alpha)
			++k;

	return double(k) / sz;
}

void Chi_sq::gen_p_levels(Hypergeom_sample* smpl, Hypergeom_distr& h0)
{
	Hypergeom_sample* generator= smpl->copy();
	generator->set_sz(sample_sz);

	double *p_vals=new double[samples_nmb];

	h0.calc_probs();
	for (int i = 0; i < samples_nmb; ++i)
	{
		generator->gen_sample();
		p_vals[i] = calc_p_val(generator, h0);
	}
		

	std::sort(p_vals, p_vals + samples_nmb);
	int j = 0, i = 1;
	while (j < samples_nmb - 1)
		if (p_vals[j + 1] < double(i) / 100)
			++j;
		else
		{
			p_lvls[i] = double(j) / samples_nmb;
			++i;
		}
	for (; i < 101; ++i)
		p_lvls[i] = 1;

	power = calc_power(p_vals, samples_nmb);

	delete[] p_vals;
	delete generator;

	return;
}

void Chi_sq::power_n_dependence(Hypergeom_sample* smpl, Hypergeom_distr& h0)
{

	Hypergeom_sample* generator = smpl->copy();

	power_n_dep = new double[steps_nmb];

	double* p_vals = new double[power_n_samples_nmb];
	h0.calc_probs();

	for (int k = 0; k < steps_nmb; ++k)
	{
		generator->set_sz(start_sz+k*step_sz);
		for (int i = 0; i < power_n_samples_nmb; ++i)
		{
			generator->gen_sample();
			p_vals[i] = calc_p_val(generator, h0);
		}
		power_n_dep[k] = calc_power(p_vals, power_n_samples_nmb);
	}
	
	delete generator;
	delete[] p_vals;

	return;
}

Hypergeom_sample::Hypergeom_sample(const Hypergeom_sample& smpl)
{
	a = smpl.a;
	b = smpl.b;
	k = smpl.k;
	sample_sz = smpl.sample_sz;

	if (smpl.sim_freq)
	{
		sim_freq = new int[sample_sz];
		memcpy_s(sim_freq, sample_sz, smpl.sim_freq, smpl.sample_sz);
	}
}

void Hypergeom_sample::set_all(int _a, int _b, int _k, int _n)
{
	a = _a;
	b = _b;
	k = _k;
	sample_sz = _n;
	delete[] sim_freq;
	sim_freq = nullptr;
}

void Hypergeom_sample::set_param(int _a, int _b, int _k)
{
	a = _a;
	b = _b;
	k = _k;
	delete[] sim_freq;
	sim_freq = nullptr;
}
