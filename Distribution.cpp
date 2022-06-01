#include "pch.h"
#include "Distribution.h"
using namespace std;

Hypergeom_distr::Hypergeom_distr(Hypergeom_distr& distr)
{
	a = distr.a;
	b = distr.b;
	k = distr.k;
	last_calc = distr.last_calc;
	if (distr.probs)
	{
		probs = new double[k+1];
		memcpy_s(probs, k+1, distr.probs, distr.k+1);
	}
	else
		probs = nullptr;
}

Hypergeom_distr::~Hypergeom_distr()
{
	if(probs)
	 delete[] probs; 
}

void Hypergeom_distr::set_param(int _a, int _b, int _k)
{
	if (_a < 1 || _b<1 || _k<1 || _k>_a || _k>_b)
		return;

	a = _a;
	b = _b;
	k = _k;
	last_calc = -1;
	if (probs)
		delete[] probs;
	probs = new double[k + 1];
}

void Hypergeom_distr::calc_all_probs()
{
	calc_probs_to(k);

	return;
}

void Hypergeom_distr::calc_next_prob()
{
	if (a < 1 || b < 1 || k < 1 || k > a || k > b||last_calc>k)
		return;

	if (last_calc > -1)
		probs[last_calc + 1] = probs[last_calc] * (a - last_calc) * (k - last_calc) / (last_calc + 1) / (b - k + last_calc + 1);
	else
	{
		int n = a + b, i = 0;
		probs[0] = 1;
		bool flg_1 = false, flg_2 = false;

		for (int j = 1; j < 2 * n - a - k + 1; ++j)
		{
			if (!flg_1)
				probs[0] *= j;
			else
				probs[0] *= j - b;

			if (!flg_2)
				probs[0] /= j;
			else
				probs[0] /= j - b + k;

			if (j == b)
				flg_1 = true;
			if (j == b - k)
				flg_2 = true;
		}
	}

	++last_calc;

	return;
}

double Hypergeom_distr::calc_probs_to(int n)
{
	if (a < 1 || b < 1 || k < 1 || k > a || k > b || last_calc>k)
		return -1;

	while (last_calc < n)
		calc_next_prob();

	return probs[n];
}

void Hypergeom_bern::gen_sample()
{
	if (sim_freq)
	{
		delete[] sim_freq;
		sim_freq = nullptr;
	}
	int a = h1.get_a();
	int b = h1.get_b();
	int k = h1.get_k();
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
	int a = h1.get_a();
	int b = h1.get_b();
	int k = h1.get_k();
	const double* probs = h1.get_probs();
	double l, x;
	int i;
	if (a < 1 || b < 1 || k < 1 || k > a || k > b)
		return;


	sim_freq = new int[k + 1];
	fill(sim_freq, sim_freq + k + 1, 0);

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> rndm(0, 1);
	for (int j = 0; j < sample_sz; ++j)
	{
		l = h1.calc_probs_to(0);
		i = 0;
		x = rndm(gen);
		while (x > l)
		{
			l += h1.calc_probs_to(i + 1);
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
	, power_n_samples_nmb(10'000)
{ 
	std::fill(p_lvls, p_lvls + 101, 0); 
}

Chi_sq::~Chi_sq()
{
	if(power_n_dep)
		delete[] power_n_dep;
}

/*!
* \param smpl Класс, генерирующий выборку с заданными H1 и методом
* \param sz Распределения соответствующее H0
*/
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

/*!
* \param _sample_sz Размер выборки, для которой получается значение p-value
* \param _samples_nmb Количество выборок, для которых считается значение p-value
* \param _alpha Заданный уровень значимости
*/
void Chi_sq::set_p_lvls(int _sample_sz, int _samples_nmb, double _alpha)
{
	sample_sz = _sample_sz;
	samples_nmb = _samples_nmb;
	alpha = _alpha;
	std::fill(p_lvls, p_lvls + 101, 0);
}

/*!
* \param _start_sz Минимальный размер выборки
* \param _steps_nmb Количество различных размеров выборок
* \param _step_sz Отличие размера следующей выборки от предыдущей
* \param _power_n_samples_nmb Количество значений p-value необходимых для получения мощности
* \param _alpha Заданный уровень значимости
*/
void Chi_sq::set_power_n(int _start_sz, int _steps_nmb, int _step_sz, int _power_n_samples_nmb, double _alpha)
{
	start_sz = _start_sz;
	steps_nmb = _steps_nmb;
	step_sz = _step_sz;
	power_n_samples_nmb = _power_n_samples_nmb;
	alpha = _alpha;
}

/*!
* \param vals Выборка из p-value
* \param sz Размер выборки
*/
double Chi_sq::calc_power(double* vals, int sz)
{
	int k = 0;

	for (int i = 0; i < sz; ++i)
		if (vals[i] < alpha)
			++k;

	return double(k) / sz;
}

/*!
* \param smpl Класс, генерирующий выборку с заданными H1 и методом
* \param sz Распределения соответствующее H0
*/
void Chi_sq::gen_p_levels(Hypergeom_sample* smpl, Hypergeom_distr& h0)
{
	Hypergeom_sample* generator= smpl->copy();
	generator->set_sz(sample_sz);

	double *p_vals=new double[samples_nmb];

	h0.calc_all_probs();
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

/*!
* \param smpl Класс, генерирующий выборку с заданными H1 и методом
* \param sz Распределения соответствующее H0
*/
void Chi_sq::power_n_dependance(Hypergeom_sample* smpl, Hypergeom_distr& h0)
{
	Hypergeom_sample* generator = smpl->copy();

	if (power_n_dep)
		delete[] power_n_dep;
	power_n_dep = new double[steps_nmb];

	double* p_vals = new double[power_n_samples_nmb];
	h0.calc_all_probs();

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

Hypergeom_sample::Hypergeom_sample(const Hypergeom_sample& smpl):h1(smpl.h1),sample_sz(smpl.sample_sz)
{
	if (smpl.sim_freq)
	{
		sim_freq = new int[sample_sz];
		memcpy_s(sim_freq, sample_sz, smpl.sim_freq, smpl.sample_sz);
	}
	else
		sim_freq = nullptr;
}

void Hypergeom_table::gen_sample()
{
	if (sim_freq)
	{
		delete[] sim_freq;
		sim_freq = nullptr;
	}
	int a = h1.get_a();
	int b = h1.get_b();
	int k = h1.get_k();
	const double* probs = h1.get_probs();
	double l, x;
	int i;
	if (a < 1 || b < 1 || k < 1 || k > a || k > b)
		return;

	h1.calc_all_probs();
	sim_freq = new int[k + 1];
	fill(sim_freq, sim_freq + k + 1, 0);

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> rndm(0, 1);
	for (int j = 0; j < sample_sz; ++j)
	{
		l = probs[0];
		i = 0;
		x = rndm(gen);
		while (x > l)
		{
			l += probs[i+1];
			++i;
		}
		++sim_freq[i];
	}

	return;
}
