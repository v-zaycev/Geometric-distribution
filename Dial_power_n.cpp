// Dial_power_n.cpp : implementation file
//

#include "pch.h"
#include "Geom_distr.h"
#include "afxdialogex.h"
#include "Dial_power_n.h"


// Dial_power_n dialog

IMPLEMENT_DYNAMIC(Dial_power_n, CDialog)

Dial_power_n::Dial_power_n(CWnd* pParent /*=nullptr*/)
	: CDialog(IDD_DIALOG3, pParent)
	, type(FALSE)
	, init_sz(0)
	, step_sz(0)
	, steps_nmb(0)
	, alpha(0)
	, a(0)
	, b(0)
	, k(0)
	, h_a(0)
	, h_b(0)
	, h_k(0)
	, sample_sz(0)
{
}

Dial_power_n::Dial_power_n(BOOL _type, int _a, int _b, int _k, int _h_a, int _h_b, int _h_k, int _init_sz, int _step_sz, int _steps_nmb, int _sample_sz, double _alpha, CWnd* pParent)
	: CDialog(IDD_DIALOG3, pParent)
	, type(_type)
	, init_sz(_init_sz)
	, step_sz(_step_sz)
	, steps_nmb(_steps_nmb)
	, alpha(_alpha)
	, a(_a)
	, b(_b)
	, k(_k)
	, h_a(_h_a)
	, h_b(_h_b)
	, h_k(_h_k)
	, sample_sz(_sample_sz)
{
}

Dial_power_n::~Dial_power_n()
{
}

void Dial_power_n::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Radio(pDX, IDC_RADIO1, type);
	DDX_Text(pDX, IDC_EDIT7, init_sz);
	DDX_Text(pDX, IDC_EDIT8, step_sz);
	DDX_Text(pDX, IDC_EDIT9, steps_nmb);
	DDX_Text(pDX, IDC_EDIT10, alpha);
	DDX_Text(pDX, IDC_EDIT1, a);
	DDX_Text(pDX, IDC_EDIT2, b);
	DDX_Text(pDX, IDC_EDIT3, k);
	DDX_Text(pDX, IDC_EDIT4, h_a);
	DDX_Text(pDX, IDC_EDIT5, h_b);
	DDX_Text(pDX, IDC_EDIT6, h_k);
	DDX_Text(pDX, IDC_EDIT11, sample_sz);
}


BEGIN_MESSAGE_MAP(Dial_power_n, CDialog)
	ON_BN_CLICKED(IDOK, &Dial_power_n::OnBnClickedOk)
END_MESSAGE_MAP()


// Dial_power_n message handlers


void Dial_power_n::OnBnClickedOk()
{
	UpdateData(TRUE);
	if (a <= 0 || b <= 0 || k <= 0 || k > a || k > b || h_a <= 0 || h_b <= 0 || h_k <= 0 || h_k > h_a || h_k > h_b)
		AfxMessageBox(L"Parameter values are incorrect. Parameters must be positive and k<=a and k<=b.");
	else if (init_sz < 50 || init_sz > 10'000)
		AfxMessageBox(L"Select the initial size from the range (50, 10.000).");
	else if (steps_nmb < 5 || steps_nmb > 100)
		AfxMessageBox(L"Select the number of steps from the range (10, 100).");
	else if (step_sz < 10 || step_sz > 1'000)
		AfxMessageBox(L"Select the step size from the range (10, 1.000).");
	else if (sample_sz < 100 || sample_sz > 100'000)
		AfxMessageBox(L"Select the sample size from the range (100, 100.000).");
	else if (alpha < 0.001 || alpha > 1)
		AfxMessageBox(L"select alpha between 0.001 and 1.");
	else
		CDialog::OnOK();
}
