// Dial_p_val.cpp : implementation file
//

#include "pch.h"
#include "Geom_distr.h"
#include "afxdialogex.h"
#include "Dial_p_val.h"


// Dial_p_val dialog

IMPLEMENT_DYNAMIC(Dial_p_val, CDialog)

Dial_p_val::Dial_p_val(CWnd* pParent)
	: CDialog(IDD_DIALOG2, pParent)
	, type(FALSE)
	, samples_nmb(10'000)
	, sample_sz(100)
	, a(10)
	, b(10)
	, k(10)
	, h_a(10)
	, h_b(10)
	, h_k(10)
	, alpha(0.05)
{
}

Dial_p_val::Dial_p_val(BOOL _type, int _sample_sz_p_val, int _sample_sz, int _a, int _b, int _k, int _h_a, int _h_b, int _h_k, double _alpha, CWnd* pParent)
	: CDialog(IDD_DIALOG2, pParent)
	, type(_type)
	, samples_nmb(_sample_sz_p_val)
	, sample_sz(_sample_sz)
	, a(_a)
	, b(_b)
	, k(_k)
	, h_a(_h_a)
	, h_b(_h_b)
	, h_k(_h_k)
	, alpha(_alpha)
{
}

Dial_p_val::~Dial_p_val()
{
}

void Dial_p_val::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Radio(pDX, IDC_RADIO1, type);
	DDX_Text(pDX, IDC_EDIT7, samples_nmb);
	DDX_Text(pDX, IDC_EDIT8, sample_sz);
	DDX_Text(pDX, IDC_EDIT1, a);
	DDX_Text(pDX, IDC_EDIT2, b);
	DDX_Text(pDX, IDC_EDIT3, k);
	DDX_Text(pDX, IDC_EDIT4, h_a);
	DDX_Text(pDX, IDC_EDIT5, h_b);
	DDX_Text(pDX, IDC_EDIT6, h_k);
	DDX_Text(pDX, IDC_EDIT9, alpha);
}


BEGIN_MESSAGE_MAP(Dial_p_val, CDialog)
	ON_BN_CLICKED(IDOK, &Dial_p_val::OnBnClickedOk)
END_MESSAGE_MAP()


// Dial_p_val message handlers

void Dial_p_val::OnBnClickedOk()
{
	UpdateData(TRUE);
	if (a <= 0 || b <= 0 || k <= 0 || k > a || k > b)
		AfxMessageBox(L"Parameter values are incorrect. Parameters must be positive and k<=a and k<=b.");
	else if (sample_sz< 50 || sample_sz>10'000)
		AfxMessageBox(L"select the number of simulated random variables between 100 and 10.000.");
	else if (samples_nmb < 1'000 || samples_nmb>1'000'000)
		AfxMessageBox(L"select the number of simulated random variables between 1.000 and 1.000.000.");
	else
		CDialog::OnOK();
}
