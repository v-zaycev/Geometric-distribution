// Dial_basic.cpp : implementation file
//

#include "pch.h"
#include "Geom_distr.h"
#include "afxdialogex.h"
#include "Dial_basic.h"


// Dial_basic dialog

IMPLEMENT_DYNAMIC(Dial_basic, CDialogEx)

Dial_basic::Dial_basic(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_DIALOG1, pParent)
	, func_type(FALSE)
	, a_val(1)
	, b_val(1)
	, k_val(1)
	, n_val(100)
{

}

Dial_basic::Dial_basic(BOOL _func_type, int _a, int _b, int _k, int _n, CWnd* pParent)
	: CDialogEx(IDD_DIALOG1, pParent)
	, func_type(_func_type)
	, a_val(_a)
	, b_val(_b)
	, k_val(_k)
	, n_val(_n)

{
}

Dial_basic::~Dial_basic()
{
}

void Dial_basic::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Radio(pDX, IDC_RADIO1, func_type);
	DDX_Text(pDX, IDC_EDIT1, a_val);
	DDX_Text(pDX, IDC_EDIT2, b_val);
	DDX_Text(pDX, IDC_EDIT3, k_val);
	DDX_Text(pDX, IDC_EDIT4, n_val);
}


BEGIN_MESSAGE_MAP(Dial_basic, CDialogEx)
	ON_BN_CLICKED(IDOK, &Dial_basic::OnBnClickedOk)
END_MESSAGE_MAP()


// Dial_basic message handlers


void Dial_basic::OnBnClickedOk()
{
	UpdateData(TRUE);
	if (a_val <= 0 || b_val <= 0 || k_val <= 0 || k_val > a_val || k_val > b_val)
		AfxMessageBox(L"Parameter values are incorrect. Parameters must be positive and k<=a and k<=b.");
	else if (n_val<100||n_val>100'000'000)
		AfxMessageBox(L"select the number of simulated random variables between 100 and 1.000.000.000.");
	else
		CDialogEx::OnOK();
}

