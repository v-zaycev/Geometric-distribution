#pragma once
#include "afxdialogex.h"


// Dial_power_n dialog

class Dial_power_n : public CDialog
{
	DECLARE_DYNAMIC(Dial_power_n)

public:
	Dial_power_n(CWnd* pParent = nullptr);   // standard constructor
	Dial_power_n(BOOL type, int a, int b, int k, int h_a, int h_b, int h_k, int init_sz, int step_sz, int steps_nmb, int sample_sz, double alpha, CWnd* pParent = nullptr);
	virtual ~Dial_power_n();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG3 };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	BOOL type;
	int init_sz;
	int step_sz;
	int steps_nmb;
	double alpha;
	int a;
	int b;
	int k;
	int h_a;
	int h_b;
	int h_k;
	afx_msg void OnBnClickedOk();
	int sample_sz;
};
