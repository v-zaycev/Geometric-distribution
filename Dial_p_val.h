#pragma once
#include "afxdialogex.h"


// Dial_p_val dialog

class Dial_p_val : public CDialog
{
	DECLARE_DYNAMIC(Dial_p_val)

public:
	Dial_p_val(CWnd* pParent = nullptr);   // standard constructor
	Dial_p_val(BOOL _type, int _sample_sz_p_val, int _sample_sz, int _a, int _b, int _k, int _h_a, int _h_b, int _h_k, double alpha, CWnd* pParent = nullptr);
	virtual ~Dial_p_val();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG2 };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	BOOL type;
	int samples_nmb;
	int sample_sz;
	int a;
	int b;
	int k;
	int h_a;
	int h_b;
	int h_k;
	double alpha;
	afx_msg void OnBnClickedOk();
};
