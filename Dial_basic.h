#pragma once
#include "afxdialogex.h"


// Dial_basic dialog

class Dial_basic : public CDialogEx
{
	DECLARE_DYNAMIC(Dial_basic)

public:
	Dial_basic(CWnd* pParent = nullptr);   // standard constructor
	Dial_basic(BOOL _func_type, int _a, int _b, int _k, int _n, CWnd* pParent = nullptr);
	virtual ~Dial_basic();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG1 };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	BOOL func_type;
	int a_val;
	int b_val;
	int k_val;
	int n_val;
	afx_msg void OnBnClickedOk();
};
