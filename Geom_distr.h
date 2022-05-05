
// Geom_distr.h : main header file for the Geom_distr application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'pch.h' before including this file for PCH"
#endif

#include "resource.h"       // main symbols


// CGeomdistrApp:
// See Geom_distr.cpp for the implementation of this class
//

class CGeomdistrApp : public CWinAppEx
{
public:
	CGeomdistrApp() noexcept;


// Overrides
public:
	virtual BOOL InitInstance();

// Implementation
	UINT  m_nAppLook;
	BOOL  m_bHiColorIcons;

	virtual void PreLoadState();
	virtual void LoadCustomState();
	virtual void SaveCustomState();

	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CGeomdistrApp theApp;
