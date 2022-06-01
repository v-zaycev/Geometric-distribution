
// Geom_distrDoc.h : interface of the CGeomdistrDoc class
//


#pragma once
#include "Distribution.h"

class CGeomdistrDoc : public CDocument
{
protected: // create from serialization only
	CGeomdistrDoc() noexcept;
	DECLARE_DYNCREATE(CGeomdistrDoc)

// Attributes
public:

// Operations
public:

// Overrides
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

// Implementation
public:
	virtual ~CGeomdistrDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// Helper function that sets search content for a Search Handler
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS

public:
	Hypergeom_distr h0;
	Hypergeom_distr h1;
	Hypergeom_sample* generator;
	Chi_sq processor;
	
	int display_cond=0;
	COLORREF th_distr = RGB(180, 180, 180);
	COLORREF mod_distr = RGB(100, 100, 100);
	COLORREF basic_line = RGB(0, 0, 0);
	COLORREF p_val_line = RGB(255, 0, 0);
	COLORREF power_n_line = RGB(0, 255, 0);
	COLORREF background = RGB(235, 235, 235);
	afx_msg void OnSetparameters();
	afx_msg void OnCalcPVal();
	afx_msg void OnRebuild();
	afx_msg void OnPvalparam();
	afx_msg void OnNPowerDep();
	afx_msg void OnPowerdependenceonn();
};
