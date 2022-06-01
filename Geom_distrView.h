
// Geom_distrView.h : interface of the CGeomdistrView class
//

#pragma once


class CGeomdistrView : public CView
{
protected: // create from serialization only
	CGeomdistrView() noexcept;
	DECLARE_DYNCREATE(CGeomdistrView)

// Attributes
public:
	CGeomdistrDoc* GetDocument() const;

// Operations
public:

// Overrides
public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	void OnDraw_basic(CDC* pDC);
	void OnDraw_p_val(CDC* pDC);
	void OnDraw_power(CDC* pDC);
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:

// Implementation
public:
	virtual ~CGeomdistrView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in Geom_distrView.cpp
inline CGeomdistrDoc* CGeomdistrView::GetDocument() const
   { return reinterpret_cast<CGeomdistrDoc*>(m_pDocument); }
#endif

