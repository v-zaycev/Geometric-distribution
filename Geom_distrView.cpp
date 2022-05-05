
// Geom_distrView.cpp : implementation of the CGeomdistrView class
//

#include "pch.h"
#include "framework.h"
#include<string>
// SHARED_HANDLERS can be defined in an ATL project implementing preview, thumbnail
// and search filter handlers and allows sharing of document code with that project.
#ifndef SHARED_HANDLERS
#include "Geom_distr.h"
#endif

#include "Geom_distrDoc.h"
#include "Geom_distrView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CGeomdistrView

IMPLEMENT_DYNCREATE(CGeomdistrView, CView)

BEGIN_MESSAGE_MAP(CGeomdistrView, CView)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
END_MESSAGE_MAP()

// CGeomdistrView construction/destruction

CGeomdistrView::CGeomdistrView() noexcept
{
	// TODO: add construction code here

}

CGeomdistrView::~CGeomdistrView()
{
}

BOOL CGeomdistrView::PreCreateWindow(CREATESTRUCT& cs)
{
	cs.cx = 1280;
	cs.cy = 1720;

	return CView::PreCreateWindow(cs);
}

// CGeomdistrView drawing

void CGeomdistrView::OnDraw(CDC* pDC)
{
	CGeomdistrDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	CRect wnd;
	CPen pen_axis, pen_p_val, pen_power_n, zero_pen;
	zero_pen.CreatePen(PS_SOLID, 0, pDoc->background);
	pen_axis.CreatePen(PS_SOLID, 2, pDoc->basic_line);
	pen_p_val.CreatePen(PS_SOLID, 2, pDoc->p_val_line);
	pen_power_n.CreatePen(PS_SOLID, 2, pDoc->power_n_line);
	CBrush brush_th;
	CBrush brush_mod;
	CBrush background;
	brush_th.CreateSolidBrush(pDoc->th_distr);
	brush_mod.CreateSolidBrush(pDoc->mod_distr);
	background.CreateSolidBrush(pDoc->background);
	GetClientRect(&wnd);
	double wdz = wnd.Width(), hgt = wnd.Height();
	double out_sp_x = wdz * 0.8, out_sp_y = hgt * 0.8;
	double x_0 = wdz * 0.05, x_1 = wdz * 0.85;
	double y_0 = hgt * 0.1, y_1 = hgt * 0.9;


	// строим оси

	switch (pDoc->display_cond)
	{
	case 1: //dstribution
	{
		if (!pDoc->generator)
			break;
		const double* th_distr = pDoc->generator->get_th_distr();
		const int* mod_distr = pDoc->generator->get_gen_distr();
		double max_val = 0;
		const int k = pDoc->generator->get_k();
		double column_1 = out_sp_x / (k + 1), column_2 = column_1 * 0.6;

		for (int i = 0; i < k+1; ++i)
		{
			if (th_distr[i] > max_val)
				max_val = th_distr[i];
			if (mod_distr[i]/double(pDoc->n) > max_val)
				max_val = mod_distr[i] / double(pDoc->n);
		}

		CBrush l1;
		CBrush l2;
		l1.CreateSolidBrush(pDoc->th_distr);
		l2.CreateSolidBrush(pDoc->mod_distr);

		pDC->SelectObject(l1);
		pDC->Rectangle(x_1 + wdz * 0.0125, y_0 + hgt * 0.05, x_1 + wdz * 0.0125 + hgt * 0.07, y_0 + hgt * 0.09);
		pDC->TextOutW(x_1 + wdz * 0.0125 + hgt * 0.09, y_0 + hgt * 0.06, L"Theoretical");

		pDC->SelectObject(l2);
		pDC->Rectangle(x_1 + wdz * 0.0125, y_0 + hgt * 0.1, x_1 + wdz * 0.0125 + hgt * 0.07, y_0 + hgt * 0.14);
		pDC->TextOutW(x_1 + wdz * 0.0125 + hgt * 0.09, y_0 + hgt * 0.11, L"Modeled");


		pDC->SelectObject(zero_pen);
		pDC->SelectObject(background);
		pDC->Rectangle(x_0, y_0, x_1, y_1);

		pDC->SelectObject(pen_axis);

		// строим оси
		pDC->MoveTo(x_0, y_1);
		pDC->LineTo(x_1, y_1);

		pDC->MoveTo(x_0, y_1);
		pDC->LineTo(x_0, y_0);

		pDC->MoveTo(x_1, y_1);
		pDC->LineTo(x_1 - 5, y_1 - 5);

		pDC->MoveTo(x_1, y_1);
		pDC->LineTo(x_1 - 5, y_1 + 5);

		pDC->MoveTo(x_0, y_0);
		pDC->LineTo(x_0 + 5, y_0 + 5);

		pDC->MoveTo(x_0, y_0);
		pDC->LineTo(x_0 - 5, y_0 + 5);

		pDC->SelectObject(brush_th);
		for (int i = 0; i < k+1; ++i)
			pDC->Rectangle(x_0 + i * column_1, y_1, x_0 + (i + 1) * column_1, y_1 - th_distr[i] * out_sp_y/ max_val);

		pDC->SelectObject(brush_mod);
		for (int i = 0; i < k+1; ++i)
			pDC->Rectangle(x_0 + (i + 0.2) * column_1, y_1, x_0 + (i + 1 - 0.2) * column_1, y_1 - mod_distr[i] / double(pDoc->n) * out_sp_y / max_val);


		for (int i = 0; i < k + 1; ++i)
			pDC->TextOutW(x_0 + (i + 0.5) * column_1 - 10, y_1 + 10, std::to_wstring(i).data());

		pDC->TextOutW(x_0-40, y_0, std::to_wstring(int(pDoc->n*max_val)).data());
		pDC->TextOutW(x_0 - 40, y_1, L"0");

		pDC->TextOutW(wdz * 0.1, hgt * 0.05, L"d.f. = ");
		pDC->TextOutW(wdz * 0.1 + 40, hgt * 0.05, std::to_wstring(pDoc->generator->get_d_f()).data());

		pDC->TextOutW(wdz * 0.1 + 140, hgt * 0.05, L"chi^2 = ");
		pDC->TextOutW(wdz * 0.1 + 190, hgt * 0.05, std::to_wstring(pDoc->generator->get_chi_sq()).substr(0,6).data());

		pDC->TextOutW(wdz * 0.1 + 290, hgt * 0.05, L"p-value = ");
		pDC->TextOutW(wdz * 0.1 + 350, hgt * 0.05, std::to_wstring(pDoc->generator->get_p_val()).substr(0, 6).data());

		break;
	}
	case 2: //p-levels
	{
		const double* p_levels = pDoc->generator->get_p_levels();
		CBrush l1;
		CBrush l2;
		l1.CreateSolidBrush(pDoc->p_val_line);
		l2.CreateSolidBrush(pDoc->basic_line);

		pDC->SelectObject(l1);
		pDC->Rectangle(x_1 + wdz * 0.0125, y_0 + hgt * 0.05, x_1 + wdz * 0.0125 + hgt * 0.07, y_0 + hgt * 0.09);
		pDC->TextOutW(x_1 + wdz * 0.0125 + hgt * 0.09, y_0 + hgt * 0.06, L"Real");

		pDC->SelectObject(l2);
		pDC->Rectangle(x_1 + wdz * 0.0125, y_0 + hgt * 0.1, x_1 + wdz * 0.0125 + hgt * 0.07, y_0 + hgt * 0.14);
		pDC->TextOutW(x_1 + wdz * 0.0125 + hgt * 0.09, y_0 + hgt * 0.11, L"Expected");

		pDC->SelectObject(zero_pen);
		pDC->SelectObject(background);
		pDC->Rectangle(x_0, y_0, x_1, y_1);

		pDC->SelectObject(pen_axis);

		// строим оси
		pDC->MoveTo(x_0, y_1);
		pDC->LineTo(x_1, y_1);

		pDC->MoveTo(x_0, y_1);
		pDC->LineTo(x_0, y_0);

		pDC->MoveTo(x_1, y_1);
		pDC->LineTo(x_1 - 5, y_1 - 5);

		pDC->MoveTo(x_1, y_1);
		pDC->LineTo(x_1 - 5, y_1 + 5);

		pDC->MoveTo(x_0, y_0);
		pDC->LineTo(x_0 + 5, y_0 + 5);

		pDC->MoveTo(x_0, y_0);
		pDC->LineTo(x_0 - 5, y_0 + 5);

		pDC->MoveTo(x_0, y_1);
		pDC->LineTo(x_1, y_0);

		for (int i = 1; i < 20; ++i)
		{
			pDC->MoveTo(x_0 + out_sp_x * double(i) / 20, y_1 - 7);
			pDC->LineTo(x_0 + out_sp_x * double(i) / 20, y_1 + 7);
			pDC->MoveTo(x_0 - 5, y_1 - out_sp_y * double(i) / 20);
			pDC->LineTo(x_0 + 5, y_1 - out_sp_y * double(i) / 20);
		}

		pDC->TextOutW(x_0 - 10, y_1 + 10, L"0");

		for (int i = 1; i < 21; ++i)
		{
			pDC->TextOutW(x_0 + out_sp_x * double(i) / 20 - 10, y_1 + 10, std::to_wstring(double(i) / 20).substr(0,4).data());
			pDC->TextOutW(x_0 - 40, y_1 - out_sp_y * double(i) / 20-5, std::to_wstring(double(i) / 20).substr(0, 4).data());
		}

		pDC->SelectObject(pen_p_val);
		pDC->MoveTo(x_0, y_1);
		for (int i = 0; i < 101; ++i)
			pDC->LineTo(x_0 + out_sp_x * double(i) / 100, y_1 - out_sp_y * p_levels[i]);
		
		pDC->TextOutW(wdz * 0.1, hgt * 0.05, L"Power = ");
		pDC->TextOutW(wdz * 0.15, hgt * 0.05, std::to_wstring(pDoc->generator->get_power()).substr(0,5).data());

		break;
	}
	case 3: //power/n dependance
	{
		if (!pDoc->generator)
			break;

		const double* power_n_dep = pDoc->generator->get_power_n();
		const int n_steps = pDoc->generator->get_steps();
		double step = out_sp_x / n_steps;
		double max_val = 0;

		for (int i = 0; i < n_steps + 1; ++i)
			if (power_n_dep[i] > max_val)
				max_val = power_n_dep[i];

		CBrush l1;
		l1.CreateSolidBrush(pDoc->power_n_line);

		pDC->SelectObject(l1);
		pDC->Rectangle(x_1 + wdz * 0.0125, y_0 + hgt * 0.05, x_1 + wdz * 0.0125 + hgt * 0.07, y_0 + hgt * 0.09);
		pDC->TextOutW(x_1 + wdz * 0.0125 + hgt * 0.09, y_0 + hgt * 0.06, L"Power");


		pDC->SelectObject(zero_pen);
		pDC->SelectObject(background);
		pDC->Rectangle(x_0, y_0, x_1, y_1);

		pDC->SelectObject(pen_axis);

		// строим оси
		pDC->MoveTo(x_0, y_1);
		pDC->LineTo(x_1, y_1);

		pDC->MoveTo(x_0, y_1);
		pDC->LineTo(x_0, y_0);

		pDC->MoveTo(x_1, y_1);
		pDC->LineTo(x_1 - 5, y_1 - 5);

		pDC->MoveTo(x_1, y_1);
		pDC->LineTo(x_1 - 5, y_1 + 5);

		pDC->MoveTo(x_0, y_0);
		pDC->LineTo(x_0 + 5, y_0 + 5);

		pDC->MoveTo(x_0, y_0);
		pDC->LineTo(x_0 - 5, y_0 + 5);

		if (pDoc->a_0 == pDoc->a_1 && pDoc->b_0 == pDoc->b_1 && pDoc->k_0 == pDoc->k_1)
		{
			pDC->MoveTo(x_0, y_1 - out_sp_y * pDoc->alpha / max_val);
			pDC->LineTo(x_1, y_1 - out_sp_y * pDoc->alpha / max_val);
			pDC->TextOutW(x_0 - 50, y_1 - out_sp_y * pDoc->alpha / max_val, std::to_wstring(pDoc->generator->get_alpha()).substr(0, 5).data());
		}

		pDC->SelectObject(pen_power_n);
		pDC->MoveTo(x_0, y_1 - out_sp_y * power_n_dep[0] / max_val);
		for (int i = 1; i < n_steps + 1; ++i)
			pDC->LineTo(x_0 + i * step, y_1 - out_sp_y * power_n_dep[i] / max_val);

		for (int i = 0; i < n_steps + 1; ++i)
			pDC->TextOutW(x_0 + i * step, y_1 + 10, std::to_wstring(pDoc->generator->get_init_pos()+i*pDoc->generator->get_step()).data());

		pDC->TextOutW(x_0 - 50, y_0, std::to_wstring(max_val).substr(0, 5).data());
		pDC->TextOutW(x_0 - 20, y_1, L"0");

		break;
	}
	}


}

void CGeomdistrView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CGeomdistrView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CGeomdistrView diagnostics

#ifdef _DEBUG
void CGeomdistrView::AssertValid() const
{
	CView::AssertValid();
}

void CGeomdistrView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CGeomdistrDoc* CGeomdistrView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CGeomdistrDoc)));
	return (CGeomdistrDoc*)m_pDocument;
}
#endif //_DEBUG


// CGeomdistrView message handlers
