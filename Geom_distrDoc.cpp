
// Geom_distrDoc.cpp : implementation of the CGeomdistrDoc class
//

#include "pch.h"
#include "framework.h"
// SHARED_HANDLERS can be defined in an ATL project implementing preview, thumbnail
// and search filter handlers and allows sharing of document code with that project.
#ifndef SHARED_HANDLERS
#include "Geom_distr.h"
#endif

#include "Geom_distrDoc.h"
#include "Dial_basic.h"
#include "Dial_p_val.h"
#include "Dial_power_n.h"
#include <propkey.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CGeomdistrDoc

IMPLEMENT_DYNCREATE(CGeomdistrDoc, CDocument)

BEGIN_MESSAGE_MAP(CGeomdistrDoc, CDocument)
	ON_COMMAND(ID_SETPARAMETERS, &CGeomdistrDoc::OnSetparameters)
	ON_COMMAND(ID_BUTTON_P_VAL, &CGeomdistrDoc::OnCalcPVal)
	ON_COMMAND(ID_REBUILD, &CGeomdistrDoc::OnRebuild)
	ON_COMMAND(ID_PVALPARAM, &CGeomdistrDoc::OnPvalparam)
	ON_COMMAND(ID_POWER_DEP, &CGeomdistrDoc::OnNPowerDep)
	ON_COMMAND(ID_POWERDEPENDENCEONN, &CGeomdistrDoc::OnPowerdependenceonn)
END_MESSAGE_MAP()


// CGeomdistrDoc construction/destruction

CGeomdistrDoc::CGeomdistrDoc() noexcept
{
	// TODO: add one-time construction code here

}

CGeomdistrDoc::~CGeomdistrDoc()
{
}

BOOL CGeomdistrDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}




// CGeomdistrDoc serialization

void CGeomdistrDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

#ifdef SHARED_HANDLERS

// Support for thumbnails
void CGeomdistrDoc::OnDrawThumbnail(CDC& dc, LPRECT lprcBounds)
{
	// Modify this code to draw the document's data
	dc.FillSolidRect(lprcBounds, RGB(255, 255, 255));

	CString strText = _T("TODO: implement thumbnail drawing here");
	LOGFONT lf;

	CFont* pDefaultGUIFont = CFont::FromHandle((HFONT) GetStockObject(DEFAULT_GUI_FONT));
	pDefaultGUIFont->GetLogFont(&lf);
	lf.lfHeight = 36;

	CFont fontDraw;
	fontDraw.CreateFontIndirect(&lf);

	CFont* pOldFont = dc.SelectObject(&fontDraw);
	dc.DrawText(strText, lprcBounds, DT_CENTER | DT_WORDBREAK);
	dc.SelectObject(pOldFont);
}

// Support for Search Handlers
void CGeomdistrDoc::InitializeSearchContent()
{
	CString strSearchContent;
	// Set search contents from document's data.
	// The content parts should be separated by ";"

	// For example:  strSearchContent = _T("point;rectangle;circle;ole object;");
	SetSearchContent(strSearchContent);
}

void CGeomdistrDoc::SetSearchContent(const CString& value)
{
	if (value.IsEmpty())
	{
		RemoveChunk(PKEY_Search_Contents.fmtid, PKEY_Search_Contents.pid);
	}
	else
	{
		CMFCFilterChunkValueImpl *pChunk = nullptr;
		ATLTRY(pChunk = new CMFCFilterChunkValueImpl);
		if (pChunk != nullptr)
		{
			pChunk->SetTextValue(PKEY_Search_Contents, value, CHUNK_TEXT);
			SetChunkValue(pChunk);
		}
	}
}

#endif // SHARED_HANDLERS

// CGeomdistrDoc diagnostics

#ifdef _DEBUG
void CGeomdistrDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CGeomdistrDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CGeomdistrDoc commands


void CGeomdistrDoc::OnSetparameters()
{
	Dial_basic d(type, a_0, b_0, k_0, n);


	if (d.DoModal() == IDOK)
	{
		a_0 = d.a_val;
		b_0 = d.b_val;
		k_0 = d.k_val;
		n = d.n_val;
			
		if (type != d.func_type||!generator)
		{
			if (generator)
			{
				delete generator;
				generator = nullptr;
			}
			if (d.func_type == 0)
				generator = new Hypergeom_distr_inv(a_0, b_0, k_0);
			else
				generator = new Hypergeom_distr_bern(a_0, b_0, k_0);

			generator->set_hypothesis(a_1, b_1, k_1, sample_sz, samples_nmb, alpha);
			generator->set_power_n_param(init_sz, steps_nmb, step_sz, power_n_sample_sz);
		}
		type = d.func_type;
		
		generator->set_param(a_0, b_0, k_0, n);
	}
	// TODO: Add your command handler code here
}

void CGeomdistrDoc::OnCalcPVal()
{
	display_cond = 2;
	if (!generator)
	{
		if (type == 0)
			generator = new Hypergeom_distr_inv(a_0, b_0, k_0);
		else
			generator = new Hypergeom_distr_bern(a_0, b_0, k_0);

		generator->set_param(a_0, b_0, k_0, n);
		generator->set_hypothesis(a_1, b_1, k_1, sample_sz, samples_nmb, alpha);
		generator->set_power_n_param(init_sz, steps_nmb, step_sz, power_n_sample_sz);
	}
		
	generator->gen_p_levels();
	UpdateAllViews(0);
	return;
}

void CGeomdistrDoc::OnRebuild()
{
	display_cond = 1;
	if (!generator)
	{
		if (type == 0)
			generator = new Hypergeom_distr_inv(a_0, b_0, k_0);
		else
			generator = new Hypergeom_distr_bern(a_0, b_0, k_0);

		generator->set_param(a_0, b_0, k_0, n);
		generator->set_hypothesis(a_1, b_1, k_1, sample_sz, samples_nmb, alpha);
		generator->set_power_n_param(init_sz, steps_nmb, step_sz, power_n_sample_sz);
	}
	generator->gen_distr(n);
	generator->calc_distr();
	generator->calc_p_value();
	UpdateAllViews(0);
	return;
}

void CGeomdistrDoc::OnPvalparam()
{
	// TODO: Add your command handler code here
	Dial_p_val d(type, sample_sz, samples_nmb, a_0, b_0, k_0, a_1, b_1, k_1, alpha);


	if (d.DoModal() == IDOK)
	{
		a_0 = d.a;
		b_0 = d.b;
		k_0 = d.k;
		a_1 = d.h_a;
		b_1 = d.h_b;
		k_1 = d.h_k;
		sample_sz = d.sample_sz_p_val;
		samples_nmb = d.sample_sz;
		k_1 = d.h_k;
		alpha = d.alpha;

		if (type != d.type||!generator)
		{
			if (generator)
			{
				delete generator;
				generator = nullptr;
			}
			if (d.type == 0)
				generator = new Hypergeom_distr_inv(a_0, b_0, k_0);
			else
				generator = new Hypergeom_distr_bern(a_0, b_0, k_0);

			generator->set_power_n_param(init_sz, steps_nmb, step_sz, sample_sz);

		}
		type = d.type;

		generator->set_param(a_0, b_0, k_0, n);
		generator->set_hypothesis(a_1, b_1, k_1,sample_sz, samples_nmb, alpha);
	}
}

void CGeomdistrDoc::OnNPowerDep()
{
	display_cond = 3;
	if (!generator)
	{
		if (type == 0)
			generator = new Hypergeom_distr_inv(a_0, b_0, k_0);
		else
			generator = new Hypergeom_distr_bern(a_0, b_0, k_0);

		generator->set_param(a_0, b_0, k_0, n);
		generator->set_hypothesis(a_1, b_1, k_1, sample_sz, samples_nmb, alpha);
		generator->set_power_n_param(init_sz, steps_nmb, step_sz, power_n_sample_sz);
	}
	
	generator->power_n_dependence();
	UpdateAllViews(0);
	return;
}

void CGeomdistrDoc::OnPowerdependenceonn()
{
	Dial_power_n d(type, a_0, b_0, k_0, a_1, b_1, k_1, init_sz,step_sz, steps_nmb, power_n_sample_sz,alpha);


	if (d.DoModal() == IDOK)
	{
		a_0 = d.a;
		b_0 = d.b;
		k_0 = d.k;
		a_1 = d.h_a;
		b_1 = d.h_b;
		k_1 = d.h_k;
		init_sz = d.init_sz;
		steps_nmb = d.steps_nmb;
		power_n_sample_sz = d.sample_sz;
		step_sz = d.step_sz;
		alpha = d.alpha;

		if (type != d.type||!generator)
		{
			if (generator)
			{
				delete generator;
				generator = nullptr;
			}
			if (d.type == 0)
				generator = new Hypergeom_distr_inv(a_0, b_0, k_0);
			else
				generator = new Hypergeom_distr_bern(a_0, b_0, k_0);
		}
		type = d.type;

		generator->set_param(a_0, b_0, k_0, n);
		generator->set_hypothesis(a_1, b_1, k_1, sample_sz, samples_nmb, alpha);
		generator->set_power_n_param(init_sz, steps_nmb, step_sz, power_n_sample_sz);
	}
}
