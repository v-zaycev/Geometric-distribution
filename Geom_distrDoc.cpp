
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
	generator = new Hypergeom_inv();
}

CGeomdistrDoc::~CGeomdistrDoc()
{
	delete generator;
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
	Dial_basic d(generator->get_type(), generator->get_a(), generator->get_b(), generator->get_k(), generator->get_n());


	if (d.DoModal() == IDOK)
	{
		if (!generator || generator->get_type() != d.func_type)
		{
			if (generator)
			{
				delete generator;
				generator = nullptr;
			}
			if (d.func_type == 0)
				generator = new Hypergeom_inv(d.a_val, d.b_val, d.k_val, d.n_val);
			else
				generator = new Hypergeom_bern(d.a_val, d.b_val, d.k_val, d.n_val);
		}
		else
			generator->set_all(d.a_val, d.b_val, d.k_val, d.n_val);
		h0.set_param(d.a_val, d.b_val, d.k_val);
	}
}

void CGeomdistrDoc::OnCalcPVal()
{
	display_cond = 2;
	if (!generator)
		generator = new Hypergeom_inv();

	processor.gen_p_levels(generator,h0);
	UpdateAllViews(0);
	return;
}

void CGeomdistrDoc::OnRebuild()
{
	display_cond = 1;
	if (!generator)
		generator = new Hypergeom_inv();

	generator->gen_sample();
	h0.calc_probs();
	processor.calc_p_val(generator, h0);
	UpdateAllViews(0);
	return;
}

void CGeomdistrDoc::OnPvalparam()
{
	Dial_p_val d(generator->get_type(), 
				 processor.get_smpls_nmb(),
				 processor.get_smpl_sz(), 
				 generator->get_a(), 
				 generator->get_b(), 
				 generator->get_k(), 
				 h0.get_a(), 
				 h0.get_b(), 
				 h0.get_k(), 
				 processor.get_alpha());


	if (d.DoModal() == IDOK)
	{
		if (!generator || generator->get_type() != d.type)
		{
			int n = 100;
			if (generator)
			{
				n = generator->get_n();
				delete generator;
				generator = nullptr;
			}
			if (d.type == 0)
				generator = new Hypergeom_inv(d.a, d.b, d.k, n);
			else
				generator = new Hypergeom_bern(d.a, d.b, d.k, n);
		}
		else
			generator->set_param(d.a, d.b, d.k);
		h0.set_param(d.h_a, d.h_b, d.h_k);
		processor.set_p_lvls( d.sample_sz, d.samples_nmb, d.alpha);
	}
}

void CGeomdistrDoc::OnNPowerDep()
{
	display_cond = 3;
	if (!generator)
		generator = new Hypergeom_inv();
	
	processor.power_n_dependence(generator, h0);
	UpdateAllViews(0);
	return;
}

void CGeomdistrDoc::OnPowerdependenceonn()
{
	Dial_power_n d(generator->get_type(), 
				   generator->get_a(), 
				   generator->get_b(), 
				   generator->get_k(), 
			       h0.get_a(), 
		           h0.get_b(), 
				   h0.get_k(), 
				   processor.get_start_pos(),
				   processor.get_step_sz(),
				   processor.get_steps_nmb(),
				   processor.get_power_n_samples_nmb(),
				   processor.get_alpha());


	if (d.DoModal() == IDOK)
	{
		if (!generator || generator->get_type() != d.type)
		{
			int n = 100;
			if (generator)
			{
				n = generator->get_n();
				delete generator;
				generator = nullptr;
			}
			if (d.type == 0)
				generator = new Hypergeom_inv(d.a, d.b, d.k, n);
			else
				generator = new Hypergeom_bern(d.a, d.b, d.k, n);
		}
		else
			generator->set_param(d.a, d.b, d.k);
		h0.set_param(d.h_a, d.h_b, d.h_k);
		processor.set_power_n(d.init_sz, d.steps_nmb, d.step_sz, d.sample_sz, d.alpha);
	}
}
