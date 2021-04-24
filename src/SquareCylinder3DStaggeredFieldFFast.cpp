#include "SquareCylinder3DStaggeredFieldFFast.h"
#include "LSSolverF.h"

#include <string>

SquareCylinder3DStaggeredFieldFFast::SquareCylinder3DStaggeredFieldFFast(float rho,
	float mu,
	float h,
	float l,
	float L,
	float H,
	float W,
	float ufreestream)
	:
	crho(rho),
	cmu(mu),
	ch(h),
	cL(L),
	cH(H),
	cW(W),

	cl(l),

	cFirstIteration(true),

	cDistRelFrenteRect(5),

	cNXU(int(cL / ch) + 1),
	cNYU(int(cW / ch) + 2),
	cNZU(int(cH / ch) + 2),

	cNXV(int(cL / ch) + 2),
	cNYV(int(cW / ch) + 1),
	cNZV(int(cH / ch) + 2),

	cNXW(int(cL / ch) + 2),
	cNYW(int(cW / ch) + 2),
	cNZW(int(cH / ch) + 1),

	cNXP(int(cL / ch) + 2),
	cNYP(int(cW / ch) + 2),
	cNZP(int(cH / ch) + 2),

	cNXCy(cl / ch),
	cNYCy(int(cNYP * 0.8)),
	cNZCy(cNXCy),

	ciCyInit(cDistRelFrenteRect* cNXCy),
	cjCyInit(cNYP / 2 - cNYCy / 2),
	ckCyInit(cNZP / 2 - cNZCy / 2),

	ciCyEnd(cDistRelFrenteRect* cNXCy + cNXCy),
	cjCyEnd((cNYP / 2 + cNYCy / 2) + (cNYCy % 2)),
	ckCyEnd((cNZP / 2 + cNZCy / 2) + (cNZCy % 2)),

	cUfreestream(ufreestream)
{
	cUn.resize(cNXU*cNYU*cNZU, cUfreestream);
	cVn.resize(cNXV*cNYV*cNZV, 0.0);
	cWn.resize(cNXW*cNYW*cNZW, 0.0);
	cPn.resize(cNXP*cNYP*cNZP, 0.0);

	cUnp = cUn;
	cVnp = cVn;
	cWnp = cWn;
	//cPnp = cPn;

	cPhi = cPn;

	cdphidx = cPn;
	cdphidy = cPn;
	cdphidz = cPn;

	cd2phidx2 = cPn;
	cd2phidy2 = cPn;
	cd2phidz2 = cPn;

	cdUdtn = cUn;
	cdUdt1 = cUn;
	cdUdt2 = cUn;
	cdUdt3 = cUn;

	cdVdtn = cVn;
	cdVdt1 = cVn;
	cdVdt2 = cVn;
	cdVdt3 = cVn;

	cdWdtn = cWn;
	cdWdt1 = cWn;
	cdWdt2 = cWn;
	cdWdt3 = cWn;
}

SquareCylinder3DStaggeredFieldFFast::~SquareCylinder3DStaggeredFieldFFast()
{}

void SquareCylinder3DStaggeredFieldFFast::SetBCFlatPlate()
{
	// freestream in all

	/*for (int j = 1; j < cNYU-1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			//cUn[0 + (k * cNXU) + (j*cNXU*cNZU)] = cUfreestream;
			cUn[1 + (k * cNXU) + (j*cNXU*cNZU)] = cUfreestream;
		}
	}*/

	for (int j = cjCyInit; j < cjCyEnd; j++)
		for (int k = ckCyInit; k < ckCyEnd; k++)
			for (int i = ciCyInit - 1; i < ciCyEnd; i++)
			{
				cUn[i + (k*cNXU) + (j*cNXU*cNZU)] = 0.0;
			}

	for (int j = cjCyInit - 1; j < cjCyEnd; j++)
		for (int k = ckCyInit; k < ckCyEnd; k++)
			for (int i = ciCyInit; i < ciCyEnd; i++)
			{
				cVn[i + (k*cNXV) + (j*cNXV*cNZV)] = 0.0;
			}

	for (int j = cjCyInit; j < cjCyEnd; j++)
		for (int k = ckCyInit - 1; k < ckCyEnd; k++)
			for (int i = ciCyInit; i < ciCyEnd; i++)
			{
				cWn[i + (k*cNXW) + (j*cNXW*cNZW)] = 0.0;
			}

	cUnp = cUn;
	cVnp = cVn;
	cWnp = cWn;
}

void SquareCylinder3DStaggeredFieldFFast::UpdateBCFlatPlate()
{

	// extrapolate outlet (i = NX - 1)

	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			{
				const int i = cNXU - 1;
				const int iu = cNXU - 2;

				cUnp[(j * cNZU * cNXU) + (k * cNXU) + iu + 1] = cUnp[(j * cNZU * cNXU) + (k * cNXU) + iu] +
					cVnp[((j - 1) * cNZV * cNXV) + (k * cNXV) + i] - cVnp[((j)* cNZV * cNXV) + (k * cNXV) + i] +
					cWnp[(j * cNZW * cNXW) + ((k - 1) * cNXW) + i] - cWnp[(j * cNZW * cNXW) + ((k)* cNXW) + i];
			}

			//cUnp[cNXU - 1 + (k * cNXU) + (j * cNXU * cNZU)] = cUnp[cNXU - 2 + (k * cNXU) + (j * cNXU * cNZU)];

		}
	}

	/*for (int j = 1; j < cNYV - 1; j++)
	{
		for (int k = 1; k < cNZV - 1; k++)
		{
			cVn[cNXV - 1 + (k * cNXV) + (j * cNXV * cNZV)] = cVn[cNXV - 2 + (k * cNXV) + (j * cNXV * cNZV)];
		}
	}

	for (int j = 1; j < cNYW - 1; j++)
	{
		for (int k = 1; k < cNZW - 1; k++)
		{
			cWn[cNXW - 1 + (k * cNXW) + (j * cNXW * cNZW)] = cWn[cNXW - 2 + (k * cNXW) + (j * cNXW * cNZW)];
		}
	}*/


	/*for(int k = 1; k < cNZU - 1; k++)
		for (int i = 1; i < cNXU - 1; i++)
		{
			cUnp[i + (k * cNXU) + ((cNYU - 1) * cNXU * cNZU)] = cUnp[i + (k * cNXU) + ((cNYU - 2) * cNXU * cNZU)];

			cUnp[i + (k * cNXU) + (0 * cNXU * cNZU)] = cUnp[i + (k * cNXU) + (1 * cNXU * cNZU)];
		}


	for (int k = 1; k < cNZW - 1; k++)
		for (int i = 1; i < cNXW - 1; i++)
		{
			cWnp[i + (k * cNXW) + ((cNYW - 1) * cNXW * cNZW)] = cWnp[i + (k * cNXW) + ((cNYW - 2) * cNXW * cNZW)];
			cWnp[i + (k * cNXW) + (0 * cNXW * cNZW)] = cWnp[i + (k * cNXW) + (1 * cNXW * cNZW)];
		}



	for (int j = 1; j < cNYU - 1; j++)
		for (int i = 1; i < cNXU - 1; i++)
		{
			cUnp[i + ((cNZU - 1) * cNXU) + (j * cNXU * cNZU)] = cUnp[i + ((cNZU - 2) * cNXU) + (j * cNXU * cNZU)];
			cUnp[i + (0 * cNXU) + (j * cNXU * cNZU)] = cUnp[i + (1 * cNXU) + (j * cNXU * cNZU)];
		}

	for (int j = 1; j < cNYV - 1; j++)
		for (int i = 1; i < cNXV - 1; i++)
		{
			cVnp[i + ((cNZV - 1) * cNXW) + (j * cNXW * cNZW)] = cVnp[i + ((cNZV - 2) * cNXW) + (j * cNXW * cNZW)];
			cVnp[i + (0 * cNXW) + (j * cNXW * cNZW)] = cVnp[i + (1 * cNXW) + (j * cNXW * cNZW)];
		}*/
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidxSystemTDUNUSED(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * phi[i + (k * NX) + (j * NX * NZ)] + 4 * phi[i + 1 + (k * NX) + (j * NX * NZ)] - phi[i + 2 + (k * NX) + (j * NX * NZ)]);
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * phi[i + (k * NX) + (j * NX * NZ)] - 4 * phi[i - 1 + (k * NX) + (j * NX * NZ)] + phi[i - 2 + (k * NX) + (j * NX * NZ)]);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*phi[i - 1 + (k * NX) + (j * NX * NZ)] - (-8 * 0.025 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (-4 * 0.025 + 1.0)*phi[i + 1 + (k * NX) + (j * NX * NZ)] - (-8 * 0.025 / 3.0 + 1.0 / 6.0)*phi[i + 2 + (k * NX) + (j * NX * NZ)] + (-2 * 0.025 / 3.0*phi[i + 3 + (k * NX) + (j * NX * NZ)]));
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*phi[i + 1 + (k * NX) + (j * NX * NZ)] - (8 * 0.09 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (4 * 0.09 + 1.0)*phi[i - 1 + (k * NX) + (j * NX * NZ)] - (8 * 0.09 / 3.0 + 1.0 / 6.0)*phi[i - 2 + (k * NX) + (j * NX * NZ)] + (2 * 0.09 / 3.0*phi[i - 3 + (k * NX) + (j * NX * NZ)]));
				}
				else
				{
					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*phi[i - 2 + (k * NX) + (j * NX * NZ)] - Eimp / 2.0*phi[i - 1 + (k * NX) + (j * NX * NZ)] + Eimp / 2.0*phi[i + 1 + (k * NX) + (j * NX * NZ)] + Fimp / 4.0*phi[i + 2 + (k * NX) + (j * NX * NZ)]);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidySystemTDUNUSED(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * phi[i + (k * NX) + (j * NX * NZ)] + 4 * phi[i + (k * NX) + ((j + 1) * NX * NZ)] - phi[i + (k * NX) + ((j + 2) * NX * NZ)]);
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * phi[i + (k * NX) + (j * NX * NZ)] - 4 * phi[i + (k * NX) + ((j - 1) * NX * NZ)] + phi[i + (k * NX) + ((j - 2) * NX * NZ)]);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*phi[i + (k * NX) + ((j - 1) * NX * NZ)] - (-8 * 0.025 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (-4 * 0.025 + 1.0)*phi[i + (k * NX) + ((j + 1)* NX * NZ)] - (-8 * 0.025 / 3.0 + 1.0 / 6.0)*phi[i + (k * NX) + ((j + 2) * NX * NZ)] + (-2 * 0.025 / 3.0*phi[i + (k * NX) + ((j + 3) * NX * NZ)]));
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*phi[i + (k * NX) + ((j + 1) * NX * NZ)] - (8 * 0.09 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (4 * 0.09 + 1.0)*phi[i + (k * NX) + ((j - 1) * NX * NZ)] - (8 * 0.09 / 3.0 + 1.0 / 6.0)*phi[i + (k * NX) + ((j - 2) * NX * NZ)] + (2 * 0.09 / 3.0*phi[i + (k * NX) + ((j - 3) * NX * NZ)]));
				}
				else
				{
					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*phi[i + (k * NX) + ((j - 2) * NX * NZ)] - Eimp / 2.0*phi[i + (k * NX) + ((j - 1) * NX * NZ)] + Eimp / 2.0*phi[i + (k * NX) + ((j + 1) * NX * NZ)] + Fimp / 4.0*phi[i + (k * NX) + ((j + 2) * NX * NZ)]);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidzSystemTDUNUSED(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * phi[i + (k * NX) + (j * NX * NZ)] + 4 * phi[i + ((k + 1) * NX) + ((j)* NX * NZ)] - phi[i + ((k + 2) * NX) + ((j)* NX * NZ)]);
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * phi[i + (k * NX) + (j * NX * NZ)] - 4 * phi[i + ((k - 1) * NX) + ((j)* NX * NZ)] + phi[i + ((k - 2) * NX) + ((j)* NX * NZ)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*phi[i + ((k - 1) * NX) + ((j)* NX * NZ)] - (-8 * 0.025 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (-4 * 0.025 + 1.0)*phi[i + ((k + 1) * NX) + ((j)* NX * NZ)] - (-8 * 0.025 / 3.0 + 1.0 / 6.0)*phi[i + ((k + 2) * NX) + ((j)* NX * NZ)] + (-2 * 0.025 / 3.0*phi[i + ((k + 3) * NX) + ((j)* NX * NZ)]));
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*phi[i + ((k + 1) * NX) + ((j)* NX * NZ)] - (8 * 0.09 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (4 * 0.09 + 1.0)*phi[i + ((k - 1) * NX) + ((j)* NX * NZ)] - (8 * 0.09 / 3.0 + 1.0 / 6.0)*phi[i + ((k - 2) * NX) + ((j)* NX * NZ)] + (2 * 0.09 / 3.0*phi[i + ((k - 3) * NX) + ((j)* NX * NZ)]));
				}
				else
				{
					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*phi[i + ((k - 2) * NX) + ((j)* NX * NZ)] - Eimp / 2.0*phi[i + ((k - 1) * NX) + ((j)* NX * NZ)] + Eimp / 2.0*phi[i + ((k + 1) * NX) + ((j)* NX * NZ)] + Fimp / 4.0*phi[i + ((k + 2) * NX) + ((j)* NX * NZ)]);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::Created2phidx2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i - 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + 1 + (k*NX) + (j*NX*NZ)])
						/(ch*ch);
				}
				else*/
				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + 1 + (k*NX) + (j*NX*NZ)] +
						4.0 * phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i + 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i - 1 + (k*NX) + (j*NX*NZ)] +
						4.0 * phi[i - 2 + (k*NX) + (j*NX*NZ)] - phi[i - 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + 2 + (k*NX) + (j*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i - 2 + (k*NX) + (j*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i - 1 + (k*NX) + (j*NX*NZ)]);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::Created2phidy2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j-1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + (k*NX) + ((j+1)*NX*NZ)])
						/ (ch*ch);
				}
				else */if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j + 1)*NX*NZ)] +
						4.0 * phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j + 3)*NX*NZ)]) / (ch*ch);
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j - 1)*NX*NZ)] +
						4.0 * phi[i + (k*NX) + ((j - 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 3)*NX*NZ)]) / (ch*ch);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + (k*NX) + ((j + 2)*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j - 2)*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j - 1)*NX*NZ)]);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::Created2phidz2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k-1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + ((k+1)*NX) + ((j)*NX*NZ)])
						/ (ch*ch);
				}
				else*/ if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
						4.0 * phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] +
						4.0 * phi[i + ((k - 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k - 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k - 2)*NX) + ((j)*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreateSparsePressureEquation(const field3D& dUdx,
	const field3D& dVdy,
	const field3D& dWdz,
	const float dt,
	std::vector<int>& ptr,
	std::vector<int>& col,
	std::vector<float>& val,
	std::vector<float>& rhs)const
{
	const int n2 = cNXP * cNYP * cNZP; // Number of points in the grid.

	ptr.clear();
	ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	col.reserve(n2 * 7);

	val.clear();
	val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2);

	for (int j = 0; j < cNYP; ++j)
	{
		for (int k = 0; k < cNZP; k++)
		{
			for (int i = 0; i < cNXP; i++)
			{
				const int index = i + (k * cNXP) + (j*cNXP*cNZP);

				if (i == cNXP - 2)
				{
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				//else if (
				//	// baixo
				//	(i == 0 && k == 0) || (i == cNXP - 1 && k == 0)
				//	|| (j == cNYP - 1 && k == 0) || (k == 0 && j == 0) ||
				//	// cima
				//	(i == 0 && k == cNZP - 1) || (i == cNXP - 1 && k == cNZP - 1)
				//	|| (j == cNYP - 1 && k == cNZP - 1) || (k == cNZP - 1 && j == 0) ||
				//	//laterais
				//	(i == 0 && j == 0) || (i == cNXP - 1 && j == 0)
				//	|| (i == 0 && j == cNYP - 1) || (i == cNXP - 1 && j == cNYP - 1)
				//	)
				//{
				//	col.push_back(index);
				//	val.push_back(1.0);
				//	rhs[index] = 0.0;
				//}
				else if (i == cNXP - 1)
				{
					// Outlet => Pressure set to zero
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (k == 0)
				{
					// surface => dPdn = 0
					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (k == cNZP - 1)
				{
					// top => dPdn = 0

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (i == 0)
				{
					// inlet => dPdn = 0

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					rhs[index] = 0.0;

				}
				else if (j == cNYP - 1)
				{
					// back => dPdn = 0

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (j == 0)
				{
					// front => dPdn = 0

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					// dentro do s�lido

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else  /*if (i == 1 || i == 2 ||  i == 3 ||
				i == cNXP - 2 || i == cNXP - 3 || i == cNXP - 4 ||
				j == 1 || j == 2 || j == 3 ||
				j == cNYP - 2 || j == cNYP - 3 || j == cNYP - 4 ||
				k == 1 || k == 2 || k == 3 ||
				k == cNZP - 2 || k == cNZP - 3 || k == cNZP - 4)*/
				{
					// reduz o stencil utilizado para 2

					float valijk = 0.0;

					// jkmi
					col.push_back(index + (-1 * cNXP) + (0 * cNXP * cNZP));
					if (k == 1 || (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jkpi
					col.push_back(index + (1 * cNXP) + (0 * cNXP * cNZP));
					if (k == cNZP - 2 || (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jkim
					col.push_back(index - 1 + (0 * cNXP) + (0 * cNXP * cNZP));
					if (i == 1 || (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}


					// jkip
					col.push_back(index + 1 + (0 * cNXP) + (0 * cNXP * cNZP));
					if (i == cNXP - 2 || (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}

					// jmki
					col.push_back(index + (0 * cNXP) + (-1 * cNXP * cNZP));
					if (j == 1 || (j == cjCyEnd && i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}


					// jpki
					col.push_back(index + (0 * cNXP) + (1 * cNXP * cNZP));
					if (j == cNYP - 2 || (j == cjCyInit - 1 && i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(valijk);

					rhs[index] = ch * ch*(crho / dt)*(
						dUdx[i - 1 + (k*(cNXU - 1)) + (j*(cNXU - 1)*cNZU)] +
						dVdy[i + (k*cNXV) + ((j - 1)*cNXV *cNZV)] +
						dWdz[i + ((k - 1)*cNXW) + (j*cNXW *(cNZW - 1))]
						);
				}
				/*else
				{
					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]

					// jkmi
					col.push_back(index + (-1 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jkmmi
					col.push_back(index + (-2 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0 / (12.0));


					// jkim
					col.push_back(index - 1 + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jkimm
					col.push_back(index - 2 + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0/(12.0));

					// jkip
					col.push_back(index + 1 + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jkipp
					col.push_back(index + 2 + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0 / (12.0));


					// jkpi
					col.push_back(index + (1 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jkppi
					col.push_back(index + (2 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0/(12.0));

					// jmki
					col.push_back(index + (0 * cNXP) + (-1 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jmmki
					col.push_back(index + (0 * cNXP) + (-2 * cNXP * cNZP));
					val.push_back(-1.0 / (12.0));


					// jpki
					col.push_back(index + (0 * cNXP) + (1 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jppki
					col.push_back(index + (0 * cNXP) + (2 * cNXP * cNZP));
					val.push_back(-1.0/(12.0));


					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-30.0*3.0 / (12.0));

					rhs[index] = ch * ch*(crho / dt)*(
						dUdx[i - 1 + (k*(cNXU - 1)) + (j*(cNXU - 1)*cNZU)] +
						dVdy[i + (k*cNXV) + ((j - 1)*cNXV *cNZV)] +
						dWdz[i + ((k - 1)*cNXW) + (j*cNXW *(cNZW - 1))]
						);
				}*/

				ptr.push_back(col.size());
			}
		}
	}

	return n2;
}

void SquareCylinder3DStaggeredFieldFFast::CalcdphidtUmom(field3D& dphidt,
	LSSolverF& lssolver,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc)
{
	int nn = 0;

	auto aa1 = aa;
	auto aa2 = aa;
	auto aa3 = aa;
	auto aa4 = aa;

	auto bb1 = bb;
	auto bb2 = bb;
	auto bb3 = bb;
	auto bb4 = bb;

	auto cc1 = cc;
	auto cc2 = cc;
	auto cc3 = cc;
	auto cc4 = cc;

	auto wnpxInterpolated = cdphidx;
	
	// interpolate v in x
	nn = CreateInterpolateXSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cdphidx);

	auto aa5 = aa;
	auto aa6 = aa;
	auto aa7 = aa;

	auto bb5 = bb;
	auto bb6 = bb;
	auto bb7 = bb;

	auto cc5 = cc;
	auto cc6 = cc;
	auto cc7 = cc;

	auto cdphidx5 = cdphidx;
	auto cdphidx6 = cdphidx;
	auto cdphidx7 = cdphidx;

	//field3D vnpInterpolated;

	//lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

	//vnpInterpolated = cdphidx;

	// interpolate w in x

	nn = CreateInterpolateXSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa1, bb1, cc1, wnpxInterpolated);

	//field3D wnpInterpolated;

	//lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);
	//wnpInterpolated = cdphidx;

	// interpolate u in x, y, z

	nn = CreateInterpolateXSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa2, bb2, cc2, cphiInterpolateX);

	//lssolver.SolveTridiagonalDestructive(cphiInterpolateX, nn, aa, bb, cc);

	nn = CreateInterpolateYSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa3, bb3, cc3, cphiInterpolateY);

	//lssolver.SolveTridiagonalDestructive(cphiInterpolateY, nn, aa, bb, cc);

	nn = CreateInterpolateZSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa4, bb4, cc4, cphiInterpolateZ);

	//lssolver.SolveTridiagonalDestructive(cphiInterpolateZ, nn, aa, bb, cc);

	std::vector<const std::vector<float>*> aaa, bbb;

	std::vector<std::vector<float>*> ccc, rrhs;

	aaa = { &aa, &aa1, &aa2 , &aa3 , &aa4 , &aa5 , &aa6 , &aa7 };
	bbb = { &bb, &bb1, &bb2 , &bb3 , &bb4 , &bb5 , &bb6 , &bb7 };
	ccc = { &cc, &cc1, &cc2 , &cc3 , &cc4 , &cc5 , &cc6 , &cc7 };

	rrhs = { &cdphidx, &wnpxInterpolated, &cphiInterpolateX, &cphiInterpolateY, &cphiInterpolateZ,
		&cdphidx5, &cdphidx6, &cdphidx7 };
	
	lssolver.SolveTridiagonalDestructiveAVX2(rrhs, nn, aaa, bbb, ccc);

	auto vnpInterpolated = cdphidx;
	auto wnpInterpolated = wnpxInterpolated;
	
	// u mom eq

	// calc duudx
	nn = CreatedphidxUConvSystemTD(aa, bb, cc, cdphidx);

	lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

	// calc duvdy
	nn = CreatedphidyUConvSystemTD(vnpInterpolated, aa, bb, cc, cdphidy);

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	// calc duwdz
	nn = CreatedphidzUConvSystemTD(wnpInterpolated, aa, bb, cc, cdphidz);

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	// calc d2udx2
	nn = Created2phidx2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, cd2phidx2);

	lssolver.SolveTridiagonalDestructive(cd2phidx2, nn, aa, bb, cc);

	// calc d2udy2
	nn = Created2phidy2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, cd2phidy2);

	lssolver.SolveTridiagonalDestructive(cd2phidy2, nn, aa, bb, cc);

	// calc d2udz2
	nn = Created2phidz2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, cd2phidz2);

	lssolver.SolveTridiagonalDestructive(cd2phidz2, nn, aa, bb, cc);

	// get dudt

	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				if (i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					dphidt[i + (k*cNXU) + (j* cNXU* cNZU)] = 0.0;
				}
				else
				{
					dphidt[i + (k*cNXU) + (j* cNXU* cNZU)] =
						-1.0*(cdphidx[i + (k* cNXU) + (j* cNXU* cNZU)] +
							cdphidy[i + (k* cNXU) + (j* cNXU* cNZU)] +
							cdphidz[i + (k* cNXU) + (j* cNXU* cNZU)])
						+ (cmu / crho)*(cd2phidx2[i + (k* cNXU) + (j* cNXU* cNZU)]
							+ cd2phidy2[i + (k* cNXU) + (j* cNXU* cNZU)] +
							cd2phidz2[i + (k* cNXU) + (j* cNXU* cNZU)]);
				}
			}
		}
	}
}

void SquareCylinder3DStaggeredFieldFFast::CalcdphidtVmom(field3D& dphidt,
	LSSolverF& lssolver,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc)
{
	int nn = 0;

	// interpolate u in y
	nn = CreateInterpolateYSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, cdphidy);

	field3D unpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	unpInterpolated = cdphidy;

	// interpolate w in y

	nn = CreateInterpolateYSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cdphidy);

	field3D wnpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	wnpInterpolated = cdphidy;

	// interpolate v in x, y, z

	nn = CreateInterpolateXSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cphiInterpolateX);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateX, nn, aa, bb, cc);

	nn = CreateInterpolateYSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cphiInterpolateY);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateY, nn, aa, bb, cc);

	nn = CreateInterpolateZSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cphiInterpolateZ);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateZ, nn, aa, bb, cc);

	// v mom eq

	// calc duvdx
	nn = CreatedphidxVConvSystemTD(unpInterpolated, aa, bb, cc, cdphidx);

	lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

	// calc dvvdy
	nn = CreatedphidyVConvSystemTD(aa, bb, cc, cdphidy);

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	// calc dvwdz
	nn = CreatedphidzVConvSystemTD(wnpInterpolated, aa, bb, cc, cdphidz);

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	// calc d2vdx2
	nn = Created2phidx2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cd2phidx2);

	lssolver.SolveTridiagonalDestructive(cd2phidx2, nn, aa, bb, cc);

	// calc d2vdy2
	nn = Created2phidy2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cd2phidy2);

	lssolver.SolveTridiagonalDestructive(cd2phidy2, nn, aa, bb, cc);

	// calc d2vdz2
	nn = Created2phidz2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cd2phidz2);

	lssolver.SolveTridiagonalDestructive(cd2phidz2, nn, aa, bb, cc);

	// get dvdt

	for (int j = 1; j < cNYV - 1; j++)
	{
		for (int k = 1; k < cNZV - 1; k++)
		{
			for (int i = 1; i < cNXV - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd  &&
					j >= cjCyInit - 1 && j < cjCyEnd)
				{
					dphidt[i + (k*cNXV) + (j* cNXV* cNZV)] = 0.0;
				}
				else
				{
					dphidt[i + (k*cNXV) + (j* cNXV* cNZV)] =
						-1.0*(cdphidx[i + (k* cNXV) + (j* cNXV* cNZV)] +
							cdphidy[i + (k* cNXV) + (j* cNXV* cNZV)] +
							cdphidz[i + (k* cNXV) + (j* cNXV* cNZV)])
						+ (cmu / crho)*(cd2phidx2[i + (k* cNXV) + (j* cNXV* cNZV)]
							+ cd2phidy2[i + (k* cNXV) + (j* cNXV* cNZV)] +
							cd2phidz2[i + (k* cNXV) + (j* cNXV* cNZV)]);
				}
			}
		}
	}
}

void SquareCylinder3DStaggeredFieldFFast::CalcdphidtWmom(field3D& dphidt,
	LSSolverF& lssolver,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc)
{
	int nn = 0;

	// interpolate u in z
	nn = CreateInterpolateZSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, cdphidz);

	field3D unpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	unpInterpolated = cdphidz;

	// interpolate v in z

	nn = CreateInterpolateZSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cdphidz);

	field3D vnpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	vnpInterpolated = cdphidz;

	// interpolate w in x, y, z

	nn = CreateInterpolateXSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cphiInterpolateX);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateX, nn, aa, bb, cc);

	nn = CreateInterpolateYSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cphiInterpolateY);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateY, nn, aa, bb, cc);

	nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cphiInterpolateZ);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateZ, nn, aa, bb, cc);

	// w mom eq

	// calc duwdx
	nn = CreatedphidxWConvSystemTD(unpInterpolated, aa, bb, cc, cdphidx);

	lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

	// calc dvwdy
	nn = CreatedphidyWConvSystemTD(vnpInterpolated, aa, bb, cc, cdphidy);

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	// calc dwwdz
	nn = CreatedphidzWConvSystemTD(aa, bb, cc, cdphidz);

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	// calc d2udx2
	nn = Created2phidx2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cd2phidx2);

	lssolver.SolveTridiagonalDestructive(cd2phidx2, nn, aa, bb, cc);

	// calc d2udy2
	nn = Created2phidy2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cd2phidy2);

	lssolver.SolveTridiagonalDestructive(cd2phidy2, nn, aa, bb, cc);

	// calc d2udz2
	nn = Created2phidz2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cd2phidz2);

	lssolver.SolveTridiagonalDestructive(cd2phidz2, nn, aa, bb, cc);

	// get dudt

	for (int j = 0; j < cNYW; j++)
	{
		for (int k = 0; k < cNZW; k++)
		{
			for (int i = 0; i < cNXW; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd &&
					j >= cjCyInit && j < cjCyEnd)
				{
					dphidt[i + (k*cNXW) + (j* cNXW* cNZW)] = 0.0;
				}
				else
				{
					dphidt[i + (k*cNXW) + (j* cNXW* cNZW)] =
						-1.0*(cdphidx[i + (k* cNXW) + (j* cNXW* cNZW)] +
							cdphidy[i + (k* cNXW) + (j* cNXW* cNZW)] +
							cdphidz[i + (k* cNXW) + (j* cNXW* cNZW)])
						+ (cmu / crho)*(cd2phidx2[i + (k* cNXW) + (j* cNXW* cNZW)]
							+ cd2phidy2[i + (k* cNXW) + (j* cNXW* cNZW)] +
							cd2phidz2[i + (k* cNXW) + (j* cNXW* cNZW)]);
				}
			}
		}
	}
}


void SquareCylinder3DStaggeredFieldFFast::SaveIKCutToCSV(const std::string& filename, int j)
{
	std::ofstream arqU(filename + "U.csv");

	for (int k = 0; k < cNZU; k++)
	{
		for (int i = 0; i < cNXU; i++)
		{
			arqU << cUnp[i + (k*cNXU) + (j * cNXU*cNZU)] << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqV(filename + "V.csv");

	for (int k = 0; k < cNZV; k++)
	{
		for (int i = 0; i < cNXV; i++)
		{
			arqV << cVnp[i + (k*cNXV) + (j * cNXV*cNZV)] << "\t";
		}

		arqV << "\n";
	}

	arqV.close();

	std::ofstream arqW(filename + "W.csv");

	for (int k = 0; k < cNZW; k++)
	{
		for (int i = 0; i < cNXW; i++)
		{
			arqW << cWnp[i + (k*cNXW) + (j * cNXW*cNZW)] << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "P.csv");

	for (int k = 0; k < cNZP; k++)
	{
		for (int i = 0; i < cNXP; i++)
		{
			arqP << cPn[i + (k*cNXP) + (j * cNXP*cNZP)] << "\t";
		}

		arqP << "\n";
	}

	arqP.close();

	/*std::ofstream arqDivU(filename + "divU.csv");

	for (int k = 0; k < cNZU; k++)
	{
		for (int i = 0; i < (cNXU-1); i++)
		{
			arqDivU << cdphidx[i + (k*(cNXU - 1)) + (j * (cNXU - 1)*cNZU)]  << "\t";
		}

		arqDivU << "\n";
	}

	arqDivU.close();

	std::ofstream arqDivV(filename + "divV.csv");

	for (int k = 0; k < cNZV; k++)
	{
		for (int i = 0; i < cNXV; i++)
		{
			arqDivV << cdphidy[i + (k* cNXV) + ((j) * cNXV *cNZV)] << "\t";
		}

		arqDivV << "\n";
	}

	arqDivV.close();

	std::ofstream arqDivW(filename + "divW.csv");

	for (int k = 0; k < cNZW-1; k++)
	{
		for (int i = 0; i < cNXW; i++)
		{
			arqDivW << cdphidz[i + (k* cNXW) + (j * cNXW *(cNZW-1))] << "\t";
		}

		arqDivW << "\n";
	}

	arqDivW.close();*/

}

void SquareCylinder3DStaggeredFieldFFast::SaveJKCutToCSV(const std::string& filename, const int i)
{
	std::ofstream arqU(filename + "TransU.csv");

	for (int j = 0; j < cNYU; j++)
	{
		for (int k = 0; k < cNZU; k++)
		{
			arqU << cUnp[i + (k*cNXU) + (j * cNXU*cNZU)] << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqV(filename + "TransV.csv");

	for (int j = 0; j < cNYV; j++)
	{
		for (int k = 0; k < cNZV; k++)
		{
			arqV << cVnp[i + (k*cNXV) + (j * cNXV*cNZV)] << "\t";
		}

		arqV << "\n";
	}

	arqV.close();

	std::ofstream arqW(filename + "TransW.csv");

	for (int j = 0; j < cNYW; j++)
	{
		for (int k = 0; k < cNZW; k++)
		{
			arqW << cWnp[i + (k*cNXW) + (j * cNXW*cNZW)] << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "TransP.csv");

	for (int j = 0; j < cNYP; j++)
	{
		for (int k = 0; k < cNZP; k++)
		{
			arqP << cPn[i + (k*cNXP) + (j * cNXP*cNZP)] << "\t";
		}

		arqP << "\n";
	}

	arqP.close();

	//std::ofstream arqDiv(filename + "Transdiv.csv");

	//for (int j = 0; j < cNYP; j++)
	//{
	//	for (int k = 0; k < cNZP; k++)
	//	{
	//		/*arqDiv << cdphidx[i + (k*cNXP) + (j * cNXP*cNZP)] + cdphidy[i + (k*cNXP) + (j * cNXP*cNZP)] +
	//			cdphidz[i + (k*cNXP) + (j * cNXP*cNZP)] << "\t";*/
	//	}

	//	arqDiv << "\n";
	//}

	//arqDiv.close();
}

void SquareCylinder3DStaggeredFieldFFast::SaveIJCutToCSV(const std::string& filename, int k)
{
	std::ofstream arqU(filename + "TransKU.csv");

	for (int j = 0; j < cNYU; j++)
	{
		for (int i = 0; i < cNXU; i++)
		{
			arqU << cUnp[i + (k*cNXU) + (j * cNXU*cNZU)] << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqV(filename + "TransKV.csv");

	for (int j = 0; j < cNYV; j++)
	{
		for (int i = 0; i < cNXV; i++)
		{
			arqV << cVnp[i + (k*cNXV) + (j * cNXV*cNZV)] << "\t";
		}

		arqV << "\n";
	}

	arqV.close();

	std::ofstream arqW(filename + "TransKW.csv");

	for (int j = 0; j < cNYW; j++)
	{
		for (int i = 0; i < cNXW; i++)
		{
			arqW << cWnp[i + (k*cNXW) + (j * cNXW*cNZW)] << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "TransKP.csv");

	for (int j = 0; j < cNYP; j++)
	{
		for (int i = 0; i < cNXP; i++)
		{
			arqP << cPn[i + (k*cNXP) + (j * cNXP*cNZP)] << "\t";
		}

		arqP << "\n";
	}

	arqP.close();
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidxUConvSystemTD(
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXU * cNYU * cNZU;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYU; j++)
		for (int k = 0; k < cNZU; k++)
			for (int i = 0; i < cNXU; i++)
			{
				const int index = i + (k * cNXU) + (j * cNXU * cNZU);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == cNXU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((i == ciCyInit - 1 || i == ciCyInit - 2 || i == ciCyEnd - 1 || i == ciCyEnd)
					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uipjk = cphiInterpolateX[i + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					const float uimjk = cphiInterpolateX[i - 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					rhs[index] = ((uipjk * uipjk)
						- (uimjk*uimjk)
						) /
						(ch);
				}
				else if (i == 1 || i == 2 || i == cNXU - 2 || i == cNXU - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uipjk = cphiInterpolateX[i + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					const float uimjk = cphiInterpolateX[i - 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					rhs[index] = ((uipjk * uipjk)
						- (uimjk*uimjk)
						) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uipjk = 0.5*(
						cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] + cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)]
						);

					const float uimjk = 0.5*(
						cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] + cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)]
						);

					rhs[index] = ((uipjk * uipjk)
						- (uimjk*uimjk)
						) /
						(ch);
				}*/
				/*else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
					const float uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const float uippjk = cUnp[i + 2 + (k * cNXU) + (j * cNXU * cNZU)];

					rhs[index] = 0.5 / ch * (-3 * uijk * uijk +
						4 * uipjk * uipjk -
						 uippjk* uippjk);
				}
				else if (i == cNXU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
					const float uimjk = cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const float uimmjk = cUnp[i - 2 + (k * cNXU) + (j * cNXU * cNZU)];

					rhs[index] = 0.5 / ch * (3 * uijk*uijk -
						4 * uimjk*uimjk +
						uimmjk*uimmjk);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uimjk = cUnp[i -1 + (k * cNXU) + (j * cNXU * cNZU)];
					const float uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
					const float uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const float uippjk = cUnp[i + 2 + (k * cNXU) + (j * cNXU * cNZU)];
					const float uipppjk = cUnp[i + 3 + (k * cNXU) + (j * cNXU * cNZU)];

					rhs[index] = (1.0 / ch) * ((-2 * 0.025 / 3 - 1.0 / 3.0)*uimjk*uimjk -
						(-8 * 0.025 / 3 + 0.5)*uijk*uijk
						+ (-4 * 0.025 + 1.0)*uipjk*uipjk -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*uippjk +
						(-2 * 0.025 / 3.0*uipppjk*uipppjk));
				}
				else if (i == cNXU - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const float uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
					const float uimjk = cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const float uimmjk = cUnp[i - 2 + (k * cNXU) + (j * cNXU * cNZU)];
					const float uimmmjk = cUnp[i - 3 + (k * cNXU) + (j * cNXU * cNZU)];

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*uipjk -
						(8 * 0.09 / 3 + 0.5)*uijk
						+ (4 * 0.09 + 1.0)*uimjk -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*uimmjk +
						(2 * 0.09 / 3.0*uimmmjk));
				}*/
				else
				{

					const float uipjk = cphiInterpolateX[i + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
					const float uippjk = cphiInterpolateX[i + 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					const float uimjk = cphiInterpolateX[i - 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
					const float uimmjk = cphiInterpolateX[i - 2 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (
						-Fimp / 6.0 * uimmjk * uimmjk
						- Eimp / 1.0*uimjk*uimjk
						+ Eimp / 1.0*uipjk*uipjk
						+ Fimp / 6.0* uippjk * uippjk);*/

					rhs[index] = bII / (3 * ch)*(uippjk*uippjk - uimmjk)
						+ aII / ch * (uipjk*uipjk - uimjk * uimjk);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidyUConvSystemTD(
	const field3D& V,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXU * cNYU * cNZU;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYU; j++)
		for (int k = 0; k < cNZU; k++)
			for (int i = 0; i < cNXU; i++)
			{
				const int index = i + (k * cNXU) + (j * cNXU * cNZU);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == cNYU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((j == cjCyInit || j == cjCyInit - 1 || j == cjCyEnd - 1 || j == cjCyEnd)
					&& k >= ckCyInit && k < ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijmk = cphiInterpolateY[i + (k *cNXU) + ((j - 1) * cNXU * cNZU)];
					const float uijpk = cphiInterpolateY[i + (k*cNXU) + ((j)* cNXU * cNZU)];

					const float vijpk = V[i + (k * (cNXV - 1)) + ((j) * (cNXV - 1) * cNZV)];
					const float vijmk = V[i + (k * (cNXV - 1)) + ((j - 1) * (cNXV - 1) * cNZV)];

					rhs[index] = ((uijpk * vijpk)
						- (uijmk * vijmk)) / ch;
				}
				else if (j == 1 || j == 2 || j == cNYU - 2 || j == cNYU - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijmk = cphiInterpolateY[i + (k *cNXU) + ((j - 1) * cNXU * cNZU)];
					const float uijpk = cphiInterpolateY[i + (k*cNXU) + ((j)* cNXU * cNZU)];

					const float vijpk = V[i + (k * (cNXV - 1)) + ((j) * (cNXV - 1) * cNZV)];
					const float vijmk = V[i + (k * (cNXV - 1)) + ((j - 1) * (cNXV - 1) * cNZV)];

					rhs[index] = ((uijpk * vijpk)
						- (uijmk * vijmk)) / ch;
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijpk = 0.5*(
						cUnp[i + (k * cNXU) + ((j) * cNXU * cNZU)] + cUnp[i + (k * cNXU) + ((j + 1) * cNXU * cNZU)]
						);

					const float uijmk = 0.5*(
						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)]
						);

					rhs[index] = ((uijpk * V[i + (k * (cNXV-1)) + ((j) * (cNXV - 1) * cNZV)])
						- (uijmk * V[i  + (k * (cNXV - 1)) + ((j-1) * (cNXV - 1) * cNZV)])) /ch;
				}
				else if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijk = cUnp[i + (k*cNXU) + (j*cNXU*cNZU)];
					const float uijpk = cUnp[i + (k*cNXU) + ((j+1)*cNXU*cNZU)];
					const float uijppk = cUnp[i + (k*cNXU) + ((j + 2)*cNXU*cNZU)];

					const float vijk = V[i + (k*(cNXV - 1)) + (j * (cNXV - 1) * cNXV)];
					const float vijpk = V[i + (k*(cNXV - 1)) + ((j+1) * (cNXV - 1) * cNXV)];
					const float vijppk = V[i + (k*(cNXV - 1)) + ((j+2) * (cNXV - 1) * cNXV)];

					rhs[index] = 0.5 / ch * (-3 * uijk * V[i + (k * cNXU) + (j * cNXU * cNZV)] +
						4 * uijpk * V[i + (k * cNXU) + ((j + 1) * cNXU * cNZV)] -
						uijppk * V[i + (k * cNXU) + ((j + 2) * cNXU * cNZV)]);
				}
				else if (j == cNYU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + ((j-1) * cNXU * cNZV)] -
						4 * cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 2) * cNXU * cNZV)] +
						cUnp[i + (k * cNXU) + ((j - 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 3) * cNXU * cNZV)]);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 1) * cNXU * cNZV)] -
						(-8 * 0.025 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + (j * cNXU * cNZV)]
						+ (-4 * 0.025 + 1.0)*cUnp[i + (k * cNXU) + ((j + 1)* cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 1)* cNXU * cNZV)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cUnp[i + (k * cNXU) + ((j + 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 2) * cNXU * cNZV)] +
						(-2 * 0.025 / 3.0*cUnp[i + (k * cNXU) + ((j + 3) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 3) * cNXU * cNZV)]));
				}
				else if (j == cNYU - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cUnp[i + (k * cNXU) + ((j + 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j) * cNXU * cNZV)] -
						(8 * 0.09 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + ((j-1) * cNXU * cNZV)]
						+ (4 * 0.09 + 1.0)*cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 2) * cNXU * cNZV)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cUnp[i + (k * cNXU) + ((j - 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 3) * cNXU * cNZV)] +
						(2 * 0.09 / 3.0*cUnp[i + (k * cNXU) + ((j - 3) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 4) * cNXU * cNZV)]));
				}*/
				else
				{

					const float uijmk = cphiInterpolateY[i + (k *cNXU) + ((j - 1) * cNXU * cNZU)];
					const float uijmmk = cphiInterpolateY[i + (k *cNXU) + ((j - 2) * cNXU * cNZU)];

					const float uijpk = cphiInterpolateY[i + (k*cNXU) + ((j)* cNXU * cNZU)];
					const float uijppk = cphiInterpolateY[i + (k*cNXU) + ((j + 1)* cNXU * cNZU)];

					const float vijpk = V[i + (k * (cNXV - 1)) + ((j) * (cNXV - 1) * cNZV)];
					const float vijppk = V[i + (k * (cNXV - 1)) + ((j + 1) * (cNXV - 1) * cNZV)];

					const float vijmk = V[i + (k * (cNXV - 1)) + ((j - 1) * (cNXV - 1) * cNZV)];
					const float vijmmk = V[i + (k * (cNXV - 1)) + ((j - 2) * (cNXV - 1) * cNZV)];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (-Fimp / 6.0*uijmmk * vijmmk -
						Eimp / 1.0*uijmk * vijmk +
						Eimp / 1.0*uijpk * vijpk +
						Fimp / 6.0*uijppk * vijppk);*/

					rhs[index] = bII / (3 * ch)*(uijppk*vijppk - uijmmk * vijmmk)
						+ aII / ch * (uijpk*vijpk - uijmk * vijmk);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidzUConvSystemTD(
	const field3D& W,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXU * cNYU * cNZU;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYU; j++)
		for (int k = 0; k < cNZU; k++)
			for (int i = 0; i < cNXU; i++)
			{
				const int index = i + (k * cNXU) + (j * cNXU * cNZU);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == cNZU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((k == ckCyInit || k == ckCyInit - 1 || k == ckCyEnd - 1 || k == ckCyEnd)
					&& i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijkp = cphiInterpolateZ[i + (k*cNXU) + (j*cNXU*(cNZU - 1))];
					const float uijkm = cphiInterpolateZ[i + ((k - 1)*cNXU) + (j*cNXU*(cNZU - 1))];

					const float wijkp = W[i + ((k) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
					const float wijkm = W[i + ((k - 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];

					rhs[index] = ((uijkp * wijkp)
						- (uijkm * wijkm)) /
						(ch);
				}
				else if (k == 1 || k == 2 || k == cNZU - 2 || k == cNZU - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijkp = cphiInterpolateZ[i + (k*cNXU) + (j*cNXU*(cNZU - 1))];
					const float uijkm = cphiInterpolateZ[i + ((k - 1)*cNXU) + (j*cNXU*(cNZU - 1))];

					const float wijkp = W[i + ((k) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
					const float wijkm = W[i + ((k - 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];

					rhs[index] = ((uijkp * wijkp)
						- (uijkm * wijkm)) /
						(ch);
				}
				/*
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float uijkp = 0.5*(
						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + ((k+1) * cNXU) + ((j) * cNXU * cNZU)]
						);

					const float uijkm = 0.5*(
						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + ((k-1) * cNXU) + ((j) * cNXU * cNZU)]
						);

					rhs[index] = ((uijkp * W[i + ((k) * (cNXW-1)) + ((j)* (cNXW - 1) * cNZW)])
						- (uijkm * W[i + ((k-1) * (cNXW - 1)) + (j  * (cNXW - 1) * cNZW)])) /
						(ch);
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + (k * cNXU) + (j * cNXU * cNZW)] +
						4 * cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZW)] -
						cUnp[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZW)]);
				}
				else if (k == cNZU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + ((k-1) * cNXU) + (j * cNXU * cNZW)] -
						4 * cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZW)] +
						cUnp[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZW)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + (k * cNXU) + (j * cNXU * cNZW)]
						+ (-4 * 0.025 + 1.0)*cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cUnp[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZW)] +
						(-2 * 0.025 / 3.0*cUnp[i + ((k + 3) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 3) * cNXU) + ((j)* cNXU * cNZW)]));
				}
				else if (k == cNZU - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k) * cNXU) + ((j)* cNXU * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + ((k-1) * cNXU) + (j * cNXU * cNZW)]
						+ (4 * 0.09 + 1.0)*cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cUnp[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZW)] +
						(2 * 0.09 / 3.0*cUnp[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 4) * cNXU) + ((j)* cNXU * cNZW)]));
				}*/
				else
				{
					const float uijkp = cphiInterpolateZ[i + (k*cNXU) + (j*cNXU*(cNZU - 1))];
					const float uijkpp = cphiInterpolateZ[i + ((k + 1)*cNXU) + (j*cNXU*(cNZU - 1))];

					const float uijkm = cphiInterpolateZ[i + ((k - 1)*cNXU) + (j*cNXU*(cNZU - 1))];
					const float uijkmm = cphiInterpolateZ[i + ((k - 2)*cNXU) + (j*cNXU*(cNZU - 1))];

					const float wijkp = W[i + ((k) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
					const float wijkpp = W[i + ((k + 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];

					const float wijkm = W[i + ((k - 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
					const float wijkmm = W[i + ((k - 2) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (
						-Fimp / 6.0* uijkmm * wijkmm
						- Eimp / 1.0* uijkm * wijkm
						+ Eimp / 1.0*uijkp * wijkp
						+ Fimp / 6.0*uijkpp * wijkpp);*/

					rhs[index] = bII / (3 * ch)*(uijkpp*wijkpp - uijkmm * wijkmm)
						+ aII / ch * (uijkp*wijkp - uijkm * wijkm);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidxVConvSystemTD(
	const field3D& U,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXV * cNYV * cNZV;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYV; j++)
		for (int k = 0; k < cNZV; k++)
			for (int i = 0; i < cNXV; i++)
			{
				const int index = i + (k * cNXV) + (j * cNXV * cNZV);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == cNXV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((i == ciCyInit || i == ciCyInit - 1 || i == ciCyEnd - 1 || i == ciCyEnd)
					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vipjk = cphiInterpolateX[i + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
					const float vimjk = cphiInterpolateX[i - 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];

					const float uipjk = U[i + (k *cNXU) + (j*cNXU*cNZU)];
					const float uimjk = U[i - 1 + (k *cNXU) + (j*cNXU*cNZU)];

					rhs[index] = ((vipjk * uipjk)
						- (vimjk * uimjk)) /
						(ch);
				}
				else if (i == 1 || i == 2 || i == cNXV - 2 || i == cNXV - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vipjk = cphiInterpolateX[i + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
					const float vimjk = cphiInterpolateX[i - 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];

					const float uipjk = U[i + (k *cNXU) + (j*cNXU*cNZU)];
					const float uimjk = U[i - 1 + (k *cNXU) + (j*cNXU*cNZU)];

					rhs[index] = ((vipjk * uipjk)
						- (vimjk * uimjk)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vipjk = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + 1 + ((k) * cNXV) + ((j)* cNXV * cNZV)]
						);

					const float vimjk = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i - 1 + ((k) * cNXV) + ((j)* cNXV * cNZV)]
						);

					rhs[index] = ((vipjk * U[i + (k * cNXU) + ((j)* cNXU  * cNZU)])
						- (vimjk * U[i  - 1 + (k * cNXU) + (j  * cNXU * cNZU)])) /
						(ch);
				}
				else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)] +
						4 * cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZU)] -
						cVnp[i + 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZU)]);
				}
				else if (i == cNXV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1+ (k * cNXU) + (j * cNXU * cNZU)] -
						4 * cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZU)] +
						cVnp[i - 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZU)]);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZU)] -
						(-8 * 0.025 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)]
						+ (-4 * 0.025 + 1.0)*cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZU)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cVnp[i + 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZU)] +
						(-2 * 0.025 / 3.0*cVnp[i + 3 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 3 + (k * cNXU) + (j * cNXU * cNZU)]));
				}
				else if (i == cNXV - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)] -
						(8 * 0.09 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZU)]
						+ (4 * 0.09 + 1.0)*cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZU)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cVnp[i - 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZU)] +
						(2 * 0.09 / 3.0*cVnp[i - 3 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 4 + (k * cNXU) + (j * cNXU * cNZU)]));
				}*/
				else
				{
					const float vipjk = cphiInterpolateX[i + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
					const float vippjk = cphiInterpolateX[i + 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];

					const float vimjk = cphiInterpolateX[i - 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
					const float vimmjk = cphiInterpolateX[i - 2 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];

					const float uipjk = U[i + (k *cNXU) + (j*cNXU*cNZU)];
					const float uippjk = U[i + 1 + (k *cNXU) + (j*cNXU*cNZU)];

					const float uimjk = U[i - 1 + (k *cNXU) + (j*cNXU*cNZU)];
					const float uimmjk = U[i - 2 + (k *cNXU) + (j*cNXU*cNZU)];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (-Fimp / 6.0*vimmjk * uimmjk -
						Eimp / 1.0*vimjk * uimjk +
						Eimp / 1.0*vipjk * uipjk +
						Fimp / 6.0*vippjk * uippjk);*/

					rhs[index] = bII / (3 * ch)*(vippjk*uippjk - vimmjk * uimmjk)
						+ aII / ch * (vipjk*uipjk - vimjk * uimjk);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidyVConvSystemTD(
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXV * cNYV * cNZV;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYV; j++)
		for (int k = 0; k < cNZV; k++)
			for (int i = 0; i < cNXV; i++)
			{
				const int index = i + (k * cNXV) + (j * cNXV * cNZV);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == cNYV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((j == cjCyInit - 1 || j == cjCyInit - 2 || j == cjCyEnd - 2 || j == cjCyEnd)
					&& k >= ckCyInit && k < ckCyEnd && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vijpk = cphiInterpolateY[i + (k * cNXV) + (j*cNXV * cNZV)];

					const float vijmk = cphiInterpolateY[i + (k * cNXV) + ((j - 1)*cNXV * cNZV)];

					rhs[index] = ((vijpk*vijpk)
						- (vijmk*vijmk)) /
						(ch);
				}
				else if (j == 1 || j == 2 || j == cNYV - 2 || j == cNYV - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vijpk = cphiInterpolateY[i + (k * cNXV) + (j*cNXV * cNZV)];

					const float vijmk = cphiInterpolateY[i + (k * cNXV) + ((j - 1)*cNXV * cNZV)];

					rhs[index] = ((vijpk*vijpk)
						- (vijmk*vijmk)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vijpk = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k)* cNXV) + ((j+1)* cNXV * cNZV)]
						);

					const float vijmk = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k) * cNXV) + ((j-1)* cNXV * cNZV)]
						);

					rhs[index] = ((vijpk*vijpk)
						- (vijmk*vijmk)) /
						(ch);
				}
				else if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] +
						4 * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] -
						cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)]);
				}
				else if (j == cNYV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] -
						4 * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] +
						cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)]);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] -
						(-8 * 0.025 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)]
						+ (-4 * 0.025 + 1.0)*cVnp[i + (k * cNXV) + ((j + 1)* cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1)* cNXV * cNZV)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] +
						(-2 * 0.025 / 3.0*cVnp[i + (k * cNXV) + ((j + 3) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 3) * cNXV * cNZV)]));
				}
				else if (j == cNYV - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] -
						(8 * 0.09 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)]
						+ (4 * 0.09 + 1.0)*cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] +
						(2 * 0.09 / 3.0*cVnp[i + (k * cNXV) + ((j - 3) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 3) * cNXV * cNZV)]));
				}*/
				else
				{
					const float vijpk = cphiInterpolateY[i + (k * cNXV) + (j*cNXV * cNZV)];
					const float vijppk = cphiInterpolateY[i + (k * cNXV) + ((j + 1)*cNXV * cNZV)];

					const float vijmk = cphiInterpolateY[i + (k * cNXV) + ((j - 1)*cNXV * cNZV)];
					const float vijmmk = cphiInterpolateY[i + (k * cNXV) + ((j - 2)*cNXV * cNZV)];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (-Fimp / 6.0*vijmmk * vijmmk -
						Eimp / 1.0*vijmk * vijmk +
						Eimp / 1.0*vijpk * vijpk +
						Fimp / 6.0*vijppk * vijppk);*/

					rhs[index] = bII / (3 * ch)*(vijppk*vijppk - vijmmk * vijmmk)
						+ aII / ch * (vijpk*vijpk - vijmk * vijmk);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidzVConvSystemTD(
	const field3D& W,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXV * cNYV * cNZV;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYV; j++)
		for (int k = 0; k < cNZV; k++)
			for (int i = 0; i < cNXV; i++)
			{
				const int index = i + (k * cNXV) + (j * cNXV * cNZV);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == cNZV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((k == ckCyInit || k == ckCyInit - 1 || k == ckCyEnd - 1 || k == ckCyEnd)
					&& i >= ciCyInit && i <= ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vijkp = cphiInterpolateZ[i + ((k)* cNXV) + (j* cNXV * (cNZV - 1))];

					const float vijkm = cphiInterpolateZ[i + ((k - 1) * cNXV) + (j* cNXV * (cNZV - 1))];

					const float wijkp = W[i + ((k)* cNXW) + (j*cNXW*cNZW)];
					const float wijkm = W[i + ((k - 1)* cNXW) + (j*cNXW*cNZW)];

					rhs[index] = ((vijkp * wijkp)
						- (vijkm * wijkm)) /
						(ch);
				}
				else if (k == 1 || k == 2 || k == cNZV - 2 || k == cNZV - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vijkp = cphiInterpolateZ[i + ((k)* cNXV) + (j* cNXV * (cNZV - 1))];

					const float vijkm = cphiInterpolateZ[i + ((k - 1) * cNXV) + (j* cNXV * (cNZV - 1))];

					const float wijkp = W[i + ((k)* cNXW) + (j*cNXW*cNZW)];
					const float wijkm = W[i + ((k - 1)* cNXW) + (j*cNXW*cNZW)];

					rhs[index] = ((vijkp * wijkp)
						- (vijkm * wijkm)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float vijkp = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k+1)* cNXV) + ((j)* cNXV * cNZV)]
						);

					const float vijkm = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)]
						);

					rhs[index] = ((vijkp * W[i + ((k) * cNXW) + ((j)* cNXW * cNZW)])
						- (vijkm * W[i + ((k - 1) * cNXW) + ((j)  * cNXW * cNZW)])) /
						(ch);
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + (k * cNXW) + (j * cNXW * cNZW)] +
						4 * cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						cVnp[i + ((k + 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == cNZV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + ((k-1) * cNXW) + (j * cNXW * cNZW)] -
						4 * cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] +
						cVnp[i + ((k - 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + (k * cNXW) + (j * cNXW * cNZW)]
						+ (-4 * 0.025 + 1.0)*cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cVnp[i + ((k + 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] +
						(-2 * 0.025 / 3.0*cVnp[i + ((k + 3) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)]));
				}
				else if (k == cNZV - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + ((k-1) * cNXW) + (j * cNXW * cNZW)]
						+ (4 * 0.09 + 1.0)*cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k-2) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cVnp[i + ((k - 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)] +
						(2 * 0.09 / 3.0*cVnp[i + ((k - 3) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 4) * cNXW) + ((j)* cNXW * cNZW)]));
				}*/
				else
				{
					const float vijkp = cphiInterpolateZ[i + ((k)* cNXV) + (j* cNXV * (cNZV - 1))];
					const float vijkpp = cphiInterpolateZ[i + ((k + 1)* cNXV) + (j* cNXV * (cNZV - 1))];

					const float vijkm = cphiInterpolateZ[i + ((k - 1) * cNXV) + (j* cNXV * (cNZV - 1))];
					const float vijkmm = cphiInterpolateZ[i + ((k - 2) * cNXV) + (j* cNXV * (cNZV - 1))];

					const float wijkp = W[i + ((k)* cNXW) + (j*cNXW*cNZW)];
					const float wijkpp = W[i + ((k + 1)* cNXW) + (j*cNXW*cNZW)];

					const float wijkm = W[i + ((k - 1)* cNXW) + (j*cNXW*cNZW)];
					const float wijkmm = W[i + ((k - 2)* cNXW) + (j*cNXW*cNZW)];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (-Fimp / 6.0*vijkmm * wijkmm
						- Eimp / 1.0* vijkm * wijkm +
						Eimp / 1.0*vijkp * wijkp +
						Fimp / 6.0*vijkpp * wijkpp);*/

					rhs[index] = bII / (3 * ch)*(vijkpp*wijkpp - vijkmm * wijkmm)
						+ aII / ch * (vijkp*wijkp - vijkm * wijkm);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidxWConvSystemTD(
	const field3D& U,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXW * cNYW * cNZW;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYW; j++)
		for (int k = 0; k < cNZW; k++)
			for (int i = 0; i < cNXW; i++)
			{
				const int index = i + (k * cNXW) + (j * cNXW * cNZW);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == cNXW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((i == ciCyInit || i == ciCyInit - 1 || i == ciCyEnd - 1 || i == ciCyEnd)
					&& k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wipjk = cphiInterpolateX[i + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const float wimjk = cphiInterpolateX[i - 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const float uipjk = U[i + (k*cNXU) + (j * cNXU * (cNZU - 1))];
					const float uimjk = U[i - 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];

					rhs[index] = ((wipjk * uipjk)
						- (wimjk * uimjk)) /
						(ch);
				}
				else if (i == 1 || i == 2 || i == cNXW - 2 || i == cNXW - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wipjk = cphiInterpolateX[i + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const float wimjk = cphiInterpolateX[i - 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const float uipjk = U[i + (k*cNXU) + (j * cNXU * (cNZU - 1))];
					const float uimjk = U[i - 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];

					rhs[index] = ((wipjk * uipjk)
						- (wimjk * uimjk)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wipjk = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i+1 + ((k)* cNXW) + ((j)* cNXW * cNZW)]
						);

					const float wimjk = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i-1 + ((k) * cNXW) + ((j)* cNXW * cNZW)]
						);

					rhs[index] = ((wipjk * U[i + ((k)* cNXU) + ((j)* cNXU  * (cNZU-1))])
						- (wimjk * U[i-1 + ((k) * cNXU) + ((j)* cNXU * (cNZU - 1))])) /
						(ch);
				}
				else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)] +
						4 * cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						cWnp[i + 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZW)]);
				}
				else if (i == cNXW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						4 * cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZW)] +
						cWnp[i - 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZW)]);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)]
						+ (-4 * 0.025 + 1.0)*cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cWnp[i + 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZW)] +
						(-2 * 0.025 / 3.0*cWnp[i + 3 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 3 + (k * cNXU) + (j * cNXU * cNZW)]));
				}
				else if (i == cNXW - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i-1 + (k * cNXU) + (j * cNXU * cNZW)]
						+ (4 * 0.09 + 1.0)*cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cWnp[i - 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i -3 + (k * cNXU) + (j * cNXU * cNZW)] +
						(2 * 0.09 / 3.0*cWnp[i - 3 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 4 + (k * cNXU) + (j * cNXU * cNZW)]));
				}*/
				else
				{
					const float wipjk = cphiInterpolateX[i + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
					const float wippjk = cphiInterpolateX[i + 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const float wimjk = cphiInterpolateX[i - 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
					const float wimmjk = cphiInterpolateX[i - 2 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const float uipjk = U[i + (k*cNXU) + (j * cNXU * (cNZU - 1))];
					const float uippjk = U[i + 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];

					const float uimjk = U[i - 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];
					const float uimmjk = U[i - 2 + (k*cNXU) + (j * cNXU * (cNZU - 1))];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (-Fimp / 6.0*wimmjk * uimmjk -
						Eimp / 1.0*wimjk * uimjk +
						Eimp / 1.0*wipjk * uipjk +
						Fimp / 6.0*wippjk * uippjk);*/

					rhs[index] = bII / (3 * ch)*(wippjk*uippjk - wimmjk * uimmjk)
						+ aII / ch * (wipjk*uipjk - wimjk * uimjk);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidyWConvSystemTD(
	const field3D& V,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXW * cNYW * cNZW;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYW; j++)
		for (int k = 0; k < cNZW; k++)
			for (int i = 0; i < cNXW; i++)
			{
				const int index = i + (k * cNXW) + (j * cNXW * cNZW);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == cNYW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((j == cjCyInit || j == cjCyInit - 1 || j == cjCyEnd - 1 || j == cjCyEnd)
					&& k >= ckCyInit - 1 && k < ckCyEnd && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wijpk = cphiInterpolateY[i + (k * cNXW) + (j * cNXW * cNZW)];

					const float wijmk = cphiInterpolateY[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)];

					const float vijpk = V[i + (k * cNXV) + (j * cNXV * (cNZV - 1))];
					const float vijmk = V[i + (k * cNXV) + ((j - 1) * cNXV * (cNZV - 1))];

					rhs[index] = ((wijpk * vijpk)
						- (wijmk * vijmk)) /
						(ch);
				}
				else if (j == 1 || j == 2 || j == cNYW - 2 || j == cNYW - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wijpk = cphiInterpolateY[i + (k * cNXW) + (j * cNXW * cNZW)];

					const float wijmk = cphiInterpolateY[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)];

					const float vijpk = V[i + (k * cNXV) + (j * cNXV * (cNZV - 1))];
					const float vijmk = V[i + (k * cNXV) + ((j - 1) * cNXV * (cNZV - 1))];

					rhs[index] = ((wijpk * vijpk)
						- (wijmk * vijmk)) /
						(ch);
				}
				/*
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wijpk = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k)* cNXW) + ((j+1)* cNXW * cNZW)]
						);

					const float wijmk = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k)* cNXW) + ((j-1)* cNXW * cNZW)]
						);

					rhs[index] = ((wijpk * V[i + ((k)* cNXV) + ((j)* cNXV * (cNZV - 1))])
						- (wijmk * V[i + ((k)* cNXV) + ((j-1)* cNXV * (cNZV-1))])) /
						(ch);
				}
				else if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + (j * cNXV * cNZW)] +
						4 * cWnp[i + (k * cNXW) + ((j + 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 1) * cNXV * cNZW)] -
						cWnp[i + (k * cNXW) + ((j + 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 2) * cNXV * cNZW)]);
				}
				else if (j == cNYW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)]* V[i + (k * cNXV) + ((j-1) * cNXV * cNZW)] -
						4 * cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 2) * cNXV * cNZW)] +
						cWnp[i + (k * cNXW) + ((j - 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 3) * cNXV * cNZW)]);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 1) * cNXV * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + (j * cNXV * cNZW)]
						+ (-4 * 0.025 + 1.0)*cWnp[i + (k * cNXW) + ((j + 1)* cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 1)* cNXV * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cWnp[i + (k * cNXW) + ((j + 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 2) * cNXV * cNZW)] +
						(-2 * 0.025 / 3.0*cWnp[i + (k * cNXW) + ((j + 3) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 3) * cNXV * cNZW)]));
				}
				else if (j == cNYW - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cWnp[i + (k * cNXW) + ((j + 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j) * cNXV * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + ((j-1) * cNXV * cNZW)]
						+ (4 * 0.09 + 1.0)*cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 2) * cNXV * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cWnp[i + (k * cNXW) + ((j - 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 3) * cNXV * cNZW)] +
						(2 * 0.09 / 3.0*cWnp[i + (k * cNXW) + ((j - 3) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 4) * cNXV * cNZW)]));
				}
				*/
				else
				{
					const float wijpk = cphiInterpolateY[i + (k * cNXW) + (j * cNXW * cNZW)];
					const float wijppk = cphiInterpolateY[i + (k * cNXW) + ((j + 1) * cNXW * cNZW)];

					const float wijmk = cphiInterpolateY[i + (k* cNXW) + ((j - 1)  * cNXW * cNZW)];
					const float wijmmk = cphiInterpolateY[i + (k * cNXW) + ((j - 2) * cNXW * cNZW)];

					const float vijpk = V[i + (k * cNXV) + (j * cNXV * (cNZV - 1))];
					const float vijppk = V[i + (k * cNXV) + ((j + 1) * cNXV * (cNZV - 1))];

					const float vijmk = V[i + (k * cNXV) + ((j - 1) * cNXV * (cNZV - 1))];
					const float vijmmk = V[i + (k * cNXV) + ((j - 2) * cNXV * (cNZV - 1))];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (-Fimp / 6.0*wijmmk * vijmmk -
						Eimp / 1.0*wijmk * vijmk +
						Eimp / 1.0*wijpk * vijpk +
						Fimp / 6.0*wijppk * vijppk);*/

					rhs[index] = bII / (3 * ch)*(wijppk*vijppk - wijmmk * vijmmk)
						+ aII / ch * (wijpk*vijpk - wijmk * vijmk);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidzWConvSystemTD(
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = cNXW * cNYW * cNZW;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYW; j++)
		for (int k = 0; k < cNZW; k++)
			for (int i = 0; i < cNXW; i++)
			{
				const int index = i + (k * cNXW) + (j * cNXW * cNZW);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == cNZW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if ((k == ckCyInit - 1 || k == ckCyInit - 2 || k == ckCyEnd - 1 || k == ckCyEnd)
					&& i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wijkp = cphiInterpolateZ[i + (k * cNXW) + (j * cNXW + (cNZW - 1))];

					const float wijkm = cphiInterpolateZ[i + ((k - 1) * cNXW) + (j * cNXW + (cNZW - 1))];

					rhs[index] = ((wijkp*wijkp)
						- (wijkm*wijkm)) /
						(ch);
				}
				else if (k == 1 || k == 2 || k == cNZW - 2 || k == cNZW - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wijkp = cphiInterpolateZ[i + (k * cNXW) + (j * cNXW + (cNZW - 1))];

					const float wijkm = cphiInterpolateZ[i + ((k - 1) * cNXW) + (j * cNXW + (cNZW - 1))];

					rhs[index] = ((wijkp*wijkp)
						- (wijkm*wijkm)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const float wijkp = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k+1)* cNXW) + ((j)* cNXW * cNZW)]
						);

					const float wijkm = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k-1)* cNXW) + ((j)* cNXW * cNZW)]
						);

					rhs[index] = ((wijkp*wijkp)
						- (wijkm*wijkm)) /
						(ch);
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] +
						4 * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == cNZW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] -
						4 * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] +
						cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)]
						+ (-4 * 0.025 + 1.0)*cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] +
						(-2 * 0.025 / 3.0*cWnp[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)]));
				}
				else if (k == cNZW - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNZW) + (j * cNXW * cNZW)]
						+ (4 * 0.09 + 1.0)*cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW *cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] +
						(2 * 0.09 / 3.0*cWnp[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)]));
				}
				*/
				else
				{
					const float wijkp = cphiInterpolateZ[i + (k * cNXW) + (j * cNXW + (cNZW - 1))];
					const float wijkpp = cphiInterpolateZ[i + ((k + 1) * cNXW) + (j * cNXW + (cNZW - 1))];

					const float wijkm = cphiInterpolateZ[i + ((k - 1) * cNXW) + (j * cNXW + (cNZW - 1))];
					const float wijkmm = cphiInterpolateZ[i + ((k - 2) * cNXW) + (j * cNXW + (cNZW - 1))];

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs[index] = 1.0 / ch * (-Fimp / 6.0*wijkmm * wijkmm -
						Eimp / 1.0*wijkm * wijkm +
						Eimp / 1.0*wijkp * wijkp +
						Fimp / 6.0*wijkpp * wijkpp);*/

					rhs[index] = bII / (3 * ch)*(wijkpp*wijkpp - wijkmm * wijkmm)
						+ aII / ch * (wijkp*wijkp - wijkm * wijkm);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreateInterpolateXSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = NXm * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = i + (k * NXm) + (j * NXm * NZ);

				if (i == 0 || i == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (i == NXm - 1 || i == NXm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]) / 2.0
						+ bI * (phi[i - 1 + (k*NX) + (j*NX*NZ)] + phi[i + 2 + (k*NX) + (j*NX*NZ)]) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreateInterpolateYSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NYm = NY - 1;
	const int n2 = NX * NYm * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + (k*NX) + ((j)*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]) / 2.0
						+ bI * (phi[i + (k*NX) + ((j - 1)*NX*NZ)] + phi[i + (k*NX) + ((j + 2)*NX*NZ)]) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreateInterpolateZSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * NY * NZm;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZm);

				if (k == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + ((k)*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]) / 2.0
						+ bI * (phi[i + ((k - 1)*NX) + (j*NX*NZ)] + phi[i + ((k + 2)*NX) + (j*NX*NZ)]) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreateInterpolateXSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = cNXP * cNYP * cNZP;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	int ii = 0;

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = i + (k * NXm) + (j * NXm * NZ);
				
				if (i == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 1 && (i == ciCyInit || i == ciCyInit - 1 || i == ciCyEnd - 1 || i == ciCyEnd)
					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 2 && (i == ciCyInit || i == ciCyInit - 1 || i == ciCyEnd - 1 || i == ciCyEnd)
					&& k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 0 && (i == ciCyInit - 1 || i == ciCyInit - 2 || i == ciCyEnd - 1 || i == ciCyEnd)
					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]) / 2.0
						+ bI * (phi[i - 1 + (k*NX) + (j*NX*NZ)] + phi[i + 2 + (k*NX) + (j*NX*NZ)]) / 2.0;
				}


				ii++;
			}

	const int maxi = cNXP * cNYP*cNZP;
	for (int i = ii; i < maxi; i++)
	{
		aa[i] = 0.0f;
		cc[i] = 0.0f;
		bb[i] = 1.0f;

		rhs[i] = 0.0f;
	}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreateInterpolateYSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NYm = NY - 1;
	const int n2 = cNXP*cNYP*cNZP;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	int ii = 0;

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else if (uvw != 1 && (j == 0 || j == 1 || j == NY - 1 || j == NY - 2))
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else if (uvw == 1 && (j == NY - 1 || j == NY - 2 || j == 0 || j == 1))
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else if (uvw == 0 && (j == cjCyInit - 1 || j == cjCyInit || j == cjCyEnd && j == cjCyEnd - 1)
					&& i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else if (uvw == 1 && (j == cjCyInit - 2 || j == cjCyInit - 1 || j == cjCyEnd && j == cjCyEnd - 1)
					&& i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else if (uvw == 2 && (j == cjCyInit - 1 || j == cjCyInit || j == cjCyEnd && j == cjCyEnd - 1)
					&& i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + (k*NX) + ((j)*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]) / 2.0
						+ bI * (phi[i + (k*NX) + ((j - 1)*NX*NZ)] + phi[i + (k*NX) + ((j + 2)*NX*NZ)]) / 2.0;
				}

				ii++;
			}

	const int maxi = cNXP * cNYP*cNZP;
	for (int i = ii; i < maxi; i++)
	{
		aa[i] = 0.0f;
		cc[i] = 0.0f;
		bb[i] = 1.0f;

		rhs[i] = 0.0f;
	}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreateInterpolateZSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = cNXP*cNYP*cNZP;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);
	
	int ii = 0;

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZm);

				if (k == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 0 && (k == ckCyInit || k == ckCyInit - 1 || k == ckCyEnd - 1 || k == ckCyEnd)
					&& i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 1 && (k == ckCyInit || k == ckCyInit - 1 || k == ckCyEnd - 1 || k == ckCyEnd)
					&& i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 2 && (k == ckCyInit - 1 || k == ckCyInit - 2 || k == ckCyEnd - 1 || k == ckCyEnd)
					&& i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + ((k)*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]) / 2.0
						+ bI * (phi[i + ((k - 1)*NX) + (j*NX*NZ)] + phi[i + ((k + 2)*NX) + (j*NX*NZ)]) / 2.0;
				}

				ii++;
			}

	const int maxi = cNXP * cNYP*cNZP;
	for (int i = ii; i < maxi; i++)
	{
		aa[i] = 0.0f;
		cc[i] = 0.0f;
		bb[i] = 1.0f;

		rhs[i] = 0.0f;
	}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidxStaggeredSystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = NXm * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = i + (k * NXm) + (j * NXm * NZ);

				if (i == 0 || i == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (i == NXm - 1 || i == NXm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidyStaggeredSystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NYm = NY - 1;
	const int n2 = NX * NYm * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0 || j == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else if (j == NYm - 1 || j == NYm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidzStaggeredSystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * NY * NZm;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZm);

				if (k == 0 || k == 1/* || true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (k == NZm - 1 || k == NZm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + ((k + 2)*NX) + (j*NX*NZ)] - phi[i + ((k - 1)*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + ((k)*NX) + (j*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidxStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const

{
	const int NXm = NX - 1;
	const int n2 = NXm * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = i + (k * NXm) + (j * NXm * NZ);

				if (i == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidyStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NYm = NY - 1;
	const int n2 = NX * NYm * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidzStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * NY * NZm;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZm);

				if (k == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + ((k + 2)*NX) + (j*NX*NZ)] - phi[i + ((k - 1)*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + ((k)*NX) + (j*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidxStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const

{
	const int NXm = NX - 1;
	const int n2 = NXm * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = i + (k * NXm) + (j * NXm * NZ);

				if (uvw == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + 1 + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				/*else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}*/
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidyStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NYm = NY - 1;
	const int n2 = NX * NYm * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (uvw == 1 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j + 1)*NX*NZ)]) / ch;
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				/*else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}*/
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::CreatedphidzStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * NY * NZm;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZm);

				if (uvw == 2 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + ((k + 1)*NX) + (j*NX*NZ)]) / ch;
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				/*else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + ((k + 2)*NX) + (j*NX*NZ)] - phi[i + ((k - 1)*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + ((k)*NX) + (j*NX*NZ)]) / ch;
				}*/
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::Created2phidx2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i - 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + 1 + (k*NX) + (j*NX*NZ)])
						/(ch*ch);
				}
				else*/
				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						// force value to be + phi at i+1;
						rhs[index] = (2.0 * -phi[i + 1 + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + 1 + (k*NX) + (j*NX*NZ)] +
							4.0 * phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i + 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + 1 + (k*NX) + (j*NX*NZ)] +
							4.0 * phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i + 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}

				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						// force value to be + phi at i-1;
						rhs[index] = (2.0 * -phi[i - 1 + (k*NX) + (j*NX*NZ)] - 5.0*phi[i - 1 + (k*NX) + (j*NX*NZ)] +
							4.0 * phi[i - 2 + (k*NX) + (j*NX*NZ)] - phi[i - 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i - 1 + (k*NX) + (j*NX*NZ)] +
							4.0 * phi[i - 2 + (k*NX) + (j*NX*NZ)] - phi[i - 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						// force value to be + phi at 0;
						rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
				}
				else if (uvw == 1 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 1 && i == ciCyInit && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i - 1 + (k*NX) + (j*NX*NZ)] - 5.0*phi[i - 1 + (k*NX) + (j*NX*NZ)] +
						4.0 * phi[i - 2 + (k*NX) + (j*NX*NZ)] - phi[i - 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 1 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 1 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + 1 + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + 1 + (k*NX) + (j*NX*NZ)] +
						4.0 * phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i + 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i - 1 + (k*NX) + (j*NX*NZ)] - 5.0*phi[i - 1 + (k*NX) + (j*NX*NZ)] +
						4.0 * phi[i - 2 + (k*NX) + (j*NX*NZ)] - phi[i - 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + 1 + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + 1 + (k*NX) + (j*NX*NZ)] +
						4.0 * phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i + 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + 2 + (k*NX) + (j*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i - 2 + (k*NX) + (j*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i - 1 + (k*NX) + (j*NX*NZ)]);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::Created2phidy2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j-1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + (k*NX) + ((j+1)*NX*NZ)])
						/ (ch*ch);
				}
				else */if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2)
					{
						rhs[index] = (2.0 * -phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j + 1)*NX*NZ)] +
							4.0 * phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j + 3)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j + 1)*NX*NZ)] +
							4.0 * phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j + 3)*NX*NZ)]) / (ch*ch);
					}

				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2)
					{
						rhs[index] = (2.0 * -phi[i + (k*NX) + ((j - 1)*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j - 1)*NX*NZ)] +
							4.0 * phi[i + (k*NX) + ((j - 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 3)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j - 1)*NX*NZ)] +
							4.0 * phi[i + (k*NX) + ((j - 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 3)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;
					if (uvw == 2)
					{
						rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2)
					{
						rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (uvw == 0 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 0 && j == cjCyInit && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + (k*NX) + ((j - 1)*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j - 1)*NX*NZ)] +
						4.0 * phi[i + (k*NX) + ((j - 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 3)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 0 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 0 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j + 1)*NX*NZ)] +
						4.0 * phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j + 3)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 2 && j == cjCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 2 && j == cjCyInit && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + (k*NX) + ((j - 1)*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j - 1)*NX*NZ)] +
						4.0 * phi[i + (k*NX) + ((j - 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 3)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 2 && j == cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 2 && j == cjCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j + 1)*NX*NZ)] +
						4.0 * phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j + 3)*NX*NZ)]) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + (k*NX) + ((j + 2)*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j - 2)*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j - 1)*NX*NZ)]);
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredFieldFFast::Created2phidz2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<float>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k-1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + ((k+1)*NX) + ((j)*NX*NZ)])
						/ (ch*ch);
				}
				else*/
				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1)
					{
						rhs[index] = (2.0 * -phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}

				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1)
					{
						rhs[index] = (2.0 * -phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] - 5.0*phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k - 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k - 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k - 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k - 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1)
					{
						rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1)
					{
						rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (uvw == 1 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (-phi[i + ((k)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 1 && k == ckCyInit && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] - 5.0*phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] +
						4.0 * phi[i + ((k - 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k - 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 1 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
						4.0 * phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 1 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}

				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (-phi[i + ((k)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] - 5.0*phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] +
						4.0 * phi[i + ((k - 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k - 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * -phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
						4.0 * phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k - 2)*NX) + ((j)*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]);
				}
			}

	return n2;
}

void SquareCylinder3DStaggeredFieldFFast::EnforceContinuityVelnp(float dt,
	std::vector<float>& aa,
	std::vector<float>& bb,
	std::vector<float>& cc,
	std::vector<int>& ptr,
	std::vector<int>& col,
	std::vector<float>& val,
	std::vector<float>& rhs,
	LSSolverF& lssolver)
{
	int nn = 0;

	if (true)
	{
		nn = CreatedphidxStaggeredSystemTDSimplified(cUnp, cNXU, cNYU, cNZU, 0,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		nn = CreatedphidyStaggeredSystemTDSimplified(cVnp, cNXV, cNYV, cNZV, 1,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		nn = CreatedphidzStaggeredSystemTDSimplified(cWnp, cNXW, cNYW, cNZW, 2,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}
	else
	{
		nn = CreatedphidxStaggeredSystemTD(cUnp, cNXU, cNYU, cNZU,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		nn = CreatedphidyStaggeredSystemTD(cVnp, cNXV, cNYV, cNZV,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		nn = CreatedphidzStaggeredSystemTD(cWnp, cNXW, cNYW, cNZW,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}

	nn = CreateSparsePressureEquation(cdphidx, cdphidy, cdphidz,
		dt, ptr, col, val, rhs);

	if (cFirstIteration)
	{
		lssolver.PrecondtionCRS(nn, ptr, col, val);

		cFirstIteration = false;
	}

	float errOutPEqt = 0.0;
	int nitPEqt = 0;

	lssolver.SolvePreconditionedCRS(std::move(cPn), nn, ptr, col, val, rhs, errOutPEqt, nitPEqt);

	std::cout << "Pressure eq solved its " << nitPEqt << " err " << errOutPEqt << std::endl;


	// project p into uvw

	if (true)
	{
		nn = CreatedphidxStaggeredSystemTDSimplified(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		nn = CreatedphidyStaggeredSystemTDSimplified(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		nn = CreatedphidzStaggeredSystemTDSimplified(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}
	else
	{
		// set dp/dn = 0
		for (int j = 0; j < cNYP; j++)
			for (int k = 0; k < cNZP; k++)
			{
				cPn[0 + (k*cNXP) + (j*cNXP*cNZP)] = cPn[1 + (k*cNXP) + (j*cNXP*cNZP)];
				cPn[cNXP - 1 + (k*cNXP) + (j*cNXP*cNZP)] = cPn[cNXP - 2 + (k*cNXP) + (j*cNXP*cNZP)];
			}

		nn = CreatedphidxStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		// set dp/dn = 0
		for (int k = 0; k < cNZP; k++)
			for (int i = 0; i < cNXP; i++)
			{
				cPn[i + (k*cNXP) + ((cNYP - 1)*cNXP*cNZP)] = cPn[i + (k*cNXP) + ((cNYP - 2)*cNXP*cNZP)];
				cPn[i + (k*cNXP) + (0 * cNXP*cNZP)] = cPn[i + (k*cNXP) + (1 * cNXP*cNZP)];
			}

		nn = CreatedphidyStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		// set dp/dn = 0
		for (int j = 0; j < cNYP; j++)
			for (int i = 0; i < cNXP; i++)
			{
				cPn[i + (0 * cNXP) + (j*cNXP*cNZP)] = cPn[i + (1 * cNXP) + (j*cNXP*cNZP)];
				cPn[i + ((cNZP - 1)*cNXP) + (j*cNXP*cNZP)] = cPn[i + ((cNZP - 2)*cNXP) + (j*cNXP*cNZP)];
			}

		nn = CreatedphidzStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}


	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				if (i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					continue;
				}
				cUnp[i + (k*cNXU) + (j*cNXU*cNZU)] -=
					(dt / crho)*(cdphidx[i + (k*(cNXP - 1)) + (j*(cNXP - 1)*cNZP)]);
			}
		}
	}

	for (int j = 1; j < cNYV - 1; j++)
	{
		for (int k = 1; k < cNZV - 1; k++)
		{
			for (int i = 1; i < cNXV - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd &&  k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					continue;
				}
				cVnp[i + (k*cNXV) + (j*cNXV*cNZV)] -=
					(dt / crho)*(cdphidy[i + (k*cNXP) + ((j)*cNXP*cNZP)]);
			}
		}
	}

	for (int j = 1; j < cNYW - 1; j++)
	{
		for (int k = 1; k < cNZW - 1; k++)
		{
			for (int i = 1; i < cNXW - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					continue;
				}
				cWnp[i + (k*cNXW) + (j*cNXW*cNZW)] -=
					(dt / crho)*(cdphidz[i + ((k)*cNXP) + (j*cNXP*(cNZP - 1))]);
			}
		}
	}
}

void SquareCylinder3DStaggeredFieldFFast::SaveField(const std::string& baseFileName) const
{
	std::ofstream uFile(baseFileName + "U.datf", std::ios::binary);

	uFile.write(reinterpret_cast<const char*>(&cNXU), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&cNYU), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&cNZU), sizeof(int64_t));

	uFile.write(reinterpret_cast<const char*>(cUnp.data()), cUnp.size() * sizeof(float));

	uFile.close();

	std::ofstream vFile(baseFileName + "V.datf", std::ios::binary);

	vFile.write(reinterpret_cast<const char*>(&cNXV), sizeof(int64_t));
	vFile.write(reinterpret_cast<const char*>(&cNYV), sizeof(int64_t));
	vFile.write(reinterpret_cast<const char*>(&cNZV), sizeof(int64_t));

	vFile.write(reinterpret_cast<const char*>(cVnp.data()), cVnp.size() * sizeof(float));

	vFile.close();

	std::ofstream wFile(baseFileName + "W.datf", std::ios::binary);

	wFile.write(reinterpret_cast<const char*>(&cNXW), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&cNYW), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&cNZW), sizeof(int64_t));

	wFile.write(reinterpret_cast<const char*>(cWnp.data()), cWnp.size() * sizeof(float));

	wFile.close();

	std::ofstream pFile(baseFileName + "P.datf", std::ios::binary);

	pFile.write(reinterpret_cast<const char*>(&cNXP), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&cNYP), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&cNZP), sizeof(int64_t));

	pFile.write(reinterpret_cast<const char*>(cPn.data()), cPn.size() * sizeof(float));

	pFile.close();
}

void SquareCylinder3DStaggeredFieldFFast::ReadField(const std::string& basefilename)
{

	std::ifstream uFile(basefilename + "U.datf", std::ios::binary);

	uFile.read(reinterpret_cast<char*>(&cNXU), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&cNYU), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&cNZU), sizeof(int64_t));

	cUn.clear();
	cUn.resize(cNXU * cNYU * cNZU);

	uFile.read(reinterpret_cast<char*>(cUn.data()), cUn.size() * sizeof(float));

	uFile.close();

	cUnp = cUn;

	std::ifstream vFile(basefilename + "V.datf", std::ios::binary);

	vFile.read(reinterpret_cast<char*>(&cNXV), sizeof(int64_t));
	vFile.read(reinterpret_cast<char*>(&cNYV), sizeof(int64_t));
	vFile.read(reinterpret_cast<char*>(&cNZV), sizeof(int64_t));

	cVn.clear();
	cVn.resize(cNXV * cNYV * cNZV);

	vFile.read(reinterpret_cast<char*>(cVn.data()), cVn.size() * sizeof(float));

	vFile.close();

	cVnp = cVn;

	std::ifstream wFile(basefilename + "W.datf", std::ios::binary);

	wFile.read(reinterpret_cast<char*>(&cNXW), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&cNYW), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&cNZW), sizeof(int64_t));

	cWn.clear();
	cWn.resize(cNXW * cNYW * cNZW);

	wFile.read(reinterpret_cast<char*>(cWn.data()), cWn.size() * sizeof(float));

	wFile.close();

	cWnp = cWn;

	std::ifstream pFile(basefilename + "P.datf", std::ios::binary);

	pFile.read(reinterpret_cast<char*>(&cNXP), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&cNYP), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&cNZP), sizeof(int64_t));

	cPn.clear();
	cPn.resize(cNXP * cNYP * cNZP);

	pFile.read(reinterpret_cast<char*>(cPn.data()), cPn.size() * sizeof(float));

	pFile.close();
}