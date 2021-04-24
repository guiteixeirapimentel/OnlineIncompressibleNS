#pragma once
#include <vector>
#include <fstream>

#include "Field3DIKJ.h"
#include "Field3DKJI.h"
#include "Field3DJIK.h"

class SquareCylinder3DStaggeredField
{
public:
	//typedef Field3DIKJ field3D;

	const double Dimp = 0.3793894912;
	const double Eimp = 1.57557379;
	const double Fimp = 0.183205192;
	const double etastar = -2.0;

	const double qm2 = (-Fimp / 4.0) + (etastar / 300.0);
	const double qm1 = (-Eimp / 2.0) + (etastar / 30.0);
	const double qp1 = (Eimp / 2.0) + (etastar / 30.0);
	const double qp2 = (Fimp / 4.0) + (etastar / 300.0);
	const double q0 = (-11.0*etastar) / 150.0;
	const double pip = Dimp + (etastar / 60.0);
	const double pim = Dimp - (etastar / 60.0);

	const double betaj2 = -0.025; //-0.025;
	const double betajm2 = 0.09;  //0.09;

	const double alfa2imp = 2.0 / 11.0;//0.26;
	const double a = (4.0 / 3.0)*(1.0 - alfa2imp);
	const double b = (1.0 / 3.0)*(10.0*alfa2imp - 1.0);

	const double alfaI = 3.0 / 10.0; // 0.42;
	const double aI = (9.0 + 10.0*alfaI) / 8.0;
	const double bI = (6 * alfaI - 1.0) / 8.0;

	const double alfaII = 0.216; //9.0 / 62.0;
	const double aII = (9.0 - 6.0*alfaII) / 8.0;
	const double bII = (22.0*alfaII - 1.0) / 8.0;

	const double velPosSUCSalfam1 = 1.0 / 2.0;
	const double velPosSUCSalfap1 = 1.0 / 6.0;
	const double velPosSUCScm2 = -1.0 / 18.0;
	const double velPosSUCScm1 = -1.0;
	const double velPosSUCSc0 = 1.0 / 2.0;
	const double velPosSUCSc1 = 5.0 / 9.0;

	const double velNegSUCSalfam1 = 1.0 / 6.0;
	const double velNegSUCSalfap1 = 1.0 / 2.0;
	const double velNegSUCScm1 = -5.0 / 9.0;
	const double velNegSUCSc0 = -1.0 / 2.0;
	const double velNegSUCSc1 = 1.0;
	const double velNegSUCSc2 = 1.0 / 18.0;

public:
	SquareCylinder3DStaggeredField(double rho, double mu, double h, double l, double L, double H, double W, double freestream);
	SquareCylinder3DStaggeredField() = delete;
	~SquareCylinder3DStaggeredField();

	void SetBCFlatPlate();
	void UpdateBCFlatPlate();

	void SaveIKCutToCSV(const std::string& filename, int j);

	void SaveIKCutToCSV(const std::string& nomeArq, const Field3D& campo, int NX, int NY, int NZ, int j);

	void SaveJKCutToCSV(const std::string& filename, int i);

	void SaveIJCutToCSV(const std::string& filename, int k);

	void SaveIJCutToCSV(const std::string& nomeArq, const Field3D& campo, int NX, int NY, int NZ, int k);

	int CreatedphidxSystemTDOUCS3(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidySystemTDOUCS3(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzSystemTDOUCS3(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreatedphidxSystemTDSUCS(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const Field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidySystemTDSUCS(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const Field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzSystemTDSUCS(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const Field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreatedphidxSystemTDUPWIND(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const Field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidySystemTDUPWIND(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const Field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzSystemTDUPWIND(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const Field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreatedphidxStaggeredSystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidyStaggeredSystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzStaggeredSystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreatedphidxStaggeredSystemTDSimplified(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidyStaggeredSystemTDSimplified(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzStaggeredSystemTDSimplified(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreatedphidxStaggeredSystemTDSimplified(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidyStaggeredSystemTDSimplified(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzStaggeredSystemTDSimplified(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int Created2phidx2SystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int Created2phidy2SystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int Created2phidz2SystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int Created2phidx2SystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int Created2phidy2SystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int Created2phidz2SystemTD(const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	/*int CreatedphidxUConvSystemTD(
		const std::vector<double>& phiInterpolateX,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidyUConvSystemTD(
		const std::vector<double>& phiInterpolateY,
		const field3D& V,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzUConvSystemTD(
		const std::vector<double>& phiInterpolateZ,
		const field3D& W,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreatedphidxVConvSystemTD(
		const std::vector<double>& phiInterpolateX,
		const field3D& U,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidyVConvSystemTD(
		const std::vector<double>& phiInterpolateY,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzVConvSystemTD(
		const std::vector<double>& phiInterpolateZ,
		const field3D& W,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreatedphidxWConvSystemTD(
		const std::vector<double>& phiInterpolateX,
		const field3D& U,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreatedphidyWConvSystemTD(
		const std::vector<double>& phiInterpolateY,
		const field3D& V,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreatedphidzWConvSystemTD(
		const std::vector<double>& phiInterpolateZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;*/

	int CreateInterpolateXPreIntSystemTD(
		const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreateInterpolateYPreIntSystemTD(
		const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreateInterpolateZPreIntSystemTD(
		const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreateInterpolateXSystemTD(
		const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DIKJ& rhs) const;

	int CreateInterpolateYSystemTD(
		const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DJIK& rhs) const;

	int CreateInterpolateZSystemTD(
		const Field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		Field3DKJI& rhs) const;

	int CreateSparsePressureEquation(const Field3D& dUdx,
		const Field3D& dVdy,
		const Field3D& dWdz,
		const double dt,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		Field3DIKJ& rhs) const;

	//void CalcdphidtUmomDiv(field3D& dphidt,
	//	std::vector<double>& aa,
	//	std::vector<double>& bb,
	//	std::vector<double>& cc);

	//void CalcdphidtVmomDiv(field3D& dphidt,
	//	std::vector<double>& aa,
	//	std::vector<double>& bb,
	//	std::vector<double>& cc);

	//void CalcdphidtWmomDiv(field3D& dphidt,
	//	std::vector<double>& aa,
	//	std::vector<double>& bb,
	//	std::vector<double>& cc);

	void CalcdphidtUmomAdv(Field3D& dphidt,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc);

	void CalcdphidtVmomAdv(Field3D& dphidt,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc);

	void CalcdphidtWmomAdv(Field3D& dphidt,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc);

	void EnforceContinuityVelnp(double dt,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		Field3DIKJ& rhs,
		bool show,
		class LSSolver& lssolver);

	void SaveField(const std::string& baseFileName) const;

	void ReadField(const std::string& basefilename);

	void InterpolateUVWInRespDirections(
		std::vector<double>& aa, std::vector<double>& bb, std::vector<double>& cc,
		std::vector<double>& aa1, std::vector<double>& bb1, std::vector<double>& cc1,
		std::vector<double>& aa2, std::vector<double>& bb2, std::vector<double>& cc2);

public:
	const double crho;
	const double cmu;
	const double ch;

	const double cL;
	const double cH;
	const double cW;

	const double cl;

	const int cDistRelFrenteRect;

	int64_t cNXU;
	int64_t cNYU;
	int64_t cNZU;
	   
	int64_t cNXV;
	int64_t cNYV;
	int64_t cNZV;
	   
	int64_t cNXW;
	int64_t cNYW;
	int64_t cNZW;
	   
	int64_t cNXP;
	int64_t cNYP;
	int64_t cNZP;
	   
	int64_t cNXCy;
	int64_t cNYCy;
	int64_t cNZCy;
	   
	int64_t ciCyInit;
	int64_t cjCyInit;
	int64_t ckCyInit;
	   
	int64_t ciCyEnd;
	int64_t cjCyEnd;
	int64_t ckCyEnd;

	double cUfreestream;

	bool cFirstIteration;

	Field3DIKJ cUn;
	Field3DIKJ cVn;
	Field3DIKJ cWn;
	Field3DIKJ cPn;

	Field3DIKJ cUnp;
	Field3DIKJ cVnp;
	Field3DIKJ cWnp;
	/*field3D cPnp;*/

	Field3DIKJ cPhi;

	Field3DIKJ cdphidx;
	Field3DJIK cdphidy;
	Field3DKJI cdphidz;

	Field3DIKJ cphiInterpolateX;
	Field3DJIK cphiInterpolateY;
	Field3DKJI cphiInterpolateZ;

	Field3DIKJ cd2phidx2;
	Field3DJIK cd2phidy2;
	Field3DKJI cd2phidz2;
	
	Field3DIKJ cdUdtn;
	Field3DIKJ cdUdt1;
	Field3DIKJ cdUdt2;
	Field3DIKJ cdUdt3;

	Field3DIKJ cdVdtn;
	Field3DIKJ cdVdt1;
	Field3DIKJ cdVdt2;
	Field3DIKJ cdVdt3;

	Field3DIKJ cdWdtn;
	Field3DIKJ cdWdt1;
	Field3DIKJ cdWdt2;
	Field3DIKJ cdWdt3;

	Field3DIKJ cUnpc;
	Field3DJIK cVnpc;
	Field3DKJI cWnpc;
};