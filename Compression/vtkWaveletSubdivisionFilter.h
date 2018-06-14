/*=========================================================================

  Program:   Mailleur 3D multi-r?olution
  Module:    vtkWaveletSubdivisionFilter.h
  Language:  C++
  Date:      2002/05
  Auteur:    Sebastien VALETTE
  This software is governed by the GPL license (see License.txt)
=========================================================================*/

#ifndef __vtkWaveletSubdivisionFilter_h
#define __vtkWaveletSubdivisionFilter_h
#include <vector>
#include <math.h>
#include <vtkObjectFactory.h>
#include <vtkDoubleArray.h> 
#include <vtkCommand.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>


#include "vtkSurface.h"
#include "vtkSparseMatrix.h"
#include "vtkArithmeticCoder.h"


#define CHILD 1
#define PARENT 2

/**
 *  A Wavelet-based Multiresolution analysis filter
 *  The vtkWaveletSubdivisionFilter is an efficient reversible mesh 
 *  simplification and reconstruction Class. It can process integer or
 *  double coordinates, handles compression.
 */
class VTK_EXPORT vtkWaveletSubdivisionFilter : public vtkObject
{
	class contour_element
	{
	public:
		contour_element *PreviousElement,*NextElement;
		vtkIdType edge;
	};

	class edge
	{
	public:
		contour_element *contour;
		vtkIdType child;
		char type;
		edge() : contour(0){};
	};

	class face
	{
	public:
		char lowrestype;
		char hirestype;
		vtkIdType vertex; //indice de sommet servant au détrompage (subdivision 1:3 ou 1:2 avec permutation d'arête)
		face() : lowrestype (-1),hirestype(-1),vertex(-1){};
	};

	class vertex
	{
	public:
		char type;
		short valence;
		vtkIdType old_to_new;
		vtkIdType new_to_old;
		vtkIdType parent1;
		vtkIdType parent2;

		vertex() {old_to_new=-1;new_to_old=-1;type=-1;};
	};

public:
   //BTX
	vtkArithmeticCoder *ArithmeticCoder;
   //ETX

	int ConnectivityOnly;
	int LiftNum;
	int LiftDen;

	int GeometryPrediction;

	static vtkWaveletSubdivisionFilter *New();
	vtkTypeMacro(vtkWaveletSubdivisionFilter,vtkObject);



	void SetDisplayEfficiency(int D)
		{this->DisplayEfficiency=D;};
	int DisplayEfficiency;

/**
 ** Sparse Matrix Alpha.
 ** Used for Lifting
 **/
	vtkSparseMatrix *Alpha;

/**
 ** Sparse Matrix P.
 ** Used for Lifting
 **/
	vtkSparseMatrix *P;

/**
 ** Geometry compression estimations
 **/
	double EstimatedGeometryBits;

/**
 ** =0 if the mesh was entirely simplified. If >0, the simplification
 ** is not reversible.
 **/
	int remaining_faces;

/**
 ** Number of 4:1 mergings
 **/
	int M4to1;

/**
 ** Number of 3:1 mergings
 **/
	int M3to1;

/**
 ** Number of 2:1 mergings
 **/
	int M2to1;

/**
 ** Number of 1:1 mergings
 **/
	int M1to1;


/**
 ** Number of edge swaps (type 1)
 **/
	int Swap1;

/**
 ** Number of edge swaps (type 2)
 **/
	int Swap2;


/**
 ** Number of raw bits needed for encoding midpoints (before entropy coding)
 **/
	int BitsMidPoints;

/**
 ** Number of raw bits set to 1 needed for encoding midpoints (before entropy coding)
 **/
	int BitsMidPoints1;

/**
 ** Number of supplementary bits for 1:3 subdivisions
 **/
	int Bits1to3;

/**
 ** Number of supplementary bits for 1:2 subdivisions
 **/
	int Bits1to2;

/**
 ** Number of supplementary bits for edge swaps
 **/
	int BitsSwap;

/**
 ** Estimation of bits needed for midpoints encoding
 **/
	double RealBitsMidPoints;

/**
 ** Quantization factor
 **/
	int Quantization;
/**
 ** Number of bits needed for encoding the low resolution mesh
 **/
	int LowresBits;

/**
 ** Data size (in bits) of the encoded connectivity
 **/
	int ConnectivityBits;

/**
 ** Data size (in bits) of the encoded geometry
 **/
	double GeometryBits;

/**
 ** Takes this->MergeInput mesh and simplifies it to a new mesh :
 ** this->MergeOutput
 **/
	void SolveInverseProblem (vtkIdType initial_face);

/**
 ** After this call, the filter will use fast 0-ring lifting
 **/
	void SetFastLiftingOn() {this->ProcessLifting=2;};

/**
 ** After this call, the filter will use lifting
 **/
	void SetLiftingOn() {this->ProcessLifting=1;};

/**
 ** After this call, the filter will use : 
 ** 0 : no lifting
 ** 1 : lifting
 ** 2 : fast 0-ring Lifting
 **/
	void SetLifting(int lift) {this->ProcessLifting=lift;};

/**
 ** After this call, the filter will not process lifting
 **/
	void SetLiftingOff() {this->ProcessLifting=0;};


/**
 ** Returns the lifting process state :
 ** 0 : no lifting
 ** 1 : lifting
 ** 2 : fast 0-ring lifting
 **/
	int IsLiftingOn() {return (this->ProcessLifting);};

/**
 ** Sets the size of the wavelet support for lifting (0-ring, 1-ring etc...).
 ** Note : 0<=Radius
 **/
	void SetLiftingRadius(int Radius) {this->LiftingRadius=Radius;};

/**
 ** Returns the size of the wavelet support for lifting (0-ring, 1-ring etc...)
 **/
	int GetLiftingRadius() {return (this->LiftingRadius);};

/**
 ** After calling this->SolveInverseProblem, this function returns the number of non merged faces
 ** Note : if the mesh was correctly simplified, the function must return 0
 **/
	int GetRemainingFaces() {return (this->remaining_faces);};

/**
 ** Sets the arithmetic used for the wavelet decomposition.
 ** 0 : double arithmetics, 1 : integer arithmetics
 ** Note : for Compression, you must use integer arithmetics
 **/
	void SetArithmeticType (int Type) {this->ArithmeticType=Type;};

/**
 ** Sets the arithmetic used for the wavelet decomposition.
 ** 0 : double arithmetics, 1 : integer arithmetics
 ** Note : for Compression, you must use integer arithmetics
 **/
	int GetArithmeticType () {return (this->ArithmeticType);};


/**
 ** Sets the subdivision type :
 ** 0 : regular subdivision (Loop).
 ** 1 : irregular subdivision.
 ** Note : by default, the type is 0, but when calling this->Solve InverseProblem, it is automatically set to 1.
 **/
	void SetSubdivisionType(int Type) {this->SubdivisionType=Type;};

/**
 ** Returns the subdivision type :
 ** 0 : regular subdivision .
 ** 1 : irregular subdivision.
 ** Note : by default, the type is 0, but when calling this->Solve InverseProblem, it is automatically set to 1.
 **/
	int GetSubdivisionType() {return (this->SubdivisionType);};


/**
 ** Array used for storing the wavelets coefficients after calling this->Approximate().
 ** Note : this array is used only when double arithmetic is used.
 **/
	vtkDoubleArray *Wavelets;

/**
 ** Array used for storing the wavelets coefficients after calling this->Approximate().
 ** Note : this array is used only when integer arithmetic is used
 **/
	vtkIntArray *IntegerWavelets;

/**
 ** Subdivides this->SubdivisionInput mesh to this->SubdivisionOutput mesh.
 ** Note : it is also used for writing and reading data to file using the IOType parameter.
 **/
	virtual void Subdivide();

/**
 ** Sets the direction of the file transfer before calling this->Execute() or this->Subdivide(): 
 ** 0 : no file transfer (for compression, you have to do it once for each filter, to build the subdivided meshes)
 ** 1 : write to file
 ** 2 : read from file
 **/
	void SetIOType(int Type) {this->IOType=Type;};

/**
 ** Returns the direction of the file transfer before calling this->Execute() or this->Subdivide(): 
 ** 0 : no file transfer (for compression, you have to do it once for each filter, to build the subdivided meshes)
 ** 1 : write to file
 ** 2 : read from file
 **/
	int GetIOType() {return (this->IOType);};


/**
 ** Returns the simplified mesh after calling this->SolveInverseProblem
 **/
	vtkSurface *GetMergeOutput() {return (this->MergeOutput);};

/**
 ** Returns the mesh to be simplified
 **/
	vtkSurface *GetMergeInput() {return (this->MergeInput);};

/**
 ** Sets the mesh to be simplified
 **/
	void SetMergeInput(vtkSurface *Input) {this->MergeInput=Input;this->Modified();Input->Register(this);};

/**
 ** Sets an array of PointIds, used for vertices renumbering
 **/
	void SetPointsIds(vtkIdList *PtIds) {this->PointsIds=PtIds;this->Modified();};

/**
 **  Sets the mesh to be subdivided
 **/
	void SetInput(vtkSurface *input);

/**
 **  Void function, written for vtkPolyDataToPolyDatafilter compatibility
 **/
	vtkPolyData * GetInput();

/**
 **   Returns the mesh to be subdivided
 **/
	vtkSurface * GetSubdivisionInput(){return (this->SubdivisionInput);};

/**
 **   Returns the result of the subdivision
 **/
	vtkSurface * GetOutput() {return (this->SubdivisionOutput);};

/**
 **   Performs the geometrical approximation of M(j) to M(j-1).
 ** NOTE : must be called only after calling this->Subdivide or this->Execute.
 **/ 
	void Approximate();

/**
 **   Performs the geometrical reconstruction of M(j) from M(j-1).
 ** NOTE : must be called only after calling this->Subdivide or this->Execute
 **/ 
	void Reconstruct();

/**
 **   Do not use, this function is only made for research purpose....
 **/ 
	void DisplayWavelet();
	void DisplayScalingFunction();
	void ClearWavelets();

/**
 **  Sets wether the filter will use a geometrical criterion WGC.
 ** 0 : no WGC; 1: use WGC
 **/ 
	void SetGeometryCriterion(int criterion) {this->GeometryCriterion=criterion;};

/**
 **  Sets the curvature threshold for the WGC.
 ** Note : 0<threshold<1, threshold=0.3 is recommended
 **/ 

	void SetCurvatureTreshold(double treshold) {this->CurvatureTreshold=treshold;};

/**
 **  Sets the Wavelet threshold for the WGC.
 ** Note : 0<threshold<1, threshold=0.25 is recommended.
 **/ 
	void SetWaveletTreshold(double treshold) {this->WaveletTreshold=treshold;};

	void Orthogonalize ();

	void SetGeometryPrediction(int p)
	{ this->GeometryPrediction=p;};

	/**
	 ** for compression purpose
	 **
	 **/
	void WriteCoefficients();
	void ReadCoefficients();

	vtkIdList *TreeFirstChildEdge;
	vtkIdList *TreeNextChildEdge;
	vtkIdList *TreeParentEdge;
	vtkIdList *EdgeMidPoints;

//BTX
	std::vector<edge> edgesvector;
	std::vector<face> faces;
	std::vector<vertex> vertices;
//ETX
	
	void InitTree();
	void AddChildEdgeToTree(vtkIdType Child, vtkIdType Parent);

	void Subdivision();

	void SaveWavelets(const char *Filename);


protected:

	vtkWaveletSubdivisionFilter():FirstElement(0), LastElement(0) ,
	FirstElementRegular(0),LastElementRegular(0){
		this->SubdivisionType=0;
		this->ProcessLifting=0;
		this->LiftingRadius=0;
		this->Alpha=vtkSparseMatrix::New();
		this->P=vtkSparseMatrix::New();

		this->SubdivisionOutput=vtkSurface::New(); 
		this->SubdivisionOutput->Register(this);
		this->SubdivisionOutput->Delete();
		this->MergeOutput=vtkSurface::New();
		this->MergeInput=0;
		this->SubdivisionInput=0;
		this->Wavelets=0;
		this->IntegerWavelets=0;
		this->IOType=0;
		this->ArithmeticType=0;
		this->EdgeMidPoints=vtkIdList::New();
		this->IOStarted=0;
		this->GeometryCriterion=0;
		this->CurvatureTreshold=0;
		this->WaveletTreshold=0;
		this->CoordinatesCoupling=0;
		this->vlist=vtkIdList::New();
		this->CellTypes=vtkIntArray::New();
		this->DisplayEfficiency=0;
		this->GeometryPrediction=0;
		this->LiftNum=0;
		this->LiftDen=100;
		this->RegularSubdivisionFlag=0;
		this->ConnectivityOnly=0;
		this->TreeFirstChildEdge=vtkIdList::New();
		this->TreeNextChildEdge=vtkIdList::New();
		this->TreeParentEdge=vtkIdList::New();
		this->ArithmeticCoder=0;
		this->FileType=0;
		}; 

  	~vtkWaveletSubdivisionFilter() {
			this->Alpha->Delete();
			this->P->Delete();
		if (this->Wavelets)
            this->Wavelets->Delete();
		if (this->IntegerWavelets)
            this->IntegerWavelets->Delete();
		this->EdgeMidPoints->Delete();
		this->TreeFirstChildEdge->Delete();
		this->TreeNextChildEdge->Delete();
		this->TreeParentEdge->Delete();
		this->vlist->Delete();
		this->CellTypes->Delete();

		if (this->SubdivisionOutput)
			this->SubdivisionOutput->UnRegister(this);
		this->MergeOutput->UnRegister(this);
		if (this->MergeInput)
				this->MergeInput->UnRegister(this);
		};

	vtkWaveletSubdivisionFilter(const vtkWaveletSubdivisionFilter&) {};
	void operator=(const vtkWaveletSubdivisionFilter&) {};


	void Execute() {};

	vtkSurface *SubdivisionInput;
	vtkSurface *SubdivisionOutput;
	void EncodeMidpoint(int code);

	int DecodeMidpoint();
	int IOStarted;


private:

	void ComputeNewVertexCoordinates(vtkIdType p1,vtkIdType p2, double *V1);
	void ComputeVertexContribution(vtkIdType p1,vtkIdType p2, double *V1);
	void ComputeVertexContribution2(vtkIdType p1,vtkIdType p2, double *V1);

	void ConquerRegularSubdivisions();
	void SolveInverseSubdivision4(vtkIdType initial_face);

	int CoordinatesCoupling; // if set to 1, the three coordinates will be encoded in the same context;

	vtkIdList *vlist;


	int FileType; //==0 : .iv 1 : .ply output
	int GeometryCriterion;
	double CurvatureTreshold;
	double WaveletTreshold;
	int SolveProcess;
	int ArithmeticType; //=0 si on utilise des flottants =1 pour arithmetique enti?e 
	int SubdivisionType;  //=0 si subdivision r?uli?e =1 si subdivision irr?uli?e

	int RegularSubdivisionFlag;


	int IOType; //=0 si on fait tout en m?oire, =1 pour ?rire sur le disque, 2 pour lire
	int ProcessLifting;  //0: pas de lifting; 1:lifting classique; 2:lifting rapide
	int LiftingRadius;  // Rayon du support des ondelettes
	vtkIdList *PointsIds;


	vtkSurface *MergeInput;
	vtkSurface *MergeOutput;

	vtkIntArray *CellTypes;

	contour_element *FirstElement;
	contour_element *LastElement;
	contour_element *FirstElementRegular;
	contour_element *LastElementRegular;
	void Merge1To1(vtkIdType v1,vtkIdType v2, vtkIdType v3, vtkIdType f1);
	void Merge2To1(vtkIdType v1,vtkIdType v2, vtkIdType v3, vtkIdType v4,vtkIdType f1,vtkIdType f2);
	void Merge3To1(vtkIdType v1,vtkIdType v2, vtkIdType v3,	vtkIdType v4, vtkIdType v5,vtkIdType f1,vtkIdType f2,vtkIdType f3);

	void EncodeFaceSwap(int code,int pos);
	void EncodeEdgeSwap(int code,int pos);
	int DecodeFaceSwap(int pos);
	int DecodeEdgeSwap(int pos);




	int PossibleParent(vtkIdType v1);
	void ConquerEdge(vtkIdType v1,vtkIdType v2,vtkIdType & face,vtkIdType &vertex);
	int AddFace(vtkIdType v1,vtkIdType v2,vtkIdType v3,vtkIdType& edge1,vtkIdType& edge2,vtkIdType& edge3);
	void switchcontour(vtkIdType edge);
	void switchcontour2(vtkIdType edge);

	int MatchParents(vtkIdType v1,vtkIdType p1,vtkIdType p2,vtkIdType f1,vtkIdType f2,int incidentfaces);
	int MatchParents2(vtkIdType v1,vtkIdType p1,vtkIdType p2,vtkIdType f1,vtkIdType f2,int incidentfaces);
	void MergeAndSwap1(vtkIdType v1,vtkIdType v2, vtkIdType v3,	vtkIdType v4,vtkIdType v5,vtkIdType f1,vtkIdType f2,vtkIdType f3);
	void MergeAndSwap2(vtkIdType v1,vtkIdType v2, vtkIdType v3,
	vtkIdType v4,vtkIdType v5,vtkIdType v6,vtkIdType f1,vtkIdType f2,vtkIdType f3, vtkIdType f4);
	vtkIdType GetEdgeMidPoint(vtkIdType v1,vtkIdType v2);
	int IsFaceSwitch(vtkIdType v1,vtkIdType v2,vtkIdType v3);
	vtkIdType GetEdgeSwap(vtkIdType v1,vtkIdType v2,vtkIdType v3);
	void MakeInMatrixWeighted(vtkSurface *ptr,vtkSparseMatrix *In);
	void MakeInMatrixUnweighted(vtkSurface *ptr, vtkSparseMatrix *In);
	void MakePMatrix(vtkSurface *ptr,vtkSurface *ptr2);

	void ComputeFastLiftingCoeff();

};

#endif
