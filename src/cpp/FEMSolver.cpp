#include "CFEMTypes_Global.h"
#include "FEMSolver.h"
#include "PhyElement.h"
#include "LAFuncs.h"
#include "PhyGlobal.h"
#include "CFEMTypes_Global.h"
//#include <string>

// in C++ do not write friend again (similar to virtual)
void FEMSolver::Input(istream& in)
{
	string buf;
    in >> buf >> dim;
    //cout << buf << '\n';
    //cout << "dim\t" << dim << '\n';
    //cout << in.good() << '\n';
	in >> buf >> ndofpn;
    //cout << buf << '\n';
    //cout << "ndofpn\t" << ndofpn << '\n';
    
	in >> buf >> buf >> nNodes;
    //cout << buf << '\n';
    //cout << "nNodes\t" << ndofpn << '\n';

	in >> buf >> buf;
	nodes.resize(nNodes);
	for (int i = 0; i < nNodes; ++i)
		nodes[i].set_nndof(ndofpn);

	int tmpi;
	for (int i = 0; i < nNodes; ++i)
	{
		in >> tmpi;
		if (tmpi != (i + 1))
		{
			THROW("incorrect id");
		}
		// use tmpi - 1 = i for id instead to simplify accessing nodes	
		// for a more robust implementation that you can have any node ids, read with arbitrary order: node 10 ... node 1000 ... node 1, you can store nodes in map instead of vector node[1000] node[10] node[1]
		nodes[i].id = i;		// id is a private member of PhyNode cannot access it. Either can have set function in PhyNode::setID(int idIN) (good practice) or for quick fix here we declare this function a friend of PhyNode (refer to PhyNode.h class PhyNode)
		nodes[i].coordinate.resize(dim);
		for (int j = 0; j < dim; ++j)
			in >> nodes[i].coordinate(j);
	}
	in >> buf >> buf >> ne;
	in >> buf >> buf >> buf >> buf >> buf;
	pes.resize(ne);
	ElementType eType;
	int matID;
	//int nNodeInElement;
	PhyElement *pe;

	for (int i = 0; i < ne; ++i)
	{
		in >> tmpi;
		if (tmpi != (i + 1))
			THROW("incorrect id");

		in >> tmpi;
		eType = (ElementType)tmpi;

		pes[i] = PhyElementFactory(eType);
		pes[i]->id = tmpi;

		// 1. we don't need to use pe instead of pes[i]. It just makes read and write simpler
		// 2. we can use pe instead pes[i] because it's a pointer (it's an address).
		// 3. OBVIOUSLY we cannot do this trick with nonpointer data

		pe = pes[i]; 
		in >> pe->matID;		// longer way which was fine: in >> pes[i]->matID;
								// another way not recommended (*pe).matID
								// ptr			*ptr	object
								// object		&object	address of the object
		
		int nNodeInElement;
		in >> nNodeInElement;
		vector<int> eNodesTmp(nNodeInElement);
		vector <PhyNode*> eNodePtrsTmp(nNodeInElement);
		for (int j = 0; j < nNodeInElement; ++j)
		{
		 in >> eNodesTmp[j];
		 --eNodesTmp[j];
		 eNodePtrsTmp[j] = &nodes[eNodesTmp[j]]; // safe here because nodes size never is going to change. If not this causes a very nasty bug to fix ... 
		}
		pe->setNodeConnectivity_Sizes(nNodeInElement, ndofpn, eNodesTmp, eNodePtrsTmp); 
	}
	in >> buf >> buf >> np;
	in >> buf >> buf >> buf ; //node	node_dof_index		value

	int nodeid, dofid;
	double value;
	PhyDof* dofPtr;
	for (int i = 0; i < np; ++i)
	{
		in >> nodeid >> dofid >> value;
		--nodeid;
		--dofid;
		// could have done the last three as the following -- and ++ after the parameter does the operator after the first operation (read here) is done
//		in >> nodeid-- >> dofid-- >> value;
		
		// good practice to with a shorter pointer rather than the full name
		dofPtr = &nodes[nodeid].ndof[dofid]; // & to get the point
		dofPtr->p = true;
		dofPtr->v = value;
		// need to 
		// A. assign these values to nodes (Step 4 in course notes)
	}


	int nnzdof;  // num of nonzero force free Dofs;
	in >> buf >> buf >> nnzdof;

	in >> buf >> buf >> buf ; //node	node_dof_index		value
	for (int i = 0; i < nnzdof; ++i)
	{
		in >> nodeid >> dofid >> value;
		--nodeid;
		--dofid;

		// good practice to with a shorter pointer rather than the full name
		dofPtr = &nodes[nodeid].ndof[dofid]; // & to get the point
//		dofPtr->p = false;	// no need for this (default is false)
		dofPtr->f = value;	// force is given
	}

	in >> buf >> buf >> nmats;
	in >> buf >> buf >> buf ; //id	numPara		Paras

	int numParas, matid;
	for (int i = 0; i < nmats; ++i)
	{
		in >> matid >> numParas;
//		--matid;
//		if (matid != i)
//			THROW("wrong material id\n");

		mats[matid].setSize(numParas);
		for (int j = 0; j < numParas; ++j)
			in >> mats[matid].paras(j);
	}

	for (int e = 0; e < ne; ++e)
	{
		pe = pes[e];
		pe->setGeometry();
		matID = pe->matID;
		pe->setInternalMaterialProperties(&mats[matID]);
	}
//	return in;
}

istream& operator>>(istream& input, FEMSolver& dat)
{
	dat.Input(input);
	return input;
}

ostream& operator<<(ostream& out, const FEMSolver& dat)
{
	out << "Nodes\n";
	out << "nNodes\t" << dat.nNodes << '\n';
	out << "id\tcrd\n";
	out << "values\nforces\n";
	if (verbose == true)
	{
		out << "position(verbose)\n";
		out << "prescribed_boolean(verbose)\n";
	}
	for (int node = 0; node < dat.nNodes; ++node)
		out << dat.nodes[node] << '\n';

	out << "Elements\n";
	out << "ne\t" << dat.ne << "\n";
	out << "id ElementType\n";
	out << "forces(verbose)\n";
	out << "specific output\n";

	for (int e = 0; e < dat.ne; ++e)
		out << (*dat.pes[e]) << '\n';
	return out;
}

FEMSolver::FEMSolver(int dimIn)
{
	dim = dimIn;
}

FEMSolver::~FEMSolver()
{
	// ne should be equal to pes.size()
	// still a better practice is
	for (int i = 0; i < pes.size(); ++i)
		delete pes[i];
}


void FEMSolver::FEMSolve(string& runName, bool verboseIn)
{
	verbose = verboseIn;
	string inputFileName;
	inputFileName = runName + ".txt";
	fstream in(inputFileName.c_str(), ios::in);
	// reading data
	Input(in);
	// can do it as
//	in >> (*this);
	in.close();
	
	/////////////////////////////////////////////////////////////////////////
	// steps

	// Step 3
	setSizes();
	// Step 4; set prescribed dofs: already done when reading the input file
	// Step 5: Set global free nodal dof: already done when reading the input file 
	// Step 6 and Step 7: dof positions; Step 7: Set F
	setPositions_F();
	// Step 8: Element dof maps Me
	// Step 9: Set element dofs ae
	setElementDofMap_ae();
	// Step 10: Compute element stifness
	Calculate_ElementStiffness_Force();

	// Step 11: Assembly from local to global system
	Assemble();
	// Step 12: Solve global (free) dof a from Ka = F
	// successful solution returns true
	if (Solve_Dofs() == false)
		THROW("Matrix solve failed\n");
	// Step 13: Assign a to nodes and elements
	Assign_dof();
	// Step 14: Compute prescribed dof forces
	UpdateFpNodalPrescribedForces();

	/////////////////////////////////////////////////////////////////////////
	// output
	string outputFileName;
	outputFileName = runName + "Output.txt";
	fstream out(outputFileName.c_str(), ios::out);
	out << (*this);
}

void FEMSolver::setSizes()
{
	ndof = nNodes * ndofpn;
	//	np already read in
	nf = ndof - np;
	K.resize(nf, nf);
	F.resize(nf);
	Fp.resize(np);
}

void FEMSolver::setPositions_F()
{
	F = 0.0;
	int cntrP = 0, cntrF = 0;	// counters for prescribed and free dofs
	PhyDof* dofPtr;
	for (int node = 0; node < nNodes; ++node)
		for (int j = 0; j < nodes[node].nndof; ++j)
		{
			dofPtr = &nodes[node].ndof[j];
			if (dofPtr->p == false) // free
			{
				dofPtr->pos = cntrF;
				// or alternatively do it all in one line:
				// dofPtr->pos = cntrF++;
				F(cntrF) = dofPtr->f;
				++cntrF;
			}
			else
				dofPtr->pos = cntrP++;
		}
	if (cntrP != np) 
		THROW("inconsistency in np (cntrP != np)\n");
	if (cntrF != nf) 
		THROW("inconsistency in nf (cntrF != nf)\n");
}

void FEMSolver::setElementDofMap_ae()
{
	PhyElement* pe;
	for (int e = 0; e < ne; ++e)
		pes[e]->setElementDofMap_ae(ndofpn);
}


void FEMSolver::Calculate_ElementStiffness_Force()
{
	PhyElement* pe;
	for (int e = 0; e < ne; ++e)
		pes[e]->Calculate_ElementStiffness_Force();
}

void FEMSolver::Assemble()
{
	PhyElement* pe;
	for (int e = 0; e < ne; ++e)
		pes[e]->AssembleStiffnessForce(K, F);
}

bool FEMSolver::Solve_Dofs()
{
	dofs.resize(nf);
	dofs = F;
	cout << "K\n" << K << endl;
	int isNonsingular = LUsolve(K, dofs);
	return (isNonsingular != 0);
}

void FEMSolver::Assign_dof()
{
	PhyNode* nPtr;
	PhyDof* dofPtr;
	for (int node = 0; node < nNodes; ++node)
	{
		nPtr = &nodes[node];
		for (int j = 0; j < nPtr->nndof; ++j)
		{
			dofPtr = &nPtr->ndof[j];
			if (dofPtr->p == false)
				dofPtr->v = dofs(dofPtr->pos);
		}
	}
	PhyElement* pe;
	int I;
	for (int e = 0; e < ne; ++e)
	{
		pe = pes[e];
		for (int j = 0; j < pe->nedof; ++j)
		{
			I = pe->dofMap[j];
			if (I >= 0)
				pe->edofs(j) = dofs(I);
		}
	}
}

void FEMSolver::UpdateFpNodalPrescribedForces()
{
	PhyElement* pe;
	for (int e = 0; e < ne; ++e)
	{
		pe = pes[e];
		pe->UpdateElementForces_GlobalFp(Fp);
	}
	for (int node = 0; node < nNodes; ++node)
		nodes[node].UpdateNodePrescribedDofForces(Fp);
}
