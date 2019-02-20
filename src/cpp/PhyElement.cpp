#include "PhyElement.h"
#include "PhyElementBar.h"
#include "PhyElementTruss.h"
#include "PhyNode.h"
#include "PhyDof.h"
#include "CFEMTypes_Global.h"
#include "PhyGlobal.h"

PhyElement* PhyElementFactory(ElementType eTypeIn)
{
	PhyElement* pePtr = NULL;
	switch (eTypeIn)
	{
//	case etDefault:
//		pePtr = new PhyElement();
//		break;
	case etBar:
		pePtr = new PhyElementBar();
		break;
	case etTruss:
		pePtr = new PhyElementTruss();
		break;
	default:
		THROW("the type is not defined");
	}
	if (pePtr != NULL)
		pePtr->eType = eTypeIn;
	return pePtr;
}

ostream& operator<<(ostream& out, const PhyElement& dat)
{
	// id ElementType
	out << dat.id << '\t' << dat.eType << '\n';
	if (verbose == true) {
		for (int i = 0; i < dat.nedof; ++i)
			out << dat.fee(i) << '\t';
	}
	out << '\n';
	dat.SpecificOutput(out);
	out << '\n';
	return out;
}

void PhyElement::setNodeConnectivity_Sizes(int nNodeInElement, int ndofpnIn, vector<int>& eNodesIn, vector <PhyNode*>& eNodePtrsIn)
{
	neNodes = nNodeInElement;
	eNodes.resize(neNodes);
	eNodes = eNodesIn;

	eNodePtrs.resize(neNodes);
	eNodePtrs = eNodePtrsIn;
	// resizing members in PhyElement

	nedof = neNodes * ndofpnIn;
	edofs.resize(nedof); 
	dofMap.resize(nedof); 
	ke.resize(nedof, nedof);
	foe.resize(nedof);	// may be better not to do it here, because the element may not have this force
	foe = 0.0;

	fde.resize(nedof);	// may be better not to do it here, because the element may not have this force
	fde = 0.0;

	fee.resize(nedof);	// may be better not to do it here, because the element may not have this force
	fee = 0.0;
}


void PhyElement::setElementDofMap_ae(int ndofpn)
{
	PhyNode* nodePtr;
	PhyDof* dofPtr;
	int cntr = 0;
	for (int node = 0; node < neNodes; ++ node)
	{
		nodePtr = eNodePtrs[node];
		for (int dof = 0; dof < ndofpn; ++dof)
		{
			dofPtr = &nodePtr->ndof[dof];
			if (dofPtr->p == true)  // subtract prescibed position by -1 so that it's distinguished from free ones (0 = -0)
			{
				edofs(cntr) = dofPtr->v;
				dofMap[cntr++] = - (dofPtr->pos + 1);
			}
			else
				dofMap[cntr++] = dofPtr->pos;
		}
	}
}


void PhyElement::AssembleStiffnessForce(MATRIX& globalK, VECTOR& globalF)
{
	fee.resize(nedof);
	if (foe.size() == nedof)
		fee = foe;
	else
		fee = 0.0;

	int I, J;
	for (int i = 0; i < nedof; ++i)
	{
		I = dofMap[i];
		if (I < 0) // prescribed dof
			continue;
		for (int j = 0; j < nedof; ++j)
		{
			J = dofMap[j];
			if (J < 0) // prescibed
				fee(i) -= ke(i, j) * edofs(j);
			else
				globalK(I, J) += ke(i, j);
		}
		globalF(I) += fee(i);
	}
}

void PhyElement::UpdateElementForces_GlobalFp(VECTOR& Fp)
{
	if (fde.size() == 0)
	{
		fde.resize(nedof);
	}
	fde = 0.0;
	if (foe.size() == 0)
	{
		foe.resize(nedof);
		foe = 0.0;
	}
	if (fee.size() == 0)
	{
		fee.resize(nedof);
	}
	fee = 0.0;

	int I;
	for (int i = 0; i < nedof; ++i)
	{
		I = dofMap[i];
		for (int j = 0; j < nedof; ++j)
				fde(i) += ke(i, j) * edofs(j);
		fee(i) = foe(i) - fde(i);
		if (I < 0) // prescribed dof
		{
			I = -1 - I;
			Fp(I) -= fee(i);
		}
	}
}
