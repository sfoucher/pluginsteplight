#ifndef STL_GRID3DBEAMVISITOR_H
#define STL_GRID3DBEAMVISITOR_H

#include "ct_itemdrawable/ct_grid3d.h"
#include "ct_itemdrawable/tools/gridtools/ct_abstractgrid3dbeamvisitor.h"

class STL_Grid3DBeamVisitor : public CT_AbstractGrid3DBeamVisitor
{
public:
    STL_Grid3DBeamVisitor(CT_Grid3D<int> *grid);
    virtual void visit(const size_t &index, const CT_Beam *beam);

private:
    CT_Grid3D<int>* _grid;
};

#endif // STL_GRID3DBEAMVISITOR_H
