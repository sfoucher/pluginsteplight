#include "stl_grid3dbeamvisitor.h"

STL_Grid3DBeamVisitor::STL_Grid3DBeamVisitor(CT_Grid3D<int> *grid) {
    _grid = grid;
}

void STL_Grid3DBeamVisitor::visit(const size_t &index, const CT_Beam *beam)
{
    Q_UNUSED(beam);

    _grid->addValueAtIndex(index, 1);
}
