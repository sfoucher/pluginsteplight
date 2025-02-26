#ifndef ST_STEPEXTRACTCURVESFROM3DGRID_H
#define ST_STEPEXTRACTCURVESFROM3DGRID_H

#include "ct_itemdrawable/ct_circle.h"
#include "ct_step/abstract/ct_abstractstep.h"

class ST_StepExtractCurvesFrom3DGrid : public CT_AbstractStep
{
    Q_OBJECT
    using SuperClass = CT_AbstractStep;

    using Vec3d = Eigen::Vector3d;
    using Vec4d = Eigen::Vector4d;
    using Vec4i = Eigen::Vector4i;

    using Circle                    = CT_Circle;
    using CirclePtr                 = Circle*;

    //using Snake     = ST_OpenActiveContours<HoughSpaceValueType>;
    //using SnakePtr  = Snake*;
public:
    ST_StepExtractCurvesFrom3DGrid();
};

#endif // ST_STEPEXTRACTCURVESFROM3DGRID_H
