#include "st_stl_step_helloworld.h"

#include <ct_log/ct_logmanager.h>

ST_STL_Step_HelloWorld::ST_STL_Step_HelloWorld() : SuperClass()
{

}

QString ST_STL_Step_HelloWorld::description() const
{
    return tr("1 - Hello World en console");
}

// Step detailled description
QString ST_STL_Step_HelloWorld::getStepDetailledDescription() const
{
    return tr("Écrire Hello World en console");
}

CT_VirtualAbstractStep* ST_STL_Step_HelloWorld::createNewInstance() const
{
    return new ST_STL_Step_HelloWorld();
}

//////////////////// PROTECTED METHODS //////////////////

void ST_STL_Step_HelloWorld::declareInputModels(CT_StepInModelStructureManager& manager)
{
    manager.addResult(_inResult, tr("Scène(s)"));
    manager.setZeroOrMoreRootGroup(_inResult, _inZeroOrMoreRootGroup);
    manager.addGroup(_inZeroOrMoreRootGroup, _inGroup);
    manager.addItem(_inGroup, _in_point_cloud, tr("Point cloud"));
}

void ST_STL_Step_HelloWorld::declareOutputModels(CT_StepOutModelStructureManager& manager)
{
    manager.addResultCopy(_inResult);
}

void ST_STL_Step_HelloWorld::fillPostInputConfigurationDialog(CT_StepConfigurableDialog* postInputConfigDialog)
{

}

void ST_STL_Step_HelloWorld::compute()
{
    using PointCloudConst       = const CT_AbstractItemDrawableWithPointCloud;
    using PointCloudConstPtr    = PointCloudConst*;

    setProgress(0);

    for (CT_StandardItemGroup* group : _inGroup.iterateOutputs(_inResult))
    {
        if( isStopped() )
        {
            return;
        }

        // -----------------------------------------------------------------------------------------------------------------
        // Get point cloud and normal cloud from computree previous step
        PointCloudConstPtr  inPointCloud = group->singularItem(_in_point_cloud);
        // -----------------------------------------------------------------------------------------------------------------

        PS_LOG->addInfoMessage(LogInterface::error, tr("Hello"));

        // -----------------------------------------------------------------------------------------------------------------

    }

    setProgress(100);
}
