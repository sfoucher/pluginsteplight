#ifndef ST_STL_STEP_HELLOWORLD_H
#define ST_STL_STEP_HELLOWORLD_H

#include "ct_itemdrawable/abstract/ct_abstractitemdrawablewithpointcloud.h"
#include "ct_step/abstract/ct_abstractstep.h"

class ST_STL_Step_HelloWorld : public CT_AbstractStep
{
    Q_OBJECT
    using SuperClass = CT_AbstractStep;

public:
    ST_STL_Step_HelloWorld();
    QString description() const override;
    QString getStepDetailledDescription() const;
    CT_VirtualAbstractStep* createNewInstance() const override;

protected:
    void declareInputModels(CT_StepInModelStructureManager& manager) final;
    void declareOutputModels(CT_StepOutModelStructureManager& manager) final;
    void fillPostInputConfigurationDialog(CT_StepConfigurableDialog* postInputConfigDialog) final;
    void compute() final;

protected:
    CT_HandleInResultGroupCopy<>                                    _inResult;
    CT_HandleInStdZeroOrMoreGroup                                   _inZeroOrMoreRootGroup;
    CT_HandleInStdGroup<>                                           _inGroup;
    CT_HandleInSingularItem<CT_AbstractItemDrawableWithPointCloud>  _in_point_cloud;

};


#endif // ST_STL_STEP_HELLOWORLD_H
