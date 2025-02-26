#ifndef ST_STL_STEPPLUGINMANAGER_H
#define ST_STL_STEPPLUGINMANAGER_H

#include "ct_abstractstepplugin.h"

class ST_STL_StepPluginManager : public CT_AbstractStepPlugin
{
public:
    ST_STL_StepPluginManager();
    ~ST_STL_StepPluginManager() override;

    QString getPluginURL() const override {return QString("No link available yet!");}

    virtual QString getPluginOfficialName() const override {return "Stepv6_ravaglia";}

    QString getPluginRISCitation() const override;
protected:

    bool loadGenericsStep() override;
    bool loadOpenFileStep() override;
    bool loadCanBeAddedFirstStep() override;
    bool loadActions() override;
    bool loadExporters() override;
    bool loadReaders() override;
};

#endif // ST_STL_STEPPLUGINMANAGER_H
